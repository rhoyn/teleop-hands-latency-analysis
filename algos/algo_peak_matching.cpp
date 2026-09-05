#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

#include "algo.hpp"

namespace {

/**
 * Envelope shaping, burst segmentation and peak pairing limits.
 *
 * Two envelopes are built from the same signal for different jobs, which is
 * why there are two smoothing constants: the sharp one at `SEG_SMOOTH_S`
 * decides where a reach starts and stops, while the blunt one at
 * `ENV_SMOOTH_S` has to leave exactly one peak per reach for the pairing to
 * be meaningful, and under-smoothing it is the fastest way to make this
 * method report nonsense. The pairing window is one-sided because the robot
 * cannot physically lead; the small backward margin only exists so a pair
 * pushed against zero is visible as a rail.
 */
constexpr double BAND_FC_HZ = 15.0;      ///< 15 Hz, one peak per reach
constexpr double ENV_SMOOTH_S = 0.060;   ///< 60 ms, blunt enough to be single
constexpr double SEG_SMOOTH_S = 0.020;   ///< 20 ms, keeps burst edges sharp
constexpr double MIN_BURST_S = 0.12;     ///< shorter is a twitch, not a reach
constexpr double MIN_REST_S = 0.12;      ///< closer than this is one reach
constexpr double SEARCH_BACK_S = 0.03;   ///< small, the robot cannot lead
constexpr double SEARCH_AHEAD_S = 0.30;  ///< 300 ms past the burst end
constexpr double TAU_MIN_MS = -20.0;     ///< negative only to expose a rail
constexpr double TAU_MAX_MS = 250.0;     ///< beyond it the pairing is wrong

/**
 * A half-open sample range covering one detected reach.
 *
 * Indices are into the full-rate signals; `end` is one past the last
 * sample.
 */
struct Burst {
  size_t begin;
  size_t end;
};

/**
 * One reach's matched pair of speed peaks, on the operator and the robot.
 *
 * The times carry the parabolic sub-sample refinement and so fall between
 * samples, while the indices are the raw integer argmaxes kept for
 * plotting, which is why both forms are stored rather than derived.
 */
struct Pair {
  double t_op;
  double t_rob;
  size_t i_op;
  size_t i_rob;
};

/**
 * Per-hand peak envelope values used to scale a panel to a common height.
 *
 * The robot attenuates the operator, so the two traces are drawn each
 * against its own maximum; the panel is about when the peaks happen, not
 * how tall they are.
 */
struct Norm {
  double rob;
  double op;
};

/**
 * Find the first sample at or after a wall-clock time.
 *
 * A linear scan, which is fine because it runs a handful of times per tile
 * on an already ordered axis, and it clamps to the last sample rather than
 * running off the end.
 *
 * @param[in] t  ascending timestamps in seconds
 * @param[in] v  time to look up
 * @returns the index of the first sample not before `v`, clamped
 * @exceptsafe basic
 */
size_t index_at(
    const std::vector<double>& t,
    double v
) {
  size_t i = 0;
  while (i + 1 < t.size() && t[i] < v) i++;
  return i;
}

/**
 * Draw both speed envelopes over a sample range, each on its own scale.
 *
 * Returns the normalisers so the caller can place markers and brackets at
 * the same relative heights as the curves; recomputing them at the call
 * site would risk the two disagreeing.
 *
 * @param[in,out] im     tile being drawn into
 * @param[in,out] ax     panel axes, range set here
 * @param[in]     in     timestamps, used for the elapsed-time axis
 * @param[in]     rob    robot speed envelope
 * @param[in]     op     operator speed envelope
 * @param[in]     ia     first sample of the range to draw
 * @param[in]     ib     one past the last sample of the range
 * @param[in]     y_top  top of the y range, headroom for the brackets
 * @returns the per-hand maxima the traces were divided by
 * @exceptsafe basic
 */
Norm draw_envelopes(
    canvas::Image& im,
    canvas::Axes& ax,
    const algo::Input& in,
    const std::vector<double>& rob,
    const std::vector<double>& op,
    size_t ia,
    size_t ib,
    double y_top
) {
  Norm nz{1e-12, 1e-12};
  for (size_t i = ia; i < ib; i++) {
    nz.rob = std::max(nz.rob, rob[i]);
    nz.op = std::max(nz.op, op[i]);
  }
  const size_t stride = std::max<size_t>(1, (ib - ia) / 1400);
  std::vector<double> xs, ys_rob, ys_op;
  for (size_t i = ia; i < ib; i += stride) {
    xs.push_back(in.t[i] - in.t[ia]);
    ys_rob.push_back(rob[i] / nz.rob);
    ys_op.push_back(op[i] / nz.op);
  }
  canvas::set_range(ax, 0.0, in.t[ib - 1] - in.t[ia], 0.0, y_top);
  canvas::grid_and_ticks(im, ax, 10, 5);
  canvas::polyline(im, ax, xs, ys_rob, canvas::BLACK);
  canvas::polyline(im, ax, xs, ys_op, canvas::VIVID);
  return nz;
}

/**
 * Write the measured gap in ms above a bracket, on its own plate.
 *
 * The filled plate is drawn first so the label stays readable where it
 * lands on top of a trace, and the text is centred on the bracket rather
 * than on either peak.
 *
 * @param[in,out] im    tile being drawn into
 * @param[in]     ax    panel axes the values are mapped through
 * @param[in]     xa_v  bracket start in data units
 * @param[in]     xb_v  bracket end in data units
 * @param[in]     y_v   bracket height in data units
 * @param[in]     ms    gap to print, in ms
 * @exceptsafe basic
 */
void gap_label(
    canvas::Image& im,
    const canvas::Axes& ax,
    double xa_v,
    double xb_v,
    double y_v,
    double ms
) {
  char lab[32];
  snprintf(lab, sizeof lab, "%.0f ms", ms);
  const int w = canvas::text_width(lab, canvas::LABEL_SCALE);
  const int h = canvas::text_height(canvas::LABEL_SCALE);
  const int mid = (canvas::sx(ax, xa_v) + canvas::sx(ax, xb_v)) / 2;
  const int x = mid - w / 2;
  const int y = canvas::sy(ax, y_v) - h - 12;
  canvas::fill_rect(im, x - 10, y - 6, w + 20, h + 14, canvas::PANEL);
  canvas::rect(im, x - 10, y - 6, w + 20, h + 14, canvas::WHITE);
  canvas::text(im, x, y, lab, canvas::WHITE, canvas::LABEL_SCALE);
}

/**
 * Draw a flat bracket with end ticks between two times.
 *
 * This is the visual the whole method is meant to sell, so it is drawn
 * three pixels thick with ticks that overshoot the bar at both ends and
 * stays legible after the tile is scaled down.
 *
 * @param[in,out] im    tile being drawn into
 * @param[in]     ax    panel axes the values are mapped through
 * @param[in]     xa_v  bracket start in data units
 * @param[in]     xb_v  bracket end in data units
 * @param[in]     y_v   bracket height in data units
 * @param[in]     c     bracket colour
 * @exceptsafe basic
 */
void bracket(
    canvas::Image& im,
    const canvas::Axes& ax,
    double xa_v,
    double xb_v,
    double y_v,
    canvas::RGB c
) {
  const int xa = canvas::sx(ax, xa_v);
  const int xb = canvas::sx(ax, xb_v);
  const int yb = canvas::sy(ax, y_v);
  for (int k = 0; k < 3; k++) canvas::line(im, xa, yb + k, xb, yb + k, c);
  for (int k = -6; k <= 8; k++) {
    canvas::put(im, xa, yb + k, c);
    canvas::put(im, xb, yb + k, c);
  }
}

/**
 * Finish an abstaining result, drawing the envelopes when there are any.
 *
 * Abstaining is a legitimate outcome and not an error: the reported latency
 * is the median over ten algorithms, so a method that could not pair two
 * peaks should stand aside instead of guessing. `Result::note` is left as
 * the caller set it and names which check failed.
 *
 * @param[in] r    partly filled result, taken by value and returned
 * @param[in] in   timestamps, or null to draw an empty tile
 * @param[in] rob  robot speed envelope, or null
 * @param[in] op   operator speed envelope, or null
 * @returns the same result with `ok` false and an explanatory tile
 * @exceptsafe basic
 */
algo::Result rejected(
    algo::Result r,
    const algo::Input* in = nullptr,
    const std::vector<double>* rob = nullptr,
    const std::vector<double>* op = nullptr
) {
  r.image = canvas::standard("PEAK MATCHING", r.tau_ms, false);
  if (in && rob && op && rob->size() > 32 && op->size() == rob->size()) {
    canvas::panel_label(
        r.image,
        0,
        "1  HAND SPEED OVER TIME, NO PEAKS COULD BE MATCHED"
    );
    canvas::Axes a0 = canvas::panel(0);
    draw_envelopes(r.image, a0, *in, *rob, *op, 0, rob->size(), 1.15);
    canvas::axis_labels(r.image, a0, "TIME (s)", "HAND SPEED (relative)");
    canvas::legend(
        r.image,
        a0,
        {{"OPERATOR", canvas::VIVID}, {"ROBOT", canvas::BLACK}},
        true
    );
  }
  return r;
}

/**
 * Moving average of half-width `half` samples, computed by prefix sums.
 *
 * The window shrinks at both ends rather than wrapping or zero-padding, so
 * a reach that starts near the beginning of a take keeps its true height.
 * Cost is linear in the signal whatever the width.
 *
 * @param[in] x     signal to smooth
 * @param[in] half  half-width in samples, so the window is 2*half + 1
 * @returns the smoothed signal, the same length as `x`
 * @exceptsafe basic
 */
std::vector<double> box_smooth(
    const std::vector<double>& x,
    size_t half
) {
  const size_t n = x.size();
  std::vector<double> p(n + 1, 0.0), out(n, 0.0);
  for (size_t i = 0; i < n; i++) p[i + 1] = p[i] + x[i];
  for (size_t i = 0; i < n; i++) {
    const size_t a = i > half ? i - half : 0;
    const size_t b = std::min(n, i + half + 1);
    out[i] = (p[b] - p[a]) / double(b - a);
  }
  return out;
}

/**
 * Low-pass in place with a one-pole filter, forwards only.
 *
 * Applied twice by the caller for a steeper roll-off. It is deliberately
 * not run backwards to make it zero phase: the resulting group delay is
 * identical on both hands and cancels out of the difference of two peak
 * times, whereas a filtfilt pass would smear each peak symmetrically and
 * blunt the very feature being timed.
 *
 * @param[in,out] x   signal filtered in place
 * @param[in]     fc  cutoff in Hz
 * @param[in]     dt  sample interval in seconds
 * @exceptsafe basic
 */
void onepole(
    std::vector<double>& x,
    double fc,
    double dt
) {
  if (x.empty()) return;
  const double a = std::exp(-2 * M_PI * fc * dt);
  double s = x[0];
  for (size_t i = 0; i < x.size(); i++) {
    s = a * s + (1 - a) * x[i];
    x[i] = s;
  }
}

/**
 * Order statistic of a copy of the samples, by partial sort.
 *
 * Takes its argument by value because `nth_element` reorders it and the
 * caller still needs the envelope. No interpolation between neighbours: it
 * only sets burst thresholds, where one sample of slop is far below the
 * noise on the envelope itself.
 *
 * @param[in] v  samples, in any order
 * @param[in] q  quantile in [0, 1]
 * @returns the sample at that quantile, or zero when `v` is empty
 * @exceptsafe basic
 */
double percentile(
    std::vector<double> v,
    double q
) {
  if (v.empty()) return 0;
  const size_t k = size_t(q * double(v.size() - 1));
  std::nth_element(v.begin(), v.begin() + k, v.end());
  return v[k];
}

/**
 * Build a smoothed speed envelope of one hand, optionally band limited.
 *
 * High-passing at 0.3 Hz first removes the gyro's slow bias drift, which
 * would otherwise survive rectification as a pedestal and lift the burst
 * thresholds. The optional low pass runs before rectifying so it shapes the
 * rate rather than the envelope, and passing `fc` of zero skips it, which
 * is how the same routine produces both the blunt envelope used for timing
 * and the sharp one used for segmentation.
 *
 * @param[in] in        both channels plus the sample step
 * @param[in] mpu       `algo::ROBOT` or `algo::OPERATOR`
 * @param[in] fc        low-pass cutoff in Hz, or zero to skip it
 * @param[in] smooth_s  averaging half-width in seconds
 * @returns the envelope at the input sample rate
 * @exceptsafe basic
 */
std::vector<double> speed_envelope(
    const algo::Input& in,
    int mpu,
    double fc,
    double smooth_s
) {
  std::vector<double> x = algo::motion_signal(in, mpu, 0.3);
  if (fc > 0) {
    onepole(x, fc, in.dt);
    onepole(x, fc, in.dt);
  }
  for (double& v : x) v = std::fabs(v);
  return box_smooth(x, size_t(smooth_s / in.dt));
}

/**
 * Segment the envelope into reaches with a hysteretic threshold.
 *
 * Thresholds are relative to the take's own 20th and 99.5th percentiles
 * rather than absolute in dps, because gain depends on how vigorously the
 * operator moved; taking the 99.5th percentile instead of the maximum stops
 * one knock raising the bar for every other reach. A run must cross the
 * high threshold to count but is followed out to the low one so the slow
 * tail stays attached, runs closer than `MIN_REST_S` are merged, and
 * anything shorter than `MIN_BURST_S` is then dropped.
 *
 * @param[in] env  sharp speed envelope of the operator hand
 * @param[in] dt   sample interval in seconds
 * @returns the accepted bursts in time order, possibly empty
 * @exceptsafe basic
 */
std::vector<Burst> find_bursts(
    const std::vector<double>& env,
    double dt
) {
  const size_t n = env.size();
  const double base = percentile(env, 0.20);
  const double span = percentile(env, 0.995) - base;
  std::vector<Burst> out;
  if (span <= 0) return out;
  const double hi = base + 0.35 * span;
  const double lo = base + 0.08 * span;
  const size_t min_rest = size_t(MIN_REST_S / dt);
  const size_t min_burst = size_t(MIN_BURST_S / dt);
  size_t i = 0;
  while (i < n) {
    if (env[i] <= lo) {
      i++;
      continue;
    }
    size_t j = i;
    double peak = 0;
    while (j < n && env[j] > lo) {
      peak = std::max(peak, env[j]);
      j++;
    }
    if (peak > hi) {
      if (!out.empty() && i - out.back().end < min_rest)
        out.back().end = j;
      else
        out.push_back({i, j});
    }
    i = j;
  }
  std::vector<Burst> keep;
  for (const Burst& b : out)
    if (b.end - b.begin >= min_burst) keep.push_back(b);
  return keep;
}

/**
 * Index of the largest sample in a half-open range.
 *
 * Ties go to the earliest sample, which matters on the robot's blunted
 * envelope where a peak can sit almost flat across several milliseconds.
 *
 * @param[in] v     samples to search
 * @param[in] from  first index of the range
 * @param[in] to    one past the last index of the range
 * @returns the index of the maximum, or `from` for an empty range
 * @exceptsafe no-throw
 */
size_t argmax_in(
    const std::vector<double>& v,
    size_t from,
    size_t to
) {
  size_t best = from;
  for (size_t i = from; i < to; i++)
    if (v[i] > v[best]) best = i;
  return best;
}

/**
 * Sub-sample offset of a peak, from a parabola through its neighbours.
 *
 * At 1 kHz a whole-sample argmax quantises the delay to 1 ms, which is a
 * visible share of the 17-21 ms seen on the fastest rigs, so the peak is
 * refined against the two samples either side. Offsets outside one sample
 * mean the vertex is not really here, usually a flat or noisy top, and are
 * dropped to zero rather than trusted.
 *
 * @param[in] v  samples containing the peak
 * @param[in] i  index of the integer argmax
 * @returns the offset in samples, within [-1, 1], zero if not resolvable
 * @exceptsafe no-throw
 */
double parabolic_offset(
    const std::vector<double>& v,
    size_t i
) {
  if (i == 0 || i + 1 >= v.size()) return 0.0;
  const double den = v[i - 1] - 2 * v[i] + v[i + 1];
  if (std::fabs(den) < 1e-12) return 0.0;
  const double s = 0.5 * (v[i - 1] - v[i + 1]) / den;
  return (s > 1.0 || s < -1.0) ? 0.0 : s;
}

/**
 * Time each reach by the gap between the two hands' speed peaks.
 *
 * Reaches are segmented on the operator hand, the fastest moment of each is
 * found on both hands as the argmax of the smoothed speed envelope, refined
 * below the sample grid with a parabola, and the difference of those two
 * times is one delay measurement; the reported value is the median over
 * reaches and the spread is their interquartile range. The robot's peak is
 * looked for from just before the operator's peak out to `SEARCH_AHEAD_S`
 * past the end of the burst, and a pair falling outside `TAU_MIN_MS` to
 * `TAU_MAX_MS` is discarded rather than reported. Comparing the timing of a
 * single landmark rather than whole waveforms is what makes this robust to
 * the robot attenuating and low-passing the operator, and it is by a wide
 * margin the most explainable method in the set: this is the one to put in
 * front of a non-technical audience, since panel 2 zooms a single reach and
 * brackets the gap being measured with the number written above it. Its
 * weakness is the same landmark, which needs one clean speed peak per
 * reach, so slow or compound motions blur the peak and widen the spread.
 * `Result::note` reports how many reaches were found and how many of them
 * yielded a peak pair inside the window.
 *
 * @param[in] in  both gyro magnitude channels at 1 kHz plus the sample step
 * @returns the delay in ms, positive when the robot lags the operator, with
 *          `ok` false when the take is tiny, no burst was found or fewer
 *          than two pairs survived the window
 * @exceptsafe basic
 */
algo::Result peak_matching(const algo::Input& in) {
  algo::Result r;
  const size_t n = in.t.size();
  if (n < 64) {
    r.note = "not enough samples";
    return rejected(r);
  }
  std::vector<double> env[algo::MPUS];
  for (int m = 0; m < algo::MPUS; m++)
    env[m] = speed_envelope(in, m, BAND_FC_HZ, ENV_SMOOTH_S);
  const std::vector<double> seg =
      speed_envelope(in, algo::OPERATOR, 0.0, SEG_SMOOTH_S);

  const std::vector<Burst> bursts = find_bursts(seg, in.dt);
  const int nbursts = (int)bursts.size();
  if (bursts.empty()) {
    r.note = "no motion bursts found on the operator hand";
    return rejected(r, &in, &env[algo::ROBOT], &env[algo::OPERATOR]);
  }

  const size_t back = size_t(SEARCH_BACK_S / in.dt);
  const size_t ahead = size_t(SEARCH_AHEAD_S / in.dt);
  std::vector<double> taus;
  std::vector<Pair> pairs;
  for (const Burst& b : bursts) {
    const size_t pk_op =
        argmax_in(env[algo::OPERATOR], b.begin, std::min(n, b.end));
    const size_t from = pk_op > back ? pk_op - back : 0;
    const size_t to = std::min(n, b.end + ahead);
    if (from + 2 >= to) continue;
    const size_t pk_rob = argmax_in(env[algo::ROBOT], from, to);
    const double t_op =
        in.t[pk_op] + parabolic_offset(env[algo::OPERATOR], pk_op) * in.dt;
    const double t_rob =
        in.t[pk_rob] + parabolic_offset(env[algo::ROBOT], pk_rob) * in.dt;
    const double tau = (t_rob - t_op) * 1000.0;
    if (tau < TAU_MIN_MS || tau > TAU_MAX_MS) continue;
    taus.push_back(tau);
    pairs.push_back({t_op, t_rob, pk_op, pk_rob});
  }

  char buf[128];
  if (taus.size() < 2) {
    snprintf(
        buf,
        sizeof buf,
        "only %d matched peaks from %d bursts",
        (int)taus.size(),
        nbursts
    );
    r.note = buf;
    return rejected(r, &in, &env[algo::ROBOT], &env[algo::OPERATOR]);
  }

  r.ok = true;
  r.tau_ms = algo::median(taus);
  r.spread_ms = algo::iqr(taus);
  snprintf(
      buf,
      sizeof buf,
      "%d reaches, %d peak pairs matched",
      nbursts,
      (int)taus.size()
  );
  r.note = buf;

  r.image = canvas::standard("PEAK MATCHING", r.tau_ms, r.ok);
  canvas::caption(r.image, "MATCHED SPEED PEAKS");

  canvas::panel_label(
      r.image,
      0,
      "1  EACH REACH, THE OPERATOR PEAK AND THE ROBOT PEAK"
  );
  {
    const double t_lo = in.t.front();
    const double t_hi = in.t.back();
    const size_t ia = index_at(in.t, t_lo);
    const size_t ib = std::min(n, index_at(in.t, t_hi) + 1);
    if (ib > ia + 8) {
      canvas::Axes a0 = canvas::panel(0);
      const Norm nz = draw_envelopes(
          r.image,
          a0,
          in,
          env[algo::ROBOT],
          env[algo::OPERATOR],
          ia,
          ib,
          1.32
      );
      const double base = in.t[ia];
      for (const Pair& p : pairs) {
        if (p.t_op < t_lo || p.t_rob > t_hi) continue;
        canvas::marker(
            r.image,
            a0,
            p.t_op - base,
            env[algo::OPERATOR][p.i_op] / nz.op,
            canvas::VIVID,
            4
        );
        canvas::marker(
            r.image,
            a0,
            p.t_rob - base,
            env[algo::ROBOT][p.i_rob] / nz.rob,
            canvas::PALE,
            4
        );
        bracket(
            r.image,
            a0,
            p.t_op - base,
            p.t_rob - base,
            1.16,
            canvas::WHITE
        );
      }
      canvas::axis_labels(r.image, a0, "TIME (s)", "HAND SPEED (relative)");
      canvas::legend(
          r.image,
          a0,
          {{"OPERATOR", canvas::VIVID},
           {"ROBOT", canvas::BLACK},
           {"DELAY GAP", canvas::WHITE, false, true}},
          true
      );
    }
  }

  canvas::panel_label(
      r.image,
      1,
      "2  ONE REACH, THE GAP BETWEEN PEAKS IS THE DELAY"
  );
  {
    size_t pick = 0;
    double closest = 1e300;
    for (size_t i = 0; i < taus.size(); i++) {
      const double d = std::fabs(taus[i] - r.tau_ms);
      if (d < closest) {
        closest = d;
        pick = i;
      }
    }
    const Pair p = pairs[pick];
    double t_lo = p.t_op - 0.22;
    double t_hi = p.t_rob + 0.28;
    if (t_lo < in.t.front()) t_lo = in.t.front();
    if (t_hi > in.t.back()) t_hi = in.t.back();
    const size_t ia = index_at(in.t, t_lo);
    const size_t ib = std::min(n, index_at(in.t, t_hi) + 1);
    if (ib > ia + 8) {
      canvas::Axes a1 = canvas::panel(1);
      const Norm nz = draw_envelopes(
          r.image,
          a1,
          in,
          env[algo::ROBOT],
          env[algo::OPERATOR],
          ia,
          ib,
          1.55
      );
      const double base = in.t[ia];
      canvas::mark_guide(r.image, a1, p.t_op - base, true, "PEAK TIME", 0);
      canvas::mark_guide(r.image, a1, p.t_rob - base, true, "PEAK TIME", 1);
      canvas::marker(
          r.image,
          a1,
          p.t_op - base,
          env[algo::OPERATOR][p.i_op] / nz.op,
          canvas::VIVID,
          6
      );
      canvas::marker(
          r.image,
          a1,
          p.t_rob - base,
          env[algo::ROBOT][p.i_rob] / nz.rob,
          canvas::PALE,
          6
      );
      bracket(r.image, a1, p.t_op - base, p.t_rob - base, 1.20, canvas::WHITE);
      gap_label(r.image, a1, p.t_op - base, p.t_rob - base, 1.20, taus[pick]);
      canvas::axis_labels(r.image, a1, "TIME (s)", "HAND SPEED (relative)");
      canvas::legend(
          r.image,
          a1,
          {{"OPERATOR", canvas::VIVID},
           {"ROBOT", canvas::BLACK},
           {"DELAY GAP", canvas::WHITE, false, true}},
          true
      );
    }
  }

  canvas::panel_label(
      r.image,
      2,
      "3  THE DELAY MEASURED ON EVERY REACH, AND THE MEDIAN"
  );
  {
    canvas::Axes a2 = canvas::panel(2);
    double lo = taus.front(), hi = taus.front();
    for (double v : taus) {
      lo = std::min(lo, v);
      hi = std::max(hi, v);
    }
    const double pad = std::max(5.0, 0.15 * (hi - lo));
    canvas::set_range(a2, 0.4, double(taus.size()) + 0.6, lo - pad, hi + pad);
    canvas::grid_and_ticks(r.image, a2, 10, 5);
    canvas::mark_guide(r.image, a2, r.tau_ms, false, "MEDIAN", 0);
    for (size_t i = 0; i < taus.size(); i++)
      canvas::marker(r.image, a2, double(i + 1), taus[i], canvas::SAGE, 5);
    canvas::axis_labels(r.image, a2, "REACH NUMBER", "DELAY (ms)");
    canvas::legend(r.image, a2, {{"ONE REACH", canvas::WHITE, true}}, true);
  }

  return r;
}

}

ALGO_REGISTER(
    "peak_matching",
    peak_matching
);
