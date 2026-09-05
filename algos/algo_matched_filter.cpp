#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

#include "algo.hpp"

namespace {

/**
 * Burst segmentation and lag search limits, tuned on reach-and-rest takes.
 *
 * Segmentation runs on the operator hand alone, because that is the hand
 * that starts a motion and the robot's copy is both weaker and later. The
 * lag window is deliberately lopsided: physically the robot can only lag,
 * and the small negative margin exists so a burst whose peak wants to sit
 * at zero is visibly at a rail rather than silently clipped.
 */
constexpr double SEG_SMOOTH_S = 0.020;  ///< 20 ms, keeps burst edges sharp
constexpr double MIN_BURST_S = 0.12;    ///< shorter is a twitch, not a reach
constexpr double MIN_REST_S = 0.12;     ///< closer than this is one reach
constexpr double BURST_PAD_S = 0.10;    ///< room for the robot to catch up
constexpr double LAG_MIN_MS = -20.0;    ///< negative only to expose a rail
constexpr double LAG_MAX_MS = 250.0;    ///< past the worst stack, 71-98 ms

/**
 * A half-open sample range covering one detected reach.
 *
 * Indices are into the full-rate operator signal; `end` is one past the
 * last sample, and both ends are widened by `BURST_PAD_S` only at the point
 * of use, not here.
 */
struct Burst {
  size_t begin;
  size_t end;
};

/**
 * Moving average of half-width `half` samples, computed by prefix sums.
 *
 * The window shrinks at both ends rather than wrapping or zero-padding, so
 * the envelope stays on scale where a burst begins near the start of a
 * recording. Cost is linear in the signal regardless of the width.
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
 * Order statistic of a copy of the samples, by partial sort.
 *
 * Takes its argument by value because `nth_element` reorders it, and the
 * caller still needs the envelope intact. No interpolation between
 * neighbours: it only sets burst thresholds, where a sample of slop is far
 * below the noise on the envelope itself.
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
 * Build a smoothed magnitude envelope of one hand's rotation rate.
 *
 * High-passing at 0.3 Hz first removes the gyro's slow bias drift, which
 * would otherwise survive the rectification as a constant pedestal and lift
 * the burst thresholds. Rectify-then-smooth gives a rough speed envelope
 * whose bumps line up with reaches; it is used only to find bursts, never
 * to measure the delay, which is done on the signed signal.
 *
 * @param[in] in        both channels plus the sample step
 * @param[in] mpu       `algo::ROBOT` or `algo::OPERATOR`
 * @param[in] smooth_s  averaging half-width in seconds
 * @returns the envelope at the input sample rate
 * @exceptsafe basic
 */
std::vector<double> speed_envelope(
    const algo::Input& in,
    int mpu,
    double smooth_s
) {
  std::vector<double> x = algo::motion_signal(in, mpu, 0.3);
  for (double& v : x) v = std::fabs(v);
  return box_smooth(x, size_t(smooth_s / in.dt));
}

/**
 * Segment the envelope into reaches with a hysteretic threshold.
 *
 * Thresholds are set relative to the take's own 20th and 99.5th percentiles
 * rather than in absolute dps, because gain varies with how vigorously the
 * operator moved; the 99.5th percentile rather than the maximum keeps one
 * knock from raising the bar for every other reach. A run must cross the
 * high threshold to count but is followed out to the low one, so the slow
 * tail of a reach stays attached, and runs separated by less than
 * `MIN_REST_S` are merged before anything shorter than `MIN_BURST_S` is
 * dropped.
 *
 * @param[in] env  smoothed speed envelope of one hand
 * @param[in] dt   sample interval in seconds
 * @returns the accepted bursts in time order, possibly empty when the take
 *          is flat
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
 * Draw both hands, z-scored, over the whole recording.
 *
 * Z-scoring is for display only: the robot is markedly weaker than the
 * operator, and at a shared scale its trace would be an almost flat line
 * next to the operator's.
 *
 * @param[in,out] im   tile being drawn into
 * @param[in,out] ax   panel axes, autoscaled here to fit both traces
 * @param[in]     in   timestamps, used for the elapsed-time axis
 * @param[in]     rob  robot signal
 * @param[in]     op   operator signal
 * @exceptsafe basic
 */
void draw_signals(
    canvas::Image& im,
    canvas::Axes& ax,
    const algo::Input& in,
    const std::vector<double>& rob,
    const std::vector<double>& op
) {
  const std::vector<double> zr = algo::znorm(rob);
  const std::vector<double> zo = algo::znorm(op);
  const size_t stride = std::max<size_t>(1, zr.size() / 1400);
  std::vector<double> xs, ys_rob, ys_op;
  for (size_t i = 0; i < zr.size() && i < in.t.size(); i += stride) {
    xs.push_back(in.t[i] - in.t.front());
    ys_rob.push_back(zr[i]);
    ys_op.push_back(zo[i]);
  }
  std::vector<double> both = ys_rob;
  both.insert(both.end(), ys_op.begin(), ys_op.end());
  canvas::autoscale(ax, xs, both);
  canvas::grid_and_ticks(im, ax, 10, 5);
  canvas::polyline(im, ax, xs, ys_rob, canvas::BLACK);
  canvas::polyline(im, ax, xs, ys_op, canvas::VIVID);
}

/**
 * Finish an abstaining result, drawing the signals when there are any.
 *
 * Abstaining is a legitimate outcome and not an error: the tool reports the
 * median over ten algorithms, so a burst-based method that found no burst
 * should stand aside rather than guess. `Result::note` is left as the
 * caller set it and explains which of the two ways it gave up.
 *
 * @param[in] r    partly filled result, taken by value and returned
 * @param[in] in   timestamps, or null to draw an empty tile
 * @param[in] rob  robot signal, or null
 * @param[in] op   operator signal, or null
 * @returns the same result with `ok` false and an explanatory tile
 * @exceptsafe basic
 */
algo::Result rejected(
    algo::Result r,
    const algo::Input* in = nullptr,
    const std::vector<double>* rob = nullptr,
    const std::vector<double>* op = nullptr
) {
  r.image = canvas::standard("MATCHED FILTER", r.tau_ms, false);
  if (in && rob && op && rob->size() > 32 && op->size() == rob->size()) {
    canvas::panel_label(
        r.image,
        0,
        "1  ROTATION RATE, NO BURST COULD BE MATCHED"
    );
    canvas::Axes a0 = canvas::panel(0);
    draw_signals(r.image, a0, *in, *rob, *op);
    canvas::axis_labels(r.image, a0, "TIME (s)", "RATE (z score)");
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
 * Match each operator burst against the robot trace to time the copy.
 *
 * Every detected reach on the operator hand, padded by `BURST_PAD_S`,
 * becomes a template that is slid over the robot signal from `LAG_MIN_MS`
 * to `LAG_MAX_MS`; the score at each lag is the normalised correlation over
 * the template window, so the robot's attenuation cancels in the
 * denominator and only the shape has to agree. The argmax of that curve,
 * refined by a parabola through its two neighbours, is one delay estimate
 * per burst, and the reported value is the median across bursts with the
 * interquartile range as the spread. Working burst by burst is what suits
 * bursty data: a single whole-recording correlation is dominated by
 * whichever reach was largest, whereas a per-burst median gives every reach
 * one vote and is unbothered by the robot low-passing the operator. Bursts
 * whose best lag lands on either end of the search window are discarded as
 * unresolved rather than reported at the rail. `Result::note` gives how
 * many bursts were matched and the median peak correlation, which is the
 * number to look at when the estimate seems unstable; panel 1 shows the
 * whole recording with the robot slid back by the measured delay, so the
 * traces should sit on top of one another if the answer is right.
 * The matched-filter formulation follows Turin, An Introduction to Matched Filters,
 * IRE Trans. Information Theory 1960.
 *
 * @param[in] in  both gyro magnitude channels at 1 kHz plus the sample step
 * @returns the delay in ms, positive when the robot lags the operator, with
 *          `ok` false when no burst was found or fewer than two survived
 *          the rail check
 * @exceptsafe basic
 */
algo::Result matched_filter(const algo::Input& in) {
  algo::Result r;
  const size_t n = in.t.size();
  const std::vector<double> op = algo::motion_signal(in, algo::OPERATOR, 0.3);
  const std::vector<double> rob = algo::motion_signal(in, algo::ROBOT, 0.3);
  const std::vector<double> seg =
      speed_envelope(in, algo::OPERATOR, SEG_SMOOTH_S);

  const std::vector<Burst> bursts = find_bursts(seg, in.dt);
  const int nbursts = (int)bursts.size();
  if (bursts.empty()) {
    r.note = "no motion bursts found on the operator hand";
    return rejected(r, &in, &rob, &op);
  }

  const int l0 = int(LAG_MIN_MS / 1000.0 / in.dt);
  const int l1 = int(LAG_MAX_MS / 1000.0 / in.dt);
  const size_t pad = size_t(BURST_PAD_S / in.dt);

  std::vector<double> tau_ms, quality, plot_curve;
  std::vector<double> plot_template, plot_match;
  double plot_best = -2;
  double plot_tau_ms = 0;
  for (const Burst& b : bursts) {
    const size_t a = b.begin > pad ? b.begin - pad : 0;
    const size_t z = std::min(n, b.end + pad);
    const size_t len = z - a;
    if (len < 32) continue;
    double mean_op = 0;
    for (size_t i = a; i < z; i++) mean_op += op[i];
    mean_op /= double(len);
    double norm_op = 0;
    for (size_t i = a; i < z; i++)
      norm_op += (op[i] - mean_op) * (op[i] - mean_op);
    norm_op = std::sqrt(norm_op);
    if (norm_op <= 0) continue;

    std::vector<double> curve;
    curve.reserve(size_t(l1 - l0 + 1));
    double best = -2;
    int best_l = l0;
    for (int l = l0; l <= l1; l++) {
      const long s = (long)a + l;
      const long e = (long)z + l;
      if (s < 0 || e > (long)n) {
        curve.push_back(0.0);
        continue;
      }
      double mean_rob = 0;
      for (long i = s; i < e; i++) mean_rob += rob[size_t(i)];
      mean_rob /= double(len);
      double num = 0, norm_rob = 0;
      for (size_t k = 0; k < len; k++) {
        const double u = op[a + k] - mean_op;
        const double v = rob[size_t(s) + k] - mean_rob;
        num += u * v;
        norm_rob += v * v;
      }
      const double den = norm_op * std::sqrt(norm_rob);
      const double v = den > 0 ? num / den : 0.0;
      curve.push_back(v);
      if (v > best) {
        best = v;
        best_l = l;
      }
    }
    if (best_l <= l0 || best_l >= l1) continue;
    const size_t i = size_t(best_l - l0);
    const double y0 = curve[i - 1], y1 = curve[i], y2 = curve[i + 1];
    const double den = y0 - 2 * y1 + y2;
    const double sub = std::fabs(den) > 1e-12 ? 0.5 * (y0 - y2) / den : 0.0;
    tau_ms.push_back((best_l + sub) * in.dt * 1000.0);
    quality.push_back(best);
    if (best > plot_best) {
      plot_best = best;
      plot_curve = curve;
      plot_tau_ms = tau_ms.back();
      plot_template.assign(op.begin() + (long)a, op.begin() + (long)z);
      plot_match.clear();
      const long s0 = (long)a + best_l;
      if (s0 >= 0 && s0 + (long)len <= (long)n)
        for (size_t k = 0; k < len; k++)
          plot_match.push_back(rob[size_t(s0) + k]);
    }
  }

  if (tau_ms.size() < 2) {
    char buf[128];
    snprintf(
        buf,
        sizeof buf,
        "only %d usable bursts of %d",
        (int)tau_ms.size(),
        nbursts
    );
    r.note = buf;
    return rejected(r, &in, &rob, &op);
  }

  r.ok = true;
  r.tau_ms = algo::median(tau_ms);
  r.spread_ms = algo::iqr(tau_ms);
  char buf[160];
  snprintf(
      buf,
      sizeof buf,
      "%d bursts matched, median peak correlation %.2f",
      (int)tau_ms.size(),
      algo::median(quality)
  );
  r.note = buf;

  r.image = canvas::standard("MATCHED FILTER", r.tau_ms, r.ok);
  canvas::caption(r.image, "BURST CORRELATION");

  canvas::panel_label(
      r.image,
      0,
      "1  WHOLE RECORDING, THE ROBOT SLID BACK BY THE DELAY"
  );
  {
    canvas::Axes a0 = canvas::panel(0);
    const long shift = std::lround(r.tau_ms / 1000.0 / in.dt);
    const size_t n = std::min(op.size(), rob.size());
    const size_t stride = std::max<size_t>(1, n / 1600);
    std::vector<double> xs, ys_op, ys_rob;
    for (size_t i = 0; i < n; i += stride) {
      const long j = long(i) + shift;
      if (j < 0 || j >= long(rob.size())) continue;
      xs.push_back(in.t[i] - in.t.front());
      ys_op.push_back(op[i]);
      ys_rob.push_back(rob[size_t(j)]);
    }
    if (xs.size() > 8) {
      std::vector<double> both = ys_op;
      both.insert(both.end(), ys_rob.begin(), ys_rob.end());
      canvas::autoscale(a0, xs, both);
      canvas::set_range(a0, xs.front(), xs.back(), a0.y0, a0.y1);
      canvas::grid_and_ticks(r.image, a0, 10, 5);
      canvas::polyline(r.image, a0, xs, ys_rob, canvas::BLACK);
      canvas::polyline(r.image, a0, xs, ys_op, canvas::VIVID);
      canvas::axis_labels(r.image, a0, "TIME (s)", "RATE (dps)");
      canvas::legend(
          r.image,
          a0,
          {{"OPERATOR", canvas::VIVID}, {"ROBOT SHIFTED", canvas::BLACK}},
          true
      );
    }
  }

  canvas::panel_label(
      r.image,
      1,
      "2  HOW WELL THEY MATCH AT EVERY LAG, THE PEAK IS THE DELAY"
  );
  if (plot_curve.size() > 8) {
    canvas::Axes a1 = canvas::panel(1);
    const size_t stride = std::max<size_t>(1, plot_curve.size() / 1400);
    std::vector<double> xs, ys;
    for (size_t i = 0; i < plot_curve.size(); i += stride) {
      xs.push_back((double(l0) + double(i)) * in.dt * 1000.0);
      ys.push_back(plot_curve[i]);
    }
    canvas::autoscale(a1, xs, ys);
    canvas::grid_and_ticks(r.image, a1, 10, 5);
    canvas::polyline(r.image, a1, xs, ys, canvas::SAGE);
    canvas::mark_guide(r.image, a1, 0.0, true, "ZERO LAG", 0);
    canvas::marker(r.image, a1, plot_tau_ms, plot_best, canvas::INK, 5);
    canvas::mark_guide(r.image, a1, plot_tau_ms, true, "BEST LAG", 1);
    canvas::axis_labels(r.image, a1, "LAG (ms)", "CORRELATION");
    canvas::legend(r.image, a1, {{"CORRELATION", canvas::SAGE}}, true);
  }

  canvas::panel_label(
      r.image,
      2,
      "3  THE DELAY FROM EACH BURST, AND THE MEDIAN"
  );
  {
    canvas::Axes a2 = canvas::panel(2);
    double lo = tau_ms.front(), hi = tau_ms.front();
    for (double v : tau_ms) {
      lo = std::min(lo, v);
      hi = std::max(hi, v);
    }
    const double pad = std::max(5.0, 0.15 * (hi - lo));
    canvas::set_range(a2, 0.4, double(tau_ms.size()) + 0.6, lo - pad, hi + pad);
    canvas::grid_and_ticks(r.image, a2, 10, 5);
    canvas::mark_guide(r.image, a2, r.tau_ms, false, "MEDIAN", 0);
    for (size_t i = 0; i < tau_ms.size(); i++)
      canvas::marker(r.image, a2, double(i + 1), tau_ms[i], canvas::SAGE, 5);
    canvas::axis_labels(r.image, a2, "BURST NUMBER", "DELAY (ms)");
    canvas::legend(r.image, a2, {{"ONE BURST", canvas::WHITE, true}}, true);
  }

  return r;
}

}

ALGO_REGISTER(
    "matched_filter",
    matched_filter
);
