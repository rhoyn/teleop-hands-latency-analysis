#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <vector>

#include "algo.hpp"

namespace {

/**
 * Tuning for the whole-recording warp.
 *
 * The band is deliberately far wider than any delay the rig could produce:
 * it exists to constrain the warp and bound the cost, not to encode a prior
 * on the answer, and narrowing it towards the expected latency would make
 * the estimate partly self-fulfilling. The smoothing and envelope lengths
 * match the burst variant on purpose, so that where the two disagree they
 * disagree about the alignment method rather than about preprocessing.
 */
constexpr double BAND_MS = 250.0;     ///< warp may stray past any real delay
constexpr double SMOOTH_MS = 8.0;     ///< strips sensor noise, keeps the motion
constexpr double ENV_MS = 60.0;       ///< envelope tracks a single reach's rise
constexpr double ACTIVE_FRAC = 0.40;  ///< 40% of the 95th pct counts as moving
constexpr double WARP_PENALTY = 0.9;  ///< without it the path collapses to runs

/**
 * Smooth a signal with a centred box filter of `w` samples.
 *
 * A prefix sum keeps this linear in the recording length, which matters
 * because the envelope is rebuilt from scratch over hundreds of thousands
 * of 1 kHz samples. The window shrinks at both edges instead of assuming
 * zeros beyond them, so the first and last samples average only what
 * actually exists and the envelope does not dip into a false rest there.
 *
 * @param[in] x  signal to smooth
 * @param[in] w  window length in samples; 0 or 1 returns `x` unchanged
 * @returns the smoothed signal, the same length as `x`
 * @exceptsafe basic
 */
std::vector<double> moving_average(
    const std::vector<double>& x,
    size_t w
) {
  const size_t n = x.size();
  if (w < 2 || n == 0) return x;
  std::vector<double> ps(n + 1, 0.0);
  for (size_t i = 0; i < n; i++) ps[i + 1] = ps[i] + x[i];
  std::vector<double> out(n, 0.0);
  const size_t h = w / 2;
  for (size_t i = 0; i < n; i++) {
    const size_t a = i > h ? i - h : 0;
    const size_t b = std::min(n, i + h + 1);
    out[i] = (ps[b] - ps[a]) / double(b - a);
  }
  return out;
}

/**
 * Compute the running RMS of a signal over `w` samples.
 *
 * This is the motion envelope that decides which samples count as active.
 * Squaring before smoothing is what makes it blind to the sign of the
 * rotation rate, so a fast oscillation reads as motion instead of averaging
 * away to zero the way a plain moving average of the high-passed signal
 * would.
 *
 * @param[in] x  high-passed rotation rate
 * @param[in] w  window length in samples
 * @returns the envelope in the same units as `x`, the same length as `x`
 * @exceptsafe basic
 */
std::vector<double> moving_rms(
    const std::vector<double>& x,
    size_t w
) {
  std::vector<double> sq(x.size());
  for (size_t i = 0; i < x.size(); i++) sq[i] = x[i] * x[i];
  std::vector<double> m = moving_average(sq, w);
  for (double& v : m) v = std::sqrt(std::max(0.0, v));
  return m;
}

/**
 * Return the value at a fractional rank of a sample set.
 *
 * Used at 0.95 to set the energy scale, deliberately not at 1.0: a single
 * knock against the table would otherwise raise the gate high enough to
 * mark the whole recording as rest. The vector is taken by value because it
 * is sorted in place, and the rank is truncated rather than interpolated,
 * which is immaterial over a recording this long.
 *
 * @param[in] v  samples, in any order
 * @param[in] p  rank in [0, 1]
 * @returns the sample at that rank, or 0 when `v` is empty
 * @exceptsafe basic
 */
double percentile(
    std::vector<double> v,
    double p
) {
  if (v.empty()) return 0;
  std::sort(v.begin(), v.end());
  return v[size_t(p * double(v.size() - 1))];
}

/**
 * Standardise a signal using only the samples flagged as motion.
 *
 * Rest dominates these recordings by sample count, so statistics taken over
 * everything would be set by the sensor noise floor and each hand would end
 * up scaled by whatever its own noise happened to be. Restricting to active
 * samples makes the two z-scores genuinely comparable, which is what allows
 * the warp to use a plain absolute difference as its local cost. Fewer than
 * eight active samples falls back to the whole signal rather than dividing
 * by a noise estimate.
 *
 * @param[in,out] x       signal standardised in place
 * @param[in]     active  1 where the envelope says the hand is moving
 * @exceptsafe no-throw
 */
void znorm_masked(
    std::vector<double>& x,
    const std::vector<char>& active
) {
  double s = 0, s2 = 0;
  size_t c = 0;
  for (size_t i = 0; i < x.size(); i++)
    if (active[i]) {
      s += x[i];
      s2 += x[i] * x[i];
      c++;
    }
  if (c < 8) {
    c = x.size();
    s = 0;
    s2 = 0;
    for (double v : x) {
      s += v;
      s2 += v * v;
    }
  }
  if (c == 0) return;
  const double m = s / double(c);
  const double sd = std::sqrt(std::max(1e-12, s2 / double(c) - m * m));
  for (double& v : x) v = (v - m) / sd;
}

/**
 * Warp two series of near-equal length against each other inside a band.
 *
 * The Sakoe-Chiba band of `band` samples does two jobs at once: it bounds
 * the cost to O(n * band) instead of O(n * m), and it caps how far the path
 * may stray from the diagonal, which is the only thing stopping the warp
 * from explaining a quiet stretch by an enormous shift. `pen` is charged on
 * every non-diagonal step and is not optional — once both series are
 * z-scored, repeating a sample costs almost nothing, so an unpenalised path
 * collapses into a handful of long many-to-one runs and the offsets read
 * off it become meaningless. Only the one-byte backtrace direction is
 * stored per band cell, never the full cost matrix, which is what keeps a
 * whole-recording warp inside memory at 1 kHz.
 *
 * @param[in]  a     first series, indexed by `pi`
 * @param[in]  b     second series, indexed by `pj`
 * @param[in]  band  half-width of the band, in samples
 * @param[in]  pen   cost added to each insertion or deletion step
 * @param[out] pi    index into `a` for each path step, ascending
 * @param[out] pj    index into `b` matched to each entry of `pi`
 * @returns false when either series is under two samples or their lengths
 *          differ by more than `band`, leaving `pi` and `pj` untouched
 * @exceptsafe basic
 */
bool banded_dtw(
    const std::vector<double>& a,
    const std::vector<double>& b,
    int band,
    double pen,
    std::vector<int>& pi,
    std::vector<int>& pj
) {
  const int n = (int)a.size();
  const int m = (int)b.size();
  if (n < 2 || m < 2 || std::abs(n - m) > band) return false;
  const int w = 2 * band + 1;
  const double inf = 1e300;
  std::vector<double> prev(m, inf), cur(m, inf);
  std::vector<uint8_t> step((size_t)n * (size_t)w, 0);
  for (int i = 0; i < n; i++) {
    const int lo = std::max(0, i - band);
    const int hi = std::min(m - 1, i + band);
    for (int j = lo; j <= hi; j++) {
      const double c = std::fabs(a[i] - b[j]);
      double best = 0;
      uint8_t s = 0;
      if (i > 0 || j > 0) {
        const double dg = (i > 0 && j > 0) ? prev[j - 1] : inf;
        const double up = (i > 0) ? prev[j] : inf;
        const double lf = (j > lo) ? cur[j - 1] : inf;
        best = dg;
        s = 1;
        if (up + pen < best) {
          best = up + pen;
          s = 2;
        }
        if (lf + pen < best) {
          best = lf + pen;
          s = 3;
        }
      }
      cur[j] = (best >= inf) ? inf : best + c;
      step[(size_t)i * (size_t)w + (size_t)(j - i + band)] = s;
    }
    prev = cur;
    std::fill(cur.begin(), cur.end(), inf);
  }
  int i = n - 1, j = m - 1;
  pi.clear();
  pj.clear();
  while (true) {
    pi.push_back(i);
    pj.push_back(j);
    if (i == 0 && j == 0) break;
    const uint8_t s = step[(size_t)i * (size_t)w + (size_t)(j - i + band)];
    if (s == 1) {
      i--;
      j--;
    } else if (s == 2) {
      i--;
    } else if (s == 3) {
      j--;
    } else {
      break;
    }
  }
  std::reverse(pi.begin(), pi.end());
  std::reverse(pj.begin(), pj.end());
  return true;
}

/**
 * Start the tile for a run that will abstain.
 *
 * Abstaining is a legitimate outcome — the reported latency is the median
 * over the algorithms that did answer — so a tile is still produced and
 * simply marked as not ok, keeping the contact sheet's layout intact.
 *
 * @param[in,out] r  result whose image is created and captioned
 * @exceptsafe basic
 */
void reject_tile(algo::Result& r) {
  r.image = canvas::standard("DTW OFFSET", r.tau_ms, r.ok);
  canvas::caption(r.image, "WARP PATH OFFSET");
}

/**
 * Draw both z-scored traces overlaid against time.
 *
 * Points are decimated to about 3000, all an 800 px panel can resolve. The
 * two share one autoscaled y axis, which is meaningful here precisely
 * because both have already been standardised over their active samples:
 * what remains visible is shape disagreement, not the robot's attenuation.
 *
 * @param[in,out] r       result whose image is drawn into
 * @param[in]     t       timestamps in seconds
 * @param[in]     robot   z-scored robot trace
 * @param[in]     operat  z-scored operator trace
 * @exceptsafe basic
 */
void draw_signals(
    algo::Result& r,
    const std::vector<double>& t,
    const std::vector<double>& robot,
    const std::vector<double>& operat
) {
  canvas::panel_label(
      r.image,
      0,
      "1  BOTH HANDS NORMALISED AND OVERLAID OVER TIME"
  );
  const size_t n = std::min(t.size(), std::min(robot.size(), operat.size()));
  if (n < 4) return;
  const size_t stride = std::max<size_t>(1, n / 3000);
  std::vector<double> ts, vr, vo;
  for (size_t i = 0; i < n; i += stride) {
    ts.push_back(t[i] - t.front());
    vr.push_back(robot[i]);
    vo.push_back(operat[i]);
  }
  std::vector<double> both = vr;
  both.insert(both.end(), vo.begin(), vo.end());
  canvas::Axes a0 = canvas::panel(0);
  canvas::autoscale(a0, ts, both);
  canvas::set_range(a0, ts.front(), ts.back(), a0.y0, a0.y1);
  canvas::grid_and_ticks(r.image, a0, 10, 5);
  canvas::polyline(r.image, a0, ts, vr, canvas::BLACK);
  canvas::polyline(r.image, a0, ts, vo, canvas::VIVID);
  canvas::axis_labels(r.image, a0, "TIME (s)", "Z SCORE");
  canvas::legend(
      r.image,
      a0,
      {{"OPERATOR", canvas::VIVID}, {"ROBOT", canvas::BLACK}},
      true
  );
}

/**
 * Estimate the delay from the offset of a single whole-recording warp.
 *
 * Both gyro magnitudes are high-passed at 0.3 Hz, box-smoothed over
 * `SMOOTH_MS` and z-scored over their active samples, then warped against
 * each other in one banded DTW of half-width `BAND_MS`. Each operator
 * sample's mean matched robot index gives an offset in samples; those are
 * converted to ms, any whose robot or operator endpoint falls in a rest
 * region is thrown away, and the median of the survivors is the delay with
 * their IQR as the spread.
 *
 * Warping everything at once is what distinguishes this from the burst
 * variant: it needs no segmentation and uses every moving sample, so it is
 * the more stable of the two on long recordings, and it copes with the fact
 * that the robot attenuates and low-passes the operator rather than
 * delaying it cleanly. The weakness is the mirror image — the path is a
 * single object, so one misaligned stretch drags a long run of it along
 * with it and the offsets stop being independent, which is why the spread
 * is the IQR of the offsets and should not be read as a fit uncertainty.
 * The per-step warp penalty was mandatory to get here at all: without it
 * the cheapest path is a few long many-to-one runs and every offset
 * collapses onto a handful of values.
 *
 * `Result::note` reports the path length, the share of offsets that
 * survived the rest gate — a low share means most of the recording was
 * still — and the band half-width. The three panels show the two
 * normalised traces, the warp offset against time so a drifting or
 * jumping path is visible, and the sorted offsets as a quantile curve with
 * the median marked, where a flat plateau around 50% is the sign of a
 * well-determined delay.
 * The band-constrained dynamic programming follows Sakoe and Chiba, Dynamic
 * programming algorithm optimization for spoken word recognition, IEEE Trans. ASSP
 * 1978; the per-step warp penalty is an addition, without which the path collapses
 * into long many-to-one runs on this data.
 *
 * @param[in] in  both 1 kHz gyro magnitudes with their shared timebase
 * @returns the delay in ms, positive when the robot lags the operator, or
 *          an abstaining `Result` when the warp could not be run or fewer
 *          than 64 offsets landed in motion on both sides
 * @exceptsafe basic
 */
algo::Result dtw_offset(const algo::Input& in) {
  algo::Result r;
  const size_t n = in.t.size();
  const size_t sw = std::max<size_t>(1, size_t(SMOOTH_MS / 1000.0 / in.dt));
  const size_t ew = std::max<size_t>(1, size_t(ENV_MS / 1000.0 / in.dt));
  std::vector<double> robot =
      moving_average(algo::motion_signal(in, algo::ROBOT, 0.3), sw);
  std::vector<double> operat =
      moving_average(algo::motion_signal(in, algo::OPERATOR, 0.3), sw);
  const std::vector<double> er = moving_rms(robot, ew);
  const std::vector<double> eo = moving_rms(operat, ew);
  const double gr = ACTIVE_FRAC * percentile(er, 0.95);
  const double go = ACTIVE_FRAC * percentile(eo, 0.95);
  std::vector<char> ar(n), ao(n);
  for (size_t i = 0; i < n; i++) {
    ar[i] = er[i] > gr ? 1 : 0;
    ao[i] = eo[i] > go ? 1 : 0;
  }
  znorm_masked(robot, ar);
  znorm_masked(operat, ao);

  const int band = std::max(8, int(BAND_MS / 1000.0 / in.dt));
  std::vector<int> pi, pj;
  if (!banded_dtw(robot, operat, band, WARP_PENALTY, pi, pj)) {
    r.note = "recording too short for a banded warp";
    reject_tile(r);
    draw_signals(r, in.t, robot, operat);
    return r;
  }

  std::vector<double> matched_sum(n, 0.0), matched_count(n, 0.0);
  for (size_t k = 0; k < pi.size(); k++) {
    matched_sum[(size_t)pj[k]] += double(pi[k]);
    matched_count[(size_t)pj[k]] += 1.0;
  }
  std::vector<double> offsets;
  for (size_t j = 0; j < n; j++) {
    if (!ao[j] || matched_count[j] < 1.0) continue;
    const long i = std::lround(matched_sum[j] / matched_count[j]);
    if (i < 0 || i >= (long)n || !ar[(size_t)i]) continue;
    offsets.push_back(double(i - (long)j) * in.dt * 1000.0);
  }
  const double kept_fraction = double(offsets.size()) / double(n);
  if (offsets.size() < 64) {
    r.note = "almost the whole warp path fell in rest regions";
    reject_tile(r);
    draw_signals(r, in.t, robot, operat);
    return r;
  }

  r.ok = true;
  r.tau_ms = algo::median(offsets);
  r.spread_ms = algo::iqr(offsets);
  char buf[128];
  snprintf(
      buf,
      sizeof buf,
      "%d path steps, %.0f%% in motion, band +/-%.0f ms",
      (int)pi.size(),
      100.0 * kept_fraction,
      BAND_MS
  );
  r.note = buf;

  r.image = canvas::standard("DTW OFFSET", r.tau_ms, r.ok);
  canvas::caption(r.image, "WARP PATH OFFSET");
  draw_signals(r, in.t, robot, operat);

  canvas::panel_label(
      r.image,
      1,
      "2  HOW FAR THE WARP SHIFTS THE HANDS OVER TIME"
  );
  const size_t stride = std::max<size_t>(1, pi.size() / 1200);
  std::vector<double> ptime, poff;
  for (size_t k = 0; k < pi.size(); k += stride) {
    ptime.push_back(in.t[(size_t)pj[k]] - in.t.front());
    poff.push_back(double(pi[k] - pj[k]) * in.dt * 1000.0);
  }
  canvas::Axes a1 = canvas::panel(1);
  canvas::autoscale(a1, ptime, poff);
  canvas::grid_and_ticks(r.image, a1, 10, 5);
  canvas::polyline(r.image, a1, ptime, poff, canvas::SLATE);
  canvas::mark_guide(r.image, a1, r.tau_ms, false, "MEDIAN", 0);
  canvas::axis_labels(r.image, a1, "TIME (s)", "OFFSET (ms)");
  canvas::legend(r.image, a1, {{"WARP OFFSET", canvas::SLATE}}, true);

  canvas::panel_label(
      r.image,
      2,
      "3  OFFSETS SORTED, THE MIDDLE VALUE IS THE DELAY"
  );
  {
    std::vector<double> sorted = offsets;
    std::sort(sorted.begin(), sorted.end());
    const size_t ostride = std::max<size_t>(1, sorted.size() / 1200);
    std::vector<double> qx, qy;
    for (size_t i = 0; i < sorted.size(); i += ostride) {
      qx.push_back(100.0 * double(i) / double(sorted.size() - 1));
      qy.push_back(sorted[i]);
    }
    canvas::Axes a2 = canvas::panel(2);
    canvas::autoscale(a2, qx, qy);
    canvas::set_range(a2, 0.0, 100.0, a2.y0, a2.y1);
    canvas::grid_and_ticks(r.image, a2, 10, 5);
    canvas::polyline(r.image, a2, qx, qy, canvas::SAGE);
    canvas::mark_guide(r.image, a2, r.tau_ms, false, "MEDIAN", 0);
    canvas::axis_labels(r.image, a2, "PERCENTILE (%)", "OFFSET (ms)");
    canvas::legend(r.image, a2, {{"SORTED OFFSETS", canvas::SAGE}}, true);
  }
  return r;
}

}

ALGO_REGISTER(
    "dtw_offset",
    dtw_offset
);
