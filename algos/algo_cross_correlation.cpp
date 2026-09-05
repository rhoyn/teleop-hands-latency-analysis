#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdio>
#include <vector>

#include "algo.hpp"

namespace {

/**
 * Windowing and rest gating for the correlation estimator.
 *
 * A window has to be long enough that its correlation peak is supported by
 * several reaches yet short enough that the delay can be treated as
 * constant across it, and three seconds is where these recordings settled.
 * The hop is far smaller than the window on purpose: neighbouring windows
 * are highly correlated, but the median only needs a population of lags,
 * not independent ones. The gate is a fraction of the loudest window rather
 * than an absolute rate, so the same value transfers between stacks whose
 * motion amplitudes differ by an order of magnitude.
 */
constexpr double WIN_S = 3.0;       ///< several reaches, delay still constant
constexpr double HOP_S = 0.5;       ///< sixfold overlap, plenty of windows
constexpr double REST_GATE = 0.05;  ///< below 5% of peak energy is noise only

/**
 * Draw both high-passed hands against time into one panel.
 *
 * The traces are decimated to at most about 2000 points each, all an 800 px
 * wide panel can resolve, and share one autoscaled y axis so the robot's
 * attenuation relative to the operator stays visible rather than being
 * normalised away.
 *
 * @param[in,out] im   tile being drawn into
 * @param[in]     idx  panel index, 0 to 2 running top to bottom
 * @param[in]     t    timestamps in seconds
 * @param[in]     a    robot trace
 * @param[in]     b    operator trace
 * @exceptsafe basic
 */
void plot_pair(
    canvas::Image& im,
    int idx,
    const std::vector<double>& t,
    const std::vector<double>& a,
    const std::vector<double>& b
) {
  const size_t n = std::min(t.size(), std::min(a.size(), b.size()));
  if (n < 4) return;
  const size_t stride = std::max<size_t>(1, n / 2000);
  std::vector<double> ts, as, bs;
  for (size_t i = 0; i < n; i += stride) {
    ts.push_back(t[i]);
    as.push_back(a[i]);
    bs.push_back(b[i]);
  }
  std::vector<double> both(as);
  both.insert(both.end(), bs.begin(), bs.end());
  canvas::Axes ax = canvas::panel(idx);
  canvas::autoscale(ax, ts, both);
  canvas::set_range(ax, ts.front(), ts.back(), ax.y0, ax.y1);
  canvas::grid_and_ticks(im, ax, 10, 5);
  canvas::polyline(im, ax, ts, as, canvas::BLACK);
  canvas::polyline(im, ax, ts, bs, canvas::VIVID);
  canvas::axis_labels(im, ax, "TIME (s)", "RATE (dps)");
  canvas::legend(
      im,
      ax,
      {{"OPERATOR", canvas::VIVID}, {"ROBOT", canvas::BLACK}},
      true
  );
}

/**
 * Build the abstaining result, still showing the two input traces.
 *
 * Abstaining is a legitimate outcome rather than an error: the reported
 * latency is the median over the ten algorithms and simply forms without
 * this one. The tile is still drawn, first panel only, so a reader can see
 * whether the recording really was too short or too quiet.
 *
 * @param[in] note  short reason recorded in `Result::note`
 * @param[in] t     timestamps in seconds
 * @param[in] x     high-passed robot trace
 * @param[in] y     high-passed operator trace
 * @returns a `Result` with `ok` false and a single populated panel
 * @exceptsafe basic
 */
algo::Result rejected(
    const char* note,
    const std::vector<double>& t,
    const std::vector<double>& x,
    const std::vector<double>& y
) {
  algo::Result r;
  r.note = note;
  r.image = canvas::standard("CROSS CORRELATION", 0.0, false);
  canvas::caption(r.image, "NO USABLE WINDOWS");
  canvas::panel_label(r.image, 0, "1  HIGH PASSED ROTATION RATE OF BOTH HANDS");
  plot_pair(r.image, 0, t, x, y);
  return r;
}

/**
 * Find the lag that maximises the cross-correlation of two windows.
 *
 * The correlation is evaluated by FFT over a length of at least twice the
 * window, so the zero padding makes it linear rather than circular and no
 * energy wraps around between the two ends of the lag axis. Only lags
 * within `algo::MAX_LAG_MS` are searched, which both bounds the answer to
 * something physically plausible and stops a distant spurious peak from
 * winning. The peak is refined by a parabola through it and its two
 * neighbours because the 1 ms grid is coarser than the differences between
 * stacks being measured. Sign follows `a` correlated against `b`, so with
 * the robot passed as `a` a positive result means the robot happens later.
 *
 * @param[in]  a       robot window
 * @param[in]  b       operator window, the same length as `a`
 * @param[in]  dt      sample period in seconds
 * @param[out] lag_ms  lag axis of the full correlation, filled only when
 *                     both `lag_ms` and `value` are non-null
 * @param[out] value   correlation per lag, scaled to peak magnitude 1
 * @returns the interpolated peak lag in ms
 * @exceptsafe basic
 */
double peak_lag_ms(
    const std::vector<double>& a,
    const std::vector<double>& b,
    double dt,
    std::vector<double>* lag_ms = nullptr,
    std::vector<double>* value = nullptr
) {
  const size_t n = a.size();
  const size_t fn = algo::next_pow2(2 * n);
  std::vector<algo::cd> A(fn, 0.0), B(fn, 0.0);
  for (size_t i = 0; i < n; i++) {
    A[i] = a[i];
    B[i] = b[i];
  }
  algo::fft(A, false);
  algo::fft(B, false);
  for (size_t i = 0; i < fn; i++) A[i] *= std::conj(B[i]);
  algo::fft(A, true);

  const int max_lag = int(algo::MAX_LAG_MS / 1000.0 / dt);
  double best = -1e300;
  int best_k = 0;
  for (int k = -max_lag; k <= max_lag; k++) {
    const size_t idx = (k >= 0) ? size_t(k) : fn - size_t(-k);
    const double v = A[idx].real();
    if (v > best) {
      best = v;
      best_k = k;
    }
  }
  const size_t i1 = (best_k >= 0) ? size_t(best_k) : fn - size_t(-best_k);
  const size_t i2 =
      (best_k + 1 >= 0) ? size_t(best_k + 1) : fn - size_t(-(best_k + 1));
  const size_t i0 =
      (best_k - 1 >= 0) ? size_t(best_k - 1) : fn - size_t(-(best_k - 1));
  const double y0 = A[i0].real(), y1 = A[i1].real(), y2 = A[i2].real();
  const double den = y0 - 2 * y1 + y2;
  const double sub = (std::fabs(den) > 1e-12) ? 0.5 * (y0 - y2) / den : 0.0;

  if (lag_ms && value) {
    double scale = 0;
    for (int k = -max_lag; k <= max_lag; k++) {
      const size_t idx = (k >= 0) ? size_t(k) : fn - size_t(-k);
      scale = std::max(scale, std::fabs(A[idx].real()));
    }
    if (scale <= 0) scale = 1.0;
    lag_ms->clear();
    value->clear();
    for (int k = -max_lag; k <= max_lag; k++) {
      const size_t idx = (k >= 0) ? size_t(k) : fn - size_t(-k);
      lag_ms->push_back(double(k) * dt * 1000.0);
      value->push_back(A[idx].real() / scale);
    }
  }
  return (best_k + sub) * dt * 1000.0;
}

/**
 * Estimate the delay from the correlation peak of overlapping windows.
 *
 * Both gyro magnitudes are high-passed at 0.3 Hz to strip gravity and slow
 * bias, then cut into `WIN_S` windows every `HOP_S`. A window carrying less
 * than `REST_GATE` of the loudest window's energy is discarded as rest,
 * since correlating two noise floors returns an essentially arbitrary lag
 * that would still enter the median with full weight. Every surviving
 * window contributes the lag of its correlation peak and the recording's
 * delay is the median of those, with the IQR as the spread.
 *
 * This is the bluntest of the ten estimators and, partly for that reason,
 * among the most robust. It assumes only that the robot resembles a shifted
 * operator, which is false in detail — the robot attenuates and low-passes
 * — but that distortion is close to symmetric in time, so it mostly lowers
 * the peak rather than moving it. Where it does struggle is amplitude:
 * windows are not normalised beyond the high-pass, so one loud burst inside
 * a window dominates the lag that window reports, and a window straddling
 * the boundary between rest and a reach reports the reach's lag with the
 * rest's noise added.
 *
 * `Result::note` reports how many windows survived the rest gate. The three
 * panels show the high-passed pair, one correlation curve with its peak
 * marked — the last accepted window's, since the arrays are overwritten
 * each time round — and every window's lag against its centre time with the
 * median drawn across.
 * This is the plain-weighted member of the generalized cross-correlation family of
 * Knapp and Carter, The Generalized Correlation Method for Estimation of Time Delay,
 * IEEE Trans. ASSP 1976.
 *
 * @param[in] in  both 1 kHz gyro magnitudes with their shared timebase
 * @returns the delay in ms, positive when the robot lags the operator, or
 *          an abstaining `Result` when the recording is shorter than one
 *          window or every window gated out as rest
 * @exceptsafe basic
 */
algo::Result cross_correlation(const algo::Input& in) {
  algo::Result r;
  const std::vector<double> x = algo::motion_signal(in, algo::ROBOT, 0.3);
  const std::vector<double> y = algo::motion_signal(in, algo::OPERATOR, 0.3);

  const size_t wn = size_t(WIN_S / in.dt);
  const size_t hop = size_t(HOP_S / in.dt);
  if (x.size() < wn)
    return rejected("recording shorter than one window", in.t, x, y);
  double energy_max = 0;
  for (size_t i = 0; i + wn <= x.size(); i += hop) {
    double e = 0;
    for (size_t k = 0; k < wn; k++)
      e += x[i + k] * x[i + k] + y[i + k] * y[i + k];
    energy_max = std::max(energy_max, e);
  }
  std::vector<double> taus;
  std::vector<double> win_t;
  std::vector<double> cc_lag_ms, cc_val;
  for (size_t i = 0; i + wn <= x.size(); i += hop) {
    double e = 0;
    for (size_t k = 0; k < wn; k++)
      e += x[i + k] * x[i + k] + y[i + k] * y[i + k];
    if (e < REST_GATE * energy_max) continue;
    const std::vector<double> xa(x.begin() + i, x.begin() + i + wn);
    const std::vector<double> ya(y.begin() + i, y.begin() + i + wn);
    taus.push_back(peak_lag_ms(xa, ya, in.dt, &cc_lag_ms, &cc_val));
    win_t.push_back(
        in.t.empty() ? 0.0 : in.t[std::min(in.t.size() - 1, i + wn / 2)]
    );
  }
  if (taus.empty())
    return rejected("every window gated out as rest", in.t, x, y);
  r.ok = true;
  r.tau_ms = algo::median(taus);
  r.spread_ms = algo::iqr(taus);
  char buf[128];
  snprintf(buf, sizeof buf, "%d windows accepted", (int)taus.size());
  r.note = buf;

  r.image = canvas::standard("CROSS CORRELATION", r.tau_ms, r.ok);
  canvas::caption(r.image, "CORR AND WINDOW TAUS");

  canvas::panel_label(r.image, 0, "1  HIGH PASSED ROTATION RATE OF BOTH HANDS");
  plot_pair(r.image, 0, in.t, x, y);

  canvas::panel_label(
      r.image,
      1,
      "2  CROSS CORRELATION, THE PEAK LAG IS THE DELAY"
  );
  if (cc_lag_ms.size() > 1) {
    canvas::Axes a1 = canvas::panel(1);
    canvas::autoscale(a1, cc_lag_ms, cc_val);
    canvas::grid_and_ticks(r.image, a1, 10, 5);
    canvas::polyline(r.image, a1, cc_lag_ms, cc_val, canvas::SAGE);
    canvas::mark_guide(r.image, a1, 0.0, true, "ZERO LAG", 0);
    size_t pk = 0;
    for (size_t i = 1; i < cc_val.size(); i++)
      if (cc_val[i] > cc_val[pk]) pk = i;
    canvas::marker(r.image, a1, cc_lag_ms[pk], cc_val[pk], canvas::WHITE, 4);
    canvas::mark_guide(r.image, a1, r.tau_ms, true, "ESTIMATE", 1);
    canvas::axis_labels(r.image, a1, "LAG (ms)", "CORRELATION");
    canvas::legend(
        r.image,
        a1,
        {{"CORRELATION", canvas::SAGE}, {"PEAK", canvas::WHITE, true}},
        true
    );
  }

  canvas::panel_label(
      r.image,
      2,
      "3  DELAY MEASURED IN EACH WINDOW AND THEIR MEDIAN"
  );
  if (!win_t.empty()) {
    canvas::Axes a2 = canvas::panel(2);
    canvas::autoscale(a2, win_t, taus);
    canvas::set_range(a2, win_t.front(), win_t.back(), a2.y0, a2.y1);
    canvas::grid_and_ticks(r.image, a2, 10, 5);
    for (size_t i = 0; i < win_t.size(); i++)
      canvas::marker(r.image, a2, win_t[i], taus[i], canvas::SLATE, 5);
    canvas::mark_guide(r.image, a2, r.tau_ms, false, "MEDIAN", 0);
    canvas::axis_labels(r.image, a2, "WINDOW CENTRE (s)", "DELAY (ms)");
    canvas::legend(r.image, a2, {{"ACCEPTED", canvas::WHITE, true}}, true);
  }
  return r;
}

}

ALGO_REGISTER(
    "cross_correlation",
    cross_correlation
);
