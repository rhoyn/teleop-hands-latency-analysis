#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <utility>
#include <vector>

#include "algo.hpp"

namespace {

/**
 * Tuning for the reported uncertainty.
 *
 * Neither number can move the MAP delay; they only shape the interval that
 * is reported as `Result::spread_ms`. The autocorrelation cap exists purely
 * to bound the cost of an O(m * kmax) sum over a residual that is hundreds
 * of thousands of samples long at 1 kHz.
 */
constexpr double CREDIBLE_MASS = 0.95;  ///< central interval reported as spread
constexpr int MAX_AUTOCORR_LAG = 1500;  ///< 1.5 s caps the correlation sum

/**
 * Estimate the integrated autocorrelation time of a residual series.
 *
 * At 1 kHz the residual is nowhere near white: the robot's own dynamics and
 * the parts of the operator's motion it fails to reproduce leave it
 * correlated over tens of milliseconds. Treating every sample as
 * independent evidence would therefore make the posterior absurdly sharp,
 * so this returns the factor the raw count is divided by. The sum is
 * truncated at the first non-positive lag, the usual initial-positive-
 * sequence rule, which keeps the noisy tail of the estimator out of the
 * total, and the result is floored at 1 so it can never inflate the count.
 *
 * @param[in] residual  robot minus the gain-scaled delayed operator
 * @returns 1 + 2 * sum of the normalised autocorrelation, at least 1, or 1
 *          for a series shorter than 32 samples or with no variance
 * @exceptsafe basic
 */
double integrated_autocorrelation_time(const std::vector<double>& residual) {
  const size_t m = residual.size();
  if (m < 32) return 1.0;
  double c0 = 0;
  for (double v : residual) c0 += v * v;
  c0 /= double(m);
  if (c0 <= 0) return 1.0;
  double sum = 0;
  const size_t kmax = std::min<size_t>(MAX_AUTOCORR_LAG, m / 4);
  for (size_t k = 1; k <= kmax; k++) {
    double c = 0;
    for (size_t i = 0; i + k < m; i++) c += residual[i] * residual[i + k];
    c /= double(m - k);
    const double rho = c / c0;
    if (rho <= 0) break;
    sum += rho;
  }
  return std::max(1.0, 1.0 + 2.0 * sum);
}

/**
 * Draw both hands z-scored and overlaid, for the abstaining tile.
 *
 * The z-scoring is for display only and is not what the estimator sees: the
 * robot's rotation rate is a fraction of the operator's, so raw traces
 * would leave one of them a flat line against the other's scale. Points are
 * decimated to about 1400, all an 800 px panel resolves.
 *
 * @param[in,out] im   tile being drawn into
 * @param[in,out] ax   panel axes, autoscaled in place
 * @param[in]     in   input, used only for its timebase
 * @param[in]     rob  robot trace
 * @param[in]     op   operator trace
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
  canvas::legend(
      im,
      ax,
      {{"OPERATOR", canvas::VIVID}, {"ROBOT", canvas::BLACK}},
      true
  );
}

/**
 * Draw the profile residual against the assumed delay.
 *
 * This is the raw evidence the posterior is built from, and reading it is
 * the quickest sanity check available: a broad, shallow trough here is what
 * a wide credible interval looks like before the likelihood is formed, and
 * a trough at the very edge of the grid means the delay ran past the 200 ms
 * search ceiling.
 *
 * @param[in,out] im       tile being drawn into
 * @param[in,out] ax       panel axes, autoscaled in place
 * @param[in]     rss      residual sum of squares per delay hypothesis
 * @param[in]     step_ms  milliseconds per hypothesis, 1 ms at 1 kHz
 * @param[in]     tau_ms   estimate to mark on the delay axis
 * @param[in]     mark     false to leave the estimate guide off entirely
 * @exceptsafe basic
 */
void draw_residual(
    canvas::Image& im,
    canvas::Axes& ax,
    const std::vector<double>& rss,
    double step_ms,
    double tau_ms,
    bool mark
) {
  std::vector<double> xs, ys;
  for (size_t i = 0; i < rss.size(); i++) {
    xs.push_back(double(i) * step_ms);
    ys.push_back(rss[i]);
  }
  if (xs.size() < 8) return;
  canvas::autoscale(ax, xs, ys);
  canvas::grid_and_ticks(im, ax, 10, 5);
  canvas::polyline(im, ax, xs, ys, canvas::SLATE);
  if (mark) canvas::mark_guide(im, ax, tau_ms, true, "ESTIMATE", 0);
  canvas::axis_labels(im, ax, "ASSUMED DELAY (ms)", "RESIDUAL");
  canvas::legend(im, ax, {{"RESIDUAL", canvas::SLATE}}, true);
}

/**
 * Finish an abstaining result, drawing whichever panels are available.
 *
 * The stages that can fail run in order, so the caller passes only the
 * intermediate data it actually reached and a null pointer simply leaves
 * that panel blank. Abstaining is a legitimate outcome and not an error:
 * the reported latency is the median over the algorithms that did answer.
 *
 * @param[in] r    partially filled result, taken by value and returned
 * @param[in] in   input for its timebase, or null
 * @param[in] rob  robot trace, or null
 * @param[in] op   operator trace, or null
 * @param[in] rss  residual per delay hypothesis, or null
 * @returns the result with `ok` left false and its tile drawn
 * @exceptsafe basic
 */
algo::Result rejected(
    algo::Result r,
    const algo::Input* in = nullptr,
    const std::vector<double>* rob = nullptr,
    const std::vector<double>* op = nullptr,
    const std::vector<double>* rss = nullptr
) {
  r.image = canvas::standard("BAYESIAN POSTERIOR", r.tau_ms, false);
  if (in && rob && op && rob->size() > 32 && op->size() == rob->size()) {
    canvas::panel_label(
        r.image,
        0,
        "1  NORMALISED ROTATION RATE OF BOTH HANDS"
    );
    canvas::Axes a0 = canvas::panel(0);
    draw_signals(r.image, a0, *in, *rob, *op);
    canvas::axis_labels(r.image, a0, "TIME (s)", "Z SCORE");
  }
  if (in && rss && rss->size() > 8) {
    canvas::panel_label(r.image, 1, "2  FIT RESIDUAL FOR EACH ASSUMED DELAY");
    canvas::Axes a1 = canvas::panel(1);
    draw_residual(r.image, a1, *rss, in->dt * 1000.0, 0.0, false);
  }
  return r;
}

/**
 * Estimate the delay and a credible interval from a Gaussian posterior.
 *
 * Both gyro magnitudes are high-passed at 0.3 Hz, and for every delay from
 * 0 to `algo::MAX_LAG_MS` the robot is regressed on the delayed operator,
 * robot[i] ~ a * operator[i - lag]. The amplitude a is profiled out
 * analytically as sxy / syy, which is what lets one scalar per hypothesis
 * absorb the robot's attenuation without adding a dimension to the search,
 * and the profile residual sxx - sxy^2 / syy is what is recorded. A
 * negative gain is rejected outright by falling back to sxx, since an
 * inverted robot is not a physical answer. Those residuals become a
 * Gaussian log likelihood, -0.5 * n_eff * log(rss / n), and a flat prior
 * over the grid turns them into a posterior; the MAP is refined by a
 * parabola on the log likelihood and the reported spread is the central
 * `CREDIBLE_MASS` interval read off the cumulative posterior.
 *
 * The effective sample size, not the raw count, is what makes that interval
 * honest, and it is the subtle part of this file. The residual at 1 kHz is
 * heavily autocorrelated, so the independent evidence actually present is
 * `samples` divided by the integrated autocorrelation time, frequently a
 * factor of tens; using the raw count would collapse the interval to a
 * fraction of a millisecond and claim a precision the recording does not
 * contain. The likelihood is scaled by the effective size while the
 * residual inside the logarithm is still divided by the raw `samples`,
 * which is harmless: that introduces the same constant offset at every
 * hypothesis and it cancels when the posterior is normalised. Where the
 * method struggles is that a single gain and a single shift cannot
 * represent a low-pass, and the robot is low-passing rather than copying,
 * so the trough it minimises tends to sit slightly late.
 *
 * `Result::note` reports the MAP delay, the 95% interval, the fitted gain —
 * well under 1, because the robot attenuates — and the effective sample
 * count; `Result::spread_ms` is the width of the interval. The three panels
 * show the robot against the gain-scaled delayed operator, the residual
 * curve over the whole grid, and the posterior density zoomed around the
 * MAP with both interval bounds marked.
 *
 * @param[in] in  both 1 kHz gyro magnitudes with their shared timebase
 * @returns the MAP delay in ms, positive when the robot lags the operator,
 *          or an abstaining `Result` when the recording is too short for a
 *          200 ms grid, the robot is silent, or the posterior underflowed
 * @exceptsafe basic
 */
algo::Result bayesian_posterior(const algo::Input& in) {
  algo::Result r;
  const std::vector<double> robot = algo::motion_signal(in, algo::ROBOT, 0.3);
  const std::vector<double> operator_hand =
      algo::motion_signal(in, algo::OPERATOR, 0.3);
  const size_t n = robot.size();
  const int max_lag = int(algo::MAX_LAG_MS / 1000.0 / in.dt);
  if (n < size_t(max_lag) * 4) {
    r.note = "recording too short for a 200 ms delay grid";
    return rejected(r, &in, &robot, &operator_hand);
  }
  const size_t first = size_t(max_lag);
  const double samples = double(n - first);

  double sxx = 0;
  for (size_t i = first; i < n; i++) sxx += robot[i] * robot[i];
  if (sxx <= 0) {
    r.note = "robot signal is silent";
    return rejected(r, &in, &robot, &operator_hand);
  }

  std::vector<double> residual_sum_squares;
  std::vector<double> amplitude;
  for (int lag = 0; lag <= max_lag; lag++) {
    double sxy = 0, syy = 0;
    for (size_t i = first; i < n; i++) {
      const double yv = operator_hand[i - size_t(lag)];
      sxy += robot[i] * yv;
      syy += yv * yv;
    }
    const double a = syy > 0 ? sxy / syy : 0.0;
    const double rss = (a > 0) ? sxx - sxy * sxy / syy : sxx;
    residual_sum_squares.push_back(std::max(rss, 1e-300));
    amplitude.push_back(a);
  }

  int map_index = 0;
  for (int i = 1; i <= max_lag; i++)
    if (residual_sum_squares[i] < residual_sum_squares[map_index])
      map_index = i;
  const double gain = amplitude[map_index];

  std::vector<double> residual;
  residual.reserve(size_t(samples));
  for (size_t i = first; i < n; i++)
    residual.push_back(robot[i] - gain * operator_hand[i - size_t(map_index)]);
  const double autocorr_time = integrated_autocorrelation_time(residual);
  const double effective_samples = std::max(4.0, samples / autocorr_time);

  std::vector<double> loglik;
  double best_ll = -1e300;
  for (int i = 0; i <= max_lag; i++) {
    const double ll =
        -0.5 * effective_samples * std::log(residual_sum_squares[i] / samples);
    loglik.push_back(ll);
    best_ll = std::max(best_ll, ll);
  }
  std::vector<double> posterior;
  double mass = 0;
  for (int i = 0; i <= max_lag; i++) {
    const double p = std::exp(loglik[i] - best_ll);
    posterior.push_back(p);
    mass += p;
  }
  if (!(mass > 0)) {
    r.note = "posterior underflowed";
    return rejected(r, &in, &robot, &operator_hand, &residual_sum_squares);
  }
  const double step_ms = in.dt * 1000.0;
  for (double& p : posterior) p /= mass * step_ms;

  double sub = 0;
  if (map_index > 0 && map_index < max_lag) {
    const double y0 = loglik[map_index - 1], y1 = loglik[map_index],
                 y2 = loglik[map_index + 1];
    const double den = y0 - 2 * y1 + y2;
    if (std::fabs(den) > 1e-15) sub = 0.5 * (y0 - y2) / den;
    if (sub < -1 || sub > 1) sub = 0;
  }
  const double map_tau = (double(map_index) + sub) * step_ms;

  const double tail = 0.5 * (1.0 - CREDIBLE_MASS);
  double cumulative = 0;
  double ci_low = 0;
  double ci_high = algo::MAX_LAG_MS;
  bool low_found = false;
  for (int i = 0; i <= max_lag; i++) {
    const double slice = posterior[i] * step_ms;
    const double before = cumulative;
    cumulative += slice;
    if (!low_found && cumulative >= tail) {
      const double frac = slice > 0 ? (tail - before) / slice : 0.0;
      ci_low = (double(i) - 0.5 + frac) * step_ms;
      low_found = true;
    }
    if (cumulative >= 1.0 - tail) {
      const double frac = slice > 0 ? (1.0 - tail - before) / slice : 0.0;
      ci_high = (double(i) - 0.5 + frac) * step_ms;
      break;
    }
  }
  ci_low = std::max(0.0, ci_low);

  r.ok = true;
  r.tau_ms = map_tau;
  r.spread_ms = ci_high - ci_low;
  char buf[192];
  snprintf(
      buf,
      sizeof buf,
      "MAP %.1f ms, 95%% CI %.1f-%.1f ms, gain %.2f, %.0f eff samples",
      map_tau,
      ci_low,
      ci_high,
      gain,
      effective_samples
  );
  r.note = buf;

  r.image = canvas::standard("BAYESIAN POSTERIOR", r.tau_ms, r.ok);
  canvas::caption(r.image, "POSTERIOR OVER TAU");

  canvas::panel_label(
      r.image,
      0,
      "1  DELAYED OPERATOR TRACE LINES UP WITH THE ROBOT"
  );
  {
    const size_t ia = first;
    const size_t ib = n;
    if (ib > ia + 16 && ib <= in.t.size()) {
      canvas::Axes a0 = canvas::panel(0);
      const size_t stride = std::max<size_t>(1, (ib - ia) / 1400);
      std::vector<double> xs, ys_rob, ys_fit;
      for (size_t i = ia; i < ib; i += stride) {
        xs.push_back(in.t[i] - in.t[ia]);
        ys_rob.push_back(robot[i]);
        ys_fit.push_back(gain * operator_hand[i - size_t(map_index)]);
      }
      std::vector<double> both = ys_rob;
      both.insert(both.end(), ys_fit.begin(), ys_fit.end());
      canvas::autoscale(a0, xs, both);
      canvas::grid_and_ticks(r.image, a0, 10, 5);
      canvas::polyline(r.image, a0, xs, ys_rob, canvas::BLACK);
      canvas::polyline(r.image, a0, xs, ys_fit, canvas::VIVID);
      canvas::axis_labels(r.image, a0, "TIME (s)", "RATE (dps)");
      canvas::legend(
          r.image,
          a0,
          {{"DELAYED OPERATOR", canvas::VIVID}, {"ROBOT", canvas::BLACK}},
          true
      );
    }
  }

  canvas::panel_label(
      r.image,
      1,
      "2  FIT RESIDUAL PER ASSUMED DELAY, LOWEST IS BEST"
  );
  {
    canvas::Axes a1 = canvas::panel(1);
    draw_residual(r.image, a1, residual_sum_squares, step_ms, r.tau_ms, true);
  }

  canvas::panel_label(
      r.image,
      2,
      "3  HOW LIKELY EACH DELAY IS, WITH 95PCT BOUNDS"
  );
  {
    canvas::Axes a2 = canvas::panel(2);
    const double half_span = std::max(3.0 * (ci_high - ci_low), 12.0);
    const double x_lo = std::max(0.0, map_tau - half_span);
    const double x_hi = std::min(algo::MAX_LAG_MS, map_tau + half_span);
    std::vector<double> xs, ys;
    double top = 0;
    for (int i = 0; i <= max_lag; i++) {
      const double x = double(i) * step_ms;
      if (x < x_lo || x > x_hi) continue;
      xs.push_back(x);
      ys.push_back(posterior[size_t(i)]);
      top = std::max(top, posterior[size_t(i)]);
    }
    if (xs.size() > 8 && top > 0) {
      canvas::set_range(a2, x_lo, x_hi, 0.0, top * 1.15);
      canvas::grid_and_ticks(r.image, a2, 10, 5);
      canvas::polyline(r.image, a2, xs, ys, canvas::SAGE);
      int slot = 0;
      if (ci_low >= x_lo && ci_low <= x_hi)
        canvas::mark_guide(r.image, a2, ci_low, true, "95% BOUNDS", slot++);
      if (ci_high >= x_lo && ci_high <= x_hi)
        canvas::mark_guide(r.image, a2, ci_high, true, "95% BOUNDS", slot++);
      canvas::mark_guide(r.image, a2, r.tau_ms, true, "ESTIMATE", slot);
      canvas::axis_labels(r.image, a2, "DELAY (ms)", "DENSITY");
      canvas::legend(r.image, a2, {{"POSTERIOR", canvas::SAGE}}, true);
    }
  }

  return r;
}

}

ALGO_REGISTER(
    "bayesian_posterior",
    bayesian_posterior
);
