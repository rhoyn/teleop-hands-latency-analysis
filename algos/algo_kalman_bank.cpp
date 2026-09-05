#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

#include <Eigen/Dense>

#include "algo.hpp"

namespace {

/**
 * Preprocessing, plant order and bank geometry, tuned on 1 kHz takes.
 *
 * These are deliberately the same values `arx_delay` uses, so the two
 * algorithms see an identical signal and their votes are correlated rather
 * than independent evidence; treat an agreement between them as one opinion
 * confirmed twice, not as two. The orders are kept small because a richer
 * plant can imitate a delay with its own poles and flatten the likelihood
 * ridge this method depends on.
 */
constexpr int DECIMATE = 2;       ///< 500 Hz state, 1 ms delay grid kept
constexpr double LP_HZ = 25.0;    ///< the robot has nothing above this
constexpr int LP_TAPS = 101;      ///< linear phase, cancels across hands
constexpr int NA = 2;             ///< two poles cover the robot roll-off
constexpr int NB = 3;             ///< three input taps, no more free delay
constexpr int NK_MAX = 200;       ///< 200 ms ceiling, so 201 filters
constexpr double Q_OVER_R = 0.1;  ///< process noise as a fraction of R
constexpr double FORGET_S = 5.0;  ///< about one reach and its rest

/**
 * Build a Hamming-windowed sinc low pass of `LP_TAPS` taps.
 *
 * The taps are normalised to unit DC gain so the filter cannot rescale the
 * two hands differently, and the odd length gives an exactly linear phase
 * whose (`LP_TAPS` - 1)/2 sample group delay is identical on the robot and
 * the operator channel, so it cancels out of the measured delay entirely.
 *
 * @param[in] fc  cutoff in Hz
 * @param[in] dt  sample interval in seconds
 * @returns the tap vector, already normalised
 * @exceptsafe basic
 */
std::vector<double> lowpass_fir(
    double fc,
    double dt
) {
  std::vector<double> h(LP_TAPS);
  const double fn = fc * dt;
  double sum = 0;
  for (int i = 0; i < LP_TAPS; i++) {
    const double n = double(i) - double(LP_TAPS - 1) / 2.0;
    const double s = (std::fabs(n) < 1e-9)
                         ? 2 * fn
                         : std::sin(2 * M_PI * fn * n) / (M_PI * n);
    const double w =
        0.54 - 0.46 * std::cos(2 * M_PI * double(i) / double(LP_TAPS - 1));
    h[i] = s * w;
    sum += h[i];
  }
  for (double& v : h) v /= sum;
  return h;
}

/**
 * Convolve with a centred kernel, returning a signal of the same length.
 *
 * Samples off either end are treated as zero, which lets the first and last
 * `h.size()`/2 samples droop; the estimator starts well past that point and
 * the burst structure sits in the interior, so the edge is never fitted.
 *
 * @param[in] x  signal to filter
 * @param[in] h  kernel, indexed about its own centre tap
 * @returns the filtered signal, the same length as `x`
 * @exceptsafe basic
 */
std::vector<double> convolve_same(
    const std::vector<double>& x,
    const std::vector<double>& h
) {
  const int n = (int)x.size(), m = (int)h.size(), c = m / 2;
  std::vector<double> y(n, 0.0);
  for (int i = 0; i < n; i++) {
    double s = 0;
    for (int k = 0; k < m; k++) {
      const int j = i + c - k;
      if (j >= 0 && j < n) s += h[k] * x[j];
    }
    y[i] = s;
  }
  return y;
}

/**
 * One integer delay hypothesis: its ARX plant, Kalman state and scores.
 *
 * The index of a hypothesis inside the bank is its delay in samples, so no
 * field records it. `loglik` accumulates over the whole batch and picks the
 * single reported delay, while `logw` is the same innovation score decayed
 * by `FORGET_S` so it can follow a delay that drifts during the take.
 */
struct Hypothesis {
  double a1 = 0, a2 = 0;
  double b[NB] = {0, 0, 0};
  double var = 1e300;
  double x0 = 0, x1 = 0;
  double p00 = 0, p01 = 0, p10 = 0, p11 = 0;
  double logw = 0;
  double loglik = 0;
};

/**
 * Least-squares fit the ARX plant for one delay hypothesis.
 *
 * Each hypothesis gets its own plant because the operator regressors are
 * read at `nk` samples of lag, so a wrong `nk` must explain the robot with
 * mismatched input and pays for it in residual variance. Regression starts
 * far enough in that every index is in range for the largest lag, and the
 * QR solve is used rather than the normal equations because the operator
 * taps are strongly collinear at 500 Hz.
 *
 * @param[in]  yf   filtered robot signal at the full 1 kHz rate
 * @param[in]  uf   filtered operator signal at the full 1 kHz rate
 * @param[in]  nk   delay hypothesis in 1 kHz samples, hence in ms
 * @param[in]  nd   number of decimated samples available
 * @param[out] hyp  plant coefficients and residual variance
 * @returns false when too few rows remain to fit, or the variance is not
 *          finite, in which case the whole bank is abandoned
 * @exceptsafe basic
 */
bool fit_plant(
    const std::vector<double>& yf,
    const std::vector<double>& uf,
    int nk,
    int nd,
    Hypothesis& hyp
) {
  const int p = std::max(NA, NB - 1) + 1;
  const int i0 = std::max(p, (nk + DECIMATE - 1) / DECIMATE + p);
  const int rows = nd - i0;
  if (rows < 4 * (NA + NB)) return false;
  Eigen::MatrixXd X(rows, NA + NB);
  Eigen::VectorXd Y(rows);
  for (int r = 0; r < rows; r++) {
    const int t = i0 + r;
    X(r, 0) = -yf[(t - 1) * DECIMATE];
    X(r, 1) = -yf[(t - 2) * DECIMATE];
    for (int j = 0; j < NB; j++) X(r, NA + j) = uf[(t - j) * DECIMATE - nk];
    Y(r) = yf[t * DECIMATE];
  }
  const Eigen::VectorXd beta = X.colPivHouseholderQr().solve(Y);
  hyp.a1 = beta(0);
  hyp.a2 = beta(1);
  for (int j = 0; j < NB; j++) hyp.b[j] = beta(NA + j);
  hyp.var = (Y - X * beta).squaredNorm() / double(rows);
  return std::isfinite(hyp.var);
}

/**
 * Start the tile for an abstaining run so the panels can still be drawn.
 *
 * @param[in,out] r  result whose `image` is replaced with a fresh tile
 * @exceptsafe basic
 */
void reject_tile(algo::Result& r) {
  r.image = canvas::standard("KALMAN BANK", r.tau_ms, r.ok);
  canvas::caption(r.image, "HYPOTHESIS LOGLIK");
}

/**
 * Draw both filtered hands against time in the first panel.
 *
 * Called on the abstaining path too, so it tolerates mismatched lengths and
 * simply draws nothing when there is too little to show. The series is
 * strided down to a few thousand points because a tile is 800 px wide.
 *
 * @param[in,out] r   result whose `image` is drawn into
 * @param[in]     t   timestamps, only the elapsed part is used
 * @param[in]     yf  filtered robot signal
 * @param[in]     uf  filtered operator signal
 * @exceptsafe basic
 */
void draw_signals(
    algo::Result& r,
    const std::vector<double>& t,
    const std::vector<double>& yf,
    const std::vector<double>& uf
) {
  canvas::panel_label(r.image, 0, "1  FILTERED ROTATION RATE OF BOTH HANDS");
  const size_t n = std::min(t.size(), std::min(yf.size(), uf.size()));
  if (n < 4) return;
  const size_t stride = std::max<size_t>(1, n / 3000);
  std::vector<double> ts, ys, us;
  for (size_t i = 0; i < n; i += stride) {
    ts.push_back(t[i] - t.front());
    ys.push_back(yf[i]);
    us.push_back(uf[i]);
  }
  std::vector<double> both = ys;
  both.insert(both.end(), us.begin(), us.end());
  canvas::Axes a0 = canvas::panel(0);
  canvas::autoscale(a0, ts, both);
  canvas::set_range(a0, ts.front(), ts.back(), a0.y0, a0.y1);
  canvas::grid_and_ticks(r.image, a0, 10, 5);
  canvas::polyline(r.image, a0, ts, ys, canvas::BLACK);
  canvas::polyline(r.image, a0, ts, us, canvas::VIVID);
  canvas::axis_labels(r.image, a0, "TIME (s)", "RATE (dps)");
  canvas::legend(
      r.image,
      a0,
      {{"OPERATOR", canvas::VIVID}, {"ROBOT", canvas::BLACK}},
      true
  );
}

/**
 * Score every integer delay with its own Kalman filter and keep the best.
 *
 * A bank of `NK_MAX` + 1 hypotheses is built, one per whole millisecond of
 * delay, each carrying an ARX plant fitted at that lag and run as a
 * second-order Kalman filter whose measurement is the robot signal. Every
 * step adds the Gaussian innovation term -0.5*(nu^2/S + log S) to that
 * filter's batch log likelihood, so the winner is the delay under which the
 * robot was least surprising; a parabola through the winner and its two
 * neighbours interpolates below the 1 ms grid, and the same score decayed
 * with a `FORGET_S` time constant gives the leading-delay track drawn in
 * panel 3, whose tail supplies the spread. Modelling the plant rather than
 * assuming a shifted copy is what suits this data, since the robot
 * attenuates and low-passes the operator; the cost is that a flexible plant
 * can partly fake a delay, which flattens the likelihood ridge on gentle
 * takes, and that the preprocessing is shared with `arx_delay` so the two
 * are not independent votes. `Result::note` gives the filter count, the
 * number of update steps and their rate, and the median tracked delay,
 * which is worth comparing against `tau_ms` as a stationarity check.
 * The bank of filters and the likelihood vote between them follow Magill's multiple-
 * model adaptive estimation, IEEE Trans. Automatic Control 1965.
 *
 * @param[in] in  both gyro magnitude channels at 1 kHz plus the sample step
 * @returns the delay in ms, positive when the robot lags the operator, with
 *          `ok` false when the take is too short for the bank or any plant
 *          fit fails
 * @exceptsafe basic
 */
algo::Result kalman_bank(const algo::Input& in) {
  algo::Result r;
  const std::vector<double> h = lowpass_fir(LP_HZ, in.dt);
  const std::vector<double> yf =
      convolve_same(algo::motion_signal(in, algo::ROBOT, 0.3), h);
  const std::vector<double> uf =
      convolve_same(algo::motion_signal(in, algo::OPERATOR, 0.3), h);
  const int n = (int)yf.size();
  const int nd = (n - NK_MAX) / DECIMATE;
  if (nd < 200) {
    r.note = "recording too short for a filter bank";
    reject_tile(r);
    draw_signals(r, in.t, yf, uf);
    return r;
  }

  std::vector<Hypothesis> bank;
  std::vector<double> vars;
  for (int nk = 0; nk <= NK_MAX; nk++) {
    Hypothesis hyp;
    if (!fit_plant(yf, uf, nk, nd, hyp)) {
      r.note = "plant fit failed";
      reject_tile(r);
      draw_signals(r, in.t, yf, uf);
      return r;
    }
    vars.push_back(hyp.var);
    bank.push_back(hyp);
  }
  const double R = algo::median(vars);
  const double Q = Q_OVER_R * R;
  const int start = (NK_MAX + DECIMATE - 1) / DECIMATE + NB;
  for (Hypothesis& hyp : bank) {
    hyp.x0 = yf[start * DECIMATE];
    hyp.x1 = yf[(start - 1) * DECIMATE];
    hyp.p00 = hyp.p11 = R;
    hyp.p01 = hyp.p10 = 0;
    hyp.logw = 0;
    hyp.loglik = 0;
  }

  const double rho = std::exp(-in.dt * DECIMATE / FORGET_S);
  const int steps = nd - start;
  std::vector<double> track_tau;

  for (int t = start; t < nd; t++) {
    const double z = yf[t * DECIMATE];
    for (size_t i = 0; i < bank.size(); i++) {
      Hypothesis& hyp = bank[i];
      double drive = 0;
      for (int j = 0; j < NB; j++)
        drive += hyp.b[j] * uf[(t - j) * DECIMATE - (int)i];
      const double xp0 = -hyp.a1 * hyp.x0 - hyp.a2 * hyp.x1 + drive;
      const double xp1 = hyp.x0;
      const double m00 = hyp.a1 * hyp.a1 * hyp.p00 +
                         hyp.a1 * hyp.a2 * (hyp.p01 + hyp.p10) +
                         hyp.a2 * hyp.a2 * hyp.p11 + Q;
      const double m01 = -hyp.a1 * hyp.p00 - hyp.a2 * hyp.p10;
      const double m10 = -hyp.a1 * hyp.p00 - hyp.a2 * hyp.p01;
      const double m11 = hyp.p00;
      const double S = m00 + R;
      const double nu = z - xp0;
      hyp.loglik += -0.5 * (nu * nu / S + std::log(S));
      hyp.logw = rho * hyp.logw - 0.5 * (nu * nu / S + std::log(S));
      const double k0 = m00 / S, k1 = m10 / S;
      hyp.x0 = xp0 + k0 * nu;
      hyp.x1 = xp1 + k1 * nu;
      hyp.p00 = m00 - k0 * m00;
      hyp.p01 = m01 - k0 * m01;
      hyp.p10 = m10 - k1 * m00;
      hyp.p11 = m11 - k1 * m01;
    }
    int best = 0;
    for (size_t i = 1; i < bank.size(); i++)
      if (bank[i].logw > bank[best].logw) best = (int)i;
    track_tau.push_back(best * in.dt * 1000.0);
  }

  int kbest = 0;
  for (size_t i = 1; i < bank.size(); i++)
    if (bank[i].loglik > bank[kbest].loglik) kbest = (int)i;
  double sub = 0;
  if (kbest > 0 && kbest + 1 < (int)bank.size()) {
    const double a = bank[kbest - 1].loglik, b = bank[kbest].loglik,
                 c = bank[kbest + 1].loglik;
    const double den = a - 2 * b + c;
    if (std::fabs(den) > 1e-18)
      sub = std::max(-1.0, std::min(1.0, 0.5 * (a - c) / den));
  }

  std::vector<double> tail(
      track_tau.begin() + track_tau.size() / 4,
      track_tau.end()
  );

  r.ok = true;
  r.tau_ms = (kbest + sub) * in.dt * 1000.0;
  r.spread_ms = algo::iqr(tail);
  char buf[128];
  snprintf(
      buf,
      sizeof buf,
      "%d filters, %d steps at %d Hz, tracked %.1f ms",
      (int)bank.size(),
      steps,
      int(1.0 / (in.dt * DECIMATE)),
      algo::median(tail)
  );
  r.note = buf;

  r.image = canvas::standard("KALMAN BANK", r.tau_ms, r.ok);
  canvas::caption(r.image, "HYPOTHESIS LOGLIK");
  draw_signals(r, in.t, yf, uf);

  canvas::panel_label(r.image, 1, "2  EVERY DELAY SCORED, THE BEST SCORE WINS");
  std::vector<double> lag_ms(bank.size()), rel(bank.size());
  double lmax = -1e300;
  for (const Hypothesis& hyp : bank) lmax = std::max(lmax, hyp.loglik);
  for (size_t i = 0; i < bank.size(); i++) {
    lag_ms[i] = double(i) * in.dt * 1000.0;
    rel[i] = bank[i].loglik - lmax;
  }
  canvas::Axes a1 = canvas::panel(1);
  canvas::autoscale(a1, lag_ms, rel);
  canvas::grid_and_ticks(r.image, a1, 10, 5);
  canvas::polyline(r.image, a1, lag_ms, rel, canvas::SAGE);
  canvas::mark_guide(r.image, a1, r.tau_ms, true, "ESTIMATE", 0);
  canvas::axis_labels(r.image, a1, "DELAY (ms)", "LOG LIKELIHOOD");
  canvas::legend(r.image, a1, {{"SCORE", canvas::SAGE}}, true);

  canvas::panel_label(
      r.image,
      2,
      "3  WHICH DELAY IS LEADING AS THE RECORDING PLAYS"
  );
  if (track_tau.size() > 1) {
    const size_t stride = std::max<size_t>(1, track_tau.size() / 2000);
    std::vector<double> tt, tv;
    for (size_t i = 0; i < track_tau.size(); i += stride) {
      tt.push_back(double((size_t(start) + i) * DECIMATE) * in.dt);
      tv.push_back(track_tau[i]);
    }
    canvas::Axes a2 = canvas::panel(2);
    canvas::set_range(
        a2,
        tt.front(),
        tt.back(),
        0.0,
        double(NK_MAX) * in.dt * 1000.0
    );
    canvas::grid_and_ticks(r.image, a2, 10, 5);
    canvas::polyline(r.image, a2, tt, tv, canvas::SLATE);
    canvas::mark_guide(r.image, a2, r.tau_ms, false, "ESTIMATE", 0);
    canvas::axis_labels(r.image, a2, "TIME (s)", "DELAY (ms)");
    canvas::legend(r.image, a2, {{"LEADING FILTER", canvas::SLATE}}, true);
  }
  return r;
}

}

ALGO_REGISTER(
    "kalman_bank",
    kalman_bank
);
