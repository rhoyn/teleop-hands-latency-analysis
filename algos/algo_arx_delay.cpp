#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

#include <Eigen/Dense>

#include "algo.hpp"

namespace {

/**
 * Tuning for the low-pass, the decimation and the delay grid.
 *
 * One-step prediction on the raw 1 kHz pair is dominated by the band above
 * roughly 10 Hz, where the two hands share no coherent content at all, so
 * the signals are low-passed and decimated before any model is fitted;
 * without that the residual is almost pure noise and the loss curve over
 * the delay grid comes out flat. The model orders are kept deliberately
 * small so that the delay, rather than a rich set of poles, has to explain
 * what is left.
 */
constexpr int DECIMATE = 2;     ///< 1 kHz down to 500 Hz, ample above LP_HZ
constexpr double LP_HZ = 25.0;  ///< above this the hands share no motion
constexpr int LP_TAPS = 101;    ///< odd, so the kernel stays linear phase
constexpr int NA = 2;           ///< two poles cover the robot's rolloff
constexpr int NB = 3;           ///< three input taps, a short shaping FIR
constexpr int NK_MAX = 200;     ///< 200 raw samples is the 200 ms ceiling
constexpr int BLOCKS = 4;       ///< four sub-ranges give an IQR to report

/**
 * Build a Hamming-windowed sinc low-pass kernel.
 *
 * The taps are normalised to unit DC gain so the filter neither rescales
 * the rotation rate nor shifts its mean, and `LP_TAPS` is odd, which keeps
 * the kernel symmetric and therefore linear phase. Phase matters far more
 * here than stopband depth: a filter that delayed one hand differently from
 * the other would bias the very quantity being measured.
 *
 * @param[in] fc  cutoff in Hz
 * @param[in] dt  sample period in seconds, 0.001 for the 1 kHz IMUs
 * @returns `LP_TAPS` coefficients summing to one
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
 * Convolve a signal with a kernel and keep the input length.
 *
 * The kernel is centred, so a symmetric `h` adds no group delay and both
 * hands come out shifted identically. Samples that would need data from
 * before the start or past the end are treated as zero, which leaves about
 * `h.size() / 2` samples attenuated at each edge; the delay grid discards
 * considerably more than that at the front anyway.
 *
 * @param[in] x  signal to filter
 * @param[in] h  kernel, assumed symmetric
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
 * One least-squares ARX fit at a single delay hypothesis.
 *
 * `loss` keeps its sentinel value when the hypothesis was never fitted,
 * which is how the caller distinguishes an unusable delay from a merely
 * poor one without a separate flag.
 */
struct Fit {
  Eigen::VectorXd beta;  ///< NA denominator then NB numerator coefficients
  double loss = 1e300;   ///< mean squared one-step residual, sentinel if unfit
  int rows = 0;          ///< regression rows used, 0 when the fit was skipped
};

/**
 * Fit an ARX model for one delay hypothesis over a range of samples.
 *
 * The regression predicts the robot from `NA` of its own past samples and
 * `NB` operator samples starting `nk` back, so the delay is carried by the
 * regressor indices rather than by any coefficient. `nk` counts raw 1 kHz
 * samples while the row index steps in decimated ones, which buys 1 ms
 * delay resolution out of a model that only runs at 500 Hz. A column-pivoted
 * QR is used because the operator regressors are strongly collinear once
 * the signal has been low-passed to 25 Hz.
 *
 * @param[in] yf    filtered robot signal at the full 1 kHz rate
 * @param[in] uf    filtered operator signal at the full 1 kHz rate
 * @param[in] nk    delay hypothesis, in raw samples
 * @param[in] from  first decimated row index of the range
 * @param[in] to    one past the last decimated row index
 * @returns the fit, or a default `Fit` whose `rows` is zero and whose
 *          `loss` is still the sentinel when the range holds too few rows
 *          to identify `NA + NB` coefficients
 * @exceptsafe basic
 */
Fit arx_fit(
    const std::vector<double>& yf,
    const std::vector<double>& uf,
    int nk,
    int from,
    int to
) {
  Fit fit;
  const int p = std::max(NA, NB - 1) + 1;
  const int i0 = std::max(from + p, (nk + DECIMATE - 1) / DECIMATE + p);
  const int i1 = to;
  const int rows = i1 - i0;
  if (rows < 4 * (NA + NB)) return fit;
  Eigen::MatrixXd X(rows, NA + NB);
  Eigen::VectorXd Y(rows);
  for (int r = 0; r < rows; r++) {
    const int t = i0 + r;
    for (int i = 0; i < NA; i++) X(r, i) = -yf[(t - 1 - i) * DECIMATE];
    for (int j = 0; j < NB; j++) X(r, NA + j) = uf[(t - j) * DECIMATE - nk];
    Y(r) = yf[t * DECIMATE];
  }
  fit.beta = X.colPivHouseholderQr().solve(Y);
  fit.loss = (Y - X * fit.beta).squaredNorm() / double(rows);
  fit.rows = rows;
  return fit;
}

/**
 * Locate the minimum of a sampled curve to sub-sample precision.
 *
 * A parabola through the minimum and its two neighbours gives a fractional
 * index, which is needed because the grid step is a whole millisecond while
 * the delays being measured are not multiples of one. The correction is
 * clamped to plus or minus one sample so a flat or noisy neighbourhood
 * cannot throw the estimate outside its bracket, and a minimum sitting on
 * either endpoint falls back to the integer index.
 *
 * @param[in]  v      loss per delay hypothesis
 * @param[out] k_out  index of the integer minimum, ignored when null
 * @returns the interpolated index of the minimum
 * @exceptsafe no-throw
 */
double refine_min(
    const std::vector<double>& v,
    int* k_out
) {
  int k = 0;
  for (size_t i = 1; i < v.size(); i++)
    if (v[i] < v[k]) k = (int)i;
  if (k_out) *k_out = k;
  if (k <= 0 || k + 1 >= (int)v.size()) return double(k);
  const double a = v[k - 1], b = v[k], c = v[k + 1];
  const double den = a - 2 * b + c;
  if (std::fabs(den) < 1e-18) return double(k);
  const double sub = 0.5 * (a - c) / den;
  return double(k) + std::max(-1.0, std::min(1.0, sub));
}

/**
 * Draw both filtered hands against time into one panel.
 *
 * The traces are decimated to at most about 2000 points each, which is all
 * an 800 px wide panel can resolve, and share one autoscaled y axis so the
 * robot's attenuation relative to the operator stays visible.
 *
 * @param[in,out] im   tile being drawn into
 * @param[in]     idx  panel index, 0 to 2 running top to bottom
 * @param[in]     t    timestamps in seconds
 * @param[in]     y    filtered robot trace
 * @param[in]     u    filtered operator trace
 * @exceptsafe basic
 */
void plot_pair(
    canvas::Image& im,
    int idx,
    const std::vector<double>& t,
    const std::vector<double>& y,
    const std::vector<double>& u
) {
  const size_t n = std::min(t.size(), std::min(y.size(), u.size()));
  if (n < 4) return;
  const size_t stride = std::max<size_t>(1, n / 2000);
  std::vector<double> ts, ys, us;
  for (size_t i = 0; i < n; i += stride) {
    ts.push_back(t[i]);
    ys.push_back(y[i]);
    us.push_back(u[i]);
  }
  std::vector<double> both(ys);
  both.insert(both.end(), us.begin(), us.end());
  canvas::Axes ax = canvas::panel(idx);
  canvas::autoscale(ax, ts, both);
  canvas::set_range(ax, ts.front(), ts.back(), ax.y0, ax.y1);
  canvas::grid_and_ticks(im, ax, 10, 5);
  canvas::polyline(im, ax, ts, ys, canvas::BLACK);
  canvas::polyline(im, ax, ts, us, canvas::VIVID);
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
 * latency is the median over the ten algorithms, and it simply forms
 * without this one. The tile is still drawn, with only the first panel
 * filled, so a reader can see what the estimator was handed.
 *
 * @param[in] note  short reason recorded in `Result::note`
 * @param[in] t     timestamps in seconds
 * @param[in] yf    filtered robot trace
 * @param[in] uf    filtered operator trace
 * @returns a `Result` with `ok` false and a single populated panel
 * @exceptsafe basic
 */
algo::Result rejected(
    const char* note,
    const std::vector<double>& t,
    const std::vector<double>& yf,
    const std::vector<double>& uf
) {
  algo::Result r;
  r.note = note;
  r.image = canvas::standard("ARX DELAY", 0.0, false);
  canvas::caption(r.image, "NO DELAY GRID");
  canvas::panel_label(
      r.image,
      0,
      "1  SMOOTHED ROTATION RATE OF BOTH HANDS OVER TIME"
  );
  plot_pair(r.image, 0, t, yf, uf);
  return r;
}

/**
 * Estimate the teleoperation delay by fitting an ARX model per hypothesis.
 *
 * Both gyro magnitudes are high-passed at 0.3 Hz to strip gravity and slow
 * bias, low-passed at `LP_HZ` and decimated to 500 Hz. For every delay nk
 * from 0 to `NK_MAX` samples the model
 * y[t] = -sum_i a_i y[t-i] + sum_j b_j u[t-nk-j] is fitted by least squares
 * and its mean squared one-step residual recorded; the nk minimising that
 * loss is the delay, refined by a parabola through the trough. Because the
 * model carries its own poles it can absorb the robot's attenuation and
 * low-pass rolloff instead of charging them to the delay, which is exactly
 * what a plain correlation cannot do. The failure mode is the mirror of
 * that strength: a too-rich model trades delay against pole placement, so
 * `NA` and `NB` stay small. Rest hurts too, since a silent stretch fits
 * equally well at every hypothesis and flattens the loss curve.
 *
 * `Result::note` reports the model orders, the rate the model was fitted
 * at, how far the loss falls from its worst hypothesis to its best — a
 * shallow drop means a poorly determined delay — and how many of the
 * `BLOCKS` sub-ranges yielded their own estimate; `Result::spread_ms` is
 * the IQR across those block estimates. The three panels show the filtered
 * pair, the loss curve over the delay grid with the grid minimum marked,
 * and a free-run simulation of the fitted model driven by the operator
 * alone. That last panel is not a one-step prediction: the model is fed its
 * own past outputs, so it shows whether the fit reproduces the robot from
 * scratch, and the run is abandoned as soon as it diverges past eight times
 * the robot's peak amplitude.
 * The ARX formulation with an explicit delay and the search over the loss curve follow
 * Ljung, System Identification: Theory for the User.
 *
 * @param[in] in  both 1 kHz gyro magnitudes with their shared timebase
 * @returns the delay in ms, positive when the robot lags the operator, or
 *          an abstaining `Result` when the recording is too short for the
 *          grid or every hypothesis failed to solve
 * @exceptsafe basic
 */
algo::Result arx_delay(const algo::Input& in) {
  algo::Result r;
  const std::vector<double> h = lowpass_fir(LP_HZ, in.dt);
  const std::vector<double> yf =
      convolve_same(algo::motion_signal(in, algo::ROBOT, 0.3), h);
  const std::vector<double> uf =
      convolve_same(algo::motion_signal(in, algo::OPERATOR, 0.3), h);
  const int n = (int)yf.size();
  const int nd = (n - NK_MAX) / DECIMATE;
  if (nd < 200)
    return rejected(
        "recording too short for a delay grid search",
        in.t,
        yf,
        uf
    );

  std::vector<double> loss;
  for (int nk = 0; nk <= NK_MAX; nk++)
    loss.push_back(arx_fit(yf, uf, nk, 0, nd).loss);

  int kmin = 0;
  const double khat = refine_min(loss, &kmin);
  if (loss[kmin] > 1e299)
    return rejected(
        "least squares failed on every delay hypothesis",
        in.t,
        yf,
        uf
    );

  const int span = nd / BLOCKS;
  std::vector<double> block_taus;
  for (int b = 0; b < BLOCKS; b++) {
    std::vector<double> bl;
    for (int nk = 0; nk <= NK_MAX; nk++)
      bl.push_back(arx_fit(yf, uf, nk, b * span, (b + 1) * span).loss);
    if (bl[0] > 1e299) continue;
    int bk = 0;
    const double bh = refine_min(bl, &bk);
    block_taus.push_back(bh * in.dt * 1000.0);
  }

  double vmin = loss[kmin], vmax = 0;
  for (double v : loss) vmax = std::max(vmax, v);

  r.ok = true;
  r.tau_ms = khat * in.dt * 1000.0;
  r.spread_ms = block_taus.size() >= 3 ? algo::iqr(block_taus) : 0.0;
  char buf[128];
  snprintf(
      buf,
      sizeof buf,
      "ARX(%d,%d) %dHz, loss drop %.0f%%, %d blocks",
      NA,
      NB,
      int(1.0 / (in.dt * DECIMATE)),
      100.0 * (1.0 - vmin / vmax),
      (int)block_taus.size()
  );
  r.note = buf;

  std::vector<double> nk_ms, nk_loss;
  for (int nk = 0; nk <= NK_MAX; nk++) {
    if (loss[size_t(nk)] > 1e299) continue;
    nk_ms.push_back(double(nk) * in.dt * 1000.0);
    nk_loss.push_back(loss[size_t(nk)]);
  }
  const Fit best = arx_fit(yf, uf, kmin, 0, nd);
  std::vector<double> pred_t, pred_y, pred_hat;
  if (best.rows > 0) {
    const int p = std::max(NA, NB - 1) + 1;
    const int i0 = std::max(p, (kmin + DECIMATE - 1) / DECIMATE + p);
    const int from = 0;
    const int to = best.rows;
    double scale = 0;
    for (int rr = from; rr < to; rr++)
      scale = std::max(scale, std::fabs(yf[size_t((i0 + rr) * DECIMATE)]));
    std::vector<double> sim;
    for (int rr = from; rr < to; rr++) {
      const int t = i0 + rr;
      double s = 0;
      for (int i = 0; i < NA; i++) {
        const int k = int(sim.size()) - 1 - i;
        const double past =
            (k >= 0) ? sim[size_t(k)] : yf[size_t((t - 1 - i) * DECIMATE)];
        s += -past * best.beta(i);
      }
      for (int j = 0; j < NB; j++)
        s += uf[size_t((t - j) * DECIMATE - kmin)] * best.beta(NA + j);
      if (!std::isfinite(s) || std::fabs(s) > 8 * scale) break;
      sim.push_back(s);
      const size_t si = size_t(t * DECIMATE);
      pred_t.push_back((si < in.t.size()) ? in.t[si] : double(si) * in.dt);
      pred_y.push_back(yf[si]);
      pred_hat.push_back(s);
    }
  }

  r.image = canvas::standard("ARX DELAY", r.tau_ms, r.ok);
  canvas::caption(r.image, "ARX LOSS VS DELAY");

  canvas::panel_label(
      r.image,
      0,
      "1  SMOOTHED ROTATION RATE OF BOTH HANDS OVER TIME"
  );
  plot_pair(r.image, 0, in.t, yf, uf);

  canvas::panel_label(
      r.image,
      1,
      "2  MODEL FIT ERROR PER ASSUMED DELAY, LOWEST WINS"
  );
  if (nk_ms.size() > 1) {
    canvas::Axes a1 = canvas::panel(1);
    canvas::autoscale(a1, nk_ms, nk_loss);
    canvas::grid_and_ticks(r.image, a1, 10, 5);
    canvas::polyline(r.image, a1, nk_ms, nk_loss, canvas::SHADE);
    canvas::marker(
        r.image,
        a1,
        double(kmin) * in.dt * 1000.0,
        vmin,
        canvas::WHITE,
        3
    );
    canvas::mark_guide(r.image, a1, r.tau_ms, true, "ESTIMATE", 0);
    canvas::axis_labels(r.image, a1, "ASSUMED DELAY (ms)", "LOSS");
    canvas::legend(
        r.image,
        a1,
        {{"FIT LOSS", canvas::SHADE}, {"GRID MINIMUM", canvas::WHITE, true}},
        true
    );
  }

  canvas::panel_label(
      r.image,
      2,
      "3  MODEL DRIVEN BY THE OPERATOR ALONE VS THE REAL ROBOT"
  );
  if (pred_t.size() > 1) {
    canvas::Axes a2 = canvas::panel(2);
    std::vector<double> both(pred_y);
    both.insert(both.end(), pred_hat.begin(), pred_hat.end());
    canvas::autoscale(a2, pred_t, both);
    canvas::set_range(a2, pred_t.front(), pred_t.back(), a2.y0, a2.y1);
    canvas::grid_and_ticks(r.image, a2, 10, 5);
    canvas::polyline(r.image, a2, pred_t, pred_y, canvas::BLACK);
    canvas::polyline(r.image, a2, pred_t, pred_hat, canvas::VIVID);
    canvas::axis_labels(r.image, a2, "TIME (s)", "RATE (dps)");
    canvas::legend(
        r.image,
        a2,
        {{"MODEL OUTPUT", canvas::VIVID}, {"ROBOT", canvas::BLACK}},
        true
    );
  }
  return r;
}

}

ALGO_REGISTER(
    "arx_delay",
    arx_delay
);
