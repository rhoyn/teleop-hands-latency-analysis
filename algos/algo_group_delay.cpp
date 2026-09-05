#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdio>
#include <vector>

#include "algo.hpp"

namespace {

/**
 * Welch segmentation and band-selection constants for 1 kHz IMU takes.
 *
 * These trade frequency resolution against the number of segments left to
 * average: the hands are only coherent over a narrow low band, so a shorter
 * segment blurs that band away while a longer one leaves too few segments
 * for the jackknife. The gates were tuned on bursty reach-and-rest takes,
 * where most of the recording is rest and would otherwise contribute
 * incoherent phase at full weight.
 */
constexpr int SEG_LONG = 2048;        ///< 2.05 s window, 0.49 Hz bins
constexpr int SEG_SHORT = 1024;       ///< fallback for a short recording
constexpr int HOP_DIV = 8;            ///< 87.5% overlap, more segments
constexpr double REST_GATE = 0.05;    ///< skip segments under 5% of peak
constexpr double COH_MIN = 0.30;      ///< below this the phase is noise
constexpr double F_LO = 0.5;          ///< below it drift and gait dominate
constexpr double F_HI = 12.0;         ///< above it the robot has nothing left
constexpr int UNWRAP_ITERS = 8;       ///< converges in about four passes
constexpr double PLOT_MAX_HZ = 40.0;  ///< plot extent only, never fitted

/**
 * Welch cross-spectral sums together with the per-segment transforms.
 *
 * `suu`, `syy` and `syu` are running sums over accepted segments rather
 * than averages, which is harmless because every consumer forms a ratio and
 * the missing 1/segments cancels. The individual segment spectra are kept
 * only so the delete-one jackknife can subtract a segment back out of those
 * same sums without recomputing any transform.
 */
struct Spectra {
  std::vector<double> f;
  std::vector<double> suu;
  std::vector<double> syy;
  std::vector<algo::cd> syu;
  std::vector<std::vector<algo::cd>> seg_u;
  std::vector<std::vector<algo::cd>> seg_y;
  int segments = 0;
  int bins = 0;
};

/**
 * Accumulate Hann-windowed cross spectra over the moving part of the take.
 *
 * Motion here is bursty, reaches separated by rest, so segments whose
 * combined energy falls under `REST_GATE` of the loudest segment are
 * dropped: rest is sensor noise whose phase is uniformly random and would
 * be averaged in at full weight. The dense 1/`HOP_DIV` hop buys segment
 * count rather than independent information, which the coherence weighting
 * downstream tolerates because it only ever compares bins to each other.
 *
 * @param[in]  y    robot gyro magnitude, high-passed, the output signal
 * @param[in]  u    operator gyro magnitude, high-passed, the input signal
 * @param[in]  dt   sample interval in seconds, 0.001 for these recordings
 * @param[in]  seg  segment length in samples
 * @param[out] sp   sums and per-segment spectra, left partly filled on a
 *                  false return
 * @returns true once at least four moving segments were accumulated, the
 *          fewest the band fit and the jackknife can work with
 * @exceptsafe basic
 */
bool welch(
    const std::vector<double>& y,
    const std::vector<double>& u,
    double dt,
    int seg,
    Spectra& sp
) {
  const int hop = seg / HOP_DIV;
  const int n = (int)y.size();
  if (n < seg + 4 * hop) return false;
  std::vector<int> starts;
  std::vector<double> energy;
  for (int i = 0; i + seg <= n; i += hop) {
    double e = 0;
    for (int k = 0; k < seg; k++)
      e += y[i + k] * y[i + k] + u[i + k] * u[i + k];
    starts.push_back(i);
    energy.push_back(e);
  }
  double emax = 0;
  for (double e : energy) emax = std::max(emax, e);
  std::vector<double> win(seg);
  for (int k = 0; k < seg; k++)
    win[k] = 0.5 - 0.5 * std::cos(2 * M_PI * double(k) / double(seg));

  sp.bins = seg / 2 + 1;
  sp.f.assign(sp.bins, 0.0);
  sp.suu.assign(sp.bins, 0.0);
  sp.syy.assign(sp.bins, 0.0);
  sp.syu.assign(sp.bins, algo::cd(0.0, 0.0));
  for (int b = 0; b < sp.bins; b++) sp.f[b] = double(b) / (double(seg) * dt);

  for (size_t s = 0; s < starts.size(); s++) {
    if (energy[s] < REST_GATE * emax) continue;
    std::vector<algo::cd> Y(seg), U(seg);
    for (int k = 0; k < seg; k++) {
      Y[k] = y[starts[s] + k] * win[k];
      U[k] = u[starts[s] + k] * win[k];
    }
    algo::fft(Y, false);
    algo::fft(U, false);
    std::vector<algo::cd> yb(sp.bins), ub(sp.bins);
    for (int b = 0; b < sp.bins; b++) {
      yb[b] = Y[b];
      ub[b] = U[b];
      sp.suu[b] += std::norm(U[b]);
      sp.syy[b] += std::norm(Y[b]);
      sp.syu[b] += Y[b] * std::conj(U[b]);
    }
    sp.seg_y.push_back(yb);
    sp.seg_u.push_back(ub);
    sp.segments++;
  }
  return sp.segments >= 4;
}

/**
 * Fit phase against frequency with the line forced through the origin.
 *
 * A free-slope fit collapses on this data: the coherent band is only about
 * one to six Hz wide, and two to five Hz of bandwidth cannot separate a
 * slope from an intercept, so the intercept quietly absorbs the delay. Zero
 * phase at zero frequency is the physically correct constraint for a pure
 * delay anyway. Each pass unwraps `phase` onto the current -omega*tau line
 * and re-solves the weighted slope, which settles well inside
 * `UNWRAP_ITERS` passes.
 *
 * @param[in] f      bin centre frequencies in Hz, ascending
 * @param[in] phase  wrapped phase of the transfer estimate, in radians
 * @param[in] w      per-bin weights, coherence/(1 - coherence)
 * @returns the delay in seconds, positive when the robot lags the operator,
 *          or zero if the weights sum to nothing
 * @exceptsafe basic
 */
double origin_slope(
    const std::vector<double>& f,
    const std::vector<double>& phase,
    const std::vector<double>& w
) {
  double tau = 0;
  std::vector<double> pu(f.size());
  for (int it = 0; it < UNWRAP_ITERS; it++) {
    double num = 0, den = 0;
    for (size_t i = 0; i < f.size(); i++) {
      const double om = 2 * M_PI * f[i];
      pu[i] =
          phase[i] + 2 * M_PI * std::round((-om * tau - phase[i]) / (2 * M_PI));
      num += w[i] * pu[i] * om;
      den += w[i] * om * om;
    }
    if (den <= 0) return 0;
    tau = -num / den;
  }
  return tau;
}

/**
 * Re-run band selection and the slope fit over a supplied set of sums.
 *
 * The jackknife has to apply the whole estimator to leave-one-out sums, so
 * band membership is recomputed rather than reused: dropping one segment
 * can push a marginal bin below `COH_MIN`, and freezing the band would
 * understate the spread that this function exists to measure.
 *
 * @param[in] sp   spectra, read only for its frequency axis and bin count
 * @param[in] suu  operator auto-spectrum sums to use instead of `sp.suu`
 * @param[in] syy  robot auto-spectrum sums to use instead of `sp.syy`
 * @param[in] syu  cross-spectrum sums to use instead of `sp.syu`
 * @returns the delay in seconds, or NaN when fewer than three bins survive
 * @exceptsafe basic
 */
double tau_from_average(
    const Spectra& sp,
    const std::vector<double>& suu,
    const std::vector<double>& syy,
    const std::vector<algo::cd>& syu
) {
  std::vector<double> bf, bp, bw;
  for (int b = 1; b < sp.bins; b++) {
    if (sp.f[b] < F_LO || sp.f[b] > F_HI) continue;
    const double denom = suu[b] * syy[b];
    if (denom <= 0) continue;
    double coh = std::norm(syu[b]) / denom;
    if (coh < COH_MIN) continue;
    coh = std::min(coh, 0.999);
    bf.push_back(sp.f[b]);
    bp.push_back(std::arg(syu[b] / suu[b]));
    bw.push_back(coh / (1.0 - coh));
  }
  if (bf.size() < 3) return std::nan("");
  return origin_slope(bf, bp, bw);
}

/**
 * Draw coherence against frequency with the fitted band bracketed.
 *
 * Shared by the accepted and the abstaining tile, which is why it tolerates
 * empty inputs and simply draws nothing.
 *
 * @param[in,out] im   tile being drawn into
 * @param[in]     idx  panel index to fill
 * @param[in]     f    frequency axis in Hz
 * @param[in]     coh  magnitude-squared coherence at those frequencies
 * @exceptsafe basic
 */
void draw_coherence(
    canvas::Image& im,
    int idx,
    const std::vector<double>& f,
    const std::vector<double>& coh
) {
  if (f.size() < 2) return;
  canvas::Axes ax = canvas::panel(idx);
  canvas::set_range(ax, 0.0, f.back(), 0.0, 1.0);
  canvas::grid_and_ticks(im, ax, 10, 5);
  canvas::mark_guide(im, ax, F_LO, true, "FITTED BAND", 0);
  canvas::mark_guide(im, ax, F_HI, true, "FITTED BAND", 1);
  canvas::mark_guide(im, ax, COH_MIN, false, "MINIMUM", 0);
  canvas::polyline(im, ax, f, coh, canvas::SLATE);
  canvas::axis_labels(im, ax, "FREQUENCY (Hz)", "COHERENCE");
  canvas::legend(im, ax, {{"COHERENCE", canvas::SLATE}}, true);
}

/**
 * Build the abstaining result, still showing why no fit was attempted.
 *
 * Abstention is a legitimate outcome rather than an error: ten algorithms
 * vote and the median is reported, so withholding a vote beats fitting a
 * slope to noise. The tile keeps the coherence curve so a reader can see
 * for themselves that no band cleared `COH_MIN`.
 *
 * @param[in] note  short reason, copied verbatim into `Result::note`
 * @param[in] f     frequency axis for the coherence panel, may be empty
 * @param[in] coh   coherence at those frequencies, may be empty
 * @returns a `Result` with `ok` false and `tau_ms` zero
 * @exceptsafe basic
 */
algo::Result rejected(
    const char* note,
    const std::vector<double>& f,
    const std::vector<double>& coh
) {
  algo::Result r;
  r.note = note;
  r.image = canvas::standard("GROUP DELAY", 0.0, false);
  canvas::caption(r.image, "NO COHERENT BAND");
  canvas::panel_label(
      r.image,
      0,
      "1  THE TWO HANDS SHARE NO USABLE FREQUENCY BAND"
  );
  draw_coherence(r.image, 0, f, coh);
  return r;
}

/**
 * Estimate the delay from the slope of the operator-to-robot phase.
 *
 * Welch averaging over the moving segments gives the H1 transfer estimate
 * H(f) = Syu(f)/Suu(f), and a pure delay tau contributes arg H =
 * -2*pi*f*tau, so the delay is the slope of the unwrapped phase against
 * frequency. That slope is fitted through the origin and weighted by
 * coherence/(1 - coherence), which is the usual variance weighting for a
 * phase estimate and makes bins the hands disagree on nearly free. Working
 * from phase alone is the point of the method on this data: the robot
 * attenuates and low-passes the operator instead of copying it, so any
 * magnitude-based estimator inherits that gain error while this one does
 * not. What it does struggle with is bandwidth, since the pair stays
 * coherent only over roughly one to six Hz, which is exactly why the fit is
 * constrained through the origin; the reported spread is the standard error
 * of a delete-one-segment jackknife over the same estimator, and
 * `Result::note` gives the segment count and segment length in ms, the
 * number of bins that entered the fit and the frequency span they cover.
 * The H1 estimator and the coherence weighting follow Bendat and Piersol, Random Data:
 * Analysis and Measurement Procedures; the origin constraint on the phase fit is a
 * departure from that treatment, forced by the narrow coherent band here.
 *
 * @param[in] in  both gyro magnitude channels at 1 kHz plus the sample step
 * @returns the delay in ms, positive when the robot lags the operator, with
 *          `ok` false when the take was too short for Welch averaging or no
 *          band cleared `COH_MIN`
 * @exceptsafe basic
 */
algo::Result group_delay(const algo::Input& in) {
  algo::Result r;
  const std::vector<double> y = algo::motion_signal(in, algo::ROBOT, 0.3);
  const std::vector<double> u = algo::motion_signal(in, algo::OPERATOR, 0.3);

  Spectra sp;
  int seg = SEG_LONG;
  if (!welch(y, u, in.dt, seg, sp)) {
    sp = Spectra();
    seg = SEG_SHORT;
    if (!welch(y, u, in.dt, seg, sp))
      return rejected("recording too short for Welch averaging", {}, {});
  }

  std::vector<double> plot_f, plot_coh, plot_gain_db;
  for (int b = 1; b < sp.bins; b++) {
    if (sp.f[b] > PLOT_MAX_HZ) break;
    const double denom = std::max(sp.suu[b] * sp.syy[b], 1e-300);
    plot_f.push_back(sp.f[b]);
    plot_coh.push_back(std::min(1.0, std::norm(sp.syu[b]) / denom));
    const double gain = std::abs(sp.syu[b]) / std::max(sp.suu[b], 1e-300);
    plot_gain_db.push_back(
        std::max(-60.0, std::min(40.0, 20.0 * std::log10(gain + 1e-12)))
    );
  }

  std::vector<double> bf, bp, bw;
  for (int b = 1; b < sp.bins; b++) {
    const double coh_raw =
        std::norm(sp.syu[b]) / std::max(sp.suu[b] * sp.syy[b], 1e-300);
    if (sp.f[b] < F_LO || sp.f[b] > F_HI) continue;
    if (coh_raw < COH_MIN) continue;
    const double coh = std::min(coh_raw, 0.999);
    bf.push_back(sp.f[b]);
    bp.push_back(std::arg(sp.syu[b] / std::max(sp.suu[b], 1e-300)));
    bw.push_back(coh / (1.0 - coh));
  }
  if (bf.size() < 3)
    return rejected(
        "no frequency band with usable coherence",
        plot_f,
        plot_coh
    );

  const double tau = origin_slope(bf, bp, bw);

  std::vector<double> jack;
  if (sp.segments >= 5) {
    for (int s = 0; s < sp.segments; s++) {
      std::vector<double> suu(sp.bins, 0.0), syy(sp.bins, 0.0);
      std::vector<algo::cd> syu(sp.bins, algo::cd(0.0, 0.0));
      for (int b = 0; b < sp.bins; b++) {
        suu[b] = sp.suu[b] - std::norm(sp.seg_u[s][b]);
        syy[b] = sp.syy[b] - std::norm(sp.seg_y[s][b]);
        syu[b] = sp.syu[b] - sp.seg_y[s][b] * std::conj(sp.seg_u[s][b]);
      }
      const double t = tau_from_average(sp, suu, syy, syu);
      if (std::isfinite(t)) jack.push_back(t * 1000.0);
    }
  }
  double spread = 0;
  if (jack.size() >= 5) {
    const double m = algo::mean(jack);
    double s2 = 0;
    for (double v : jack) s2 += (v - m) * (v - m);
    spread = std::sqrt(s2 * double(jack.size() - 1) / double(jack.size()));
  }

  r.ok = true;
  r.tau_ms = tau * 1000.0;
  r.spread_ms = spread;
  char buf[128];
  snprintf(
      buf,
      sizeof buf,
      "%d segs of %d ms, %d bins, %.1f-%.1f Hz",
      sp.segments,
      int(seg * in.dt * 1000.0),
      (int)bf.size(),
      bf.front(),
      bf.back()
  );
  r.note = buf;

  std::vector<double> unwrapped(bf.size()), fitted(bf.size());
  for (size_t i = 0; i < bf.size(); i++) {
    const double om = 2 * M_PI * bf[i];
    unwrapped[i] =
        bp[i] + 2 * M_PI * std::round((-om * tau - bp[i]) / (2 * M_PI));
    fitted[i] = -om * tau;
  }
  std::vector<double> all_x(bf), all_y(unwrapped);
  all_x.insert(all_x.end(), bf.begin(), bf.end());
  all_y.insert(all_y.end(), fitted.begin(), fitted.end());
  r.image = canvas::standard("GROUP DELAY", r.tau_ms, r.ok);
  canvas::caption(r.image, "PHASE SLOPE FIT");

  canvas::panel_label(
      r.image,
      0,
      "1  HOW WELL THE HANDS AGREE AT EACH FREQUENCY"
  );
  draw_coherence(r.image, 0, plot_f, plot_coh);

  canvas::panel_label(
      r.image,
      1,
      "2  PHASE FALLS WITH FREQUENCY, THE SLOPE IS THE DELAY"
  );
  {
    canvas::Axes a1 = canvas::panel(1);
    canvas::autoscale(a1, all_x, all_y);
    canvas::grid_and_ticks(r.image, a1, 10, 5);
    canvas::mark_guide(r.image, a1, 0.0, false, "ZERO PHASE", 0);
    canvas::polyline(r.image, a1, bf, fitted, canvas::BLACK);
    canvas::polyline(r.image, a1, bf, unwrapped, canvas::VIVID);
    for (size_t i = 0; i < bf.size(); i++)
      canvas::marker(r.image, a1, bf[i], unwrapped[i], canvas::VIVID, 4);
    canvas::axis_labels(r.image, a1, "FREQUENCY (Hz)", "PHASE (rad)");
    canvas::legend(
        r.image,
        a1,
        {{"MEASURED", canvas::VIVID, true}, {"FITTED SLOPE", canvas::BLACK}},
        true
    );
  }

  canvas::panel_label(
      r.image,
      2,
      "3  GAIN FOR CONTEXT, THIS METHOD ONLY USES PHASE"
  );
  if (plot_f.size() > 1) {
    canvas::Axes a2 = canvas::panel(2);
    canvas::autoscale(a2, plot_f, plot_gain_db);
    canvas::set_range(a2, 0.0, plot_f.back(), a2.y0, a2.y1);
    canvas::grid_and_ticks(r.image, a2, 10, 5);
    canvas::mark_guide(r.image, a2, F_LO, true, "FITTED BAND", 0);
    canvas::mark_guide(r.image, a2, F_HI, true, "FITTED BAND", 1);
    canvas::mark_guide(r.image, a2, 0.0, false, "UNITY GAIN", 0);
    canvas::polyline(r.image, a2, plot_f, plot_gain_db, canvas::SHADE);
    canvas::axis_labels(r.image, a2, "FREQUENCY (Hz)", "GAIN (dB)");
    canvas::legend(r.image, a2, {{"GAIN", canvas::SHADE}}, true);
  }
  return r;
}

}

ALGO_REGISTER(
    "group_delay",
    group_delay
);
