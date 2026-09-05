#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdio>
#include <cstring>
#include <vector>

#include "algo.hpp"

namespace lat {
using cd = std::complex<double>;

/**
 * Everything one estimation run produced, both the answer and the plots.
 *
 * This struct predates `algo::Result` and is kept because the estimator
 * below is the original standalone tool, imported whole; the wrapper at the
 * bottom of the file translates it. `tau_ms` is the median over accepted
 * windows and `tau_full_ms` the single whole-recording peak, kept separately
 * so the two can be compared, and the plotting vectors are carried out of
 * the estimator rather than recomputed because the panels must show exactly
 * the data the estimate was made from.
 */
struct Result {
  bool ok = false;
  char msg[192] = "";
  double tau_ms = 0;
  double corr = 0;
  double iqr_ms = 0;
  double tau_full_ms = 0;
  int nwin = 0;
  std::vector<double> cc_lag_ms, cc_val;
  std::vector<double> ov_t;
  std::vector<double> ov_ahp, ov_bhp;
  std::vector<double> win_t, win_tau;
  std::vector<int> win_pass;
};

/**
 * Round up to the next power of two.
 *
 * Named apart from `algo::next_pow2` because this file predates the shared
 * header and still carries its own copy of the transform code.
 *
 * @param[in] n  length to round up
 * @returns the smallest power of two not less than `n`
 * @exceptsafe no-throw
 */
static size_t lat_next_pow2(size_t n) {
  size_t p = 1;
  while (p < n) p <<= 1;
  return p;
}

/**
 * Transform in place with an iterative radix-2 Cooley-Tukey pass.
 *
 * The size must already be a power of two; the bit-reversal permutation
 * runs first and the butterflies then work in place. The inverse scales by
 * 1/n, so a forward pass followed by an inverse round-trips exactly up to
 * rounding.
 *
 * @param[in,out] a    samples, transformed in place, size a power of two
 * @param[in]     inv  true for the inverse transform
 * @exceptsafe basic
 */
static void fft(
    std::vector<cd>& a,
    bool inv
) {
  size_t n = a.size();
  for (size_t i = 1, j = 0; i < n; i++) {
    size_t bit = n >> 1;
    for (; j & bit; bit >>= 1) j ^= bit;
    j ^= bit;
    if (i < j) std::swap(a[i], a[j]);
  }
  for (size_t len = 2; len <= n; len <<= 1) {
    double ang = 2 * M_PI / (double)len * (inv ? 1 : -1);
    cd wlen(std::cos(ang), std::sin(ang));
    for (size_t i = 0; i < n; i += len) {
      cd w(1);
      for (size_t k = 0; k < len / 2; k++) {
        cd u = a[i + k], v = a[i + k + len / 2] * w;
        a[i + k] = u + v;
        a[i + k + len / 2] = u - v;
        w *= wlen;
      }
    }
  }
  if (inv) {
    for (auto& x : a) x /= (double)n;
  }
}

/**
 * A raw circular cross-correlation and the norms needed to scale it.
 *
 * `z` is unnormalised and `nf` samples long, with negative lags wrapped
 * around to the top end, which is why every read goes through a helper
 * rather than indexing it directly. The two energies are kept so a lag can
 * be turned into a correlation coefficient without a second pass.
 */
struct XC {
  std::vector<double> z;
  size_t nf;
  double nx, ny;
};
/**
 * Cross-correlate two equal-length signals by FFT, means removed.
 *
 * Both signals are mean-subtracted first so a residual DC offset cannot
 * dominate the peak, then zero padded to at least twice their length so the
 * circular correlation the FFT computes agrees with the linear one over
 * every lag that will actually be read. The sign convention matters and is
 * set here: `z[m]` sums `x[i + m] * y[i]`, so a positive lag means `x`
 * carries the event later than `y`, and the caller passes the robot as `x`
 * and the operator as `y` to make a positive delay mean the robot lags.
 *
 * @param[in] x  first signal, the robot in every call here
 * @param[in] y  second signal, the operator in every call here
 * @returns the unnormalised correlation with both signal energies
 * @exceptsafe basic
 */
static XC xcorr(
    const std::vector<double>& x,
    const std::vector<double>& y
) {
  size_t n = x.size();
  size_t nf = lat_next_pow2(2 * n);
  std::vector<cd> X(nf, 0), Y(nf, 0);
  double mx = 0, my = 0;
  for (double v : x) mx += v;
  for (double v : y) my += v;
  mx /= (double)n;
  my /= (double)n;
  double nx = 0, ny = 0;
  for (size_t i = 0; i < n; i++) {
    double a = x[i] - mx, b = y[i] - my;
    X[i] = a;
    Y[i] = b;
    nx += a * a;
    ny += b * b;
  }
  fft(X, false);
  fft(Y, false);
  for (size_t i = 0; i < nf; i++) X[i] *= std::conj(Y[i]);
  fft(X, true);
  XC r;
  r.nf = nf;
  r.nx = nx;
  r.ny = ny;
  r.z.resize(nf);
  for (size_t i = 0; i < nf; i++) r.z[i] = X[i].real();
  return r;
}

/**
 * One correlation peak: its refined lag, its height and its integer bin.
 *
 * `lag_samp` carries the sub-sample refinement while `best_lag` is the raw
 * argmax, kept because the caller checks it against the search limits.
 */
struct Est {
  double lag_samp;
  double corr;
  long best_lag;
};
/**
 * Find the correlation peak within a lag limit and refine it sub-sample.
 *
 * At 1 kHz a whole-sample argmax quantises the answer to 1 ms, a visible
 * share of the 17-21 ms seen on the quickest rigs, so the peak is refined
 * against its two neighbours. The fit is done on logs when all three
 * samples are positive, which finds the centre of a Gaussian-shaped peak
 * exactly and suits a correlation lobe far better than a parabola; the
 * plain parabola is the fallback for a peak sitting on negative values. An
 * offset outside one sample means the vertex is not really here and is
 * discarded rather than trusted.
 *
 * @param[in] xc      correlation to search, with negative lags wrapped
 * @param[in] maxlag  largest magnitude of lag to consider, in samples
 * @returns the refined lag in samples, the correlation coefficient at the
 *          peak, and the integer bin it came from
 * @exceptsafe basic
 */
static Est peak(
    const XC& xc,
    long maxlag
) {
  double norm = std::sqrt(xc.nx * xc.ny) + 1e-30;
  auto val = [&](long L) {
    size_t idx = (L >= 0) ? (size_t)L : (size_t)((long)xc.nf + L);
    return xc.z[idx];
  };
  long best = 0;
  double bestv = -1e300;
  for (long L = -maxlag; L <= maxlag; L++) {
    double v = val(L);
    if (v > bestv) {
      bestv = v;
      best = L;
    }
  }
  double y0 = val(best - 1), y1 = val(best), y2 = val(best + 1), off = 0;
  if (y0 > 0 && y1 > 0 && y2 > 0) {
    double l0 = std::log(y0), l1 = std::log(y1), l2 = std::log(y2);
    double den = l0 - 2 * l1 + l2;
    if (std::fabs(den) > 1e-12) off = 0.5 * (l0 - l2) / den;
  } else {
    double den = y0 - 2 * y1 + y2;
    if (std::fabs(den) > 1e-12) off = 0.5 * (y0 - y2) / den;
  }
  if (off > 1 || off < -1) off = 0;
  Est e;
  e.lag_samp = (double)best + off;
  e.corr = bestv / norm;
  e.best_lag = best;
  return e;
}

/**
 * Linearly interpolate a signal onto a uniform time grid.
 *
 * The two IMUs are nominally 1 kHz but arrive on their own timestamps, and
 * cross-correlation assumes a common uniform grid, so both are resampled
 * onto one before anything else happens. The scan pointer only moves
 * forwards, which makes this linear but also means `t` must be ascending.
 * Past the end of the input the last sample is held rather than
 * extrapolated, so a short channel produces a flat tail instead of a ramp.
 *
 * @param[in] t   source timestamps in seconds, ascending
 * @param[in] y   source samples, paired with `t` by index
 * @param[in] t0  first output time in seconds
 * @param[in] dt  output sample interval in seconds
 * @param[in] n   number of output samples
 * @returns the resampled signal, `n` samples long
 * @exceptsafe basic
 */
static std::vector<double> resample(
    const std::vector<double>& t,
    const std::vector<double>& y,
    double t0,
    double dt,
    size_t n
) {
  std::vector<double> out(n);
  size_t j = 0;
  for (size_t i = 0; i < n; i++) {
    double tt = t0 + (double)i * dt;
    while (j + 1 < t.size() && t[j + 1] < tt) j++;
    if (j + 1 >= t.size()) {
      out[i] = y.empty() ? 0.0 : y.back();
      continue;
    }
    double f = (tt - t[j]) / (t[j + 1] - t[j] + 1e-30);
    out[i] = y[j] + f * (y[j + 1] - y[j]);
  }
  return out;
}

/**
 * One-pole high pass in place, forwards only.
 *
 * Duplicated from `algo::highpass` because this file is the imported
 * original; the coefficient and the state update are the same, so changing
 * one without the other would silently split the two.
 *
 * @param[in,out] x   signal filtered in place
 * @param[in]     fc  corner frequency in Hz
 * @param[in]     dt  sample interval in seconds
 * @exceptsafe basic
 */
static void highpass(
    std::vector<double>& x,
    double fc,
    double dt
) {
  if (x.empty()) return;
  double a = std::exp(-2 * M_PI * fc * dt);
  double px = x[0], py = 0;
  for (size_t i = 0; i < x.size(); i++) {
    double xi = x[i];
    py = a * (py + xi - px);
    px = xi;
    x[i] = py;
  }
}

/**
 * Estimate the delay by cross-correlating three-second sliding windows.
 *
 * Both hands are resampled onto a uniform 1 kHz grid and high-passed at
 * 0.3 Hz to strip gyro bias drift, then correlated once over the whole
 * recording for the overview curve and again inside three-second windows
 * hopped every half second. Motion is bursty, so each window is gated
 * twice: windows whose geometric-mean energy falls below `efloor` are rest
 * and are recorded as rejected without being correlated at all, and windows
 * whose peak correlation misses `corr_gate` are drawn but not counted. The
 * reported delay is the median of the accepted windows, with their
 * interquartile range as the spread; a run with no accepted window falls
 * back to the whole-recording peak but will not be marked ok. Window and
 * hop are three and half a second because a window has to contain a whole
 * reach and still leave several windows in a short take, and the search is
 * capped at 1.5 s overall and 1 s per window, far wider than any real
 * teleoperation delay, so that a spurious match is visible rather than
 * folded into the range.
 *
 * @param[in] t   timestamps in seconds, ascending, both hands share them
 * @param[in] g1  robot gyro magnitude, so a positive result means it lags
 * @param[in] g2  operator gyro magnitude
 * @returns the estimate, its diagnostics and the plotting series, with `ok`
 *          false and `msg` set when the take was too short, too still, or
 *          the two hands were not clearly related
 * @exceptsafe basic
 */
static Result estimate(
    const std::vector<double>& t,
    const std::vector<double>& g1,
    const std::vector<double>& g2
) {
  Result r;
  if (t.size() < 100) {
    snprintf(r.msg, sizeof r.msg, "not enough samples");
    return r;
  }
  double dur = t.back() - t.front();
  double fs = 1000.0, dt = 1.0 / fs;
  size_t n = (size_t)(dur / dt);
  if (n < 512) {
    snprintf(r.msg, sizeof r.msg, "span too short (%.2fs)", dur);
    return r;
  }
  double t0 = t.front();

  auto a_raw = resample(t, g1, t0, dt, n);
  auto b_raw = resample(t, g2, t0, dt, n);
  auto a = a_raw, b = b_raw;
  highpass(a, 0.3, dt);
  highpass(b, 0.3, dt);

  long maxlag = (long)(1.5 * fs);  ///< far past any real delay
  if (maxlag > (long)n - 1) maxlag = (long)n - 1;
  XC xc = xcorr(a, b);
  Est ef = peak(xc, maxlag);
  r.tau_full_ms = ef.lag_samp / fs * 1000.0;
  r.corr = ef.corr;

  long span = std::min<long>((long)(0.8 * fs), maxlag);
  double norm = std::sqrt(xc.nx * xc.ny) + 1e-30;
  for (long L = -span; L <= span; L++) {
    size_t idx = (L >= 0) ? (size_t)L : (size_t)((long)xc.nf + L);
    r.cc_lag_ms.push_back((double)L / fs * 1000.0);
    r.cc_val.push_back(xc.z[idx] / norm);
  }

  size_t win = (size_t)(3.0 * fs), hop = (size_t)(0.5 * fs);
  const double efloor = 15.0;    ///< below this dps^2 the window is rest
  const double corr_gate = 0.6;  ///< weaker matches are drawn, not used
  std::vector<double> ests;
  for (size_t s = 0; win >= 1 && s + win <= n; s += hop) {
    double tc = t0 + (double)(s + win / 2) * dt;
    std::vector<double> aw(a.begin() + s, a.begin() + s + win);
    std::vector<double> bw(b.begin() + s, b.begin() + s + win);
    double ea = 0, eb = 0;
    for (double v : aw) ea += v * v;
    for (double v : bw) eb += v * v;
    ea /= (double)win;
    eb /= (double)win;
    if (std::sqrt(ea * eb) < efloor) {
      r.win_t.push_back(tc);
      r.win_tau.push_back(NAN);
      r.win_pass.push_back(0);
      continue;
    }
    XC xw = xcorr(aw, bw);
    Est ew = peak(xw, (long)(1.0 * fs));
    double tw = ew.lag_samp / fs * 1000.0;
    int pass = (ew.corr >= corr_gate) ? 1 : 0;
    r.win_t.push_back(tc);
    r.win_tau.push_back(tw);
    r.win_pass.push_back(pass);
    if (pass) ests.push_back(tw);
  }
  r.nwin = (int)ests.size();
  if (!ests.empty()) {
    std::sort(ests.begin(), ests.end());
    auto pct = [&](double p) {
      double x = p * (double)(ests.size() - 1);
      size_t i = (size_t)x;
      double f = x - (double)i;
      return (i + 1 < ests.size()) ? ests[i] * (1 - f) + ests[i + 1] * f
                                   : ests[i];
    };
    double med = pct(0.5);
    r.iqr_ms = pct(0.75) - pct(0.25);
    r.tau_ms = med;
  } else {
    r.tau_ms = r.tau_full_ms;
  }

  size_t stride = std::max<size_t>(1, n / 4000);
  for (size_t i = 0; i < n; i += stride) {
    r.ov_t.push_back(t0 + (double)i * dt);
    r.ov_ahp.push_back(a[i]);
    r.ov_bhp.push_back(b[i]);
  }

  r.ok = (r.nwin >= 3 && std::fabs(r.corr) > 0.3);
  if (!r.ok) {
    if (r.nwin < 3)
      snprintf(
          r.msg,
          sizeof r.msg,
          "insufficient shared motion (only %d gated window(s)) - move both "
          "sensors together",
          r.nwin
      );
    else
      snprintf(
          r.msg,
          sizeof r.msg,
          "low correlation (rho=%.2f) - signals not clearly related",
          r.corr
      );
  }
  return r;
}

}
namespace {

/**
 * Time the robot against the operator by windowed cross-correlation.
 *
 * This is the method the whole tool was built around and is now simply one
 * of the ten that vote, which is worth remembering when it disagrees with
 * the median: it has no special standing. The estimator high-passes the
 * gyro magnitude of both hands at 0.3 Hz, cross-correlates them inside
 * three-second windows, refines each window's peak to a fraction of a
 * sample, and reports the median over the windows that passed the rest and
 * correlation gates. Gyro magnitude is used because it is rotation
 * invariant, so the two IMUs never need a shared orientation, and the
 * windowing is what makes it work on bursty data, since a single
 * whole-recording correlation would be dominated by the largest reach and
 * would silently average over a delay that drifted. Its weakness is the
 * robot's low-pass behaviour, which broadens the correlation lobe and makes
 * the peak location less sharp than a landmark-based method's.
 * `Result::note` reports the whole-recording correlation and how many
 * windows were accepted, and on the abstaining path it carries the
 * estimator's own message instead. Unlike most of the other algorithms this
 * one fills its tile through `canvas::draw_panel` rather than drawing into
 * panel axes directly.
 * The generalized cross-correlation framing follows Knapp and Carter, The Generalized
 * Correlation Method for Estimation of Time Delay, IEEE Trans. ASSP 1976, using plain
 * weighting rather than the PHAT variant that paper also introduces.
 *
 * @param[in] in  both gyro magnitude channels plus their timestamps
 * @returns the delay in ms, positive when the robot lags the operator, with
 *          `ok` false when fewer than three windows were accepted or the
 *          whole-recording correlation stayed under 0.3
 * @exceptsafe basic
 */
algo::Result windowed_gcc(const algo::Input& in) {
  const lat::Result L = lat::estimate(in.t, in.gmag[0], in.gmag[1]);
  algo::Result r;
  r.ok = L.ok;
  r.tau_ms = L.tau_ms;
  r.spread_ms = L.iqr_ms;
  if (L.ok) {
    char buf[128];
    snprintf(buf, sizeof buf, "rho %.2f, %d windows", L.corr, L.nwin);
    r.note = buf;
  } else {
    r.note = L.msg;
  }

  r.image = canvas::tile(
      "WINDOWED GCC",
      r.tau_ms,
      r.ok,
      "CORRELATION PEAK PER WINDOW"
  );

  {
    canvas::Panel p;
    p.explain = "1  BOTH HANDS AFTER HIGH PASS, THE ROBOT TRAILS";
    p.xlabel = "TIME (s)";
    p.ylabel = "RATE (dps)";
    p.series = {
        {"OPERATOR", canvas::VIVID, &L.ov_t, &L.ov_bhp, canvas::LINE, 5},
        {"ROBOT", canvas::BLACK, &L.ov_t, &L.ov_ahp, canvas::LINE, 5}
    };
    canvas::draw_panel(r.image, 0, p);
  }
  {
    canvas::Panel p;
    p.explain = "2  CROSS CORRELATION, ITS PEAK IS THE DELAY";
    p.xlabel = "LAG (ms)";
    p.ylabel = "CORRELATION";
    p.series = {
        {"R(lag)", canvas::SAGE, &L.cc_lag_ms, &L.cc_val, canvas::LINE, 5}
    };
    p.fixed_range = true;
    p.x0 = -120;
    p.x1 = 200;
    p.y0 = -1.0;
    p.y1 = 1.0;
    p.guides = {{0.0, canvas::DIM, "ZERO LAG", true}};
    if (r.ok) p.guides.push_back({r.tau_ms, canvas::INK, "ESTIMATE", true});
    canvas::draw_panel(r.image, 1, p);
  }
  {
    std::vector<double> pass_t, pass_v, fail_t, fail_v;
    for (size_t i = 0; i < L.win_t.size(); i++) {
      if (std::isnan(L.win_tau[i])) continue;
      if (L.win_pass[i]) {
        pass_t.push_back(L.win_t[i]);
        pass_v.push_back(L.win_tau[i]);
      } else {
        fail_t.push_back(L.win_t[i]);
        fail_v.push_back(L.win_tau[i]);
      }
    }
    canvas::Panel p;
    p.explain = "3  ONE DELAY PER WINDOW, THE MEDIAN IS REPORTED";
    p.xlabel = "WINDOW CENTRE (s)";
    p.ylabel = "DELAY (ms)";
    p.fixed_range = true;
    p.x0 = L.win_t.empty() ? 0 : L.win_t.front();
    p.x1 = L.win_t.empty() ? 1 : L.win_t.back();
    p.y0 = 0;
    p.y1 = 150;
    p.series = {
        {"ACCEPTED", canvas::VIVID, &pass_t, &pass_v, canvas::MARKERS, 6},
        {"REJECTED", canvas::BLACK, &fail_t, &fail_v, canvas::MARKERS, 6}
    };
    if (r.ok) p.guides = {{r.tau_ms, canvas::INK, "MEDIAN", false}};
    canvas::draw_panel(r.image, 2, p);
  }

  return r;
}

}

ALGO_REGISTER(
    "windowed_gcc",
    windowed_gcc
);
