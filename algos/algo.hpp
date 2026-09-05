#pragma once

#include <algorithm>
#include <cmath>
#include <complex>
#include <string>
#include <vector>

#include "canvas.hpp"

namespace algo {

using cd = std::complex<double>;

/**
 * Number of MPU6050 sensors a recording carries.
 *
 * The rig is fixed at two, one strapped to each hand, and every per-sensor
 * array is sized by this. It is a description of the hardware rather than a
 * tunable: raising it would not by itself teach `load_csv` to parse a third
 * sensor's six columns.
 */
inline constexpr int MPUS = 2;

/**
 * Index of the robot-hand sensor, MPU #1 in the CSV column order.
 *
 * This is the following hand, so its motion is the delayed copy. It is the
 * one wired to `I2C_A` and read by `loop()`, so its six columns arrive
 * first in every row. The robot is taken to follow the operator always,
 * which is the only direction a working teleoperation link can produce, so
 * `run_all` reports the magnitude of every estimate and the tool never
 * emits a negative lag.
 */
inline constexpr int ROBOT = 0;

/**
 * Index of the operator-hand sensor, MPU #2 in the CSV column order.
 *
 * This is the reference, the hand that moves first. Estimators search for
 * the shift that maps this signal onto `ROBOT`. The two indices are not
 * interchangeable: several estimators treat the reference and the delayed
 * copy differently, so swapping them changes the magnitude of the answer,
 * not just its sign.
 */
inline constexpr int OPERATOR = 1;

/**
 * Widest lag, in milliseconds, that any estimator is allowed to consider.
 *
 * Search windows are clamped to plus or minus this value, which bounds cost
 * and, more importantly, stops a repetitive hand motion from being matched
 * a whole cycle late. A link slower than 200 ms is outside what this tool
 * is built to measure, and an estimate pinned against the clamp is a sign
 * the two recordings do not describe the same motion.
 */
inline constexpr double MAX_LAG_MS = 200.0;

/**
 * The signals every algorithm sees, identical across all ten of them.
 *
 * Only the rotation-invariant gyro magnitudes are handed over, never the
 * raw axes: the two sensors sit on hands with no shared orientation, so
 * anything axis-wise would need an alignment that was never measured. All
 * vectors have the same length, one entry per accepted CSV row.
 */
struct Input {
  std::vector<double> t;           ///< seconds since the first sample
  std::vector<double> gmag[MPUS];  ///< gyro magnitude in dps, per sensor
  double dt = 0.001;               ///< measured mean sample period, seconds
};

/**
 * What one algorithm concluded about one recording.
 *
 * `ok` is a status, not an error flag. False means the algorithm inspected
 * the recording and declined to answer; the run continues, the tile renders
 * as "rejected" and the tau is left out of the median. Genuine failures are
 * thrown as exceptions and never encoded here, so any `Result` that comes
 * back is one of those two honest outcomes.
 */
struct Result {
  std::string name;      ///< stamped by `run_all` from the registry entry
  bool ok = false;       ///< accepted, not error-free: false is abstention
  double tau_ms = 0;     ///< lag magnitude in ms, robot behind operator
  double spread_ms = 0;  ///< this algorithm's own uncertainty, in ms
  std::string note;      ///< one line of reasoning drawn onto the tile
  canvas::Image image;   ///< the 800x800 tile this algorithm rendered
};

/**
 * Signature every algorithm implements: one recording in, one verdict out.
 *
 * A plain function pointer rather than `std::function` because every
 * algorithm is a free function and the registry is populated before `main`,
 * where a type with a non-trivial constructor would be one more thing to
 * order correctly.
 */
using Fn = Result (*)(const Input&);

/**
 * One registry slot: an algorithm's display name and its entry point.
 *
 * The name is held here rather than inside the algorithm so `run_all` can
 * stamp it onto every `Result` after the call. That keeps the results table
 * and the tiles labelled consistently even if an algorithm neglects to name
 * itself.
 */
struct Entry {
  std::string name;
  Fn fn;
};

/**
 * Return the process-wide list of registered algorithms.
 *
 * A function rather than a global object so the vector is constructed on
 * first use. That is what lets registrations run safely during static
 * initialisation, in whatever order the linker chose.
 *
 * @returns the mutable registry, in the order entries were added
 * @exceptsafe basic
 */
std::vector<Entry>& registry();

/**
 * Self-registering handle that adds one algorithm to the registry.
 *
 * A namespace-scope instance runs its constructor during static
 * initialisation, before `main`, so merely linking a translation unit is
 * what enrols its algorithm. Declare one through `ALGO_REGISTER` rather
 * than by hand.
 */
struct Register {
  /**
   * Append one algorithm to the registry.
   *
   * Reaches the list through `registry()` rather than through a global
   * object, because at this point in start-up a global vector in another
   * translation unit may not have been constructed yet.
   *
   * @param[in] name  label shown in the results table and on the tile
   * @param[in] fn    entry point invoked once per recording
   * @exceptsafe basic
   */
  Register(
      const char* name,
      Fn fn
  ) {
    registry().push_back({name, fn});
  }
};

/**
 * Enrol an algorithm in the registry from its own translation unit.
 *
 * Expands to a namespace-scope `Register` object named after `FN`, whose
 * constructor runs before `main` and appends the name and pointer to
 * `registry()`. Nothing central has to be edited to add an algorithm:
 * dropping a .cpp into `algos/` suffices, since the Makefile globs the
 * directory and the constructor does the enrolling. Invoke it at namespace
 * scope, once per algorithm, and terminate the invocation with a semicolon.
 *
 * @param[in] NAME  display name, a string literal
 * @param[in] FN    the `Fn` to register; must be a bare identifier, since
 *                  it is pasted into the generated object's name
 */
#define ALGO_REGISTER(NAME, FN) \
  static ::algo::Register FN##_registration(NAME, FN)

/**
 * Run every registered algorithm over one recording.
 *
 * The ten estimators are deliberately independent and none may veto
 * another, so nothing is filtered here: abstentions come back alongside
 * accepted results and are sorted out downstream.
 *
 * @param[in] in  the shared signals and sample period
 * @returns one `Result` per registered algorithm, in registration order
 * @throws std::runtime_error if an algorithm throws, since a throw is a
 *         defect rather than an abstention
 * @exceptsafe basic
 */
std::vector<Result> run_all(const Input& in);

/**
 * Median lag across the algorithms that accepted the recording.
 *
 * This is the figure the tool reports. Abstentions are dropped before the
 * median is taken rather than counted as zero, because `tau_ms` carries no
 * meaning on a rejected result.
 *
 * @param[in] r  every result from `run_all`
 * @returns the median lag in milliseconds, or 0
 *          if no algorithm accepted
 * @exceptsafe basic
 */
double median_tau_ms(const std::vector<Result>& r);

/**
 * Median absolute deviation of the accepted lags.
 *
 * Quoted as the spread beside the headline figure. It is the raw MAD, with
 * no 1.4826 scaling to a normal-equivalent sigma, so it must not be read as
 * a standard deviation.
 *
 * @param[in] r  every result from `run_all`
 * @returns the MAD in milliseconds, or 0 if no algorithm accepted
 * @exceptsafe basic
 */
double mad_tau_ms(const std::vector<Result>& r);

/**
 * Count the algorithms that produced a usable estimate.
 *
 * Reported as "n of 10" so a reader can weigh the median: the same number
 * backed by three estimators deserves less trust than one backed by ten,
 * and nothing else in the output distinguishes the two.
 *
 * @param[in] r  every result from `run_all`
 * @returns how many results carry `ok` true
 * @exceptsafe no-throw
 */
int accepted(const std::vector<Result>& r);

/**
 * High-pass a signal in place with a one-pole filter.
 *
 * The pole sits at exp(-2*pi*fc*dt), the usual DC-blocker mapping, accurate
 * while `fc` stays far below the 500 Hz Nyquist of a 1 kHz recording and
 * increasingly warped as it approaches it. Filtering this way is what
 * strips per-sensor gyro bias and slow drift while leaving the hand motion,
 * and because it is causal and applied identically to both channels its own
 * group delay cancels out of any lag estimate. State starts at zero, so the
 * first few time constants of output are a settling transient.
 *
 * @param[in,out] x   signal to filter, replaced by its high-passed self
 * @param[in]     fc  cutoff frequency in Hz
 * @param[in]     dt  sample period in seconds
 * @exceptsafe no-throw
 */
inline void highpass(
    std::vector<double>& x,
    double fc,
    double dt
) {
  if (x.empty()) return;
  const double a = std::exp(-2 * M_PI * fc * dt);
  double px = x[0], py = 0;
  for (size_t i = 0; i < x.size(); i++) {
    const double xi = x[i];
    py = a * (py + xi - px);
    px = xi;
    x[i] = py;
  }
}

/**
 * Round a length up to the next power of two, at least one.
 *
 * The radix-2 `fft` accepts nothing else, and padding a correlation out to
 * such a length is also what keeps circular wrap-around from folding one
 * end of the recording onto the other. Returns 1 for 0 and for 1; an `n`
 * above the largest representable power of two would loop forever rather
 * than saturate, a size no caller in this tool can reach.
 *
 * @param[in] n  minimum length required
 * @returns the smallest power of two that is at least `n`
 * @exceptsafe no-throw
 */
inline size_t next_pow2(size_t n) {
  size_t p = 1;
  while (p < n) p <<= 1;
  return p;
}

/**
 * Transform a sequence in place, forward or inverse.
 *
 * Iterative radix-2 Cooley-Tukey: a bit-reversal permutation followed by
 * log2(n) butterfly passes. The length must be a power of two and this is
 * not checked, so a wrongly sized input yields silent garbage instead of an
 * error; size it with `next_pow2` first. Twiddle factors are advanced by
 * repeated multiplication inside each block, which is cheap but lets a
 * little phase error accumulate across long transforms.
 *
 * @param[in,out] a    sequence to transform, overwritten by its transform
 * @param[in]     inv  true for the inverse, which also divides by n so that
 *                     a forward and inverse pair is the identity
 * @exceptsafe no-throw
 */
inline void fft(
    std::vector<cd>& a,
    bool inv
) {
  const size_t n = a.size();
  for (size_t i = 1, j = 0; i < n; i++) {
    size_t bit = n >> 1;
    for (; j & bit; bit >>= 1) j ^= bit;
    j ^= bit;
    if (i < j) std::swap(a[i], a[j]);
  }
  for (size_t len = 2; len <= n; len <<= 1) {
    const double ang = 2 * M_PI / double(len) * (inv ? 1 : -1);
    const cd wlen(std::cos(ang), std::sin(ang));
    for (size_t i = 0; i < n; i += len) {
      cd w(1);
      for (size_t k = 0; k < len / 2; k++) {
        const cd u = a[i + k], v = a[i + k + len / 2] * w;
        a[i + k] = u + v;
        a[i + k + len / 2] = u - v;
        w *= wlen;
      }
    }
  }
  if (inv)
    for (cd& z : a) z /= double(n);
}

/**
 * Return the median of a sample.
 *
 * Takes its argument by value on purpose: callers pass vectors they still
 * need in the original order, and sorting in place would quietly corrupt
 * them. Even-sized samples average the two central values. An empty sample
 * gives 0, which is indistinguishable from a genuine zero median, so a
 * caller that cares must test for emptiness itself.
 *
 * @param[in] v  sample, copied and sorted internally
 * @returns the median, or 0 if `v` is empty
 * @exceptsafe basic
 */
inline double median(std::vector<double> v) {
  if (v.empty()) return 0;
  std::sort(v.begin(), v.end());
  const size_t n = v.size();
  return (n % 2) ? v[n / 2] : 0.5 * (v[n / 2 - 1] + v[n / 2]);
}

/**
 * Return the interquartile range of a sample.
 *
 * Quartiles are taken by nearest rank with no interpolation, so on small
 * samples the value steps rather than moving smoothly as points are added.
 * Fewer than four points cannot bracket a quartile at all and yield 0, the
 * same answer a genuinely constant sample gives.
 *
 * @param[in] v  sample, copied and sorted internally
 * @returns the gap between the upper and lower quartile, or 0 if fewer than
 *          four points were given
 * @exceptsafe basic
 */
inline double iqr(std::vector<double> v) {
  if (v.size() < 4) return 0;
  std::sort(v.begin(), v.end());
  return v[(3 * v.size()) / 4] - v[v.size() / 4];
}

/**
 * Return the arithmetic mean of a sample.
 *
 * Plain running summation with no compensation: across the roughly 10^5
 * samples of one recording the accumulated rounding stays orders of
 * magnitude below the sensor noise, so a Kahan sum would buy nothing here.
 *
 * @param[in] v  sample to average
 * @returns the mean, or 0 if `v` is empty
 * @exceptsafe no-throw
 */
inline double mean(const std::vector<double>& v) {
  if (v.empty()) return 0;
  double s = 0;
  for (double x : v) s += x;
  return s / double(v.size());
}

/**
 * Build the comparable motion signal for one sensor.
 *
 * The gyro magnitude is used instead of individual axes because it is
 * rotation invariant: the two MPU6050s ride on hands at unrelated
 * attitudes and were never aligned to a common frame, and a magnitude needs
 * no such alignment. High-passing at `hp_fc` then removes the per-sensor
 * bias and slow drift and leaves the deliberate hand motion the estimators
 * line up; both channels receive the same filter, so its delay cancels.
 *
 * @param[in] in     the loaded recording
 * @param[in] mpu    `ROBOT` or `OPERATOR`
 * @param[in] hp_fc  high-pass cutoff in Hz, chosen per algorithm
 * @returns the filtered gyro magnitude in dps, one value per sample
 * @exceptsafe basic
 */
inline std::vector<double> motion_signal(
    const Input& in,
    int mpu,
    double hp_fc
) {
  std::vector<double> x = in.gmag[mpu];
  highpass(x, hp_fc, in.dt);
  return x;
}

/**
 * Rescale a signal to zero mean and unit standard deviation.
 *
 * Lets the two hands be compared on shape alone. The operator generally
 * moves harder than the robot manages to repeat, and without this the
 * livelier channel would dominate any squared-error or dot-product
 * comparison. The deviation is the population form, divided by n; a spread
 * below 1e-12 is replaced by 1 so a constant input flattens to zeros rather
 * than turning into NaN and poisoning everything downstream.
 *
 * @param[in] x  signal, taken by value and returned normalised
 * @returns the normalised copy, or the input unchanged if it was empty
 * @exceptsafe basic
 */
inline std::vector<double> znorm(std::vector<double> x) {
  if (x.empty()) return x;
  const double m = mean(x);
  double v = 0;
  for (double a : x) v += (a - m) * (a - m);
  v = std::sqrt(v / double(x.size()));
  if (v < 1e-12) v = 1.0;
  for (double& a : x) a = (a - m) / v;
  return x;
}

}
