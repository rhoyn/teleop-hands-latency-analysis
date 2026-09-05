#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <vector>

#include "algo.hpp"

namespace {

/**
 * Tuning for burst segmentation and the banded warp.
 *
 * Every one of these was settled empirically on real recordings. The energy
 * gates are fractions of each hand's own 95th percentile rather than
 * absolute rotation rates, which is what lets the same numbers work across
 * the robot and the operator despite the robot moving far less; taking the
 * 95th percentile rather than the maximum stops one knock against the table
 * from raising the gate above every genuine burst. The burst limits come
 * from how the recordings are actually made: reaches lasting a few hundred
 * milliseconds, sometimes run together into trains with no rest between.
 */
constexpr double BAND_MS = 250.0;      ///< warp may stray past any real delay
constexpr double PAD_MS = 250.0;       ///< robot context each side of a burst
constexpr double SMOOTH_MS = 8.0;      ///< strips sensor noise, keeps motion
constexpr double ENV_MS = 60.0;        ///< envelope tracks a reach's rise
constexpr double ACTIVE_FRAC = 0.40;   ///< 40% of the 95th pct counts as moving
constexpr double BURST_FRAC = 0.15;    ///< lower gate, so burst edges survive
constexpr double WARP_PENALTY = 0.9;   ///< without it the path collapses
constexpr double MERGE_GAP_MS = 60.0;  ///< rests below this stay in one burst
constexpr double MIN_BURST_MS = 200.0;   ///< shorter cannot pin a 200 ms lag
constexpr double MAX_BURST_MS = 1000.0;  ///< longer runs are split into pieces

/**
 * One contiguous stretch of operator motion, as inclusive sample indices.
 *
 * Indices are into the full 1 kHz recording, not into any sub-array, so a
 * burst can be used to slice either hand without a further offset.
 */
struct Burst {
  size_t a;  ///< first sample of the burst, inclusive
  size_t b;  ///< last sample of the burst, inclusive
};

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
 * This is the motion envelope that the burst detector thresholds. Squaring
 * before smoothing is what makes it blind to the sign of the rotation rate,
 * so a fast oscillation reads as motion instead of averaging away to zero
 * the way a plain moving average of the high-passed signal would.
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
 * swallow every real burst. The vector is taken by value because it is
 * sorted in place, and the rank is truncated rather than interpolated,
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
 * the warp below to use a plain absolute difference as its local cost.
 * Fewer than eight active samples falls back to the whole signal rather
 * than dividing by a noise estimate.
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
 * Return the weighted median of a set of values.
 *
 * Each burst contributes one delay, but bursts differ wildly in how many
 * aligned samples backed that delay, so a long confident burst ought to
 * outweigh a brief one. A median rather than a weighted mean is used
 * because a single badly warped burst should not move the answer at all,
 * only fail to hold it back.
 *
 * @param[in] v  values, one per burst
 * @param[in] w  non-negative weight for each entry of `v`, same length
 * @returns the first value whose cumulative weight reaches half the total,
 *          or 0 when `v` is empty
 * @exceptsafe basic
 */
double weighted_median(
    const std::vector<double>& v,
    const std::vector<double>& w
) {
  if (v.empty()) return 0;
  std::vector<size_t> order(v.size());
  for (size_t i = 0; i < v.size(); i++) order[i] = i;
  std::sort(order.begin(), order.end(), [&](size_t a, size_t b) {
    return v[a] < v[b];
  });
  double total = 0;
  for (double x : w) total += x;
  double run = 0;
  for (size_t k : order) {
    run += w[k];
    if (run >= 0.5 * total) return v[k];
  }
  return v[order.back()];
}

/**
 * Cut bursts longer than `max_len` into equal contiguous pieces.
 *
 * The recordings are often trains of movements with no rest between them,
 * so energy gating alone can return a single burst tens of seconds long;
 * aligning that gives one delay, no population to take a median over and a
 * spread of exactly zero. Splitting restores a set of independent
 * estimates, at the cost of cut points that fall mid-motion and therefore
 * of pieces whose first and last samples have no matching context.
 *
 * @param[in] in       bursts as found by the energy gate
 * @param[in] max_len  longest burst kept intact, in samples
 * @returns the bursts with over-long ones replaced by their pieces, in
 *          time order
 * @exceptsafe basic
 */
std::vector<Burst> split_long(
    const std::vector<Burst>& in,
    size_t max_len
) {
  std::vector<Burst> out;
  for (const Burst& b : in) {
    const size_t len = b.b - b.a + 1;
    if (len <= max_len) {
      out.push_back(b);
      continue;
    }
    const size_t parts = (len + max_len - 1) / max_len;
    const size_t step = len / parts;
    for (size_t k = 0; k < parts; k++) {
      const size_t a = b.a + k * step;
      const size_t e = (k + 1 == parts) ? b.b : a + step - 1;
      out.push_back({a, e});
    }
  }
  return out;
}

/**
 * Group active samples into bursts, bridging short inactive gaps.
 *
 * A single reach makes the envelope dip below the gate several times on its
 * way through, so runs separated by fewer than `gap` samples are merged
 * instead of being reported as separate movements. Runs shorter than
 * `min_len` are dropped outright: they cannot support a delay measurement
 * that is itself allowed to reach 200 ms.
 *
 * @param[in] active   1 where the operator envelope is above the gate
 * @param[in] gap      inactive samples tolerated inside a single burst
 * @param[in] min_len  shortest burst kept, in samples
 * @returns inclusive index ranges, in time order and non-overlapping
 * @exceptsafe basic
 */
std::vector<Burst> find_bursts(
    const std::vector<char>& active,
    size_t gap,
    size_t min_len
) {
  std::vector<Burst> out;
  const size_t n = active.size();
  size_t i = 0;
  while (i < n) {
    if (!active[i]) {
      i++;
      continue;
    }
    size_t last = i, j = i;
    while (j < n && j - last <= gap) {
      if (active[j]) last = j;
      j++;
    }
    if (last - i + 1 >= min_len) out.push_back({i, last});
    i = last + 1;
  }
  return out;
}

/**
 * Align a query burst against a longer reference with banded DTW.
 *
 * The query must be consumed whole while the reference may be entered and
 * left anywhere, which is what lets an operator burst be located inside a
 * padded stretch of robot signal without knowing the delay in advance. The
 * band is centred on `diag`, the reference index the query would start at
 * with zero delay, so the cost is linear in the burst length rather than
 * quadratic and the path can never explain a quiet stretch with an absurd
 * shift. `pen` is charged on every non-diagonal step: without it the
 * cheapest path is a long many-to-one run that smears a single sample
 * across the burst and reports a meaningless offset.
 *
 * @param[in]  query  z-scored operator burst
 * @param[in]  ref    z-scored robot stretch, padded on both sides
 * @param[in]  diag   reference index `query[0]` sits at with zero delay
 * @param[in]  band   half-width of the search band, in samples
 * @param[in]  pen    cost added to each insertion or deletion step
 * @param[out] qi     query index of each path step, ascending
 * @param[out] rj     reference index matched to each entry of `qi`
 * @returns false when either series is under four samples, when the band
 *          leaves no admissible cell on some row, when every predecessor
 *          was unreachable, or when the backtrace runs off the reference;
 *          `qi` and `rj` are then not usable
 * @exceptsafe basic
 */
bool subsequence_dtw(
    const std::vector<double>& query,
    const std::vector<double>& ref,
    int diag,
    int band,
    double pen,
    std::vector<int>& qi,
    std::vector<int>& rj
) {
  const int m = (int)query.size();
  const int l = (int)ref.size();
  if (m < 4 || l < 4) return false;
  const int w = 2 * band + 1;
  const double inf = 1e300;
  std::vector<double> prev(l, inf), cur(l, inf);
  std::vector<uint8_t> step((size_t)m * (size_t)w, 0);
  for (int i = 0; i < m; i++) {
    const int lo = std::max(0, i + diag - band);
    const int hi = std::min(l - 1, i + diag + band);
    if (lo > hi) return false;
    for (int j = lo; j <= hi; j++) {
      const double c = std::fabs(query[(size_t)i] - ref[(size_t)j]);
      double best = 0;
      uint8_t s = 0;
      if (i > 0) {
        const double dg = (j > 0) ? prev[(size_t)j - 1] : inf;
        const double up = prev[(size_t)j];
        const double lf = (j > lo) ? cur[(size_t)j - 1] : inf;
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
        if (best >= inf) return false;
      }
      cur[(size_t)j] = best + c;
      step[(size_t)i * (size_t)w + (size_t)(j - i - diag + band)] = s;
    }
    prev = cur;
    std::fill(cur.begin(), cur.end(), inf);
  }
  const int lo = std::max(0, m - 1 + diag - band);
  const int hi = std::min(l - 1, m - 1 + diag + band);
  int bestj = lo;
  double bestv = inf;
  for (int j = lo; j <= hi; j++)
    if (prev[(size_t)j] < bestv) {
      bestv = prev[(size_t)j];
      bestj = j;
    }
  qi.clear();
  rj.clear();
  int i = m - 1, j = bestj;
  while (i >= 0) {
    qi.push_back(i);
    rj.push_back(j);
    if (i == 0) break;
    const uint8_t s =
        step[(size_t)i * (size_t)w + (size_t)(j - i - diag + band)];
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
    if (j < 0) return false;
  }
  std::reverse(qi.begin(), qi.end());
  std::reverse(rj.begin(), rj.end());
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
  r.image = canvas::standard("BURST DTW", r.tau_ms, r.ok);
  canvas::caption(r.image, "LAG PER BURST");
}

/**
 * Draw both motion envelopes with the detected burst edges marked.
 *
 * The traces are decimated to about 3000 points, all an 800 px panel can
 * resolve. Burst starts and ends are dropped onto the operator envelope,
 * since that is the hand the gate was applied to, so a reader can judge at
 * a glance whether the segmentation split the movements sensibly.
 *
 * @param[in,out] r       result whose image is drawn into
 * @param[in]     t       timestamps in seconds
 * @param[in]     er      robot envelope
 * @param[in]     eo      operator envelope
 * @param[in]     bursts  ranges found by the gate, may be empty
 * @exceptsafe basic
 */
void draw_energy(
    algo::Result& r,
    const std::vector<double>& t,
    const std::vector<double>& er,
    const std::vector<double>& eo,
    const std::vector<Burst>& bursts
) {
  canvas::panel_label(
      r.image,
      0,
      "1  MOTION ENERGY WITH THE BURSTS THAT WERE FOUND"
  );
  const size_t n = std::min(t.size(), std::min(er.size(), eo.size()));
  if (n < 4) return;
  const size_t stride = std::max<size_t>(1, n / 3000);
  std::vector<double> ts, vr, vo;
  for (size_t i = 0; i < n; i += stride) {
    ts.push_back(t[i] - t.front());
    vr.push_back(er[i]);
    vo.push_back(eo[i]);
  }
  std::vector<double> both = vr;
  both.insert(both.end(), vo.begin(), vo.end());
  canvas::Axes a0 = canvas::panel(0);
  canvas::autoscale(a0, ts, both);
  canvas::set_range(a0, ts.front(), ts.back(), 0.0, a0.y1);
  canvas::grid_and_ticks(r.image, a0, 10, 5);
  canvas::polyline(r.image, a0, ts, vr, canvas::BLACK);
  canvas::polyline(r.image, a0, ts, vo, canvas::VIVID);
  for (const Burst& b : bursts) {
    if (b.a >= n || b.b >= n) continue;
    canvas::marker(r.image, a0, t[b.a] - t.front(), eo[b.a], canvas::SHADE, 4);
    canvas::marker(r.image, a0, t[b.b] - t.front(), eo[b.b], canvas::BLACK, 4);
  }
  canvas::axis_labels(r.image, a0, "TIME (s)", "ENERGY (dps)");
  canvas::legend(
      r.image,
      a0,
      {{"OPERATOR", canvas::VIVID},
       {"ROBOT", canvas::BLACK},
       {"BURST START", canvas::WHITE, true},
       {"BURST END", canvas::BLACK, true}},
      true
  );
}

/**
 * Estimate the delay by aligning each motion burst separately with DTW.
 *
 * Both gyro magnitudes are high-passed at 0.3 Hz, box-smoothed over
 * `SMOOTH_MS` and z-scored over their active samples only. An RMS envelope
 * gates the operator into bursts, which are merged across short dips and
 * then split at `MAX_BURST_MS`; each burst is aligned against a padded
 * stretch of robot signal by banded subsequence DTW. Every query sample's
 * mean matched reference index yields one offset in ms, the burst's delay
 * is the median of the offsets whose endpoints both fall in motion, and the
 * recording's delay is the weighted median across bursts, weighted by how
 * many samples backed each one.
 *
 * Working burst by burst suits this data specifically because the robot
 * attenuates and low-passes the operator rather than copying it: DTW is
 * free to warp within a burst and so tolerates a distorted response where a
 * single fixed shift would not, and segmenting first stops a long rest from
 * contributing anything. The costs are real too. The delay is only ever
 * measured where there is motion, a burst the robot barely responded to is
 * silently dropped rather than counted as disagreement, and the split
 * points inside a train of movements are arbitrary, so the per-burst spread
 * understates how uncertain a single burst really is.
 *
 * `Result::note` reports how many of the detected bursts produced a usable
 * alignment and the band half-width; `Result::spread_ms` is the IQR of the
 * per-burst delays. The three panels show the two envelopes with the burst
 * edges marked, the best-supported single burst with the robot shifted by
 * the estimate so the reader can judge the fit by eye, and every burst's
 * delay against its centre time with the median drawn across.
 * The band-constrained dynamic programming follows Sakoe and Chiba, Dynamic
 * programming algorithm optimization for spoken word recognition, IEEE Trans. ASSP
 * 1978, and the open-boundary subsequence variant follows Mueller, Information
 * Retrieval for Music and Motion, chapter 4.
 *
 * @param[in] in  both 1 kHz gyro magnitudes with their shared timebase
 * @returns the delay in ms, positive when the robot lags the operator, or
 *          an abstaining `Result` when the recording is under 64 samples,
 *          no burst passed the energy gate, or none of them aligned
 * @exceptsafe basic
 */
algo::Result burst_dtw(const algo::Input& in) {
  algo::Result r;
  const size_t n = in.t.size();
  if (n < 64) {
    r.note = "not enough samples";
    reject_tile(r);
    return r;
  }
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

  const size_t gap = std::max<size_t>(1, size_t(MERGE_GAP_MS / 1000.0 / in.dt));
  const size_t minlen =
      std::max<size_t>(8, size_t(MIN_BURST_MS / 1000.0 / in.dt));
  std::vector<char> ao_burst(n);
  const double gb = BURST_FRAC * percentile(eo, 0.95);
  for (size_t i = 0; i < n; i++) ao_burst[i] = eo[i] > gb ? 1 : 0;
  const size_t maxlen =
      std::max<size_t>(minlen * 2, size_t(MAX_BURST_MS / 1000.0 / in.dt));
  const std::vector<Burst> bursts =
      split_long(find_bursts(ao_burst, gap, minlen), maxlen);
  if (bursts.empty()) {
    r.note = "no motion burst passed the energy threshold";
    reject_tile(r);
    draw_energy(r, in.t, er, eo, bursts);
    return r;
  }

  const int band = std::max(8, int(BAND_MS / 1000.0 / in.dt));
  const int pad = std::max(band, int(PAD_MS / 1000.0 / in.dt));
  std::vector<double> taus, weights, centres;
  std::vector<double> rep_ms, rep_op, rep_rob;
  double rep_weight = -1;
  for (const Burst& b : bursts) {
    const size_t r0 = b.a > (size_t)pad ? b.a - (size_t)pad : 0;
    const size_t r1 = std::min(n - 1, b.b + (size_t)pad);
    const std::vector<double> query(
        operat.begin() + (long)b.a,
        operat.begin() + (long)b.b + 1
    );
    const std::vector<double> ref(
        robot.begin() + (long)r0,
        robot.begin() + (long)r1 + 1
    );
    std::vector<int> qi, rj;
    if (!subsequence_dtw(
            query,
            ref,
            (int)(b.a - r0),
            band,
            WARP_PENALTY,
            qi,
            rj
        ))
      continue;
    const size_t m = query.size();
    std::vector<double> acc(m, 0.0), cnt(m, 0.0);
    for (size_t k = 0; k < qi.size(); k++) {
      acc[(size_t)qi[k]] += double(rj[k]);
      cnt[(size_t)qi[k]] += 1.0;
    }
    std::vector<double> local;
    for (size_t i = 0; i < m; i++) {
      if (cnt[i] < 1.0) continue;
      const long gi = (long)r0 + std::lround(acc[i] / cnt[i]);
      if (gi < 0 || gi >= (long)n || !ar[(size_t)gi] || !ao[b.a + i]) continue;
      local.push_back(double(gi - (long)(b.a + i)) * in.dt * 1000.0);
    }
    if (local.size() < minlen / 4) continue;
    const double burst_tau = algo::median(local);
    taus.push_back(burst_tau);
    weights.push_back(double(local.size()));
    centres.push_back(0.5 * (in.t[b.a] + in.t[b.b]) - in.t.front());
    if (double(local.size()) > rep_weight) {
      rep_weight = double(local.size());
      const long shift = std::lround(burst_tau / 1000.0 / in.dt);
      const size_t rstride = std::max<size_t>(1, m / 800);
      rep_ms.clear();
      rep_op.clear();
      rep_rob.clear();
      for (size_t k = 0; k < m; k += rstride) {
        const long gi = (long)(b.a + k) + shift;
        if (gi < 0 || gi >= (long)n) continue;
        rep_ms.push_back((in.t[b.a + k] - in.t[b.a]) * 1000.0);
        rep_op.push_back(operat[b.a + k]);
        rep_rob.push_back(robot[(size_t)gi]);
      }
    }
  }

  if (taus.empty()) {
    r.note = "no burst produced a usable alignment";
    reject_tile(r);
    draw_energy(r, in.t, er, eo, bursts);
    return r;
  }
  r.ok = true;
  r.tau_ms = weighted_median(taus, weights);
  r.spread_ms = algo::iqr(taus);
  char buf[128];
  snprintf(
      buf,
      sizeof buf,
      "%d/%d bursts aligned, band +/-%.0f ms",
      (int)taus.size(),
      (int)bursts.size(),
      BAND_MS
  );
  r.note = buf;

  r.image = canvas::standard("BURST DTW", r.tau_ms, r.ok);
  canvas::caption(r.image, "LAG PER BURST");
  draw_energy(r, in.t, er, eo, bursts);

  canvas::panel_label(
      r.image,
      1,
      "2  ONE BURST WITH THE ROBOT SHIFTED BY THE DELAY"
  );
  if (rep_ms.size() > 1) {
    std::vector<double> both = rep_op;
    both.insert(both.end(), rep_rob.begin(), rep_rob.end());
    canvas::Axes a1 = canvas::panel(1);
    canvas::autoscale(a1, rep_ms, both);
    canvas::set_range(a1, rep_ms.front(), rep_ms.back(), a1.y0, a1.y1);
    canvas::grid_and_ticks(r.image, a1, 10, 5);
    canvas::polyline(r.image, a1, rep_ms, rep_rob, canvas::BLACK);
    canvas::polyline(r.image, a1, rep_ms, rep_op, canvas::VIVID);
    canvas::axis_labels(r.image, a1, "TIME IN BURST (ms)", "Z SCORE");
    canvas::legend(
        r.image,
        a1,
        {{"OPERATOR", canvas::VIVID}, {"SHIFTED ROBOT", canvas::BLACK}},
        true
    );
  }

  canvas::panel_label(
      r.image,
      2,
      "3  DELAY MEASURED IN EACH BURST AND THEIR MEDIAN"
  );
  canvas::Axes a2 = canvas::panel(2);
  canvas::autoscale(a2, centres, taus);
  canvas::set_range(a2, 0.0, in.t.back() - in.t.front(), a2.y0, a2.y1);
  canvas::grid_and_ticks(r.image, a2, 10, 5);
  for (size_t i = 0; i < taus.size(); i++)
    canvas::marker(r.image, a2, centres[i], taus[i], canvas::SLATE, 5);
  canvas::mark_guide(r.image, a2, r.tau_ms, false, "MEDIAN", 0);
  canvas::axis_labels(r.image, a2, "BURST CENTRE (s)", "DELAY (ms)");
  canvas::legend(r.image, a2, {{"BURST DELAY", canvas::WHITE, true}}, true);
  return r;
}

}

ALGO_REGISTER(
    "burst_dtw",
    burst_dtw
);
