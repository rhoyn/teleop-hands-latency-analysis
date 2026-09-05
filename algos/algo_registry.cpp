#include <cstdio>
#include <filesystem>
#include <fstream>
#include <limits.h>
#include <unistd.h>

#include "algo.hpp"

#define STB_TRUETYPE_IMPLEMENTATION
#include <stb_truetype.h>

namespace canvas {

/**
 * Load the bundled monospace font once and hand back the shared instance.
 *
 * The file is looked for beside the executable, then one directory up, then
 * in the working directory, so the tool behaves the same run from `build/`
 * as from the repository root. The search happens exactly once: if it
 * fails, a warning goes to stderr and every later call returns the same
 * unusable `Font`, so tiles render without text rather than aborting a run
 * whose measurements have already been made. Text is decoration here; the
 * numbers also go to stdout.
 *
 * @returns the process-wide font, with `Font::ok` false if none was loaded
 * @exceptsafe basic
 */
Font& font() {
  static Font f;
  static bool tried = false;
  if (tried) return f;
  tried = true;
  char self[PATH_MAX];
  const ssize_t n = readlink("/proc/self/exe", self, sizeof(self) - 1);
  std::filesystem::path base =
      std::filesystem::path(n > 0 ? std::string(self, (size_t)n) : ".")
          .parent_path();
  for (const std::filesystem::path& root :
       {base, base.parent_path(), std::filesystem::path(".")}) {
    const std::filesystem::path p = root / "assets" / "JetBrainsMono.ttf";
    if (!std::filesystem::exists(p)) continue;
    std::ifstream in(p, std::ios::binary);
    f.data.assign(
        std::istreambuf_iterator<char>(in),
        std::istreambuf_iterator<char>()
    );
    if (f.data.empty()) continue;
    if (stbtt_InitFont(
            &f.info,
            f.data.data(),
            stbtt_GetFontOffsetForIndex(f.data.data(), 0)
        ))
      f.ok = true;
    break;
  }
  if (!f.ok)
    fprintf(
        stderr,
        "warning: assets/JetBrainsMono.ttf not found, tiles will have no text\n"
    );
  return f;
}

}

namespace algo {

/**
 * Return the process-wide list of registered algorithms.
 *
 * The vector is a function-local static, so it is built on first use rather
 * than at load time. That ordering is what makes `ALGO_REGISTER` sound:
 * registrations fire during static initialisation in an order the standard
 * leaves unspecified across translation units, and each one reaches the
 * list through this call instead of touching an object that might not be
 * alive yet.
 *
 * @returns the mutable registry, in the order entries were added
 * @exceptsafe basic
 */
std::vector<Entry>& registry() {
  static std::vector<Entry> entries;
  return entries;
}

/**
 * Run every registered algorithm over the same recording.
 *
 * The ten estimators are independent by construction and none may veto
 * another, so this loop filters nothing: an abstention arrives as a
 * `Result` with `ok` false and is carried through to the report, where it
 * renders as "rejected". The name is stamped from the registry entry after
 * the call, so an algorithm cannot mislabel itself in the results table.
 *
 * @param[in] in  the shared signals and sample period
 * @returns one `Result` per registered algorithm, in registration order
 * @throws std::runtime_error if an algorithm throws; failures propagate
 *         because a throw is a defect, whereas abstention is a status
 * @exceptsafe basic
 */
std::vector<Result> run_all(const Input& in) {
  std::vector<Result> out;
  for (const Entry& e : registry()) {
    Result r = e.fn(in);
    r.name = e.name;
    r.tau_ms = std::fabs(r.tau_ms);
    out.push_back(r);
  }
  return out;
}

/**
 * Take the median lag over the algorithms that accepted.
 *
 * This is the number the tool reports, and taking it across ten
 * independent methods is the whole design: no single technique decides the
 * answer, and one estimator that went astray cannot pull a median far.
 * Abstentions are dropped rather than counted as zero, since `tau_ms` means
 * nothing on a rejected result.
 *
 * @param[in] r  every result from `run_all`
 * @returns the median lag in milliseconds, or 0 if no algorithm accepted
 * @exceptsafe basic
 */
double median_tau_ms(const std::vector<Result>& r) {
  std::vector<double> v;
  for (const Result& x : r)
    if (x.ok) v.push_back(x.tau_ms);
  return median(v);
}

/**
 * Measure the spread of the accepted lags as a median absolute deviation.
 *
 * The MAD is used rather than a standard deviation because the sample is
 * ten values that may include a wild one, and a squared measure would let
 * that single value describe the whole set. It is returned unscaled, with
 * no 1.4826 conversion to a normal-equivalent sigma, so it is not
 * comparable to one.
 *
 * @param[in] r  every result from `run_all`
 * @returns the MAD in milliseconds, or 0 if no algorithm accepted
 * @exceptsafe basic
 */
double mad_tau_ms(const std::vector<Result>& r) {
  std::vector<double> v;
  for (const Result& x : r)
    if (x.ok) v.push_back(x.tau_ms);
  if (v.empty()) return 0;
  const double m = median(v);
  std::vector<double> d;
  for (double x : v) d.push_back(std::fabs(x - m));
  return median(d);
}

/**
 * Count the algorithms that produced a usable estimate.
 *
 * Printed as "n of 10" beside the median so a reader can weigh it: the same
 * figure backed by three estimators deserves less confidence than one
 * backed by ten, and the median alone does not show the difference.
 *
 * @param[in] r  every result from `run_all`
 * @returns how many results carry `ok` true
 * @exceptsafe no-throw
 */
int accepted(const std::vector<Result>& r) {
  int n = 0;
  for (const Result& x : r)
    if (x.ok) n++;
  return n;
}

}
