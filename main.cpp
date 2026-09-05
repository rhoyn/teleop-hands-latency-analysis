
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "algos/algo.hpp"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include <stb_image_write.h>

/**
 * Word drawn across the branding tile that opens the output grid.
 *
 * The tile exists so that a screenshot of the report, detached from any
 * filename or directory, still says where it came from. It is drawn one
 * glyph at a time rather than as a single string, because the logo wants
 * letter tracking that `canvas::text` does not provide.
 */
static const char* LOGO_TEXT = "rhoyn";

/**
 * Extra space inserted between logo glyphs, as a fraction of glyph height.
 *
 * Expressed relative to the glyph box rather than in pixels so the logo
 * keeps its proportions if `LOGO_SCALE` changes. `make_logo_tile` measures
 * the word with the same rule and then subtracts one tracking step, since
 * the layout loop also adds a gap after the final letter.
 */
static constexpr float LOGO_TRACKING = 0.18f;

/**
 * Glyph scale for the logo, in `canvas::PX_PER_SCALE` units.
 *
 * At nine pixels per unit this gives a 189 px tall word inside the 800 px
 * tile, about the largest size that still leaves margin once tracking has
 * widened the word.
 */
static constexpr int LOGO_SCALE = 21;

/**
 * Glyph scale for the headline latency figure on the results tile.
 *
 * Four times the ordinary title scale, so the one number a reader looks for
 * first dominates its tile and survives the whole grid being viewed shrunk
 * to fit a screen.
 */
static constexpr int RESULT_SCALE = canvas::TITLE_SCALE * 4;

/**
 * JPEG quality used when writing the stitched report image.
 *
 * The tiles are flat colour crossed by thin text and hairline plots, which
 * is precisely what JPEG handles worst; 92 keeps label edges clean while
 * still encoding a multi-megapixel grid into a file worth attaching to a
 * message.
 */
static constexpr int JPEG_QUALITY = 92;

/**
 * One recording as `load_csv` lays it out in memory.
 *
 * Both the raw columns and their magnitudes are kept: the raw form is what
 * the CSV holds and what makes a parsing problem visible, while the
 * estimators only ever consume the rotation-invariant magnitudes. Every
 * vector here has the same length, one entry per accepted row.
 */
struct Data {
  std::vector<double> t;        ///< seconds since the first sample
  std::vector<double> s[12];    ///< raw columns, 6 per MPU: ax,ay,az,gx,gy,gz
  std::vector<double> gmag[2];  ///< gyro magnitude in dps, per MPU
};

/**
 * Print the command line synopsis.
 *
 * Goes to stderr, where every other diagnostic in this tool goes, so that
 * redirecting stdout to capture the measurements still leaves the complaint
 * visible on the terminal.
 *
 * @param[in] prog  program name as invoked, normally `argv[0]`
 * @exceptsafe no-throw
 */
static void usage(const char* prog) {
  fprintf(stderr, "usage: %s INPUT.csv OUTPUT.jpg\n", prog);
}

/**
 * Read a recording from CSV into per-channel arrays.
 *
 * A row holds thirteen comma separated fields: a 32-bit microsecond
 * timestamp followed by six channels for each MPU6050, ax,ay,az in g then
 * gx,gy,gz in degrees per second. That counter wraps about every 71
 * minutes, so a raw value below its predecessor is taken as a wrap and adds
 * 2^32 microseconds to a running base; times are then rebased on the first
 * sample and converted to seconds, which keeps the doubles small enough
 * that microsecond resolution is never at risk. Lines not starting with a
 * digit or minus, and rows with fewer than thirteen parsable fields, are
 * skipped in silence, so a header line or a recording truncated mid-row
 * costs one sample rather than the run.
 *
 * @param[in] path  CSV file to read
 * @returns the recording, with times in seconds from the first sample and
 *          the accelerometer and gyro magnitudes precomputed
 * @throws std::runtime_error if the file cannot be opened, or if it holds
 *         no usable row
 * @exceptsafe basic
 */
static Data load_csv(const std::string& path) {
  Data d;
  std::ifstream f(path);
  if (!f) throw std::runtime_error("cannot open " + path);
  std::string line;
  bool first = true;
  uint64_t base = 0, prev = 0, t0 = 0;
  long rows = 0;
  while (std::getline(f, line)) {
    if (line.empty()) continue;
    char c = line[0];
    if (c != '-' && (c < '0' || c > '9')) continue;
    double v[13];
    const char* p = line.c_str();
    char* end;
    int i = 0;
    for (; i < 13; i++) {
      v[i] = strtod(p, &end);
      if (end == p) break;
      p = end;
      if (*p == ',') p++;
    }
    if (i < 13) continue;
    uint32_t raw = (uint32_t)(int64_t)v[0];
    if (first) {
      prev = raw;
      t0 = raw;
      base = 0;
      first = false;
    } else if (raw < prev)
      base += (1ULL << 32);
    prev = raw;
    uint64_t abs_us = base + raw;
    double ts = (double)(abs_us - t0) / 1e6;
    d.t.push_back(ts);
    for (int k = 0; k < 12; k++) d.s[k].push_back(v[1 + k]);
    rows++;
  }
  if (d.t.empty()) throw std::runtime_error("no usable rows in " + path);

  size_t n = d.t.size();
  for (int m = 0; m < algo::MPUS; m++) {
    d.gmag[m].resize(n);
    int b = m * 6;
    for (size_t i = 0; i < n; i++) {
      double gx = d.s[b + 3][i], gy = d.s[b + 4][i], gz = d.s[b + 5][i];
      d.gmag[m][i] = std::sqrt(gx * gx + gy * gy + gz * gz);
    }
  }
  fprintf(
      stderr,
      "loaded %ld rows, [%.3f .. %.3f] s\n",
      rows,
      d.t.front(),
      d.t.back()
  );
  return d;
}

/**
 * Draw the branding tile that fills the first cell of the grid.
 *
 * Laying the word out letter by letter is what allows tracking, which
 * `canvas::text` cannot apply; the total width is measured by the same loop
 * beforehand so the result can be centred exactly. One tracking step is
 * removed from that total because the loop appends a gap after the last
 * letter as well as between letters.
 *
 * @returns an opaque tile of `canvas::TILE` square
 * @exceptsafe basic
 */
static canvas::Image make_logo_tile() {
  canvas::Image im;
  im.px.assign(size_t(canvas::TILE) * canvas::TILE * 3, 0);
  canvas::fill_rect(im, 0, 0, canvas::TILE, canvas::TILE, canvas::BG);
  canvas::border(im);

  const std::string logo = LOGO_TEXT;
  const int glyph_h = LOGO_SCALE * canvas::PX_PER_SCALE;
  const int tracking = int(LOGO_TRACKING * float(glyph_h));

  int total = -tracking;
  for (char c : logo)
    total += canvas::text_width(std::string(1, c), LOGO_SCALE) + tracking;

  int x = (canvas::TILE - total) / 2;
  const int y = (canvas::TILE - glyph_h) / 2;
  for (char c : logo) {
    const std::string one(1, c);
    canvas::text(im, x, y, one, canvas::WHITE, LOGO_SCALE);
    x += canvas::text_width(one, LOGO_SCALE) + tracking;
  }
  return im;
}

/**
 * Draw the results tile: the headline figure over the per-algorithm table.
 *
 * Rows are ordered with the accepted estimates first, ascending by lag, so
 * abstentions collect at the bottom where they read as a footnote rather
 * than as a wall of failures. Every row is set in the same colour: an
 * abstention is already legible as one because its value column reads
 * "rejected", and the tile makes no claim about which estimates disagree
 * with the consensus. The MAD under the headline figure is where a reader
 * judges the spread.
 *
 * @param[in] ar   every algorithm result, accepted or not
 * @param[in] med  median lag across the accepted results, in milliseconds
 * @param[in] mad  spread of those results, in milliseconds
 * @returns the rendered tile
 * @exceptsafe basic
 */
static canvas::Image make_table_tile(
    const std::vector<algo::Result>& ar,
    double med,
    double mad
) {
  canvas::Image im;
  im.px.assign(size_t(canvas::TILE) * canvas::TILE * 3, 0);
  canvas::fill_rect(im, 0, 0, canvas::TILE, canvas::TILE, canvas::BG);
  canvas::border(im);
  canvas::text(im, 16, 12, "RESULTS", canvas::WHITE, canvas::TITLE_SCALE);

  char head[32];
  snprintf(head, sizeof head, "%.1f ms", med);
  const int hw = canvas::text_width(head, RESULT_SCALE);
  canvas::text(
      im,
      (canvas::TILE - hw) / 2,
      56,
      head,
      canvas::WHITE,
      RESULT_SCALE
  );

  char sub[80];
  snprintf(
      sub,
      sizeof sub,
      "MEDIAN OF %d ALGORITHMS, MAD %.1f ms",
      (int)ar.size(),
      mad
  );
  canvas::text(
      im,
      (canvas::TILE - canvas::text_width(sub, canvas::LABEL_SCALE)) / 2,
      220,
      sub,
      canvas::DIM,
      canvas::LABEL_SCALE
  );

  std::vector<const algo::Result*> rows;
  for (const algo::Result& a : ar) rows.push_back(&a);
  std::sort(
      rows.begin(),
      rows.end(),
      [](const algo::Result* a, const algo::Result* b) {
        if (a->ok != b->ok) return a->ok;
        return a->tau_ms < b->tau_ms;
      }
  );

  const int top = 253;
  const int bottom = canvas::TILE - 27;
  const int x = canvas::PAD_L;
  const int w = canvas::TILE - 2 * canvas::PAD_L;
  canvas::fill_rect(im, x, top, w, bottom - top, canvas::PANEL);
  canvas::rect(im, x, top, w, bottom - top, canvas::GRID);

  const int step = (bottom - top) / int(rows.size() + 2);
  int y = top + step;
  canvas::text(im, x + 20, y, "ALGORITHM", canvas::DIM, canvas::LABEL_SCALE);
  {
    const std::string h = "LATENCY";
    canvas::text(
        im,
        x + w - canvas::text_width(h, canvas::LABEL_SCALE) - 20,
        y,
        h,
        canvas::DIM,
        canvas::LABEL_SCALE
    );
  }
  y += step;
  for (size_t i = 0; i < rows.size(); i++) {
    const algo::Result* a = rows[i];
    char n[80];
    snprintf(n, sizeof n, "%2d   %s", int(i + 1), a->name.c_str());
    canvas::text(im, x + 20, y, n, canvas::WHITE, canvas::LABEL_SCALE);
    char v[32];
    if (a->ok)
      snprintf(v, sizeof v, "%.1f ms", a->tau_ms);
    else
      snprintf(v, sizeof v, "rejected");
    canvas::text(
        im,
        x + w - canvas::text_width(v, canvas::LABEL_SCALE) - 20,
        y,
        v,
        canvas::WHITE,
        canvas::LABEL_SCALE
    );
    y += step;
  }
  return im;
}

/**
 * Composite every tile into one image and write it out as a JPEG.
 *
 * The grid is two columns wide: the logo tile, then the results table tile,
 * then one 800x800 tile per algorithm in registration order, so a reader
 * meets the answer before the ten arguments for it. A tile whose pixel
 * buffer is not exactly one `canvas::TILE` square is skipped instead of
 * copied, so an algorithm that returned an unrendered image leaves a
 * background-coloured hole rather than shifting or corrupting its
 * neighbours.
 *
 * @param[in] ar    algorithm results supplying every tile after the first
 *                  two
 * @param[in] med   median lag in milliseconds, for the results tile
 * @param[in] mad   spread of the accepted results, in milliseconds
 * @param[in] path  destination JPEG path
 * @throws std::runtime_error if the JPEG cannot be written
 * @exceptsafe basic
 */
static void stitch_algo_tiles(
    const std::vector<algo::Result>& ar,
    double med,
    double mad,
    const std::string& path
) {
  std::vector<canvas::Image> tiles;
  tiles.push_back(make_logo_tile());
  tiles.push_back(make_table_tile(ar, med, mad));
  for (const algo::Result& a : ar) tiles.push_back(a.image);

  const int cols = 2;
  const int rows = (int(tiles.size()) + cols - 1) / cols;
  const int W = cols * canvas::TILE;
  const int H = rows * canvas::TILE;
  std::vector<unsigned char> out((size_t)W * H * 3);
  for (size_t i = 0; i + 2 < out.size(); i += 3) {
    out[i] = canvas::BG.r;
    out[i + 1] = canvas::BG.g;
    out[i + 2] = canvas::BG.b;
  }
  for (size_t k = 0; k < tiles.size(); k++) {
    const int cx = int(k % cols) * canvas::TILE;
    const int cy = int(k / cols) * canvas::TILE;
    const canvas::Image& im = tiles[k];
    if (im.px.size() != size_t(canvas::TILE) * canvas::TILE * 3) continue;
    for (int y = 0; y < canvas::TILE; y++)
      memcpy(
          &out[((size_t)(cy + y) * W + cx) * 3],
          &im.px[(size_t)y * canvas::TILE * 3],
          (size_t)canvas::TILE * 3
      );
  }
  stbi_flip_vertically_on_write(0);
  if (stbi_write_jpg(path.c_str(), W, H, 3, out.data(), JPEG_QUALITY))
    fprintf(
        stderr,
        "wrote %s (%dx%d, %d tiles)\n",
        path.c_str(),
        W,
        H,
        (int)tiles.size()
    );
  else
    throw std::runtime_error("cannot write " + path);
}

/**
 * Measure one recording end to end and write its report image.
 *
 * The sample period handed to the estimators is the mean over the whole
 * recording rather than the nominal 1 ms, because a sensor running slightly
 * off its 1 kHz target would otherwise bias every lag in proportion. All
 * ten algorithms then measure the same input independently and the reported
 * figure is the median of those that accepted, with the MAD as its spread:
 * no single method decides the answer, and one outlier cannot move a median
 * far. A positive figure means the robot hand lags the operator hand.
 *
 * @param[in] path  input CSV recording
 * @param[in] out   destination JPEG path for the stitched report
 * @throws std::runtime_error if the recording cannot be read or the report
 *         cannot be written
 * @exceptsafe basic
 */
static void analyze(
    const std::string& path,
    const std::string& out
) {
  const Data d = load_csv(path);

  algo::Input ain;
  ain.t = d.t;
  for (int m = 0; m < algo::MPUS; m++) {
    ain.gmag[m] = d.gmag[m];
  }
  ain.dt = (d.t.size() > 1)
               ? (d.t.back() - d.t.front()) / double(d.t.size() - 1)
               : 0.001;
  const std::vector<algo::Result> AR = algo::run_all(ain);
  const double algo_median_ms = algo::median_tau_ms(AR);
  const double algo_mad_ms = algo::mad_tau_ms(AR);
  printf(
      "median over %d/%d algorithms: %.1f ms (MAD %.1f ms)\n",
      algo::accepted(AR),
      (int)AR.size(),
      algo_median_ms,
      algo_mad_ms
  );
  for (const algo::Result& a : AR)
    printf(
        "  %-26s %8.1f ms  %-8s %s\n",
        a.name.c_str(),
        a.tau_ms,
        a.ok ? "ok" : "REJECT",
        a.note.c_str()
    );
  fflush(stdout);
  {
    stitch_algo_tiles(AR, algo_median_ms, algo_mad_ms, out);
  }
}

/**
 * Parse the command line, run the analysis, and report any failure.
 *
 * Written as a function-try-block so that every `std::runtime_error` raised
 * anywhere below arrives at one handler, which prints the message to stderr
 * and exits 1. Nothing is caught deeper down, and that is deliberate: it is
 * what keeps the measurement code free of status plumbing and leaves
 * `Result::ok` free to mean abstention rather than error.
 *
 * @param[in] argc  argument count, which must be exactly 3
 * @param[in] argv  program name, then input CSV, then output JPEG
 * @returns 0 on success, 1 on wrong usage or any failure
 * @exceptsafe no-throw
 */
int main(
    int argc,
    char** argv
) try {
  if (argc != 3) {
    usage(argv[0]);
    return 1;
  }
  analyze(argv[1], argv[2]);
  return 0;
} catch (const std::exception& failure) {
  fprintf(stderr, "%s\n", failure.what());
  return 1;
}
