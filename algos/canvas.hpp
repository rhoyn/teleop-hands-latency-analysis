#pragma once

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

#include <stb_truetype.h>

namespace canvas {

/**
 * Fix the tile geometry, text scaling and the two-colour palette.
 *
 * Every algorithm renders into one 800x800 px tile: an 80 px header followed
 * by exactly `PANELS` stacked plots, each one a 20 px label strip above a
 * 181 px plot box with a 29 px gap beneath it. Those numbers are chosen so
 * that three panels plus the header fill the tile with a hair of slack at the
 * bottom, and `panel(i)` is the only thing that should ever add them up.
 * Text is sized in integer steps of `PX_PER_SCALE` px so glyphs land on whole
 * pixels and survive the halving a contact sheet applies.
 *
 * The palette is deliberately two colours: the green background `BG` and
 * white. `GRID`, `DIM`, `PALE`, `SAGE`, `VIVID`, `SLATE`, `SHADE` and `INK`
 * are all white aliases, kept only so a call site can still name its intent;
 * giving one of them its own value breaks the rule. `BLACK` is the single
 * exception and exists so the second trace of a paired plot can be told
 * apart from the first.
 */
inline constexpr int PX_PER_SCALE =
    9;  ///< glyph height in px for one scale step

/**
 * Hold the memory-resident TrueType face every glyph in a tile comes from.
 *
 * `info` points into `data`, so the byte vector has to outlive it and must
 * not be reallocated; the face is therefore only ever handed out by `font()`
 * as a reference to a function-local static. When the asset cannot be found
 * `ok` stays false and all text drawing degrades to nothing at all, which is
 * preferred over aborting a long analysis run.
 */
struct Font {
  std::vector<unsigned char> data;  ///< raw .ttf bytes that `info` indexes into
  stbtt_fontinfo info;              ///< stb_truetype face parsed from `data`
  bool ok = false;                  ///< false when the font could not be loaded
};

/**
 * Return the process-wide font, loading assets/JetBrainsMono.ttf on first use.
 *
 * The asset path is resolved from /proc/self/exe, then from the executable's
 * parent directory and the working directory, so tiles look the same whether
 * the binary runs from the build tree or from an install prefix. The load is
 * attempted exactly once: on failure a warning goes to stderr and the
 * returned `Font` keeps `ok` false forever after.
 *
 * @returns the shared `Font`, always a valid object but possibly not `ok`
 * @exceptsafe basic
 */
Font& font();

inline constexpr int TILE =
    800;  ///< tile edge in px; tiles are square by contract
inline constexpr int PAD_L = 100;  ///< left gutter, sized for rotated y labels
inline constexpr int PAD_R = PAD_L;  ///< mirrored so the plot box stays centred
inline constexpr int HEADER_H = 80;  ///< title and latency band above panel 0
inline constexpr int PANEL_LABEL_H = 20;  ///< prose strip above each plot box
inline constexpr int PANEL_H =
    181;  ///< plot box height; 3 of them fill the tile
inline constexpr int PANEL_GAP = 29;  ///< room below a panel's x tick labels
inline constexpr int PANELS =
    3;  ///< panels per tile; `panel(i)` assumes i < this
inline constexpr int TITLE_SCALE = 4;  ///< 36 px header text
inline constexpr int LABEL_SCALE =
    2;  ///< 18 px tick, axis, legend and strip text

/**
 * Carry one opaque 8-bit-per-channel colour.
 *
 * There is no alpha channel: transparency only ever enters through the
 * coverage byte handed to `blend`, which comes from a rasterised glyph.
 * Channel order matches the interleaved layout of `Image::px`.
 */
struct RGB {
  unsigned char r, g, b;
};

inline constexpr RGB BG = {
    46,
    111,
    64
};  ///< tile background; the only non-white ink
inline constexpr RGB PANEL = {
    46,
    111,
    64
};  ///< plot box fill; same as `BG` by design
inline constexpr RGB GRID = {
    255,
    255,
    255
};  ///< gridlines, ticks and every border
inline constexpr RGB WHITE = {255, 255, 255};  ///< the palette's second colour
inline constexpr RGB DIM = {
    255,
    255,
    255
};  ///< white alias for de-emphasised labels
inline constexpr RGB PALE = {
    255,
    255,
    255
};  ///< white alias for call-site intent
inline constexpr RGB SAGE = {
    255,
    255,
    255
};  ///< white alias for call-site intent
inline constexpr RGB VIVID = {
    255,
    255,
    255
};  ///< white alias for call-site intent
inline constexpr RGB SLATE = {
    255,
    255,
    255
};  ///< white alias for call-site intent
inline constexpr RGB SHADE = {
    255,
    255,
    255
};  ///< white alias; not actually dimmed
inline constexpr RGB INK = {
    255,
    255,
    255
};  ///< white alias used by guide defaults
inline constexpr RGB BLACK = {
    0,
    0,
    0
};  ///< only non-palette ink: paired second trace

/**
 * Own the RGB pixel buffer a tile is composed into.
 *
 * `px` is interleaved 8-bit RGB, row major, `w * h * 3` bytes, and has to be
 * sized by the caller before any drawing happens; `standard()` is what
 * normally does that. The clip box lives on the image instead of being an
 * argument to every primitive so that `Clip` can narrow it for the length of
 * a scope. `put` honours it, but `blend` does not, so text is never cropped
 * by a panel it was drawn inside.
 */
struct Image {
  int w = TILE, h = TILE;         ///< pixel size; defaults to a whole tile
  std::vector<unsigned char> px;  ///< interleaved RGB, `w * h * 3` bytes
  int clip_x0 = 0, clip_y0 = 0, clip_x1 = TILE,
      clip_y1 = TILE;  ///< clip box in px
};

/**
 * Narrow an `Image`'s clip box for the current scope, then restore it.
 *
 * Plot primitives emit whole shapes and lean on clipping to keep them inside
 * their panel, so this guard is constructed on entry to `polyline`, `marker`,
 * `vline` and `hline`. It saves the previous box by value, which makes
 * nesting safe, and it holds a reference to the image, so it must never
 * outlive the image it was built from.
 */
struct Clip {
  Image& im;           ///< image whose clip box is being borrowed
  int x0, y0, x1, y1;  ///< the clip box to put back on destruction
  /**
   * Install a new clip box, remembering the one it replaces.
   *
   * The box is half open on both axes: a pixel survives when its coordinate is
   * at least the low bound and strictly below the high bound. Nothing is
   * validated, so an inverted box simply clips everything away.
   *
   * @param[in,out] i    image whose clip box is being narrowed
   * @param[in]     cx0  left edge in tile pixels, inclusive
   * @param[in]     cy0  top edge in tile pixels, inclusive
   * @param[in]     cx1  right edge in tile pixels, exclusive
   * @param[in]     cy1  bottom edge in tile pixels, exclusive
   * @exceptsafe no-throw
   */
  Clip(
      Image& i,
      int cx0,
      int cy0,
      int cx1,
      int cy1
  )
      : im(i),
        x0(i.clip_x0),
        y0(i.clip_y0),
        x1(i.clip_x1),
        y1(i.clip_y1) {
    im.clip_x0 = cx0;
    im.clip_y0 = cy0;
    im.clip_x1 = cx1;
    im.clip_y1 = cy1;
  }
  /**
   * Restore the clip box that was in force before this guard was built.
   *
   * @exceptsafe no-throw
   */
  ~Clip() {
    im.clip_x0 = x0;
    im.clip_y0 = y0;
    im.clip_x1 = x1;
    im.clip_y1 = y1;
  }
};

/**
 * Map data coordinates onto one panel's pixel box.
 *
 * The pixel fields describe the plot box alone, not the tile; `sx` and `sy`
 * combine them with the data range to place a sample, and `sy` flips the y
 * axis so that larger values sit higher. The defaults describe panel 0 over
 * a unit data range, so a fresh `Axes` is usable before `set_range` runs,
 * and every panel shares the same width and gutters, which is what keeps the
 * three plots of a tile visually aligned.
 */
struct Axes {
  int px = PAD_L,
      py = HEADER_H + PANEL_LABEL_H;  ///< plot box top-left in tile px
  int pw = TILE - PAD_L - PAD_R, ph = PANEL_H;  ///< plot box size in pixels
  double x0 = 0, x1 = 1, y0 = 0, y1 = 1;  ///< data range mapped onto the box
};

/**
 * Return the top of panel `i`'s label strip, in tile pixels.
 *
 * Panels stack at a fixed pitch below the header, label strip first, so the
 * plot box of the same index starts `PANEL_LABEL_H` further down. The index
 * is not range checked; anything at or above `PANELS` yields a row that falls
 * off the bottom of the tile.
 *
 * @param[in] i  zero-based panel index
 * @returns the y of the panel's label strip in tile pixels
 * @exceptsafe no-throw
 */
inline int panel_top(int i) {
  return HEADER_H + i * (PANEL_LABEL_H + PANEL_H + PANEL_GAP);
}

/**
 * Build the `Axes` covering panel `i`'s plot box.
 *
 * Only the vertical origin varies between panels; the width, the height and
 * the gutters are shared, which is what makes a tile's three y axes line up.
 * The data range is left at the default unit square, so a caller has to
 * follow up with `set_range`, `autoscale`, or simply use `draw_panel`.
 *
 * @param[in] i  zero-based panel index, expected below `PANELS`
 * @returns axes positioned over panel `i`, with a meaningless data range
 * @exceptsafe basic
 */
inline Axes panel(int i) {
  Axes ax;
  ax.py = panel_top(i) + PANEL_LABEL_H;
  return ax;
}

/**
 * Write one opaque pixel, honouring the image bounds and the clip box.
 *
 * This is the only primitive that stores opaque colour into `Image::px`, so
 * every shape in this header inherits its clipping for free. Writes outside
 * the image or the clip box are dropped rather than wrapped, which is what
 * lets shape generators emit geometry without pre-clipping it.
 *
 * @param[in,out] im  target image
 * @param[in]     x   column in tile pixels
 * @param[in]     y   row in tile pixels
 * @param[in]     c   colour to store
 * @exceptsafe no-throw
 */
inline void put(
    Image& im,
    int x,
    int y,
    RGB c
) {
  if (x < 0 || y < 0 || x >= im.w || y >= im.h) return;
  if (x < im.clip_x0 || y < im.clip_y0 || x >= im.clip_x1 || y >= im.clip_y1)
    return;
  const size_t i = (size_t(y) * im.w + x) * 3;
  im.px[i] = c.r;
  im.px[i + 1] = c.g;
  im.px[i + 2] = c.b;
}

/**
 * Fill a solid axis-aligned rectangle.
 *
 * The rectangle covers `[x, x + w)` by `[y, y + h)`, so a zero or negative
 * extent draws nothing. Every pixel goes through `put`, so the clip box
 * applies and the rectangle may end up drawn only in part.
 *
 * @param[in,out] im  target image
 * @param[in]     x   left edge in tile pixels
 * @param[in]     y   top edge in tile pixels
 * @param[in]     w   width in pixels
 * @param[in]     h   height in pixels
 * @param[in]     c   fill colour
 * @exceptsafe no-throw
 */
inline void fill_rect(
    Image& im,
    int x,
    int y,
    int w,
    int h,
    RGB c
) {
  for (int j = y; j < y + h; j++)
    for (int i = x; i < x + w; i++) put(im, i, j, c);
}

/**
 * Stroke the one-pixel outline of an axis-aligned rectangle.
 *
 * The outline lies inside the same half-open box `fill_rect` would fill, so
 * an outline and a fill with identical arguments never disagree by a pixel.
 * Corner pixels are written twice, which is harmless because drawing here is
 * opaque rather than blended.
 *
 * @param[in,out] im  target image
 * @param[in]     x   left edge in tile pixels
 * @param[in]     y   top edge in tile pixels
 * @param[in]     w   width in pixels
 * @param[in]     h   height in pixels
 * @param[in]     c   stroke colour
 * @exceptsafe no-throw
 */
inline void rect(
    Image& im,
    int x,
    int y,
    int w,
    int h,
    RGB c
) {
  for (int i = x; i < x + w; i++) {
    put(im, i, y, c);
    put(im, i, y + h - 1, c);
  }
  for (int j = y; j < y + h; j++) {
    put(im, x, j, c);
    put(im, x + w - 1, j, c);
  }
}

/**
 * Stroke a one-pixel Bresenham line between two points.
 *
 * The integer error formulation avoids floating point entirely, so a line is
 * reproducible bit for bit between machines, which matters because tiles from
 * different runs get compared. Both endpoints are drawn, nothing is
 * anti-aliased, and clipping happens per pixel inside `put`.
 *
 * @param[in,out] im  target image
 * @param[in]     x0  start column in tile pixels
 * @param[in]     y0  start row in tile pixels
 * @param[in]     x1  end column in tile pixels
 * @param[in]     y1  end row in tile pixels
 * @param[in]     c   stroke colour
 * @exceptsafe no-throw
 */
inline void line(
    Image& im,
    int x0,
    int y0,
    int x1,
    int y1,
    RGB c
) {
  int dx = std::abs(x1 - x0), sx = x0 < x1 ? 1 : -1;
  int dy = -std::abs(y1 - y0), sy = y0 < y1 ? 1 : -1;
  int err = dx + dy;
  for (;;) {
    put(im, x0, y0, c);
    if (x0 == x1 && y0 == y1) break;
    const int e2 = 2 * err;
    if (e2 >= dy) {
      err += dy;
      x0 += sx;
    }
    if (e2 <= dx) {
      err += dx;
      y0 += sy;
    }
  }
}

/**
 * Rasterise a string into an 8-bit coverage mask via stb_truetype.
 *
 * Glyphs are laid out along a single baseline using the face's advance widths
 * and composited with a max, so overlapping glyphs cannot darken each other.
 * The nominal cell height is `scale * PX_PER_SCALE` px and one pixel of slack
 * is added on each axis for glyphs that overhang their advance. Only the low
 * byte of each `char` is looked up, so input is effectively Latin-1 and
 * multi-byte UTF-8 comes out as mojibake.
 *
 * @param[in]  str       string to rasterise
 * @param[in]  scale     integer text scale, in `PX_PER_SCALE` steps
 * @param[out] w         mask width in px, zero when nothing was rasterised
 * @param[out] h         mask height in px, zero when nothing was rasterised
 * @param[out] baseline  baseline row inside the mask, measured from its top
 * @returns a row-major `w * h` coverage mask, empty when the string is empty
 *          or no font could be loaded
 * @exceptsafe basic
 */
inline std::vector<unsigned char> text_mask(
    const std::string& str,
    int scale,
    int& w,
    int& h,
    int& baseline
) {
  std::vector<unsigned char> mask;
  Font& f = font();
  w = h = baseline = 0;
  if (!f.ok || str.empty()) return mask;
  const float px = float(scale * PX_PER_SCALE);
  const float sf = stbtt_ScaleForPixelHeight(&f.info, px);
  int ascent = 0, descent = 0, gap = 0;
  stbtt_GetFontVMetrics(&f.info, &ascent, &descent, &gap);
  baseline = int(float(ascent) * sf);
  h = int(float(ascent - descent) * sf) + 2;
  float total = 0;
  for (char ch : str) {
    int aw = 0, lsb = 0;
    stbtt_GetCodepointHMetrics(&f.info, (unsigned char)ch, &aw, &lsb);
    total += float(aw) * sf;
  }
  w = int(total) + 2;
  if (w <= 0 || h <= 0) return mask;
  mask.assign(size_t(w) * h, 0);
  float xpos = 0;
  for (char ch : str) {
    const int c = (unsigned char)ch;
    int aw = 0, lsb = 0;
    stbtt_GetCodepointHMetrics(&f.info, c, &aw, &lsb);
    int gw = 0, gh = 0, gx = 0, gy = 0;
    unsigned char* bmp =
        stbtt_GetCodepointBitmap(&f.info, sf, sf, c, &gw, &gh, &gx, &gy);
    if (bmp) {
      for (int j = 0; j < gh; j++)
        for (int i = 0; i < gw; i++) {
          const int dx = int(xpos) + gx + i;
          const int dy = baseline + gy + j;
          if (dx < 0 || dy < 0 || dx >= w || dy >= h) continue;
          unsigned char& d = mask[size_t(dy) * w + dx];
          d = std::max(d, bmp[j * gw + i]);
        }
      stbtt_FreeBitmap(bmp, nullptr);
    }
    xpos += float(aw) * sf;
  }
  return mask;
}

/**
 * Composite a colour over one pixel using an 8-bit coverage value.
 *
 * Coverage below 16 is discarded, which stops the faintest glyph fringes from
 * muddying the flat background at the cost of very slightly harder text
 * edges. Unlike `put` this ignores the image's clip box and checks only the
 * image bounds, so text drawn inside a narrowed scope still escapes it.
 *
 * @param[in,out] im  target image
 * @param[in]     x   column in tile pixels
 * @param[in]     y   row in tile pixels
 * @param[in]     a   coverage, 0 for transparent through 255 for opaque
 * @param[in]     c   colour to composite over the existing pixel
 * @exceptsafe no-throw
 */
inline void blend(
    Image& im,
    int x,
    int y,
    unsigned char a,
    RGB c
) {
  if (a < 16 || x < 0 || y < 0 || x >= im.w || y >= im.h) return;
  if (x < im.clip_x0 || y < im.clip_y0 || x >= im.clip_x1 || y >= im.clip_y1)
    return;
  const size_t o = (size_t(y) * im.w + x) * 3;
  const float t = float(a) / 255.0f;
  im.px[o] = (unsigned char)(im.px[o] * (1 - t) + c.r * t);
  im.px[o + 1] = (unsigned char)(im.px[o + 1] * (1 - t) + c.g * t);
  im.px[o + 2] = (unsigned char)(im.px[o + 2] * (1 - t) + c.b * t);
}

/**
 * Draw a string with the top-left corner of its mask at a given point.
 *
 * The anchor is the corner of the rasterised block rather than the baseline,
 * so vertical centring is done by subtracting half of `text_height`. Every
 * call rasterises afresh because there is no glyph cache, which is affordable
 * given a tile carries only a few dozen short labels.
 *
 * @param[in,out] im     target image
 * @param[in]     x      left edge of the text block in tile pixels
 * @param[in]     y      top edge of the text block in tile pixels
 * @param[in]     str    string to draw
 * @param[in]     c      text colour
 * @param[in]     scale  integer text scale, in `PX_PER_SCALE` steps
 * @exceptsafe basic
 */
inline void draw_text(
    Image& im,
    int x,
    int y,
    const std::string& str,
    RGB c,
    int scale
) {
  int w = 0, h = 0, bl = 0;
  const std::vector<unsigned char> m = text_mask(str, scale, w, h, bl);
  for (int j = 0; j < h; j++)
    for (int i = 0; i < w; i++)
      blend(im, x + i, y + j, m[size_t(j) * w + i], c);
}

/**
 * Blit a string rotated a quarter turn anticlockwise, reading upwards.
 *
 * This is what places a y-axis title alongside a panel. The mask from
 * `text_mask` is sampled transposed, so `x` is the left edge of the resulting
 * column and `y` is its bottom edge: the text grows upwards from `y` and
 * rightwards from `x`. A caller must therefore reserve `text_width` of
 * vertical room and `text_height` of horizontal room, which is exactly what
 * `axis_labels` does.
 *
 * @param[in,out] im     target image
 * @param[in]     x      left edge of the rotated block in tile pixels
 * @param[in]     y      bottom edge of the rotated block in tile pixels
 * @param[in]     str    string to draw
 * @param[in]     c      text colour
 * @param[in]     scale  integer text scale, in `PX_PER_SCALE` steps
 * @exceptsafe basic
 */
inline void draw_text_rot90(
    Image& im,
    int x,
    int y,
    const std::string& str,
    RGB c,
    int scale
) {
  int w = 0, h = 0, bl = 0;
  const std::vector<unsigned char> m = text_mask(str, scale, w, h, bl);
  for (int j = 0; j < h; j++)
    for (int i = 0; i < w; i++)
      blend(im, x + j, y - i, m[size_t(j) * w + i], c);
}

/**
 * Measure a string's advance width in pixels without rasterising it.
 *
 * Layout code centres and right-aligns labels with this, so it has to agree
 * with `text_mask`, which sums the same advances; it omits only that
 * function's one pixel of slack. With no font loaded it falls back to half an
 * em per character so that a text-free tile still lays out plausibly.
 *
 * @param[in] str    string to measure
 * @param[in] scale  integer text scale, in `PX_PER_SCALE` steps
 * @returns the advance width in pixels, rounded down
 * @exceptsafe basic
 */
inline int text_width(
    const std::string& str,
    int scale
) {
  Font& f = font();
  if (!f.ok) return int(str.size()) * scale * PX_PER_SCALE / 2;
  const float sf =
      stbtt_ScaleForPixelHeight(&f.info, float(scale * PX_PER_SCALE));
  float w = 0;
  for (char ch : str) {
    int aw = 0, lsb = 0;
    stbtt_GetCodepointHMetrics(&f.info, (unsigned char)ch, &aw, &lsb);
    w += float(aw) * sf;
  }
  return int(w);
}

/**
 * Draw a string at a point; the short spelling used throughout this header.
 *
 * Kept as a thin alias so that plotting code reads as a run of drawing verbs
 * rather than of implementation names, and so the rasterising back end can be
 * swapped in one place.
 *
 * @param[in,out] im   target image
 * @param[in]     x    left edge of the text block in tile pixels
 * @param[in]     y    top edge of the text block in tile pixels
 * @param[in]     str  string to draw
 * @param[in]     c    text colour
 * @param[in]     s    integer text scale, in `PX_PER_SCALE` steps
 * @exceptsafe basic
 */
inline void text(
    Image& im,
    int x,
    int y,
    const std::string& str,
    RGB c,
    int s
) {
  draw_text(im, x, y, str, c, s);
}

/**
 * Draw a string rotated a quarter turn; the short spelling for axis titles.
 *
 * The alias exists for the same reason as `text`, and it inherits the anchor
 * convention of `draw_text_rot90`: `y` is the bottom of the column and the
 * glyphs climb from there.
 *
 * @param[in,out] im   target image
 * @param[in]     x    left edge of the rotated block in tile pixels
 * @param[in]     y    bottom edge of the rotated block in tile pixels
 * @param[in]     str  string to draw
 * @param[in]     c    text colour
 * @param[in]     s    integer text scale, in `PX_PER_SCALE` steps
 * @exceptsafe basic
 */
inline void text_rot(
    Image& im,
    int x,
    int y,
    const std::string& str,
    RGB c,
    int s
) {
  draw_text_rot90(im, x, y, str, c, s);
}

/**
 * Return the nominal line height of a text scale, in pixels.
 *
 * This is the requested cell height, not the height `text_mask` reports for
 * one particular string; layout uses it so that stacked rows of labels keep
 * an even pitch whether or not their glyphs happen to have descenders.
 *
 * @param[in] scale  integer text scale, in `PX_PER_SCALE` steps
 * @returns the line height in pixels
 * @exceptsafe no-throw
 */
inline int text_height(int scale) { return scale * PX_PER_SCALE; }

/**
 * Draw the tile's two-pixel outer border.
 *
 * Two nested one-pixel rectangles are used instead of one thicker stroke
 * because the border has to survive the halving a contact sheet applies,
 * where a single pixel would alias away on some rows.
 *
 * @param[in,out] im  tile-sized image to frame
 * @exceptsafe no-throw
 */
inline void border(Image& im) {
  rect(im, 0, 0, TILE, TILE, GRID);
  rect(im, 1, 1, TILE - 2, TILE - 2, GRID);
}

/**
 * Build a blank tile: background, border, header text and empty panel boxes.
 *
 * This is where every algorithm's figure starts. It allocates the pixel
 * buffer, paints the background, writes the title and the headline latency,
 * and outlines all three plot boxes so that the layout is identical across
 * tiles even when a panel is later left empty. The title is truncated at 30
 * characters to keep it inside the header, and a rejected estimate prints the
 * word "rejected" rather than a misleading number.
 *
 * @param[in] title   algorithm name, truncated to 30 characters
 * @param[in] tau_ms  estimated latency in ms, printed to a tenth
 * @param[in] ok      false when the estimate was rejected and must not print
 * @returns a fully painted tile that carries no data yet
 * @exceptsafe basic
 */
inline Image standard(
    const std::string& title,
    double tau_ms,
    bool ok
) {
  Image im;
  im.px.assign(size_t(TILE) * TILE * 3, 0);
  fill_rect(im, 0, 0, TILE, TILE, BG);
  border(im);
  text(im, 16, 12, title.substr(0, 30), WHITE, TITLE_SCALE);
  char buf[48];
  if (ok)
    snprintf(buf, sizeof buf, "%.1f ms", tau_ms);
  else
    snprintf(buf, sizeof buf, "rejected");
  text(im, 16, 44, buf, ok ? WHITE : SHADE, TITLE_SCALE);
  for (int i = 0; i < PANELS; i++) {
    const Axes ax = panel(i);
    fill_rect(im, ax.px, ax.py, ax.pw, ax.ph, PANEL);
    rect(im, ax.px, ax.py, ax.pw, ax.ph, GRID);
  }
  return im;
}

/**
 * Write the explanatory strip above a panel's plot box.
 *
 * The strip says what the panel shows and why, so it carries prose rather
 * than a terse title. It is truncated at 60 characters, which is roughly what
 * fits across the tile at `LABEL_SCALE`.
 *
 * @param[in,out] im  target tile
 * @param[in]     i   zero-based panel index
 * @param[in]     s   strip text, truncated to 60 characters
 * @exceptsafe basic
 */
inline void panel_label(
    Image& im,
    int i,
    const std::string& s
) {
  text(im, PAD_L, panel_top(i) + 2, s.substr(0, 60), DIM, LABEL_SCALE);
}

/**
 * Set the data range a panel's axes map onto its pixel box.
 *
 * Degenerate and inverted spans are widened to one unit instead of being
 * rejected, because a constant signal is an ordinary input here and a zero
 * span would make `sx` and `sy` divide by zero. Axes are never flipped: `x0`
 * and `y0` always end up at the left edge and the bottom edge.
 *
 * @param[in,out] ax  axes to modify
 * @param[in]     x0  lower x bound in data units
 * @param[in]     x1  upper x bound, pushed above `x0` when it is not
 * @param[in]     y0  lower y bound in data units
 * @param[in]     y1  upper y bound, pushed above `y0` when it is not
 * @exceptsafe no-throw
 */
inline void set_range(
    Axes& ax,
    double x0,
    double x1,
    double y0,
    double y1
) {
  if (x1 <= x0) x1 = x0 + 1;
  if (y1 <= y0) y1 = y0 + 1;
  ax.x0 = x0;
  ax.x1 = x1;
  ax.y0 = y0;
  ax.y1 = y1;
}

/**
 * Fit the axes to one series, padding the y range only.
 *
 * The x range is taken tight because it is nearly always time, where an empty
 * margin reads as missing data; y gets six percent of its span added at each
 * end so that extremes do not sit on the panel border. Either vector being
 * empty leaves the axes untouched rather than inventing a range, and the two
 * vectors are scanned independently, so a length mismatch is not detected.
 *
 * @param[in,out] ax  axes to fit
 * @param[in]     x   x samples, scanned for their extremes
 * @param[in]     y   y samples, scanned for their extremes
 * @exceptsafe basic
 */
inline void autoscale(
    Axes& ax,
    const std::vector<double>& x,
    const std::vector<double>& y
) {
  if (x.empty() || y.empty()) return;
  const auto xm = std::minmax_element(x.begin(), x.end());
  const auto ym = std::minmax_element(y.begin(), y.end());
  const double pad = 0.06 * (*ym.second - *ym.first);
  set_range(ax, *xm.first, *xm.second, *ym.first - pad, *ym.second + pad);
}

/**
 * Project a data x value onto a column of the panel's plot box.
 *
 * The mapping spans `pw - 1` pixels so that `x0` and `x1` land exactly on the
 * first and the last column of the box. Values outside the range map outside
 * the box and are not clamped; the clip box is what keeps them off the tile.
 * The conversion truncates towards zero, so up to a pixel of jitter between
 * neighbouring samples is expected and normal.
 *
 * @param[in] ax  axes carrying the pixel box and the data range
 * @param[in] x   value in data units
 * @returns the column in tile pixels
 * @exceptsafe no-throw
 */
inline int sx(
    const Axes& ax,
    double x
) {
  return ax.px + int((x - ax.x0) / (ax.x1 - ax.x0) * (ax.pw - 1));
}

/**
 * Project a data y value onto a row of the panel's plot box.
 *
 * The row is measured up from the bottom of the box, so larger values sit
 * higher, which is the reverse of the image's own row order; this is the only
 * place that flip happens. As with `sx` the span is `ph - 1` pixels, values
 * outside the range are not clamped, and the conversion truncates towards
 * zero rather than rounding.
 *
 * @param[in] ax  axes carrying the pixel box and the data range
 * @param[in] y   value in data units
 * @returns the row in tile pixels
 * @exceptsafe no-throw
 */
inline int sy(
    const Axes& ax,
    double y
) {
  return ax.py + ax.ph - 1 - int((y - ax.y0) / (ax.y1 - ax.y0) * (ax.ph - 1));
}

/**
 * Draw a series as straight segments between consecutive samples.
 *
 * The trace is clipped to the inside of the plot box so it can never
 * overwrite the border. Long series are decimated to roughly four samples per
 * column, which is invisible at tile resolution and keeps a million-sample
 * trace cheap; because the stride simply steps by a constant, the final
 * sample can be dropped. Fewer than two samples draw nothing, and the shorter
 * of the two vectors decides the length.
 *
 * @param[in,out] im  target tile
 * @param[in]     ax  axes defining the plot box and the data range
 * @param[in]     x   x samples
 * @param[in]     y   y samples, paired with `x` by index
 * @param[in]     c   stroke colour
 * @exceptsafe basic
 */
inline void polyline(
    Image& im,
    const Axes& ax,
    const std::vector<double>& x,
    const std::vector<double>& y,
    RGB c
) {
  const Clip
      clip(im, ax.px + 1, ax.py + 1, ax.px + ax.pw - 1, ax.py + ax.ph - 1);
  const size_t n = std::min(x.size(), y.size());
  if (n < 2) return;
  const size_t step = (n > size_t(ax.pw) * 4) ? n / (size_t(ax.pw) * 4) : 1;
  int lx = sx(ax, x[0]), ly = sy(ax, y[0]);
  for (size_t i = step; i < n; i += step) {
    const int cx = sx(ax, x[i]), cy = sy(ax, y[i]);
    line(im, lx, ly, cx, cy, c);
    lx = cx;
    ly = cy;
  }
  const int ex = sx(ax, x[n - 1]), ey = sy(ax, y[n - 1]);
  line(im, lx, ly, ex, ey, c);
}

/**
 * Draw a filled circular marker at one data point.
 *
 * The disc is rasterised by a squared-radius test per pixel, so it is not
 * anti-aliased and a small radius looks a little square; that is accepted
 * because the two-colour palette has no intermediate tone to soften an edge
 * with. Unlike `vline` and `hline`, this honours the colour it is handed,
 * which is how the second trace of a paired plot gets drawn in `BLACK`.
 * Clipping to the plot box means a marker near an edge appears as a partial
 * disc rather than being dropped.
 *
 * @param[in,out] im  target tile
 * @param[in]     ax  axes defining the plot box and the data range
 * @param[in]     x   marker centre x in data units
 * @param[in]     y   marker centre y in data units
 * @param[in]     c   fill colour
 * @param[in]     r   radius in pixels
 * @exceptsafe basic
 */
inline void marker(
    Image& im,
    const Axes& ax,
    double x,
    double y,
    RGB c,
    int r
) {
  const Clip
      clip(im, ax.px + 1, ax.py + 1, ax.px + ax.pw - 1, ax.py + ax.ph - 1);
  const int cx = sx(ax, x), cy = sy(ax, y);
  for (int j = -r; j <= r; j++)
    for (int i = -r; i <= r; i++)
      if (i * i + j * j <= r * r) put(im, cx + i, cy + j, c);
}

/**
 * Draw a full-height vertical rule at a data x value.
 *
 * The colour argument is overwritten with white before anything is drawn.
 * That is intentional, and it surprises people: a rule is structure rather
 * than data, and the two-colour palette gives structure the same single ink
 * as everything else. Pass whichever constant reads best at the call site,
 * because it cannot change the output. The rule spans the full panel height
 * and is clipped to the inside of the plot box.
 *
 * @param[in,out] im  target tile
 * @param[in]     ax  axes defining the plot box and the data range
 * @param[in]     x   position in data units
 * @param[in]     c   ignored; the rule is always white
 * @exceptsafe basic
 */
inline void vline(
    Image& im,
    const Axes& ax,
    double x,
    RGB c
) {
  c = WHITE;
  const Clip
      clip(im, ax.px + 1, ax.py + 1, ax.px + ax.pw - 1, ax.py + ax.ph - 1);
  const int cx = sx(ax, x);
  for (int j = ax.py; j < ax.py + ax.ph; j++) put(im, cx, j, c);
}

/**
 * Draw a full-width horizontal rule at a data y value.
 *
 * As with `vline`, the colour argument is overwritten with white before the
 * rule is drawn; that is deliberate and is the house palette rule, not an
 * oversight. The rule spans the full panel width and is clipped to the inside
 * of the plot box so it never lands on the border.
 *
 * @param[in,out] im  target tile
 * @param[in]     ax  axes defining the plot box and the data range
 * @param[in]     y   position in data units
 * @param[in]     c   ignored; the rule is always white
 * @exceptsafe basic
 */
inline void hline(
    Image& im,
    const Axes& ax,
    double y,
    RGB c
) {
  c = WHITE;
  const Clip
      clip(im, ax.px + 1, ax.py + 1, ax.px + ax.pw - 1, ax.py + ax.ph - 1);
  const int cy = sy(ax, y);
  for (int i = ax.px; i < ax.px + ax.pw; i++) put(im, i, cy, c);
}

/**
 * Write the header subtitle to the right of the latency readout.
 *
 * The x offset is measured from a worst-case "0000.0 ms" readout rather than
 * from the digits actually printed, so the subtitle does not shuffle sideways
 * between tiles with different latencies. It is truncated at 34 characters,
 * which is what is left of the header width.
 *
 * @param[in,out] im  target tile
 * @param[in]     s   subtitle text, truncated to 34 characters
 * @exceptsafe basic
 */
inline void caption(
    Image& im,
    const std::string& s
) {
  text(
      im,
      16 + text_width("0000.0 ms", TITLE_SCALE) + 16,
      44,
      s.substr(0, 34),
      DIM,
      LABEL_SCALE
  );
}

/**
 * Format a tick value at the precision its step implies.
 *
 * The decimal count is derived from the tick spacing so that neighbouring
 * labels differ in print, then capped at three digits to keep them inside the
 * gutter. Values that would otherwise need many digits switch to scientific
 * notation with the exponent stripped of its sign padding and leading zeros,
 * so 1e-05 comes out as 1.0e-5. A value that rounds to zero from below prints
 * as "0" rather than "-0", which would look like a different tick.
 *
 * @param[in] v     value to print, in data units
 * @param[in] step  spacing between ticks, which sets the precision
 * @returns the formatted tick label
 * @exceptsafe basic
 */
inline std::string fmt_tick(
    double v,
    double step
) {
  char buf[32];
  const double m = std::max(std::fabs(v), std::fabs(step));
  if (m >= 1e5 || (m > 0 && m < 1e-3)) {
    snprintf(buf, sizeof buf, "%.1e", v);
    std::string s(buf);
    const size_t e = s.find('e');
    if (e != std::string::npos) {
      std::string mant = s.substr(0, e);
      std::string ex = s.substr(e + 1);
      const bool neg = !ex.empty() && ex[0] == '-';
      if (!ex.empty() && (ex[0] == '+' || ex[0] == '-')) ex = ex.substr(1);
      while (ex.size() > 1 && ex[0] == '0') ex = ex.substr(1);
      s = mant + "e" + (neg ? "-" : "") + ex;
    }
    return s;
  }
  const int dp = std::max(0, int(std::ceil(-std::log10(step) + 0.2)));
  snprintf(buf, sizeof buf, "%.*f", std::min(dp, 3), v);
  std::string s(buf);
  if (s == "-0" || s == "-0.0" || s == "-0.00" || s == "-0.000")
    s = s.substr(1);
  return s;
}

/**
 * Draw a panel's dotted gridlines, tick marks and tick labels.
 *
 * Ticks sit at evenly spaced fractions of the data range rather than at
 * rounded values, so their labels carry whatever precision the range needs
 * and `fmt_tick` works that out from the spacing. Gridlines are dotted every
 * fourth pixel because a solid white line at full contrast would compete with
 * the data for attention. Call this before drawing any series so that the
 * grid ends up underneath them.
 *
 * @param[in,out] im  target tile
 * @param[in]     ax  axes defining the plot box and the data range
 * @param[in]     nx  x ticks including both ends; must be greater than one
 * @param[in]     ny  y ticks including both ends; must be greater than one
 * @exceptsafe basic
 */
inline void grid_and_ticks(
    Image& im,
    const Axes& ax,
    int nx,
    int ny
) {
  const int ts = LABEL_SCALE;
  const double xs = (ax.x1 - ax.x0) / double(nx - 1);
  for (int i = 0; i < nx; i++) {
    const double v = ax.x0 + xs * double(i);
    const int x = sx(ax, v);
    for (int y = ax.py; y < ax.py + ax.ph; y += 4) put(im, x, y, GRID);
    for (int k = 0; k < 4; k++) put(im, x, ax.py + ax.ph + k, DIM);
    const std::string t = fmt_tick(v, xs);
    text(im, x - text_width(t, ts) / 2, ax.py + ax.ph + 1, t, DIM, ts);
  }
  const double ys = (ax.y1 - ax.y0) / double(ny - 1);
  for (int i = 0; i < ny; i++) {
    const double v = ax.y0 + ys * double(i);
    const int y = sy(ax, v);
    for (int x = ax.px; x < ax.px + ax.pw; x += 4) put(im, x, y, GRID);
    for (int k = 1; k <= 4; k++) put(im, ax.px - k, y, DIM);
    const std::string t = fmt_tick(v, ys);
    text(
        im,
        ax.px - 1 - text_width(t, ts),
        y - text_height(ts) / 2,
        t,
        DIM,
        ts
    );
  }
}

/**
 * Measure the widest y tick label a panel would print.
 *
 * The rotated y-axis title has to clear the tick labels, whose width depends
 * on the data range, so this regenerates exactly the strings
 * `grid_and_ticks` would draw instead of guessing from a digit count. The
 * tick count passed here must match the one used for drawing, or the gutter
 * reserved will be the wrong size.
 *
 * @param[in] ax  axes carrying the data range
 * @param[in] ny  y ticks including both ends; must be greater than one
 * @returns the widest label width in pixels at `LABEL_SCALE`
 * @exceptsafe basic
 */
inline int y_tick_width(
    const Axes& ax,
    int ny
) {
  const double ys = (ax.y1 - ax.y0) / double(ny - 1);
  int w = 0;
  for (int i = 0; i < ny; i++)
    w = std::max(
        w,
        text_width(fmt_tick(ax.y0 + ys * double(i), ys), LABEL_SCALE)
    );
  return w;
}

/**
 * Write a panel's x and y axis titles.
 *
 * The x title is centred beneath the tick labels; the y title is rotated a
 * quarter turn and pushed left of the widest tick label, measured for five
 * ticks so that it agrees with what `draw_panel` actually draws. Call it
 * after `grid_and_ticks` so the measured gutter and the drawn one match.
 * Units belong in these strings, since nothing else on the tile carries them.
 *
 * @param[in,out] im    target tile
 * @param[in]     ax    axes defining the plot box and the data range
 * @param[in]     xlab  x axis title, units included
 * @param[in]     ylab  y axis title, units included
 * @exceptsafe basic
 */
inline void axis_labels(
    Image& im,
    const Axes& ax,
    const std::string& xlab,
    const std::string& ylab
) {
  const int s = LABEL_SCALE;
  text(
      im,
      ax.px + (ax.pw - text_width(xlab, s)) / 2,
      ax.py + ax.ph + 1 + text_height(s),
      xlab,
      DIM,
      s
  );
  text_rot(
      im,
      ax.px - y_tick_width(ax, 5) - text_height(s) - 2,
      ax.py + (ax.ph + text_width(ylab, s)) / 2,
      ylab,
      DIM,
      s
  );
}

/**
 * Describe one row of a panel legend.
 *
 * Only data series earn a legend row; a straight guide line is labelled where
 * it is drawn instead, which keeps the key short enough to overlay the plot.
 * The two flags pick the swatch: a dot for a marker series, a dimension line
 * with end caps for a bracket, and a plain bar when neither is set. They are
 * not mutually exclusive in the type, but `dot` wins if both are set.
 */
struct LegendItem {
  std::string label;     ///< row text, always drawn white
  RGB color = WHITE;     ///< swatch colour, matching the series it stands for
  bool dot = false;      ///< draw a disc swatch, for a marker series
  bool bracket = false;  ///< draw a dimension-line swatch with end caps
};

/**
 * Draw a boxed key overlaying a corner of the plot box.
 *
 * The box is filled with the panel colour and outlined so that it stays
 * readable over data, and it is sized from the widest label plus room for a
 * swatch. It sits inside the plot box and therefore can hide data, so the
 * caller picks whichever side is emptier. An empty item list draws nothing
 * rather than an empty frame. Label text is always white even when the swatch
 * is not, so a `BLACK` trace still gets a legible name.
 *
 * @param[in,out] im     target tile
 * @param[in]     ax     axes defining the plot box
 * @param[in]     items  rows, drawn top to bottom in the order given
 * @param[in]     right  true to pin the box to the right edge, else the left
 * @exceptsafe basic
 */
inline void legend(
    Image& im,
    const Axes& ax,
    const std::vector<LegendItem>& items,
    bool right
) {
  if (items.empty()) return;
  const int s = LABEL_SCALE;
  const int lh = text_height(s) + 8;
  int w = 0;
  for (const LegendItem& it : items) w = std::max(w, text_width(it.label, s));
  w += 44;
  const int h = int(items.size()) * lh + 8;
  const int x = right ? ax.px + ax.pw - w - 8 : ax.px + 8;
  const int y = ax.py + 8;
  fill_rect(im, x, y, w, h, PANEL);
  rect(im, x, y, w, h, GRID);
  for (size_t i = 0; i < items.size(); i++) {
    const int ly = y + 4 + int(i) * lh;
    const int cy = ly + text_height(s) / 2;
    if (items[i].dot) {
      const int cx = x + 16;
      for (int j = -5; j <= 5; j++)
        for (int k = -5; k <= 5; k++)
          if (j * j + k * k <= 25) put(im, cx + k, cy + j, items[i].color);
    } else if (items[i].bracket) {
      fill_rect(im, x + 6, cy - 1, 22, 3, items[i].color);
      fill_rect(im, x + 6, cy - 7, 3, 15, items[i].color);
      fill_rect(im, x + 25, cy - 7, 3, 15, items[i].color);
    } else {
      fill_rect(im, x + 6, cy - 1, 22, 3, items[i].color);
    }
    text(im, x + 36, ly, items[i].label, WHITE, s);
  }
}

/**
 * Draw a guide line and label it on the plot rather than in the legend.
 *
 * Guides mark structure such as a threshold or an estimated delay, and they
 * are named where they are drawn so that a reader never has to match a colour
 * against a key; with two colours the key could not tell them apart anyway.
 * A vertical label flips to the far side of its line when it would overrun
 * the right edge of the plot box. `slot` staggers labels so that several
 * guides in one panel do not overprint each other, and the caller is expected
 * to count vertical and horizontal guides separately.
 *
 * @param[in,out] im        target tile
 * @param[in]     ax        axes defining the plot box and the data range
 * @param[in]     value     position in data units, on x or y per `vertical`
 * @param[in]     vertical  true for a vertical rule, false for a horizontal
 * @param[in]     label     text beside the line; empty draws the line bare
 * @param[in]     slot      zero-based stacking row among guides of this
 *                          orientation, which offsets the label
 * @exceptsafe basic
 */
inline void mark_guide(
    Image& im,
    const Axes& ax,
    double value,
    bool vertical,
    const std::string& label,
    int slot
) {
  if (vertical) {
    vline(im, ax, value, WHITE);
    if (label.empty()) return;
    const int cx = sx(ax, value);
    const int tw = text_width(label, LABEL_SCALE);
    int tx = cx + 6;
    if (tx + tw > ax.px + ax.pw - 4) tx = cx - 6 - tw;
    const int row = (text_height(LABEL_SCALE) + 6) * (slot + 1);
    text(im, tx, ax.py + ax.ph - row - 2, label, WHITE, LABEL_SCALE);
  } else {
    hline(im, ax, value, WHITE);
    if (label.empty()) return;
    const int cy = sy(ax, value);
    const int off = (text_height(LABEL_SCALE) + 4) * slot;
    text(
        im,
        ax.px + 8,
        cy - text_height(LABEL_SCALE) - 4 - off,
        label,
        WHITE,
        LABEL_SCALE
    );
  }
}

/**
 * Select how a `Series` is rendered inside a panel.
 *
 * `LINE` joins consecutive samples and suits a densely and regularly sampled
 * signal; `MARKERS` draws an independent disc per sample and is what to use
 * for sparse or unordered estimates, where a joining line would imply a
 * continuity that is not in the data. The choice also picks the legend
 * swatch, a bar or a dot.
 */
enum SeriesKind { LINE, MARKERS };

/**
 * Describe one data trace to be drawn in a panel.
 *
 * The sample vectors are borrowed rather than copied, so they must outlive
 * the `draw_panel` call; a null or empty pair is skipped in silence, which
 * lets a caller declare an optional trace unconditionally. Samples pair by
 * index and the shorter vector wins. Declaration order matters twice: it is
 * the legend order, and it is reversed at draw time so that the first series
 * declared ends up on top of the others.
 */
struct Series {
  std::string label;  ///< legend text; empty keeps it out of the legend
  RGB color = WHITE;  ///< `WHITE`, or `BLACK` for a paired second trace
  const std::vector<double>* x =
      nullptr;  ///< borrowed x samples; null skips series
  const std::vector<double>* y =
      nullptr;             ///< borrowed y samples, paired by index
  SeriesKind kind = LINE;  ///< joined line or one disc per sample
  int radius = 5;  ///< marker radius in px, unused when `kind` is `LINE`
};

/**
 * Describe a straight reference line drawn across a panel.
 *
 * A guide is annotation rather than data: `mark_guide` labels it on the plot
 * and it never appears in the legend. Its value is folded into an autoscaled
 * range, so a guide that falls outside the data widens the panel instead of
 * being clipped out of sight.
 */
struct Guide {
  double value = 0;      ///< position in data units, on x or y per `vertical`
  RGB color = INK;       ///< ignored; guide lines are always drawn white
  std::string label;     ///< text beside the line; empty draws the line bare
  bool vertical = true;  ///< true for a vertical rule, false for horizontal
};

/**
 * Declare everything one panel of a tile should show.
 *
 * This is the description `draw_panel` renders, so a caller assembles a panel
 * as data and never issues drawing calls of its own; that indirection is what
 * keeps every tile in house style. The range is derived from the series and
 * the guides unless `fixed_range` is set, which is how panels that have to be
 * comparable across algorithms pin their axes to the same numbers.
 */
struct Panel {
  std::string explain;  ///< prose strip above the plot box, cut at 60 chars
  std::string xlabel;   ///< x axis title, units included
  std::string ylabel;   ///< y axis title, units included
  std::vector<Series>
      series;  ///< traces; the first declared is drawn last, on top
  std::vector<Guide>
      guides;  ///< reference lines, labelled on the plot not the key
  bool fixed_range =
      false;  ///< true to use the bounds below instead of autoscaling
  double x0 = 0, x1 = 1, y0 = 0,
         y1 = 1;             ///< bounds used only when `fixed_range`
  bool legend_right = true;  ///< pin the key to the right edge of the plot box
};

/**
 * Render a whole panel from its description, in house style.
 *
 * This is the shared entry point every tile goes through, and the order of
 * work is the style: the range is either taken verbatim or derived from the
 * series and the guides together, then ten x ticks and five y ticks with
 * their gridlines go down first so that they sit under the data, then the
 * guides with their on-plot labels, then the series in reverse declaration
 * order so the first declared trace finishes on top, then the axis titles,
 * and finally a legend assembled from whichever series carry a label.
 *
 * A panel with no usable series and no fixed range returns straight after the
 * label strip, leaving an empty box rather than inventing an axis; note that
 * its guides are then not drawn either.
 *
 * @param[in,out] im     target tile
 * @param[in]     index  zero-based panel index, expected below `PANELS`
 * @param[in]     p      panel contents; the vectors its series borrow must
 *                       stay alive for the duration of the call
 * @exceptsafe basic
 */
inline void draw_panel(
    Image& im,
    int index,
    const Panel& p
) {
  panel_label(im, index, p.explain);
  Axes ax = panel(index);

  if (p.fixed_range) {
    set_range(ax, p.x0, p.x1, p.y0, p.y1);
  } else {
    bool first = true;
    double xa = 0, xb = 0, ya = 0, yb = 0;
    for (const Series& s : p.series) {
      if (!s.x || !s.y || s.x->empty() || s.y->empty()) continue;
      const auto xm = std::minmax_element(s.x->begin(), s.x->end());
      const auto ym = std::minmax_element(s.y->begin(), s.y->end());
      if (first) {
        xa = *xm.first;
        xb = *xm.second;
        ya = *ym.first;
        yb = *ym.second;
        first = false;
      } else {
        xa = std::min(xa, *xm.first);
        xb = std::max(xb, *xm.second);
        ya = std::min(ya, *ym.first);
        yb = std::max(yb, *ym.second);
      }
    }
    if (first) {
      for (const Guide& g : p.guides)
        mark_guide(im, ax, g.value, g.vertical, g.label, 0);
      return;
    }
    for (const Guide& g : p.guides) {
      if (g.vertical) {
        xa = std::min(xa, g.value);
        xb = std::max(xb, g.value);
      } else {
        ya = std::min(ya, g.value);
        yb = std::max(yb, g.value);
      }
    }
    const double pad = 0.06 * (yb - ya);
    set_range(ax, xa, xb, ya - pad, yb + pad);
  }

  grid_and_ticks(im, ax, 10, 5);

  int vslot = 0, hslot = 0;
  for (const Guide& g : p.guides)
    mark_guide(
        im,
        ax,
        g.value,
        g.vertical,
        g.label,
        g.vertical ? vslot++ : hslot++
    );

  for (auto it = p.series.rbegin(); it != p.series.rend(); ++it) {
    const Series& s = *it;
    if (!s.x || !s.y || s.x->empty() || s.y->empty()) continue;
    if (s.kind == LINE) {
      polyline(im, ax, *s.x, *s.y, s.color);
    } else {
      const size_t n = std::min(s.x->size(), s.y->size());
      for (size_t i = 0; i < n; i++)
        marker(im, ax, (*s.x)[i], (*s.y)[i], s.color, s.radius);
    }
  }

  axis_labels(im, ax, p.xlabel, p.ylabel);

  std::vector<LegendItem> keys;
  for (const Series& s : p.series)
    if (!s.label.empty())
      keys.push_back({s.label, s.color, s.kind == MARKERS, false});
  if (!keys.empty()) legend(im, ax, keys, p.legend_right);
}

/**
 * Build a finished blank tile: header, subtitle, border and panel boxes.
 *
 * This is the one-call front door an algorithm uses before filling panels
 * with `draw_panel`. It exists so that header composition lives in this
 * header rather than being reassembled, slightly differently, at every call
 * site.
 *
 * @param[in] name      algorithm name, truncated to 30 characters
 * @param[in] tau_ms    estimated latency in milliseconds
 * @param[in] ok        false when the estimate was rejected
 * @param[in] subtitle  header subtitle, truncated to 34 characters
 * @returns a painted tile whose three panels are still empty
 * @exceptsafe basic
 */
inline Image tile(
    const std::string& name,
    double tau_ms,
    bool ok,
    const std::string& subtitle
) {
  Image im = standard(name, tau_ms, ok);
  caption(im, subtitle);
  return im;
}

}
