# Visualization Improvements — Milestone 2 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add optional rasterization of dense read-panel layers, a per-CpG fraction-methylated summary track, and a per-read mean-methylation sidebar to `plot_methylation()`.

**Architecture:** All three features are argument-gated and default to preserving current output. Rasterization wraps only the heavy geom layers via `ggrastr` (a Suggests dependency) with a graceful fallback when absent. The fraction track and sidebar are additional optional panels/layers built from existing `methylation_data` fields.

**Tech Stack:** R, ggplot2, patchwork, testthat (edition 3), devtools, roxygen2, ggrastr (Suggests).

## Global Constraints

- R package; tests use `testthat` edition 3 via `devtools::test()`.
- All new behaviour argument-gated; existing default output unchanged.
- `ggrastr` is a **Suggests** dependency only — never `library()` it unconditionally; guard with `requireNamespace("ggrastr", quietly = TRUE)`.
- Internal helpers prefixed `.`, reached in tests via `ggmethylation:::`.
- Exported functions use roxygen2 (`markdown = TRUE`); run `devtools::document()` after doc edits.
- Assumes Milestone 1 is merged (relies on `theme_ggmethylation()`, `.PROB_GRADIENT`, and the extended panel-assembly block in `plot_methylation()`).

---

### Task 1: Optional rasterization of dense read-panel layers

**Files:**
- Modify: `DESCRIPTION` (add `ggrastr` to `Suggests`)
- Modify: `R/build_read_panel.R` (`build_read_panel()` signature; layer construction)
- Modify: `R/plot_methylation.R` (thread `rasterize`, `raster_dpi`)
- Create: `tests/testthat/test-rasterize.R`

**Interfaces:**
- Consumes: nothing new.
- Produces: `plot_methylation(rasterize = FALSE, raster_dpi = 300)`. Internal helper `.maybe_rasterize(layer, rasterize, dpi)` returns the layer wrapped in `ggrastr::rasterise()` when enabled and available, else the layer unchanged.

- [ ] **Step 1: Add `ggrastr` to Suggests**

Edit `DESCRIPTION`, appending to the `Suggests:` block:

```
    ggrastr,
```

- [ ] **Step 2: Write the failing test**

```r
# tests/testthat/test-rasterize.R
test_that("maybe_rasterize returns layer unchanged when disabled", {
  lyr <- ggplot2::geom_point()
  out <- ggmethylation:::.maybe_rasterize(lyr, rasterize = FALSE, dpi = 300)
  expect_identical(out, lyr)
})

test_that("maybe_rasterize falls back gracefully when ggrastr absent", {
  # When rasterize=TRUE but package missing, must warn once and return layer.
  lyr <- ggplot2::geom_point()
  if (!requireNamespace("ggrastr", quietly = TRUE)) {
    expect_warning(
      out <- ggmethylation:::.maybe_rasterize(lyr, rasterize = TRUE, dpi = 300),
      "ggrastr"
    )
    expect_identical(out, lyr)
  } else {
    out <- ggmethylation:::.maybe_rasterize(lyr, rasterize = TRUE, dpi = 300)
    expect_false(identical(out, lyr))
  }
})
```

- [ ] **Step 3: Run test to verify it fails**

Run: `Rscript -e 'devtools::test_file("tests/testthat/test-rasterize.R")'`
Expected: FAIL — `.maybe_rasterize` not found.

- [ ] **Step 4: Implement `.maybe_rasterize()`**

Add to `R/build_read_panel.R`:

```r
# Wrap a ggplot2 layer in ggrastr::rasterise when rasterization is requested
# and ggrastr is installed. Falls back to the raw layer (with a one-time
# warning) otherwise.
.maybe_rasterize <- function(layer, rasterize, dpi = 300) {
  if (!isTRUE(rasterize)) return(layer)
  if (!requireNamespace("ggrastr", quietly = TRUE)) {
    warning("`rasterize = TRUE` requires the 'ggrastr' package; drawing vector layers instead.",
            call. = FALSE)
    return(layer)
  }
  ggrastr::rasterise(layer, dpi = dpi)
}
```

- [ ] **Step 5: Apply to heavy layers**

Add `rasterize = FALSE`, `raster_dpi = 300` to the `build_read_panel()`
signature. Wrap the read-polygon `geom_polygon` layers and the mod-prob segment
layer. Because `.add_mod_prob_segments()` returns `p + geom + scale`, the
simplest approach is to build the heavy layers as objects and add the wrapped
form. For the read polygons, replace e.g.:

```r
      ggplot2::geom_polygon(
        data = read_polys,
        ggplot2::aes(x = .data$x, y = .data$y, group = .data$polygon_id,
                     fill = .data$group),
        colour = NA
      )
```

with:

```r
      .maybe_rasterize(
        ggplot2::geom_polygon(
          data = read_polys,
          ggplot2::aes(x = .data$x, y = .data$y, group = .data$polygon_id,
                       fill = .data$group),
          colour = NA
        ),
        rasterize, raster_dpi
      )
```

Apply the same wrapping to the strand and plain `geom_polygon` branches. For the
mod-prob segments, refactor `.add_mod_prob_segments()` to accept `rasterize`,
`raster_dpi` and wrap its `geom_segment`:

```r
.add_mod_prob_segments <- function(p, sites_plot, half_height, line_width,
                                    colour_low, colour_high,
                                    call_mode = "continuous",
                                    call_threshold = 0.5,
                                    call_ambiguous = NULL,
                                    rasterize = FALSE, raster_dpi = 300) {
  # ... build `seg` (the geom_segment) and `scl` (the scale) per call_mode ...
  p + .maybe_rasterize(seg, rasterize, raster_dpi) + scl
}
```

Refactor the binary/continuous branches so each computes a `seg` layer and a
`scl` scale, then share the final `p + .maybe_rasterize(seg, ...) + scl` return.

- [ ] **Step 6: Thread through `plot_methylation()`**

Add `rasterize = FALSE`, `raster_dpi = 300` to `plot_methylation()` and
`.plot_multi_methylation()`, and pass into every `build_read_panel()` call.

- [ ] **Step 7: Rendering smoke test**

```r
test_that("plot_methylation accepts rasterize argument", {
  reads <- data.frame(
    read_name = "r1", start = 1000L, end = 2000L, strand = "+",
    lane = 0L, mean_mod_prob = 0.5, clip_side = NA_character_,
    sa_chrom = NA_character_, stringsAsFactors = FALSE
  )
  sites <- data.frame(
    read_name = "r1", position = seq(1100, 1900, length.out = 5),
    mod_prob = c(0.1,0.3,0.5,0.7,0.9), mod_code = "m", stringsAsFactors = FALSE
  )
  md <- make_test_data(reads, sites)
  expect_no_error(
    suppressWarnings(
      ggmethylation::plot_methylation(md, rasterize = TRUE, show_supplementary = FALSE)
    )
  )
})
```

- [ ] **Step 8: Run tests, document, commit**

Run: `Rscript -e 'devtools::test_file("tests/testthat/test-rasterize.R"); devtools::test_file("tests/testthat/test-build_read_panel.R")'`
Expected: PASS.

Add `@param rasterize`, `@param raster_dpi` roxygen; NEWS bullet; `devtools::document()`.

```bash
git add DESCRIPTION R/build_read_panel.R R/plot_methylation.R tests/testthat/test-rasterize.R NEWS.md man/
git commit -m "feat: optional ggrastr rasterization of dense read layers"
```

---

### Task 2: Per-CpG fraction-methylated summary track

**Files:**
- Create: `R/fraction_track.R`
- Create: `tests/testthat/test-fraction_track.R`
- Modify: `R/plot_methylation.R` (`show_fraction` arg; panel assembly)

**Interfaces:**
- Consumes: `theme_ggmethylation()`, `.PROB_GRADIENT` (Milestone 1); `data$sites`; `call_threshold` (Milestone 1 arg).
- Produces: `.compute_site_fractions(sites, threshold, group_col = NULL)` → data.frame `position`, `frac`, `coverage`, and `group` when grouped. `.build_fraction_panel(frac_df, region_start, region_end, grouped, group_colours)` → a ggplot.

- [ ] **Step 1: Write the failing test**

```r
# tests/testthat/test-fraction_track.R
test_that("compute_site_fractions computes methylated fraction and coverage", {
  sites <- data.frame(
    position = c(100, 100, 100, 200, 200),
    mod_prob = c(0.9, 0.8, 0.1, 0.2, 0.3),
    stringsAsFactors = FALSE
  )
  res <- ggmethylation:::.compute_site_fractions(sites, threshold = 0.5)
  expect_equal(res$coverage[res$position == 100], 3L)
  # 2 of 3 calls >= 0.5 at position 100
  expect_equal(res$frac[res$position == 100], 2 / 3)
  expect_equal(res$frac[res$position == 200], 0)
})

test_that("compute_site_fractions splits by group when requested", {
  sites <- data.frame(
    position = c(100, 100, 100, 100),
    mod_prob = c(0.9, 0.1, 0.9, 0.9),
    group    = c("1", "1", "2", "2"),
    stringsAsFactors = FALSE
  )
  res <- ggmethylation:::.compute_site_fractions(sites, threshold = 0.5,
                                                 group_col = "group")
  expect_equal(res$frac[res$group == "1"], 0.5)
  expect_equal(res$frac[res$group == "2"], 1)
})
```

- [ ] **Step 2: Run test to verify it fails**

Run: `Rscript -e 'devtools::test_file("tests/testthat/test-fraction_track.R")'`
Expected: FAIL — `.compute_site_fractions` not found.

- [ ] **Step 3: Implement `.compute_site_fractions()`**

```r
# R/fraction_track.R

# Per-position methylated fraction (calls >= threshold) and coverage.
# When group_col is given, computes per group.
.compute_site_fractions <- function(sites, threshold = 0.5, group_col = NULL) {
  if (is.null(sites) || nrow(sites) == 0L) {
    base <- data.frame(position = numeric(0), frac = numeric(0),
                       coverage = integer(0), stringsAsFactors = FALSE)
    if (!is.null(group_col)) base[[group_col]] <- character(0)
    return(base)
  }
  sites$.meth <- as.integer(sites$mod_prob >= threshold)
  keys <- if (is.null(group_col)) "position" else c("position", group_col)
  agg_frac <- stats::aggregate(sites$.meth,
                               by = sites[keys], FUN = mean)
  agg_cov  <- stats::aggregate(sites$.meth,
                               by = sites[keys], FUN = length)
  out <- merge(agg_frac, agg_cov, by = keys)
  names(out)[(length(keys) + 1L):(length(keys) + 2L)] <- c("frac", "coverage")
  out$coverage <- as.integer(out$coverage)
  out[order(out$position), , drop = FALSE]
}
```

- [ ] **Step 4: Implement `.build_fraction_panel()`**

```r
# Point track of methylated fraction per site, alpha scaled by coverage.
.build_fraction_panel <- function(frac_df, region_start, region_end,
                                  grouped = FALSE, group_colours = NULL) {
  aes_base <- if (grouped) {
    ggplot2::aes(x = .data$position, y = .data$frac,
                 alpha = .data$coverage, colour = .data$group)
  } else {
    ggplot2::aes(x = .data$position, y = .data$frac, alpha = .data$coverage)
  }
  p <- ggplot2::ggplot(frac_df, aes_base) +
    ggplot2::geom_point(size = 1, colour = if (grouped) NULL else .PROB_GRADIENT$high) +
    ggplot2::scale_alpha_continuous(range = c(0.2, 1), name = "Coverage") +
    ggplot2::scale_y_continuous(limits = c(0, 1), name = "Fraction\nmethylated") +
    ggplot2::scale_x_continuous(labels = scales::comma_format()) +
    ggplot2::coord_cartesian(xlim = c(region_start, region_end)) +
    ggplot2::labs(x = "Genomic position (bp)") +
    theme_ggmethylation()
  if (grouped && !is.null(group_colours)) {
    p <- p + ggplot2::scale_colour_manual(values = group_colours,
                                          na.value = "grey50", name = "Group")
  }
  p
}
```

Note: `geom_point(colour = NULL)` is invalid; when `grouped` the colour comes
from the aesthetic, so build the geom conditionally:

```r
  geom <- if (grouped) ggplot2::geom_point(size = 1) else
          ggplot2::geom_point(size = 1, colour = .PROB_GRADIENT$high)
  p <- ggplot2::ggplot(frac_df, aes_base) + geom + ...
```

- [ ] **Step 5: Wire `show_fraction` into `plot_methylation()`**

Add `show_fraction = FALSE`. Build the panel after `p_bottom`:

```r
  p_fraction <- NULL
  if (isTRUE(show_fraction)) {
    grouped <- !is.null(data$group_tag)
    frac_df <- .compute_site_fractions(
      data$sites, threshold = call_threshold,
      group_col = if (grouped) "group" else NULL
    )
    if (nrow(frac_df) > 0L) {
      p_fraction <- .build_fraction_panel(frac_df, region_start, region_end,
                                          grouped = grouped,
                                          group_colours = group_colours)
    }
  }
```

Append `p_fraction` to the `panels` list (after `p_bottom`, before/after
`p_delta` — place it directly after `p_bottom`) and add a default height of
`0.25`. Extend the `panel_heights` length checks accordingly. When both
`p_delta` and `p_fraction` are present the bottom-most panel keeps the x-axis
title; hide it on the others (reuse the `axis.title.x = element_blank()` pattern
from Milestone 1's delta wiring).

- [ ] **Step 6: Smoke test + run**

```r
test_that("plot_methylation adds a fraction panel when requested", {
  reads <- data.frame(
    read_name = c("r1","r2"), start = 1000L, end = 2000L, strand = "+",
    lane = 0L, mean_mod_prob = 0.5, clip_side = NA_character_,
    sa_chrom = NA_character_, stringsAsFactors = FALSE
  )
  sites <- data.frame(
    read_name = rep(c("r1","r2"), each = 4),
    position = rep(seq(1100, 1900, length.out = 4), 2),
    mod_prob = c(0.9,0.1,0.8,0.2, 0.7,0.3,0.6,0.4), mod_code = "m",
    stringsAsFactors = FALSE
  )
  md <- make_test_data(reads, sites)
  p <- ggmethylation::plot_methylation(md, show_fraction = TRUE, show_supplementary = FALSE)
  expect_s3_class(p, "patchwork")
})
```

Run: `Rscript -e 'devtools::test_file("tests/testthat/test-fraction_track.R")'`
Expected: PASS.

- [ ] **Step 7: Document + commit**

Add `@param show_fraction` roxygen (note it uses `call_threshold`); NEWS bullet;
`devtools::document()`.

```bash
git add R/fraction_track.R R/plot_methylation.R tests/testthat/test-fraction_track.R NEWS.md man/
git commit -m "feat: add per-CpG fraction-methylated summary track (show_fraction)"
```

---

### Task 3: Per-read mean-methylation sidebar

**Files:**
- Modify: `R/build_read_panel.R` (`show_read_summary` arg; sidebar layer)
- Modify: `R/plot_methylation.R` (thread `show_read_summary`)
- Create: `tests/testthat/test-read_sidebar.R`

**Interfaces:**
- Consumes: `data$reads` (`lane`, `mean_mod_prob`), `.PROB_GRADIENT`, arrow geometry (`region_start`, `region_end`).
- Produces: `.make_read_summary_tiles(reads, region_start, region_end, half_height, width_frac = 0.02)` → data.frame `xmin`, `xmax`, `ymin`, `ymax`, `mean_mod_prob`, `lane`. `build_read_panel(show_read_summary = FALSE)` draws a marginal tile column just left of `region_start`.

- [ ] **Step 1: Write the failing test**

```r
# tests/testthat/test-read_sidebar.R
test_that("make_read_summary_tiles places tiles left of region start", {
  reads <- data.frame(
    read_name = c("r1", "r2"), lane = c(0L, 1L),
    mean_mod_prob = c(0.2, 0.8), stringsAsFactors = FALSE
  )
  tiles <- ggmethylation:::.make_read_summary_tiles(
    reads, region_start = 1000, region_end = 2000, half_height = 0.35
  )
  expect_equal(nrow(tiles), 2L)
  # All tiles sit at or left of region_start
  expect_true(all(tiles$xmax <= 1000))
  # ymin/ymax straddle the lane by half_height
  expect_equal(tiles$ymin[tiles$lane == 0L], -0.35)
  expect_equal(tiles$ymax[tiles$lane == 0L],  0.35)
})
```

- [ ] **Step 2: Run test to verify it fails**

Run: `Rscript -e 'devtools::test_file("tests/testthat/test-read_sidebar.R")'`
Expected: FAIL — `.make_read_summary_tiles` not found.

- [ ] **Step 3: Implement the tile builder**

```r
# Build one marginal tile per read, coloured by mean_mod_prob, positioned just
# left of the region start. width_frac is the tile width as a fraction of the
# region span.
.make_read_summary_tiles <- function(reads, region_start, region_end,
                                     half_height, width_frac = 0.02) {
  if (nrow(reads) == 0L) {
    return(data.frame(xmin = numeric(0), xmax = numeric(0),
                      ymin = numeric(0), ymax = numeric(0),
                      mean_mod_prob = numeric(0), lane = numeric(0),
                      stringsAsFactors = FALSE))
  }
  span  <- region_end - region_start
  w     <- span * width_frac
  gap   <- span * 0.005
  xmax  <- region_start - gap
  xmin  <- xmax - w
  data.frame(
    xmin = xmin, xmax = xmax,
    ymin = reads$lane - half_height,
    ymax = reads$lane + half_height,
    mean_mod_prob = reads$mean_mod_prob,
    lane = reads$lane,
    stringsAsFactors = FALSE
  )
}
```

- [ ] **Step 4: Draw the sidebar in `build_read_panel()`**

Add `show_read_summary = FALSE` to the signature. After the read polygons are
added and before the theme block, add:

```r
  if (isTRUE(show_read_summary)) {
    tiles <- .make_read_summary_tiles(data$reads, region_start, region_end,
                                      half_height)
    if (nrow(tiles) > 0L) {
      p <- p +
        ggnewscale::new_scale_fill() +
        ggplot2::geom_rect(
          data = tiles,
          ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax,
                       ymin = .data$ymin, ymax = .data$ymax,
                       fill = .data$mean_mod_prob),
          inherit.aes = FALSE
        ) +
        ggplot2::scale_fill_gradient(
          low = .PROB_GRADIENT$low, high = .PROB_GRADIENT$high,
          limits = c(0, 1), name = "Read mean\nmod. prob."
        )
    }
  }
```

The `coord_cartesian(xlim = c(region_start, region_end))` clips the tiles, so
widen the left limit when the sidebar is on:

```r
  x_left <- if (isTRUE(show_read_summary)) {
    region_start - (region_end - region_start) * 0.03
  } else region_start
  # ... use xlim = c(x_left, region_end) in coord_cartesian ...
```

Replace the existing `coord_cartesian(xlim = c(region_start, region_end))` in
the theme block with `coord_cartesian(xlim = c(x_left, region_end))`.

- [ ] **Step 5: Thread through `plot_methylation()`**

Add `show_read_summary = FALSE` to `plot_methylation()` and
`.plot_multi_methylation()`; pass into every `build_read_panel()` call. Note:
the smooth/delta/fraction panels keep `xlim = c(region_start, region_end)`, so a
small x-axis misalignment with the read panel is expected when the sidebar is
on — document this as a known limitation.

- [ ] **Step 6: Smoke test + run**

```r
test_that("plot_methylation renders with read summary sidebar", {
  reads <- data.frame(
    read_name = "r1", start = 1000L, end = 2000L, strand = "+",
    lane = 0L, mean_mod_prob = 0.5, clip_side = NA_character_,
    sa_chrom = NA_character_, stringsAsFactors = FALSE
  )
  sites <- data.frame(
    read_name = "r1", position = seq(1100, 1900, length.out = 5),
    mod_prob = c(0.1,0.3,0.5,0.7,0.9), mod_code = "m", stringsAsFactors = FALSE
  )
  md <- make_test_data(reads, sites)
  expect_no_error(
    ggmethylation::plot_methylation(md, show_read_summary = TRUE, show_supplementary = FALSE)
  )
})
```

Run: `Rscript -e 'devtools::test_file("tests/testthat/test-read_sidebar.R"); devtools::test()'`
Expected: PASS.

- [ ] **Step 7: Document + commit**

Add `@param show_read_summary` roxygen (note the known x-axis-alignment
limitation across panels); NEWS bullet; `devtools::document()`.

```bash
git add R/build_read_panel.R R/plot_methylation.R tests/testthat/test-read_sidebar.R NEWS.md man/
git commit -m "feat: add per-read mean-methylation sidebar (show_read_summary)"
```

---

## Self-Review Notes

- **Spec coverage:** C1 → Task 1; A3 → Task 2; B2 → Task 3. All Milestone-2 features covered.
- **Type consistency:** `.compute_site_fractions()` produces `position/frac/coverage[/group]`, consumed by `.build_fraction_panel()`. `.make_read_summary_tiles()` produces `xmin/xmax/ymin/ymax/mean_mod_prob/lane`, consumed by the `geom_rect` in Task 3.
- **Dependency:** `ggrastr` guarded via `requireNamespace`; tests branch on its availability so the suite passes with or without it installed.
- **Known limitations flagged inline:** sidebar x-axis alignment vs lower panels (Task 3 Step 5); `geom_point(colour=NULL)` corrected to conditional geom (Task 2 Step 4).
- **Ordering:** Tasks are independent of each other but all assume Milestone 1's `theme_ggmethylation()`, `.PROB_GRADIENT`, `call_threshold` arg, and extended panel-assembly block are present.
```