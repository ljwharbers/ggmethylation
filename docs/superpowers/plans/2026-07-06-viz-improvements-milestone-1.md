# Visualization Improvements — Milestone 1 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a group-difference delta track, confidence ribbons on the smooth panel, a binarized call mode, and centralized theming/colorblind-safe palettes to `plot_methylation()`.

**Architecture:** Each feature is gated behind a new argument that defaults to preserving current output (except `show_ci`, which defaults on). New arguments thread from `plot_methylation()` down through `build_read_panel()` / the smooth-panel builders. A new `R/palettes.R` centralizes colour constants and the exported `theme_ggmethylation()`. The delta track computes both group loess curves on a shared position grid and renders a diverging area panel appended to the patchwork composite.

**Tech Stack:** R, ggplot2, patchwork, stats::loess, testthat (edition 3), devtools, roxygen2.

## Global Constraints

- R package; tests use `testthat` edition 3 via `devtools::test()`.
- All new behaviour must be argument-gated; existing default output must not change except deliberate colorblind-safety upgrades documented in `NEWS.md`.
- Internal helpers are prefixed `.` and reached in tests via `ggmethylation:::`.
- Exported functions use roxygen2 with `markdown = TRUE`; run `devtools::document()` after editing docs.
- Integration tests must skip gracefully when the BAM fixture is absent (follow `test-integration.R`).
- No new hard dependencies in `Imports`; Milestone 1 uses only existing deps (`ggplot2`, `stats`, `patchwork`, `scales`).

---

### Task 1: Centralized palettes and `theme_ggmethylation()`

Foundational. Extracts the theme block duplicated in `build_read_panel()` and `.smooth_panel_base()` into one exported function, and gathers colour constants (including the diverging palette used by Task 6) into one file.

**Files:**
- Create: `R/palettes.R`
- Create: `tests/testthat/test-palettes.R`
- Modify: `R/build_read_panel.R` (theme block near line 540–550)
- Modify: `R/plot_methylation.R` (`.smooth_panel_base()` near line 196–209)
- Modify: `NEWS.md` (create if absent)

**Interfaces:**
- Produces: `theme_ggmethylation()` → a ggplot2 theme object (exported). Internal constants `.OKABE_ITO`, `.GROUP_PALETTE_DEFAULT` (named char vector), `.DELTA_DIVERGING` (list with `neg`, `pos`, `zero` colours), `.PROB_GRADIENT` (list with `low`, `high`).

- [ ] **Step 1: Write the failing test**

```r
# tests/testthat/test-palettes.R
test_that("theme_ggmethylation returns a ggplot2 theme", {
  th <- ggmethylation::theme_ggmethylation()
  expect_s3_class(th, "theme")
})

test_that("palette constants have expected shape", {
  expect_true(is.character(ggmethylation:::.OKABE_ITO))
  expect_gte(length(ggmethylation:::.OKABE_ITO), 8L)
  gp <- ggmethylation:::.GROUP_PALETTE_DEFAULT
  expect_true(all(c("1", "2") %in% names(gp)))
  dv <- ggmethylation:::.DELTA_DIVERGING
  expect_true(all(c("neg", "pos", "zero") %in% names(dv)))
})
```

- [ ] **Step 2: Run test to verify it fails**

Run: `Rscript -e 'devtools::test_file("tests/testthat/test-palettes.R")'`
Expected: FAIL — object `theme_ggmethylation` / `.OKABE_ITO` not found.

- [ ] **Step 3: Write minimal implementation**

```r
# R/palettes.R

# Okabe-Ito colorblind-safe qualitative palette
.OKABE_ITO <- c(
  "#E69F00", "#56B4E9", "#009E73", "#F0E442",
  "#0072B2", "#D55E00", "#CC79A7", "#000000"
)

# Default group palette. Named "1"/"2" to match HP haplotype tag output.
# Uses two well-separated Okabe-Ito hues (blue / orange).
.GROUP_PALETTE_DEFAULT <- c("1" = "#0072B2", "2" = "#E69F00")

# Probability gradient endpoints (colorblind-safe grey -> red)
.PROB_GRADIENT <- list(low = "#BDBDBD", high = "#C62828")

# Diverging palette for the delta track
.DELTA_DIVERGING <- list(neg = "#0072B2", pos = "#D55E00", zero = "grey60")

#' ggmethylation plot theme
#'
#' Shared minimal theme used across `ggmethylation` panels. Provides the base
#' `theme_minimal()` plus compact legend styling and suppressed minor grid.
#'
#' @return A [ggplot2::theme] object.
#' @export
theme_ggmethylation <- function() {
  ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      legend.text      = ggplot2::element_text(size = ggplot2::rel(0.75)),
      legend.title     = ggplot2::element_text(size = ggplot2::rel(0.75)),
      legend.key.size  = ggplot2::unit(0.4, "cm")
    )
}
```

- [ ] **Step 4: Run test to verify it passes**

Run: `Rscript -e 'devtools::document(); devtools::test_file("tests/testthat/test-palettes.R")'`
Expected: PASS.

- [ ] **Step 5: Refactor call sites to use the shared theme**

In `R/plot_methylation.R`, replace the `theme_minimal()` + `theme(...)` portion of `.smooth_panel_base()` with `theme_ggmethylation()`:

```r
.smooth_panel_base <- function(region_start, region_end) {
  list(
    ggplot2::scale_y_continuous(limits = c(0, 1), name = "Mean modification\nprobability"),
    ggplot2::scale_x_continuous(labels = scales::comma_format()),
    ggplot2::coord_cartesian(xlim = c(region_start, region_end)),
    theme_ggmethylation()
  )
}
```

In `R/build_read_panel.R`, replace the `ggplot2::theme_minimal() + ggplot2::theme(...)` block (near line 540) so the shared theme is applied first, then the read-panel-specific overrides layer on top:

```r
  p <- p +
    ggplot2::scale_y_reverse() +
    ggplot2::coord_cartesian(xlim = c(region_start, region_end)) +
    theme_ggmethylation() +
    ggplot2::theme(
      axis.text.y        = ggplot2::element_blank(),
      axis.ticks.y       = ggplot2::element_blank(),
      axis.title.y       = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank()
    ) +
    ggplot2::labs(x = NULL)
```

- [ ] **Step 6: Update `group_colours` default to the centralized palette**

In `R/plot_methylation.R`, the `plot_methylation()` signature currently has
`group_colours = c("1" = "#95babc", "2" = "#efbb76")`. Change the default to
reference the constant so there is one source of truth:

```r
                             group_colours = .GROUP_PALETTE_DEFAULT,
```

Do the same for the `.plot_multi_methylation()` internal call path if it hard-codes a default (it inherits from the caller, so no change needed there).

- [ ] **Step 7: Record the palette change in NEWS**

Create/append `NEWS.md`:

```markdown
# ggmethylation (development version)

## New features
* `theme_ggmethylation()` exposes the shared panel theme.

## Changes
* Default `group_colours` updated to a colorblind-safe Okabe-Ito blue/orange
  pair. Pass an explicit `group_colours` vector to restore prior colours.
```

- [ ] **Step 8: Run the full suite to confirm no regressions**

Run: `Rscript -e 'devtools::test()'`
Expected: PASS. Update any existing snapshot that legitimately changed only due to the deliberate palette default (review the diff before accepting).

- [ ] **Step 9: Commit**

```bash
git add R/palettes.R tests/testthat/test-palettes.R R/build_read_panel.R R/plot_methylation.R NEWS.md NAMESPACE man/
git commit -m "feat: add theme_ggmethylation() and centralize colorblind-safe palettes"
```

---

### Task 2: Confidence-interval columns in `smooth_methylation()`

**Files:**
- Modify: `R/smooth_methylation.R`
- Modify: `tests/testthat/test-smooth_methylation.R`

**Interfaces:**
- Consumes: nothing new.
- Produces: `smooth_methylation()` return gains numeric columns `lower` and `upper` (present in every return path). On the loess path they equal `fit ± 1.96 * se.fit` clamped to `[0, 1]`; on the raw-means fallback and empty paths they are `NA_real_`.

- [ ] **Step 1: Write the failing test**

```r
test_that("smooth_methylation returns clamped lower/upper on loess path", {
  sites <- data.frame(
    position = 1:5 * 100,
    mod_prob = c(0.1, 0.2, 0.3, 0.4, 0.5),
    group = "A",
    stringsAsFactors = FALSE
  )
  result <- ggmethylation:::smooth_methylation(sites)
  expect_true(all(c("lower", "upper") %in% names(result)))
  ok <- !is.na(result$lower)
  expect_true(all(result$lower[ok] >= 0 & result$lower[ok] <= 1))
  expect_true(all(result$upper[ok] >= 0 & result$upper[ok] <= 1))
  expect_true(all(result$upper[ok] >= result$lower[ok]))
})

test_that("smooth_methylation sets NA CI on raw-means fallback", {
  sites <- data.frame(
    position = c(1, 1, 2, 3),
    mod_prob = c(0.8, 0.6, 0.5, 0.2),
    group = "A",
    stringsAsFactors = FALSE
  )
  result <- ggmethylation:::smooth_methylation(sites)
  expect_true(all(is.na(result$lower)))
  expect_true(all(is.na(result$upper)))
})

test_that("smooth_methylation empty input includes lower/upper cols", {
  result <- ggmethylation:::smooth_methylation(NULL)
  expect_true(all(c("lower", "upper") %in% names(result)))
})
```

- [ ] **Step 2: Run test to verify it fails**

Run: `Rscript -e 'devtools::test_file("tests/testthat/test-smooth_methylation.R")'`
Expected: FAIL — `lower`/`upper` not in names.

- [ ] **Step 3: Write the implementation**

Edit `R/smooth_methylation.R`. Update `out_cols`, the empty-return frame, both
per-group branches, and the composite-key split so `lower`/`upper` survive.

```r
  out_cols <- c("position", "mean_prob", "lower", "upper", effective_group_col)

  if (is.null(sites) || nrow(sites) == 0L) {
    out <- data.frame(
      position  = numeric(0L),
      mean_prob = numeric(0L),
      lower     = numeric(0L),
      upper     = numeric(0L),
      group     = character(0L),
      stringsAsFactors = FALSE
    )
    names(out)[5L] <- effective_group_col
    return(out)
  }
```

In the per-group loop, replace the raw-means and loess branches:

```r
    if (nrow(agg) < 4L) {
      # Too few unique positions for loess; return raw means, no CI
      df <- agg
      df$lower <- NA_real_
      df$upper <- NA_real_
    } else {
      df <- tryCatch(
        suppressWarnings({
          fit  <- stats::loess(mean_prob ~ position, data = agg, span = effective_span)
          grid <- seq(min(agg$position, na.rm = TRUE),
                      max(agg$position, na.rm = TRUE),
                      length.out = 200L)
          pr   <- stats::predict(fit, newdata = data.frame(position = grid), se = TRUE)
          lower <- pmin(pmax(pr$fit - 1.96 * pr$se.fit, 0), 1)
          upper <- pmin(pmax(pr$fit + 1.96 * pr$se.fit, 0), 1)
          data.frame(position = grid, mean_prob = pr$fit,
                     lower = lower, upper = upper)
        }),
        error = function(e) {
          agg$lower <- NA_real_
          agg$upper <- NA_real_
          agg
        }
      )
    }
```

The final `out <- out[, out_cols, drop = FALSE]` now carries `lower`/`upper`.
The composite-key split block (mod_code) is unaffected because it only rewrites
the group column.

- [ ] **Step 4: Run test to verify it passes**

Run: `Rscript -e 'devtools::test_file("tests/testthat/test-smooth_methylation.R")'`
Expected: PASS. Also run `test-smooth_deletion_breaks.R` to confirm the deletion-break helper still works with the wider frame.

Run: `Rscript -e 'devtools::test_file("tests/testthat/test-smooth_deletion_breaks.R")'`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add R/smooth_methylation.R tests/testthat/test-smooth_methylation.R
git commit -m "feat: add loess confidence-interval columns to smooth_methylation"
```

---

### Task 3: Render confidence ribbon on the smooth panel

**Files:**
- Modify: `R/plot_methylation.R` (smooth-panel branches ~ lines 456–553; `.insert_deletion_breaks` ~ line 109)
- Modify: `tests/testthat/test-smooth_deletion_breaks.R` or add `tests/testthat/test-smooth_ci_ribbon.R`

**Interfaces:**
- Consumes: `smooth_methylation()` `lower`/`upper` columns (Task 2).
- Produces: `plot_methylation(show_ci = TRUE)` adds a `geom_ribbon` behind each smooth line. `.insert_deletion_breaks()` nulls `lower`/`upper` alongside `mean_prob`.

- [ ] **Step 1: Write the failing test**

```r
# tests/testthat/test-smooth_ci_ribbon.R
test_that("insert_deletion_breaks nulls CI columns inside deletions", {
  smoothed <- data.frame(
    position  = c(100, 150, 200, 250),
    mean_prob = c(0.5, 0.5, 0.5, 0.5),
    lower     = c(0.4, 0.4, 0.4, 0.4),
    upper     = c(0.6, 0.6, 0.6, 0.6),
    group     = "A",
    stringsAsFactors = FALSE
  )
  ranges <- data.frame(del_start = 140, del_end = 210,
                       group = "A", stringsAsFactors = FALSE)
  out <- ggmethylation:::.insert_deletion_breaks(smoothed, ranges, "group")
  masked <- out$position >= 140 & out$position <= 210 & !is.na(out$position)
  expect_true(all(is.na(out$mean_prob[masked])))
  expect_true(all(is.na(out$lower[masked])))
  expect_true(all(is.na(out$upper[masked])))
})
```

- [ ] **Step 2: Run test to verify it fails**

Run: `Rscript -e 'devtools::test_file("tests/testthat/test-smooth_ci_ribbon.R")'`
Expected: FAIL — `lower`/`upper` not nulled (only `mean_prob` currently masked).

- [ ] **Step 3: Update `.insert_deletion_breaks()` to null CI columns**

In `R/plot_methylation.R`, inside `.insert_deletion_breaks()`, after the line
`smoothed$mean_prob[in_grp & in_del] <- NA_real_`, add:

```r
    if ("lower" %in% names(smoothed)) smoothed$lower[in_grp & in_del] <- NA_real_
    if ("upper" %in% names(smoothed)) smoothed$upper[in_grp & in_del] <- NA_real_
```

Note: `id_cols <- setdiff(names(smoothed), c("position", "mean_prob"))` currently
treats `lower`/`upper` as identity columns. Change it to also exclude them:

```r
  id_cols <- setdiff(names(smoothed), c("position", "mean_prob", "lower", "upper"))
```

And when building sentinel rows, set `s1$lower <- NA_real_; s1$upper <- NA_real_`
(and likewise `s2`) so the rebind aligns columns:

```r
      s1$position  <- del_start - 0.5
      s1$mean_prob <- NA_real_
      if ("lower" %in% names(s1)) { s1$lower <- NA_real_; s1$upper <- NA_real_ }
      s2$position  <- del_end + 0.5
      s2$mean_prob <- NA_real_
      if ("lower" %in% names(s2)) { s2$lower <- NA_real_; s2$upper <- NA_real_ }
```

- [ ] **Step 4: Add `show_ci` arg and ribbon layers**

Add `show_ci = TRUE` to the `plot_methylation()` signature (and thread it into
`.plot_multi_methylation()` with the same default). Define a small helper near
the top of `R/plot_methylation.R`:

```r
# Add a CI ribbon behind the smooth line(s) when requested and columns exist.
.add_ci_ribbon <- function(p, smoothed, show_ci, fill_aes = NULL) {
  if (!isTRUE(show_ci)) return(p)
  if (!all(c("lower", "upper") %in% names(smoothed))) return(p)
  rib <- smoothed[!is.na(smoothed$lower) & !is.na(smoothed$upper), , drop = FALSE]
  if (nrow(rib) == 0L) return(p)
  aes_args <- list(x = quote(.data$position),
                   ymin = quote(.data$lower),
                   ymax = quote(.data$upper))
  if (!is.null(fill_aes)) aes_args$fill <- fill_aes
  ribbon <- ggplot2::geom_ribbon(
    data = rib,
    mapping = do.call(ggplot2::aes, aes_args),
    alpha = 0.2, colour = NA,
    inherit.aes = FALSE,
    show.legend = FALSE
  )
  # Insert ribbon *before* existing line layers so it renders behind them.
  p$layers <- c(list(ribbon), p$layers)
  p
}
```

In each smooth-panel branch, after constructing `p_bottom` (before returning),
call the helper. For the ungrouped single-code branch:

```r
      p_bottom <- .add_ci_ribbon(p_bottom, smoothed, show_ci)
```

For grouped-by-group and multi-code branches, pass the matching fill aesthetic
so ribbon colour tracks the line:

```r
      p_bottom <- .add_ci_ribbon(p_bottom, smoothed, show_ci,
                                 fill_aes = quote(.data$group))
```

For a grouped ribbon, also add `ggplot2::scale_fill_manual(values = group_colours,
na.value = "grey50", guide = "none")` when `group_colours` is non-NULL, using
`ggnewscale::new_scale_fill()` is **not** needed because the line uses `colour`,
not `fill`.

- [ ] **Step 5: Add a rendering smoke test**

```r
test_that("plot_methylation smooth panel gains a ribbon layer when grouped", {
  # Build minimal grouped methylation_data (reuse make_test_data pattern)
  reads <- data.frame(
    read_name = c("r1", "r2"), start = c(1000L, 1000L), end = c(2000L, 2000L),
    strand = c("+", "+"), lane = 0L, mean_mod_prob = 0.5,
    group = c("1", "2"), clip_side = NA_character_, sa_chrom = NA_character_,
    stringsAsFactors = FALSE
  )
  sites <- data.frame(
    read_name = rep(c("r1", "r2"), each = 5),
    position = rep(seq(1100, 1900, length.out = 5), 2),
    mod_prob = c(0.1,0.3,0.5,0.7,0.9, 0.2,0.4,0.6,0.8,0.95),
    mod_code = "m", group = rep(c("1","2"), each = 5),
    stringsAsFactors = FALSE
  )
  md <- make_test_data(reads, sites)
  md$group_tag <- "HP"
  p <- ggmethylation::plot_methylation(md, show_ci = TRUE, show_supplementary = FALSE)
  expect_s3_class(p, "patchwork")
})
```

Add `make_test_data` to this test file or source it; if the helper lives only in
`test-build_read_panel.R`, copy the minimal constructor into the new test file.

- [ ] **Step 6: Run tests**

Run: `Rscript -e 'devtools::test_file("tests/testthat/test-smooth_ci_ribbon.R")'`
Expected: PASS.

- [ ] **Step 7: Update roxygen and NEWS**

Add `@param show_ci` documentation to `plot_methylation()` and a NEWS bullet.
Run `Rscript -e 'devtools::document()'`.

- [ ] **Step 8: Commit**

```bash
git add R/plot_methylation.R tests/testthat/test-smooth_ci_ribbon.R NEWS.md man/
git commit -m "feat: draw loess confidence ribbon on smooth panel (show_ci)"
```

---

### Task 4: Binarized call mode in the read panel

**Files:**
- Modify: `R/plot_methylation.R` (signature + pass-through to `build_read_panel`)
- Modify: `R/build_read_panel.R` (`build_read_panel()` signature; `.add_mod_prob_segments()`)
- Create: `tests/testthat/test-binary_call_mode.R`

**Interfaces:**
- Consumes: nothing new.
- Produces: `plot_methylation(call_mode = "binary", call_threshold = 0.5, call_ambiguous = NULL)`. `.add_mod_prob_segments()` gains `call_mode`, `call_threshold`, `call_ambiguous` params and branches its colour aesthetic/scale.

- [ ] **Step 1: Write the failing test**

```r
# tests/testthat/test-binary_call_mode.R
test_that("binary call classification splits at threshold", {
  probs <- c(0.1, 0.49, 0.5, 0.9)
  cls <- ggmethylation:::.classify_calls(probs, threshold = 0.5, ambiguous = NULL)
  expect_equal(cls, c("unmethylated", "unmethylated", "methylated", "methylated"))
})

test_that("ambiguous band labels near-threshold calls", {
  probs <- c(0.1, 0.45, 0.5, 0.55, 0.9)
  cls <- ggmethylation:::.classify_calls(probs, threshold = 0.5, ambiguous = 0.1)
  expect_equal(cls, c("unmethylated", "ambiguous", "methylated", "ambiguous", "methylated"))
})
```

- [ ] **Step 2: Run test to verify it fails**

Run: `Rscript -e 'devtools::test_file("tests/testthat/test-binary_call_mode.R")'`
Expected: FAIL — `.classify_calls` not found.

- [ ] **Step 3: Implement the classifier in `R/build_read_panel.R`**

```r
# Classify continuous modification probabilities into discrete calls.
# ambiguous: NULL for a hard threshold, or a numeric half-width; probs within
# [threshold - w, threshold + w) that are not exactly >= threshold-only cases
# are labelled "ambiguous". At/above threshold and outside the band -> methylated.
.classify_calls <- function(probs, threshold = 0.5, ambiguous = NULL) {
  cls <- ifelse(probs >= threshold, "methylated", "unmethylated")
  if (!is.null(ambiguous) && ambiguous > 0) {
    band <- abs(probs - threshold) < ambiguous
    cls[band] <- "ambiguous"
  }
  cls
}
```

- [ ] **Step 4: Branch `.add_mod_prob_segments()` on call mode**

Replace `.add_mod_prob_segments()` with a version that accepts the new params:

```r
.add_mod_prob_segments <- function(p, sites_plot, half_height, line_width,
                                    colour_low, colour_high,
                                    call_mode = "continuous",
                                    call_threshold = 0.5,
                                    call_ambiguous = NULL) {
  if (identical(call_mode, "binary") && nrow(sites_plot) > 0L) {
    sites_plot$.call <- .classify_calls(sites_plot$mod_prob,
                                        call_threshold, call_ambiguous)
    vals <- c(unmethylated = colour_low, methylated = colour_high,
              ambiguous = "grey75")
    return(
      p +
        ggplot2::geom_segment(
          data = sites_plot,
          ggplot2::aes(
            x = .data$position, xend = .data$position,
            y = .data$lane - half_height, yend = .data$lane + half_height,
            colour = .data$.call
          ),
          linewidth = line_width
        ) +
        ggplot2::scale_colour_manual(
          values = vals, name = "Call",
          breaks = c("unmethylated", "methylated", "ambiguous")
        )
    )
  }
  p +
    ggplot2::geom_segment(
      data = sites_plot,
      ggplot2::aes(
        x = .data$position, xend = .data$position,
        y = .data$lane - half_height, yend = .data$lane + half_height,
        colour = .data$mod_prob
      ),
      linewidth = line_width
    ) +
    ggplot2::scale_colour_gradient(
      low = colour_low, high = colour_high,
      limits = c(0, 1),
      name = "Modification\nprobability"
    )
}
```

- [ ] **Step 5: Thread params through `build_read_panel()`**

Add `call_mode = "continuous"`, `call_threshold = 0.5`, `call_ambiguous = NULL`
to the `build_read_panel()` signature, and pass them into every
`.add_mod_prob_segments(...)` call site (three branches: grouped, strand,
plain):

```r
    p <- .add_mod_prob_segments(p, sites_plot, half_height, line_width,
                                colour_low, colour_high,
                                call_mode, call_threshold, call_ambiguous)
```

- [ ] **Step 6: Thread params through `plot_methylation()`**

Add `call_mode = c("continuous", "binary")`, `call_threshold = 0.5`,
`call_ambiguous = NULL` to `plot_methylation()`. At the top of the body:

```r
  call_mode <- match.arg(call_mode)
```

Pass all three into both `build_read_panel()` calls in `plot_methylation()` and
into the `.plot_multi_methylation()` signature + its `build_read_panel()` call.

- [ ] **Step 7: Rendering smoke test**

```r
test_that("plot_methylation renders in binary mode", {
  reads <- data.frame(
    read_name = "r1", start = 1000L, end = 2000L, strand = "+",
    lane = 0L, mean_mod_prob = 0.5, clip_side = NA_character_,
    sa_chrom = NA_character_, stringsAsFactors = FALSE
  )
  sites <- data.frame(
    read_name = "r1", position = seq(1100, 1900, length.out = 5),
    mod_prob = c(0.1, 0.3, 0.5, 0.7, 0.9), mod_code = "m",
    stringsAsFactors = FALSE
  )
  md <- make_test_data(reads, sites)
  p <- ggmethylation::plot_methylation(md, call_mode = "binary",
                                       show_supplementary = FALSE)
  expect_true(inherits(p, "ggplot") || inherits(p, "patchwork"))
})
```

- [ ] **Step 8: Run tests**

Run: `Rscript -e 'devtools::test_file("tests/testthat/test-binary_call_mode.R")'`
Expected: PASS. Then run `test-build_read_panel.R` to confirm continuous mode unchanged.

- [ ] **Step 9: Document + commit**

Add `@param call_mode`, `@param call_threshold`, `@param call_ambiguous` to
`plot_methylation()` roxygen; add NEWS bullet; `devtools::document()`.

```bash
git add R/build_read_panel.R R/plot_methylation.R tests/testthat/test-binary_call_mode.R NEWS.md man/
git commit -m "feat: add binarized call mode to read panel (call_mode)"
```

---

### Task 5: Shared-grid delta computation

**Files:**
- Create: `R/delta_track.R`
- Create: `tests/testthat/test-delta_track.R`

**Interfaces:**
- Consumes: `smooth_methylation()` (fits per group).
- Produces: `.compute_group_delta(sites, group_col, span, n_grid = 200)` → either `NULL` (with a `message()`) when the number of non-NA groups ≠ 2, or a data.frame with columns `position`, `delta`, `sign` (`"pos"`/`"neg"`/`"zero"`), where `delta = value(group2) - value(group1)` on a shared grid, `NA` where either group lacks support.

- [ ] **Step 1: Write the failing test**

```r
# tests/testthat/test-delta_track.R
test_that("compute_group_delta returns NULL and messages for non-2 groups", {
  sites <- data.frame(
    position = 1:5 * 100, mod_prob = seq(0.1, 0.5, length.out = 5),
    group = "A", stringsAsFactors = FALSE
  )
  expect_message(
    res <- ggmethylation:::.compute_group_delta(sites, "group", span = NULL),
    "exactly two groups"
  )
  expect_null(res)
})

test_that("compute_group_delta returns signed delta on shared grid for 2 groups", {
  sites <- data.frame(
    position = rep(1:6 * 100, 2),
    mod_prob = c(rep(0.2, 6), rep(0.8, 6)),
    group    = rep(c("1", "2"), each = 6),
    stringsAsFactors = FALSE
  )
  res <- ggmethylation:::.compute_group_delta(sites, "group", span = NULL)
  expect_true(all(c("position", "delta", "sign") %in% names(res)))
  # group2 (0.8) - group1 (0.2) ~ +0.6 where both supported
  ok <- !is.na(res$delta)
  expect_true(mean(res$delta[ok]) > 0)
  expect_true(all(res$sign[ok & res$delta > 0] == "pos"))
})
```

- [ ] **Step 2: Run test to verify it fails**

Run: `Rscript -e 'devtools::test_file("tests/testthat/test-delta_track.R")'`
Expected: FAIL — `.compute_group_delta` not found.

- [ ] **Step 3: Implement `.compute_group_delta()`**

```r
# R/delta_track.R

# Compute a signed methylation difference between exactly two groups over a
# shared position grid. Returns NULL (with a message) when != 2 groups.
.compute_group_delta <- function(sites, group_col, span = NULL, n_grid = 200L) {
  groups <- unique(sites[[group_col]])
  groups <- sort(groups[!is.na(groups)])
  if (length(groups) != 2L) {
    message("Delta track requires exactly two groups; skipping.")
    return(NULL)
  }

  fit_one <- function(grp, grid) {
    sub <- sites[!is.na(sites[[group_col]]) & sites[[group_col]] == grp, , drop = FALSE]
    agg <- stats::aggregate(mod_prob ~ position, data = sub, FUN = mean)
    if (nrow(agg) < 4L) {
      # Not enough points for loess: interpolate raw means, NA outside support
      return(stats::approx(agg$position, agg$mean_prob <- agg$mod_prob,
                           xout = grid, rule = 1)$y)
    }
    eff_span <- if (is.null(span)) max(0.15, min(0.75, 15 / nrow(agg))) else span
    fit <- tryCatch(
      suppressWarnings(stats::loess(mod_prob ~ position, data = agg, span = eff_span)),
      error = function(e) NULL
    )
    if (is.null(fit)) return(rep(NA_real_, length(grid)))
    pred <- suppressWarnings(stats::predict(fit, newdata = data.frame(position = grid)))
    pred[grid < min(agg$position) | grid > max(agg$position)] <- NA_real_
    pmin(pmax(pred, 0), 1)
  }

  all_pos <- sites$position[!is.na(sites$position)]
  grid <- seq(min(all_pos), max(all_pos), length.out = n_grid)

  v1 <- fit_one(groups[1L], grid)
  v2 <- fit_one(groups[2L], grid)
  delta <- v2 - v1
  sign <- ifelse(is.na(delta), NA_character_,
                 ifelse(delta > 0, "pos", ifelse(delta < 0, "neg", "zero")))
  data.frame(position = grid, delta = delta, sign = sign,
             stringsAsFactors = FALSE)
}
```

Note the `agg$mean_prob <- agg$mod_prob` inline assignment in the `approx` call
is fragile; write it as two statements instead:

```r
    if (nrow(agg) < 4L) {
      return(stats::approx(agg$position, agg$mod_prob, xout = grid, rule = 1)$y)
    }
```

- [ ] **Step 4: Run test to verify it passes**

Run: `Rscript -e 'devtools::test_file("tests/testthat/test-delta_track.R")'`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add R/delta_track.R tests/testthat/test-delta_track.R
git commit -m "feat: add shared-grid group delta computation"
```

---

### Task 6: Delta panel rendering and layout integration

**Files:**
- Modify: `R/delta_track.R` (add `.build_delta_panel()`)
- Modify: `R/plot_methylation.R` (`show_delta` arg; panel assembly ~ lines 578–602)
- Modify: `tests/testthat/test-delta_track.R`

**Interfaces:**
- Consumes: `.compute_group_delta()` (Task 5), `.DELTA_DIVERGING` (Task 1), `theme_ggmethylation()` (Task 1).
- Produces: `.build_delta_panel(delta_df, region_start, region_end)` → a ggplot; `plot_methylation(show_delta = FALSE)` appends the delta panel below the smooth panel when grouping yields exactly two groups.

- [ ] **Step 1: Write the failing test**

```r
test_that("build_delta_panel returns a ggplot", {
  df <- data.frame(position = seq(1000, 2000, length.out = 50),
                   delta = sin(seq(0, 3, length.out = 50)) * 0.4,
                   stringsAsFactors = FALSE)
  df$sign <- ifelse(df$delta >= 0, "pos", "neg")
  p <- ggmethylation:::.build_delta_panel(df, 1000, 2000)
  expect_s3_class(p, "ggplot")
})

test_that("plot_methylation adds delta panel for 2 groups when show_delta", {
  reads <- data.frame(
    read_name = c("r1","r2"), start = 1000L, end = 2000L, strand = "+",
    lane = 0L, mean_mod_prob = 0.5, group = c("1","2"),
    clip_side = NA_character_, sa_chrom = NA_character_, stringsAsFactors = FALSE
  )
  sites <- data.frame(
    read_name = rep(c("r1","r2"), each = 6),
    position = rep(seq(1100, 1900, length.out = 6), 2),
    mod_prob = c(rep(0.2,6), rep(0.8,6)), mod_code = "m",
    group = rep(c("1","2"), each = 6), stringsAsFactors = FALSE
  )
  md <- make_test_data(reads, sites); md$group_tag <- "HP"
  p_no  <- ggmethylation::plot_methylation(md, show_delta = FALSE, show_supplementary = FALSE)
  p_yes <- ggmethylation::plot_methylation(md, show_delta = TRUE,  show_supplementary = FALSE)
  expect_gt(length(p_yes$patches$plots), length(p_no$patches$plots))
})
```

(If `$patches$plots` introspection is brittle across patchwork versions, assert
`expect_s3_class(p_yes, "patchwork")` instead and rely on the panel-count being
covered by the message/skip test below.)

- [ ] **Step 2: Run test to verify it fails**

Run: `Rscript -e 'devtools::test_file("tests/testthat/test-delta_track.R")'`
Expected: FAIL — `.build_delta_panel` not found / `show_delta` arg unused.

- [ ] **Step 3: Implement `.build_delta_panel()`**

```r
# Render the signed delta as a diverging area around a zero baseline.
.build_delta_panel <- function(delta_df, region_start, region_end) {
  df <- delta_df[!is.na(delta_df$delta), , drop = FALSE]
  ggplot2::ggplot(df, ggplot2::aes(x = .data$position, y = .data$delta)) +
    ggplot2::geom_area(
      ggplot2::aes(fill = .data$sign),
      alpha = 0.85, na.rm = TRUE
    ) +
    ggplot2::geom_hline(yintercept = 0, colour = .DELTA_DIVERGING$zero,
                        linewidth = 0.4) +
    ggplot2::scale_fill_manual(
      values = c(pos = .DELTA_DIVERGING$pos, neg = .DELTA_DIVERGING$neg,
                 zero = .DELTA_DIVERGING$zero),
      guide = "none"
    ) +
    ggplot2::scale_x_continuous(labels = scales::comma_format()) +
    ggplot2::coord_cartesian(xlim = c(region_start, region_end)) +
    ggplot2::labs(x = "Genomic position (bp)", y = "Δ methylation\n(group2 - group1)") +
    theme_ggmethylation()
}
```

Note: `geom_area` with a `fill` grouping can leave gaps at sign changes; this is
acceptable for v1. If a filled-to-zero look is required later, switch to two
`geom_ribbon` layers (one clamped `ymin=0, ymax=pmax(delta,0)`, one
`ymin=pmin(delta,0), ymax=0`).

- [ ] **Step 4: Wire `show_delta` into `plot_methylation()`**

Add `show_delta = FALSE` to the signature (and to `.plot_multi_methylation()`;
for multi-sample, delta is out of scope — accept the arg and `message()` +
ignore). In the single-sample panel-assembly section (after `p_bottom` and the
optional `p_gene` are built, before `wrap_plots`):

```r
  p_delta <- NULL
  if (isTRUE(show_delta) && !is.null(data$group_tag)) {
    delta_df <- .compute_group_delta(data$sites, "group", span = smooth_span)
    if (!is.null(delta_df)) {
      # Break the delta over consensus deletions too, for visual consistency
      delta_df$group <- "delta"  # single-line identity for the break helper
      delta_df <- .apply_deletion_breaks(
        stats::setNames(delta_df[c("position", "delta", "group")],
                        c("position", "mean_prob", "group")),
        data$cigar_features, data$reads, "group", min_indel_size, show_cigar
      )
      names(delta_df)[names(delta_df) == "mean_prob"] <- "delta"
      delta_df$sign <- ifelse(is.na(delta_df$delta), NA_character_,
                              ifelse(delta_df$delta >= 0, "pos", "neg"))
      p_delta <- .build_delta_panel(delta_df, region_start, region_end)
      # The smooth panel is no longer the bottom-most; hide its x-axis title/labels
      p_bottom <- p_bottom +
        ggplot2::theme(axis.title.x = ggplot2::element_blank())
    }
  }
```

- [ ] **Step 5: Extend panel assembly and heights**

Replace the panel-assembly block so the delta panel is appended last:

```r
  panels   <- list(p_top, p_bottom)
  if (!is.null(p_gene)) panels <- c(list(p_gene), panels)
  if (!is.null(p_delta)) panels <- c(panels, list(p_delta))
  n_panels <- length(panels)

  if (is.null(panel_heights)) {
    heights <- c(
      if (!is.null(p_gene)) 0.08 else NULL,
      1,                       # reads
      0.25,                    # smooth
      if (!is.null(p_delta)) 0.2 else NULL
    )
  } else {
    if (length(panel_heights) != n_panels) {
      stop(sprintf(
        "`panel_heights` has length %d but there are %d panels.",
        length(panel_heights), n_panels
      ), call. = FALSE)
    }
    heights <- panel_heights
  }

  patchwork::wrap_plots(panels, ncol = 1, heights = heights)
```

- [ ] **Step 6: Run tests**

Run: `Rscript -e 'devtools::test_file("tests/testthat/test-delta_track.R")'`
Expected: PASS. Then run the full suite:

Run: `Rscript -e 'devtools::test()'`
Expected: PASS (confirm `panel_heights` length tests elsewhere still hold; update any that hard-code panel counts).

- [ ] **Step 7: Document + commit**

Add `@param show_delta` to `plot_methylation()` roxygen (note the exactly-two-groups
requirement and that it appends a bottom panel; update the `@param panel_heights`
note to mention the optional extra delta panel). Add NEWS bullet.
Run `devtools::document()`.

```bash
git add R/delta_track.R R/plot_methylation.R tests/testthat/test-delta_track.R NEWS.md man/
git commit -m "feat: add group-difference delta track (show_delta)"
```

---

## Self-Review Notes

- **Spec coverage:** A1 → Tasks 5–6; A2 → Tasks 2–3; B1 → Task 4; C2 → Task 1. All Milestone-1 features covered.
- **Type consistency:** `.compute_group_delta()` produces `position/delta/sign`; `.build_delta_panel()` consumes the same. `smooth_methylation()` `lower/upper` produced in Task 2, consumed in Task 3 and by `.add_ci_ribbon()`. `.classify_calls()` produced and consumed within Task 4.
- **Known fragilities flagged inline:** patchwork panel-count introspection (Task 6 Step 1), `geom_area` sign-change gaps (Task 6 Step 3), the delta deletion-break column renaming round-trip (Task 6 Step 4 — verify column order after `.apply_deletion_breaks`).
- **Ordering:** Task 1 must land first (palette constants + theme are consumed by Tasks 3 and 6). Tasks 2→3 are ordered. Task 4 is independent. Tasks 5→6 are ordered and depend on Task 1.
