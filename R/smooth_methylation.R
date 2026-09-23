# Convert continuous mod_prob values to 0/1 calls for binary-mode aggregation.
#
# Reuses .classify_calls() so the threshold semantics are defined in exactly one
# place (notably: a probability exactly equal to `threshold` is "methylated").
# Ambiguous sites and NA probabilities are dropped, so a downstream per-position
# mean() becomes "fraction methylated among confident calls".
#
# With `ambiguous = NULL` (the default, matching `call_ambiguous`) no site is
# ever labelled ambiguous, so nothing is dropped and the denominator equals
# full coverage — a pure hard threshold.
.binarize_sites = function(sites, threshold = 0.5, ambiguous = NULL) {
  if (is.null(sites) || nrow(sites) == 0L) return(sites)
  cls = .classify_calls(sites$mod_prob, threshold, ambiguous)
  keep = !is.na(cls) & cls != "ambiguous"
  sites = sites[keep, , drop = FALSE]
  sites$mod_prob = as.numeric(cls[keep] == "methylated")
  rownames(sites) = NULL
  sites
}

# Adaptive loess span targeting ~15 data points per local fit.
.adaptive_span = function(n_positions) max(0.15, min(0.75, 15 / n_positions))

# Loess fit of `mean_prob ~ position`, predicted at `grid` with a 95% band
# (+/- 1.96 se) clamped to [0, 1]. With `mask_outside`, grid points outside the
# data's own range are set to NA instead of being extrapolated.
.loess_predict = function(agg, grid, span, mask_outside = FALSE) {
  fit = stats::loess(mean_prob ~ position, data = agg, span = span)
  pr  = stats::predict(fit, newdata = data.frame(position = grid), se = TRUE)
  out = data.frame(
    position  = grid,
    mean_prob = pr$fit,
    lower     = pmin(pmax(pr$fit - 1.96 * pr$se.fit, 0), 1),
    upper     = pmin(pmax(pr$fit + 1.96 * pr$se.fit, 0), 1)
  )
  if (mask_outside) {
    outside = grid < min(agg$position) | grid > max(agg$position)
    out[outside, c("mean_prob", "lower", "upper")] = NA_real_
  }
  out
}

# Smooth one group's sites: per-position means, then loess on `grid` (or on a
# 200-point grid over the group's own range when `grid` is NULL). Groups with
# fewer than 4 unique positions return raw means (own grid) or a linear
# interpolation (shared grid), without a confidence band.
.smooth_group = function(sites, span = NULL, grid = NULL) {
  agg = stats::aggregate(mod_prob ~ position, data = sites, FUN = mean)
  names(agg) = c("position", "mean_prob")
  if (is.null(span)) span = .adaptive_span(nrow(agg))
  own_grid = is.null(grid)

  if (nrow(agg) < 4L) {
    if (own_grid) return(cbind(agg, lower = NA_real_, upper = NA_real_))
    interp = stats::approx(agg$position, agg$mean_prob, xout = grid, rule = 1)$y
    return(data.frame(position = grid, mean_prob = interp,
                      lower = NA_real_, upper = NA_real_))
  }
  if (own_grid) grid = seq(min(agg$position), max(agg$position), length.out = 200L)

  tryCatch(
    suppressWarnings(.loess_predict(agg, grid, span, mask_outside = !own_grid)),
    error = function(e) {
      if (own_grid) cbind(agg, lower = NA_real_, upper = NA_real_) else
        data.frame(position = grid, mean_prob = NA_real_, lower = NA_real_, upper = NA_real_)
    }
  )
}

#' Fit a loess curve to x/y data and predict on a 200-point grid
#'
#' Single-group shortcut for [smooth_methylation()], used by
#' [plot_insertion_locus()].
#'
#' @param x Numeric vector of positions.
#' @param y Numeric vector of values (same length as `x`).
#' @param span Numeric or NULL. Loess span. When NULL, uses adaptive rule
#'   `max(0.15, min(0.75, 15 / n_unique_x))`.
#'
#' @return A data.frame with columns `position` and `mean_prob`. Returns
#'   raw per-x means when fewer than 4 unique x values are present.
#'
#' @keywords internal
.smooth_xy = function(x, y, span = NULL) {
  .smooth_group(data.frame(position = x, mod_prob = y), span)[c("position", "mean_prob")]
}

#' Smooth methylation probabilities per group using loess
#'
#' Computes loess-smoothed mean modification probability per group across
#' genomic positions. Per-site means are calculated first, then a loess curve
#' is fit and predicted on a regular grid of 200 points spanning the range
#' of input positions.
#'
#' @param sites A data.frame with columns `position`, `mod_prob`, and a
#'   grouping column (specified by `group_col`).
#' @param group_col Name of the grouping column (default `"group"`). Sites
#'   whose group is `NA` are dropped.
#' @param mod_code_col Optional name of a modification-code column; when set,
#'   one line is fitted per group and code.
#' @param span Loess smoothing span. When `NULL` (default), an adaptive span is
#'   computed per group as `max(0.15, min(0.75, 15 / n_unique_sites))`, targeting
#'   approximately 15 data points per local fit regardless of region size or CpG
#'   density. Pass an explicit numeric value (e.g. `0.3`) to use a fixed span.
#' @param grid Optional numeric vector of shared evaluation positions. When
#'   `NULL` (default), each group is evaluated on its own 200-point grid
#'   spanning that group's own position range (the original behavior). When
#'   supplied, every group is instead evaluated at exactly these positions
#'   (in this order), which is required to compare/subtract groups index-wise
#'   (e.g. for a delta track). Predicted `mean_prob`/`lower`/`upper` values
#'   are set to `NA_real_` wherever a grid position falls outside that
#'   group's own position range, since loess `predict()` would otherwise
#'   extrapolate. Groups with fewer than 4 unique positions are interpolated
#'   via `stats::approx(..., rule = 1)` instead (also `NA` outside range),
#'   with `lower`/`upper` set to `NA_real_`.
#'
#' @return A data.frame with columns `position`, `mean_prob`, `lower`,
#'   `upper`, and the grouping column. When `grid` is `NULL`, positions are a
#'   regular grid of 200 points per group, spanning that group's own position
#'   range, with loess-predicted values. `lower`/`upper` are the loess fit
#'   `± 1.96 * se.fit`, clamped to `[0, 1]`. Groups with fewer than 4 unique
#'   positions return raw per-site means instead, with `lower`/`upper` set
#'   to `NA_real_`. When `grid` is supplied, positions are the shared `grid`
#'   values for every group instead (see `grid` above).
#'
#' @keywords internal
smooth_methylation = function(sites, group_col = "group",
                               mod_code_col = NULL, span = NULL, grid = NULL) {
  # One line per group, or per (group, mod_code) pair when mod_code_col is set.
  # Sites without a group are not smoothed.
  id_cols = group_col
  if (!is.null(mod_code_col) && mod_code_col %in% names(sites)) {
    id_cols = c(group_col, mod_code_col)
  }
  out_cols = c("position", "mean_prob", "lower", "upper", id_cols)
  if (!is.null(sites)) sites = sites[!is.na(sites[[group_col]]), , drop = FALSE]

  if (is.null(sites) || nrow(sites) == 0L) {
    out = data.frame(position = numeric(0L), mean_prob = numeric(0L),
                     lower = numeric(0L), upper = numeric(0L))
    for (col in id_cols) out[[col]] = character(0L)
    return(out)
  }

  keys = unique(sites[id_cols])
  pieces = lapply(seq_len(nrow(keys)), function(k) {
    in_key = Reduce(`&`, lapply(id_cols, function(col) sites[[col]] == keys[[col]][k]))
    df = .smooth_group(sites[in_key, , drop = FALSE], span, grid)
    for (col in id_cols) df[[col]] = keys[[col]][k]
    df
  })
  out = do.call(rbind, pieces)
  rownames(out) = NULL
  out[, out_cols, drop = FALSE]
}
