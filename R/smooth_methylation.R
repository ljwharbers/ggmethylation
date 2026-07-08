#' Smooth methylation probabilities per group using loess
#'
#' Computes loess-smoothed mean modification probability per group across
#' genomic positions. Per-site means are calculated first, then a loess curve
#' is fit and predicted on a regular grid of 200 points spanning the range
#' of input positions.
#'
#' @param sites A data.frame with columns `position`, `mod_prob`, and a
#'   grouping column (specified by `group_col`).
#' @param group_col Name of the grouping column (default `"group"`).
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
smooth_methylation <- function(sites, group_col = "group",
                               mod_code_col = NULL, span = NULL, grid = NULL) {
  if (!is.null(mod_code_col) && mod_code_col %in% names(sites)) {
    sites$.smooth_group <- paste(sites[[group_col]], sites[[mod_code_col]], sep = ":::")
    effective_group_col <- ".smooth_group"
  } else {
    effective_group_col <- group_col
  }

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

  groups <- unique(sites[[effective_group_col]])
  groups <- groups[!is.na(groups)]
  result_list <- vector("list", length(groups))

  for (k in seq_along(groups)) {
    grp <- groups[k]
    sub <- sites[sites[[effective_group_col]] == grp, , drop = FALSE]

    # Compute per-site mean modification probability
    agg <- stats::aggregate(
      mod_prob ~ position,
      data = sub,
      FUN = mean
    )
    names(agg) <- c("position", "mean_prob")

    effective_span <- if (is.null(span)) {
      max(0.15, min(0.75, 15 / nrow(agg)))
    } else {
      span
    }

    if (is.null(grid)) {
      if (nrow(agg) < 4L) {
        # Too few unique positions for loess; return raw means, no CI
        df <- agg
        df$lower <- NA_real_
        df$upper <- NA_real_
      } else {
        df <- tryCatch(
          suppressWarnings({
            fit      <- stats::loess(mean_prob ~ position, data = agg, span = effective_span)
            own_grid <- seq(min(agg$position, na.rm = TRUE),
                        max(agg$position, na.rm = TRUE),
                        length.out = 200L)
            pr   <- stats::predict(fit, newdata = data.frame(position = own_grid), se = TRUE)
            lower <- pmin(pmax(pr$fit - 1.96 * pr$se.fit, 0), 1)
            upper <- pmin(pmax(pr$fit + 1.96 * pr$se.fit, 0), 1)
            data.frame(position = own_grid, mean_prob = pr$fit,
                       lower = lower, upper = upper)
          }),
          error = function(e) {
            agg$lower <- NA_real_
            agg$upper <- NA_real_
            agg
          }
        )
      }
    } else {
      # Shared external grid: evaluate every group at the same positions so
      # they can be compared/subtracted index-wise (e.g. delta track).
      if (nrow(agg) < 4L) {
        mean_prob <- stats::approx(agg$position, agg$mean_prob, xout = grid, rule = 1)$y
        df <- data.frame(
          position  = grid,
          mean_prob = mean_prob,
          lower     = NA_real_,
          upper     = NA_real_
        )
      } else {
        df <- tryCatch(
          suppressWarnings({
            fit <- stats::loess(mean_prob ~ position, data = agg, span = effective_span)
            pr  <- stats::predict(fit, newdata = data.frame(position = grid), se = TRUE)
            mean_prob <- pr$fit
            lower <- pmin(pmax(pr$fit - 1.96 * pr$se.fit, 0), 1)
            upper <- pmin(pmax(pr$fit + 1.96 * pr$se.fit, 0), 1)
            out_of_support <- grid < min(agg$position, na.rm = TRUE) |
              grid > max(agg$position, na.rm = TRUE)
            mean_prob[out_of_support] <- NA_real_
            lower[out_of_support] <- NA_real_
            upper[out_of_support] <- NA_real_
            data.frame(position = grid, mean_prob = mean_prob,
                       lower = lower, upper = upper)
          }),
          error = function(e) {
            data.frame(
              position  = grid,
              mean_prob = NA_real_,
              lower     = NA_real_,
              upper     = NA_real_
            )
          }
        )
      }
    }

    df[[effective_group_col]] <- grp
    result_list[[k]] <- df
  }

  out <- do.call(rbind, result_list)
  rownames(out) <- NULL
  out <- out[, out_cols, drop = FALSE]

  # Split composite key back into original columns
  if (!is.null(mod_code_col) && mod_code_col %in% names(sites)) {
    parts <- strsplit(out[[".smooth_group"]], ":::", fixed = TRUE)
    out[[group_col]]    <- vapply(parts, `[[`, character(1L), 1L)
    out[[mod_code_col]] <- vapply(parts, `[[`, character(1L), 2L)
    out[[".smooth_group"]] <- NULL
  }

  out
}
