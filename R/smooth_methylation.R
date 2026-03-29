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
#' @param span Loess smoothing span (default 0.3).
#'
#' @return A data.frame with columns `position`, `mean_prob`, and the
#'   grouping column. Positions are a regular grid of 200 points with
#'   loess-predicted values. Groups with fewer than 4 unique positions
#'   return raw per-site means instead.
#'
#' @keywords internal
smooth_methylation <- function(sites, group_col = "group", span = 0.3) {
  out_cols <- c("position", "mean_prob", group_col)

  if (is.null(sites) || nrow(sites) == 0L) {
    out <- data.frame(
      position = numeric(0L),
      mean_prob = numeric(0L),
      group = character(0L),
      stringsAsFactors = FALSE
    )
    names(out)[3L] <- group_col
    return(out)
  }

  groups <- unique(sites[[group_col]])
  groups <- groups[!is.na(groups)]
  result_list <- vector("list", length(groups))

  for (k in seq_along(groups)) {
    grp <- groups[k]
    sub <- sites[sites[[group_col]] == grp, , drop = FALSE]

    # Compute per-site mean modification probability
    agg <- stats::aggregate(
      mod_prob ~ position,
      data = sub,
      FUN = mean
    )
    names(agg) <- c("position", "mean_prob")

    if (nrow(agg) < 4L) {
      # Too few unique positions for loess; return raw means
      df <- agg
    } else {
      fit <- stats::loess(mean_prob ~ position, data = agg, span = span)
      grid <- seq(min(agg$position), max(agg$position), length.out = 200L)
      pred <- stats::predict(fit, newdata = data.frame(position = grid))
      df <- data.frame(position = grid, mean_prob = pred)
    }

    df[[group_col]] <- grp
    result_list[[k]] <- df
  }

  out <- do.call(rbind, result_list)
  rownames(out) <- NULL
  out[, out_cols, drop = FALSE]
}
