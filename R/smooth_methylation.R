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
smooth_methylation <- function(sites, group_col = "group",
                               mod_code_col = NULL, span = 0.3) {
  if (!is.null(mod_code_col) && mod_code_col %in% names(sites)) {
    sites$.smooth_group <- paste(sites[[group_col]], sites[[mod_code_col]], sep = ":::")
    effective_group_col <- ".smooth_group"
  } else {
    effective_group_col <- group_col
  }

  out_cols <- c("position", "mean_prob", effective_group_col)

  if (is.null(sites) || nrow(sites) == 0L) {
    out <- data.frame(
      position  = numeric(0L),
      mean_prob = numeric(0L),
      group     = character(0L),
      stringsAsFactors = FALSE
    )
    names(out)[3L] <- effective_group_col
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

    if (nrow(agg) < 4L) {
      # Too few unique positions for loess; return raw means
      df <- agg
    } else {
      df <- tryCatch(
        suppressWarnings({
          fit  <- stats::loess(mean_prob ~ position, data = agg, span = span)
          grid <- seq(min(agg$position, na.rm = TRUE),
                      max(agg$position, na.rm = TRUE),
                      length.out = 200L)
          pred <- stats::predict(fit, newdata = data.frame(position = grid))
          data.frame(position = grid, mean_prob = pred)
        }),
        error = function(e) agg
      )
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
