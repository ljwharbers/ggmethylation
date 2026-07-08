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
      return(stats::approx(agg$position, agg$mod_prob, xout = grid, rule = 1)$y)
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
