# Compute a signed methylation difference between exactly two groups over a
# shared position grid. Returns NULL (with a message) when != 2 groups.
.compute_group_delta <- function(sites, group_col, span = NULL, n_grid = 200L) {
  groups <- unique(sites[[group_col]])
  groups <- sort(groups[!is.na(groups)])
  if (length(groups) != 2L) {
    warning(
      "Delta track requires exactly two groups; found ", length(groups),
      ". Skipping the delta panel.",
      call. = FALSE
    )
    return(NULL)
  }

  all_pos <- sites$position[!is.na(sites$position)]
  grid <- seq(min(all_pos), max(all_pos), length.out = n_grid)

  smoothed <- smooth_methylation(sites, group_col = group_col, span = span, grid = grid)
  v1 <- smoothed$mean_prob[smoothed[[group_col]] == groups[1L]]
  v2 <- smoothed$mean_prob[smoothed[[group_col]] == groups[2L]]
  delta <- v2 - v1
  sign <- ifelse(is.na(delta), NA_character_,
                 ifelse(delta > 0, "pos", ifelse(delta < 0, "neg", "zero")))
  out <- data.frame(position = grid, delta = delta, sign = sign,
                    stringsAsFactors = FALSE)
  # The sorted group names, so callers can colour "pos"/"neg" by the group each
  # sign belongs to (delta = groups[2] - groups[1]).
  attr(out, "groups") <- groups
  out
}

# Render the signed delta as a diverging area around a zero baseline.
#
# `fill_pos`/`fill_neg` default to the standalone diverging palette; callers
# with a group palette in hand pass the two group colours instead, so the delta
# area matches the groups in the panels above (pos = group 2, neg = group 1).
.build_delta_panel <- function(delta_df, region_start, region_end,
                               y_label = "\u0394 mod.\nprobability",
                               fill_pos = .DELTA_DIVERGING$pos,
                               fill_neg = .DELTA_DIVERGING$neg) {
  df <- delta_df[!is.na(delta_df$delta), , drop = FALSE]
  ggplot2::ggplot(df, ggplot2::aes(x = .data$position, y = .data$delta)) +
    ggplot2::geom_area(
      ggplot2::aes(fill = .data$sign),
      alpha = 0.85, na.rm = TRUE
    ) +
    ggplot2::geom_hline(yintercept = 0, colour = .DELTA_DIVERGING$zero,
                        linewidth = 0.4) +
    ggplot2::scale_fill_manual(
      values = c(pos = fill_pos, neg = fill_neg,
                 zero = .DELTA_DIVERGING$zero),
      guide = "none"
    ) +
    ggplot2::scale_x_continuous(labels = scales::comma_format()) +
    ggplot2::coord_cartesian(xlim = c(region_start, region_end)) +
    ggplot2::labs(x = "Genomic position (bp)", y = y_label) +
    theme_ggmethylation() +
    .compact_y_title()
}
