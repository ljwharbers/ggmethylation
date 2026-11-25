# Compute a signed methylation difference between exactly two groups over a
# shared position grid. Returns NULL (with a message) when != 2 groups.
.compute_group_delta = function(sites, group_col, span = NULL, n_grid = 200L) {
  groups = unique(sites[[group_col]])
  groups = sort(groups[!is.na(groups)])
  if (length(groups) != 2L) {
    message("Delta track requires exactly two groups; skipping.")
    return(NULL)
  }

  all_pos = sites$position[!is.na(sites$position)]
  grid = seq(min(all_pos), max(all_pos), length.out = n_grid)

  smoothed = smooth_methylation(sites, group_col = group_col, span = span, grid = grid)
  v1 = smoothed$mean_prob[smoothed[[group_col]] == groups[1L]]
  v2 = smoothed$mean_prob[smoothed[[group_col]] == groups[2L]]
  delta = v2 - v1
  out = data.frame(position = grid, delta = delta, sign = .delta_sign(delta),
                    stringsAsFactors = FALSE)
  # The sorted group names, so callers can colour "pos"/"neg" by the group each
  # sign belongs to (delta = groups[2] - groups[1]).
  attr(out, "groups") = groups
  out
}

# "pos" / "neg" / "zero" (NA stays NA), used to fill the delta area.
.delta_sign = function(delta) {
  ifelse(is.na(delta), NA_character_,
         ifelse(delta > 0, "pos", ifelse(delta < 0, "neg", "zero")))
}

# Ribbon data for the delta area: one "pos" band (0 .. max(delta, 0)) and one
# "neg" band (min(delta, 0) .. 0) per stretch of non-NA delta values. The NA
# rows that .insert_deletion_breaks() adds split the stretches, so the area is
# interrupted over consensus deletions exactly like the smooth lines. Points
# where the delta changes sign are added at the interpolated zero crossing so
# the two bands meet without a gap, and bands that are zero-height throughout
# are dropped. Each row keeps the signed `delta` it was built from.
.delta_ribbon_data = function(delta_df) {
  df = data.frame(position = delta_df$position, delta = delta_df$delta)
  df = df[order(df$position), , drop = FALSE]
  n = nrow(df)
  if (n > 1L) {
    d = df$delta
    x = df$position
    i = which(!is.na(d[-n]) & !is.na(d[-1L]) & d[-n] * d[-1L] < 0)
    if (length(i) > 0L) {
      x0 = x[i] - d[i] * (x[i + 1L] - x[i]) / (d[i + 1L] - d[i])
      df = rbind(df, data.frame(position = x0, delta = 0))
      df = df[order(df$position, is.na(df$delta)), , drop = FALSE]
    }
  }
  piece = cumsum(is.na(df$delta))
  keep  = !is.na(df$delta)
  df    = df[keep, , drop = FALSE]
  piece = piece[keep]

  bands = rbind(
    data.frame(position = df$position, delta = df$delta, ymin = 0,
               ymax = pmax(df$delta, 0), sign = rep("pos", nrow(df)),
               piece = piece, stringsAsFactors = FALSE),
    data.frame(position = df$position, delta = df$delta,
               ymin = pmin(df$delta, 0), ymax = 0, sign = rep("neg", nrow(df)),
               piece = piece, stringsAsFactors = FALSE)
  )
  bands$run = paste(bands$sign, bands$piece, sep = "_")
  if (nrow(bands) == 0L) return(bands)
  height = tapply(bands$ymax - bands$ymin, bands$run, max)
  bands = bands[bands$run %in% names(height)[height > 0], , drop = FALSE]
  rownames(bands) = NULL
  bands
}

# Render the signed delta as a diverging area around a zero baseline.
#
# `fill_pos`/`fill_neg` default to the standalone diverging palette; callers
# with a group palette in hand pass the two group colours instead, so the delta
# area matches the groups in the panels above (pos = group 2, neg = group 1).
.build_delta_panel = function(delta_df, region_start, region_end,
                               y_label = "\u0394 methylation",
                               fill_pos = .DELTA_DIVERGING$pos,
                               fill_neg = .DELTA_DIVERGING$neg) {
  bands = .delta_ribbon_data(delta_df)
  ggplot2::ggplot(bands, ggplot2::aes(x = .data$position)) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = .data$ymin, ymax = .data$ymax,
                   fill = .data$sign, group = .data$run),
      alpha = 0.85, colour = NA
    ) +
    ggplot2::geom_hline(yintercept = 0, colour = .DELTA_DIVERGING$zero,
                        linewidth = 0.4) +
    ggplot2::scale_fill_manual(
      values = c(pos = fill_pos, neg = fill_neg),
      guide = "none"
    ) +
    ggplot2::scale_x_continuous(labels = scales::comma_format()) +
    ggplot2::coord_cartesian(xlim = c(region_start, region_end)) +
    ggplot2::labs(x = "Genomic position (bp)", y = y_label) +
    theme_ggmethylation()
}
