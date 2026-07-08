# Compute a signed methylation difference between exactly two groups over a
# shared position grid. Returns NULL (with a message) when != 2 groups.
.compute_group_delta <- function(sites, group_col, span = NULL, n_grid = 200L) {
  groups <- unique(sites[[group_col]])
  groups <- sort(groups[!is.na(groups)])
  if (length(groups) != 2L) {
    message("Delta track requires exactly two groups; skipping.")
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
  data.frame(position = grid, delta = delta, sign = sign, stringsAsFactors = FALSE)
}
