# Tests for CI ribbon rendering (show_ci) and .insert_deletion_breaks() CI
# column handling in R/plot_methylation.R

# Local minimal methylation_data constructor (mirrors
# test-build_read_panel.R::make_test_data; not shared via a helper-*.R file,
# so duplicated here per testthat edition-3 auto-sourcing rules).
make_test_data <- function(reads_df, sites_df,
                            region_start = 1000L, region_end = 2000L) {
  gr <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges   = IRanges::IRanges(start = region_start, end = region_end)
  )
  structure(
    list(
      reads        = reads_df,
      sites        = sites_df,
      region       = gr,
      mod_code     = "m",
      group_tag    = NULL,
      cigar_features = data.frame(
        read_name  = character(0),
        type       = character(0),
        ref_start  = integer(0),
        ref_end    = integer(0),
        length     = integer(0),
        stringsAsFactors = FALSE
      )
    ),
    class = "methylation_data"
  )
}

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

test_that("insert_deletion_breaks sentinel rows have NA lower/upper", {
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
  sentinels <- out[out$position %in% c(139.5, 210.5), ]
  expect_equal(nrow(sentinels), 2L)
  expect_true(all(is.na(sentinels$lower)))
  expect_true(all(is.na(sentinels$upper)))
})

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

test_that("plot_methylation omits ribbon layer when show_ci = FALSE", {
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
  p <- ggmethylation::plot_methylation(md, show_ci = FALSE, show_supplementary = FALSE)
  expect_s3_class(p, "patchwork")
  p_bottom <- p[[2]]
  geom_classes <- vapply(p_bottom$layers, function(l) class(l$geom)[1], character(1L))
  expect_false("GeomRibbon" %in% geom_classes)
})

test_that("grouped multi-code ribbon keeps (group, mod_code) polygons separate", {
  # Regression test for the group_aes fix in the grouped + multi-mod-code
  # branch of plot_methylation() (colour = group, linetype = mod_code).
  # Without passing group_aes = interaction(group, mod_code) to
  # .add_ci_ribbon(), ggplot2's default grouping collapses to fill = group
  # only, so the "m" and "h" ribbon bands for a given group get drawn as a
  # single self-crossing polygon (x jumps backward where one band's grid
  # ends and the other's begins). 2 groups x 2 mod codes, 7 unique positions
  # each, so every (group, mod_code) combination gets a real loess CI
  # (smooth_methylation() requires >= 4 unique positions per sub-group).
  # mod_prob values are deliberately non-monotonic (not near-perfectly
  # linear) so the loess fit has non-degenerate residual variance and
  # se.fit/lower/upper don't come out NaN for any combination.
  reads <- data.frame(
    read_name = c("r1", "r2"), start = c(1000L, 1000L), end = c(2000L, 2000L),
    strand = c("+", "+"), lane = 0L, mean_mod_prob = 0.5,
    group = c("1", "2"), clip_side = NA_character_, sa_chrom = NA_character_,
    stringsAsFactors = FALSE
  )
  positions <- seq(1100, 1900, length.out = 7)
  sites <- data.frame(
    read_name = rep(c("r1", "r2"), each = 14),
    position = rep(rep(positions, 2), 2),
    mod_prob = c(0.1, 0.3, 0.2, 0.5, 0.4, 0.7, 0.6,        # r1, mod_code "m"
                 0.2, 0.5, 0.3, 0.6, 0.4, 0.8, 0.5,        # r1, mod_code "h"
                 0.15, 0.4, 0.25, 0.55, 0.35, 0.75, 0.5,   # r2, mod_code "m"
                 0.25, 0.55, 0.3, 0.65, 0.4, 0.85, 0.55),  # r2, mod_code "h"
    mod_code = rep(rep(c("m", "h"), each = 7), 2),
    group = rep(c("1", "2"), each = 14),
    stringsAsFactors = FALSE
  )
  md <- make_test_data(reads, sites)
  md$group_tag <- "HP"

  p <- ggmethylation::plot_methylation(md, show_ci = TRUE, show_supplementary = FALSE)
  p_bottom <- p[[2]]
  geom_classes <- vapply(p_bottom$layers, function(l) class(l$geom)[1], character(1L))
  ribbon_idx <- which(geom_classes == "GeomRibbon")
  expect_length(ribbon_idx, 1L)

  built <- ggplot2::ggplot_build(p_bottom)
  ribbon_data <- built$data[[ribbon_idx]]

  # With the fix, ggplot2 assigns 4 distinct built groups (one per
  # group x mod_code combination). Without it, this collapses to 2.
  built_groups <- unique(ribbon_data$group)
  expect_length(built_groups, 4L)

  # Within each built group, x (as laid out by ggplot2) must be
  # non-decreasing; a backward jump indicates a self-crossing polygon.
  for (g in built_groups) {
    x <- ribbon_data$x[ribbon_data$group == g]
    expect_true(all(diff(x) >= 0),
               info = paste("backward x jump found in built group", g))
  }
})
