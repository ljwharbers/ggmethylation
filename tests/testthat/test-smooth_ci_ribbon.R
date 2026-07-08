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
