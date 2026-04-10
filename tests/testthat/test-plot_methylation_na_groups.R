# Tests for grouped plotting when reads with missing group labels are retained

make_grouped_plot_data <- function() {
  gr <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges   = IRanges::IRanges(start = 100L, end = 300L)
  )

  reads <- data.frame(
    read_name = c("r1", "r2", "r3"),
    start     = c(100L, 120L, 180L),
    end       = c(170L, 190L, 260L),
    strand    = c("+", "+", "-"),
    group     = c("1", NA, "2"),
    stringsAsFactors = FALSE
  )

  sites <- data.frame(
    read_name = c("r1", "r1", "r2", "r2", "r3", "r3"),
    position  = c(105L, 130L, 125L, 150L, 200L, 230L),
    mod_prob  = c(0.1, 0.3, 0.4, 0.6, 0.5, 0.7),
    mod_code  = rep("m", 6L),
    group     = c("1", "1", NA, NA, "2", "2"),
    stringsAsFactors = FALSE
  )

  cigar_features <- data.frame(
    type        = c("D", "D", "D"),
    ref_start   = c(135L, 140L, 210L),
    ref_end     = c(145L, 150L, 220L),
    query_start = c(NA_integer_, NA_integer_, NA_integer_),
    query_end   = c(NA_integer_, NA_integer_, NA_integer_),
    length      = c(11L, 11L, 11L),
    read_name   = c("r1", "r2", "r3"),
    stringsAsFactors = FALSE
  )

  structure(
    list(
      reads          = reads,
      sites          = sites,
      region         = gr,
      mod_code       = "m",
      group_tag      = "HP",
      snv_position   = NULL,
      sequences      = setNames(character(0), character(0)),
      cigars         = setNames(character(0), character(0)),
      cigar_features = cigar_features
    ),
    class = "methylation_data"
  )
}

test_that("plot_methylation handles retained NA groups when showing CIGAR deletions", {
  md <- make_grouped_plot_data()

  p <- expect_no_error(
    plot_methylation(md, show_cigar = TRUE, min_indel_size = 1L)
  )

  expect_s3_class(p, "patchwork")

  top_panel <- p$patches$plots[[1L]]
  poly_data <- ggplot2::layer_data(top_panel, 1L)

  expect_true(nrow(poly_data) > 0L)
  expect_true("grey50" %in% poly_data$fill)
})