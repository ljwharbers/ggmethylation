# Tests for binarized call mode (call_mode = "binary") in the read panel:
# .classify_calls() and plot_methylation(call_mode = "binary", ...).

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

test_that("binary call classification splits at threshold", {
  probs <- c(0.1, 0.49, 0.5, 0.9)
  cls <- ggmethylation:::.classify_calls(probs, threshold = 0.5, ambiguous = NULL)
  expect_equal(cls, c("unmethylated", "unmethylated", "methylated", "methylated"))
})

test_that("ambiguous band labels near-threshold calls", {
  probs <- c(0.1, 0.45, 0.5, 0.55, 0.9)
  cls <- ggmethylation:::.classify_calls(probs, threshold = 0.5, ambiguous = 0.1)
  expect_equal(cls, c("unmethylated", "ambiguous", "methylated", "ambiguous", "methylated"))
})

test_that("plot_methylation renders in binary mode", {
  reads <- data.frame(
    read_name = "r1", start = 1000L, end = 2000L, strand = "+",
    lane = 0L, mean_mod_prob = 0.5, clip_side = NA_character_,
    sa_chrom = NA_character_, stringsAsFactors = FALSE
  )
  sites <- data.frame(
    read_name = "r1", position = seq(1100, 1900, length.out = 5),
    mod_prob = c(0.1, 0.3, 0.5, 0.7, 0.9), mod_code = "m",
    stringsAsFactors = FALSE
  )
  md <- make_test_data(reads, sites)
  p <- ggmethylation::plot_methylation(md, call_mode = "binary",
                                       show_supplementary = FALSE)
  expect_true(inherits(p, "ggplot") || inherits(p, "patchwork"))
})
