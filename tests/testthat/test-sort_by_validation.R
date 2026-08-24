# Regression tests for .validate_sort_by().
#
# An unknown `sort_by` column used to yield NULL from `data$reads[[col]]`;
# order(NULL) returns integer(0), which silently dropped every read and
# rendered an empty plot with no error at all.

make_sort_md <- function() {
  reads <- data.frame(
    read_name = c("r1", "r2", "r3"),
    start     = c(1000L, 1200L, 1100L),
    end       = c(1800L, 2000L, 1900L),
    strand    = c("+", "-", "+"),
    stringsAsFactors = FALSE
  )
  sites <- data.frame(
    read_name = rep(c("r1", "r2", "r3"), each = 2L),
    position  = c(1100L, 1200L, 1300L, 1400L, 1150L, 1250L),
    mod_prob  = c(0.9, 0.1, 0.4, 0.6, 0.2, 0.8),
    mod_code  = "m",
    stringsAsFactors = FALSE
  )
  structure(
    list(
      reads           = reads,
      sites           = sites,
      insertion_sites = sites[0, , drop = FALSE],
      cigar_features  = NULL,
      region          = GenomicRanges::GRanges(
        "chr1", IRanges::IRanges(1000, 2000)
      ),
      mod_code        = "m",
      group_tag       = NULL
    ),
    class = "methylation_data"
  )
}

test_that("validate_sort_by accepts NULL and valid columns", {
  reads <- make_sort_md()$reads

  expect_silent(ggmethylation:::.validate_sort_by(NULL, reads))
  expect_silent(ggmethylation:::.validate_sort_by("start", reads))
  expect_silent(ggmethylation:::.validate_sort_by(c("start", "strand"), reads))
})

test_that("validate_sort_by rejects an unknown column and names it", {
  reads <- make_sort_md()$reads

  expect_error(
    ggmethylation:::.validate_sort_by("position", reads),
    'Unknown `sort_by` column: "position"'
  )
  # The message should list what IS available, so the fix is obvious.
  expect_error(
    ggmethylation:::.validate_sort_by("position", reads),
    "Available columns:"
  )
})

test_that("validate_sort_by reports every unknown column at once", {
  reads <- make_sort_md()$reads

  expect_error(
    ggmethylation:::.validate_sort_by(c("start", "nope", "alsonope"), reads),
    "Unknown `sort_by` columns:"
  )
})

test_that("validate_sort_by rejects non-character input", {
  reads <- make_sort_md()$reads

  expect_error(ggmethylation:::.validate_sort_by(1L, reads),
               "must be a character vector")
})

test_that("plot_methylation errors on an unknown sort_by instead of plotting nothing", {
  md <- make_sort_md()

  # "position" is a column of $sites, not $reads -- an easy and previously
  # silent mistake (the package vignette itself used to make it).
  expect_error(
    plot_methylation(md, sort_by = "position", show_supplementary = FALSE),
    "Unknown `sort_by` column"
  )
})

test_that("plot_methylation still accepts valid sort_by values", {
  md <- make_sort_md()

  for (col in c("start", "end", "strand", "read_name")) {
    p <- plot_methylation(md, sort_by = col, show_supplementary = FALSE)
    expect_true(inherits(p, "gg") || inherits(p, "patchwork"))
  }
})

test_that("mean_mod_prob is sortable even though it is computed, not stored", {
  md <- make_sort_md()
  expect_false("mean_mod_prob" %in% names(md$reads))

  # plot_methylation() adds it before sorting, so it must validate cleanly.
  p <- plot_methylation(md, sort_by = "mean_mod_prob",
                        show_supplementary = FALSE)
  expect_true(inherits(p, "gg") || inherits(p, "patchwork"))
})
