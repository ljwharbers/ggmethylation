# Tests for plot_insertion_locus().
#
# The exported renderer was only exercised incidentally from
# test-insertion_loci.R (a list_insertion_loci() test file). Nothing targeted
# the stitched-coordinate panel, `flank`, `include_noncarriers` or
# `show_smoothed`.
#
# These run against the committed MEG3 fixture, which carries 54 insertion
# sites -- see data-raw/build_test_fixture.sh.

bam_path <- testthat::test_path("fixtures", "hg002_fiberseq_MEG3.bam")
region   <- "chr14:95055000-95070000"

meg3_md <- function() {
  read_methylation(bam_path, region, mod_code = "m")
}

first_locus <- function(md) {
  loci <- list_insertion_loci(md, tol_pos = 10L, tol_len = 0.20,
                              min_reads = 2L)
  testthat::skip_if(nrow(loci) == 0L, "no insertion loci in the fixture")
  loci$locus_id[1L]
}

test_that("the fixture yields insertion loci to plot", {
  md   <- meg3_md()
  loci <- list_insertion_loci(md, tol_pos = 10L, tol_len = 0.20,
                              min_reads = 2L)

  expect_true(is.data.frame(loci))
  expect_gt(nrow(loci), 0L)
  expect_true(all(c("locus_id", "anchor_pos", "n_carriers",
                    "n_noncarriers", "median_length") %in% names(loci)))
  expect_true(all(loci$n_carriers >= 2L))
})

test_that("plot_insertion_locus returns a patchwork by default", {
  md <- meg3_md()
  p  <- plot_insertion_locus(md, first_locus(md))

  expect_true(inherits(p, "patchwork") || inherits(p, "gg"))
})

test_that("show_smoothed = FALSE returns a single ggplot, TRUE adds a panel", {
  md  <- meg3_md()
  lid <- first_locus(md)

  p_bare <- plot_insertion_locus(md, lid, show_smoothed = FALSE)
  expect_s3_class(p_bare, "ggplot")
  expect_false(inherits(p_bare, "patchwork"))

  p_full <- plot_insertion_locus(md, lid, show_smoothed = TRUE)
  expect_s3_class(p_full, "patchwork")
})

# The read polygons are layer 1; each row is one stitched sub-segment and
# carries the carrier_status of the read it came from.
polygon_data <- function(p) p$layers[[1L]]$data

test_that("include_noncarriers = FALSE drops the non-carrier reads", {
  md  <- meg3_md()
  lid <- first_locus(md)

  p_with    <- plot_insertion_locus(md, lid, include_noncarriers = TRUE,
                                    show_smoothed = FALSE)
  p_without <- plot_insertion_locus(md, lid, include_noncarriers = FALSE,
                                    show_smoothed = FALSE)

  expect_true("non-carrier" %in% polygon_data(p_with)$carrier_status)
  expect_false("non-carrier" %in% polygon_data(p_without)$carrier_status)

  # Carriers are unaffected by the switch.
  n_carrier <- function(p) sum(polygon_data(p)$carrier_status == "carrier")
  expect_gt(n_carrier(p_without), 0L)
  expect_identical(n_carrier(p_with), n_carrier(p_without))

  expect_lt(nrow(polygon_data(p_without)), nrow(polygon_data(p_with)))
})

test_that("flank widens the stitched x-axis by twice the added flank", {
  md  <- meg3_md()
  lid <- first_locus(md)

  span <- function(fl) {
    p <- plot_insertion_locus(md, lid, flank = fl, show_smoothed = FALSE)
    diff(ggplot2::ggplot_build(p)$layout$panel_params[[1L]]$x.range)
  }

  s100 <- span(100L)
  s400 <- span(400L)

  expect_gt(s400, s100)
  # flank is applied on both sides, so +300 bp per side => +600 bp total.
  expect_equal(s400 - s100, 600, tolerance = 1)
})

test_that("an unknown locus_id errors and points at list_insertion_loci", {
  md <- meg3_md()

  expect_error(
    plot_insertion_locus(md, "INS_chr14_1_1bp"),
    "not found"
  )
  expect_error(
    plot_insertion_locus(md, "INS_chr14_1_1bp"),
    "list_insertion_loci"
  )
})

test_that("a non-methylation_data input is rejected", {
  expect_error(plot_insertion_locus(list(), "whatever"),
               "methylation_data")
})

# Pull the continuous probability ramp out of a panel and evaluate it, so the
# comparison is on rendered colours rather than on how the scale was declared.
prob_ramp <- function(p, at = c(0, 0.5, 1)) {
  sc <- Filter(function(s) "colour" %in% s$aesthetics, p$scales$scales)
  testthat::expect_length(sc, 1L)
  sc[[1L]]$palette(at)
}

test_that("locus colours match plot_methylation's probability ramp (issue #18)", {
  # Issue #18 asked for the locus plot to use the same colours as the main
  # plot. Compare the two ramps directly instead of hard-coding hex values in
  # two places.
  md  <- meg3_md()
  lid <- first_locus(md)

  p_locus <- plot_insertion_locus(md, lid, show_smoothed = FALSE)
  p_main  <- plot_methylation(md, show_supplementary = FALSE)

  expect_identical(prob_ramp(p_locus), prob_ramp(p_main[[1L]]))

  # And they are the documented defaults.
  expect_identical(prob_ramp(p_locus, c(0, 1)), c("#BDBDBD", "#C62828"))
})

test_that("colour_low/colour_high override the locus ramp", {
  md  <- meg3_md()
  lid <- first_locus(md)

  p <- plot_insertion_locus(md, lid, show_smoothed = FALSE,
                            colour_low = "white", colour_high = "black")

  expect_identical(prob_ramp(p, c(0, 1)), c("#FFFFFF", "#000000"))
})
