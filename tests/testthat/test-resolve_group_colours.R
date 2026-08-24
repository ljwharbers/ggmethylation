# Regression tests for .resolve_group_colours().
#
# The default group palette is named "1"/"2" because that is what an HP
# haplotype tag yields. Any other grouping -- SNV genotype ("REF"/"ALT"), a
# custom BAM tag -- shared no names with it, so scale_*_manual() matched
# nothing: every group fell through to na.value grey and ggplot2 warned
# "No shared levels found ...". SNV grouping is a documented feature, so it
# lost its colours by default.

test_that("an exactly matching palette is returned unchanged", {
  pal <- c("1" = "#0072B2", "2" = "#E69F00")
  expect_identical(
    ggmethylation:::.resolve_group_colours(pal, c("1", "2", "1")),
    pal
  )
})

test_that("the default palette adapts to non-haplotype group names", {
  res <- ggmethylation:::.resolve_group_colours(
    ggmethylation:::.GROUP_PALETTE_DEFAULT,
    c("REF", "ALT", "REF")
  )

  expect_identical(sort(names(res)), c("ALT", "REF"))
  # Assigned over .ordered_plot_groups() order, i.e. sorted: ALT then REF.
  expect_identical(unname(res), unname(ggmethylation:::.GROUP_PALETTE_DEFAULT))
  expect_identical(res[["ALT"]], "#0072B2")
  expect_identical(res[["REF"]], "#E69F00")
})

test_that("adapting the default palette is silent", {
  expect_silent(
    ggmethylation:::.resolve_group_colours(
      ggmethylation:::.GROUP_PALETTE_DEFAULT, c("REF", "ALT")
    )
  )
})

test_that("a user-supplied palette that does not fit warns before adapting", {
  expect_warning(
    res <- ggmethylation:::.resolve_group_colours(
      c("1" = "red", "2" = "blue"), c("REF", "ALT")
    ),
    "do not cover the groups present"
  )
  expect_identical(sort(names(res)), c("ALT", "REF"))
})

test_that("an unnamed palette is assigned positionally without warning", {
  expect_silent(
    res <- ggmethylation:::.resolve_group_colours(
      c("red", "blue"), c("REF", "ALT")
    )
  )
  expect_identical(res[["ALT"]], "red")
  expect_identical(res[["REF"]], "blue")
})

test_that("extra groups are topped up with distinct colours, not recycled", {
  # Recycling would give two groups the same colour, which makes the figure
  # unreadable; top up from Okabe-Ito instead.
  res <- ggmethylation:::.resolve_group_colours(
    c("red", "blue"), c("a", "b", "c")
  )
  expect_length(res, 3L)
  expect_identical(res[["a"]], "red")
  expect_identical(res[["b"]], "blue")
  expect_false(anyDuplicated(res) > 0L)
  expect_true(res[["c"]] %in% ggmethylation:::.OKABE_ITO)
})

test_that("a partially-named palette keeps its match and fills the rest", {
  expect_warning(
    res <- ggmethylation:::.resolve_group_colours(
      c("1" = "darkgreen"), c("1", "B")
    ),
    "do not cover the groups present"
  )
  expect_identical(res[["1"]], "darkgreen")
  expect_false(is.na(res[["B"]]))
  expect_false(anyDuplicated(res) > 0L)
})

test_that("NULL passes through and NA-only groups are left alone", {
  expect_null(ggmethylation:::.resolve_group_colours(NULL, c("a", "b")))

  pal <- c("1" = "red")
  expect_identical(
    ggmethylation:::.resolve_group_colours(pal, c(NA_character_, NA_character_)),
    pal
  )
})

test_that("NA groups do not consume a colour", {
  res <- ggmethylation:::.resolve_group_colours(
    c("red", "blue"), c("REF", "ALT", NA_character_)
  )
  # Only the two real groups get colours; NA is handled by na.value.
  expect_length(res, 2L)
  expect_false(any(is.na(names(res))))
})

test_that("plot_methylation colours SNV groups without warning", {
  reads <- data.frame(
    read_name = c("r1", "r2"),
    start     = c(1000L, 1200L),
    end       = c(1800L, 2000L),
    strand    = c("+", "-"),
    group     = c("REF", "ALT"),
    stringsAsFactors = FALSE
  )
  sites <- data.frame(
    read_name = rep(c("r1", "r2"), each = 2L),
    position  = c(1100L, 1200L, 1300L, 1400L),
    mod_prob  = c(0.9, 0.1, 0.4, 0.6),
    mod_code  = "m",
    group     = rep(c("REF", "ALT"), each = 2L),
    stringsAsFactors = FALSE
  )
  md <- structure(
    list(
      reads = reads, sites = sites,
      insertion_sites = sites[0, , drop = FALSE], cigar_features = NULL,
      region = GenomicRanges::GRanges("chr1", IRanges::IRanges(1000, 2000)),
      mod_code = "m", group_tag = "SNV"
    ),
    class = "methylation_data"
  )

  p <- plot_methylation(md, show_supplementary = FALSE)

  # Printing is what triggers scale resolution, and therefore the old warning.
  png(tempfile(fileext = ".png"))
  on.exit(dev.off(), add = TRUE)
  expect_no_warning(print(p))
})
