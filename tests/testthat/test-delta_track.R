test_that("compute_group_delta returns NULL and messages for non-2 groups", {
  sites <- data.frame(
    position = 1:5 * 100, mod_prob = seq(0.1, 0.5, length.out = 5),
    group = "A", stringsAsFactors = FALSE
  )
  expect_message(
    res <- ggmethylation:::.compute_group_delta(sites, "group", span = NULL),
    "exactly two groups"
  )
  expect_null(res)
})

test_that("compute_group_delta returns signed delta on shared grid for 2 groups", {
  sites <- data.frame(
    position = rep(1:6 * 100, 2),
    mod_prob = c(rep(0.2, 6), rep(0.8, 6)),
    group    = rep(c("1", "2"), each = 6),
    stringsAsFactors = FALSE
  )
  res <- ggmethylation:::.compute_group_delta(sites, "group", span = NULL)
  expect_true(all(c("position", "delta", "sign") %in% names(res)))
  # group2 (0.8) - group1 (0.2) ~ +0.6 where both supported
  ok <- !is.na(res$delta)
  expect_true(mean(res$delta[ok]) > 0)
  expect_true(all(res$sign[ok & res$delta > 0] == "pos"))
})

test_that("compute_group_delta handles both groups with < 4 unique positions (approx fallback)", {
  # Both groups have only 3 unique positions -> smooth_methylation() falls
  # back to stats::approx(rule = 1) for each group on the shared grid.
  sites <- data.frame(
    position = c(100, 200, 300, 150, 250, 350),
    mod_prob = c(0.1, 0.2, 0.3, 0.9, 0.8, 0.7),
    group    = c("1", "1", "1", "2", "2", "2"),
    stringsAsFactors = FALSE
  )
  res <- ggmethylation:::.compute_group_delta(sites, "group", span = NULL, n_grid = 5L)

  # grid = seq(100, 350, length.out = 5) = 100, 162.5, 225, 287.5, 350
  expect_equal(res$position, c(100, 162.5, 225, 287.5, 350))

  # Hand-computed via stats::approx(rule = 1):
  # group "1": approx(c(100,200,300), c(0.1,0.2,0.3), xout = grid)
  #   -> 0.1, 0.1625, 0.225, 0.2875, NA (350 is outside range)
  # group "2": approx(c(150,250,350), c(0.9,0.8,0.7), xout = grid)
  #   -> NA (100 is outside range), 0.8875, 0.825, 0.7625, 0.7
  # delta = group2 - group1
  expected_delta <- c(NA, 0.8875 - 0.1625, 0.825 - 0.225, 0.7625 - 0.2875, NA)
  expect_equal(res$delta, expected_delta, tolerance = 1e-8)

  # Outside the overlapping support of both groups: NA delta and NA sign
  expect_true(is.na(res$delta[res$position == 100]))
  expect_true(is.na(res$delta[res$position == 350]))
  expect_true(is.na(res$sign[res$position == 100]))
  expect_true(is.na(res$sign[res$position == 350]))

  # Within overlapping support: non-NA, positive (group 2 > group 1 everywhere)
  inside <- res$position %in% c(162.5, 225, 287.5)
  expect_false(any(is.na(res$delta[inside])))
  expect_true(all(res$sign[inside] == "pos"))
})

test_that("build_delta_panel returns a ggplot", {
  df <- data.frame(position = seq(1000, 2000, length.out = 50),
                   delta = sin(seq(0, 3, length.out = 50)) * 0.4,
                   stringsAsFactors = FALSE)
  df$sign <- ifelse(df$delta >= 0, "pos", "neg")
  p <- ggmethylation:::.build_delta_panel(df, 1000, 2000)
  expect_s3_class(p, "ggplot")
})

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

test_that("plot_methylation adds delta panel for 2 groups when show_delta", {
  reads <- data.frame(
    read_name = c("r1","r2"), start = 1000L, end = 2000L, strand = "+",
    lane = 0L, mean_mod_prob = 0.5, group = c("1","2"),
    clip_side = NA_character_, sa_chrom = NA_character_, stringsAsFactors = FALSE
  )
  sites <- data.frame(
    read_name = rep(c("r1","r2"), each = 6),
    position = rep(seq(1100, 1900, length.out = 6), 2),
    mod_prob = c(rep(0.2,6), rep(0.8,6)), mod_code = "m",
    group = rep(c("1","2"), each = 6), stringsAsFactors = FALSE
  )
  md <- make_test_data(reads, sites); md$group_tag <- "HP"
  p_no  <- ggmethylation::plot_methylation(md, show_delta = FALSE, show_supplementary = FALSE)
  p_yes <- ggmethylation::plot_methylation(md, show_delta = TRUE,  show_supplementary = FALSE)
  expect_gt(length(p_yes$patches$plots), length(p_no$patches$plots))
})
