# Tests for smooth_methylation in R/smooth_methylation.R
# Signature: smooth_methylation(sites, group_col = "group", span = NULL)
# sites: data.frame with columns position, mod_prob, and a grouping column
# Returns: data.frame with columns position, mean_prob, and the grouping column

test_that("smooth_methylation returns empty df on NULL input", {
  result <- ggmethylation:::smooth_methylation(NULL)
  expect_equal(nrow(result), 0L)
  expect_true(all(c("position", "mean_prob", "group") %in% names(result)))
})

test_that("smooth_methylation returns empty df on 0-row input", {
  sites <- data.frame(
    position = numeric(0),
    mod_prob = numeric(0),
    group = character(0),
    stringsAsFactors = FALSE
  )
  result <- ggmethylation:::smooth_methylation(sites)
  expect_equal(nrow(result), 0L)
})

test_that("smooth_methylation returns raw means when < 4 unique positions", {
  sites <- data.frame(
    position = c(1, 1, 2, 3),
    mod_prob = c(0.8, 0.6, 0.5, 0.2),
    group = "A",
    stringsAsFactors = FALSE
  )
  result <- ggmethylation:::smooth_methylation(sites)
  # 3 unique positions -> raw means returned
  expect_equal(nrow(result), 3L)
  # mean of 0.8 and 0.6 = 0.7 at position 1
  pos1_prob <- result$mean_prob[result$position == 1]
  expect_equal(pos1_prob, 0.7)
})

test_that("smooth_methylation fits loess when >= 4 unique positions", {
  sites <- data.frame(
    position = 1:5 * 100,
    mod_prob = c(0.1, 0.2, 0.3, 0.4, 0.5),
    group = "A",
    stringsAsFactors = FALSE
  )
  result <- ggmethylation:::smooth_methylation(sites)
  # Loess produces 200-point grid
  expect_true(nrow(result) > 5L)
  # Predicted values should be clamped to [0, 1] by the loess fit on this data
  valid_probs <- result$mean_prob[!is.na(result$mean_prob)]
  expect_true(all(valid_probs >= 0 & valid_probs <= 1))
})

test_that("smooth_methylation handles multiple groups independently", {
  sites <- data.frame(
    position = rep(1:5 * 100, 2),
    mod_prob = rep(c(0.1, 0.2, 0.3, 0.4, 0.5), 2),
    group = rep(c("A", "B"), each = 5),
    stringsAsFactors = FALSE
  )
  result <- ggmethylation:::smooth_methylation(sites)
  expect_true("A" %in% result$group)
  expect_true("B" %in% result$group)
})

test_that("smooth_methylation respects custom group_col name", {
  sites <- data.frame(
    position = 1:5 * 100,
    mod_prob = c(0.1, 0.2, 0.3, 0.4, 0.5),
    haplotype = "H1",
    stringsAsFactors = FALSE
  )
  result <- ggmethylation:::smooth_methylation(sites, group_col = "haplotype")
  expect_true("haplotype" %in% names(result))
  expect_false("group" %in% names(result))
})

test_that("smooth_methylation result has expected columns", {
  sites <- data.frame(
    position = 1:5 * 100,
    mod_prob = c(0.1, 0.2, 0.3, 0.4, 0.5),
    group = "A",
    stringsAsFactors = FALSE
  )
  result <- ggmethylation:::smooth_methylation(sites)
  expect_true(all(c("position", "mean_prob", "group") %in% names(result)))
})

test_that("smooth_methylation drops NA group values", {
  sites <- data.frame(
    position = 1:6 * 100,
    mod_prob = rep(0.5, 6),
    group = c("A", "A", "A", "A", "A", NA),
    stringsAsFactors = FALSE
  )
  result <- ggmethylation:::smooth_methylation(sites)
  expect_false(any(is.na(result$group)))
})

test_that("smooth_methylation loess grid has 200 points per group", {
  sites <- data.frame(
    position = 1:5 * 100,
    mod_prob = c(0.1, 0.2, 0.3, 0.4, 0.5),
    group = "A",
    stringsAsFactors = FALSE
  )
  result <- ggmethylation:::smooth_methylation(sites)
  expect_equal(nrow(result), 200L)
})

test_that("smooth_methylation adaptive span (span = NULL) produces 200-point grid", {
  sites <- data.frame(
    position = 1:20 * 100,
    mod_prob = rep(0.5, 20),
    group = "A",
    stringsAsFactors = FALSE
  )
  result <- ggmethylation:::smooth_methylation(sites, span = NULL)
  expect_equal(nrow(result), 200L)
  expect_true(all(c("position", "mean_prob", "group") %in% names(result)))
})

test_that("smooth_methylation adaptive span scales with site count", {
  # 100 sites -> span = max(0.15, min(0.75, 15/100)) = 0.15
  sites_dense <- data.frame(
    position = 1:100 * 10,
    mod_prob = rep(0.5, 100),
    group = "A",
    stringsAsFactors = FALSE
  )
  # 10 sites -> span = max(0.15, min(0.75, 15/10)) = 0.75
  sites_sparse <- data.frame(
    position = 1:10 * 100,
    mod_prob = rep(0.5, 10),
    group = "A",
    stringsAsFactors = FALSE
  )
  # Both should complete without error and return 200-point grids
  expect_equal(nrow(ggmethylation:::smooth_methylation(sites_dense)), 200L)
  expect_equal(nrow(ggmethylation:::smooth_methylation(sites_sparse)), 200L)
})

test_that("smooth_methylation explicit numeric span overrides adaptive default", {
  # Use non-trivial data with a clear local peak so different spans produce
  # visibly different fits
  sites <- data.frame(
    position = c(100, 200, 300, 400, 500, 600, 700, 800, 900, 1000,
                 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000),
    mod_prob = c(0.1, 0.1, 0.9, 0.9, 0.9, 0.1, 0.1, 0.8, 0.8, 0.1,
                 0.1, 0.1, 0.7, 0.7, 0.7, 0.1, 0.1, 0.1, 0.1, 0.1),
    group = "A",
    stringsAsFactors = FALSE
  )
  # Explicit span = 0.2 (tight) vs adaptive span = 0.75 (20 sites, max-clamped)
  result_fixed <- ggmethylation:::smooth_methylation(sites, span = 0.2)
  result_adapt <- ggmethylation:::smooth_methylation(sites, span = NULL)
  # Both return valid 200-point grids
  expect_equal(nrow(result_fixed), 200L)
  expect_equal(nrow(result_adapt), 200L)
  # The tight span tracks local peaks; the adaptive (wide) span smooths them out
  valid_fixed <- result_fixed$mean_prob[!is.na(result_fixed$mean_prob)]
  valid_adapt <- result_adapt$mean_prob[!is.na(result_adapt$mean_prob)]
  expect_false(isTRUE(all.equal(valid_fixed, valid_adapt[seq_along(valid_fixed)])))
})

test_that("smooth_methylation returns clamped lower/upper on loess path", {
  sites <- data.frame(
    position = 1:5 * 100,
    mod_prob = c(0.1, 0.2, 0.3, 0.4, 0.5),
    group = "A",
    stringsAsFactors = FALSE
  )
  result <- ggmethylation:::smooth_methylation(sites)
  expect_true(all(c("lower", "upper") %in% names(result)))
  ok <- !is.na(result$lower)
  expect_true(all(result$lower[ok] >= 0 & result$lower[ok] <= 1))
  expect_true(all(result$upper[ok] >= 0 & result$upper[ok] <= 1))
  expect_true(all(result$upper[ok] >= result$lower[ok]))
})

test_that("smooth_methylation sets NA CI on raw-means fallback", {
  sites <- data.frame(
    position = c(1, 1, 2, 3),
    mod_prob = c(0.8, 0.6, 0.5, 0.2),
    group = "A",
    stringsAsFactors = FALSE
  )
  result <- ggmethylation:::smooth_methylation(sites)
  expect_true(all(is.na(result$lower)))
  expect_true(all(is.na(result$upper)))
})

test_that("smooth_methylation empty input includes lower/upper cols", {
  result <- ggmethylation:::smooth_methylation(NULL)
  expect_true(all(c("lower", "upper") %in% names(result)))
})

test_that("smooth_methylation grid path: >= 4 points evaluates on shared grid and masks out-of-support", {
  sites <- data.frame(
    position = 1:5 * 100,  # 100, 200, 300, 400, 500
    mod_prob = c(0.1, 0.2, 0.3, 0.4, 0.5),
    group = "A",
    stringsAsFactors = FALSE
  )
  grid <- c(0, 100, 300, 500, 600)
  result <- ggmethylation:::smooth_methylation(sites, grid = grid)
  # Output must align row-for-row with the supplied grid
  expect_equal(nrow(result), length(grid))
  expect_equal(result$position, grid)
  # Outside this group's own position range (100-500): NA in mean_prob/lower/upper
  expect_true(is.na(result$mean_prob[result$position == 0]))
  expect_true(is.na(result$mean_prob[result$position == 600]))
  expect_true(is.na(result$lower[result$position == 0]))
  expect_true(is.na(result$upper[result$position == 600]))
  # Inside the range: not NA
  inside <- result$position %in% c(100, 300, 500)
  expect_false(any(is.na(result$mean_prob[inside])))
})

test_that("smooth_methylation grid path: < 4 unique positions interpolates via approx(rule = 1)", {
  sites <- data.frame(
    position = c(1, 2, 3),
    mod_prob = c(0.8, 0.5, 0.2),
    group = "A",
    stringsAsFactors = FALSE
  )
  grid <- c(0, 1, 1.5, 2, 3, 4)
  result <- ggmethylation:::smooth_methylation(sites, grid = grid)

  expect_equal(nrow(result), length(grid))
  expect_equal(result$position, grid)

  # Hand-computed stats::approx(x = c(1,2,3), y = c(0.8,0.5,0.2), rule = 1) values
  expected <- c(NA, 0.8, 0.65, 0.5, 0.2, NA)
  expect_equal(result$mean_prob, expected, tolerance = 1e-8)

  # lower/upper are always NA on the fallback path
  expect_true(all(is.na(result$lower)))
  expect_true(all(is.na(result$upper)))
})

test_that("smooth_methylation grid path: two groups produce row-aligned output on the same grid", {
  sites <- data.frame(
    position = c(rep(1:5 * 100, 2)),
    mod_prob = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.9, 0.8, 0.7, 0.6, 0.5),
    group    = rep(c("A", "B"), each = 5),
    stringsAsFactors = FALSE
  )
  grid <- seq(100, 500, length.out = 7)
  result <- ggmethylation:::smooth_methylation(sites, grid = grid)

  pos_a <- result$position[result$group == "A"]
  pos_b <- result$position[result$group == "B"]
  expect_equal(pos_a, grid)
  expect_equal(pos_b, grid)
})
