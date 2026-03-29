# Tests for smooth_methylation in R/smooth_methylation.R
# Signature: smooth_methylation(sites, group_col = "group", span = 0.3)
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
