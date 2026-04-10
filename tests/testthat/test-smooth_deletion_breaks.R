# Tests for .consensus_deletion_ranges() and .insert_deletion_breaks()
# in R/plot_methylation.R

# --- Helpers to build minimal test data ---

make_cigar_features <- function(read_names, starts, ends, lengths = NULL) {
  n <- length(read_names)
  if (is.null(lengths)) lengths <- ends - starts + 1L
  data.frame(
    type        = rep("D", n),
    ref_start   = as.integer(starts),
    ref_end     = as.integer(ends),
    query_start = NA_integer_,
    query_end   = NA_integer_,
    length      = as.integer(lengths),
    read_name   = read_names,
    stringsAsFactors = FALSE
  )
}

make_reads <- function(read_names, group_values, group_col = "group") {
  df <- data.frame(
    read_name = read_names,
    start     = rep(100L, length(read_names)),
    end       = rep(1000L, length(read_names)),
    strand    = rep("+", length(read_names)),
    stringsAsFactors = FALSE
  )
  df[[group_col]] <- group_values
  df
}

make_smoothed <- function(positions, mean_probs, group_col = "group",
                          group_val = "all", extra_cols = list()) {
  df <- data.frame(
    position  = positions,
    mean_prob = mean_probs,
    stringsAsFactors = FALSE
  )
  df[[group_col]] <- group_val
  for (nm in names(extra_cols)) df[[nm]] <- extra_cols[[nm]]
  df
}

# =============================================================================
# .consensus_deletion_ranges() tests
# =============================================================================

test_that(".consensus_deletion_ranges returns empty df when no deletions", {
  cf <- data.frame(
    type = character(0L), ref_start = integer(0L), ref_end = integer(0L),
    query_start = integer(0L), query_end = integer(0L),
    length = integer(0L), read_name = character(0L),
    stringsAsFactors = FALSE
  )
  reads <- make_reads(c("r1", "r2"), c("all", "all"))
  result <- ggmethylation:::.consensus_deletion_ranges(cf, reads, "group")
  expect_equal(nrow(result), 0L)
  expect_true(all(c("del_start", "del_end", "group") %in% names(result)))
})

test_that(".consensus_deletion_ranges keeps deletion at exactly the threshold", {
  # 3/4 = 0.75 >= 0.75 -> keep
  cf <- make_cigar_features(c("r1", "r2", "r3"), c(500, 500, 500), c(600, 600, 600))
  reads <- make_reads(c("r1", "r2", "r3", "r4"), rep("all", 4))
  result <- ggmethylation:::.consensus_deletion_ranges(cf, reads, "group")
  expect_equal(nrow(result), 1L)
  expect_equal(result$del_start, 500L)
  expect_equal(result$del_end, 600L)
})

test_that(".consensus_deletion_ranges drops deletion below threshold", {
  # 2/4 = 0.50 < 0.75 -> drop
  cf <- make_cigar_features(c("r1", "r2"), c(500, 500), c(600, 600))
  reads <- make_reads(c("r1", "r2", "r3", "r4"), rep("all", 4))
  result <- ggmethylation:::.consensus_deletion_ranges(cf, reads, "group")
  expect_equal(nrow(result), 0L)
})

test_that(".consensus_deletion_ranges handles multiple groups correctly", {
  # Group "1": 3/3 = 100% - keep; Group "2": 1/3 = 33% - drop
  cf <- make_cigar_features(
    c("r1", "r2", "r3", "r4"),
    c(500,  500,  500,  500),
    c(600,  600,  600,  600)
  )
  reads <- make_reads(
    c("r1", "r2", "r3", "r4", "r5", "r6"),
    c("1",  "1",  "1",  "2",  "2",  "2")
  )
  result <- ggmethylation:::.consensus_deletion_ranges(cf, reads, "group")
  expect_equal(nrow(result), 1L)
  expect_equal(result$group, "1")
})

test_that(".consensus_deletion_ranges ignores NA-group deletions", {
  cf <- make_cigar_features(
    c("r1", "r2", "r3"),
    c(500, 520, 700),
    c(600, 620, 800)
  )
  reads <- make_reads(
    c("r1", "r2", "r3"),
    c("1", NA, "2")
  )

  result <- ggmethylation:::.consensus_deletion_ranges(cf, reads, "group")

  expect_false(any(is.na(result$group)))
  expect_setequal(result$group, c("1", "2"))
})

test_that(".consensus_deletion_ranges merges overlapping deletions", {
  # Two overlapping deletions from different reads should merge
  cf <- make_cigar_features(
    c("r1",  "r2",  "r3"),
    c(500,   580,   500),
    c(620,   700,   700)
  )
  reads <- make_reads(c("r1", "r2", "r3"), rep("all", 3))
  result <- ggmethylation:::.consensus_deletion_ranges(cf, reads, "group")
  # All 3 reads out of 3 carry deletions covering [500, 700] -> one merged range
  expect_equal(nrow(result), 1L)
  expect_equal(result$del_start, 500L)
  expect_equal(result$del_end, 700L)
})

test_that(".consensus_deletion_ranges keeps non-overlapping deletions separate", {
  # Two non-overlapping intervals, each supported by 3/3 reads
  cf <- make_cigar_features(
    c("r1",  "r2",  "r3",  "r1",  "r2",  "r3"),
    c(500,   500,   500,   800,   800,   800),
    c(600,   600,   600,   900,   900,   900)
  )
  reads <- make_reads(c("r1", "r2", "r3"), rep("all", 3))
  result <- ggmethylation:::.consensus_deletion_ranges(cf, reads, "group")
  expect_equal(nrow(result), 2L)
  expect_equal(sort(result$del_start), c(500L, 800L))
})

test_that(".consensus_deletion_ranges handles single read in group", {
  # 1/1 = 100% >= 75% -> keep
  cf <- make_cigar_features("r1", 500, 600)
  reads <- make_reads("r1", "all")
  result <- ggmethylation:::.consensus_deletion_ranges(cf, reads, "group")
  expect_equal(nrow(result), 1L)
})

# =============================================================================
# .insert_deletion_breaks() tests
# =============================================================================

test_that(".insert_deletion_breaks returns smoothed unchanged when no ranges", {
  sm <- make_smoothed(1:10, rep(0.5, 10))
  ranges <- data.frame(group = character(0L), del_start = integer(0L),
                       del_end = integer(0L), stringsAsFactors = FALSE)
  result <- ggmethylation:::.insert_deletion_breaks(sm, ranges, "group")
  expect_equal(result$position,  sm$position)
  expect_equal(result$mean_prob, sm$mean_prob)
})

test_that(".insert_deletion_breaks sets grid points within deletion to NA", {
  positions <- c(100, 200, 300, 400, 500, 600, 700, 800, 900, 1000)
  sm <- make_smoothed(positions, rep(0.5, 10))
  ranges <- data.frame(group = "all", del_start = 300L, del_end = 500L,
                       stringsAsFactors = FALSE)
  result <- ggmethylation:::.insert_deletion_breaks(sm, ranges, "group")
  in_del <- result$position >= 300 & result$position <= 500 &
    !is.na(result$position)
  expect_true(all(is.na(result$mean_prob[in_del])))
})

test_that(".insert_deletion_breaks inserts sentinel rows at half-integer positions", {
  positions <- c(100, 200, 400, 500, 600)
  sm <- make_smoothed(positions, rep(0.5, 5))
  ranges <- data.frame(group = "all", del_start = 300L, del_end = 350L,
                       stringsAsFactors = FALSE)
  result <- ggmethylation:::.insert_deletion_breaks(sm, ranges, "group")
  expect_true(299.5 %in% result$position)
  expect_true(350.5 %in% result$position)
  expect_true(is.na(result$mean_prob[result$position == 299.5]))
  expect_true(is.na(result$mean_prob[result$position == 350.5]))
})

test_that(".insert_deletion_breaks preserves sort order", {
  positions <- c(100, 300, 500, 700, 900)
  sm <- make_smoothed(positions, rep(0.5, 5))
  ranges <- data.frame(group = "all", del_start = 400L, del_end = 450L,
                       stringsAsFactors = FALSE)
  result <- ggmethylation:::.insert_deletion_breaks(sm, ranges, "group")
  expect_equal(result$position, sort(result$position))
})

test_that(".insert_deletion_breaks only affects the matching group", {
  sm <- rbind(
    make_smoothed(c(300, 500, 700), c(0.5, 0.5, 0.5), group_val = "A"),
    make_smoothed(c(300, 500, 700), c(0.5, 0.5, 0.5), group_val = "B")
  )
  ranges <- data.frame(group = "A", del_start = 400L, del_end = 600L,
                       stringsAsFactors = FALSE)
  result <- ggmethylation:::.insert_deletion_breaks(sm, ranges, "group")
  b_rows <- result[result$group == "B", ]
  expect_true(all(!is.na(b_rows$mean_prob)))
})

test_that(".insert_deletion_breaks creates sentinels per mod_code in multi-code path", {
  # smoothed has group + mod_code columns
  sm <- data.frame(
    position  = c(100, 300, 500, 700, 100, 300, 500, 700),
    mean_prob = rep(0.5, 8),
    group     = rep("all", 8),
    mod_code  = c(rep("m", 4), rep("h", 4)),
    stringsAsFactors = FALSE
  )
  ranges <- data.frame(group = "all", del_start = 250L, del_end = 350L,
                       stringsAsFactors = FALSE)
  result <- ggmethylation:::.insert_deletion_breaks(sm, ranges, "group")
  # Sentinel at 249.5 should appear for both mod_codes
  sents <- result[result$position == 249.5, ]
  expect_equal(nrow(sents), 2L)
  expect_setequal(sents$mod_code, c("m", "h"))
})
