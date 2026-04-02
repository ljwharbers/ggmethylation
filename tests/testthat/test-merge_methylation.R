# Tests for merge_methylation() in R/merge_methylation.R
# merge_methylation(..., .list = NULL) -> multi_methylation_data

# Helper: build a minimal methylation_data without reading a BAM
make_methylation_data <- function(chrom = "chr1", start = 1000L, end = 2000L,
                                  mod_code = "m", n_reads = 3L) {
  gr <- GenomicRanges::GRanges(
    seqnames = chrom,
    ranges   = IRanges::IRanges(start = start, end = end)
  )

  reads <- data.frame(
    read_name = paste0("r", seq_len(n_reads)),
    start     = rep(start, n_reads),
    end       = rep(end,   n_reads),
    strand    = rep("+",   n_reads),
    stringsAsFactors = FALSE
  )

  sites <- data.frame(
    position  = integer(0L),
    mod_prob  = numeric(0L),
    read_name = character(0L),
    mod_code  = character(0L),
    stringsAsFactors = FALSE
  )

  structure(
    list(
      reads        = reads,
      sites        = sites,
      region       = gr,
      mod_code     = mod_code,
      group_tag    = NULL,
      snv_position = NULL,
      sequences    = setNames(character(0L), character(0L)),
      cigars       = setNames(character(0L), character(0L))
    ),
    class = "methylation_data"
  )
}

# --- Test 1: merging two valid objects returns correct structure ---

test_that("merge_methylation returns multi_methylation_data with correct structure", {
  md1 <- make_methylation_data(n_reads = 3L)
  md2 <- make_methylation_data(n_reads = 5L)

  merged <- merge_methylation(sampleA = md1, sampleB = md2)

  expect_s3_class(merged, "multi_methylation_data")
  expect_named(merged, c("samples", "region", "mod_code"))
  expect_named(merged$samples, c("sampleA", "sampleB"))
  expect_equal(merged$mod_code, "m")
})

# --- Test 2: per-sample read counts are correct ---

test_that("merge_methylation preserves per-sample read counts", {
  md1 <- make_methylation_data(n_reads = 3L)
  md2 <- make_methylation_data(n_reads = 7L)

  merged <- merge_methylation(s1 = md1, s2 = md2)

  expect_equal(nrow(merged$samples$s1$reads), 3L)
  expect_equal(nrow(merged$samples$s2$reads), 7L)
})

# --- Test 3: error when samples have different regions ---

test_that("merge_methylation errors on different regions", {
  md1 <- make_methylation_data(chrom = "chr1", start = 1000L, end = 2000L)
  md2 <- make_methylation_data(chrom = "chr2", start = 1000L, end = 2000L)

  expect_error(
    merge_methylation(a = md1, b = md2),
    regexp = "same genomic region"
  )
})

test_that("merge_methylation errors when start positions differ", {
  md1 <- make_methylation_data(start = 1000L, end = 2000L)
  md2 <- make_methylation_data(start = 3000L, end = 4000L)

  expect_error(
    merge_methylation(a = md1, b = md2),
    regexp = "same genomic region"
  )
})

# --- Test 4: error when samples have different mod_code ---

test_that("merge_methylation errors on different mod_code", {
  md1 <- make_methylation_data(mod_code = "m")
  md2 <- make_methylation_data(mod_code = "h")

  expect_error(
    merge_methylation(a = md1, b = md2),
    regexp = "mod_code"
  )
})

# --- Test 5: error when samples are unnamed ---

test_that("merge_methylation errors when samples are unnamed", {
  md1 <- make_methylation_data()
  md2 <- make_methylation_data()

  expect_error(
    merge_methylation(md1, md2),
    regexp = "must be named"
  )
})

test_that("merge_methylation errors when some names are empty strings", {
  md1 <- make_methylation_data()
  md2 <- make_methylation_data()
  lst <- list(md1, md2)
  names(lst) <- c("a", "")

  expect_error(
    merge_methylation(.list = lst),
    regexp = "must be named"
  )
})

# --- Test 6: print() method runs without error ---

test_that("print.multi_methylation_data runs without error", {
  md1 <- make_methylation_data(n_reads = 2L)
  md2 <- make_methylation_data(n_reads = 4L)
  merged <- merge_methylation(x = md1, y = md2)

  expect_no_error(print(merged))
  expect_invisible(print(merged))
})
