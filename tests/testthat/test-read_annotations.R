# Tests for read_annotations() in R/read_annotations.R

# Skip all tests if GenomicFeatures is not available
skip_if_not_installed("GenomicFeatures")

# Build a minimal synthetic TxDb once for all tests in this file
txdb <- GenomicFeatures::makeTxDb(
  transcripts = data.frame(
    tx_id     = 1L,
    tx_name   = "tx1",
    tx_chrom  = "chr1",
    tx_strand = "+",
    tx_start  = 1000L,
    tx_end    = 2000L
  ),
  splicings = data.frame(
    tx_id      = 1L,
    exon_rank  = 1L,
    exon_start = 1000L,
    exon_end   = 1200L
  ),
  genes = data.frame(
    tx_id   = 1L,
    gene_id = "GENE1"
  )
)

# --- Test 1: overlapping region returns correct transcript and exon counts ---

test_that("read_annotations returns 1 transcript and 1 exon for overlapping region", {
  ann <- read_annotations(txdb = txdb, region = "chr1:900-2100")

  expect_s3_class(ann, "gene_annotations")
  expect_named(ann, c("transcripts", "exons", "cds", "utr5", "utr3", "region"))
  expect_equal(nrow(ann$transcripts), 1L)
  expect_equal(nrow(ann$exons), 1L)
  expect_equal(ann$transcripts$tx_name, "tx1")
})

# --- Test: $cds, $utr5, $utr3 fields exist in returned object ---

test_that("read_annotations includes cds, utr5, utr3 fields in returned object", {
  ann <- read_annotations(txdb = txdb, region = "chr1:900-2100")

  expect_true("cds"  %in% names(ann))
  expect_true("utr5" %in% names(ann))
  expect_true("utr3" %in% names(ann))

  # For this minimal TxDb there is no CDS/UTR data; fields should be data.frames
  expect_s3_class(ann$cds,  "data.frame")
  expect_s3_class(ann$utr5, "data.frame")
  expect_s3_class(ann$utr3, "data.frame")

  # Verify expected column names
  expect_named(ann$cds,  c("tx_id", "cds_start", "cds_end"))
  expect_named(ann$utr5, c("tx_id", "utr_start", "utr_end"))
  expect_named(ann$utr3, c("tx_id", "utr_start", "utr_end"))
})

# --- Test 2: non-overlapping region returns 0 transcripts gracefully ---

test_that("read_annotations returns 0 transcripts for non-overlapping region", {
  ann <- read_annotations(txdb = txdb, region = "chr2:1-100")

  expect_s3_class(ann, "gene_annotations")
  expect_equal(nrow(ann$transcripts), 0L)
  expect_equal(nrow(ann$exons), 0L)
  # Empty result should still contain cds/utr fields
  expect_true("cds"  %in% names(ann))
  expect_true("utr5" %in% names(ann))
  expect_true("utr3" %in% names(ann))
})

# --- Test 3: error when both txdb and gtf are provided ---

test_that("read_annotations errors when both txdb and gtf are provided", {
  expect_error(
    read_annotations(txdb = txdb, gtf = "some.gtf", region = "chr1:900-2100"),
    regexp = "exactly one"
  )
})

# --- Test 4: error when neither txdb nor gtf is provided ---

test_that("read_annotations errors when neither txdb nor gtf is provided", {
  expect_error(
    read_annotations(region = "chr1:900-2100"),
    regexp = "must be provided"
  )
})

# --- Test 5: print() method runs without error ---

test_that("print.gene_annotations runs without error", {
  ann <- read_annotations(txdb = txdb, region = "chr1:900-2100")

  expect_no_error(print(ann))
  expect_invisible(print(ann))
})

test_that("print.gene_annotations works for empty annotations", {
  ann <- read_annotations(txdb = txdb, region = "chr2:1-100")

  expect_no_error(print(ann))
  expect_invisible(print(ann))
})

# --- Test 6: clear_annotation_cache() empties the cache ---

test_that("clear_annotation_cache() returns NULL invisibly and leaves cache empty", {
  # Seed the cache with a dummy entry
  assign("dummy_key", txdb, envir = ggmethylation:::.annotation_cache)

  result <- clear_annotation_cache()

  expect_null(result)
  expect_equal(length(ls(ggmethylation:::.annotation_cache)), 0L)
})
