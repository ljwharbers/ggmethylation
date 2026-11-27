# Tests for read_annotations() in R/read_annotations.R

# Skip all tests if GenomicFeatures is not available
skip_if_not_installed("GenomicFeatures")

# Build a minimal synthetic TxDb once for all tests in this file.
# makeTxDb() has moved to txdbmaker (deprecated in GenomicFeatures). The
# "without CDS information" and "chromosome lengths" warnings are expected
# for this deliberately minimal TxDb and are muffled.
make_txdb = if (requireNamespace("txdbmaker", quietly = TRUE)) {
  txdbmaker::makeTxDb
} else {
  getExportedValue("GenomicFeatures", "makeTxDb")
}
txdb = withCallingHandlers(make_txdb(
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
), warning = function(w) {
  if (grepl("without CDS information|chromosome lengths", conditionMessage(w))) {
    invokeRestart("muffleWarning")
  }
})

# --- Test 1: overlapping region returns correct transcript and exon counts ---

test_that("read_annotations returns 1 transcript and 1 exon for overlapping region", {
  ann = read_annotations(txdb = txdb, region = "chr1:900-2100")

  expect_s3_class(ann, "gene_annotations")
  expect_named(ann, c("transcripts", "exons", "cds", "utr5", "utr3", "region"))
  expect_equal(nrow(ann$transcripts), 1L)
  expect_equal(nrow(ann$exons), 1L)
  expect_equal(ann$transcripts$tx_name, "tx1")
})

# --- Test: $cds, $utr5, $utr3 fields exist in returned object ---

test_that("read_annotations includes cds, utr5, utr3 fields in returned object", {
  ann = read_annotations(txdb = txdb, region = "chr1:900-2100")

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
  ann = read_annotations(txdb = txdb, region = "chr1:5000-6000")
  expect_equal(nrow(ann$transcripts), 0L)
  expect_equal(nrow(ann$exons), 0L)

  # A chromosome the TxDb does not know warns (instead of leaking
  # GenomicRanges' "no sequence levels in common" warning)
  expect_warning(
    ann <- read_annotations(txdb = txdb, region = "chr2:1-100"),
    "Chromosome 'chr2' of `region` is not among the annotation's sequences \\(chr1\\)"
  )

  expect_s3_class(ann, "gene_annotations")
  expect_equal(nrow(ann$transcripts), 0L)
  expect_equal(nrow(ann$exons), 0L)
  # Empty result should still contain cds/utr fields
  expect_true("cds"  %in% names(ann))
  expect_true("utr5" %in% names(ann))
  expect_true("utr3" %in% names(ann))
})

test_that("read_annotations warns on a chr-prefix naming mismatch", {
  expect_warning(
    ann <- read_annotations(txdb = txdb, region = "1:900-2100"),
    "'chr1' vs '1'"
  )
  expect_s3_class(ann, "gene_annotations")
  expect_equal(nrow(ann$transcripts), 0L)
  expect_named(ann, c("transcripts", "exons", "cds", "utr5", "utr3", "region"))
})

test_that("read_annotations builds a TxDb from a GTF without deprecation warnings", {
  skip_if_not_installed("rtracklayer")
  skip_if_not_installed("txdbmaker")
  gtf = tempfile(fileext = ".gtf")
  attrs = 'gene_id "G1"; transcript_id "T1"; gene_name "GENEA";'
  writeLines(c(
    paste("chr1", "test", "transcript", 1000, 2000, ".", "+", ".", attrs, sep = "\t"),
    paste("chr1", "test", "exon",       1000, 1200, ".", "+", ".", attrs, sep = "\t"),
    paste("chr1", "test", "exon",       1800, 2000, ".", "+", ".", attrs, sep = "\t")
  ), gtf)
  on.exit(clear_annotation_cache(), add = TRUE)
  msgs = character(0)
  ann = withCallingHandlers(
    suppressMessages(read_annotations(gtf = gtf, region = "chr1:900-2100")),
    warning = function(w) {
      msgs <<- c(msgs, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_false(any(grepl("has moved to the txdbmaker package", msgs)))
  expect_equal(ann$transcripts$gene_name, "GENEA")
  expect_equal(nrow(ann$exons), 2L)
})

# A TxDb with gene GENE2 carrying two transcripts of different span
# (tx2: 3000-3400, tx3: 3000-3800) next to the one-transcript GENE1.
txdb2 = withCallingHandlers(make_txdb(
  transcripts = data.frame(
    tx_id     = 1:3,
    tx_name   = c("tx1", "tx2", "tx3"),
    tx_chrom  = "chr1",
    tx_strand = "+",
    tx_start  = c(1000L, 3000L, 3000L),
    tx_end    = c(2000L, 3400L, 3800L)
  ),
  splicings = data.frame(
    tx_id      = c(1L, 2L, 3L, 3L),
    exon_rank  = c(1L, 1L, 1L, 2L),
    exon_start = c(1000L, 3000L, 3000L, 3700L),
    exon_end   = c(1200L, 3400L, 3200L, 3800L)
  ),
  genes = data.frame(
    tx_id   = 1:3,
    gene_id = c("GENE1", "GENE2", "GENE2")
  )
), warning = function(w) {
  if (grepl("without CDS information|chromosome lengths", conditionMessage(w))) {
    invokeRestart("muffleWarning")
  }
})

test_that("collapse_transcripts keeps one transcript per gene, or all of them", {
  skip_if_not_installed("AnnotationDbi")
  ann = read_annotations(txdb = txdb2, region = "chr1:900-4000")
  expect_identical(ann$transcripts$tx_name, c("tx1", "tx3"))  # longest span
  expect_identical(ann$transcripts$gene_name, c("GENE1", "GENE2"))
  expect_setequal(ann$exons$tx_id, c(1L, 3L))

  ann_all = read_annotations(txdb = txdb2, region = "chr1:900-4000",
                             collapse_transcripts = FALSE)
  expect_setequal(ann_all$transcripts$tx_name, c("tx1", "tx2", "tx3"))
  expect_identical(sort(ann_all$transcripts$gene_name),
                   c("GENE1", "GENE2", "GENE2"))
  expect_setequal(ann_all$exons$tx_id, 1:3)
})

test_that("collapse_transcripts is validated", {
  for (bad in list(NA, "yes", c(TRUE, FALSE), NULL)) {
    expect_error(
      read_annotations(txdb = txdb2, region = "chr1:900-4000",
                       collapse_transcripts = bad),
      "`collapse_transcripts` must be TRUE or FALSE"
    )
  }
})

test_that("without AnnotationDbi, transcripts are named after themselves", {
  local_mocked_bindings(.has_annotationdbi = function() FALSE)
  ann = read_annotations(txdb = txdb2, region = "chr1:900-4000")
  # Every transcript is its own gene, so collapsing keeps all three
  expect_setequal(ann$transcripts$tx_name, c("tx1", "tx2", "tx3"))
  expect_identical(ann$transcripts$gene_name, ann$transcripts$tx_name)
})

test_that("a GTF without transcript rows still gives gene names and collapses", {
  skip_if_not_installed("rtracklayer")
  skip_if_not_installed("txdbmaker")
  gtf = tempfile(fileext = ".gtf")
  a1 = 'gene_id "G1"; transcript_id "T1"; gene_name "GENEA";'
  a2 = 'gene_id "G1"; transcript_id "T2"; gene_name "GENEA";'
  a3 = 'gene_id "G2"; transcript_id "T3";'
  # Exon rows only, as in UCSC Table Browser GTF exports
  writeLines(c(
    paste("chr1", "test", "exon", 1000, 1200, ".", "+", ".", a1, sep = "\t"),
    paste("chr1", "test", "exon", 1800, 2000, ".", "+", ".", a1, sep = "\t"),
    paste("chr1", "test", "exon", 1000, 1200, ".", "+", ".", a2, sep = "\t"),
    paste("chr1", "test", "exon", 1500, 1600, ".", "+", ".", a2, sep = "\t"),
    paste("chr1", "test", "exon", 3000, 3500, ".", "-", ".", a3, sep = "\t")
  ), gtf)
  on.exit(clear_annotation_cache(), add = TRUE)

  ann = suppressMessages(read_annotations(gtf = gtf, region = "chr1:900-4000"))
  expect_identical(ann$transcripts$tx_name, c("T1", "T3"))
  expect_identical(ann$transcripts$gene_name, c("GENEA", "G2"))  # gene_id fallback

  ann_all = suppressMessages(read_annotations(gtf = gtf, region = "chr1:900-4000",
                                              collapse_transcripts = FALSE))
  expect_setequal(ann_all$transcripts$tx_name, c("T1", "T2", "T3"))
  expect_false(anyNA(ann_all$transcripts$gene_name))
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
  ann = read_annotations(txdb = txdb, region = "chr1:900-2100")

  expect_no_error(print(ann))
  expect_invisible(print(ann))
})

test_that("print.gene_annotations works for empty annotations", {
  ann = read_annotations(txdb = txdb, region = "chr1:5000-6000")

  expect_no_error(print(ann))
  expect_invisible(print(ann))
})

# --- Test 6: clear_annotation_cache() empties the cache ---

test_that("clear_annotation_cache() returns NULL invisibly and leaves cache empty", {
  # Seed the cache with a dummy entry
  assign("dummy_key", txdb, envir = ggmethylation:::.annotation_cache)

  result = clear_annotation_cache()

  expect_null(result)
  expect_equal(length(ls(ggmethylation:::.annotation_cache)), 0L)
})
