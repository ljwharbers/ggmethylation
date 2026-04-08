# Tests for build_gene_panel() in R/build_gene_panel.R

skip_if_not_installed("ggplot2")
skip_if_not_installed("GenomicRanges")
skip_if_not_installed("IRanges")

test_that("build_gene_panel renders without CDS (fallback)", {
  # Create minimal annotations with no cds
  ann <- structure(
    list(
      transcripts = data.frame(
        tx_id = 1L, tx_name = "tx1", gene_name = "GENE1",
        strand = "+", tx_start = 1000L, tx_end = 5000L,
        stringsAsFactors = FALSE
      ),
      exons = data.frame(
        tx_id = 1L, exon_start = c(1000L, 3000L), exon_end = c(2000L, 5000L),
        stringsAsFactors = FALSE
      ),
      cds  = data.frame(tx_id = integer(0), cds_start = integer(0), cds_end = integer(0)),
      utr5 = data.frame(tx_id = integer(0), utr_start = integer(0), utr_end = integer(0)),
      utr3 = data.frame(tx_id = integer(0), utr_start = integer(0), utr_end = integer(0)),
      region = GenomicRanges::GRanges("chr1", IRanges::IRanges(1000, 5000))
    ),
    class = "gene_annotations"
  )
  p <- build_gene_panel(ann, region_start = 1000L, region_end = 5000L)
  expect_s3_class(p, "gg")
})

test_that("build_gene_panel renders with CDS (IGV style)", {
  ann <- structure(
    list(
      transcripts = data.frame(
        tx_id = 1L, tx_name = "tx1", gene_name = "GENE1",
        strand = "+", tx_start = 1000L, tx_end = 5000L,
        stringsAsFactors = FALSE
      ),
      exons = data.frame(
        tx_id = 1L, exon_start = c(1000L, 3000L), exon_end = c(2000L, 5000L),
        stringsAsFactors = FALSE
      ),
      cds  = data.frame(tx_id = 1L, cds_start = c(1200L, 3000L), cds_end = c(2000L, 4800L)),
      utr5 = data.frame(tx_id = 1L, utr_start = 1000L, utr_end = 1199L),
      utr3 = data.frame(tx_id = 1L, utr_start = 4801L, utr_end = 5000L),
      region = GenomicRanges::GRanges("chr1", IRanges::IRanges(1000, 5000))
    ),
    class = "gene_annotations"
  )
  p <- build_gene_panel(ann, region_start = 1000L, region_end = 5000L)
  expect_s3_class(p, "gg")
})

test_that("build_gene_panel handles empty annotations gracefully", {
  ann <- structure(
    list(
      transcripts = data.frame(
        tx_id = integer(0), tx_name = character(0), gene_name = character(0),
        strand = character(0), tx_start = integer(0), tx_end = integer(0),
        stringsAsFactors = FALSE
      ),
      exons = data.frame(
        tx_id = integer(0), exon_start = integer(0), exon_end = integer(0),
        stringsAsFactors = FALSE
      ),
      cds  = data.frame(tx_id = integer(0), cds_start = integer(0), cds_end = integer(0)),
      utr5 = data.frame(tx_id = integer(0), utr_start = integer(0), utr_end = integer(0)),
      utr3 = data.frame(tx_id = integer(0), utr_start = integer(0), utr_end = integer(0)),
      region = GenomicRanges::GRanges("chr1", IRanges::IRanges(1000, 5000))
    ),
    class = "gene_annotations"
  )
  p <- build_gene_panel(ann, region_start = 1000L, region_end = 5000L)
  expect_s3_class(p, "gg")
})
