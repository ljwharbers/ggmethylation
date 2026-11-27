# Every label the package draws must work on the default pdf() device, which
# cannot draw characters outside its single-byte encoding (a "conversion
# failure ... in 'mbcsToSbcs'" error or warning). ggsave("x.pdf") uses it.

on_pdf = function(p) {
  f = tempfile(fileext = ".pdf")
  grDevices::pdf(f)
  on.exit({ grDevices::dev.off(); unlink(f) })
  print(p)
  invisible(TRUE)
}

test_that("the delta panel draws on pdf(), in continuous and binary mode", {
  md = suppressWarnings(read_methylation(fixture_bam(), "chr1:1000-7000",
                                         group_tag = "HP"))
  expect_no_warning(on_pdf(plot_methylation(md, show_delta = TRUE)))
  expect_no_warning(on_pdf(plot_methylation(md, show_delta = TRUE,
                                            call_threshold = 0.5)))
})

test_that("the variant overlay draws on pdf(), with an SV cut at the region edge", {
  skip_if_not_installed("VariantAnnotation")
  skip_if_not_installed("SummarizedExperiment")
  vars = read_variants(fixture_vcf(), "chr1:500-7000")
  expect_true(all(c("SNV", "DEL", "BND") %in% vars$variants$type))
  for (region in c("chr1:1000-7000", "chr1:3950-4050")) {
    md = suppressWarnings(read_methylation(fixture_bam(), region))
    expect_no_warning(on_pdf(plot_methylation(md, variants = vars)))
  }
})

test_that("the gene track draws on pdf()", {
  md = suppressWarnings(read_methylation(fixture_bam(), "chr1:1000-7000"))
  no_rows = function(prefix) {
    out = data.frame(tx_id = integer(0), start = integer(0), end = integer(0))
    names(out)[2:3] = paste0(prefix, c("_start", "_end"))
    out
  }
  ann = structure(
    list(
      transcripts = data.frame(
        tx_id = 1:2, tx_name = c("tx1", "tx2"), gene_name = c("GENE1", "GENE2"),
        strand = c("+", "-"), tx_start = c(1500L, 5200L), tx_end = c(5000L, 6800L),
        stringsAsFactors = FALSE
      ),
      exons = data.frame(
        tx_id = c(1L, 1L, 2L), exon_start = c(1500L, 3000L, 5200L),
        exon_end = c(2000L, 5000L, 6800L), stringsAsFactors = FALSE
      ),
      cds = no_rows("cds"), utr5 = no_rows("utr"), utr3 = no_rows("utr"),
      region = GenomicRanges::GRanges("chr1", IRanges::IRanges(1000L, 7000L))
    ),
    class = "gene_annotations"
  )
  expect_no_warning(on_pdf(plot_methylation(md, annotations = ann)))
})

test_that("the insertion locus plot draws on pdf()", {
  md = suppressWarnings(read_methylation(fixture_bam(), "chr1:1000-7000",
                                         mod_code = c("m", "h")))
  loci = list_insertion_loci(md, min_length = 50L)
  expect_gt(nrow(loci), 0L)
  expect_no_warning(on_pdf(plot_insertion_locus(md, loci[1L, ])))
})
