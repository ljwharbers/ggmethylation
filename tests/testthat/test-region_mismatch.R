# Annotations / variants read for a different region than the plotted data
# would silently give an empty track; build_gene_panel(),
# build_variant_overlay() and plot_methylation() (for the annotations'
# chromosome) warn instead.

.vd_for = function(chrom, start, end) {
  structure(
    list(
      variants = data.frame(
        position = integer(0), ref = character(0), alt = character(0),
        type = character(0), end = integer(0), mate_chrom = character(0),
        mate_pos = integer(0), stringsAsFactors = FALSE
      ),
      region = GenomicRanges::GRanges(chrom, IRanges::IRanges(start, end)),
      bnd_match_tol = 50L
    ),
    class = "variant_data"
  )
}

.md_chr1 = structure(
  list(
    reads = data.frame(
      read_name = character(0), start = integer(0), end = integer(0),
      strand = character(0), lane = integer(0), stringsAsFactors = FALSE
    ),
    sequences = NULL, cigars = NULL,
    region = GenomicRanges::GRanges("chr1", IRanges::IRanges(1000, 7000))
  ),
  class = "methylation_data"
)

test_that("build_variant_overlay warns when variants are for another chromosome", {
  expect_warning(
    ggmethylation:::build_variant_overlay(.md_chr1, .vd_for("chr2", 1000, 7000)),
    "chr2:1000-7000, which does not overlap the plotted region chr1:1000-7000"
  )
})

test_that("build_variant_overlay warns when variants do not overlap the region", {
  expect_warning(
    ggmethylation:::build_variant_overlay(.md_chr1, .vd_for("chr1", 8000, 9000)),
    "does not overlap"
  )
})

test_that("build_variant_overlay is silent for an overlapping or wider region", {
  expect_silent(ggmethylation:::build_variant_overlay(.md_chr1, .vd_for("chr1", 1000, 7000)))
  expect_silent(ggmethylation:::build_variant_overlay(.md_chr1, .vd_for("chr1", 1, 100000)))
})

test_that("build_variant_overlay ignores NA-typed and MNV rows", {
  # NA types (missing alleles) must not be turned into all-NA SNV / BND rows,
  # which previously triggered the 'requires sequences' warning and an empty
  # BND label layer. MNVs are not drawn.
  vd = .vd_for("chr1", 1000, 7000)
  vd$variants = data.frame(
    position = c(2000L, 3000L), ref = c(NA, "GT"), alt = c(NA, "AC"),
    type = c(NA, "MNV"), end = c(2000L, 3000L),
    mate_chrom = c(NA_character_, NA_character_),
    mate_pos = c(NA_integer_, NA_integer_), stringsAsFactors = FALSE
  )
  expect_silent(ov <- ggmethylation:::build_variant_overlay(.md_chr1, vd))
  expect_null(ov$snv)
  expect_null(ov$sv)
  expect_null(ov$bnd)
})

.ann_at = function(start, end) {
  structure(
    list(
      transcripts = data.frame(
        tx_id = 1L, tx_name = "tx1", gene_name = "GENE1", strand = "+",
        tx_start = start, tx_end = end, stringsAsFactors = FALSE
      ),
      exons = data.frame(tx_id = 1L, exon_start = start, exon_end = end),
      cds = NULL, utr5 = NULL, utr3 = NULL,
      region = GenomicRanges::GRanges("chr1", IRanges::IRanges(start, end))
    ),
    class = "gene_annotations"
  )
}

test_that("build_gene_panel warns when annotations do not overlap the range", {
  expect_warning(
    build_gene_panel(.ann_at(20000L, 30000L), region_start = 1000L, region_end = 7000L),
    "chr1:20000-30000, which does not overlap the plotted range 1000-7000"
  )
  expect_silent(build_gene_panel(.ann_at(2000L, 3000L), 1000L, 7000L))
})

test_that("plot_methylation warns when annotations are on another chromosome", {
  md = suppressWarnings(read_methylation(fixture_bam(), "chr1:1000-7000"))
  ann = structure(
    list(
      transcripts = data.frame(
        tx_id = 1L, tx_name = "tx1", gene_name = "GENE1",
        strand = "+", tx_start = 1500L, tx_end = 5000L,
        stringsAsFactors = FALSE
      ),
      exons = data.frame(
        tx_id = 1L, exon_start = c(1500L, 3000L), exon_end = c(2000L, 5000L),
        stringsAsFactors = FALSE
      ),
      cds  = data.frame(tx_id = integer(0), cds_start = integer(0),
                        cds_end = integer(0)),
      utr5 = data.frame(tx_id = integer(0), utr_start = integer(0),
                        utr_end = integer(0)),
      utr3 = data.frame(tx_id = integer(0), utr_start = integer(0),
                        utr_end = integer(0)),
      region = GenomicRanges::GRanges("chr2", IRanges::IRanges(1000L, 7000L))
    ),
    class = "gene_annotations"
  )
  expect_warning(plot_methylation(md, annotations = ann), "chromosome 'chr2'")

  ann$region = GenomicRanges::GRanges("chr1", IRanges::IRanges(1000L, 7000L))
  expect_no_warning(plot_methylation(md, annotations = ann))
})
