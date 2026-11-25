# The gene track must use the same x coordinate system as the read, smooth and
# delta panels (coord_cartesian(xlim = region) with default expansion), so that
# gene features line up with the reads below them.

# Hand-built annotation with one two-exon gene that starts left of the region.
.ann_overhang = function(region_start = 1000L, region_end = 7000L) {
  structure(
    list(
      transcripts = data.frame(
        tx_id = 1L, tx_name = "tx1", gene_name = "GENE1",
        strand = "+", tx_start = 500L, tx_end = 5000L,
        stringsAsFactors = FALSE
      ),
      exons = data.frame(
        tx_id = 1L, exon_start = c(500L, 3000L), exon_end = c(2000L, 5000L),
        stringsAsFactors = FALSE
      ),
      cds  = data.frame(tx_id = integer(0), cds_start = integer(0), cds_end = integer(0)),
      utr5 = data.frame(tx_id = integer(0), utr_start = integer(0), utr_end = integer(0)),
      utr3 = data.frame(tx_id = integer(0), utr_start = integer(0), utr_end = integer(0)),
      region = GenomicRanges::GRanges("chr1", IRanges::IRanges(region_start, region_end))
    ),
    class = "gene_annotations"
  )
}

.x_range = function(p) ggplot2::ggplot_build(p)$layout$panel_params[[1L]]$x.range

test_that("gene panel shares its x range with the read and smooth panels", {
  bam = testthat::test_path("fixtures", "synthetic.bam")
  md  = read_methylation(bam, "chr1:1000-7000")
  p   = plot_methylation(md, annotations = .ann_overhang())

  # patchwork: the last panel is the plot itself, the others are in $patches
  x_ranges = c(lapply(p$patches$plots, .x_range), list(.x_range(p)))
  expect_length(x_ranges, 3L)  # gene, read, smooth
  # coord_cartesian(xlim = c(1000, 7000)) with the default 5% expansion:
  # 1000 - 300 = 700 and 7000 + 300 = 7300
  for (xr in x_ranges) expect_equal(xr, c(700, 7300))
})

test_that("gene panel clips (not drops) features overhanging the region", {
  p  = build_gene_panel(.ann_overhang(), region_start = 1000L, region_end = 7000L)
  pd = ggplot2::ggplot_build(p)
  rects = do.call(rbind, lapply(
    Filter(function(d) all(c("xmin", "xmax") %in% names(d)), pd$data),
    function(d) data.frame(xmin = d$xmin, xmax = d$xmax)
  ))
  # exon 500-2000 is drawn clipped to 1000-2000; exon 3000-5000 unchanged
  expect_true(any(rects$xmin == 1000 & rects$xmax == 2000))
  expect_true(any(rects$xmin == 3000 & rects$xmax == 5000))
  expect_equal(.x_range(p), c(700, 7300))
})

test_that("strand arrows near the region edge are drawn and kept inside it", {
  # Exons 4000-6699 and 6990-9000 (clipped to 6990-7000): the intron is
  # 6700..6989, its arrow starts at the midpoint 6844.5 and is
  # max(200, 289 * 0.05) = 200 bp long, i.e. would end at 7044.5. It must be
  # drawn, with its end clamped to the region end 7000.
  ann = .ann_overhang()
  ann$transcripts$tx_start = 4000L
  ann$transcripts$tx_end   = 9000L
  ann$exons = data.frame(tx_id = 1L, exon_start = c(4000L, 6990L),
                         exon_end = c(6699L, 9000L))
  p  = build_gene_panel(ann, region_start = 1000L, region_end = 7000L)
  pd = suppressWarnings(ggplot2::ggplot_build(p))
  segs = Filter(function(d) all(c("x", "xend") %in% names(d)), pd$data)
  arrows = do.call(rbind, lapply(segs, function(d) d[, c("x", "xend")]))
  arrows = arrows[arrows$x == 6844.5, , drop = FALSE]
  expect_equal(nrow(arrows), 1L)
  expect_equal(arrows$xend, 7000)
})

test_that("empty gene panel uses the same x range as the other panels", {
  ann = .ann_overhang()
  ann$transcripts = ann$transcripts[0, ]
  ann$exons       = ann$exons[0, ]
  p = build_gene_panel(ann, region_start = 1000L, region_end = 7000L)
  expect_equal(.x_range(p), c(700, 7300))
})
