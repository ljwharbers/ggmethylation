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

test_that("non-overlapping genes are packed into the same row", {
  # Two transcripts that do not overlap — packing should assign both to row 1.
  ann <- structure(
    list(
      transcripts = data.frame(
        tx_id = 1:2, tx_name = c("tx1", "tx2"), gene_name = c("GENE1", "GENE2"),
        strand = c("+", "+"), tx_start = c(1000L, 8000L), tx_end = c(3000L, 10000L),
        stringsAsFactors = FALSE
      ),
      exons = data.frame(
        tx_id = 1:2,
        exon_start = c(1000L, 8000L),
        exon_end   = c(3000L, 10000L),
        stringsAsFactors = FALSE
      ),
      cds  = data.frame(tx_id = integer(0), cds_start = integer(0), cds_end = integer(0)),
      utr5 = data.frame(tx_id = integer(0), utr_start = integer(0), utr_end = integer(0)),
      utr3 = data.frame(tx_id = integer(0), utr_start = integer(0), utr_end = integer(0)),
      region = GenomicRanges::GRanges("chr1", IRanges::IRanges(1000, 10000))
    ),
    class = "gene_annotations"
  )
  p <- build_gene_panel(ann, region_start = 1000L, region_end = 10000L)
  expect_s3_class(p, "gg")
  # Both genes should be packed onto the same display row (y-coordinate).
  # Absolute value depends on the internal row_spacing multiplier, so assert
  # equality rather than a specific number.
  pd <- ggplot2::ggplot_build(p)
  seg_y <- pd$data[[1L]]$y
  expect_equal(length(unique(seg_y)), 1L)
})

test_that("genes clipped to region edge still show CDS boxes", {
  # A transcript that extends beyond region_end: without pre-clipping its CDS
  # boxes would be censored, leaving only a thin intron line.
  ann <- structure(
    list(
      transcripts = data.frame(
        tx_id = 1L, tx_name = "tx1", gene_name = "GENE1",
        strand = "+", tx_start = 1000L, tx_end = 8000L,
        stringsAsFactors = FALSE
      ),
      exons = data.frame(
        tx_id = 1L, exon_start = c(1000L, 4000L), exon_end = c(2000L, 8000L),
        stringsAsFactors = FALSE
      ),
      cds  = data.frame(tx_id = 1L, cds_start = c(1200L, 4000L), cds_end = c(2000L, 8000L)),
      utr5 = data.frame(tx_id = 1L, utr_start = 1000L, utr_end = 1199L),
      utr3 = data.frame(tx_id = integer(0), utr_start = integer(0), utr_end = integer(0)),
      region = GenomicRanges::GRanges("chr1", IRanges::IRanges(1000, 5000))
    ),
    class = "gene_annotations"
  )
  # Region only shows 1000-5000 but transcript goes to 8000
  p <- build_gene_panel(ann, region_start = 1000L, region_end = 5000L)
  expect_s3_class(p, "gg")
  pd <- ggplot2::ggplot_build(p)
  # CDS rect layer: xmax values should be clamped to 5000 (region_end), not 8000
  # The CDS geom_rect is layer 3 (intron, utr, cds); find it by checking for xmax
  cds_layer <- pd$data[[3L]]
  expect_true(all(cds_layer$xmax <= 5000L))
  expect_true(any(cds_layer$xmax == 5000L))
})

test_that("gene name labels are positioned below the gene body row", {
  ann <- structure(
    list(
      transcripts = data.frame(
        tx_id = 1L, tx_name = "tx1", gene_name = "GENE1",
        strand = "+", tx_start = 1000L, tx_end = 5000L,
        stringsAsFactors = FALSE
      ),
      exons = data.frame(
        tx_id = 1L, exon_start = 1000L, exon_end = 5000L,
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
  pd <- ggplot2::ggplot_build(p)
  # Find the geom_text layer by looking for the one that has a 'label' column.
  # Gene body row = 1; label y should be below row 1 (i.e. < 1).
  text_idx   <- which(vapply(pd$data, function(d) "label" %in% names(d), logical(1L)))
  text_layer <- pd$data[[text_idx]]
  expect_true(all(text_layer$y < 1))
})

test_that("non-coding transcripts get exon boxes when region also has coding tx", {
  # Regression: when any transcript in the region has CDS, the renderer
  # previously took the IGV three-tier branch and only drew rectangles for
  # tx_ids present in `cds`/`utr*`. Non-coding transcripts (lncRNAs, many
  # LOC* entries, miRNAs) have exons but no CDS/UTR, so they rendered as
  # only a thin intron line. They must now receive uniform exon boxes.
  ann <- structure(
    list(
      transcripts = data.frame(
        tx_id = 1:2, tx_name = c("tx_coding", "tx_noncoding"),
        gene_name = c("GENE1", "LOC123"),
        strand = c("+", "+"),
        tx_start = c(1000L, 6000L), tx_end = c(5000L, 10000L),
        stringsAsFactors = FALSE
      ),
      exons = data.frame(
        tx_id      = c(1L, 1L, 2L, 2L),
        exon_start = c(1000L, 3000L, 6000L, 8500L),
        exon_end   = c(2000L, 5000L, 7500L, 10000L),
        stringsAsFactors = FALSE
      ),
      cds  = data.frame(tx_id = c(1L, 1L),
                        cds_start = c(1200L, 3000L),
                        cds_end   = c(2000L, 4800L)),
      utr5 = data.frame(tx_id = 1L, utr_start = 1000L, utr_end = 1199L),
      utr3 = data.frame(tx_id = 1L, utr_start = 4801L, utr_end = 5000L),
      region = GenomicRanges::GRanges("chr1", IRanges::IRanges(1000, 10000))
    ),
    class = "gene_annotations"
  )
  p <- build_gene_panel(ann, region_start = 1000L, region_end = 10000L)
  pd <- ggplot2::ggplot_build(p)

  # Collect all rectangle layers and look for geometry matching the
  # non-coding transcript's exon coordinates (6000-7500 and 8500-10000).
  rect_layers <- Filter(function(d) all(c("xmin", "xmax") %in% names(d)), pd$data)
  xranges <- do.call(rbind, lapply(rect_layers, function(d) {
    data.frame(xmin = d$xmin, xmax = d$xmax)
  }))
  found_exon1 <- any(xranges$xmin == 6000L & xranges$xmax == 7500L)
  found_exon2 <- any(xranges$xmin == 8500L & xranges$xmax == 10000L)
  expect_true(found_exon1)
  expect_true(found_exon2)
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
