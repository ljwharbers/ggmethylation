# plot_methylation() argument handling, on the synthetic fixture BAM.

fixture_md = function(...) {
  set.seed(1)
  suppressWarnings(read_methylation(fixture_bam(), "chr1:1000-7000", ...))
}

# Draw every panel on a throwaway PNG device: scale and conversion problems
# surface here, not in plot_methylation().
render = function(p) {
  f = tempfile(fileext = ".png")
  grDevices::png(f)
  on.exit({ grDevices::dev.off(); unlink(f) })
  invisible(patchwork::patchworkGrob(p))
}

layer_geoms = function(p) {
  vapply(p$layers, function(l) class(l$geom)[1L], character(1L))
}

# One row per read of the read panel's bar layer (layer 1): name, span, strand,
# group and lane, in the order the reads were packed.
packed_reads = function(rp) {
  d = rp$layers[[1L]]$data
  d = d[!startsWith(d$polygon_id, ".legend_"), , drop = FALSE]
  d[!duplicated(d$read_name), intersect(c("read_name", "start", "end", "strand",
                                          "group", "lane"), names(d))]
}

no_lane_overlaps = function(reads) {
  all(vapply(split(reads, reads$lane), function(l) {
    l = l[order(l$start), , drop = FALSE]
    nrow(l) < 2L || all(utils::tail(l$start, -1L) > utils::head(l$end, -1L))
  }, logical(1L)))
}

test_that("sort_by rejects columns that don't exist in data$reads", {
  # order() on a NULL column returns integer(0), which used to silently drop
  # every read (the vignette once passed sort_by = "position").
  md = fixture_md(group_tag = "HP")
  expect_error(plot_methylation(md, sort_by = "position"), "position")
  expect_error(
    plot_methylation(merge_methylation(A = md, B = md), sort_by = "position"),
    "position"
  )
  expect_s3_class(plot_methylation(md, sort_by = "mean_mod_prob"), "patchwork")
})

test_that("a custom sort_by changes the read order and still packs without overlaps", {
  md = fixture_md()
  default = packed_reads(plot_methylation(md, show_cigar = FALSE)[[1L]])
  expect_false(is.unsorted(default$start))

  by_strand = packed_reads(plot_methylation(md, show_cigar = FALSE,
                                            sort_by = "strand")[[1L]])
  expect_setequal(by_strand$read_name, md$reads$read_name)
  expect_false(identical(by_strand$read_name, default$read_name))
  expect_identical(by_strand$strand, sort(by_strand$strand))
  expect_true(no_lane_overlaps(by_strand))

  grouped = fixture_md(group_tag = "HP")
  default_grp = packed_reads(plot_methylation(grouped, show_cigar = FALSE)[[1L]])
  by_group = packed_reads(plot_methylation(grouped, show_cigar = FALSE,
                                           sort_by = c("group", "start"))[[1L]])
  expect_setequal(by_group$read_name, grouped$reads$read_name)
  expect_false(identical(by_group$read_name, default_grp$read_name))
  expect_identical(by_group$read_name,
                   by_group$read_name[order(by_group$group, by_group$start)])
  expect_true(no_lane_overlaps(by_group))
})

test_that("multi-sample plots mix grouped and ungrouped samples", {
  # Regression: combining the samples' site tables used to fail in rbind()
  # because only the grouped sample carried a `group` column.
  grouped   = fixture_md(group_tag = "HP")
  ungrouped = fixture_md()
  p = plot_methylation(merge_methylation(A = grouped, B = ungrouped))
  expect_s3_class(p, "patchwork")
  smooth = ggplot2::ggplot_build(p[[length(p)]])$data
  expect_length(smooth, 2L)  # CI ribbon + lines
  expect_gt(nrow(smooth[[2L]]), 0L)
})

test_that("multi-sample read panels use the same lane packing as single-sample", {
  # The multi-sample path used to skip clip-aware packing.
  md = fixture_md()
  single = plot_methylation(md)
  multi = plot_methylation(merge_methylation(A = md, B = md))
  lanes = function(p) sort(unique(ggplot2::ggplot_build(p)$data[[1L]]$y))
  expect_equal(lanes(multi[[1L]]), lanes(single[[1L]]))
})

# --- panel_heights -------------------------------------------------------------

test_that("panel_heights sets the relative panel heights", {
  md = fixture_md(group_tag = "HP")
  p = plot_methylation(md, panel_heights = c(3, 1))
  expect_identical(p$patches$layout$heights, c(3, 1))
  expect_error(plot_methylation(md, panel_heights = c(3, 1, 1)),
               "`panel_heights` has length 3 but there are 2 panels")
  expect_error(plot_methylation(md, panel_heights = 1),
               "`panel_heights` has length 1 but there are 2 panels")
})

test_that("panel_heights counts the delta panel", {
  md = fixture_md(group_tag = "HP")
  p = plot_methylation(md, show_delta = TRUE, panel_heights = c(1, 0.4, 0.3))
  expect_length(p, 3L)
  expect_identical(p$patches$layout$heights, c(1, 0.4, 0.3))
  expect_error(plot_methylation(md, show_delta = TRUE, panel_heights = c(3, 1)),
               "has length 2 but there are 3 panels")
})

test_that("panel_heights takes one read panel per sample", {
  md = fixture_md(group_tag = "HP")
  mm = merge_methylation(A = md, B = md)
  p = plot_methylation(mm, panel_heights = c(2, 2, 1))
  expect_identical(p$patches$layout$heights, c(2, 2, 1))
  expect_error(plot_methylation(mm, panel_heights = c(3, 1)),
               "has length 2 but there are 3 panels")
})

# --- variants ---------------------------------------------------------------------

test_that("variants are overlaid on the read panel", {
  skip_if_not_installed("VariantAnnotation")
  md = fixture_md(group_tag = "HP")
  vd = read_variants(fixture_vcf(), "chr1:1000-7000")
  plain = plot_methylation(md)[[1L]]
  p = plot_methylation(md, variants = vd)
  rp = p[[1L]]
  expect_gt(length(rp$layers), length(plain$layers))
  geoms = layer_geoms(rp)
  colours = vapply(rp$layers, function(l) {
    col = l$aes_params$colour
    if (is.null(col)) NA_character_ else as.character(col)
  }, character(1L))

  # SNV (chr1:3000 A>G): asterisks on the ALT reads, at the SNV position
  snv = rp$layers[[which(geoms == "GeomPoint")]]$data
  expect_gt(nrow(snv), 0L)
  expect_true(all(snv$position == 3000L))
  # <DEL> 3899-4050: dark red line plus two end ticks
  expect_equal(sum(geoms == "GeomSegment" & colours %in% "#B71C1C"), 3L)
  # BND to chr7:5000: one label
  bnd = rp$layers[[which(geoms == "GeomText")]]$data
  expect_match(bnd$mate_label, "chr7:5,000")
  # SA reads matching the BND mate get the validated border
  expect_true(any(geoms == "GeomPolygon" & colours %in% "#880E4F"))
  # nothing is drawn on the smooth panel
  expect_identical(layer_geoms(p[[2L]]), layer_geoms(plot_methylation(md)[[2L]]))
  expect_no_warning(render(p))
})

# --- colours --------------------------------------------------------------------

test_that("colour_low and colour_high set the probability gradient", {
  reads = data.frame(read_name = "r1", start = 1000L, end = 2000L, strand = "+",
                     clip_side = NA_character_, sa_chrom = NA_character_,
                     stringsAsFactors = FALSE)
  sites = data.frame(read_name = "r1", position = c(1200, 1500, 1800),
                     mod_prob = c(0, 0.5, 1), mod_code = "m",
                     stringsAsFactors = FALSE)
  md = make_test_data(reads, sites)
  p = plot_methylation(md, colour_low = "#112233", colour_high = "#DDEEFF",
                       show_supplementary = FALSE)
  rp = p[[1L]]
  seg = ggplot2::ggplot_build(rp)$data[[which(layer_geoms(rp) == "GeomSegment")[1L]]]
  seg = seg[order(seg$x), ]
  expect_identical(toupper(seg$colour[c(1L, 3L)]), c("#112233", "#DDEEFF"))
  expect_false(seg$colour[2L] %in% c("#112233", "#DDEEFF"))
  expect_no_warning(render(p))

  # On the fixture: binary calls are drawn in exactly the given colours
  fx = fixture_md()
  p = plot_methylation(fx, colour_low = "navy", colour_high = "gold",
                       colour_ambiguous = "pink", call_threshold = 0.5,
                       call_ambiguous = 0.2)
  rp = p[[1L]]
  calls = ggplot2::ggplot_build(rp)$data[[which(layer_geoms(rp) == "GeomSegment")[1L]]]
  expect_setequal(unique(calls$colour), c("navy", "gold", "pink"))
  expect_no_warning(render(p))
})

test_that("colour_strand fills read bars with strand_colours", {
  md = fixture_md()
  p = plot_methylation(md, colour_strand = TRUE,
                       strand_colours = c("+" = "orchid", "-" = "darkgreen"))
  rp = p[[1L]]
  bars = rp$layers[[1L]]
  expect_identical(rlang::as_label(bars$mapping$fill), "strand")
  built = ggplot2::ggplot_build(rp)
  fill_col = grep("^fill", names(built$data[[1L]]), value = TRUE)[1L]
  expect_setequal(unique(built$data[[1L]][[fill_col]]), c("orchid", "darkgreen"))
  strand_scale = Filter(function(s) identical(s$name, "Strand"), built$plot$scales$scales)
  expect_length(strand_scale, 1L)
  expect_identical(unname(strand_scale[[1L]]$map(c("+", "-"))), c("orchid", "darkgreen"))
  expect_no_warning(render(p))

  # Without colour_strand the bars are a fixed grey
  plain = ggplot2::ggplot_build(plot_methylation(md)[[1L]])$data[[1L]]
  expect_identical(unique(plain$fill), "#B0BEC5")

  # Grouped data: group colour wins, with a message
  grouped = fixture_md(group_tag = "HP")
  expect_message(pg <- plot_methylation(grouped, colour_strand = TRUE,
                                        strand_colours = c("+" = "orchid", "-" = "darkgreen")),
                 "colour_strand")
  gfill = ggplot2::ggplot_build(pg[[1L]])$data[[1L]]
  gfill = gfill[[grep("^fill", names(gfill), value = TRUE)[1L]]]
  expect_false(any(c("orchid", "darkgreen") %in% gfill))
  expect_true(all(c("#0072B2", "#E69F00") %in% gfill))
})

# --- NA groups ------------------------------------------------------------------

# Grouped data whose reads include one with a missing group label.
make_grouped_plot_data = function() {
  gr = GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges   = IRanges::IRanges(start = 100L, end = 300L)
  )

  reads = data.frame(
    read_name = c("r1", "r2", "r3"),
    start     = c(100L, 120L, 180L),
    end       = c(170L, 190L, 260L),
    strand    = c("+", "+", "-"),
    group     = c("1", NA, "2"),
    stringsAsFactors = FALSE
  )

  sites = data.frame(
    read_name = c("r1", "r1", "r2", "r2", "r3", "r3"),
    position  = c(105L, 130L, 125L, 150L, 200L, 230L),
    mod_prob  = c(0.1, 0.3, 0.4, 0.6, 0.5, 0.7),
    mod_code  = rep("m", 6L),
    group     = c("1", "1", NA, NA, "2", "2"),
    stringsAsFactors = FALSE
  )

  cigar_features = data.frame(
    type        = c("D", "D", "D"),
    ref_start   = c(135L, 140L, 210L),
    ref_end     = c(145L, 150L, 220L),
    query_start = c(NA_integer_, NA_integer_, NA_integer_),
    query_end   = c(NA_integer_, NA_integer_, NA_integer_),
    length      = c(11L, 11L, 11L),
    read_name   = c("r1", "r2", "r3"),
    stringsAsFactors = FALSE
  )

  structure(
    list(
      reads          = reads,
      sites          = sites,
      region         = gr,
      mod_code       = "m",
      group_tag      = "HP",
      snv_position   = NULL,
      sequences      = setNames(character(0), character(0)),
      cigars         = setNames(character(0), character(0)),
      cigar_features = cigar_features
    ),
    class = "methylation_data"
  )
}

test_that("plot_methylation handles retained NA groups when showing CIGAR deletions", {
  md = make_grouped_plot_data()

  p = expect_no_error(
    plot_methylation(md, show_cigar = TRUE, min_indel_size = 1L)
  )

  expect_s3_class(p, "patchwork")

  top_panel = p$patches$plots[[1L]]
  poly_data = ggplot2::layer_data(top_panel, 1L)

  expect_true(nrow(poly_data) > 0L)
  expect_true("grey50" %in% poly_data$fill)
})
