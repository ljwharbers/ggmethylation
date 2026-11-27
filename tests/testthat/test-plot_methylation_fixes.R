# Regression tests for the 0.5.0 main-plot audit (M1-M13): colours, gaps,
# binary-mode clamping, mod codes, SA markers, legends and argument checks.

fixture_md = function(...) {
  set.seed(1)
  suppressWarnings(read_methylation(
    testthat::test_path("fixtures", "synthetic.bam"), "chr1:1000-7000", ...
  ))
}

# Minimal methylation_data (testthat edition 3 has no shared helper file here).
make_md = function(reads, sites, region_start = 1000L, region_end = 2000L,
                   cigar_features = NULL, group_tag = NULL, mod_code = "m") {
  if (is.null(cigar_features)) {
    cigar_features = data.frame(read_name = character(0), type = character(0),
                                ref_start = integer(0), ref_end = integer(0),
                                length = integer(0), stringsAsFactors = FALSE)
  }
  structure(
    list(
      reads = reads, sites = sites,
      region = GenomicRanges::GRanges("chr1", IRanges::IRanges(region_start, region_end)),
      mod_code = mod_code, group_tag = group_tag, cigar_features = cigar_features
    ),
    class = "methylation_data"
  )
}

make_reads = function(read_name, start = 1000L, end = 2000L, group = NULL,
                      strand = "+", clip_side = NA_character_,
                      sa_chrom = NA_character_) {
  df = data.frame(read_name = read_name, start = as.integer(start),
                  end = as.integer(end), strand = strand, clip_side = clip_side,
                  sa_chrom = sa_chrom, stringsAsFactors = FALSE)
  if (!is.null(group)) df$group = group
  df
}

# Both helpers draw on a throwaway PNG device (no stray Rplots.pdf). A bitmap
# device is used because pdf(NULL) failed to measure the delta panel's
# "Δ" label under R CMD check ("conversion failure ... in 'mbcsToSbcs'").

# Render every panel to grobs: scale and "Removed rows" warnings fire here,
# not in plot_methylation().
render = function(p) {
  f = tempfile(fileext = ".png")
  grDevices::png(f)
  on.exit({ grDevices::dev.off(); unlink(f) })
  invisible(patchwork::patchworkGrob(p))
}

# Polygon pieces drawn by one layer (each NA-separated stretch is one child).
grob_pieces = function(p, layer) {
  f = tempfile(fileext = ".png")
  grDevices::png(f)
  on.exit({ grDevices::dev.off(); unlink(f) })
  ggplot2::layer_grob(p, layer)[[1L]]$children
}

layer_of = function(p, geom_class) {
  which(vapply(p$layers, function(l) inherits(l$geom, geom_class), logical(1L)))
}

# --- M1: groups the palette does not name ------------------------------------

test_that("SNV groups get distinct colours in every panel, without scale warnings", {
  md = fixture_md(snv = list(position = 3000, ref = "A", alt = "G"))
  expect_setequal(unique(md$reads$group), c("ALT", "REF"))
  p = plot_methylation(md)
  expect_no_warning(render(p))

  ext = ggmethylation:::.GROUP_PALETTE_EXTENSION
  bars = ggplot2::ggplot_build(p[[1L]])$data[[1L]]
  expect_setequal(unique(bars$fill), ext[1:2])
  smooth = p[[2L]]
  lines = ggplot2::ggplot_build(smooth)$data[[layer_of(smooth, "GeomLine")]]
  expect_setequal(unique(lines$colour), ext[1:2])
})

test_that("a third haplotype gets its own colour; 1/2 keep the defaults", {
  reads = make_reads(c("r1", "r2", "r3", "r4"), group = c("1", "2", "3", NA))
  pos = seq(1100, 1900, by = 100)
  sites = data.frame(read_name = rep(reads$read_name, each = length(pos)),
                     position = rep(pos, 4), mod_prob = 0.5, mod_code = "m",
                     group = rep(reads$group, each = length(pos)),
                     stringsAsFactors = FALSE)
  md = make_md(reads, sites, group_tag = "HP")
  p = plot_methylation(md, show_supplementary = FALSE)
  expect_no_warning(render(p))
  bars = ggplot2::ggplot_build(p[[1L]])$data[[1L]]
  expect_setequal(unique(bars$fill),
                  c("#0072B2", "#E69F00", "#009E73", ggmethylation:::.GROUP_NA_COLOUR))
})

test_that("a grouped object with zero sites renders without scale warnings", {
  reads = make_reads(c("r1", "r2"), group = c("1", "2"))
  sites = data.frame(read_name = character(0), position = integer(0),
                     mod_prob = numeric(0), mod_code = character(0),
                     group = character(0), stringsAsFactors = FALSE)
  md = make_md(reads, sites, group_tag = "HP")
  expect_no_warning(render(plot_methylation(md, show_supplementary = FALSE)))
})

# --- M3: area and ribbon split at consensus deletions ------------------------

# Group "1" (one read) carries a 1400-1500 deletion, so it has no sites there;
# group "2" covers everything.
deletion_md = function() {
  reads = make_reads(c("r1", "r2"), group = c("1", "2"))
  pos1 = c(1100, 1200, 1300, 1600, 1700, 1800, 1900)
  pos2 = seq(1100, 1900, by = 100)
  sites = data.frame(
    read_name = c(rep("r1", length(pos1)), rep("r2", length(pos2))),
    position = c(pos1, pos2),
    mod_prob = c(rep(0.2, length(pos1)), rep(0.8, length(pos2))),
    mod_code = "m",
    group = c(rep("1", length(pos1)), rep("2", length(pos2))),
    stringsAsFactors = FALSE
  )
  cf = data.frame(read_name = "r1", type = "D", ref_start = 1400L,
                  ref_end = 1500L, length = 101L, stringsAsFactors = FALSE)
  make_md(reads, sites, cigar_features = cf, group_tag = "HP")
}

test_that("the CI ribbon is split at a consensus deletion", {
  p = plot_methylation(deletion_md(), show_supplementary = FALSE, show_cigar = TRUE)
  expect_no_warning(render(p))
  smooth = p[[2L]]
  rib = layer_of(smooth, "GeomRibbon")
  expect_length(rib, 1L)
  pieces = grob_pieces(smooth, rib)
  # group 1: two pieces either side of the deletion; group 2: one piece
  expect_length(pieces, 3L)
})

test_that("the delta area is split at a consensus deletion", {
  p = plot_methylation(deletion_md(), show_supplementary = FALSE, show_cigar = TRUE,
                       show_delta = TRUE)
  expect_no_warning(render(p))
  delta = p[[3L]]
  pieces = grob_pieces(delta, 1L)
  expect_length(pieces, 2L)
  x = ggplot2::ggplot_build(delta)$data[[1L]]$x
  expect_false(any(x > 1400 & x < 1500))
})

test_that("delta bands meet at the zero crossing and drop all-zero bands", {
  df = data.frame(position = c(1, 2, 3, 4), delta = c(-0.2, 0.2, NA, 0.4))
  bands = ggmethylation:::.delta_ribbon_data(df)
  # crossing at x = 1.5 is present in both bands
  expect_true(all(c("pos_0", "neg_0") %in% bands$run))
  expect_true(1.5 %in% bands$position[bands$run == "pos_0"])
  expect_true(1.5 %in% bands$position[bands$run == "neg_0"])
  # after the NA: a new piece, positive only
  expect_true("pos_1" %in% bands$run)
  expect_false("neg_1" %in% bands$run)
})

# --- M4: binary-mode fit stays inside [0, 1] ---------------------------------

test_that("binary mode drops no smooth-line rows and raises no warning", {
  md = fixture_md(group_tag = "HP")
  p = plot_methylation(md, call_threshold = 0.5, call_ambiguous = 0.2)
  expect_no_warning(render(p))
  smooth = p[[2L]]
  expect_true(all(smooth$data$mean_prob >= 0 & smooth$data$mean_prob <= 1,
                  na.rm = TRUE))
  built = ggplot2::ggplot_build(smooth)
  line = built$data[[layer_of(smooth, "GeomLine")]]
  # every fitted value survives the y scale
  expect_equal(sum(!is.na(line$y)), sum(!is.na(smooth$data$mean_prob)))
})

# --- M5: several modification codes -----------------------------------------

test_that("multi-sample smooth keeps one line per sample, group and code", {
  md = fixture_md(group_tag = "HP", mod_code = c("m", "h"))
  p = plot_methylation(merge_methylation(A = md, B = md))
  smooth = p[[length(p)]]
  expect_setequal(unique(smooth$data$mod_code), c("m", "h"))
  line = ggplot2::ggplot_build(smooth)$data[[layer_of(smooth, "GeomLine")]]
  # 2 samples x 2 groups x 2 codes
  expect_equal(length(unique(line$group)), 8L)
})

test_that("ungrouped multi-sample smooth splits by code too", {
  md = fixture_md(mod_code = c("m", "h"))
  p = plot_methylation(merge_methylation(A = md, B = md))
  smooth = p[[length(p)]]
  line = ggplot2::ggplot_build(smooth)$data[[layer_of(smooth, "GeomLine")]]
  expect_equal(length(unique(line$group)), 4L)
})

test_that("the delta panel is skipped with a message for several codes", {
  md = fixture_md(group_tag = "HP", mod_code = c("m", "h"))
  expect_message(p <- plot_methylation(md, show_delta = TRUE),
                 "single modification code")
  expect_length(p, 2L)
})

test_that("calls for different codes at the same position do not overplot", {
  reads = make_reads("r1")
  sites = data.frame(read_name = "r1", position = c(1500, 1500),
                     mod_prob = c(0.9, 0.1), mod_code = c("m", "h"),
                     stringsAsFactors = FALSE)
  md = make_md(reads, sites, mod_code = c("m", "h"))
  p = plot_methylation(md, show_supplementary = FALSE, show_cigar = FALSE)
  rp = p[[1L]]
  seg = ggplot2::ggplot_build(rp)$data[[layer_of(rp, "GeomSegment")[1L]]]
  expect_equal(nrow(seg), 2L)
  # disjoint strips: one ends where the other starts
  spans = lapply(seq_len(2L), function(i) sort(c(seg$y[i], seg$yend[i])))
  expect_false(isTRUE(all.equal(spans[[1L]], spans[[2L]])))
  expect_match(rp$labels$y, "m, h")
})

# --- M6 / M7: supplementary markers ------------------------------------------

test_that("a split SA read gets one marker per clipped end", {
  reads = make_reads("r1", clip_side = "right", sa_chrom = "chr7")
  sites = data.frame(read_name = "r1", position = c(1100, 1900), mod_prob = 0.5,
                     mod_code = "m", stringsAsFactors = FALSE)
  cf = data.frame(read_name = "r1", type = "D", ref_start = 1400L,
                  ref_end = 1500L, length = 101L, stringsAsFactors = FALSE)
  md = make_md(reads, sites, cigar_features = cf)
  rp = plot_methylation(md, show_cigar = TRUE)[[1L]]
  polys = ggplot2::ggplot_build(rp)$data[layer_of(rp, "GeomPolygon")]
  sa = polys[[2L]]  # bars, then SA markers
  expect_equal(length(unique(sa$group)), 1L)
  expect_true(all(sa$x >= 2000))  # at the clipped right end, not at 1400
})

test_that("a site next to a clipped end is drawn", {
  reads = make_reads("r1", start = 1000L, end = 3000L, clip_side = "left",
                     sa_chrom = "chr3")
  sites = data.frame(read_name = "r1", position = c(1005, 2000), mod_prob = 0.5,
                     mod_code = "m", stringsAsFactors = FALSE)
  md = make_md(reads, sites, region_end = 5000L)
  rp = plot_methylation(md, show_cigar = FALSE)[[1L]]
  seg = ggplot2::ggplot_build(rp)$data[[layer_of(rp, "GeomSegment")[1L]]]
  expect_true(1005 %in% seg$x)
  # and the SA marker lies entirely left of the read start, not over the site
  sa = ggplot2::ggplot_build(rp)$data[[layer_of(rp, "GeomPolygon")[2L]]]
  expect_true(all(sa$x <= 1000))
})

# --- M10: one ribbon per line in a mixed merge -------------------------------

test_that("ungrouped samples in a mixed merge each get their own ribbon", {
  grouped = fixture_md(group_tag = "HP")
  ungrouped = fixture_md()
  p = plot_methylation(merge_methylation(A = grouped, B = ungrouped, C = ungrouped))
  smooth = p[[length(p)]]
  built = ggplot2::ggplot_build(smooth)
  n_rib  = length(unique(built$data[[layer_of(smooth, "GeomRibbon")]]$group))
  n_line = length(unique(built$data[[layer_of(smooth, "GeomLine")]]$group))
  expect_equal(n_line, 4L)
  expect_equal(n_rib, n_line)
})

# --- M11: consensus deletions use covering reads -----------------------------

test_that("consensus deletion fraction counts only reads covering the deletion", {
  cf = data.frame(type = "D", ref_start = c(1400L, 1400L), ref_end = c(1500L, 1500L),
                  length = 101L, read_name = c("r1", "r2"), stringsAsFactors = FALSE)
  # r1/r2 carry the deletion; r3-r8 are elsewhere in a wide window.
  reads = data.frame(read_name = paste0("r", 1:8),
                     start = c(1000L, 1000L, rep(5000L, 6)),
                     end   = c(2000L, 2000L, rep(6000L, 6)),
                     group = "1", stringsAsFactors = FALSE)
  res = ggmethylation:::.consensus_deletion_ranges(cf, reads, "group")
  expect_equal(nrow(res), 1L)
  # a covering read without the deletion still counts: 2 of 3 < 0.75
  reads$start[3] = 1300L
  reads$end[3] = 1600L
  res = ggmethylation:::.consensus_deletion_ranges(cf, reads, "group")
  expect_equal(nrow(res), 0L)
})

# --- M12: argument validation -------------------------------------------------

test_that("plot_methylation validates its numeric and colour arguments", {
  md = fixture_md(group_tag = "HP")
  expect_error(plot_methylation(md, call_threshold = NA_real_), "call_threshold")
  expect_error(plot_methylation(md, call_threshold = NA), "call_threshold")
  expect_error(plot_methylation(md, call_threshold = 0.5, call_ambiguous = "a"),
               "call_ambiguous")
  expect_error(plot_methylation(md, call_threshold = 0.5, call_ambiguous = -0.1),
               "call_ambiguous")
  expect_error(plot_methylation(md, call_threshold = 0.5, call_ambiguous = c(0.1, 0.2)),
               "call_ambiguous")
  expect_error(plot_methylation(md, smooth_span = -1), "smooth_span")
  expect_error(plot_methylation(md, smooth_span = "a"), "smooth_span")
  expect_error(plot_methylation(md, smooth_span = NA_real_), "smooth_span")
  expect_error(plot_methylation(md, group_colours = 3), "group_colours")
  expect_error(plot_methylation(md, group_colours = c("1" = "notacolour")),
               "invalid colour")
  expect_error(plot_methylation(md, group_colours = c(a = "red", a = "blue")),
               "group_colours")
  expect_error(plot_methylation(md, group_colours = "red"), "1 unnamed colour")
  expect_error(plot_methylation(md, line_width = -1), "line_width")
  expect_error(plot_methylation(md, line_width = "thin"), "line_width")
  expect_error(plot_methylation(md, min_indel_size = "x"), "min_indel_size")
  expect_error(plot_methylation(md, min_indel_size = NA_integer_), "min_indel_size")
  expect_error(plot_methylation(md, show_delta = NA), "show_delta")
  expect_error(plot_methylation(md, smooth_span = 0), "smooth_span")
  expect_error(plot_methylation(md, line_width = 0), "line_width")
  expect_error(plot_methylation(md, call_threshold = 1.5), "call_threshold")
  for (flag in c("colour_strand", "show_cigar", "show_supplementary", "show_ci")) {
    args = stats::setNames(list(md, "yes"), c("data", flag))
    expect_error(do.call(plot_methylation, args), sprintf("`%s` must be TRUE or FALSE", flag))
  }
})

test_that("plot_methylation rejects invalid colours before drawing anything", {
  md = fixture_md(group_tag = "HP")
  for (arg in c("colour_low", "colour_high", "colour_ambiguous")) {
    bad = stats::setNames(list(md, "notacolour"), c("data", arg))
    expect_error(do.call(plot_methylation, bad),
                 sprintf("`%s` contains invalid colour", arg))
    two = stats::setNames(list(md, c("red", "blue")), c("data", arg))
    expect_error(do.call(plot_methylation, two), sprintf("`%s` must be a single colour", arg))
    na = stats::setNames(list(md, NA_character_), c("data", arg))
    expect_error(do.call(plot_methylation, na), arg)
  }
  # strand_colours is checked (names and values) only when it is used
  mu = fixture_md()
  bad_strand = c("+" = "notacolour", "-" = "red")
  expect_error(plot_methylation(mu, colour_strand = TRUE, strand_colours = bad_strand),
               "`strand_colours` contains invalid colour")
  expect_error(plot_methylation(mu, colour_strand = TRUE, strand_colours = c("red", "blue")),
               "`strand_colours` must be a named vector")
  expect_s3_class(plot_methylation(mu, strand_colours = bad_strand), "patchwork")
})

test_that("plot_methylation validates panel_heights values and sort_by type", {
  md = fixture_md(group_tag = "HP")
  for (bad in list(c(-1, 1), c(0, 1), c(NA, 1), c(Inf, 1), c("1", "1"), numeric(0))) {
    expect_error(plot_methylation(md, panel_heights = bad),
                 "`panel_heights` must be NULL or a numeric vector")
  }
  for (bad in list(1, NA_character_, c("start", NA), character(0), TRUE)) {
    expect_error(plot_methylation(md, sort_by = bad),
                 "`sort_by` must be NULL or a character vector")
  }
})

test_that("plot_methylation's literal colour defaults match the internal constants", {
  f = formals(plot_methylation)
  expect_identical(eval(f$colour_low), ggmethylation:::.PROB_GRADIENT$low)
  expect_identical(eval(f$colour_high), ggmethylation:::.PROB_GRADIENT$high)
  expect_identical(eval(f$colour_ambiguous), ggmethylation:::.CALL_AMBIGUOUS_DEFAULT)
  expect_identical(eval(f$group_colours), ggmethylation:::.GROUP_PALETTE_DEFAULT)
})

test_that("the default group palette fills in unnamed groups silently", {
  # SNV groups REF/ALT are not named by the default HP palette: no message for
  # the default, one for a user-supplied palette that misses them.
  md = fixture_md(snv = list(position = 3000, ref = "A", alt = "G"))
  expect_no_message(plot_methylation(md))
  expect_message(plot_methylation(md, group_colours = c("1" = "red", "2" = "blue")),
                 "does not name group")
})

test_that("show_delta on ungrouped data says why there is no delta panel", {
  md = fixture_md()
  expect_message(p <- plot_methylation(md, show_delta = TRUE), "grouped data")
  expect_length(p, 2L)
})

# --- M13: legends and line types ---------------------------------------------

legend_titles = function(p) {
  g = render(p)
  boxes = g$grobs[grepl("^guide-box", g$layout$name)]
  titles = unlist(lapply(boxes, function(b) {
    unlist(lapply(b$grobs, function(leg) {
      if (!inherits(leg, "gtable")) return(NULL)
      ttl = leg$grobs[grepl("^title", leg$layout$name)]
      vapply(ttl, function(t) {
        txt = unlist(lapply(t$children, function(ch) ch$label))
        if (length(txt) == 0L) NA_character_ else paste(txt, collapse = "")
      }, character(1L))
    }))
  }))
  titles[!is.na(titles)]
}

test_that("shared legends appear once in a multi-sample figure", {
  grouped = fixture_md(group_tag = "HP")
  p = plot_methylation(merge_methylation(A = grouped, B = grouped))
  titles = legend_titles(p)
  expect_equal(sum(titles == "Group"), 1L)
  expect_equal(sum(titles == "SA partner"), 1L)
  expect_equal(sum(grepl("^Modification", titles)), 1L)
})

test_that("the smooth panel does not repeat the read panel's group legend", {
  p = plot_methylation(fixture_md(group_tag = "HP"))
  titles = legend_titles(p)
  expect_equal(sum(titles == "Group"), 1L)
})

test_that("SA partner colours are the same in every sample", {
  # Sample B has only one of A's two partner chromosomes; per-panel hue
  # palettes used to give chr7 a different colour in each panel.
  sites = function(rn) data.frame(read_name = rn, position = 1500, mod_prob = 0.5,
                                  mod_code = "m", stringsAsFactors = FALSE)
  a = make_md(make_reads(c("a1", "a2"), clip_side = "right",
                         sa_chrom = c("chr3", "chr7")), sites(c("a1", "a2")))
  b = make_md(make_reads("b1", clip_side = "right", sa_chrom = "chr7"), sites("b1"))
  p = plot_methylation(merge_methylation(A = a, B = b))
  sa_fills = function(rp) {
    d = ggplot2::ggplot_build(rp)$data[[layer_of(rp, "GeomPolygon")[2L]]]
    # drop zero-area legend placeholders
    width = stats::ave(d$x, d$group, FUN = function(v) diff(range(v)))
    unique(d$fill[width > 0])
  }
  hue = scales::hue_pal()(2L)  # chr3, chr7
  expect_setequal(sa_fills(p[[1L]]), hue)
  expect_equal(sa_fills(p[[2L]]), hue[2L])
  expect_equal(sum(legend_titles(p) == "SA partner"), 1L)
})

test_that("more than four modification codes get distinct line types", {
  codes = c("a", "b", "c", "d", "e")
  reads = make_reads(c("r1", "r2"), group = c("1", "2"))
  pos = seq(1100, 1900, by = 100)
  grid = expand.grid(position = pos, mod_code = codes, read_name = c("r1", "r2"),
                     stringsAsFactors = FALSE)
  grid$mod_prob = 0.5
  grid$group = ifelse(grid$read_name == "r1", "1", "2")
  md = make_md(reads, grid, group_tag = "HP", mod_code = codes)
  p = plot_methylation(md, show_supplementary = FALSE)
  expect_no_warning(render(p))
  smooth = p[[2L]]
  line = ggplot2::ggplot_build(smooth)$data[[layer_of(smooth, "GeomLine")]]
  expect_false(anyNA(line$linetype))
  expect_equal(length(unique(line$linetype)), 5L)
})

test_that("too many line types stop with a clear message", {
  expect_error(ggmethylation:::.linetype_values(letters[1:14], "samples"),
               "Too many samples")
})
