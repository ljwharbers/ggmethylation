# Plot-content tests for plot_insertion_locus().
#
# Stitched x convention (hand-derived, see R/plot_insertion_locus.R): with
# anchor A (the reference base right after the insertion) and drawn insertion
# width W, reference base p < A is drawn at x = p, inserted base k (1-based) at
# x = A + k - 1, and reference base p >= A at x = p + W. The boundaries sit
# half-way between bases, at A - 0.5 and A + W - 0.5.

# Synthetic locus at A = 1000: carriers c1, c2 (50 bp) and c3 (55 bp, anchored
# at 1001), so W = median(50, 50, 55) = 50 and c3 overflows; non-carriers n1, n2
# span the locus. All reads cover 800-1200.
pil_md = function(codes = "m") {
  rn = c("c1", "c2", "c3", "n1", "n2")
  reads = data.frame(
    read_name = rn, start = 800L, end = 1200L, bam_pos = 800L, strand = "+",
    is_supplementary = FALSE, sa_chrom = NA_character_, sa_pos = NA_integer_,
    clip_side = NA_character_, stringsAsFactors = FALSE
  )
  cf = data.frame(
    type = "I", ref_start = c(1000L, 1000L, 1001L), ref_end = NA_integer_,
    query_start = 201L, query_end = 200L + c(50L, 50L, 55L),
    length = c(50L, 50L, 55L), read_name = c("c1", "c2", "c3"),
    stringsAsFactors = FALSE
  )
  sites = expand.grid(position = c(990L, 999L, 1000L, 1010L), read_name = rn,
                      mod_code = codes, stringsAsFactors = FALSE)
  sites$mod_prob = 0.5
  ins = data.frame(
    read_name  = c("c1", "c1", "c2", "c3", "c3"),
    ref_anchor = c(1000L, 1000L, 1000L, 1001L, 1001L),
    ins_offset = c(1L, 50L, 25L, 1L, 53L),
    ins_length = c(50L, 50L, 50L, 55L, 55L),
    mod_prob   = 0.5,
    mod_code   = codes[1L],
    stringsAsFactors = FALSE
  )
  ins$query_pos = 200L + ins$ins_offset
  structure(
    list(reads = reads, sites = sites, insertion_sites = ins,
         region = GenomicRanges::GRanges("chr1", IRanges::IRanges(1L, 2000L)),
         mod_code = codes, group_tag = NULL, cigar_features = cf),
    class = "methylation_data"
  )
}

pil_locus = function(md) {
  loci = list_insertion_loci(md, tol_pos = 10L, tol_len = 0.2, min_reads = 2L)
  stopifnot(nrow(loci) == 1L)
  loci[1L, ]
}

layer_of = function(p, geom) {
  i = which(vapply(p$layers, function(l) inherits(l$geom, geom), logical(1L)))
  ggplot2::layer_data(p, i[1L])
}

test_that("synthetic locus is what the tests assume", {
  locus = pil_locus(pil_md())
  expect_equal(locus$anchor_pos, 1000L)
  expect_equal(locus$median_length, 50L)
  expect_equal(locus$carriers[[1L]], c(c1 = 50L, c2 = 50L, c3 = 55L))
  expect_equal(locus$noncarriers[[1L]], c("n1", "n2"))
})

# --- I3: reference base at the anchor is not drawn on the first inserted base ---

test_that("reference base at the anchor goes to the right flank (I3)", {
  md = pil_md()
  p = plot_insertion_locus(md, pil_locus(md), show_smoothed = FALSE)
  sites = layer_of(p, "GeomSegment")
  # 999 stays at 999, inserted base 1 is at 1000, reference 1000 moves to 1050
  # and reference 1010 to 1060.
  expect_true(all(c(999, 1000, 1050, 1060) %in% sites$x))
  # No two calls of the same read share a column.
  expect_false(any(duplicated(sites[, c("x", "y")])))
})

# --- I4: over-length inserted bases are clipped at the insertion boundary ---

test_that("inserted bases beyond the drawn width are not drawn in the right flank (I4)", {
  md = pil_md()
  p = plot_insertion_locus(md, pil_locus(md), show_smoothed = FALSE)
  sites = layer_of(p, "GeomSegment")
  # Expected columns, derived by hand: flanks 990, 999, 1050, 1060; inserted
  # offsets 1, 25, 50 -> 1000, 1024, 1049. c3's offset 53 (x = 1052) is clipped.
  expect_setequal(unique(sites$x), c(990, 999, 1000, 1024, 1049, 1050, 1060))
})

test_that("each smooth line stays inside one region (I3, I4)", {
  md = pil_md()
  p = plot_insertion_locus(md, pil_locus(md))
  sm = layer_of(p[[2]], "GeomLine")
  sm = sm[!is.na(sm$y), ]
  region = cut(sm$x, c(799.5, 999.5, 1049.5, 1249.5), labels = FALSE)
  expect_false(anyNA(region))
  expect_true(all(tapply(region, sm$group, function(r) length(unique(r))) == 1L))
  # The insertion's own line spans offsets 1..50 only.
  mid = sm[region == 2L, ]
  expect_equal(range(mid$x), c(1000, 1049))
})

# --- I5: non-carriers below carriers, centred separator, no lane numbers ---

test_that("non-carriers are drawn below carriers with a centred separator (I5)", {
  md = pil_md()
  p = plot_insertion_locus(md, pil_locus(md), show_smoothed = FALSE)
  polys = layer_of(p, "GeomPolygon")
  is_carrier = polys$fill == "#0072B2"
  expect_true(any(!is_carrier))
  # Larger built y is higher on the page whatever the scale direction.
  expect_gt(min(polys$y[is_carrier]), max(polys$y[!is_carrier]))
  # The separator sits in the middle of the gap between the two blocks.
  c_y = range(polys$y[is_carrier])
  n_y = range(polys$y[!is_carrier])
  gap = if (c_y[1L] > n_y[2L]) c(n_y[2L], c_y[1L]) else c(c_y[2L], n_y[1L])
  expect_equal(layer_of(p, "GeomHline")$yintercept, mean(gap))
  # Lane numbers are not shown.
  built_theme = ggplot2::calc_element("axis.text.y.left", ggplot2::theme_get() + p$theme)
  expect_s3_class(built_theme, "element_blank")
})

# --- I6: identical x ranges in the stacked panels ---

test_that("read and smooth panels share the same x range (I6)", {
  md = pil_md()
  p = plot_insertion_locus(md, pil_locus(md))
  xr = function(q) ggplot2::ggplot_build(q)$layout$panel_params[[1L]]$x.range
  expect_equal(xr(p[[1L]]), xr(p[[2L]]))
  # Window: 200 reference bases per flank plus the 50 bp insertion.
  expect_equal(xr(p[[1L]]), c(799.5, 1249.5))
})

# --- I7: stitched axis labels, mod codes, argument validation ---

test_that("x axis labels give reference coordinates in the flanks (I7)", {
  md = pil_md()
  check_axis = function(q) {
    pp = ggplot2::ggplot_build(q)$layout$panel_params[[1L]]
    br = pp$x$get_breaks()
    lab = pp$x$get_labels()
    keep = !is.na(br)
    br = br[keep]
    lab = lab[keep]
    left  = br < 999.5
    right = br > 1049.5
    mid   = !left & !right
    expect_true(any(left) && any(right) && sum(mid) == 1L)
    expect_equal(as.numeric(gsub(",", "", lab[left])), br[left])
    expect_equal(as.numeric(gsub(",", "", lab[right])), br[right] - 50)
    expect_match(lab[mid], "+1", fixed = TRUE)
    expect_match(lab[mid], "+50 bp", fixed = TRUE)
    el = ggplot2::calc_element("axis.text.x.bottom", ggplot2::theme_get() + q$theme)
    expect_false(inherits(el, "element_blank"))
  }
  p = plot_insertion_locus(md, pil_locus(md))
  check_axis(p[[2L]])
  check_axis(plot_insertion_locus(md, pil_locus(md), show_smoothed = FALSE))
})

test_that("different modification codes are not averaged together (I7)", {
  md = pil_md(codes = c("m", "h"))
  md$sites$mod_prob = ifelse(md$sites$mod_code == "m", 0.9, 0.1)
  md$insertion_sites$mod_prob = 0.9  # all inserted-base calls are code "m"
  p =plot_insertion_locus(md, pil_locus(md))
  sm = layer_of(p[[2L]], "GeomLine")
  sm = sm[!is.na(sm$y), ]
  expect_true(all(abs(sm$y - 0.9) < 1e-8 | abs(sm$y - 0.1) < 1e-8))
  expect_true(any(abs(sm$y - 0.1) < 1e-8))
  # One linetype per code.
  expect_equal(length(unique(sm$linetype)), 2L)
})

test_that("plot_insertion_locus() validates its arguments (I7)", {
  md = pil_md()
  locus = pil_locus(md)
  for (bad in list(NA, "a", -1, 0, c(10, 20))) {
    expect_error(plot_insertion_locus(md, locus, flank = bad), "flank")
  }
  for (bad in list(NA, "yes", c(TRUE, FALSE))) {
    expect_error(plot_insertion_locus(md, locus, include_noncarriers = bad),
                 "include_noncarriers")
    expect_error(plot_insertion_locus(md, locus, show_smoothed = bad),
                 "show_smoothed")
  }
})

# --- Calls drawn as in plot_methylation() ---

test_that("calls of different codes occupy disjoint strips of a lane", {
  md = pil_md(codes = c("m", "h"))
  md$sites$mod_prob = ifelse(md$sites$mod_code == "m", 1, 0)
  p = plot_insertion_locus(md, pil_locus(md), show_smoothed = FALSE)
  seg = layer_of(p, "GeomSegment")
  # Flank columns only (inserted-base calls are all code "m" here).
  seg = seg[seg$x %in% c(990, 999, 1050, 1060), ]
  seg$code = ifelse(seg$colour == "#C62828", "m", "h")
  expect_setequal(unique(seg$code), c("m", "h"))
  # Lane of each segment: its midpoint, rounded (sign-agnostic for the
  # reversed y scale).
  seg$lane = round(abs(seg$y + seg$yend) / 2)
  for (ln in unique(seg$lane)) {
    s = seg[seg$lane == ln, ]
    lo = function(cd) range(abs(c(s$y[s$code == cd], s$yend[s$code == cd])))
    rm = lo("m")
    rh = lo("h")
    expect_true(rm[2L] <= rh[1L] || rh[2L] <= rm[1L])
  }
  # Two calls at one column of one read no longer overplot each other.
  expect_false(any(duplicated(seg[, c("x", "y", "yend")])))
  # The strips are named on the y axis, as in plot_methylation().
  expect_match(p$labels$y, "m, h", fixed = TRUE)
})

test_that("call_threshold gives binary calls and a fraction-methylated panel", {
  md = pil_md()
  md$sites$mod_prob = rep(c(0.9, 0.1), length.out = nrow(md$sites))
  md$insertion_sites$mod_prob = 0.9
  locus = pil_locus(md)
  p = plot_insertion_locus(md, locus, call_threshold = 0.5)
  sc = p[[1L]]$scales$get_scales("colour")
  expect_s3_class(sc, "ScaleDiscrete")
  seg = layer_of(p[[1L]], "GeomSegment")
  expect_setequal(unique(seg$colour), c("#C62828", "#BDBDBD"))

  sm = layer_of(p[[2L]], "GeomLine")
  sm = sm[!is.na(sm$y), ]
  expect_true(all(sm$y >= 0 & sm$y <= 1))
  # Inserted bases are all called methylated: their line is exactly 1.
  mid = sm[sm$x >= 1000 & sm$x <= 1049, ]
  expect_true(nrow(mid) > 0L && all(abs(mid$y - 1) < 1e-8))
  expect_equal(p[[2L]]$scales$get_scales("y")$name, "Fraction methylated")

  # The continuous panel averages the raw probabilities instead.
  p_cont = plot_insertion_locus(md, locus)
  expect_equal(p_cont[[2L]]$scales$get_scales("y")$name, "Mean mod. prob.")
  expect_s3_class(p_cont[[1L]]$scales$get_scales("colour"), "ScaleContinuous")

  # Ambiguous calls get their own colour.
  md$sites$mod_prob[1:4] = 0.52
  p_amb = plot_insertion_locus(md, locus, call_threshold = 0.5,
                               call_ambiguous = 0.1, colour_ambiguous = "black",
                               show_smoothed = FALSE)
  expect_true("black" %in% layer_of(p_amb, "GeomSegment")$colour)
})

test_that("smooth-panel linetypes follow .linetype_values()", {
  codes = c("m", "h", "a")
  md = pil_md(codes = codes)
  md$sites$mod_prob = unname(c(m = 0.9, h = 0.5, a = 0.1)[md$sites$mod_code])
  md$insertion_sites$mod_prob = 0.9  # inserted-base calls are code "m"
  p =plot_insertion_locus(md, pil_locus(md))
  sm = layer_of(p[[2L]], "GeomLine")
  sm = sm[!is.na(sm$y), ]
  code_of_y = c("0.9" = "m", "0.5" = "h", "0.1" = "a")
  sm$code = code_of_y[format(round(sm$y, 1))]
  expected = .linetype_values(codes, "modification codes")
  got = vapply(codes, function(cd) unique(sm$linetype[sm$code == cd]), character(1L))
  expect_equal(got, expected[codes])
})

test_that("colour defaults are written out and match .PROB_GRADIENT", {
  f = formals(plot_insertion_locus)
  expect_type(f$colour_low, "character")
  expect_type(f$colour_high, "character")
  expect_equal(f$colour_low, .PROB_GRADIENT$low)
  expect_equal(f$colour_high, .PROB_GRADIENT$high)
})

test_that("plot_insertion_locus() validates the call and colour arguments", {
  md = pil_md()
  locus = pil_locus(md)
  for (bad in list(-0.1, 1.5, NA, "a", c(0.2, 0.4))) {
    expect_error(plot_insertion_locus(md, locus, call_threshold = bad),
                 "`call_threshold`")
  }
  expect_error(plot_insertion_locus(md, locus, call_ambiguous = 0.1),
               "`call_ambiguous`.*`call_threshold`")
  for (bad in list(-0.1, NA, "a")) {
    expect_error(plot_insertion_locus(md, locus, call_threshold = 0.5,
                                      call_ambiguous = bad), "`call_ambiguous`")
  }
  for (bad in list(0, -1, NA, "a")) {
    expect_error(plot_insertion_locus(md, locus, line_width = bad), "`line_width`")
  }
  expect_error(plot_insertion_locus(md, locus, colour_low = "notacolour"), "`colour_low`")
  expect_error(plot_insertion_locus(md, locus, colour_high = c("red", "blue")),
               "`colour_high`")
  expect_error(plot_insertion_locus(md, locus, colour_ambiguous = NA_character_),
               "`colour_ambiguous`")
})

test_that("list_insertion_loci() validates its arguments (I7)", {
  md = pil_md()
  for (bad in list(NA, "a", -1, c(1, 2))) {
    expect_error(list_insertion_loci(md, tol_pos = bad), "tol_pos")
    expect_error(list_insertion_loci(md, tol_len = bad), "tol_len")
  }
  for (bad in list(NA, "a", 0, -1)) {
    expect_error(list_insertion_loci(md, min_reads = bad), "min_reads")
    expect_error(list_insertion_loci(md, min_length = bad), "min_length")
  }
})
