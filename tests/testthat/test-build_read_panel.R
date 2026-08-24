# Tests for build_read_panel and related helpers in R/build_read_panel.R

# Helper: build a minimal methylation_data object
make_test_data <- function(reads_df, sites_df,
                            region_start = 1000L, region_end = 2000L) {
  gr <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges   = IRanges::IRanges(start = region_start, end = region_end)
  )
  structure(
    list(
      reads        = reads_df,
      sites        = sites_df,
      region       = gr,
      mod_code     = "m",
      group_tag    = NULL,
      cigar_features = data.frame(
        read_name  = character(0),
        type       = character(0),
        ref_start  = integer(0),
        ref_end    = integer(0),
        length     = integer(0),
        stringsAsFactors = FALSE
      )
    ),
    class = "methylation_data"
  )
}

# --- .make_read_polygons ---

test_that(".make_read_polygons produces arrow for + strand", {
  reads <- data.frame(
    read_name       = "r1",
    start           = 1000L,
    end             = 1500L,
    strand          = "+",
    lane            = 1L,
    is_first_segment = TRUE,
    is_last_segment  = TRUE,
    stringsAsFactors = FALSE
  )
  polys <- ggmethylation:::.make_read_polygons(reads, arrow_w = 10, half_height = 0.35)
  # + strand arrow: 5 vertices (tip extends beyond end)
  expect_equal(nrow(polys), 5L)
  # Arrow tip x is end + arrow_w
  expect_true(any(polys$x > 1500L))
})

test_that(".make_read_polygons produces arrow for - strand", {
  reads <- data.frame(
    read_name        = "r1",
    start            = 1000L,
    end              = 1500L,
    strand           = "-",
    lane             = 1L,
    is_first_segment = TRUE,
    is_last_segment  = TRUE,
    stringsAsFactors = FALSE
  )
  polys <- ggmethylation:::.make_read_polygons(reads, arrow_w = 10, half_height = 0.35)
  # - strand arrow: 5 vertices (tip extends before start)
  expect_equal(nrow(polys), 5L)
  expect_true(any(polys$x < 1000L))
})

# --- .make_sa_overlay_polygons ---

# Helper: one SA-annotated read row, with the columns the overlay reads.
sa_read_row <- function(sa_side, clip_side = "both", strand = "+",
                        start = 1000L, end = 1500L, ...) {
  base <- data.frame(
    read_name        = "r1",
    start            = start,
    end              = end,
    strand           = strand,
    lane             = 1L,
    clip_side        = clip_side,
    sa_chrom         = "chr7",
    sa_side          = sa_side,
    is_first_segment = TRUE,
    is_last_segment  = TRUE,
    stringsAsFactors = FALSE
  )
  extra <- list(...)
  for (nm in names(extra)) base[[nm]] <- extra[[nm]]
  base
}

test_that(".make_sa_overlay_polygons marks only the breakpoint flank", {
  polys <- ggmethylation:::.make_sa_overlay_polygons(
    sa_read_row("right"), arrow_w = 10, half_height = 0.35,
    region_start = 1000L, region_end = 2000L
  )
  expect_equal(length(unique(polys$polygon_id)), 1L)
  expect_true(all(polys$x >= 1500L))
})

test_that(".make_sa_overlay_polygons marks the left flank when that is the side", {
  polys <- ggmethylation:::.make_sa_overlay_polygons(
    sa_read_row("left"), arrow_w = 10, half_height = 0.35,
    region_start = 1000L, region_end = 2000L
  )
  expect_equal(length(unique(polys$polygon_id)), 1L)
  expect_true(all(polys$x <= 1000L))
})

test_that(".make_sa_overlay_polygons ignores clip_side entirely", {
  # The old behaviour drew on both ends here, because adapter trimming makes
  # clip_side == "both" for nearly every long read.
  polys <- ggmethylation:::.make_sa_overlay_polygons(
    sa_read_row(NA_character_, clip_side = "both"),
    arrow_w = 10, half_height = 0.35,
    region_start = 1000L, region_end = 2000L
  )
  expect_equal(nrow(polys), 0L)
})

test_that(".make_sa_overlay_polygons draws two flanks only when both are real", {
  polys <- ggmethylation:::.make_sa_overlay_polygons(
    sa_read_row("both"), arrow_w = 10, half_height = 0.35,
    region_start = 1000L, region_end = 2000L
  )
  expect_equal(length(unique(polys$polygon_id)), 2L)
  expect_true(any(polys$x <= 1000L))
  expect_true(any(polys$x >= 1500L))
})

test_that(".make_sa_overlay_polygons colours each flank by its own partner", {
  polys <- ggmethylation:::.make_sa_overlay_polygons(
    sa_read_row("both", sa_chrom_left = "chr3", sa_chrom_right = "chr9"),
    arrow_w = 10, half_height = 0.35,
    region_start = 1000L, region_end = 2000L
  )
  left  <- polys[polys$x <= 1000L, , drop = FALSE]
  right <- polys[polys$x >= 1500L, , drop = FALSE]
  expect_equal(unique(left$sa_chrom), "chr3")
  expect_equal(unique(right$sa_chrom), "chr9")
})

test_that(".make_sa_overlay_polygons skips interior edges of a split read", {
  # A read split on two large deletions arrives as three rows; only the last
  # one carries the read's true right end.
  segs <- do.call(rbind, list(
    sa_read_row("right", start = 1000L, end = 1100L,
                is_first_segment = TRUE,  is_last_segment = FALSE),
    sa_read_row("right", start = 1200L, end = 1300L,
                is_first_segment = FALSE, is_last_segment = FALSE),
    sa_read_row("right", start = 1400L, end = 1500L,
                is_first_segment = FALSE, is_last_segment = TRUE)
  ))
  polys <- ggmethylation:::.make_sa_overlay_polygons(
    segs, arrow_w = 10, half_height = 0.35,
    region_start = 1000L, region_end = 2000L
  )
  expect_equal(length(unique(polys$polygon_id)), 1L)
  expect_true(all(polys$x >= 1500L))
})

test_that(".make_sa_overlay_polygons returns no rows without an sa_side column", {
  reads <- sa_read_row("right")
  reads$sa_side <- NULL
  polys <- ggmethylation:::.make_sa_overlay_polygons(
    reads, arrow_w = 10, half_height = 0.35,
    region_start = 1000L, region_end = 2000L
  )
  expect_equal(nrow(polys), 0L)
})

# --- .ensure_sa_side (back-compatibility) ---

test_that(".ensure_sa_side back-fills legacy objects and warns", {
  legacy <- sa_read_row("right")
  legacy$sa_side <- NULL
  expect_warning(out <- ggmethylation:::.ensure_sa_side(legacy),
                 "predates per-side")
  expect_equal(out$sa_side, "both")
})

test_that(".ensure_sa_side leaves a current object untouched and silent", {
  current <- sa_read_row("right")
  expect_silent(out <- ggmethylation:::.ensure_sa_side(current))
  expect_equal(out$sa_side, "right")
})

# --- build_read_panel: dots at region boundary ---

test_that("build_read_panel does not warn about removed polygon rows", {
  # Read whose end == region_end: arrowhead tip was previously censored by
  # scale_x_continuous(limits=...) causing the polygon to vanish.
  region_start <- 1000L
  region_end   <- 2000L

  reads <- data.frame(
    read_name = "r1",
    start     = 1000L,
    end       = 2000L,  # right at region boundary
    strand    = "+",
    lane      = 1L,
    stringsAsFactors = FALSE
  )
  sites <- data.frame(
    read_name = "r1",
    position  = 1500L,
    mod_prob  = 0.8,
    mod_code  = "m",
    stringsAsFactors = FALSE
  )
  data <- make_test_data(reads, sites, region_start, region_end)

  # Should produce a plot without warnings about removed rows
  expect_no_warning(
    p <- ggmethylation:::build_read_panel(
      data            = data,
      separator_lanes = numeric(0),
      region_start    = region_start,
      region_end      = region_end,
      colour_low      = "blue",
      colour_high     = "red",
      line_width      = 1,
      colour_strand   = FALSE,
      strand_colours  = c("+" = "grey60", "-" = "grey60"),
      group_colours   = NULL,
      mod_code_shapes = c(m = 16L)
    )
  )
  expect_s3_class(p, "gg")
})

# --- build_read_panel: dots in deletion gaps excluded when show_cigar=TRUE ---

test_that("build_read_panel excludes dots in deletion gaps with show_cigar=TRUE", {
  region_start <- 1000L
  region_end   <- 3000L

  # Read with a large deletion in the middle
  reads <- data.frame(
    read_name = "r1",
    start     = 1000L,
    end       = 2500L,
    strand    = "+",
    lane      = 1L,
    stringsAsFactors = FALSE
  )
  # Site at position 1800 — inside the deletion gap (1200:2100)
  sites <- data.frame(
    read_name = "r1",
    position  = 1800L,
    mod_prob  = 0.5,
    mod_code  = "m",
    stringsAsFactors = FALSE
  )
  cigar_features <- data.frame(
    read_name = "r1",
    type      = "D",
    ref_start = 1200L,
    ref_end   = 2100L,
    length    = 900L,
    stringsAsFactors = FALSE
  )
  data <- make_test_data(reads, sites, region_start, region_end)
  data$cigar_features <- cigar_features

  # Extract sites_plot by building the panel; the dot at 1800 should be absent
  p <- ggmethylation:::build_read_panel(
    data            = data,
    separator_lanes = numeric(0),
    region_start    = region_start,
    region_end      = region_end,
    colour_low      = "blue",
    colour_high     = "red",
    line_width      = 1,
    colour_strand   = FALSE,
    strand_colours  = c("+" = "grey60", "-" = "grey60"),
    group_colours   = NULL,
    mod_code_shapes = c(m = 16L),
    show_cigar      = TRUE,
    cigar_features  = cigar_features,
    min_indel_size  = 50L
  )
  # Extract the geom_point layer data
  point_data <- ggplot2::layer_data(p, i = 2L)
  expect_equal(nrow(point_data), 0L)
})

# --- build_read_panel: supplementary indicators end to end ---

# Build a panel for one chimeric read and return its layer data.
sa_panel <- function(sa_side, sites) {
  reads <- data.frame(
    read_name        = "r1",
    start            = 1000L,
    end              = 2000L,
    strand           = "+",
    lane             = 1L,
    clip_side        = "both",     # what adapter trimming leaves behind
    sa_chrom         = "chr7",
    sa_side          = sa_side,
    stringsAsFactors = FALSE
  )
  data <- make_test_data(reads, sites, 1000L, 3000L)
  ggmethylation:::build_read_panel(
    data               = data,
    separator_lanes    = numeric(0),
    region_start       = 1000L,
    region_end         = 3000L,
    colour_low         = "blue",
    colour_high        = "red",
    line_width         = 1,
    colour_strand      = FALSE,
    strand_colours     = c("+" = "grey60", "-" = "grey60"),
    group_colours      = NULL,
    mod_code_shapes    = c(m = 16L),
    show_supplementary = TRUE
  )
}

test_that("build_read_panel draws one indicator, on the sa_side flank", {
  sites <- data.frame(
    read_name = "r1", position = 1500L, mod_prob = 0.8, mod_code = "m",
    stringsAsFactors = FALSE
  )
  p <- sa_panel("right", sites)

  # Layer 1 is the read bar; layer 2 is the SA overlay.
  sa_layer <- ggplot2::layer_data(p, i = 2L)
  expect_equal(length(unique(sa_layer$group)), 1L)
  expect_true(all(sa_layer$x >= 2000L))
})

test_that("build_read_panel suppresses dots only under the drawn indicator", {
  # With the old clip_side-driven logic both ends were marked, so the dot at
  # 1001 was dropped as well.
  sites <- data.frame(
    read_name = "r1",
    position  = c(1001L, 1500L, 1999L),
    mod_prob  = 0.8,
    mod_code  = "m",
    stringsAsFactors = FALSE
  )
  p <- sa_panel("right", sites)
  dots <- ggplot2::layer_data(p, i = 3L)
  expect_equal(sort(dots$x), c(1001, 1500))
})
