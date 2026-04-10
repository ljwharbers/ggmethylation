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
      dot_size        = 1,
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
    dot_size        = 1,
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
