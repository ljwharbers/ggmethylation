# Tests for binarized call mode (call_mode = "binary"):
# .classify_calls() and .binarize_sites(), the read panel, and the
# fraction-methylated smooth/delta panels in plot_methylation().

# Local minimal methylation_data constructor (mirrors
# test-build_read_panel.R::make_test_data; not shared via a helper-*.R file,
# so duplicated here per testthat edition-3 auto-sourcing rules).
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

test_that("binary call classification splits at threshold", {
  probs <- c(0.1, 0.49, 0.5, 0.9)
  cls <- ggmethylation:::.classify_calls(probs, threshold = 0.5, ambiguous = NULL)
  expect_equal(cls, c("unmethylated", "unmethylated", "methylated", "methylated"))
})

test_that("ambiguous band labels near-threshold calls", {
  probs <- c(0.1, 0.45, 0.5, 0.55, 0.9)
  cls <- ggmethylation:::.classify_calls(probs, threshold = 0.5, ambiguous = 0.1)
  expect_equal(cls, c("unmethylated", "ambiguous", "methylated", "ambiguous", "methylated"))
})

test_that("plot_methylation renders in binary mode", {
  reads <- data.frame(
    read_name = "r1", start = 1000L, end = 2000L, strand = "+",
    lane = 0L, mean_mod_prob = 0.5, clip_side = NA_character_,
    sa_chrom = NA_character_, stringsAsFactors = FALSE
  )
  sites <- data.frame(
    read_name = "r1", position = seq(1100, 1900, length.out = 5),
    mod_prob = c(0.1, 0.3, 0.5, 0.7, 0.9), mod_code = "m",
    stringsAsFactors = FALSE
  )
  md <- make_test_data(reads, sites)
  p <- ggmethylation::plot_methylation(md, call_mode = "binary",
                                       show_supplementary = FALSE)
  expect_true(inherits(p, "ggplot") || inherits(p, "patchwork"))
})


# --- .binarize_sites() ------------------------------------------------------

test_that("binarize_sites maps probabilities to 0/1 at the threshold", {
  sites <- data.frame(
    read_name = "r1", position = 1:4,
    mod_prob = c(0.1, 0.49, 0.5, 0.9), mod_code = "m",
    stringsAsFactors = FALSE
  )
  out <- ggmethylation:::.binarize_sites(sites, threshold = 0.5)
  expect_equal(out$mod_prob, c(0, 0, 1, 1))
})

test_that("binarize_sites honours a non-default threshold", {
  sites <- data.frame(
    read_name = "r1", position = 1:3,
    mod_prob = c(0.7, 0.8, 0.85), mod_code = "m",
    stringsAsFactors = FALSE
  )
  out <- ggmethylation:::.binarize_sites(sites, threshold = 0.8)
  expect_equal(out$mod_prob, c(0, 1, 1))
})

test_that("binarize_sites drops no rows when ambiguous is NULL", {
  sites <- data.frame(
    read_name = "r1", position = 1:5,
    mod_prob = c(0.45, 0.48, 0.5, 0.52, 0.55), mod_code = "m",
    stringsAsFactors = FALSE
  )
  out <- ggmethylation:::.binarize_sites(sites, threshold = 0.5, ambiguous = NULL)
  expect_equal(nrow(out), 5L)
  expect_equal(out$mod_prob, c(0, 0, 1, 1, 1))
})

test_that("binarize_sites drops ambiguous rows from the denominator", {
  sites <- data.frame(
    read_name = "r1", position = 1:5,
    mod_prob = c(0.05, 0.45, 0.5, 0.55, 0.95), mod_code = "m",
    stringsAsFactors = FALSE
  )
  out <- ggmethylation:::.binarize_sites(sites, threshold = 0.5, ambiguous = 0.1)
  # 0.45 and 0.55 are ambiguous; 0.5 is always methylated
  expect_equal(nrow(out), 3L)
  expect_equal(out$position, c(1L, 3L, 5L))
  expect_equal(mean(out$mod_prob), 2 / 3)
})

test_that("binarize_sites preserves other columns and drops NA probabilities", {
  sites <- data.frame(
    read_name = c("r1", "r2", "r3"), position = 1:3,
    mod_prob = c(0.9, NA_real_, 0.1), mod_code = "m",
    group = c("A", "A", "B"),
    stringsAsFactors = FALSE
  )
  out <- ggmethylation:::.binarize_sites(sites, threshold = 0.5)
  expect_equal(nrow(out), 2L)
  expect_equal(out$read_name, c("r1", "r3"))
  expect_equal(out$group, c("A", "B"))
  expect_false(anyNA(out$mod_prob))
})

test_that("binarize_sites returns 0-row input unchanged", {
  sites <- data.frame(
    read_name = character(0), position = numeric(0),
    mod_prob = numeric(0), mod_code = character(0),
    stringsAsFactors = FALSE
  )
  out <- ggmethylation:::.binarize_sites(sites, threshold = 0.5)
  expect_equal(nrow(out), 0L)
})


# --- Smooth panel reflects the binary calls ---------------------------------

# Pull the y-scale name off a ggplot (the smooth panel sets it via
# scale_y_continuous(name = ...), not labs()).
y_scale_name <- function(p) {
  for (s in p$scales$scales) {
    if (any(s$aesthetics == "y")) return(s$name)
  }
  NA_character_
}

# Reads/sites where every call is confidently methylated at prob 0.8: the mean
# probability is 0.8 but the methylated fraction is 1.0.
make_uniform_md <- function(prob = 0.8, n_sites = 8L) {
  reads <- data.frame(
    read_name = "r1", start = 1000L, end = 2000L, strand = "+",
    lane = 0L, mean_mod_prob = prob, clip_side = NA_character_,
    sa_chrom = NA_character_, stringsAsFactors = FALSE
  )
  sites <- data.frame(
    read_name = "r1",
    position  = seq(1100, 1900, length.out = n_sites),
    mod_prob  = rep(prob, n_sites), mod_code = "m",
    stringsAsFactors = FALSE
  )
  make_test_data(reads, sites)
}

test_that("binary smooth panel plots fraction methylated, not mean probability", {
  md <- make_uniform_md(prob = 0.8)

  p_cont <- ggmethylation::plot_methylation(md, show_supplementary = FALSE)
  p_bin  <- ggmethylation::plot_methylation(md, call_mode = "binary",
                                            show_supplementary = FALSE)

  cont_vals <- p_cont[[2]]$data$mean_prob
  bin_vals  <- p_bin[[2]]$data$mean_prob

  expect_equal(mean(cont_vals, na.rm = TRUE), 0.8, tolerance = 1e-6)
  expect_equal(mean(bin_vals, na.rm = TRUE), 1.0, tolerance = 1e-6)
})

test_that("binary smooth panel relabels the y axis", {
  md <- make_uniform_md()
  p_cont <- ggmethylation::plot_methylation(md, show_supplementary = FALSE)
  p_bin  <- ggmethylation::plot_methylation(md, call_mode = "binary",
                                            show_supplementary = FALSE)
  expect_match(y_scale_name(p_cont[[2]]), "Mean modification")
  expect_match(y_scale_name(p_bin[[2]]), "Fraction")
})

test_that("continuous mode smooth panel is unchanged by call_threshold", {
  md <- make_uniform_md(prob = 0.8)
  p_a <- ggmethylation::plot_methylation(md, show_supplementary = FALSE)
  p_b <- ggmethylation::plot_methylation(md, call_threshold = 0.9,
                                         call_ambiguous = 0.3,
                                         show_supplementary = FALSE)
  expect_equal(p_a[[2]]$data$mean_prob, p_b[[2]]$data$mean_prob)
})

test_that("ambiguous sites are excluded from the smoothed fraction", {
  # Two reads per position: one confidently methylated, one ambiguous.
  # With the ambiguous call dropped the fraction is 1.0, not 0.5.
  positions <- seq(1100, 1900, length.out = 8L)
  reads <- data.frame(
    read_name = c("r1", "r2"), start = 1000L, end = 2000L, strand = "+",
    lane = 0:1, mean_mod_prob = c(0.95, 0.5), clip_side = NA_character_,
    sa_chrom = NA_character_, stringsAsFactors = FALSE
  )
  sites <- data.frame(
    read_name = rep(c("r1", "r2"), each = 8L),
    position  = rep(positions, times = 2L),
    mod_prob  = c(rep(0.95, 8L), rep(0.45, 8L)), mod_code = "m",
    stringsAsFactors = FALSE
  )
  md <- make_test_data(reads, sites)

  p_hard <- ggmethylation::plot_methylation(md, call_mode = "binary",
                                            show_supplementary = FALSE)
  p_amb  <- ggmethylation::plot_methylation(md, call_mode = "binary",
                                            call_ambiguous = 0.1,
                                            show_supplementary = FALSE)
  # Hard threshold: r2's 0.45 counts as unmethylated -> fraction 0.5
  expect_equal(mean(p_hard[[2]]$data$mean_prob, na.rm = TRUE), 0.5,
               tolerance = 1e-6)
  # Ambiguous band: r2's 0.45 is dropped entirely -> fraction 1.0
  expect_equal(mean(p_amb[[2]]$data$mean_prob, na.rm = TRUE), 1.0,
               tolerance = 1e-6)
})

test_that("all-ambiguous data still builds a plot", {
  positions <- seq(1100, 1900, length.out = 8L)
  reads <- data.frame(
    read_name = "r1", start = 1000L, end = 2000L, strand = "+",
    lane = 0L, mean_mod_prob = 0.5, clip_side = NA_character_,
    sa_chrom = NA_character_, stringsAsFactors = FALSE
  )
  sites <- data.frame(
    read_name = "r1", position = positions,
    mod_prob = rep(0.52, 8L), mod_code = "m",
    stringsAsFactors = FALSE
  )
  md <- make_test_data(reads, sites)

  expect_equal(
    nrow(ggmethylation:::.binarize_sites(md$sites, 0.5, ambiguous = 0.2)),
    0L
  )
  expect_silent(
    p <- ggmethylation::plot_methylation(md, call_mode = "binary",
                                         call_threshold = 0.5,
                                         call_ambiguous = 0.2,
                                         show_supplementary = FALSE)
  )
  expect_s3_class(p, "patchwork")
})


# --- Per-read mean_mod_prob (read sorting) ----------------------------------

test_that("mean_mod_prob becomes a fraction of methylated calls in binary mode", {
  reads <- data.frame(
    read_name = "r1", start = 1000L, end = 2000L, strand = "+",
    lane = 0L, mean_mod_prob = 0, clip_side = NA_character_,
    sa_chrom = NA_character_, stringsAsFactors = FALSE
  )
  sites <- data.frame(
    read_name = "r1", position = c(1100, 1200, 1300),
    mod_prob = c(0.9, 0.9, 0.1), mod_code = "m",
    stringsAsFactors = FALSE
  )
  md <- make_test_data(reads, sites)

  # mean_mod_prob is computed inside plot_methylation(); read it back off the
  # read-bar layer data of the top panel.
  read_mean <- function(p) {
    unique(p[[1]]$layers[[1]]$data$mean_mod_prob)
  }

  p_cont <- ggmethylation::plot_methylation(md, show_supplementary = FALSE)
  p_bin  <- ggmethylation::plot_methylation(md, call_mode = "binary",
                                            show_supplementary = FALSE)
  expect_equal(read_mean(p_cont), mean(c(0.9, 0.9, 0.1)), tolerance = 1e-6)
  expect_equal(read_mean(p_bin), 2 / 3, tolerance = 1e-6)
})


# --- Grouped + delta panel -------------------------------------------------

test_that("binary mode grouped delta panel builds and is relabelled", {
  positions <- seq(1100, 1900, length.out = 8L)
  reads <- data.frame(
    read_name = c("r1", "r2"), start = 1000L, end = 2000L, strand = "+",
    lane = 0:1, mean_mod_prob = c(0.9, 0.1), clip_side = NA_character_,
    sa_chrom = NA_character_, group = c("1", "2"),
    stringsAsFactors = FALSE
  )
  sites <- data.frame(
    read_name = rep(c("r1", "r2"), each = 8L),
    position  = rep(positions, times = 2L),
    mod_prob  = c(rep(0.9, 8L), rep(0.1, 8L)), mod_code = "m",
    group     = rep(c("1", "2"), each = 8L),
    stringsAsFactors = FALSE
  )
  md <- make_test_data(reads, sites)
  md$group_tag <- "HP"

  p <- ggmethylation::plot_methylation(md, call_mode = "binary",
                                       show_delta = TRUE,
                                       show_supplementary = FALSE)
  expect_s3_class(p, "patchwork")
  expect_match(p[[3]]$labels$y, "fraction methylated")
  # Group "1" is fully methylated, group "2" fully unmethylated -> delta = -1
  expect_equal(mean(p[[3]]$data$delta, na.rm = TRUE), -1, tolerance = 1e-6)
})
