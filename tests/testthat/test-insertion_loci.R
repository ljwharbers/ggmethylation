# Tests for insertion_loci.R

# Helper: build a minimal methylation_data-like list for unit testing
make_test_md <- function(reads_df, cigar_features_df, ins_sites_df = NULL,
                         region = GenomicRanges::GRanges(
                           "chr1", IRanges::IRanges(1L, 1000L))) {
  if (is.null(ins_sites_df)) {
    ins_sites_df <- data.frame(
      read_name  = character(0L),
      ref_anchor = integer(0L),
      query_pos  = integer(0L),
      ins_offset = integer(0L),
      ins_length = integer(0L),
      mod_prob   = numeric(0L),
      mod_code   = character(0L),
      stringsAsFactors = FALSE
    )
  }
  structure(
    list(
      reads           = reads_df,
      sites           = data.frame(position = integer(0), mod_prob = numeric(0),
                                   read_name = character(0), mod_code = character(0),
                                   stringsAsFactors = FALSE),
      insertion_sites = ins_sites_df,
      region          = region,
      mod_code        = "m",
      group_tag       = NULL,
      snv_position    = NULL,
      sequences       = setNames(character(0), character(0)),
      cigars          = setNames(character(0), character(0)),
      cigar_features  = cigar_features_df
    ),
    class = "methylation_data"
  )
}

make_reads <- function(names, starts, ends) {
  data.frame(
    read_name        = names,
    start            = as.integer(starts),
    end              = as.integer(ends),
    bam_pos          = as.integer(starts),
    strand           = rep("+", length(names)),
    is_supplementary = rep(FALSE, length(names)),
    sa_chrom         = rep(NA_character_, length(names)),
    sa_pos           = rep(NA_integer_, length(names)),
    clip_side        = rep(NA_character_, length(names)),
    stringsAsFactors = FALSE
  )
}

make_cf <- function(read_names, ref_starts, lengths) {
  data.frame(
    type        = "I",
    ref_start   = as.integer(ref_starts),
    ref_end     = NA_integer_,
    query_start = 1L,
    query_end   = as.integer(lengths),
    length      = as.integer(lengths),
    read_name   = read_names,
    stringsAsFactors = FALSE
  )
}

# --- insertion_sites accessor ---

test_that("insertion_sites() returns data$insertion_sites", {
  ins <- data.frame(read_name = "r1", ref_anchor = 100L, query_pos = 5L,
                    ins_offset = 1L, ins_length = 10L, mod_prob = 0.8,
                    mod_code = "m", stringsAsFactors = FALSE)
  md <- make_test_md(
    reads_df        = make_reads("r1", 1L, 200L),
    cigar_features_df = make_cf("r1", 100L, 10L),
    ins_sites_df    = ins
  )
  result <- ggmethylation::insertion_sites(md)
  expect_equal(nrow(result), 1L)
  expect_equal(result$ref_anchor, 100L)
})

test_that("insertion_sites() errors on non-methylation_data input", {
  expect_error(ggmethylation::insertion_sites(list()), "methylation_data")
})

# --- list_insertion_loci ---

test_that("list_insertion_loci() returns empty df when no insertions", {
  md <- make_test_md(
    reads_df          = make_reads(character(0), integer(0), integer(0)),
    cigar_features_df = data.frame(
      type = character(0), ref_start = integer(0), ref_end = integer(0),
      query_start = integer(0), query_end = integer(0), length = integer(0),
      read_name = character(0), stringsAsFactors = FALSE
    )
  )
  out <- ggmethylation::list_insertion_loci(md)
  expect_equal(nrow(out), 0L)
  expect_true("locus_id" %in% names(out))
})

test_that("list_insertion_loci() clusters reads with close positions and similar lengths", {
  # Two reads with insertions at pos 100 and 105 (within tol_pos=10) of length 50
  # One read with insertion at pos 500, length 50 (separate locus)
  reads <- make_reads(c("r1", "r2", "r3"), c(1, 1, 400), c(300, 300, 700))
  cf    <- make_cf(c("r1", "r2", "r3"), c(100L, 105L, 500L), c(50L, 50L, 50L))
  md    <- make_test_md(reads, cf)

  out <- ggmethylation::list_insertion_loci(md, tol_pos = 10L, tol_len = 0.20,
                                              min_reads = 2L)
  # Only the first cluster (r1+r2) meets min_reads=2; r3 is singleton
  expect_equal(nrow(out), 1L)
  expect_equal(out$n_carriers, 2L)
})

test_that("list_insertion_loci() separates clusters by length tolerance", {
  # Three reads at ref_start 100, but lengths 50, 50, 200 (200 > 20% off median 50)
  reads <- make_reads(c("r1", "r2", "r3"), c(1, 1, 1), c(500, 500, 500))
  cf    <- make_cf(c("r1", "r2", "r3"), c(100L, 100L, 100L), c(50L, 50L, 200L))
  md    <- make_test_md(reads, cf)

  out <- ggmethylation::list_insertion_loci(md, tol_pos = 10L, tol_len = 0.20,
                                              min_reads = 2L)
  # r1+r2 form one locus; r3 is singleton -> only 1 locus passes min_reads=2
  expect_equal(nrow(out), 1L)
  expect_equal(out$median_length, 50L)
})

test_that("list_insertion_loci() locus_id is deterministic", {
  reads <- make_reads(c("r1", "r2"), c(1, 1), c(300, 300))
  cf    <- make_cf(c("r1", "r2"), c(100L, 104L), c(60L, 60L))
  md    <- make_test_md(reads, cf)

  out1 <- ggmethylation::list_insertion_loci(md)
  out2 <- ggmethylation::list_insertion_loci(md)
  expect_equal(out1$locus_id, out2$locus_id)
  expect_match(out1$locus_id, "^INS_chr1_\\d+_\\d+bp$")
})

test_that("list_insertion_loci() counts non-carriers correctly", {
  # r1 and r2 carry the insertion; r3 spans the locus but does not
  reads <- make_reads(c("r1", "r2", "r3"), c(1, 1, 1), c(300, 300, 300))
  cf <- rbind(
    make_cf(c("r1", "r2"), c(100L, 102L), c(50L, 50L))
  )
  md <- make_test_md(reads, cf)

  out <- ggmethylation::list_insertion_loci(md, min_reads = 2L)
  expect_equal(out$n_carriers, 2L)
  expect_equal(out$n_noncarriers, 1L)
})

test_that("list_insertion_loci() populates mean_ins_mod_prob from insertion_sites", {
  reads <- make_reads(c("r1", "r2"), c(1, 1), c(300, 300))
  cf    <- make_cf(c("r1", "r2"), c(100L, 103L), c(50L, 50L))
  ins_sites <- data.frame(
    read_name  = c("r1", "r2"),
    ref_anchor = c(100L, 103L),
    query_pos  = c(5L, 5L),
    ins_offset = c(1L, 1L),
    ins_length = c(50L, 50L),
    mod_prob   = c(0.8, 0.4),
    mod_code   = c("m", "m"),
    stringsAsFactors = FALSE
  )
  md  <- make_test_md(reads, cf, ins_sites_df = ins_sites)
  out <- ggmethylation::list_insertion_loci(md, min_reads = 2L)

  expect_equal(out$mean_ins_mod_prob, mean(c(0.8, 0.4)))
})

# --- plot_insertion_locus smoke tests ---

test_that("plot_insertion_locus() returns a ggplot object (show_smoothed = FALSE)", {
  reads <- make_reads(c("r1", "r2", "r3"), c(1L, 1L, 1L), c(300L, 300L, 300L))
  cf    <- make_cf(c("r1", "r2"), c(100L, 103L), c(50L, 50L))
  md    <- make_test_md(reads, cf)

  loci <- ggmethylation::list_insertion_loci(md, min_reads = 2L)
  expect_equal(nrow(loci), 1L)

  p <- ggmethylation::plot_insertion_locus(md, loci[1L, ], show_smoothed = FALSE)
  expect_true(inherits(p, c("gg", "patchwork")))
})

test_that("plot_insertion_locus() returns patchwork when show_smoothed = TRUE", {
  ins_sites <- data.frame(
    read_name  = c("r1", "r2"),
    ref_anchor = c(100L, 103L),
    query_pos  = c(55L, 55L),
    ins_offset = c(1L, 1L),
    ins_length = c(50L, 50L),
    mod_prob   = c(0.8, 0.6),
    mod_code   = c("m", "m"),
    stringsAsFactors = FALSE
  )
  reads <- make_reads(c("r1", "r2", "r3"), c(1L, 1L, 1L), c(300L, 300L, 300L))
  cf    <- make_cf(c("r1", "r2"), c(100L, 103L), c(50L, 50L))
  md    <- make_test_md(reads, cf, ins_sites_df = ins_sites)

  loci <- ggmethylation::list_insertion_loci(md, min_reads = 2L)
  p <- ggmethylation::plot_insertion_locus(md, loci[1L, ], show_smoothed = TRUE)
  expect_true(inherits(p, c("gg", "patchwork")))
})

test_that("plot_insertion_locus() requires one row of list_insertion_loci()", {
  reads <- make_reads(c("r1", "r2"), c(1L, 1L), c(300L, 300L))
  cf    <- make_cf(c("r1", "r2"), c(100L, 103L), c(50L, 50L))
  md    <- make_test_md(reads, cf)
  loci  <- ggmethylation::list_insertion_loci(md, min_reads = 2L)

  expect_error(ggmethylation::plot_insertion_locus(md, "INS_chr1_99999_100bp"),
               "list_insertion_loci")
  expect_error(ggmethylation::plot_insertion_locus(md, loci[0L, ]),
               "list_insertion_loci")
})

test_that("list_insertion_loci() returns each locus' reads and tolerances", {
  reads <- make_reads(c("r1", "r2", "r3", "r4"), 1L, 300L)
  cf    <- make_cf(c("r1", "r2"), c(100L, 103L), c(50L, 52L))
  md    <- make_test_md(reads, cf)

  loci <- ggmethylation::list_insertion_loci(md, tol_pos = 12L, tol_len = 0.1,
                                             min_reads = 2L)
  expect_equal(loci$carriers[[1L]], c(r1 = 50L, r2 = 52L))
  expect_equal(loci$noncarriers[[1L]], c("r3", "r4"))
  expect_equal(loci$n_noncarriers, 2L)
  expect_equal(loci$tol_pos, 12L)
  expect_equal(loci$tol_len, 0.1)

  empty <- ggmethylation::list_insertion_loci(make_test_md(reads, cf[0L, ]))
  expect_true(all(c("carriers", "noncarriers", "tol_pos", "tol_len") %in% names(empty)))
})

test_that("plot_insertion_locus() respects include_noncarriers = FALSE", {
  reads <- make_reads(c("r1", "r2", "r3"), c(1L, 1L, 1L), c(300L, 300L, 300L))
  cf    <- make_cf(c("r1", "r2"), c(100L, 103L), c(50L, 50L))
  md    <- make_test_md(reads, cf)

  loci <- ggmethylation::list_insertion_loci(md, min_reads = 2L)
  p <- ggmethylation::plot_insertion_locus(md, loci[1L, ],
                                            show_smoothed = FALSE,
                                            include_noncarriers = FALSE)
  expect_true(inherits(p, c("gg", "patchwork")))
})

test_that("plot_insertion_locus() uses the shared package colour palette", {
  # Defaults should match plot_methylation()'s shared grey -> red gradient.
  fun_args <- formals(ggmethylation::plot_insertion_locus)
  expect_equal(eval(fun_args$colour_low), "#BDBDBD")
  expect_equal(eval(fun_args$colour_high), "#C62828")

  # Enough per-read sites (>= 4 per group/region) are needed for the smoothed
  # panel to actually build a patchwork (rather than falling back to the read
  # panel alone) -- see smooth_piece()'s `length(x_vals) < 4L` guard.
  ins_sites <- data.frame(
    read_name  = c("r1", "r2"),
    ref_anchor = c(100L, 103L),
    query_pos  = c(55L, 55L),
    ins_offset = c(1L, 1L),
    ins_length = c(50L, 50L),
    mod_prob   = c(0.8, 0.6),
    mod_code   = c("m", "m"),
    stringsAsFactors = FALSE
  )
  sites <- data.frame(
    read_name = rep(c("r1", "r2"), each = 4L),
    position  = c(60L, 65L, 70L, 75L, 62L, 67L, 72L, 77L),
    mod_prob  = c(0.1, 0.3, 0.5, 0.7, 0.2, 0.4, 0.6, 0.8),
    mod_code  = "m",
    stringsAsFactors = FALSE
  )
  reads <- make_reads(c("r1", "r2", "r3"), c(1L, 1L, 1L), c(300L, 300L, 300L))
  cf    <- make_cf(c("r1", "r2"), c(100L, 103L), c(50L, 50L))
  md    <- make_test_md(reads, cf, ins_sites_df = ins_sites)
  md$sites <- sites

  loci <- ggmethylation::list_insertion_loci(md, min_reads = 2L)
  p <- ggmethylation::plot_insertion_locus(md, loci[1L, ], show_smoothed = TRUE)
  expect_true(inherits(p, "patchwork"))

  # Read-panel carrier fill should be the Okabe-Ito blue, non-carrier neutral grey.
  poly_layer_data <- ggplot2::layer_data(p[[1]], 1L)
  expect_true("#0072B2" %in% poly_layer_data$fill)
  expect_true("#999999" %in% poly_layer_data$fill)

  # Smooth-panel line colours should match the same palette.
  smooth_layer_data <- ggplot2::layer_data(p[[2]], 1L)
  expect_true(all(smooth_layer_data$colour %in% c("#0072B2", "#999999")))
})
