test_that("compute_group_delta returns NULL and warns for non-2 groups", {
  # This is a warning rather than a message on purpose: a message is invisible
  # in scripts and knitr, so the delta panel used to vanish silently.
  sites <- data.frame(
    position = 1:5 * 100, mod_prob = seq(0.1, 0.5, length.out = 5),
    group = "A", stringsAsFactors = FALSE
  )
  expect_warning(
    res <- ggmethylation:::.compute_group_delta(sites, "group", span = NULL),
    "exactly two groups"
  )
  expect_null(res)
})

test_that("compute_group_delta returns signed delta on shared grid for 2 groups", {
  sites <- data.frame(
    position = rep(1:6 * 100, 2),
    mod_prob = c(rep(0.2, 6), rep(0.8, 6)),
    group    = rep(c("1", "2"), each = 6),
    stringsAsFactors = FALSE
  )
  res <- ggmethylation:::.compute_group_delta(sites, "group", span = NULL)
  expect_true(all(c("position", "delta", "sign") %in% names(res)))
  # group2 (0.8) - group1 (0.2) ~ +0.6 where both supported
  ok <- !is.na(res$delta)
  expect_true(mean(res$delta[ok]) > 0)
  expect_true(all(res$sign[ok & res$delta > 0] == "pos"))
})

test_that("compute_group_delta handles both groups with < 4 unique positions (approx fallback)", {
  # Both groups have only 3 unique positions -> smooth_methylation() falls
  # back to stats::approx(rule = 1) for each group on the shared grid.
  sites <- data.frame(
    position = c(100, 200, 300, 150, 250, 350),
    mod_prob = c(0.1, 0.2, 0.3, 0.9, 0.8, 0.7),
    group    = c("1", "1", "1", "2", "2", "2"),
    stringsAsFactors = FALSE
  )
  res <- ggmethylation:::.compute_group_delta(sites, "group", span = NULL, n_grid = 5L)

  # grid = seq(100, 350, length.out = 5) = 100, 162.5, 225, 287.5, 350
  expect_equal(res$position, c(100, 162.5, 225, 287.5, 350))

  # Hand-computed via stats::approx(rule = 1):
  # group "1": approx(c(100,200,300), c(0.1,0.2,0.3), xout = grid)
  #   -> 0.1, 0.1625, 0.225, 0.2875, NA (350 is outside range)
  # group "2": approx(c(150,250,350), c(0.9,0.8,0.7), xout = grid)
  #   -> NA (100 is outside range), 0.8875, 0.825, 0.7625, 0.7
  # delta = group2 - group1
  expected_delta <- c(NA, 0.8875 - 0.1625, 0.825 - 0.225, 0.7625 - 0.2875, NA)
  expect_equal(res$delta, expected_delta, tolerance = 1e-8)

  # Outside the overlapping support of both groups: NA delta and NA sign
  expect_true(is.na(res$delta[res$position == 100]))
  expect_true(is.na(res$delta[res$position == 350]))
  expect_true(is.na(res$sign[res$position == 100]))
  expect_true(is.na(res$sign[res$position == 350]))

  # Within overlapping support: non-NA, positive (group 2 > group 1 everywhere)
  inside <- res$position %in% c(162.5, 225, 287.5)
  expect_false(any(is.na(res$delta[inside])))
  expect_true(all(res$sign[inside] == "pos"))
})

test_that("compute_group_delta attaches the sorted group names", {
  sites <- data.frame(
    position = rep(1:6 * 100, 2),
    mod_prob = c(rep(0.2, 6), rep(0.8, 6)),
    group    = rep(c("2", "1"), each = 6),
    stringsAsFactors = FALSE
  )
  res <- ggmethylation:::.compute_group_delta(sites, "group", span = NULL)
  # Sorted, so that delta = groups[2] - groups[1] regardless of input order.
  expect_equal(attr(res, "groups"), c("1", "2"))
})

test_that("build_delta_panel honours explicit pos/neg fills", {
  df_pos <- data.frame(position = seq(1000, 2000, length.out = 50),
                       delta = 0.3, sign = "pos", stringsAsFactors = FALSE)
  df_neg <- data.frame(position = seq(1000, 2000, length.out = 50),
                       delta = -0.3, sign = "neg", stringsAsFactors = FALSE)

  p_pos <- ggmethylation:::.build_delta_panel(
    df_pos, 1000, 2000, "y", fill_pos = "purple", fill_neg = "darkgreen"
  )
  p_neg <- ggmethylation:::.build_delta_panel(
    df_neg, 1000, 2000, "y", fill_pos = "purple", fill_neg = "darkgreen"
  )

  expect_equal(unique(ggplot2::ggplot_build(p_pos)$data[[1]]$fill), "purple")
  expect_equal(unique(ggplot2::ggplot_build(p_neg)$data[[1]]$fill), "darkgreen")
})

test_that("build_delta_panel defaults to the standalone diverging palette", {
  df <- data.frame(position = seq(1000, 2000, length.out = 50),
                   delta = 0.3, sign = "pos", stringsAsFactors = FALSE)
  p <- ggmethylation:::.build_delta_panel(df, 1000, 2000)
  expect_equal(unique(ggplot2::ggplot_build(p)$data[[1]]$fill),
               ggmethylation:::.DELTA_DIVERGING$pos)
})

test_that("build_delta_panel returns a ggplot", {
  df <- data.frame(position = seq(1000, 2000, length.out = 50),
                   delta = sin(seq(0, 3, length.out = 50)) * 0.4,
                   stringsAsFactors = FALSE)
  df$sign <- ifelse(df$delta >= 0, "pos", "neg")
  p <- ggmethylation:::.build_delta_panel(df, 1000, 2000)
  expect_s3_class(p, "ggplot")
})

# Local minimal methylation_data constructor (mirrors
# test-build_read_panel.R::make_test_data; not shared via a helper-*.R file,
# so duplicated here per testthat edition-3 auto-sourcing rules).
make_test_data <- function(reads_df, sites_df,
                            region_start = 1000L, region_end = 2000L,
                            cigar_features = data.frame(
                              read_name  = character(0),
                              type       = character(0),
                              ref_start  = integer(0),
                              ref_end    = integer(0),
                              length     = integer(0),
                              stringsAsFactors = FALSE
                            )) {
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
      cigar_features = cigar_features
    ),
    class = "methylation_data"
  )
}

test_that("plot_methylation adds delta panel for 2 groups when show_delta", {
  reads <- data.frame(
    read_name = c("r1","r2"), start = 1000L, end = 2000L, strand = "+",
    lane = 0L, mean_mod_prob = 0.5, group = c("1","2"),
    clip_side = NA_character_, sa_chrom = NA_character_, stringsAsFactors = FALSE
  )
  sites <- data.frame(
    read_name = rep(c("r1","r2"), each = 6),
    position = rep(seq(1100, 1900, length.out = 6), 2),
    mod_prob = c(rep(0.2,6), rep(0.8,6)), mod_code = "m",
    group = rep(c("1","2"), each = 6), stringsAsFactors = FALSE
  )
  md <- make_test_data(reads, sites); md$group_tag <- "HP"
  p_no  <- ggmethylation::plot_methylation(md, show_delta = FALSE, show_supplementary = FALSE)
  p_yes <- ggmethylation::plot_methylation(md, show_delta = TRUE,  show_supplementary = FALSE)
  expect_gt(length(p_yes$patches$plots), length(p_no$patches$plots))
})

test_that("plot_methylation delta panel breaks over consensus deletions using real group values", {
  # Regression test for the fix where .consensus_deletion_ranges() must be
  # computed against the REAL per-read group values ("1"/"2" from
  # data$reads) *before* the resulting ranges' `group` column is relabeled
  # to the delta line's placeholder identity "delta". A naive (buggy)
  # rewrite that relabels the delta frame's group to "delta" first and then
  # feeds the *unrelabeled* ranges (keyed by "1"/"2") into
  # .insert_deletion_breaks() would compare "delta" == "1"/"2" (always
  # FALSE), silently never masking anything - the delta line would look
  # identical whether or not show_cigar masking applied.
  reads <- data.frame(
    read_name = c("r1", "r2"), start = 1000L, end = 2000L, strand = "+",
    lane = 0L, mean_mod_prob = 0.5, group = c("1", "2"),
    clip_side = NA_character_, sa_chrom = NA_character_, stringsAsFactors = FALSE
  )
  sites <- data.frame(
    read_name = rep(c("r1", "r2"), each = 6),
    position = rep(seq(1100, 1900, length.out = 6), 2),
    mod_prob = c(rep(0.2, 6), rep(0.8, 6)), mod_code = "m",
    group = rep(c("1", "2"), each = 6), stringsAsFactors = FALSE
  )
  # Read r1 (group "1") carries a 100bp consensus deletion at [1400, 1500],
  # well above the default min_indel_size = 50L.
  cigar_features <- data.frame(
    read_name = "r1", type = "D", ref_start = 1400L, ref_end = 1500L,
    length = 100L, stringsAsFactors = FALSE
  )
  md <- make_test_data(reads, sites, cigar_features = cigar_features)
  md$group_tag <- "HP"

  # patchwork::wrap_plots() builds the composite by folding the panel list
  # with `+`; the LAST panel supplied (the delta panel, per the ordering
  # comment above `panels <- c(panels, list(p_delta))`) ends up as the
  # top-level ggplot object itself, with the earlier panels tucked away in
  # `$patches$plots`. So the delta panel's own data is simply `p$data`.
  extract_delta_data <- function(p) p$data

  p_cigar_off <- ggmethylation::plot_methylation(
    md, show_delta = TRUE, show_cigar = FALSE, show_supplementary = FALSE
  )
  p_cigar_on <- ggmethylation::plot_methylation(
    md, show_delta = TRUE, show_cigar = TRUE, show_supplementary = FALSE
  )

  d_off <- extract_delta_data(p_cigar_off)
  d_on  <- extract_delta_data(p_cigar_on)

  in_del_off <- d_off$position >= 1400 & d_off$position <= 1500
  in_del_on  <- d_on$position  >= 1400 & d_on$position  <= 1500

  # Sanity: without cigar masking, the grid actually covers the deletion
  # interval (otherwise the assertions below would be vacuously true).
  expect_true(any(in_del_off))

  # With the fix, show_cigar = TRUE must remove/blank all delta points
  # inside the consensus-deletion interval, and thus emit strictly fewer
  # rows than the unmasked (show_cigar = FALSE) delta data.
  expect_false(any(in_del_on))
  expect_lt(nrow(d_on), nrow(d_off))
})


# --- Delta panel fills follow the group palette ------------------------------

# The delta panel is the last panel handed to patchwork::wrap_plots(), so it
# ends up as the top-level ggplot object (see the note above
# extract_delta_data()). Stripping the patchwork class lets ggplot_build()
# render that panel alone rather than the whole composite.
delta_fills <- function(p) {
  # Drop ONLY the patchwork class. Assigning c("gg", "ggplot") wholesale used to
  # work when ggplot objects were plain S3, but ggplot2 >= 4.0 makes them S7
  # ("ggplot2::ggplot" / "S7_object"); overwriting the class vector destroys S7
  # dispatch and ggplot_build() then fails to find a method.
  class(p) <- setdiff(class(p), "patchwork")
  unique(ggplot2::ggplot_build(p)$data[[1]]$fill)
}

# Two groups whose methylation levels are constant, so the delta has a single
# sign throughout and the panel resolves to exactly one fill colour.
# `higher` names the group given prob 0.8; the other gets 0.2.
make_delta_md <- function(group_names = c("1", "2"), higher = group_names[2L]) {
  probs <- ifelse(group_names == higher, 0.8, 0.2)
  reads <- data.frame(
    read_name = c("r1", "r2"), start = 1000L, end = 2000L, strand = "+",
    lane = 0:1, mean_mod_prob = probs, group = group_names,
    clip_side = NA_character_, sa_chrom = NA_character_,
    stringsAsFactors = FALSE
  )
  sites <- data.frame(
    read_name = rep(c("r1", "r2"), each = 6),
    position  = rep(seq(1100, 1900, length.out = 6), 2),
    mod_prob  = rep(probs, each = 6), mod_code = "m",
    group     = rep(group_names, each = 6),
    stringsAsFactors = FALSE
  )
  md <- make_test_data(reads, sites)
  md$group_tag <- "HP"
  md
}

test_that("delta panel takes its fills from group_colours", {
  cols <- c("1" = "darkgreen", "2" = "purple")

  # Group "2" higher -> delta positive -> group "2"'s colour.
  p_pos <- ggmethylation::plot_methylation(
    make_delta_md(higher = "2"), show_delta = TRUE,
    show_supplementary = FALSE, group_colours = cols
  )
  expect_equal(delta_fills(p_pos), "purple")

  # Group "1" higher -> delta negative -> group "1"'s colour.
  p_neg <- ggmethylation::plot_methylation(
    make_delta_md(higher = "1"), show_delta = TRUE,
    show_supplementary = FALSE, group_colours = cols
  )
  expect_equal(delta_fills(p_neg), "darkgreen")
})

test_that("delta panel falls back to the diverging palette for unnamed groups", {
  # Groups "A"/"B" are not named in the default group_colours, so neither half
  # can be resolved; the panel must fall back rather than emit NA fills.
  p <- ggmethylation::plot_methylation(
    make_delta_md(group_names = c("A", "B"), higher = "B"),
    show_delta = TRUE, show_supplementary = FALSE
  )
  expect_equal(delta_fills(p), ggmethylation:::.DELTA_DIVERGING$pos)
})

test_that("delta panel falls back when only one group is named", {
  # A partial match would leave one half group-coloured and the other not -
  # exactly the half-matching state this feature removes. Require both.
  p <- ggmethylation::plot_methylation(
    make_delta_md(group_names = c("1", "B"), higher = "B"),
    show_delta = TRUE, show_supplementary = FALSE,
    group_colours = c("1" = "darkgreen")
  )
  expect_equal(delta_fills(p), ggmethylation:::.DELTA_DIVERGING$pos)
})

test_that("group_colours still reach the delta panel under cigar masking", {
  # Regression: .compute_group_delta() returns the group names as an attribute,
  # but the show_cigar branch rebuilds delta_df with data.frame(), which drops
  # attributes. Reading the attribute after that rebuild would yield NULL and
  # silently fall back to the diverging palette whenever masking applied.
  md <- make_delta_md(higher = "2")
  md$cigar_features <- data.frame(
    read_name = "r1", type = "D", ref_start = 1400L, ref_end = 1500L,
    length = 100L, stringsAsFactors = FALSE
  )
  p <- ggmethylation::plot_methylation(
    md, show_delta = TRUE, show_cigar = TRUE, show_supplementary = FALSE,
    group_colours = c("1" = "darkgreen", "2" = "purple")
  )
  expect_equal(delta_fills(p), "purple")
})
