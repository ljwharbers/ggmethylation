# .split_reads_on_deletions() and .make_read_polygons() are vectorised. The
# reference implementations below are the earlier per-read loops, kept
# verbatim: the vectorised versions must return identical data frames
# (column order and types, row order, values, row names).

ref_split_reads_on_deletions = function(reads, del_df) {
  if (is.null(del_df) || nrow(del_df) == 0L) return(reads)

  result_list = vector("list", nrow(reads))

  for (i in seq_len(nrow(reads))) {
    rn      = reads$read_name[i]
    r_start = reads$start[i]
    r_end   = reads$end[i]

    rdels = del_df[del_df$read_name == rn, , drop = FALSE]

    if (nrow(rdels) == 0L) {
      result_list[[i]] = reads[i, , drop = FALSE]
      next
    }

    rdels = rdels[order(rdels$ref_start), , drop = FALSE]

    segments = list()
    cur_start = r_start

    for (d in seq_len(nrow(rdels))) {
      d_start = max(rdels$ref_start[d], r_start)
      d_end   = min(rdels$ref_end[d],   r_end)

      if (d_start > cur_start) {
        seg        = reads[i, , drop = FALSE]
        seg$start  = cur_start
        seg$end    = d_start - 1L
        segments   = c(segments, list(seg))
      }
      cur_start = d_end + 1L
    }

    if (cur_start <= r_end) {
      seg        = reads[i, , drop = FALSE]
      seg$start  = cur_start
      seg$end    = r_end
      segments   = c(segments, list(seg))
    }

    result_list[[i]] = if (length(segments) > 0L)
      do.call(rbind, segments)
    else
      reads[i, , drop = FALSE]
  }

  out = do.call(rbind, result_list)
  rownames(out) = NULL
  out
}

ref_make_read_polygons = function(reads, arrow_w, half_height) {
  if (nrow(reads) == 0L) {
    return(data.frame(
      x = numeric(0), y = numeric(0), polygon_id = character(0),
      stringsAsFactors = FALSE
    ))
  }

  poly_list = vector("list", nrow(reads))

  for (i in seq_len(nrow(reads))) {
    r = reads[i, , drop = FALSE]
    s   = r$start
    e   = r$end
    ln  = r$lane
    hh  = half_height
    aw  = arrow_w
    st  = r$strand

    is_first = if ("is_first_segment" %in% names(r)) r$is_first_segment else TRUE
    is_last  = if ("is_last_segment" %in% names(r)) r$is_last_segment else TRUE

    read_len = e - s
    aw = min(aw, read_len)

    draw_right_arrow = (st == "+") && is_last
    draw_left_arrow  = (st == "-") && is_first

    if (draw_right_arrow) {
      xs = c(s, e, e + aw, e, s)
      ys = c(ln - hh, ln - hh, ln, ln + hh, ln + hh)
    } else if (draw_left_arrow) {
      xs = c(s - aw, s, e, e, s)
      ys = c(ln, ln - hh, ln - hh, ln + hh, ln + hh)
    } else {
      xs = c(s, e, e, s)
      ys = c(ln - hh, ln - hh, ln + hh, ln + hh)
    }

    pid = paste0("read_", i)
    n_verts = length(xs)

    verts = r[rep(1L, n_verts), , drop = FALSE]
    verts$x = xs
    verts$y = ys
    verts$polygon_id = pid

    poly_list[[i]] = verts
  }

  out = do.call(rbind, poly_list)
  rownames(out) = NULL
  out
}

# The segment flags exactly as build_read_panel() adds them.
add_segment_flags = function(reads_plot) {
  seg_start = reads_plot$start
  reads_plot$is_first_segment = seg_start == stats::ave(seg_start, reads_plot$read_name, FUN = min)
  reads_plot$is_last_segment  = seg_start == stats::ave(seg_start, reads_plot$read_name, FUN = max)
  reads_plot
}

expect_same_read_geometry = function(reads, dels, arrow_w = 18, half_height = 0.35) {
  split_new = ggmethylation:::.split_reads_on_deletions(reads, dels)
  expect_identical(split_new, ref_split_reads_on_deletions(reads, dels))
  flagged = add_segment_flags(split_new)
  expect_identical(ggmethylation:::.make_read_polygons(flagged, arrow_w, half_height),
                   ref_make_read_polygons(flagged, arrow_w, half_height))
  # without the flags (every segment may carry an arrowhead)
  expect_identical(ggmethylation:::.make_read_polygons(split_new, arrow_w, half_height),
                   ref_make_read_polygons(split_new, arrow_w, half_height))
}

test_that("vectorised read geometry matches the per-read loops on the fixture", {
  for (args in list(list(), list(group_tag = "HP"),
                    list(snv = list(position = 3000, ref = "A", alt = "G")))) {
    md = suppressWarnings(do.call(read_methylation,
                                  c(list(fixture_bam(), "chr1:1000-7000"), args)))
    reads = ggmethylation:::.prepare_reads(md, md$sites)$data$reads
    for (min_size in c(50L, 1L)) {
      expect_same_read_geometry(reads, ggmethylation:::.large_deletions(md$cigar_features, min_size))
    }
    expect_same_read_geometry(reads, NULL)
    expect_same_read_geometry(reads, md$cigar_features[0L, ])
  }
  # a narrow window: reads clipped to the region, deletions crossing its edges
  md = suppressWarnings(read_methylation(fixture_bam(), "chr1:3950-4100", group_tag = "HP"))
  reads = ggmethylation:::.prepare_reads(md, md$sites)$data$reads
  expect_same_read_geometry(reads, ggmethylation:::.large_deletions(md$cigar_features, 1L),
                            arrow_w = 0.45)
})

test_that("vectorised read geometry matches the per-read loops on edge cases", {
  reads = data.frame(
    read_name = paste0("r", 1:7),
    start  = c(1000L, 1000L, 1200L, 1500L, 1600L, 1700L, 1800L),
    end    = c(2000L, 1500L, 1900L, 1505L, 1990L, 1800L, 1802L),
    strand = c("+", "-", "+", "-", "*", "+", "-"),
    lane   = c(1L, 2L, 3L, 2L, 4L, 5L, 6L),
    clip_side = NA_character_, sa_chrom = NA_character_,
    stringsAsFactors = FALSE
  )
  sites = data.frame(read_name = "r1", position = 1100, mod_prob = 0.5,
                     mod_code = "m", stringsAsFactors = FALSE)
  cf = data.frame(
    read_name = c("r1", "r1", "r1", "r2", "r3", "r4", "r5", "r6"),
    type      = "D",
    # r1: two overlapping deletions plus one touching the read end; r2: one at
    # the read start; r3: one past the read end; r4: covers the whole read;
    # r5: unsorted input; r6: shares one base with the read start
    ref_start = c(1400L, 1300L, 1950L, 1000L, 2100L, 1490L, 1700L, 1650L),
    ref_end   = c(1500L, 1450L, 2000L, 1100L, 2200L, 1510L, 1750L, 1700L),
    length    = 100L, stringsAsFactors = FALSE
  )
  md = make_test_data(reads, sites, cigar_features = cf)
  expect_same_read_geometry(md$reads, md$cigar_features)
  expect_same_read_geometry(md$reads, md$cigar_features, arrow_w = 2)

  # double coordinates and lanes; a factor column carried through
  dbl = reads
  dbl$start = as.numeric(dbl$start)
  dbl$lane = as.numeric(dbl$lane)
  dbl$tag = factor(c("a", "b", "a", "b", "a", "b", "a"))
  cf_dbl = cf
  cf_dbl$ref_end = as.numeric(cf_dbl$ref_end)
  expect_same_read_geometry(dbl, cf_dbl)
  expect_same_read_geometry(dbl, cf)
  expect_same_read_geometry(reads, cf_dbl)

  # rectangles only: integer x is kept
  rect = reads
  rect$strand = "*"
  expect_same_read_geometry(rect, cf)
  expect_type(ggmethylation:::.make_read_polygons(rect, 18, 0.35)$x, "integer")

  # deletions of reads that are not in `reads`, and non-default row names
  other = cf
  other$read_name = paste0("x", seq_len(nrow(other)))
  expect_same_read_geometry(reads, other)
  shuffled = reads[c(3, 1, 2, 7, 5, 4, 6), ]
  expect_same_read_geometry(shuffled, cf)

  # a single read and zero reads
  expect_same_read_geometry(reads[1L, ], cf)
  expect_identical(ggmethylation:::.make_read_polygons(reads[0L, ], 18, 0.35),
                   ref_make_read_polygons(reads[0L, ], 18, 0.35))
})
