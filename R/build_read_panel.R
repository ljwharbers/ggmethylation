# Internal function — not exported

# Split reads into non-deletion sub-segments for IGV-style gap display.
# Returns a data.frame with the same columns as `reads` but potentially more
# rows: each row is a contiguous non-deleted stretch of a read.
#
# Per read, with its deletions sorted by ref_start and each clipped to the read
# (d_start = max(ref_start, start), d_end = min(ref_end, end)): a segment runs
# from the read start, or from the previous deletion's d_end + 1, up to
# d_start - 1 when that is not empty, and a last one from the last d_end + 1
# to the read end. A read without any segment is kept whole. Vectorised over
# all deletions at once; segments are collected in typed pieces so that the
# start/end column types are the same as row-binding the segments would give.
.split_reads_on_deletions = function(reads, del_df) {
  if (is.null(del_df) || nrow(del_df) == 0L) return(reads)

  # Deletion rows of every read row, matched by name (each alignment of a read
  # gets all of that read's deletions), in ref_start order within a read.
  per_read = split(seq_len(nrow(del_df)), del_df$read_name)[as.character(reads$read_name)]
  ri = rep(seq_len(nrow(reads)), lengths(per_read))
  di = unlist(per_read, use.names = FALSE)
  o  = order(ri, del_df$ref_start[di])
  ri = ri[o]
  di = di[o]

  r_start = reads$start[ri]
  r_end   = reads$end[ri]
  d_start = pmax(del_df$ref_start[di], r_start)
  d_end   = pmin(del_df$ref_end[di], r_end)
  after   = d_end + 1L
  first   = !duplicated(ri)
  last    = !duplicated(ri, fromLast = TRUE)
  prev    = seq_along(ri) - 1L  # index of the previous deletion (if !first)

  before = r_start
  before[!first] = after[prev[!first]]
  gap      = d_start > before        # a segment before this deletion
  has_tail = last & after <= r_end   # a segment after the read's last deletion
  b1 = which(gap & first)
  b2 = which(gap & !first)
  cl = which(has_tail)
  split_rows = unique(ri[c(b1, b2, cl)])
  whole = setdiff(seq_len(nrow(reads)), split_rows)

  # Order within a read: segments by deletion, then the last one.
  pos = seq_along(ri)
  pieces = list(
    list(i = whole,  k = rep(0, length(whole)),
         start = reads$start[whole], end = reads$end[whole]),
    list(i = ri[b1], k = pos[b1], start = r_start[b1],        end = d_start[b1] - 1L),
    list(i = ri[b2], k = pos[b2], start = after[prev[b2]],    end = d_start[b2] - 1L),
    list(i = ri[cl], k = rep(Inf, length(cl)), start = after[cl], end = r_end[cl])
  )
  pieces = Filter(function(p) length(p$i) > 0L, pieces)
  i   = unlist(lapply(pieces, `[[`, "i"))
  k   = unlist(lapply(pieces, `[[`, "k"))
  out_order = order(i, k)

  out = reads[i[out_order], , drop = FALSE]
  out$start = do.call(c, lapply(pieces, `[[`, "start"))[out_order]
  out$end   = do.call(c, lapply(pieces, `[[`, "end"))[out_order]
  rownames(out) = NULL
  out
}

# Build arrow-tipped polygon vertices for each read segment.
#
# Each read becomes a 5-vertex (arrow) or 4-vertex (rectangle) polygon.
# + strand: arrowhead at right end.  - strand: arrowhead at left end.
# For split reads, only the terminal segment in the read direction gets
# the arrowhead (controlled by is_first_segment / is_last_segment).
#
# @param reads data.frame with columns: start, end, strand, lane,
#   read_name, and optionally is_first_segment, is_last_segment plus any
#   extra columns that should be carried through (group, etc.).
# @param arrow_w Numeric. Width of the arrowhead in genomic coordinates.
# @param half_height Numeric. Half the height of the read rectangle in
#   lane units.
# @return data.frame with columns x, y, polygon_id, plus all original
#   read columns (repeated per vertex).
.make_read_polygons = function(reads, arrow_w, half_height) {
  if (nrow(reads) == 0L) {
    return(data.frame(
      x = numeric(0), y = numeric(0), polygon_id = character(0),
      stringsAsFactors = FALSE
    ))
  }

  n  = nrow(reads)
  s  = reads$start
  e  = reads$end
  ln = reads$lane
  hh = half_height
  is_first = if ("is_first_segment" %in% names(reads)) reads$is_first_segment else rep(TRUE, n)
  is_last  = if ("is_last_segment" %in% names(reads)) reads$is_last_segment else rep(TRUE, n)
  aw = pmin(arrow_w, e - s)

  # Which segments get an arrowhead
  right = (reads$strand == "+" & is_last) %in% TRUE
  left  = !right & (reads$strand == "-" & is_first) %in% TRUE
  arrow = right | left

  # One row of vertices per read: right arrow (5), left arrow (5), or plain
  # rectangle (the first 4).
  lo = ln - hh
  hi = ln + hh
  xs = cbind(s, e, e, s, s)
  ys = cbind(lo, lo, hi, hi, hi)
  xs[right, ] = cbind(s, e, e + aw, e, s)[right, ]
  ys[right, ] = cbind(lo, lo, ln, hi, hi)[right, ]
  xs[left, ]  = cbind(s - aw, s, e, e, s)[left, ]
  ys[left, ]  = cbind(ln, lo, lo, hi, hi)[left, ]
  n_verts = ifelse(arrow, 5L, 4L)
  keep = t(col(xs) <= n_verts)

  # Vertex order is read by read; types as c() of each read's vertices gives.
  x = t(xs)[keep]
  storage.mode(x) = typeof(c(s[0L], e[0L], if (any(arrow)) aw[0L]))
  y = t(ys)[keep]
  storage.mode(y) = typeof(c(ln[0L], lo[0L]))

  idx = rep(seq_len(n), n_verts)
  out = reads[idx, , drop = FALSE]
  out$x = x
  out$y = y
  out$polygon_id = paste0("read_", idx)
  rownames(out) = NULL
  out
}

# Build colored indicator polygons at breakpoint ends of supplementary reads.
#
# For each read with a non-NA sa_chrom, a small indicator is placed at the
# clipped end: a triangle (matching the arrowhead) when the clip side
# coincides with the arrowhead direction, or a small rectangle tab otherwise.
# Both are drawn outward from the read end, never over the bar itself.
#
# `reads` may be the deletion-split segments (see .split_reads_on_deletions()):
# the left marker then goes only on a read's first segment and the right
# marker only on its last, so an internal deletion edge never gets one.
#
# @param reads data.frame with SA-annotated reads (sa_chrom, clip_side columns,
#   optionally is_first_segment / is_last_segment).
# @param arrow_w Numeric.
# @param half_height Numeric.
# @return data.frame with x, y, polygon_id, sa_chrom, lane columns,
#   or a 0-row data.frame if no SA reads.
.make_sa_overlay_polygons = function(reads, arrow_w, half_height) {
  sa_reads = reads[!is.na(reads$sa_chrom), , drop = FALSE]
  if (nrow(sa_reads) == 0L) {
    return(data.frame(
      x = numeric(0), y = numeric(0), polygon_id = character(0),
      sa_chrom = character(0), lane = numeric(0),
      vcf_validated = logical(0),
      stringsAsFactors = FALSE
    ))
  }

  poly_list = vector("list", 2L * nrow(sa_reads))
  pid_counter = 0L

  for (i in seq_len(nrow(sa_reads))) {
    r  = sa_reads[i, , drop = FALSE]
    s  = r$start
    e  = r$end
    ln = r$lane
    hh = half_height
    st = r$strand
    cs = r$clip_side

    read_len = e - s
    aw = min(arrow_w, read_len)

    # Determine which sides to overlay
    sides = character(0)
    if (is.na(cs)) {
      if (st == "+") sides = "right" else sides = "left"
    } else if (cs == "both") {
      sides = c("left", "right")
    } else {
      sides = cs
    }
    if ("is_first_segment" %in% names(r) && !isTRUE(r$is_first_segment)) {
      sides = setdiff(sides, "left")
    }
    if ("is_last_segment" %in% names(r) && !isTRUE(r$is_last_segment)) {
      sides = setdiff(sides, "right")
    }

    # Arrowhead side depends on strand
    arrow_side = if (st == "+") "right" else "left"

    for (side in sides) {
      pid_counter = pid_counter + 1L
      pid = paste0("sa_", pid_counter)

      if (side == arrow_side) {
        # Triangle matching the arrowhead shape
        if (side == "right") {
          xs = c(e, e + aw, e)
          ys = c(ln - hh, ln, ln + hh)
        } else {
          xs = c(s - aw, s, s)
          ys = c(ln, ln - hh, ln + hh)
        }
      } else {
        # Small rectangle tab on the non-arrowhead side, extending outward
        if (side == "left") {
          xs = c(s - aw, s, s, s - aw)
          ys = c(ln - hh, ln - hh, ln + hh, ln + hh)
        } else {
          xs = c(e, e + aw, e + aw, e)
          ys = c(ln - hh, ln - hh, ln + hh, ln + hh)
        }
      }

      poly_list[[pid_counter]] = data.frame(
        x             = xs,
        y             = ys,
        polygon_id    = pid,
        sa_chrom      = r$sa_chrom,
        lane          = ln,
        vcf_validated = if ("vcf_validated" %in% names(r)) isTRUE(r$vcf_validated) else FALSE,
        stringsAsFactors = FALSE
      )
    }
  }

  out = do.call(rbind, poly_list[seq_len(pid_counter)])
  rownames(out) = NULL
  out
}

#' Build a ggplot2 read-level methylation panel
#'
#' Constructs the read panel (horizontal read bars + modification probability
#' lines) for a single `methylation_data` object that has already been sorted
#' and lane-packed. Used by [plot_methylation()] for every sample it draws.
#'
#' Read bars are filled by group when the data is grouped, by strand when
#' `colour_strand = TRUE`, and in a fixed grey otherwise.
#'
#' With several modification codes in `data$sites`, each read bar is split
#' into horizontal strips, one per code in the order of `data$mod_code` (top
#' to bottom), so calls for different codes at the same position do not
#' overplot; the y-axis title lists the order.
#'
#' @param data A `methylation_data` object. `$reads` must already contain a
#'   `lane` column (from `pack_reads()`) and a `mean_mod_prob` column.
#'   `$cigar_features` supplies the indels drawn when `show_cigar = TRUE`.
#' @param separator_lanes Numeric vector of y-positions where horizontal
#'   dashed separator lines should be drawn (between groups). Pass
#'   `numeric(0)` when no separators are needed.
#' @param region_start Integer. Left boundary of the x-axis.
#' @param region_end Integer. Right boundary of the x-axis.
#' @param colour_low Colour for low modification probability.
#' @param colour_high Colour for high modification probability.
#' @param colour_ambiguous Colour for ambiguous calls. Only used when
#'   `call_threshold` and `call_ambiguous` are set.
#' @param line_width Linewidth of modification site markers.
#' @param colour_strand Logical. Colour read bars by strand when ungrouped.
#' @param strand_colours Named character vector with `"+"` and `"-"` entries.
#' @param group_colours Colours per group as accepted by [plot_methylation()]
#'   (named, unnamed or NULL); resolved over the group levels by the same
#'   rule as the smooth and delta panels.
#' @param show_x_axis Logical. When `FALSE` (default), x-axis text and ticks
#'   are hidden. Set to `TRUE` for the bottom-most read panel.
#' @param variant_overlay A list returned by [build_variant_overlay()], or `NULL`.
#'   When non-NULL, SNV point marks, SV spans, and BND position markers are drawn
#'   on the read panel.
#' @param show_cigar Logical. When `TRUE`, structural variants from CIGAR
#'   strings are overlaid on reads. Default `FALSE`.
#' @param min_indel_size Integer. Minimum size (in bp) for insertions and
#'   deletions to be displayed. Features smaller than this threshold are
#'   suppressed. Default `50`.
#' @param show_supplementary Logical. Draw supplementary-alignment indicators.
#' @param call_threshold,call_ambiguous Binary-call settings (`NULL` threshold
#'   for continuous colouring); see [plot_methylation()].
#' @param legend_levels `NULL`, or a list with `group` (group levels),
#'   `group_na` (logical: show an NA key) and `sa` (SA partner chromosomes).
#'   [plot_methylation()] passes the union over all samples so every read
#'   panel has identical legends and colours, which patchwork then merges.
#'   `NULL` uses the levels present in `data`.
#'
#' @return A [ggplot2::ggplot] object.
#'
#' @keywords internal
build_read_panel = function(data,
                             separator_lanes,
                             region_start,
                             region_end,
                             colour_low,
                             colour_high,
                             colour_ambiguous = .CALL_AMBIGUOUS_DEFAULT,
                             line_width,
                             colour_strand,
                             strand_colours,
                             group_colours,
                             show_x_axis        = FALSE,
                             variant_overlay    = NULL,
                             show_cigar         = FALSE,
                             min_indel_size     = 50L,
                             show_supplementary = FALSE,
                             call_threshold     = NULL,
                             call_ambiguous     = NULL,
                             legend_levels      = NULL) {
  if (is.null(legend_levels)) {
    legend_levels = list(
      group    = .group_levels(data$reads$group),
      group_na = any(is.na(data$reads$group)),
      sa       = .group_levels(data$reads$sa_chrom)
    )
  }
  cigar_features = if (isTRUE(show_cigar)) data$cigar_features
  has_cigar = !is.null(cigar_features) && nrow(cigar_features) > 0L

  # When show_cigar is TRUE, split reads on large deletions so the thick read
  # bar has IGV-style gaps instead of running through deletion regions.
  reads_plot = data$reads
  if (has_cigar) {
    reads_plot = .split_reads_on_deletions(
      data$reads, .large_deletions(cigar_features, min_indel_size)
    )
  }

  # Only the outermost segments of a split read carry an arrowhead: flag the
  # first and last segment of each read in genomic order (unsplit reads get
  # both flags).
  seg_start = reads_plot$start
  reads_plot$is_first_segment = seg_start == stats::ave(seg_start, reads_plot$read_name, FUN = min)
  reads_plot$is_last_segment  = seg_start == stats::ave(seg_start, reads_plot$read_name, FUN = max)

  # Merge vcf_validated flag from variant_overlay$sa_reads into reads_plot
  # Only needed when supplementary overlays will actually be drawn.
  if (isTRUE(show_supplementary) &&
      !is.null(variant_overlay) && !is.null(variant_overlay$sa_reads)) {
    validated_flag = variant_overlay$sa_reads[, c("read_name", "vcf_validated"), drop = FALSE]
    reads_plot = merge(reads_plot, validated_flag, by = "read_name", all.x = TRUE)
    reads_plot$vcf_validated[is.na(reads_plot$vcf_validated)] = FALSE
  }

  # Merge lane info and segment extents into sites.  Using reads_plot (the
  # deletion-split version when show_cigar=TRUE) means dots that fall in
  # large-deletion gaps are excluded along with any remaining soft-clipped
  # positions.  When show_cigar=FALSE, reads_plot == data$reads so behaviour
  # is unchanged.
  sites_plot = merge(
    data$sites,
    reads_plot[, c("read_name", "lane", "start", "end"), drop = FALSE],
    by = "read_name"
  )
  # Keep only dots that lie within a visible read segment
  sites_plot = sites_plot[
    sites_plot$position >= sites_plot$start &
      sites_plot$position <= sites_plot$end,
    , drop = FALSE
  ]
  sites_plot$start = NULL
  sites_plot$end   = NULL

  # --- Arrow geometry parameters ---
  arrow_w     = (region_end - region_start) * 0.003
  half_height = 0.35

  # Sites next to a clipped end are kept: the SA marker is drawn outward from
  # the read end (see .make_sa_overlay_polygons()), so it never covers them.

  # --- Read bars: filled by group, by strand, or a fixed grey ---
  grouped  = !is.null(data$group_tag)
  fill_var = if (grouped) "group" else if (isTRUE(colour_strand)) "strand"
  read_polys = .make_read_polygons(reads_plot, arrow_w, half_height)
  bar_aes = list(x = quote(.data$x), y = quote(.data$y), group = quote(.data$polygon_id))

  p = ggplot2::ggplot()
  if (is.null(fill_var)) {
    p = p + ggplot2::geom_polygon(data = read_polys, do.call(ggplot2::aes, bar_aes),
                                   fill = "#B0BEC5", colour = NA)
  } else {
    bar_aes$fill = .data_col(fill_var)
    fill_scale = if (!grouped) {
      ggplot2::scale_fill_manual(values = strand_colours, name = "Strand")
    } else {
      colours = .resolve_group_colours(group_colours,
                                       c(legend_levels$group, data$reads$group))
      ggplot2::scale_fill_manual(
        values = colours,
        limits = c(names(colours), if (isTRUE(legend_levels$group_na)) NA),
        na.value = .GROUP_NA_COLOUR, name = "Group"
      )
    }
    if (grouped) {
      # Every legend key must be drawn (also for groups absent from this
      # sample) so that all samples' legends match and patchwork merges them.
      read_polys = .add_legend_placeholders(
        read_polys, "group",
        c(legend_levels$group, if (isTRUE(legend_levels$group_na)) NA)
      )
    }
    p = p + ggplot2::geom_polygon(data = read_polys, do.call(ggplot2::aes, bar_aes),
                                   colour = NA) + fill_scale
  }

  if (isTRUE(show_supplementary)) {
    # A bar fill scale is already on the plot unless the bars are plain grey.
    p = .add_sa_overlay(p, reads_plot, if (grouped) arrow_w * 1.5 else arrow_w,
                         half_height, variant_overlay,
                         needs_new_scale = !is.null(fill_var),
                         sa_levels = legend_levels$sa)
  }
  mod_codes = .present_codes(data$mod_code, sites_plot)
  p = .add_mod_prob_segments(p, sites_plot, half_height, line_width,
                              colour_low, colour_high, colour_ambiguous,
                              call_threshold, call_ambiguous,
                              mod_codes = mod_codes)

  if (length(separator_lanes) > 0L) {
    p = p +
      ggplot2::geom_hline(
        yintercept = separator_lanes,
        linetype = "dashed", colour = "grey40", linewidth = 0.4
      )
  }

  # SNV indicator line
  if (!is.null(data$snv_position)) {
    p = p +
      ggplot2::geom_vline(
        xintercept = data$snv_position,
        linetype = "dashed", colour = "black", linewidth = 0.5
      )
  }

  # CIGAR structural variant overlays
  if (has_cigar) {
    cf = merge(
      cigar_features,
      data$reads[, c("read_name", "lane", "start", "end"), drop = FALSE],
      by = "read_name"
    )
    cf = cf[cf$type %in% c("I", "D") & cf$length >= min_indel_size, , drop = FALSE]

    # Deletions: draw black line segments
    del_df = cf[cf$type == "D", , drop = FALSE]
    if (nrow(del_df) > 0L) {
      p = p + ggplot2::geom_segment(
        data = del_df,
        ggplot2::aes(
          x = .data$ref_start, xend = .data$ref_end,
          y = .data$lane, yend = .data$lane
        ),
        linewidth = 0.3, colour = "black", linetype = "solid",
        inherit.aes = FALSE
      )
    }

    # Insertions: draw purple I-beam markers spanning the read bar
    ins_df = cf[cf$type == "I", , drop = FALSE]
    if (nrow(ins_df) > 0L) {
      p = p + ggplot2::geom_errorbar(
        data = ins_df,
        ggplot2::aes(
          x    = .data$ref_start,
          ymin = .data$lane - half_height,
          ymax = .data$lane + half_height
        ),
        colour    = "#7B1FA2",
        linewidth = 0.4,
        width     = arrow_w * 2,
        inherit.aes = FALSE
      )
    }
  }

  # --- Variant overlay layers ---
  if (!is.null(variant_overlay)) {
    for (lyr in c(variant_overlay$snv, variant_overlay$sv, variant_overlay$bnd)) {
      p = p + lyr
    }
  }

  p = p +
    ggplot2::scale_y_reverse() +
    ggplot2::coord_cartesian(xlim = c(region_start, region_end)) +
    theme_ggmethylation() +
    ggplot2::theme(
      axis.text.y        = ggplot2::element_blank(),
      axis.ticks.y       = ggplot2::element_blank(),
      axis.title.y       = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank()
    ) +
    ggplot2::labs(x = NULL)
  if (length(mod_codes) > 1L) {
    # Name the per-code strips (top to bottom) on the otherwise blank y axis.
    p = p +
      ggplot2::labs(y = paste("Strips (top to bottom):", paste(mod_codes, collapse = ", "))) +
      ggplot2::theme(axis.title.y = ggplot2::element_text(size = ggplot2::rel(0.7),
                                                           angle = 90))
  }

  if (!show_x_axis) {
    p = p +
      ggplot2::theme(
        axis.text.x  = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank()
      )
  }

  p
}

# Add a zero-area polygon (three identical vertices on an existing vertex, so
# nothing is drawn and no scale range changes) for each of `levels` missing
# from `polys[[col]]`. ggplot2 (>= 3.5) only draws legend keys for values
# present in a layer's data; the placeholders make every panel draw the same
# keys, so identical legends from several samples can be merged by patchwork.
# (show.legend = TRUE would do the same but breaks keys of layers that
# ggnewscale later renames.)
.add_legend_placeholders = function(polys, col, levels) {
  if (nrow(polys) == 0L || length(levels) == 0L) return(polys)
  present = unique(polys[[col]])
  missing = levels[!levels %in% present]
  if (length(missing) == 0L) return(polys)
  extra = polys[rep(1L, 3L * length(missing)), , drop = FALSE]
  extra[[col]] = rep(missing, each = 3L)
  extra$polygon_id = rep(paste0(".legend_", seq_along(missing)), each = 3L)
  if ("vcf_validated" %in% names(extra)) extra$vcf_validated = FALSE
  out = rbind(polys, extra)
  rownames(out) = NULL
  out
}

# Add SA supplementary-alignment overlay polygons to an existing ggplot.
#
# When `needs_new_scale = TRUE` (grouped / strand-coloured branches), a
# `ggnewscale::new_scale_fill()` call is prepended so the SA fill scale does not
# replace the read-bar fill scale already on the plot.  The plain branch omits
# it because the read bars are drawn with a fixed colour, not a scale.
#
# `sa_levels` fixes the legend keys and colours (hue palette over the sorted
# levels), so the same partner chromosome has the same colour in every sample.
.add_sa_overlay = function(p, reads_plot, arrow_w, half_height,
                             variant_overlay,
                             needs_new_scale = TRUE,
                             sa_levels = NULL) {
  if (!("sa_chrom" %in% names(reads_plot))) return(p)

  sa_polys = .make_sa_overlay_polygons(reads_plot, arrow_w, half_height)
  if (nrow(sa_polys) == 0L) return(p)
  sa_levels = .group_levels(c(sa_levels, sa_polys$sa_chrom))
  sa_colours = stats::setNames(scales::hue_pal()(length(sa_levels)), sa_levels)
  sa_polys = .add_legend_placeholders(sa_polys, "sa_chrom", sa_levels)

  if (needs_new_scale) p = p + ggnewscale::new_scale_fill()

  p = p +
    ggplot2::geom_polygon(
      data = sa_polys,
      ggplot2::aes(
        x     = .data$x,
        y     = .data$y,
        group = .data$polygon_id,
        fill  = .data$sa_chrom
      ),
      colour = NA,
      inherit.aes = FALSE
    ) +
    ggplot2::scale_fill_manual(values = sa_colours, limits = sa_levels,
                               name = "SA partner")

  # VCF-validated SA border
  if (!is.null(variant_overlay) &&
      "vcf_validated" %in% names(sa_polys) &&
      any(sa_polys$vcf_validated, na.rm = TRUE)) {
    p = p + ggplot2::geom_polygon(
      data = sa_polys[sa_polys$vcf_validated %in% TRUE, , drop = FALSE],
      ggplot2::aes(x = .data$x, y = .data$y, group = .data$polygon_id),
      fill = NA, colour = "#880E4F", linewidth = 0.6, inherit.aes = FALSE
    )
  }

  p
}

# Classify continuous modification probabilities into discrete calls.
# ambiguous: NULL for a hard threshold, or a numeric half-width; probs within
# (threshold - w, threshold + w), excluding the threshold value itself, are
# labelled "ambiguous". A prob exactly equal to threshold is always
# "methylated" (the >= threshold-only case); everything else at/above
# threshold and outside the band is also "methylated".
.classify_calls = function(probs, threshold = 0.5, ambiguous = NULL) {
  cls = ifelse(probs >= threshold, "methylated", "unmethylated")
  if (!is.null(ambiguous) && ambiguous > 0) {
    band = abs(probs - threshold) < ambiguous & probs != threshold
    cls[band] = "ambiguous"
  }
  cls
}

# Add mod-prob segment layer and colour scale to an existing ggplot.
#
# With a non-NULL `call_threshold`, sites are classified into discrete calls
# (methylated/unmethylated/ambiguous) via .classify_calls() and coloured with
# a manual discrete scale instead of the continuous gradient.
#
# With more than one entry in `mod_codes`, the bar height is divided into one
# horizontal strip per code (first code on top: the y scale is reversed), so
# calls of different codes at the same read/position are both visible.
.add_mod_prob_segments = function(p, sites_plot, half_height, line_width,
                                    colour_low, colour_high,
                                    colour_ambiguous = .CALL_AMBIGUOUS_DEFAULT,
                                    call_threshold = NULL,
                                    call_ambiguous = NULL,
                                    mod_codes = NULL) {
  sites_plot$.y0 = sites_plot$lane - half_height
  sites_plot$.y1 = sites_plot$lane + half_height
  if (length(mod_codes) > 1L && nrow(sites_plot) > 0L) {
    k = match(sites_plot$mod_code, mod_codes)
    strip = 2 * half_height / length(mod_codes)
    sites_plot$.y0 = sites_plot$lane - half_height + (k - 1L) * strip
    sites_plot$.y1 = sites_plot$.y0 + strip
  }
  if (!is.null(call_threshold) && nrow(sites_plot) > 0L) {
    sites_plot$.call = .classify_calls(sites_plot$mod_prob,
                                        call_threshold, call_ambiguous)
    vals = c(unmethylated = colour_low, methylated = colour_high,
              ambiguous = colour_ambiguous)
    return(
      p +
        ggplot2::geom_segment(
          data = sites_plot,
          ggplot2::aes(
            x = .data$position, xend = .data$position,
            y = .data$.y0, yend = .data$.y1,
            colour = .data$.call
          ),
          linewidth = line_width
        ) +
        ggplot2::scale_colour_manual(
          values = vals, name = "Call",
          breaks = c("unmethylated", "methylated", "ambiguous")
        )
    )
  }
  p +
    ggplot2::geom_segment(
      data = sites_plot,
      ggplot2::aes(
        x      = .data$position,
        xend   = .data$position,
        y      = .data$.y0,
        yend   = .data$.y1,
        colour = .data$mod_prob
      ),
      linewidth = line_width
    ) +
    ggplot2::scale_colour_gradient(
      low = colour_low, high = colour_high,
      limits = c(0, 1),
      name = "Modification\nprobability"
    )
}
