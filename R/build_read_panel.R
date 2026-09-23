# Internal function — not exported

# Split reads into non-deletion sub-segments for IGV-style gap display.
# Returns a data.frame with the same columns as `reads` but potentially more
# rows: each row is a contiguous non-deleted stretch of a read.
.split_reads_on_deletions <- function(reads, del_df) {
  if (is.null(del_df) || nrow(del_df) == 0L) return(reads)

  result_list <- vector("list", nrow(reads))

  for (i in seq_len(nrow(reads))) {
    rn      <- reads$read_name[i]
    r_start <- reads$start[i]
    r_end   <- reads$end[i]

    rdels <- del_df[del_df$read_name == rn, , drop = FALSE]

    if (nrow(rdels) == 0L) {
      result_list[[i]] <- reads[i, , drop = FALSE]
      next
    }

    rdels <- rdels[order(rdels$ref_start), , drop = FALSE]

    segments <- list()
    cur_start <- r_start

    for (d in seq_len(nrow(rdels))) {
      d_start <- max(rdels$ref_start[d], r_start)
      d_end   <- min(rdels$ref_end[d],   r_end)

      if (d_start > cur_start) {
        seg        <- reads[i, , drop = FALSE]
        seg$start  <- cur_start
        seg$end    <- d_start - 1L
        segments   <- c(segments, list(seg))
      }
      cur_start <- d_end + 1L
    }

    if (cur_start <= r_end) {
      seg        <- reads[i, , drop = FALSE]
      seg$start  <- cur_start
      seg$end    <- r_end
      segments   <- c(segments, list(seg))
    }

    result_list[[i]] <- if (length(segments) > 0L)
      do.call(rbind, segments)
    else
      reads[i, , drop = FALSE]
  }

  out <- do.call(rbind, result_list)
  rownames(out) <- NULL
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
.make_read_polygons <- function(reads, arrow_w, half_height) {
  if (nrow(reads) == 0L) {
    return(data.frame(
      x = numeric(0), y = numeric(0), polygon_id = character(0),
      stringsAsFactors = FALSE
    ))
  }

  poly_list <- vector("list", nrow(reads))

  for (i in seq_len(nrow(reads))) {
    r <- reads[i, , drop = FALSE]
    s   <- r$start
    e   <- r$end
    ln  <- r$lane
    hh  <- half_height
    aw  <- arrow_w
    st  <- r$strand

    is_first <- if ("is_first_segment" %in% names(r)) r$is_first_segment else TRUE
    is_last  <- if ("is_last_segment" %in% names(r)) r$is_last_segment else TRUE

    read_len <- e - s
    aw <- min(aw, read_len)

    # Determine whether this segment gets an arrowhead
    draw_right_arrow <- (st == "+") && is_last
    draw_left_arrow  <- (st == "-") && is_first

    if (draw_right_arrow) {
      # Arrow pointing right: 5 vertices
      xs <- c(s, e, e + aw, e, s)
      ys <- c(ln - hh, ln - hh, ln, ln + hh, ln + hh)
    } else if (draw_left_arrow) {
      # Arrow pointing left: 5 vertices
      xs <- c(s - aw, s, e, e, s)
      ys <- c(ln, ln - hh, ln - hh, ln + hh, ln + hh)
    } else {
      # Plain rectangle: 4 vertices
      xs <- c(s, e, e, s)
      ys <- c(ln - hh, ln - hh, ln + hh, ln + hh)
    }

    pid <- paste0("read_", i)
    n_verts <- length(xs)

    # Replicate the row data for each vertex
    verts <- r[rep(1L, n_verts), , drop = FALSE]
    verts$x <- xs
    verts$y <- ys
    verts$polygon_id <- pid

    poly_list[[i]] <- verts
  }

  out <- do.call(rbind, poly_list)
  rownames(out) <- NULL
  out
}

# Build colored indicator polygons at breakpoint ends of supplementary reads.
#
# For each read with a non-NA sa_chrom, a small indicator is placed at the
# clipped end: a triangle (matching the arrowhead) when the clip side
# coincides with the arrowhead direction, or a small rectangle tab otherwise.
#
# @param reads data.frame with SA-annotated reads (sa_chrom, clip_side columns).
# @param arrow_w Numeric.
# @param half_height Numeric.
# @return data.frame with x, y, polygon_id, sa_chrom, lane columns,
#   or a 0-row data.frame if no SA reads.
.make_sa_overlay_polygons <- function(reads, arrow_w, half_height) {
  sa_reads <- reads[!is.na(reads$sa_chrom), , drop = FALSE]
  if (nrow(sa_reads) == 0L) {
    return(data.frame(
      x = numeric(0), y = numeric(0), polygon_id = character(0),
      sa_chrom = character(0), lane = numeric(0),
      vcf_validated = logical(0),
      stringsAsFactors = FALSE
    ))
  }

  poly_list <- vector("list", 2L * nrow(sa_reads))
  pid_counter <- 0L

  for (i in seq_len(nrow(sa_reads))) {
    r  <- sa_reads[i, , drop = FALSE]
    s  <- r$start
    e  <- r$end
    ln <- r$lane
    hh <- half_height
    st <- r$strand
    cs <- r$clip_side

    read_len <- e - s
    aw <- min(arrow_w, read_len)

    # Determine which sides to overlay
    sides <- character(0)
    if (is.na(cs)) {
      if (st == "+") sides <- "right" else sides <- "left"
    } else if (cs == "both") {
      sides <- c("left", "right")
    } else {
      sides <- cs
    }

    # Arrowhead side depends on strand
    arrow_side <- if (st == "+") "right" else "left"

    for (side in sides) {
      pid_counter <- pid_counter + 1L
      pid <- paste0("sa_", pid_counter)

      if (side == arrow_side) {
        # Triangle matching the arrowhead shape
        if (side == "right") {
          xs <- c(e, e + aw, e)
          ys <- c(ln - hh, ln, ln + hh)
        } else {
          xs <- c(s - aw, s, s)
          ys <- c(ln, ln - hh, ln + hh)
        }
      } else {
        # Small rectangle tab on the non-arrowhead side, extending outward
        if (side == "left") {
          xs <- c(s - aw, s, s, s - aw)
          ys <- c(ln - hh, ln - hh, ln + hh, ln + hh)
        } else {
          xs <- c(e, e + aw, e + aw, e)
          ys <- c(ln - hh, ln - hh, ln + hh, ln + hh)
        }
      }

      poly_list[[pid_counter]] <- data.frame(
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

  out <- do.call(rbind, poly_list[seq_len(pid_counter)])
  rownames(out) <- NULL
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
#' @param group_colours Named character vector of colours per group, or NULL.
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
#'
#' @return A [ggplot2::ggplot] object.
#'
#' @keywords internal
build_read_panel <- function(data,
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
                             call_ambiguous     = NULL) {
  cigar_features = if (isTRUE(show_cigar)) data$cigar_features
  has_cigar = !is.null(cigar_features) && nrow(cigar_features) > 0L

  # When show_cigar is TRUE, split reads on large deletions so the thick read
  # bar has IGV-style gaps instead of running through deletion regions.
  reads_plot <- data$reads
  if (has_cigar) {
    reads_plot <- .split_reads_on_deletions(
      data$reads, .large_deletions(cigar_features, min_indel_size)
    )
  }

  # Only the outermost segments of a split read carry an arrowhead: flag the
  # first and last segment of each read in genomic order (unsplit reads get
  # both flags).
  seg_start <- reads_plot$start
  reads_plot$is_first_segment <- seg_start == stats::ave(seg_start, reads_plot$read_name, FUN = min)
  reads_plot$is_last_segment  <- seg_start == stats::ave(seg_start, reads_plot$read_name, FUN = max)

  # Merge vcf_validated flag from variant_overlay$sa_reads into reads_plot
  # Only needed when supplementary overlays will actually be drawn.
  if (isTRUE(show_supplementary) &&
      !is.null(variant_overlay) && !is.null(variant_overlay$sa_reads)) {
    validated_flag <- variant_overlay$sa_reads[, c("read_name", "vcf_validated"), drop = FALSE]
    reads_plot <- merge(reads_plot, validated_flag, by = "read_name", all.x = TRUE)
    reads_plot$vcf_validated[is.na(reads_plot$vcf_validated)] <- FALSE
  }

  # Merge lane info and segment extents into sites.  Using reads_plot (the
  # deletion-split version when show_cigar=TRUE) means dots that fall in
  # large-deletion gaps are excluded along with any remaining soft-clipped
  # positions.  When show_cigar=FALSE, reads_plot == data$reads so behaviour
  # is unchanged.
  sites_plot <- merge(
    data$sites,
    reads_plot[, c("read_name", "lane", "start", "end"), drop = FALSE],
    by = "read_name"
  )
  # Keep only dots that lie within a visible read segment
  sites_plot <- sites_plot[
    sites_plot$position >= sites_plot$start &
      sites_plot$position <= sites_plot$end,
    , drop = FALSE
  ]
  sites_plot$start <- NULL
  sites_plot$end   <- NULL

  # --- Arrow geometry parameters ---
  arrow_w     <- (region_end - region_start) * 0.003
  half_height <- 0.35

  # --- Suppress modification dots within SA indicator regions ---
  if (isTRUE(show_supplementary) &&
      "sa_chrom"  %in% names(data$reads) &&
      "clip_side" %in% names(data$reads)) {
    sa_reads <- data$reads[!is.na(data$reads$sa_chrom) & !is.na(data$reads$clip_side), , drop = FALSE]
    keep <- rep(TRUE, nrow(sites_plot))
    for (i in seq_len(nrow(sa_reads))) {
      r      <- sa_reads[i, ]
      sa_ext <- min(arrow_w, r$end - r$start)
      idx    <- sites_plot$read_name == r$read_name
      if (r$clip_side %in% c("left",  "both")) keep <- keep & !(idx & sites_plot$position <= r$start + sa_ext)
      if (r$clip_side %in% c("right", "both")) keep <- keep & !(idx & sites_plot$position >= r$end - sa_ext)
    }
    sites_plot <- sites_plot[keep, , drop = FALSE]
  }

  # --- Read bars: filled by group, by strand, or a fixed grey ---
  grouped  <- !is.null(data$group_tag)
  fill_var <- if (grouped) "group" else if (isTRUE(colour_strand)) "strand"
  read_polys <- .make_read_polygons(reads_plot, arrow_w, half_height)
  bar_aes <- list(x = quote(.data$x), y = quote(.data$y), group = quote(.data$polygon_id))

  p <- ggplot2::ggplot()
  if (is.null(fill_var)) {
    p <- p + ggplot2::geom_polygon(data = read_polys, do.call(ggplot2::aes, bar_aes),
                                   fill = "#B0BEC5", colour = NA)
  } else {
    bar_aes$fill <- .data_col(fill_var)
    fill_scale <- if (!grouped) {
      ggplot2::scale_fill_manual(values = strand_colours, name = "Strand")
    } else if (!is.null(group_colours)) {
      ggplot2::scale_fill_manual(values = group_colours, na.value = "grey50", name = "Group")
    } else {
      ggplot2::scale_fill_discrete(name = "Group")
    }
    p <- p + ggplot2::geom_polygon(data = read_polys, do.call(ggplot2::aes, bar_aes),
                                   colour = NA) + fill_scale
  }

  if (isTRUE(show_supplementary)) {
    # A bar fill scale is already on the plot unless the bars are plain grey.
    p <- .add_sa_overlay(p, reads_plot, if (grouped) arrow_w * 1.5 else arrow_w,
                         half_height, variant_overlay,
                         needs_new_scale = !is.null(fill_var))
  }
  p <- .add_mod_prob_segments(p, sites_plot, half_height, line_width,
                              colour_low, colour_high, colour_ambiguous,
                              call_threshold, call_ambiguous)

  if (length(separator_lanes) > 0L) {
    p <- p +
      ggplot2::geom_hline(
        yintercept = separator_lanes,
        linetype = "dashed", colour = "grey40", linewidth = 0.4
      )
  }

  # SNV indicator line
  if (!is.null(data$snv_position)) {
    p <- p +
      ggplot2::geom_vline(
        xintercept = data$snv_position,
        linetype = "dashed", colour = "black", linewidth = 0.5
      )
  }

  # CIGAR structural variant overlays
  if (has_cigar) {
    cf <- merge(
      cigar_features,
      data$reads[, c("read_name", "lane", "start", "end"), drop = FALSE],
      by = "read_name"
    )
    cf <- cf[cf$type %in% c("I", "D") & cf$length >= min_indel_size, , drop = FALSE]

    # Deletions: draw black line segments
    del_df <- cf[cf$type == "D", , drop = FALSE]
    if (nrow(del_df) > 0L) {
      p <- p + ggplot2::geom_segment(
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
    ins_df <- cf[cf$type == "I", , drop = FALSE]
    if (nrow(ins_df) > 0L) {
      p <- p + ggplot2::geom_errorbar(
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
      p <- p + lyr
    }
  }

  p <- p +
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

  if (!show_x_axis) {
    p <- p +
      ggplot2::theme(
        axis.text.x  = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank()
      )
  }

  p
}

# Add SA supplementary-alignment overlay polygons to an existing ggplot.
#
# When `needs_new_scale = TRUE` (grouped / strand-coloured branches), a
# `ggnewscale::new_scale_fill()` call is prepended so the SA fill scale does not
# replace the read-bar fill scale already on the plot.  The plain branch omits
# it because the read bars are drawn with a fixed colour, not a scale.
.add_sa_overlay <- function(p, reads_plot, arrow_w, half_height,
                             variant_overlay,
                             needs_new_scale = TRUE) {
  if (!("sa_chrom" %in% names(reads_plot))) return(p)

  sa_polys <- .make_sa_overlay_polygons(reads_plot, arrow_w, half_height)
  if (nrow(sa_polys) == 0L) return(p)

  if (needs_new_scale) p <- p + ggnewscale::new_scale_fill()

  p <- p +
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
    ggplot2::scale_fill_hue(name = "SA partner")

  # VCF-validated SA border
  if (!is.null(variant_overlay) &&
      "vcf_validated" %in% names(sa_polys) &&
      any(sa_polys$vcf_validated, na.rm = TRUE)) {
    p <- p + ggplot2::geom_polygon(
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
.classify_calls <- function(probs, threshold = 0.5, ambiguous = NULL) {
  cls <- ifelse(probs >= threshold, "methylated", "unmethylated")
  if (!is.null(ambiguous) && ambiguous > 0) {
    band <- abs(probs - threshold) < ambiguous & probs != threshold
    cls[band] <- "ambiguous"
  }
  cls
}

# Add mod-prob segment layer and colour scale to an existing ggplot.
#
# With a non-NULL `call_threshold`, sites are classified into discrete calls
# (methylated/unmethylated/ambiguous) via .classify_calls() and coloured with
# a manual discrete scale instead of the continuous gradient.
.add_mod_prob_segments <- function(p, sites_plot, half_height, line_width,
                                    colour_low, colour_high,
                                    colour_ambiguous = .CALL_AMBIGUOUS_DEFAULT,
                                    call_threshold = NULL,
                                    call_ambiguous = NULL) {
  if (!is.null(call_threshold) && nrow(sites_plot) > 0L) {
    sites_plot$.call <- .classify_calls(sites_plot$mod_prob,
                                        call_threshold, call_ambiguous)
    vals <- c(unmethylated = colour_low, methylated = colour_high,
              ambiguous = colour_ambiguous)
    return(
      p +
        ggplot2::geom_segment(
          data = sites_plot,
          ggplot2::aes(
            x = .data$position, xend = .data$position,
            y = .data$lane - half_height, yend = .data$lane + half_height,
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
        y      = .data$lane - half_height,
        yend   = .data$lane + half_height,
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
