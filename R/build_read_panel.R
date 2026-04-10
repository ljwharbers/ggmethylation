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
# @param region_start Integer.
# @param region_end Integer.
# @return data.frame with x, y, polygon_id, sa_chrom, lane columns,
#   or a 0-row data.frame if no SA reads.
.make_sa_overlay_polygons <- function(reads, arrow_w, half_height,
                                      region_start, region_end) {
  sa_reads <- reads[!is.na(reads$sa_chrom), , drop = FALSE]
  if (nrow(sa_reads) == 0L) {
    return(data.frame(
      x = numeric(0), y = numeric(0), polygon_id = character(0),
      sa_chrom = character(0), lane = numeric(0),
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
        x          = xs,
        y          = ys,
        polygon_id = pid,
        sa_chrom   = r$sa_chrom,
        lane       = ln,
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
#' Constructs the top read panel (horizontal read bars + modification
#' probability lines) for a single `methylation_data` object that has already
#' been sorted and lane-packed. This function encapsulates the grouped,
#' strand-coloured, and plain rendering branches so that both
#' `plot_methylation()` (single-sample) and `.plot_multi_methylation()`
#' (multi-sample) can share the same drawing code.
#'
#' @param data A `methylation_data` object. `$reads` must already contain a
#'   `lane` column (from `pack_reads()`) and a `mean_mod_prob` column.
#' @param separator_lanes Numeric vector of y-positions where horizontal
#'   dashed separator lines should be drawn (between groups). Pass
#'   `numeric(0)` when no separators are needed.
#' @param region_start Integer. Left boundary of the x-axis.
#' @param region_end Integer. Right boundary of the x-axis.
#' @param colour_low Colour for low modification probability.
#' @param colour_high Colour for high modification probability.
#' @param line_width Linewidth of modification site markers.
#' @param colour_strand Logical. Colour read bars by strand when ungrouped.
#' @param strand_colours Named character vector with `"+"` and `"-"` entries.
#' @param group_colours Named character vector of colours per group, or NULL.
#' @param mod_code_shapes Named integer vector mapping mod codes to point
#'   shapes.
#' @param show_x_axis Logical. When `FALSE` (default), x-axis text and ticks
#'   are hidden. Set to `TRUE` for the bottom-most read panel.
#' @param variant_bases Data.frame returned by [extract_variant_bases()], or
#'   `NULL`. When non-NULL, per-read base letters are drawn at variant positions
#'   and coloured by `variant_class` (ref/alt/del/other).
#' @param variant_positions Numeric vector of genomic positions at which to draw
#'   vertical dashed marker lines, or `NULL`.
#' @param show_cigar Logical. When `TRUE`, structural variants from CIGAR
#'   strings are overlaid on reads. Default `FALSE`.
#' @param cigar_features Data.frame of CIGAR features (from
#'   `methylation_data$cigar_features`), or `NULL`.
#' @param min_indel_size Integer. Minimum size (in bp) for insertions and
#'   deletions to be displayed. Features smaller than this threshold are
#'   suppressed. Default `50`.
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
                             line_width,
                             colour_strand,
                             strand_colours,
                             group_colours,
                             mod_code_shapes,
                             show_x_axis        = FALSE,
                             variant_bases      = NULL,
                             variant_positions  = NULL,
                             show_cigar         = FALSE,
                             cigar_features     = NULL,
                             min_indel_size     = 50L,
                             show_supplementary = FALSE) {
  codes      <- unique(data$sites$mod_code)
  multi_code <- length(codes) > 1L

  # When show_cigar is TRUE, split reads on large deletions so the thick read
  # bar has IGV-style gaps instead of running through deletion regions.
  reads_plot <- data$reads
  if (isTRUE(show_cigar) && !is.null(cigar_features) && nrow(cigar_features) > 0L) {
    large_dels <- cigar_features[
      cigar_features$type == "D" & cigar_features$length >= min_indel_size,
      , drop = FALSE
    ]
    if (nrow(large_dels) > 0L) {
      reads_plot <- .split_reads_on_deletions(data$reads, large_dels)
    }
  }

  # --- Annotate split-read segments with first/last flags ---
  # For each read_name, segments ordered by start position.  The first segment
  # in genomic order is marked is_first_segment = TRUE; the last is
  # is_last_segment = TRUE.  Unsplit reads get both TRUE.
  reads_plot$is_first_segment <- TRUE
  reads_plot$is_last_segment  <- TRUE

  rn_tab <- table(reads_plot$read_name)
  split_names <- names(rn_tab[rn_tab > 1L])
  if (length(split_names) > 0L) {
    for (rn in split_names) {
      idx <- which(reads_plot$read_name == rn)
      ord <- order(reads_plot$start[idx])
      sorted_idx <- idx[ord]
      reads_plot$is_first_segment[sorted_idx] <- FALSE
      reads_plot$is_last_segment[sorted_idx]  <- FALSE
      reads_plot$is_first_segment[sorted_idx[1L]] <- TRUE
      reads_plot$is_last_segment[sorted_idx[length(sorted_idx)]] <- TRUE
    }
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
    sa_reads <- data$reads[!is.na(data$reads$sa_chrom), , drop = FALSE]
    if (nrow(sa_reads) > 0L) {
      keep <- rep(TRUE, nrow(sites_plot))
      for (i in seq_len(nrow(sa_reads))) {
        r        <- sa_reads[i, ]
        rn       <- r$read_name
        s        <- r$start
        e        <- r$end
        cs       <- r$clip_side
        if (is.na(cs)) next
        read_len <- e - s
        sa_ext   <- min(arrow_w, read_len)
        idx      <- sites_plot$read_name == rn
        if (cs %in% c("left",  "both")) keep <- keep & !(idx & sites_plot$position <= s + sa_ext)
        if (cs %in% c("right", "both")) keep <- keep & !(idx & sites_plot$position >= e - sa_ext)
      }
      sites_plot <- sites_plot[keep, , drop = FALSE]
    }
  }

  # --- Build read polygons ---
  read_polys <- .make_read_polygons(reads_plot, arrow_w, half_height)

  p <- ggplot2::ggplot()

  # --- Read bar colouring ---
  if (!is.null(data$group_tag)) {
    # Grouped: fill read polygons by group, then colour scale for dots
    p <- p +
      ggplot2::geom_polygon(
        data = read_polys,
        ggplot2::aes(
          x     = .data$x,
          y     = .data$y,
          group = .data$polygon_id,
          fill  = .data$group
        ),
        colour = NA
      )

    if (!is.null(group_colours)) {
      p <- p + ggplot2::scale_fill_manual(
        values   = group_colours,
        na.value = "grey50",
        name     = "Group"
      )
    } else {
      p <- p + ggplot2::scale_fill_discrete(name = "Group")
    }

    # --- SA overlay (grouped) ---
    if (isTRUE(show_supplementary) && "sa_chrom" %in% names(reads_plot)) {
      sa_polys <- .make_sa_overlay_polygons(
        reads_plot, arrow_w * 1.5, half_height, region_start, region_end
      )
      if (nrow(sa_polys) > 0L) {
        p <- p +
          ggnewscale::new_scale_fill() +
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

      }
    }

    p <- p +
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

    if (length(separator_lanes) > 0L) {
      p <- p +
        ggplot2::geom_hline(
          yintercept = separator_lanes,
          linetype = "dashed", colour = "grey40", linewidth = 0.4
        )
    }
  } else {
    # Ungrouped path
    if (isTRUE(colour_strand)) {
      # Strand-coloured: fill polygons by strand
      p <- p +
        ggplot2::geom_polygon(
          data = read_polys,
          ggplot2::aes(
            x     = .data$x,
            y     = .data$y,
            group = .data$polygon_id,
            fill  = .data$strand
          ),
          colour = NA
        ) +
        ggplot2::scale_fill_manual(values = strand_colours, name = "Strand")

      # --- SA overlay (strand) ---
      if (isTRUE(show_supplementary) && "sa_chrom" %in% names(reads_plot)) {
        sa_polys <- .make_sa_overlay_polygons(
          reads_plot, arrow_w, half_height, region_start, region_end
        )
        if (nrow(sa_polys) > 0L) {
          p <- p +
            ggnewscale::new_scale_fill() +
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
        }
      }

      p <- p +
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
    } else {
      # Plain (no grouping, no strand colouring)
      p <- p +
        ggplot2::geom_polygon(
          data = read_polys,
          ggplot2::aes(
            x     = .data$x,
            y     = .data$y,
            group = .data$polygon_id
          ),
          fill = "#B0BEC5",
          colour = NA
        )

      # --- SA overlay (plain) ---
      if (isTRUE(show_supplementary) && "sa_chrom" %in% names(reads_plot)) {
        sa_polys <- .make_sa_overlay_polygons(
          reads_plot, arrow_w, half_height, region_start, region_end
        )
        if (nrow(sa_polys) > 0L) {
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
        }
      }

      p <- p +
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
  if (isTRUE(show_cigar) && !is.null(cigar_features) && nrow(cigar_features) > 0L) {
    # Merge lane info
    cf <- merge(
      cigar_features,
      data$reads[, c("read_name", "lane", "start", "end"), drop = FALSE],
      by = "read_name"
    )

    # Filter small indels
    cf <- cf[!(cf$type %in% c("I", "D") & cf$length < min_indel_size), , drop = FALSE]

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

    # Insertions: draw purple "I" markers
    ins_df <- cf[cf$type == "I", , drop = FALSE]
    if (nrow(ins_df) > 0L) {
      p <- p + ggplot2::geom_text(
        data = ins_df,
        ggplot2::aes(x = .data$ref_start, y = .data$lane, label = "I"),
        colour = "#7B1FA2", size = 2, fontface = "bold",
        inherit.aes = FALSE
      )
    }

  }

  # Variant base letters
  if (!is.null(variant_bases) && nrow(variant_bases) > 0L) {
    p <- p +
      ggnewscale::new_scale_colour() +
      ggplot2::geom_text(
        data = variant_bases,
        ggplot2::aes(x = .data$position, y = .data$lane,
                     label = .data$base, colour = .data$variant_class),
        size = 2, fontface = "bold", show.legend = FALSE,
        inherit.aes = FALSE
      ) +
      ggplot2::scale_colour_manual(
        values = c(ref = "#9E9E9E", alt = "#E53935", other = "#FDD835", del = "#212121"),
        guide  = "none"
      )
  }

  # Vertical lines at variant positions
  if (!is.null(variant_positions) && length(variant_positions) > 0L) {
    p <- p +
      ggplot2::geom_vline(
        xintercept = variant_positions,
        linetype = "dashed", colour = "#424242", linewidth = 0.4, alpha = 0.7
      )
  }

  p <- p +
    ggplot2::scale_y_reverse() +
    ggplot2::coord_cartesian(xlim = c(region_start, region_end)) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.y        = ggplot2::element_blank(),
      axis.ticks.y       = ggplot2::element_blank(),
      axis.title.y       = ggplot2::element_blank(),
      panel.grid.minor   = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      legend.text        = ggplot2::element_text(size = ggplot2::rel(0.75)),
      legend.title       = ggplot2::element_text(size = ggplot2::rel(0.75)),
      legend.key.size    = ggplot2::unit(0.4, "cm")
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
