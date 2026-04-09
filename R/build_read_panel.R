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
    too_short <- read_len < 2 * aw

    # Determine whether this segment gets an arrowhead
    draw_right_arrow <- (st == "+") && is_last && !too_short
    draw_left_arrow  <- (st == "-") && is_first && !too_short

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

# Build colored overlay polygons at breakpoint ends of supplementary reads.
#
# For each read with a non-NA sa_chrom, an overlay polygon is placed at the
# clipped end of the read to indicate the SA partner chromosome.
#
# @param reads data.frame with SA-annotated reads (sa_chrom, clip_side columns).
# @param arrow_w Numeric.
# @param half_height Numeric.
# @param region_start Integer.
# @param region_end Integer.
# @return data.frame with x, y, polygon_id, sa_chrom, label_x columns,
#   or a 0-row data.frame if no SA reads.
.make_sa_overlay_polygons <- function(reads, arrow_w, half_height,
                                      region_start, region_end) {
  sa_reads <- reads[!is.na(reads$sa_chrom), , drop = FALSE]
  if (nrow(sa_reads) == 0L) {
    return(data.frame(
      x = numeric(0), y = numeric(0), polygon_id = character(0),
      sa_chrom = character(0), label_x = numeric(0),
      lane = numeric(0),
      stringsAsFactors = FALSE
    ))
  }

  region_span <- region_end - region_start
  poly_list <- vector("list", 2L * nrow(sa_reads))  # upper bound: 2 overlays per read for clip_side=="both"
  pid_counter <- 0L

  for (i in seq_len(nrow(sa_reads))) {
    r  <- sa_reads[i, , drop = FALSE]
    s  <- r$start
    e  <- r$end
    ln <- r$lane
    hh <- half_height
    aw <- arrow_w
    st <- r$strand
    cs <- r$clip_side

    read_len <- e - s
    sa_extent <- max(aw * 4, 0.03 * region_span)
    sa_extent <- min(sa_extent, 0.5 * read_len)

    # Determine which sides to overlay
    sides <- character(0)
    if (is.na(cs)) {
      # Fallback: put overlay at arrowhead end
      if (st == "+") sides <- "right" else sides <- "left"
    } else if (cs == "both") {
      sides <- c("left", "right")
    } else {
      sides <- cs
    }

    for (side in sides) {
      pid_counter <- pid_counter + 1L
      pid <- paste0("sa_", pid_counter)

      if (side == "right") {
        ov_start <- e - sa_extent
        ov_end   <- e
        # Include arrowhead if + strand and this is last segment
        is_last <- if ("is_last_segment" %in% names(r)) r$is_last_segment else TRUE
        if (st == "+" && is_last) {
          xs <- c(ov_start, ov_end, ov_end + aw, ov_end, ov_start)
          ys <- c(ln - hh, ln - hh, ln, ln + hh, ln + hh)
        } else {
          xs <- c(ov_start, ov_end, ov_end, ov_start)
          ys <- c(ln - hh, ln - hh, ln + hh, ln + hh)
        }
        lx <- (ov_start + ov_end) / 2
      } else {
        # side == "left"
        ov_start <- s
        ov_end   <- s + sa_extent
        is_first <- if ("is_first_segment" %in% names(r)) r$is_first_segment else TRUE
        if (st == "-" && is_first) {
          xs <- c(ov_start - aw, ov_start, ov_end, ov_end, ov_start)
          ys <- c(ln, ln - hh, ln - hh, ln + hh, ln + hh)
        } else {
          xs <- c(ov_start, ov_end, ov_end, ov_start)
          ys <- c(ln - hh, ln - hh, ln + hh, ln + hh)
        }
        lx <- (ov_start + ov_end) / 2
      }

      n_verts <- length(xs)
      poly_list[[pid_counter]] <- data.frame(
        x          = xs,
        y          = ys,
        polygon_id = pid,
        sa_chrom   = r$sa_chrom,
        label_x    = lx,
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
#' probability dots) for a single `methylation_data` object that has already
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
#' @param dot_size Size of modification dots.
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
                             dot_size,
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

  # Merge lane info into sites
  sites_plot <- merge(
    data$sites,
    data$reads[, c("read_name", "lane"), drop = FALSE],
    by = "read_name"
  )

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

  # --- Arrow geometry parameters ---
  arrow_w     <- (region_end - region_start) * 0.003
  half_height <- 0.35

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

        # Labels
        sa_label_df <- unique(sa_polys[, c("polygon_id", "label_x", "lane", "sa_chrom"),
                                       drop = FALSE])
        p <- p +
          ggplot2::geom_text(
            data = sa_label_df,
            ggplot2::aes(
              x = .data$label_x, y = .data$lane, label = .data$sa_chrom
            ),
            size = 1.5, fontface = "bold", colour = "white",
            hjust = 0.5, vjust = 0.5, inherit.aes = FALSE
          )
      }
    }

    p <- p +
      ggplot2::geom_point(
        data = sites_plot,
        ggplot2::aes(
          x      = .data$position,
          y      = .data$lane,
          colour = .data$mod_prob,
          shape  = if (multi_code) .data$mod_code else NULL
        ),
        size = dot_size
      ) +
      ggplot2::scale_colour_gradient(
        low = colour_low, high = colour_high,
        limits = c(0, 1),
        name = "Modification\nprobability"
      )

    if (multi_code) {
      p <- p +
        ggplot2::scale_shape_manual(
          values = mod_code_shapes,
          name   = "Modification"
        )
    }

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

          sa_label_df <- unique(sa_polys[, c("polygon_id", "label_x", "lane", "sa_chrom"),
                                         drop = FALSE])
          p <- p +
            ggplot2::geom_text(
              data = sa_label_df,
              ggplot2::aes(
                x = .data$label_x, y = .data$lane, label = .data$sa_chrom
              ),
              size = 1.5, fontface = "bold", colour = "white",
              hjust = 0.5, vjust = 0.5, inherit.aes = FALSE
            )
        }
      }

      p <- p +
        ggplot2::geom_point(
          data = sites_plot,
          ggplot2::aes(
            x      = .data$position,
            y      = .data$lane,
            colour = .data$mod_prob,
            shape  = if (multi_code) .data$mod_code else NULL
          ),
          size = dot_size
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

          sa_label_df <- unique(sa_polys[, c("polygon_id", "label_x", "lane", "sa_chrom"),
                                         drop = FALSE])
          p <- p +
            ggplot2::geom_text(
              data = sa_label_df,
              ggplot2::aes(
                x = .data$label_x, y = .data$lane, label = .data$sa_chrom
              ),
              size = 1.5, fontface = "bold", colour = "white",
              hjust = 0.5, vjust = 0.5, inherit.aes = FALSE
            )
        }
      }

      p <- p +
        ggplot2::geom_point(
          data = sites_plot,
          ggplot2::aes(
            x      = .data$position,
            y      = .data$lane,
            colour = .data$mod_prob,
            shape  = if (multi_code) .data$mod_code else NULL
          ),
          size = dot_size
        ) +
        ggplot2::scale_colour_gradient(
          low = colour_low, high = colour_high,
          limits = c(0, 1),
          name = "Modification\nprobability"
        )
    }

    if (multi_code) {
      p <- p +
        ggplot2::scale_shape_manual(
          values = mod_code_shapes,
          name   = "Modification"
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
    ggplot2::scale_x_continuous(limits = c(region_start, region_end)) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.y        = ggplot2::element_blank(),
      axis.ticks.y       = ggplot2::element_blank(),
      axis.title.y       = ggplot2::element_blank(),
      panel.grid.minor   = ggplot2::element_blank(),
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
