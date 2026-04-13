# Internal function — not exported

#' Build ggplot2 layers for SV spans on reads
#'
#' For each SV in `sv_df` (type DEL/DUP/INV), intersects with overlapping reads
#' and returns a flat list of ggplot2 layer objects (preceded by
#' `ggnewscale::new_scale_colour()`) that draw the SV spans on the read panel.
#'
#' @param reads Data.frame with columns `read_name`, `start`, `end`, `lane`
#'   (the `$reads` element of a `methylation_data` object).
#' @param sv_df Data.frame with columns `position` (int), `end` (int), `type`
#'   (chr: `"DEL"`, `"DUP"`, or `"INV"`). Only the SV rows from
#'   `variant_data$variants`. May be `NULL` or zero-row.
#' @param region_start Integer. Left boundary of the visible plot window.
#' @param region_end Integer. Right boundary of the visible plot window.
#'
#' @return A flat list of ggplot2 layer objects (`ggnewscale::new_scale_colour()`
#'   plus `geom_segment` and `geom_text` layers), or `NULL` when there are no
#'   SV–read overlaps to draw.
#'
#' @keywords internal
build_sv_layer <- function(reads, sv_df, region_start, region_end) {
  # Return NULL for empty/NULL input
  if (is.null(sv_df) || nrow(sv_df) == 0L) return(NULL)

  # --- Build per-read per-SV segment data frame ---
  result_list <- vector("list", nrow(reads) * nrow(sv_df))
  result_idx  <- 0L

  for (v in seq_len(nrow(sv_df))) {
    sv <- sv_df[v, , drop = FALSE]
    sv_pos  <- sv$position
    sv_end  <- sv$end
    sv_type <- sv$type

    # Clip SV span to visible plot window
    span_left  <- max(sv_pos, region_start)
    span_right <- min(sv_end, region_end)

    # Skip SVs entirely outside the plot window
    if (span_right < span_left) next

    clipped_left  <- sv_pos < region_start
    clipped_right <- sv_end > region_end

    for (i in seq_len(nrow(reads))) {
      r_start <- reads$start[i]
      r_end   <- reads$end[i]

      # Check read overlaps the SV span
      if (r_start > sv_end || r_end < sv_pos) next

      # Clip to the read
      seg_left  <- max(span_left,  r_start)
      seg_right <- min(span_right, r_end)

      if (seg_right < seg_left) next

      result_idx <- result_idx + 1L
      result_list[[result_idx]] <- data.frame(
        read_name     = reads$read_name[i],
        lane          = reads$lane[i],
        seg_left      = seg_left,
        seg_right     = seg_right,
        type          = sv_type,
        clipped_left  = clipped_left,
        clipped_right = clipped_right,
        stringsAsFactors = FALSE
      )
    }
  }

  if (result_idx == 0L) return(NULL)

  seg_df <- do.call(rbind, result_list[seq_len(result_idx)])
  rownames(seg_df) <- NULL

  # Colour lookup per type
  sv_colours <- c(DEL = "#B71C1C", DUP = "#00695C", INV = "#4A148C")

  # Attach colour column for use in manual scale
  seg_df$sv_colour <- sv_colours[seg_df$type]

  # --- Subset by type ---
  del_df <- seg_df[seg_df$type == "DEL", , drop = FALSE]
  dup_df <- seg_df[seg_df$type == "DUP", , drop = FALSE]
  inv_df <- seg_df[seg_df$type == "INV", , drop = FALSE]

  # Chevron labels (clipping indicators)
  left_clips  <- seg_df[seg_df$clipped_left,  , drop = FALSE]
  right_clips <- seg_df[seg_df$clipped_right, , drop = FALSE]

  layers <- list()

  # --- DEL: dark red bracket-style line + tick marks at ends ---
  if (nrow(del_df) > 0L) {
    # Main horizontal segment
    layers <- c(layers, list(
      ggplot2::geom_segment(
        data = del_df,
        ggplot2::aes(
          x    = .data$seg_left,
          xend = .data$seg_right,
          y    = .data$lane,
          yend = .data$lane
        ),
        colour      = "#B71C1C",
        linewidth   = 1.2,
        inherit.aes = FALSE
      ),
      # Left tick (downward: lane-0.3 to lane+0.3)
      ggplot2::geom_segment(
        data = del_df,
        ggplot2::aes(
          x    = .data$seg_left,
          xend = .data$seg_left,
          y    = .data$lane - 0.3,
          yend = .data$lane + 0.3
        ),
        colour      = "#B71C1C",
        linewidth   = 1.2,
        inherit.aes = FALSE
      ),
      # Right tick
      ggplot2::geom_segment(
        data = del_df,
        ggplot2::aes(
          x    = .data$seg_right,
          xend = .data$seg_right,
          y    = .data$lane - 0.3,
          yend = .data$lane + 0.3
        ),
        colour      = "#B71C1C",
        linewidth   = 1.2,
        inherit.aes = FALSE
      )
    ))
  }

  # --- DUP: teal double parallel lines offset ±0.12 from lane centre ---
  if (nrow(dup_df) > 0L) {
    dup_upper <- dup_df
    dup_lower <- dup_df
    dup_upper$lane_offset <- dup_df$lane - 0.12
    dup_lower$lane_offset <- dup_df$lane + 0.12

    layers <- c(layers, list(
      ggplot2::geom_segment(
        data = dup_upper,
        ggplot2::aes(
          x    = .data$seg_left,
          xend = .data$seg_right,
          y    = .data$lane_offset,
          yend = .data$lane_offset
        ),
        colour      = "#00695C",
        linewidth   = 0.8,
        inherit.aes = FALSE
      ),
      ggplot2::geom_segment(
        data = dup_lower,
        ggplot2::aes(
          x    = .data$seg_left,
          xend = .data$seg_right,
          y    = .data$lane_offset,
          yend = .data$lane_offset
        ),
        colour      = "#00695C",
        linewidth   = 0.8,
        inherit.aes = FALSE
      )
    ))
  }

  # --- INV: purple dashed single line ---
  if (nrow(inv_df) > 0L) {
    layers <- c(layers, list(
      ggplot2::geom_segment(
        data = inv_df,
        ggplot2::aes(
          x    = .data$seg_left,
          xend = .data$seg_right,
          y    = .data$lane,
          yend = .data$lane
        ),
        colour      = "#4A148C",
        linewidth   = 1.2,
        linetype    = "dashed",
        inherit.aes = FALSE
      )
    ))
  }

  # --- Clipping chevrons ---
  has_left_clips  <- nrow(left_clips)  > 0L
  has_right_clips <- nrow(right_clips) > 0L

  if (has_left_clips || has_right_clips) {
    layers <- c(layers, list(ggnewscale::new_scale_colour()))

    if (has_left_clips) {
      layers <- c(layers, list(
        ggplot2::geom_text(
          data = left_clips,
          ggplot2::aes(
            x      = .data$seg_left,
            y      = .data$lane,
            label  = "\u00AB",
            colour = .data$sv_colour
          ),
          size        = 2,
          vjust       = -0.5,
          inherit.aes = FALSE,
          show.legend = FALSE
        )
      ))
    }

    if (has_right_clips) {
      layers <- c(layers, list(
        ggplot2::geom_text(
          data = right_clips,
          ggplot2::aes(
            x      = .data$seg_right,
            y      = .data$lane,
            label  = "\u00BB",
            colour = .data$sv_colour
          ),
          size        = 2,
          vjust       = -0.5,
          inherit.aes = FALSE,
          show.legend = FALSE
        )
      ))
    }

    layers <- c(layers, list(ggplot2::scale_colour_identity(guide = "none")))
  }

  layers
}
