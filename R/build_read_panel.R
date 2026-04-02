# Internal function — not exported

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
                             show_x_axis = FALSE) {
  codes      <- unique(data$sites$mod_code)
  multi_code <- length(codes) > 1L

  # Merge lane info into sites
  sites_plot <- merge(
    data$sites,
    data$reads[, c("read_name", "lane"), drop = FALSE],
    by = "read_name"
  )

  if (!is.null(data$group_tag)) {
    # Grouped: colour read bars by group, then second colour scale for dots
    p <- ggplot2::ggplot() +
      ggplot2::geom_segment(
        data = data$reads,
        ggplot2::aes(
          x = .data$start, xend = .data$end,
          y = .data$lane, yend = .data$lane,
          colour = .data$group
        ),
        linewidth = 2, lineend = "round"
      )

    if (!is.null(group_colours)) {
      p <- p + ggplot2::scale_colour_manual(values = group_colours, name = "Group")
    } else {
      p <- p + ggplot2::scale_colour_discrete(name = "Group")
    }

    p <- p +
      ggnewscale::new_scale_colour() +
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
      p <- ggplot2::ggplot() +
        ggplot2::geom_segment(
          data = data$reads,
          ggplot2::aes(
            x = .data$start, xend = .data$end,
            y = .data$lane, yend = .data$lane,
            colour = .data$strand
          ),
          linewidth = 2, lineend = "round"
        ) +
        ggplot2::scale_colour_manual(values = strand_colours, name = "Strand") +
        ggnewscale::new_scale_colour() +
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
      p <- ggplot2::ggplot() +
        ggplot2::geom_segment(
          data = data$reads,
          ggplot2::aes(
            x = .data$start, xend = .data$end,
            y = .data$lane, yend = .data$lane
          ),
          linewidth = 2, colour = "#B0BEC5", lineend = "round"
        ) +
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
