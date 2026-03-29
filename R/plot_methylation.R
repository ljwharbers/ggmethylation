#' Plot read-level methylation data
#'
#' Creates a ggplot2 visualisation of read-level base modification data. When
#' no grouping is present, produces a single panel showing reads as grey bars
#' with coloured modification dots. When groups are present, adds a bottom
#' panel with loess-smoothed mean modification probability per group and
#' combines the panels using patchwork.
#'
#' @param data A `methylation_data` object returned by [read_methylation()].
#' @param sort_by Character vector of column names from `data$reads` used
#'   to sort reads before packing into lanes. Default NULL uses `c("start")`
#'   when ungrouped or `c("group", "start", "mean_mod_prob")` when grouped.
#' @param colour_low Colour for low modification probability (default
#'   `"#BDBDBD"`).
#' @param colour_high Colour for high modification probability (default
#'   `"#C62828"`).
#' @param dot_size Size of modification dots (default 1.5).
#' @param group_colours Named character vector of colours per group, or NULL
#'   to use ggplot2 defaults.
#' @param smooth_span Loess smoothing span for the bottom panel (default 0.3).
#' @param panel_heights Numeric vector of length 2 giving relative heights of
#'   the top (reads) and bottom (smooth) panels (default `c(3, 1)`).
#'
#' @return A [ggplot2::ggplot] object (ungrouped) or a
#'   [patchwork::patchwork] composite (grouped).
#'
#' @examples
#' \dontrun{
#' md <- read_methylation("sample.bam", "chr1:1000-2000")
#' plot_methylation(md)
#'
#' md <- read_methylation("sample.bam", "chr1:1000-2000", group_tag = "HP")
#' plot_methylation(md, group_colours = c("1" = "steelblue", "2" = "coral"))
#' }
#'
#' @export
plot_methylation <- function(data, sort_by = NULL,
                             colour_low = "#BDBDBD",
                             colour_high = "#C62828",
                             dot_size = 1.5,
                             group_colours = NULL,
                             smooth_span = 0.3,
                             panel_heights = c(3, 1)) {
  # --- 1. Validate input ---
  if (!inherits(data, "methylation_data")) {
    stop(
      "'data' must be a 'methylation_data' object from read_methylation().",
      call. = FALSE
    )
  }

  region_start <- GenomicRanges::start(data$region)
  region_end <- GenomicRanges::end(data$region)

  # --- 2. Handle empty data ---
  if (nrow(data$reads) == 0L) {
    p <- ggplot2::ggplot() +
      ggplot2::annotate(
        "text", x = 0.5, y = 0.5,
        label = "No reads in region", size = 5
      ) +
      ggplot2::theme_void()
    return(p)
  }

  # --- 3. Compute mean mod prob per read (for sorting) ---
  read_means <- stats::setNames(
    tapply(data$sites$mod_prob, data$sites$read_name, mean),
    NULL
  )
  data$reads$mean_mod_prob <- as.numeric(
    read_means[data$reads$read_name]
  )
  data$reads$mean_mod_prob[is.na(data$reads$mean_mod_prob)] <- 0

  # --- 4. Sort reads ---
  if (is.null(sort_by)) {
    if (is.null(data$group_tag)) {
      sort_by <- "start"
    } else {
      sort_by <- c("start", "group", "mean_mod_prob")
    }
  }

  sort_args <- lapply(sort_by, function(col) data$reads[[col]])
  ord <- do.call(order, sort_args)
  data$reads <- data$reads[ord, , drop = FALSE]
  rownames(data$reads) <- NULL

  # --- 5. Pack reads into lanes ---
  data$reads$lane <- integer(nrow(data$reads))
  separator_lanes <- numeric(0)

  if (!is.null(data$group_tag)) {
    groups_ordered <- sort(unique(data$reads$group[!is.na(data$reads$group)]))
    lane_offset <- 0L
    for (grp in groups_ordered) {
      idx <- which(data$reads$group == grp)
      data$reads$lane[idx] <- pack_reads(data$reads[idx, ]) + lane_offset
      lane_offset <- max(data$reads$lane[idx]) + 2L
      separator_lanes <- c(separator_lanes, lane_offset - 1L)
    }
    # Drop the trailing separator (after the last group)
    separator_lanes <- separator_lanes[-length(separator_lanes)]
  } else {
    data$reads$lane <- pack_reads(data$reads)
  }

  # --- 6. Build top panel ---
  # Merge lane info into sites
  sites_plot <- merge(
    data$sites,
    data$reads[, c("read_name", "lane"), drop = FALSE],
    by = "read_name"
  )

  if (!is.null(data$group_tag)) {
    # Grouped: colour read bars by group, then add a second colour scale for dots
    p_top <- ggplot2::ggplot() +
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
      p_top <- p_top +
        ggplot2::scale_colour_manual(values = group_colours, name = "Group")
    } else {
      p_top <- p_top + ggplot2::scale_colour_discrete(name = "Group")
    }

    p_top <- p_top +
      ggnewscale::new_scale_colour() +
      ggplot2::geom_point(
        data = sites_plot,
        ggplot2::aes(
          x = .data$position, y = .data$lane,
          colour = .data$mod_prob
        ),
        size = dot_size
      ) +
      ggplot2::scale_colour_gradient(
        low = colour_low, high = colour_high,
        limits = c(0, 1),
        name = "Modification\nprobability"
      )

    if (length(separator_lanes) > 0L) {
      p_top <- p_top +
        ggplot2::geom_hline(
          yintercept = separator_lanes,
          linetype = "dashed", colour = "grey40", linewidth = 0.4
        )
    }
  } else {
    # Ungrouped: fixed bar colour, single colour scale for dots
    p_top <- ggplot2::ggplot() +
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
          x = .data$position, y = .data$lane,
          colour = .data$mod_prob
        ),
        size = dot_size
      ) +
      ggplot2::scale_colour_gradient(
        low = colour_low, high = colour_high,
        limits = c(0, 1),
        name = "Modification\nprobability"
      )
  }

  if (!is.null(data$snv_position)) {
    p_top <- p_top +
      ggplot2::geom_vline(
        xintercept = data$snv_position,
        linetype = "dashed", colour = "black", linewidth = 0.5
      )
  }

  p_top <- p_top +
    ggplot2::scale_y_reverse() +
    ggplot2::scale_x_continuous(limits = c(region_start, region_end)) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank()
    ) +
    ggplot2::labs(x = NULL)

  # --- 7. Build bottom panel ---
  if (is.null(data$group_tag)) {
    # Ungrouped: smooth all sites together as one line
    sites_smooth <- data$sites
    sites_smooth$group <- "all"
    smoothed <- smooth_methylation(sites_smooth, group_col = "group",
                                   span = smooth_span)

    p_bottom <- ggplot2::ggplot(
      smoothed,
      ggplot2::aes(x = .data$position, y = .data$mean_prob)
    ) +
      ggplot2::geom_line(linewidth = 1, colour = "#C62828") +
      ggplot2::scale_y_continuous(
        limits = c(0, 1),
        name = "Mean modification\nprobability"
      ) +
      ggplot2::scale_x_continuous(limits = c(region_start, region_end)) +
      ggplot2::theme_minimal() +
      ggplot2::theme(panel.grid.minor = ggplot2::element_blank()) +
      ggplot2::labs(x = "Genomic position (bp)")
  } else {
    # Grouped: one smoothed line per group
    smoothed <- smooth_methylation(
      data$sites,
      group_col = "group",
      span = smooth_span
    )

    p_bottom <- ggplot2::ggplot(
      smoothed,
      ggplot2::aes(
        x = .data$position, y = .data$mean_prob,
        colour = .data$group
      )
    ) +
      ggplot2::geom_line(linewidth = 1) +
      ggplot2::scale_y_continuous(
        limits = c(0, 1),
        name = "Mean modification\nprobability"
      ) +
      ggplot2::scale_x_continuous(limits = c(region_start, region_end)) +
      ggplot2::theme_minimal() +
      ggplot2::theme(panel.grid.minor = ggplot2::element_blank()) +
      ggplot2::labs(x = "Genomic position (bp)", colour = "Group")

    if (!is.null(group_colours)) {
      p_bottom <- p_bottom +
        ggplot2::scale_colour_manual(values = group_colours)
    }
  }

  if (!is.null(data$snv_position)) {
    p_bottom <- p_bottom +
      ggplot2::geom_vline(
        xintercept = data$snv_position,
        linetype = "dashed", colour = "black", linewidth = 0.5
      )
  }

  # --- 9. Combine with patchwork ---
  patchwork::wrap_plots(p_top, p_bottom, ncol = 1, heights = panel_heights)
}
