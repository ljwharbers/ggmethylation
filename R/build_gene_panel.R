# Internal function — not exported

#' Build a ggplot2 gene annotation panel
#'
#' Renders transcript models (intron lines, exon boxes, strand arrows, gene
#' labels) for a `gene_annotations` object into a single ggplot2 panel
#' suitable for assembly with patchwork.
#'
#' @param annotations A `gene_annotations` object from [read_annotations()].
#' @param region_start Integer. Start coordinate of the plot x-axis.
#' @param region_end Integer. End coordinate of the plot x-axis.
#' @param bottom_label Logical. When `TRUE` (default), the x-axis label
#'   "Genomic position (bp)" is shown. When `FALSE`, the x-axis is hidden
#'   (useful when this panel sits above another panel that carries the label).
#'
#' @return A [ggplot2::ggplot] object.
#'
#' @keywords internal
build_gene_panel <- function(annotations, region_start, region_end,
                             bottom_label = TRUE) {
  transcripts_df <- annotations$transcripts

  # --- Handle empty annotations ---
  if (nrow(transcripts_df) == 0L) {
    p <- ggplot2::ggplot() +
      ggplot2::annotate(
        "text", x = (region_start + region_end) / 2, y = 0.5,
        label = "No gene annotations in region", size = 3.5, colour = "grey50"
      ) +
      ggplot2::scale_x_continuous(limits = c(region_start, region_end),
                                  expand = c(0, 0)) +
      ggplot2::scale_y_continuous(limits = c(0, 1)) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        axis.text.y  = ggplot2::element_blank(),
        axis.ticks.y = ggplot2::element_blank(),
        axis.title.y = ggplot2::element_blank(),
        panel.grid   = ggplot2::element_blank()
      )

    if (bottom_label) {
      p <- p + ggplot2::labs(x = "Genomic position (bp)")
    } else {
      p <- p + ggplot2::labs(x = NULL) +
        ggplot2::theme(
          axis.text.x  = ggplot2::element_blank(),
          axis.ticks.x = ggplot2::element_blank()
        )
    }
    return(p)
  }

  # --- Assign row numbers to transcripts ---
  transcripts_df$row <- seq_len(nrow(transcripts_df))

  # --- Join row numbers onto exons ---
  exons_df     <- annotations$exons
  exons_joined <- merge(
    exons_df,
    transcripts_df[, c("tx_id", "row"), drop = FALSE],
    by = "tx_id"
  )

  # --- Build strand arrow data ---
  arrow_df <- do.call(rbind, lapply(seq_len(nrow(transcripts_df)), function(i) {
    tx    <- transcripts_df[i, ]
    width <- tx$tx_end - tx$tx_start
    if (width < 200L) return(NULL)
    positions <- tx$tx_start + width * c(0.25, 0.50, 0.75)
    arrow_len <- width * 0.02
    if (tx$strand == "+") {
      data.frame(
        x = positions, xend = positions + arrow_len,
        y = tx$row,    yend = tx$row
      )
    } else {
      data.frame(
        x = positions, xend = positions - arrow_len,
        y = tx$row,    yend = tx$row
      )
    }
  }))

  n_rows <- nrow(transcripts_df)

  # --- Assemble plot ---
  p <- ggplot2::ggplot() +
    # Intron/body line
    ggplot2::geom_segment(
      data = transcripts_df,
      ggplot2::aes(
        x = .data$tx_start, xend = .data$tx_end,
        y = .data$row,      yend = .data$row
      ),
      colour = "black", linewidth = 0.4
    ) +
    # Exon boxes
    ggplot2::geom_rect(
      data = exons_joined,
      ggplot2::aes(
        xmin = .data$exon_start, xmax = .data$exon_end,
        ymin = .data$row - 0.3,  ymax = .data$row + 0.3
      ),
      fill = "black", colour = NA
    ) +
    # Gene name labels at left margin
    ggplot2::geom_text(
      data = transcripts_df,
      ggplot2::aes(
        x = region_start, y = .data$row, label = .data$gene_name
      ),
      hjust = 1.1, size = 2.5, colour = "black"
    ) +
    ggplot2::scale_x_continuous(
      limits = c(region_start, region_end), expand = c(0, 0)
    ) +
    ggplot2::scale_y_continuous(limits = c(0.5, n_rows + 0.5)) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.y  = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      panel.grid   = ggplot2::element_blank()
    )

  # --- Add strand arrows ---
  if (!is.null(arrow_df) && nrow(arrow_df) > 0L) {
    p <- p + ggplot2::geom_segment(
      data = arrow_df,
      ggplot2::aes(
        x = .data$x, xend = .data$xend,
        y = .data$y, yend = .data$yend
      ),
      arrow    = ggplot2::arrow(
        length = ggplot2::unit(0.15, "cm"), type = "open"
      ),
      colour    = "black",
      linewidth = 0.4
    )
  }

  # --- X-axis label ---
  if (bottom_label) {
    p <- p + ggplot2::labs(x = "Genomic position (bp)")
  } else {
    p <- p + ggplot2::labs(x = NULL) +
      ggplot2::theme(
        axis.text.x  = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank()
      )
  }

  p
}
