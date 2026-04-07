# Internal function — not exported

#' Build a ggplot2 gene annotation panel
#'
#' Renders transcript models (intron lines, exon/CDS/UTR boxes, strand arrows,
#' gene labels) for a `gene_annotations` object into a single ggplot2 panel
#' suitable for assembly with patchwork.
#'
#' When the annotations object contains CDS data (`$cds`), an IGV-style
#' three-tier rendering is used: thin intron line, medium-height UTR boxes,
#' and full-height CDS boxes. Without CDS data a uniform exon box is drawn.
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

  # --- Check for CDS/UTR data ---
  has_cds  <- !is.null(annotations$cds)  && nrow(annotations$cds)  > 0L
  has_utr5 <- !is.null(annotations$utr5) && nrow(annotations$utr5) > 0L
  has_utr3 <- !is.null(annotations$utr3) && nrow(annotations$utr3) > 0L

  if (has_cds) {
    # Join CDS data with row numbers
    cds_joined <- merge(
      annotations$cds,
      transcripts_df[, c("tx_id", "row"), drop = FALSE],
      by = "tx_id"
    )
  }

  utr_df <- NULL
  if (has_utr5 || has_utr3) {
    utr_parts <- list()
    if (has_utr5) {
      utr5_rows <- merge(annotations$utr5, transcripts_df[, c("tx_id", "row"), drop = FALSE], by = "tx_id")
      utr_parts[["utr5"]] <- utr5_rows
    }
    if (has_utr3) {
      utr3_rows <- merge(annotations$utr3, transcripts_df[, c("tx_id", "row"), drop = FALSE], by = "tx_id")
      utr_parts[["utr3"]] <- utr3_rows
    }
    utr_df <- do.call(rbind, utr_parts)
  }

  # --- Build strand arrow data (placed at intron midpoints) ---
  arrow_df <- do.call(rbind, lapply(seq_len(nrow(transcripts_df)), function(i) {
    tx       <- transcripts_df[i, ]
    tx_exons <- exons_joined[exons_joined$tx_id == tx$tx_id, , drop = FALSE]
    if (nrow(tx_exons) == 0L) return(NULL)
    tx_exons <- tx_exons[order(tx_exons$exon_start), , drop = FALSE]

    # Build intron intervals from gaps between exons
    intron_starts <- tx_exons$exon_end[-nrow(tx_exons)] + 1L
    intron_ends   <- tx_exons$exon_start[-1L] - 1L
    intron_widths <- intron_ends - intron_starts

    # Also consider full span if single exon
    if (length(intron_starts) == 0L) {
      width <- tx$tx_end - tx$tx_start
      if (width < 200L) return(NULL)
      intron_starts <- tx$tx_start
      intron_ends   <- tx$tx_end
      intron_widths <- width
    }

    wide_introns <- intron_widths > 200L
    if (!any(wide_introns)) return(NULL)

    positions <- (intron_starts[wide_introns] + intron_ends[wide_introns]) / 2L
    arrow_len <- pmax(200L, (intron_ends[wide_introns] - intron_starts[wide_introns]) * 0.05)

    if (tx$strand == "+") {
      data.frame(x = positions, xend = positions + arrow_len, y = tx$row, yend = tx$row)
    } else {
      data.frame(x = positions, xend = positions - arrow_len, y = tx$row, yend = tx$row)
    }
  }))

  n_rows <- nrow(transcripts_df)

  # --- Assemble plot ---
  p <- ggplot2::ggplot() +
    # Thin intron/body line
    ggplot2::geom_segment(
      data = transcripts_df,
      ggplot2::aes(x = .data$tx_start, xend = .data$tx_end, y = .data$row, yend = .data$row),
      colour = "black", linewidth = 0.3
    )

  if (has_cds) {
    # UTR regions (medium height)
    if (!is.null(utr_df) && nrow(utr_df) > 0L) {
      p <- p + ggplot2::geom_rect(
        data = utr_df,
        ggplot2::aes(
          xmin = .data$utr_start, xmax = .data$utr_end,
          ymin = .data$row - 0.15, ymax = .data$row + 0.15
        ),
        fill = "black", colour = NA
      )
    }
    # CDS regions (full height)
    p <- p + ggplot2::geom_rect(
      data = cds_joined,
      ggplot2::aes(
        xmin = .data$cds_start, xmax = .data$cds_end,
        ymin = .data$row - 0.3,  ymax = .data$row + 0.3
      ),
      fill = "black", colour = NA
    )
  } else {
    # Fallback: uniform exon boxes
    p <- p + ggplot2::geom_rect(
      data = exons_joined,
      ggplot2::aes(
        xmin = .data$exon_start, xmax = .data$exon_end,
        ymin = .data$row - 0.3,  ymax = .data$row + 0.3
      ),
      fill = "black", colour = NA
    )
  }

  # Gene labels at left margin
  p <- p + ggplot2::geom_text(
    data = transcripts_df,
    ggplot2::aes(x = region_start, y = .data$row, label = .data$gene_name),
    hjust = 1.1, size = 2.5, colour = "black"
  ) +
  ggplot2::coord_cartesian(clip = "off") +
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
