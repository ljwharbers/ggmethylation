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
#' Gene names are rendered underneath each gene body, centred on the visible
#' portion of the transcript. Non-overlapping genes are packed into shared rows
#' using the same greedy interval-scheduling algorithm as the read panel.
#' Geometry extending past the visible region is pre-clipped so that partial
#' genes still show their exon/CDS/UTR boxes rather than just a thin intron
#' line.
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

  # Helper: clamp start/end columns to [region_start, region_end] and drop
  # zero-width rows (i.e. geometry that falls completely outside the view).
  clip_range <- function(df, start_col, end_col) {
    if (is.null(df) || nrow(df) == 0L) return(df)
    df[[start_col]] <- pmax(df[[start_col]], region_start)
    df[[end_col]]   <- pmin(df[[end_col]],   region_end)
    df[df[[end_col]] > df[[start_col]], , drop = FALSE]
  }

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

  # --- Check for CDS/UTR data ---
  has_cds  <- !is.null(annotations$cds)  && nrow(annotations$cds)  > 0L
  has_utr5 <- !is.null(annotations$utr5) && nrow(annotations$utr5) > 0L
  has_utr3 <- !is.null(annotations$utr3) && nrow(annotations$utr3) > 0L

  # --- Pre-clip all geometry to the visible region ---
  # Transcripts: clip tx_start/tx_end so the intron line is always visible
  # for partial genes; remove genes entirely outside the view.
  transcripts_df <- clip_range(transcripts_df, "tx_start", "tx_end")

  # If clipping left nothing, fall back to the empty-annotations panel.
  if (nrow(transcripts_df) == 0L) {
    return(build_gene_panel(
      structure(
        list(
          transcripts = transcripts_df,
          exons       = annotations$exons[integer(0), , drop = FALSE],
          cds         = NULL,
          utr5        = NULL,
          utr3        = NULL,
          region      = annotations$region
        ),
        class = "gene_annotations"
      ),
      region_start = region_start,
      region_end   = region_end,
      bottom_label = bottom_label
    ))
  }

  # --- Sort by tx_start for stable packing ---
  transcripts_df <- transcripts_df[order(transcripts_df$tx_start), , drop = FALSE]

  # --- Assign rows via greedy interval packing ---
  gap_bp <- max(50L, as.integer((region_end - region_start) * 0.02))
  pack_input <- data.frame(start = transcripts_df$tx_start,
                           end   = transcripts_df$tx_end)
  transcripts_df$row <- pack_reads(pack_input, gap = gap_bp)

  n_rows <- max(transcripts_df$row)

  # --- Apply vertical spacing between rows ---
  # Multiplier > 1 increases the gap between stacked transcripts so that
  # labels sit clear of the row beneath. Box heights remain ±0.3, so the
  # visible gap between adjacent boxes grows from (1 - 0.6) to
  # (row_spacing - 0.6).
  row_spacing        <- 1.3
  transcripts_df$row <- transcripts_df$row * row_spacing

  # --- Join row numbers onto exons ---
  exons_df     <- annotations$exons
  exons_joined <- merge(
    exons_df,
    transcripts_df[, c("tx_id", "row"), drop = FALSE],
    by = "tx_id"
  )
  exons_joined <- clip_range(exons_joined, "exon_start", "exon_end")

  if (has_cds) {
    cds_joined <- merge(
      annotations$cds,
      transcripts_df[, c("tx_id", "row"), drop = FALSE],
      by = "tx_id"
    )
    cds_joined <- clip_range(cds_joined, "cds_start", "cds_end")
  }

  utr_df <- NULL
  if (has_utr5 || has_utr3) {
    utr_parts <- list()
    if (has_utr5) {
      utr5_rows <- merge(annotations$utr5, transcripts_df[, c("tx_id", "row"), drop = FALSE], by = "tx_id")
      utr_parts[["utr5"]] <- clip_range(utr5_rows, "utr_start", "utr_end")
    }
    if (has_utr3) {
      utr3_rows <- merge(annotations$utr3, transcripts_df[, c("tx_id", "row"), drop = FALSE], by = "tx_id")
      utr_parts[["utr3"]] <- clip_range(utr3_rows, "utr_start", "utr_end")
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
    # Filter arrows whose midpoints are outside the view
    in_view   <- positions >= region_start & positions <= region_end
    if (!any(in_view)) return(NULL)
    positions  <- positions[in_view]
    arrow_len  <- pmax(200L, (intron_ends[wide_introns][in_view] -
                              intron_starts[wide_introns][in_view]) * 0.05)

    if (tx$strand == "+") {
      data.frame(x = positions, xend = positions + arrow_len, y = tx$row, yend = tx$row)
    } else {
      data.frame(x = positions, xend = positions - arrow_len, y = tx$row, yend = tx$row)
    }
  }))

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
    if (nrow(cds_joined) > 0L) {
      p <- p + ggplot2::geom_rect(
        data = cds_joined,
        ggplot2::aes(
          xmin = .data$cds_start, xmax = .data$cds_end,
          ymin = .data$row - 0.3,  ymax = .data$row + 0.3
        ),
        fill = "black", colour = NA
      )
    }
    # Non-coding transcripts (e.g. lncRNAs, many LOC* entries, miRNAs): they
    # have exons but no CDS / UTR rows, so they would otherwise render as a
    # bare intron line. Draw uniform full-height exon boxes for them, matching
    # IGV's convention for non-coding RNAs.
    noncoding_exons <- exons_joined[
      !exons_joined$tx_id %in% cds_joined$tx_id, , drop = FALSE
    ]
    if (nrow(noncoding_exons) > 0L) {
      p <- p + ggplot2::geom_rect(
        data = noncoding_exons,
        ggplot2::aes(
          xmin = .data$exon_start, xmax = .data$exon_end,
          ymin = .data$row - 0.3,  ymax = .data$row + 0.3
        ),
        fill = "black", colour = NA
      )
    }
  } else {
    # Fallback: uniform exon boxes
    if (nrow(exons_joined) > 0L) {
      p <- p + ggplot2::geom_rect(
        data = exons_joined,
        ggplot2::aes(
          xmin = .data$exon_start, xmax = .data$exon_end,
          ymin = .data$row - 0.3,  ymax = .data$row + 0.3
        ),
        fill = "black", colour = NA
      )
    }
  }

  # Gene labels centred underneath each gene body (on the visible span)
  label_df           <- transcripts_df
  label_df$label_x   <- (pmax(label_df$tx_start, region_start) +
                          pmin(label_df$tx_end,   region_end)) / 2
  label_df$label_y   <- label_df$row - 0.45

  p <- p +
    ggplot2::geom_text(
      data  = label_df,
      ggplot2::aes(x = .data$label_x, y = .data$label_y, label = .data$gene_name),
      vjust = 1, size = 2.5, colour = "black"
    ) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::scale_x_continuous(
      limits = c(region_start, region_end), expand = c(0, 0)
    ) +
    ggplot2::scale_y_continuous(
      limits = c(row_spacing - 0.8, n_rows * row_spacing + 0.5)
    ) +
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
