# Per-locus insertion visualization for ggmethylation

#' Plot modification probabilities at an insertion locus
#'
#' Produces a composite ggplot / patchwork figure for a single insertion locus,
#' showing carrier and non-carrier reads in a stitched coordinate system
#' (left reference flank | insertion bases in literal bp | right reference flank),
#' with an optional loess-smoothed comparison panel below.
#'
#' Use [list_insertion_loci()] to find loci; pass one row of its output.
#'
#' @param data A `methylation_data` object produced by [read_methylation()].
#' @param locus One row of the data.frame returned by [list_insertion_loci()]
#'   for `data`, e.g. `loci[1, ]` or `loci[loci$locus_id == id, ]`. It carries
#'   the locus' carrier and non-carrier reads and the clustering tolerances.
#' @param flank Integer. Number of reference base pairs to show on each side of
#'   the insertion (default 200).
#' @param show_smoothed Logical. Whether to add a loess-smoothed modification
#'   probability panel below the read panel (default `TRUE`).
#' @param include_noncarriers Logical. Whether to include reads that span the
#'   locus but do not carry the insertion (default `TRUE`). Non-carriers are
#'   shown below carriers with an empty middle region.
#' @param colour_low Character. Colour for low modification probability
#'   (default `"#BDBDBD"`).
#' @param colour_high Character. Colour for high modification probability
#'   (default `"#C62828"`).
#'
#' @return A `patchwork` object (or a single ggplot when `show_smoothed = FALSE`).
#'
#' @examples
#' \dontrun{
#' md   <- read_methylation("sample.bam", "chr21:34500000-34510000")
#' loci <- list_insertion_loci(md, tol_pos = 10L, tol_len = 0.20, min_reads = 2L)
#' plot_insertion_locus(md, loci[1, ])
#' }
#'
#' @export
plot_insertion_locus <- function(data, locus,
                                  flank              = 200L,
                                  show_smoothed      = TRUE,
                                  include_noncarriers = TRUE,
                                  colour_low         = .PROB_GRADIENT$low,
                                  colour_high        = .PROB_GRADIENT$high) {
  if (!inherits(data, "methylation_data")) {
    stop("'data' must be a methylation_data object.", call. = FALSE)
  }
  locus_cols = c("locus_id", "anchor_pos", "median_length", "carriers",
                 "noncarriers", "tol_pos")
  if (!is.data.frame(locus) || nrow(locus) != 1L || !all(locus_cols %in% names(locus))) {
    stop("`locus` must be one row of list_insertion_loci() output, e.g. loci[1, ].",
         call. = FALSE)
  }

  # --- 1. Carrier / non-carrier reads of the locus ---
  locus_id    = locus$locus_id
  anchor      = locus$anchor_pos
  mid_width   = locus$median_length
  tol_pos     = locus$tol_pos
  flank       = as.integer(flank)

  carrier_ins_len  = locus$carriers[[1L]]   # insertion length, named by read
  carrier_names    = names(carrier_ins_len)
  noncarrier_names = locus$noncarriers[[1L]]

  reads_all <- data$reads
  all_names <- reads_all$read_name
  carrier_reads    <- reads_all[all_names %in% carrier_names, , drop = FALSE]
  noncarrier_reads <- reads_all[include_noncarriers & all_names %in% noncarrier_names, , drop = FALSE]

  if (nrow(carrier_reads) == 0L) {
    stop("No carrier reads found for locus '", locus_id, "'.", call. = FALSE)
  }

  # --- 2. Assign lanes via pack_reads ---
  carrier_reads$carrier_status    <- rep("carrier",     nrow(carrier_reads))
  noncarrier_reads$carrier_status <- rep("non-carrier", nrow(noncarrier_reads))

  # Pack carriers and non-carriers separately so carriers come first (lower lane)
  carrier_reads$lane <- pack_reads(carrier_reads)
  n_carrier_lanes    <- max(carrier_reads$lane)

  if (nrow(noncarrier_reads) > 0L) {
    noncarrier_reads$lane <- pack_reads(noncarrier_reads) + n_carrier_lanes + 1L
    reads_packed          <- rbind(carrier_reads, noncarrier_reads)
    separator_lane        <- n_carrier_lanes + 0.5
  } else {
    reads_packed   <- carrier_reads
    separator_lane <- NULL
  }

  # --- 3. Build stitched sub-segments for read polygons ---
  # Each read becomes 2 (non-carrier) or 3 (carrier) sub-segment rows.
  # Stitched x-axis:
  #   Left flank:  ref coords [anchor - flank, anchor]  (unchanged)
  #   Middle:      insertion-internal coords -> stitched x = anchor + (ins_offset - 1)
  #                middle zone spans [anchor, anchor + mid_width - 1]
  #   Right flank: ref coords [anchor, anchor + flank] shifted by mid_width
  #                -> stitched x = ref_x + mid_width

  seg_list <- vector("list", nrow(reads_packed))

  for (i in seq_len(nrow(reads_packed))) {
    r <- reads_packed[i, , drop = FALSE]
    rn <- r$read_name
    ln <- r$lane
    st <- r$strand

    r_start <- pmin(r$bam_pos, r$start)
    r_end   <- r$end

    is_carrier <- r$carrier_status == "carrier"
    ins_len    <- if (is_carrier) carrier_ins_len[[rn]] else 0L

    # Left flank sub-segment
    lf_start <- max(r_start, anchor - flank)
    lf_end   <- min(r_end, anchor)
    if (lf_start > lf_end) lf_end <- lf_start  # zero-width if read starts after anchor

    # Right flank sub-segment (stitched coords)
    rf_ref_end <- min(r_end, anchor + flank)
    rf_start_st <- anchor + mid_width
    rf_end_st   <- anchor + mid_width + max(0L, rf_ref_end - anchor)

    if (is_carrier) {
      mid_end_st <- anchor + min(ins_len, mid_width)
      overflow   <- ins_len > mid_width

      segs <- data.frame(
        read_name         = rn,
        start             = c(lf_start,  anchor,     rf_start_st),
        end               = c(lf_end,    mid_end_st, rf_end_st),
        lane              = ln,
        strand            = st,
        carrier_status    = "carrier",
        is_first_segment  = c(r_start >= (anchor - flank), FALSE, FALSE),
        is_last_segment   = c(FALSE, FALSE, r_end <= (anchor + flank)),
        overflow          = c(FALSE, overflow, FALSE),
        stringsAsFactors  = FALSE
      )
    } else {
      segs <- data.frame(
        read_name         = rn,
        start             = c(lf_start,  rf_start_st),
        end               = c(lf_end,    rf_end_st),
        lane              = ln,
        strand            = st,
        carrier_status    = "non-carrier",
        is_first_segment  = c(r_start >= (anchor - flank), FALSE),
        is_last_segment   = c(FALSE, r_end <= (anchor + flank)),
        overflow          = c(FALSE, FALSE),
        stringsAsFactors  = FALSE
      )
    }
    seg_list[[i]] <- segs
  }

  segs_all <- do.call(rbind, seg_list)
  rownames(segs_all) <- NULL

  # Build polygon data via the existing helper
  half_height <- 0.4
  arrow_w     <- 0.03 * (2 * flank + mid_width)
  poly_data   <- .make_read_polygons(segs_all, arrow_w, half_height)

  # Fill colour by carrier_status
  carrier_fill    <- "#0072B2"
  noncarrier_fill <- "#999999"
  poly_data$fill_col <- ifelse(
    poly_data$carrier_status == "carrier", carrier_fill, noncarrier_fill
  )

  # --- 4. Build mod-site data in stitched coords ---
  sites <- data$sites
  lane_map <- reads_packed[, c("read_name", "lane"), drop = FALSE]

  # Left flank sites: ref positions in [anchor - flank, anchor]
  lf_sites <- sites[sites$position >= (anchor - flank) &
                    sites$position <= anchor, , drop = FALSE]
  lf_sites  <- merge(lf_sites, lane_map, by = "read_name")
  lf_sites$stitch_x <- lf_sites$position  # unchanged

  # Right flank sites: ref positions in (anchor, anchor + flank], shifted
  rf_sites <- sites[sites$position > anchor &
                    sites$position <= (anchor + flank), , drop = FALSE]
  rf_sites  <- merge(rf_sites, lane_map, by = "read_name")
  rf_sites$stitch_x <- rf_sites$position + mid_width

  # Insertion sites: carriers at this locus
  ins_sites <- data$insertion_sites
  ins_locus <- if (nrow(ins_sites) > 0L) {
    ins_sites[ins_sites$read_name %in% carrier_names &
              abs(ins_sites$ref_anchor - anchor) <= tol_pos, , drop = FALSE]
  } else {
    ins_sites[integer(0), , drop = FALSE]
  }
  ins_locus <- merge(ins_locus, lane_map, by = "read_name")
  ins_locus$stitch_x <- anchor + ins_locus$ins_offset - 1L

  # Combine all sites
  stitched = function(df) {
    data.frame(stitch_x = df$stitch_x, mod_prob = df$mod_prob, lane = df$lane,
               stringsAsFactors = FALSE)
  }
  all_sites <- rbind(stitched(lf_sites), stitched(rf_sites), stitched(ins_locus))

  # --- 5. Build stitched x-axis break labels ---
  sep1 <- anchor          # left flank / middle boundary
  sep2 <- anchor + mid_width  # middle / right flank boundary

  # --- 6. Build read panel plot ---
  n_lanes <- max(reads_packed$lane)
  y_limits <- c(0.5, n_lanes + 0.5)

  p_reads <- ggplot2::ggplot() +
    ggplot2::geom_polygon(
      data = poly_data,
      ggplot2::aes(x = .data$x, y = .data$y,
                   group = .data$polygon_id, fill = .data$fill_col),
      colour = NA
    ) +
    ggplot2::scale_fill_identity() +
    ggplot2::geom_vline(xintercept = c(sep1, sep2),
                        linetype = "dashed", colour = "grey50", linewidth = 0.5) +
    ggplot2::coord_cartesian(ylim = y_limits, expand = FALSE) +
    theme_ggmethylation() +
    ggplot2::theme(
      axis.title  = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      panel.grid  = ggplot2::element_blank()
    )

  if (nrow(all_sites) > 0L) {
    p_reads <- p_reads +
      ggplot2::geom_segment(
        data = all_sites,
        ggplot2::aes(
          x      = .data$stitch_x,
          xend   = .data$stitch_x,
          y      = .data$lane - half_height,
          yend   = .data$lane + half_height,
          colour = .data$mod_prob
        ),
        linewidth = 0.5
      ) +
      ggplot2::scale_colour_gradient(
        low = colour_low, high = colour_high, limits = c(0, 1),
        name = "Mod. prob."
      )
  }

  # Overflow markers for carriers with ins_length > mid_width
  overflow_segs <- segs_all[segs_all$overflow, , drop = FALSE]
  if (nrow(overflow_segs) > 0L) {
    overflow_pts <- data.frame(
      stitch_x = overflow_segs$end,
      lane     = overflow_segs$lane
    )
    p_reads <- p_reads +
      ggplot2::geom_text(
        data = overflow_pts,
        ggplot2::aes(x = .data$stitch_x, y = .data$lane),
        label = ">", size = 3, colour = "black"
      )
  }

  # Separator line between carrier and non-carrier blocks
  if (!is.null(separator_lane)) {
    p_reads <- p_reads +
      ggplot2::geom_hline(yintercept = separator_lane,
                          linetype = "solid", colour = "grey70", linewidth = 0.5)
  }

  if (!show_smoothed) return(p_reads)

  # --- 7. Smoothed bottom panel ---
  # Compute smooth curves for carriers and non-carriers in left/right flanks,
  # and for carriers in the insertion middle.

  smooth_piece <- function(x_vals, y_vals, group_label, region) {
    if (length(x_vals) < 4L) return(NULL)
    df <- .smooth_xy(x_vals, y_vals)
    df$group  <- group_label
    df$region <- region
    df
  }
  carriers_in    <- function(df) df[df$read_name %in% carrier_names, , drop = FALSE]
  noncarriers_in <- function(df) df[df$read_name %in% noncarrier_names, , drop = FALSE]
  lf_c <- carriers_in(lf_sites)
  rf_c <- carriers_in(rf_sites)
  lf_nc <- noncarriers_in(lf_sites)
  rf_nc <- noncarriers_in(rf_sites)

  # Carriers and non-carriers in each flank; carriers only in the insertion
  # (at literal ins_offset, shifted onto the stitched axis).
  smooth_list <- Filter(Negate(is.null), list(
    smooth_piece(lf_c$position, lf_c$mod_prob, "carrier", "left"),
    if (include_noncarriers) smooth_piece(lf_nc$position, lf_nc$mod_prob, "non-carrier", "left"),
    smooth_piece(anchor + ins_locus$ins_offset - 1L, ins_locus$mod_prob, "carrier", "middle"),
    smooth_piece(rf_c$position + mid_width, rf_c$mod_prob, "carrier", "right"),
    if (include_noncarriers) smooth_piece(rf_nc$position + mid_width, rf_nc$mod_prob, "non-carrier", "right")
  ))

  if (length(smooth_list) == 0L) {
    # No smoothing data; return read panel only
    return(p_reads)
  }

  smooth_df <- do.call(rbind, smooth_list)
  rownames(smooth_df) <- NULL

  # Insert NA break rows at region boundaries to prevent line bridging
  groups_present <- unique(smooth_df$group)
  breaks_list <- lapply(groups_present, function(g) {
    data.frame(
      position  = c(sep1, sep2),
      mean_prob = c(NA_real_, NA_real_),
      group     = g,
      region    = c("break", "break"),
      stringsAsFactors = FALSE
    )
  })
  smooth_df <- rbind(smooth_df, do.call(rbind, breaks_list))
  smooth_df <- smooth_df[order(smooth_df$group, smooth_df$position), , drop = FALSE]

  group_colours <- c("carrier" = "#0072B2", "non-carrier" = "#999999")

  p_smooth <- ggplot2::ggplot(
    smooth_df,
    ggplot2::aes(x = .data$position, y = .data$mean_prob,
                 colour = .data$group, group = interaction(.data$group, .data$region))
  ) +
    ggplot2::geom_line(linewidth = 0.8, na.rm = TRUE) +
    ggplot2::geom_vline(xintercept = c(sep1, sep2),
                        linetype = "dashed", colour = "grey50", linewidth = 0.5) +
    ggplot2::scale_colour_manual(values = group_colours, name = NULL) +
    ggplot2::scale_y_continuous(limits = c(0, 1), name = "Mean mod. prob.") +
    theme_ggmethylation() +
    ggplot2::theme(
      axis.title.x  = ggplot2::element_blank(),
      axis.text.x   = ggplot2::element_blank(),
      axis.ticks.x  = ggplot2::element_blank(),
      panel.grid    = ggplot2::element_blank()
    )

  patchwork::wrap_plots(p_reads, p_smooth, ncol = 1L,
                        heights = c(3L, 1L)) &
    ggplot2::theme(legend.position = "right")
}
