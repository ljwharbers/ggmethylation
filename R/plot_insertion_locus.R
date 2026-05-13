# Per-locus insertion visualization for ggmethylation

#' Plot modification probabilities at an insertion locus
#'
#' Produces a composite ggplot / patchwork figure for a single insertion locus,
#' showing carrier and non-carrier reads in a stitched coordinate system
#' (left reference flank | insertion bases in literal bp | right reference flank),
#' with an optional loess-smoothed comparison panel below.
#'
#' Use [list_insertion_loci()] to discover available `locus_id` values.
#'
#' @param m A `methylation_data` object produced by [read_methylation()].
#' @param locus_id Character. A `locus_id` value from [list_insertion_loci()].
#' @param flank Integer. Number of reference base pairs to show on each side of
#'   the insertion (default 200).
#' @param show_smoothed Logical. Whether to add a loess-smoothed modification
#'   probability panel below the read panel (default `TRUE`).
#' @param include_noncarriers Logical. Whether to include reads that span the
#'   locus but do not carry the insertion (default `TRUE`). Non-carriers are
#'   shown below carriers with an empty middle region.
#' @param tol_pos Integer. Position tolerance used for re-clustering (must match
#'   the value used when calling [list_insertion_loci()]; default 10).
#' @param tol_len Numeric. Length tolerance used for re-clustering (default 0.20).
#' @param colour_low Character. Colour for low modification probability
#'   (default `"#313695"`).
#' @param colour_high Character. Colour for high modification probability
#'   (default `"#A50026"`).
#'
#' @return A `patchwork` object (or a single ggplot when `show_smoothed = FALSE`).
#'
#' @export
plot_insertion_locus <- function(m, locus_id,
                                  flank              = 200L,
                                  show_smoothed      = TRUE,
                                  include_noncarriers = TRUE,
                                  tol_pos            = 10L,
                                  tol_len            = 0.20,
                                  colour_low         = "#313695",
                                  colour_high        = "#A50026") {
  if (!inherits(m, "methylation_data")) {
    stop("'m' must be a methylation_data object.", call. = FALSE)
  }

  # --- 1. Find the requested locus ---
  loci <- list_insertion_loci(m, tol_pos = tol_pos, tol_len = tol_len,
                               min_reads = 1L)
  locus_row <- loci[loci$locus_id == locus_id, , drop = FALSE]
  if (nrow(locus_row) == 0L) {
    stop("locus_id '", locus_id, "' not found. ",
         "Run list_insertion_loci() to see available IDs.", call. = FALSE)
  }

  anchor      <- locus_row$anchor_pos
  mid_width   <- locus_row$median_length
  flank       <- as.integer(flank)

  # Carrier read names: rows in cigar_features[type=="I"] whose ref_start is
  # within tol_pos of anchor and whose length is within tol_len of mid_width
  cf_ins <- m$cigar_features[m$cigar_features$type == "I" &
                              !is.na(m$cigar_features$ref_start), , drop = FALSE]
  carrier_mask <- abs(cf_ins$ref_start - anchor) <= tol_pos &
                  abs(cf_ins$length - mid_width) / pmax(mid_width, 1L) <= tol_len
  carrier_cf       <- cf_ins[carrier_mask, , drop = FALSE]
  carrier_names    <- unique(carrier_cf$read_name)

  # Build per-carrier lookup: read_name -> insertion length (use first I-row per read)
  carrier_ins_len <- setNames(
    vapply(carrier_names, function(rn) {
      rows <- carrier_cf[carrier_cf$read_name == rn, , drop = FALSE]
      rows$length[1L]
    }, integer(1L)),
    carrier_names
  )

  reads_all <- m$reads
  all_names <- reads_all$read_name

  # Non-carrier reads: span [anchor - tol_pos, anchor + tol_pos], not in carrier set
  read_start <- pmin(reads_all$bam_pos, reads_all$start)
  read_end   <- reads_all$end
  spans_locus <- read_start <= (anchor + tol_pos) & read_end >= (anchor - tol_pos)
  noncarrier_mask_all <- spans_locus & !(all_names %in% carrier_names)
  noncarrier_names <- all_names[noncarrier_mask_all]

  carrier_reads    <- reads_all[all_names %in% carrier_names, , drop = FALSE]
  noncarrier_reads <- if (include_noncarriers) {
    reads_all[noncarrier_mask_all, , drop = FALSE]
  } else {
    reads_all[integer(0), , drop = FALSE]
  }

  if (nrow(carrier_reads) == 0L) {
    stop("No carrier reads found for locus '", locus_id, "'.", call. = FALSE)
  }

  # --- 2. Assign lanes via pack_reads ---
  carrier_reads$carrier_status    <- "carrier"
  noncarrier_reads$carrier_status <- "non-carrier"

  # Pack carriers and non-carriers separately so carriers come first (lower lane)
  carrier_packed    <- pack_reads(carrier_reads)
  n_carrier_lanes   <- max(carrier_packed$lane)

  if (nrow(noncarrier_reads) > 0L) {
    nc_packed       <- pack_reads(noncarrier_reads)
    nc_packed$lane  <- nc_packed$lane + n_carrier_lanes + 1L
    reads_packed    <- rbind(carrier_packed, nc_packed)
    separator_lane  <- n_carrier_lanes + 0.5
  } else {
    reads_packed    <- carrier_packed
    separator_lane  <- NULL
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
  carrier_fill    <- "#5E81AC"
  noncarrier_fill <- "#AAAAAA"
  poly_data$fill_col <- ifelse(
    poly_data$carrier_status == "carrier", carrier_fill, noncarrier_fill
  )

  # --- 4. Build mod-site data in stitched coords ---
  sites <- m$sites
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
  ins_sites <- m$insertion_sites
  ins_locus <- if (nrow(ins_sites) > 0L) {
    ins_sites[ins_sites$read_name %in% carrier_names &
              abs(ins_sites$ref_anchor - anchor) <= tol_pos, , drop = FALSE]
  } else {
    ins_sites[integer(0), , drop = FALSE]
  }
  ins_locus <- merge(ins_locus, lane_map, by = "read_name")
  ins_locus$stitch_x <- anchor + ins_locus$ins_offset - 1L

  # Combine all sites
  all_sites <- rbind(
    if (nrow(lf_sites) > 0L)
      data.frame(stitch_x = lf_sites$stitch_x, mod_prob = lf_sites$mod_prob,
                 lane = lf_sites$lane, stringsAsFactors = FALSE)
    else
      data.frame(stitch_x = numeric(0), mod_prob = numeric(0),
                 lane = numeric(0), stringsAsFactors = FALSE),
    if (nrow(rf_sites) > 0L)
      data.frame(stitch_x = rf_sites$stitch_x, mod_prob = rf_sites$mod_prob,
                 lane = rf_sites$lane, stringsAsFactors = FALSE)
    else
      data.frame(stitch_x = numeric(0), mod_prob = numeric(0),
                 lane = numeric(0), stringsAsFactors = FALSE),
    if (nrow(ins_locus) > 0L)
      data.frame(stitch_x = ins_locus$stitch_x, mod_prob = ins_locus$mod_prob,
                 lane = ins_locus$lane, stringsAsFactors = FALSE)
    else
      data.frame(stitch_x = numeric(0), mod_prob = numeric(0),
                 lane = numeric(0), stringsAsFactors = FALSE)
  )

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
    ggplot2::theme_minimal() +
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

  smooth_list <- list()

  .add_smooth <- function(x_vals, y_vals, group_label, region) {
    if (length(x_vals) < 4L) return(invisible(NULL))
    df <- .smooth_xy(x_vals, y_vals)
    df$group  <- group_label
    df$region <- region
    smooth_list[[length(smooth_list) + 1L]] <<- df
  }

  carrier_read_names    <- carrier_names
  noncarrier_read_names <- noncarrier_names

  # Left flank — carriers
  lf_c <- lf_sites[lf_sites$read_name %in% carrier_read_names, , drop = FALSE]
  if (nrow(lf_c) > 0L) .add_smooth(lf_c$position, lf_c$mod_prob, "carrier", "left")

  # Left flank — non-carriers
  if (include_noncarriers) {
    lf_nc <- lf_sites[lf_sites$read_name %in% noncarrier_read_names, , drop = FALSE]
    if (nrow(lf_nc) > 0L) .add_smooth(lf_nc$position, lf_nc$mod_prob, "non-carrier", "left")
  }

  # Insertion middle — carriers only
  if (nrow(ins_locus) > 0L) {
    # Use literal ins_offset as x (1..ins_length), shifted to stitched axis
    ins_c <- ins_locus
    .add_smooth(anchor + ins_c$ins_offset - 1L, ins_c$mod_prob, "carrier", "middle")
  }

  # Right flank — carriers
  rf_c <- rf_sites[rf_sites$read_name %in% carrier_read_names, , drop = FALSE]
  if (nrow(rf_c) > 0L)
    .add_smooth(rf_c$position + mid_width, rf_c$mod_prob, "carrier", "right")

  # Right flank — non-carriers
  if (include_noncarriers) {
    rf_nc <- rf_sites[rf_sites$read_name %in% noncarrier_read_names, , drop = FALSE]
    if (nrow(rf_nc) > 0L)
      .add_smooth(rf_nc$position + mid_width, rf_nc$mod_prob, "non-carrier", "right")
  }

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

  group_colours <- c("carrier" = "#5E81AC", "non-carrier" = "#AAAAAA")

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
    ggplot2::theme_minimal() +
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
