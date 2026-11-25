# Per-locus insertion visualization for ggmethylation

# Read fill (read panel) and line colour (smooth panel) per carrier status.
.CARRIER_COLOURS = c("carrier" = "#0072B2", "non-carrier" = "#999999")

#' Plot modification probabilities at an insertion locus
#'
#' Produces a composite ggplot / patchwork figure for a single insertion locus,
#' showing carrier and non-carrier reads in a stitched coordinate system
#' (left reference flank | insertion bases in literal bp | right reference flank),
#' with an optional loess-smoothed comparison panel below.
#'
#' With anchor `A` (`locus$anchor_pos`, the first reference base after the
#' insertion) and `W = locus$median_length`, the stitched x axis has one
#' column per base: reference bases `A - flank` to `A - 1` keep their
#' coordinate, inserted bases 1 to `W` sit at `A` to `A + W - 1`, and reference
#' bases `A` to `A + flank - 1` are shifted right by `W`. The axis labels give
#' reference coordinates in the flanks and the inserted-base range in the
#' middle; dashed lines mark the two insertion boundaries. Inserted bases
#' beyond `W` (carriers with a longer insertion than the median) are not drawn;
#' such carriers are marked with `>` at the right boundary.
#'
#' Use [list_insertion_loci()] to find loci; pass one row of its output.
#'
#' @param data A `methylation_data` object produced by [read_methylation()].
#' @param locus One row of the data.frame returned by [list_insertion_loci()]
#'   for `data`, e.g. `loci[1, ]` or `loci[loci$locus_id == id, ]`. It carries
#'   the locus' carrier and non-carrier reads and the clustering tolerances.
#' @param flank Positive integer. Number of reference base pairs to show on
#'   each side of the insertion (default 200).
#' @param show_smoothed Logical. Whether to add a loess-smoothed modification
#'   probability panel below the read panel (default `TRUE`). It has one line
#'   per carrier status and region, and per modification code (as linetype)
#'   when `data` holds more than one.
#' @param include_noncarriers Logical. Whether to include reads that span the
#'   locus but do not carry the insertion (default `TRUE`). Non-carriers are
#'   shown below carriers, separated by a line, with an empty middle region.
#' @param colour_low Character. Colour for low modification probability
#'   (default `"#BDBDBD"`).
#' @param colour_high Character. Colour for high modification probability
#'   (default `"#C62828"`).
#'
#' @return A `patchwork` object (or a single ggplot when `show_smoothed = FALSE`
#'   or there are too few calls to smooth).
#'
#' @examples
#' \dontrun{
#' md   = read_methylation("sample.bam", "chr21:34500000-34510000")
#' loci = list_insertion_loci(md, tol_pos = 10L, tol_len = 0.20, min_reads = 2L)
#' plot_insertion_locus(md, loci[1, ])
#' }
#'
#' @export
plot_insertion_locus = function(data, locus,
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
  flank               = .ins_check_number(flank, "flank", min = 1, whole = TRUE)
  show_smoothed       = .ins_check_flag(show_smoothed, "show_smoothed")
  include_noncarriers = .ins_check_flag(include_noncarriers, "include_noncarriers")

  # --- 1. Carrier / non-carrier reads of the locus ---
  locus_id    = locus$locus_id
  anchor      = locus$anchor_pos
  mid_width   = locus$median_length
  tol_pos     = locus$tol_pos

  carrier_ins_len  = locus$carriers[[1L]]   # insertion length, named by read
  carrier_names    = names(carrier_ins_len)
  noncarrier_names = if (include_noncarriers) locus$noncarriers[[1L]] else character(0L)

  reads_all = data$reads
  all_names = reads_all$read_name
  carrier_reads    = reads_all[all_names %in% carrier_names, , drop = FALSE]
  noncarrier_reads = reads_all[all_names %in% noncarrier_names, , drop = FALSE]

  if (nrow(carrier_reads) == 0L) {
    stop("No carrier reads found for locus '", locus_id, "'.", call. = FALSE)
  }

  # --- 2. Stitched coordinates ---
  # One x column per base:
  #   left flank   reference A - flank .. A - 1  ->  x = position
  #   insertion    inserted base 1 .. W          ->  x = A + ins_offset - 1
  #   right flank  reference A .. A + flank - 1  ->  x = position + W
  # where A = anchor (first reference base after the insertion) and
  # W = mid_width. Region boundaries lie half-way between columns.
  sep1  = anchor - 0.5               # left flank | insertion
  sep2  = anchor + mid_width - 0.5   # insertion | right flank
  x_lim = c(anchor - flank - 0.5, anchor + mid_width + flank - 0.5)

  # --- 3. Lanes: carriers on top, then a blank lane, then non-carriers ---
  carrier_reads$carrier_status    = rep("carrier",     nrow(carrier_reads))
  noncarrier_reads$carrier_status = rep("non-carrier", nrow(noncarrier_reads))
  carrier_reads$lane = pack_reads(carrier_reads)
  n_carrier_lanes    = max(carrier_reads$lane)
  separator_lane     = NULL
  reads_packed       = carrier_reads
  if (nrow(noncarrier_reads) > 0L) {
    noncarrier_reads$lane = pack_reads(noncarrier_reads) + n_carrier_lanes + 1L
    reads_packed          = rbind(carrier_reads, noncarrier_reads)
    separator_lane        = n_carrier_lanes + 1L
  }
  n_lanes = max(reads_packed$lane)

  # --- 4. Read polygons in stitched coordinates ---
  segs_all = do.call(rbind, lapply(seq_len(nrow(reads_packed)), function(i) {
    .stitched_segments(reads_packed[i, , drop = FALSE], carrier_ins_len,
                       anchor, mid_width, flank)
  }))
  rownames(segs_all) = NULL

  half_height = 0.4
  arrow_w     = 0.03 * (2 * flank + mid_width)
  poly_data   = .make_read_polygons(segs_all, arrow_w, half_height)

  # --- 5. Modification calls in stitched coordinates ---
  lane_map = reads_packed[, c("read_name", "lane", "carrier_status"), drop = FALSE]
  sites = data$sites

  lf_sites = sites[sites$position >= (anchor - flank) & sites$position < anchor, ,
                   drop = FALSE]
  lf_sites = merge(lf_sites, lane_map, by = "read_name")
  lf_sites$stitch_x = lf_sites$position
  lf_sites$region   = rep("left", nrow(lf_sites))

  rf_sites = sites[sites$position >= anchor & sites$position < (anchor + flank), ,
                   drop = FALSE]
  rf_sites = merge(rf_sites, lane_map, by = "read_name")
  rf_sites$stitch_x = rf_sites$position + mid_width
  rf_sites$region   = rep("right", nrow(rf_sites))

  # The carriers' own inserted bases; those beyond the drawn width are clipped.
  ins_locus = .carrier_insertion_sites(data$insertion_sites, carrier_ins_len,
                                       anchor, tol_pos)
  ins_locus = ins_locus[ins_locus$ins_offset <= mid_width, , drop = FALSE]
  ins_locus = merge(ins_locus, lane_map, by = "read_name")
  ins_locus$stitch_x = anchor + ins_locus$ins_offset - 1L
  ins_locus$region   = rep("middle", nrow(ins_locus))

  site_cols = c("stitch_x", "mod_prob", "mod_code", "lane", "carrier_status", "region")
  all_sites = rbind(lf_sites[site_cols], rf_sites[site_cols], ins_locus[site_cols])

  # --- 6. Shared x axis ---
  # Both panels get the same scale and xlim (with no expansion), so the
  # boundary lines line up when the panels are stacked.
  axis_breaks = .stitched_axis(anchor, mid_width, flank)
  x_axis = list(
    ggplot2::scale_x_continuous(breaks = axis_breaks$breaks,
                                labels = axis_breaks$labels,
                                expand = c(0, 0)),
    ggplot2::geom_vline(xintercept = c(sep1, sep2),
                        linetype = "dashed", colour = "grey50", linewidth = 0.5)
  )

  # --- 7. Read panel ---
  p_reads = ggplot2::ggplot() +
    ggplot2::geom_polygon(
      data = poly_data,
      ggplot2::aes(x = .data$x, y = .data$y,
                   group = .data$polygon_id, fill = .data$carrier_status),
      colour = NA
    ) +
    ggplot2::scale_fill_manual(values = .CARRIER_COLOURS, guide = "none") +
    x_axis +
    ggplot2::scale_y_reverse(expand = c(0, 0)) +
    ggplot2::coord_cartesian(xlim = x_lim, ylim = c(0.5, n_lanes + 0.5)) +
    theme_ggmethylation() +
    ggplot2::theme(
      axis.title   = ggplot2::element_blank(),
      axis.text.y  = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      panel.grid   = ggplot2::element_blank()
    )

  if (nrow(all_sites) > 0L) {
    p_reads = p_reads +
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

  # Overflow markers for carriers with an insertion longer than the drawn
  # width, just inside the insertion's right boundary.
  overflow_lanes = segs_all$lane[segs_all$overflow]
  if (length(overflow_lanes) > 0L) {
    p_reads = p_reads +
      ggplot2::geom_text(
        data = data.frame(x = sep2, lane = overflow_lanes),
        ggplot2::aes(x = .data$x, y = .data$lane),
        label = ">", hjust = 1, size = 3, colour = "black"
      )
  }

  # Separator line, centred in the blank lane between carriers and non-carriers
  if (!is.null(separator_lane)) {
    p_reads = p_reads +
      ggplot2::geom_hline(yintercept = separator_lane,
                          linetype = "solid", colour = "grey70", linewidth = 0.5)
  }

  if (!show_smoothed) return(p_reads)

  # --- 8. Smoothed bottom panel ---
  # One loess line per carrier status x region x mod code; non-carriers have
  # no middle region. Lines never cross a region boundary.
  pieces = split(all_sites, list(all_sites$carrier_status, all_sites$region,
                                 all_sites$mod_code), drop = TRUE)
  smooth_list = lapply(pieces, function(df) {
    if (nrow(df) < 4L) return(NULL)
    out = .smooth_xy(df$stitch_x, df$mod_prob)
    out$carrier_status = df$carrier_status[1L]
    out$region         = df$region[1L]
    out$mod_code       = df$mod_code[1L]
    out
  })
  smooth_list = Filter(Negate(is.null), smooth_list)

  if (length(smooth_list) == 0L) return(p_reads)

  smooth_df = do.call(rbind, smooth_list)
  rownames(smooth_df) = NULL

  # Codes in data order, as plot_methylation() assigns linetypes.
  codes = unique(c(data$sites$mod_code, all_sites$mod_code))
  codes = codes[codes %in% all_sites$mod_code]
  multi_code = length(codes) > 1L
  p_smooth = ggplot2::ggplot(smooth_df, ggplot2::aes(
    x = .data$position, y = .data$mean_prob, colour = .data$carrier_status,
    group = interaction(.data$carrier_status, .data$region, .data$mod_code)
  )) +
    ggplot2::geom_line(
      mapping = if (multi_code) ggplot2::aes(linetype = .data$mod_code),
      linewidth = 0.8, na.rm = TRUE
    ) +
    ggplot2::scale_colour_manual(values = .CARRIER_COLOURS, name = NULL) +
    x_axis +
    ggplot2::scale_y_continuous(limits = c(0, 1), name = "Mean mod. prob.") +
    ggplot2::coord_cartesian(xlim = x_lim) +
    theme_ggmethylation() +
    ggplot2::theme(
      axis.title.x  = ggplot2::element_blank(),
      panel.grid    = ggplot2::element_blank()
    )
  if (multi_code) {
    p_smooth = p_smooth + ggplot2::scale_linetype_manual(
      values = stats::setNames(
        rep_len(c("solid", "dashed", "dotdash", "dotted"), length(codes)), codes),
      name = "Modification"
    )
  }

  p_reads = p_reads + ggplot2::theme(
    axis.text.x  = ggplot2::element_blank(),
    axis.ticks.x = ggplot2::element_blank()
  )

  patchwork::wrap_plots(p_reads, p_smooth, ncol = 1L,
                        heights = c(3L, 1L)) &
    ggplot2::theme(legend.position = "right")
}

# Sub-segments of one read in stitched coordinates (see plot_insertion_locus()
# for the coordinate convention). Columns as expected by .make_read_polygons(),
# plus `overflow` (TRUE on a carrier's insertion segment when its insertion is
# longer than the drawn width).
.stitched_segments = function(r, carrier_ins_len, anchor, mid_width, flank) {
  r_start = min(r$bam_pos, r$start)
  r_end   = r$end
  is_carrier = r$carrier_status == "carrier"
  ins_len    = if (is_carrier) carrier_ins_len[[r$read_name]] else 0L

  # Reference bases lo..hi drawn from lo - 0.5 to hi + 0.5, shifted by `shift`.
  ref_seg = function(lo, hi, shift) {
    if (lo > hi) return(NULL)
    c(lo - 0.5 + shift, hi + 0.5 + shift)
  }
  segs = list(
    ref_seg(max(r_start, anchor - flank), min(r_end, anchor - 1L), 0),
    if (is_carrier) c(anchor - 0.5, anchor - 0.5 + min(ins_len, mid_width)),
    ref_seg(max(r_start, anchor), min(r_end, anchor + flank - 1L), mid_width)
  )
  overflow = c(FALSE, is_carrier && ins_len > mid_width, FALSE)
  keep = !vapply(segs, is.null, logical(1L))
  segs = segs[keep]
  overflow = overflow[keep]
  n = length(segs)
  if (n == 0L) return(NULL)

  data.frame(
    read_name        = r$read_name,
    start            = vapply(segs, `[`, numeric(1L), 1L),
    end              = vapply(segs, `[`, numeric(1L), 2L),
    lane             = r$lane,
    strand           = r$strand,
    carrier_status   = r$carrier_status,
    # Arrowheads only where the read itself ends inside the window.
    is_first_segment = seq_len(n) == 1L & r_start >= (anchor - flank),
    is_last_segment  = seq_len(n) == n & r_end <= (anchor + flank - 1L),
    overflow         = overflow,
    stringsAsFactors = FALSE
  )
}

# Breaks and labels for the stitched x axis: reference coordinates in each
# flank and one "+1...+W bp" label in the middle of the insertion. Flank breaks
# too close to the insertion label, or to the window edges (where the panel has
# no expansion and the label would be cut off), are dropped.
.stitched_axis = function(anchor, mid_width, flank) {
  total    = 2 * flank + mid_width
  edge_gap = 0.03 * total
  flank_breaks = function(lo, hi) {
    b = scales::extended_breaks(n = 4L)(c(lo, hi))
    b[b >= lo & b <= hi]
  }
  left  = flank_breaks(anchor - flank, anchor - 1L)
  right = flank_breaks(anchor, anchor + flank - 1L)
  left  = left[left - (anchor - flank) >= edge_gap]
  right = right[(anchor + flank - 1L) - right >= edge_gap]
  mid   = anchor + (mid_width - 1) / 2
  min_gap = 0.15 * total
  left  = left[abs(left - mid) >= min_gap]
  right = right[abs(right + mid_width - mid) >= min_gap]
  list(
    breaks = c(left, mid, right + mid_width),
    labels = c(scales::comma(left),
               sprintf("+1\u2026+%d bp", as.integer(mid_width)),
               scales::comma(right))
  )
}
