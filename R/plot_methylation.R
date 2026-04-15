# Internal helpers for interrupting the smooth line at consensus deletion regions

.ordered_plot_groups <- function(groups) {
  ordered <- sort(unique(groups[!is.na(groups)]))
  if (any(is.na(groups))) {
    c(ordered, NA_character_)
  } else {
    ordered
  }
}

.match_plot_group <- function(values, group) {
  if (is.na(group)) {
    is.na(values)
  } else {
    !is.na(values) & values == group
  }
}

# For each group, returns genomic intervals [del_start, del_end] where at least
# `threshold` fraction of reads in that group carry a deletion.
# cigar_features: pre-filtered to type == "D" and min_indel_size by caller
# reads: data.frame with read_name and the group column
.consensus_deletion_ranges <- function(cigar_features, reads, group_col,
                                       threshold = 0.75) {
  empty <- data.frame(
    del_start = integer(0L),
    del_end   = integer(0L),
    stringsAsFactors = FALSE
  )
  empty[[group_col]] <- character(0L)

  dels <- cigar_features[cigar_features$type == "D", , drop = FALSE]
  if (nrow(dels) == 0L) return(empty)

  # Attach group membership to each deletion row
  dels <- merge(
    dels,
    reads[, c("read_name", group_col), drop = FALSE],
    by = "read_name", all.x = FALSE
  )
  if (nrow(dels) == 0L) return(empty)

  groups <- unique(dels[[group_col]])
  groups <- groups[!is.na(groups)]
  result_list <- vector("list", length(groups))

  for (k in seq_along(groups)) {
    grp  <- groups[k]
    grp_mask <- .match_plot_group(dels[[group_col]], grp)
    sub  <- dels[grp_mask, , drop = FALSE]
    n_reads <- sum(.match_plot_group(reads[[group_col]], grp))
    if (n_reads == 0L) next

    sub <- sub[order(sub$ref_start), , drop = FALSE]

    # Sweep-merge overlapping deletion intervals, accumulating read names
    merged_starts <- integer(0L)
    merged_ends   <- integer(0L)
    merged_reads  <- list()

    cur_start <- sub$ref_start[1L]
    cur_end   <- sub$ref_end[1L]
    cur_rds   <- sub$read_name[1L]

    for (i in seq_len(nrow(sub))[-1L]) {
      if (sub$ref_start[i] <= cur_end + 1L) {
        cur_end <- max(cur_end, sub$ref_end[i])
        cur_rds <- c(cur_rds, sub$read_name[i])
      } else {
        merged_starts <- c(merged_starts, cur_start)
        merged_ends   <- c(merged_ends,   cur_end)
        merged_reads  <- c(merged_reads,  list(unique(cur_rds)))
        cur_start <- sub$ref_start[i]
        cur_end   <- sub$ref_end[i]
        cur_rds   <- sub$read_name[i]
      }
    }
    merged_starts <- c(merged_starts, cur_start)
    merged_ends   <- c(merged_ends,   cur_end)
    merged_reads  <- c(merged_reads,  list(unique(cur_rds)))

    keep <- vapply(
      merged_reads,
      function(rds) length(rds) / n_reads >= threshold,
      logical(1L)
    )

    if (!any(keep)) next

    df <- data.frame(
      del_start = merged_starts[keep],
      del_end   = merged_ends[keep],
      stringsAsFactors = FALSE
    )
    df[[group_col]] <- grp
    result_list[[k]] <- df
  }

  out <- do.call(rbind, Filter(Negate(is.null), result_list))
  if (is.null(out)) return(empty)
  rownames(out) <- NULL
  out
}

# Insert NA breaks into a smoothed data frame at each consensus deletion interval.
# Sentinel NA rows are added at del_start - 0.5 and del_end + 0.5 to guarantee
# a visible gap even when no grid point falls inside the deletion.
.insert_deletion_breaks <- function(smoothed, deletion_ranges, group_col) {
  if (nrow(deletion_ranges) == 0L) return(smoothed)

  id_cols <- setdiff(names(smoothed), c("position", "mean_prob"))
  sentinel_list <- vector("list", nrow(deletion_ranges))

  for (i in seq_len(nrow(deletion_ranges))) {
    grp_val   <- deletion_ranges[[group_col]][i]
    del_start <- deletion_ranges$del_start[i]
    del_end   <- deletion_ranges$del_end[i]

    # Mask grid points within the deletion
    in_grp <- smoothed[[group_col]] == grp_val
    in_del <- smoothed$position >= del_start & smoothed$position <= del_end
    smoothed$mean_prob[in_grp & in_del] <- NA_real_

    # Build sentinel rows — one pair per unique line-identity combo in this group
    grp_rows  <- smoothed[in_grp, id_cols, drop = FALSE]
    templates <- unique(grp_rows)

    if (nrow(templates) == 0L) next

    sentinels <- vector("list", nrow(templates) * 2L)
    for (j in seq_len(nrow(templates))) {
      s1 <- templates[j, , drop = FALSE]
      s1$position  <- del_start - 0.5
      s1$mean_prob <- NA_real_
      s2 <- templates[j, , drop = FALSE]
      s2$position  <- del_end + 0.5
      s2$mean_prob <- NA_real_
      sentinels[[2L * j - 1L]] <- s1
      sentinels[[2L * j]]      <- s2
    }
    sentinel_list[[i]] <- do.call(rbind, sentinels)
  }

  all_sentinels <- do.call(rbind, Filter(Negate(is.null), sentinel_list))
  if (!is.null(all_sentinels) && nrow(all_sentinels) > 0L) {
    smoothed <- rbind(smoothed, all_sentinels[, names(smoothed), drop = FALSE])
  }

  smoothed <- smoothed[order(smoothed[[group_col]], smoothed$position), , drop = FALSE]
  rownames(smoothed) <- NULL
  smoothed
}

# Convenience wrapper: apply deletion breaks when show_cigar is TRUE.
.apply_deletion_breaks <- function(smoothed, cigar_features, reads, group_col,
                                   min_indel_size, show_cigar) {
  if (!isTRUE(show_cigar)) return(smoothed)
  if (is.null(cigar_features) || nrow(cigar_features) == 0L) return(smoothed)
  dels <- cigar_features[
    cigar_features$type == "D" & cigar_features$length >= min_indel_size,
    , drop = FALSE
  ]
  if (nrow(dels) == 0L) return(smoothed)
  ranges <- .consensus_deletion_ranges(dels, reads, group_col)
  .insert_deletion_breaks(smoothed, ranges, group_col)
}

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
#' @param line_width Linewidth of modification site lines (default 0.2).
#' @param colour_strand Logical. When `TRUE`, read bars are coloured by strand
#'   (`"+"` = forward, `"-"` = reverse). Ignored when data is grouped; group
#'   colour takes precedence. Default `FALSE`.
#' @param strand_colours Named character vector with `"+"` and `"-"` entries
#'   giving bar colours for each strand. Only used when `colour_strand = TRUE`
#'   and data is ungrouped. Default `c("+" = "#4393C3", "-" = "#D6604D")`.
#' @param group_colours Named character vector of colours per group. Defaults
#'   to `c("1" = "#95babc", "2" = "#efbb76")`, matching the typical HP
#'   haplotype tag output. Pass `NULL` to use ggplot2 defaults, or supply a
#'   fully named vector for other group names.
#' @param smooth_span Loess smoothing span for the bottom panel. When `NULL`
#'   (default), an adaptive span is computed per group as
#'   `max(0.15, min(0.75, 15 / n_unique_sites))`, targeting ~15 data points
#'   per local fit regardless of region size or CpG density. Pass an explicit
#'   numeric value (e.g. `0.3`) to use a fixed span.
#' @param panel_heights Numeric vector giving relative heights of the panels.
#'   When `annotations = NULL`, expects length 2 (`reads`, `smooth`); when
#'   annotations are provided, expects length 3 (`gene track`, `reads`,
#'   `smooth`) — the gene panel is placed at the top, reads in the middle, and
#'   the smoothed probability curve at the bottom. Pass `NULL` to use the
#'   built-in defaults.
#' @param annotations A `gene_annotations` object returned by
#'   [read_annotations()], or `NULL` (default). When provided, a gene
#'   annotation track is placed at the top of the composite figure (above the
#'   read panel and the smooth panel).
#' @param variants A `variant_data` object returned by [read_variants()], or
#'   `NULL` (default). When provided, per-read base letters are drawn at
#'   variant positions (coloured by match to ref/alt), and vertical dashed
#'   lines are added at every variant position across all panels. Requires
#'   that `data` was produced by a current version of [read_methylation()] that
#'   stores sequences and CIGARs in the object.
#' @param mod_code_shapes Named integer vector mapping modification codes to
#'   point shapes used in the legend when multiple modification codes are
#'   present. Names must match the mod codes in `data$sites$mod_code`. If
#'   `NULL` (default), shapes are assigned automatically (circle, square,
#'   triangle, diamond for up to four codes).
#' @param show_cigar Logical. When `TRUE` (default), structural variants from
#'   CIGAR strings (insertions, deletions) are displayed on reads. Insertions
#'   appear as bold purple "I" markers; deletions open a gap in the read bar
#'   and draw a thin black line. Set to `FALSE` to hide these features.
#' @param min_indel_size Integer. Minimum size (in bp) for insertions and
#'   deletions to be displayed when `show_cigar = TRUE`. Indels smaller than
#'   this threshold are hidden to reduce visual clutter from common small
#'   indels. Default `50`.
#' @param show_supplementary Logical. When `TRUE` (default), a coloured halo
#'   is drawn around read bars indicating the chromosome of the supplementary
#'   alignment partner (from the SA BAM tag). The original bar colouring
#'   (group, strand, or default grey) is preserved inside the halo. Reads with
#'   no supplementary alignment have no halo.
#' @param bnd_match_tol Integer. Position tolerance (bp) for matching
#'   supplementary-alignment breakpoints to VCF BND calls. The SA matching runs
#'   whenever `variants` is supplied and reads carry SA tags; the visual border
#'   marking requires `show_supplementary = TRUE`. Default 50.
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
                             line_width = 0.2,
                             colour_strand = FALSE,
                             strand_colours = c("+" = "#4393C3", "-" = "#D6604D"),
                             group_colours = c("1" = "#95babc", "2" = "#efbb76"),
                             mod_code_shapes = NULL,
                             smooth_span = NULL,
                             panel_heights = NULL,
                             annotations = NULL,
                             variants = NULL,
                             show_cigar = TRUE,
                             min_indel_size = 50L,
                             show_supplementary = TRUE,
                             bnd_match_tol = 50L) {
  # --- 1. Validate input ---
  if (inherits(data, "multi_methylation_data")) {
    return(.plot_multi_methylation(
      data               = data,
      sort_by            = sort_by,
      colour_low         = colour_low,
      colour_high        = colour_high,
      line_width         = line_width,
      colour_strand      = colour_strand,
      strand_colours     = strand_colours,
      group_colours      = group_colours,
      mod_code_shapes    = mod_code_shapes,
      smooth_span        = smooth_span,
      panel_heights      = panel_heights,
      annotations        = annotations,
      variants           = variants,
      show_cigar         = show_cigar,
      min_indel_size     = min_indel_size,
      show_supplementary = show_supplementary,
      bnd_match_tol      = bnd_match_tol
    ))
  }

  if (!inherits(data, "methylation_data") && !inherits(data, "multi_methylation_data")) {
    stop(
      "`data` must be a `methylation_data` or `multi_methylation_data` object.",
      call. = FALSE
    )
  }
  if (!is.logical(colour_strand) || length(colour_strand) != 1L) {
    stop("'colour_strand' must be a single logical value.", call. = FALSE)
  }
  if (isTRUE(colour_strand) && !all(c("+", "-") %in% names(strand_colours))) {
    stop("'strand_colours' must be a named vector with '+' and '-' entries.", call. = FALSE)
  }
  if (isTRUE(colour_strand) && !is.null(data$group_tag)) {
    message("'colour_strand = TRUE' is ignored when reads are grouped; group colour takes precedence.")
  }

  region_start <- GenomicRanges::start(data$region)
  region_end <- GenomicRanges::end(data$region)

  # --- 1b. Resolve mod_code shape mapping ---
  codes <- unique(data$sites$mod_code)
  default_shapes <- c(16L, 15L, 17L, 18L)   # circle, square, triangle, diamond

  if (length(codes) > length(default_shapes)) {
    stop("More than 4 mod codes present; supply mod_code_shapes explicitly.", call. = FALSE)
  }

  if (is.null(mod_code_shapes)) {
    mod_code_shapes <- stats::setNames(
      default_shapes[seq_along(codes)],
      codes
    )
  } else {
    # Validate that all codes present in data have an entry
    missing_codes <- setdiff(codes, names(mod_code_shapes))
    if (length(missing_codes) > 0L) {
      warning(
        "mod_code_shapes has no entry for code(s): ",
        paste(missing_codes, collapse = ", "),
        ". Using default shapes.", call. = FALSE
      )
      extras <- stats::setNames(
        default_shapes[seq_along(missing_codes)],
        missing_codes
      )
      mod_code_shapes <- c(mod_code_shapes, extras)
    }
  }

  multi_code <- length(codes) > 1L

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
    groups_ordered <- .ordered_plot_groups(data$reads$group)
    lane_offset <- 0L
    for (grp in groups_ordered) {
      idx <- which(.match_plot_group(data$reads$group, grp))
      data$reads$lane[idx] <- pack_reads(data$reads[idx, ], clip_side = data$reads$clip_side[idx]) + lane_offset
      lane_offset <- max(data$reads$lane[idx]) + 2L
      separator_lanes <- c(separator_lanes, lane_offset - 1L)
    }
    # Drop the trailing separator (after the last group)
    separator_lanes <- separator_lanes[-length(separator_lanes)]
  } else {
    data$reads$lane <- pack_reads(data$reads, clip_side = data$reads$clip_side)
  }

  # --- 5b. Build variant overlay ---
  if (!is.null(variants) && !inherits(variants, "variant_data")) {
    stop("`variants` must be a `variant_data` object returned by read_variants().",
         call. = FALSE)
  }
  variant_ov <- if (!is.null(variants)) {
    build_variant_overlay(data, variants, bnd_match_tol = bnd_match_tol)
  } else {
    NULL
  }

  # --- 6. Build top panel ---
  p_top <- build_read_panel(
    data               = data,
    separator_lanes    = separator_lanes,
    region_start       = region_start,
    region_end         = region_end,
    colour_low         = colour_low,
    colour_high        = colour_high,
    line_width         = line_width,
    colour_strand      = colour_strand,
    strand_colours     = strand_colours,
    group_colours      = group_colours,
    mod_code_shapes    = mod_code_shapes,
    show_x_axis        = FALSE,
    variant_overlay    = variant_ov,
    show_cigar         = show_cigar,
    cigar_features     = if (isTRUE(show_cigar)) data$cigar_features else NULL,
    min_indel_size     = min_indel_size,
    show_supplementary = show_supplementary
  )

  # --- 7. Build bottom panel ---
  default_linetypes <- c("solid", "dashed", "dotdash", "dotted")

  if (is.null(data$group_tag)) {
    # Ungrouped smooth panel
    sites_smooth <- data$sites
    sites_smooth$group <- "all"

    reads_grouped <- data$reads
    reads_grouped$group <- "all"

    if (!multi_code) {
      # Single code: one line, fixed colour
      smoothed <- smooth_methylation(sites_smooth, group_col = "group",
                                     span = smooth_span)
      smoothed <- .apply_deletion_breaks(smoothed, data$cigar_features,
                                         reads_grouped, "group",
                                         min_indel_size, show_cigar)

      p_bottom <- ggplot2::ggplot(
        smoothed,
        ggplot2::aes(x = .data$position, y = .data$mean_prob)
      ) +
        ggplot2::geom_line(linewidth = 1, colour = "#C62828") +
        ggplot2::scale_y_continuous(
          limits = c(0, 1),
          name = "Mean modification\nprobability"
        ) +
        ggplot2::scale_x_continuous(labels = scales::comma_format()) +
        ggplot2::coord_cartesian(xlim = c(region_start, region_end)) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
          panel.grid.minor = ggplot2::element_blank(),
          legend.text      = ggplot2::element_text(size = ggplot2::rel(0.75)),
          legend.title     = ggplot2::element_text(size = ggplot2::rel(0.75)),
          legend.key.size  = ggplot2::unit(0.4, "cm")
        ) +
        ggplot2::labs(x = "Genomic position (bp)")
    } else {
      # Multi-code: one line per code, colour by mod_code
      smoothed <- smooth_methylation(sites_smooth, group_col = "group",
                                     mod_code_col = "mod_code", span = smooth_span)
      smoothed <- .apply_deletion_breaks(smoothed, data$cigar_features,
                                         reads_grouped, "group",
                                         min_indel_size, show_cigar)

      p_bottom <- ggplot2::ggplot(
        smoothed,
        ggplot2::aes(
          x      = .data$position,
          y      = .data$mean_prob,
          colour = .data$mod_code
        )
      ) +
        ggplot2::geom_line(linewidth = 1) +
        ggplot2::scale_y_continuous(
          limits = c(0, 1),
          name = "Mean modification\nprobability"
        ) +
        ggplot2::scale_x_continuous(labels = scales::comma_format()) +
        ggplot2::coord_cartesian(xlim = c(region_start, region_end)) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
          panel.grid.minor = ggplot2::element_blank(),
          legend.text      = ggplot2::element_text(size = ggplot2::rel(0.75)),
          legend.title     = ggplot2::element_text(size = ggplot2::rel(0.75)),
          legend.key.size  = ggplot2::unit(0.4, "cm")
        ) +
        ggplot2::labs(x = "Genomic position (bp)", colour = "Modification")
    }
  } else {
    # Grouped smooth panel
    if (!multi_code) {
      # Single code: one line per group, colour by group
      smoothed <- smooth_methylation(
        data$sites,
        group_col = "group",
        span = smooth_span
      )
      smoothed <- .apply_deletion_breaks(smoothed, data$cigar_features,
                                         data$reads, "group",
                                         min_indel_size, show_cigar)

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
        ggplot2::scale_x_continuous(labels = scales::comma_format()) +
        ggplot2::coord_cartesian(xlim = c(region_start, region_end)) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
          panel.grid.minor = ggplot2::element_blank(),
          legend.text      = ggplot2::element_text(size = ggplot2::rel(0.75)),
          legend.title     = ggplot2::element_text(size = ggplot2::rel(0.75)),
          legend.key.size  = ggplot2::unit(0.4, "cm")
        ) +
        ggplot2::labs(x = "Genomic position (bp)", colour = "Group")

      if (!is.null(group_colours)) {
        p_bottom <- p_bottom +
          ggplot2::scale_colour_manual(values = group_colours, na.value = "grey50")
      }
    } else {
      # Multi-code + grouped: colour by group, linetype by mod_code
      smoothed <- smooth_methylation(
        data$sites,
        group_col    = "group",
        mod_code_col = "mod_code",
        span         = smooth_span
      )
      smoothed <- .apply_deletion_breaks(smoothed, data$cigar_features,
                                         data$reads, "group",
                                         min_indel_size, show_cigar)

      p_bottom <- ggplot2::ggplot(
        smoothed,
        ggplot2::aes(
          x        = .data$position,
          y        = .data$mean_prob,
          colour   = .data$group,
          linetype = .data$mod_code
        )
      ) +
        ggplot2::geom_line(linewidth = 1) +
        ggplot2::scale_y_continuous(
          limits = c(0, 1),
          name = "Mean modification\nprobability"
        ) +
        ggplot2::scale_x_continuous(labels = scales::comma_format()) +
        ggplot2::coord_cartesian(xlim = c(region_start, region_end)) +
        ggplot2::scale_linetype_manual(
          values = stats::setNames(
            default_linetypes[seq_along(codes)],
            codes
          ),
          name = "Modification"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
          panel.grid.minor = ggplot2::element_blank(),
          legend.text      = ggplot2::element_text(size = ggplot2::rel(0.75)),
          legend.title     = ggplot2::element_text(size = ggplot2::rel(0.75)),
          legend.key.size  = ggplot2::unit(0.4, "cm")
        ) +
        ggplot2::labs(x = "Genomic position (bp)", colour = "Group")

      if (!is.null(group_colours)) {
        p_bottom <- p_bottom +
          ggplot2::scale_colour_manual(values = group_colours, na.value = "grey50")
      }
    }
  }

  if (!is.null(data$snv_position)) {
    p_bottom <- p_bottom +
      ggplot2::geom_vline(
        xintercept = data$snv_position,
        linetype = "dashed", colour = "black", linewidth = 0.5
      )
  }

  # --- 8b. Build gene annotation panel (optional) ---
  p_gene <- NULL
  if (!is.null(annotations)) {
    if (!inherits(annotations, "gene_annotations")) {
      stop("`annotations` must be a `gene_annotations` object from read_annotations().",
           call. = FALSE)
    }
    p_gene <- build_gene_panel(
      annotations  = annotations,
      region_start = region_start,
      region_end   = region_end,
      bottom_label = FALSE
    )
  }

  # --- 9. Combine with patchwork ---
  # Panel order (when gene panel present): gene (top), reads, smooth (bottom
  # with genomic coordinates). When no gene panel: reads, smooth.
  panels   <- list(p_top, p_bottom)
  n_panels <- 2L

  if (!is.null(p_gene)) {
    panels   <- c(list(p_gene), panels)
    n_panels <- 3L
  }

  # Compute heights
  if (is.null(panel_heights)) {
    heights <- c(if (!is.null(p_gene)) 0.08 else NULL, 1, 0.25)
  } else {
    if (length(panel_heights) != n_panels) {
      stop(sprintf(
        "`panel_heights` has length %d but there are %d panels.",
        length(panel_heights), n_panels
      ), call. = FALSE)
    }
    heights <- panel_heights
  }

  patchwork::wrap_plots(panels, ncol = 1, heights = heights)
}

# Internal multi-sample renderer
# Not exported — called by plot_methylation() when data is multi_methylation_data.

.plot_multi_methylation <- function(data, sort_by, colour_low, colour_high,
                                    line_width, colour_strand, strand_colours,
                                    group_colours, mod_code_shapes,
                                    smooth_span, panel_heights, annotations,
                                    variants, show_cigar = FALSE,
                                    min_indel_size = 50L,
                                    show_supplementary = FALSE,
                                    bnd_match_tol = 50L) {

  region_start <- GenomicRanges::start(data$region)
  region_end   <- GenomicRanges::end(data$region)

  # Resolve mod_code shapes from combined codes across all samples
  all_codes <- unique(unlist(lapply(data$samples, function(s) {
    unique(s$sites$mod_code)
  })))
  default_shapes <- c(16L, 15L, 17L, 18L)

  if (is.null(mod_code_shapes)) {
    mod_code_shapes <- stats::setNames(
      default_shapes[seq_along(all_codes)],
      all_codes
    )
  } else {
    missing_codes <- setdiff(all_codes, names(mod_code_shapes))
    if (length(missing_codes) > 0L) {
      warning(
        "mod_code_shapes has no entry for code(s): ",
        paste(missing_codes, collapse = ", "),
        ". Using default shapes.", call. = FALSE
      )
      extras <- stats::setNames(
        default_shapes[seq_along(missing_codes)],
        missing_codes
      )
      mod_code_shapes <- c(mod_code_shapes, extras)
    }
  }

  # --- 1. Build per-sample read panels ---
  sample_names  <- names(data$samples)
  n_samples     <- length(sample_names)
  sample_panels <- vector("list", n_samples)

  # Combined sites list for shared smooth panel (built alongside)
  combined_sites_list <- vector("list", n_samples)

  for (i in seq_len(n_samples)) {
    nm <- sample_names[i]
    s  <- data$samples[[nm]]

    # Skip empty samples gracefully
    if (nrow(s$reads) == 0L) {
      sample_panels[[i]] <- ggplot2::ggplot() +
        ggplot2::annotate(
          "text", x = 0.5, y = 0.5,
          label = sprintf("No reads (%s)", nm), size = 4
        ) +
        ggplot2::theme_void() +
        ggplot2::labs(title = nm)
      combined_sites_list[[i]] <- NULL
      next
    }

    # Compute mean_mod_prob per read
    read_means <- stats::setNames(
      tapply(s$sites$mod_prob, s$sites$read_name, mean),
      NULL
    )
    s$reads$mean_mod_prob <- as.numeric(read_means[s$reads$read_name])
    s$reads$mean_mod_prob[is.na(s$reads$mean_mod_prob)] <- 0

    # Sort reads
    sort_by_s <- sort_by
    if (is.null(sort_by_s)) {
      if (is.null(s$group_tag)) {
        sort_by_s <- "start"
      } else {
        sort_by_s <- c("start", "group", "mean_mod_prob")
      }
    }
    sort_args <- lapply(sort_by_s, function(col) s$reads[[col]])
    ord <- do.call(order, sort_args)
    s$reads <- s$reads[ord, , drop = FALSE]
    rownames(s$reads) <- NULL

    # Pack reads into lanes
    s$reads$lane  <- integer(nrow(s$reads))
    separator_lanes <- numeric(0)

    if (!is.null(s$group_tag)) {
      groups_ordered <- .ordered_plot_groups(s$reads$group)
      lane_offset <- 0L
      for (grp in groups_ordered) {
        idx <- which(.match_plot_group(s$reads$group, grp))
        s$reads$lane[idx] <- pack_reads(s$reads[idx, ]) + lane_offset
        lane_offset <- max(s$reads$lane[idx]) + 2L
        separator_lanes <- c(separator_lanes, lane_offset - 1L)
      }
      if (length(separator_lanes) > 0L) {
        separator_lanes <- separator_lanes[-length(separator_lanes)]
      }
    } else {
      s$reads$lane <- pack_reads(s$reads)
    }

    # Build per-sample variant overlay
    if (!is.null(variants) && !inherits(variants, "variant_data")) {
      stop("`variants` must be a `variant_data` object returned by read_variants().",
           call. = FALSE)
    }
    s_variant_ov <- if (!is.null(variants)) {
      build_variant_overlay(s, variants, bnd_match_tol = bnd_match_tol)
    } else {
      NULL
    }

    # Build read panel
    p_reads <- build_read_panel(
      data               = s,
      separator_lanes    = separator_lanes,
      region_start       = region_start,
      region_end         = region_end,
      colour_low         = colour_low,
      colour_high        = colour_high,
      line_width         = line_width,
      colour_strand      = colour_strand,
      strand_colours     = strand_colours,
      group_colours      = group_colours,
      mod_code_shapes    = mod_code_shapes,
      show_x_axis        = FALSE,
      variant_overlay    = s_variant_ov,
      show_cigar         = show_cigar,
      cigar_features     = if (isTRUE(show_cigar)) s$cigar_features else NULL,
      min_indel_size     = min_indel_size,
      show_supplementary = show_supplementary
    )
    p_reads <- p_reads + ggplot2::labs(title = nm)
    sample_panels[[i]] <- p_reads

    # Collect sites for shared smooth panel, tagged with sample name
    sites_tagged <- s$sites
    n_sites <- nrow(sites_tagged)
    sites_tagged$sample <- if (n_sites > 0L) rep(nm, n_sites) else character(0L)
    if (!is.null(s$group_tag) && "group" %in% names(sites_tagged) && n_sites > 0L) {
      sites_tagged$smooth_key <- paste(nm, sites_tagged$group, sep = ":")
    } else {
      sites_tagged$smooth_key <- if (n_sites > 0L) rep(nm, n_sites) else character(0L)
    }
    combined_sites_list[[i]] <- if (n_sites > 0L) sites_tagged else NULL
  }

  # --- 2. Build shared smooth panel ---
  combined_sites <- do.call(rbind, Filter(Negate(is.null), combined_sites_list))
  rownames(combined_sites) <- NULL

  if (!is.null(combined_sites) && nrow(combined_sites) > 0L) {
    smoothed <- smooth_methylation(
      combined_sites,
      group_col = "smooth_key",
      span      = smooth_span
    )

    if (isTRUE(show_cigar)) {
      # Collect deletions and reads from all samples tagged with their smooth_key
      all_dels_list  <- vector("list", n_samples)
      all_reads_list <- vector("list", n_samples)
      for (i in seq_len(n_samples)) {
        nm <- sample_names[i]
        s  <- data$samples[[nm]]
        if (!is.null(s$cigar_features) && nrow(s$cigar_features) > 0L) {
          s_dels <- s$cigar_features[
            s$cigar_features$type == "D" &
              s$cigar_features$length >= min_indel_size, , drop = FALSE
          ]
          if (nrow(s_dels) > 0L) all_dels_list[[i]] <- s_dels
        }
        if (nrow(s$reads) > 0L) {
          s_reads <- s$reads
          if (!is.null(s$group_tag) && "group" %in% names(s_reads)) {
            s_reads$smooth_key <- paste(nm, s_reads$group, sep = ":")
          } else {
            s_reads$smooth_key <- rep(nm, nrow(s_reads))
          }
          all_reads_list[[i]] <- s_reads
        }
      }
      all_dels  <- do.call(rbind, Filter(Negate(is.null), all_dels_list))
      all_reads <- do.call(rbind, Filter(Negate(is.null), all_reads_list))
      if (!is.null(all_dels) && nrow(all_dels) > 0L &&
          !is.null(all_reads) && nrow(all_reads) > 0L) {
        ranges  <- .consensus_deletion_ranges(all_dels, all_reads, "smooth_key")
        smoothed <- .insert_deletion_breaks(smoothed, ranges, "smooth_key")
      }
    }

    # Determine if any sample has grouping by checking for ":" in smooth_key
    multi_has_groups <- any(grepl(":", smoothed$smooth_key, fixed = TRUE))

    if (multi_has_groups) {
      # Split "sampleName:groupValue" into separate columns
      parts <- strsplit(smoothed$smooth_key, ":", fixed = TRUE)
      smoothed$smooth_sample <- vapply(parts, `[[`, character(1L), 1L)
      smoothed$smooth_group  <- vapply(
        parts, function(x) paste(x[-1L], collapse = ":"), character(1L)
      )
      p_smooth <- ggplot2::ggplot(
        smoothed,
        ggplot2::aes(
          x        = .data$position,
          y        = .data$mean_prob,
          colour   = .data$smooth_group,
          linetype = .data$smooth_sample
        )
      ) +
        ggplot2::geom_line(linewidth = 1) +
        ggplot2::scale_y_continuous(
          limits = c(0, 1),
          name   = "Mean modification\nprobability"
        ) +
        ggplot2::scale_x_continuous(labels = scales::comma_format()) +
        ggplot2::coord_cartesian(xlim = c(region_start, region_end)) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
          panel.grid.minor = ggplot2::element_blank(),
          legend.text      = ggplot2::element_text(size = ggplot2::rel(0.75)),
          legend.title     = ggplot2::element_text(size = ggplot2::rel(0.75)),
          legend.key.size  = ggplot2::unit(0.4, "cm")
        ) +
        ggplot2::labs(x = "Genomic position (bp)", colour = "Group",
                      linetype = "Sample")

      if (!is.null(group_colours)) {
        p_smooth <- p_smooth +
          ggplot2::scale_colour_manual(values = group_colours, na.value = "grey50")
      }
    } else {
      # No grouping: colour lines by sample name
      p_smooth <- ggplot2::ggplot(
        smoothed,
        ggplot2::aes(
          x      = .data$position,
          y      = .data$mean_prob,
          colour = .data$smooth_key
        )
      ) +
        ggplot2::geom_line(linewidth = 1) +
        ggplot2::scale_y_continuous(
          limits = c(0, 1),
          name   = "Mean modification\nprobability"
        ) +
        ggplot2::scale_x_continuous(labels = scales::comma_format()) +
        ggplot2::coord_cartesian(xlim = c(region_start, region_end)) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
          panel.grid.minor = ggplot2::element_blank(),
          legend.text      = ggplot2::element_text(size = ggplot2::rel(0.75)),
          legend.title     = ggplot2::element_text(size = ggplot2::rel(0.75)),
          legend.key.size  = ggplot2::unit(0.4, "cm")
        ) +
        ggplot2::labs(x = "Genomic position (bp)", colour = "Sample")
    }

  } else {
    p_smooth <- ggplot2::ggplot() +
      ggplot2::annotate(
        "text", x = 0.5, y = 0.5,
        label = "No methylation sites", size = 4
      ) +
      ggplot2::theme_void()
  }

  # --- 3. Build gene annotation panel (optional) ---
  p_gene <- NULL
  if (!is.null(annotations)) {
    if (!inherits(annotations, "gene_annotations")) {
      stop("`annotations` must be a `gene_annotations` object from read_annotations().",
           call. = FALSE)
    }
    p_gene <- build_gene_panel(
      annotations  = annotations,
      region_start = region_start,
      region_end   = region_end,
      bottom_label = FALSE
    )
  }

  # --- 4. Assemble panels ---
  # Panel order (when gene panel present): gene (top), sample reads..., smooth
  # (bottom with genomic coordinates). When no gene panel: sample reads, smooth.
  if (!is.null(p_gene)) {
    all_panels <- c(list(p_gene), sample_panels, list(p_smooth))
    n_panels   <- n_samples + 2L
  } else {
    all_panels <- c(sample_panels, list(p_smooth))
    n_panels   <- n_samples + 1L
  }

  # Compute heights
  if (is.null(panel_heights)) {
    heights <- c(if (!is.null(p_gene)) 0.6 else NULL, rep(3, n_samples), 1)
  } else {
    if (length(panel_heights) != n_panels) {
      stop(sprintf(
        "`panel_heights` has length %d but there are %d panels.",
        length(panel_heights), n_panels
      ), call. = FALSE)
    }
    heights <- panel_heights
  }

  patchwork::wrap_plots(all_panels, ncol = 1, heights = heights)
}
