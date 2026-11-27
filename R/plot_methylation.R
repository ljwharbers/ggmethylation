# plot_methylation(): argument validation, read sorting and lane packing, and assembly of the patchwork composite.

# Stop early on a `sort_by` column that doesn't exist: order() on a NULL column
# returns integer(0) and would silently drop every read.
.check_sort_by = function(sort_by, reads) {
  missing_cols = setdiff(sort_by, names(reads))
  if (length(missing_cols) > 0L) {
    stop("`sort_by` column(s) not found in `data$reads`: ",
         paste(missing_cols, collapse = ", "), call. = FALSE)
  }
}

# Centred text placeholder used for empty panels.
.message_plot = function(label, size) {
  ggplot2::ggplot() +
    ggplot2::annotate("text", x = 0.5, y = 0.5, label = label, size = size) +
    ggplot2::theme_void()
}

# Site frame used for every aggregate (smooth panel, delta panel, per-read
# means). In binary mode these aggregate 0/1 calls instead of raw
# probabilities, so they report a fraction of methylated calls. The raw
# `data$sites` is still what build_read_panel() colours, so it can keep the
# ambiguous category.
.sites_for_aggregation = function(sites, call_threshold, call_ambiguous) {
  if (is.null(call_threshold)) sites else .binarize_sites(sites, call_threshold, call_ambiguous)
}

# Sort reads and pack them into lanes (per group when grouped, with a blank
# separator lane between groups). Adds `mean_mod_prob` and `lane` to
# `data$reads`; returns list(data, separator_lanes).
.prepare_reads = function(data, sites_agg, sort_by = NULL) {
  reads = data$reads
  reads$mean_mod_prob = .read_mean_mod_prob(sites_agg, reads$read_name)
  reads$mean_mod_prob[is.na(reads$mean_mod_prob)] = 0

  if (is.null(sort_by)) {
    sort_by = if (is.null(data$group_tag)) "start" else c("start", "group", "mean_mod_prob")
  }
  .check_sort_by(sort_by, reads)
  reads = reads[do.call(order, lapply(sort_by, function(col) reads[[col]])), , drop = FALSE]
  rownames(reads) = NULL

  separator_lanes = numeric(0)
  if (is.null(data$group_tag)) {
    reads$lane = pack_reads(reads, clip_side = reads$clip_side)
  } else {
    reads$lane = integer(nrow(reads))
    lane_offset = 0L
    for (grp in .ordered_plot_groups(reads$group)) {
      idx = which(.match_plot_group(reads$group, grp))
      reads$lane[idx] = pack_reads(reads[idx, ], clip_side = reads$clip_side[idx]) + lane_offset
      lane_offset = max(reads$lane[idx]) + 2L
      separator_lanes = c(separator_lanes, lane_offset - 1L)
    }
    # No separator after the last group
    separator_lanes = utils::head(separator_lanes, -1L)
  }

  data$reads = reads
  list(data = data, separator_lanes = separator_lanes)
}

#' Plot read-level methylation data
#'
#' Creates a ggplot2 visualisation of read-level base modification data,
#' returned as a [patchwork::patchwork] composite: an optional gene track, one
#' read panel per sample (reads as bars with a coloured mark per modification
#' site), a panel with the loess-smoothed mean modification probability (one
#' line per group when the data is grouped, otherwise one overall line), and an
#' optional delta panel (`show_delta`). Setting `call_threshold` switches to
#' binary calls: sites are coloured methylated/unmethylated and the smooth
#' panel plots the fraction of methylated calls instead of the mean probability.
#'
#' Groups are coloured consistently in every panel (read bars, smooth lines,
#' delta area) by one palette resolved from `group_colours`; reads whose group
#' is `NA` (e.g. unphased) are drawn in grey. Legends shared by several panels
#' are collected once at the right of the figure.
#'
#' When the data holds several modification codes (e.g.
#' `mod_code = c("m", "h")`), each read bar is split into horizontal strips,
#' one per code, in the order of `data$mod_code` from top to bottom (the read
#' panel's y-axis title lists the order), and the smooth panel draws one line
#' per code (by line type).
#'
#' @param data A `methylation_data` object returned by [read_methylation()],
#'   or a `multi_methylation_data` object from [merge_methylation()].
#' @param sort_by Character vector of column names from `data$reads` used
#'   to sort reads before packing into lanes (reads are always packed per
#'   group, so group order is fixed). Default NULL uses `"start"` when
#'   ungrouped or `c("start", "group", "mean_mod_prob")` when grouped.
#'   `mean_mod_prob` is the per-read mean modification probability, or the
#'   per-read fraction of methylated calls when `call_threshold` is set.
#' @param colour_low Colour for low modification probability (default
#'   `"#BDBDBD"`).
#' @param colour_high Colour for high modification probability (default
#'   `"#C62828"`).
#' @param colour_ambiguous Colour for ambiguous calls (default `"#78909C"`, a
#'   slate blue-grey deliberately off the `colour_low`/`colour_high` ramp so it
#'   reads as "no confident call" rather than "intermediate methylation"). Only
#'   used when `call_threshold` and `call_ambiguous` are set.
#' @param line_width Positive number: linewidth of modification site lines
#'   (default 0.2).
#' @param colour_strand Logical. When `TRUE`, read bars are coloured by strand
#'   (`"+"` = forward, `"-"` = reverse). Ignored when data is grouped; group
#'   colour takes precedence. Default `FALSE`.
#' @param strand_colours Named character vector with `"+"` and `"-"` entries
#'   giving bar colours for each strand. Only used when `colour_strand = TRUE`
#'   and data is ungrouped. Default `c("+" = "#4393C3", "-" = "#D6604D")`.
#' @param group_colours Character vector of colours per group, used for the
#'   read bars, the smooth lines and the delta area alike. Defaults to
#'   `c("1" = "#0072B2", "2" = "#E69F00")` (the Okabe-Ito colorblind-safe
#'   blue/orange pair), matching HP haplotype tags. Groups the vector does not
#'   name (e.g. `"REF"`/`"ALT"` from SNV grouping, or a third haplotype) get
#'   further Okabe-Ito colours (with a message when you supplied the vector).
#'   An unnamed vector is matched to the sorted group names by position and
#'   needs at least one colour per group. Pass `NULL` to use ggplot2's default
#'   hue palette.
#' @param smooth_span `NULL` or a positive number: the loess smoothing span for
#'   the smooth panel. When `NULL` (default), an adaptive span is computed per
#'   group as `max(0.15, min(0.75, 15 / n_unique_sites))`, targeting ~15 data
#'   points per local fit regardless of region size or CpG density. Pass an
#'   explicit value (e.g. `0.3`) to use a fixed span.
#' @param panel_heights Numeric vector giving relative heights of the panels,
#'   top to bottom: the gene track (when `annotations` is given), one read
#'   panel per sample, the smooth panel, and the delta panel (when shown).
#'   Its length must equal the number of panels. Pass `NULL` (default) to use
#'   the built-in defaults.
#' @param annotations A `gene_annotations` object returned by
#'   [read_annotations()], or `NULL` (default). When provided, a gene
#'   annotation track is placed at the top of the composite figure (above the
#'   read panel and the smooth panel).
#' @param variants A `variant_data` object returned by [read_variants()], or
#'   `NULL` (default). When provided, VCF records are overlaid on the read
#'   panel(s):
#'   * **SNV**: a small black asterisk at the variant position on each read
#'     whose base there equals the ALT allele. Reads with the REF or another
#'     base are not marked; no base letters are drawn.
#'   * **DEL, DUP, INV** (symbolic SVs with an `END`): a styled segment on
#'     every read overlapping the SV span, clipped to the read and the region:
#'     DEL a dark red line with end ticks, DUP two teal parallel lines, INV a
#'     purple dashed line. Chevrons (« ») mark a span cut at the region
#'     edge.
#'   * **BND**: one text label `"BND mate_chrom:mate_pos"` at the bottom of
#'     the read panel, centred on the breakend position. When
#'     `show_supplementary = TRUE`, the SA marker of each read whose
#'     supplementary alignment lies on the BND mate chromosome within the
#'     object's `bnd_match_tol` of the mate position gets a dark magenta
#'     border.
#'   * Short insertions, deletions and MNVs, and symbolic types other than
#'     DEL/DUP/INV (e.g. `<CNV>`), are not drawn.
#'
#'   No vertical lines are drawn for variants. Requires that `data` was
#'   produced by a current version of [read_methylation()], which stores
#'   sequences and CIGARs.
#' @param show_cigar Logical. When `TRUE` (default), structural variants from
#'   CIGAR strings (insertions, deletions) are displayed on reads. Insertions
#'   appear as purple I-beam markers spanning the full height of the read bar;
#'   deletions open a gap in the read bar and draw a thin black line. Smooth
#'   lines (and the delta area) are interrupted where at least 75% of a
#'   group's reads spanning the whole deletion carry it. Set to `FALSE` to hide these
#'   features.
#' @param min_indel_size Non-negative number. Minimum size (in bp) for
#'   insertions and deletions to be displayed when `show_cigar = TRUE`. Indels
#'   smaller than this threshold are hidden to reduce visual clutter from
#'   common small indels. Default `50`.
#' @param show_supplementary Logical. When `TRUE` (default), reads with a
#'   supplementary alignment (SA BAM tag) get a small marker at each clipped
#'   end, filled by the chromosome of the SA partner. The marker is drawn
#'   outside the read bar, so no modification site is hidden.
#' @param show_ci Logical. When `TRUE` (default), a shaded ribbon showing the
#'   loess confidence interval (`lower`/`upper` from [smooth_methylation()])
#'   is drawn behind each smooth line in the smooth panel. Has no effect when
#'   fewer than 4 unique positions are available for a group/code (no CI is
#'   computed in that case). Set to `FALSE` to hide the ribbon.
#' @param call_threshold `NULL` (default) or a number in `[0, 1]`. With
#'   `NULL`, modification sites are coloured by a continuous probability
#'   gradient (`colour_low` to `colour_high`) and the smooth panel plots the
#'   mean modification probability per position. A number switches to
#'   **binary calls**: sites with probability at or above the threshold are
#'   methylated, the rest unmethylated (or ambiguous, see `call_ambiguous`),
#'   coloured with a discrete scale. The same split drives the read panel and
#'   every aggregate, so the two always agree.
#'
#'   Binary calls change what every aggregate reports: the smooth panel, the
#'   delta panel (`show_delta`) and the per-read `mean_mod_prob` used for
#'   sorting all aggregate 0/1 calls rather than raw probabilities, so the
#'   smooth panel plots the **fraction of calls that are methylated** and its
#'   y-axis is relabelled accordingly. The two quantities genuinely differ: a
#'   position where every read reports probability 0.8 has a mean probability of
#'   0.8 but a methylated fraction of 1.0. The fraction is the conventional
#'   percent-methylation / beta value; the continuous mean is systematically
#'   pulled toward 0.5 by basecaller uncertainty.
#'
#'   Because binarised values are 0 or 1, positions with thin coverage
#'   contribute coarser values than with continuous probabilities, so the
#'   binary curve is noisier where few reads overlap. The loess span averages
#'   over neighbouring positions. Note that the loess fit weights every
#'   position equally regardless of how many reads cover it, so the confidence
#'   ribbon (`show_ci`) reflects the scatter and density of positions, not read
#'   coverage: it does not widen where coverage is thin. Fitted values are
#'   clamped to `[0, 1]`.
#' @param call_ambiguous `NULL` (default) for a hard threshold, or a single
#'   non-negative number `w`: sites whose probability lies in the open band
#'   `(call_threshold - w, call_threshold + w)` are labelled "ambiguous"
#'   rather than methylated/unmethylated, except a probability exactly equal
#'   to `call_threshold`, which is always methylated. Ambiguous sites are
#'   excluded from both the numerator and the denominator of the methylated
#'   fraction, so the smooth panel reports the fraction methylated *among
#'   confident calls*. The default `NULL` (or `0`) discards nothing — every
#'   site is called at `call_threshold`. Requires `call_threshold`.
#' @param show_delta Logical. When `TRUE`, and grouping yields exactly two
#'   groups, an additional panel is appended below the smooth panel showing
#'   the signed difference in loess-smoothed modification probability between
#'   the two groups (group2 - group1), filled with the colour of whichever
#'   group is higher. Requires `data` to be grouped (`group_tag` or `snv` in
#'   [read_methylation()]) with exactly two distinct group values and a single
#'   modification code; otherwise the panel is skipped with a message. Not
#'   supported for multi-sample data (`multi_methylation_data`); the argument
#'   is accepted there but ignored (with a message). Default `FALSE`.
#'
#' @return A [patchwork::patchwork] composite (also for ungrouped data), or,
#'   when a single-sample `data` has no reads, a placeholder [ggplot2::ggplot]
#'   saying so.
#'
#' @examples
#' bam = system.file("extdata", "synthetic.bam", package = "ggmethylation")
#'
#' # Reads grouped by haplotype: one smooth line per haplotype. Binary calls
#' # at probability 0.5, plus the difference between the haplotypes
#' md = read_methylation(bam, "chr1:1000-7000", group_tag = "HP",
#'                       max_reads = 30, seed = 1)
#' plot_methylation(md, call_threshold = 0.5, show_delta = TRUE)
#'
#' # 5mC and 5hmC: one strip per code on each read, one line type per code
#' md2 = read_methylation(bam, "chr1:1000-7000", mod_code = c("m", "h"),
#'                        group_tag = "HP", max_reads = 30, seed = 1)
#' plot_methylation(md2, group_colours = c("1" = "steelblue", "2" = "coral"))
#'
#' @examplesIf requireNamespace("VariantAnnotation", quietly = TRUE) && requireNamespace("SummarizedExperiment", quietly = TRUE)
#' # Overlay the calls of a VCF (an SNV, a deletion and a breakend)
#' vcf = system.file("extdata", "synthetic.vcf.gz", package = "ggmethylation")
#' vd = read_variants(vcf, "chr1:1000-7000")
#' plot_methylation(md, variants = vd)
#'
#' @export
plot_methylation = function(data, sort_by = NULL,
                            colour_low = "#BDBDBD",
                            colour_high = "#C62828",
                            colour_ambiguous = "#78909C",
                            line_width = 0.2,
                            colour_strand = FALSE,
                            strand_colours = c("+" = "#4393C3", "-" = "#D6604D"),
                            group_colours = c("1" = "#0072B2", "2" = "#E69F00"),
                            smooth_span = NULL,
                            panel_heights = NULL,
                            annotations = NULL,
                            variants = NULL,
                            show_cigar = TRUE,
                            min_indel_size = 50L,
                            show_supplementary = TRUE,
                            show_ci = TRUE,
                            call_threshold = NULL,
                            call_ambiguous = NULL,
                            show_delta = FALSE) {
  # --- 1. Validate input ---
  multi = inherits(data, "multi_methylation_data")
  if (!multi && !inherits(data, "methylation_data")) {
    stop(
      "`data` must be a `methylation_data` or `multi_methylation_data` object.",
      call. = FALSE
    )
  }
  .validate_plot_args(
    sort_by = sort_by, colour_low = colour_low, colour_high = colour_high,
    colour_ambiguous = colour_ambiguous,
    colour_strand = colour_strand, strand_colours = strand_colours,
    call_threshold = call_threshold, call_ambiguous = call_ambiguous,
    smooth_span = smooth_span, group_colours = group_colours,
    line_width = line_width, min_indel_size = min_indel_size,
    panel_heights = panel_heights,
    show_cigar = show_cigar, show_supplementary = show_supplementary,
    show_ci = show_ci, show_delta = show_delta
  )
  if (!is.null(variants) && !inherits(variants, "variant_data")) {
    stop("`variants` must be a `variant_data` object returned by read_variants().",
         call. = FALSE)
  }
  if (!is.null(annotations) && !inherits(annotations, "gene_annotations")) {
    stop("`annotations` must be a `gene_annotations` object from read_annotations().",
         call. = FALSE)
  }

  samples = if (multi) data$samples else list(data)
  grouped_samples = !vapply(samples, function(s) is.null(s$group_tag), logical(1L))
  if (isTRUE(colour_strand) && any(grouped_samples)) {
    message("`colour_strand = TRUE` is ignored when reads are grouped; group colour takes precedence.")
  }
  if (multi && isTRUE(show_delta)) {
    message("`show_delta` is not supported for multi-sample data; ignoring.")
  }
  if (!multi && isTRUE(show_delta) && is.null(data$group_tag)) {
    message("`show_delta` needs grouped data (`group_tag` or `snv` in read_methylation()); skipping the delta panel.")
  }
  if (!multi && nrow(data$reads) == 0L) {
    return(.message_plot("No reads in region", size = 5))
  }

  region_start = GenomicRanges::start(data$region)
  region_end   = GenomicRanges::end(data$region)

  # The gene panel only sees coordinates, so the chromosome is checked here.
  ann_region = annotations$region
  if (inherits(ann_region, "GRanges") && length(ann_region) == 1L) {
    ann_chrom  = as.character(GenomicRanges::seqnames(ann_region))
    data_chrom = as.character(GenomicRanges::seqnames(data$region))
    if (!identical(ann_chrom, data_chrom)) {
      warning(sprintf(
        "`annotations` were read for chromosome '%s' but the data are on '%s'.",
        ann_chrom, data_chrom
      ), call. = FALSE)
    }
  }

  # One group palette for every panel and sample, resolved over the union of
  # the groups present in any sample's reads.
  all_groups = unlist(lapply(samples[grouped_samples], function(s) s$reads$group))
  resolved_colours = .resolve_group_colours(group_colours, all_groups)
  # Legend keys shared by every read panel, so patchwork can merge them.
  legend_levels = list(
    group = names(resolved_colours),
    group_na = any(is.na(all_groups)),
    sa = .group_levels(unlist(lapply(samples, function(s) s$reads$sa_chrom)))
  )
  opts = list(smooth_span = smooth_span, min_indel_size = min_indel_size,
              show_cigar = show_cigar, show_ci = show_ci, binary = !is.null(call_threshold),
              group_colours = resolved_colours)

  # --- 2. Sort and pack each sample's reads ---
  prepared = lapply(samples, function(s) {
    if (nrow(s$reads) == 0L) return(NULL)
    sites_agg = .sites_for_aggregation(s$sites, call_threshold, call_ambiguous)
    c(.prepare_reads(s, sites_agg, sort_by), list(sites_agg = sites_agg))
  })

  # --- 3. Read panel(s) ---
  sample_names = if (multi) names(samples) else ""
  read_panels = Map(function(prep, nm) {
    if (is.null(prep)) {
      return(.message_plot(sprintf("No reads (%s)", nm), size = 4) + ggplot2::labs(title = nm))
    }
    s = prep$data
    p = build_read_panel(
      data               = s,
      separator_lanes    = prep$separator_lanes,
      region_start       = region_start,
      region_end         = region_end,
      colour_low         = colour_low,
      colour_high        = colour_high,
      colour_ambiguous   = colour_ambiguous,
      line_width         = line_width,
      colour_strand      = colour_strand,
      strand_colours     = strand_colours,
      group_colours      = resolved_colours,
      show_x_axis        = FALSE,
      variant_overlay    = build_variant_overlay(s, variants),
      show_cigar         = show_cigar,
      min_indel_size     = min_indel_size,
      show_supplementary = show_supplementary,
      call_threshold     = call_threshold,
      call_ambiguous     = call_ambiguous,
      legend_levels      = legend_levels
    )
    if (multi) p + ggplot2::labs(title = nm) else p
  }, prepared, sample_names)
  names(read_panels) = NULL

  # --- 4. Smooth panel, plus the delta panel for two groups ---
  p_delta = NULL
  if (multi) {
    p_smooth = .multi_smooth_panel(samples, prepared, region_start, region_end, opts)
  } else {
    prep = prepared[[1L]]
    p_smooth = .single_smooth_panel(prep$data, prep$sites_agg, region_start, region_end, opts)
    if (!is.null(data$snv_position)) {
      p_smooth = p_smooth +
        ggplot2::geom_vline(
          xintercept = data$snv_position,
          linetype = "dashed", colour = "black", linewidth = 0.5
        )
    }
    if (isTRUE(show_delta) && !is.null(data$group_tag)) {
      p_delta = .delta_panel(prep$data, prep$sites_agg, region_start, region_end, opts)
      if (!is.null(p_delta)) {
        # The smooth panel is no longer the bottom-most; hide its x-axis title
        p_smooth = p_smooth + ggplot2::theme(axis.title.x = ggplot2::element_blank())
      }
    }
  }

  # --- 5. Gene track on top, then combine ---
  p_gene = if (!is.null(annotations)) {
    build_gene_panel(annotations, region_start, region_end, bottom_label = FALSE)
  }
  panels = c(if (!is.null(p_gene)) list(p_gene), read_panels, list(p_smooth),
             if (!is.null(p_delta)) list(p_delta))

  heights = panel_heights
  if (is.null(heights)) {
    heights = if (multi) {
      c(if (!is.null(p_gene)) 0.6, rep(3, length(samples)), 1)
    } else {
      # With a delta panel the two summary panels get more room, or their
      # y-axis titles run into each other at ordinary figure heights.
      if (is.null(p_delta)) {
        c(if (!is.null(p_gene)) 0.08, 1, 0.25)
      } else {
        c(if (!is.null(p_gene)) 0.08, 1, 0.35, 0.3)
      }
    }
  } else if (length(heights) != length(panels)) {
    stop(sprintf(
      "`panel_heights` has length %d but there are %d panels.",
      length(heights), length(panels)
    ), call. = FALSE)
  }

  # Identical legends (modification scale, group, SA partner) from several
  # panels are merged into one.
  patchwork::wrap_plots(panels, ncol = 1, heights = heights, guides = "collect")
}

# Argument checks for plot_methylation(), each with a message naming the
# argument (instead of an R-level "missing value where TRUE/FALSE needed" or
# a silent fallback further down).
.validate_plot_args = function(sort_by, colour_low, colour_high, colour_ambiguous,
                               colour_strand, strand_colours, call_threshold,
                               call_ambiguous, smooth_span, group_colours,
                               line_width, min_indel_size, panel_heights,
                               show_cigar, show_supplementary, show_ci,
                               show_delta) {
  .check_flag(colour_strand, "colour_strand")
  .check_flag(show_cigar, "show_cigar")
  .check_flag(show_supplementary, "show_supplementary")
  .check_flag(show_ci, "show_ci")
  .check_flag(show_delta, "show_delta")

  # Checked here, not by .check_sort_by(), so the type is caught before the
  # "No reads" early return. An empty vector would make order() drop every read.
  if (!is.null(sort_by) &&
      !(is.character(sort_by) && length(sort_by) > 0L && !anyNA(sort_by))) {
    stop("`sort_by` must be NULL or a character vector of column names without NAs.",
         call. = FALSE)
  }

  .check_colours(colour_low, "colour_low", n = 1L)
  .check_colours(colour_high, "colour_high", n = 1L)
  .check_colours(colour_ambiguous, "colour_ambiguous", n = 1L)
  if (isTRUE(colour_strand)) {
    if (!all(c("+", "-") %in% names(strand_colours))) {
      stop("`strand_colours` must be a named vector with \"+\" and \"-\" entries.",
           call. = FALSE)
    }
    .check_colours(strand_colours, "strand_colours")
  }

  if (!is.null(call_threshold)) {
    .check_number(call_threshold, "call_threshold", min = 0, max = 1)
  }
  if (!is.null(call_ambiguous)) {
    if (is.null(call_threshold)) {
      stop("`call_ambiguous` only applies to binary calls; set `call_threshold` too.",
           call. = FALSE)
    }
    .check_number(call_ambiguous, "call_ambiguous", min = 0)
  }
  # .check_number()'s bounds are inclusive; these two must be strictly positive.
  if (!is.null(smooth_span)) {
    .check_number(smooth_span, "smooth_span", min = 0)
    if (smooth_span == 0) {
      stop("`smooth_span` must be NULL or a single positive number (e.g. 0.3).",
           call. = FALSE)
    }
  }
  .check_number(line_width, "line_width", min = 0)
  if (line_width == 0) {
    stop("`line_width` must be a single positive number.", call. = FALSE)
  }
  .check_number(min_indel_size, "min_indel_size", min = 0)

  # Only the values are checked here; their number is checked against the
  # panels once those are known.
  if (!is.null(panel_heights) &&
      !(is.numeric(panel_heights) && length(panel_heights) > 0L &&
          all(is.finite(panel_heights)) && all(panel_heights > 0))) {
    stop("`panel_heights` must be NULL or a numeric vector of finite positive numbers.",
         call. = FALSE)
  }

  if (!is.null(group_colours)) {
    .check_colours(group_colours, "group_colours")
    nms = names(group_colours)
    if (!is.null(nms) && (any(is.na(nms) | !nzchar(nms)) || anyDuplicated(nms))) {
      stop("`group_colours` must be either unnamed or named with unique, non-empty group names.",
           call. = FALSE)
    }
  }
  invisible(TRUE)
}
