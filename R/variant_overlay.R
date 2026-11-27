# Internal functions for VCF variant overlay — not exported

#' Extract per-read bases at variant positions
#'
#' For each combination of (read, variant position) where the read's clipped
#' span (`start`..`end`) contains the variant, maps the position through the
#' read's CIGAR to its query base and classifies it as `"ref"`, `"alt"`,
#' `"del"` (no aligned query base there: a deletion, skip or uncovered
#' position), or `"other"`. Bases are compared case-insensitively; positions
#' past the end of the stored sequence are skipped. Each read's CIGAR is
#' decomposed once and all its variants are mapped in one vectorised step.
#' Rows are ordered by read (in `reads` order), then by variant (in `variants`
#' order).
#'
#' @param reads Data.frame. The `$reads` element of a `methylation_data` object.
#'   Must contain columns `read_name`, `start`, `end`, and `lane`.
#' @param sequences Named character vector. The `$sequences` element of a
#'   `methylation_data` object, keyed by `read_name`.
#' @param cigars Named character vector. The `$cigars` element of a
#'   `methylation_data` object, keyed by `read_name`.
#' @param variants Data.frame with columns `position`, `ref`, `alt`, `type`
#'   (as returned by [read_variants()]).
#'
#' @return A data.frame with columns `read_name` (character), `position`
#'   (integer), `base` (character), `variant_class` (character), and
#'   `lane` (integer).
#'
#' @keywords internal
extract_variant_bases = function(reads, sequences, cigars, variants) {
  empty = data.frame(
    read_name     = character(0L),
    position      = integer(0L),
    base          = character(0L),
    variant_class = character(0L),
    lane          = integer(0L),
    stringsAsFactors = FALSE
  )
  if (is.null(variants) || nrow(variants) == 0L) return(empty)

  var_pos = variants$position
  ref_uc  = toupper(variants$ref)
  alt_uc  = toupper(variants$alt)

  # Reads without a stored sequence or CIGAR are skipped (first match by name)
  seq_idx = match(reads$read_name, names(sequences))
  cig_idx = match(reads$read_name, names(cigars))
  usable  = which(!is.na(seq_idx) & !is.na(cig_idx) &
                    !is.na(unname(sequences)[seq_idx]) &
                    !is.na(unname(cigars)[cig_idx]))

  per_read = vector("list", nrow(reads))
  for (i in usable) {
    # Variants inside this read's clipped span, in `variants` order
    vi = which(var_pos >= reads$start[i] & var_pos <= reads$end[i])
    if (length(vi) == 0L) next

    # The CIGAR is walked from bam_pos (the unclipped alignment start), not
    # the region-clipped reads$start.
    q = .ref_to_query(cigars[[cig_idx[i]]], reads$bam_pos[i], var_pos[vi])

    in_del = is.na(q)  # no aligned query base: deletion, skip or uncovered
    base   = rep("-", length(vi))
    if (!all(in_del)) {
      base[!in_del] = substring(sequences[[seq_idx[i]]], q[!in_del], q[!in_del])
    }
    keep = in_del | nzchar(base)  # query position past the stored sequence
    if (!any(keep)) next
    vi     = vi[keep]
    base   = base[keep]
    in_del = in_del[keep]

    ub  = toupper(base)
    cls = ifelse(!is.na(ref_uc[vi]) & ub == ref_uc[vi], "ref",
                 ifelse(!is.na(alt_uc[vi]) & ub == alt_uc[vi], "alt", "other"))
    cls[in_del] = "del"
    per_read[[i]] = list(ri = rep(i, length(vi)), vi = vi, base = base, cls = cls)
  }

  per_read = per_read[!vapply(per_read, is.null, logical(1L))]
  if (length(per_read) == 0L) return(empty)

  ri = unlist(lapply(per_read, `[[`, "ri"), use.names = FALSE)
  vi = unlist(lapply(per_read, `[[`, "vi"), use.names = FALSE)
  data.frame(
    read_name     = reads$read_name[ri],
    position      = var_pos[vi],
    base          = unlist(lapply(per_read, `[[`, "base"), use.names = FALSE),
    variant_class = unlist(lapply(per_read, `[[`, "cls"), use.names = FALSE),
    lane          = reads$lane[ri],
    stringsAsFactors = FALSE
  )
}

# Query positions of the reference positions `targets` in one read, from its
# CIGAR and alignment start `pos`: vectorised ref_to_seq(). NA where a target
# has no aligned query base (inside D/N, or outside the M/=/X blocks).
.ref_to_query = function(cigar, pos, targets) {
  dc = decompose_cigar(cigar, pos)
  dc = dc[dc$type %in% c("M", "=", "X"), , drop = FALSE]
  q  = rep(NA_integer_, length(targets))
  if (nrow(dc) == 0L) return(q)
  k  = findInterval(targets, dc$ref_start)
  ok = k > 0L
  ok[ok] = targets[ok] <= dc$ref_end[k[ok]]
  q[ok] = dc$query_start[k[ok]] + as.integer(targets[ok] - dc$ref_start[k[ok]])
  q
}

# Warn when the region a variant_data object was read for is on another
# chromosome than, or does not overlap, the plotted methylation_data region:
# the overlay would then be silently empty. Mere inequality (e.g. variants read
# for a wider window) is fine.
.check_variant_region = function(data_region, variants_region) {
  if (!inherits(data_region, "GRanges") || !inherits(variants_region, "GRanges")) {
    return(invisible(NULL))
  }
  fmt = function(gr) {
    sprintf("%s:%d-%d", as.character(GenomicRanges::seqnames(gr)),
            GenomicRanges::start(gr), GenomicRanges::end(gr))
  }
  same_chrom = as.character(GenomicRanges::seqnames(data_region)) ==
    as.character(GenomicRanges::seqnames(variants_region))
  overlaps = same_chrom &&
    GenomicRanges::start(variants_region) <= GenomicRanges::end(data_region) &&
    GenomicRanges::end(variants_region) >= GenomicRanges::start(data_region)
  if (!overlaps) {
    warning(
      "`variants` were read for region ", fmt(variants_region),
      ", which does not overlap the plotted region ", fmt(data_region),
      "; no variants will be shown.",
      call. = FALSE
    )
  }
  invisible(NULL)
}

#' Compose all variant overlay layers for a single methylation_data object
#'
#' Dispatches to [build_snv_layer()], [build_sv_layer()], [build_bnd_layer()],
#' and [match_sa_to_vcf_bnd()] based on which variant types are present in
#' `variants$variants`. Returns a named list suitable for consumption by
#' [build_read_panel()].
#'
#' @param data A `methylation_data` object.
#' @param variants A `variant_data` object returned by [read_variants()], or
#'   `NULL`. Its `bnd_match_tol` sets the SA-to-BND matching tolerance
#'   (50 bp when absent). A warning is issued when its `region` is on another
#'   chromosome than, or does not overlap, `data$region`.
#'
#' @return A named list with elements `snv`, `sv`, `bnd`, and `sa_reads`
#'   (each is a list of ggplot2 layer objects or `NULL`), or `NULL` when
#'   `variants` is `NULL` or not a `variant_data` object.
#'
#' @keywords internal
build_variant_overlay = function(data, variants) {
  if (is.null(variants) || !inherits(variants, "variant_data")) {
    return(NULL)
  }

  .check_variant_region(data$region, variants$region)

  vdf = variants$variants

  # Partition by type (%in% so that NA types are dropped, not turned into
  # all-NA rows).
  # Only single-nucleotide variants are supported by the asterisk overlay.
  # Insertions, deletions and MNVs cannot be resolved to a single reference
  # base, so routing them through extract_variant_bases() would classify every
  # spanning read as "other" and mark all reads (issue #20). They, and
  # symbolic types other than DEL/DUP/INV (e.g. <CNV>), are not drawn.
  snv_rows = vdf[vdf$type %in% "SNV", , drop = FALSE]
  sv_rows  = vdf[vdf$type %in% c("DEL", "DUP", "INV"), , drop = FALSE]
  bnd_rows = vdf[vdf$type %in% "BND", , drop = FALSE]

  # --- SNV layer ---
  snv_layers = NULL
  if (nrow(snv_rows) > 0L) {
    if (!is.null(data$sequences) && length(data$sequences) > 0L) {
      vbases     = extract_variant_bases(data$reads, data$sequences, data$cigars, snv_rows)
      snv_layers = build_snv_layer(vbases)
    } else {
      warning(
        "VCF overlay requires sequences stored in methylation_data. ",
        "Re-run read_methylation() to get a current object with sequences.",
        call. = FALSE
      )
    }
  }

  # --- SV layer ---
  region_start = GenomicRanges::start(data$region)
  region_end   = GenomicRanges::end(data$region)
  sv_layers    = build_sv_layer(data$reads, sv_rows, region_start, region_end)

  # --- BND layer and SA matching ---
  bnd_layers         = build_bnd_layer(bnd_rows)
  sa_reads_validated = NULL
  if ("sa_chrom" %in% names(data$reads)) {
    tol = if (is.null(variants$bnd_match_tol)) 50L else variants$bnd_match_tol
    sa_reads_validated = match_sa_to_vcf_bnd(data$reads, bnd_rows, tol = tol)
  }

  list(
    snv      = snv_layers,
    sv       = sv_layers,
    bnd      = bnd_layers,
    sa_reads = sa_reads_validated
  )
}
