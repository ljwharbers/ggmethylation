# Internal functions for VCF variant overlay — not exported

#' Extract per-read bases at variant positions
#'
#' For each combination of (read, variant position) where the read spans the
#' variant, walks the CIGAR string to find the query base and classifies it as
#' `"ref"`, `"alt"`, `"del"` (position falls in a deletion), or `"other"`.
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
extract_variant_bases <- function(reads, sequences, cigars, variants) {
  # Handle empty variants immediately
  if (is.null(variants) || nrow(variants) == 0L) {
    return(data.frame(
      read_name     = character(0L),
      position      = integer(0L),
      base          = character(0L),
      variant_class = character(0L),
      lane          = integer(0L),
      stringsAsFactors = FALSE
    ))
  }

  result_list <- vector("list", nrow(reads) * nrow(variants))
  result_idx  <- 0L

  for (i in seq_len(nrow(reads))) {
    read_name <- reads$read_name[i]

    # Skip reads without a stored sequence
    if (!read_name %in% names(sequences)) next
    if (!read_name %in% names(cigars))    next

    seq_str <- sequences[[read_name]]
    cigar   <- cigars[[read_name]]
    r_start <- reads$start[i]
    r_end   <- reads$end[i]
    lane    <- reads$lane[i]

    for (v in seq_len(nrow(variants))) {
      var_pos <- variants$position[v]

      # Skip if variant is outside this read's span
      if (var_pos < r_start || var_pos > r_end) next

      # Map genomic position to query position via CIGAR
      q_pos <- ref_to_seq(
        cigar          = cigar,
        ref_start      = r_start,
        target_ref_pos = var_pos
      )

      if (is.na(q_pos)) {
        # Position falls in a deletion
        base          <- "-"
        variant_class <- "del"
      } else {
        base <- substr(seq_str, q_pos, q_pos)
        if (nchar(base) == 0L || base == "") {
          next  # q_pos out of sequence length — skip
        }
        variant_class <- if (base == variants$ref[v]) {
          "ref"
        } else if (base == variants$alt[v]) {
          "alt"
        } else {
          "other"
        }
      }

      result_idx <- result_idx + 1L
      result_list[[result_idx]] <- data.frame(
        read_name     = read_name,
        position      = var_pos,
        base          = base,
        variant_class = variant_class,
        lane          = lane,
        stringsAsFactors = FALSE
      )
    }
  }

  if (result_idx == 0L) {
    return(data.frame(
      read_name     = character(0L),
      position      = integer(0L),
      base          = character(0L),
      variant_class = character(0L),
      lane          = integer(0L),
      stringsAsFactors = FALSE
    ))
  }

  do.call(rbind, result_list[seq_len(result_idx)])
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
#'   `NULL`.
#' @param bnd_match_tol Integer. Position tolerance (bp) for matching SA
#'   breakpoints to VCF BND calls. Default `50L`.
#'
#' @return A named list with elements `snv`, `sv`, `bnd`, and `sa_reads`
#'   (each is a list of ggplot2 layer objects or `NULL`), or `NULL` when
#'   `variants` is `NULL` or not a `variant_data` object.
#'
#' @keywords internal
build_variant_overlay <- function(data, variants, bnd_match_tol = 50L) {
  if (is.null(variants) || !inherits(variants, "variant_data")) {
    return(NULL)
  }

  vdf <- variants$variants

  # Partition by type
  snv_rows <- vdf[vdf$type %in% c("SNV", "insertion", "deletion"), , drop = FALSE]
  sv_rows  <- vdf[vdf$type %in% c("DEL", "DUP", "INV"),           , drop = FALSE]
  bnd_rows <- vdf[vdf$type == "BND",                               , drop = FALSE]

  # --- SNV layer ---
  snv_layers <- NULL
  if (nrow(snv_rows) > 0L) {
    if (!is.null(data$sequences) && length(data$sequences) > 0L) {
      vbases     <- extract_variant_bases(data$reads, data$sequences, data$cigars, snv_rows)
      snv_layers <- build_snv_layer(vbases)
    } else {
      warning(
        "VCF overlay requires sequences stored in methylation_data. ",
        "Re-run read_methylation() to get a current object with sequences.",
        call. = FALSE
      )
    }
  }

  # --- SV layer ---
  region_start <- GenomicRanges::start(data$region)
  region_end   <- GenomicRanges::end(data$region)
  sv_layers    <- build_sv_layer(data$reads, sv_rows, region_start, region_end)

  # --- BND layer and SA matching ---
  bnd_layers         <- build_bnd_layer(bnd_rows)
  sa_reads_validated <- NULL
  if ("sa_chrom" %in% names(data$reads)) {
    sa_reads_validated <- match_sa_to_vcf_bnd(data$reads, bnd_rows, tol = bnd_match_tol)
  }

  list(
    snv      = snv_layers,
    sv       = sv_layers,
    bnd      = bnd_layers,
    sa_reads = sa_reads_validated
  )
}
