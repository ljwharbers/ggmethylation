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
