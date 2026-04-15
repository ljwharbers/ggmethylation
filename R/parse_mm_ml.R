# MM/ML tag parser for ggmethylation

#' Convert query (read) positions to reference (genomic) positions via CIGAR
#'
#' Walks through a CIGAR string to build a mapping from 1-based query positions
#' to 1-based reference positions. Positions that fall within soft clips or
#' insertions return `NA`.
#'
#' @param cigar Character. A CIGAR string (e.g., `"5S10M2I3M4D5M"`).
#' @param pos Integer. The 1-based leftmost reference position of the
#'   alignment (the POS field from BAM, which corresponds to the first
#'   M/=/X base after any leading S/H).
#' @param query_positions Integer vector. 1-based positions within the query
#'   sequence to map.
#'
#' @return Integer vector of the same length as `query_positions`, with
#'   reference positions or `NA` for unmappable positions.
#'
#' @keywords internal
seq_to_ref <- function(cigar, pos, query_positions) {
  # Parse CIGAR into operations and lengths
  parsed_cigar <- split_cigar(cigar)
  ops  <- parsed_cigar$ops
  lens <- parsed_cigar$lens

  # Pre-allocate a mapping vector from query offset to ref position.
  # Compute total query consumption to size the vector.
  query_consumers <- ops %in% c("M", "I", "S", "=", "X")
  total_query <- sum(lens[query_consumers])
  query_to_ref <- rep(NA_integer_, total_query)

  query_offset <- 0L
  ref_offset <- 0L

  for (i in seq_along(ops)) {
    op <- ops[i]
    len <- lens[i]

    if (op %in% c("M", "=", "X")) {
      # Consumes both query and reference
      idx <- seq_len(len)
      query_to_ref[query_offset + idx] <- pos + ref_offset + idx - 1L
      query_offset <- query_offset + len
      ref_offset <- ref_offset + len
    } else if (op == "I") {
      # Consumes query only -- positions map to NA (already NA)
      query_offset <- query_offset + len
    } else if (op == "S") {
      # Soft clip: consumes query only -- positions map to NA (already NA)
      query_offset <- query_offset + len
    } else if (op %in% c("D", "N")) {
      # Consumes reference only
      ref_offset <- ref_offset + len
    } else if (op == "H") {
      # Hard clip: consumes neither
    }
  }

  query_to_ref[query_positions]
}

#' Parse MM and ML tags for a single read
#'
#' Extracts modification probabilities and their genomic positions from the
#' MM and ML tags of a single BAM read. Supports the SAM specification for
#' base modification tags.
#'
#' @param seq Character. The read sequence from the BAM SEQ field.
#' @param mm_tag Character. The MM tag value, e.g. `"C+m,0,1,3;C+h,2;"`.
#' @param ml_tag Integer vector. ML tag values (0--255), concatenated for all
#'   modifications listed in the MM tag.
#' @param mod_code Character. The modification code to extract (e.g., `"m"`
#'   for 5mC, `"h"` for 5hmC).
#' @param strand Character. `"+"` or `"-"`, the strand of the alignment.
#' @param cigar Character. The CIGAR string for the alignment.
#' @param pos Integer. 1-based leftmost alignment position (BAM POS field).
#'
#' @return A `data.frame` with columns:
#'   \describe{
#'     \item{position}{Integer. 1-based genomic position of the modification.}
#'     \item{mod_prob}{Numeric. Probability of modification (0--1), computed
#'       as ML value / 255.}
#'   }
#'   Returns an empty data.frame (zero rows, same columns) when the
#'   requested modification is not present.
#'
#' @keywords internal
parse_mm_ml <- function(seq, mm_tag, ml_tag, mod_code, strand, cigar, pos) {
  empty_result <- data.frame(
    position = integer(0),
    mod_prob = numeric(0),
    stringsAsFactors = FALSE
  )


  # --- Guard clauses ---
  if (is.null(mm_tag) || is.na(mm_tag) || !nzchar(mm_tag)) {
    return(empty_result)
  }
  if (is.null(ml_tag) || length(ml_tag) == 0) {
    return(empty_result)
  }


  # --- 1. Parse MM tag string ---
  # Remove trailing semicolon if present and split entries
  mm_clean <- sub(";$", "", mm_tag)
  entries <- strsplit(mm_clean, ";")[[1]]
  entries <- trimws(entries)

  # Parse each entry to find the target mod_code
  target_idx <- NA_integer_
  parsed_entries <- vector("list", length(entries))

  for (j in seq_along(entries)) {
    entry <- entries[j]
    # Strip '?' and '.' flags (implicit unmodified bases indicators)
    entry <- gsub("[?.]", "", entry)

    # Split on comma: first element is like "C+m", rest are deltas
    parts <- strsplit(entry, ",")[[1]]
    spec <- parts[1]
    deltas <- if (length(parts) > 1) as.integer(parts[2:length(parts)]) else integer(0)

    # Parse the spec: canonical_base + strand_char + code
    # e.g., "C+m" -> base="C", mm_strand="+", code="m"
    canonical_base <- substr(spec, 1, 1)
    mm_strand <- substr(spec, 2, 2)
    code <- substring(spec, 3)

    parsed_entries[[j]] <- list(
      canonical_base = canonical_base,
      mm_strand = mm_strand,
      code = code,
      deltas = deltas,
      n_values = length(deltas)
    )

    if (code == mod_code) {
      target_idx <- j
    }
  }

  if (is.na(target_idx)) {
    return(empty_result)
  }

  target <- parsed_entries[[target_idx]]
  deltas <- target$deltas

  if (length(deltas) == 0) {
    return(empty_result)
  }


  # --- 2. Extract ML probabilities for this modification ---
  # Count ML values consumed by entries before the target
  ml_offset <- 0L
  if (target_idx > 1) {
    for (j in seq_len(target_idx - 1)) {
      ml_offset <- ml_offset + parsed_entries[[j]]$n_values
    }
  }
  ml_values <- ml_tag[ml_offset + seq_len(target$n_values)]


  # --- 3. Find canonical base positions in the read sequence ---
  canonical_base <- target$canonical_base
  mm_strand <- target$mm_strand

  # Determine search base and scan direction based on both the MM strand

  # field (+/-) and the alignment strand.
  # MM strand "+" = modification on same strand as read
  # MM strand "-" = modification on opposite strand
  if (mm_strand == "+") {
    if (strand == "+") {
      search_base <- canonical_base
      reverse_scan <- FALSE
    } else {
      search_base <- complement_base(canonical_base)
      reverse_scan <- TRUE
    }
  } else {
    # MM strand "-": modification on complementary strand
    if (strand == "+") {
      search_base <- complement_base(canonical_base)
      reverse_scan <- TRUE
    } else {
      search_base <- canonical_base
      reverse_scan <- FALSE
    }
  }

  # Find all 1-based positions of the search base in seq
  seq_chars <- strsplit(seq, "")[[1]]
  canonical_positions <- which(toupper(seq_chars) == toupper(search_base))

  if (reverse_scan) {
    # Reverse the order: MM deltas are 5'->3' of original molecule
    canonical_positions <- rev(canonical_positions)
  }


  # --- 4. Apply delta offsets to find modified base indices ---
  current <- 0L
  modified_indices <- integer(length(deltas))
  for (i in seq_along(deltas)) {
    current <- current + deltas[i] + 1L
    modified_indices[i] <- current
  }

  # Guard against out-of-bounds indices
  valid <- modified_indices >= 1L & modified_indices <= length(canonical_positions)
  if (!any(valid)) {
    return(empty_result)
  }

  modified_indices <- modified_indices[valid]
  ml_values <- ml_values[valid]

  modified_seq_positions <- canonical_positions[modified_indices]


  # --- 5. Convert sequence positions to genomic positions ---
  ref_positions <- seq_to_ref(cigar, pos, modified_seq_positions)


  # --- 6. Build result, filtering out NA positions ---
  keep <- !is.na(ref_positions)
  if (!any(keep)) {
    return(empty_result)
  }

  data.frame(
    position = ref_positions[keep],
    mod_prob = ml_values[keep] / 255,
    stringsAsFactors = FALSE
  )
}
