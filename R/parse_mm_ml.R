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
#' @return A named `list` with two `data.frame` elements:
#'   \describe{
#'     \item{sites}{Reference-aligned modifications. Columns `position`
#'       (integer, 1-based genomic position) and `mod_prob` (numeric,
#'       ML value / 255). When the MM entry uses the `.` implicit-unmodified
#'       flag, also includes rows for unlisted canonical positions with
#'       `mod_prob = 0`.}
#'     \item{insertion_sites}{Modifications whose query base falls inside
#'       a CIGAR `I` interval (no reference position). Columns `query_pos`
#'       (integer, 1-based position in the read sequence) and `mod_prob`
#'       (numeric). Modifications inside soft/hard clips are still dropped.
#'       Implicit-zero positions (`.` flag) inside insertions are not emitted.}
#'   }
#'   Returns both frames empty (zero rows, same columns) when the
#'   requested modification is not present.
#'
#' @keywords internal
parse_mm_ml <- function(seq, mm_tag, ml_tag, mod_code, strand, cigar, pos) {
  empty_sites <- data.frame(
    position = integer(0),
    mod_prob = numeric(0),
    stringsAsFactors = FALSE
  )
  empty_insertion_sites <- data.frame(
    query_pos = integer(0),
    mod_prob  = numeric(0),
    stringsAsFactors = FALSE
  )
  empty_result <- list(sites = empty_sites, insertion_sites = empty_insertion_sites)


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

    # Split on comma: first element is the spec (e.g. "C+m", "C+m?", "C+m."),
    # rest are skip-count deltas.
    parts  <- strsplit(entry, ",")[[1]]
    spec   <- parts[1]
    deltas <- if (length(parts) > 1) as.integer(parts[2:length(parts)]) else integer(0)

    # Extract the optional implicit-base flag ('?' or '.') from the end of the
    # spec. Per SAM spec: '.' means unlisted canonical bases are implicitly
    # unmodified (mod_prob = 0); '?' or absent means no information.
    last_char <- if (nchar(spec) > 0L) substr(spec, nchar(spec), nchar(spec)) else ""
    flag_char <- if (last_char %in% c("?", ".")) last_char else ""
    if (nzchar(flag_char)) spec <- substr(spec, 1L, nchar(spec) - 1L)

    # Parse the spec: canonical_base + strand_char + code
    # e.g., "C+m" -> base="C", mm_strand="+", code="m"
    canonical_base <- substr(spec, 1, 1)
    mm_strand <- substr(spec, 2, 2)
    code <- substring(spec, 3)

    parsed_entries[[j]] <- list(
      canonical_base = canonical_base,
      mm_strand = mm_strand,
      code = code,
      flag = flag_char,
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

  # With '.' flag an empty delta list means every canonical base is implicitly
  # unmodified; we still need to emit them. Only bail out early for '?' / none.
  if (length(deltas) == 0 && !identical(target$flag, ".")) {
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
    # With '.' flag, even if no listed modifications are in-bounds we can still
    # emit implicit-zero rows for all canonical positions.
    if (!identical(target$flag, ".") || length(canonical_positions) == 0L) {
      return(empty_result)
    }
    modified_indices <- integer(0L)
    ml_values        <- integer(0L)
  } else {
    modified_indices <- modified_indices[valid]
    ml_values        <- ml_values[valid]
  }

  modified_seq_positions <- canonical_positions[modified_indices]

  # --- 4b. Unwalked canonical positions for '.' (implicit-unmodified) flag ---
  # When the MM entry uses '.', any canonical base NOT listed in the delta walk
  # is implicitly unmodified (mod_prob = 0). Compute those indices now; they
  # will be emitted as explicit zeros in Section 6.
  target_flag <- target$flag
  unwalked_indices <- if (identical(target_flag, ".") && length(canonical_positions) > 0L) {
    setdiff(seq_along(canonical_positions), modified_indices)
  } else {
    integer(0L)
  }

  # --- 5. Convert sequence positions to genomic positions ---
  ref_positions <- seq_to_ref(cigar, pos, modified_seq_positions)

  # Map implicit-zero positions to reference coordinates (NA for I/S positions)
  zero_ref_positions <- if (length(unwalked_indices) > 0L) {
    seq_to_ref(cigar, pos, canonical_positions[unwalked_indices])
  } else {
    integer(0L)
  }


  # --- 6. Build result: ref-aligned sites and insertion sites ---
  ref_mask <- !is.na(ref_positions)

  sites_df <- if (any(ref_mask)) {
    data.frame(
      position = ref_positions[ref_mask],
      mod_prob = ml_values[ref_mask] / 255,
      stringsAsFactors = FALSE
    )
  } else {
    empty_sites
  }

  # Append implicit-zero rows for '.' flag positions that map to reference.
  # Positions inside insertions/clips return NA from seq_to_ref and are dropped.
  if (length(zero_ref_positions) > 0L) {
    zero_ref_mask <- !is.na(zero_ref_positions)
    if (any(zero_ref_mask)) {
      zero_df <- data.frame(
        position = zero_ref_positions[zero_ref_mask],
        mod_prob = 0,
        stringsAsFactors = FALSE
      )
      sites_df <- if (nrow(sites_df) > 0L) rbind(sites_df, zero_df) else zero_df
    }
  }

  # Classify NA positions: insertion (CIGAR I) vs soft/hard clip (dropped)
  insertion_sites_df <- empty_insertion_sites
  na_mask <- !ref_mask
  if (any(na_mask)) {
    cf <- decompose_cigar(cigar, pos)
    i_rows <- cf[cf$type == "I", , drop = FALSE]
    if (nrow(i_rows) > 0L) {
      na_qpos <- modified_seq_positions[na_mask]
      na_ml   <- ml_values[na_mask]
      in_ins  <- vapply(na_qpos, function(qp) {
        any(qp >= i_rows$query_start & qp <= i_rows$query_end)
      }, logical(1L))
      if (any(in_ins)) {
        insertion_sites_df <- data.frame(
          query_pos = na_qpos[in_ins],
          mod_prob  = na_ml[in_ins] / 255,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  list(sites = sites_df, insertion_sites = insertion_sites_df)
}
