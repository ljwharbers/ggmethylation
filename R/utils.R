# Internal utility functions for ggmethylation

#' Parse a genomic region string
#'
#' Parses a region string of the form `"chr1:1000-2000"` into its components.
#' Commas in numeric positions (e.g., `"chr1:1,000-2,000"`) are stripped
#' before parsing.
#'
#' @param region Character string specifying a genomic region, e.g.
#'   `"chr1:1000-2000"`.
#'
#' @return A named list with elements:
#'   \describe{
#'     \item{chrom}{Character. The chromosome name.}
#'     \item{start}{Integer. The start coordinate.}
#'     \item{end}{Integer. The end coordinate.}
#'   }
#'
#' @keywords internal
parse_region <- function(region) {
  # Strip commas from the string before validation
  clean <- gsub(",", "", region)

  if (!grepl("^[\\w.]+:\\d+-\\d+$", clean, perl = TRUE)) {
    stop(
      "Invalid region format: '", region, "'. ",
      "Expected format: 'chr:start-end' (e.g., 'chr1:1000-2000').",
      call. = FALSE
    )
  }

  parts <- strsplit(clean, "[:-]")[[1]]
  list(
    chrom = parts[1],
    start = as.integer(parts[2]),
    end   = as.integer(parts[3])
  )
}

#' Validate that a BAM file and its index exist
#'
#' Checks that the BAM file exists and that a corresponding `.bai` index file
#' can be found. Looks for the index at both `{bam}.bai` and the path obtained
#' by replacing the `.bam` extension with `.bam.bai`.
#'
#' @param bam Character. Path to a BAM file.
#'
#' @return The path to the BAI index file (invisibly).
#'
#' @keywords internal
validate_bam_index <- function(bam) {
  if (!file.exists(bam)) {
    stop("BAM file does not exist: '", bam, "'.", call. = FALSE)
  }

  bai_candidates <- c(
    paste0(bam, ".bai"),
    sub("\\.bam$", ".bai", bam)
  )

  for (bai in bai_candidates) {
    if (file.exists(bai)) {
      return(bai)
    }
  }

  stop(
    "BAM index not found for '", bam, "'. ",
    "Please create one with: samtools index ", bam,
    call. = FALSE
  )
}

#' Return the complement of a DNA base
#'
#' Maps A to T, T to A, C to G, and G to C.
#'
#' @param base Single character. One of `"A"`, `"T"`, `"C"`, or `"G"`.
#'
#' @return Single character. The complementary base.
#'
#' @keywords internal
#' Map a reference position to a query (read) position via CIGAR
#'
#' Walks the CIGAR string to find the 1-based query position that corresponds
#' to a given 1-based reference position. Returns `NA` if the position falls
#' in a deletion, is not covered by the read, or is in a soft/hard clip.
#'
#' @param cigar Character. A CIGAR string (e.g., `"5S10M2D3M"`).
#' @param ref_start Integer. The 1-based reference position of the first
#'   aligned base (the BAM POS field).
#' @param target_ref_pos Integer. The 1-based reference position to look up.
#'
#' @return Integer. The 1-based query position, or `NA_integer_`.
#'
#' @keywords internal
ref_to_seq <- function(cigar, ref_start, target_ref_pos) {
  ops  <- regmatches(cigar, gregexpr("[A-Z=]", cigar))[[1]]
  lens <- as.integer(regmatches(cigar, gregexpr("\\d+", cigar))[[1]])

  q_pos <- 1L
  r_pos <- ref_start

  for (i in seq_along(ops)) {
    op  <- ops[i]
    len <- lens[i]

    if (op %in% c("M", "=", "X")) {
      if (target_ref_pos >= r_pos && target_ref_pos < r_pos + len) {
        return(q_pos + (target_ref_pos - r_pos))
      }
      q_pos <- q_pos + len
      r_pos <- r_pos + len
    } else if (op %in% c("D", "N")) {
      if (target_ref_pos >= r_pos && target_ref_pos < r_pos + len) {
        return(NA_integer_)
      }
      r_pos <- r_pos + len
    } else if (op %in% c("I", "S")) {
      q_pos <- q_pos + len
    }
    # H: hard clip — consumes neither
  }

  NA_integer_
}

#' Return the complement of a DNA base
#'
#' Maps A to T, T to A, C to G, and G to C.
#'
#' @param base Single character. One of `"A"`, `"T"`, `"C"`, or `"G"`.
#'
#' @return Single character. The complementary base.
#'
#' @keywords internal
complement_base <- function(base) {
  comp <- c(A = "T", T = "A", C = "G", G = "C")
  result <- comp[toupper(base)]
  if (is.na(result)) {
    stop("Unknown base: '", base, "'.", call. = FALSE)
  }
  as.character(result)
}

#' Decompose a CIGAR string into individual operations
#'
#' Walks the CIGAR string and returns a data.frame describing each operation
#' with reference and query coordinate ranges.
#'
#' @param cigar Character. A CIGAR string (e.g., `"5S10M2I3D5M"`).
#' @param pos Integer. The 1-based reference position of the first aligned
#'   base (the BAM POS field).
#'
#' @return A data.frame with columns:
#'   \describe{
#'     \item{type}{Character. CIGAR operation: `"M"`, `"I"`, `"D"`, `"S"`,
#'       `"H"`, `"N"`, `"="`, or `"X"`.}
#'     \item{ref_start}{Integer. Reference start position (`NA` for I/S/H).}
#'     \item{ref_end}{Integer. Reference end position (`NA` for I/S/H).}
#'     \item{query_start}{Integer. Query start position (`NA` for D/N/H).}
#'     \item{query_end}{Integer. Query end position (`NA` for D/N/H).}
#'     \item{length}{Integer. Operation length.}
#'   }
#'
#' @keywords internal
#' Parse an SA (Supplementary Alignment) BAM auxiliary tag
#'
#' Splits a semicolon-delimited SA tag string into a data.frame, one row per
#' supplementary alignment entry.  The SA tag format is:
#' `"rname,pos,strand,CIGAR,mapQ,NM;"` (trailing semicolon).
#'
#' @param sa_string Character. The SA tag value (e.g.
#'   `"chr5,45000,+,50M,60,0;"`), or `NULL`/`NA`.
#'
#' @return A data.frame with columns `rname` (character), `pos` (integer),
#'   `strand` (character), `cigar` (character), `mapq` (integer), `nm`
#'   (integer).  Returns a zero-row data.frame for NULL/NA/empty input.
#'
#' @keywords internal
parse_sa_tag <- function(sa_string) {
  empty <- data.frame(
    rname  = character(0L),
    pos    = integer(0L),
    strand = character(0L),
    cigar  = character(0L),
    mapq   = integer(0L),
    nm     = integer(0L),
    stringsAsFactors = FALSE
  )
  if (is.null(sa_string) || length(sa_string) == 0L ||
        is.na(sa_string) || !nzchar(sa_string)) {
    return(empty)
  }
  entries <- strsplit(sa_string, ";", fixed = TRUE)[[1L]]
  entries <- entries[nzchar(entries)]
  if (length(entries) == 0L) return(empty)

  rows <- lapply(entries, function(e) {
    fields <- strsplit(e, ",", fixed = TRUE)[[1L]]
    if (length(fields) < 6L) return(NULL)
    data.frame(
      rname  = fields[1L],
      pos    = as.integer(fields[2L]),
      strand = fields[3L],
      cigar  = fields[4L],
      mapq   = as.integer(fields[5L]),
      nm     = as.integer(fields[6L]),
      stringsAsFactors = FALSE
    )
  })
  rows <- Filter(Negate(is.null), rows)
  if (length(rows) == 0L) return(empty)
  do.call(rbind, rows)
}

decompose_cigar <- function(cigar, pos) {
  ops  <- regmatches(cigar, gregexpr("[A-Z=]", cigar))[[1]]
  lens <- as.integer(regmatches(cigar, gregexpr("\\d+", cigar))[[1]])

  n <- length(ops)
  type        <- character(n)
  ref_start   <- rep(NA_integer_, n)
  ref_end     <- rep(NA_integer_, n)
  query_start <- rep(NA_integer_, n)
  query_end   <- rep(NA_integer_, n)
  op_length   <- integer(n)

  ref_offset   <- as.integer(pos)
  query_offset <- 1L

  for (i in seq_len(n)) {
    op  <- ops[i]
    len <- lens[i]
    type[i]      <- op
    op_length[i] <- len

    if (op %in% c("M", "=", "X")) {
      ref_start[i]   <- ref_offset
      ref_end[i]     <- ref_offset + len - 1L
      query_start[i] <- query_offset
      query_end[i]   <- query_offset + len - 1L
      ref_offset     <- ref_offset + len
      query_offset   <- query_offset + len
    } else if (op == "I") {
      ref_start[i]   <- ref_offset  # insertion point
      query_start[i] <- query_offset
      query_end[i]   <- query_offset + len - 1L
      query_offset   <- query_offset + len
    } else if (op %in% c("D", "N")) {
      ref_start[i] <- ref_offset
      ref_end[i]   <- ref_offset + len - 1L
      ref_offset   <- ref_offset + len
    } else if (op == "S") {
      query_start[i] <- query_offset
      query_end[i]   <- query_offset + len - 1L
      query_offset   <- query_offset + len
    } else if (op == "H") {
      # Hard clip: consumes neither reference nor query
    }
  }

  data.frame(
    type        = type,
    ref_start   = ref_start,
    ref_end     = ref_end,
    query_start = query_start,
    query_end   = query_end,
    length      = op_length,
    stringsAsFactors = FALSE
  )
}

#' Compute reference-space width from a CIGAR string
#'
#' Vectorised over a character vector of CIGAR strings.  For each CIGAR,
#' sums the lengths of operations that consume reference bases
#' (M, =, X, D, N).
#'
#' @param cigar Character vector of CIGAR strings.
#'
#' @return Integer vector of reference-space widths.
#'
#' @keywords internal
cigar_ref_width <- function(cigar) {
  vapply(cigar, function(cig) {
    if (is.na(cig) || cig == "*") return(0L)
    ops  <- regmatches(cig, gregexpr("[A-Z=]", cig))[[1]]
    lens <- as.integer(regmatches(cig, gregexpr("\\d+", cig))[[1]])
    sum(lens[ops %in% c("M", "=", "X", "D", "N")])
  }, integer(1L), USE.NAMES = FALSE)
}

#' Detect which side(s) of a read are soft/hard clipped
#'
#' Vectorised over a character vector of CIGAR strings.  Checks whether
#' the first and/or last CIGAR operation is S or H.
#'
#' @param cigar Character vector of CIGAR strings.
#'
#' @return Character vector with values `"left"`, `"right"`, `"both"`,
#'   or `NA_character_`.
#'
#' @keywords internal
detect_clip_side <- function(cigar) {
  vapply(cigar, function(cig) {
    if (is.na(cig) || cig == "*") return(NA_character_)
    ops <- regmatches(cig, gregexpr("[A-Z=]", cig))[[1]]
    left  <- ops[1] %in% c("S", "H")
    right <- ops[length(ops)] %in% c("S", "H")
    if (left && right) "both"
    else if (left)      "left"
    else if (right)     "right"
    else                NA_character_
  }, character(1L), USE.NAMES = FALSE)
}
