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
