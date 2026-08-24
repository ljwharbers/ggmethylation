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

#' Parse a CIGAR string into operation codes and lengths
#'
#' @param cigar Character. A single CIGAR string.
#'
#' @return A named list with elements `ops` (character vector of operation
#'   codes) and `lens` (integer vector of lengths).
#'
#' @keywords internal
split_cigar <- function(cigar) {
  list(
    ops  = regmatches(cigar, gregexpr("[A-Z=]", cigar))[[1]],
    lens = as.integer(regmatches(cigar, gregexpr("\\d+", cigar))[[1]])
  )
}

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
  parsed_cigar <- split_cigar(cigar)
  ops  <- parsed_cigar$ops
  lens <- parsed_cigar$lens

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
decompose_cigar <- function(cigar, pos) {
  parsed_cigar <- split_cigar(cigar)
  ops  <- parsed_cigar$ops
  lens <- parsed_cigar$lens

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
    parsed_cigar <- split_cigar(cig)
    sum(parsed_cigar$lens[parsed_cigar$ops %in% c("M", "=", "X", "D", "N")])
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
    ops   <- split_cigar(cig)$ops
    left  <- ops[1] %in% c("S", "H")
    right <- ops[length(ops)] %in% c("S", "H")
    if (left && right) "both"
    else if (left)      "left"
    else if (right)     "right"
    else                NA_character_
  }, character(1L), USE.NAMES = FALSE)
}

#' Query-space extent of an alignment in original read orientation
#'
#' Computes where an alignment sits along the read as it came off the
#' sequencer, so that a primary alignment and the entries of its `SA` tag can
#' be compared on a common axis.
#'
#' Hard-clipped bases are counted towards the full read length, which is what
#' makes the comparison valid: a primary alignment normally soft-clips the
#' portion aligned elsewhere while its supplementary counterpart hard-clips the
#' portion aligned here, yet both CIGARs then describe the same total read
#' length.  A reverse-strand alignment is stored reverse-complemented relative
#' to the read, so its *trailing* clip is the read-5' offset.
#'
#' @param cigar Character. A single CIGAR string.
#' @param strand Character. `"+"` or `"-"`.
#'
#' @return A named list with `start`, `end` (0-based, inclusive, in original
#'   read orientation) and `qlen` (full read length), or `NULL` when the CIGAR
#'   is missing, unaligned, or has no aligned query bases.
#'
#' @keywords internal
.query_extent <- function(cigar, strand) {
  if (length(cigar) != 1L || is.na(cigar) || cigar == "*") return(NULL)
  parsed <- split_cigar(cigar)
  ops    <- parsed$ops
  lens   <- parsed$lens
  n      <- length(ops)
  if (n == 0L || n != length(lens) || anyNA(lens)) return(NULL)

  lead  <- if (ops[1L] %in% c("S", "H")) lens[1L] else 0L
  trail <- if (n > 1L && ops[n] %in% c("S", "H")) lens[n] else 0L
  qlen  <- sum(lens[ops %in% c("M", "I", "S", "=", "X", "H")])
  aln   <- qlen - lead - trail
  if (aln <= 0L) return(NULL)

  start <- if (identical(strand, "-")) trail else lead
  list(start = start, end = start + aln - 1L, qlen = qlen)
}

#' Build the supplementary-alignment columns for a set of alignments
#'
#' Vectorised over alignments.  For each one, [sa_partner_sides()] decides which
#' reference flank every `SA` partner joins; the best partner (by MAPQ) is then
#' recorded overall and per flank.
#'
#' `sa_chrom` / `sa_pos` are populated from the best partner regardless of
#' whether its flank could be determined, so BND matching still works for reads
#' whose junction side is undecidable.
#'
#' @param cigar Character vector of CIGAR strings, one per alignment.
#' @param strand Character vector of strands, one per alignment.
#' @param sa_tags List or character vector of raw `SA` tag values, one per
#'   alignment, or `NULL` when the BAM carries no `SA` tag at all.
#'
#' @return A data.frame with `length(cigar)` rows and columns `sa_chrom`,
#'   `sa_pos`, `sa_side`, `sa_chrom_left`, `sa_pos_left`, `sa_chrom_right`,
#'   `sa_pos_right`. All `NA` for alignments with no usable `SA` entry.
#'
#' @keywords internal
sa_columns <- function(cigar, strand, sa_tags) {
  n <- length(cigar)
  out <- data.frame(
    sa_chrom       = rep(NA_character_, n),
    sa_pos         = rep(NA_integer_,   n),
    sa_side        = rep(NA_character_, n),
    sa_chrom_left  = rep(NA_character_, n),
    sa_pos_left    = rep(NA_integer_,   n),
    sa_chrom_right = rep(NA_character_, n),
    sa_pos_right   = rep(NA_integer_,   n),
    stringsAsFactors = FALSE
  )
  if (is.null(sa_tags) || n == 0L) return(out)

  for (i in seq_len(n)) {
    sided <- sa_partner_sides(cigar[i], strand[i], sa_tags[[i]])
    if (nrow(sided) == 0L) next

    best <- .best_by_mapq(sided)
    out$sa_chrom[i] <- sided$rname[best]
    out$sa_pos[i]   <- sided$pos[best]

    for (side in c("left", "right")) {
      rows <- which(!is.na(sided$side) & sided$side == side)
      if (length(rows) == 0L) next
      b <- rows[.best_by_mapq(sided[rows, , drop = FALSE])]
      out[[paste0("sa_chrom_", side)]][i] <- sided$rname[b]
      out[[paste0("sa_pos_",   side)]][i] <- sided$pos[b]
    }

    has_left  <- !is.na(out$sa_chrom_left[i])
    has_right <- !is.na(out$sa_chrom_right[i])
    out$sa_side[i] <- if (has_left && has_right) {
      "both"
    } else if (has_left) {
      "left"
    } else if (has_right) {
      "right"
    } else {
      NA_character_
    }
  }

  out
}

#' Index of the highest-MAPQ row of an SA entry table
#'
#' Ties resolve to the first row; missing MAPQ sorts last.
#'
#' @param entries A data.frame with a `mapq` column and at least one row.
#'
#' @return Integer row index.
#'
#' @keywords internal
.best_by_mapq <- function(entries) {
  which.max(ifelse(is.na(entries$mapq), -1L, entries$mapq))
}

#' Locate supplementary-alignment partners relative to a primary alignment
#'
#' For a chimeric read, each `SA` entry joins the primary alignment at exactly
#' one point in read coordinates.  This helper works out, for every entry,
#' which *reference* flank of the primary alignment that junction falls on.
#'
#' The clipped side of a CIGAR is not a usable proxy: long ONT/PacBio reads are
#' routinely soft-clipped at both ends by adapter and quality trimming, so
#' [detect_clip_side()] reports `"both"` for reads that have only a single
#' supplementary partner.  Instead, the entry's position along the read is
#' compared with the primary's via `.query_extent()`: an entry lying downstream
#' of the primary in read space joins at the primary's read-3' end, which is
#' its reference right edge on the forward strand and its left edge on the
#' reverse strand.
#'
#' @param cigar Character. CIGAR of the primary (this) alignment.
#' @param strand Character. Strand of the primary alignment, `"+"` or `"-"`.
#' @param sa_string Character. The raw `SA` tag value, or `NA`.
#'
#' @return A data.frame with one row per parseable `SA` entry and columns
#'   `rname`, `pos`, `mapq`, and `side` (`"left"`, `"right"`, or `NA` when the
#'   junction side cannot be determined). Zero rows when there are no entries.
#'
#' @keywords internal
sa_partner_sides <- function(cigar, strand, sa_string) {
  entries <- parse_sa_tag(sa_string)
  out <- data.frame(
    rname = entries$rname,
    pos   = entries$pos,
    mapq  = entries$mapq,
    side  = rep(NA_character_, nrow(entries)),
    stringsAsFactors = FALSE
  )
  if (nrow(out) == 0L) return(out)

  prim <- .query_extent(cigar, strand)
  if (is.null(prim)) return(out)

  # The primary's read-3' end maps to this reference edge.
  three_prime_side <- if (identical(strand, "-")) "left" else "right"
  five_prime_side  <- if (identical(three_prime_side, "right")) "left" else "right"

  for (i in seq_len(nrow(entries))) {
    ext <- .query_extent(entries$cigar[i], entries$strand[i])
    if (is.null(ext)) next
    # An exact tie carries no directional information.
    if (ext$start == prim$start) next
    out$side[i] <- if (ext$start > prim$start) three_prime_side else five_prime_side
  }

  out
}

#' Convert a region string to a GRanges object
#'
#' Wraps [parse_region()] and constructs a [GenomicRanges::GRanges] from the
#' result.
#'
#' @param region Character. A genomic region string, e.g. `"chr1:1000-2000"`.
#'
#' @return A [GenomicRanges::GRanges] object with one range.
#'
#' @keywords internal
region_to_granges <- function(region) {
  parsed <- parse_region(region)
  GenomicRanges::GRanges(
    seqnames = parsed$chrom,
    ranges   = IRanges::IRanges(start = parsed$start, end = parsed$end)
  )
}


#' Validate `sort_by` against the available read columns
#'
#' `plot_methylation()` sorts reads with `order()` over columns pulled out of
#' `$reads` by name. An unknown name yields `NULL`, and `order(NULL)` returns
#' `integer(0)` -- which silently drops every read and produces an empty plot
#' rather than an error. This helper turns that into an explicit failure.
#'
#' @param sort_by Character vector of column names to sort by.
#' @param reads The `$reads` data frame the names must exist in.
#'
#' @return `sort_by`, invisibly, when every name is valid.
#'
#' @keywords internal
.validate_sort_by <- function(sort_by, reads) {
  if (is.null(sort_by)) {
    return(invisible(sort_by))
  }
  if (!is.character(sort_by)) {
    stop("`sort_by` must be a character vector of column names.", call. = FALSE)
  }

  missing <- setdiff(sort_by, names(reads))
  if (length(missing) > 0L) {
    stop(
      "Unknown `sort_by` column", if (length(missing) > 1L) "s" else "", ": ",
      paste0("\"", missing, "\"", collapse = ", "), ".\n",
      "Available columns: ",
      paste0("\"", names(reads), "\"", collapse = ", "), ".",
      call. = FALSE
    )
  }

  invisible(sort_by)
}
