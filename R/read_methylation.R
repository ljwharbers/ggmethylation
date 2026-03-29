#' Read methylation data from a modBAM file
#'
#' Parses a BAM file with MM/ML tags and extracts base modification data for
#' a given genomic region. This is the main entry point for data extraction
#' (Layer 1 of the package).
#'
#' @param bam Character. Path to a BAM file containing MM/ML modification tags.
#' @param region Character. Genomic region string, e.g. `"chr1:1000-2000"`.
#' @param mod_code Character. Modification code from the MM tag (default `"m"`
#'   for 5-methylcytosine).
#' @param group_tag Character or NULL. BAM tag used to group reads (e.g.
#'   `"HP"` for haplotype, `"RG"` for read group). NULL disables grouping.
#' @param max_reads Integer. Maximum number of reads to return (default 200).
#'   If more reads overlap the region, a random subset is kept.
#'
#' @return A `methylation_data` object (S3 list) with elements:
#'   \describe{
#'     \item{reads}{Data.frame with columns `read_name`, `start`, `end`,
#'       `strand`, and optionally `group`.}
#'     \item{sites}{Data.frame with columns `position`, `mod_prob`,
#'       `read_name`, `mod_code`, and optionally `group`.}
#'     \item{region}{A [GenomicRanges::GRanges] object for the queried region.}
#'     \item{mod_code}{The modification code used.}
#'     \item{group_tag}{The grouping tag used, or NULL.}
#'   }
#'
#' @examples
#' \dontrun{
#' md <- read_methylation("sample.bam", "chr1:1000-2000")
#' md <- read_methylation("sample.bam", "chr1:1000-2000",
#'   group_tag = "HP", max_reads = 100
#' )
#' }
#'
#' @export
read_methylation <- function(bam, region, mod_code = "m", group_tag = NULL,
                             snv_position = NULL, ref_base = NULL,
                             alt_base = NULL, max_reads = 200L) {
  # --- 1. Validate inputs ---
  mod_code <- as.character(mod_code)          # coerce in case of factor
  if (length(mod_code) == 0L || any(!nzchar(mod_code))) {
    stop("'mod_code' must be a non-empty character vector of modification codes.",
         call. = FALSE)
  }
  mod_code <- unique(mod_code)                # drop accidental duplicates

  if (!is.null(group_tag) && !is.null(snv_position)) {
    stop("'group_tag' and 'snv_position' are mutually exclusive.", call. = FALSE)
  }
  if (!is.null(snv_position) && (is.null(ref_base) || is.null(alt_base))) {
    stop("'ref_base' and 'alt_base' must be provided when 'snv_position' is set.",
         call. = FALSE)
  }
  validate_bam_index(bam)
  parsed <- parse_region(region)

  # --- 2. Set up Rsamtools query ---
  gr <- GenomicRanges::GRanges(
    seqnames = parsed$chrom,
    ranges = IRanges::IRanges(start = parsed$start, end = parsed$end)
  )

  what <- c("qname", "pos", "cigar", "strand", "seq")

  tag_names <- c("MM", "ML")
  if (!is.null(group_tag)) tag_names <- c(tag_names, group_tag)

  param <- Rsamtools::ScanBamParam(
    which = gr,
    what = what,
    tag = tag_names
  )

  bam_data <- Rsamtools::scanBam(bam, param = param)[[1]]

  # --- 3. Check for reads ---
  if (is.null(bam_data$qname) || length(bam_data$qname) == 0L) {
    warning("No reads found in region '", region, "'.", call. = FALSE)
    return(empty_methylation_data(gr, mod_code, group_tag))
  }

  # --- 4. Build reads data.frame ---
  ref_widths <- GenomicAlignments::cigarWidthAlongReferenceSpace(bam_data$cigar)

  reads <- data.frame(
    read_name = bam_data$qname,
    start = bam_data$pos,
    end = bam_data$pos + ref_widths - 1L,
    strand = as.character(bam_data$strand),
    stringsAsFactors = FALSE
  )

  # Add group column
  if (!is.null(group_tag)) {
    group_values <- bam_data$tag[[group_tag]]
    if (is.null(group_values) || all(is.na(group_values))) {
      warning(
        "Group tag '", group_tag, "' not found in any read. ",
        "Ignoring grouping.", call. = FALSE
      )
      reads$group <- NA_character_
      group_tag <- NULL
    } else {
      reads$group <- as.character(group_values)
    }
  }

  # --- 5. SNV-based grouping ---
  # bam_indices maps filtered-reads rows back to original bam_data positions
  bam_indices <- seq_len(nrow(reads))

  if (!is.null(snv_position)) {
    n_all <- nrow(reads)
    snv_bases <- character(n_all)
    for (j in seq_len(n_all)) {
      seq_str <- as.character(bam_data$seq[[j]])
      q_pos <- ref_to_seq(bam_data$cigar[j], bam_data$pos[j], snv_position)
      snv_bases[j] <- if (!is.na(q_pos) && q_pos >= 1L && q_pos <= nchar(seq_str)) {
        toupper(substr(seq_str, q_pos, q_pos))
      } else {
        NA_character_
      }
    }
    reads$group <- ifelse(
      snv_bases == toupper(ref_base), "REF",
      ifelse(snv_bases == toupper(alt_base), "ALT", NA_character_)
    )
    snv_keep <- which(!is.na(reads$group))
    if (length(snv_keep) == 0L) {
      warning("No reads carry REF or ALT at snv_position ", snv_position, ".",
              call. = FALSE)
      return(empty_methylation_data(gr, mod_code, "SNV", snv_position))
    }
    bam_indices <- snv_keep
    reads <- reads[snv_keep, , drop = FALSE]
    rownames(reads) <- NULL
    group_tag <- "SNV"
  }

  # --- 6. Downsample if needed ---
  n_reads <- nrow(reads)
  keep_local <- seq_len(n_reads)

  if (n_reads > max_reads) {
    keep_local <- sort(sample(n_reads, max_reads))
    reads <- reads[keep_local, , drop = FALSE]
    rownames(reads) <- NULL
  }

  # Translate to original bam_data indices
  keep <- bam_indices[keep_local]

  # --- 7. Parse MM/ML tags for each read ---
  sites_list <- vector("list", length(keep))

  for (j in seq_along(keep)) {
    idx       <- keep[j]
    seq_str   <- as.character(bam_data$seq[[idx]])
    mm        <- bam_data$tag$MM[idx]
    ml        <- bam_data$tag$ML[[idx]]

    code_results <- vector("list", length(mod_code))
    for (ci in seq_along(mod_code)) {
      result <- parse_mm_ml(
        seq      = seq_str,
        mm_tag   = mm,
        ml_tag   = ml,
        mod_code = mod_code[ci],
        strand   = reads$strand[j],
        cigar    = bam_data$cigar[idx],
        pos      = bam_data$pos[idx]
      )
      result$read_name <- if (nrow(result) > 0L) reads$read_name[j] else character(0L)
      result$mod_code  <- if (nrow(result) > 0L) mod_code[ci]       else character(0L)
      code_results[[ci]] <- result
    }
    sites_list[[j]] <- do.call(rbind, code_results)
  }

  sites <- do.call(rbind, sites_list)

  if (is.null(sites) || nrow(sites) == 0L) {
    sites <- data.frame(
      position  = integer(0L),
      mod_prob  = numeric(0L),
      read_name = character(0L),
      mod_code  = character(0L),
      stringsAsFactors = FALSE
    )
  }

  # Add group column to sites by joining on read_name
  if (!is.null(group_tag)) {
    sites$group <- reads$group[match(sites$read_name, reads$read_name)]
  }

  # --- 7. Clip reads to region ---
  reads$start <- pmax(reads$start, parsed$start)
  reads$end <- pmin(reads$end, parsed$end)

  # --- 8. Filter sites to region ---
  sites <- sites[
    sites$position >= parsed$start & sites$position <= parsed$end,
    ,
    drop = FALSE
  ]
  rownames(sites) <- NULL

  # --- 9. Return methylation_data object ---
  # (mod_code column is already populated per-code in the inner loop above)
  structure(
    list(
      reads = reads,
      sites = sites,
      region = gr,
      mod_code = mod_code,
      group_tag = group_tag,
      snv_position = snv_position
    ),
    class = "methylation_data"
  )
}

#' Create an empty methylation_data object
#'
#' @param gr A [GenomicRanges::GRanges] object.
#' @param mod_code Character. Modification code.
#' @param group_tag Character or NULL. Grouping tag.
#'
#' @return An empty `methylation_data` object.
#'
#' @keywords internal
empty_methylation_data <- function(gr, mod_code, group_tag, snv_position = NULL) {
  reads <- data.frame(
    read_name = character(0L),
    start = integer(0L),
    end = integer(0L),
    strand = character(0L),
    stringsAsFactors = FALSE
  )

  sites <- data.frame(
    position = integer(0L),
    mod_prob = numeric(0L),
    read_name = character(0L),
    mod_code = character(0L),
    stringsAsFactors = FALSE
  )

  if (!is.null(group_tag)) {
    reads$group <- character(0L)
    sites$group <- character(0L)
  }

  structure(
    list(
      reads = reads,
      sites = sites,
      region = gr,
      mod_code = mod_code,
      group_tag = group_tag,
      snv_position = snv_position
    ),
    class = "methylation_data"
  )
}

#' Print a methylation_data object
#'
#' @param x A `methylation_data` object.
#' @param ... Additional arguments (ignored).
#'
#' @return `x`, invisibly.
#'
#' @export
print.methylation_data <- function(x, ...) {
  chrom <- as.character(GenomicRanges::seqnames(x$region))
  start <- GenomicRanges::start(x$region)
  end <- GenomicRanges::end(x$region)

  cat("methylation_data object\n")
  cat(sprintf("Region: %s:%d-%d\n", chrom, start, end))
  cat(sprintf("Reads: %d\n", nrow(x$reads)))
  cat(sprintf("Modification sites: %d\n", nrow(x$sites)))
  cat(sprintf("Modification code(s): %s\n", paste(x$mod_code, collapse = ", ")))

  if (!is.null(x$group_tag)) {
    n_groups <- length(unique(x$reads$group[!is.na(x$reads$group)]))
    cat(sprintf("Group tag: %s (%d groups)\n", x$group_tag, n_groups))
  }

  invisible(x)
}
