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
#' @param per_group_downsample Logical. When `TRUE` and grouping is active
#'   (via `group_tag` or `snv_position`), the `max_reads` cap is applied
#'   independently per group. When `FALSE` (default), the existing global
#'   cap behaviour is unchanged.
#' @param min_mapq Integer. Minimum mapping quality (MAPQ) threshold
#'   (default \code{0L}, no filtering). Reads with MAPQ below this value,
#'   or with missing MAPQ, are excluded before MM/ML parsing.
#' @param strand_filter Character vector. Which strands to include. One or
#'   both of \code{"+"} and \code{"-"} (default \code{c("+", "-")}, both
#'   strands). Use \code{"+"} or \code{"-"} alone to restrict to a single
#'   strand.
#' @param min_read_length Integer. Minimum read length in reference-space
#'   base pairs (default \code{0L}, no filtering). Reads shorter than this
#'   value are excluded before MM/ML parsing.
#' @param drop_na_group Logical. When `TRUE` and `group_tag` is set, reads
#'   where the tag is absent (i.e. group is `NA`) are removed before
#'   downsampling. Default `FALSE` preserves the existing behaviour of keeping
#'   unphased reads in the data.
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
                             alt_base = NULL, max_reads = 200L,
                             per_group_downsample = FALSE,
                             min_mapq = 0L,
                             strand_filter = c("+", "-"),
                             min_read_length = 0L,
                             drop_na_group = FALSE) {
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
  strand_filter <- match.arg(strand_filter, choices = c("+", "-"), several.ok = TRUE)
  if (!is.integer(min_mapq))        min_mapq <- as.integer(min_mapq)
  if (!is.integer(min_read_length)) min_read_length <- as.integer(min_read_length)
  validate_bam_index(bam)
  parsed <- parse_region(region)

  # --- 2. Set up Rsamtools query ---
  gr <- GenomicRanges::GRanges(
    seqnames = parsed$chrom,
    ranges = IRanges::IRanges(start = parsed$start, end = parsed$end)
  )

  what <- c("qname", "pos", "cigar", "strand", "seq", "mapq")

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

  # bam_indices maps filtered-reads rows back to original bam_data positions
  bam_indices <- seq_len(nrow(reads))

  # Drop reads with NA group if requested
  if (isTRUE(drop_na_group) && !is.null(group_tag) && "group" %in% names(reads)) {
    na_mask <- is.na(reads$group)
    if (any(na_mask)) {
      reads       <- reads[!na_mask, , drop = FALSE]
      bam_indices <- bam_indices[!na_mask]
      rownames(reads) <- NULL
    }
  }

  # --- 4b. Apply read-level filters ---
  n_before <- nrow(reads)
  filter_mask <- rep(TRUE, n_before)

  # MAPQ filter
  if (min_mapq > 0L) {
    mapq_vec <- bam_data$mapq
    filter_mask <- filter_mask & (!is.na(mapq_vec) & mapq_vec >= min_mapq)
  }

  # Strand filter
  if (!identical(sort(strand_filter), c("+", "-"))) {
    filter_mask <- filter_mask & (reads$strand %in% strand_filter)
  }

  # Read length filter (ref_widths computed earlier from CIGAR)
  if (min_read_length > 0L) {
    filter_mask <- filter_mask & (ref_widths >= min_read_length)
  }

  if (!all(filter_mask)) {
    n_removed <- sum(!filter_mask)
    frac_removed <- n_removed / n_before
    if (frac_removed > 0.5) {
      warning(sprintf(
        "%.0f%% of reads (%d/%d) were removed by filters (min_mapq=%d, strand_filter=c(%s), min_read_length=%d).",
        frac_removed * 100, n_removed, n_before, min_mapq,
        paste(sprintf('"%s"', strand_filter), collapse = ", "),
        min_read_length
      ), call. = FALSE)
    }
    reads       <- reads[filter_mask, , drop = FALSE]
    bam_indices <- bam_indices[filter_mask]
    rownames(reads) <- NULL
  }

  if (nrow(reads) == 0L) {
    warning("No reads remain after applying filters.", call. = FALSE)
    return(empty_methylation_data(gr, mod_code, group_tag))
  }

  # --- 5. SNV-based grouping ---

  if (!is.null(snv_position)) {
    n_all <- nrow(reads)
    snv_bases <- character(n_all)
    for (j in seq_len(n_all)) {
      idx_j <- bam_indices[j]
      seq_str <- as.character(bam_data$seq[[idx_j]])
      q_pos <- ref_to_seq(bam_data$cigar[idx_j], bam_data$pos[idx_j], snv_position)
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
    bam_indices <- bam_indices[snv_keep]
    reads <- reads[snv_keep, , drop = FALSE]
    rownames(reads) <- NULL
    group_tag <- "SNV"
  }

  # --- 6. Downsample if needed ---
  n_reads <- nrow(reads)
  keep_local <- seq_len(n_reads)

  active_grouping <- per_group_downsample &&
                     !is.null(group_tag) && "group" %in% names(reads)

  if (active_grouping) {
    groups <- unique(reads$group[!is.na(reads$group)])
    keep_local <- unlist(lapply(groups, function(g) {
      idx <- which(reads$group == g)
      if (length(idx) > max_reads) idx <- sort(sample(idx, max_reads))
      idx
    }), use.names = FALSE)
    keep_local <- sort(keep_local)
    reads <- reads[keep_local, , drop = FALSE]
    rownames(reads) <- NULL
  } else if (n_reads > max_reads) {
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
  # Retain sequences and CIGARs for kept reads (no additional I/O needed)
  sequences <- setNames(as.character(bam_data$seq[keep]), reads$read_name)
  cigars    <- setNames(bam_data$cigar[keep], reads$read_name)

  # (mod_code column is already populated per-code in the inner loop above)
  structure(
    list(
      reads = reads,
      sites = sites,
      region = gr,
      mod_code = mod_code,
      group_tag = group_tag,
      snv_position = snv_position,
      sequences = sequences,
      cigars = cigars
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
      snv_position = snv_position,
      sequences = setNames(character(0), character(0)),
      cigars    = setNames(character(0), character(0))
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

  n_plus  <- sum(x$reads$strand == "+", na.rm = TRUE)
  n_minus <- sum(x$reads$strand == "-", na.rm = TRUE)

  read_lengths <- x$reads$end - x$reads$start + 1L
  med_len <- if (length(read_lengths) > 0L) as.integer(median(read_lengths)) else NA_integer_

  cat("methylation_data object\n")
  cat(sprintf("Region: %s:%d-%d\n", chrom, start, end))
  cat(sprintf("Reads: %d  (+ strand: %d,  - strand: %d)\n",
              nrow(x$reads), n_plus, n_minus))
  if (!is.na(med_len))
    cat(sprintf("Median read length: %d bp\n", med_len))
  cat(sprintf("Modification sites: %d\n", nrow(x$sites)))
  cat(sprintf("Modification code(s): %s\n", paste(x$mod_code, collapse = ", ")))

  if (!is.null(x$group_tag)) {
    n_groups <- length(unique(x$reads$group[!is.na(x$reads$group)]))
    cat(sprintf("Group tag: %s (%d groups)\n", x$group_tag, n_groups))
    groups <- sort(unique(x$reads$group[!is.na(x$reads$group)]))
    for (g in groups) {
      n_g   <- sum(x$reads$group == g, na.rm = TRUE)
      mods  <- x$sites$mod_prob[!is.na(x$sites$group) & x$sites$group == g]
      mean_g <- if (length(mods) > 0L) round(mean(mods, na.rm = TRUE), 2L) else NA_real_
      cat(sprintf("  %s: %d reads, mean methylation %.2f\n", g, n_g, mean_g))
    }
  }

  invisible(x)
}

#' Summary of a methylation_data object
#'
#' @param object A \code{methylation_data} object.
#' @param ... Unused.
#' @return A named list with summary statistics (invisibly).
#' @export
summary.methylation_data <- function(object, ...) {
  chrom <- as.character(GenomicRanges::seqnames(object$region))
  reg_str <- sprintf("%s:%d-%d", chrom,
                     GenomicRanges::start(object$region),
                     GenomicRanges::end(object$region))

  strand_df <- as.data.frame(table(strand = object$reads$strand),
                              stringsAsFactors = FALSE)
  names(strand_df)[2] <- "n_reads"

  lens <- object$reads$end - object$reads$start + 1L
  rl <- if (length(lens) > 0L)
    list(median = as.integer(median(lens)), min = min(lens), max = max(lens))
  else
    list(median = NA_integer_, min = NA_integer_, max = NA_integer_)

  overall_mean <- if (nrow(object$sites) > 0L)
    round(mean(object$sites$mod_prob, na.rm = TRUE), 4L)
  else NA_real_

  groups_df <- NULL
  if (!is.null(object$group_tag)) {
    g_levels <- sort(unique(object$reads$group[!is.na(object$reads$group)]))
    groups_df <- do.call(rbind, lapply(g_levels, function(g) {
      n_g  <- sum(object$reads$group == g, na.rm = TRUE)
      mods <- object$sites$mod_prob[
        !is.na(object$sites$group) & object$sites$group == g]
      data.frame(
        group           = g,
        n_reads         = n_g,
        mean_mod_prob   = if (length(mods) > 0L) round(mean(mods,   na.rm = TRUE), 4L) else NA_real_,
        median_mod_prob = if (length(mods) > 0L) round(median(mods, na.rm = TRUE), 4L) else NA_real_,
        stringsAsFactors = FALSE
      )
    }))
  }

  out <- list(
    region                = reg_str,
    n_reads               = nrow(object$reads),
    n_sites               = nrow(object$sites),
    mod_code              = object$mod_code,
    group_tag             = object$group_tag,
    strand                = strand_df,
    read_length           = rl,
    groups                = groups_df,
    overall_mean_mod_prob = overall_mean
  )

  cat("methylation_data summary\n========================\n")
  cat(sprintf("Region:            %s\n", reg_str))
  cat(sprintf("Reads:             %d\n", out$n_reads))
  for (i in seq_len(nrow(strand_df)))
    cat(sprintf("  %s strand:        %d\n", strand_df$strand[i], strand_df$n_reads[i]))
  cat(sprintf("  Read length:     min=%d, median=%d, max=%d bp\n",
              rl$min, rl$median, rl$max))
  cat(sprintf("Sites:             %d\n", out$n_sites))
  cat(sprintf("Modification:      %s\n", paste(out$mod_code, collapse = ", ")))
  cat(sprintf("Overall mean mod:  %.2f\n", overall_mean))
  if (!is.null(groups_df)) {
    cat(sprintf("\nGroup breakdown (%s):\n", object$group_tag))
    print(groups_df, row.names = FALSE, digits = 4)
  }

  invisible(out)
}
