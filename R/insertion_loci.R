# Insertion locus clustering and accessor for ggmethylation

#' Access inserted-base modification calls from a methylation_data object
#'
#' Convenience accessor for `methylation_data$insertion_sites`.
#'
#' @param m A `methylation_data` object produced by [read_methylation()].
#'
#' @return A data.frame with columns `read_name`, `ref_anchor`, `query_pos`,
#'   `ins_offset`, `ins_length`, `mod_prob`, `mod_code`, and optionally
#'   `group`. Zero rows when no inserted-base modifications were found.
#'
#' @export
insertion_sites <- function(m) {
  if (!inherits(m, "methylation_data")) {
    stop("'m' must be a methylation_data object.", call. = FALSE)
  }
  m$insertion_sites
}

#' List insertion loci from a methylation_data object
#'
#' Clusters insertion events across reads into loci using a greedy single-pass
#' algorithm. Insertions are grouped when their reference anchor positions and
#' lengths are within the specified tolerances of each other. Returns a summary
#' tibble sorted by anchor position; pass a `locus_id` value to
#' [plot_insertion_locus()] to visualise a specific locus.
#'
#' @param m A `methylation_data` object produced by [read_methylation()].
#' @param tol_pos Integer. Maximum distance in base pairs between consecutive
#'   insertion anchor positions to be merged into the same locus (default 10).
#' @param tol_len Numeric. Maximum fractional difference in insertion length
#'   relative to the running cluster median to be merged into the same locus
#'   (default 0.20, i.e., 20 percent).
#' @param min_reads Integer. Minimum number of carrier reads required for a
#'   locus to be reported (default 2).
#'
#' @return A data.frame (tibble-compatible) with one row per locus and columns:
#'   \describe{
#'     \item{locus_id}{Character. Deterministic identifier
#'       `"INS_<chrom>_<anchor_pos>_<median_length>bp"`.}
#'     \item{chrom}{Character. Chromosome name.}
#'     \item{anchor_pos}{Integer. Median reference anchor position across
#'       carrier reads in the locus.}
#'     \item{n_carriers}{Integer. Number of reads carrying an insertion at this
#'       locus.}
#'     \item{n_noncarriers}{Integer. Number of reads that span the locus
#'       (within `tol_pos` of `anchor_pos`) but do not carry an insertion.}
#'     \item{median_length}{Integer. Median insertion length across carriers.}
#'     \item{length_min}{Integer. Minimum insertion length in the locus.}
#'     \item{length_max}{Integer. Maximum insertion length in the locus.}
#'     \item{mean_ins_mod_prob}{Numeric or NA. Mean modification probability
#'       across all insertion-site calls at this locus. NA when
#'       `m$insertion_sites` has no rows for this locus.}
#'   }
#'   Sorted by `anchor_pos`. Returns a zero-row data.frame when no insertions
#'   pass the `min_reads` threshold.
#'
#' @export
list_insertion_loci <- function(m, tol_pos = 10L, tol_len = 0.20,
                                min_reads = 2L) {
  if (!inherits(m, "methylation_data")) {
    stop("'m' must be a methylation_data object.", call. = FALSE)
  }
  tol_pos   <- as.integer(tol_pos)
  min_reads <- as.integer(min_reads)

  cf <- m$cigar_features
  ins <- cf[cf$type == "I" & !is.na(cf$ref_start), , drop = FALSE]

  empty_out <- data.frame(
    locus_id         = character(0L),
    chrom            = character(0L),
    anchor_pos       = integer(0L),
    n_carriers       = integer(0L),
    n_noncarriers    = integer(0L),
    median_length    = integer(0L),
    length_min       = integer(0L),
    length_max       = integer(0L),
    mean_ins_mod_prob = numeric(0L),
    stringsAsFactors = FALSE
  )

  if (nrow(ins) == 0L) return(empty_out)

  chrom <- as.character(GenomicRanges::seqnames(m$region))

  # Sort by ref_start for single-pass greedy clustering
  ord <- order(ins$ref_start)
  ins <- ins[ord, , drop = FALSE]

  cluster_id  <- integer(nrow(ins))
  current_id  <- 1L
  cluster_id[1L] <- current_id
  prev_anchor  <- ins$ref_start[1L]
  cluster_lens <- ins$length[1L]

  for (i in seq_len(nrow(ins))[-1L]) {
    cur_pos <- ins$ref_start[i]
    cur_len <- ins$length[i]
    med_len <- median(cluster_lens)

    frac_diff <- if (med_len > 0) abs(cur_len - med_len) / med_len else 0

    if ((cur_pos - prev_anchor) > tol_pos || frac_diff > tol_len) {
      current_id  <- current_id + 1L
      cluster_lens <- cur_len
    } else {
      cluster_lens <- c(cluster_lens, cur_len)
    }
    cluster_id[i] <- current_id
    prev_anchor   <- cur_pos
  }

  n_clusters <- max(cluster_id)
  result_list <- vector("list", n_clusters)

  all_read_names <- m$reads$read_name
  ins_sites <- m$insertion_sites

  for (cid in seq_len(n_clusters)) {
    rows <- ins[cluster_id == cid, , drop = FALSE]

    if (nrow(rows) < min_reads) next

    med_len <- as.integer(median(rows$length))

    # Post-filter to match plot_insertion_locus() carrier definition: within
    # tol_len of the final cluster median, not the running merge median.
    carrier_rows  <- rows[abs(rows$length - med_len) / pmax(med_len, 1L) <= tol_len,
                          , drop = FALSE]
    carrier_names <- unique(carrier_rows$read_name)
    n_carriers    <- length(carrier_names)
    if (n_carriers < min_reads) next

    anchor  <- as.integer(median(carrier_rows$ref_start))
    med_len <- as.integer(median(carrier_rows$length))

    # Non-carriers: reads spanning anchor_pos +/- tol_pos, not in carrier set
    read_start <- pmin(m$reads$bam_pos, m$reads$start)
    read_end   <- m$reads$end
    spans_locus <- read_start <= (anchor + tol_pos) & read_end >= (anchor - tol_pos)
    noncarrier_mask <- spans_locus & !(all_read_names %in% carrier_names)
    n_noncarriers <- sum(noncarrier_mask)

    # Mean mod_prob from insertion_sites for carrier reads at this locus
    mean_mod <- NA_real_
    if (!is.null(ins_sites) && nrow(ins_sites) > 0L) {
      ins_sub <- ins_sites[
        ins_sites$read_name %in% carrier_names &
        ins_sites$ref_anchor >= (anchor - tol_pos) &
        ins_sites$ref_anchor <= (anchor + tol_pos), , drop = FALSE
      ]
      if (nrow(ins_sub) > 0L) {
        mean_mod <- mean(ins_sub$mod_prob, na.rm = TRUE)
      }
    }

    locus_id <- sprintf("INS_%s_%d_%dbp", chrom, anchor, med_len)

    result_list[[cid]] <- data.frame(
      locus_id          = locus_id,
      chrom             = chrom,
      anchor_pos        = anchor,
      n_carriers        = n_carriers,
      n_noncarriers     = n_noncarriers,
      median_length     = med_len,
      length_min        = min(carrier_rows$length),
      length_max        = max(carrier_rows$length),
      mean_ins_mod_prob = mean_mod,
      stringsAsFactors  = FALSE
    )
  }

  result_list <- Filter(Negate(is.null), result_list)
  if (length(result_list) == 0L) return(empty_out)

  out <- do.call(rbind, result_list)
  rownames(out) <- NULL
  out[order(out$anchor_pos), , drop = FALSE]
}
