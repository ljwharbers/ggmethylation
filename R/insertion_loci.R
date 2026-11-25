# Insertion locus clustering and accessor for ggmethylation

#' Access inserted-base modification calls from a methylation_data object
#'
#' Convenience accessor for `methylation_data$insertion_sites`.
#'
#' @param data A `methylation_data` object produced by [read_methylation()].
#'
#' @return A data.frame with columns `read_name`, `ref_anchor`, `query_pos`,
#'   `ins_offset`, `ins_length`, `mod_prob`, `mod_code`, and optionally
#'   `group`. Zero rows when no inserted-base modifications were found.
#'
#' @export
insertion_sites = function(data) {
  if (!inherits(data, "methylation_data")) {
    stop("'data' must be a methylation_data object.", call. = FALSE)
  }
  data$insertion_sites
}

#' List insertion loci from a methylation_data object
#'
#' Clusters insertion events (CIGAR `I` operations) across reads into loci.
#' Insertions are first grouped by position: sorted by reference anchor, a new
#' group starts wherever two consecutive anchors are more than `tol_pos` apart.
#' Within each position group, insertions are then sorted by length and split
#' wherever an insertion differs from the running median length of its
#' sub-group by more than `tol_len`; each length sub-group is finally split by
#' position again, so that two same-length insertions joined only by a chain
#' of other-length insertions are not merged. An insertion of a different
#' length at the same place therefore never splits a locus. Returns a summary
#' data.frame sorted by anchor position; pass one of its rows to
#' [plot_insertion_locus()] to visualise that locus.
#'
#' A read is a **carrier** when it has an insertion in the locus whose length is
#' within `tol_len` of the locus median length. A read is a **non-carrier** when
#' its alignment spans the locus on both sides (it starts at or before
#' `anchor_pos - tol_pos` and ends at or after `anchor_pos + tol_pos`) and it
#' has no insertion of at least `min_length` within `tol_pos` of the anchor.
#' Reads with an insertion of a different length at the locus (another allele)
#' are neither carriers nor non-carriers, and neither are reads that end or
#' start within `tol_pos` of the anchor (often soft-clipped at the insertion
#' breakpoint).
#'
#' @param data A `methylation_data` object produced by [read_methylation()].
#' @param tol_pos Non-negative integer. Maximum distance in base pairs between
#'   consecutive insertion anchor positions to be merged into the same locus,
#'   and the margin a non-carrier must extend past the anchor on each side
#'   (default 10).
#' @param tol_len Non-negative number. Maximum fractional difference in
#'   insertion length relative to the running cluster median to be merged into
#'   the same locus (default 0.20, i.e., 20 percent).
#' @param min_reads Positive integer. Minimum number of carrier reads required
#'   for a locus to be reported (default 2).
#' @param min_length Positive integer. Insertions shorter than this many bases
#'   are ignored entirely: they form no locus, and a read carrying only such
#'   insertions near a locus still counts as a non-carrier. The default, 1,
#'   keeps every insertion; raise it (e.g. to 50 for structural-variant-sized
#'   insertions) to ignore small sequencing-error insertions.
#'
#' @return A data.frame with one row per locus and columns:
#'   \describe{
#'     \item{locus_id}{Character. Deterministic identifier
#'       `"INS_<chrom>_<anchor_pos>_<median_length>bp"`.}
#'     \item{chrom}{Character. Chromosome name.}
#'     \item{anchor_pos}{Integer. Median reference anchor position across
#'       carrier reads in the locus. Like `ref_start` of an `I` operation in
#'       `data$cigar_features`, this is the first reference base after the
#'       insertion: the inserted bases lie between `anchor_pos - 1` and
#'       `anchor_pos`.}
#'     \item{n_carriers}{Integer. Number of reads carrying an insertion at this
#'       locus.}
#'     \item{n_noncarriers}{Integer. Number of reads that span the locus
#'       (from `anchor_pos - tol_pos` to `anchor_pos + tol_pos`) and carry no
#'       insertion there.}
#'     \item{median_length}{Integer. Median insertion length across carriers.}
#'     \item{length_min}{Integer. Minimum insertion length in the locus.}
#'     \item{length_max}{Integer. Maximum insertion length in the locus.}
#'     \item{mean_ins_mod_prob}{Numeric or NA. Mean modification probability
#'       across the carriers' inserted-base calls at this locus (all
#'       modification codes pooled). NA when `data$insertion_sites` has no
#'       rows for this locus.}
#'     \item{carriers}{List-column. Per locus, a named integer vector of the
#'       carrier reads' insertion lengths, named by `read_name`.}
#'     \item{noncarriers}{List-column. Per locus, the `read_name`s of the
#'       non-carrier reads.}
#'     \item{tol_pos, tol_len}{The tolerances used, so [plot_insertion_locus()]
#'       selects insertion calls exactly as they were clustered.}
#'   }
#'   Sorted by `anchor_pos`. Returns a zero-row data.frame when no insertions
#'   pass the `min_reads` threshold.
#'
#' @examples
#' \dontrun{
#' md   = read_methylation("sample.bam", "chr21:34500000-34510000")
#' loci = list_insertion_loci(md, tol_pos = 10L, tol_len = 0.20, min_reads = 2L,
#'                            min_length = 50L)
#' loci$locus_id
#' plot_insertion_locus(md, loci[1, ])
#' }
#'
#' @export
list_insertion_loci = function(data, tol_pos = 10L, tol_len = 0.20,
                                min_reads = 2L, min_length = 1L) {
  if (!inherits(data, "methylation_data")) {
    stop("'data' must be a methylation_data object.", call. = FALSE)
  }
  tol_pos    = .ins_check_number(tol_pos, "tol_pos", min = 0, whole = TRUE)
  tol_len    = .ins_check_number(tol_len, "tol_len", min = 0)
  min_reads  = .ins_check_number(min_reads, "min_reads", min = 1, whole = TRUE)
  min_length = .ins_check_number(min_length, "min_length", min = 1, whole = TRUE)

  cf = data$cigar_features
  ins = cf[cf$type == "I" & !is.na(cf$ref_start) & cf$length >= min_length, ,
           drop = FALSE]

  empty_out = data.frame(
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
  empty_out$carriers    = list()
  empty_out$noncarriers = list()
  empty_out$tol_pos     = integer(0L)
  empty_out$tol_len     = numeric(0L)

  if (nrow(ins) == 0L) return(empty_out)

  chrom = as.character(GenomicRanges::seqnames(data$region))

  clusters = .cluster_insertions(ins$ref_start, ins$length, tol_pos, tol_len)
  result_list = vector("list", length(clusters))

  all_read_names = data$reads$read_name
  read_start = pmin(data$reads$bam_pos, data$reads$start)
  # `end` is clipped to the queried region; the unclipped end decides whether a
  # read spans a locus near the region edge (hand-built objects may lack it).
  read_end   = if (is.null(data$reads$bam_end)) {
    data$reads$end
  } else {
    pmax(data$reads$bam_end, data$reads$end)
  }
  ins_sites  = data$insertion_sites

  for (cid in seq_along(clusters)) {
    rows = ins[clusters[[cid]], , drop = FALSE]
    if (length(unique(rows$read_name)) < min_reads) next

    med_len = as.integer(stats::median(rows$length))

    # Post-filter to match plot_insertion_locus() carrier definition: within
    # tol_len of the final cluster median, not the running merge median.
    carrier_rows  = rows[abs(rows$length - med_len) / pmax(med_len, 1L) <= tol_len,
                          , drop = FALSE]
    carrier_names = unique(carrier_rows$read_name)
    n_carriers    = length(carrier_names)
    if (n_carriers < min_reads) next

    anchor  = as.integer(stats::median(carrier_rows$ref_start))
    med_len = as.integer(stats::median(carrier_rows$length))

    # Non-carriers: reads spanning anchor_pos +/- tol_pos on both sides, with
    # no insertion (of at least min_length) within tol_pos of the anchor.
    # Reads with an other-length insertion there belong to neither set.
    near = abs(ins$ref_start - anchor) <= tol_pos
    other_allele  = setdiff(unique(ins$read_name[near]), carrier_names)
    spans_locus   = read_start <= (anchor - tol_pos) & read_end >= (anchor + tol_pos)
    noncarrier_mask = spans_locus & !(all_read_names %in% c(carrier_names, other_allele))
    noncarriers   = all_read_names[noncarrier_mask]

    # Insertion length per carrier read (its first matching insertion)
    first = !duplicated(carrier_rows$read_name)
    carrier_len = stats::setNames(carrier_rows$length[first],
                                  carrier_rows$read_name[first])

    # Mean mod_prob over the carriers' own inserted bases at this locus
    mean_mod = NA_real_
    ins_sub = .carrier_insertion_sites(ins_sites, carrier_len, anchor, tol_pos)
    if (nrow(ins_sub) > 0L) {
      mean_mod = mean(ins_sub$mod_prob, na.rm = TRUE)
    }

    locus_id = sprintf("INS_%s_%d_%dbp", chrom, anchor, med_len)

    row = data.frame(
      locus_id          = locus_id,
      chrom             = chrom,
      anchor_pos        = anchor,
      n_carriers        = n_carriers,
      n_noncarriers     = length(noncarriers),
      median_length     = med_len,
      length_min        = min(carrier_rows$length),
      length_max        = max(carrier_rows$length),
      mean_ins_mod_prob = mean_mod,
      stringsAsFactors  = FALSE
    )
    row$carriers    = list(carrier_len)
    row$noncarriers = list(noncarriers)
    row$tol_pos     = tol_pos
    row$tol_len     = tol_len
    result_list[[cid]] = row
  }

  result_list = Filter(Negate(is.null), result_list)
  if (length(result_list) == 0L) return(empty_out)

  out = do.call(rbind, result_list)
  out = out[order(out$anchor_pos, out$median_length), , drop = FALSE]
  rownames(out) = NULL
  out
}

# Cluster insertions by position, then by length within each position group,
# then by position again within each length group. Returns a list of integer
# index vectors into `pos`/`len`.
.cluster_insertions = function(pos, len, tol_pos, tol_len) {
  out = list()
  for (p_idx in .split_by_gap(pos, seq_along(pos), tol_pos)) {
    for (l_idx in .split_by_length(len, p_idx, tol_len)) {
      out = c(out, .split_by_gap(pos, l_idx, tol_pos))
    }
  }
  out
}

# Split `idx` (indices into `pos`) wherever consecutive sorted positions are
# more than `tol` apart. Each group is ordered by position, ties by index.
.split_by_gap = function(pos, idx, tol) {
  idx = idx[order(pos[idx], idx)]
  new_group = c(TRUE, diff(pos[idx]) > tol)
  unname(split(idx, cumsum(new_group)))
}

# Split `idx` (indices into `len`) by length: sorted ascending, an insertion
# starts a new group when its fractional difference to the running median of
# the current group exceeds `tol`.
.split_by_length = function(len, idx, tol) {
  idx = idx[order(len[idx])]
  groups = list()
  current = idx[1L]
  for (i in idx[-1L]) {
    med = stats::median(len[current])
    if (abs(len[i] - med) / med > tol) {
      groups = c(groups, list(current))
      current = i
    } else {
      current = c(current, i)
    }
  }
  c(groups, list(current))
}

# Inserted-base calls belonging to a locus: the carrier reads' calls anchored
# within `tol_pos` of `anchor`, from the insertion whose length is the one
# recorded for that carrier (so a second, other insertion nearby on the same
# read is not mixed in).
.carrier_insertion_sites = function(ins_sites, carrier_len, anchor, tol_pos) {
  if (is.null(ins_sites) || nrow(ins_sites) == 0L) {
    return(ins_sites[integer(0L), , drop = FALSE])
  }
  keep = ins_sites$read_name %in% names(carrier_len) &
    abs(ins_sites$ref_anchor - anchor) <= tol_pos
  keep[keep] = ins_sites$ins_length[keep] ==
    carrier_len[ins_sites$read_name[keep]]
  ins_sites[keep, , drop = FALSE]
}

# Validate a length-1 finite number >= `min` (optionally a whole number) and
# return it (as integer when `whole`).
.ins_check_number = function(x, name, min = -Inf, whole = FALSE) {
  ok = is.numeric(x) && length(x) == 1L && !is.na(x) && is.finite(x) && x >= min
  if (ok && whole) ok = x == round(x)
  if (!ok) {
    what = if (whole) "a single whole number" else "a single number"
    stop(sprintf("`%s` must be %s >= %s.", name, what, format(min)), call. = FALSE)
  }
  if (whole) as.integer(x) else x
}

# Validate a single TRUE/FALSE.
.ins_check_flag = function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop(sprintf("`%s` must be TRUE or FALSE.", name), call. = FALSE)
  }
  x
}
