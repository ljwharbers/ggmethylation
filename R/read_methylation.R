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
#' @param snv `NULL`, or a list `list(position =, ref =, alt =)` to group reads
#'   by the base they carry at an SNV instead of by a BAM tag: reads with the
#'   `ref` base form group `"REF"`, reads with `alt` form `"ALT"`, and all
#'   other reads are dropped. Cannot be combined with `group_tag`.
#' @param max_reads Integer. Maximum number of reads to return (default 200).
#'   If more reads overlap the region, a random subset is kept.
#' @param per_group_downsample Logical. When `TRUE` and grouping is active
#'   (via `group_tag` or `snv`), the `max_reads` cap is applied
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
#'       `bam_pos` (original unclipped alignment start, equal to the BAM POS
#'       field), `strand`, `is_supplementary` (logical), `sa_chrom` (chr or
#'       NA), `sa_pos` (int or NA), `clip_side` (chr: `"left"`, `"right"`,
#'       `"both"`, or NA), and optionally `group`. `start` and `end` are
#'       clipped to the queried region; `bam_pos` retains the original start.}
#'     \item{sites}{Data.frame with columns `position`, `mod_prob`,
#'       `read_name`, `mod_code`, and optionally `group`.}
#'     \item{region}{A [GenomicRanges::GRanges] object for the queried region.}
#'     \item{mod_code}{The modification code(s) used.}
#'     \item{group_tag}{The grouping tag used, or `NULL`.}
#'     \item{snv_position}{The SNV position used for grouping, or `NULL`.}
#'     \item{sequences}{Named character vector of read sequences, keyed by
#'       `read_name`. Used internally by [plot_methylation()] for variant
#'       overlay.}
#'     \item{cigars}{Named character vector of CIGAR strings, keyed by
#'       `read_name`. Used internally for CIGAR feature decomposition.}
#'     \item{cigar_features}{Data.frame of structural CIGAR features (type,
#'       ref_start, ref_end, query_start, query_end, length, read_name) for
#'       insertions (`"I"`), deletions (`"D"`), and skipped regions (`"N"`).}
#'     \item{insertion_sites}{Data.frame of modification calls that fall on
#'       inserted bases (no reference coordinate). Columns: `read_name`,
#'       `ref_anchor` (reference position the insertion sits before),
#'       `query_pos` (1-based position in the read sequence), `ins_offset`
#'       (1-based position within the insertion), `ins_length`, `mod_prob`,
#'       `mod_code`, and optionally `group`. Zero rows when no insertions
#'       carry modification calls.}
#'   }
#'
#' @examples
#' \dontrun{
#' md = read_methylation("sample.bam", "chr1:1000-2000")
#' md = read_methylation("sample.bam", "chr1:1000-2000",
#'   group_tag = "HP", max_reads = 100
#' )
#' md = read_methylation("sample.bam", "chr1:1000-2000",
#'   snv = list(position = 1500, ref = "C", alt = "T")
#' )
#' }
#'
#' @export
read_methylation = function(bam, region, mod_code = "m", group_tag = NULL,
                             snv = NULL, max_reads = 200L,
                             per_group_downsample = FALSE,
                             min_mapq = 0L,
                             strand_filter = c("+", "-"),
                             min_read_length = 0L,
                             drop_na_group = FALSE) {
  # --- 1. Validate inputs ---
  mod_code = as.character(mod_code)          # coerce in case of factor
  if (length(mod_code) == 0L || any(!nzchar(mod_code))) {
    stop("'mod_code' must be a non-empty character vector of modification codes.",
         call. = FALSE)
  }
  mod_code = unique(mod_code)                # drop accidental duplicates

  if (!is.null(group_tag) && !is.null(snv)) {
    stop("'group_tag' and 'snv' are mutually exclusive.", call. = FALSE)
  }
  snv = .check_snv(snv)
  strand_filter = match.arg(strand_filter, choices = c("+", "-"), several.ok = TRUE)
  min_mapq = as.integer(min_mapq)
  min_read_length = as.integer(min_read_length)
  validate_bam_index(bam)

  # --- 2. Query the BAM ---
  gr = region_to_granges(region)
  param = Rsamtools::ScanBamParam(
    which = gr,
    what  = c("qname", "flag", "pos", "cigar", "strand", "seq", "mapq"),
    tag   = c("MM", "ML", "SA", group_tag)
  )
  bam_data = Rsamtools::scanBam(bam, param = param)[[1]]
  if (length(bam_data$qname) == 0L) {
    warning("No reads found in region '", region, "'.", call. = FALSE)
    return(empty_methylation_data(gr, mod_code, group_tag))
  }

  # `reads$.idx` maps each row back to its record in `bam_data`, so reads can
  # be filtered freely; it is dropped before returning.
  reads = .bam_reads(bam_data)

  # --- 3. Grouping by BAM tag ---
  if (!is.null(group_tag)) {
    group_values = bam_data$tag[[group_tag]]
    if (is.null(group_values) || all(is.na(group_values))) {
      warning(
        "Group tag '", group_tag, "' not found in any read. ",
        "Ignoring grouping.", call. = FALSE
      )
      reads$group = NA_character_
      group_tag = NULL
    } else {
      reads$group = as.character(group_values)
      if (isTRUE(drop_na_group)) reads = reads[!is.na(reads$group), , drop = FALSE]
    }
  }

  # --- 4. Read-level filters ---
  reads = .filter_reads(reads, bam_data, min_mapq, strand_filter, min_read_length)
  if (nrow(reads) == 0L) {
    warning("No reads remain after applying filters.", call. = FALSE)
    return(empty_methylation_data(gr, mod_code, group_tag))
  }

  # --- 5. Grouping by SNV genotype ---
  if (!is.null(snv)) {
    reads$group = .snv_genotype(bam_data, reads$.idx, snv$position, snv$ref, snv$alt)
    reads = reads[!is.na(reads$group), , drop = FALSE]
    if (nrow(reads) == 0L) {
      warning("No reads carry REF or ALT at SNV position ", snv$position, ".",
              call. = FALSE)
      return(empty_methylation_data(gr, mod_code, "SNV", snv$position))
    }
    group_tag = "SNV"
  }

  # --- 6. Downsample ---
  reads = .downsample_reads(reads, max_reads,
                            per_group = per_group_downsample && !is.null(group_tag))
  rownames(reads) = NULL

  # --- 7. Parse MM/ML tags and CIGARs ---
  parsed = .parse_read_mods(bam_data, reads, mod_code)
  sites = parsed$sites
  insertion_sites = parsed$insertion_sites
  cigar_features = parsed$cigar_features
  if (!is.null(group_tag)) {
    sites$group = reads$group[match(sites$read_name, reads$read_name)]
    insertion_sites$group = reads$group[match(insertion_sites$read_name, reads$read_name)]
  }

  # --- 8. Restrict everything to the region ---
  region_start = GenomicRanges::start(gr)
  region_end   = GenomicRanges::end(gr)
  reads$start = pmax(reads$start, region_start)
  reads$end   = pmin(reads$end, region_end)

  in_region = function(df, keep) {
    df = df[keep, , drop = FALSE]
    rownames(df) = NULL
    df
  }
  sites = in_region(sites, sites$position >= region_start & sites$position <= region_end)
  insertion_sites = in_region(
    insertion_sites,
    insertion_sites$ref_anchor >= region_start & insertion_sites$ref_anchor <= region_end
  )
  # Features overlapping the region; insertions (no ref_end) are a point at
  # ref_start.
  feature_end = ifelse(is.na(cigar_features$ref_end), cigar_features$ref_start,
                       cigar_features$ref_end)
  cigar_features = in_region(
    cigar_features,
    is.na(cigar_features$ref_start) |
      (cigar_features$ref_start <= region_end & feature_end >= region_start)
  )

  # --- 9. Return methylation_data object ---
  # Sequences and CIGARs are kept for the variant overlay.
  sequences = stats::setNames(as.character(bam_data$seq[reads$.idx]), reads$read_name)
  cigars    = stats::setNames(bam_data$cigar[reads$.idx], reads$read_name)
  reads$.idx = NULL

  .new_methylation_data(
    reads = reads, sites = sites, insertion_sites = insertion_sites, region = gr,
    mod_code = mod_code, group_tag = group_tag, snv_position = snv$position,
    sequences = sequences, cigars = cigars, cigar_features = cigar_features
  )
}

# Validate the `snv` argument: NULL, or list(position =, ref =, alt =) with a
# single position and single-base alleles. Returns it with an integer position.
.check_snv = function(snv) {
  if (is.null(snv)) return(NULL)
  snv = as.list(snv)
  ok = all(c("position", "ref", "alt") %in% names(snv)) &&
    all(lengths(snv[c("position", "ref", "alt")]) == 1L) &&
    !is.na(suppressWarnings(as.integer(snv$position))) &&
    all(nchar(c(snv$ref, snv$alt)) == 1L)
  if (!ok) {
    stop("`snv` must be list(position = <integer>, ref = <base>, alt = <base>), ",
         "e.g. list(position = 1500, ref = \"C\", alt = \"T\").", call. = FALSE)
  }
  list(position = as.integer(snv$position), ref = snv$ref, alt = snv$alt)
}

# Build the per-read table from a scanBam() result: coordinates, strand,
# supplementary flag, first SA partner and clip side, plus `.idx` (row index
# into `bam_data`).
.bam_reads = function(bam_data) {
  reads = data.frame(
    read_name = bam_data$qname,
    start     = bam_data$pos,
    end       = bam_data$pos + cigar_ref_width(bam_data$cigar) - 1L,
    bam_pos   = bam_data$pos,
    strand    = as.character(bam_data$strand),
    stringsAsFactors = FALSE
  )
  reads$is_supplementary = bitwAnd(bam_data$flag, 0x800L) > 0L

  sa_tags = bam_data$tag[["SA"]]
  if (!is.null(sa_tags)) {
    sa_parsed = lapply(sa_tags, parse_sa_tag)
    reads$sa_chrom = vapply(sa_parsed, function(x) {
      if (nrow(x) == 0L) NA_character_ else x$rname[1L]
    }, character(1L))
    reads$sa_pos = vapply(sa_parsed, function(x) {
      if (nrow(x) == 0L) NA_integer_ else x$pos[1L]
    }, integer(1L))
  } else {
    reads$sa_chrom = NA_character_
    reads$sa_pos   = NA_integer_
  }
  reads$clip_side = detect_clip_side(bam_data$cigar)

  # When a primary and its supplementary alignment(s) both fall in the region
  # the qname repeats: the primary keeps its name, supplementary copies become
  # "<qname>_supp1", "<qname>_supp2", ...
  renamed = reads$is_supplementary &
    reads$read_name %in% reads$read_name[duplicated(reads$read_name)]
  if (any(renamed)) {
    nth = stats::ave(seq_along(reads$read_name)[renamed], reads$read_name[renamed],
                     FUN = seq_along)
    reads$read_name[renamed] = paste0(reads$read_name[renamed], "_supp", nth)
  }

  reads$.idx = seq_len(nrow(reads))
  reads
}

# Apply the MAPQ, strand and reference-length filters, warning when they
# remove more than half of the reads.
.filter_reads = function(reads, bam_data, min_mapq, strand_filter, min_read_length) {
  keep = rep(TRUE, nrow(reads))
  if (min_mapq > 0L) {
    mapq = bam_data$mapq[reads$.idx]
    keep = keep & !is.na(mapq) & mapq >= min_mapq
  }
  if (!setequal(strand_filter, c("+", "-"))) {
    keep = keep & reads$strand %in% strand_filter
  }
  if (min_read_length > 0L) {
    keep = keep & cigar_ref_width(bam_data$cigar[reads$.idx]) >= min_read_length
  }
  if (all(keep)) return(reads)

  n_removed = sum(!keep)
  if (n_removed / length(keep) > 0.5) {
    warning(sprintf(
      "%.0f%% of reads (%d/%d) were removed by filters (min_mapq=%d, strand_filter=c(%s), min_read_length=%d).",
      100 * n_removed / length(keep), n_removed, length(keep), min_mapq,
      paste(sprintf('"%s"', strand_filter), collapse = ", "),
      min_read_length
    ), call. = FALSE)
  }
  reads[keep, , drop = FALSE]
}

# "REF" / "ALT" / NA per read, from the base each read carries at `position`.
.snv_genotype = function(bam_data, idx, position, ref_base, alt_base) {
  bases = vapply(idx, function(i) {
    seq_str = as.character(bam_data$seq[[i]])
    q_pos = ref_to_seq(bam_data$cigar[i], bam_data$pos[i], position)
    if (!is.na(q_pos) && q_pos >= 1L && q_pos <= nchar(seq_str)) {
      toupper(substr(seq_str, q_pos, q_pos))
    } else {
      NA_character_
    }
  }, character(1L))
  ifelse(bases == toupper(ref_base), "REF",
         ifelse(bases == toupper(alt_base), "ALT", NA_character_))
}

# Keep at most `max_reads` reads, sampled at random (in original order). With
# `per_group`, the cap applies to each non-NA group separately and reads
# without a group are dropped.
.downsample_reads = function(reads, max_reads, per_group = FALSE) {
  if (per_group) {
    groups = unique(reads$group[!is.na(reads$group)])
    keep = unlist(lapply(groups, function(g) {
      idx = which(reads$group == g)
      if (length(idx) > max_reads) idx = sort(sample(idx, max_reads))
      idx
    }), use.names = FALSE)
    return(reads[sort(keep), , drop = FALSE])
  }
  if (nrow(reads) > max_reads) {
    return(reads[sort(sample(nrow(reads), max_reads)), , drop = FALSE])
  }
  reads
}

# Parse every read's MM/ML tags (per modification code) and CIGAR.
# Returns list(sites, insertion_sites, cigar_features), not yet restricted to
# the region.
.parse_read_mods = function(bam_data, reads, mod_code) {
  n = nrow(reads)
  sites_list = ins_list = cigar_list = vector("list", n)

  for (j in seq_len(n)) {
    i         = reads$.idx[j]
    read_name = reads$read_name[j]
    seq_str   = as.character(bam_data$seq[[i]])
    dc        = decompose_cigar(bam_data$cigar[i], bam_data$pos[i])

    per_code = lapply(mod_code, function(code) {
      parsed = parse_mm_ml(
        seq      = seq_str,
        mm_tag   = bam_data$tag$MM[i],
        ml_tag   = bam_data$tag$ML[[i]],
        mod_code = code,
        strand   = reads$strand[j],
        cigar    = bam_data$cigar[i],
        pos      = bam_data$pos[i]
      )
      s = parsed$sites
      s$read_name = rep(read_name, nrow(s))
      s$mod_code  = rep(code, nrow(s))
      list(sites = s,
           ins = .anchor_insertion_sites(parsed$insertion_sites, dc, read_name, code))
    })
    sites_list[[j]] = do.call(rbind, lapply(per_code, `[[`, "sites"))
    ins_list[[j]]   = do.call(rbind, lapply(per_code, `[[`, "ins"))

    dc = dc[dc$type %in% c("I", "D", "N"), , drop = FALSE]
    dc$read_name = rep(read_name, nrow(dc))
    cigar_list[[j]] = dc
  }

  bind = function(pieces, empty) {
    out = do.call(rbind, c(list(empty), pieces))
    rownames(out) = NULL
    out
  }
  list(
    sites           = bind(sites_list, .empty_sites()),
    insertion_sites = bind(ins_list, .empty_insertion_sites()),
    cigar_features  = bind(cigar_list, .empty_cigar_features())
  )
}

# Attach each inserted-base call to the CIGAR `I` operation it sits in:
# the reference anchor, insertion length and 1-based offset within it.
.anchor_insertion_sites = function(ins, dc, read_name, code) {
  i_rows = dc[dc$type == "I", , drop = FALSE]
  if (nrow(ins) == 0L || nrow(i_rows) == 0L) return(.empty_insertion_sites())
  op = vapply(ins$query_pos, function(qp) {
    m = which(qp >= i_rows$query_start & qp <= i_rows$query_end)
    if (length(m) == 0L) NA_integer_ else m[1L]
  }, integer(1L))
  data.frame(
    read_name  = read_name,
    ref_anchor = i_rows$ref_start[op],
    query_pos  = ins$query_pos,
    ins_offset = ins$query_pos - i_rows$query_start[op] + 1L,
    ins_length = i_rows$length[op],
    mod_prob   = ins$mod_prob,
    mod_code   = code,
    stringsAsFactors = FALSE
  )
}

# Zero-row templates for the methylation_data tables.
.empty_reads = function() {
  data.frame(
    read_name        = character(0L),
    start            = integer(0L),
    end              = integer(0L),
    bam_pos          = integer(0L),
    strand           = character(0L),
    is_supplementary = logical(0L),
    sa_chrom         = character(0L),
    sa_pos           = integer(0L),
    clip_side        = character(0L),
    stringsAsFactors = FALSE
  )
}

.empty_sites = function() {
  data.frame(
    position  = integer(0L),
    mod_prob  = numeric(0L),
    read_name = character(0L),
    mod_code  = character(0L),
    stringsAsFactors = FALSE
  )
}

.empty_insertion_sites = function() {
  data.frame(
    read_name  = character(0L),
    ref_anchor = integer(0L),
    query_pos  = integer(0L),
    ins_offset = integer(0L),
    ins_length = integer(0L),
    mod_prob   = numeric(0L),
    mod_code   = character(0L),
    stringsAsFactors = FALSE
  )
}

.empty_cigar_features = function() {
  data.frame(
    type        = character(0L),
    ref_start   = integer(0L),
    ref_end     = integer(0L),
    query_start = integer(0L),
    query_end   = integer(0L),
    length      = integer(0L),
    read_name   = character(0L),
    stringsAsFactors = FALSE
  )
}

.new_methylation_data = function(reads, sites, insertion_sites, region, mod_code,
                                 group_tag, snv_position, sequences, cigars,
                                 cigar_features) {
  structure(
    list(
      reads           = reads,
      sites           = sites,
      insertion_sites = insertion_sites,
      region          = region,
      mod_code        = mod_code,
      group_tag       = group_tag,
      snv_position    = snv_position,
      sequences       = sequences,
      cigars          = cigars,
      cigar_features  = cigar_features
    ),
    class = "methylation_data"
  )
}

#' Create an empty methylation_data object
#'
#' @param gr A [GenomicRanges::GRanges] object.
#' @param mod_code Character. Modification code.
#' @param group_tag Character or NULL. Grouping tag.
#' @param snv_position Integer or NULL. SNV position used for grouping.
#'
#' @return An empty `methylation_data` object.
#'
#' @keywords internal
empty_methylation_data = function(gr, mod_code, group_tag, snv_position = NULL) {
  reads = .empty_reads()
  sites = .empty_sites()
  insertion_sites = .empty_insertion_sites()
  if (!is.null(group_tag)) {
    reads$group           = character(0L)
    sites$group           = character(0L)
    insertion_sites$group = character(0L)
  }
  no_reads = stats::setNames(character(0), character(0))
  .new_methylation_data(
    reads = reads, sites = sites, insertion_sites = insertion_sites, region = gr,
    mod_code = mod_code, group_tag = group_tag, snv_position = snv_position,
    sequences = no_reads, cigars = no_reads, cigar_features = .empty_cigar_features()
  )
}

# Per-group read counts and mean/median modification probability, one row per
# non-NA group in sorted order (NULL when there are no groups).
.group_stats = function(x) {
  groups = sort(unique(x$reads$group[!is.na(x$reads$group)]))
  do.call(rbind, lapply(groups, function(g) {
    mods = x$sites$mod_prob[!is.na(x$sites$group) & x$sites$group == g]
    data.frame(
      group           = g,
      n_reads         = sum(x$reads$group == g, na.rm = TRUE),
      mean_mod_prob   = if (length(mods) > 0L) mean(mods, na.rm = TRUE) else NA_real_,
      median_mod_prob = if (length(mods) > 0L) stats::median(mods, na.rm = TRUE) else NA_real_,
      stringsAsFactors = FALSE
    )
  }))
}

#' Print a methylation_data object
#'
#' @param x A `methylation_data` object.
#' @param ... Additional arguments (ignored).
#'
#' @return `x`, invisibly.
#'
#' @export
print.methylation_data = function(x, ...) {
  chrom = as.character(GenomicRanges::seqnames(x$region))
  start = GenomicRanges::start(x$region)
  end = GenomicRanges::end(x$region)

  n_plus  = sum(x$reads$strand == "+", na.rm = TRUE)
  n_minus = sum(x$reads$strand == "-", na.rm = TRUE)

  read_lengths = x$reads$end - x$reads$start + 1L
  med_len = if (length(read_lengths) > 0L) as.integer(stats::median(read_lengths)) else NA_integer_

  cat("methylation_data object\n")
  cat(sprintf("Region: %s:%d-%d\n", chrom, start, end))
  cat(sprintf("Reads: %d  (+ strand: %d,  - strand: %d)\n",
              nrow(x$reads), n_plus, n_minus))
  if (!is.na(med_len))
    cat(sprintf("Median read length: %d bp\n", med_len))
  cat(sprintf("Modification sites: %d\n", nrow(x$sites)))
  cat(sprintf("Modification code(s): %s\n", paste(x$mod_code, collapse = ", ")))

  if (!is.null(x$group_tag)) {
    stats_df = .group_stats(x)
    n_groups = if (is.null(stats_df)) 0L else nrow(stats_df)
    cat(sprintf("Group tag: %s (%d groups)\n", x$group_tag, n_groups))
    for (i in seq_len(n_groups)) {
      cat(sprintf("  %s: %d reads, mean methylation %.2f\n", stats_df$group[i],
                  stats_df$n_reads[i], round(stats_df$mean_mod_prob[i], 2L)))
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
summary.methylation_data = function(object, ...) {
  chrom = as.character(GenomicRanges::seqnames(object$region))
  reg_str = sprintf("%s:%d-%d", chrom,
                     GenomicRanges::start(object$region),
                     GenomicRanges::end(object$region))

  strand_df = as.data.frame(table(strand = object$reads$strand),
                              stringsAsFactors = FALSE)
  names(strand_df)[2] = "n_reads"

  lens = object$reads$end - object$reads$start + 1L
  rl = if (length(lens) > 0L)
    list(median = as.integer(stats::median(lens)), min = min(lens), max = max(lens))
  else
    list(median = NA_integer_, min = NA_integer_, max = NA_integer_)

  overall_mean = if (nrow(object$sites) > 0L)
    round(mean(object$sites$mod_prob, na.rm = TRUE), 4L)
  else NA_real_

  groups_df = NULL
  if (!is.null(object$group_tag)) {
    groups_df = .group_stats(object)
    groups_df$mean_mod_prob   = round(groups_df$mean_mod_prob, 4L)
    groups_df$median_mod_prob = round(groups_df$median_mod_prob, 4L)
  }

  out = list(
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
