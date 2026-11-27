#' Write methylation data to disk
#'
#' Saves the `$reads` and/or `$sites` data frames from a `methylation_data`
#' object to tab-separated (TSV) or BED files. Files can optionally be
#' gzip-compressed. A `multi_methylation_data` object gets one set of files
#' per sample.
#'
#' @param data A `methylation_data` object as returned by [read_methylation()],
#'   or a `multi_methylation_data` object as returned by
#'   [merge_methylation()]. For the latter, each sample is written with the
#'   prefix `{prefix}_{sample}`; sample names must not contain a path
#'   separator (`/` or a backslash).
#' @param prefix Character. Path prefix for output files. The appropriate
#'   suffix (`_reads.tsv`, `_sites.tsv`, or `_sites.bed`) is appended
#'   automatically. Intermediate directories are created if they do not exist.
#' @param format Character. Output format for the sites file: `"tsv"` (default)
#'   or `"bed"` (6-column BED, 0-based coordinates, no header).
#' @param reads Logical. Write the `$reads` data frame? Default `TRUE`.
#' @param sites Logical. Write the `$sites` data frame? Default `TRUE`.
#' @param gzip Logical. Compress output files with gzip? Default `FALSE`.
#'   When `TRUE`, `.gz` is appended to each file name.
#'
#' @return `data`, invisibly. Allows use in a pipeline.
#'
#' @details
#' The file names below are for a `methylation_data` object. For a
#' `multi_methylation_data` object, `{prefix}` becomes `{prefix}_{sample}`
#' for each sample, e.g. `out/run_A_reads.tsv` and `out/run_B_reads.tsv` for
#' `prefix = "out/run"` and samples `A` and `B`.
#'
#' **Reads file** (`{prefix}_reads.tsv[.gz]`)
#'
#' Always written as TSV with a header line. Columns: `read_name`, `start`,
#' `end`, `strand`, `group` (only when the object carries a group tag),
#' `mean_mod_prob`. `mean_mod_prob` is the per-read mean of `mod_prob` over
#' all of the read's rows in `$sites`, pooling every modification code in the
#' object: with `mod_code = c("m", "h")` it averages 5mC and 5hmC
#' probabilities together. For a per-code mean, keep one code in
#' `data$sites` before writing (e.g.
#' `data$sites = data$sites[data$sites$mod_code == "m", ]`). Reads with no
#' called sites receive `NA`.
#'
#' **Sites file — TSV** (`{prefix}_sites.tsv[.gz]`)
#'
#' Columns: `position`, `mod_prob`, `read_name`, `mod_code`, `group`
#' (if present). Coordinates are 1-based (as stored in the object).
#'
#' **Sites file — BED** (`{prefix}_sites.bed[.gz]`)
#'
#' 6-column BED (no header), sorted by chromosome and start. Coordinates are
#' converted to 0-based half-open intervals (`chromStart = position - 1`,
#' `chromEnd = position`). `name` is `<read_name>:<mod_code>` (split on the
#' last colon to recover both). `score` is the integer
#' `round(mod_prob * 1000)`, in `0`--`1000` as BED requires. `strand` is the
#' strand of the read's alignment (`.` if unknown).
#'
#' @examples
#' \dontrun{
#' md = read_methylation("sample.bam", "chr1:1000-2000")
#'
#' # Write both TSV files to the current directory
#' write_methylation(md, prefix = "output/sample_chr1")
#'
#' # BED format for sites, gzip-compressed
#' write_methylation(md, prefix = "output/sample_chr1",
#'                   format = "bed", gzip = TRUE)
#'
#' # Reads only
#' write_methylation(md, prefix = "output/sample_chr1",
#'                   sites = FALSE)
#'
#' # Several samples: output/cohort_A_reads.tsv, output/cohort_B_reads.tsv, ...
#' md2 = read_methylation("sample2.bam", "chr1:1000-2000")
#' write_methylation(merge_methylation(A = md, B = md2), prefix = "output/cohort")
#'
#' # Pipe-friendly usage
#' md |>
#'   write_methylation(prefix = "output/sample_chr1") |>
#'   plot_methylation()
#' }
#'
#' @seealso [read_methylation()]
#' @export
write_methylation = function(data,
                              prefix,
                              format = c("tsv", "bed"),
                              reads  = TRUE,
                              sites  = TRUE,
                              gzip   = FALSE) {

  # --- Input validation ---
  multi = inherits(data, "multi_methylation_data")
  if (!multi && !inherits(data, "methylation_data")) {
    stop("`data` must be a `methylation_data` or `multi_methylation_data` object.",
         call. = FALSE)
  }
  .check_string(prefix, "prefix")
  format = match.arg(format)
  .check_flag(reads, "reads")
  .check_flag(sites, "sites")
  .check_flag(gzip, "gzip")
  if (multi) {
    samples = names(data$samples)
    bad = samples[grepl("[/\\\\]", samples)]
    if (length(bad) > 0L) {
      stop(sprintf(
        "Sample names are used in file names and must not contain a path separator: %s.",
        paste(sprintf("'%s'", bad), collapse = ", ")
      ), call. = FALSE)
    }
  }
  if (!reads && !sites) {
    warning("Nothing to write: both `reads` and `sites` are FALSE.", call. = FALSE)
    return(invisible(data))
  }

  if (multi) {
    for (s in samples) {
      .write_methylation_one(data$samples[[s]], paste0(prefix, "_", s),
                             format, reads, sites, gzip)
    }
  } else {
    .write_methylation_one(data, prefix, format, reads, sites, gzip)
  }
  invisible(data)
}

# Internal helpers --------------------------------------------------------

# Write the files of one methylation_data object (arguments already checked).
.write_methylation_one = function(data, prefix, format, reads, sites, gzip) {
  # --- Create output directory if needed ---
  out_dir = dirname(prefix)
  if (!is.na(out_dir) && out_dir != "." && !dir.exists(out_dir))
    dir.create(out_dir, recursive = TRUE)

  # --- Compute mean_mod_prob per read ---
  data$reads$mean_mod_prob = .read_mean_mod_prob(data$sites, data$reads$read_name)

  # --- Write reads ---
  if (reads) {
    read_cols = c("read_name", "start", "end", "strand")
    if (!is.null(data$group_tag) && "group" %in% names(data$reads))
      read_cols = c(read_cols, "group")
    read_cols = c(read_cols, "mean_mod_prob")
    .write_table_maybe_gz(data$reads[, read_cols, drop = FALSE],
                          paste0(prefix, "_reads.tsv"), gzip)
  }

  # --- Write sites ---
  if (sites) {
    if (format == "tsv") {
      site_cols = c("position", "mod_prob", "read_name", "mod_code")
      if ("group" %in% names(data$sites))
        site_cols = c(site_cols, "group")
      .write_table_maybe_gz(data$sites[, site_cols, drop = FALSE],
                            paste0(prefix, "_sites.tsv"), gzip)
    } else {
      # BED format
      chrom = as.character(GenomicRanges::seqnames(data$region))
      strand = data$reads$strand[match(data$sites$read_name, data$reads$read_name)]
      strand[is.na(strand) | !strand %in% c("+", "-")] = "."
      bed = data.frame(
        chrom      = rep(chrom, nrow(data$sites)),
        chromStart = as.integer(data$sites$position - 1L),
        chromEnd   = as.integer(data$sites$position),
        name       = sprintf("%s:%s", data$sites$read_name, data$sites$mod_code),
        score      = as.integer(round(pmin(pmax(data$sites$mod_prob, 0), 1) * 1000)),
        strand     = strand,
        stringsAsFactors = FALSE
      )
      bed = bed[order(bed$chrom, bed$chromStart, bed$name), , drop = FALSE]
      .write_table_maybe_gz(bed, paste0(prefix, "_sites.bed"), gzip,
                            header = FALSE)
    }
  }

  invisible(NULL)
}

#' Write a tab-separated table, optionally gzip-compressed
#'
#' @param df Data frame to write.
#' @param path Output file path (without `.gz`; `.gz` appended when `gzip = TRUE`).
#' @param gzip Logical. Compress the file?
#' @param header Logical. Write a header line? `FALSE` for BED output.
#'
#' @return `path` (with `.gz` suffix if applicable), invisibly.
#'
#' @keywords internal
.write_table_maybe_gz = function(df, path, gzip, header = TRUE) {
  if (gzip) path = paste0(path, ".gz")
  con = if (gzip) gzfile(path, "w") else file(path, "w")
  on.exit(close(con))
  utils::write.table(df, con, sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = header)
  invisible(path)
}
