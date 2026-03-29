#' Write methylation data to disk
#'
#' Saves the `$reads` and/or `$sites` data frames from a `methylation_data`
#' object to tab-separated (TSV) or BED files. Files can optionally be
#' gzip-compressed.
#'
#' @param data A `methylation_data` object as returned by [read_methylation()].
#' @param prefix Character. Path prefix for output files. The appropriate
#'   suffix (`_reads.tsv`, `_sites.tsv`, or `_sites.bed`) is appended
#'   automatically. Intermediate directories are created if they do not exist.
#' @param format Character. Output format for the sites file: `"tsv"` (default)
#'   or `"bed"` (6-column BED, 0-based coordinates, no header).
#' @param reads Logical. Write the `$reads` data frame? Default `TRUE`.
#' @param sites Logical. Write the `$sites` data frame? Default `TRUE`.
#' @param gzip Logical. Compress output files with gzip? Default `FALSE`.
#'   When `TRUE`, `.gz` is appended to each file name.
#' @param ... Reserved for future use; currently ignored.
#'
#' @return `data`, invisibly. Allows use in a pipeline.
#'
#' @details
#' **Reads file** (`{prefix}_reads.tsv[.gz]`)
#'
#' Always written as TSV with a header line. Columns: `read_name`, `start`,
#' `end`, `strand`, `group` (only when the object carries a group tag),
#' `mean_mod_prob`. `mean_mod_prob` is the per-read mean of `mod_prob` values
#' across all sites; reads with no called sites receive `NA`.
#'
#' **Sites file — TSV** (`{prefix}_sites.tsv[.gz]`)
#'
#' Columns: `position`, `mod_prob`, `read_name`, `mod_code`, `group`
#' (if present). Coordinates are 1-based (as stored in the object).
#'
#' **Sites file — BED** (`{prefix}_sites.bed[.gz]`)
#'
#' 6-column BED (no header). Coordinates are converted to 0-based half-open
#' intervals (`chromStart = position - 1`, `chromEnd = position`). The `score`
#' column contains `mod_prob` (float in [0, 1]).
#'
#' @examples
#' \dontrun{
#' md <- read_methylation("sample.bam", "chr1:1000-2000")
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
#' # Pipe-friendly usage
#' md |>
#'   write_methylation(prefix = "output/sample_chr1") |>
#'   plot_methylation()
#' }
#'
#' @seealso [read_methylation()]
#' @export
write_methylation <- function(data,
                              prefix,
                              format = c("tsv", "bed"),
                              reads  = TRUE,
                              sites  = TRUE,
                              gzip   = FALSE,
                              ...) {

  # --- Input validation ---
  if (!inherits(data, "methylation_data"))
    stop("'data' must be a methylation_data object.", call. = FALSE)
  if (!is.character(prefix) || length(prefix) != 1L || nchar(prefix) == 0L)
    stop("'prefix' must be a non-empty character string.", call. = FALSE)
  format <- match.arg(format)
  if (!reads && !sites) {
    warning("Nothing to write: both 'reads' and 'sites' are FALSE.")
    return(invisible(data))
  }

  # --- Create output directory if needed ---
  out_dir <- dirname(prefix)
  if (!is.na(out_dir) && out_dir != "." && !dir.exists(out_dir))
    dir.create(out_dir, recursive = TRUE)

  # --- Compute mean_mod_prob per read ---
  if (nrow(data$sites) > 0L) {
    mp <- tapply(data$sites$mod_prob, data$sites$read_name, mean, na.rm = TRUE)
    data$reads$mean_mod_prob <- as.numeric(mp[data$reads$read_name])
  } else {
    data$reads$mean_mod_prob <- NA_real_
  }

  # --- Write reads ---
  if (reads) {
    read_cols <- c("read_name", "start", "end", "strand")
    if (!is.null(data$group_tag) && "group" %in% names(data$reads))
      read_cols <- c(read_cols, "group")
    read_cols <- c(read_cols, "mean_mod_prob")
    .write_tsv_maybe_gz(data$reads[, read_cols, drop = FALSE],
                        paste0(prefix, "_reads.tsv"), gzip)
  }

  # --- Write sites ---
  if (sites) {
    if (format == "tsv") {
      site_cols <- c("position", "mod_prob", "read_name", "mod_code")
      if ("group" %in% names(data$sites))
        site_cols <- c(site_cols, "group")
      .write_tsv_maybe_gz(data$sites[, site_cols, drop = FALSE],
                          paste0(prefix, "_sites.tsv"), gzip)
    } else {
      # BED format
      chrom <- as.character(GenomicRanges::seqnames(data$region))
      strand_map <- data$reads$strand
      names(strand_map) <- data$reads$read_name
      bed <- data.frame(
        chrom      = chrom,
        chromStart = data$sites$position - 1L,
        chromEnd   = data$sites$position,
        name       = data$sites$read_name,
        score      = data$sites$mod_prob,
        strand     = strand_map[data$sites$read_name],
        stringsAsFactors = FALSE
      )
      .write_bed_maybe_gz(bed, paste0(prefix, "_sites.bed"), gzip)
    }
  }

  invisible(data)
}

# Internal helpers --------------------------------------------------------

#' Write a data frame as TSV, optionally gzip-compressed
#'
#' @param df Data frame to write.
#' @param path Output file path (without `.gz`; `.gz` appended when `gzip = TRUE`).
#' @param gzip Logical. Compress the file?
#'
#' @return `path` (with `.gz` suffix if applicable), invisibly.
#'
#' @keywords internal
.write_tsv_maybe_gz <- function(df, path, gzip) {
  if (gzip) path <- paste0(path, ".gz")
  con <- if (gzip) gzfile(path, "w") else file(path, "w")
  on.exit(close(con))
  write.table(df, con, sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = TRUE)
  invisible(path)
}

#' Write a data frame as BED (no header), optionally gzip-compressed
#'
#' @param df Data frame to write.
#' @param path Output file path (without `.gz`; `.gz` appended when `gzip = TRUE`).
#' @param gzip Logical. Compress the file?
#'
#' @return `path` (with `.gz` suffix if applicable), invisibly.
#'
#' @keywords internal
.write_bed_maybe_gz <- function(df, path, gzip) {
  if (gzip) path <- paste0(path, ".gz")
  con <- if (gzip) gzfile(path, "w") else file(path, "w")
  on.exit(close(con))
  write.table(df, con, sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = FALSE)
  invisible(path)
}
