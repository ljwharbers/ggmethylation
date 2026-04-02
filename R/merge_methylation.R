#' Merge multiple methylation data objects into a multi-sample object
#'
#' Combines two or more [`methylation_data`][read_methylation()] objects
#' into a single `multi_methylation_data` object suitable for multi-sample
#' plotting with [plot_methylation()].
#'
#' All samples must cover the same genomic region and use the same
#' modification code(s). Every sample must be supplied with a unique name,
#' either via R's named-argument syntax (`merge_methylation(A = md1, B = md2)`)
#' or via the `.list` parameter (`merge_methylation(.list = list(A = md1, B = md2))`).
#'
#' @param ... Named `methylation_data` objects. Every argument **must** be
#'   named; an error is raised if any name is missing or empty.
#' @param .list A named list of `methylation_data` objects. When provided,
#'   `...` is ignored. This is useful for programmatic construction of the
#'   sample list.
#'
#' @return A `multi_methylation_data` object (S3 list) with elements:
#'   \describe{
#'     \item{samples}{A named list of `methylation_data` objects, one per
#'       sample.}
#'     \item{region}{A [GenomicRanges::GRanges] object for the shared
#'       genomic region.}
#'     \item{mod_code}{Character vector of modification code(s) shared by
#'       all samples.}
#'   }
#'
#' @examples
#' \dontrun{
#' md1 <- read_methylation("sample1.bam", "chr1:1000-2000")
#' md2 <- read_methylation("sample2.bam", "chr1:1000-2000")
#' merged <- merge_methylation(sampleA = md1, sampleB = md2)
#' print(merged)
#' plot_methylation(merged)
#' }
#'
#' @export
merge_methylation <- function(..., .list = NULL) {
  # --- 1. Collect samples ---
  if (!is.null(.list)) {
    samples_list <- .list
  } else {
    samples_list <- list(...)
  }

  # --- 2. Validate names ---
  nms <- names(samples_list)
  if (is.null(nms) || any(!nzchar(nms))) {
    stop(
      "All samples passed to `merge_methylation()` must be named.",
      call. = FALSE
    )
  }

  # --- 3. Validate all are methylation_data ---
  for (nm in nms) {
    if (!inherits(samples_list[[nm]], "methylation_data")) {
      stop(
        sprintf(
          "Sample '%s' is not a `methylation_data` object. ",
          nm
        ),
        "Use `read_methylation()` to create samples.",
        call. = FALSE
      )
    }
  }

  # --- 4. Validate shared region ---
  chroms <- vapply(
    samples_list,
    function(x) as.character(GenomicRanges::seqnames(x$region)),
    character(1L)
  )
  starts <- vapply(
    samples_list,
    function(x) GenomicRanges::start(x$region),
    integer(1L)
  )
  ends <- vapply(
    samples_list,
    function(x) GenomicRanges::end(x$region),
    integer(1L)
  )

  if (length(unique(chroms)) > 1L || length(unique(starts)) > 1L ||
      length(unique(ends)) > 1L) {
    region_strs <- sprintf(
      "%s:%d-%d", chroms, starts, ends
    )
    stop(
      "All samples must cover the same genomic region. Found:\n",
      paste(sprintf("  %s: %s", nms, region_strs), collapse = "\n"),
      call. = FALSE
    )
  }

  # --- 5. Validate shared mod_code ---
  first_code <- samples_list[[1L]]$mod_code
  for (nm in nms[-1L]) {
    if (!identical(samples_list[[nm]]$mod_code, first_code)) {
      stop(
        sprintf(
          "Sample '%s' has mod_code '%s' but sample '%s' has mod_code '%s'. ",
          nm,
          paste(samples_list[[nm]]$mod_code, collapse = ", "),
          nms[1L],
          paste(first_code, collapse = ", ")
        ),
        "All samples must use the same modification code(s).",
        call. = FALSE
      )
    }
  }

  # --- 6. Return multi_methylation_data ---
  structure(
    list(
      samples  = samples_list,
      region   = samples_list[[1L]]$region,
      mod_code = samples_list[[1L]]$mod_code
    ),
    class = "multi_methylation_data"
  )
}

#' Print a multi_methylation_data object
#'
#' @param x A `multi_methylation_data` object.
#' @param ... Additional arguments (ignored).
#'
#' @return `x`, invisibly.
#'
#' @export
print.multi_methylation_data <- function(x, ...) {
  chrom <- as.character(GenomicRanges::seqnames(x$region))
  start <- GenomicRanges::start(x$region)
  end   <- GenomicRanges::end(x$region)

  cat("multi_methylation_data object\n")
  cat(sprintf("Region:   %s:%d-%d\n", chrom, start, end))
  cat(sprintf("Samples:  %d\n", length(x$samples)))
  cat(sprintf("Mod code: %s\n", paste(x$mod_code, collapse = ", ")))
  cat("\nPer-sample summary:\n")

  for (nm in names(x$samples)) {
    s <- x$samples[[nm]]
    n_reads <- nrow(s$reads)

    if (!is.null(s$group_tag) && "group" %in% names(s$reads)) {
      grps <- sort(unique(s$reads$group[!is.na(s$reads$group)]))
      grp_counts <- vapply(grps, function(g) {
        sum(s$reads$group == g, na.rm = TRUE)
      }, integer(1L))
      grp_str <- paste(
        sprintf("%s=%d", grps, grp_counts),
        collapse = ", "
      )
      cat(sprintf("  %s: %d reads [%s]\n", nm, n_reads, grp_str))
    } else {
      cat(sprintf("  %s: %d reads\n", nm, n_reads))
    }
  }

  invisible(x)
}

#' Summary of a multi_methylation_data object
#'
#' @param object A `multi_methylation_data` object.
#' @param ... Unused.
#'
#' @return A named list of per-sample summaries (invisibly). Each element is
#'   the list returned by [summary.methylation_data()].
#'
#' @export
summary.multi_methylation_data <- function(object, ...) {
  out <- lapply(names(object$samples), function(nm) {
    cat(sprintf("\n=== Sample: %s ===\n", nm))
    summary(object$samples[[nm]])
  })
  names(out) <- names(object$samples)
  invisible(out)
}
