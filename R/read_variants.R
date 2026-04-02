#' Read variants from a VCF file
#'
#' Loads variant calls from a bgzipped, tabix-indexed VCF file for a given
#' genomic region. Returns a `variant_data` object containing a data.frame of
#' variants and the queried region as a `GRanges` object.
#'
#' Requires the `VariantAnnotation` and `SummarizedExperiment` packages from
#' Bioconductor. Install them with:
#' ```r
#' BiocManager::install(c("VariantAnnotation", "SummarizedExperiment"))
#' ```
#'
#' @param vcf Character. Path to a bgzipped VCF file (`.vcf.gz`) with a
#'   corresponding tabix index (`.vcf.gz.tbi`).
#' @param region Character. Genomic region string, e.g. `"chr1:1000-2000"`.
#'   The same format accepted by [read_methylation()].
#'
#' @return A `variant_data` object (S3 list) with elements:
#'   \describe{
#'     \item{variants}{Data.frame with columns `position` (integer),
#'       `ref` (character), `alt` (character), and `type` (character:
#'       `"SNV"`, `"insertion"`, or `"deletion"`).}
#'     \item{region}{A [GenomicRanges::GRanges] object for the queried region.}
#'   }
#'
#' @examples
#' \dontrun{
#' vd <- read_variants("calls.vcf.gz", "chr1:1000000-1050000")
#' print(vd)
#' }
#'
#' @export
read_variants <- function(vcf, region) {
  if (!requireNamespace("VariantAnnotation", quietly = TRUE)) {
    stop(
      "Package 'VariantAnnotation' is required to read VCF files. ",
      "Install it with: BiocManager::install('VariantAnnotation')",
      call. = FALSE
    )
  }
  if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    stop(
      "Package 'SummarizedExperiment' is required to read VCF files. ",
      "Install it with: BiocManager::install('SummarizedExperiment')",
      call. = FALSE
    )
  }

  # Parse region string
  parsed <- parse_region(region)
  region_gr <- GenomicRanges::GRanges(
    seqnames = parsed$chrom,
    ranges   = IRanges::IRanges(start = parsed$start, end = parsed$end)
  )

  # Read VCF for the requested region
  vcf_obj <- tryCatch(
    VariantAnnotation::readVcf(
      vcf,
      param = VariantAnnotation::ScanVcfParam(which = region_gr)
    ),
    error = function(e) {
      stop(
        "Failed to read VCF file '", vcf, "': ", conditionMessage(e),
        call. = FALSE
      )
    }
  )

  # Handle empty result (no variants in region)
  if (length(vcf_obj) == 0L) {
    variants_df <- data.frame(
      position = integer(0L),
      ref      = character(0L),
      alt      = character(0L),
      type     = character(0L),
      stringsAsFactors = FALSE
    )
    return(structure(
      list(variants = variants_df, region = region_gr),
      class = "variant_data"
    ))
  }

  # Extract per-variant data
  rr      <- SummarizedExperiment::rowRanges(vcf_obj)
  pos     <- GenomicRanges::start(rr)
  ref_vec <- as.character(VariantAnnotation::ref(vcf_obj))
  alt_list <- VariantAnnotation::alt(vcf_obj)
  alt_vec  <- as.character(unlist(lapply(alt_list, function(x) as.character(x)[1L])))

  # Classify variant type
  type <- ifelse(
    nchar(ref_vec) == 1L & nchar(alt_vec) == 1L, "SNV",
    ifelse(nchar(ref_vec) < nchar(alt_vec), "insertion", "deletion")
  )

  variants_df <- data.frame(
    position = pos,
    ref      = ref_vec,
    alt      = alt_vec,
    type     = type,
    stringsAsFactors = FALSE
  )
  rownames(variants_df) <- NULL

  structure(
    list(variants = variants_df, region = region_gr),
    class = "variant_data"
  )
}

#' Print a variant_data object
#'
#' @param x A `variant_data` object returned by [read_variants()].
#' @param ... Additional arguments (ignored).
#'
#' @return `x`, invisibly.
#'
#' @export
print.variant_data <- function(x, ...) {
  chrom <- as.character(GenomicRanges::seqnames(x$region))
  start <- GenomicRanges::start(x$region)
  end   <- GenomicRanges::end(x$region)

  cat("variant_data object\n")
  cat(sprintf("Region: %s:%d-%d\n", chrom, start, end))
  cat(sprintf("Variants: %d total\n", nrow(x$variants)))

  if (nrow(x$variants) > 0L) {
    type_counts <- table(x$variants$type)
    for (tp in c("SNV", "insertion", "deletion")) {
      n <- if (tp %in% names(type_counts)) type_counts[[tp]] else 0L
      cat(sprintf("  %s: %d\n", tp, n))
    }
  }

  invisible(x)
}
