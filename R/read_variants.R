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
#'       `ref` (character), `alt` (character), `type` (character:
#'       `"SNV"`, `"insertion"`, `"deletion"`, `"BND"`, `"DEL"`, `"DUP"`,
#'       or `"INV"`), `end` (integer), `mate_chrom` (character), and
#'       `mate_pos` (integer).}
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
      position   = integer(0L),
      ref        = character(0L),
      alt        = character(0L),
      type       = character(0L),
      end        = integer(0L),
      mate_chrom = character(0L),
      mate_pos   = integer(0L),
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

  # Extract INFO fields (may be NULL if absent in VCF)
  info_df   <- VariantAnnotation::info(vcf_obj)
  svtype_vec <- if ("SVTYPE" %in% names(info_df)) {
    sv <- info_df[["SVTYPE"]]
    if (is.list(sv)) {
      vapply(sv, function(x) if (length(x) > 0L) as.character(x[[1L]]) else NA_character_, character(1L))
    } else {
      as.character(sv)
    }
  } else {
    NULL
  }
  end_info   <- if ("END" %in% names(info_df)) {
    ei <- info_df[["END"]]
    if (is.list(ei)) {
      vapply(ei, function(x) if (length(x) > 0L) as.integer(x[[1L]]) else NA_integer_, integer(1L))
    } else {
      as.integer(ei)
    }
  } else {
    NULL
  }
  svlen_info <- if ("SVLEN"  %in% names(info_df)) info_df[["SVLEN"]]  else NULL

  # Classify each variant
  n <- length(pos)
  classified <- vector("list", n)
  for (i in seq_len(n)) {
    svtype_i <- if (!is.null(svtype_vec)) svtype_vec[i] else NA_character_
    end_i    <- if (!is.null(end_info))   end_info[i]   else NA_integer_
    svlen_i  <- if (!is.null(svlen_info)) svlen_info[i] else NA_real_
    classified[[i]] <- classify_variant_row(
      alt    = alt_vec[i],
      ref    = ref_vec[i],
      pos    = pos[i],
      svtype = svtype_i,
      end    = end_i,
      svlen  = svlen_i
    )
  }

  type       <- vapply(classified, `[[`, character(1L), "type")
  end_vec    <- vapply(classified, `[[`, integer(1L),   "end")
  mate_chrom <- vapply(classified, `[[`, character(1L), "mate_chrom")
  mate_pos   <- vapply(classified, `[[`, integer(1L),   "mate_pos")

  variants_df <- data.frame(
    position   = pos,
    ref        = ref_vec,
    alt        = alt_vec,
    type       = type,
    end        = end_vec,
    mate_chrom = mate_chrom,
    mate_pos   = mate_pos,
    stringsAsFactors = FALSE
  )
  rownames(variants_df) <- NULL

  structure(
    list(variants = variants_df, region = region_gr),
    class = "variant_data"
  )
}

#' Classify a single variant record
#'
#' Internal helper used by [read_variants()] to classify one VCF row into a
#' type and derive `end`, `mate_chrom`, and `mate_pos`.
#'
#' @param alt    Character. The ALT allele string.
#' @param ref    Character. The REF allele string.
#' @param pos    Integer. The 1-based start position of the variant.
#' @param svtype Character or NA. The INFO/SVTYPE value.
#' @param end    Integer or NA. The INFO/END value.
#' @param svlen  Numeric or NA (may be a list). The INFO/SVLEN value.
#'
#' @return A named list with four elements: `type` (character), `end`
#'   (integer), `mate_chrom` (character), `mate_pos` (integer).
#'
#' @keywords internal
classify_variant_row <- function(alt, ref, pos, svtype = NA_character_,
                                 end = NA_integer_, svlen = NA_real_) {
  # --- Symbolic SVs: <DEL>, <DUP>, <INV> ---
  if (!is.na(alt) && grepl("^<(DEL|DUP|INV)>$", alt)) {
    sv_symbol <- sub("^<(.+)>$", "\\1", alt)

    # Determine end: prefer INFO/END, fall back to pos + abs(SVLEN)
    end_i <- NA_integer_
    if (length(end) == 1L && !is.na(end)) {
      end_i <- as.integer(end)
    }
    if (is.na(end_i)) {
      # SVLEN may be a list (list-type INFO field) — unlist safely
      sl <- svlen
      if (is.list(sl)) sl <- sl[[1L]]
      if (length(sl) >= 1L && !is.na(sl[1L])) {
        end_i <- as.integer(pos + abs(sl[1L]))
      }
    }
    if (is.na(end_i)) end_i <- as.integer(pos)

    return(list(
      type       = sv_symbol,
      end        = end_i,
      mate_chrom = NA_character_,
      mate_pos   = NA_integer_
    ))
  }

  # --- BND: SVTYPE == "BND" or ALT matches a breakend bracket pattern ---
  is_bnd <- (!is.na(svtype) && svtype == "BND") ||
            (!is.na(alt) && grepl("[\\[\\]]", alt, perl = TRUE))

  if (is_bnd) {
    parsed_bnd <- parse_bnd_alt(alt)
    return(list(
      type       = "BND",
      end        = as.integer(pos),
      mate_chrom = parsed_bnd$mate_chrom,
      mate_pos   = parsed_bnd$mate_pos
    ))
  }

  # --- Standard SNV / insertion / deletion ---
  tp <- ifelse(
    nchar(ref) == 1L & nchar(alt) == 1L, "SNV",
    ifelse(nchar(ref) < nchar(alt), "insertion", "deletion")
  )
  list(
    type       = tp,
    end        = as.integer(pos),
    mate_chrom = NA_character_,
    mate_pos   = NA_integer_
  )
}

#' Parse a BND ALT string to extract mate chromosome and position
#'
#' Handles all four VCF breakend notation forms:
#' - `N[chr:pos[`
#' - `N]chr:pos]`
#' - `[chr:pos[N`
#' - `]chr:pos]N`
#'
#' where `N` is any nucleotide (A/T/G/C) or a dot.
#'
#' @param alt_str Character scalar. The ALT field value from a VCF BND record.
#'
#' @return A list with two elements:
#'   \describe{
#'     \item{mate_chrom}{Character: the mate chromosome, or `NA_character_` if
#'       `alt_str` does not match a BND pattern.}
#'     \item{mate_pos}{Integer: the mate position, or `NA_integer_` if
#'       `alt_str` does not match a BND pattern.}
#'   }
#'
#' @keywords internal
parse_bnd_alt <- function(alt_str) {
  if (is.na(alt_str) || !grepl("[\\[\\]]", alt_str, perl = TRUE)) {
    return(list(mate_chrom = NA_character_, mate_pos = NA_integer_))
  }

  # Match bracket-enclosed chr:pos in any of the four orientations.
  # Pattern: one or more characters inside [ ] that contain at least one colon.
  m <- regmatches(alt_str, regexpr("[\\[\\]]([^\\[\\]]+:[0-9]+)[\\[\\]]",
                                   alt_str, perl = TRUE))
  if (length(m) == 0L || nchar(m) == 0L) {
    return(list(mate_chrom = NA_character_, mate_pos = NA_integer_))
  }

  # Strip the surrounding brackets to get "chr:pos"
  inner <- sub("^[\\[\\]]", "", sub("[\\[\\]]$", "", m, perl = TRUE), perl = TRUE)
  parts <- strsplit(inner, ":", fixed = TRUE)[[1L]]

  if (length(parts) < 2L) {
    return(list(mate_chrom = NA_character_, mate_pos = NA_integer_))
  }

  # chromosome is everything except the last part (handles "chr1", "chr_Un_gl000220" etc.)
  mate_chrom <- paste(parts[-length(parts)], collapse = ":")
  mate_pos   <- suppressWarnings(as.integer(parts[length(parts)]))

  if (is.na(mate_pos)) {
    return(list(mate_chrom = NA_character_, mate_pos = NA_integer_))
  }

  list(mate_chrom = mate_chrom, mate_pos = mate_pos)
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
    for (tp in c("SNV", "insertion", "deletion", "BND", "DEL", "DUP", "INV")) {
      n <- if (tp %in% names(type_counts)) type_counts[[tp]] else 0L
      cat(sprintf("  %s: %d\n", tp, n))
    }
  }

  invisible(x)
}
