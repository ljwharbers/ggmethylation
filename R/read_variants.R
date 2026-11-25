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
#' Only the INFO fields `SVTYPE`, `END` and `SVLEN` are read; genotypes are
#' not loaded. Fields declared with `Number=.` (as written by pbsv, SVIM or
#' cuteSV) are supported; their first value is used.
#'
#' Limitations:
#' * Only the **first ALT allele** of a multi-allelic record is used.
#' * The **FILTER column is ignored**: records failing filters are returned
#'   like `PASS` records. Pre-filter the VCF (e.g. `bcftools view -f PASS`)
#'   if needed.
#' * Records whose first ALT is a spanning deletion (`*`) are dropped.
#'
#' If the region's chromosome is neither declared in the VCF header nor
#' present in its tabix index (e.g. `"chr1"` vs `"1"`), a warning is issued
#' and an empty `variant_data` object is returned.
#'
#' @param vcf Character. Path to a bgzipped VCF file (`.vcf.gz`) with a
#'   corresponding tabix index (`.vcf.gz.tbi`).
#' @param region Character. Genomic region string, e.g. `"chr1:1000-2000"`.
#'   The same format accepted by [read_methylation()].
#' @param bnd_match_tol Single non-negative number. Position tolerance (bp)
#'   used by [plot_methylation()] when matching reads' supplementary-alignment
#'   (SA) position to the mate position of the BND calls here; matched reads
#'   get a border on their SA indicator (requires `show_supplementary = TRUE`).
#'   Default `50`.
#'
#' @return A `variant_data` object (S3 list) with elements:
#'   \describe{
#'     \item{variants}{Data.frame with columns `position` (integer),
#'       `ref` (character), `alt` (character), `type` (character:
#'       `"SNV"`, `"MNV"`, `"insertion"`, `"deletion"`, `"BND"`, `"DEL"`,
#'       `"DUP"`, `"INV"`, another symbolic type such as `"CNV"`, or `NA`
#'       when the alleles are missing), `end` (integer), `mate_chrom`
#'       (character), and `mate_pos` (integer).}
#'     \item{region}{A [GenomicRanges::GRanges] object for the queried region.}
#'     \item{bnd_match_tol}{The `bnd_match_tol` value.}
#'   }
#'
#' @examples
#' \dontrun{
#' vd = read_variants("calls.vcf.gz", "chr1:1000000-1050000")
#' print(vd)
#' }
#'
#' @export
read_variants = function(vcf, region, bnd_match_tol = 50L) {
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

  if (!is.numeric(bnd_match_tol) || length(bnd_match_tol) != 1L ||
      is.na(bnd_match_tol) || bnd_match_tol < 0) {
    stop("`bnd_match_tol` must be a single non-negative number.", call. = FALSE)
  }

  # Parse region string
  region_gr = region_to_granges(region)
  chrom     = as.character(GenomicRanges::seqnames(region_gr))

  # Check the region's chromosome against the VCF's contigs (header ##contig
  # lines and the tabix index). A naming mismatch such as "chr1" vs "1"
  # otherwise surfaces as a cryptic readVcf() error or a silently empty result.
  header = tryCatch(
    VariantAnnotation::scanVcfHeader(vcf),
    error = function(e) {
      stop("Failed to read VCF file '", vcf, "': ", conditionMessage(e),
           call. = FALSE)
    }
  )
  header_contigs = as.character(
    GenomicRanges::seqnames(GenomicRanges::seqinfo(header))
  )
  tabix_contigs = tryCatch(Rsamtools::seqnamesTabix(vcf),
                           error = function(e) character(0L))
  if (!chrom %in% c(header_contigs, tabix_contigs)) {
    warning(
      "Chromosome '", chrom, "' of `region` is not among the VCF's contigs (",
      .format_contigs(unique(c(header_contigs, tabix_contigs))),
      "); returning no variants. ",
      "Check the chromosome naming (e.g. 'chr1' vs '1').",
      call. = FALSE
    )
    return(.new_variant_data(.empty_variants_df(), region_gr, bnd_match_tol))
  }
  if (length(tabix_contigs) > 0L && !chrom %in% tabix_contigs) {
    # Declared in the header but no records on it (scanTabix would error)
    return(.new_variant_data(.empty_variants_df(), region_gr, bnd_match_tol))
  }

  # Read VCF for the requested region: only the INFO fields used below and no
  # genotypes. Requesting INFO fields the header does not declare warns, so
  # restrict to those present.
  info_fields = intersect(c("SVTYPE", "END", "SVLEN"),
                          rownames(VariantAnnotation::info(header)))
  vcf_obj = tryCatch(
    VariantAnnotation::readVcf(
      vcf,
      # A genome label (rather than the header's Seqinfo) lets VCFs without
      # ##contig lines be queried; the chromosome was checked above.
      genome = "unknown",
      param = VariantAnnotation::ScanVcfParam(
        which = region_gr,
        info  = if (length(info_fields) > 0L) info_fields else NA,
        geno  = NA
      )
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
    return(.new_variant_data(.empty_variants_df(), region_gr, bnd_match_tol))
  }

  # Extract per-variant data
  rr      = SummarizedExperiment::rowRanges(vcf_obj)
  pos     = GenomicRanges::start(rr)
  ref_vec = as.character(VariantAnnotation::ref(vcf_obj))
  alt_list = VariantAnnotation::alt(vcf_obj)
  alt_vec  = as.character(unlist(lapply(alt_list, function(x) as.character(x)[1L])))

  # Extract INFO fields (NULL if absent in VCF). Fields declared Number=. or
  # Number=A (pbsv, SVIM, cuteSV, ...) come back as list-like columns
  # (e.g. CompressedIntegerList); keep their first value.
  info_df    = VariantAnnotation::info(vcf_obj)
  info_col   = function(nm) if (nm %in% names(info_df)) info_df[[nm]] else NULL
  svtype_vec = .info_first(info_col("SVTYPE"), as.character, NA_character_)
  end_info   = .info_first(info_col("END"),    as.integer,   NA_integer_)
  svlen_info = .info_first(info_col("SVLEN"),  as.numeric,   NA_real_)

  # Classify each variant
  n = length(pos)
  classified = vector("list", n)
  for (i in seq_len(n)) {
    svtype_i = if (!is.null(svtype_vec)) svtype_vec[i] else NA_character_
    end_i    = if (!is.null(end_info))   end_info[i]   else NA_integer_
    svlen_i  = if (!is.null(svlen_info)) svlen_info[i] else NA_real_
    classified[[i]] = classify_variant_row(
      alt    = alt_vec[i],
      ref    = ref_vec[i],
      pos    = pos[i],
      svtype = svtype_i,
      end    = end_i,
      svlen  = svlen_i
    )
  }

  type       = vapply(classified, `[[`, character(1L), "type")
  end_vec    = vapply(classified, `[[`, integer(1L),   "end")
  mate_chrom = vapply(classified, `[[`, character(1L), "mate_chrom")
  mate_pos   = vapply(classified, `[[`, integer(1L),   "mate_pos")

  variants_df = data.frame(
    position   = pos,
    ref        = ref_vec,
    alt        = alt_vec,
    type       = type,
    end        = end_vec,
    mate_chrom = mate_chrom,
    mate_pos   = mate_pos,
    stringsAsFactors = FALSE
  )
  # A spanning-deletion ALT ("*") only says the site is covered by a deletion
  # described in another record; it is not a variant of its own.
  variants_df = variants_df[is.na(variants_df$alt) | variants_df$alt != "*", ,
                            drop = FALSE]
  rownames(variants_df) = NULL

  .new_variant_data(variants_df, region_gr, bnd_match_tol)
}

.empty_variants_df = function() {
  data.frame(
    position   = integer(0L),
    ref        = character(0L),
    alt        = character(0L),
    type       = character(0L),
    end        = integer(0L),
    mate_chrom = character(0L),
    mate_pos   = integer(0L),
    stringsAsFactors = FALSE
  )
}

# Up to ten contig names for a message, e.g. "chr1, chr2, ... (25 total)"
.format_contigs = function(contigs) {
  if (length(contigs) == 0L) return("none found")
  shown = paste(utils::head(contigs, 10L), collapse = ", ")
  if (length(contigs) > 10L) {
    shown = paste0(shown, ", ... (", length(contigs), " total)")
  }
  shown
}

# Flatten one INFO column to a plain vector with one value per record.
# VariantAnnotation returns fields declared with Number=. / A / R as list-like
# columns: base lists or S4 Lists such as CompressedIntegerList /
# CompressedCharacterList, for which is.list() is FALSE. Keep the first value
# of each record (NA when empty), coerced with `cast`.
.info_first = function(x, cast, na) {
  if (is.null(x)) return(NULL)
  if (is.atomic(x)) return(cast(x))
  vapply(as.list(x), function(e) if (length(e) > 0L) cast(e[[1L]]) else na, na)
}

.new_variant_data = function(variants, region, bnd_match_tol) {
  structure(
    list(variants = variants, region = region,
         bnd_match_tol = as.integer(bnd_match_tol)),
    class = "variant_data"
  )
}

#' Classify a single variant record
#'
#' Internal helper used by [read_variants()] to classify one VCF row into a
#' type and derive `end`, `mate_chrom`, and `mate_pos`.
#'
#' Short variants are typed from the allele lengths: `"SNV"` (1 vs 1 base),
#' `"MNV"` (equal lengths > 1), `"insertion"` (REF shorter than ALT) or
#' `"deletion"` (REF longer than ALT). A spanning-deletion ALT (`"*"`), like an
#' `NA` allele, gives type `NA`.
#'
#' @param alt    Character. The ALT allele string.
#' @param ref    Character. The REF allele string.
#' @param pos    Integer. The 1-based start position of the variant.
#' @param svtype Character or NA. The INFO/SVTYPE value.
#' @param end    Integer or NA. The INFO/END value.
#' @param svlen  Numeric or NA (may be a list or S4 List). The INFO/SVLEN value.
#'
#' @return A named list with four elements: `type` (character), `end`
#'   (integer), `mate_chrom` (character), `mate_pos` (integer).
#'
#' @keywords internal
classify_variant_row = function(alt, ref, pos, svtype = NA_character_,
                                 end = NA_integer_, svlen = NA_real_) {
  # --- Symbolic SVs: <DEL>, <DUP>, <INV> ---
  if (!is.na(alt) && grepl("^<(DEL|DUP|INV)>$", alt)) {
    sv_symbol = sub("^<(.+)>$", "\\1", alt)

    # Determine end: prefer INFO/END, fall back to pos + abs(SVLEN)
    end_i = NA_integer_
    if (length(end) == 1L && !is.na(end)) {
      end_i = as.integer(end)
    }
    if (is.na(end_i)) {
      # SVLEN may be list-like (base list or S4 List) — take its first record
      sl = svlen
      if (!is.atomic(sl)) sl = as.list(sl)[[1L]]
      if (length(sl) >= 1L && !is.na(sl[1L])) {
        end_i = as.integer(pos + abs(sl[1L]))
      }
    }
    if (is.na(end_i)) end_i = as.integer(pos)

    return(list(
      type       = sv_symbol,
      end        = end_i,
      mate_chrom = NA_character_,
      mate_pos   = NA_integer_
    ))
  }

  # --- Unknown symbolic SV (e.g., <CNV>, <DUP:TANDEM>) ---
  # This must come after the DEL/DUP/INV branch and before BND.
  if (!is.na(alt) && grepl("^<[^>]+>$", alt, perl = TRUE)) {
    return(list(
      type       = sub("^<(.+)>$", "\\1", alt),
      end        = as.integer(if (length(end) == 1L && !is.na(end)) end else pos),
      mate_chrom = NA_character_,
      mate_pos   = NA_integer_
    ))
  }

  # --- BND: SVTYPE == "BND" or ALT matches a breakend bracket pattern ---
  is_bnd = (!is.na(svtype) && svtype == "BND") ||
            (!is.na(alt) && grepl("[\\[\\]]", alt, perl = TRUE))

  if (is_bnd) {
    parsed_bnd = parse_bnd_alt(alt)
    return(list(
      type       = "BND",
      end        = as.integer(pos),
      mate_chrom = parsed_bnd$mate_chrom,
      mate_pos   = parsed_bnd$mate_pos
    ))
  }

  # --- NA guard: alt or ref is NA, or a spanning-deletion "*" allele ---
  if (is.na(alt) || is.na(ref) || alt == "*") {
    return(list(
      type       = NA_character_,
      end        = as.integer(pos),
      mate_chrom = NA_character_,
      mate_pos   = NA_integer_
    ))
  }

  # --- Standard SNV / MNV / insertion / deletion ---
  n_ref = nchar(ref)
  n_alt = nchar(alt)
  tp = if (n_ref == 1L && n_alt == 1L) {
    "SNV"
  } else if (n_ref == n_alt) {
    "MNV"
  } else if (n_ref < n_alt) {
    "insertion"
  } else {
    "deletion"
  }
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
parse_bnd_alt = function(alt_str) {
  if (is.na(alt_str) || !grepl("[\\[\\]]", alt_str, perl = TRUE)) {
    return(list(mate_chrom = NA_character_, mate_pos = NA_integer_))
  }

  # Match bracket-enclosed chr:pos in any of the four orientations.
  # Pattern: one or more characters inside [ ] that contain at least one colon.
  m = regmatches(alt_str, regexpr("[\\[\\]]([^\\[\\]]+:[0-9]+)[\\[\\]]",
                                   alt_str, perl = TRUE))
  if (length(m) == 0L || nchar(m) == 0L) {
    return(list(mate_chrom = NA_character_, mate_pos = NA_integer_))
  }

  # Strip the surrounding brackets to get "chr:pos"
  inner = sub("^[\\[\\]]", "", sub("[\\[\\]]$", "", m, perl = TRUE), perl = TRUE)
  parts = strsplit(inner, ":", fixed = TRUE)[[1L]]

  if (length(parts) < 2L) {
    return(list(mate_chrom = NA_character_, mate_pos = NA_integer_))
  }

  # chromosome is everything except the last part (handles "chr1", "chr_Un_gl000220" etc.)
  mate_chrom = paste(parts[-length(parts)], collapse = ":")
  mate_pos   = suppressWarnings(as.integer(parts[length(parts)]))

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
print.variant_data = function(x, ...) {
  chrom      = as.character(GenomicRanges::seqnames(x$region))
  start      = GenomicRanges::start(x$region)
  region_end = GenomicRanges::end(x$region)

  cat("variant_data object\n")
  cat(sprintf("Region: %s:%d-%d\n", chrom, start, region_end))
  cat(sprintf("Variants: %d total\n", nrow(x$variants)))

  if (nrow(x$variants) > 0L) {
    type_counts = table(x$variants$type)
    for (tp in c("SNV", "MNV", "insertion", "deletion", "BND", "DEL", "DUP", "INV")) {
      n = if (tp %in% names(type_counts)) type_counts[[tp]] else 0L
      cat(sprintf("  %s: %d\n", tp, n))
    }
  }

  invisible(x)
}
