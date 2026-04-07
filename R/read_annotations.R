# Package-level cache for TxDb objects built from GTF/GFF files.
# Keyed by normalised file path; avoids re-parsing large files within a session.
.annotation_cache <- new.env(parent = emptyenv())

# Internal helper: return list(txdb, granges) for `gtf`, building and caching on first use.
.get_or_build_txdb <- function(gtf) {
  if (!requireNamespace("GenomicFeatures", quietly = TRUE)) {
    stop(
      "The 'GenomicFeatures' package is required but not installed. ",
      "Install it with: BiocManager::install('GenomicFeatures')",
      call. = FALSE
    )
  }
  if (!requireNamespace("rtracklayer", quietly = TRUE)) {
    stop(
      "The 'rtracklayer' package is required to import GTF/GFF files. ",
      "Install it with: BiocManager::install('rtracklayer')",
      call. = FALSE
    )
  }

  key <- normalizePath(gtf, mustWork = TRUE)

  if (exists(key, envir = .annotation_cache, inherits = FALSE)) {
    message("Using cached TxDb for: ", key)
    return(.annotation_cache[[key]])
  }

  message("Importing GTF/GFF file — this may take a moment for large files ...")
  gr_gtf <- rtracklayer::import(gtf)
  txdb   <- GenomicFeatures::makeTxDbFromGRanges(gr_gtf)
  result <- list(txdb = txdb, granges = gr_gtf)
  assign(key, result, envir = .annotation_cache)
  result
}

# Internal helper: resolve gene names for transcripts.
# Uses GTF gene_name attribute when available, falls back to AnnotationDbi lookup,
# then to tx_name.
.resolve_gene_names <- function(txs, txdb, gtf_granges = NULL) {
  if (!is.null(gtf_granges)) {
    # Filter to transcript rows only
    tx_rows <- gtf_granges[gtf_granges$type == "transcript"]
    # Match on transcript_id metadata column
    tx_ids_gtf <- tx_rows$transcript_id
    idx <- match(as.character(txs$tx_name), tx_ids_gtf)
    # Extract gene_name, fall back to gene_id
    if (!is.null(tx_rows$gene_name)) {
      gene_names <- tx_rows$gene_name[idx]
      gene_fallback <- if (!is.null(tx_rows$gene_id)) tx_rows$gene_id[idx] else as.character(txs$tx_name)
      gene_name_vec <- ifelse(is.na(gene_names) | gene_names == "", gene_fallback, gene_names)
    } else if (!is.null(tx_rows$gene_id)) {
      gene_name_vec <- tx_rows$gene_id[idx]
      gene_name_vec <- ifelse(is.na(gene_name_vec), as.character(txs$tx_name), gene_name_vec)
    } else {
      gene_name_vec <- as.character(txs$tx_name)
    }
    return(as.character(gene_name_vec))
  }

  # No GTF — try AnnotationDbi lookup
  gene_name_vec <- as.character(txs$tx_name)  # fallback

  gene_map <- tryCatch(
    AnnotationDbi::select(
      txdb,
      keys    = as.character(txs$tx_name),
      keytype = "TXNAME",
      columns = "GENEID"
    ),
    error = function(e) NULL
  )

  if (!is.null(gene_map) && "GENEID" %in% names(gene_map) &&
      any(!is.na(gene_map$GENEID))) {
    idx <- match(as.character(txs$tx_name), gene_map$TXNAME)
    mapped_gene <- gene_map$GENEID[idx]
    candidate <- ifelse(is.na(mapped_gene), as.character(txs$tx_name), mapped_gene)

    # If results look like Entrez IDs (all non-NA values are numeric), try mapIds for SYMBOL
    non_na <- candidate[!is.na(mapped_gene)]
    if (length(non_na) > 0L && all(grepl("^[0-9]+$", non_na))) {
      symbol_map <- tryCatch(
        AnnotationDbi::mapIds(
          txdb,
          keys    = as.character(txs$tx_name),
          keytype = "TXNAME",
          column  = "SYMBOL"
        ),
        error = function(e) NULL
      )
      if (!is.null(symbol_map) && any(!is.na(symbol_map))) {
        sym_idx <- match(as.character(txs$tx_name), names(symbol_map))
        sym_vals <- symbol_map[sym_idx]
        candidate <- ifelse(is.na(sym_vals) | sym_vals == "", candidate, sym_vals)
      }
    }
    gene_name_vec <- candidate
  }

  gene_name_vec
}

# Internal helper: collapse to one canonical transcript per gene.
# Prefers transcript with longest total CDS; ties broken by longest transcript span.
.select_canonical_transcripts <- function(transcripts_df, exons_df, cds_df = NULL) {
  if (nrow(transcripts_df) == 0L) {
    return(list(transcripts = transcripts_df, exons = exons_df))
  }

  if (!is.null(cds_df) && nrow(cds_df) > 0L) {
    # Compute total CDS length per tx_id
    cds_lengths <- tapply(
      cds_df$cds_end - cds_df$cds_start + 1L,
      cds_df$tx_id,
      sum
    )
    transcripts_df$cds_length <- cds_lengths[as.character(transcripts_df$tx_id)]
    transcripts_df$cds_length[is.na(transcripts_df$cds_length)] <- 0L
    transcripts_df$tx_span <- transcripts_df$tx_end - transcripts_df$tx_start

    # Per gene: pick tx with max cds_length; tie-break by tx_span
    selected_ids <- by(transcripts_df, transcripts_df$gene_name, function(grp) {
      grp <- grp[order(-grp$cds_length, -grp$tx_span), ]
      grp$tx_id[1L]
    })
  } else {
    # No CDS: pick tx with longest span per gene
    transcripts_df$tx_span <- transcripts_df$tx_end - transcripts_df$tx_start
    selected_ids <- by(transcripts_df, transcripts_df$gene_name, function(grp) {
      grp <- grp[order(-grp$tx_span), ]
      grp$tx_id[1L]
    })
  }

  selected_ids <- as.integer(unlist(selected_ids))

  # Clean up helper columns
  transcripts_df$tx_span    <- NULL
  if (!is.null(transcripts_df$cds_length)) transcripts_df$cds_length <- NULL

  transcripts_df <- transcripts_df[transcripts_df$tx_id %in% selected_ids, , drop = FALSE]
  exons_df       <- exons_df[exons_df$tx_id %in% selected_ids, , drop = FALSE]

  list(transcripts = transcripts_df, exons = exons_df)
}

#' Clear the internal TxDb annotation cache
#'
#' Removes all `TxDb` objects cached from GTF/GFF imports. Call this to free
#' memory or to force re-parsing of a GTF file that has changed on disk.
#'
#' @return `NULL`, invisibly.
#'
#' @examples
#' clear_annotation_cache()
#'
#' @export
clear_annotation_cache <- function() {
  rm(list = ls(.annotation_cache), envir = .annotation_cache)
  invisible(NULL)
}

#' Load gene annotations for a genomic region
#'
#' Extracts transcript and exon models overlapping a given genomic region from
#' either a pre-built `TxDb` object or a GTF/GFF file. The result is used as
#' the `annotations` argument to [plot_methylation()] to add a gene track.
#'
#' When a GTF/GFF path is supplied the file is parsed once and the resulting
#' `TxDb` is cached internally for the rest of the R session. Subsequent calls
#' with the same file path reuse the cache, so it is safe and efficient to call
#' `read_annotations()` for multiple regions without manually pre-loading.
#'
#' @param txdb A `TxDb` object (from the `GenomicFeatures` package), or `NULL`.
#'   Exactly one of `txdb` or `gtf` must be provided.
#' @param gtf Character. Path to a GTF or GFF annotation file, or `NULL`.
#'   Exactly one of `txdb` or `gtf` must be provided. On the first call the
#'   file is imported with `rtracklayer::import()` and converted to a `TxDb`
#'   via `GenomicFeatures::makeTxDbFromGRanges()`; subsequent calls with the
#'   same path reuse the cached `TxDb`. Use [clear_annotation_cache()] to free
#'   the cached object or force re-parsing.
#' @param region Character. Genomic region string in the format
#'   `"chr:start-end"` (e.g., `"chr1:1000000-2000000"`). The same format
#'   accepted by [read_methylation()].
#' @param collapse_transcripts Logical. When `TRUE` (default), only one
#'   canonical transcript per gene is retained (longest CDS, or longest span
#'   if no CDS data are available).
#'
#' @return A `gene_annotations` S3 object (a list) with elements:
#'   \describe{
#'     \item{transcripts}{Data frame with one row per transcript overlapping
#'       the region. Columns: `tx_id`, `tx_name`, `gene_name`, `strand`,
#'       `tx_start`, `tx_end`.}
#'     \item{exons}{Data frame with one row per exon. Columns: `tx_id`,
#'       `exon_start`, `exon_end`.}
#'     \item{cds}{Data frame with CDS intervals. Columns: `tx_id`,
#'       `cds_start`, `cds_end`.}
#'     \item{utr5}{Data frame with 5' UTR intervals. Columns: `tx_id`,
#'       `utr_start`, `utr_end`.}
#'     \item{utr3}{Data frame with 3' UTR intervals. Columns: `tx_id`,
#'       `utr_start`, `utr_end`.}
#'     \item{region}{A [GenomicRanges::GRanges] object for the queried region.}
#'   }
#'
#' @examples
#' \dontrun{
#' library(TxDb.Hsapiens.UCSC.hg38.knownGene)
#' txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
#' ann <- read_annotations(txdb = txdb, region = "chr1:1000000-2000000")
#'
#' # From a GTF file — parsed once, cached for subsequent calls
#' ann1 <- read_annotations(gtf = "Homo_sapiens.GRCh38.gtf", region = "chr1:1000000-2000000")
#' ann2 <- read_annotations(gtf = "Homo_sapiens.GRCh38.gtf", region = "chr2:5000000-6000000")
#' }
#'
#' @export
read_annotations <- function(txdb = NULL, gtf = NULL, region,
                             collapse_transcripts = TRUE) {
  # --- 1. Validate exactly one source ---
  has_txdb <- !is.null(txdb)
  has_gtf  <- !is.null(gtf)

  if (has_txdb && has_gtf) {
    stop("Provide exactly one of `txdb` or `gtf`, not both.", call. = FALSE)
  }
  if (!has_txdb && !has_gtf) {
    stop("One of `txdb` or `gtf` must be provided.", call. = FALSE)
  }

  # --- 2. Check GenomicFeatures availability ---
  if (!requireNamespace("GenomicFeatures", quietly = TRUE)) {
    stop(
      "The 'GenomicFeatures' package is required but not installed. ",
      "Install it with: BiocManager::install('GenomicFeatures')",
      call. = FALSE
    )
  }

  # --- 3. Build TxDb from GTF if needed (cached) ---
  gtf_granges <- NULL
  if (has_gtf) {
    result      <- .get_or_build_txdb(gtf)
    txdb        <- result$txdb
    gtf_granges <- result$granges
  }

  # --- 4. Parse region ---
  parsed    <- parse_region(region)
  region_gr <- GenomicRanges::GRanges(
    seqnames = parsed$chrom,
    ranges   = IRanges::IRanges(start = parsed$start, end = parsed$end)
  )

  # --- 5. Extract overlapping transcripts ---
  txs <- GenomicFeatures::transcriptsByOverlaps(txdb, region_gr)

  # Empty data frames for CDS/UTR (used in empty result and as fallback)
  empty_cds  <- data.frame(tx_id = integer(0L), cds_start = integer(0L),
                           cds_end = integer(0L), stringsAsFactors = FALSE)
  empty_utr  <- data.frame(tx_id = integer(0L), utr_start = integer(0L),
                           utr_end = integer(0L), stringsAsFactors = FALSE)

  # --- 6. Handle empty result ---
  if (length(txs) == 0L) {
    transcripts_df <- data.frame(
      tx_id     = integer(0L),
      tx_name   = character(0L),
      gene_name = character(0L),
      strand    = character(0L),
      tx_start  = integer(0L),
      tx_end    = integer(0L),
      stringsAsFactors = FALSE
    )
    exons_df <- data.frame(
      tx_id      = integer(0L),
      exon_start = integer(0L),
      exon_end   = integer(0L),
      stringsAsFactors = FALSE
    )
    return(structure(
      list(
        transcripts = transcripts_df,
        exons       = exons_df,
        cds         = empty_cds,
        utr5        = empty_utr,
        utr3        = empty_utr,
        region      = region_gr
      ),
      class = "gene_annotations"
    ))
  }

  # --- 7. Extract exons for those transcripts ---
  tx_ids       <- txs$tx_id   # integer IDs
  exons_by_tx  <- GenomicFeatures::exonsBy(txdb, by = "tx")
  exons_by_tx  <- exons_by_tx[names(exons_by_tx) %in% as.character(tx_ids)]
  exons_unlisted <- unlist(exons_by_tx)
  exons_raw    <- as.data.frame(exons_unlisted)
  # The tx_id is carried in names(exons_unlisted), not a column (column
  # 'group_name' was removed in newer Bioconductor versions)
  exons_df <- data.frame(
    tx_id      = as.integer(names(exons_unlisted)),
    exon_start = exons_raw$start,
    exon_end   = exons_raw$end,
    stringsAsFactors = FALSE
  )

  # --- 8. Extract CDS and UTR regions ---
  cds_df <- tryCatch({
    cds_by_tx <- GenomicFeatures::cdsBy(txdb, by = "tx")
    cds_by_tx <- cds_by_tx[names(cds_by_tx) %in% as.character(tx_ids)]
    if (length(cds_by_tx) == 0L) {
      empty_cds
    } else {
      cds_unlisted <- unlist(cds_by_tx)
      data.frame(
        tx_id     = as.integer(names(cds_unlisted)),
        cds_start = GenomicRanges::start(cds_unlisted),
        cds_end   = GenomicRanges::end(cds_unlisted),
        stringsAsFactors = FALSE
      )
    }
  }, error = function(e) empty_cds)

  utr5_df <- tryCatch({
    utr5_by_tx <- GenomicFeatures::fiveUTRsByTranscript(txdb)
    utr5_by_tx <- utr5_by_tx[names(utr5_by_tx) %in% as.character(tx_ids)]
    if (length(utr5_by_tx) == 0L) {
      empty_utr
    } else {
      utr5_unlisted <- unlist(utr5_by_tx)
      data.frame(
        tx_id     = as.integer(names(utr5_unlisted)),
        utr_start = GenomicRanges::start(utr5_unlisted),
        utr_end   = GenomicRanges::end(utr5_unlisted),
        stringsAsFactors = FALSE
      )
    }
  }, error = function(e) empty_utr)

  utr3_df <- tryCatch({
    utr3_by_tx <- GenomicFeatures::threeUTRsByTranscript(txdb)
    utr3_by_tx <- utr3_by_tx[names(utr3_by_tx) %in% as.character(tx_ids)]
    if (length(utr3_by_tx) == 0L) {
      empty_utr
    } else {
      utr3_unlisted <- unlist(utr3_by_tx)
      data.frame(
        tx_id     = as.integer(names(utr3_unlisted)),
        utr_start = GenomicRanges::start(utr3_unlisted),
        utr_end   = GenomicRanges::end(utr3_unlisted),
        stringsAsFactors = FALSE
      )
    }
  }, error = function(e) empty_utr)

  # --- 9. Resolve gene names ---
  gene_name_vec <- .resolve_gene_names(txs, txdb, if (has_gtf) gtf_granges else NULL)

  # --- 10. Build transcripts data.frame ---
  transcripts_df <- data.frame(
    tx_id     = as.integer(txs$tx_id),
    tx_name   = as.character(txs$tx_name),
    gene_name = gene_name_vec,
    strand    = as.character(BiocGenerics::strand(txs)),
    tx_start  = GenomicRanges::start(txs),
    tx_end    = GenomicRanges::end(txs),
    stringsAsFactors = FALSE
  )

  # --- 11. Optionally collapse to canonical transcripts ---
  if (collapse_transcripts) {
    canonical      <- .select_canonical_transcripts(transcripts_df, exons_df, cds_df)
    transcripts_df <- canonical$transcripts
    exons_df       <- canonical$exons
    # Filter CDS/UTR to retained tx_ids
    kept_ids <- transcripts_df$tx_id
    cds_df   <- cds_df[cds_df$tx_id   %in% kept_ids, , drop = FALSE]
    utr5_df  <- utr5_df[utr5_df$tx_id %in% kept_ids, , drop = FALSE]
    utr3_df  <- utr3_df[utr3_df$tx_id %in% kept_ids, , drop = FALSE]
  }

  # --- 12. Return S3 object ---
  structure(
    list(
      transcripts = transcripts_df,
      exons       = exons_df,
      cds         = cds_df,
      utr5        = utr5_df,
      utr3        = utr3_df,
      region      = region_gr
    ),
    class = "gene_annotations"
  )
}

#' Print a gene_annotations object
#'
#' @param x A `gene_annotations` object.
#' @param ... Additional arguments (ignored).
#'
#' @return `x`, invisibly.
#'
#' @export
print.gene_annotations <- function(x, ...) {
  chrom <- as.character(GenomicRanges::seqnames(x$region))
  start <- GenomicRanges::start(x$region)
  end   <- GenomicRanges::end(x$region)

  n_tx    <- nrow(x$transcripts)
  genes   <- if (n_tx > 0L) unique(x$transcripts$gene_name) else character(0L)
  n_genes <- length(genes)

  cat("gene_annotations object\n")
  cat(sprintf("Region: %s:%d-%d\n", chrom, start, end))
  cat(sprintf("Transcripts: %d\n", n_tx))
  cat(sprintf("Unique gene names: %d\n", n_genes))
  if (n_genes > 0L) {
    shown <- if (n_genes <= 10L) genes else c(genes[seq_len(10L)], "...")
    cat(sprintf("  %s\n", paste(shown, collapse = ", ")))
  }
  invisible(x)
}
