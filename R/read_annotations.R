#' Load gene annotations for a genomic region
#'
#' Extracts transcript and exon models overlapping a given genomic region from
#' either a pre-built `TxDb` object or a GTF/GFF file. The result is used as
#' the `annotations` argument to [plot_methylation()] to add a gene track.
#'
#' @param txdb A `TxDb` object (from the `GenomicFeatures` package), or `NULL`.
#'   Exactly one of `txdb` or `gtf` must be provided.
#' @param gtf Character. Path to a GTF or GFF annotation file, or `NULL`.
#'   Exactly one of `txdb` or `gtf` must be provided. When supplied, the file
#'   is imported with `rtracklayer::import()` and converted to a `TxDb` via
#'   `GenomicFeatures::makeTxDbFromGRanges()`.
#' @param region Character. Genomic region string in the format
#'   `"chr:start-end"` (e.g., `"chr1:1000000-2000000"`). The same format
#'   accepted by [read_methylation()].
#'
#' @return A `gene_annotations` S3 object (a list) with elements:
#'   \describe{
#'     \item{transcripts}{Data frame with one row per transcript overlapping
#'       the region. Columns: `tx_id`, `tx_name`, `gene_name`, `strand`,
#'       `tx_start`, `tx_end`.}
#'     \item{exons}{Data frame with one row per exon. Columns: `tx_id`,
#'       `exon_start`, `exon_end`.}
#'     \item{region}{A [GenomicRanges::GRanges] object for the queried region.}
#'   }
#'
#' @examples
#' \dontrun{
#' library(TxDb.Hsapiens.UCSC.hg38.knownGene)
#' txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
#' ann <- read_annotations(txdb = txdb, region = "chr1:1000000-2000000")
#'
#' # From a GTF file
#' ann <- read_annotations(gtf = "Homo_sapiens.GRCh38.gtf", region = "1:1000000-2000000")
#' }
#'
#' @export
read_annotations <- function(txdb = NULL, gtf = NULL, region) {
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

  # --- 3. Build TxDb from GTF if needed ---
  if (has_gtf) {
    if (!requireNamespace("rtracklayer", quietly = TRUE)) {
      stop(
        "The 'rtracklayer' package is required to import GTF/GFF files. ",
        "Install it with: BiocManager::install('rtracklayer')",
        call. = FALSE
      )
    }
    gr_gtf <- rtracklayer::import(gtf)
    txdb   <- GenomicFeatures::makeTxDbFromGRanges(gr_gtf)
  }

  # --- 4. Parse region ---
  parsed    <- parse_region(region)
  region_gr <- GenomicRanges::GRanges(
    seqnames = parsed$chrom,
    ranges   = IRanges::IRanges(start = parsed$start, end = parsed$end)
  )

  # --- 5. Extract overlapping transcripts ---
  txs <- GenomicFeatures::transcriptsByOverlaps(txdb, region_gr)

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
      list(transcripts = transcripts_df, exons = exons_df, region = region_gr),
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

  # --- 8. Try to get gene names ---
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
    # Match gene IDs back to transcript order
    idx <- match(as.character(txs$tx_name), gene_map$TXNAME)
    mapped_gene <- gene_map$GENEID[idx]
    gene_name_vec <- ifelse(is.na(mapped_gene), as.character(txs$tx_name), mapped_gene)
  }

  # --- 9. Build transcripts data.frame ---
  transcripts_df <- data.frame(
    tx_id     = as.integer(txs$tx_id),
    tx_name   = as.character(txs$tx_name),
    gene_name = gene_name_vec,
    strand    = as.character(BiocGenerics::strand(txs)),
    tx_start  = GenomicRanges::start(txs),
    tx_end    = GenomicRanges::end(txs),
    stringsAsFactors = FALSE
  )

  # --- 10. Return S3 object ---
  structure(
    list(
      transcripts = transcripts_df,
      exons       = exons_df,
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
