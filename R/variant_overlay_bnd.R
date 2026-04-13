# Internal functions — not exported

#' Flag SA-tagged reads whose breakpoint matches a VCF BND call
#'
#' For each read carrying an SA (supplementary alignment) tag, checks whether
#' the supplementary chromosome and position match any BND record in `bnd_df`
#' within a position tolerance. Both the chromosome and position conditions must
#' be satisfied simultaneously.
#'
#' @param reads Data.frame — the `$reads` element of a `methylation_data`
#'   object. Relevant columns: `read_name` (chr), `sa_chrom` (chr, `NA` for
#'   non-SA reads), `sa_pos` (int, `NA` for non-SA reads).
#' @param bnd_df Data.frame — BND rows from `variant_data$variants`, with
#'   columns `position` (int), `mate_chrom` (chr), `mate_pos` (int). May be
#'   `NULL` or zero-row.
#' @param tol Integer, default `50L` — position tolerance in bp. A read's SA
#'   position is considered matching when `abs(bnd_df$position - sa_pos) <= tol`.
#'
#' @return `reads` with one new logical column `vcf_validated`:
#'   - `TRUE`  if `sa_chrom == bnd_df$mate_chrom` AND
#'             `abs(bnd_df$position - sa_pos) <= tol` for ANY BND row.
#'   - `FALSE` otherwise (including reads with no SA tag, i.e., `is.na(sa_chrom)`).
#'
#' @keywords internal
match_sa_to_vcf_bnd <- function(reads, bnd_df, tol = 50L) {
  # Default: all reads are not validated
  reads$vcf_validated <- FALSE

  # Early exit when there are no BND records to match against
  if (is.null(bnd_df) || nrow(bnd_df) == 0L) {
    return(reads)
  }

  # Only consider reads that have an SA tag (non-NA sa_chrom and sa_pos)
  has_sa <- !is.na(reads$sa_chrom) & !is.na(reads$sa_pos)

  if (!any(has_sa)) {
    return(reads)
  }

  # For each SA-carrying read, check whether ANY BND row satisfies both
  # the chromosome and position conditions simultaneously.
  sa_indices <- which(has_sa)

  for (i in sa_indices) {
    r_chrom <- reads$sa_chrom[i]
    r_pos   <- reads$sa_pos[i]

    chrom_match <- bnd_df$mate_chrom == r_chrom
    pos_match   <- abs(bnd_df$position - r_pos) <= tol

    if (any(chrom_match & pos_match)) {
      reads$vcf_validated[i] <- TRUE
    }
  }

  reads
}


#' Build ggplot2 layers marking VCF BND positions on the read panel
#'
#' Returns a flat list of ggplot2 layer objects (a vertical line and a text
#' label) for each BND position. The vertical line marks the BND call;
#' the label shows the mate chromosome using the join symbol (⋈).
#'
#' Because the read panel uses `scale_y_reverse()`, `y = -Inf` places the
#' text at the *top* of the panel (the intended position).
#'
#' @param bnd_df Data.frame — BND rows from `variant_data$variants` (columns:
#'   `position` (int), `mate_chrom` (chr), `mate_pos` (int)). May be `NULL`
#'   or zero-row.
#'
#' @return A flat list with two ggplot2 layer objects per BND
#'   (`geom_vline` + `geom_text`), or `NULL` when `bnd_df` is `NULL` or
#'   zero-row.
#'
#' @keywords internal
build_bnd_layer <- function(bnd_df) {
  if (is.null(bnd_df) || nrow(bnd_df) == 0L) return(NULL)

  list(
    ggplot2::geom_vline(
      xintercept = bnd_df$position,
      linetype   = "solid",
      colour     = "#880E4F",
      linewidth  = 0.6,
      alpha      = 0.8
    ),
    ggplot2::geom_text(
      data        = bnd_df,
      mapping     = ggplot2::aes(
        x     = .data$position,
        label = paste0("\u22c8 ", .data$mate_chrom)
      ),
      y           = -Inf,
      vjust       = -0.3,
      hjust       = 0.5,
      size        = 1.8,
      colour      = "#880E4F",
      inherit.aes = FALSE
    )
  )
}
