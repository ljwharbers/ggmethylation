# Internal function — not exported

#' Build ggplot2 layers for SNV marks on reads
#'
#' Takes pre-extracted per-read base calls (from [extract_variant_bases()]) and
#' returns a list of ggplot2 layer objects that draw shape-encoded points at the
#' position/lane of each non-ref call. The layers can be added directly to a
#' ggplot with `+`.
#'
#' @param variant_bases A data.frame with columns `read_name`, `position`,
#'   `lane`, and `variant_class` (values: `"ref"`, `"alt"`, `"other"`,
#'   `"del"`), as returned by [extract_variant_bases()]. May be `NULL` or
#'   zero-row.
#'
#' @return A list with one ggplot2 `geom_point` layer that can be appended to a
#'   ggplot with `+`, or `NULL` when there are no non-ref calls to draw. All
#'   variant classes are drawn as a small red asterisk (shape 8).
#'
#' @keywords internal
build_snv_layer <- function(variant_bases) {
  # Return NULL for empty/NULL input
  if (is.null(variant_bases) || nrow(variant_bases) == 0L) return(NULL)

  # Filter to non-ref only
  non_ref <- variant_bases[variant_bases$variant_class != "ref", , drop = FALSE]

  if (nrow(non_ref) == 0L) return(NULL)

  list(
    ggplot2::geom_point(
      data = non_ref,
      ggplot2::aes(
        x = .data$position,
        y = .data$lane
      ),
      colour      = "#000000",
      shape       = 8L,
      size        = 0.8,
      stroke      = 0.3,
      inherit.aes = FALSE
    )
  )
}
