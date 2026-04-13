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
#' @return A list of four ggplot2 layer objects
#'   (`ggnewscale::new_scale_colour()`, `geom_point`, `scale_colour_manual`,
#'   `scale_shape_manual`) that can be appended to a ggplot with `+`, or
#'   `NULL` when there are no non-ref calls to draw.
#'
#' @keywords internal
build_snv_layer <- function(variant_bases) {
  # Return NULL for empty/NULL input
  if (is.null(variant_bases) || nrow(variant_bases) == 0L) return(NULL)

  # Filter to non-ref only
  non_ref <- variant_bases[variant_bases$variant_class != "ref", , drop = FALSE]

  if (nrow(non_ref) == 0L) return(NULL)

  list(
    ggnewscale::new_scale_colour(),
    ggplot2::geom_point(
      data = non_ref,
      ggplot2::aes(
        x      = .data$position,
        y      = .data$lane,
        colour = .data$variant_class,
        shape  = .data$variant_class
      ),
      size        = 1.5,
      inherit.aes = FALSE
    ),
    ggplot2::scale_colour_manual(
      values = c(alt = "#E53935", other = "#FDD835", del = "#212121"),
      # #212121 = soft black, consistent with existing variant colour palette in build_read_panel.R
      guide  = "none"
    ),
    ggplot2::scale_shape_manual(
      values = c(alt = 19L, other = 18L, del = 1L),
      guide  = "none"
    )
  )
}
