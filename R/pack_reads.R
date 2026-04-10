#' Pack reads into horizontal lanes for compact visualization
#'
#' Assigns reads to horizontal lanes (rows) using greedy interval scheduling,
#' similar to how a genome browser packs aligned reads. Each lane contains
#' non-overlapping reads separated by at least `gap` base pairs.
#'
#' @param reads A data.frame with columns `read_name`, `start`, and `end`.
#'   Rows should already be sorted in the desired display order.
#' @param gap Minimum gap in base pairs between reads on the same lane
#'   (default 10).
#' @param clip_side Character vector (same length as `nrow(reads)`) indicating
#'   soft-clip status of each read: `"left"`, `"right"`, `"both"`, or `NA`.
#'   When non-NULL, adjacent clipped reads receive a wider gap so that
#'   supplementary-alignment indicators do not overlap. Default `NULL` (no
#'   clip-aware packing).
#' @param clip_gap Minimum gap in base pairs between reads when at least one
#'   side is clipped (default 100).
#'
#' @return An integer vector of lane assignments (same length as
#'   `nrow(reads)`), where lane 1 is the top.
#'
#' @keywords internal
pack_reads <- function(reads, gap = 10, clip_side = NULL, clip_gap = 100L) {
  n <- nrow(reads)
  if (n == 0L) {
    return(integer(0L))
  }

  lanes <- list()
  result <- integer(n)

  for (i in seq_len(n)) {
    placed <- FALSE
    for (j in seq_along(lanes)) {
      # Determine the effective gap for this pair of reads
      if (!is.null(clip_side)) {
        prev_clip <- lanes[[j]]$clip
        curr_clip <- clip_side[i]
        if ((!is.na(prev_clip) && prev_clip %in% c("right", "both")) ||
            (!is.na(curr_clip) && curr_clip %in% c("left", "both"))) {
          effective_gap <- clip_gap
        } else {
          effective_gap <- gap
        }
      } else {
        effective_gap <- gap
      }

      if (reads$start[i] > lanes[[j]]$end + effective_gap) {
        result[i] <- j
        lanes[[j]] <- list(
          end  = reads$end[i],
          clip = if (!is.null(clip_side)) clip_side[i] else NA
        )
        placed <- TRUE
        break
      }
    }
    if (!placed) {
      lanes <- c(lanes, list(list(
        end  = reads$end[i],
        clip = if (!is.null(clip_side)) clip_side[i] else NA
      )))
      result[i] <- length(lanes)
    }
  }

  result
}
