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
#'
#' @return An integer vector of lane assignments (same length as
#'   `nrow(reads)`), where lane 1 is the top.
#'
#' @keywords internal
pack_reads <- function(reads, gap = 10) {
  n <- nrow(reads)
  if (n == 0L) {
    return(integer(0L))
  }

  lanes <- list()
  result <- integer(n)

  for (i in seq_len(n)) {
    placed <- FALSE
    for (j in seq_along(lanes)) {
      if (reads$start[i] > lanes[[j]] + gap) {
        result[i] <- j
        lanes[[j]] <- reads$end[i]
        placed <- TRUE
        break
      }
    }
    if (!placed) {
      lanes <- c(lanes, list(reads$end[i]))
      result[i] <- length(lanes)
    }
  }

  result
}
