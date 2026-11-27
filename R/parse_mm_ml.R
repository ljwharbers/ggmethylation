# MM/ML tag parser for ggmethylation

#' Parse MM and ML tags for a single read
#'
#' Extracts modification probabilities and their genomic positions from the
#' MM and ML tags of a single BAM read, following the SAM specification for
#' base modification tags (SAMtags, section "Base modifications").
#'
#' Bases are counted along the original read, 5' to 3', on its top strand.
#' Which base is looked for in `seq` and the direction of the scan depend
#' only on the alignment strand: for a reverse-strand alignment (`strand =
#' "-"`, FLAG 0x10) the stored SEQ is reverse-complemented, so the complement
#' of the MM base is searched from the end of `seq`. The `+`/`-` field of an
#' MM entry only says which strand of the molecule carries the modification;
#' it does not change which bases are counted. Canonical base `N` matches any
#' base.
#'
#' An entry may carry several codes (e.g. `C+mh`); its ML values are then
#' interleaved per base in code order. A numeric code (a ChEBI identifier
#' such as `C+76792`) is one code. When several entries carry the requested
#' code (e.g. `C+m` and `G-m`), calls from all of them are returned.
#'
#' The `.` flag means unlisted canonical bases are unmodified: they are
#' returned with `mod_prob = 0`. The `?` flag means unlisted bases carry no
#' information and are not returned. An entry with **no** flag is treated
#' like `?`. This deliberately deviates from the SAM specification (which
#' treats a missing flag as `.`): older basecallers that write no flag list
#' only CpG cytosines, so following the specification would report every
#' non-CpG C as unmodified.
#'
#' A tag that cannot be interpreted is rejected as a whole, with a warning,
#' and no calls are returned: an unparseable entry, a missing or negative
#' skip count, an ML length that differs from the number of calls in MM, an
#' entry that lists more bases than SEQ contains, an empty SEQ (`*`), an
#' `mn` value that differs from the length of `seq` (a hard-clipped record
#' whose MM/ML still describe the full read), or a hard-clipped record (an
#' `H` operation in `cigar`) without `mn`. In the last case there is no way
#' to tell whether MM/ML were rewritten for the clipped SEQ, and if they were
#' not, every call would be placed on the wrong base. When `mn` is given and
#' equals the length of `seq`, a hard-clipped record is trusted.
#'
#' @param seq Character. The read sequence from the BAM SEQ field.
#' @param mm_tag Character. The MM tag value, e.g. `"C+m,0,1,3;C+h,2;"`.
#' @param ml_tag Integer vector. ML tag values (0--255), concatenated for all
#'   modifications listed in the MM tag.
#' @param mod_code Character. The modification code to extract (e.g., `"m"`
#'   for 5mC, `"h"` for 5hmC).
#' @param strand Character. `"+"` or `"-"`, the strand of the alignment.
#' @param cigar Character. The CIGAR string for the alignment.
#' @param pos Integer. 1-based leftmost alignment position (BAM POS field).
#' @param mn Integer or `NULL`. The value of the optional `MN` tag (the SEQ
#'   length the MM/ML tags were written for). `NULL` or `NA` means the tag is
#'   absent: a record without hard clips is then trusted, a hard-clipped one
#'   is rejected.
#'
#' @return A named `list` with two `data.frame` elements:
#'   \describe{
#'     \item{sites}{Reference-aligned modifications. Columns `position`
#'       (integer, 1-based genomic position) and `mod_prob` (numeric,
#'       ML value / 255). When the MM entry uses the `.` implicit-unmodified
#'       flag, also includes rows for unlisted canonical positions with
#'       `mod_prob = 0`.}
#'     \item{insertion_sites}{Modifications whose query base falls inside
#'       a CIGAR `I` interval (no reference position). Columns `query_pos`
#'       (integer, 1-based position in the read sequence) and `mod_prob`
#'       (numeric). Modifications inside soft/hard clips are still dropped.
#'       Implicit-zero positions (`.` flag) inside insertions are not emitted.}
#'   }
#'   Returns both frames empty (zero rows, same columns) when the
#'   requested modification is not present or the tags are malformed.
#'
#' @keywords internal
parse_mm_ml = function(seq, mm_tag, ml_tag, mod_code, strand, cigar, pos, mn = NULL) {
  parsed = tryCatch(
    .parse_mods(seq, mm_tag, ml_tag, mod_code, strand, cigar, pos, mn = mn),
    ggmethylation_malformed_mm = function(e) {
      warning("Malformed MM/ML tags, no calls returned: ", conditionMessage(e),
              call. = FALSE)
      NULL
    }
  )
  if (is.null(parsed)) parsed = .empty_mods()
  sites = parsed$sites[, c("position", "mod_prob"), drop = FALSE]
  ins   = parsed$insertion_sites[, c("query_pos", "mod_prob"), drop = FALSE]
  rownames(sites) = NULL
  rownames(ins)   = NULL
  list(sites = sites, insertion_sites = ins)
}

# Zero-row result of .parse_mods().
.empty_mods = function() {
  list(
    sites = data.frame(position = integer(0), mod_prob = numeric(0),
                       mod_code = character(0), stringsAsFactors = FALSE),
    insertion_sites = data.frame(query_pos = integer(0), mod_prob = numeric(0),
                                 mod_code = character(0), stringsAsFactors = FALSE)
  )
}

# Signal a malformed or inconsistent MM/ML tag. Callers catch the
# `ggmethylation_malformed_mm` class: parse_mm_ml() turns it into a warning,
# .parse_read_mods() counts it and warns once for all reads.
.mm_malformed = function(...) {
  stop(structure(
    class = c("ggmethylation_malformed_mm", "error", "condition"),
    list(message = paste0(...), call = NULL)
  ))
}

# Split an MM tag into entries. Each entry is a list with `base` (canonical
# base, upper case), `mod_strand` ("+"/"-"), `codes` (character vector: one
# element per single-letter code, or one ChEBI number), `flag` ("", "?" or
# "."), `deltas` (integer skip counts) and `n_ml` (ML values the entry uses:
# one per code per listed base). Signals a malformed tag on any entry it
# cannot parse.
.parse_mm_tag = function(mm_tag) {
  entries = trimws(strsplit(mm_tag, ";", fixed = TRUE)[[1L]])
  entries = entries[nzchar(entries)]
  lapply(entries, function(entry) {
    parts = trimws(strsplit(entry, ",", fixed = TRUE)[[1L]])
    # strsplit() drops a trailing empty field; keep it so "C+m,0," is caught.
    if (grepl(",$", entry)) parts = c(parts, "")
    m = regmatches(parts[1L], regexec("^([ACGTUN])([+-])([A-Za-z]+|[0-9]+)([.?]?)$",
                                      parts[1L]))[[1L]]
    if (length(m) == 0L) .mm_malformed("cannot parse MM entry '", entry, "'")
    delta_str = parts[-1L]
    if (!all(grepl("^[0-9]+$", delta_str))) {
      .mm_malformed("MM entry '", entry, "' has a missing or non-numeric skip count (delta)")
    }
    deltas = as.integer(delta_str)
    if (anyNA(deltas)) .mm_malformed("MM entry '", entry, "' has an out-of-range delta")
    codes = if (grepl("^[0-9]+$", m[4L])) m[4L] else strsplit(m[4L], "")[[1L]]
    list(entry = entry, base = m[2L], mod_strand = m[3L], codes = codes, flag = m[5L],
         deltas = deltas, n_ml = length(deltas) * length(codes))
  })
}

# Reference position and insertion operation for every query position of a
# read, from a decompose_cigar() table. `ref` is NA outside M/=/X; `ins_op` is
# the row of the `I` operation a query position sits in (NA elsewhere).
.query_ref_map = function(dc, seq_len) {
  qlen = max(c(seq_len, dc$query_end), na.rm = TRUE)
  ref = ins_op = rep(NA_integer_, qlen)
  for (k in which(dc$type %in% c("M", "=", "X"))) {
    ref[dc$query_start[k]:dc$query_end[k]] = dc$ref_start[k]:dc$ref_end[k]
  }
  for (k in which(dc$type == "I")) {
    ins_op[dc$query_start[k]:dc$query_end[k]] = k
  }
  list(ref = ref, ins_op = ins_op)
}

# Core MM/ML parser: parse the tag, sequence and CIGAR of one read once and
# extract every code in `mod_codes`. Returns list(sites = data.frame(position,
# mod_prob, mod_code), insertion_sites = data.frame(query_pos, mod_prob,
# mod_code)), rows grouped by code in `mod_codes` order. Signals
# `ggmethylation_malformed_mm` when the tags cannot be interpreted (see
# parse_mm_ml()). `dc` may pass a precomputed decompose_cigar(cigar, pos).
.parse_mods = function(seq, mm_tag, ml_tag, mod_codes, strand, cigar, pos,
                       mn = NULL, dc = NULL) {
  if (is.null(mm_tag) || length(mm_tag) == 0L || is.na(mm_tag) || !nzchar(mm_tag)) {
    return(.empty_mods())
  }
  entries = .parse_mm_tag(mm_tag)
  wanted = vapply(entries, function(e) any(e$codes %in% mod_codes), logical(1L))
  if (!any(wanted)) return(.empty_mods())

  # --- Consistency of MM with ML, MN and SEQ ---
  if (is.null(ml_tag) || all(is.na(ml_tag))) ml_tag = integer(0)
  n_ml = sum(vapply(entries, `[[`, integer(1L), "n_ml"))
  if (length(ml_tag) != n_ml) {
    .mm_malformed("ML has ", length(ml_tag), " values but the MM tag lists ", n_ml, " calls")
  }
  seq = toupper(as.character(seq))
  if (identical(seq, "*")) seq = ""
  seq_len = nchar(seq)
  has_mn = !is.null(mn) && length(mn) == 1L && !is.na(mn)
  if (has_mn && mn != seq_len) {
    .mm_malformed("MN tag (", mn, ") does not match the SEQ length (", seq_len, ")")
  }
  # A hard clip removed bases from SEQ. Without MN there is no way to tell
  # whether MM/ML were rewritten for the clipped SEQ or still describe the
  # full read (then every call would land on the wrong base), so the record
  # is not trusted.
  if (!has_mn && isTRUE(grepl("H", cigar, fixed = TRUE))) {
    .mm_malformed("CIGAR ", cigar, " is hard-clipped and there is no MN tag to ",
                  "confirm that MM/ML describe the clipped SEQ")
  }
  if (seq_len == 0L) {
    if (any(vapply(entries[wanted], function(e) length(e$deltas) > 0L, logical(1L)))) {
      .mm_malformed("SEQ is empty ('*') but the MM tag lists calls")
    }
    return(.empty_mods())
  }

  # --- Sequence and CIGAR context, computed once ---
  seq_chars = strsplit(seq, "")[[1L]]
  reverse = identical(strand, "-")
  complement = c(A = "T", C = "G", G = "C", T = "A", U = "A")
  base_positions = function(base) {
    # 1-based SEQ positions of `base`, in the order of the original read.
    hits = if (base == "N") {
      seq_along(seq_chars)
    } else {
      which(seq_chars == if (reverse) complement[[base]] else if (base == "U") "T" else base)
    }
    if (reverse) rev(hits) else hits
  }
  if (is.null(dc)) dc = decompose_cigar(cigar, pos)
  qmap = .query_ref_map(dc, seq_len)

  # --- Walk each entry; collect calls per code ---
  site_pieces = ins_pieces = list()
  ml_offset = 0L
  for (k in seq_along(entries)) {
    e = entries[[k]]
    if (wanted[k]) {
      canon = base_positions(e$base)
      listed = cumsum(e$deltas + 1L)
      if (length(listed) > 0L && listed[length(listed)] > length(canon)) {
        .mm_malformed("MM entry '", e$entry, "' lists more ", e$base,
                      " bases than SEQ contains (", length(canon), ")")
      }
      q_listed = canon[listed]
      ref_listed = qmap$ref[q_listed]
      in_ins = is.na(ref_listed) & !is.na(qmap$ins_op[q_listed])
      q_zero = if (e$flag == ".") canon[-listed] else integer(0)
      if (length(listed) == 0L && e$flag == ".") q_zero = canon
      ref_zero = qmap$ref[q_zero]
      ref_zero = ref_zero[!is.na(ref_zero)]

      for (ci in which(e$codes %in% mod_codes)) {
        code = e$codes[ci]
        ml = ml_tag[ml_offset + (seq_along(listed) - 1L) * length(e$codes) + ci]
        position = c(ref_listed[!is.na(ref_listed)], ref_zero)
        site_pieces[[length(site_pieces) + 1L]] = data.frame(
          position = position,
          mod_prob = c(ml[!is.na(ref_listed)] / 255, rep(0, length(ref_zero))),
          mod_code = rep(code, length(position)),
          stringsAsFactors = FALSE
        )
        ins_pieces[[length(ins_pieces) + 1L]] = data.frame(
          query_pos = q_listed[in_ins],
          mod_prob  = ml[in_ins] / 255,
          mod_code  = rep(code, sum(in_ins)),
          stringsAsFactors = FALSE
        )
      }
    }
    ml_offset = ml_offset + e$n_ml
  }

  empty = .empty_mods()
  bind = function(pieces, template) {
    pieces = Filter(function(p) nrow(p) > 0L, pieces)
    if (length(pieces) == 0L) return(template)
    out = do.call(rbind, pieces)
    out = out[order(match(out$mod_code, mod_codes)), , drop = FALSE]
    rownames(out) = NULL
    out
  }
  list(sites = bind(site_pieces, empty$sites),
       insertion_sites = bind(ins_pieces, empty$insertion_sites))
}
