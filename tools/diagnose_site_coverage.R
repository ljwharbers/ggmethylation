#!/usr/bin/env Rscript
# Diagnose per-read methylation site coverage.
#
# Investigates why different reads have different sets of called sites at
# the same reference positions. Classifies each (read, position) pair into
# one of five categories and prints a summary table.
#
# Usage (from the repository root):
#   Rscript tools/diagnose_site_coverage.R [bam] [region]
#
# Defaults to the standard test BAM and a 3 kb region on chr22.

library(ggmethylation)
library(Rsamtools)

# ---- Configuration -----------------------------------------------------------

BAM_PATH <- "/staging/leuven/stg_00096/home/lharbers/repositories/ggmethylation/data/PTCL8_PB_tumor_chr21_22_subset.bam"
REGION   <- "chr22:23500000-23503000"
MOD_CODES <- c("m")  # extend with "a" etc. if present

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1) BAM_PATH <- args[1]
if (length(args) >= 2) REGION   <- args[2]

if (!file.exists(BAM_PATH)) {
  stop("BAM not found: ", BAM_PATH, "\n",
       "Run from the repo root or pass a path as the first argument.")
}

# ---- 1. Detect MM tag flags --------------------------------------------------

cat("=== MM tag flag distribution ===\n")

gr <- ggmethylation:::region_to_granges(REGION)
param_raw <- ScanBamParam(
  which = gr,
  what  = c("qname", "strand"),
  tag   = c("MM")
)
raw <- scanBam(BAM_PATH, param = param_raw)[[1]]
mm_tags <- raw$tag$MM
mm_tags  <- mm_tags[!is.na(mm_tags) & nzchar(mm_tags)]

# Parse all MM entries across all reads, extract per-(mod_code, flag) counts
parse_mm_flags <- function(mm_tag) {
  mm_clean <- sub(";$", "", mm_tag)
  entries  <- strsplit(mm_clean, ";")[[1]]
  lapply(entries, function(e) {
    parts <- strsplit(e, ",")[[1]]
    spec  <- parts[1]
    last  <- substr(spec, nchar(spec), nchar(spec))
    flag  <- if (last %in% c("?", ".")) last else "(none)"
    if (flag != "(none)") spec <- substr(spec, 1, nchar(spec) - 1)
    code  <- substring(spec, 3)
    list(code = code, flag = flag)
  })
}

flag_rows <- do.call(rbind, lapply(mm_tags, function(mm) {
  entries <- parse_mm_flags(mm)
  data.frame(
    mod_code = vapply(entries, `[[`, character(1), "code"),
    flag     = vapply(entries, `[[`, character(1), "flag"),
    stringsAsFactors = FALSE
  )
}))

flag_summary <- aggregate(n ~ mod_code + flag,
                          data = cbind(flag_rows, n = 1L),
                          FUN  = sum)
flag_summary <- flag_summary[order(flag_summary$mod_code, flag_summary$flag), ]
print(flag_summary, row.names = FALSE)
cat("\n")

# ---- 2. Load methylation data ------------------------------------------------

cat("=== Loading methylation data ===\n")
md_list <- lapply(MOD_CODES, function(mc) {
  tryCatch(
    read_methylation(BAM_PATH, REGION, mod_code = mc, max_reads = 50L),
    error = function(e) { message("Skipping mod_code '", mc, "': ", e$message); NULL }
  )
})
names(md_list) <- MOD_CODES
md_list <- Filter(Negate(is.null), md_list)

if (length(md_list) == 0L) {
  stop("No methylation data could be loaded. Check BAM and region.")
}

# ---- 3. Classify (read, ref_position) pairs ----------------------------------
# For each position covered by >=2 reads, classify every overlapping read.
# Categories:
#   call              — position present in $sites for this read
#   canonical-no-call — read has canonical base at this ref pos but no $sites entry
#   non-canonical     — read has a non-canonical (SNV/error) base at this ref pos
#   inside-insertion  — position is inside a CIGAR I interval for this read
#   deletion-outside  — position inside D/N or outside the read's alignment span

classify_read_pos <- function(md) {
  reads     <- md$reads
  sites     <- md$sites
  seqs      <- md$sequences
  cigars    <- md$cigars
  cf        <- md$cigar_features
  ins_sites <- md$insertion_sites

  # Collect all ref positions that appear in $sites
  all_pos  <- sort(unique(sites$position))

  # Only keep positions covered by >= 2 reads
  pos_by_read <- lapply(seq_len(nrow(reads)), function(i) {
    rname  <- reads$read_name[i]
    rstart <- reads$start[i]
    rend   <- reads$end[i]
    all_pos[all_pos >= rstart & all_pos <= rend]
  })
  pos_counts  <- table(unlist(pos_by_read))
  multi_pos   <- as.integer(names(pos_counts[pos_counts >= 2L]))

  if (length(multi_pos) == 0L) {
    message("No reference positions covered by >=2 reads in this region/dataset.")
    return(data.frame(category = character(0), n = integer(0),
                      stringsAsFactors = FALSE))
  }

  # Determine canonical base per mod_code  (simple: C for m/h, A for a)
  canonical_for_code <- function(mc) {
    switch(mc, m = "C", h = "C", a = "A", f = "A", "C")
  }
  canon_base <- canonical_for_code(md$mod_code[1])

  rows <- lapply(seq_len(nrow(reads)), function(i) {
    rname  <- reads$read_name[i]
    rstart <- reads$start[i]
    rend   <- reads$end[i]
    cigar  <- cigars[rname]
    bam_pos <- reads$bam_pos[i]
    seq_str <- seqs[rname]

    # Positions to classify for this read
    pos_here <- multi_pos[multi_pos >= rstart & multi_pos <= rend]
    if (length(pos_here) == 0L) return(NULL)

    # Sites present in $sites for this read
    called_pos <- sites$position[sites$read_name == rname]

    # Insertions for this read (ref_start = anchor before insertion)
    ins_anchors <- if (nrow(cf) > 0L) {
      cf_read <- cf[cf$read_name == rname & cf$type == "I", , drop = FALSE]
      if (nrow(cf_read) > 0L) cf_read$ref_start else integer(0)
    } else integer(0)

    # Deletions / skipped regions for this read
    del_rows <- if (nrow(cf) > 0L) {
      cf[cf$read_name == rname & cf$type %in% c("D", "N"), , drop = FALSE]
    } else cf[integer(0), , drop = FALSE]

    vapply(pos_here, function(p) {
      if (p %in% called_pos) return("call")

      # Check if inside a deletion or skipped region
      if (nrow(del_rows) > 0L &&
          any(p >= del_rows$ref_start & p <= del_rows$ref_end, na.rm = TRUE))
        return("deletion-outside")

      # Check if position is an insertion anchor (ref pos just before insertion)
      if (p %in% ins_anchors) return("inside-insertion")

      # Map ref position to query position to inspect the base
      q_pos <- ggmethylation:::ref_to_seq(cigar, bam_pos, p)
      if (is.na(q_pos) || q_pos < 1L || q_pos > nchar(seq_str))
        return("deletion-outside")

      read_base <- toupper(substr(seq_str, q_pos, q_pos))
      if (read_base == canon_base) "canonical-no-call" else "non-canonical"
    }, character(1L))
  })

  all_cats <- unlist(rows)
  if (length(all_cats) == 0L) return(data.frame(category = character(0),
                                                  n = integer(0),
                                                  stringsAsFactors = FALSE))
  tab <- as.data.frame(table(category = all_cats), stringsAsFactors = FALSE)
  names(tab)[2] <- "n"
  tab$pct <- round(100 * tab$n / sum(tab$n), 1)
  tab[order(-tab$n), ]
}

# ---- 4. Print results --------------------------------------------------------

all_results <- lapply(names(md_list), function(mc) {
  cat(sprintf("=== mod_code '%s': site-coverage classification ===\n", mc))
  tab <- classify_read_pos(md_list[[mc]])
  print(tab, row.names = FALSE)
  cat(sprintf("Total (read x position) pairs assessed: %d\n\n",
              sum(tab$n)))
  cbind(mod_code = mc, tab)
})

combined <- do.call(rbind, all_results)

# ---- 5. Save TSV -------------------------------------------------------------

out_tsv <- file.path("tools", "site_coverage_classification.tsv")
write.table(combined, out_tsv, sep = "\t", row.names = FALSE, quote = FALSE)
cat("Per-category table saved to:", out_tsv, "\n")

# ---- 6. Interpretation guide -------------------------------------------------

cat("\nInterpretation:\n")
cat("  call              — position present in $sites (basecaller listed it)\n")
cat("  canonical-no-call — read has the canonical base but basecaller did not list it\n")
cat("                      (basecaller omission; fix: basecaller should use '.' flag)\n")
cat("  non-canonical     — read has a different base at this position (SNV/error)\n")
cat("  inside-insertion  — position is inside a CIGAR I interval for this read\n")
cat("                      (present in $insertion_sites, not $sites — by design)\n")
cat("  deletion-outside  — position inside D/N or outside the read's alignment span\n")
cat("                      (read genuinely has no base here — not a bug)\n")
