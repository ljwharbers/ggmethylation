# Build the synthetic modBAM + VCF used by tests/testthat/fixtures/.
#
# ~70 reads on a random 8 kb chr1 carrying HP tags (some reads untagged),
# 5mC (m) and 5hmC (h) MM/ML calls, a 150 bp deletion shared by HP1 reads,
# a ~120 bp insertion locus, soft clips with SA tags, low-MAPQ reads, one
# supplementary alignment sharing a primary's name, and an SNV at chr1:3000
# carried by HP2 reads. The VCF holds that SNV, a <DEL> and a BND.
#
# Regenerate (needs samtools, bgzip and tabix on PATH):
#   Rscript data-raw/make_test_fixture.R $TMPDIR/fx
#   samtools sort -o tests/testthat/fixtures/synthetic.bam $TMPDIR/fx/reads.sam
#   samtools index tests/testthat/fixtures/synthetic.bam
#   bgzip -c $TMPDIR/fx/vars.vcf > tests/testthat/fixtures/synthetic.vcf.gz
#   tabix -p vcf tests/testthat/fixtures/synthetic.vcf.gz
set.seed(42)
out = commandArgs(trailingOnly = TRUE)[1]
dir.create(out, showWarnings = FALSE, recursive = TRUE)

ref_len = 8000L
bases = c("A", "C", "G", "T")
ref = sample(bases, ref_len, replace = TRUE, prob = c(0.3, 0.2, 0.2, 0.3))
# sprinkle CpGs
cg = sort(sample(seq(10, ref_len - 10), 500))
ref[cg] = "C"; ref[cg + 1L] = "G"
ref[3000] = "A"  # SNV site
ref_str = paste(ref, collapse = "")
writeLines(c(">chr1", ref_str), file.path(out, "ref.fa"))

rand_seq = function(n) sample(bases, n, replace = TRUE, prob = c(0.25, 0.3, 0.2, 0.25))
comp = c(A = "T", C = "G", G = "C", T = "A")

prob_at = function(refpos, grp, code) {
  base = if (identical(grp, "1")) 0.8 else if (identical(grp, "2")) 0.2 else 0.5
  if (!is.na(refpos) && refpos > 4500 && refpos < 5500) base = 1 - base
  if (code == "h") base = base / 3
  p = if (is.na(refpos)) 0.9 else base + stats::rnorm(1, 0, 0.15)
  min(max(p, 0), 1)
}

make_read = function(i, grp, strand, start, len, del = NULL, ins = NULL,
                     lclip = 0L, rclip = 0L, sa = NULL, flag_extra = 0L,
                     name = sprintf("read%03d", i), mapq = 60L, snv_alt = FALSE) {
  # Walk building query seq and cigar
  ops = character(0); seq = character(0); qref = integer(0)
  if (lclip > 0) { ops = c(ops, sprintf("%dS", lclip)); s = rand_seq(lclip); seq = c(seq, s); qref = c(qref, rep(NA, lclip)) }
  pos = start; end = start + len - 1L
  breaks = list()
  if (!is.null(del)) breaks[[length(breaks) + 1]] = list(type = "D", at = del[1], n = del[2])
  if (!is.null(ins)) breaks[[length(breaks) + 1]] = list(type = "I", at = ins[1], n = ins[2])
  breaks = breaks[order(vapply(breaks, `[[`, numeric(1), "at"))]
  for (b in breaks) {
    m = b$at - pos
    if (m > 0) {
      ops = c(ops, sprintf("%dM", m)); seq = c(seq, ref[pos:(b$at - 1L)]); qref = c(qref, pos:(b$at - 1L)); pos = b$at
    }
    if (b$type == "D") { ops = c(ops, sprintf("%dD", b$n)); pos = pos + b$n }
    else { ops = c(ops, sprintf("%dI", b$n)); s = rand_seq(b$n); seq = c(seq, s); qref = c(qref, rep(NA, b$n)) }
  }
  m = end - pos + 1L
  ops = c(ops, sprintf("%dM", m)); seq = c(seq, ref[pos:end]); qref = c(qref, pos:end)
  if (rclip > 0) { ops = c(ops, sprintf("%dS", rclip)); seq = c(seq, rand_seq(rclip)); qref = c(qref, rep(NA, rclip)) }
  if (snv_alt) seq[which(qref == 3000L)] = "G"

  # MM/ML in original read orientation
  if (strand == "+") {
    cpos = which(seq == "C")
  } else {
    cpos = rev(which(seq == "G"))
  }
  mm = character(0); ml = integer(0)
  for (code in c("m", "h")) {
    probs = vapply(cpos, function(k) prob_at(qref[k], grp, code), numeric(1))
    mm = c(mm, paste0("C+", code, "?", paste0(",", rep(0L, length(cpos)), collapse = "")))
    ml = c(ml, as.integer(round(probs * 255)))
  }
  mm_tag = paste0(paste(mm, collapse = ";"), ";")
  flag = (if (strand == "-") 16L else 0L) + flag_extra
  tags = c(sprintf("MM:Z:%s", mm_tag), sprintf("ML:B:C,%s", paste(ml, collapse = ",")))
  if (!is.null(grp)) tags = c(tags, sprintf("HP:i:%s", grp))
  if (!is.null(sa)) tags = c(tags, sprintf("SA:Z:%s", sa))
  paste(c(name, flag, "chr1", start, mapq, paste(ops, collapse = ""), "*", 0, 0,
          paste(seq, collapse = ""), "*", tags), collapse = "\t")
}

lines = character(0)
i = 0L
for (k in 1:70) {
  i = i + 1L
  grp = if (k %% 7 == 0) NULL else as.character((k %% 2) + 1)
  strand = if (runif(1) < 0.5) "+" else "-"
  start = sample(200:5000, 1); len = sample(800:2800, 1)
  if (start + len > ref_len - 300) len = ref_len - 300 - start
  del = NULL; ins = NULL; lclip = 0L; rclip = 0L; sa = NULL
  if (identical(grp, "1") && start < 3800 && start + len > 4300) del = c(3900, 150)
  if (start < 6000 && start + len > 6200 && k %% 3 != 0) ins = c(6100, 120 + sample(-8:8, 1))
  if (k %% 11 == 0) { rclip = 400L; sa = "chr7,5000,+,400M,60,0;" }
  if (k %% 13 == 0) { lclip = 300L; sa = "chr3,12000,-,300M,60,0;" }
  mapq = if (k %% 9 == 0) 5L else 60L
  snv_alt = identical(grp, "2")
  lines = c(lines, make_read(i, grp, strand, start, len, del, ins, lclip, rclip, sa,
                             mapq = mapq, snv_alt = snv_alt))
}
# a supplementary alignment sharing a primary's name
lines = c(lines, make_read(999, "1", "+", 5200, 600, flag_extra = 2048L, name = "read011",
                           sa = "chr1,4000,+,600M,60,0;"))

hdr = c("@HD\tVN:1.6\tSO:unsorted", sprintf("@SQ\tSN:chr1\tLN:%d", ref_len),
        "@SQ\tSN:chr3\tLN:100000", "@SQ\tSN:chr7\tLN:100000")
writeLines(c(hdr, lines), file.path(out, "reads.sam"))

vcf = c(
  "##fileformat=VCFv4.2",
  sprintf("##contig=<ID=chr1,length=%d>", ref_len),
  "##contig=<ID=chr7,length=100000>",
  '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="SV type">',
  '##INFO=<ID=END,Number=1,Type=Integer,Description="End">',
  '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Len">',
  "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
  "chr1\t3000\tsnv1\tA\tG\t50\tPASS\t.",
  sprintf("chr1\t3899\tdel1\t%s\t<DEL>\t50\tPASS\tSVTYPE=DEL;END=4050;SVLEN=-150", ref[3899]),
  sprintf("chr1\t6000\tbnd1\t%s\t%s[chr7:5000[\t50\tPASS\tSVTYPE=BND", ref[6000], ref[6000])
)
writeLines(vcf, file.path(out, "vars.vcf"))
