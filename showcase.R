# =============================================================================
# ggmethylation — interactive showcase script
#
# Demonstrates all major use cases of the package.
# Requires access to the staging BAM/VCF files on the HPC cluster.
#
# Run sections interactively; do not source() the whole file at once.
# =============================================================================

devtools::load_all()  # or: install.packages(".", repos = NULL, type = "source")

library(ggplot2)
library(patchwork)

# -- Shared paths -------------------------------------------------------------

BAM = "/staging/leuven/stg_00096/home/averham/LR_SOMATIC_T2T/BL1/bamfiles/BL1_tumor.bam"
VCF_SVs = "/staging/leuven/stg_00096/home/averham/LR_SOMATIC_T2T/BL1/variants/severus/somatic_SVs/severus_somatic.vcf.gz"
VCF_SNVs = "/staging/leuven/stg_00096/home/averham/LR_SOMATIC_T2T/BL1/variants/clairsto/somatic.vcf.gz"

REGION     = "chr8:128861685-128872405"

# =============================================================================
# 1. Basic single-sample, ungrouped
# =============================================================================

meth = read_methylation(BAM, REGION, mod_code = "m")

print(meth)
summary(meth)

plot_methylation(meth)


# =============================================================================
# 2. Grouped by haplotype tag (HP)
# =============================================================================

meth_hp = read_methylation(BAM, REGION,
                             mod_code    = "m",
                             group_tag   = "HP",
                             drop_na_group = FALSE)

print(meth_hp)

plot_methylation(meth_hp)


# =============================================================================
# 3. Grouped by SNV genotype
# =============================================================================

meth_snv = read_methylation(BAM, REGION,
                              mod_code     = "m",
                              snv_position = 127740000L,
                              ref_base     = "A",
                              alt_base     = "T")

plot_methylation(meth_snv)


# =============================================================================
# 4. Read filters — MAPQ, strand, length, downsampling
# =============================================================================

meth_filtered = read_methylation(BAM, REGION,
                                   mod_code        = "m",
                                   min_mapq        = 20L,
                                   strand_filter   = "+",
                                   min_read_length = 5000L,
                                   max_reads       = 50L)

plot_methylation(meth_filtered)

# Per-group downsampling: apply max_reads cap independently within each group
# (useful when one haplotype has far more reads than the other)
meth_hp_balanced = read_methylation(BAM, REGION,
                                      mod_code            = "m",
                                      group_tag           = "HP",
                                      drop_na_group       = TRUE,
                                      max_reads           = 50L,
                                      per_group_downsample = TRUE)

plot_methylation(meth_hp_balanced)


# =============================================================================
# 5. sort_by options
# =============================================================================

plot_methylation(meth_hp, sort_by = "group")
plot_methylation(meth_hp, sort_by = "start")
plot_methylation(meth_hp, sort_by = "mean_mod_prob")


# =============================================================================
# 6. Colour and aesthetics customisation
# =============================================================================

plot_methylation(meth,
                 colour_low  = "#FFFFFF",
                 colour_high = "#1A237E",
                 line_width  = 0.5)

# Colour reads by strand instead of a single colour
plot_methylation(meth,
                 colour_strand  = TRUE,
                 strand_colours = c("+" = "#2196F3", "-" = "#FF5722"))

# Custom group colours
plot_methylation(meth_hp,
                 group_colours = c("1" = "#388E3C", "2" = "#7B1FA2"))

# Adjust panel height ratio (reads : smooth)
plot_methylation(meth_hp, panel_heights = c(4, 1))

# Custom shape mapping when multiple modification codes are present
# (e.g. 5mC "m" and 5hmC "h" from a dual-modification BAM)
# meth_dual = read_methylation(BAM, REGION, mod_code = c("m", "h"))
# plot_methylation(meth_dual, mod_code_shapes = c(m = 16L, h = 17L))


# =============================================================================
# 7. Gene annotations — from UCSC ncbiRefSeq (canonical gene symbols)
# =============================================================================

# Use genome = "hg38" or "chm13" to automatically download and cache the UCSC
# ncbiRefSeq GTF. Gene names are canonical symbols (e.g. "MYC").
# The GTF is cached locally; subsequent calls skip the download.
#
# Supported genomes:
#   "hg38"  — hg38.ncbiRefSeq.gtf.gz
#   "chm13" — hs1.ncbiRefSeq.gtf.gz  (T2T CHM13)
#
# To change the cache location (e.g. to scratch on HPC):
#   options(ggmethylation.cache_dir = "$VSC_SCRATCH/ggmethylation_cache")
annot = read_annotations(genome = "chm13", region = REGION)

print(annot)

# Second call reuses both the downloaded GTF and the TxDb cache
annot2 = read_annotations(genome = "chm13", region = "chr8:127700000-127800000")

# Opt out of transcript collapsing to show all isoforms
annot_all_isoforms = read_annotations(genome = "chm13", region = REGION,
                                       collapse_transcripts = FALSE)
print(annot_all_isoforms)

# Clear TxDb cache to free memory or force re-parse (does not delete downloaded GTF)
clear_annotation_cache()

# CHM13 / T2T genome
# annot_chm13 = read_annotations(genome = "chm13", region = REGION)

# From a local GTF file (e.g. custom annotation)
# annot_gtf = read_annotations(gtf = "/path/to/custom.gtf", region = REGION)

# From a pre-built TxDb (e.g., Bioconductor annotation package)
# library(TxDb.Hsapiens.UCSC.hg38.knownGene)
# annot_txdb = read_annotations(txdb = TxDb.Hsapiens.UCSC.hg38.knownGene,
#                                region = REGION)

# IGV-style gene track: thin intron lines, medium UTR boxes, full CDS boxes
plot_methylation(meth_hp, annotations = annot)

# Combined with all isoforms visible
plot_methylation(meth_hp, annotations = annot_all_isoforms)


# =============================================================================
# 8. CIGAR structural variant visualisation
# =============================================================================

# show_cigar = TRUE is the default. It overlays per-read structural features
# derived from the CIGAR string:
#   - Purple bold "I" at insertion positions
#   - Thin black line across deletion spans (with IGV-style gap in read bar)
# Only indels >= min_indel_size (default 50 bp) are shown.

plot_methylation(meth)                          # CIGAR on by default
plot_methylation(meth, show_cigar = FALSE)      # disable CIGAR overlay
plot_methylation(meth, min_indel_size = 20L)    # show smaller indels (>= 20 bp)

# Works with grouping and annotations too
plot_methylation(meth_hp, annotations = annot)


# =============================================================================
# 8b. Supplementary alignment visualisation
# =============================================================================

# show_supplementary = TRUE is the default. Reads with a SA (supplementary
# alignment) tag get a coloured halo at their clipped end(s), coloured by the
# mate chromosome. When variants = <vcf BND data>, SA reads whose breakpoint
# matches a BND call are additionally outlined in dark red.

plot_methylation(meth_hp)                             # SA halos on by default
plot_methylation(meth_hp, show_supplementary = FALSE) # disable SA reads

# Combined: CIGAR off, SA halos on
plot_methylation(meth_hp, show_cigar = FALSE, show_supplementary = TRUE)




# =============================================================================
# 9. Variant overlay (VCF)
# =============================================================================

vars = read_variants(VCF_SVs, REGION)

print(vars)
plot_methylation(meth_hp, variants = vars)

# BND matching tolerance (bp) for VCF-validated SA border
plot_methylation(meth_hp, variants = vars, bnd_match_tol = 100L)



# Combined: annotations + variants
plot_methylation(meth_hp,
                 annotations = annot,
                 variants    = vars)


# =============================================================================
# 10. Multi-sample comparison
# =============================================================================

merged = merge_methylation(sample1 = meth_hp,
                           sample2 = meth_hp)

# Equivalent using .list
merged = merge_methylation(.list = list(sample1 = meth_hp, sample2 = meth_hp))

print(merged)
summary(merged)

plot_methylation(merged)

plot_methylation(merged, annotations = annot, variants = vars,
                 smooth_span = 0.1)


# =============================================================================
# 11. Write methylation data to disk
# =============================================================================

write_methylation(meth, prefix = "/tmp/ptcl8_meth", format = "tsv")
write_methylation(meth, prefix = "/tmp/ptcl8_meth", format = "bed", gzip = TRUE)

# Pipe-friendly (returns data invisibly)
meth |>
  write_methylation(prefix = "/tmp/ptcl8_pipe", format = "tsv") |>
  plot_methylation()


# =============================================================================
# 12. Manipulating patchwork output
# =============================================================================

# Panel order depends on which panels are present:
#   No annotations:  [[1]] reads,  [[2]] smooth
#   With annotations: [[1]] gene,  [[2]] reads, [[3]] smooth

p_no_annot = plot_methylation(meth_hp)
reads_panel_2  = p_no_annot[[1]]
smooth_panel_2 = p_no_annot[[2]]

p = plot_methylation(meth_hp, annotations = annot)
gene_panel   = p[[1]]
reads_panel  = p[[2]]
smooth_panel = p[[3]]

# Modify a panel and reassemble
reads_panel + labs(y = "Read lane") + theme_classic()

# Add a global title
p + plot_annotation(
  title    = REGION,
  subtitle = "CpG methylation — PTCL8 tumour (haplotype-grouped)"
)

# Apply a theme across all panels
p & theme(text = element_text(size = 14))
p & theme_bw()

# Mark a position of interest
reads_panel + geom_vline(xintercept = 127740000, colour = "blue",
                          linetype = "dashed", linewidth = 0.5)
