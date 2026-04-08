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

BAM_PTCL8  <- "/staging/leuven/stg_00096/home/averham/LR_SOMATIC_hg38/PTCL8_PB/bamfiles/PTCL8_PB_tumor.bam"
BAM_AITL4  <- "/staging/leuven/stg_00096/home/averham/LR_SOMATIC_hg38/AITL4/bamfiles/AITL4_tumor.bam"
VCF_AITL4  <- "/staging/leuven/stg_00096/home/averham/LR_SOMATIC_hg38/AITL4/variants/clairsto/somatic.vcf.gz"

REGION     <- "chr8:127733434-127744951"
REGION     <- "chr8:127598491-127599663"

# =============================================================================
# 1. Basic single-sample, ungrouped
# =============================================================================

meth <- read_methylation(BAM_PTCL8, REGION, mod_code = "m")

print(meth)
summary(meth)

plot_methylation(meth)


# =============================================================================
# 2. Grouped by haplotype tag (HP)
# =============================================================================

meth_hp <- read_methylation(BAM_PTCL8, REGION,
                             mod_code    = "m",
                             group_tag   = "HP",
                             drop_na_group = TRUE)

print(meth_hp)

plot_methylation(meth_hp)


# =============================================================================
# 3. Grouped by SNV genotype
# =============================================================================

meth_snv <- read_methylation(BAM_PTCL8, REGION,
                              mod_code     = "m",
                              snv_position = 127740000L,
                              ref_base     = "A",
                              alt_base     = "T")

plot_methylation(meth_snv)


# =============================================================================
# 4. Read filters — MAPQ, strand, length, downsampling
# =============================================================================

meth_filtered <- read_methylation(BAM_PTCL8, REGION,
                                   mod_code        = "m",
                                   min_mapq        = 20L,
                                   strand_filter   = "+",
                                   min_read_length = 5000L,
                                   max_reads       = 50L)

plot_methylation(meth_filtered)


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
                 dot_size    = 2)

# Colour reads by strand instead of a single colour
plot_methylation(meth,
                 colour_strand  = TRUE,
                 strand_colours = c("+" = "#2196F3", "-" = "#FF5722"))

# Custom group colours
plot_methylation(meth_hp,
                 group_colours = c("1" = "#388E3C", "2" = "#7B1FA2"))

# Adjust panel height ratio (reads : smooth)
plot_methylation(meth_hp, panel_heights = c(4, 1))


# =============================================================================
# 7. Gene annotations — from UCSC ncbiRefSeq (canonical gene symbols)
# =============================================================================

# Use genome = "hg38" to automatically download and cache the UCSC ncbiRefSeq
# GTF (hs1.ncbiRefSeq.gtf.gz). Gene names are canonical symbols (e.g. "MYC").
# The GTF is cached locally; subsequent calls skip the download.
#
# To change the cache location (e.g. to scratch on HPC):
#   options(ggmethylation.cache_dir = "$VSC_SCRATCH/ggmethylation_cache")
annot <- read_annotations(genome = "hg38", region = REGION)

print(annot)

# Second call reuses both the downloaded GTF and the TxDb cache
annot2 <- read_annotations(genome = "hg38", region = "chr8:127700000-127800000")

# Opt out of transcript collapsing to show all isoforms
annot_all_isoforms <- read_annotations(genome = "hg38", region = REGION,
                                       collapse_transcripts = FALSE)
print(annot_all_isoforms)

# Clear TxDb cache to free memory or force re-parse (does not delete downloaded GTF)
clear_annotation_cache()

# hg38 genome is also supported
# annot_hg38 <- read_annotations(genome = "hg38", region = REGION)

# From a local GTF file (e.g. custom annotation)
# annot_gtf <- read_annotations(gtf = "/path/to/custom.gtf", region = REGION)

# From a pre-built TxDb (e.g., Bioconductor annotation package)
# library(TxDb.Hsapiens.UCSC.hg38.knownGene)
# annot_txdb <- read_annotations(txdb = TxDb.Hsapiens.UCSC.hg38.knownGene,
#                                region = REGION)

# IGV-style gene track: thin intron lines, medium UTR boxes, full CDS boxes
plot_methylation(meth_hp, annotations = annot)

# Combined with all isoforms visible
plot_methylation(meth_hp, annotations = annot_all_isoforms)


# =============================================================================
# 8. CIGAR structural variant visualisation
# =============================================================================

# show_cigar = TRUE overlays per-read structural features derived from the
# CIGAR string:
#   - Purple bold "I" at insertion positions
#   - Thin black line across deletion spans
#   - Orange triangle markers at soft/hard-clipped read ends

plot_methylation(meth, show_cigar = TRUE)

# Works with grouping and annotations too
plot_methylation(meth_hp, annotations = annot, show_cigar = TRUE)


# =============================================================================
# 9. Variant overlay (VCF)
# =============================================================================

vars <- read_variants(VCF_AITL4, REGION)

print(vars)

meth_aitl4 <- read_methylation(BAM_AITL4, REGION,
                                mod_code  = "m",
                                group_tag = "HP")

plot_methylation(meth_aitl4, variants = vars)

# Combined: annotations + variants
plot_methylation(meth_aitl4,
                 annotations = annot,
                 variants    = vars)


# =============================================================================
# 10. Multi-sample comparison
# =============================================================================

merged <- merge_methylation(PTCL8 = meth_hp,
                             AITL4 = meth_aitl4)

# Equivalent using .list
merged <- merge_methylation(.list = list(PTCL8 = meth_hp, AITL4 = meth_aitl4))

print(merged)
summary(merged)

plot_methylation(merged)

plot_methylation(merged, annotations = annot, variants = vars)


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

p <- plot_methylation(meth_hp, annotations = annot)

# Access individual panels
reads_panel  <- p[[1]]
smooth_panel <- p[[2]]
gene_panel   <- p[[3]]

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
