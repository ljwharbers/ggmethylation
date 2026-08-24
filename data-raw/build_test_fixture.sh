#!/usr/bin/env bash
#
# Rebuild the integration-test BAM fixture.
#
# Run this manually on the VSC cluster from the package root:
#   bash data-raw/build_test_fixture.sh
#
# ---------------------------------------------------------------------------
# Why this locus
# ---------------------------------------------------------------------------
# MEG3 is an imprinted lncRNA whose promoter DMR carries allele-specific
# methylation.  With haplotype-tagged reads the two haplotypes therefore differ
# for a real biological reason, which makes this window a genuine test of the
# grouped smooth, the delta track and the CI ribbon rather than a test against
# noise.
#
# MEG3 gene span in CHM13 (hs1.ncbiRefSeq): chr14:95060996-95095916 (+)
# DLK1, immediately upstream:               chr14:94961805-94973139 (+)
# The window below brackets the MEG3 TSS / promoter DMR.
#
# ---------------------------------------------------------------------------
# Fixture contents (verified after building)
# ---------------------------------------------------------------------------
#   52 reads, 969 KB
#   MM and ML tags on 52/52 reads          -> MM/ML parsing
#   HP on 51/52, split 30 (hap1) / 21 (hap2) -> group_tag = "HP", delta track
#   exactly 1 read without HP              -> drop_na_group path
#   both mod codes on every read: A+a (6mA) and C+m (5mC)
#                                          -> multi-modification shape aesthetic
#   MM entries carry no ?/. flag           -> "omit unlisted canonical bases"
#
# If you change the region or the sampling, re-run the verification block at
# the bottom and update tests/testthat/test-integration.R accordingly.
# ---------------------------------------------------------------------------

set -euo pipefail

SAMTOOLS="${SAMTOOLS:-$VSC_DATA/software/micromamba/bin/samtools}"

SRC="/staging/leuven/stg_00096/home/meftyc/projects/phd/DeepFiberNet/data/output/py_output/longphase/phased.bam"
REGION="chr14:95055000-95070000"
OUT="tests/testthat/fixtures/hg002_fiberseq_MEG3.bam"

# The source is a 130 GB longphase-phased HG002 PacBio Fiber-seq BAM aligned to
# CHM13.  Its reads carry ~40 aux tags each; keeping only the six the package
# actually needs takes the subset from ~3.5 MB of SAM text down to under 1 MB.
mkdir -p "$(dirname "$OUT")"

"$SAMTOOLS" view -b \
    --keep-tag MM,ML,HP,PS,NM \
    -o "$OUT" \
    "$SRC" "$REGION"

"$SAMTOOLS" index "$OUT"

# --- verification -----------------------------------------------------------
echo "--- $OUT ($(du -h "$OUT" | cut -f1)) ---"
echo "reads:      $("$SAMTOOLS" view -c "$OUT")"
echo "with MM:    $("$SAMTOOLS" view "$OUT" | grep -c 'MM:Z:')"
echo "with ML:    $("$SAMTOOLS" view "$OUT" | grep -c 'ML:B:')"
echo "HP split:"
"$SAMTOOLS" view "$OUT" | grep -o 'HP:i:[0-9]*' | sort | uniq -c
echo "reads without HP: $("$SAMTOOLS" view "$OUT" | grep -vc 'HP:i:')"
echo "with C+m:   $("$SAMTOOLS" view "$OUT" | grep -c 'C+m')"
echo "with A+a:   $("$SAMTOOLS" view "$OUT" | grep -c 'A+a')"
