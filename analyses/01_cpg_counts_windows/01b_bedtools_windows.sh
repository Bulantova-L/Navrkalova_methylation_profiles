#!/usr/bin/env bash
# Bins per-CpG methylation (output/{sample}_cov_m.bedgraph, from
# 01_count_methylated_cpgs.sh) into fixed-size genome windows using bedtools,
# for the PCA step.
#
# Writes:
#   output/windows.bed           chrom,start,end (genome-wide bins, made once)
#   output/{sample}_windows.bed  chrom,start,end,mean_pct,n_cpgs (bedtools map output)
#
# Requires bedtools on PATH.

set -euo pipefail

# ---- config ---------------------------------------------------------------
SAMPLES=(cfDNA_BAR cfDNA_LOS cfDNA_MUD cfDNA_Pilot1 cfDNA_VN)
OUT_DIR="analyses/01_cpg_counts_windows/output"   # spoustet z methylation_profiles/
CHROM_SIZES="data/hg38.chrom.sizes"
WINDOW_SIZE=100000
SORT_MEM="2G"      # cap sort's in-memory buffer so it doesn't get OOM-killed
SORT_TMP="$OUT_DIR/tmp"  # disk-backed scratch dir for sort (NOT /tmp - that's
                          # often a small RAM-backed tmpfs on pods/containers,
                          # which just moves the OOM problem instead of fixing it)
# -----------------------------------------------------------------------------

WINDOWS_BED="$OUT_DIR/windows.bed"
CHROM_SIZES_CHR="$OUT_DIR/chrom.sizes.chr"

mkdir -p "$OUT_DIR" "$SORT_TMP"

# modkit.bed uses "chr1", "chr2", ... - normalize chrom.sizes to match in case
# it uses bare names ("1", "2", ...), which is what causes bedtools map to
# report windows.bed as missing chromosomes that the bedgraph does have.
awk 'BEGIN{OFS="\t"} {chrom=$1; if (chrom !~ /^chr/) chrom="chr"chrom; print chrom, $2}' \
    "$CHROM_SIZES" > "$CHROM_SIZES_CHR"

echo "[$(date '+%F %T')] creating $WINDOWS_BED ..." >&2
bedtools makewindows -g "$CHROM_SIZES_CHR" -w "$WINDOW_SIZE" | LC_ALL=C sort -k1,1 -k2,2n > "$WINDOWS_BED"

for sample in "${SAMPLES[@]}"; do
    bedgraph="$OUT_DIR/${sample}_cov_m.bedgraph"
    bedgraph_sorted="$OUT_DIR/${sample}_cov_m.sorted.bedgraph"
    win_out="$OUT_DIR/${sample}_windows.bed"

    if [[ ! -f "$bedgraph" ]]; then
        echo "WARNING: missing $bedgraph, skipping $sample" >&2
        continue
    fi

    # bedtools map requires -a and -b in the same strict lexicographic order
    # (chr1, chr10, chr11, ... chr2, ...) - modkit's own output order doesn't
    # match that, so both files must be explicitly (re)sorted this way.
    echo "[$(date '+%F %T')] sorting $bedgraph ..." >&2
    LC_ALL=C sort -k1,1 -k2,2n -S "$SORT_MEM" -T "$SORT_TMP" "$bedgraph" > "$bedgraph_sorted"

    echo "[$(date '+%F %T')] mapping $sample onto windows ..." >&2
    bedtools map -a "$WINDOWS_BED" -b "$bedgraph_sorted" -c 4 -o mean,count > "$win_out"
done

echo "Done. Next: Rscript analyses/01_cpg_counts_windows/02_analysis.R"
