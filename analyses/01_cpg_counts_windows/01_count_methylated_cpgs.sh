#!/usr/bin/env bash
# Single-pass extraction from modkit pileup bed files, per sample.
#
# This is the one unavoidable full-genome scan (each modkit.bed is ~180GB).
# Everything downstream (windowing, Venn, PCA) works off the much smaller
# files written here - no need to touch the raw modkit.bed again.
#
# For each sample writes to analyses/01_cpg_counts_windows/output/:
#   {sample}_meth_cpgs.bed   chrom,start,end of CpGs called methylated
#                            (mod_code == m, coverage >= MIN_COV, pct >= MIN_PCT)
#   {sample}_cov_m.bedgraph  chrom,start,end,pct_modified for ALL m-code CpGs
#                            with coverage >= MIN_COV (input for bedtools map)
#   methylated_cpg_counts.tsv  sample<TAB>n_methylated_cpgs (one row per sample)
#
# Only needs awk - no bedtools/conda dependency for this step.

set -euo pipefail

# ---- config ---------------------------------------------------------------
RESULTS_DIR="results"                 # contains {sample}/methylation/{sample}_modkit.bed
SAMPLES=(cfDNA_BAR cfDNA_LOS cfDNA_MUD cfDNA_Pilot1 cfDNA_VN)
OUT_DIR="analyses/01_cpg_counts_windows/output"   # spoustet z methylation_profiles/
MIN_COV=5        # minimum valid coverage (modkit col 10) to consider a site
MIN_PCT=50       # minimum percent-modified (modkit col 11) to call "methylated"
# -----------------------------------------------------------------------------

mkdir -p "$OUT_DIR"
COUNTS_FILE="$OUT_DIR/methylated_cpg_counts.tsv"
echo -e "sample\tn_methylated_cpgs" > "$COUNTS_FILE"

for sample in "${SAMPLES[@]}"; do
    bed="$RESULTS_DIR/$sample/methylation/${sample}_modkit.bed"
    meth_out="$OUT_DIR/${sample}_meth_cpgs.bed"
    bedgraph_out="$OUT_DIR/${sample}_cov_m.bedgraph"
    count_out="$OUT_DIR/${sample}_meth_count.txt"

    if [[ ! -f "$bed" ]]; then
        echo "WARNING: missing $bed, skipping $sample" >&2
        continue
    fi

    echo "[$(date '+%F %T')] processing $sample ..." >&2

    awk -v OFS='\t' \
        -v methout="$meth_out" \
        -v bedgraphout="$bedgraph_out" \
        -v countout="$count_out" \
        -v mincov="$MIN_COV" \
        -v minpct="$MIN_PCT" '
    $4 == "m" && $10 >= mincov {
        print $1, $2, $3, $11 > bedgraphout
        if ($11 >= minpct) {
            print $1, $2, $3 > methout
            meth_count++
        }
    }
    END {
        print meth_count+0 > countout
    }
    ' "$bed"

    n=$(cat "$count_out")
    echo -e "${sample}\t${n}" >> "$COUNTS_FILE"
    echo "[$(date '+%F %T')] $sample: $n methylated CpGs" >&2
done

echo "Done. Summary: $COUNTS_FILE"
echo "Next: bash analyses/01_cpg_counts_windows/01b_bedtools_windows.sh"
