#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"   # analyses/07_gse54200_dlbcl

[[ -s reference/hg17ToHg38.over.chain.gz ]] || wget -P reference https://hgdownload.soe.ucsc.edu/goldenPath/hg17/liftOver/hg17ToHg38.over.chain.gz

# liftover to hg38 is done in R (DLBCL_metylace.R), which writes reference/output_hg38_new.bed
Rscript DLBCL_metylace.R

# sort so it matches the *_modkit_sorted.bed chrom/pos order for bedtools -sorted/map
# reference/output_hg38_sorted.bed = GSE54200_MANIFEST in the Snakefile
sort -k1,1 -k2,2n reference/output_hg38_new.bed > reference/output_hg38_sorted.bed

# NOTE: the loop below was later moved into the Snakefile (rule bedtools_gse54200_map ->
# results/<sample>/methylation/<sample>_gse54200_mapped.bed), which is what DLBCL_metylace.R reads.
# Kept for reference; output/cfDNA_*_in_GSE54200loci.bed are from this older version.

# bedtools itself is single-threaded (no -threads/-T option on any subcommand);
# the 4 samples are independent, so run them concurrently in the background instead
for s in cfDNA_BAR cfDNA_MUD cfDNA_Pilot1 cfDNA_VN; do
  (
    modkit_bed=../../results/${s}/methylation/${s}_modkit_sorted.bed
    # m and h rows for the same CpG share the same Nvalid_cov (col10) -- it's the
    # total coverage at that position, not per-mod-code -- so first collapse each
    # position to a single row: Nmod = Nmod_m + Nmod_h, Nvalid_cov taken once.
    # (still one pass, still streamed -- no separate huge file written to disk)
    bedtools map -a reference/output_hg38_sorted.bed \
      -b <(awk -F'\t' -v OFS='\t' '
            $4 == "m" || $4 == "h" {
              key = $1 OFS $2
              if (key != prevkey) {
                if (prevkey != "") print prevchrom, prevstart, prevend, nmod, ncov
                nmod = 0; ncov = $10
              }
              nmod += $12
              prevkey = key; prevchrom = $1; prevstart = $2; prevend = $3
            }
            END { if (prevkey != "") print prevchrom, prevstart, prevend, nmod, ncov }
          ' "$modkit_bed") \
      -c 4,5 -o sum,sum -null 0 \
    | awk -F'\t' 'BEGIN{OFS="\t"} {pct=($NF>0)?100*$(NF-1)/$NF:"NA"; print $0, pct}' \
    > output/${s}_5mC5hmC_in_GSE54200loci.bed
  ) &
done
wait

