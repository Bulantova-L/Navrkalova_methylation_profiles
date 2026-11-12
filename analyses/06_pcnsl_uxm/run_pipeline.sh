#!/usr/bin/env bash
# run_pipeline.sh
#
# Orchestrátor - ukázka, jak jednotlivé kroky navazují na sebe pro JEDEN
# vzorek. Pro 5 vzorků (4x DLBCL + 1x zdravé cfDNA) tuto smyčku spusťte
# pro každý vzorek zvlášť (BAM soubory), pak teprve agregujte výstupy
# (05_run_uxm.sh bere všechny pat.gz najednou; 07_marker_methylation_table.py
# bere všechny bedmethyl najednou).
#
# NEPOUŠTĚJTE naslepo na celý dataset - projděte README.md, projděte si
# caveaty v hlavičkách jednotlivých skriptů a nejdřív vyzkoušejte na
# jednom malém vzorku / jednom regionu.
set -euo pipefail

REF_DIR="reference/hg38"
TOOLS_DIR="../../../tools"   # sdilena slozka Navrkalova_cfDNA/tools (UXM_deconv, wgbs_tools, ichorCNA); spoustet z analyses/06_pcnsl_uxm
UXM_DIR="${TOOLS_DIR}/UXM_deconv"
MARKERS_BED="pcnsl_markers.hg38.bed"

# --- -1) UXM_deconv (git clone + svazkový atlas, jen jednou) ---
scripts/00_install_uxm.sh "$TOOLS_DIR"
export PATH="$(cd "$UXM_DIR" && pwd):$PATH"

# --- 0) referenční genom + CpG index (jen jednou, sdílené pro všechny vzorky)
scripts/01_prepare_reference.sh "$REF_DIR"

# --- 1) PCNSL markery -> BED (jen jednou)
python3 scripts/06_build_pcnsl_marker_bed.py \
  --markers markers/pcnsl_markers_450k.tsv \
  --out "$MARKERS_BED" \
  --flank 250

declare -A BAMS=(
  [dlbcl1]="/cesta/k/dlbcl1.modbam"
  [dlbcl2]="/cesta/k/dlbcl2.modbam"
  [dlbcl3]="/cesta/k/dlbcl3.modbam"
  [dlbcl4]="/cesta/k/dlbcl4.modbam"
  [healthy1]="/cesta/k/healthy1.modbam"
)

for sample in "${!BAMS[@]}"; do
  bam="${BAMS[$sample]}"
  echo "=== Vzorek: $sample ==="

  # --- 2) bedMethyl (region-level profil, pro markery i QC)
  scripts/02_modkit_pileup.sh "$bam" "$REF_DIR/genome.fa" "out/${sample}"

  # --- 3) per-read CpG volání (omezeno na markerové regiony +/- okolí,
  #        ať to netrvá celý den - pro celogenomovou UXM dekonvoluci
  #        vynechte poslední argument, ale počítejte s výrazně delším
  #        během a větší pamětí)
  scripts/03_modkit_extract_calls.sh "$bam" "$REF_DIR/genome.fa" \
    "out/${sample}.calls.tsv.gz" "$MARKERS_BED"

  # --- 4) konverze na PAT (fragment-level, pro UXM)
  python3 scripts/04_calls_to_pat.py \
    --calls "out/${sample}.calls.tsv.gz" \
    --cpg-index "$REF_DIR/wgbstools_ref/CpG.bed.gz" \
    --out "out/${sample}.pat.gz"
done

# --- 5) UXM dekonvoluce přes všechny vzorky najednou (svazkový U25 atlas)
scripts/05_run_uxm.sh "$UXM_DIR" "out" "out/deconv_result" U25

# --- 6) metylace na PCNSL markerech přes všechny vzorky najednou
sample_args=()
for sample in "${!BAMS[@]}"; do
  sample_args+=("${sample}=out/${sample}.sorted.bedmethyl.gz")
done
python3 scripts/07_marker_methylation_table.py \
  --markers "$MARKERS_BED" \
  --out "out/pcnsl_marker_matrix.tsv" \
  "${sample_args[@]}"

echo "[HOTOVO] Zkontrolujte out/deconv_result.deconv.csv a out/pcnsl_marker_matrix.tsv"
