#!/usr/bin/env bash
# 02_modkit_pileup.sh
#
# Z Dorado modBAM vytvoří per-sample bedMethyl (region/site-level metylace).
# Tohle je "hrubý" metylační profil vzorku - použije se jak pro QC/vizualizaci
# v IGV, tak jako vstup pro 07_marker_methylation_table.py (PCNSL markery).
#
# Použití:
#   ./02_modkit_pileup.sh <sample.modbam> <ref.fa> <output_prefix>
#
# Ověřte si aktuální flagy: `modkit pileup --help` (verze modkitu se mění).
set -euo pipefail

BAM="$1"
REF="$2"
OUT="$3"

if ! command -v modkit >/dev/null 2>&1; then
  echo "[CHYBA] modkit nenalezen v PATH (cargo install / ONT release binary)." >&2
  exit 1
fi

if [ ! -f "${BAM}.bai" ] && [ ! -f "${BAM%.bam}.bai" ]; then
  samtools index -@4 "$BAM"
fi

modkit pileup "$BAM" "${OUT}.bedmethyl" \
  --ref "$REF" \
  --cpg \
  --combine-strands \
  --only-tabs \
  --threads 8

sort -k1,1 -k2,2n "${OUT}.bedmethyl" > "${OUT}.sorted.bedmethyl"
rm -f "${OUT}.bedmethyl"
bgzip -f "${OUT}.sorted.bedmethyl"
tabix -p bed "${OUT}.sorted.bedmethyl.gz"

echo "[OK] ${OUT}.sorted.bedmethyl.gz (+.tbi) hotovo"
echo "     Sloupce (dle modkit dokumentace, ověřte 'modkit pileup --help'):"
echo "     1-3 chrom/start/end, 4 mod kód, 5 skóre, 6 strand, 10 Nvalid_cov,"
echo "     11 %modifikace, 12 Nmod, 13 Ncanonical, 14 Nother_mod, ..."
