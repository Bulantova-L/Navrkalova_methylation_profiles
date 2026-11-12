#!/usr/bin/env bash
# 05_run_uxm.sh
#
# Spustí UXM dekonvoluci (Loyfer et al. 2023, Nature - "A DNA methylation
# atlas of normal human cell types") na sadě pat.gz souborů z kroku 04,
# proti referenčnímu atlasu, který je SOUČÁSTÍ UXM_deconv repozitáře -
# viz scripts/00_install_uxm.sh, který ho zpřístupní ve složce
# <UXM_deconv>/supplemental/.
#
# (Oprava oproti dřívější verzi tohoto skriptu: dřív jsem tu měl
# instrukce hledat/stahovat atlas externě, protože se mi to při rešerši
# nepodařilo ověřit. Po skutečném klonu repozitáře se ukázalo, že atlas
# je prostě přímo uvnitř - žádné externí hledání není potřeba.)
#
# Dvě varianty atlasu (obě hg38):
#   Atlas.U25.l4.hg38.tsv   - 25 širokých buněčných typů (rychlejší, hrubší)
#   Atlas.U250.l4.hg38.tsv  - ~250 jemnějších subtypů (podrobnější, pomalejší)
# Pro první orientační odhad ctDNA/hematopoetické frakce doporučuji začít
# s U25 variantou.
#
# Použití:
#   ./05_run_uxm.sh <cesta_k_UXM_deconv_klonu> <pat_dir> <output_prefix> [U25|U250]
set -euo pipefail

UXM_DIR="$1"
PAT_DIR="$2"
OUT="$3"
ATLAS_RESOLUTION="${4:-U25}"

if ! command -v uxm >/dev/null 2>&1; then
  echo "[CHYBA] 'uxm' nenalezen v PATH." >&2
  echo "        Spusťte nejdřív: scripts/00_install_uxm.sh a přidejte" >&2
  echo "        vypsanou cestu do PATH (export PATH=\"...UXM_deconv:\$PATH\")." >&2
  exit 1
fi

ATLAS_FILE="${UXM_DIR}/supplemental/Atlas.${ATLAS_RESOLUTION}.l4.hg38.tsv"
if [ ! -f "$ATLAS_FILE" ]; then
  echo "[CHYBA] Atlas nenalezen: $ATLAS_FILE" >&2
  echo "        Zkontrolujte cestu k UXM_deconv klonu a název souboru" >&2
  echo "        (ls \$UXM_DIR/supplemental/)." >&2
  exit 1
fi

echo "[INFO] Používám svazkový atlas: $ATLAS_FILE"

uxm deconv "$PAT_DIR"/*.pat.gz \
  --atlas "$ATLAS_FILE" \
  -o "${OUT}.deconv.csv"

echo "[OK] ${OUT}.deconv.csv - odhadované podíly buněčných typů / nádorové frakce"
echo "     Zkontrolujte, že atlas obsahuje referenci relevantní pro B-lymfocyty/"
echo "     lymfom (např. 'B-cell' typ) - bez toho nebude dekonvoluce tumor"
echo "     frakce u DLBCL smysluplná. Jako fallback/srovnání zvažte i"
echo "     meth_atlas (Moss 2018, https://github.com/nloyfer/meth_atlas) -"
echo "     jednodušší region-level NNLS dekonvoluce s vlastním atlasem."
