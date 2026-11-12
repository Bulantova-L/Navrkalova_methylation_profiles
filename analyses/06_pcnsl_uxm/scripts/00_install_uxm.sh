#!/usr/bin/env bash
# 00_install_uxm.sh
#
# UXM_deconv (https://github.com/nloyfer/UXM_deconv) nejde nainstalovat
# přes pip/conda - nemá setup.py ani pyproject.toml, je to jen sbírka
# skriptů. Oficiální README doporučuje přesně tohle: naklonovat repo a
# přidat ho do PATH. Ověřeno funkční.
#
# Tenhle skript navíc rovnou ukáže, kde v klonu najdete referenční atlas
# (Loyfer 2023) - je součástí repozitáře, není potřeba ho hledat/stahovat
# odjinud.
#
# Použití:
#   ./00_install_uxm.sh [cílový_adresář, default ./tools]
#
# Po doběhnutí si vypsaný `export PATH=...` řádek přidejte do svého
# ~/.bashrc / aktivačního skriptu conda prostředí, ať se `uxm` najde
# i v nových shellech.
set -euo pipefail

TARGET_DIR="${1:-./tools}"
mkdir -p "$TARGET_DIR"
cd "$TARGET_DIR"

if [ -d "UXM_deconv" ]; then
  echo "[INFO] $TARGET_DIR/UXM_deconv už existuje, aktualizuji (git pull)..."
  git -C UXM_deconv pull --ff-only
else
  echo "[INFO] Klonuji UXM_deconv do $TARGET_DIR/UXM_deconv ..."
  git clone https://github.com/nloyfer/UXM_deconv.git
fi

chmod +x UXM_deconv/uxm UXM_deconv/src/*.py 2>/dev/null || true

UXM_ABS_PATH="$(cd UXM_deconv && pwd)"

echo
echo "[OK] UXM_deconv připraven v: $UXM_ABS_PATH"
echo "     Referenční atlasy (Loyfer 2023, hg38) jsou přímo tady - není"
echo "     potřeba nic stahovat externě:"
echo "       $UXM_ABS_PATH/supplemental/Atlas.U25.l4.hg38.tsv   (25 širokých buněčných typů)"
echo "       $UXM_ABS_PATH/supplemental/Atlas.U250.l4.hg38.tsv  (~250 jemnějších subtypů)"
echo
echo "     Přidejte do PATH (a do ~/.bashrc, aby to zůstalo trvale):"
echo "       export PATH=\"$UXM_ABS_PATH:\$PATH\""
echo
echo "     Ověření: uxm --help"
