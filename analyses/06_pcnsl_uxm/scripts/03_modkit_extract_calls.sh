#!/usr/bin/env bash
# 03_modkit_extract_calls.sh
#
# Extrahuje per-read, per-CpG modifikační volání z modBAM - to je vstup
# pro fragment-level konverzi do PAT formátu (skript 04_calls_to_pat.py).
# Bez tohoto kroku nejde dělat UXM dekonvoluci (ta potřebuje vzory na
# úrovni jednotlivých molekul/readů, ne jen agregovanou % metylaci).
#
# DOPORUČENÍ: pokud vás zajímá primárně bod 2 (PCNSL markery), omezte
# extrakci na regiony markerů (+- flank) přes --regions <bed> - ušetří
# to obrovské množství času/paměti oproti celogenomové extrakci.
#
# Použití:
#   ./03_modkit_extract_calls.sh <sample.modbam> <ref.fa> <output.tsv.gz> [regions.bed]
set -euo pipefail

BAM="$1"
REF="$2"
OUT="$3"
REGIONS="${4:-}"

if ! command -v modkit >/dev/null 2>&1; then
  echo "[CHYBA] modkit nenalezen v PATH." >&2
  exit 1
fi

RAW_OUT="${OUT%.gz}"

EXTRA_ARGS=()
if [ -n "$REGIONS" ]; then
  EXTRA_ARGS+=(--include-bed "$REGIONS")
fi

# POZOR: přesný název subpříkazu/flagů se mezi verzemi modkitu lišil
# (novější verze: `modkit extract calls`, starší: `modkit extract`).
# Ověřte: modkit extract --help  /  modkit extract calls --help
modkit extract calls "$BAM" "$RAW_OUT" \
  --ref "$REF" \
  --mapped-only \
  --threads 8 \
  "${EXTRA_ARGS[@]}"

gzip -f "$RAW_OUT"

echo "[OK] $OUT hotovo (per-read CpG volání)"
echo "     Očekávané klíčové sloupce (dle dokumentace, ověřte hlavičku souboru!):"
echo "     read_id, ref_position, chrom, call_prob, call_code (m/h/-), ..."
