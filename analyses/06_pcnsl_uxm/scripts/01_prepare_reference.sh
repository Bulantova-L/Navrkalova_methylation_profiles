#!/usr/bin/env bash
# 01_prepare_reference.sh
#
# Připraví hg38 referenční FASTA + CpG index potřebný pro wgbstools/UXM
# (a pro konverzi modkit -> PAT ve skriptu 04_calls_to_pat.py).
#
# POZOR: přesné jméno subpříkazu/flagy wgbstools se mezi verzemi mění.
# Než spustíte na produkčních datech, ověřte:
#   wgbstools --help
#   wgbstools init_genome --help
#
# Použití:
#   ./01_prepare_reference.sh [výstupní_adresář]
set -euo pipefail

REF_DIR="${1:-./reference/hg38}"
mkdir -p "$REF_DIR"
cd "$REF_DIR"

if ! command -v wgbstools >/dev/null 2>&1; then
  echo "[CHYBA] wgbstools nenalezen v PATH." >&2
  echo "        Instalace: git clone https://github.com/nloyfer/wgbs_tools && cd wgbs_tools && pip install -e ." >&2
  exit 1
fi

if [ ! -f genome.fa ]; then
  echo "[INFO] Stahuji GRCh38/hg38 referenci (UCSC, no-alt analysis set)..."
  wget -O genome.fa.gz \
    https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/latest/hg38.fa.gz
  gunzip genome.fa.gz
fi

if ! command -v samtools >/dev/null 2>&1; then
  echo "[CHYBA] samtools nenalezen v PATH (conda install -c bioconda samtools)." >&2
  exit 1
fi
samtools faidx genome.fa

echo "[INFO] Inicializuji wgbstools genome referenci pro hg38 (staví CpG index)..."
# Tenhle krok vytvoří CpG.bed.gz (mapování genomová pozice <-> CpG index),
# které skript 04_calls_to_pat.py používá k převodu souřadnic z modkitu
# do PAT formátu. Název výstupního adresáře/parametrů si ověřte přes --help.
wgbstools init_genome genome.fa --name hg38 --out_dir "$(pwd)/wgbstools_ref"

echo "[OK] Reference + CpG index by měly být v: $REF_DIR/wgbstools_ref"
echo "     Zkontrolujte, že tam je CpG.bed.gz (nebo obdobný soubor) a chrom.sizes."
echo "     Pokud struktura/jména souborů neodpovídají, spusťte 'wgbstools init_genome --help'"
echo "     a uprav skript 04_calls_to_pat.py podle skutečného výstupu."
