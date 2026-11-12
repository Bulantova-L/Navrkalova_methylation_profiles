#!/usr/bin/env bash
# Predzpracovani pro metylaci genovych oblasti (gen + UPSTREAM bp pred zacatkem genu)
#
#   1) GTF -> BED genovych oblasti: [start - UPSTREAM, end] na +, [start, end + UPSTREAM] na - retezci
#   2) FASTA -> BED vsech CpG (C na + retezci, G na - retezci); modkit pileup bezel bez --cpg,
#      proto je nutne CpG kontext dofiltrovat
#   3) Pro kazdy vzorek z modkit bedMethyl: jen 5mC ("m"), CpG, N_valid >= MIN_COV
#      -> na kazdou oblast sum(N_valid), sum(N_mod), pocet CpG volani
#
# Vystup (OUT_DIR): <sample>_genes_up<UPSTREAM>_cov<MIN_COV>.tsv
#   chr start end gene_id gene_name strand gene_biotype n_valid n_mod n_cpg   ("." = bez dat)
# Grafy pak: Rscript analyses/04_genic_methylation/genic_methylation_pca_diff.R
#
# Spusteni (pracovni adresar = methylation_profiles, conda env s bedtools, napr. envs/modkit.yaml):
#   bash analyses/04_genic_methylation/genic_methylation_preprocess.sh
# Hotove kroky se preskakuji (soubor existuje). Pro prepocet smazat prislusny soubor v OUT_DIR.

set -euo pipefail

# ---- Parametry ----
GTF="/mnt/data/ceitec_cfg2/710000-CEITEC/713000-cmm/713016-bioit/resources/references/homo_sapiens/GRCh38/annot/r115/GRCh38.gtf"
GENOME_FA="../data/hg38_with_chr.fa"           # reference s "chr" prefixem; vyzaduje .fai
RESULTS="results"                               # results/<sample>/methylation/<sample><MODKIT_SUFFIX>
MODKIT_SUFFIX="_modkit_sorted.bed"              # LC_ALL=C sort -k1,1 -k2,2n (rule sort_modkit_bed)
SAMPLES=(cfDNA_BAR cfDNA_MUD cfDNA_Pilot1 cfDNA_VN)

UPSTREAM=3000              # bp pred zacatkem genu
BIOTYPES="protein_coding"  # carkou oddeleny seznam, "all" = vsechny geny
MIN_COV=5                  # min. N_valid na CpG pozici (jeden retezec)
N_PARALLEL=2               # kolik vzorku zpracovavat soubezne (I/O narocne)

OUT_DIR="analyses/04_genic_methylation/output"
BEDTOOLS="${BEDTOOLS:-$(command -v bedtools || ls .snakemake/conda/*/bin/bedtools 2>/dev/null | head -n1 || true)}"

# ---- Kontroly ----
[[ -n "$BEDTOOLS" && -x "$BEDTOOLS" ]] || { echo "bedtools nenalezen - aktivuj conda env nebo BEDTOOLS=/cesta/bedtools" >&2; exit 1; }
[[ -f "$GTF" ]]              || { echo "Chybi GTF: $GTF" >&2; exit 1; }
[[ -f "$GENOME_FA.fai" ]]    || { echo "Chybi $GENOME_FA.fai (samtools faidx $GENOME_FA)" >&2; exit 1; }
for s in "${SAMPLES[@]}"; do
  f="$RESULTS/$s/methylation/${s}${MODKIT_SUFFIX}"
  [[ -f "$f" ]] || { echo "Chybi $f" >&2; exit 1; }
done
echo "bedtools: $BEDTOOLS"

mkdir -p "$OUT_DIR"
export LC_ALL=C
MAIN_CHR_RE='^chr([0-9]+|X|Y)$'

GENOME_TXT="$OUT_DIR/genome.txt"
GENES_BED="$OUT_DIR/genes_up${UPSTREAM}.bed"
CPG_BED="$OUT_DIR/cpg_sites_stranded.bed"

log() { echo "[$(date '+%F %T')] $*"; }

# ---- 0) Genome file (hlavni chromozomy, C-locale poradi = poradi v modkit_sorted.bed) ----
if [[ ! -s "$GENOME_TXT" ]]; then
  awk -v re="$MAIN_CHR_RE" -v OFS='\t' '$1 ~ re {print $1, $2}' "$GENOME_FA.fai" | sort -k1,1 > "$GENOME_TXT"
  log "genome.txt: $(wc -l < "$GENOME_TXT") chromozomu"
fi

# ---- 1) Genove oblasti z GTF ----
if [[ ! -s "$GENES_BED" ]]; then
  log "Genove oblasti z GTF ..."
  awk -F'\t' -v OFS='\t' -v up="$UPSTREAM" -v bt="$BIOTYPES" -v re="$MAIN_CHR_RE" '
    function attr(str, key) {
      if (match(str, key " \"[^\"]+\"")) return substr(str, RSTART + length(key) + 2, RLENGTH - length(key) - 3)
      return ""
    }
    BEGIN { n = split(bt, b, ","); for (i = 1; i <= n; i++) B[b[i]] = 1 }
    NR == FNR { len[$1] = $2; next }                  # genome.txt
    $3 != "gene" { next }
    {
      chr = ($1 ~ /^chr/) ? $1 : "chr" $1
      if (chr !~ re || !(chr in len)) next
      type = attr($9, "gene_biotype")
      if (bt != "all" && !(type in B)) next
      id = attr($9, "gene_id"); name = attr($9, "gene_name"); if (name == "") name = id
      # GTF 1-based vcetne -> BED 0-based polootevreny
      if ($7 == "+") { s = $4 - 1 - up; e = $5 } else { s = $4 - 1; e = $5 + up }
      if (s < 0) s = 0
      if (e > len[chr]) e = len[chr]
      print chr, s, e, id, name, $7, type
    }' "$GENOME_TXT" "$GTF" | sort -k1,1 -k2,2n -k3,3n > "$GENES_BED.tmp"
  mv "$GENES_BED.tmp" "$GENES_BED"
  log "Genovych oblasti: $(wc -l < "$GENES_BED")"
fi

# ---- 2) CpG pozice z reference ----
if [[ ! -s "$CPG_BED" ]]; then
  log "Hledam CpG v referenci (muze trvat ~10-20 min) ..."
  CPG_TMP="$OUT_DIR/cpg_tmp"
  rm -rf "$CPG_TMP"; mkdir -p "$CPG_TMP"
  # Pozice v chromozomu jsou vzestupne -> staci zapsat po chromozomech a spojit v C-locale poradi
  awk -v re="$MAIN_CHR_RE" -v dir="$CPG_TMP" -v OFS='\t' '
    function emit(p) {               # p = 1-based pozice C v CpG
      print chr, p - 1, p,     ".", 0, "+" > out   # modkit: + retezec hlasi C
      print chr, p,     p + 1, ".", 0, "-" > out   #         - retezec hlasi G
    }
    /^>/ {
      if (out != "") close(out)
      chr = substr($1, 2); use = (chr ~ re); pos = 0; prev = ""
      out = use ? dir "/" chr ".bed" : ""
      next
    }
    use {
      s = toupper($0); L = length(s)
      if (prev == "C" && substr(s, 1, 1) == "G") emit(pos)   # CpG pres zalomeni radku
      i = 1
      while ((j = index(substr(s, i), "CG")) > 0) { k = i + j - 1; emit(pos + k); i = k + 1 }
      pos += L; prev = substr(s, L, 1)
    }' "$GENOME_FA"
  cut -f1 "$GENOME_TXT" | while read -r c; do
    [[ -f "$CPG_TMP/$c.bed" ]] || { echo "V referenci chybi $c" >&2; exit 1; }
    cat "$CPG_TMP/$c.bed"
  done > "$CPG_BED.tmp"
  mv "$CPG_BED.tmp" "$CPG_BED"
  rm -rf "$CPG_TMP"
  log "CpG pozic (oba retezce): $(wc -l < "$CPG_BED")"
fi

# ---- 3) modkit -> geny, pro kazdy vzorek ----
aggregate_sample() {
  local s="$1"
  local inp="$RESULTS/$s/methylation/${s}${MODKIT_SUFFIX}"
  local out="$OUT_DIR/${s}_genes_up${UPSTREAM}_cov${MIN_COV}.tsv"
  if [[ -s "$out" ]]; then log "$s: hotovo drive, preskakuji"; return 0; fi
  # vlastni .tmp pro kazdy beh - stare .tmp (jiny vlastnik / prerusene behy) nevadi
  rm -f "$out".tmp* 2>/dev/null || true
  local tmp="$out.tmp.$$"
  log "$s: start"
  # b: chr start end . N_valid strand N_mod  -> map: sum N_valid, sum N_mod, pocet CpG volani
  awk -F'\t' -v OFS='\t' -v mc="$MIN_COV" -v re="$MAIN_CHR_RE" \
      '$4 == "m" && $10 >= mc && $1 ~ re {print $1, $2, $3, ".", $10, $6, $12}' "$inp" \
    | "$BEDTOOLS" intersect -a stdin -b "$CPG_BED" -s -u -sorted -g "$GENOME_TXT" \
    | "$BEDTOOLS" map -a "$GENES_BED" -b stdin -c 5,7,5 -o sum,sum,count -g "$GENOME_TXT" \
    > "$tmp"
  mv -f "$tmp" "$out"
  log "$s: hotovo -> $out"
}

log "Agreguji ${#SAMPLES[@]} vzorku ($N_PARALLEL soubezne) ..."
for s in "${SAMPLES[@]}"; do
  aggregate_sample "$s" &
  while (( $(jobs -rp | wc -l) >= N_PARALLEL )); do wait -n || true; done
done
wait || true

# Vystup vznika pres .tmp + mv, takze existujici soubor = uspesne dokonceny vzorek
missing=()
for s in "${SAMPLES[@]}"; do
  [[ -s "$OUT_DIR/${s}_genes_up${UPSTREAM}_cov${MIN_COV}.tsv" ]] || missing+=("$s")
done
if (( ${#missing[@]} )); then
  echo "CHYBA - nedokoncene vzorky: ${missing[*]}" >&2
  exit 1
fi
log "Hotovo. Dalsi krok: Rscript analyses/04_genic_methylation/genic_methylation_pca_diff.R"
