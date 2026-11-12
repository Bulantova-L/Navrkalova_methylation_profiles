#!/usr/bin/env bash
# Odhad podilu nadorove DNA (ctDNA) v cfDNA z kopiovych zmen (CNV) - ichorCNA
#
# Princip: nadorova DNA nese zisky/ztraty chromozomovych ramen, DNA z normalnich krvinek je plocha.
# ichorCNA spocita ready v 1Mb binech, opravi GC a mapovatelnost, normalizuje proti zdravym vzorkum
# (panel of normals, PoN) a z amplitudy CNV odhadne tumor fraction (TF). Nezavisi na metylacnich
# volanich ani na basecalling modifikacnim modelu.
#
# Kroky:
#   1) samtools bedcov -c: BAM -> WIG (pocet readu s MAPQ >= MIN_MAPQ v binech BIN_SIZE)
#      (misto HMMcopy readCounter - jeho bamtools parser pada na @PG DS: tagu z dorado)
#   2) (volitelne) vlastni PoN ze zdravych vzorku (CONTROLS)
#   3) runIchorCNA.R pro kazdy vzorek se zabudovanym hg38 PoN (zdravi darci, Illumina ULP-WGS)
#      a pokud existuje vlastni PoN, tak i s nim (srovnani)
#   4) tabulka TF + ploidie pro vsechny vzorky
#
# Rezim kratkych fragmentu (in silico size selection, napr. 90-150 bp):
#   ctDNA je kratsi nez cfDNA z krvinek (~167 bp) -> pocitat jen kratke ready obohati nadorovy signal.
#   U ONT read typicky pokryva cely fragment; delka fragmentu = delka zarovnani na referenci (z CIGAR).
#   FRAG_MIN=90 FRAG_MAX=150 bash analyses/05_tumor_fraction_ichorcna/tumor_fraction_ichorcna.sh
#   -> vystup v analyses/05_tumor_fraction_ichorcna/output/ichorcna_frag90-150/ (vlastni PoN se stavi z kratkych fragmentu kontrol),
#      navic fraglen/<sample>.hist.tsv (distribuce delek) a fragment_summary.tsv (podil kratkych readu)
#
# Priprava (jednou):
#   conda env create -f envs/ichorcna.yaml
#   git clone https://github.com/broadinstitute/ichorCNA.git ../tools/ichorCNA   # skripty + referencni soubory
# Spusteni (pracovni adresar = methylation_profiles):
#   conda activate ichorcna
#   bash analyses/05_tumor_fraction_ichorcna/tumor_fraction_ichorcna.sh
#   s dalsim hotovym PoN (.rds, 1 Mb biny, hg38, UCSC nazvy chromozomu):
#   (externi PoN nejdriv prevest: Rscript analyses/05_tumor_fraction_ichorcna/convert_pon_to_ichorcna.R ...)
#   EXTRA_PON=analyses/05_tumor_fraction_ichorcna/pon/ONT_PoN_1control_1Mb.GR.rds EXTRA_PON_NAME=pon_ONT_PoN_1control_1Mb bash analyses/05_tumor_fraction_ichorcna/tumor_fraction_ichorcna.sh
#   -> vystup v analyses/05_tumor_fraction_ichorcna/output/ichorcna/pon_ONT_PoN_1control_1Mb/<sample>/
#
# Vystup: analyses/05_tumor_fraction_ichorcna/output/ichorcna/<pon>/<sample>/ (<sample>.params.txt, genome-wide grafy *_genomeWide.pdf)
#         analyses/05_tumor_fraction_ichorcna/output/ichorcna/tumor_fraction_summary.tsv

set -euo pipefail

# ---- Parametry ----
ICHOR_DIR="../tools/ichorCNA"                  # klon https://github.com/broadinstitute/ichorCNA
RESULTS="results"
SAMPLES=(cfDNA_BAR cfDNA_LOS cfDNA_MUD cfDNA_Pilot1 cfDNA_VN)
CONTROLS=(cfDNA_VN)                            # zdrave kontroly pro vlastni PoN (cim vic, tim lip)
BIN_SIZE=1000000                               # 1 Mb - jen pro tuto velikost jsou v ichorCNA hg38 GC/map soubory
MIN_MAPQ=20
N_PARALLEL=2                                   # kolik BAM pocitat soubezne (I/O narocne)
FRAG_MIN="${FRAG_MIN:-}"                       # prazdne = vsechny ready; jinak delka fragmentu (bp) od-do
FRAG_MAX="${FRAG_MAX:-}"
EXTRA_PON="${EXTRA_PON:-}"                     # volitelne: dalsi hotovy PoN (.rds), napr. analyses/05_tumor_fraction_ichorcna/pon/ONT_PoN_1control_1Mb.GR.rds
EXTRA_PON_NAME="${EXTRA_PON_NAME:-}"           # nazev podadresare vystupu; prazdne = podle jmena souboru
if [[ -n "$EXTRA_PON" ]]; then
  [[ -s "$EXTRA_PON" ]] || { echo "Chybi EXTRA_PON: $EXTRA_PON" >&2; exit 1; }
  [[ "$EXTRA_PON" == *.rds ]] || { echo "EXTRA_PON musi byt .rds (ichorCNA cte PoN pres readRDS)" >&2; exit 1; }
  [[ -n "$EXTRA_PON_NAME" ]] || EXTRA_PON_NAME="pon_$(basename "$EXTRA_PON" .rds)"
fi
if [[ -n "$FRAG_MIN" || -n "$FRAG_MAX" ]]; then
  [[ "$FRAG_MIN" =~ ^[0-9]+$ && "$FRAG_MAX" =~ ^[0-9]+$ ]] || { echo "Zadej obe FRAG_MIN a FRAG_MAX (cela cisla)" >&2; exit 1; }
  OUT_DIR="analyses/05_tumor_fraction_ichorcna/output/ichorcna_frag${FRAG_MIN}-${FRAG_MAX}"
else
  OUT_DIR="analyses/05_tumor_fraction_ichorcna/output/ichorcna"
fi

EXT="$ICHOR_DIR/inst/extdata"
GC_WIG="$EXT/gc_hg38_1000kb.wig"
MAP_WIG="$EXT/map_hg38_1000kb.wig"
CENTROMERE="$EXT/GRCh38.GCA_000001405.2_centromere_acen.txt"
BUILTIN_PON="$EXT/HD_ULP_PoN_hg38_1Mb_median_normAutosome_median.rds"
RUN_ICHOR="$ICHOR_DIR/scripts/runIchorCNA.R"
CREATE_PON="$ICHOR_DIR/scripts/createPanelOfNormals.R"

log() { echo "[$(date '+%F %T')] $*"; }

# ---- Kontroly ----
for t in samtools Rscript; do command -v "$t" >/dev/null || { echo "Chybi $t - conda activate ichorcna" >&2; exit 1; }; done
for f in "$GC_WIG" "$MAP_WIG" "$CENTROMERE" "$BUILTIN_PON" "$RUN_ICHOR" "$CREATE_PON"; do
  [[ -f "$f" ]] || { echo "Chybi $f - naklonuj ichorCNA do $ICHOR_DIR" >&2; exit 1; }
done
bam_of() { echo "$RESULTS/$1/methylation/$1.sorted.bam"; }
for s in "${SAMPLES[@]}"; do
  b="$(bam_of "$s")"
  [[ -f "$b" ]] || { echo "Chybi $b" >&2; exit 1; }
  [[ -f "$b.bai" || -f "${b%.bam}.bai" ]] || { echo "Chybi index $b.bai (samtools index $b)" >&2; exit 1; }
done
{ samtools bedcov 2>&1 || true; } | grep -q -- " -c " || { echo "samtools bedcov neumi -c (potreba samtools >= 1.11)" >&2; exit 1; }
mkdir -p "$OUT_DIR/wig"

# ---- 1) BAM -> WIG ----
# Biny presne podle GC wig z ichorCNA (stejny pocet a poloha binu), orezane na delku chromozomu z BAM
BINS_BED="$OUT_DIR/bins_${BIN_SIZE}.bed"
if [[ ! -s "$BINS_BED" ]]; then
  samtools view -H "$(bam_of "${SAMPLES[0]}")" \
    | awk -F'\t' '$1 == "@SQ" { for (i = 2; i <= NF; i++) { split($i, kv, ":"); if (kv[1] == "SN") sn = kv[2]; if (kv[1] == "LN") ln = kv[2] } print sn "\t" ln }' \
    > "$OUT_DIR/chrom_sizes.txt"
  awk -v OFS='\t' '
    NR == FNR { len[$1] = $2; next }
    /^fixedStep/ {
      for (i = 2; i <= NF; i++) { split($i, kv, "="); p[kv[1]] = kv[2] }
      chr = (p["chrom"] ~ /^chr/) ? p["chrom"] : "chr" p["chrom"]
      start = p["start"] - 1; step = p["step"]; k = 0
      use = (chr ~ /^chr([0-9]+|X|Y)$/) && (chr in len)
      next
    }
    use {
      s = start + k * step; e = s + step; k++
      if (s >= len[chr]) next
      if (e > len[chr]) e = len[chr]
      print chr, s, e
    }' "$OUT_DIR/chrom_sizes.txt" "$GC_WIG" > "$BINS_BED.tmp.$$"
  mv -f "$BINS_BED.tmp.$$" "$BINS_BED"
  log "Binu: $(wc -l < "$BINS_BED")"
fi

count_bins() {
  local s="$1" wig="$OUT_DIR/wig/$1.${BIN_SIZE}.wig"
  [[ -s "$wig" ]] && { log "$s: WIG hotovo drive"; return 0; }
  log "$s: pocitam ready v binech ..."
  if [[ -z "$FRAG_MIN" ]]; then
    # -c: pocet readu v binu (posledni sloupec); -G 0xF04: bez unmapped, secondary, QC-fail, duplikatu, supplementary
    samtools bedcov -c -Q "$MIN_MAPQ" -G 0xF04 "$BINS_BED" "$(bam_of "$s")" \
      | awk -v step="$BIN_SIZE" '
          $1 != chr { chr = $1; printf "fixedStep chrom=%s start=%d step=%d span=%d\n", chr, $2 + 1, step, step }
          { print $NF }' > "$wig.tmp.$$"
  else
    # Delka fragmentu = delka zarovnani na referenci (CIGAR operace M, D, N, =, X).
    # Read se pocita do binu podle zacatku zarovnani; zaroven se uklada histogram delek vsech readu.
    mkdir -p "$OUT_DIR/fraglen"
    samtools view -q "$MIN_MAPQ" -F 0xF04 "$(bam_of "$s")" \
      | cut -f3,4,6 \
      | awk -F'\t' -v lo="$FRAG_MIN" -v hi="$FRAG_MAX" -v bs="$BIN_SIZE" -v bins="$BINS_BED" \
            -v hist="$OUT_DIR/fraglen/$s.hist.tsv.tmp" '
          $1 !~ /^chr([0-9]+|X|Y)$/ || $3 == "*" { next }
          {
            n = split($3, num, /[MIDNSHP=X]/); split($3, op, /[0-9]+/)   # op[i+1] patri k num[i]
            span = 0
            for (i = 1; i < n; i++) if (op[i + 1] ~ /^[MDN=X]$/) span += num[i]
            h[(span > 1000) ? 1000 : span]++
            if (span >= lo && span <= hi) cnt[$1 SUBSEP int(($2 - 1) / bs)]++
          }
          END {
            while ((getline line < bins) > 0) {
              split(line, f, "\t")
              if (f[1] != chr) { chr = f[1]; printf "fixedStep chrom=%s start=%d step=%d span=%d\n", chr, f[2] + 1, bs, bs }
              print cnt[chr SUBSEP int(f[2] / bs)] + 0
            }
            print "length\treads" > hist
            for (l = 0; l <= 1000; l++) if (l in h) print l "\t" h[l] > hist
          }' > "$wig.tmp.$$"
    mv -f "$OUT_DIR/fraglen/$s.hist.tsv.tmp" "$OUT_DIR/fraglen/$s.hist.tsv"
  fi
  mv -f "$wig.tmp.$$" "$wig"
  log "$s: WIG hotovo"
}

for s in "${SAMPLES[@]}"; do
  count_bins "$s" &
  while (( $(jobs -rp | wc -l) >= N_PARALLEL )); do wait -n || true; done
done
wait || true
for s in "${SAMPLES[@]}"; do
  [[ -s "$OUT_DIR/wig/$s.${BIN_SIZE}.wig" ]] || { echo "CHYBA: chybi WIG pro $s" >&2; exit 1; }
done

if [[ -n "$FRAG_MIN" ]]; then
  # Podil kratkych readu a median delky - vyssi podil kratkych ~ vice ctDNA (nezavisle na CNV)
  FRAG_SUMMARY="$OUT_DIR/fragment_summary.tsv"
  printf 'sample\treads_total\treads_%s_%s\tfraction_short\tmedian_length\tfraction_100_150\tfraction_151_220\n' "$FRAG_MIN" "$FRAG_MAX" > "$FRAG_SUMMARY"
  for s in "${SAMPLES[@]}"; do
    awk -F'\t' -v s="$s" -v lo="$FRAG_MIN" -v hi="$FRAG_MAX" '
      NR == 1 { next }
      { L[NR] = $1; C[NR] = $2; tot += $2
        if ($1 >= lo && $1 <= hi) sh += $2
        if ($1 >= 100 && $1 <= 150) a += $2
        if ($1 >= 151 && $1 <= 220) b += $2 }
      END {
        for (i = 2; i <= NR; i++) { cum += C[i]; if (cum >= tot / 2) { med = L[i]; break } }
        printf "%s\t%d\t%d\t%.4f\t%d\t%.4f\t%.4f\n", s, tot, sh, sh / tot, med, a / tot, b / tot
      }' "$OUT_DIR/fraglen/$s.hist.tsv" >> "$FRAG_SUMMARY"
  done
  log "Fragmenty:"; cat "$FRAG_SUMMARY"
fi

# ---- 2) Vlastni PoN ze zdravych kontrol ----
OWN_PON_PREFIX="$OUT_DIR/pon_own"
OWN_PON="${OWN_PON_PREFIX}_median.rds"
if [[ ! -s "$OWN_PON" ]]; then
  printf '%s\n' "${CONTROLS[@]/#/$OUT_DIR/wig/}" | sed "s/\$/.${BIN_SIZE}.wig/" > "$OUT_DIR/pon_own_files.txt"
  log "Vlastni PoN z: ${CONTROLS[*]}"
  Rscript "$CREATE_PON" --filelist "$OUT_DIR/pon_own_files.txt" \
    --gcWig "$GC_WIG" --mapWig "$MAP_WIG" --centromere "$CENTROMERE" \
    --chrs "c(1:22, 'X')" \
    --genomeStyle UCSC --method median --outfile "$OWN_PON_PREFIX" \
    || log "VAROVANI: vlastni PoN se nepodarilo vytvorit - pokracuji jen se zabudovanym"
fi

# ---- 3) ichorCNA ----
# Nastaveni pro nizke TF (doporuceni autoru ichorCNA pro cfDNA s TF < 5-10 %)
# --libdir: pouzit R kod z klonovaneho repa (stejna verze jako runIchorCNA.R); conda r-ichorcna je starsi
# a jeho HMMsegment/runEM s timto skriptem pada ("lambdas[, , i] ... incorrect number of subscripts")
run_ichor() {
  local s="$1" pon_name="$2" pon_file="$3"
  local out="$OUT_DIR/$pon_name/$s"
  [[ -s "$out/$s.params.txt" ]] && { log "$s [$pon_name]: hotovo drive"; return 0; }
  mkdir -p "$out"
  log "$s [$pon_name]: ichorCNA ..."
  Rscript "$RUN_ICHOR" --id "$s" \
    --WIG "$OUT_DIR/wig/$s.${BIN_SIZE}.wig" \
    --gcWig "$GC_WIG" --mapWig "$MAP_WIG" --centromere "$CENTROMERE" \
    --normalPanel "$pon_file" \
    --ploidy "c(2,3)" \
    --normal "c(0.5,0.6,0.7,0.8,0.9,0.95,0.99,0.995,0.999)" \
    --maxCN 5 --includeHOMD False \
    --estimateNormal True --estimatePloidy True \
    --estimateScPrevalence False --scStates "c()" \
    --chrs "c(1:22, 'X')" --chrTrain "c(1:22)" \
    --txnE 0.9999 --txnStrength 10000 \
    --genomeStyle UCSC --genomeBuild hg38 \
    --libdir "$ICHOR_DIR" \
    --outDir "$out" > "$out/$s.ichorcna.log" 2>&1 \
    || { log "$s [$pon_name]: CHYBA, viz $out/$s.ichorcna.log"; return 1; }
}

for s in "${SAMPLES[@]}"; do
  run_ichor "$s" pon_builtin "$BUILTIN_PON" || true
  if [[ -s "$OWN_PON" ]]; then
    # kontrola v PoN by se normalizovala sama proti sobe -> pro ni jen zabudovany PoN
    [[ " ${CONTROLS[*]} " == *" $s "* ]] || run_ichor "$s" pon_own "$OWN_PON" || true
  fi
  # EXTRA_PON: kontroly se spousti taky (nevime, z ceho byl PoN postaven); pokud obsahuje kontrolu, jeji TF ~ 0 neni informativni
  [[ -z "$EXTRA_PON" ]] || run_ichor "$s" "$EXTRA_PON_NAME" "$EXTRA_PON" || true
done

# ---- 4) Souhrn ----
SUMMARY="$OUT_DIR/tumor_fraction_summary.tsv"
printf 'sample\tpon\ttumor_fraction\tploidy\tsubclone_fraction\tgc_map_correction_MAD\n' > "$SUMMARY"
for pon in pon_builtin pon_own ${EXTRA_PON_NAME:+"$EXTRA_PON_NAME"}; do
  for s in "${SAMPLES[@]}"; do
    p="$OUT_DIR/$pon/$s/$s.params.txt"
    [[ -s "$p" ]] || continue
    # 2. radek params.txt: Sample Tumor_Fraction Ploidy ...
    awk -F'\t' -v s="$s" -v pon="$pon" '
      NR == 2 { tf = $2; pl = $3 }
      /^Subclone Fraction:/ { split($0, a, ":"); sc = a[2] }
      /^GC-Map correction MAD:/ { split($0, a, ":"); mad = a[2] }
      END { gsub(/^[ \t]+/, "", sc); gsub(/^[ \t]+/, "", mad); print s, pon, tf, pl, sc, mad }' OFS='\t' "$p" >> "$SUMMARY"
  done
done
log "Hotovo:"; cat "$SUMMARY"
log "Genome-wide grafy: $OUT_DIR/<pon>/<sample>/<sample>/<sample>_genomeWide.pdf"
