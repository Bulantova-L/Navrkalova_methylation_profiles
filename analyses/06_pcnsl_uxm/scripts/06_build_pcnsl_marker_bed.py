#!/usr/bin/env python3
"""
06_build_pcnsl_marker_bed.py

Namapuje kurátorský seznam 450K/EPIC CpG probe ID (viz
markers/pcnsl_markers_450k.tsv) na hg38 genomické souřadnice a zapíše
BED soubor markerových regionů (probe pozice +/- flank).

*** ZÁSADNÍ INTERPRETAČNÍ CAVEAT ***
Tyto markery (Clinical Epigenetics 2021, 13:172, DOI:
10.1186/s13148-021-01091-9, "Methylated markers accurately distinguish
primary central nervous system lymphomas (PCNSL) from other CNS
tumors") byly vybrány, aby odlišily PRIMÁRNÍ CNS lymfom od gliomů
(GBM/LGG) - tedy DIAGNOSTICKÁ úloha. NENÍ to validovaný prediktor
rizika CNS relapsu u systémového DLBCL, což je vaše skutečná otázka
(bod 2 zadání). Použijte tento panel jako explorativní výchozí bod
(je to nejkonkrétnější publikovaný CpG-level signál, co jsem při
rešerši našel), ne jako hotové řešení - ideálně doplnit o vlastní
literární rešerši k CNS-homing/BBB-transit/imunitní únik genům
zmíněným kolegyní (např. CXCR4, PD-L1/PD-L2 lokus 9p24.1, MGMT,
DAPK1 promotorová hypermetylace), pro které genomické souřadnice
tento skript zatím NEobsahuje - je potřeba je doplnit stejným
způsobem (probe ID nebo přímo souřadnice genu/promotoru).

Souřadnice se řeší přes hg38-remapovaný Infinium manifest (Zhou lab,
https://zwdzwd.github.io/InfiniumAnnotation) - tato varianta je
doporučovaná i GDC/TCGA namísto původního Illumina hg19 manifestu.
URL níže jsou KANDIDÁTI - ověřte, že se dají stáhnout, případně
aktualizujte podle aktuální struktury webu Zhou lab.

Použití:
  python 06_build_pcnsl_marker_bed.py \\
      --markers markers/pcnsl_markers_450k.tsv \\
      --out pcnsl_markers.hg38.bed \\
      --flank 250
"""
import argparse
import csv
import gzip
import sys
import urllib.error
import urllib.request

MANIFEST_URL_CANDIDATES = [
    "https://zwdzwd.github.io/InfiniumAnnotation/current/EPIC/EPIC.hg38.manifest.tsv.gz",
    "https://zwdzwd.github.io/InfiniumAnnotation/current/HM450/HM450.hg38.manifest.tsv.gz",
]


def download_manifest(dest_path):
    last_err = None
    for url in MANIFEST_URL_CANDIDATES:
        try:
            print(f"[INFO] Zkouším stáhnout manifest: {url}", file=sys.stderr)
            urllib.request.urlretrieve(url, dest_path)
            print(f"[OK] Manifest stažen do {dest_path}", file=sys.stderr)
            return dest_path
        except (urllib.error.URLError, urllib.error.HTTPError) as e:
            print(f"[VAROVÁNÍ] Nepodařilo se ({e}), zkouším další URL...", file=sys.stderr)
            last_err = e
    raise SystemExit(
        "[CHYBA] Nepodařilo se stáhnout žádný kandidátní manifest.\n"
        f"        Poslední chyba: {last_err}\n"
        "        Ruční alternativy:\n"
        "        - https://zwdzwd.github.io/InfiniumAnnotation (najděte aktuální odkaz)\n"
        "        - Bioconductor: IlluminaHumanMethylationEPICanno.ilm10b5.hg38 (R)\n"
        "        - GDC probe annotation: https://gdc.cancer.gov (improved DNA methylation "
        "array probe annotation)\n"
        "        Stažený soubor pak předejte přes --manifest <cesta>."
    )


def load_manifest_coords(manifest_path, wanted_probes):
    """Vrátí dict probe_id -> (chrom, start, end). Očekává sloupce typu
    probeID/Probe_ID a CpG_chrm/CpG_beg/CpG_end (Zhou lab layout) -
    pokud se liší, uprav COLUMN_CANDIDATES."""
    coords = {}
    opener = gzip.open if manifest_path.endswith(".gz") else open
    with opener(manifest_path, "rt") as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader)
        header_lower = [h.lower() for h in header]

        def find_col(candidates):
            for c in candidates:
                if c.lower() in header_lower:
                    return header_lower.index(c.lower())
            return None

        i_probe = find_col(["probeid", "probe_id", "name", "ilmnid"])
        i_chrom = find_col(["cpg_chrm", "chr", "chrom", "seqnames"])
        i_start = find_col(["cpg_beg", "start", "pos"])
        i_end = find_col(["cpg_end", "end"])

        if i_probe is None or i_chrom is None or i_start is None:
            raise SystemExit(
                f"[CHYBA] Nerozpoznal jsem sloupce manifestu: {header}\n"
                "        Uprav load_manifest_coords() v 06_build_pcnsl_marker_bed.py."
            )

        for row in reader:
            if not row:
                continue
            probe = row[i_probe]
            if probe in wanted_probes:
                chrom = row[i_chrom]
                start = int(row[i_start])
                end = int(row[i_end]) if i_end is not None and row[i_end] else start + 1
                coords[probe] = (chrom, start, end)
            if len(coords) == len(wanted_probes):
                break
    return coords


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--markers", required=True, help="TSV se sloupci probe_id, gene, note")
    ap.add_argument("--manifest", default=None, help="lokální cesta k manifestu (jinak se zkusí stáhnout)")
    ap.add_argument("--manifest-cache", default="manifest.hg38.tsv.gz", help="kam manifest uložit, pokud se stahuje")
    ap.add_argument("--out", required=True, help="výstupní BED (hg38)")
    ap.add_argument("--flank", type=int, default=250, help="+/- bp kolem probe pozice (default 250)")
    args = ap.parse_args()

    markers = []
    with open(args.markers) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            markers.append(row)
    wanted = {m["probe_id"] for m in markers}
    print(f"[INFO] Načteno {len(wanted)} probe ID z {args.markers}", file=sys.stderr)

    manifest_path = args.manifest or download_manifest(args.manifest_cache)
    coords = load_manifest_coords(manifest_path, wanted)

    missing = wanted - set(coords)
    if missing:
        print(f"[VAROVÁNÍ] Pro {len(missing)} probe(s) se nenašly souřadnice: "
              f"{sorted(missing)} - zkontrolujte manifest/probe ID ručně "
              f"(např. https://www.ncbi.nlm.nih.gov/probe/ nebo Illumina manifest).",
              file=sys.stderr)

    rows = []
    for m in markers:
        probe = m["probe_id"]
        if probe not in coords:
            continue
        chrom, start, end = coords[probe]
        b_start = max(0, start - args.flank)
        b_end = end + args.flank
        label = f"{probe}|{m.get('gene', 'NA')}"
        rows.append((chrom, b_start, b_end, label))

    rows.sort(key=lambda r: (r[0], r[1]))
    with open(args.out, "w") as fh:
        for chrom, start, end, label in rows:
            fh.write(f"{chrom}\t{start}\t{end}\t{label}\n")

    print(f"[OK] Zapsáno {len(rows)} markerových regionů do {args.out}", file=sys.stderr)
    print("[DALŠÍ KROK] Doplňte manuálně regiony pro CXCR4, PD-L1/PD-L2 (9p24.1), "
          "MGMT, DAPK1 promotor apod. - najděte hg38 souřadnice promotoru/TSS "
          "(např. UCSC Genome Browser / Ensembl) a přidejte řádky do BED ručně, "
          "nebo je nejprve zapište do markers/pcnsl_markers_450k.tsv jako "
          "probe_id (pokud existuje odpovídající array sonda) a skript spusťte znovu.",
          file=sys.stderr)


if __name__ == "__main__":
    main()
