#!/usr/bin/env python3
"""
07_marker_methylation_table.py

Pro sadu vzorků (bedMethyl ze skriptu 02) spočítá průměrnou % metylaci
(váženou pokrytím) na PCNSL markerových regionech (BED ze skriptu 06)
a sestaví tabulku vzorky x markery - přímo srovnatelnou (po vydělení
100) s beta-values z EPIC array publikací.

Vyžaduje nástroj `bedtools` v PATH (bedtools intersect).

Použití:
  python 07_marker_methylation_table.py \\
      --markers pcnsl_markers.hg38.bed \\
      --out marker_methylation_matrix.tsv \\
      sample1=sample1.sorted.bedmethyl.gz \\
      sample2=sample2.sorted.bedmethyl.gz \\
      ...
"""
import argparse
import gzip
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import pandas as pd

# Sloupce bedMethyl (modkit pileup), 0-based indexy - viz README/skript 02.
# Pokud vaše verze modkitu má jiné pořadí, upravte zde.
COL_CHROM = 0
COL_START = 1
COL_END = 2
COL_VALID_COV = 9
COL_PCT_MOD = 10


def parse_sample_args(pairs):
    samples = {}
    for p in pairs:
        if "=" not in p:
            raise SystemExit(f"[CHYBA] Očekávám 'jmeno=cesta_k_bedmethyl.gz', dostal jsem: {p}")
        name, path = p.split("=", 1)
        samples[name] = path
    return samples


def intersect_bedmethyl(bedmethyl_gz, markers_bed, tmp_dir):
    """bedtools intersect markerů s bedMethyl -> DataFrame se sloupci
    marker_label, valid_cov, pct_mod."""
    tmp_bm = Path(tmp_dir) / "tmp.bedmethyl"
    with gzip.open(bedmethyl_gz, "rt") as fin, open(tmp_bm, "w") as fout:
        shutil.copyfileobj(fin, fout)

    cmd = ["bedtools", "intersect", "-a", markers_bed, "-b", str(tmp_bm), "-wa", "-wb"]
    result = subprocess.run(cmd, capture_output=True, text=True, check=True)

    rows = []
    for line in result.stdout.splitlines():
        f = line.split("\t")
        # -a (marker BED): chrom,start,end,label (4 sloupce)
        # -b (bedmethyl): dalších >=11 sloupců, offset 4
        label = f[3]
        valid_cov = float(f[4 + COL_VALID_COV])
        pct_mod = float(f[4 + COL_PCT_MOD])
        rows.append((label, valid_cov, pct_mod))

    return pd.DataFrame(rows, columns=["marker", "valid_cov", "pct_mod"])


def weighted_mean_per_marker(df):
    """Váženě zprůměruje % metylace přes všechny CpG pozice patřící
    danému markeru (váha = valid coverage dané pozice)."""
    def agg(g):
        w = g["valid_cov"]
        if w.sum() == 0:
            return pd.Series({"pct_methylation": float("nan"), "total_coverage": 0})
        weighted = (g["pct_mod"] * w).sum() / w.sum()
        return pd.Series({"pct_methylation": weighted, "total_coverage": w.sum()})

    return df.groupby("marker").apply(agg).reset_index()


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--markers", required=True, help="BED se markerovými regiony (skript 06)")
    ap.add_argument("--out", required=True, help="výstupní TSV matice vzorky x markery")
    ap.add_argument("samples", nargs="+", help="jmeno=cesta_k_bedmethyl.gz, opakovaně pro každý vzorek")
    args = ap.parse_args()

    if shutil.which("bedtools") is None:
        raise SystemExit("[CHYBA] 'bedtools' nenalezen v PATH (conda install -c bioconda bedtools).")

    samples = parse_sample_args(args.samples)

    pct_table = {}
    cov_table = {}
    with tempfile.TemporaryDirectory() as tmp_dir:
        for name, path in samples.items():
            print(f"[INFO] Zpracovávám vzorek {name} ({path})...", file=sys.stderr)
            df = intersect_bedmethyl(path, args.markers, tmp_dir)
            if df.empty:
                print(f"[VAROVÁNÍ] Vzorek {name}: žádný průnik s markery - "
                      f"zkontrolujte chromozomová jména (chr1 vs 1) a genome build.",
                      file=sys.stderr)
                continue
            agg = weighted_mean_per_marker(df)
            pct_table[name] = agg.set_index("marker")["pct_methylation"]
            cov_table[name] = agg.set_index("marker")["total_coverage"]

    pct_df = pd.DataFrame(pct_table)
    pct_df["mean_beta_equivalent"] = pct_df.mean(axis=1) / 100.0
    pct_df.to_csv(args.out, sep="\t")

    cov_out = str(Path(args.out).with_suffix("")) + ".coverage.tsv"
    pd.DataFrame(cov_table).to_csv(cov_out, sep="\t")

    print(f"[OK] Metylační matice: {args.out}", file=sys.stderr)
    print(f"[OK] Tabulka pokrytí:  {cov_out}", file=sys.stderr)
    print("[POZOR] U vzorků s nízkým total_coverage (pár readů na marker) "
          "berte pct_methylation jen orientačně - u n=4 nádorových vzorků "
          "nejde dělat statistické testy, jde o deskriptivní srovnání.",
          file=sys.stderr)


if __name__ == "__main__":
    main()
