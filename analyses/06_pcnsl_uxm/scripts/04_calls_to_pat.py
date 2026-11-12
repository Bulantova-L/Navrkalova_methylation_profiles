#!/usr/bin/env python3
"""
04_calls_to_pat.py

Konvertuje per-read/per-CpG volání z `modkit extract calls` do PAT formátu
(pat.gz), který používá wgbs_tools / UXM_deconv pro fragment-level
dekonvoluci (Loyfer 2023).

DŮLEŽITÉ - toto NENÍ oficiální konvertor od nloyfer/wgbs_tools ani ONT.
Modkit i wgbs_tools jsou aktivně vyvíjené nástroje a jejich přesné
sloupce/formáty se mění mezi verzemi. Než tento skript pustíte na celý
dataset:

  1. Spusťte ho na malém testovacím regionu (např. jeden gen, pár tisíc
     readů) a výstupní pat.gz zkontrolujte příkazem `wgbstools view`
     nebo `wgbstools vis` - ověřte, že se vzory (C/T řetězce) a CpG
     indexy shodují s tím, co byste čekali z IGV/bedMethyl na stejném
     regionu.
  2. Ověřte hlavičku vašeho `modkit extract calls` výstupu
     (`zcat calls.tsv.gz | head -1`) - pokud se jména sloupců liší od
     očekávaných níže, upravte COLUMN_ALIASES.
  3. Ověřte, že váš `wgbstools init_genome` výstup skutečně obsahuje
     genomickou pozici <-> CpG-index mapování ve formátu, který zde
     předpokládáme (BED: chrom, start, end, cpg_index). Pokud se liší,
     uprav funkci `load_cpg_index`.

Princip PAT formátu (wgbs_tools): každý řádek = jeden unikátní
"metylační fragment", tj. (chromozom, index první CpG v readu,
řetězec vzoru přes po sobě jdoucí CpG, počet readů se stejným vzorem
na stejné pozici). Znak vzoru: 'C' = metylovaná CpG, 'T' = nemetylovaná
CpG (konvence odvozená z bisulfitové terminologie, přebíraná i pro
long-read data - viz --meth-char/--unmeth-char níže, pokud vaše verze
UXM/wgbstools očekává jiné znaky, např. 'M'/'U').

Použití:
  python 04_calls_to_pat.py \\
      --calls sample1.calls.tsv.gz \\
      --cpg-index reference/hg38/wgbstools_ref/CpG.bed.gz \\
      --out sample1.pat.gz \\
      --min-prob 0.8 \\
      --mod-code m
"""
import argparse
import gzip
import sys
from collections import defaultdict

import numpy as np
import pandas as pd

# Možná jména sloupců v modkit extract-calls výstupu (mění se mezi verzemi
# modkitu) -> interní jednotné jméno, které skript dál používá.
COLUMN_ALIASES = {
    "read_id": ["read_id", "read_name"],
    "chrom": ["chrom", "chromosome", "contig", "ref_chrom"],
    "ref_position": ["ref_position", "reference_position", "pos", "ref_pos"],
    "call_code": ["call_code", "mod_code", "modification_code", "base_mod_code"],
    "call_prob": ["call_prob", "mod_qual", "modification_probability", "mod_prob"],
}


def resolve_columns(df_columns):
    resolved = {}
    for canonical, aliases in COLUMN_ALIASES.items():
        found = next((a for a in aliases if a in df_columns), None)
        if found is None:
            raise SystemExit(
                f"[CHYBA] Nenašel jsem sloupec pro '{canonical}' mezi {list(df_columns)}.\n"
                f"        Upravte COLUMN_ALIASES v 04_calls_to_pat.py podle skutečné hlavičky."
            )
        resolved[canonical] = found
    return resolved


def load_cpg_index(cpg_index_path):
    """
    Načte CpG index vytvořený `wgbstools init_genome` (očekávaný formát:
    BED-like soubor se sloupci chrom, start, end, cpg_index - 1-based,
    sekvenční napříč genomem). Vrací dict: chrom -> (pozice[], indexy[])
    seřazené vzestupně, pro rychlé vyhledávání přes np.searchsorted.

    Pokud váš skutečný CpG.bed.gz má jiný sloupcový layout, uprav tuto
    funkci - je to jediné místo, kam je potřeba sáhnout.
    """
    print(f"[INFO] Načítám CpG index: {cpg_index_path}", file=sys.stderr)
    chrom_pos = defaultdict(list)
    chrom_idx = defaultdict(list)
    opener = gzip.open if cpg_index_path.endswith(".gz") else open
    with opener(cpg_index_path, "rt") as fh:
        for i, line in enumerate(fh):
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            chrom = parts[0]
            start = int(parts[1])
            # cpg_index = 4. sloupec, pokud existuje, jinak pořadové číslo
            # v rámci souboru (fallback - MÉNĚ SPOLEHLIVÉ, ověřte!)
            cpg_index = int(parts[3]) if len(parts) >= 4 and parts[3].isdigit() else i + 1
            chrom_pos[chrom].append(start)
            chrom_idx[chrom].append(cpg_index)

    index_map = {}
    for chrom in chrom_pos:
        pos_arr = np.array(chrom_pos[chrom], dtype=np.int64)
        idx_arr = np.array(chrom_idx[chrom], dtype=np.int64)
        order = np.argsort(pos_arr)
        index_map[chrom] = (pos_arr[order], idx_arr[order])
    print(f"[INFO] CpG index načten pro {len(index_map)} kontigů.", file=sys.stderr)
    return index_map


def positions_to_cpg_index(chrom, positions, index_map):
    if chrom not in index_map:
        return np.full(len(positions), -1, dtype=np.int64)
    pos_arr, idx_arr = index_map[chrom]
    where = np.searchsorted(pos_arr, positions)
    where = np.clip(where, 0, len(pos_arr) - 1)
    matched = pos_arr[where] == positions
    result = np.where(matched, idx_arr[where], -1)
    return result


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--calls", required=True, help="modkit extract calls TSV(.gz)")
    ap.add_argument("--cpg-index", required=True, help="CpG.bed.gz z wgbstools init_genome")
    ap.add_argument("--out", required=True, help="výstupní pat.gz")
    ap.add_argument("--min-prob", type=float, default=0.8, help="práh pravděpodobnosti volání (default 0.8)")
    ap.add_argument("--mod-code", default="m", help="modkit kód modifikace k zpracování (default 'm' = 5mC; '-' = kanonická)")
    ap.add_argument("--meth-char", default="C", help="znak pro metylovanou CpG v pattern stringu (default 'C', dle wgbs_tools konvence)")
    ap.add_argument("--unmeth-char", default="T", help="znak pro nemetylovanou CpG (default 'T')")
    ap.add_argument("--min-cpgs", type=int, default=1, help="min. počet CpG na read, aby byl zahrnut (default 1)")
    ap.add_argument("--chunksize", type=int, default=2_000_000, help="velikost chunku při čtení TSV (default 2e6 řádků)")
    args = ap.parse_args()

    index_map = load_cpg_index(args.cpg_index)

    # read_id -> list of (cpg_index, is_methylated)
    read_calls = defaultdict(list)

    n_rows = 0
    n_kept = 0
    n_unmapped = 0

    reader = pd.read_csv(args.calls, sep="\t", chunksize=args.chunksize, dtype=str)
    cols = None
    for chunk in reader:
        if cols is None:
            cols = resolve_columns(chunk.columns)
        n_rows += len(chunk)

        chunk = chunk.rename(columns={v: k for k, v in cols.items()})
        chunk["call_prob"] = pd.to_numeric(chunk["call_prob"], errors="coerce")
        chunk["ref_position"] = pd.to_numeric(chunk["ref_position"], errors="coerce")

        # ponecháme jen volání s dostatečnou jistotou a buď modifikovaný
        # kód, který nás zajímá (--mod-code), nebo kanonické volání ('-')
        keep = (chunk["call_prob"] >= args.min_prob) & (
            (chunk["call_code"] == args.mod_code) | (chunk["call_code"] == "-")
        )
        chunk = chunk[keep].dropna(subset=["ref_position"])
        if chunk.empty:
            continue

        for chrom, sub in chunk.groupby("chrom"):
            positions = sub["ref_position"].astype(np.int64).to_numpy()
            cpg_idx = positions_to_cpg_index(chrom, positions, index_map)
            sub = sub.assign(cpg_index=cpg_idx)
            unmapped = sub["cpg_index"] < 0
            n_unmapped += int(unmapped.sum())
            sub = sub[~unmapped]
            n_kept += len(sub)

            for read_id, code, cpg_i in zip(sub["read_id"], sub["call_code"], sub["cpg_index"]):
                is_meth = code == args.mod_code
                read_calls[read_id].append((int(cpg_i), is_meth))

    if n_rows == 0:
        raise SystemExit("[CHYBA] Vstupní soubor je prázdný nebo se nepodařilo načíst.")

    pct_unmapped = 100 * n_unmapped / max(1, n_unmapped + n_kept)
    print(f"[INFO] Zpracováno řádků: {n_rows}, ponecháno po filtru: {n_kept}, "
          f"nenamapováno na CpG index: {n_unmapped} ({pct_unmapped:.2f}%)", file=sys.stderr)
    if pct_unmapped > 5:
        print("[VAROVÁNÍ] Vysoký podíl nenamapovaných pozic - zkontrolujte, že "
              "referenční FASTA a CpG index odpovídají stejné genomové sestavě "
              "(hg38) a stejným chromozomovým názvům (chr1 vs 1).", file=sys.stderr)

    # Sestavení fragmentů (chrom_of_first_cpg neznáme přímo - dopočítáme
    # zpětně z indexu jen pro řazení výstupu; PAT formát chrom bere z
    # kontextu run - proto si chrom uchováváme zvlášť per read).
    # Pozn.: protože jsme groupovali podle chrom uvnitř chunku, čtení
    # rozdělená přes hranici chromozomu (nemělo by nastat u CpG) igno-
    # rujeme implicitně tím, že bereme chrom prvního záznamu readu.
    read_chrom = {}
    reader2 = pd.read_csv(args.calls, sep="\t", usecols=[cols["read_id"], cols["chrom"]], dtype=str)
    reader2 = reader2.rename(columns={cols["read_id"]: "read_id", cols["chrom"]: "chrom"})
    for read_id, chrom in zip(reader2["read_id"], reader2["chrom"]):
        if read_id not in read_chrom:
            read_chrom[read_id] = chrom

    fragment_counts = defaultdict(int)  # (chrom, first_index, pattern) -> count
    n_reads_total = 0
    n_reads_kept = 0
    for read_id, calls in read_calls.items():
        n_reads_total += 1
        if len(calls) < args.min_cpgs:
            continue
        calls_sorted = sorted(calls, key=lambda x: x[0])
        first_index = calls_sorted[0][0]
        pattern = "".join(args.meth_char if is_meth else args.unmeth_char for _, is_meth in calls_sorted)
        chrom = read_chrom.get(read_id, "NA")
        fragment_counts[(chrom, first_index, pattern)] += 1
        n_reads_kept += 1

    print(f"[INFO] Readů celkem: {n_reads_total}, readů zahrnutých do PAT "
          f"(>= {args.min_cpgs} CpG): {n_reads_kept}", file=sys.stderr)

    rows = sorted(fragment_counts.items(), key=lambda kv: (kv[0][0], kv[0][1]))

    out_path = args.out
    opener = gzip.open if out_path.endswith(".gz") else open
    with opener(out_path, "wt") as fh:
        for (chrom, first_index, pattern), count in rows:
            fh.write(f"{chrom}\t{first_index}\t{pattern}\t{count}\n")

    print(f"[OK] Zapsáno {len(rows)} unikátních fragmentů do {out_path}", file=sys.stderr)
    print("[DALŠÍ KROK] Než použijete v uxm deconv, zkontrolujte výstup: "
          "'wgbstools index' + 'wgbstools view' na malém regionu, případně "
          "srovnejte agregovanou % metylaci z tohoto pat.gz s bedMethyl "
          "výstupem ze skriptu 02 pro pár náhodných CpG - měly by se shodovat.",
          file=sys.stderr)


if __name__ == "__main__":
    main()
