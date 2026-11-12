# Pipeline: modkit → UXM konverze, referenční atlas, PCNSL markery na hg38

Sada skriptů pro tři úlohy z vaší poznámky (preliminary data: 4x DLBCL +
1x zdravé cfDNA, Nanopore, Dorado + Modkit):

1. konverze `modkit extract calls` výstupu do PAT formátu (`wgbs_tools` /
   `UXM_deconv`) pro fragment-level dekonvoluci ctDNA vs. hematopoetická cfDNA,
2. výběr/stažení referenčního metylačního atlasu a spuštění `uxm deconv`,
3. mapování publikovaného panelu PCNSL metylačních markerů na hg38 a
   extrakce jejich metylace z vašich vzorků.

## Než začnete – čtěte prosím

Tohle je **první funkční draft**, ne otestovaný produkční pipeline. `modkit`,
`wgbs_tools`/`UXM_deconv` i manifesty array sond jsou aktivně vyvíjené
nástroje třetích stran a jejich přesné rozhraní (názvy sloupců, CLI flagy,
URL ke stažení referencí) se mění mezi verzemi. Na několika místech jsem
narazil při rešerši na omezení nástrojů, které jsem používal k ověřování
(timeouty, blokace), takže konkrétní detaily (viz níže) **je potřeba ověřit
u vás na infrastruktuře před spuštěním na celém datasetu** – každý skript má
v hlavičce komentář, co přesně zkontrolovat a příkazem `--help`.

Doporučený postup: nejdřív spustit na jednom malém regionu / jednom vzorku,
zkontrolovat mezivýsledky (viz "Kontrolní body" níže), teprve pak pustit na
celý dataset.

## Přehled kroků

| # | Skript | Co dělá | Vstup | Výstup |
|---|--------|---------|-------|--------|
| - | `scripts/00_install_uxm.sh` | naklonuje UXM_deconv + ukáže svazkový atlas (viz caveat a) | - | `Navrkalova_cfDNA/tools/UXM_deconv/` |
| 0 | `scripts/01_prepare_reference.sh` | stáhne hg38 FASTA, spustí `wgbstools init_genome` (staví CpG index) | - | `genome.fa`, `CpG.bed.gz` |
| 1 | `scripts/02_modkit_pileup.sh` | region-level bedMethyl (agregovaná % metylace na CpG) | modBAM | `*.sorted.bedmethyl.gz` |
| 2 | `scripts/03_modkit_extract_calls.sh` | per-read CpG volání (nutné pro fragment-level PAT) | modBAM | `*.calls.tsv.gz` |
| 3 | `scripts/04_calls_to_pat.py` | konverze na PAT formát (`chrom, cpg_index, pattern, count`) | calls.tsv.gz + CpG index | `*.pat.gz` |
| 4 | `scripts/05_run_uxm.sh` | dekonvoluce ctDNA/hematopoetické frakce přes UXM (svazkový atlas) | pat.gz + atlas | `deconv_result.deconv.csv` |
| 5 | `scripts/06_build_pcnsl_marker_bed.py` | 8 CpG markerů (Clin Epigenetics 2021) → hg38 BED | probe ID + manifest | `pcnsl_markers.hg38.bed` |
| 6 | `scripts/07_marker_methylation_table.py` | % metylace na markerech napříč vzorky | bedmethyl + BED | `pcnsl_marker_matrix.tsv` |

`run_pipeline.sh` ukazuje, jak kroky navazují na sebe pro 5 vzorků.

## Prerekvizity (nainstalovat/ověřit dostupnost v PATH)

- `modkit` (ONT) – https://github.com/nanoporetech/modkit
- `wgbstools` – https://github.com/nloyfer/wgbs_tools
- `uxm` (UXM_deconv) – https://github.com/nloyfer/UXM_deconv
- `samtools`, `bedtools`, `bgzip`/`tabix` (htslib)
- Python 3.9+ s `pandas`, `numpy`

`modkit`, `wgbstools`, samtools/bedtools/htslib a Python stack se založí
jedním příkazem z přiloženého `environment.yml`. `uxm` (UXM_deconv) na to
navazuje samostatným krokem, protože nejde nainstalovat přes pip/conda
(viz caveat a níže):

```bash
conda env create -f environment.yml   # nebo: mamba env create -f environment.yml
conda activate pcnsl-methylation

# nastroje jsou ve sdilene slozce Navrkalova_cfDNA/tools (mimo git repozitar methylation_profiles)
./scripts/00_install_uxm.sh ../../../tools     # naklonuje UXM_deconv (pokud tam jeste neni) + ukáže cestu k atlasu
export PATH="$(cd ../../../tools/UXM_deconv && pwd):$PATH"   # a přidejte tenhle řádek i do ~/.bashrc

modkit --version && wgbstools --help && uxm --help   # ověření, že vše je v PATH
```

## Kontrolní body (validace mezivýsledků)

- **Po kroku 0**: zkontrolujte, že `wgbstools_ref/` obsahuje soubor s CpG
  indexem (očekávaně `CpG.bed.gz`) – pokud ne, uprav `load_cpg_index()` ve
  skriptu 04 podle skutečného výstupu vaší verze `wgbstools`.
- **Po kroku 3 (PAT)**: na malém regionu porovnejte agregovanou % metylaci
  spočtenou z `pat.gz` (počet 'C' vs 'T' na dané CpG pozici napříč
  fragmenty) s hodnotou ve stejné pozici v `bedmethyl.gz` z kroku 1 – měly
  by se přibližně shodovat. Výrazný nesoulad = něco je špatně namapované
  (často: rozdílná konvence chromozomových jmen `chr1` vs `1`, nebo posun
  souřadnic o 1 bp).
- **Po kroku 4 (UXM)**: zkontrolujte, že atlas skutečně obsahuje
  buněčný typ relevantní pro B-lymfocyty/lymfom – bez toho je odhad tumor
  frakce u DLBCL nesmyslný.
- **Krok 5–6 (PCNSL markery)**: viz caveat níže – jde o explorativní
  srovnání, ne diagnostický nástroj.

## Důležité věcné caveaty

**a) UXM instalace a atlas – opraveno a ověřeno.** V dřívější verzi tohoto
souboru jsem měl `uxm` jako `pip install git+...` v `environment.yml` a
psal, že se mi nepodařilo ověřit URL externě hostovaného atlasu. Po reálném
vyzkoušení (a chybě, kterou to vyvolalo při `conda env create`) se ukázalo:
`UXM_deconv` nemá `setup.py`/`pyproject.toml` (je to jen sbírka skriptů),
takže pip instalace nutně selže – proto je teď v samostatném kroku
`scripts/00_install_uxm.sh` (git clone + PATH, přesně jak radí oficiální
README repozitáře). A referenční atlas (Loyfer 2023) žádné externí
stahování nepotřebuje – je přímo součástí klonu repozitáře, ve složce
`supplemental/` (`Atlas.U25.l4.hg38.tsv` = 25 širokých typů,
`Atlas.U250.l4.hg38.tsv` = ~250 jemnějších subtypů, oba hg38). `05_run_uxm.sh`
na něj teď rovnou odkazuje. Jako alternativní srovnání pořád zvažte i
**`meth_atlas`** (Moss 2018, https://github.com/nloyfer/meth_atlas) –
jednodušší region-level NNLS dekonvoluce s vlastním menším atlasem.

**b) PCNSL markery jsou diagnostické, ne prognostické.** Panel 8 CpG
markerů (`markers/pcnsl_markers_450k.tsv`) pochází z práce, která odlišuje
**primární CNS lymfom od gliomu** (GBM/LGG) – ne z práce, která by
predikovala **CNS relaps u systémového DLBCL** (což je vaše skutečná
otázka). Je to nejkonkrétnější publikovaný CpG-level signál, který jsem při
rešerši našel, ale berte ho jako explorativní výchozí bod, který je potřeba
doplnit vlastní rešerší (CXCR4, PD-L1/PD-L2 lokus 9p24.1, MGMT a DAPK1
promotorová hypermetylace jsou v literatuře zmiňované kandidáty – jejich
přesné hg38 souřadnice ale skript zatím neobsahuje, je potřeba je doplnit
stejným postupem).

**c) n = 4 nádorové vzorky.** Na takhle malém vzorku nejde dělat
statisticky podloženou klasifikaci rizika CNS relapsu – smysluplný výstup
této fáze je deskriptivní srovnání s publikovanými signaturami
(hypothesis-generating), ne validovaný prediktor.

**d) Dvě CpG z markerového panelu mají v `pcnsl_markers_450k.tsv` nejasné
jméno genu** (`cg05491001`, `cg25567674`) – při extrakci textu ze zdrojové
publikace se nepodařilo spolehlivě získat název genu; probe ID samotné je
ale správné, doplňte gen ručně z tabulky v plném textu publikace.

## Zdroje použité při návrhu pipeline

- modkit bedMethyl formát: https://software-docs.nanoporetech.com/modkit/intro_pileup.html
- modkit extract calls: https://nanoporetech.github.io/modkit/intro_extract.html
- modkit GitHub: https://github.com/nanoporetech/modkit
- wgbs_tools: https://github.com/nloyfer/wgbs_tools
- UXM_deconv: https://github.com/nloyfer/UXM_deconv
- meth_atlas (Moss 2018): https://github.com/nloyfer/meth_atlas
- PCNSL markerový panel: Clinical Epigenetics (2021) 13:172,
  "Methylated markers accurately distinguish primary central nervous
  system lymphomas (PCNSL) from other CNS tumors",
  https://doi.org/10.1186/s13148-021-01091-9
- hg38 remapovaný Infinium manifest (Zhou lab): https://zwdzwd.github.io/InfiniumAnnotation
- GDC o zlepšené anotaci array sond: https://gdc.cancer.gov/content/improved-dna-methylation-array-probe-annotation
