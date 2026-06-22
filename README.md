

```markdown
# Synthetic GBM cfDNA Generation Pipeline
## Overview
This repository contains a biologically-grounded, citation-backed pipeline for generating synthetic cell-free DNA (cfDNA) datasets targeted at glioma subtype classification. 
Due to the lack of public glioblastoma (GBM) plasma cfDNA datasets, this methodology reconstructs realistic cfDNA fragments by mapping microarray beta values to the reference genome and simulating the physical dynamics of post-cfMeDIP-seq enrichment. 
## Pipeline Architecture
The generation pipeline follows strict biological constraints derived from fragmentomics literature:
* **Data Sourcing & Baseline:** Utilizes TCGA-GBM, TCGA-LGG, and GEO datasets for tumor methylation patterns, mapped against the hg38 reference genome and a healthy plasma methylation atlas (GSE122126).
* **DMR-First Windowing:** Extracts sequences exclusively from the Top 300 Differentially Methylated Regions (DMRs) identified via `limma-trend`, restricting extraction to biologically relevant areas like CpG Islands, Shores, Shelves, and FANTOM5 Enhancers.
* **Beta Value Tokenization:** Converts 450K array beta values into discrete methylation tokens (`<m>`, `<um>`) and fuses them with A/T/G/C nucleotides to create Pleiades-compatible text sequences.
* **Nucleosome-Aware Fragmentation:** Fragments DNA strings based on nuclease cleavage preferences and nucleosome footprints. Simulates the physiological bimodal distribution: ~166 bp for healthy fragments and ~134-144 bp for tumor-derived fragments[cite: 1].
* **cfMeDIP-seq Simulation:** Enforces a physical enrichment filter, requiring 2 to 8 methylated CpGs per fragment, simulating the actual capture dynamics of immunoprecipitation assays[cite: 1].
* **Fractional Mixing:** Accurately mixes simulated tumor fractions with healthy background fractions (WBCs, RBC progenitors, endothelium)[cite: 1].
## Output Format
The pipeline outputs serialized JSON files directly compatible with the Pleiades foundational model architecture[cite: 1]. 
```json
{
  "label": 1,
  "subtype": "GBM",
  "ctdna_fraction": 0.000031,
  "regions": {
    "region_0": {
      "chromosome": "chr2",
      "genomic_start": 89768000,
      "fragments": [
        {
          "tokens": ["<cfdna>", "A", "T", "<m>", "C", "</cfdna>"],
          "fragment_length": 141,
          "cpg_count": 4,
          "strand": "+"
        }
      ]
    }
  }
}
```

---

## How to Run

### Prerequisites

Install dependencies:

```bash
pip install numpy pandas scipy tqdm pyfaidx
```

You also need the hg38 reference FASTA. On Kaggle, add it as a dataset. Locally, download from UCSC:

```bash
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz
```

### Required Data

| Dataset | Source | Label |
|---|---|---|
| TCGA-GBM methylation (SeSAMe Level 3) | [GDC Portal](https://portal.gdc.cancer.gov) | 1 |
| TCGA-LGG methylation (SeSAMe Level 3) | [GDC Portal](https://portal.gdc.cancer.gov) | 2 |
| GSE161944 (DMG H3K27M, 850K) | [NCBI GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE161944) | 3 |
| GSE40279 (healthy controls, 450K) | [NCBI GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE40279) | 0 |
| hg38 reference FASTA | [UCSC](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/) | — |
| Illumina 450K manifest (optional) | [Illumina](https://support.illumina.com) | — |

> The Illumina manifest is optional — the pipeline fetches probe coordinates from NCBI GEO (GPL13534) automatically if the local CSV is missing.

### File Structure

Place all `.py` files in the same directory:

```
m5_pipeline/
├── main.py
├── pipeline.py
├── fragments.py
├── dmr_finder.py
├── data_loader.py
├── augmentation.py
├── serializer.py
├── models.py
└── constants.py
```

### Configuration

Open `main.py` and set the paths at the top of the file:

```python
HG38_FASTA       = "/path/to/hg38.fa"
MANIFEST_PATH    = "/path/to/HumanMethylation450_15017482_v1-2.csv"  # optional
TCGA_GBM_DIR     = "/path/to/tcga-gbm-methylation/"
TCGA_LGG_DIR     = "/path/to/tcga-lgg-methylation/"
GEO_DMG_PATH     = "/path/to/GSE161944_series_matrix.txt"
GEO_HEALTHY_PATH = "/path/to/GSE40279_series_matrix.txt"
OUTPUT_DIR       = "/path/to/output/"
```

TCGA directories should contain UUID subfolders in GDC format:
```
tcga-gbm-methylation/
└── <UUID>/
    └── <UUID>.methylation_array.sesame.level3betas.txt
```

### Run

```bash
python main.py
```

The pipeline will log progress at each stage. A full run across all four classes typically takes 2–4 hours depending on TCGA cohort size and available RAM.

### Output

One JSON file per patient, organized by label:

```
output/
├── label_0_Healthy/
│   ├── GSM123456.json
│   └── ...
├── label_1_GBM/
│   ├── TCGA-06-0124.json
│   └── ...
├── label_2_LGG/
│   └── ...
└── label_3_DMG_H3K27M/
    └── ...
```

Each JSON contains up to 750 fragments per patient by default (`N_FRAGS_TARGET` in `main.py`).

### Kaggle

If running on Kaggle, the default paths in `main.py` already point to `/kaggle/input/` and `/kaggle/working/`. Add your datasets via **Add Data** and match the dataset slugs to the paths, or update them to match your dataset names.
```
