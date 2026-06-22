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
