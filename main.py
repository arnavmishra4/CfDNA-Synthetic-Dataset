"""
main.py — Entrypoint for NeuroSight M5 synthetic cfDNA pipeline.

Usage (Kaggle or local):
    python main.py

Set paths in the CONFIG block below before running.
"""

import logging
import sys

import numpy as np

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)],
)
log = logging.getLogger(__name__)


# ===========================================================================
# CONFIG — set these paths before running
# ===========================================================================

# hg38 reference FASTA (needs .fai index alongside it, or pyfaidx will build one)
HG38_FASTA      = "/kaggle/input/hg38-reference/hg38.fa"

# Illumina 450K manifest CSV (optional — pipeline fetches from NCBI if missing)
MANIFEST_PATH   = "/kaggle/input/illumina-manifest/HumanMethylation450_15017482_v1-2.csv"

# TCGA-GBM per-patient methylation files (label=1)
# Each file: <UUID>/<UUID>.methylation_array.sesame.level3betas.txt
TCGA_GBM_DIR    = "/kaggle/input/tcga-gbm-methylation"

# TCGA-LGG per-patient methylation files (label=2)
TCGA_LGG_DIR    = "/kaggle/input/tcga-lgg-methylation"

# GEO matrix for DMG H3K27M (GSE161944, label=3)
GEO_DMG_PATH    = "/kaggle/input/gse161944/GSE161944_series_matrix.txt"

# GEO matrix for healthy controls (GSE40279, label=0)
GEO_HEALTHY_PATH = "/kaggle/input/gse40279/GSE40279_series_matrix.txt"

# Where to write the output JSON files
OUTPUT_DIR      = "/kaggle/working/m5_synthetic_cfdna"

# Pipeline settings
SEED            = 42
N_FRAGS_TARGET  = 750   # fragments per patient in final JSON
N_DMRS          = 300   # top DMRs to keep per class
DMG_AUGMENT_N   = 200   # synthetic DMG patients to generate via Beta augmentation

# ===========================================================================


def main():
    from data_loader import DataLoader, compute_common_cpg_space
    from augmentation import DMGAugmentor
    from dmr_finder import DMRFinder
    from pipeline import SyntheticCfDNAPipeline

    loader = DataLoader()

    # ── 1. Load methylation matrices ──────────────────────────────────────
    log.info("=== Loading methylation matrices ===")
    gbm_matrix     = loader.load_tcga_per_patient(TCGA_GBM_DIR,  label=1)
    lgg_matrix     = loader.load_tcga_per_patient(TCGA_LGG_DIR,  label=2)
    dmg_matrix_raw = loader.load_geo_matrix(GEO_DMG_PATH)
    healthy_matrix = loader.load_geo_matrix(GEO_HEALTHY_PATH)

    # ── 2. Augment small DMG cohort ───────────────────────────────────────
    log.info("=== Augmenting DMG cohort ===")
    augmentor = DMGAugmentor()
    augmentor.fit(dmg_matrix_raw)
    dmg_matrix = augmentor.sample(n=DMG_AUGMENT_N)

    # ── 3. Align to common CpG space ─────────────────────────────────────
    log.info("=== Computing common CpG space ===")
    common_cpgs = compute_common_cpg_space(
        gbm_matrix, lgg_matrix, dmg_matrix, healthy_matrix
    )
    gbm_matrix     = gbm_matrix.loc[common_cpgs]
    lgg_matrix     = lgg_matrix.loc[common_cpgs]
    dmg_matrix     = dmg_matrix.loc[common_cpgs]
    healthy_matrix = healthy_matrix.loc[common_cpgs]

    matrices = {
        1: gbm_matrix,
        2: lgg_matrix,
        3: dmg_matrix,
    }

    # ── 4. Find DMRs ──────────────────────────────────────────────────────
    log.info("=== Finding DMRs ===")
    from data_loader import load_illumina_manifest
    manifest    = load_illumina_manifest(MANIFEST_PATH)
    # Filter manifest to common CpG space for efficiency
    manifest    = manifest.loc[manifest.index.intersection(common_cpgs)]

    dmr_finder  = DMRFinder(cpg_manifest=manifest)
    dmrs_by_label = dmr_finder.find_dmrs(
        matrices={**matrices, 0: healthy_matrix},
        n_dmrs=N_DMRS,
    )
    log.info(f"DMRs found per label: { {k: len(v) for k, v in dmrs_by_label.items()} }")

    # ── 5. Run fragment generation pipeline ───────────────────────────────
    log.info("=== Running fragment pipeline ===")
    pipeline = SyntheticCfDNAPipeline(
        hg38_fasta=HG38_FASTA,
        manifest_path=MANIFEST_PATH,
        output_dir=OUTPUT_DIR,
        seed=SEED,
    )
    pipeline.run(
        matrices=matrices,
        dmrs_by_label=dmrs_by_label,
        healthy_matrix=healthy_matrix,
        n_frags_target=N_FRAGS_TARGET,
    )

    log.info(f"Done. Output written to: {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
