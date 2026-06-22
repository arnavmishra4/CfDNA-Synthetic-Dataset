"""
constants.py — Pipeline-wide constants for NeuroSight M5 synthetic cfDNA pipeline.

Rules implemented per: Synthetic_GBM_cfDNA_Rulebook.pdf (Arnav Mishra, 2025)
"""

# ---------------------------------------------------------------------------
# Methylation thresholds
# ---------------------------------------------------------------------------

BETA_HIGH = 0.70
BETA_LOW  = 0.30

# ---------------------------------------------------------------------------
# DMR identification
# ---------------------------------------------------------------------------

N_DMRS      = 300
DELTA_BETA  = 0.30
MIN_CPGS    = 2

# ---------------------------------------------------------------------------
# Fragment length model (Rule 5)
# ---------------------------------------------------------------------------

FRAG_HEALTHY_MONO   = 166
FRAG_TUMOR_MONO     = 139
FRAG_DINU           = 320
FRAG_SD_HEALTHY     = 15
FRAG_SD_TUMOR       = 12
FRAG_MONO_RATIO     = 0.75
FRAG_MIN            = 100
FRAG_MAX            = 400
FRAG_TUMOR_MONO_MIN = 134
FRAG_TUMOR_MONO_MAX = 160

# ---------------------------------------------------------------------------
# CpG per fragment (Rule 8)
# ---------------------------------------------------------------------------

MIN_CPG_PER_FRAG = 2
MAX_CPG_PER_FRAG = 8

# ---------------------------------------------------------------------------
# ctDNA fraction (Rule 7)
# ---------------------------------------------------------------------------

CTDNA_MEDIAN       = 3.1e-5
CTDNA_RANGE_MAX    = 1e-4
ENRICHED_CTDNA_MIN = 0.05
ENRICHED_CTDNA_MAX = 0.30

# ---------------------------------------------------------------------------
# Dataset sizing
# ---------------------------------------------------------------------------

FRAGS_PER_PATIENT_MIN = 500
FRAGS_PER_PATIENT_MAX = 1000

WINDOW_BP = 300

DMG_AUGMENT_TARGET = 200

# ---------------------------------------------------------------------------
# 4-class label map — no 5th class
# ---------------------------------------------------------------------------

LABEL_MAP = {
    0: "Healthy",
    1: "GBM",
    2: "LGG",
    3: "DMG_H3K27M",
}

# ---------------------------------------------------------------------------
# Healthy cell-type fractions
# ---------------------------------------------------------------------------

HEALTHY_CELL_FRACTIONS = {
    "WBC":         0.55,
    "RBC_prog":    0.30,
    "Endothelium": 0.10,
    "Hepatocyte":  0.01,
    "Other":       0.04,
}
