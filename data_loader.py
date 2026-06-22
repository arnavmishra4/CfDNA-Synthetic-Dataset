"""
data_loader.py — Data loading utilities for NeuroSight M5 synthetic cfDNA pipeline.

Handles TCGA per-patient files, GEO matrix files, and the Illumina 450K manifest.
"""

import os
import glob
import gzip
import logging
import urllib.request
import io as _io
from pathlib import Path
from typing import List

import numpy as np
import pandas as pd
from tqdm import tqdm

log = logging.getLogger(__name__)


class DataLoader:

    def load_tcga_per_patient(
        self,
        root_dir: str,
        label: int,
        na_threshold: float = 0.70,
    ) -> pd.DataFrame:
        """
        Load TCGA per-patient files (one .txt per UUID subfolder).
        Each file: no header, tab-sep, col0=cg_id, col1=beta_value.
        NAs ~67% per file are expected for SeSAMe Level 3 processed data.
        """
        log.info(f"Loading TCGA data from {root_dir} (label={label})")
        txt_files = glob.glob(
            os.path.join(root_dir, "**", "*.methylation_array.sesame.level3betas.txt"),
            recursive=True,
        )
        if not txt_files:
            raise FileNotFoundError(
                f"No TCGA methylation files found under {root_dir}. "
                "Expected: <UUID>/<UUID>.methylation_array.sesame.level3betas.txt"
            )

        frames = []
        for fpath in tqdm(txt_files, desc=f"TCGA label={label}"):
            patient_id = Path(fpath).stem
            df = pd.read_csv(
                fpath, sep="\t", header=None,
                names=["cg_id", patient_id], index_col="cg_id",
            )
            frames.append(df)

        matrix = pd.concat(frames, axis=1)
        log.info(f"  Raw shape: {matrix.shape}")

        na_rate = matrix.isna().mean(axis=1)
        matrix  = matrix.loc[na_rate <= na_threshold]
        log.info(
            f"  After NA filter ({na_threshold:.0%}): {matrix.shape} "
            f"({(na_rate > na_threshold).sum()} probes dropped)"
        )
        matrix = matrix.apply(lambda col: col.fillna(col.median()), axis=1)
        log.info(f"  Beta mean={matrix.values.mean():.3f}  std={matrix.values.std():.3f}")
        return matrix

    def load_geo_matrix(
        self,
        filepath: str,
        id_col: str = "ID_REF",
    ) -> pd.DataFrame:
        """
        Load a GEO matrix file (probes x samples).
        Handles GSE161944 (850K DMG) and GSE40279 (450K healthy).
        """
        log.info(f"Loading GEO matrix: {filepath}")
        df = pd.read_csv(filepath, sep="\t", index_col=id_col, low_memory=False)
        df = df.apply(pd.to_numeric, errors="coerce")
        na_rate = df.isna().mean(axis=1)
        df = df.loc[na_rate < 0.80]
        df = df.apply(lambda col: col.fillna(col.median()), axis=1)
        log.info(f"  Shape: {df.shape}  Beta mean={df.values.mean():.3f}")
        return df


def compute_common_cpg_space(*matrices: pd.DataFrame) -> List[str]:
    """Return the sorted intersection of probe IDs across all input matrices."""
    common = set(matrices[0].index)
    for m in matrices[1:]:
        common &= set(m.index)
    common = sorted(common)
    log.info(f"Common CpG space: {len(common):,} probes")
    return common


def load_illumina_manifest(manifest_path: str) -> pd.DataFrame:
    """
    Load probe coordinates for the Illumina HumanMethylation450 array.

    Tries sources in order:
      1. Local CSV file at manifest_path (if it exists)
      2. GPL13534 platform table fetched directly from NCBI GEO (no download needed)

    Returns:
        DataFrame indexed by cg_id with columns ['chromosome', 'start', 'end'].
    """
    # Option 1: local CSV
    if os.path.exists(manifest_path):
        log.info(f"Loading 450K manifest from local file: {manifest_path}")
        df = pd.read_csv(manifest_path, skiprows=7, low_memory=False)
        df = df[["IlmnID", "CHR", "MAPINFO"]].dropna()
        df.columns = ["cg_id", "chromosome", "start"]
        df["start"] = df["start"].astype(int)
        df["end"]   = df["start"] + 1
        return df.set_index("cg_id")

    # Option 2: fetch GPL13534 from NCBI GEO
    log.info(
        "Manifest CSV not found — fetching GPL13534 probe coordinates "
        "from NCBI GEO (one-time ~10 MB download) ..."
    )
    url = "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL13nnn/GPL13534/soft/GPL13534_family.soft.gz"
    try:
        with urllib.request.urlopen(url, timeout=120) as resp:
            raw = resp.read()
    except Exception as e:
        raise RuntimeError(
            f"Could not fetch GPL13534 from NCBI GEO: {e}\n"
            "Please download HumanMethylation450_15017482_v1-2.csv from Illumina "
            "and add it as a Kaggle dataset at /kaggle/input/illumina-manifest/."
        )

    log.info("Parsing GPL13534 platform table ...")
    records  = []
    in_table = False
    idx_id = idx_chr = idx_pos = None

    with gzip.open(_io.BytesIO(raw), "rt", encoding="latin-1") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("!platform_table_begin"):
                in_table   = True
                header     = next(fh).rstrip("\n").split("\t")
                cols_lower = [c.lower() for c in header]
                idx_id  = cols_lower.index("id")
                idx_chr = next(i for i, c in enumerate(cols_lower) if "chr" in c)
                idx_pos = next(i for i, c in enumerate(cols_lower) if "mapinfo" in c)
                continue
            if line.startswith("!platform_table_end"):
                break
            if not in_table:
                continue
            parts = line.split("\t")
            if len(parts) <= max(idx_id, idx_chr, idx_pos):
                continue
            cg_id = parts[idx_id].strip()
            chrom = parts[idx_chr].strip()
            pos   = parts[idx_pos].strip()
            if not cg_id.startswith("cg") or not chrom or not pos:
                continue
            try:
                records.append((cg_id, chrom, int(float(pos))))
            except ValueError:
                continue

    if not records:
        raise RuntimeError(
            "GPL13534 parsed 0 records — platform table format may have changed. "
            "Please supply the Illumina manifest CSV manually."
        )

    df = pd.DataFrame(records, columns=["cg_id", "chromosome", "start"])
    df["end"] = df["start"] + 1
    df = df.set_index("cg_id")
    log.info(f"GPL13534 manifest loaded: {len(df):,} probes")
    return df
