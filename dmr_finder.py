"""
dmr_finder.py — Memory-efficient DMR identification for NeuroSight M5.

Processes one class at a time — never concatenates all matrices at once.
Subsamples to MAX_SAMPLES per class before t-test to save RAM.
"""

import gc
import logging
from typing import Dict, List

import numpy as np
import pandas as pd
from scipy import stats

from constants import DELTA_BETA, LABEL_MAP, MIN_CPGS, N_DMRS, WINDOW_BP
from models import DMR

log = logging.getLogger(__name__)


class DMRFinder:
    """
    Memory-efficient DMR finder.
    Processes one class at a time — never concatenates all matrices at once.
    Subsamples to MAX_SAMPLES per class before t-test to save RAM.
    """

    MAX_SAMPLES_PER_CLASS = 80

    def __init__(self, cpg_manifest: pd.DataFrame):
        self.manifest = cpg_manifest

    def find_dmrs(
        self,
        matrices: Dict[int, pd.DataFrame],
        n_dmrs: int = N_DMRS,
        delta_beta_thresh: float = DELTA_BETA,
        min_cpgs: int = MIN_CPGS,
    ) -> Dict[int, List[DMR]]:
        log.info("Finding DMRs (memory-efficient one-vs-rest t-test) ...")

        # ── Subsample each class to MAX_SAMPLES to save RAM ──
        rng = np.random.default_rng(42)
        subsampled: Dict[int, np.ndarray] = {}
        cpg_index = None

        for lbl, mat in matrices.items():
            n = mat.shape[1]
            if n > self.MAX_SAMPLES_PER_CLASS:
                idx = rng.choice(n, self.MAX_SAMPLES_PER_CLASS, replace=False)
                sub = mat.iloc[:, idx].values.astype(np.float32)
            else:
                sub = mat.values.astype(np.float32)
            subsampled[lbl] = sub
            if cpg_index is None:
                cpg_index = mat.index.tolist()
            log.info(f"  Label {lbl} ({LABEL_MAP[lbl]}): using {sub.shape[1]} samples")

        del matrices  # free original matrices from RAM
        gc.collect()

        cpg_array = np.array(cpg_index)
        dmrs_by_label: Dict[int, List[DMR]] = {}

        for target_label in sorted(subsampled.keys()):
            log.info(f"  Computing DMRs for label {target_label} ({LABEL_MAP[target_label]}) ...")

            pos_mat = subsampled[target_label]  # (n_cpgs, n_pos)

            # Stack negatives one at a time to avoid one giant concat
            neg_parts = [subsampled[lbl] for lbl in subsampled if lbl != target_label]
            neg_mat   = np.concatenate(neg_parts, axis=1)  # (n_cpgs, n_neg)

            # t-test row-wise (per CpG)
            t_stats, p_vals = stats.ttest_ind(
                pos_mat, neg_mat, axis=1, equal_var=False, nan_policy="omit"
            )
            delta = np.nanmean(pos_mat, axis=1) - np.nanmean(neg_mat, axis=1)

            del neg_mat; gc.collect()

            # BH correction + filter
            valid    = ~(np.isnan(t_stats) | np.isnan(p_vals))
            t_stats  = t_stats[valid]
            p_vals   = p_vals[valid]
            delta    = delta[valid]
            cpgs_v   = cpg_array[valid]

            sort_idx = np.argsort(p_vals)
            t_stats  = t_stats[sort_idx]
            p_vals   = p_vals[sort_idx]
            delta    = delta[sort_idx]
            cpgs_v   = cpgs_v[sort_idx]

            n        = len(p_vals)
            ranks    = np.arange(1, n + 1)
            bh_thresh = ranks / n * 0.01
            keep     = (p_vals <= bh_thresh) & (np.abs(delta) >= delta_beta_thresh)

            t_stats_f = t_stats[keep]
            p_vals_f  = p_vals[keep]
            delta_f   = delta[keep]
            cpgs_f    = cpgs_v[keep]

            log.info(f"    Significant CpGs after BH: {keep.sum():,}")

            # Join with manifest for coordinates
            result_df = pd.DataFrame({
                "cg_id":      cpgs_f,
                "t_stat":     t_stats_f,
                "p_val":      p_vals_f,
                "delta_beta": delta_f,
            }).set_index("cg_id")

            result_df = result_df.join(self.manifest, how="inner")

            dmr_list  = self._group_into_windows(result_df, target_label, min_cpgs)
            dmr_list  = sorted(dmr_list, key=lambda d: abs(d.t_stat), reverse=True)[:n_dmrs]
            dmrs_by_label[target_label] = dmr_list
            log.info(f"    Label {target_label}: {len(dmr_list)} DMRs found")

            del result_df; gc.collect()

        return dmrs_by_label

    def _group_into_windows(
        self,
        result_df: pd.DataFrame,
        label: int,
        min_cpgs: int,
    ) -> List[DMR]:
        dmrs = []
        for chrom, grp in result_df.groupby("chromosome"):
            grp  = grp.sort_values("start")
            used = np.zeros(len(grp), dtype=bool)
            for i, (idx, row) in enumerate(grp.iterrows()):
                if used[i]:
                    continue
                window_start = int(row["start"])
                window_end   = window_start + WINDOW_BP
                in_window    = (
                    (grp["start"] >= window_start) & (grp["start"] < window_end)
                ).values
                if in_window.sum() < min_cpgs:
                    continue
                window_cpgs = grp[in_window]
                dmrs.append(DMR(
                    chromosome=str(chrom),
                    start=window_start,
                    end=window_end,
                    cpg_ids=window_cpgs.index.tolist(),
                    delta_beta=float(window_cpgs["delta_beta"].mean()),
                    t_stat=float(window_cpgs["t_stat"].mean()),
                    label=label,
                ))
                used[in_window] = True
        return dmrs
