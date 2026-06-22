"""
augmentation.py — DMG augmentation via Beta distribution fitting.

Augments the small DMG H3K27M cohort (GSE161944, ~20 real samples) to
DMG_AUGMENT_TARGET by fitting independent Beta distributions per probe.
"""

import logging
from typing import Dict, Tuple

import numpy as np
import pandas as pd
from scipy.stats import beta as beta_dist
from tqdm import tqdm

from constants import DMG_AUGMENT_TARGET

log = logging.getLogger(__name__)


class DMGAugmentor:
    """
    Augments the small DMG H3K27M cohort (GSE161944, ~20 real samples) to
    DMG_AUGMENT_TARGET by fitting independent Beta distributions per probe.
    """

    def fit(self, matrix: pd.DataFrame) -> None:
        log.info(
            f"Fitting Beta distributions for DMG augmentation "
            f"({matrix.shape[1]} real samples, GSE161944 only) ..."
        )
        self._params: Dict[str, Tuple[float, float]] = {}
        for probe_id, row in tqdm(matrix.iterrows(), total=len(matrix), desc="Beta fit"):
            vals = row.dropna().values.clip(1e-6, 1 - 1e-6)
            if len(vals) < 3:
                mu = float(np.mean(vals)) if len(vals) else 0.5
                self._params[probe_id] = (mu * 10, (1 - mu) * 10)
                continue
            try:
                a, b, _, _ = beta_dist.fit(vals, floc=0, fscale=1)
                self._params[probe_id] = (a, b)
            except Exception:
                mu = float(np.mean(vals))
                self._params[probe_id] = (mu * 10, (1 - mu) * 10)

    def sample(self, n: int = DMG_AUGMENT_TARGET) -> pd.DataFrame:
        log.info(f"Sampling {n} synthetic DMG patients ...")
        data = {
            probe_id: beta_dist.rvs(a, b, size=n).clip(0, 1)
            for probe_id, (a, b) in self._params.items()
        }
        df = pd.DataFrame(data).T
        df.columns = [f"DMG_synthetic_{i:04d}" for i in range(n)]
        log.info(f"  Augmented DMG matrix: {df.shape}")
        return df
