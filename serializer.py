"""
serializer.py — JSON serialization for synthetic cfDNA patient records.
"""

import logging
from typing import Dict, List, Optional

import numpy as np

from constants import LABEL_MAP, WINDOW_BP
from models import Fragment

log = logging.getLogger(__name__)


def serialize_patient(
    label: int,
    patient_id: str,
    fragments: List[Fragment],
    raw_ctdna_fraction: float,
    enriched_ctdna_fraction: float,
    n_fragments_target: int = 750,
    rng: Optional[np.random.Generator] = None,
) -> Dict:
    """
    Serialize a patient's fragments to the Pleiades-compatible JSON format.

    Stores both raw_ctdna_fraction (biological ground truth) and
    enriched_ctdna_fraction (actual mixing ratio post-cfMeDIP, 5-30%).
    """
    if rng is None:
        rng = np.random.default_rng()

    if len(fragments) > n_fragments_target:
        idxs      = rng.choice(len(fragments), size=n_fragments_target, replace=False)
        fragments = [fragments[i] for i in idxs]

    regions: Dict[str, Dict] = {}
    region_counter = 0
    chrom_groups: Dict[str, List[Fragment]] = {}
    for frag in fragments:
        chrom_groups.setdefault(frag.chromosome, []).append(frag)

    for chrom, chrom_frags in sorted(chrom_groups.items()):
        chrom_frags  = sorted(chrom_frags, key=lambda f: f.genomic_start)
        window_start = chrom_frags[0].genomic_start if chrom_frags else 0
        window_frags: List[Fragment] = []

        for frag in chrom_frags:
            if frag.genomic_start - window_start > WINDOW_BP and window_frags:
                regions[f"region_{region_counter}"] = _make_region(chrom, window_start, window_frags)
                region_counter += 1
                window_start = frag.genomic_start
                window_frags = []
            window_frags.append(frag)

        if window_frags:
            regions[f"region_{region_counter}"] = _make_region(chrom, window_start, window_frags)
            region_counter += 1

    return {
        "label":                   label,
        "subtype":                 LABEL_MAP[label],
        "patient_id":              patient_id,
        "raw_ctdna_fraction":      round(raw_ctdna_fraction, 9),
        "enriched_ctdna_fraction": round(enriched_ctdna_fraction, 4),
        "n_fragments":             len(fragments),
        "regions":                 regions,
    }


def _make_region(chrom: str, window_start: int, frags: List[Fragment]) -> Dict:
    return {
        "chromosome": int(chrom.lstrip("chr")) if chrom.lstrip("chr").isdigit() else chrom,
        "genomic_start": window_start,
        "fragments": [_frag_to_dict(f) for f in frags],
    }


def _frag_to_dict(f: Fragment) -> Dict:
    return {
        "tokens":          f.tokens,
        "fragment_length": f.fragment_length,
        "cpg_count":       f.cpg_count,
        "strand":          f.strand,
    }
