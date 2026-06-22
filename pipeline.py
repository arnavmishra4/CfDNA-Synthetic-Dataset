"""
pipeline.py — Top-level orchestrator for the NeuroSight M5 synthetic cfDNA pipeline.

Ties together data loading, DMG augmentation, DMR finding, fragment generation,
ctDNA mixing, and JSON serialization into a single end-to-end run.
"""

import json
import logging
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
from tqdm import tqdm

from constants import (
    FRAGS_PER_PATIENT_MIN, FRAGS_PER_PATIENT_MAX,
    LABEL_MAP, MIN_CPGS,
)
from data_loader import load_illumina_manifest
from fragments import (
    FragmentGenerator, SequenceExtractor,
    build_methylated_string, find_cpg_positions, mix_ctdna_fragments,
)
from models import DMR, Fragment
from serializer import serialize_patient

log = logging.getLogger(__name__)


class SyntheticCfDNAPipeline:

    def __init__(self, hg38_fasta: str, manifest_path: str, output_dir: str, seed: int = 42):
        self.seq_extractor = SequenceExtractor(hg38_fasta)
        self.manifest      = load_illumina_manifest(manifest_path)
        self.output_dir    = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.rng = np.random.default_rng(seed)
        log.info(f"Output directory: {self.output_dir}")

    def run(
        self,
        matrices: Dict[int, "pd.DataFrame"],
        dmrs_by_label: Dict[int, List[DMR]],
        healthy_matrix: "pd.DataFrame",
        n_frags_target: int = 750,
    ) -> None:
        log.info("Pre-generating healthy fragment pool ...")
        healthy_frag_pool = self._generate_healthy_pool(healthy_matrix, dmrs_by_label)

        self._process_label_0(healthy_matrix, healthy_frag_pool, n_frags_target)

        for label, matrix in matrices.items():
            if label == 0:
                continue
            dmrs = dmrs_by_label.get(label, [])
            if not dmrs:
                log.warning(f"No DMRs for label {label}, skipping.")
                continue
            log.info(
                f"\n{'='*60}\nLabel={label} ({LABEL_MAP[label]}) — "
                f"{matrix.shape[1]} patients, {len(dmrs)} DMRs\n{'='*60}"
            )
            self._process_label(label, matrix, dmrs, healthy_frag_pool, n_frags_target)

        log.info("Pipeline complete.")

    def _generate_healthy_pool(
        self,
        healthy_matrix: "pd.DataFrame",
        dmrs_by_label: Dict[int, List[DMR]],
    ) -> List[Fragment]:
        """
        FIX-3: Per-probe beta values looked up by cg_id (not column-wide mean).
        """
        all_dmrs = [d for dl in dmrs_by_label.values() for d in dl]
        unique_dmrs: Dict[Tuple, DMR] = {}
        for d in all_dmrs:
            key = (d.chromosome, d.start, d.end)
            if key not in unique_dmrs:
                unique_dmrs[key] = d

        n_sample    = min(50, healthy_matrix.shape[1])
        sample_cols = self.rng.choice(healthy_matrix.columns, size=n_sample, replace=False)

        pool: List[Fragment] = []
        gen = FragmentGenerator(is_tumor=False, rng=self.rng)

        for (chrom, start, end), dmr in tqdm(unique_dmrs.items(), desc="Healthy pool"):
            seq = self.seq_extractor.get_sequence(chrom, start, end)
            if not seq:
                continue
            cpg_pos = find_cpg_positions(seq)
            if len(cpg_pos) < MIN_CPGS:
                continue

            for col in sample_cols:
                patient_betas = healthy_matrix[col]
                # FIX-3: look up actual per-probe beta by cg_id
                betas: Dict[int, float] = {}
                for k, cg_id in enumerate(dmr.cpg_ids):
                    if cg_id in patient_betas.index and k < len(cpg_pos):
                        val = patient_betas[cg_id]
                        if not np.isnan(val):
                            betas[cpg_pos[k]] = float(val)

                tokens = build_methylated_string(seq, cpg_pos, betas, self.rng)
                frags  = gen.fragment_token_stream(tokens, chrom, start)
                pool.extend(FragmentGenerator.filter_by_cpg(frags))

        log.info(f"  Healthy pool: {len(pool):,} fragments")
        return pool

    def _process_label_0(
        self,
        healthy_matrix: "pd.DataFrame",
        healthy_pool: List[Fragment],
        n_frags_target: int,
    ) -> None:
        log.info(f"Generating label=0 (Healthy) — {healthy_matrix.shape[1]} patients")
        label_dir = self.output_dir / "label_0_Healthy"
        label_dir.mkdir(exist_ok=True)

        for patient_id in tqdm(healthy_matrix.columns, desc="Healthy patients"):
            n_frags = int(self.rng.integers(FRAGS_PER_PATIENT_MIN, FRAGS_PER_PATIENT_MAX))
            idxs    = self.rng.choice(len(healthy_pool), size=min(n_frags, len(healthy_pool)), replace=True)
            frags   = [healthy_pool[i] for i in idxs]
            patient_json = serialize_patient(
                label=0, patient_id=str(patient_id), fragments=frags,
                raw_ctdna_fraction=0.0, enriched_ctdna_fraction=0.0,
                n_fragments_target=n_frags_target, rng=self.rng,
            )
            with open(label_dir / f"{patient_id}.json", "w") as f:
                json.dump(patient_json, f, separators=(",", ":"))

    def _process_label(
        self,
        label: int,
        matrix: "pd.DataFrame",
        dmrs: List[DMR],
        healthy_pool: List[Fragment],
        n_frags_target: int,
    ) -> None:
        label_dir = self.output_dir / f"label_{label}_{LABEL_MAP[label]}"
        label_dir.mkdir(exist_ok=True)
        gen = FragmentGenerator(is_tumor=True, rng=self.rng)

        for patient_id in tqdm(matrix.columns, desc=f"Label {label}"):
            patient_betas = matrix[patient_id]
            tumor_frags: List[Fragment] = []

            for dmr in dmrs:
                seq = self.seq_extractor.get_sequence(dmr.chromosome, dmr.start, dmr.end)
                if not seq:
                    continue
                cpg_pos = find_cpg_positions(seq)
                if len(cpg_pos) < MIN_CPGS:
                    continue

                betas: Dict[int, float] = {}
                for k, cg_id in enumerate(dmr.cpg_ids):
                    if cg_id in patient_betas.index and k < len(cpg_pos):
                        val = patient_betas[cg_id]
                        if not np.isnan(val):
                            betas[cpg_pos[k]] = float(val)

                tokens = build_methylated_string(seq, cpg_pos, betas, self.rng)
                frags  = gen.fragment_token_stream(tokens, dmr.chromosome, dmr.start)
                tumor_frags.extend(FragmentGenerator.filter_by_cpg(frags))

            # FIX-2: unpack all three return values
            mixed, raw_ctdna, enriched_ctdna = mix_ctdna_fragments(
                tumor_frags, healthy_pool, self.rng, enriched=True
            )
            patient_json = serialize_patient(
                label=label, patient_id=str(patient_id), fragments=mixed,
                raw_ctdna_fraction=raw_ctdna, enriched_ctdna_fraction=enriched_ctdna,
                n_fragments_target=n_frags_target, rng=self.rng,
            )
            with open(label_dir / f"{patient_id}.json", "w") as f:
                json.dump(patient_json, f, separators=(",", ":"))
