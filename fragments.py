"""
fragments.py — Sequence extraction, token conversion, and fragment generation
for the NeuroSight M5 synthetic cfDNA pipeline.

Integrates three rulebook fixes from the original patch cell:
  FIX-1 (fragment_token_stream) — tail fragments that are truncated by the
         window boundary are merged into their predecessor instead of being
         emitted as distribution-violating remainder fragments.
  FIX-2 (filter_by_cpg)         — counts total CpG loci (<m> + <um>), not
         just methylated ones, matching Rule 8's definition of "CpG site."
  FIX-3 (mix_ctdna_fragments)   — healthy/tumor counts are derived FROM the
         tumor pool so enriched_frac is the realized ratio, not an upper
         bound drowned out by a large healthy pool.
"""

import logging
from typing import Dict, List, Tuple, Optional

import numpy as np
from pyfaidx import Fasta

from constants import (
    BETA_HIGH, BETA_LOW,
    FRAG_DINU, FRAG_MIN, FRAG_MAX,
    FRAG_MONO_RATIO,
    FRAG_SD_HEALTHY, FRAG_SD_TUMOR,
    FRAG_HEALTHY_MONO, FRAG_TUMOR_MONO,
    FRAG_TUMOR_MONO_MIN, FRAG_TUMOR_MONO_MAX,
    MIN_CPG_PER_FRAG, MAX_CPG_PER_FRAG,
    CTDNA_MEDIAN, CTDNA_RANGE_MAX,
    ENRICHED_CTDNA_MIN, ENRICHED_CTDNA_MAX,
)
from models import Fragment

log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Sequence Extraction
# ---------------------------------------------------------------------------

class SequenceExtractor:

    def __init__(self, hg38_fasta: str):
        log.info(f"Loading hg38 reference: {hg38_fasta}")
        self.fasta = Fasta(hg38_fasta, rebuild=True)

    def get_sequence(self, chrom: str, start: int, end: int) -> str:
        chrom_key = chrom if chrom in self.fasta else f"chr{chrom}"
        try:
            return str(self.fasta[chrom_key][start:end]).upper()
        except Exception:
            return ""


# ---------------------------------------------------------------------------
# Token / CpG helpers
# ---------------------------------------------------------------------------

def beta_to_token(beta: float, rng: np.random.Generator) -> str:
    if np.isnan(beta):
        return "<um>"
    if beta > BETA_HIGH:
        return "<m>"
    if beta < BETA_LOW:
        return "<um>"
    return "<m>" if rng.random() < beta else "<um>"


def build_methylated_string(
    sequence: str,
    cpg_positions: List[int],
    beta_values: Dict[int, float],
    rng: np.random.Generator,
) -> List[str]:
    tokens  = []
    cpg_set = set(cpg_positions)
    for i, nuc in enumerate(sequence):
        if i in cpg_set:
            tokens.append(beta_to_token(beta_values.get(i, np.nan), rng))
            tokens.append("G")
        elif nuc in "ATCG":
            tokens.append(nuc)
        else:
            tokens.append("N")
    return tokens


def find_cpg_positions(sequence: str) -> List[int]:
    seq, positions = sequence.upper(), []
    for i in range(len(seq) - 1):
        if seq[i] == "C" and seq[i + 1] == "G":
            positions.append(i)
    return positions


# ---------------------------------------------------------------------------
# Fragment Generation
# ---------------------------------------------------------------------------

class FragmentGenerator:

    def __init__(self, is_tumor: bool, rng: np.random.Generator):
        self.is_tumor   = is_tumor
        self.rng        = rng
        self._mono_mean = FRAG_TUMOR_MONO   if is_tumor else FRAG_HEALTHY_MONO
        self._mono_sd   = FRAG_SD_TUMOR     if is_tumor else FRAG_SD_HEALTHY

    def sample_fragment_length(self) -> int:
        if self.rng.random() < FRAG_MONO_RATIO:
            length = int(self.rng.normal(self._mono_mean, self._mono_sd))
            if self.is_tumor:
                length = int(np.clip(length, FRAG_TUMOR_MONO_MIN, FRAG_TUMOR_MONO_MAX))
            else:
                length = int(np.clip(length, FRAG_MIN, FRAG_MAX))
        else:
            length = int(np.clip(self.rng.normal(FRAG_DINU, 15), FRAG_MIN, FRAG_MAX))
        return length

    def fragment_token_stream(
        self,
        all_tokens: List[str],
        chrom: str,
        window_start: int,
    ) -> List[Fragment]:
        """
        FIX-1: tail fragments truncated by the window boundary are merged
        into their predecessor (or dropped if no predecessor exists) rather
        than being emitted as distribution-violating remainder fragments.
        """
        fragments, pos, total = [], 0, len(all_tokens)
        while pos < total:
            sampled_len = self.sample_fragment_length()
            target_end  = pos + sampled_len
            actual_end  = min(self._find_gc_cut(all_tokens, target_end, window=10), total)
            frag_tokens = all_tokens[pos:actual_end]
            if len(frag_tokens) < FRAG_MIN:
                break
            nt_before     = sum(1 for t in all_tokens[:pos] if len(t) == 1)
            genomic_start = window_start + nt_before

            is_window_truncated = (actual_end == total) and (len(frag_tokens) < sampled_len * 0.8)

            if is_window_truncated and fragments:
                # Merge the truncated remainder into the previous fragment
                # instead of emitting it as a separate, distribution-violating
                # fragment. Strip the previous </cfdna> close tag before
                # appending the remainder + a fresh close tag.
                prev = fragments[-1]
                merged_tokens = prev.tokens[:-1] + frag_tokens + ["</cfdna>"]
                fragments[-1] = Fragment(
                    tokens=merged_tokens,
                    fragment_length=prev.fragment_length + len(frag_tokens),
                    cpg_count=prev.cpg_count + frag_tokens.count("<m>"),
                    strand=prev.strand,
                    genomic_start=prev.genomic_start,
                    chromosome=prev.chromosome,
                )
                pos = actual_end
                continue
            elif is_window_truncated and not fragments:
                # No prior fragment to merge into; drop the malformed remainder
                # rather than emit a fragment that violates the length model.
                break

            fragments.append(Fragment(
                tokens=["<cfdna>"] + frag_tokens + ["</cfdna>"],
                fragment_length=len(frag_tokens),
                cpg_count=frag_tokens.count("<m>"),
                strand="+" if self.rng.random() > 0.5 else "-",
                genomic_start=genomic_start,
                chromosome=chrom,
            ))
            pos = actual_end
        return fragments

    def _find_gc_cut(self, tokens: List[str], target: int, window: int = 10) -> int:
        lo, hi = max(0, target - window), min(len(tokens), target + window)
        if lo >= hi:
            return target
        gc_tokens  = {"G", "C", "<m>", "<um>"}
        best_pos, best_score = target, -1
        for i in range(lo, hi):
            score = sum(1 for t in tokens[max(0, i - 3): min(len(tokens), i + 3)]
                        if t in gc_tokens)
            if score > best_score:
                best_score, best_pos = score, i
        return best_pos

    @staticmethod
    def filter_by_cpg(fragments: List[Fragment]) -> List[Fragment]:
        """
        FIX-2: counts total CpG loci (<m> + <um>), not just methylated ones,
        matching Rule 8's definition of "CpG site," not "methylated CpG site."
        """
        def total_cpg_loci(frag: Fragment) -> int:
            return frag.tokens.count("<m>") + frag.tokens.count("<um>")

        return [f for f in fragments
                if MIN_CPG_PER_FRAG <= total_cpg_loci(f) <= MAX_CPG_PER_FRAG]


# ---------------------------------------------------------------------------
# ctDNA Signal Mixing
# ---------------------------------------------------------------------------

def mix_ctdna_fragments(
    tumor_fragments: List[Fragment],
    healthy_fragments: List[Fragment],
    rng: np.random.Generator,
    enriched: bool = True,
) -> Tuple[List[Fragment], float, float]:
    """
    FIX-3: Healthy/tumor counts are derived FROM the tumor pool, not from
    the combined pool size, so enriched_frac is the actual realized ratio
    rather than an upper bound that healthy fragments can drown out.

    Returns (mixed, raw_ctdna_fraction, enriched_ctdna_fraction).
    Both values are stored in the JSON so metadata matches actual data.
    """
    if not tumor_fragments:
        return healthy_fragments, 0.0, 0.0

    raw_ctdna     = float(rng.uniform(CTDNA_MEDIAN, CTDNA_RANGE_MAX))
    enriched_frac = float(rng.uniform(ENRICHED_CTDNA_MIN, ENRICHED_CTDNA_MAX)) \
                    if enriched else raw_ctdna

    n_tumor   = len(tumor_fragments)
    n_healthy = max(1, int(round(n_tumor * (1.0 - enriched_frac) / enriched_frac)))
    n_healthy = min(n_healthy, max(len(healthy_fragments), 1))

    sel_tumor   = rng.choice(n_tumor, size=n_tumor, replace=False)
    sel_healthy = rng.choice(
        len(healthy_fragments),
        size=min(n_healthy, len(healthy_fragments)),
        replace=len(healthy_fragments) < n_healthy,
    )

    mixed = [tumor_fragments[i] for i in sel_tumor] + \
            [healthy_fragments[i] for i in sel_healthy]
    rng.shuffle(mixed)

    realized_frac = n_tumor / len(mixed) if mixed else 0.0
    return mixed, raw_ctdna, realized_frac
