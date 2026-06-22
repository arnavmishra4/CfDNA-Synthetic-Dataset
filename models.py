"""
models.py — Data structures for NeuroSight M5 synthetic cfDNA pipeline.
"""

from dataclasses import dataclass
from typing import List


@dataclass
class DMR:
    chromosome: str
    start: int
    end: int
    cpg_ids: List[str]
    delta_beta: float
    t_stat: float
    label: int


@dataclass
class Fragment:
    tokens: List[str]
    fragment_length: int
    cpg_count: int
    strand: str
    genomic_start: int
    chromosome: str
