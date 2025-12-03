"""
Default configuration values used across the TS pipeline.
"""

from dataclasses import dataclass, field
from typing import List, Sequence


def _split_keywords(keywords: str | Sequence[str]) -> List[str]:
    """
    Normalize keyword configuration into a list of plain tokens without the leading '!'.
    """
    if isinstance(keywords, str):
        keywords = keywords.replace("!", "").split()
    return [kw.strip() for kw in keywords if kw.strip()]


@dataclass
class OpiSettings:
    """Settings for running ORCA through OPI."""

    orca_path: str | None = None  # Optional override for the ORCA binary path
    mpi_path: str | None = None  # Optional OpenMPI path for Runner
    n_cores: int = 54
    max_core_mb: int = 8000
    method_keywords: List[str] = field(
        default_factory=lambda: _split_keywords("! B3LYP D3BJ def2-SVP")
    )
    geom_maxiter: int | None = 200  # None means use ORCA default (max(3N,50))
    ts_keywords: List[str] = field(
        default_factory=lambda: _split_keywords("OptTS Freq")
    )
    ts_recalc_hess: int | None = 5  # steps between Hessian rebuilds on TS retries; None to disable
    ts_max_restarts: int = 1  # TS optimization restarts
    opt_keywords: List[str] = field(default_factory=lambda: _split_keywords("Opt Freq"))
    irc_keywords: List[str] = field(default_factory=lambda: _split_keywords("IRC"))
    irc_maxiter: int = 50
    irc_max_retries: int = 1  # single-direction IRC retries on SMILES mismatch


@dataclass
class SOAPSettings:
    """
    Parameters for SOAP descriptor generation.

    average_mode options:
    - "off": per-atom descriptors (shape n_atoms x n_features), used by current dedup logic.
    - "outer": average local power spectra to a single global vector.
    - "inner": average neighbor densities first, then build a single global vector.
    - True/False: aliases for "outer"/"off".
    """

    r_cut: float = 6.0
    n_max: int = 8
    l_max: int = 6
    average_mode: str = "off"  # per-atom descriptors
    kernel_gamma: float | None = None  # if None, set to 1 / feature_dim
    threshold_similarity: float = 0.99


@dataclass
class FrequencyCheck:
    """Thresholds for saddle point verification."""

    min_imag_threshold: float = 20.0  # cm^-1 absolute value
    expected_imag_count: int = 1


@dataclass
class PipelineConfig:
    """Aggregate config."""

    opi: OpiSettings = field(default_factory=OpiSettings)
    soap: SOAPSettings = field(default_factory=SOAPSettings)
    freq: FrequencyCheck = field(default_factory=FrequencyCheck)
    max_workers: int = 1  # parallelism for processing multiple TS inputs
