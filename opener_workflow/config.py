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
    ts_keywords: List[str] = field(
        default_factory=lambda: _split_keywords("! B3LYP D3BJ def2-SVP OptTS Freq")
    )
    ts_restart_keywords: List[str] = field(
        default_factory=lambda: _split_keywords("TightSCF XQC")
    )
    ts_restart_blocks: List[str] = field(
        default_factory=lambda: ["%method SpecialGridAtoms 1:Grid7 end"]
    )
    opt_keywords: List[str] = field(
        default_factory=lambda: _split_keywords("! B3LYP D3BJ def2-SVP Opt Freq")
    )
    irc_keywords: List[str] = field(
        default_factory=lambda: _split_keywords("! B3LYP D3BJ def2-SVP IRC")
    )
    geom_maxiter: int = 200
    irc_maxiter: int = 50
    max_restarts: int = 2


@dataclass
class SOAPSettings:
    """Parameters for SOAP descriptor generation."""

    r_cut: float = 10.0
    n_max: int = 6
    l_max: int = 4
    average_mode: str = "off"  # per-atom descriptors
    kernel_gamma: float | None = None  # if None, set to 1 / feature_dim
    threshold_similarity: float = 0.9


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
