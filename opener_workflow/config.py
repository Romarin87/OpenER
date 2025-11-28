"""
Default configuration values used across the TS pipeline.
"""

from dataclasses import dataclass, field
from typing import List


@dataclass
class OrcaSettings:
    """Settings for ORCA submissions."""

    executable: str = "orca"
    ts_keywords: str = "! B3LYP D3BJ def2-SVP OptTS Freq"
    common_resources: List[str] = field(
        default_factory=lambda: [
            "%pal nprocs 4 end",
            "%maxcore 4000",
        ]
    )
    geom_block: List[str] = field(default_factory=lambda: ["%geom MaxIter 200 end"])
    ts_restart_blocks: List[str] = field(
        default_factory=lambda: [
            "%scf Convergence Tight end",
            "%scf XQC true end",
            "%geom MaxIter 200 ReStart true end",
            "%method SpecialGridAtoms 1:Grid7 end",
        ]
    )
    irc_keywords: str = "! B3LYP D3BJ def2-SVP IRC"
    irc_block: str = "\n".join(
        [
            "%irc",
            "  Direction Both",
            "  MaxIter 50",
            "  StepSize 0.1",
            "end",
        ]
    )
    irc_block_lqa: str = "\n".join(
        [
            "%irc",
            "  Direction Both",
            "  MaxIter 50",
            "  StepSize 0.1",
            "  Method LQA",
            "end",
        ]
    )
    irc_block_hpc: str = "\n".join(
        [
            "%irc",
            "  Direction Both",
            "  MaxIter 50",
            "  StepSize 0.1",
            "  Method HPC",
            "end",
        ]
    )
    opt_keywords: str = "! B3LYP D3BJ def2-SVP Opt Freq"


@dataclass
class SOAPSettings:
    """Parameters for SOAP descriptor generation."""

    r_cut: float = 6.0
    n_max: int = 8
    l_max: int = 6
    kernel_gamma: float = 1.0
    threshold_similarity: float = 0.999


@dataclass
class FrequencyCheck:
    """Thresholds for saddle point verification."""

    min_imag_threshold: float = 20.0  # cm^-1 absolute value
    expected_imag_count: int = 1


@dataclass
class PipelineConfig:
    """Aggregate config."""

    orca: OrcaSettings = field(default_factory=OrcaSettings)
    soap: SOAPSettings = field(default_factory=SOAPSettings)
    freq: FrequencyCheck = field(default_factory=FrequencyCheck)
