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
class CsvInputSettings:
    """Column names used when reading CSV inputs."""

    react_col: str = "React"  # default: React; reactant path column
    prod_col: str = "Prod"  # default: Prod; product path column
    ts_col: str = "TS"  # default: TS; TS guess path column
    label_col: str = "Label"  # default: Label; optional label column
    charge_col: str = "Charge"  # default: Charge; optional per-row charge column
    mult_col: str = "Mult"  # default: Mult; optional per-row multiplicity column


@dataclass
class OpiSettings:
    """Settings for running ORCA through OPI."""

    # Environment
    orca_path: str | None = None  # default: None; override ORCA binary path
    mpi_path: str | None = None  # default: None; override OpenMPI path

    # Resources
    n_cores: int = 1  # default: 1; cores per job
    max_core_mb: int = 8000  # default: 8000; memory per core in MB

    # Method and geometry
    method_keywords: List[str] = field(
        default_factory=lambda: _split_keywords("! GFN2-xTB")
    )  # default: GFN2-xTB
    geom_maxiter: int | None = 200  # default: 200; None uses ORCA default (max(3N,50))

    # TS optimization
    ts_keywords: List[str] = field(
        default_factory=lambda: _split_keywords("OptTS Freq")
    )  # default: OptTS Freq
    ts_recalc_hess: int | None = 5  # default: 5; Hessian rebuild interval on TS retries; None disables
    ts_max_restarts: int = 1  # default: 1; TS optimization restarts

    # Endpoint optimization
    opt_keywords: List[str] = field(default_factory=lambda: _split_keywords("Opt Freq"))  # default: Opt Freq

    # IRC
    irc_keywords: List[str] = field(default_factory=lambda: _split_keywords("IRC"))  # default: IRC
    irc_maxiter: int = 50  # default: 50; IRC steps
    irc_max_retries: int = 1  # default: 1; single-direction IRC retries on SMILES mismatch
    irc_recalc_hess: int | None = 5  # default: 5; Hessian rebuild interval during IRC; None disables

    # CINEB (NEB-TS)
    enable_cineb: bool = True  # default: True; enable CINEB pre-step
    cineb_keywords: List[str] = field(default_factory=lambda: _split_keywords("neb-ts"))  # default: neb-ts
    cineb_nimages: int | None = None  # default (ORCA 6.1): 8
    cineb_maxiter: int | None = None  # default (ORCA 6.1): 500 (LBFGS), 1000 (VPO/FIRE)
    cineb_interpolation: str | None = None  # default (ORCA 6.1): IDPP
    cineb_springconst: float | None = None  # default (ORCA 6.1): 0.01 (Eh/Bohr)
    cineb_use_ts_guess: bool = True  # default: True; pass TS guess into NEB if provided


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

    r_cut: float = 6.0  # default: 6.0; SOAP cutoff radius
    n_max: int = 8  # default: 8; radial basis size
    l_max: int = 6  # default: 6; angular basis size
    average_mode: str = "off"  # default: off; per-atom descriptors
    kernel_gamma: float | None = None  # default: None; auto 1/feature_dim
    threshold_similarity: float = 0.99  # default: 0.99; duplicate cutoff


@dataclass
class FrequencyCheck:
    """Thresholds for saddle point verification."""

    min_imag_threshold: float = 20.0  # default: 20.0; cm^-1 absolute value
    expected_imag_count: int = 1  # default: 1


@dataclass
class PipelineConfig:
    """Aggregate config."""

    opi: OpiSettings = field(default_factory=OpiSettings)  # default: OpiSettings()
    soap: SOAPSettings = field(default_factory=SOAPSettings)  # default: SOAPSettings()
    freq: FrequencyCheck = field(default_factory=FrequencyCheck)  # default: FrequencyCheck()
    input: CsvInputSettings = field(default_factory=CsvInputSettings)  # default: CsvInputSettings()
    max_workers: int = 1  # default: 1; parallelism for processing multiple inputs
    enable_dedup: bool = False  # default: False; enable SOAP dedup
