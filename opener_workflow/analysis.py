"""
Analysis helpers: vibrational frequency parsing and saddle point checks.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import List

from .config import FrequencyCheck


def read_frequencies(output_path: Path) -> List[float]:
    """Parse vibrational frequencies (cm^-1) from an ORCA output."""
    freqs: List[float] = []
    lines = output_path.read_text(errors="ignore").splitlines()
    reading = False
    for line in lines:
        if "VIBRATIONAL FREQUENCIES" in line.upper():
            reading = True
            continue
        if reading:
            if not line.strip():
                if freqs:
                    break
                continue
            if "NORMAL MODES" in line.upper():
                break
            matches = re.findall(r"[-+]?\d+\.\d+|\d+", line)
            if not matches:
                continue
            try:
                freq = float(matches[-1])
            except ValueError:
                continue
            freqs.append(freq)
    return freqs


def is_valid_saddle_point(freqs: List[float], cfg: FrequencyCheck) -> bool:
    """
    Check whether the frequency list corresponds to a first-order saddle point.
    """
    significant_imag = [f for f in freqs if f < -cfg.min_imag_threshold]
    if len(significant_imag) != cfg.expected_imag_count:
        return False
    if freqs and freqs[0] >= -cfg.min_imag_threshold:
        return False
    return True


def is_minimum(freqs: List[float], cfg: FrequencyCheck) -> bool:
    """Return True if no imaginary modes larger than threshold are present."""
    return not any(f < -cfg.min_imag_threshold for f in freqs)
