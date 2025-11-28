"""
Helpers to extract IRC endpoints from ORCA outputs.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import List, Tuple

from ase import Atoms

from .orca_runner import parse_orca_geometries


def extract_energies_and_geometries(output_path: Path) -> List[Tuple[float, Atoms]]:
    """
    Greedily pair energies with subsequent geometry blocks.

    This is a heuristic but works for typical ORCA IRC outputs where each step
    prints a 'Total Energy       :    -xxx' line followed by Cartesian
    coordinates.
    """
    text = output_path.read_text(errors="ignore")
    energy_matches = re.findall(r"Total Energy\s*[:=]\s*([-+]?\d+\.\d+)", text, flags=re.IGNORECASE)
    geometries = parse_orca_geometries(output_path)
    pairs: List[Tuple[float, Atoms]] = []
    for e, geom in zip(energy_matches, geometries):
        try:
            energy_val = float(e)
        except Exception:
            continue
        pairs.append((energy_val, geom))
    # If energies are missing, still return geometries with placeholder zeros.
    if not pairs and geometries:
        pairs = [(0.0, g) for g in geometries]
    return pairs


def select_irc_endpoints(output_path: Path) -> Tuple[Atoms, Atoms]:
    """
    Choose representative reactant/product guesses from IRC output.

    Strategy:
    - Split the sequence into two halves (reverse/forward).
    - Pick the lowest-energy geometry in each half.
    """
    pairs = extract_energies_and_geometries(output_path)
    if not pairs:
        raise ValueError(f"No IRC geometries parsed from {output_path}")
    mid = max(1, len(pairs) // 2)
    first_half = pairs[:mid]
    second_half = pairs[mid:] if len(pairs) > 1 else pairs

    def _lowest(seq: List[Tuple[float, Atoms]]) -> Atoms:
        return sorted(seq, key=lambda x: x[0])[0][1]

    reactant_guess = _lowest(first_half)
    product_guess = _lowest(second_half)
    return reactant_guess, product_guess
