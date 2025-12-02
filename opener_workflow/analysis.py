"""
Analysis helpers: vibrational frequency parsing, saddle point checks, and SMILES comparisons.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import List, Tuple

import io
from ase import Atoms
from ase.io import write
from openbabel import pybel

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
            # Lines with actual frequencies contain "cm". Skip other info (e.g., scaling factor).
            match = re.search(r"([-+]?\d*\.\d+|[-+]?\d+)\s*cm", line, flags=re.IGNORECASE)
            if not match:
                continue
            try:
                freq = float(match.group(1))
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
    if freqs:
        min_freq = min(freqs)
        if min_freq >= -cfg.min_imag_threshold:
            return False
    return True


def is_minimum(freqs: List[float], cfg: FrequencyCheck) -> bool:
    """Return True if no imaginary modes larger than threshold are present."""
    return not any(f < -cfg.min_imag_threshold for f in freqs)


def atoms_to_smiles(atoms: Atoms, isomeric: bool = True) -> str:
    """Convert an ASE Atoms to a canonical SMILES string."""
    buf = io.StringIO()
    write(buf, atoms, format="xyz")
    xyz = buf.getvalue()
    mol = pybel.readstring("xyz", xyz)
    # pybel write options expect a dict in newer OpenBabel versions.
    opts = {"c": None}
    if isomeric:
        opts["i"] = None
    smiles = mol.write("smi", opt=opts).strip()
    return smiles


def compare_endpoints(
    irc_reactant: Atoms,
    irc_product: Atoms,
    opt_reactant: Atoms,
    opt_product: Atoms,
    isomeric: bool = True,
) -> Tuple[bool, str, str, str, str]:
    """
    Return whether SMILES match between IRC endpoints and optimized minima.
    """
    smiles_irc_r = atoms_to_smiles(irc_reactant, isomeric=isomeric)
    smiles_irc_p = atoms_to_smiles(irc_product, isomeric=isomeric)
    smiles_opt_r = atoms_to_smiles(opt_reactant, isomeric=isomeric)
    smiles_opt_p = atoms_to_smiles(opt_product, isomeric=isomeric)
    ok = smiles_irc_r == smiles_opt_r and smiles_irc_p == smiles_opt_p
    return ok, smiles_irc_r, smiles_irc_p, smiles_opt_r, smiles_opt_p
