"""
Lightweight IO helpers for ORCA jobs.
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import Iterable, List, Sequence

from ase import Atoms
from ase.io import read, write


def read_xyz_frames(path: str | os.PathLike) -> List[Atoms]:
    """Read one or multiple structures from an .xyz file."""
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(path)
    # ":" reads all frames if multiple structures are concatenated
    atoms_list = read(p, index=":")
    if isinstance(atoms_list, Atoms):
        atoms_list = [atoms_list]
    return atoms_list


def write_xyz(atoms: Atoms, path: str | os.PathLike) -> None:
    """Write a single geometry to xyz."""
    write(path, atoms, format="xyz")


def atoms_to_orca_input(
    atoms: Atoms,
    keywords: str,
    blocks: Sequence[str],
    charge: int = 0,
    mult: int = 1,
) -> str:
    """
    Build the contents of an ORCA input file for a given geometry.

    Geometry is embedded as Cartesian coordinates in Angstrom.
    """
    lines: List[str] = [keywords, ""]
    lines.extend(blocks)
    lines.append("")
    lines.append(f"* xyz {charge} {mult}")
    for sym, pos in zip(atoms.get_chemical_symbols(), atoms.positions):
        lines.append(f"{sym:2s} {pos[0]:15.8f} {pos[1]:15.8f} {pos[2]:15.8f}")
    lines.append("*")
    return "\n".join(lines) + "\n"


def read_orca_input_geometry(path: str | os.PathLike) -> Atoms:
    """Extract geometry from an ORCA .inp file."""
    lines = Path(path).read_text().splitlines()
    start = None
    for i, line in enumerate(lines):
        if line.strip().lower().startswith("* xyz"):
            start = i + 1
            break
    if start is None:
        raise ValueError(f"No '* xyz' section found in {path}")
    coords = []
    for line in lines[start:]:
        if line.strip().startswith("*"):
            break
        if not line.strip():
            continue
        parts = line.split()
        if len(parts) < 4:
            continue
        sym = parts[0]
        x, y, z = map(float, parts[-3:])
        coords.append((sym, (x, y, z)))
    if not coords:
        raise ValueError(f"No coordinates parsed from {path}")
    symbols = [c[0] for c in coords]
    positions = [c[1] for c in coords]
    return Atoms(symbols=symbols, positions=positions)


def ensure_dir(path: str | os.PathLike) -> Path:
    """Create directory if missing and return Path."""
    p = Path(path)
    p.mkdir(parents=True, exist_ok=True)
    return p


def chunked(iterable: Iterable, n: int) -> Iterable[List]:
    """Yield successive n-sized chunks."""
    chunk: List = []
    for item in iterable:
        chunk.append(item)
        if len(chunk) == n:
            yield chunk
            chunk = []
    if chunk:
        yield chunk
