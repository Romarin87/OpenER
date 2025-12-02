"""
Lightweight IO helpers for geometry files and conversions.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path
from typing import Iterable, List, Sequence

from ase import Atoms
from ase.io import read, write

_OPI_SRC = Path(__file__).resolve().parents[1] / "opi" / "src"
if _OPI_SRC.exists() and str(_OPI_SRC) not in sys.path:
    sys.path.insert(0, str(_OPI_SRC))

try:  # pragma: no cover - optional dependency made available at runtime
    from opi.input.structures import Structure  # type: ignore
except Exception:  # noqa: BLE001
    Structure = None  # type: ignore


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


def read_last_xyz_frame(path: str | os.PathLike) -> Atoms:
    """Read the last frame from an xyz file (ORCA writes single-frame xyz by default)."""
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(path)
    atoms = read(p, index=-1)
    # Some readers may return a list if index handling changes; guard for safety.
    if isinstance(atoms, list):
        atoms = atoms[-1]
    return atoms


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


def structure_to_atoms(structure: Structure) -> Atoms:
    """Convert an OPI Structure into an ASE Atoms object."""
    if Structure is None:
        raise ImportError("OPI is not available; cannot convert Structure to Atoms")
    symbols: List[str] = []
    coords: List[Sequence[float]] = []
    for atom in structure.atoms:
        elem = getattr(atom, "element", None)
        if elem is None:
            raise ValueError("Structure atom missing element information")
        symbols.append(getattr(elem, "symbol", str(elem)))
        coord_obj = getattr(atom, "coordinates", None)
        if coord_obj is None:
            raise ValueError("Structure atom missing coordinates")
        coords.append(coord_obj.to_list())
    atoms = Atoms(symbols=symbols, positions=coords)
    atoms.info["charge"] = getattr(structure, "charge", 0)
    atoms.info["multiplicity"] = getattr(structure, "multiplicity", 1)
    return atoms


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
