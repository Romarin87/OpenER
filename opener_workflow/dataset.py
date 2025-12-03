"""
Dataset loading utilities with filtering and XYZ export.
"""

from __future__ import annotations

import argparse
import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator, List, Optional, Sequence
import json
import sqlite3

from ase import Atoms
from tqdm import tqdm

from .dedup import composition_key
from .io_utils import ensure_dir, read_last_xyz_frame, read_orca_input_geometry, read_xyz_frames, write_xyz

logger = logging.getLogger("opener.dataset")


@dataclass
class DatasetFilters:
    """Simple filters applied to structures when loading datasets."""

    include_elements: Optional[Sequence[str]] = None
    exclude_elements: Optional[Sequence[str]] = None
    min_atoms: int = 1
    max_atoms: Optional[int] = None
    max_heavy_atoms: Optional[int] = None  # counts non-hydrogen atoms
    compositions: Optional[Sequence[str]] = None  # restrict to specific composition keys

    def __post_init__(self) -> None:
        self.include_elements = set(self.include_elements) if self.include_elements else None
        self.exclude_elements = set(self.exclude_elements) if self.exclude_elements else None
        self.compositions = set(self.compositions) if self.compositions else None

    def accept(self, atoms: Atoms) -> bool:
        symbols = atoms.get_chemical_symbols()
        n_atoms = len(symbols)
        if n_atoms < self.min_atoms:
            return False
        if self.max_atoms is not None and n_atoms > self.max_atoms:
            return False
        if self.include_elements and not set(symbols).issubset(self.include_elements):
            return False
        if self.exclude_elements and set(symbols).intersection(self.exclude_elements):
            return False
        if self.max_heavy_atoms is not None:
            heavy = [s for s in symbols if s != "H"]
            if len(heavy) > self.max_heavy_atoms:
                return False
        if self.compositions and composition_key(atoms) not in self.compositions:
            return False
        return True


@dataclass
class DatasetEntry:
    """A single structure entry from the dataset."""

    label: str
    atoms: Atoms
    source: str


def _read_structures(path: Path) -> List[Atoms]:
    """Read structures from a single file, supporting xyz and ORCA inp."""
    if path.suffix.lower() == ".xyz":
        return read_xyz_frames(path)
    if path.suffix.lower() == ".inp":
        return [read_orca_input_geometry(path)]
    raise ValueError(f"Unsupported file type: {path.suffix}")


def iter_structures(
    dataset_root: str | Path,
    filters: Optional[DatasetFilters] = None,
    suffixes: Sequence[str] = (".xyz", ".inp"),
    progress: bool = False,
) -> Iterator[DatasetEntry]:
    """
    Yield filtered structures from a dataset directory or single file.
    """
    root = Path(dataset_root)
    filters = filters or DatasetFilters()

    if root.is_file():
        files = [root]
    else:
        files: List[Path] = []
        for suf in suffixes:
            files.extend(sorted(root.glob(f"**/*{suf}")))
    iterator = files
    if progress:
        iterator = tqdm(files, desc="Loading structures")

    for path in iterator:
        try:
            structures = _read_structures(path)
        except Exception as exc:  # noqa: BLE001
            logger.warning("Skip %s: %s", path, exc)
            continue
        multi = len(structures) > 1
        for idx, atoms in enumerate(structures):
            label = f"{path.stem}_frame{idx}" if multi else path.stem
            if filters.accept(atoms):
                yield DatasetEntry(label=label, atoms=atoms, source=str(path))


def save_as_xyz(entries: Iterable[DatasetEntry], out_dir: str | Path, overwrite: bool = False) -> List[Path]:
    """Persist dataset entries to xyz files; returns written paths."""
    out_root = ensure_dir(out_dir)
    written: List[Path] = []
    for entry in entries:
        base = entry.label
        out_path = out_root / f"{base}.xyz"
        if out_path.exists() and not overwrite:
            suffix = 1
            while (out_root / f"{base}_{suffix}.xyz").exists():
                suffix += 1
            out_path = out_root / f"{base}_{suffix}.xyz"
        comment = f"source: {entry.source}"
        write_xyz(entry.atoms, out_path, comment=comment)
        written.append(out_path)
    return written


def _ts_xyz_from_job_dir(job_dir: Path) -> Optional[Path]:
    """
    Locate the TS xyz for a pipeline job directory.

    Prefer result.json -> outputs.ts_xyz; otherwise require TS_opt/ts_opt.xyz.
    """
    res_path = job_dir / "result.json"
    if res_path.exists():
        try:
            data = json.loads(res_path.read_text())
            ts_xyz_val = data.get("outputs", {}).get("ts_xyz")
            if ts_xyz_val:
                ts_xyz_path = Path(ts_xyz_val)
                if not ts_xyz_path.is_absolute():
                    ts_xyz_path = job_dir / ts_xyz_path
                if ts_xyz_path.exists():
                    return ts_xyz_path
        except Exception:  # noqa: BLE001
            pass

    ts_dir = job_dir / "TS_opt"
    if ts_dir.exists():
        preferred = ts_dir / "ts_opt.xyz"
        if preferred.exists():
            return preferred
    return None


def iter_db_structures(
    db_path: str | Path,
    filters: Optional[DatasetFilters] = None,
    compositions: Optional[Sequence[str]] = None,
    limit: Optional[int] = None,
    progress: bool = False,
) -> Iterator[DatasetEntry]:
    """
    Yield TS structures stored in the SOAP SQLite DB by following source paths to xyz files.
    """
    filters = filters or DatasetFilters()
    query = "SELECT id, composition, source_path FROM soap_entries"
    params: List = []
    if compositions:
        placeholders = ",".join("?" for _ in compositions)
        query += f" WHERE composition IN ({placeholders})"
        params.extend(compositions)
    query += " ORDER BY id"
    if limit is not None:
        query += " LIMIT ?"
        params.append(limit)

    rows: List[tuple] = []
    with sqlite3.connect(db_path) as conn:
        cur = conn.execute(query, params)
        rows = cur.fetchall()

    iterator = rows
    if progress:
        iterator = tqdm(rows, desc="Reading DB entries")

    for row_id, comp, source_path in iterator:
        if not source_path:
            logger.warning("DB entry %s missing source_path; skipping", row_id)
            continue
        job_dir = Path(source_path)
        ts_xyz = _ts_xyz_from_job_dir(job_dir)
        if ts_xyz is None:
            logger.warning("DB entry %s: TS xyz not found under %s", row_id, job_dir)
            continue
        try:
            atoms = read_last_xyz_frame(ts_xyz)
        except Exception as exc:  # noqa: BLE001
            logger.warning("DB entry %s: failed to read %s (%s)", row_id, ts_xyz, exc)
            continue
        if not filters.accept(atoms):
            continue
        label = job_dir.name or f"entry{row_id}"
        # Ensure stable label; append composition to avoid collisions.
        label = f"{label}_{comp}"
        yield DatasetEntry(label=label, atoms=atoms, source=str(ts_xyz))


def _cli() -> None:
    parser = argparse.ArgumentParser(description="Load and filter structures, save as XYZ files.")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--data", help="Dataset directory or single xyz/inp file")
    group.add_argument("--db", help="SOAP SQLite database path (soap_db.sqlite)")
    parser.add_argument("--out", required=True, help="Output directory for filtered xyz files")
    parser.add_argument("--include-elements", nargs="*", default=None, help="Only keep structures using these elements")
    parser.add_argument("--exclude-elements", nargs="*", default=None, help="Drop structures containing these elements")
    parser.add_argument("--min-atoms", type=int, default=1, help="Minimum atom count to keep")
    parser.add_argument("--max-atoms", type=int, default=None, help="Maximum atom count to keep")
    parser.add_argument("--max-heavy-atoms", type=int, default=None, help="Maximum heavy (non-H) atom count")
    parser.add_argument(
        "--compositions",
        nargs="*",
        default=None,
        help="Optional list of composition keys (e.g., C2H6O1) to keep",
    )
    parser.add_argument(
        "--suffixes",
        nargs="*",
        default=[".xyz", ".inp"],
        help="File suffixes to scan under the dataset directory",
    )
    parser.add_argument(
        "--db-compositions",
        nargs="*",
        default=None,
        help="Optional composition keys (e.g., C2H6O1) to select from DB",
    )
    parser.add_argument("--db-limit", type=int, default=None, help="Limit number of DB entries")
    parser.add_argument("--overwrite", action="store_true", help="Overwrite existing xyz files if present")
    parser.add_argument("--no-progress", action="store_true", help="Disable progress bar")
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")

    filt = DatasetFilters(
        include_elements=args.include_elements,
        exclude_elements=args.exclude_elements,
        min_atoms=args.min_atoms,
        max_atoms=args.max_atoms,
        max_heavy_atoms=args.max_heavy_atoms,
        compositions=args.compositions,
    )

    if args.db:
        entries = list(
            iter_db_structures(
                args.db,
                filters=filt,
                compositions=args.db_compositions,
                limit=args.db_limit,
                progress=not args.no_progress,
            )
        )
    else:
        entries = list(
            iter_structures(
                args.data,
                filters=filt,
                suffixes=args.suffixes,
                progress=not args.no_progress,
            )
        )
    logger.info("Loaded %d structures after filtering", len(entries))
    written = save_as_xyz(entries, args.out, overwrite=args.overwrite)
    logger.info("Saved %d xyz files to %s", len(written), args.out)


if __name__ == "__main__":
    _cli()
