"""
SOAP-based deduplication for optimized transition states.
"""

from __future__ import annotations

import json
import sqlite3
from collections import Counter
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
from ase import Atoms
from dscribe.descriptors import SOAP
from scipy.spatial import cKDTree

from .config import SOAPSettings


def composition_key(atoms: Atoms) -> str:
    """Return a canonical composition string, e.g., C2H6O1."""
    counts = Counter(atoms.get_chemical_symbols())
    parts = [f"{el}{counts[el]}" for el in sorted(counts)]
    return "".join(parts)


class SOAPDeduplicator:
    """Persist SOAP fingerprints and quickly flag duplicates."""

    def __init__(self, settings: SOAPSettings, db_path: Path):
        self.settings = settings
        self.db_path = Path(db_path)
        self._init_db()
        self.soap_cache: Dict[str, SOAP] = {}

    def _init_db(self) -> None:
        self.db_path.parent.mkdir(parents=True, exist_ok=True)
        with sqlite3.connect(self.db_path) as conn:
            conn.execute(
                """
                CREATE TABLE IF NOT EXISTS soap_entries (
                    id INTEGER PRIMARY KEY,
                    composition TEXT NOT NULL,
                    dim INTEGER NOT NULL,
                    fingerprint BLOB NOT NULL,
                    source_path TEXT,
                    metadata TEXT
                )
                """
            )
            conn.execute("CREATE INDEX IF NOT EXISTS idx_comp ON soap_entries(composition)")

    def _descriptor(self, atoms: Atoms) -> SOAP:
        comp = composition_key(atoms)
        if comp in self.soap_cache:
            return self.soap_cache[comp]
        species = sorted(set(atoms.get_chemical_symbols()))
        soap = SOAP(
            species=species,
            rcut=self.settings.r_cut,
            nmax=self.settings.n_max,
            lmax=self.settings.l_max,
            average=True,
        )
        self.soap_cache[comp] = soap
        return soap

    def fingerprint(self, atoms: Atoms) -> np.ndarray:
        soap = self._descriptor(atoms)
        vec = soap.create(atoms, n_jobs=1)
        # Ensure 1D vector
        return np.array(vec, dtype=np.float32).ravel()

    def _load_records(self, composition: str) -> List[Dict]:
        with sqlite3.connect(self.db_path) as conn:
            cur = conn.execute(
                "SELECT id, dim, fingerprint, source_path, metadata FROM soap_entries WHERE composition=?",
                (composition,),
            )
            rows = cur.fetchall()
        records = []
        for rid, dim, blob, source, meta_json in rows:
            vec = np.frombuffer(blob, dtype=np.float32, count=dim)
            meta = json.loads(meta_json) if meta_json else {}
            records.append({"id": rid, "vector": vec, "source": source, "metadata": meta})
        return records

    def _store_record(
        self, composition: str, vec: np.ndarray, source: Optional[str], metadata: Dict
    ) -> None:
        with sqlite3.connect(self.db_path) as conn:
            conn.execute(
                "INSERT INTO soap_entries(composition, dim, fingerprint, source_path, metadata) "
                "VALUES (?, ?, ?, ?, ?)",
                (
                    composition,
                    int(vec.size),
                    sqlite3.Binary(vec.astype(np.float32).tobytes()),
                    source,
                    json.dumps(metadata),
                ),
            )

    def check_duplicate(self, atoms: Atoms) -> Tuple[bool, Optional[str]]:
        """Return (is_duplicate, matched_source_path)."""
        comp = composition_key(atoms)
        vec = self.fingerprint(atoms)
        existing = self._load_records(comp)
        if not existing:
            return False, None

        matrix = np.vstack([r["vector"] for r in existing])
        tree = cKDTree(matrix)
        dist, idx = tree.query(vec, k=min(5, len(existing)))
        if np.isscalar(dist):
            dist = [dist]
            idx = [idx]
        for d, i in zip(dist, idx):
            if np.isfinite(d) and d < self.settings.threshold_distance:
                return True, existing[int(i)]["source"]

        norms = np.linalg.norm(matrix, axis=1) * np.linalg.norm(vec)
        with np.errstate(divide="ignore", invalid="ignore"):
            sims = matrix.dot(vec) / norms
        best_idx = int(np.nanargmax(sims))
        best_sim = sims[best_idx]
        if best_sim > self.settings.threshold_similarity:
            return True, existing[best_idx]["source"]
        return False, None

    def register(self, atoms: Atoms, source: Optional[str] = None, metadata: Optional[Dict] = None) -> None:
        """Persist a new fingerprint after successful verification."""
        meta = metadata or {}
        comp = composition_key(atoms)
        vec = self.fingerprint(atoms)
        self._store_record(comp, vec, source, meta)
