"""
SOAP-based deduplication for optimized transition states.
"""

from __future__ import annotations

import json
import sqlite3
from collections import Counter
from contextlib import contextmanager
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import threading

import numpy as np
from ase import Atoms
from dscribe.descriptors import SOAP
from dscribe.kernels import AverageKernel

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
        self._lock_path = self.db_path.with_suffix(self.db_path.suffix + ".lock")
        self._soap_lock = threading.Lock()
        self._init_db()
        self.soap_cache: Dict[str, SOAP] = {}

    @contextmanager
    def _db_lock(self):
        """Simple file lock to serialize DB access across processes."""
        try:
            import fcntl
        except ImportError:
            # Non-POSIX; no locking available, yield directly.
            yield
            return
        self._lock_path.parent.mkdir(parents=True, exist_ok=True)
        with open(self._lock_path, "a+") as lock_fh:
            fcntl.flock(lock_fh, fcntl.LOCK_EX)
            try:
                yield
            finally:
                fcntl.flock(lock_fh, fcntl.LOCK_UN)

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
        with self._soap_lock:
            if comp in self.soap_cache:
                return self.soap_cache[comp]
            species = sorted(set(atoms.get_chemical_symbols()))
            soap = SOAP(
                species=species,
                r_cut=self.settings.r_cut,
                n_max=self.settings.n_max,
                l_max=self.settings.l_max,
                average=self.settings.average_mode,
            )
            self.soap_cache[comp] = soap
            return soap

    def fingerprint(self, atoms: Atoms) -> np.ndarray:
        soap = self._descriptor(atoms)
        vec = soap.create(atoms, n_jobs=1)
        # Ensure 1D vector
        return np.array(vec, dtype=np.float32).ravel()

    def _load_records(self, composition: str) -> List[Dict]:
        with self._db_lock():
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
        with self._db_lock():
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

    def check_duplicate(
        self, atoms: Atoms
    ) -> Tuple[bool, Optional[str], Optional[float], int, Optional[str]]:
        """Return (is_duplicate, matched_source_path, best_similarity, existing_count, debug_message)."""
        comp = composition_key(atoms)
        vec = self.fingerprint(atoms).ravel()
        existing = self._load_records(comp)
        existing_count = len(existing)
        if not existing:
            return False, None, None, existing_count, None

        existing_vecs = [r["vector"].ravel() for r in existing]
        # Dimension mismatch guard: if vectors differ in length, skip similarity to avoid crashes.
        if any(v.shape != vec.shape for v in existing_vecs):
            return False, "__incompatible__", None, existing_count, "fingerprint dimension mismatch"

        matrix = np.atleast_2d(np.vstack(existing_vecs))
        vec_2d = vec.reshape(1, -1)

        # Laplacian average kernel for similarity scoring
        kernel = AverageKernel(metric="laplacian", gamma=self.settings.kernel_gamma)
        try:
            sims = kernel.create(vec_2d, matrix)[0]
        except Exception as exc:
            return False, None, None, existing_count, f"similarity computation failed: {exc}"
        best_idx = int(np.nanargmax(sims))
        best_sim = sims[best_idx]
        if best_sim > self.settings.threshold_similarity:
            return True, existing[best_idx]["source"], float(best_sim), existing_count, None
        return False, existing[best_idx]["source"], float(best_sim), existing_count, None

    def register(self, atoms: Atoms, source: Optional[str] = None, metadata: Optional[Dict] = None) -> int:
        """Persist a new fingerprint after successful verification. Returns total count for this composition."""
        meta = metadata or {}
        comp = composition_key(atoms)
        vec = self.fingerprint(atoms)
        self._store_record(comp, vec, source, meta)
        return len(self._load_records(comp))
