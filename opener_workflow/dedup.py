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

    def _descriptor_matrix(self, atoms: Atoms) -> np.ndarray:
        """
        Return per-atom SOAP descriptors (n_atoms, n_features) as float32.

        The default config uses average='off' to match the provided Laplacian kernel
        logic; if averaging is enabled, we still force a 2D shape for robustness.
        """
        soap = self._descriptor(atoms)
        desc = soap.create(atoms, n_jobs=1)
        arr = np.array(desc, dtype=np.float32)
        if arr.ndim == 1:
            arr = arr.reshape(1, -1)
        return arr

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

    @staticmethod
    def _reshape_descriptor(vec: np.ndarray, meta: Dict) -> Optional[np.ndarray]:
        """Recover (n_atoms, n_features) matrix from stored vector and metadata."""
        shape = meta.get("soap_shape")
        if isinstance(shape, (list, tuple)) and len(shape) == 2:
            n_atoms, n_feat = shape
            try:
                n_atoms_i = int(n_atoms)
                n_feat_i = int(n_feat)
            except Exception:
                return None
            if n_atoms_i * n_feat_i != vec.size or n_atoms_i <= 0 or n_feat_i <= 0:
                return None
            return vec.reshape((n_atoms_i, n_feat_i))
        return None

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

    @staticmethod
    def _laplacian_local_kernel(a: np.ndarray, b: np.ndarray, gamma: float) -> np.ndarray:
        """
        Pairwise Laplacian kernel between atomic environments.

        Args:
            a: (n_atoms_a, n_feat)
            b: (n_atoms_b, n_feat)
        """
        # Broadcast differences, accumulate L1 distance per atom pair, then exponentiate.
        diff = np.abs(a[:, None, :] - b[None, :, :])  # (na, nb, n_feat)
        dist = np.sum(diff, axis=2)  # (na, nb)
        return np.exp(-gamma * dist).astype(np.float32)

    def _self_similarity(self, desc: np.ndarray, gamma: float) -> float:
        """Return sqrt of mean self-kernel to match provided normalization."""
        local = self._laplacian_local_kernel(desc, desc, gamma)
        return float(np.sqrt(np.mean(local)))

    def _normalized_similarity(self, cand: np.ndarray, ref: np.ndarray, gamma: float) -> Optional[float]:
        """
        Compute normalized average Laplacian similarity between two structures.
        """
        try:
            cross = float(np.mean(self._laplacian_local_kernel(cand, ref, gamma)))
            self_cand = self._self_similarity(cand, gamma)
            self_ref = self._self_similarity(ref, gamma)
            denom = self_cand * self_ref
            if denom <= 0 or not np.isfinite(denom):
                return None
            sim = cross / denom
            if not np.isfinite(sim):
                return None
            return sim
        except Exception:
            return None

    def check_duplicate(
        self, atoms: Atoms
    ) -> Tuple[bool, Optional[str], Optional[float], int, Optional[str]]:
        """Return (is_duplicate, matched_source_path, best_similarity, existing_count, debug_message)."""
        comp = composition_key(atoms)
        desc = self._descriptor_matrix(atoms)
        if desc.size == 0 or desc.ndim != 2:
            return False, "__incompatible__", None, 0, "empty SOAP descriptor"

        n_feat = desc.shape[1]
        gamma = self.settings.kernel_gamma if self.settings.kernel_gamma is not None else (1.0 / float(n_feat))
        if gamma <= 0 or not np.isfinite(gamma):
            return False, "__incompatible__", None, 0, "invalid gamma for SOAP kernel"

        existing = self._load_records(comp)
        existing_count = len(existing)
        if not existing:
            return False, None, None, existing_count, None

        sims: List[Tuple[float, Dict]] = []
        skipped = 0
        for rec in existing:
            reshaped = self._reshape_descriptor(rec["vector"], rec["metadata"])
            if reshaped is None or reshaped.shape[1] != n_feat:
                skipped += 1
                continue
            sim = self._normalized_similarity(desc, reshaped, gamma)
            if sim is None:
                skipped += 1
                continue
            sims.append((sim, rec))

        if not sims:
            msg = "no compatible SOAP descriptors"
            if skipped and existing_count:
                msg += f" (skipped {skipped} of {existing_count})"
            return False, None, None, existing_count, msg

        best_sim, best_rec = max(sims, key=lambda x: x[0])
        if best_sim > self.settings.threshold_similarity:
            return True, best_rec["source"], float(best_sim), existing_count, None
        return False, best_rec["source"], float(best_sim), existing_count, None

    def register(self, atoms: Atoms, source: Optional[str] = None, metadata: Optional[Dict] = None) -> int:
        """Persist a new fingerprint after successful verification. Returns total count for this composition."""
        meta = dict(metadata or {})
        comp = composition_key(atoms)
        desc = self._descriptor_matrix(atoms)
        meta["soap_shape"] = [int(dim) for dim in desc.shape]
        self._store_record(comp, desc.ravel(), source, meta)
        return len(self._load_records(comp))
