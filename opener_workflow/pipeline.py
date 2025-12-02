"""
End-to-end workflow orchestration for TS optimization and verification.
"""

from __future__ import annotations

import argparse
import json
import logging
from dataclasses import asdict, dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional
from concurrent.futures import ThreadPoolExecutor, as_completed

from ase import Atoms
from tqdm import tqdm

from .analysis import compare_endpoints, is_minimum, is_valid_saddle_point, read_frequencies
from .config import PipelineConfig
from .dedup import SOAPDeduplicator
from .io_utils import (
    ensure_dir,
    read_orca_input_geometry,
    read_xyz_frames,
)
from .orca_runner import OpiJobError, OpiRunner

logger = logging.getLogger("opener.workflow")


@dataclass
class PipelineResult:
    label: str
    status: str
    detail: str
    outputs: Dict[str, str] = field(default_factory=dict)
    metadata: Dict = field(default_factory=dict)


def _timestamped_workdir(base: str | Path = "runs") -> Path:
    """Return a timestamped working directory path."""
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    return Path(base) / stamp


class TransitionStatePipeline:
    """Coordinates all workflow stages."""

    def __init__(
        self,
        config: Optional[PipelineConfig] = None,
        workdir: str | Path | None = None,
        db_path: str | Path = "data/soap_db.sqlite",
        isomeric_smiles: bool = True,
    ):
        self.cfg = config or PipelineConfig()
        self.workdir = ensure_dir(workdir or _timestamped_workdir())
        self.dedup = SOAPDeduplicator(self.cfg.soap, Path(db_path))
        self.runner = OpiRunner(self.cfg.opi)
        self.isomeric_smiles = isomeric_smiles
        self.max_workers = max(1, self.cfg.max_workers)

    def _process_atoms(
        self, atoms: Atoms, label: str, charge: int = 0, mult: int = 1
    ) -> PipelineResult:
        job_dir = ensure_dir(self.workdir / label)
        ts_dir = ensure_dir(job_dir / "TS_opt")
        irc_dir = ensure_dir(job_dir / "IRC")
        rp_dir = ensure_dir(job_dir / "RP_opt")
        outputs: Dict[str, str] = {}
        metadata: Dict = {}
        logs: List[str] = []

        def _log(msg: str) -> None:
            logger.info("%s: %s", label, msg)
            logs.append(msg)

        def _finalize(res: PipelineResult) -> PipelineResult:
            (job_dir / "run.log").write_text("\n".join(logs))
            (job_dir / "result.json").write_text(json.dumps(asdict(res), indent=2))
            return res

        try:
            _log("Starting TS optimization")
            ts_output, ts_xyz, ts_atoms = self.runner.optimize_ts(
                atoms, job_name="ts_opt", workdir=ts_dir, charge=charge, mult=mult
            )
            outputs.update({"ts_out": str(ts_output.get_outfile()), "ts_xyz": str(ts_xyz)})
        except OpiJobError as exc:
            _log(f"TS optimization failed: {exc}")
            return _finalize(PipelineResult(label, "failed", f"TS optimization failed: {exc}", outputs))

        freqs = read_frequencies(ts_output)
        metadata["ts_freqs"] = freqs
        if not is_valid_saddle_point(freqs, self.cfg.freq):
            _log("Failed saddle-point check")
            return _finalize(
                PipelineResult(label, "not_saddle", "Failed saddle-point check", outputs, metadata)
            )

        is_dup, match = self.dedup.check_duplicate(ts_atoms)
        if is_dup:
            msg = f"Duplicate of {match}" if match else "Duplicate structure"
            _log(msg)
            return _finalize(PipelineResult(label, "duplicate", msg, outputs, metadata))

        try:
            _log("Running IRC")
            irc_output, back_xyz, forward_xyz, irc_reactant, irc_product = self.runner.run_irc(
                ts_atoms, job_name="irc", workdir=irc_dir, charge=charge, mult=mult
            )
            outputs.update(
                {
                    "irc_out": str(irc_output.get_outfile()),
                    "irc_backward_xyz": str(back_xyz),
                    "irc_forward_xyz": str(forward_xyz),
                }
            )
        except Exception as exc:  # noqa: BLE001
            _log(f"IRC failed: {exc}")
            return _finalize(PipelineResult(label, "failed", f"IRC failed: {exc}", outputs, metadata))

        try:
            _log("Optimizing IRC endpoints")
            opt_r_output, opt_r_xyz, opt_r_atoms = self.runner.optimize_minimum(
                irc_reactant, job_name="reactant", workdir=rp_dir, charge=charge, mult=mult
            )
            opt_p_output, opt_p_xyz, opt_p_atoms = self.runner.optimize_minimum(
                irc_product, job_name="product", workdir=rp_dir, charge=charge, mult=mult
            )
            outputs.update(
                {
                    "reactant_out": str(opt_r_output.get_outfile()),
                    "product_out": str(opt_p_output.get_outfile()),
                    "reactant_xyz": str(opt_r_xyz),
                    "product_xyz": str(opt_p_xyz),
                }
            )
        except OpiJobError as exc:
            _log(f"Endpoint optimization failed: {exc}")
            return _finalize(
                PipelineResult(label, "failed", f"Endpoint optimization failed: {exc}", outputs, metadata)
            )

        freqs_r = read_frequencies(opt_r_output)
        freqs_p = read_frequencies(opt_p_output)
        metadata["reactant_freqs"] = freqs_r
        metadata["product_freqs"] = freqs_p
        if not is_minimum(freqs_r, self.cfg.freq) or not is_minimum(freqs_p, self.cfg.freq):
            _log("Endpoint minima have imaginary modes")
            return _finalize(
                PipelineResult(
                    label,
                    "failed",
                    "Endpoint minima have imaginary modes",
                    outputs,
                    metadata,
                )
            )

        ok, s_irc_r, s_irc_p, s_opt_r, s_opt_p = compare_endpoints(
            irc_reactant, irc_product, opt_r_atoms, opt_p_atoms, isomeric=self.isomeric_smiles
        )
        metadata.update(
            {
                "smiles_irc_reactant": s_irc_r,
                "smiles_irc_product": s_irc_p,
                "smiles_opt_reactant": s_opt_r,
                "smiles_opt_product": s_opt_p,
            }
        )
        if not ok:
            _log("SMILES mismatch between IRC endpoints and optimized minima")
            return _finalize(
                PipelineResult(label, "failed", "SMILES mismatch between IRC and minima", outputs, metadata)
            )

        self.dedup.register(ts_atoms, source=str(job_dir), metadata=metadata)
        _log("Pipeline completed successfully")
        return _finalize(PipelineResult(label, "success", "Completed TS pipeline", outputs, metadata))

    def run_directory(self, ts_dir: str | Path, charge: int = 0, mult: int = 1) -> List[PipelineResult]:
        ts_dir = Path(ts_dir)
        inputs = sorted(list(ts_dir.glob("*.xyz")) + list(ts_dir.glob("*.inp")))
        tasks: List[tuple[str, Atoms]] = []
        results: List[PipelineResult] = []

        for path in inputs:
            try:
                atoms_list = read_xyz_frames(path) if path.suffix == ".xyz" else [read_orca_input_geometry(path)]
            except Exception as exc:  # noqa: BLE001
                results.append(PipelineResult(path.stem, "failed", f"Input parse error: {exc}"))
                continue
            for idx, atoms in enumerate(atoms_list):
                label = f"{path.stem}_frame{idx}"
                tasks.append((label, atoms))

        if self.max_workers <= 1 or len(tasks) <= 1:
            for label, atoms in tqdm(tasks, desc="TS structures"):
                results.append(self._process_atoms(atoms, label, charge=charge, mult=mult))
            return results

        with ThreadPoolExecutor(max_workers=self.max_workers) as executor:
            future_map = {
                executor.submit(self._process_atoms, atoms, label, charge=charge, mult=mult): label
                for label, atoms in tasks
            }
            for fut in tqdm(as_completed(future_map), total=len(future_map), desc="TS structures"):
                results.append(fut.result())
        return results


def _cli() -> None:
    parser = argparse.ArgumentParser(description="TS workflow driver")
    parser.add_argument("--ts-dir", required=True, help="Folder containing TS initial xyz/inp files")
    parser.add_argument(
        "--workdir",
        default=None,
        help="Working directory for ORCA/OPI jobs; default uses runs/<timestamp>",
    )
    parser.add_argument("--db", default="data/soap_db.sqlite", help="SQLite database for SOAP fingerprints")
    parser.add_argument("--charge", type=int, default=0, help="Total molecular charge")
    parser.add_argument("--mult", type=int, default=1, help="Spin multiplicity")
    parser.add_argument("--non-isomeric", action="store_true", help="Ignore stereochemistry in SMILES comparison")
    parser.add_argument("--json", default=None, help="Optional JSON file to write summary results")
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")
    pipeline = TransitionStatePipeline(
        workdir=args.workdir,
        db_path=args.db,
        isomeric_smiles=not args.non_isomeric,
    )
    logger.info("Working directory: %s", pipeline.workdir)
    results = pipeline.run_directory(args.ts_dir, charge=args.charge, mult=args.mult)
    for res in results:
        logger.info("%s: %s - %s", res.label, res.status, res.detail)
    if args.json:
        Path(args.json).write_text(json.dumps([asdict(r) for r in results], indent=2))


if __name__ == "__main__":
    _cli()
