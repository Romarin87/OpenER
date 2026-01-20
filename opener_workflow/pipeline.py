"""
End-to-end workflow orchestration for TS optimization and verification.
"""

from __future__ import annotations

import argparse
import csv
import json
import logging
import time
from dataclasses import asdict, dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Callable, Dict, List, Optional
from concurrent.futures import ThreadPoolExecutor, as_completed

from ase import Atoms
from tqdm import tqdm

from .analysis import compare_endpoints, is_minimum, is_valid_saddle_point, read_frequencies
from .config import PipelineConfig
from .dedup import SOAPDeduplicator, composition_key
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


@dataclass
class InputRow:
    label: str
    ts_path: Optional[Path]
    reactant_path: Optional[Path]
    product_path: Optional[Path]
    charge: int
    mult: int
    row_index: int


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
        enable_dedup: bool | None = None,
        enable_cineb: bool | None = None,
        isomeric_smiles: bool = True,
    ):
        self.cfg = config or PipelineConfig()
        self.workdir = ensure_dir(workdir or _timestamped_workdir())
        self.enable_dedup = (
            self.cfg.enable_dedup if enable_dedup is None else bool(enable_dedup)
        )
        self.enable_cineb = (
            self.cfg.opi.enable_cineb if enable_cineb is None else bool(enable_cineb)
        )
        self.dedup = (
            SOAPDeduplicator(self.cfg.soap, Path(db_path)) if self.enable_dedup else None
        )
        self.runner = OpiRunner(self.cfg.opi)
        self.isomeric_smiles = isomeric_smiles
        self.max_workers = max(1, self.cfg.max_workers)

    @staticmethod
    def _path_from_cell(value: Optional[str]) -> Optional[Path]:
        if value is None:
            return None
        value = value.strip()
        if not value:
            return None
        path = Path(value)
        if not path.is_absolute():
            path = Path.cwd() / path
        return path

    @staticmethod
    def _label_from_row(
        row: Dict[str, str],
        row_index: int,
        label_col: str,
        ts_path: Optional[Path],
        reactant_path: Optional[Path],
        product_path: Optional[Path],
    ) -> str:
        label_val = row.get(label_col) if label_col else None
        if label_val:
            label_val = label_val.strip()
            if label_val:
                return label_val
        for path in (ts_path, reactant_path, product_path):
            if path is not None:
                return path.stem
        return f"row{row_index}"

    @staticmethod
    def _read_single_structure(path: Path) -> Atoms:
        if path.suffix.lower() == ".xyz":
            frames = read_xyz_frames(path)
            if len(frames) != 1:
                raise ValueError(f"Expected single frame in {path}")
            return frames[0]
        if path.suffix.lower() == ".inp":
            return read_orca_input_geometry(path)
        raise ValueError(f"Unsupported structure file: {path}")

    def _prepare_ts_guess(
        self,
        ts_guess: Optional[Atoms],
        reactant_atoms: Optional[Atoms],
        product_atoms: Optional[Atoms],
        job_dir: Path,
        outputs: Dict[str, str],
        metadata: Dict,
        log: Callable[[str], None],
        charge: int = 0,
        mult: int = 1,
    ) -> Atoms:
        """Optionally run CINEB to refine the TS initial guess."""
        if not self.enable_cineb:
            if ts_guess is None:
                raise OpiJobError("TS guess is required when CINEB is disabled")
            return ts_guess

        if reactant_atoms is None or product_atoms is None:
            raise OpiJobError("CINEB enabled but reactant/product inputs are missing")

        ts_guess_for_cineb = ts_guess
        if ts_guess is not None and not self.cfg.opi.cineb_use_ts_guess:
            log("CINEB TS guess ignored because cineb_use_ts_guess=False")
            ts_guess_for_cineb = None

        cineb_dir = ensure_dir(job_dir / "CINEB")
        log("Starting CINEB refinement")
        cineb_start = time.perf_counter()
        output, ts_xyz, ts_atoms, selected_image, selected_energy = self.runner.run_cineb(
            reactant_atoms,
            product_atoms,
            ts_guess=ts_guess_for_cineb,
            job_name="cineb",
            workdir=cineb_dir,
            charge=charge,
            mult=mult,
            use_ts_guess=self.cfg.opi.cineb_use_ts_guess,
        )
        try:
            cineb_out = output.get_outfile()
        except Exception:  # noqa: BLE001
            cineb_out = cineb_dir / "cineb.out"
        outputs.update(
            {
                "cineb_out": str(cineb_out),
                "cineb_ts_xyz": str(ts_xyz),
            }
        )
        metadata["cineb_selected_image"] = selected_image
        if selected_energy is not None and selected_energy == selected_energy:
            metadata["cineb_selected_energy"] = selected_energy
            log(f"CINEB selected image {selected_image} (E={selected_energy:.8f})")
        else:
            log(f"CINEB selected image {selected_image}")
        log(f"CINEB completed in {time.perf_counter() - cineb_start:.1f}s")
        return ts_atoms

    def _process_row(self, row: InputRow) -> PipelineResult:
        label = row.label
        job_dir = ensure_dir(self.workdir / label)
        ts_dir = ensure_dir(job_dir / "TS_opt")
        irc_dir = ensure_dir(job_dir / "IRC")
        rp_dir = ensure_dir(job_dir / "RP_opt")
        outputs: Dict[str, str] = {}
        metadata: Dict = {
            "csv_row": row.row_index,
            "input": {
                "ts": str(row.ts_path) if row.ts_path else "",
                "reactant": str(row.reactant_path) if row.reactant_path else "",
                "product": str(row.product_path) if row.product_path else "",
            },
        }
        logs: List[str] = []

        def _log(msg: str) -> None:
            logger.info("%s: %s", label, msg)
            logs.append(msg)

        def _finalize(res: PipelineResult) -> PipelineResult:
            (job_dir / "run.log").write_text("\n".join(logs))
            (job_dir / "result.json").write_text(json.dumps(asdict(res), indent=2))
            return res

        try:
            ts_guess = None
            if row.ts_path and (not self.enable_cineb or self.cfg.opi.cineb_use_ts_guess):
                ts_guess = self._read_single_structure(row.ts_path)
            reactant_atoms = (
                self._read_single_structure(row.reactant_path) if row.reactant_path else None
            )
            product_atoms = (
                self._read_single_structure(row.product_path) if row.product_path else None
            )
        except Exception as exc:  # noqa: BLE001
            _log(f"Input parse error: {exc}")
            return _finalize(PipelineResult(label, "failed", f"Input parse error: {exc}", outputs, metadata))

        try:
            ts_atoms = self._prepare_ts_guess(
                ts_guess,
                reactant_atoms,
                product_atoms,
                job_dir,
                outputs,
                metadata,
                _log,
                charge=row.charge,
                mult=row.mult,
            )
        except OpiJobError as exc:
            _log(f"CINEB failed: {exc}")
            return _finalize(PipelineResult(label, "failed", f"CINEB failed: {exc}", outputs, metadata))
        except Exception as exc:  # noqa: BLE001
            _log(f"CINEB failed: {exc}")
            return _finalize(PipelineResult(label, "failed", f"CINEB failed: {exc}", outputs, metadata))

        try:
            ts_start = time.perf_counter()
            _log("Starting TS optimization")
            ts_output, ts_xyz, ts_atoms = self.runner.optimize_ts(
                ts_atoms, job_name="ts_opt", workdir=ts_dir, charge=row.charge, mult=row.mult
            )
            outputs.update({"ts_out": str(ts_output.get_outfile()), "ts_xyz": str(ts_xyz)})
            _log(f"TS optimization completed in {time.perf_counter() - ts_start:.1f}s")
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

        comp = composition_key(ts_atoms)
        if self.enable_dedup and self.dedup:
            dedup_start = time.perf_counter()
            is_dup, match, best_sim, existing_count, sim_debug = self.dedup.check_duplicate(ts_atoms)
            if match == "__incompatible__":
                msg = (
                    "SOAP fingerprint dimension mismatch with DB; "
                    "current settings appear incompatible with existing entries. "
                    "Skipping registration to avoid corruption."
                )
                _log(msg)
                return _finalize(PipelineResult(label, "failed", msg, outputs, metadata))
            if best_sim is not None:
                _log(
                    f"SOAP dedup check in {time.perf_counter() - dedup_start:.2f}s; "
                    f"comp={comp}, existing={existing_count}, max similarity {best_sim:.4f}, "
                    f"threshold {self.cfg.soap.threshold_similarity}"
                    + (f" vs {match}" if match else "")
                )
            else:
                _log(
                    f"SOAP dedup check in {time.perf_counter() - dedup_start:.2f}s; "
                    f"comp={comp}, existing={existing_count}, similarity unavailable"
                    + (f" ({sim_debug})" if sim_debug else "")
                )
            metadata["soap_best_similarity"] = best_sim
            if is_dup:
                msg = (
                    f"Duplicate of {match} (sim={best_sim:.4f})"
                    if match and best_sim is not None
                    else "Duplicate structure"
                )
                _log(msg)
                return _finalize(PipelineResult(label, "duplicate", msg, outputs, metadata))
            if best_sim is not None:
                _log(f"SOAP max similarity {best_sim:.4f}" + (f" vs {match}" if match else ""))
        else:
            _log("Dedup disabled; skipping SOAP check")

        def _run_irc_and_endpoints(irc_maxiter: int):
            irc_start = time.perf_counter()
            _log(f"Running IRC (maxiter={irc_maxiter}, direction=both)")
            irc_output, back_xyz, forward_xyz, irc_reactant, irc_product = self.runner.run_irc(
                ts_atoms,
                job_name="irc",
                workdir=irc_dir,
                charge=row.charge,
                mult=row.mult,
                irc_maxiter=irc_maxiter,
                recalc_hess=self.cfg.opi.irc_recalc_hess,
            )
            outputs.update(
                {
                    "irc_out": str(irc_output.get_outfile()),
                    "irc_backward_xyz": str(back_xyz) if back_xyz else "",
                    "irc_forward_xyz": str(forward_xyz) if forward_xyz else "",
                }
            )
            _log(f"IRC completed in {time.perf_counter() - irc_start:.1f}s")

            opt_start = time.perf_counter()
            _log("Optimizing IRC endpoints")
            opt_r_output, opt_r_xyz, opt_r_atoms = self.runner.optimize_minimum(
                irc_reactant, job_name="reactant", workdir=rp_dir, charge=row.charge, mult=row.mult
            )
            opt_p_output, opt_p_xyz, opt_p_atoms = self.runner.optimize_minimum(
                irc_product, job_name="product", workdir=rp_dir, charge=row.charge, mult=row.mult
            )
            outputs.update(
                {
                    "reactant_out": str(opt_r_output.get_outfile()),
                    "product_out": str(opt_p_output.get_outfile()),
                    "reactant_xyz": str(opt_r_xyz),
                    "product_xyz": str(opt_p_xyz),
                }
            )
            _log(f"Endpoint optimizations completed in {time.perf_counter() - opt_start:.1f}s")

            freqs_r = read_frequencies(opt_r_output)
            freqs_p = read_frequencies(opt_p_output)
            metadata["reactant_freqs"] = freqs_r
            metadata["product_freqs"] = freqs_p
            if not is_minimum(freqs_r, self.cfg.freq) or not is_minimum(freqs_p, self.cfg.freq):
                return False, "Endpoint minima have imaginary modes", None, irc_reactant, irc_product, opt_r_atoms, opt_p_atoms

            ok_smiles, s_irc_r, s_irc_p, s_opt_r, s_opt_p = compare_endpoints(
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
            if not ok_smiles:
                return False, "SMILES mismatch between IRC and minima", "smiles", irc_reactant, irc_product, opt_r_atoms, opt_p_atoms
            return True, "ok", None, irc_reactant, irc_product, opt_r_atoms, opt_p_atoms

        # First attempt: both directions
        try:
            success, detail, tag, irc_reactant, irc_product, opt_r_atoms, opt_p_atoms = _run_irc_and_endpoints(
                self.cfg.opi.irc_maxiter
            )
        except OpiJobError as exc:
            _log(f"Endpoint optimization failed: {exc}")
            return _finalize(PipelineResult(label, "failed", f"Endpoint optimization failed: {exc}", outputs, metadata))
        except Exception as exc:  # noqa: BLE001
            _log(f"IRC failed: {exc}")
            return _finalize(PipelineResult(label, "failed", f"IRC failed: {exc}", outputs, metadata))

        if not success and tag == "smiles":
            # Check which side mismatched and rerun single-direction IRC with higher maxiter for that side.
            mismatch_reactant = metadata.get("smiles_irc_reactant") != metadata.get("smiles_opt_reactant")
            mismatch_product = metadata.get("smiles_irc_product") != metadata.get("smiles_opt_product")
            last_retry_maxiter: float | None = None

            if mismatch_reactant:
                for retry_idx in range(self.cfg.opi.irc_max_retries):
                    irc_maxiter_retry = max(self.cfg.opi.irc_maxiter * (retry_idx + 2), self.cfg.opi.irc_maxiter + 1)
                    last_retry_maxiter = irc_maxiter_retry
                    try:
                        _log(
                            f"SMILES mismatch on reactant; rerunning backward IRC (attempt {retry_idx + 1}/{self.cfg.opi.irc_max_retries}) with maxiter={irc_maxiter_retry}"
                        )
                        _, back_xyz, _, irc_reactant_new, _ = self.runner.run_irc(
                            ts_atoms,
                            job_name=f"irc_retry_back{retry_idx + 1}",
                            workdir=irc_dir,
                            charge=row.charge,
                            mult=row.mult,
                            irc_maxiter=irc_maxiter_retry,
                            recalc_hess=self.cfg.opi.irc_recalc_hess,
                            direction="backward",
                        )
                        outputs["irc_backward_xyz_retry"] = str(back_xyz) if back_xyz else ""
                        opt_r_output, opt_r_xyz, opt_r_atoms = self.runner.optimize_minimum(
                            irc_reactant_new, job_name=f"reactant_retry{retry_idx + 1}", workdir=rp_dir, charge=row.charge, mult=row.mult
                        )
                        outputs["reactant_out_retry"] = str(opt_r_output.get_outfile())
                        outputs["reactant_xyz_retry"] = str(opt_r_xyz)
                        freqs_r = read_frequencies(opt_r_output)
                        metadata["reactant_freqs"] = freqs_r
                        if not is_minimum(freqs_r, self.cfg.freq):
                            _log("Endpoint minima have imaginary modes after backward retry")
                            return _finalize(
                                PipelineResult(
                                    label,
                                    "failed",
                                    "Endpoint minima have imaginary modes after backward retry",
                                    outputs,
                                    metadata,
                                )
                            )
                        irc_reactant = irc_reactant_new
                        break
                    except Exception as exc:  # noqa: BLE001
                        _log(f"Backward IRC retry attempt {retry_idx + 1} failed: {exc}")
                        if retry_idx + 1 >= self.cfg.opi.irc_max_retries:
                            return _finalize(PipelineResult(label, "failed", f"Backward IRC retry failed: {exc}", outputs, metadata))

            if mismatch_product:
                for retry_idx in range(self.cfg.opi.irc_max_retries):
                    irc_maxiter_retry = max(self.cfg.opi.irc_maxiter * (retry_idx + 2), self.cfg.opi.irc_maxiter + 1)
                    last_retry_maxiter = irc_maxiter_retry
                    try:
                        _log(
                            f"SMILES mismatch on product; rerunning forward IRC (attempt {retry_idx + 1}/{self.cfg.opi.irc_max_retries}) with maxiter={irc_maxiter_retry}"
                        )
                        _, _, forward_xyz, _, irc_product_new = self.runner.run_irc(
                            ts_atoms,
                            job_name=f"irc_retry_forward{retry_idx + 1}",
                            workdir=irc_dir,
                            charge=row.charge,
                            mult=row.mult,
                            irc_maxiter=irc_maxiter_retry,
                            recalc_hess=self.cfg.opi.irc_recalc_hess,
                            direction="forward",
                        )
                        outputs["irc_forward_xyz_retry"] = str(forward_xyz) if forward_xyz else ""
                        opt_p_output, opt_p_xyz, opt_p_atoms = self.runner.optimize_minimum(
                            irc_product_new, job_name=f"product_retry{retry_idx + 1}", workdir=rp_dir, charge=row.charge, mult=row.mult
                        )
                        outputs["product_out_retry"] = str(opt_p_output.get_outfile())
                        outputs["product_xyz_retry"] = str(opt_p_xyz)
                        freqs_p = read_frequencies(opt_p_output)
                        metadata["product_freqs"] = freqs_p
                        if not is_minimum(freqs_p, self.cfg.freq):
                            _log("Endpoint minima have imaginary modes after forward retry")
                            return _finalize(
                                PipelineResult(
                                    label,
                                    "failed",
                                    "Endpoint minima have imaginary modes after forward retry",
                                    outputs,
                                    metadata,
                                )
                            )
                        irc_product = irc_product_new
                        break
                    except Exception as exc:  # noqa: BLE001
                        _log(f"Forward IRC retry attempt {retry_idx + 1} failed: {exc}")
                        if retry_idx + 1 >= self.cfg.opi.irc_max_retries:
                            return _finalize(PipelineResult(label, "failed", f"Forward IRC retry failed: {exc}", outputs, metadata))

            # Re-evaluate SMILES with updated endpoints
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
                _log("SMILES mismatch between IRC and minima after directional retries")
                # Cache this failed TS in the dedup DB to skip reruns of the same structure.
                if self.enable_dedup and self.dedup:
                    try:
                        meta_fail = dict(metadata)
                        meta_fail["status"] = "irc_smiles_mismatch"
                        meta_fail["irc_retry_maxiter"] = last_retry_maxiter
                        stored_total = self.dedup.register(
                            ts_atoms, source=str(job_dir), metadata=meta_fail
                        )
                        _log(f"Registered failed TS (smiles mismatch) for comp={comp} (total {stored_total})")
                    except Exception as exc:  # noqa: BLE001
                        _log(f"Failed to register mismatched TS in dedup DB: {exc}")
                return _finalize(
                    PipelineResult(label, "failed", "SMILES mismatch between IRC and minima after directional retries", outputs, metadata)
                )

        elif not success:
            if detail == "Endpoint minima have imaginary modes":
                _log(detail)
                return _finalize(PipelineResult(label, "failed", detail, outputs, metadata))
            return _finalize(PipelineResult(label, "failed", detail, outputs, metadata))

        if self.enable_dedup and self.dedup:
            stored_total = self.dedup.register(ts_atoms, source=str(job_dir), metadata=metadata)
            _log(f"Pipeline completed successfully; stored SOAP entry for comp={comp} (total {stored_total})")
        else:
            _log(f"Pipeline completed successfully; dedup disabled for comp={comp}")
        return _finalize(PipelineResult(label, "success", "Completed TS pipeline", outputs, metadata))

    def run_csv(
        self,
        csv_path: str | Path,
        default_charge: int = 0,
        default_mult: int = 1,
    ) -> List[PipelineResult]:
        csv_path = Path(csv_path)
        if not csv_path.exists():
            raise FileNotFoundError(csv_path)

        cfg = self.cfg.input
        results: List[PipelineResult] = []
        tasks: List[InputRow] = []

        with csv_path.open(newline="") as fh:
            reader = csv.DictReader(fh)
            if not reader.fieldnames:
                raise ValueError("CSV header is missing")
            header = {name.strip() for name in reader.fieldnames if name}

            react_col = cfg.react_col
            prod_col = cfg.prod_col
            ts_col = cfg.ts_col
            label_col = cfg.label_col
            charge_col = cfg.charge_col
            mult_col = cfg.mult_col

            required_cols = [ts_col]
            if self.enable_cineb:
                required_cols = [react_col, prod_col]
            missing_required = [col for col in required_cols if col not in header]
            if missing_required:
                raise ValueError(f"CSV missing required columns: {', '.join(missing_required)}")

            optional_cols = [label_col, charge_col, mult_col]
            if self.enable_cineb:
                optional_cols.append(ts_col)
            missing_optional = [col for col in optional_cols if col not in header]
            for col in missing_optional:
                logger.warning("CSV missing optional column %s", col)

            for row_index, row in enumerate(reader, start=1):
                normalized_row = {
                    (key.strip() if key else ""): (value or "")
                    for key, value in row.items()
                    if key is not None
                }
                react_path = (
                    self._path_from_cell(normalized_row.get(react_col))
                    if react_col in header
                    else None
                )
                prod_path = (
                    self._path_from_cell(normalized_row.get(prod_col))
                    if prod_col in header
                    else None
                )
                ts_path = (
                    self._path_from_cell(normalized_row.get(ts_col)) if ts_col in header else None
                )
                label = self._label_from_row(
                    normalized_row, row_index, label_col, ts_path, react_path, prod_path
                )

                def _fail(detail: str) -> None:
                    results.append(PipelineResult(label, "failed", f"CSV input error: {detail}"))

                if self.enable_cineb:
                    if react_path is None or prod_path is None:
                        _fail("React/Prod paths are required when CINEB is enabled")
                        continue
                    if ts_path is not None and not self.cfg.opi.cineb_use_ts_guess:
                        logger.warning(
                            "%s: TS path provided but cineb_use_ts_guess=False; ignoring TS input",
                            label,
                        )
                        ts_path = None
                else:
                    if ts_path is None:
                        _fail("TS path is required when CINEB is disabled")
                        continue

                for path in (react_path, prod_path, ts_path):
                    if path is not None and not path.exists():
                        _fail(f"Input file not found: {path}")
                        break
                else:
                    charge = default_charge
                    if charge_col in header:
                        val = (normalized_row.get(charge_col) or "").strip()
                        if val:
                            try:
                                charge = int(val)
                            except ValueError:
                                _fail(f"Invalid charge value: {val}")
                                continue

                    mult = default_mult
                    if mult_col in header:
                        val = (normalized_row.get(mult_col) or "").strip()
                        if val:
                            try:
                                mult = int(val)
                            except ValueError:
                                _fail(f"Invalid mult value: {val}")
                                continue

                    tasks.append(
                        InputRow(
                            label=label,
                            ts_path=ts_path,
                            reactant_path=react_path,
                            product_path=prod_path,
                            charge=charge,
                            mult=mult,
                            row_index=row_index,
                        )
                    )

        if self.max_workers <= 1 or len(tasks) <= 1:
            for row in tqdm(tasks, desc="CSV inputs"):
                results.append(self._process_row(row))
            return results

        with ThreadPoolExecutor(max_workers=self.max_workers) as executor:
            future_map = {executor.submit(self._process_row, row): row.label for row in tasks}
            for fut in tqdm(as_completed(future_map), total=len(future_map), desc="CSV inputs"):
                results.append(fut.result())
        return results


def _cli() -> None:
    parser = argparse.ArgumentParser(description="TS workflow driver")
    parser.add_argument("--input-csv", required=True, help="CSV file containing input paths")
    parser.add_argument("--col-react", default=None, help="CSV column name for reactant paths")
    parser.add_argument("--col-prod", default=None, help="CSV column name for product paths")
    parser.add_argument("--col-ts", default=None, help="CSV column name for TS guess paths")
    parser.add_argument("--col-label", default=None, help="CSV column name for labels")
    parser.add_argument("--col-charge", default=None, help="CSV column name for per-row charge")
    parser.add_argument("--col-mult", default=None, help="CSV column name for per-row multiplicity")
    parser.add_argument(
        "--workdir",
        default=None,
        help="Working directory for ORCA/OPI jobs; default uses runs/<timestamp>",
    )
    parser.add_argument("--db", default="data/soap_db.sqlite", help="SQLite database for SOAP fingerprints")
    parser.add_argument("--charge", type=int, default=0, help="Total molecular charge")
    parser.add_argument("--mult", type=int, default=1, help="Spin multiplicity")
    parser.add_argument("--non-isomeric", action="store_true", help="Ignore stereochemistry in SMILES comparison")
    parser.add_argument("--cineb", action="store_true", help="Enable CINEB refinement before TS optimization")
    parser.add_argument("--no-dedup", action="store_true", help="Disable SOAP deduplication")
    parser.add_argument("--json", default=None, help="Optional JSON file to write summary results")
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")
    config = PipelineConfig()
    if args.col_react:
        config.input.react_col = args.col_react
    if args.col_prod:
        config.input.prod_col = args.col_prod
    if args.col_ts:
        config.input.ts_col = args.col_ts
    if args.col_label:
        config.input.label_col = args.col_label
    if args.col_charge:
        config.input.charge_col = args.col_charge
    if args.col_mult:
        config.input.mult_col = args.col_mult

    pipeline = TransitionStatePipeline(
        config=config,
        workdir=args.workdir,
        db_path=args.db,
        enable_dedup=not args.no_dedup,
        enable_cineb=args.cineb,
        isomeric_smiles=not args.non_isomeric,
    )
    logger.info("Working directory: %s", pipeline.workdir)
    results = pipeline.run_csv(args.input_csv, default_charge=args.charge, default_mult=args.mult)
    for res in results:
        logger.info("%s: %s - %s", res.label, res.status, res.detail)
    if args.json:
        Path(args.json).write_text(json.dumps([asdict(r) for r in results], indent=2))


if __name__ == "__main__":
    _cli()
