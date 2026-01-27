"""
OPI-backed helpers to run ORCA calculations used in the workflow.
"""

from __future__ import annotations

import logging
import os
import sys
from pathlib import Path
from typing import Iterable, List, Sequence

from ase import Atoms

from .config import OpiSettings
from .io_utils import ensure_dir, read_last_xyz_frame, structure_to_atoms, write_xyz

# Ensure the bundled OPI source is importable when the package is not installed.
_OPI_SRC = Path(__file__).resolve().parents[1] / "opi" / "src"
if _OPI_SRC.exists() and str(_OPI_SRC) not in sys.path:
    sys.path.insert(0, str(_OPI_SRC))

from opi.core import Calculator  
from opi.input.arbitrary_string import ArbitraryStringPos  
from opi.input.blocks.block_neb import BlockNeb  
from opi.input.blocks.block_irc import BlockIrc  
from opi.input.structures import Structure  
from opi.output.core import Output  

logger = logging.getLogger("opener.orca_runner")


class OpiJobError(RuntimeError):
    """Raised when an ORCA job fails using the OPI wrapper."""


def _normalize_keywords(keywords: Sequence[str] | str) -> List[str]:
    """Return a flat list of keywords without leading exclamation marks."""
    tokens: List[str] = []
    parts: Iterable[str]
    if isinstance(keywords, str):
        parts = keywords.replace("!", " ").split()
    else:
        collected: List[str] = []
        for kw in keywords:
            collected.extend(str(kw).replace("!", " ").split())
        parts = collected

    for kw in parts:
        kw = kw.strip()
        if kw:
            tokens.append(kw)
    return tokens


def _neb_mode_from_keywords(keywords: Sequence[str]) -> str:
    """Infer NEB mode from keywords (neb-ts vs neb-ci)."""
    text = " ".join(_normalize_keywords(keywords)).lower().replace("_", "-")
    if "neb-ts" in text or "nebts" in text:
        return "neb-ts"
    if "neb-ci" in text or "nebci" in text or "cineb" in text:
        return "neb-ci"
    raise OpiJobError("NEB mode keywords must include neb-ts or neb-ci")


def _geom_block_text(maxiter: int | None, restart: bool = False, recalc_hess: int | None = None) -> str:
    """Build a %geom block text; restart flag is kept for compatibility."""
    lines = ["%geom"]
    if maxiter is not None:
        lines.append(f"  MaxIter {maxiter}")
    if restart:
        lines.append("  ReStart true")
    if recalc_hess and recalc_hess > 0:
        lines.append(f"  Recalc_Hess {int(recalc_hess)}")
    lines.append("end")
    return "\n".join(lines)


def _geometry_converged_from_outfile(outfile: Path) -> bool:
    """
    Fallback convergence detection by scanning the ORCA .out file when OPI parsing is incomplete.
    """
    try:
        text = Path(outfile).read_text(errors="ignore").upper()
    except Exception:
        return False
    markers = (
        "GEOMETRY OPTIMIZATION CONVERGED",
        "THE OPTIMIZATION HAS CONVERGED",
        "OPTIMIZATION CONVERGED",
    )
    return any(marker in text for marker in markers)


class OpiRunner:
    """Run ORCA calculations via the OPI Calculator wrapper."""

    def __init__(self, settings: OpiSettings):
        self.settings = settings
        self._configure_environment()

    def _configure_environment(self) -> None:
        """Set environment variables expected by OPI for binary discovery."""
        if self.settings.orca_path:
            os.environ["OPI_ORCA"] = self.settings.orca_path
        if self.settings.mpi_path:
            os.environ["OPI_MPI"] = self.settings.mpi_path

    def _run_calculation(
        self,
        atoms: Atoms,
        job_name: str,
        workdir: Path,
        keywords: Sequence[str],
        *,
        block_strings: Sequence[str] | None = None,
        blocks: Sequence[object] | None = None,
        charge: int = 0,
        mult: int = 1,
    ) -> Output:
        """
        Create an OPI calculator, write the input, run ORCA, and return the parsed Output.
        """
        workdir = ensure_dir(workdir)
        calc = Calculator(basename=job_name, working_dir=workdir)
        calc.structure = Structure.from_ase(atoms, charge=charge, multiplicity=mult)
        calc.input.ncores = self.settings.n_cores
        calc.input.memory = self.settings.max_core_mb

        full_keywords = list(self.settings.method_keywords) + list(keywords)
        for kw in _normalize_keywords(full_keywords):
            calc.input.add_simple_keywords(kw)

        for block_str in block_strings or ():
            calc.input.add_arbitrary_string(block_str, pos=ArbitraryStringPos.TOP)

        for block in blocks or ():
            calc.input.add_blocks(block)

        try:
            calc.write_input()
            ok = calc.run()
            output = calc.get_output()
        except Exception as exc:  # noqa: BLE001
            # Continue to surface failures where ORCA itself did not finish.
            raise OpiJobError(f"Failed to run ORCA job {job_name}: {exc}") from exc

        try:
            outfile = output.get_outfile()
        except Exception:  # noqa: BLE001
            outfile = workdir / f"{job_name}.out"

        if not ok or not output.terminated_normally():
            raise OpiJobError(f"ORCA job {job_name} failed; inspect {Path(outfile).name}")

        # Try to ensure JSONs exist for downstream parsing; tolerate failures.
        try:
            if output.results_properties is None:
                output.parse(do_create_property_json=True, do_create_gbw_json=False)
        except FileNotFoundError:
            try:
                output.parse(do_create_property_json=True, do_create_gbw_json=True)
            except Exception:  # noqa: BLE001
                short_path = Path(outfile).name
                logger.warning("Parse JSON for %s failed; fallback to text output %s", job_name, short_path)
        except Exception:  # noqa: BLE001
            short_path = Path(outfile).name
            logger.warning("Parse JSON for %s failed; fallback to text output %s", job_name, short_path)
        return output

    def _final_atoms_from_output(self, output: Output, workdir: Path, job_label: str) -> tuple[Path, Atoms]:
        """Return final atoms and ensure an xyz exists for downstream steps."""
        xyz_path = Path(workdir) / f"{job_label}.xyz"
        if xyz_path.exists():
            atoms = read_last_xyz_frame(xyz_path)
        else:
            structure = output.get_structure()
            if structure is None:
                raise OpiJobError(f"No final structure found for {job_label}; see {output.get_outfile()}")
            atoms = structure_to_atoms(structure)
            write_xyz(atoms, xyz_path)
        return xyz_path, atoms

    def run_cineb(
        self,
        reactant: Atoms,
        product: Atoms,
        ts_guess: Atoms | None,
        job_name: str,
        workdir: Path,
        charge: int = 0,
        mult: int = 1,
        nimages: int | None = None,
        maxiter: int | None = None,
        interpolation: str | None = None,
        springconst: float | None = None,
        use_ts_guess: bool | None = None,
    ) -> tuple[Output, Path, Atoms, int, float | None, Path]:
        """
        Run NEB-TS/CI and select a TS guess from the converged NEB output.
        """
        if len(reactant) != len(product):
            raise OpiJobError("Reactant/product atom counts differ; NEB requires matching atoms")
        workdir = ensure_dir(workdir)
        end_xyz = workdir / f"{job_name}_product.xyz"
        write_xyz(product, end_xyz)

        block_kwargs: dict = {"neb_end_xyzfile": end_xyz.name}
        nimages = self.settings.cineb_nimages if nimages is None else nimages
        if nimages is not None:
            block_kwargs["nimages"] = nimages
        maxiter = self.settings.cineb_maxiter if maxiter is None else maxiter
        if maxiter is not None:
            block_kwargs["maxiter"] = maxiter
        interpolation = self.settings.cineb_interpolation if interpolation is None else interpolation
        if interpolation:
            block_kwargs["interpolation"] = interpolation
        springconst = self.settings.cineb_springconst if springconst is None else springconst
        if springconst is not None:
            block_kwargs["springconst"] = springconst
        use_ts_guess = self.settings.cineb_use_ts_guess if use_ts_guess is None else use_ts_guess
        if ts_guess is not None and use_ts_guess:
            ts_input = workdir / f"{job_name}_ts_input.xyz"
            write_xyz(ts_guess, ts_input)
            block_kwargs["ts"] = ts_input.name

        neb_block = BlockNeb(**block_kwargs)
        output = self._run_calculation(
            reactant,
            job_name,
            workdir,
            self.settings.cineb_keywords,
            blocks=[neb_block],
            charge=charge,
            mult=mult,
        )

        if output.results_properties is None:
            try:
                output.parse(do_create_property_json=True, do_create_gbw_json=False)
            except Exception:  # noqa: BLE001
                pass

        neb_mode = _neb_mode_from_keywords(self.settings.cineb_keywords)
        if neb_mode == "neb-ts":
            converged = workdir / f"{job_name}_NEB-TS_converged.xyz"
            label = "NEB-TS"
        else:
            converged = workdir / f"{job_name}_NEB-CI_converged.xyz"
            label = "NEB-CI"

        if not converged.exists():
            raise OpiJobError(f"{label} converged xyz not found: {converged.name}")

        selected_atoms = read_last_xyz_frame(converged)
        ts_guess_path = workdir / f"{job_name}_ts_guess.xyz"
        write_xyz(selected_atoms, ts_guess_path)
        return output, ts_guess_path, selected_atoms, -1, None, converged

    def optimize_ts(
        self,
        atoms: Atoms,
        job_name: str,
        workdir: Path,
        charge: int = 0,
        mult: int = 1,
    ) -> tuple[Output, Path, Atoms]:
        """
        Run TS optimization with restart strategies when SCF/geometry convergence fails.

        Returns the Output object, xyz path, and the final geometry.
        """
        attempt = 0
        last_outfile: Path | None = None
        while attempt <= self.settings.ts_max_restarts:
            suffix = "" if attempt == 0 else f"_retry{attempt}"
            job_label = job_name + suffix
            recalc = self.settings.ts_recalc_hess if attempt > 0 else None
            block_strings = [_geom_block_text(self.settings.geom_maxiter, restart=False, recalc_hess=recalc)]
            keywords = self.settings.ts_keywords

            output = self._run_calculation(
                atoms,
                job_label,
                workdir,
                keywords,
                block_strings=block_strings,
                charge=charge,
                mult=mult,
            )
            try:
                last_outfile = Path(output.get_outfile())
            except Exception:
                last_outfile = Path(workdir) / f"{job_label}.out"

            converged = False
            try:
                converged = bool(output.geometry_optimization_converged())
            except Exception:
                converged = False
            if not converged and last_outfile.exists():
                converged = _geometry_converged_from_outfile(last_outfile)

            if converged:
                xyz_path, final_atoms = self._final_atoms_from_output(output, workdir, job_label)
                return output, xyz_path, final_atoms

            attempt += 1

        raise OpiJobError(
            f"TS optimization failed after {self.settings.ts_max_restarts} restarts; see {Path(last_outfile).name if last_outfile else 'output'}"
        )

    def run_irc(
        self,
        atoms: Atoms,
        job_name: str,
        workdir: Path,
        charge: int = 0,
        mult: int = 1,
        irc_maxiter: int | None = None,
        recalc_hess: int | None = None,
        direction: str = "both",
    ) -> tuple[Output, Path | None, Path | None, Atoms | None, Atoms | None]:
        """
        Submit IRC calculation and return available endpoints.

        direction: "both" (default), "forward", or "backward".
        """
        maxiter = irc_maxiter if irc_maxiter is not None else self.settings.irc_maxiter
        recalc = self.settings.irc_recalc_hess if recalc_hess is None else recalc_hess
        if direction not in {"both", "forward", "backward"}:
            raise ValueError(f"Unsupported IRC direction: {direction}")
        irc_block = BlockIrc(direction=direction, maxiter=maxiter)
        geom_block = _geom_block_text(self.settings.geom_maxiter, restart=False, recalc_hess=recalc)
        output = self._run_calculation(
            atoms,
            job_name,
            workdir,
            self.settings.irc_keywords,
            block_strings=[geom_block],
            blocks=[irc_block],
            charge=charge,
            mult=mult,
        )
        back_xyz = workdir / f"{job_name}_IRC_B.xyz"
        forward_xyz = workdir / f"{job_name}_IRC_F.xyz"

        back_atoms = forward_atoms = None
        if direction in {"both", "backward"}:
            if not back_xyz.exists():
                raise OpiJobError(f"IRC backward xyz not found: {back_xyz.name}")
            back_atoms = read_last_xyz_frame(back_xyz)
        if direction in {"both", "forward"}:
            if not forward_xyz.exists():
                raise OpiJobError(f"IRC forward xyz not found: {forward_xyz.name}")
            forward_atoms = read_last_xyz_frame(forward_xyz)
        return output, (back_xyz if back_atoms is not None else None), (forward_xyz if forward_atoms is not None else None), back_atoms, forward_atoms

    def optimize_minimum(
        self,
        atoms: Atoms,
        job_name: str,
        workdir: Path,
        charge: int = 0,
        mult: int = 1,
    ) -> tuple[Output, Path, Atoms]:
        """Standard Opt+Freq at same theory level as TS."""
        block_strings = [_geom_block_text(self.settings.geom_maxiter, restart=False, recalc_hess=None)]
        output = self._run_calculation(
            atoms,
            job_name,
            workdir,
            self.settings.opt_keywords,
            block_strings=block_strings,
            charge=charge,
            mult=mult,
        )
        converged = False
        try:
            converged = bool(output.geometry_optimization_converged())
        except Exception:
            converged = False
        outfile: Path | None = None
        if not converged:
            try:
                outfile = Path(output.get_outfile())
            except Exception:
                outfile = workdir / f"{job_name}.out"
            if outfile.exists():
                converged = _geometry_converged_from_outfile(outfile)

        if not converged:
            out_name = outfile.name if outfile else f"{job_name}.out"
            raise OpiJobError(f"Minima optimization failed; inspect {out_name}")
        xyz_path, final_atoms = self._final_atoms_from_output(output, workdir, job_name)
        return output, xyz_path, final_atoms


# Backwards compatibility for existing imports
OrcaJobError = OpiJobError
OrcaRunner = OpiRunner
