"""
OPI-backed helpers to run ORCA calculations used in the workflow.
"""

from __future__ import annotations

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

from opi.core import Calculator  # type: ignore  # noqa: E402
from opi.input.arbitrary_string import ArbitraryStringPos  # type: ignore  # noqa: E402
from opi.input.blocks.block_irc import BlockIrc  # type: ignore  # noqa: E402
from opi.input.structures import Structure  # type: ignore  # noqa: E402
from opi.output.core import Output  # type: ignore  # noqa: E402


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


def _geom_block_text(maxiter: int, restart: bool) -> str:
    """Build a %geom block text with optional restart flag."""
    lines = ["%geom", f"  MaxIter {maxiter}"]
    if restart:
        lines.append("  ReStart true")
    lines.append("end")
    return "\n".join(lines)


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

        for kw in _normalize_keywords(keywords):
            calc.input.add_simple_keywords(kw)

        for block_str in block_strings or ():
            calc.input.add_arbitrary_string(block_str, pos=ArbitraryStringPos.TOP)

        for block in blocks or ():
            calc.input.add_blocks(block)

        try:
            calc.write_input()
            ok = calc.run()
            output = calc.get_output()
            # Ensure JSON files exist so we can parse properties like frequencies.
            try:
                output.parse(do_create_property_json=True, do_create_gbw_json=False)
            except FileNotFoundError:
                output.parse(do_create_property_json=True, do_create_gbw_json=True)
        except Exception as exc:  # noqa: BLE001
            raise OpiJobError(f"Failed to run ORCA job {job_name}: {exc}") from exc

        try:
            outfile = output.get_outfile()
        except Exception:  # noqa: BLE001
            outfile = workdir / f"{job_name}.out"

        if not ok or not output.terminated_normally():
            raise OpiJobError(f"ORCA job {job_name} failed; inspect {outfile}")
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

    def _restart_geometry(self, output: Output, workdir: Path, job_label: str, fallback: Atoms) -> Atoms:
        """Pick the best available geometry to seed a restart."""
        structure = output.get_structure()
        if structure is not None:
            return structure_to_atoms(structure)
        xyz_path = Path(workdir) / f"{job_label}.xyz"
        if xyz_path.exists():
            return read_last_xyz_frame(xyz_path)
        return fallback

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
        current_atoms = atoms
        last_outfile: Path | None = None
        while attempt <= self.settings.max_restarts:
            restart = attempt > 0
            suffix = "" if attempt == 0 else f"_retry{attempt}"
            job_label = job_name + suffix
            block_strings = [_geom_block_text(self.settings.geom_maxiter, restart=restart)]
            if restart:
                block_strings.extend(self.settings.ts_restart_blocks)
                keywords = list(self.settings.ts_keywords) + list(self.settings.ts_restart_keywords)
            else:
                keywords = self.settings.ts_keywords

            output = self._run_calculation(
                current_atoms,
                job_label,
                workdir,
                keywords,
                block_strings=block_strings,
                charge=charge,
                mult=mult,
            )
            last_outfile = output.get_outfile()
            if output.geometry_optimization_converged():
                xyz_path, final_atoms = self._final_atoms_from_output(output, workdir, job_label)
                return output, xyz_path, final_atoms

            current_atoms = self._restart_geometry(output, workdir, job_label, fallback=current_atoms)
            attempt += 1

        raise OpiJobError(
            f"TS optimization failed after {self.settings.max_restarts} restarts; see {last_outfile}"
        )

    def run_irc(
        self,
        atoms: Atoms,
        job_name: str,
        workdir: Path,
        charge: int = 0,
        mult: int = 1,
    ) -> tuple[Output, Path, Path, Atoms, Atoms]:
        """
        Submit IRC calculation in both directions and return endpoints.
        """
        irc_block = BlockIrc(direction="both", maxiter=self.settings.irc_maxiter)
        output = self._run_calculation(
            atoms,
            job_name,
            workdir,
            self.settings.irc_keywords,
            blocks=[irc_block],
            charge=charge,
            mult=mult,
        )
        back_xyz = workdir / f"{job_name}_IRC_B.xyz"
        forward_xyz = workdir / f"{job_name}_IRC_F.xyz"
        if not back_xyz.exists() or not forward_xyz.exists():
            raise OpiJobError(f"IRC endpoint xyz files not found: {back_xyz} and/or {forward_xyz}")
        back_atoms = read_last_xyz_frame(back_xyz)
        forward_atoms = read_last_xyz_frame(forward_xyz)
        return output, back_xyz, forward_xyz, back_atoms, forward_atoms

    def optimize_minimum(
        self,
        atoms: Atoms,
        job_name: str,
        workdir: Path,
        charge: int = 0,
        mult: int = 1,
    ) -> tuple[Output, Path, Atoms]:
        """Standard Opt+Freq at same theory level as TS."""
        block_strings = [_geom_block_text(self.settings.geom_maxiter, restart=False)]
        output = self._run_calculation(
            atoms,
            job_name,
            workdir,
            self.settings.opt_keywords,
            block_strings=block_strings,
            charge=charge,
            mult=mult,
        )
        if not output.geometry_optimization_converged():
            raise OpiJobError(f"Minima optimization failed; inspect {output.get_outfile()}")
        xyz_path, final_atoms = self._final_atoms_from_output(output, workdir, job_name)
        return output, xyz_path, final_atoms


# Backwards compatibility for existing imports
OrcaJobError = OpiJobError
OrcaRunner = OpiRunner
