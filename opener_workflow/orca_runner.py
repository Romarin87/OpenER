"""
Utilities to build and submit ORCA jobs with simple monitoring/restart logic.
"""

from __future__ import annotations

import subprocess
import textwrap
from pathlib import Path
from typing import Optional, Sequence

from ase import Atoms

from .config import OrcaSettings
from .io_utils import atoms_to_orca_input, ensure_dir, write_xyz


class OrcaJobError(RuntimeError):
    """Raised when an ORCA job fails irrecoverably."""


def run_orca_input(
    input_text: str,
    job_name: str,
    workdir: Path,
    executable: str,
) -> Path:
    """Write an ORCA input and execute it, returning the output path."""
    ensure_dir(workdir)
    inp = workdir / f"{job_name}.inp"
    out = workdir / f"{job_name}.out"
    inp.write_text(input_text)
    try:
        with open(out, "w") as fout:
            subprocess.run(
                [executable, inp.name],
                cwd=workdir,
                check=True,
                stdout=fout,
                stderr=subprocess.STDOUT,
            )
    except subprocess.CalledProcessError as exc:
        raise OrcaJobError(f"ORCA job {job_name} failed; inspect {out}") from exc
    return out


def parse_orca_geometries(output_path: Path) -> list[Atoms]:
    """
    Parse every 'CARTESIAN COORDINATES (ANGSTROEM)' block into ASE Atoms.
    """
    lines = output_path.read_text(errors="ignore").splitlines()
    blocks: list[list[tuple[str, tuple[float, float, float]]]] = []
    coords: list[tuple[str, tuple[float, float, float]]] = []
    reading = False
    for line in lines:
        if "CARTESIAN COORDINATES (ANGSTROEM)" in line.upper():
            if coords:
                blocks.append(coords)
            coords = []
            reading = False
            continue
        if "----" in line and not reading:
            reading = True
            continue
        if reading:
            if not line.strip():
                if coords:
                    blocks.append(coords)
                coords = []
                reading = False
                continue
            parts = line.split()
            if len(parts) < 4:
                continue
            symbol = parts[0]
            try:
                x, y, z = map(float, parts[-3:])
            except ValueError:
                continue
            coords.append((symbol, (x, y, z)))
    if coords:
        blocks.append(coords)

    geoms: list[Atoms] = []
    for block in blocks:
        symbols = [c[0] for c in block]
        positions = [c[1] for c in block]
        geoms.append(Atoms(symbols=symbols, positions=positions))
    return geoms


def parse_orca_geometry(output_path: Path) -> Optional[Atoms]:
    """
    Parse the last 'CARTESIAN COORDINATES (ANGSTROEM)' block into an ASE Atoms.
    """
    geoms = parse_orca_geometries(output_path)
    return geoms[-1] if geoms else None


def orca_optimization_converged(output_path: Path) -> tuple[bool, str]:
    """Detect whether ORCA reports convergence and return (ok, message)."""
    text = output_path.read_text(errors="ignore")
    lower = text.lower()
    if "the optimization has converged" in lower or "optimization converged" in lower:
        return True, "Optimization converged"
    if "normal termination" in lower:
        return True, "Normal termination"
    if "scf convergence failure" in lower:
        return False, "SCF did not converge"
    if "maximum number of optimization cycles" in lower or "too many optimization cycles" in lower:
        return False, "Max optimization cycles exceeded"
    if "error" in lower:
        return False, "Error reported by ORCA"
    return False, "Convergence flag not found"


class OrcaRunner:
    """Thin wrapper around ORCA execution."""

    def __init__(self, settings: OrcaSettings):
        self.settings = settings

    def _run_job(
        self,
        atoms: Atoms,
        job_name: str,
        workdir: Path,
        keywords: str,
        blocks: Sequence[str],
        charge: int = 0,
        mult: int = 1,
    ) -> Path:
        text = atoms_to_orca_input(
            atoms=atoms, keywords=keywords, blocks=blocks, charge=charge, mult=mult
        )
        return run_orca_input(
            input_text=text,
            job_name=job_name,
            workdir=workdir,
            executable=self.settings.executable,
        )

    def optimize_ts(
        self,
        atoms: Atoms,
        job_name: str,
        workdir: Path,
        charge: int = 0,
        mult: int = 1,
        max_restarts: int = 2,
    ) -> tuple[Path, Atoms]:
        """
        Run TS optimization with restart strategies when SCF/geometry convergence fails.

        Returns the output path and the final geometry.
        """
        attempt = 0
        current_atoms = atoms
        while attempt <= max_restarts:
            suffix = "" if attempt == 0 else f"_retry{attempt}"
            out = self._run_job(
                current_atoms,
                job_name + suffix,
                workdir,
                keywords=self.settings.ts_keywords,
                blocks=self.settings.ts_blocks if attempt == 0 else self.settings.restart_blocks,
                charge=charge,
                mult=mult,
            )
            ok, msg = orca_optimization_converged(out)
            if ok:
                final_atoms = parse_orca_geometry(out)
                if final_atoms is None:
                    raise OrcaJobError(f"Could not parse final geometry from {out}")
                xyz_path = workdir / f"{job_name}_opt.xyz"
                write_xyz(final_atoms, xyz_path)
                return out, final_atoms
            attempt += 1
            next_geom = parse_orca_geometry(out)
            if next_geom is not None:
                current_atoms = next_geom
            if attempt > max_restarts:
                raise OrcaJobError(
                    f"TS optimization failed after {max_restarts} restarts ({msg}); see {out}"
                )

    def run_irc(
        self,
        atoms: Atoms,
        job_name: str,
        workdir: Path,
        charge: int = 0,
        mult: int = 1,
    ) -> Path:
        """Submit IRC calculation in both directions."""
        blocks = list(self.settings.ts_blocks) + [self.settings.irc_block]
        return self._run_job(
            atoms,
            job_name,
            workdir,
            keywords=self.settings.irc_keywords,
            blocks=blocks,
            charge=charge,
            mult=mult,
        )

    def optimize_minimum(
        self,
        atoms: Atoms,
        job_name: str,
        workdir: Path,
        charge: int = 0,
        mult: int = 1,
    ) -> tuple[Path, Atoms]:
        """Standard Opt+Freq at same theory level as TS."""
        out = self._run_job(
            atoms,
            job_name,
            workdir,
            keywords=self.settings.opt_keywords,
            blocks=self.settings.ts_blocks,
            charge=charge,
            mult=mult,
        )
        ok, _ = orca_optimization_converged(out)
        if not ok:
            raise OrcaJobError(f"Minima optimization failed; inspect {out}")
        final_atoms = parse_orca_geometry(out)
        if final_atoms is None:
            raise OrcaJobError(f"Could not parse final geometry from {out}")
        xyz_path = workdir / f"{job_name}_opt.xyz"
        write_xyz(final_atoms, xyz_path)
        return out, final_atoms


def summarize_job_template(settings: OrcaSettings) -> str:
    """Helper to display chosen ORCA keywords."""
    return textwrap.dedent(
        f"""
        TS keywords: {settings.ts_keywords}
        TS blocks:
        {chr(10).join(settings.ts_blocks)}

        Restart blocks:
        {chr(10).join(settings.restart_blocks)}

        IRC keywords: {settings.irc_keywords}
        IRC block:
        {settings.irc_block}

        Opt keywords: {settings.opt_keywords}
        """
    ).strip()
