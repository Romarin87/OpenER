"""
Canonical SMILES utilities using OpenBabel.
"""

from __future__ import annotations

import io
from typing import Tuple

from ase import Atoms
from ase.io import write
from openbabel import pybel


def atoms_to_smiles(atoms: Atoms, isomeric: bool = True) -> str:
    """Convert an ASE Atoms to a canonical SMILES string."""
    buf = io.StringIO()
    write(buf, atoms, format="xyz")
    xyz = buf.getvalue()
    mol = pybel.readstring("xyz", xyz)
    # pybel write options: 'c' canonical, 'i' include isomerism
    opts = "ci" if isomeric else "c"
    smiles = mol.write("smi", opt=opts).strip()
    return smiles


def compare_endpoints(
    irc_reactant: Atoms,
    irc_product: Atoms,
    opt_reactant: Atoms,
    opt_product: Atoms,
    isomeric: bool = True,
) -> Tuple[bool, str, str, str, str]:
    """
    Return whether SMILES match between IRC endpoints and optimized minima.
    """
    smiles_irc_r = atoms_to_smiles(irc_reactant, isomeric=isomeric)
    smiles_irc_p = atoms_to_smiles(irc_product, isomeric=isomeric)
    smiles_opt_r = atoms_to_smiles(opt_reactant, isomeric=isomeric)
    smiles_opt_p = atoms_to_smiles(opt_product, isomeric=isomeric)
    ok = smiles_irc_r == smiles_opt_r and smiles_irc_p == smiles_opt_p
    return ok, smiles_irc_r, smiles_irc_p, smiles_opt_r, smiles_opt_p
