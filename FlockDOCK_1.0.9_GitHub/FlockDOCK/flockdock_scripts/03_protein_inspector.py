from Bio.PDB import PDBParser
from collections import Counter
from pathlib import Path
import sys

# ============================================================
# Check input
# ============================================================

if len(sys.argv) != 2:

    print("\nUsage:")
    print("python 03_protein_inspector.py ALOX12\n")
    sys.exit()

protein = sys.argv[1]

# ============================================================
# Project Paths
# ============================================================

PROJECT = Path(os.environ.get("UNIDOCK_PROJECT", os.environ.get("ONEDOCK_PROJECT", os.getcwd())))

PDB_FILE = PROJECT / "04_proteins" / protein / f"{protein}.pdb"

if not PDB_FILE.exists():

    print(f"\nPDB not found:\n{PDB_FILE}")
    sys.exit()

# ============================================================
# Load Structure
# ============================================================

parser = PDBParser(QUIET=True)

structure = parser.get_structure(protein, str(PDB_FILE))

# ============================================================
# Constants
# ============================================================

METALS = {

    "FE", "ZN", "MG", "MN", "CA",
    "CU", "CO", "NI", "NA", "K",
    "CD", "HG", "SR", "CS",
    "RB", "BA", "AL", "YB"

}

COMMON_SOLVENTS = {

    "HOH",
    "SO4",
    "PO4",
    "GOL",
    "PEG",
    "EDO",
    "ACT",
    "FMT",
    "EOH",
    "DMS"

}

# ============================================================
# Statistics
# ============================================================

chains = set()

protein_residues = 0
protein_atoms = 0

waters = 0

het_counter = Counter()
metal_counter = Counter()
ligand_counter = Counter()

# ============================================================
# Parse Structure
# ============================================================

for model in structure:

    for chain in model:

        chains.add(chain.id)

        for residue in chain:

            hetflag = residue.id[0]
            resname = residue.get_resname().strip()

            # ------------------------------------------------
            # Protein
            # ------------------------------------------------

            if hetflag == " ":

                protein_residues += 1
                protein_atoms += len(residue)

            # ------------------------------------------------
            # Water
            # ------------------------------------------------

            elif resname == "HOH":

                waters += 1

            # ------------------------------------------------
            # HETATM
            # ------------------------------------------------

            else:

                het_counter[resname] += 1

                if resname in METALS:

                    metal_counter[resname] += 1

                elif resname not in COMMON_SOLVENTS:

                    ligand_counter[resname] += 1

# ============================================================
# Report
# ============================================================

print("=" * 70)
print(f"Protein Inspection Report : {protein}")
print("=" * 70)

print(f"PDB File            : {PDB_FILE.name}")
print(f"Chains              : {', '.join(sorted(chains))}")
print(f"Protein Residues    : {protein_residues}")
print(f"Protein Atoms       : {protein_atoms}")
print(f"Water Molecules     : {waters}")
print(f"Total HETATM        : {sum(het_counter.values())}")

print()

# ============================================================
# Metal Ions
# ============================================================

print("Metal Ions")
print("-" * 35)

if metal_counter:

    for metal in sorted(metal_counter):

        print(f"{metal:<10}{metal_counter[metal]}")

else:

    print("None")

print()

# ============================================================
# Likely Ligands
# ============================================================

print("Likely Ligands")
print("-" * 35)

if ligand_counter:

    for ligand in sorted(ligand_counter):

        print(f"{ligand:<10}{ligand_counter[ligand]}")

else:

    print("None")

print()

# ============================================================
# All HETATM Residues
# ============================================================

print("All HETATM Molecules")
print("-" * 35)

if het_counter:

    for het in sorted(het_counter):

        print(f"{het:<10}{het_counter[het]}")

else:

    print("None")

print("=" * 70)

print("\nNOTE:")
print("Protein completeness (missing residues/atoms) should be")
print("checked with PDBFixer before Molecular Dynamics if required.")