"""
05_extract_protein.py

Select the correct chain for the TARGET protein, clean it for docking, and
(if the crystal already has a ligand bound in the active site) emit that
ligand's box so downstream docking can use the KNOWN pocket instead of a
blind P2Rank prediction.

Key behaviour:
  1. Chain selection is UniProt-aware. The target's accession is matched
     against the PDB's DBREF records, so for a complex (e.g. GC / vitamin-D
     binding protein co-crystallised with actin) we pick the target chain,
     not whichever chain happens to carry a cofactor.
  2. If a genuine ligand is bound to the target chain, its bounding box is
     written to <GENE>_ligand_box.csv  ->  use this as the docking box.
  3. Crystallisation additives, buffers, sugars and waters are stripped;
     catalytic metals and modified residues are kept as part of the protein.

Usage:
    python 05_extract_protein.py GENE
    python 05_extract_protein.py GENE custom_structure.pdb
"""

from Bio.PDB import PDBParser, PDBIO, Select
from pathlib import Path
from collections import defaultdict
import os
import sys
import numpy as np

# Shared with 04_download_protein.py and 08_extract_pocket.py so that all three
# steps agree on what counts as a pocket-defining ligand.
from ligand_filter import classify, MIN_LIGAND_ATOMS, MAX_LIGAND_ATOMS

try:
    import requests
except ImportError:
    requests = None

# ============================================================
# Input
# ============================================================

DRUGS_ONLY = "--drugs-only" in sys.argv
_argv = [a for a in sys.argv[1:] if a != "--drugs-only"]

if len(_argv) not in [1, 2]:
    print("\nUsage:")
    print("  python 05_extract_protein.py GENE [custom_structure.pdb] "
          "[--drugs-only]\n")
    sys.exit()

PROTEIN = _argv[0]

PROJECT = Path(os.environ.get(
    "UNIDOCK_PROJECT",
    os.environ.get("ONEDOCK_PROJECT", os.getcwd())))
OUT_DIR = PROJECT / "04_proteins" / PROTEIN

# ---- docking-box settings (used when a co-crystal ligand is found) ----
BOX_PADDING = 5.0      # A added on each side of the ligand
BOX_MIN     = 20.0     # A minimum box edge

# ============================================================
# Which PDB should be used?
# ============================================================

if len(_argv) == 2:
    pdb_filename = _argv[1]
else:
    # Prefer downloaded experimental candidates. If none exist, use the
    # AlphaFold model selected by step 04. The downloader also copies the
    # chosen model to <GENE>.pdb as the canonical structure.
    pdb_files = sorted(OUT_DIR.glob(f"{PROTEIN}_rank*.pdb"))
    if not pdb_files:
        pdb_files = sorted(OUT_DIR.glob(f"{PROTEIN}_AF_*.pdb"))
    if not pdb_files and (OUT_DIR / f"{PROTEIN}.pdb").is_file():
        pdb_files = [OUT_DIR / f"{PROTEIN}.pdb"]
    if not pdb_files:
        print(f"\nNo experimental or AlphaFold structures found for {PROTEIN}")
        sys.exit(1)
    pdb_filename = pdb_files[0].name

PDB_FILE = OUT_DIR / pdb_filename

# ============================================================
# Constants
# ============================================================

METALS = {
    "FE", "FE2", "FE3", "ZN", "ZN2", "MG", "MN", "CA", "CU", "CO",
    "NI", "CD", "HG", "SR", "CS", "RB", "BA", "AL", "YB",
}

# Additive / buffer / cryo / glycan / lipid rejection now lives in
# ligand_filter.ADDITIVES, shared with steps 04 and 08. The local REMOVE set
# that used to live here did not contain SCN, CN, AZI, MLI, URE or any lipid,
# and applied no heavy-atom floor, so a 3-atom thiocyanate could be selected as
# the "active-site ligand".

# Modified residues -> treat as part of the protein (keep, never count as ligand).
MODIFIED_RES = {
    "MSE", "CSO", "CSD", "CSX", "CME", "OCS", "SMC", "SEP", "TPO", "PTR",
    "KCX", "LLP", "MLY", "M3L", "ALY", "HYP", "PCA", "SAC", "CGU", "CAS",
    "CAF", "YCM", "NEP", "HIC", "MHO", "FME",
}


def heavy_atom_count(residue):
    return sum(1 for a in residue.get_atoms() if a.element != "H")


def ligand_verdict(residue):
    """
    'drug', 'cofactor' or 'reject' for a HETATM residue.

    Metals and modified residues are handled separately as part of the
    protein, so they are rejected here before the shared filter is consulted.
    """
    resname = residue.get_resname().strip().upper()
    if resname in MODIFIED_RES:
        return "reject", "modified residue, kept as protein"
    if resname in METALS:
        return "reject", "metal ion, kept as protein"
    return classify(resname, heavy_atom_count(residue), DRUGS_ONLY)


def is_real_ligand(resname):
    """Retained for the selectors: component-level check, no size test."""
    verdict, _ = classify(resname, None, DRUGS_ONLY)
    return (verdict != "reject"
            and resname not in METALS
            and resname not in MODIFIED_RES)


# ============================================================
# Target UniProt  (for UniProt-aware chain selection)
# ============================================================

def read_cached_uniprot():
    f = OUT_DIR / f"{PROTEIN}_uniprot.txt"
    if f.exists():
        acc = f.read_text().strip()
        return acc or None
    return None


def resolve_uniprot(gene):
    if requests is None:
        return None
    params = {
        "query": f"gene_exact:{gene} AND organism_id:9606 AND reviewed:true",
        "fields": "accession", "format": "json", "size": 1,
    }
    try:
        r = requests.get("https://rest.uniprot.org/uniprotkb/search", params=params, timeout=30)
        r.raise_for_status()
        results = r.json().get("results", [])
        return results[0].get("primaryAccession") if results else None
    except Exception:
        return None


def parse_chain_uniprot(pdb_path):
    """Map chain_id -> {UniProt accessions} from DBREF / DBREF1+DBREF2 records."""
    mapping = defaultdict(set)
    pending = None
    with open(pdb_path) as fh:
        for line in fh:
            rec = line[:6].strip()
            if rec == "DBREF":
                try:
                    chain = line[12].strip()
                    db = line[26:32].strip().upper()
                    acc = line[33:41].strip()
                except IndexError:
                    continue
                if db in ("UNP", "UNIPROT") and acc:
                    mapping[chain].add(acc)
            elif rec == "DBREF1":
                try:
                    pending = (line[12].strip(), line[26:32].strip().upper())
                except IndexError:
                    pending = None
            elif rec == "DBREF2":
                if pending and pending[1] in ("UNP", "UNIPROT"):
                    acc = line[18:40].strip().split()[0] if line[18:40].strip() else ""
                    if acc:
                        mapping[pending[0]].add(acc)
                pending = None
    return mapping


target_acc = read_cached_uniprot() or resolve_uniprot(PROTEIN)
if target_acc:
    (OUT_DIR / f"{PROTEIN}_uniprot.txt").write_text(target_acc)   # cache for reuse
chain2acc = parse_chain_uniprot(PDB_FILE)

# ============================================================
# Load Structure
# ============================================================

parser = PDBParser(QUIET=True)
structure = parser.get_structure(PDB_FILE.stem, str(PDB_FILE))
model = structure[0]

# ============================================================
# Summarise chains
# ============================================================

chain_summary = defaultdict(lambda: {"protein": 0, "waters": 0, "ligands": [],
                                     "drugs": [], "cofactors": [],
                                     "ignored": [], "metals": []})

for chain in model:
    for residue in chain:
        resname = residue.get_resname().strip()
        hetflag = residue.id[0]
        if hetflag == " " or resname in MODIFIED_RES:
            chain_summary[chain.id]["protein"] += 1
            continue
        if resname in ("HOH", "DOD"):
            chain_summary[chain.id]["waters"] += 1
            continue
        if resname in METALS:
            chain_summary[chain.id]["metals"].append(resname)
            continue
        verdict, reason = ligand_verdict(residue)
        if verdict == "drug":
            chain_summary[chain.id]["drugs"].append(resname)
            chain_summary[chain.id]["ligands"].append(resname)
        elif verdict == "cofactor":
            chain_summary[chain.id]["cofactors"].append(resname)
            chain_summary[chain.id]["ligands"].append(resname)
        else:
            chain_summary[chain.id]["ignored"].append(
                f"{resname} ({heavy_atom_count(residue)} atoms: {reason})")

# ============================================================
# Select the TARGET chain
# ============================================================

protein_chains = [c for c, s in chain_summary.items() if s["protein"] > 0]

matched = [c for c in protein_chains if target_acc and target_acc in chain2acc.get(c, set())]

if matched:
    candidates = matched
    selection_mode = f"UniProt-matched ({target_acc})"
else:
    candidates = protein_chains
    if target_acc:
        selection_mode = f"FALLBACK - no DBREF chain matched {target_acc}; using longest chain"
    else:
        selection_mode = "FALLBACK - target UniProt unknown; using longest chain"

# Among candidates: prefer a chain that carries a real ligand (defines the
# pocket), then the most complete chain.
def chain_key(cid):
    s = chain_summary[cid]
    return (2 if s["drugs"] else (1 if s["cofactors"] else 0), s["protein"])

best_chain = sorted(candidates, key=chain_key, reverse=True)[0]
info = chain_summary[best_chain]

print()
print("=" * 70)
print("Selected Chain")
print("=" * 70)
print(f"Selection mode     : {selection_mode}")
print(f"Chain              : {best_chain}")
print(f"Chain -> UniProt    : {sorted(chain2acc.get(best_chain, [])) or 'n/a'}")
print(f"Protein Residues   : {info['protein']}")
print(f"Waters             : {info['waters']}")
print(f"Metal Ions         : {sorted(set(info['metals']))}")
print(f"Drug-like ligands  : {sorted(set(info['drugs'])) or 'none'}")
print(f"Cofactors/substrate: {sorted(set(info['cofactors'])) or 'none'}")
if info["ignored"]:
    print(f"Ignored components : {sorted(set(info['ignored']))}")
if "FALLBACK" in selection_mode:
    print("WARNING: chain not verified against target UniProt - inspect before trusting results.")

# ============================================================
# Co-crystal ligand pocket  (the KNOWN active site)
# ============================================================

def best_pocket_ligand(chain_obj):
    """
    Pocket-defining ligand on the chain.

    A drug-like ligand always beats a cofactor or substrate, and only then does
    size decide. Picking purely by size, as this did before, hands the box to
    whichever lipid or cryoprotectant happens to be largest: on ATP2A2 the
    phosphatidylcholine PCW is bigger than the catalytic ATP.
    """
    ranked = []
    for residue in chain_obj:
        if residue.id[0] == " ":
            continue
        resname = residue.get_resname().strip()
        if resname in ("HOH", "DOD"):
            continue
        verdict, reason = ligand_verdict(residue)
        if verdict == "reject":
            continue
        ranked.append((0 if verdict == "drug" else 1,
                       -heavy_atom_count(residue), residue, verdict))
    if not ranked:
        return None, None
    ranked.sort(key=lambda t: (t[0], t[1]))
    return ranked[0][2], ranked[0][3]


ligand_res, ligand_verdict_label = best_pocket_ligand(model[best_chain])
ligand_box_file = None

if ligand_res is not None:
    coords = np.array([a.coord for a in ligand_res.get_atoms() if a.element != "H"])
    mn, mx = coords.min(axis=0), coords.max(axis=0)
    center = (mn + mx) / 2.0
    size = np.maximum((mx - mn) + 2 * BOX_PADDING, BOX_MIN)

    lig_name = ligand_res.get_resname().strip()
    ligand_box_file = OUT_DIR / f"{PROTEIN}_ligand_box.csv"
    with open(ligand_box_file, "w") as f:
        f.write("Protein,Source,Ligand,Ligand_Verdict,Heavy_Atoms,"
                "Center_X,Center_Y,Center_Z,Size_X,Size_Y,Size_Z\n")
        f.write(f"{PROTEIN},cocrystal_ligand,{lig_name},"
                f"{ligand_verdict_label},{len(coords)},"
                f"{center[0]:.3f},{center[1]:.3f},{center[2]:.3f},"
                f"{size[0]:.3f},{size[1]:.3f},{size[2]:.3f}\n")

    print()
    print("=" * 70)
    print("CO-CRYSTAL POCKET FOUND  ->  use this box for docking")
    print("=" * 70)
    print(f"Active-site ligand : {lig_name}  ({len(coords)} heavy atoms, "
          f"{ligand_verdict_label})")
    print(f"Box center         : {center[0]:.3f}  {center[1]:.3f}  {center[2]:.3f}")
    print(f"Box size           : {size[0]:.3f}  {size[1]:.3f}  {size[2]:.3f}")
    print(f"Saved              : {ligand_box_file}")
else:
    print()
    print("No pocket-defining ligand on the target chain "
          f"(floor {MIN_LIGAND_ATOMS}-{MAX_LIGAND_ATOMS} heavy atoms) - "
          "step 08 will fall through to P2Rank.")

# ============================================================
# Selection report
# ============================================================

report = OUT_DIR / f"{PROTEIN}_selection.txt"
with open(report, "w") as f:
    f.write(f"Protein : {PROTEIN}\n")
    f.write(f"Structure : {PDB_FILE.name}\n")
    f.write(f"Target UniProt : {target_acc or 'unknown'}\n")
    f.write(f"Selection mode : {selection_mode}\n")
    f.write(f"Selected Chain : {best_chain}\n")
    f.write(f"Chain UniProt : {sorted(chain2acc.get(best_chain, [])) or 'n/a'}\n\n")
    f.write(f"Protein Residues : {info['protein']}\n")
    f.write(f"Waters : {info['waters']}\n\n")
    f.write("Metals\n----------------------\n")
    f.write(("\n".join(sorted(set(info["metals"]))) or "None") + "\n")
    f.write("\nDrug-like ligands\n----------------------\n")
    f.write(("\n".join(sorted(set(info["drugs"]))) or "None") + "\n")
    f.write("\nCofactors / substrates\n----------------------\n")
    f.write(("\n".join(sorted(set(info["cofactors"]))) or "None") + "\n")
    f.write("\nIgnored components\n----------------------\n")
    f.write(("\n".join(sorted(set(info["ignored"]))) or "None") + "\n")
    if ligand_res is not None:
        f.write(f"\nPocket source : co-crystal ligand "
                f"{ligand_res.get_resname().strip()} ({ligand_verdict_label})\n")
        f.write(f"Box file : {ligand_box_file.name}\n")
    else:
        f.write("\nPocket source : P2Rank (no co-crystal ligand)\n")

# ============================================================
# Selectors
# ============================================================

class ProteinSelect(Select):
    def accept_chain(self, chain):
        return chain.id == best_chain

    def accept_residue(self, residue):
        resname = residue.get_resname().strip()
        hetflag = residue.id[0]
        if hetflag == " ":
            return True
        if resname in MODIFIED_RES:      # keep modified residues as protein
            return True
        if resname in METALS:            # keep catalytic metals
            return True
        return False                     # drop waters, additives, and ligands


class LigandSelect(Select):
    """
    Only the selected pocket ligand, not every ligand on the chain.

    Step 08 counts heavy atoms in this file when validating the box, so mixing
    several ligands here would inflate that count.
    """

    def accept_chain(self, chain):
        return chain.id == best_chain

    def accept_residue(self, residue):
        return ligand_res is not None and residue is ligand_res


# ============================================================
# Save protein + ligand
# ============================================================

io = PDBIO()
io.set_structure(structure)

protein_out = OUT_DIR / f"{PROTEIN}_clean.pdb"
io.save(str(protein_out), ProteinSelect())

ligand_out = None
if ligand_res is not None:
    ligand_out = OUT_DIR / f"{PROTEIN}_ligand.pdb"
    io.set_structure(structure)
    io.save(str(ligand_out), LigandSelect())
    print("\nLigand Saved")
    print(ligand_out)
else:
    print("\nNo pocket-defining ligand detected.")

# ============================================================
# Summary
# ============================================================

print()
print("=" * 70)
print("Output Files")
print("=" * 70)
print(f"Clean Protein : {protein_out}")
if ligand_out:
    print(f"Ligand        : {ligand_out}")
if ligand_box_file:
    print(f"Ligand Box    : {ligand_box_file}")
print(f"Report        : {report}")
print("\nDone.")