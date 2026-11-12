"""Validate ligand SDFs and prepare PDBQT files with hash-based cache invalidation."""
from __future__ import annotations
import hashlib
import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

PROJECT_DIR = Path(os.environ.get("UNIDOCK_PROJECT", os.environ.get("GENDOCK_PROJECT", os.environ.get("ONEDOCK_PROJECT", os.getcwd())))).resolve()
SDF_DIR, PDBQT_DIR, LOG_DIR = PROJECT_DIR / "02_sdf", PROJECT_DIR / "03_pdbqt", PROJECT_DIR / "08_logs"
PDBQT_DIR.mkdir(parents=True, exist_ok=True)
LOG_DIR.mkdir(parents=True, exist_ok=True)
if not SDF_DIR.is_dir():
    raise FileNotFoundError(f"SDF directory not found: {SDF_DIR}")

SUPPORTED_ELEMENTS = {"H", "C", "N", "O", "F", "P", "S", "Cl", "Br", "I"}
SUPPORTED_PDBQT_TYPES = {"H", "HD", "C", "A", "N", "NA", "OA", "SA", "S", "P", "F", "Cl", "Br", "I"}
MAX_ROTATABLE_BONDS = 32
MEEKO_TIMEOUT_SECONDS = int(os.environ.get("GENDOCK_MEEKO_TIMEOUT", "180"))


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def validate_pdbqt(path: Path):
    atom_count = 0
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith(("ATOM  ", "HETATM")):
                atom_count += 1
                fields = line.split()
                if len(fields) < 2:
                    return False, "Malformed atom record"
                atom_type = fields[-1]
                if atom_type not in SUPPORTED_PDBQT_TYPES:
                    return False, f"Unsupported AutoDock atom type: {atom_type}"
                try:
                    xyz = [float(line[a:b]) for a, b in ((30, 38), (38, 46), (46, 54))]
                    if not all(map(lambda x: abs(x) < 1e6, xyz)):
                        return False, "Non-finite or unreasonable atom coordinates"
                except Exception:
                    return False, "Invalid atom coordinates"
    if atom_count == 0:
        return False, "No ATOM/HETATM records in PDBQT"
    return True, ""


def read_sdf(path: Path):
    supplier = Chem.SDMolSupplier(str(path), removeHs=False)
    if len(supplier) < 1:
        return None
    return supplier[0]

sdf_files = sorted(SDF_DIR.glob("*.sdf"))
print(f"FlockDOCK 1.0.9 | SDF files: {len(sdf_files):,}", flush=True)
if not sdf_files:
    raise RuntimeError("No SDF files found. Run 01_generate_sdf.py first.")

success = skipped = failed = 0
failed_rows, rejected_rows = [], []
active_stems = {p.stem for p in sdf_files}

# Remove obsolete PDBQTs that no longer have an input SDF.
for stale in PDBQT_DIR.glob("*.pdbqt"):
    if stale.stem not in active_stems:
        stale.unlink(missing_ok=True)
        Path(str(stale) + ".sha256").unlink(missing_ok=True)

for idx, sdf in enumerate(sdf_files, start=1):
    pdbqt = PDBQT_DIR / f"{sdf.stem}.pdbqt"
    sidecar = Path(str(pdbqt) + ".sha256")
    try:
        source_hash = sha256_file(sdf)
        # Cached conversion is accepted only when its exact input SDF hash matches and the output validates.
        if pdbqt.is_file() and sidecar.is_file() and sidecar.read_text(encoding="utf-8").strip() == source_hash:
            valid, reason = validate_pdbqt(pdbqt)
            if valid:
                skipped += 1
                continue
        pdbqt.unlink(missing_ok=True)
        sidecar.unlink(missing_ok=True)
        mol = read_sdf(sdf)
        if mol is None:
            rejected_rows.append({"Ligand": sdf.stem, "Reason": "RDKit could not read SDF or no molecule record"})
            failed += 1
            continue
        unsupported = sorted({atom.GetSymbol() for atom in mol.GetAtoms() if atom.GetSymbol() not in SUPPORTED_ELEMENTS})
        if unsupported:
            rejected_rows.append({"Ligand": sdf.stem, "Reason": f"Unsupported element(s): {', '.join(unsupported)}"})
            failed += 1
            continue
        torsions = rdMolDescriptors.CalcNumRotatableBonds(mol)
        if torsions > MAX_ROTATABLE_BONDS:
            rejected_rows.append({"Ligand": sdf.stem, "Reason": f"{torsions} rotatable bonds exceeds limit {MAX_ROTATABLE_BONDS}"})
            failed += 1
            continue
        fd, tmp_name = tempfile.mkstemp(prefix=f".{sdf.stem}_", suffix=".pdbqt", dir=PDBQT_DIR)
        os.close(fd)
        tmp_path = Path(tmp_name)
        tmp_path.unlink(missing_ok=True)
        cmd = [sys.executable, "-m", "meeko.cli.mk_prepare_ligand", "-v", "-i", str(sdf),
               "-o", str(tmp_path), "--rename_atoms", "--add_index_map"]
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, encoding="utf-8", errors="replace",
                                    timeout=MEEKO_TIMEOUT_SECONDS, check=False)
        except subprocess.TimeoutExpired:
            failed_rows.append({"Ligand": sdf.stem, "Error": f"Meeko timed out after {MEEKO_TIMEOUT_SECONDS} seconds"})
            failed += 1
            tmp_path.unlink(missing_ok=True)
            continue
        if result.returncode != 0 or not tmp_path.is_file():
            failed_rows.append({"Ligand": sdf.stem, "Return_Code": result.returncode,
                                "STDOUT": (result.stdout or "")[-4000:], "STDERR": (result.stderr or "")[-4000:],
                                "Error": "Meeko returned non-zero or did not create PDBQT"})
            failed += 1
            tmp_path.unlink(missing_ok=True)
            continue
        valid, reason = validate_pdbqt(tmp_path)
        if not valid:
            failed_rows.append({"Ligand": sdf.stem, "STDOUT": (result.stdout or "")[-2000:],
                                "STDERR": (result.stderr or "")[-2000:], "Error": reason})
            failed += 1
            tmp_path.unlink(missing_ok=True)
            continue
        os.replace(tmp_path, pdbqt)
        side_tmp = Path(str(sidecar) + ".tmp")
        side_tmp.write_text(source_hash + "\n", encoding="utf-8")
        os.replace(side_tmp, sidecar)
        success += 1
    except Exception as exc:
        failed += 1
        failed_rows.append({"Ligand": sdf.stem, "Error": f"{type(exc).__name__}: {exc}"})
    if idx % 100 == 0 or idx == len(sdf_files):
        print(f"Processed {idx:,}/{len(sdf_files):,} | prepared {success:,} | cached {skipped:,} | failed {failed:,}", flush=True)

# Always rewrite audit files so an old failure list cannot be mistaken for this run's failures.
pd.DataFrame(failed_rows, columns=["Ligand", "Return_Code", "STDOUT", "STDERR", "Error"]).to_csv(LOG_DIR / "failed_pdbqt_generation.csv", index=False)
pd.DataFrame(rejected_rows, columns=["Ligand", "Reason"]).to_csv(LOG_DIR / "rejected_ligands.csv", index=False)
valid_outputs = [p for p in PDBQT_DIR.glob("*.pdbqt") if p.stem in active_stems and validate_pdbqt(p)[0]]
print("\nPDBQT preparation summary", flush=True)
print(f"Prepared : {success:,}", flush=True)
print(f"Cached   : {skipped:,}", flush=True)
print(f"Failed   : {failed:,}", flush=True)
print(f"Usable outputs: {len(valid_outputs):,}/{len(sdf_files):,}", flush=True)
if not valid_outputs:
    raise RuntimeError("No valid ligand PDBQT files were produced. Review failed_pdbqt_generation.csv and rejected_ligands.csv.")
