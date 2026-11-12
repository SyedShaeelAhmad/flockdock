"""Generate deterministic, validated 3D ligand SDF files with safe caching."""
from __future__ import annotations
import hashlib
import json
import os
import re
import tempfile
from pathlib import Path

import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem

PROJECT_DIR = Path(os.environ.get("UNIDOCK_PROJECT", os.environ.get("GENDOCK_PROJECT", os.environ.get("ONEDOCK_PROJECT", os.getcwd())))).resolve()
INPUT_CSV = PROJECT_DIR / "01_csv" / "Docking_Library.csv"
OUTPUT_DIR = PROJECT_DIR / "02_sdf"
LOG_DIR = PROJECT_DIR / "08_logs"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
LOG_DIR.mkdir(parents=True, exist_ok=True)

if not INPUT_CSV.is_file():
    raise FileNotFoundError(f"Compound CSV not found: {INPUT_CSV}")

df = pd.read_csv(INPUT_CSV, dtype=str, keep_default_na=False, low_memory=False)
df.columns = [str(c).strip() for c in df.columns]

def resolve_col(frame, aliases):
    lower = {str(c).lower().strip(): c for c in frame.columns}
    for alias in aliases:
        if alias.lower() in lower:
            return lower[alias.lower()]
    for alias in aliases:
        for col in frame.columns:
            if alias.lower() in str(col).lower():
                return col
    return None

COLS = {
    "ligand_id": resolve_col(df, ["LOTUS_ID", "Ligand_LOTUS_ID", "ligand_id", "compound_id", "molecule_id", "id"]),
    "smiles": resolve_col(df, ["Canonical_SMILES", "Ligand_SMILES", "SMILES", "smiles", "canonical smiles"]),
    "compound_name": resolve_col(df, ["compound_name", "Ligand_name", "Compound_Name", "name"]),
    "gene": resolve_col(df, ["Gene", "Target_Gene", "Original_Target_Gene"]),
    "target_name": resolve_col(df, ["Target_Name", "LOTUS_target_name", "target_name"]),
    "activity_type": resolve_col(df, ["standard_type", "Activity_Type", "activity_type"]),
    "activity_value": resolve_col(df, ["standard_value", "Activity_Value", "activity_value", "pchembl_value"]),
}
for key in ("ligand_id", "smiles"):
    if COLS[key] is None:
        raise KeyError(f"Could not find required {key} column. Available columns: {list(df.columns)}")


def cell(row, key):
    col = COLS[key]
    return str(row.get(col, "")).strip() if col else ""


def safe_filename(name):
    return "".join(c if (c.isalnum() or c in "._-") else "_" for c in str(name)).strip("._")


def normalize_smiles(value):
    value = re.sub(r"\\+([@()\[\]+-])", r"\1", str(value)).strip()
    if len(value) >= 2 and value[0] == value[-1] and value[0] in {"'", '"'}:
        value = value[1:-1].strip()
    return value


def row_hash(row, normalized_smiles):
    payload = {str(k): str(v) for k, v in row.items()}
    payload["__normalized_smiles__"] = normalized_smiles
    return hashlib.sha256(json.dumps(payload, sort_keys=True, ensure_ascii=False).encode("utf-8")).hexdigest()


def valid_cached_sdf(path, expected_hash):
    try:
        supplier = Chem.SDMolSupplier(str(path), removeHs=False)
        mol = supplier[0] if len(supplier) else None
        return bool(mol is not None and mol.HasProp("_FLOCKDOCK_SOURCE_HASH")
                    and mol.GetProp("_FLOCKDOCK_SOURCE_HASH") == expected_hash
                    and mol.GetNumConformers() > 0)
    except Exception:
        return False

# Deduplicate exact compound IDs while retaining first occurrence; reject unsafe collisions.
df[COLS["ligand_id"]] = df[COLS["ligand_id"]].astype(str).str.strip()
df = df[df[COLS["ligand_id"]] != ""].drop_duplicates(subset=COLS["ligand_id"], keep="first").reset_index(drop=True)
name_to_id = {}
for ligand_id in df[COLS["ligand_id"]]:
    safe = safe_filename(ligand_id)
    if not safe:
        raise ValueError(f"Compound ID cannot be converted to a safe filename: {ligand_id!r}")
    if safe in name_to_id and name_to_id[safe] != ligand_id:
        raise ValueError(f"Compound IDs {name_to_id[safe]!r} and {ligand_id!r} map to the same filename {safe!r}.")
    name_to_id[safe] = ligand_id

print(f"FlockDOCK 1.0.9 | unique compound IDs: {len(df):,}", flush=True)
print(f"Input CSV: {INPUT_CSV}", flush=True)

success = skipped = failed = 0
failed_list = []
retained = set()

for idx, row in df.iterrows():
    ligand_id = cell(row, "ligand_id")
    raw_smiles = cell(row, "smiles")
    smiles = normalize_smiles(raw_smiles)
    outfile = OUTPUT_DIR / f"{safe_filename(ligand_id)}.sdf"
    if not smiles:
        failed += 1
        failed_list.append({"LOTUS_ID": ligand_id, "SMILES": raw_smiles, "Error": "Missing SMILES"})
        outfile.unlink(missing_ok=True)
        continue
    source_hash = row_hash(row, smiles)
    if outfile.is_file() and valid_cached_sdf(outfile, source_hash):
        retained.add(outfile.name)
        skipped += 1
        continue
    outfile.unlink(missing_ok=True)
    tmp = None
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError("RDKit rejected normalized SMILES")
        mol = Chem.AddHs(mol)
        params = AllChem.ETKDGv3()
        params.randomSeed = 42
        params.numThreads = 1
        status = AllChem.EmbedMolecule(mol, params)
        if status != 0:
            params.useRandomCoords = True
            status = AllChem.EmbedMolecule(mol, params)
        if status != 0:
            raise ValueError("3D embedding failed (ETKDG and random-coordinate retry)")
        if AllChem.MMFFHasAllMoleculeParams(mol):
            AllChem.MMFFOptimizeMolecule(mol, maxIters=1000)
            force_field = "MMFF"
        elif AllChem.UFFHasAllMoleculeParams(mol):
            AllChem.UFFOptimizeMolecule(mol, maxIters=1000)
            force_field = "UFF"
        else:
            raise ValueError("No complete MMFF or UFF parameters for this molecule")
        conf = mol.GetConformer()
        for atom_idx in range(mol.GetNumAtoms()):
            p = conf.GetAtomPosition(atom_idx)
            if not all(map(lambda x: abs(float(x)) < 1e6, (p.x, p.y, p.z))):
                raise ValueError("Non-finite or unreasonable 3D coordinates generated")
        mol.SetProp("_Name", ligand_id)
        mol.SetProp("LOTUS_ID", ligand_id)
        mol.SetProp("Compound_Name", cell(row, "compound_name"))
        mol.SetProp("Canonical_SMILES", Chem.MolToSmiles(Chem.RemoveHs(mol), isomericSmiles=True))
        mol.SetProp("Original_SMILES", raw_smiles)
        mol.SetProp("Original_Target_Gene", cell(row, "gene"))
        mol.SetProp("Original_Target_Name", cell(row, "target_name"))
        mol.SetProp("Activity_Type", cell(row, "activity_type"))
        mol.SetProp("Activity_Value", cell(row, "activity_value"))
        mol.SetProp("_FLOCKDOCK_SOURCE_HASH", source_hash)
        mol.SetProp("_FLOCKDOCK_FORCE_FIELD", force_field)
        fd, tmp = tempfile.mkstemp(prefix=f".{outfile.stem}_", suffix=".sdf", dir=OUTPUT_DIR)
        os.close(fd)
        writer = Chem.SDWriter(tmp)
        try:
            writer.write(mol)
        finally:
            writer.close()
        # Validate the serialized file before replacing a previous cache.
        if not valid_cached_sdf(Path(tmp), source_hash):
            raise ValueError("Generated SDF failed post-write validation")
        os.replace(tmp, outfile)
        tmp = None
        retained.add(outfile.name)
        success += 1
    except Exception as exc:
        failed += 1
        outfile.unlink(missing_ok=True)
        failed_list.append({"LOTUS_ID": ligand_id, "Compound_Name": cell(row, "compound_name"),
                            "SMILES": raw_smiles, "Normalized_SMILES": smiles, "Error": str(exc)})
    finally:
        if tmp and os.path.exists(tmp):
            try: os.remove(tmp)
            except OSError: pass
    if (idx + 1) % 100 == 0:
        print(f"Processed {idx + 1:,}/{len(df):,} rows | generated {success:,} | cached {skipped:,} | failed {failed:,}", flush=True)

# Remove stale SDFs from libraries that were edited or replaced; never pass obsolete ligands downstream.
removed = 0
for path in OUTPUT_DIR.glob("*.sdf"):
    if path.name not in retained:
        path.unlink(missing_ok=True)
        removed += 1

pd.DataFrame(failed_list, columns=["LOTUS_ID", "Compound_Name", "SMILES", "Normalized_SMILES", "Error"]).to_csv(
    LOG_DIR / "failed_sdf_generation.csv", index=False)
print("\nSDF preparation summary", flush=True)
print(f"Unique IDs : {len(df):,}", flush=True)
print(f"Generated  : {success:,}", flush=True)
print(f"Cached     : {skipped:,}", flush=True)
print(f"Failed     : {failed:,}", flush=True)
print(f"Stale SDFs removed: {removed:,}", flush=True)
if not retained:
    raise RuntimeError("No valid SDF ligands remain after preparation. See 08_logs/failed_sdf_generation.csv.")
