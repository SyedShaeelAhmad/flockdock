"""Collect and quality-check UniDock pose scores for one target."""
from __future__ import annotations
import math
import os
import sys
from pathlib import Path

import pandas as pd

PROJECT_DIR = Path(os.environ.get("UNIDOCK_PROJECT", os.environ.get("GENDOCK_PROJECT", os.environ.get("ONEDOCK_PROJECT", os.getcwd())))).resolve()
RESULT_DIR = PROJECT_DIR / "07_results"
LIGAND_DIR = PROJECT_DIR / "03_pdbqt"

if len(sys.argv) != 2:
    raise SystemExit("Usage: python 10_collect_results.py PROTEIN")
protein = sys.argv[1]
protein_dir = RESULT_DIR / protein
if not protein_dir.is_dir():
    raise SystemExit(f"Result folder not found: {protein_dir}")

expected = {p.stem for p in LIGAND_DIR.glob("*.pdbqt") if p.is_file() and p.stat().st_size > 0}
if not expected:
    raise SystemExit("No valid input ligand PDBQT files were found in 03_pdbqt.")
files = sorted(protein_dir.glob("*_out.pdbqt"))
if not files:
    raise SystemExit(f"No current UniDock pose files found in {protein_dir}")

observed = {p.name[:-len("_out.pdbqt")] for p in files}
missing = sorted(expected - observed)
extra = sorted(observed - expected)
results, failed = [], []

for path in files:
    ligand = path.name[:-len("_out.pdbqt")]
    affinities = []
    try:
        with path.open("r", encoding="utf-8", errors="replace") as handle:
            for line in handle:
                if line.startswith("REMARK VINA RESULT:"):
                    tokens = line.split()
                    try:
                        value = float(tokens[3])
                        if math.isfinite(value):
                            affinities.append(value)
                    except (IndexError, ValueError):
                        continue
        if not affinities:
            raise ValueError("No finite REMARK VINA RESULT score found")
        results.append({"Ligand": ligand, "Best_Affinity": min(affinities), "Pose_Count": len(affinities)})
    except Exception as exc:
        failed.append({"Ligand": ligand, "Reason": f"{type(exc).__name__}: {exc}"})

if not results:
    pd.DataFrame(failed, columns=["Ligand", "Reason"]).to_csv(protein_dir / "Failed_Parse.csv", index=False)
    raise SystemExit("No valid docking scores could be parsed. See Failed_Parse.csv.")

scores = pd.DataFrame(results).drop_duplicates(subset="Ligand").sort_values(
    "Best_Affinity", ascending=True, kind="stable").reset_index(drop=True)
scores.to_csv(protein_dir / "Docking_Scores.csv", index=False)
pd.DataFrame(failed, columns=["Ligand", "Reason"]).to_csv(protein_dir / "Failed_Parse.csv", index=False)
pd.DataFrame({"Ligand": missing}).to_csv(protein_dir / "Missing_Docking.csv", index=False)
pd.DataFrame({"Output": extra}).to_csv(protein_dir / "Unexpected_Docking_Output.csv", index=False)
coverage = len(set(scores["Ligand"]) & expected) / len(expected) * 100.0
pd.DataFrame([{
    "Protein": protein,
    "Expected_Ligands": len(expected),
    "Output_Files": len(files),
    "Valid_Scored_Ligands": len(scores),
    "Missing_Output_Files": len(missing),
    "Parse_Failures": len(failed),
    "Unexpected_Output_Files": len(extra),
    "Scored_Coverage_Percent": round(coverage, 2),
}]).to_csv(protein_dir / "Docking_Coverage.csv", index=False)

print("=" * 70)
print(f"Protein                 : {protein}")
print(f"Input ligands           : {len(expected)}")
print(f"Output pose files       : {len(files)}")
print(f"Valid scored ligands    : {len(scores)}")
print(f"Missing output files    : {len(missing)}")
print(f"Score parse failures    : {len(failed)}")
print(f"Scored coverage         : {coverage:.2f}%")
print(f"Best score              : {scores.iloc[0]['Best_Affinity']:.3f} kcal/mol")
print("=" * 70)
print(f"Saved: {protein_dir / 'Docking_Scores.csv'}")
print(f"QC   : {protein_dir / 'Docking_Coverage.csv'}")
if missing or failed or extra:
    print("WARNING: results are incomplete; inspect Missing_Docking.csv, Failed_Parse.csv and Unexpected_Docking_Output.csv.")
    raise SystemExit(2)
