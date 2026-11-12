"""FlockDOCK orchestration entry point (1.0.9)."""
from __future__ import annotations
import csv
import json
import math
import os
import re
import shutil
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

APP_VERSION = "1.0.9"
PROJECT = Path(os.environ.get("UNIDOCK_PROJECT", os.environ.get("GENDOCK_PROJECT", os.environ.get("ONEDOCK_PROJECT", os.getcwd())))).expanduser().resolve()
SCRIPTS = Path(__file__).resolve().parent
CONFIG = PROJECT / "flockdock_project.json"


def fail(message):
    raise RuntimeError(message)


if not CONFIG.is_file():
    fail(f"Project configuration not found: {CONFIG}. Start the run from the FlockDOCK GUI.")
try:
    cfg = json.loads(CONFIG.read_text(encoding="utf-8"))
except Exception as exc:
    fail(f"Could not read project configuration {CONFIG}: {exc}")

proteins = cfg.get("proteins")
if not isinstance(proteins, list) or not proteins:
    fail("No target proteins are configured. Add at least one target in FlockDOCK.")
proteins = list(dict.fromkeys(str(x).strip().upper() for x in proteins if str(x).strip()))
if not proteins:
    fail("Target list is empty after validation.")
for protein in proteins:
    if not re.fullmatch(r"[A-Z0-9][A-Z0-9_.-]{0,31}", protein):
        fail(f"Unsafe or invalid target identifier: {protein!r}. Use a gene symbol/accession without path characters.")

ligand = Path(str(cfg.get("ligand_file", ""))).expanduser()
if not ligand.is_absolute():
    ligand = (PROJECT / ligand).resolve()
if not ligand.is_file():
    fail(f"Compound library not found: {ligand}")
settings = cfg.get("settings") or {}
search_mode = str(settings.get("search_mode", "fast")).lower()
if search_mode not in {"fast", "balance", "detail", "balanced", "slow"}:
    fail(f"Unsupported search mode: {search_mode}")
try:
    exhaustiveness = int(settings.get("exhaustiveness", 8))
    num_modes = int(settings.get("num_modes", 9))
except (TypeError, ValueError):
    fail("Exhaustiveness and binding modes must be integers.")
if not 1 <= exhaustiveness <= 128:
    fail("Exhaustiveness must be between 1 and 128.")
if not 1 <= num_modes <= 20:
    fail("Binding modes must be between 1 and 20.")

try:
    df = pd.read_csv(ligand, dtype=str, keep_default_na=False, low_memory=False)
except Exception as exc:
    fail(f"Could not read compound CSV: {exc}")
df.columns = [str(c).strip() for c in df.columns]
lower = {c.lower(): c for c in df.columns}

def resolve(aliases):
    for alias in aliases:
        if alias.lower() in lower:
            return lower[alias.lower()]
    for alias in aliases:
        for col in df.columns:
            if alias.lower() in col.lower():
                return col
    return None

id_col = resolve(["LOTUS_ID", "Ligand_LOTUS_ID", "ligand_id", "compound_id", "molecule_id", "id"])
smiles_col = resolve(["Canonical_SMILES", "Ligand_SMILES", "SMILES", "smiles", "canonical smiles"])
if id_col is None or smiles_col is None:
    fail(f"CSV must contain compound ID and SMILES columns. Found: {list(df.columns)}")
valid = df[id_col].astype(str).str.strip().ne("") & df[smiles_col].astype(str).str.strip().ne("")
if not valid.any():
    fail("The compound CSV contains no rows with both a compound ID and SMILES.")

input_dir = PROJECT / "01_csv"
input_dir.mkdir(parents=True, exist_ok=True)
master = input_dir / "Docking_Library.csv"
if ligand.resolve() != master.resolve():
    import shutil
    shutil.copy2(ligand, master)

run_started = datetime.now(timezone.utc)
run_id = run_started.strftime("%Y%m%dT%H%M%SZ")
start_clock = time.monotonic()
print(f"FlockDOCK {APP_VERSION} | project config version: {cfg.get('version', 'legacy/unspecified')} | Python: {sys.executable}", flush=True)
print(f"Project: {PROJECT}", flush=True)
print(f"Run ID: {run_id}", flush=True)
print(f"Compound library: {int(valid.sum()):,} usable rows / {len(df):,} total", flush=True)
print(f"Targets: {', '.join(proteins)}", flush=True)
print("Design: all valid compounds are screened against every selected target.", flush=True)

# The project manifest is a compact reproducibility record for this run.
manifest = {
    "flockdock_version": APP_VERSION,
    "project_config_version": cfg.get("version"),
    "run_id": run_id,
    "started_utc": run_started.isoformat(),
    "project_dir": str(PROJECT),
    "compound_csv": str(ligand),
    "compound_rows": int(len(df)),
    "usable_compound_rows": int(valid.sum()),
    "targets": proteins,
    "settings": {"search_mode": search_mode, "exhaustiveness": exhaustiveness, "num_modes": num_modes,
                 "drugs_only": bool(settings.get("drugs_only", False))},
}
logs = PROJECT / "08_logs"
logs.mkdir(parents=True, exist_ok=True)
manifest_path = logs / f"run_manifest_{run_id}.json"
manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")

# Version and run identity are passed to all stages through environment variables.
env = os.environ.copy()
env.update({
    "FLOCKDOCK_PROJECT": str(PROJECT), "UNIDOCK_PROJECT": str(PROJECT),
    "GENDOCK_PROJECT": str(PROJECT), "ONEDOCK_PROJECT": str(PROJECT),
    "FLOCKDOCK_SEARCH_MODE": search_mode, "GENDOCK_SEARCH_MODE": search_mode,
    "ONEDOCK_SEARCH_MODE": search_mode,
    "GENDOCK_EXHAUSTIVENESS": str(exhaustiveness), "GENDOCK_NUM_MODES": str(num_modes),
    "ONEDOCK_NUM_MODES": str(num_modes), "GENDOCK_VERSION": APP_VERSION, "GENDOCK_RUN_ID": run_id,
})


def run(args, label, allow_failure=False, allowed_codes=()):
    print("STAGE: " + label, flush=True)
    print("$ " + " ".join(map(str, args)), flush=True)
    stage_start = time.monotonic()
    try:
        result = subprocess.run([sys.executable, "-u", *map(str, args)], cwd=str(PROJECT), env=env, check=False)
    except OSError as exc:
        if allow_failure:
            print(f"WARNING: could not start {label}: {exc}", flush=True)
            return 127
        raise RuntimeError(f"Could not start {label}: {exc}") from exc
    duration = time.monotonic() - stage_start
    print(f"Stage elapsed: {duration:.1f} s | exit code: {result.returncode}", flush=True)
    if result.returncode != 0 and result.returncode not in set(allowed_codes):
        if allow_failure:
            print(f"WARNING: {label} returned exit code {result.returncode}; trying the next validated pocket source.", flush=True)
            return result.returncode
        raise RuntimeError(f"{label} failed (exit code {result.returncode}). Check 08_logs for the detailed stage log.")
    return 0

# Shared ligand preparation is a prerequisite for all targets.
run([SCRIPTS / "01_generate_sdf.py"], "Ligand 3D preparation")
run([SCRIPTS / "02_generate_pdbqt.py"], "Ligand QC and PDBQT preparation")

uniprot_ids = cfg.get("uniprot_ids") or {}
target_status = []
for protein in proteins:
    target_start = time.monotonic()
    current_stage = "Target initialization"
    quality_warnings = []
    try:
        print("\n" + "=" * 72 + "\nTARGET: " + protein + "\n" + "=" * 72, flush=True)
        # Preserve the previous run's result table, but keep the live result
        # directory clean so stale poses can never be interpreted as new output.
        previous_results = PROJECT / "07_results" / protein
        if previous_results.is_dir() and any(previous_results.iterdir()):
            archive_root = PROJECT / "07_results_archive"
            archive_root.mkdir(parents=True, exist_ok=True)
            archive_path = archive_root / f"{protein}_{run_id}"
            suffix = 1
            while archive_path.exists():
                suffix += 1
                archive_path = archive_root / f"{protein}_{run_id}_{suffix}"
            shutil.move(str(previous_results), str(archive_path))
            print(f"Archived previous results: {archive_path}", flush=True)
        accession = str(uniprot_ids.get(protein, "")).strip().upper()
        download_args = [SCRIPTS / "04_download_protein.py", protein]
        if accession:
            download_args.append(accession)
            print(f"Using supplied UniProt accession for {protein}: {accession}", flush=True)
        current_stage = "Protein structure download"
        run(download_args, "Protein structure: " + protein)
        cleanup = [SCRIPTS / "05_extract_protein.py", protein]
        if settings.get("drugs_only"):
            cleanup.append("--drugs-only")
        current_stage = "Protein cleanup"
        run(cleanup, "Protein cleanup: " + protein)
        current_stage = "Receptor preparation"
        run([SCRIPTS / "06_prepare_receptor.py", protein], "Receptor preparation: " + protein)
        receptor_report = PROJECT / "05_protein_pdbqt" / f"{protein}_report.txt"
        if receptor_report.is_file():
            receptor_text = receptor_report.read_text(encoding="utf-8", errors="replace")
            if "Unsupported Ni/Cu atom(s) were removed" in receptor_text:
                quality_warnings.append("Ni/Cu was removed from temporary Meeko input; metal-dependent binding may be misrepresented")
            if "Selenium/SEC/MSE was normalized" in receptor_text:
                quality_warnings.append("Selenium/SEC/MSE was normalized for Meeko; selenoprotein active-site chemistry may be altered")

        selection_report = PROJECT / "04_proteins" / protein / f"{protein}_selection.txt"
        if selection_report.is_file():
            selection_text = selection_report.read_text(encoding="utf-8", errors="replace")
            structure_line = next((line for line in selection_text.splitlines() if line.lower().startswith("structure")), "")
            if "_AF_" in structure_line or "alphafold" in structure_line.lower():
                quality_warnings.append("AlphaFold-predicted structure used; local confidence and binding-site plausibility require review")

        pocket_cfg = (cfg.get("pocket_configs") or {}).get(protein, {"mode": "automatic"}) or {"mode": "automatic"}
        mode = str(pocket_cfg.get("mode", "automatic")).lower()
        if mode == "automatic":
            # Important: P2Rank is not a prerequisite for a pocket. Step 08 first
            # checks experimental co-crystal ligands and step-05 ligand evidence,
            # then falls back to P2Rank. This also allows salvage when P2Rank
            # legitimately returns zero pockets or Java is unavailable.
            current_stage = "Pocket prediction"
            p2rank_code = run([SCRIPTS / "07_predict_pocket.py", protein],
                              "Pocket prediction (optional fallback): " + protein, allow_failure=True)
            current_stage = "Pocket selection"
            pocket_args = [SCRIPTS / "08_extract_pocket.py", protein]
            if settings.get("drugs_only"):
                pocket_args.append("--drugs-only")
        elif mode == "coordinates":
            center, size = pocket_cfg.get("center", []), pocket_cfg.get("size", [])
            if len(center) != 3 or len(size) != 3:
                raise RuntimeError(f"Invalid coordinate pocket for {protein}: center and size need 3 values each.")
            vals = [float(v) for v in [*center, *size]]
            if not all(math.isfinite(v) for v in vals) or any(v <= 0 or v > 100 for v in vals[3:]):
                raise RuntimeError(f"Invalid/non-finite pocket coordinates or box size for {protein}.")
            pocket_args = [SCRIPTS / "08_extract_pocket.py", protein, "--custom-center", *map(str, vals[:3]),
                           "--custom-size", *map(str, vals[3:])]
            current_stage = "User-defined coordinate pocket"
        elif mode == "residues":
            residues = str(pocket_cfg.get("residues", "")).strip()
            chain = str(pocket_cfg.get("chain", "")).strip()
            pad = float(pocket_cfg.get("pad", 5.0))
            if not residues or not chain:
                raise RuntimeError(f"Custom residue pocket for {protein} requires a chain and residue list.")
            if not math.isfinite(pad) or not 0 <= pad <= 30:
                raise RuntimeError(f"Pocket padding for {protein} must be between 0 and 30 Å.")
            pocket_args = [SCRIPTS / "08_extract_pocket.py", protein, "--custom-residues", residues,
                           "--custom-chain", chain, "--custom-pad", str(pad)]
            current_stage = "User-defined residue pocket"
        else:
            raise RuntimeError(f"Unknown pocket mode for {protein}: {mode}")
        run(pocket_args, "Pocket selection: " + protein)
        pocket_table = PROJECT / "06_grid" / protein / "Pocket_Config.csv"
        if pocket_table.is_file():
            try:
                selected_pocket = pd.read_csv(pocket_table)
                if not selected_pocket.empty:
                    pocket_source = str(selected_pocket.iloc[0].get("Pocket_Source", ""))
                    probability = pd.to_numeric(pd.Series([selected_pocket.iloc[0].get("Probability", float("nan"))]), errors="coerce").iloc[0]
                    if pocket_source.lower() == "p2rank" and pd.notna(probability) and float(probability) < 0.30:
                        quality_warnings.append(f"Low-confidence P2Rank pocket (probability {float(probability):.3f})")
            except Exception as exc:
                quality_warnings.append(f"Pocket QC could not be fully read: {exc}")
        current_stage = "GPU docking"
        run([SCRIPTS / "09_run_unidock.py", protein], "GPU docking: " + protein)
        current_stage = "Result collection and QC"
        collection_code = run([SCRIPTS / "10_collect_results.py", protein], "Collect results: " + protein, allowed_codes={2})
        if collection_code == 2:
            quality_warnings.append("Docking coverage incomplete or one or more pose scores failed parsing; inspect result QC files")
        current_stage = "Result annotation"
        run([SCRIPTS / "11_add_compound_names.py", protein], "Annotate results: " + protein)
        elapsed = time.monotonic() - target_start
        target_state = "COMPLETED_WITH_WARNINGS" if quality_warnings else "COMPLETED"
        target_reason = "; ".join(quality_warnings)
        target_status.append({"Protein": protein, "Status": target_state, "Failed_Stage": "", "Reason": target_reason,
                              "Duration_s": round(elapsed, 2)})
        print(f"TARGET {target_state}: {protein} ({elapsed / 60:.1f} min)", flush=True)
    except Exception as exc:
        elapsed = time.monotonic() - target_start
        reason = f"{current_stage}: {type(exc).__name__}: {exc}"
        target_status.append({"Protein": protein, "Status": "SKIPPED", "Failed_Stage": current_stage,
                              "Reason": reason, "Duration_s": round(elapsed, 2)})
        print(f"\nTARGET SKIPPED: {protein}", flush=True)
        print(f"Reason: {reason}", flush=True)
        print("The pipeline will continue with the next target.", flush=True)

status_path = logs / "target_status.csv"
status_df = pd.DataFrame(target_status, columns=["Protein", "Status", "Failed_Stage", "Reason", "Duration_s"])
status_df.to_csv(status_path, index=False)
completed = int(status_df["Status"].astype(str).str.startswith("COMPLETED").sum()) if not status_df.empty else 0
skipped = len(target_status) - completed
elapsed_total = time.monotonic() - start_clock
manifest.update({"finished_utc": datetime.now(timezone.utc).isoformat(), "duration_s": round(elapsed_total, 2),
                 "completed_targets": completed, "skipped_targets": skipped,
                 "target_status": target_status, "status_csv": str(status_path)})
summary_path = logs / f"run_summary_{run_id}.json"
summary_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")

print("\nTARGET SUMMARY", flush=True)
print(f"Completed: {completed}/{len(target_status)}", flush=True)
print(f"Skipped  : {skipped}/{len(target_status)}", flush=True)
print(f"Elapsed  : {elapsed_total / 60:.1f} min", flush=True)
print(f"Status file: {status_path}", flush=True)
print(f"Run summary: {summary_path}", flush=True)
if completed == 0:
    fail("No targets completed successfully. Review target_status.csv and individual stage logs.")
print("STAGE: Complete", flush=True)
