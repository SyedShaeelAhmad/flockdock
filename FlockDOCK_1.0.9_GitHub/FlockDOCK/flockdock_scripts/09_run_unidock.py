"""Stage inputs into WSL2's Linux filesystem and run UniDock GPU docking."""
from __future__ import annotations
import glob
import io
import math
import os
import shlex
import shutil
import subprocess
import sys
import tarfile
import time
import uuid
from pathlib import Path

import pandas as pd

PROJECT_DIR = Path(os.environ.get("UNIDOCK_PROJECT", os.environ.get("GENDOCK_PROJECT", os.environ.get("ONEDOCK_PROJECT", os.getcwd())))).resolve()
LIGAND_DIR = PROJECT_DIR / "03_pdbqt"
RECEPTOR_DIR = PROJECT_DIR / "05_protein_pdbqt"
GRID_DIR = PROJECT_DIR / "06_grid"
RESULT_DIR = PROJECT_DIR / "07_results"
LOG_DIR = PROJECT_DIR / "08_logs"
LOG_DIR.mkdir(parents=True, exist_ok=True)

SEARCH_MODE = os.environ.get("GENDOCK_SEARCH_MODE", os.environ.get("ONEDOCK_SEARCH_MODE", "fast"))
NUM_MODES = int(os.environ.get("GENDOCK_NUM_MODES", os.environ.get("ONEDOCK_NUM_MODES", "9")))
EXHAUSTIVENESS = int(os.environ.get("GENDOCK_EXHAUSTIVENESS", "8"))
DEFAULT_BOX_SIZE = 25.0
WSL_DISTRO = os.environ.get("ONEDOCK_WSL_DISTRO", "Ubuntu")

if len(sys.argv) != 2:
    raise SystemExit("Usage: python 09_run_unidock.py PROTEIN")
protein = sys.argv[1]
receptor = RECEPTOR_DIR / f"{protein}.pdbqt"
config = GRID_DIR / protein / "Pocket_Config.csv"
output_folder = RESULT_DIR / protein

if not receptor.is_file() or receptor.stat().st_size < 500:
    raise SystemExit(f"ERROR: receptor not found or appears empty: {receptor}")
if not config.is_file():
    raise SystemExit(f"ERROR: pocket configuration not found: {config}")
if not 1 <= EXHAUSTIVENESS <= 128 or not 1 <= NUM_MODES <= 20:
    raise SystemExit("ERROR: exhaustiveness must be 1-128 and num_modes must be 1-20.")
if SEARCH_MODE not in {"fast", "balance", "detail", "balanced", "slow"}:
    raise SystemExit(f"ERROR: unsupported UniDock search mode: {SEARCH_MODE}")

try:
    cfg = pd.read_csv(config)
except Exception as exc:
    raise SystemExit(f"ERROR: could not read Pocket_Config.csv: {exc}")
if cfg.empty:
    raise SystemExit("ERROR: Pocket_Config.csv contains no pocket rows.")
for column in ("Center_X", "Center_Y", "Center_Z"):
    if column not in cfg.columns:
        raise SystemExit(f"ERROR: Pocket_Config.csv is missing column: {column}")
try:
    center = [float(cfg.loc[0, c]) for c in ("Center_X", "Center_Y", "Center_Z")]
except (TypeError, ValueError) as exc:
    raise SystemExit(f"ERROR: invalid pocket center: {exc}")
if not all(math.isfinite(v) and abs(v) < 1e5 for v in center):
    raise SystemExit("ERROR: pocket center contains non-finite or unreasonable coordinates.")

def read_size(column):
    if column in cfg.columns:
        try:
            value = float(cfg.loc[0, column])
            if math.isfinite(value) and 0 < value <= 100:
                return value
        except (TypeError, ValueError):
            pass
    return DEFAULT_BOX_SIZE

size = [read_size(c) for c in ("Size_X", "Size_Y", "Size_Z")]
pocket_source = str(cfg.loc[0, "Pocket_Source"]) if "Pocket_Source" in cfg.columns else "unknown"
ligands = sorted(p for p in LIGAND_DIR.glob("*.pdbqt") if p.is_file() and p.stat().st_size > 0)
if not ligands:
    raise SystemExit("ERROR: no non-empty ligand PDBQT files found.")

wsl_exe = shutil.which("wsl.exe") or shutil.which("wsl")
if not wsl_exe:
    raise SystemExit("ERROR: WSL2 is not installed or wsl.exe is not on PATH. Run setup_FlockDOCK.bat.")

def wsl(*args, capture=False, timeout=60, input_data=None):
    return subprocess.run([wsl_exe, "-d", WSL_DISTRO, "--", *map(str, args)],
                          input=input_data, capture_output=capture,
                          text=input_data is None, encoding="utf-8" if input_data is None else None,
                          errors="replace" if input_data is None else None, timeout=timeout)

def bash(script, capture=False, timeout=60, input_data=None):
    return wsl("bash", "-lc", script, capture=capture, timeout=timeout, input_data=input_data)

def wsl_path(path):
    result = bash("wslpath -a -u " + shlex.quote(str(Path(path).resolve())), capture=True)
    if result.returncode != 0 or not result.stdout.strip():
        raise RuntimeError(f"Could not convert Windows path to WSL path: {path}\n{result.stderr}")
    return result.stdout.strip()

def shellq(value):
    return shlex.quote(str(value))

# Do not let old docking outputs survive a new attempt. The pipeline archives
# the previous target result directory before the new target run begins.
if output_folder.exists():
    shutil.rmtree(output_folder)
output_folder.mkdir(parents=True, exist_ok=True)

try:
    gpu = bash("command -v nvidia-smi >/dev/null 2>&1 && nvidia-smi --query-gpu=name,driver_version --format=csv,noheader", capture=True)
except subprocess.TimeoutExpired:
    raise SystemExit("ERROR: timed out while checking the WSL2 GPU. Restart WSL with 'wsl --shutdown' and retry.")
if gpu.returncode != 0 or not gpu.stdout.strip():
    raise SystemExit("ERROR: NVIDIA GPU is not visible inside WSL2. Update the Windows NVIDIA driver with WSL support and verify nvidia-smi inside Ubuntu.")
print("WSL GPU:", gpu.stdout.strip().replace("\n", " | "), flush=True)

home_result = bash('printf %s "$HOME"', capture=True)
linux_home = home_result.stdout.strip()
if home_result.returncode != 0 or not linux_home.startswith("/"):
    raise SystemExit("ERROR: could not resolve the WSL Linux home directory.")
micromamba = f"{linux_home}/.onedock/micromamba/micromamba"
mamba_root = f"{linux_home}/.onedock/mamba"
mamba_prefix = f"MAMBA_ROOT_PREFIX={shellq(mamba_root)}"
probe = bash(f"{mamba_prefix} {shellq(micromamba)} run -n onedock unidock --version", capture=True)
if probe.returncode != 0:
    print("UniDock Linux environment missing; installing UniDock 1.2.0...", flush=True)
    install = f'''set -euo pipefail
mkdir -p {shellq(str(Path(micromamba).parent))}
if [ ! -x {shellq(micromamba)} ]; then
  curl -L --fail --silent --show-error https://micro.mamba.pm/api/micromamba/linux-64/latest -o /tmp/micromamba.tar.bz2
  tar -xjf /tmp/micromamba.tar.bz2 -C {shellq(str(Path(micromamba).parent))} --strip-components=1 bin/micromamba
fi
export MAMBA_ROOT_PREFIX={shellq(mamba_root)}
{shellq(micromamba)} create -y -n onedock -c conda-forge unidock=1.2.0
{shellq(micromamba)} run -n onedock unidock --version'''
    install_result = bash(install, capture=True, timeout=1800)
    print(install_result.stdout or "", flush=True)
    if install_result.stderr:
        print(install_result.stderr, flush=True)
    if install_result.returncode != 0:
        raise SystemExit("ERROR: automatic UniDock installation failed. See the setup output above.")

# Stage files on the Linux filesystem: UniDock can be unreliable on /mnt/c paths,
# and WSL may not be able to enumerate OneDrive's cloud-backed project folder.
# Stream a tar archive over stdin so staging does not depend on WSL mounting the
# Windows project directory at all.
run_id = f"{protein}_{os.getpid()}_{uuid.uuid4().hex[:8]}"
stage = f"{linux_home}/.onedock/runs/{run_id}"
stage_q = shellq(stage)
archive = io.BytesIO()
with tarfile.open(fileobj=archive, mode="w:gz") as tar:
    tar.add(receptor, arcname="receptor.pdbqt")
    for ligand in ligands:
        tar.add(ligand, arcname=f"ligands/{ligand.name}")
prep_script = f'''set -euo pipefail
rm -rf {stage_q}
mkdir -p {stage_q}/results
tar -xzf - -C {stage_q}
find {stage_q}/ligands -maxdepth 1 -type f -name '*.pdbqt' -print | sort > {stage_q}/ligand_index.txt
count=$(wc -l < {stage_q}/ligand_index.txt)
if [ "$count" -eq 0 ]; then echo "No PDBQT ligands were transferred to WSL2" >&2; exit 3; fi
printf 'Staged ligand count: %s\\n' "$count"
'''
prep = bash(prep_script, capture=True, timeout=300, input_data=archive.getvalue())
if prep.returncode != 0:
    raise SystemExit(f"ERROR: could not stage docking inputs in WSL2.\n{prep.stdout}\n{prep.stderr}")

print("\n" + "=" * 70)
print(f"Protein      : {protein}")
print(f"Ligands      : {len(ligands)}")
print(f"Pocket       : {pocket_source}")
print(f"Box center   : {'  '.join(f'{v:.3f}' for v in center)}")
print(f"Box size     : {'  '.join(f'{v:.2f}' for v in size)}")
print(f"Search Mode  : {SEARCH_MODE}")
print("GPU backend  : WSL2 / Linux UniDock 1.2.0")
print("WSL staging  : Linux filesystem (~/.onedock/runs/...)")
print("=" * 70, flush=True)

parts = [micromamba, "run", "-n", "onedock", "unidock", "--receptor", f"{stage}/receptor.pdbqt",
         "--ligand_index", f"{stage}/ligand_index.txt", "--center_x", str(center[0]), "--center_y", str(center[1]),
         "--center_z", str(center[2]), "--size_x", str(round(size[0], 3)), "--size_y", str(round(size[1], 3)),
         "--size_z", str(round(size[2], 3)), "--search_mode", SEARCH_MODE, "--exhaustiveness", str(EXHAUSTIVENESS),
         "--num_modes", str(NUM_MODES), "--dir", f"{stage}/results", "--verbosity", "2"]
run_command = mamba_prefix + " exec " + " ".join(shellq(part) for part in parts)
cmd = [wsl_exe, "-d", WSL_DISTRO, "--", "bash", "-lc", run_command]
log_file = LOG_DIR / f"{protein}_unidock.log"
print("\nLaunching UniDock GPU inside WSL2...\n", flush=True)
start = time.monotonic()
with log_file.open("w", encoding="utf-8") as log:
    log.write("=== UniDock command ===\n" + " ".join(shellq(x) for x in parts) + "\n\n")
    log.write(f"Linux staging: {stage}\nExpected ligands: {len(ligands)}\n")
    log.flush()
    proc = subprocess.Popen(cmd, stdout=log, stderr=subprocess.STDOUT, text=True)
    previous = -1
    try:
        while proc.poll() is None:
            probe_count = bash(f"find {stage_q}/results -maxdepth 1 -type f -name '*_out.pdbqt' | wc -l", capture=True, timeout=30)
            try:
                done = int(probe_count.stdout.strip()) if probe_count.returncode == 0 else 0
            except ValueError:
                done = 0
            if done != previous:
                print(f"Docking progress: {done}/{len(ligands)}", flush=True)
                previous = done
            time.sleep(1)
        proc.wait()
    except BaseException:
        # Do not leave a GPU job running in WSL if GUI cancellation, a timeout,
        # or a progress-probe error interrupts the Python monitor.
        if proc.poll() is None:
            try:
                proc.terminate()
                proc.wait(timeout=15)
            except Exception:
                try:
                    proc.kill()
                    proc.wait(timeout=10)
                except Exception:
                    pass
        raise

# Copy any generated poses back even on non-zero exit, so partial results can be inspected.
copy_script = f'''set -euo pipefail
mkdir -p {shellq(wsl_path(output_folder))}
shopt -s nullglob
files=({stage_q}/results/*_out.pdbqt)
if ((${{#files[@]}} > 0)); then cp -f -- "${{files[@]}}" {shellq(wsl_path(output_folder))}/; fi
printf '%s\\n' "${{#files[@]}}"
'''
copy_back = bash(copy_script, capture=True, timeout=300)
if copy_back.returncode != 0:
    raise SystemExit(f"ERROR: failed to copy UniDock result files back to the project.\n{copy_back.stdout}\n{copy_back.stderr}\nStaging retained at: {stage}")

produced = sorted(output_folder.glob("*_out.pdbqt"))
expected_names = {p.stem for p in ligands}
observed_names = {p.name[:-len("_out.pdbqt")] for p in produced}
missing = sorted(expected_names - observed_names)
extra = sorted(observed_names - expected_names)
if proc.returncode != 0:
    raise SystemExit(f"ERROR: UniDock failed inside WSL2 (exit code {proc.returncode}). Partial poses: {len(produced)}/{len(ligands)}. See {log_file}. WSL staging retained at {stage}.")
if not produced:
    raise SystemExit(f"ERROR: UniDock returned success but produced no poses. See {log_file}. WSL staging retained at {stage}.")

pd.DataFrame([{"Expected_Ligands": len(ligands), "Produced_Poses": len(produced),
               "Missing_Ligands": len(missing), "Extra_Output_Files": len(extra),
               "Coverage_Percent": round(100 * len(observed_names & expected_names) / len(expected_names), 2),
               "UniDock_Exit_Code": proc.returncode, "Runtime_s": round(time.monotonic() - start, 2)}]).to_csv(
    output_folder / "Docking_Run_QC.csv", index=False)
if missing:
    pd.DataFrame({"Ligand": missing}).to_csv(output_folder / "Missing_Docking.csv", index=False)
if extra:
    pd.DataFrame({"Output": extra}).to_csv(output_folder / "Unexpected_Docking_Output.csv", index=False)

# Keep stage files when UniDock exited cleanly but returned partial output; they are useful for diagnosis.
if not missing and not extra:
    bash(f"rm -rf {stage_q}", capture=True, timeout=120)

print("\n" + "=" * 70)
print("GPU Batch Docking Finished" if not missing else "GPU Batch Docking Finished with Missing Ligands")
print("=" * 70)
print(f"Runtime      : {time.monotonic() - start:.1f} s")
print(f"Poses        : {len(produced)}/{len(ligands)}")
print(f"Missing      : {len(missing)}")
print(f"Results saved: {output_folder}")
print(f"Detailed log : {log_file}")
if missing:
    print("WARNING: docking coverage is incomplete; inspect Missing_Docking.csv and Docking_Run_QC.csv.")
