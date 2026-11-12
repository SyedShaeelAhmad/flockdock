import os
import sys
import subprocess
import shlex
import shutil
import tarfile
import urllib.request
from pathlib import Path

# ============================================================
# FlockDOCK 1.0.9 - P2Rank pocket prediction
# Robust Windows launcher + automatic P2Rank installation.
# ============================================================

PROJECT_DIR = os.environ.get(
    "UNIDOCK_PROJECT",
    os.environ.get("ONEDOCK_PROJECT", os.getcwd())
)
PROTEIN_DIR = os.path.join(PROJECT_DIR, "04_proteins")
OUTPUT_DIR = os.path.join(PROJECT_DIR, "06_grid")
LOG_DIR = os.path.join(PROJECT_DIR, "08_logs")
P2RANK_DIR = os.environ.get(
    "ONEDOCK_P2RANK_DIR",
    os.path.join(PROJECT_DIR, "tools", "p2rank")
)
P2RANK_OVERRIDE = os.environ.get("ONEDOCK_P2RANK", "").strip()
P2RANK_VERSION = "2.5.1"
P2RANK_URL = f"https://github.com/rdk/p2rank/releases/download/{P2RANK_VERSION}/p2rank_{P2RANK_VERSION}.tar.gz"

os.makedirs(OUTPUT_DIR, exist_ok=True)
os.makedirs(LOG_DIR, exist_ok=True)

def java_version(java):
    try:
        r = subprocess.run([java, "-version"], capture_output=True, text=True, timeout=15)
        text = (r.stderr or "") + "\n" + (r.stdout or "")
        import re
        m = re.search(r'version "([^" ]+)"', text)
        if not m:
            return None
        v = m.group(1)
        if v.startswith("1."):
            return int(v.split(".")[1])
        return int(v.split(".")[0])
    except Exception:
        return None

def find_java17():
    candidates=[]
    configured=os.environ.get("ONEDOCK_JAVA", "").strip()
    if configured:
        candidates.append(configured)
    jh=os.environ.get("JAVA_HOME", "").strip()
    if jh:
        candidates.append(os.path.join(jh, "bin", "java.exe"))
    # Project-local Java override written by setup_FlockDOCK.bat
    marker=Path(PROJECT_DIR)/"tools"/"java17_path.txt"
    try:
        if marker.is_file():
            candidates.append(marker.read_text(encoding="utf-8").strip())
    except Exception:
        pass
    for name in ("java.exe", "java"):
        try:
            candidates.extend(subprocess.check_output(["where.exe", name], text=True, stderr=subprocess.DEVNULL).splitlines())
        except Exception:
            pass
    # Common Temurin/Java install locations; recursively inspect only Java executables.
    for base in [r"C:\Program Files\Eclipse Adoptium", r"C:\Program Files\Java", r"C:\Program Files\Microsoft\jdk"]:
        if os.path.isdir(base):
            for dp, _, fs in os.walk(base):
                if "java.exe" in fs:
                    candidates.append(os.path.join(dp, "java.exe"))
    seen=set()
    for c in candidates:
        c=os.path.abspath(os.path.expandvars(c.strip('"')))
        if not c or c in seen or not os.path.isfile(c):
            continue
        seen.add(c)
        v=java_version(c)
        if v is not None and v >= 17:
            return c, v
    return None, None

def java_check():
    java, version = find_java17()
    if not java:
        return False, "Java 17+ was not found. A Java 8/11 installation is not sufficient for P2Rank 2.5.1."
    return True, f"Java {version} at {java}"

def git_bash():
    # Never mistake WSL's C:\Windows\System32\bash.exe for Git Bash.
    # P2Rank's Windows launcher is a Bash script and must run under Git for Windows.
    candidates = [
        os.environ.get("ONEDOCK_BASH", ""),
        r"C:\Program Files\Git\bin\bash.exe",
        r"C:\Program Files\Git\usr\bin\bash.exe",
        r"C:\Program Files (x86)\Git\bin\bash.exe",
        r"C:\Program Files (x86)\Git\usr\bin\bash.exe",
    ]
    git = shutil.which("git.exe") or shutil.which("git")
    if git:
        git_root = Path(git).resolve().parent.parent
        candidates += [str(git_root / "bin" / "bash.exe"), str(git_root / "usr" / "bin" / "bash.exe")]
    for c in candidates:
        if not c or not os.path.isfile(c):
            continue
        try:
            r = subprocess.run([c, "--version"], capture_output=True, text=True, timeout=10)
            text = ((r.stdout or "") + " " + (r.stderr or "")).lower()
            if "mingw" in text or "msys" in text or "git for windows" in text:
                return c
        except Exception:
            continue
    return None

def candidate_paths():
    candidates = []
    if P2RANK_OVERRIDE:
        candidates.append(P2RANK_OVERRIDE)
    for base in [
        P2RANK_DIR,
        os.path.join(P2RANK_DIR, f"p2rank_{P2RANK_VERSION}"),
    ]:
        # P2Rank officially recommends Git Bash on Windows because its
        # native cmd/PowerShell launchers can suffer from command-line
        # formatting/path issues. Prefer the Unix launcher when Git Bash is
        # available; keep the .bat launcher only as a fallback.
        if os.name == "nt" and git_bash():
            candidates += [
                os.path.join(base, "prank"),
                os.path.join(base, "distro", "prank"),
                os.path.join(base, "prank.exe"),
                os.path.join(base, "prank.bat"),
                os.path.join(base, "distro", "prank.bat"),
            ]
        else:
            candidates += [
                os.path.join(base, "prank.exe"),
                os.path.join(base, "prank.bat"),
                os.path.join(base, "prank"),
                os.path.join(base, "distro", "prank.bat"),
                os.path.join(base, "distro", "prank"),
            ]
    if os.path.isdir(P2RANK_DIR):
        for base, _, files in os.walk(P2RANK_DIR):
            for fn in files:
                if fn.lower() in ("prank.bat", "prank.exe", "prank"):
                    candidates.append(os.path.join(base, fn))
    seen=set()
    for c in candidates:
        if not c: continue
        c=os.path.abspath(c)
        if c in seen: continue
        seen.add(c)
        if os.path.isfile(c):
            return c
    return None

def safe_extract(tf, path):
    base=os.path.abspath(path)
    for member in tf.getmembers():
        target=os.path.abspath(os.path.join(path, member.name))
        if os.path.commonpath([base,target]) != base:
            raise RuntimeError(f"Unsafe P2Rank archive member: {member.name}")
    tf.extractall(path)

def install_p2rank():
    os.makedirs(P2RANK_DIR, exist_ok=True)
    archive=os.path.join(P2RANK_DIR, f"p2rank_{P2RANK_VERSION}.tar.gz")
    print(f"\nP2Rank {P2RANK_VERSION} not found.")
    print("Downloading:", P2RANK_URL)
    try:
        req=urllib.request.Request(P2RANK_URL, headers={"User-Agent":"FlockDOCK/1.0.9"})
        with urllib.request.urlopen(req, timeout=120) as r, open(archive,"wb") as f:
            shutil.copyfileobj(r,f)
        print("Extracting P2Rank...")
        with tarfile.open(archive,"r:gz") as tf:
            safe_extract(tf,P2RANK_DIR)
        os.remove(archive)
    except Exception as e:
        try:
            if os.path.exists(archive): os.remove(archive)
        except Exception: pass
        print("P2Rank installation failed:",e)
        return None
    found=candidate_paths()
    if found: print("P2Rank installed:",found)
    return found

def ensure_p2rank():
    return candidate_paths() or install_p2rank()

def run_p2rank(executable, input_pdb, output_folder):
    # IMPORTANT on Windows: P2Rank's prank.bat uses files relative to its own
    # installation directory. Running it from the FlockDOCK project directory
    # causes the Windows error "The system cannot find the path specified."
    # Always launch the batch file from its directory.
    args=["predict","-f",input_pdb,"-o",output_folder]
    low=executable.lower()
    if os.name=="nt" and low.endswith(".bat"):
        exe_dir=os.path.dirname(os.path.abspath(executable))
        # Use CALL so cmd.exe executes the batch wrapper correctly even when the
        # project path contains spaces (e.g. OneDrive\Desktop).
        exe=os.path.abspath(executable)
        # Batch files must be launched through cmd.exe as a single command string.
        # Passing the /c command as a separate argv element can make Windows
        # preserve escaped quotes (\"...\"), causing:
        #   '\"C:\\...\\prank.bat\"' is not recognized...
        # Use shell=True here deliberately: on Windows this routes the command
        # through cmd.exe and preserves normal Windows quoting for paths with spaces.
        command_line='call "'+exe+'" ' + ' '.join('"'+str(a)+'"' for a in args)
        return subprocess.run(command_line,capture_output=True,text=True,timeout=600,cwd=exe_dir,shell=True)
    elif os.name=="nt" and os.path.basename(executable).lower()=="prank":
        bash=git_bash()
        if not bash:
            raise RuntimeError("P2Rank's Unix launcher 'prank' was found, but Git Bash was not found. Install Git for Windows or set ONEDOCK_BASH to bash.exe.")
        # P2Rank's own documentation recommends Git Bash on Windows.
        # Run from the P2Rank installation directory so its relative paths
        # (models/, config/, etc.) resolve correctly. Convert all Windows
        # paths with cygpath rather than relying on shell path heuristics.
        def bash_path(p):
            r=subprocess.run(
                [bash, "-lc", "cygpath -u -- " + shlex.quote(os.path.abspath(p))],
                capture_output=True, text=True, timeout=10
            )
            if r.returncode==0 and r.stdout.strip():
                return r.stdout.strip()
            return os.path.abspath(p).replace("\\","/")
        exe=bash_path(executable)
        inp=bash_path(input_pdb)
        out=bash_path(output_folder)
        exe_dir=bash_path(os.path.dirname(os.path.abspath(executable)))
        # Use a shell command so ./prank is interpreted by Bash, not by
        # Windows CreateProcess. shlex.quote protects OneDrive/space paths.
        java17, java_ver = find_java17()
        if not java17:
            raise RuntimeError("P2Rank 2.5.1 requires Java 17+. No compatible Java runtime was found.")
        # P2Rank 2.5.1 is compiled for Java 17 (class version 61).
        # Explicitly set JAVA_HOME so an older Java 8/11 on PATH can never be selected.
        java_home_win = os.path.dirname(os.path.dirname(java17))
        java_home = bash_path(java_home_win)
        command_line = (
            "export JAVA_HOME=" + shlex.quote(java_home) +
            " && export PATH=\"$JAVA_HOME/bin:$PATH\""
            " && cd " + shlex.quote(exe_dir) +
            " && ./prank predict -f " + shlex.quote(inp) +
            " -o " + shlex.quote(out)
        )
        command=[bash,"-lc",command_line]
    else:
        command=[executable]+args
    print("\nExecuting:\n"+" ".join(f'"{x}"' if " " in str(x) else str(x) for x in command)+"\n")
    return subprocess.run(command,capture_output=True,text=True,timeout=600,cwd=PROJECT_DIR)

def locate_prediction(output_folder, protein):
    candidates=[
        os.path.join(output_folder,f"{protein}_clean.pdb_predictions.csv"),
        os.path.join(output_folder,"predictions.csv"),
    ]
    candidates += [
        os.path.join(output_folder,f)
        for f in os.listdir(output_folder) if f.lower().endswith("_predictions.csv")
    ] if os.path.isdir(output_folder) else []
    for c in candidates:
        if os.path.isfile(c) and os.path.getsize(c)>0:
            return c
    return None

def valid_prediction_rows(path):
    """Count only finite, plausible P2Rank pocket rows."""
    import csv
    import math
    try:
        with open(path, newline="", encoding="utf-8-sig") as handle:
            reader = csv.DictReader(handle)
            required = {"rank", "center_x", "center_y", "center_z", "score", "probability"}
            if not reader.fieldnames or not required.issubset({x.strip() for x in reader.fieldnames}):
                return 0
            count = 0
            for row in reader:
                try:
                    values = {str(k).strip(): v for k, v in row.items() if k}
                    nums = {k: float(values[k]) for k in required}
                    if not all(math.isfinite(v) for v in nums.values()):
                        continue
                    if nums["rank"] < 1 or not nums["rank"].is_integer():
                        continue
                    if any(abs(nums[k]) > 100000 for k in ("center_x", "center_y", "center_z")):
                        continue
                    if not 0.0 <= nums["probability"] <= 1.0:
                        continue
                    count += 1
                except (TypeError, ValueError, KeyError):
                    continue
            return count
    except Exception:
        return 0

if len(sys.argv)!=2:
    print("Usage: python 07_predict_pocket.py PROTEIN"); sys.exit(1)
protein=sys.argv[1]
input_pdb=os.path.join(PROTEIN_DIR,protein,f"{protein}_clean.pdb")
output_folder=os.path.join(OUTPUT_DIR,protein)
log_file=os.path.join(LOG_DIR,f"{protein}_p2rank.log")
report_file=os.path.join(output_folder,f"{protein}_p2rank_report.txt")
os.makedirs(output_folder,exist_ok=True)
# Invalidate previous predictions before *any* prerequisite checks. If Java or
# P2Rank is missing, step 08 must never fall through to an old pocket CSV.
for old_name in os.listdir(output_folder):
    if old_name.lower().endswith("_predictions.csv") or old_name.lower() == "predictions.csv":
        try:
            os.remove(os.path.join(output_folder, old_name))
        except OSError:
            pass
with open(log_file, "w", encoding="utf-8") as f:
    f.write(f"FlockDOCK 1.0.9 P2Rank attempt\nProtein: {protein}\nInput: {input_pdb}\n")
if not os.path.exists(input_pdb):
    print("ERROR: Protein not found:",input_pdb); sys.exit(1)
java_ok,java_info=java_check()
if not java_ok:
    with open(log_file, "a", encoding="utf-8") as f:
        f.write("\nPREFLIGHT FAILED\n" + java_info + "\n")
    print("\nERROR: P2Rank requires Java 17+.\n"+java_info); sys.exit(1)
P2RANK=ensure_p2rank()
if not P2RANK:
    with open(log_file, "a", encoding="utf-8") as f:
        f.write("\nPREFLIGHT FAILED\nP2Rank executable not found or installation failed.\n")
    print("\nERROR: P2Rank could not be installed/found.")
    print("Expected under:",P2RANK_DIR)
    sys.exit(1)
print("="*60); print("Predicting binding pocket:",protein); print("="*60)
print("P2Rank executable:",P2RANK)
try:
    result=run_p2rank(P2RANK,input_pdb,output_folder)
except Exception as e:
    with open(log_file, "a", encoding="utf-8") as f:
        f.write(f"\nLAUNCH FAILED\n{type(e).__name__}: {e}\n")
    print("\nERROR launching P2Rank:",e); sys.exit(1)
with open(log_file,"w",encoding="utf-8") as f:
    f.write("COMMAND\n=======\n"+" ".join(map(str,result.args))+"\n\nSTDOUT\n======\n"+(result.stdout or "")+"\n\nSTDERR\n======\n"+(result.stderr or ""))
prediction=locate_prediction(output_folder,protein)
if result.returncode==0 and prediction:
    n_pockets = valid_prediction_rows(prediction)
    if n_pockets < 1:
        with open(log_file, "a", encoding="utf-8") as f:
            f.write("\\n\\nVALIDATION\\n==========\\nP2Rank completed but returned zero usable pocket rows.\\n")
        print("\\nP2Rank ran but found no usable binding pockets for this structure.")
        print("This target will be marked as skipped; choose a different structure or define a pocket manually.")
        print("P2Rank log:", log_file)
        sys.exit(2)
    canonical=os.path.join(output_folder,f"{protein}_clean.pdb_predictions.csv")
    if os.path.abspath(prediction)!=os.path.abspath(canonical):
        shutil.copy2(prediction,canonical)
        prediction=canonical
    with open(report_file,"w",encoding="utf-8") as f:
        f.write("FlockDOCK P2Rank Prediction Report\n\n")
        f.write(f"Protein : {protein}\nInput   : {input_pdb}\nP2Rank  : {P2RANK}\nOutput  : {prediction}\nStatus  : SUCCESS\n")
    print("\n"+"="*60+"\nSUCCESS\n"+"="*60)
    print("Usable pockets:", n_pockets)
    print("Prediction:",prediction)
    print("Report:",report_file)
    print("Log:",log_file)
    sys.exit(0)
print("\n"+"="*60+"\nFAILED\n"+"="*60)
print("Log:",log_file)
if result.stderr.strip(): print(result.stderr)
sys.exit(1)
