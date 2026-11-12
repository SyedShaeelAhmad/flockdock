import os
import sys
import subprocess

# ============================================================
# Prepare Receptor using Meeko
# ============================================================

PROJECT_DIR = os.environ.get("UNIDOCK_PROJECT", os.environ.get("ONEDOCK_PROJECT", os.getcwd()))

PROTEIN_DIR = os.path.join(PROJECT_DIR, "04_proteins")
OUTPUT_DIR = os.path.join(PROJECT_DIR, "05_protein_pdbqt")
LOG_DIR = os.path.join(PROJECT_DIR, "08_logs")

os.makedirs(OUTPUT_DIR, exist_ok=True)
os.makedirs(LOG_DIR, exist_ok=True)

# ============================================================
# Input
# ============================================================

if len(sys.argv) != 2:

    print("Usage:")
    print("python 06_prepare_receptor.py PROTEIN")
    sys.exit(1)

protein = sys.argv[1]

input_pdb = os.path.join(
    PROTEIN_DIR,
    protein,
    f"{protein}_clean.pdb"
)

output_base = os.path.join(
    OUTPUT_DIR,
    protein
)

output_pdbqt = output_base + ".pdbqt"

log_file = os.path.join(
    LOG_DIR,
    f"{protein}_meeko.log"
)

report_file = os.path.join(
    OUTPUT_DIR,
    f"{protein}_report.txt"
)

if not os.path.exists(input_pdb):

    print(f"\nERROR\nInput PDB not found:\n{input_pdb}")
    sys.exit(1)

# Remove old log

if os.path.exists(log_file):
    os.remove(log_file)

# ============================================================
# Meeko-compatible PDB sanitizer
# ============================================================
# Some structures (notably GPX4) contain selenocysteine (SEC) and/or
# selenium atoms.  The installed Meeko/RDKit combination does not support
# elemental Se and Meeko also expects one residue name per residue key.
# For receptor preparation only, create a temporary compatibility PDB:
#   SEC -> CYS
#   MSE -> MET
#   Se/SE atom -> S with the standard CYS/MET atom name (SG/SD)
# The original clean PDB is never modified.
#
# IMPORTANT: do the residue-name conversion on EVERY ATOM/HETATM record,
# not only on selenium-containing records.  GPX4 can contain mixed SEC/CYS
# records at the same residue number; converting only the Se line leaves
# Meeko with the fatal "{'SEC','CYS'}" residue-name conflict.

def make_meeko_compatible_pdb(src):

    out = os.path.join(
        os.path.dirname(src),
        f"{protein}_clean_meeko_compatible.pdb"
    )

    changed_residues = 0
    changed_atoms = 0
    removed_unsupported_metals = []
    seen_residue_conversion = set()

    # Meeko/RDKit currently raises during intermolecular-bond perception
    # when a Ni or Cu atom is present ("doesn't have an implemented covalent
    # radius").  Ni is removed ONLY from this temporary Meeko input copy;
    # the original *_clean.pdb remains untouched.  The report/log explicitly
    # records the removal because metal-dependent docking needs special care.
    unsupported_meeko_metals = {"NI", "CU"}

    with open(src, "r", errors="replace") as fin, open(out, "w", newline="") as fout:

        for line in fin:

            record = line[:6].strip() if len(line) >= 6 else ""

            if record in ("ATOM", "HETATM") and len(line) >= 20:

                resname = line[17:20].strip().upper()
                atomname = line[12:16].strip().upper()
                element = line[76:78].strip().upper() if len(line) >= 78 else ""

                # Infer element from atom name if the PDB element column is blank.
                if not element:
                    letters = "".join(ch for ch in atomname if ch.isalpha())
                    element = letters[:2].upper() if letters[:2].upper() in unsupported_meeko_metals else (letters[:1].upper() if letters else "")

                if element in unsupported_meeko_metals:
                    atom_id = line[6:11].strip()
                    residue_id = line[17:27].strip()
                    removed_unsupported_metals.append(f"{element} atom {atom_id} ({residue_id})")
                    continue

                # Normalize modified amino-acid residue names on all atoms.
                if resname == "SEC":
                    line = line[:17] + "CYS" + line[20:]
                    resname = "CYS"
                    key = line[21:27]
                    if key not in seen_residue_conversion:
                        seen_residue_conversion.add(key)
                        changed_residues += 1

                elif resname == "MSE":
                    line = line[:17] + "MET" + line[20:]
                    resname = "MET"
                    key = line[21:27]
                    if key not in seen_residue_conversion:
                        seen_residue_conversion.add(key)
                        changed_residues += 1

                # Normalize selenium itself to sulfur for Meeko compatibility.
                if element == "SE" or atomname == "SE":
                    new_atom = " SG " if resname == "CYS" else (" SD " if resname == "MET" else " S  ")
                    line = line[:12] + new_atom + line[16:76]
                    tail = line[78:] if len(line) > 78 else "\n"
                    line = line[:76] + " S" + tail
                    changed_atoms += 1

            fout.write(line)

    return out, changed_atoms, changed_residues, removed_unsupported_metals

meeko_input_pdb = input_pdb

# ============================================================
# Run Meeko
# ============================================================

def run_meeko(extra_args, description):

    if os.path.exists(output_pdbqt):
        os.remove(output_pdbqt)

    command = [

        sys.executable,

        "-m",

        "meeko.cli.mk_prepare_receptor",

        "--read_pdb", meeko_input_pdb,

        "-o", output_base,

        "-p"

    ]

    command.extend(extra_args)

    print()
    print("=" * 70)
    print(description)
    print("=" * 70)
    print()
    print("Executing:")
    print()
    print(" ".join(command))
    print()

    try:

        result = subprocess.run(

            command,

            capture_output=True,

            text=True,

            timeout=300

        )

    except subprocess.TimeoutExpired:

        print("ERROR : Meeko timed out.")

        with open(log_file, "a") as f:

            f.write("=" * 80 + "\n")
            f.write(description + "\n")
            f.write("TIMEOUT\n\n")

        return False

    # --------------------------------------------------------
    # Save complete log
    # --------------------------------------------------------

    with open(log_file, "a") as f:

        f.write("=" * 80 + "\n")
        f.write(description + "\n\n")
        f.write("COMMAND:\n")
        f.write(" ".join(command) + "\n\n")
        f.write("STDOUT:\n")
        f.write(result.stdout)
        f.write("\n\nSTDERR:\n")
        f.write(result.stderr)
        f.write("\n\n")

    if result.stdout.strip():
        print(result.stdout)

    if result.stderr.strip():
        print(result.stderr)

    success = (

        result.returncode == 0

        and

        os.path.exists(output_pdbqt)

        and

        os.path.getsize(output_pdbqt) > 100

    )

    return success

# ============================================================
# Preparation Strategies
# ============================================================

strategies = [

    (
        [],
        "Attempt 1 : Standard Meeko"
    ),

    (
        ["--allow_bad_res"],
        "Attempt 2 : allow_bad_res"
    ),

    (
        ["--default_altloc", "A"],
        "Attempt 3 : default_altloc A"
    ),

    (
        [
            "--allow_bad_res",
            "--default_altloc", "A"
        ],
        "Attempt 4 : allow_bad_res + default_altloc A"
    )

]

# ============================================================
# Try All Strategies
# ============================================================

print()
print("=" * 70)
print(f"Preparing receptor : {protein}")
print("=" * 70)

success = False
successful_strategy = None

for args, title in strategies:

    if run_meeko(args, title):

        success = True
        successful_strategy = title
        break

# ============================================================
# Selenium fallback
# ============================================================
# If all normal Meeko strategies fail, retry once with a sanitized copy
# for structures containing selenium. This specifically handles GPX4-like
# selenoprotein structures without changing the user's source PDB.

if not success:

    try:
        compatible_pdb, changed_atoms, changed_residues, removed_metals = make_meeko_compatible_pdb(input_pdb)

        if changed_atoms > 0 or removed_metals:
            meeko_input_pdb = compatible_pdb

            print()
            print("=" * 70)
            print(f"Compatibility fallback : {changed_atoms} Se atom(s) converted; {changed_residues} SEC/MSE residue(s) normalized")
            if removed_metals:
                print(f"Removed unsupported Ni/Cu atom(s) from temporary Meeko input: {len(removed_metals)}")
                for item in removed_metals:
                    print(f"  - {item}")
                print("WARNING: Ni/Cu is absent from the prepared receptor. Validate metal-dependent binding before interpreting docking results.")
            print("=" * 70)
            print(f"Temporary receptor : {compatible_pdb}")

            selenium_strategies = [
                (
                    ["--default_altloc", "A"],
                    "Attempt 5 : selenium-compatible PDB + default_altloc A"
                ),
                (
                    ["--allow_bad_res", "--default_altloc", "A"],
                    "Attempt 6 : selenium-compatible PDB + allow_bad_res + default_altloc A"
                )
            ]

            for args, title in selenium_strategies:
                if run_meeko(args, title):
                    success = True
                    successful_strategy = title
                    break

    except Exception as exc:
        print(f"Selenium compatibility fallback failed: {exc}")

# ============================================================
# Result
# ============================================================

if not success:

    print()
    print("=" * 70)
    print("FAILED")
    print("=" * 70)
    print(f"Could not prepare receptor : {protein}")
    print(f"Log file : {log_file}")
    sys.exit(1)

# ============================================================
# Save Report
# ============================================================

with open(report_file, "w") as f:

    f.write("=" * 60 + "\n")
    f.write("Receptor Preparation Report\n")
    f.write("=" * 60 + "\n\n")

    f.write(f"Protein : {protein}\n")
    f.write(f"Input   : {input_pdb}\n")
    f.write(f"Output  : {output_pdbqt}\n\n")
    if 'removed_metals' in locals() and removed_metals:
        f.write("WARNING: Unsupported Ni/Cu atom(s) were removed from the temporary Meeko input only; original clean PDB was preserved.\n")
        for item in removed_metals:
            f.write(f"  - {item}\n")
        f.write("Metal-dependent binding should be reviewed before interpreting scores.\n\n")
    if 'changed_atoms' in locals() and (changed_atoms > 0 or changed_residues > 0):
        f.write("WARNING: Selenium/SEC/MSE was normalized in the temporary Meeko input for parameterization compatibility. The original clean PDB was preserved.\n")
        f.write(f"  Selenium atoms converted: {changed_atoms}\n")
        f.write(f"  Modified residues normalized: {changed_residues}\n")
        f.write("This changes receptor chemistry; selenoprotein/active-site interpretation requires manual review.\n\n")

    f.write(f"Strategy : {successful_strategy}\n")
    f.write("Status   : SUCCESS\n")

# ============================================================
# Final Summary
# ============================================================

print()
print("=" * 70)
print("SUCCESS")
print("=" * 70)

print(f"Protein : {protein}")
print(f"Saved   : {output_pdbqt}")
print(f"Strategy: {successful_strategy}")
print(f"Report  : {report_file}")
print(f"Log     : {log_file}")