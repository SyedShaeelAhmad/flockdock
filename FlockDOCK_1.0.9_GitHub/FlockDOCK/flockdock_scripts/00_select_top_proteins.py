import os
import tempfile
import pandas as pd

# ============================================================
# Select ALL Proteins from the CSV
# ============================================================

PROJECT_DIR = os.environ.get("UNIDOCK_PROJECT", os.environ.get("ONEDOCK_PROJECT", os.getcwd()))

INPUT_CSV = os.path.join(
    PROJECT_DIR,
    "01_csv",
    "Docking_Library.csv"
)

OUTPUT_CSV = os.path.join(
    PROJECT_DIR,
    "proteins.csv"
)

# ============================================================
# Check Input
# ============================================================

if not os.path.exists(INPUT_CSV):

    raise FileNotFoundError(
        f"Input file not found:\n{INPUT_CSV}"
    )

os.makedirs(PROJECT_DIR, exist_ok=True)

# ============================================================
# Load Dataset
# ============================================================

df = pd.read_csv(INPUT_CSV)

# ============================================================
# Validate Columns
# ============================================================

required = ["Gene"]

missing = [c for c in required if c not in df.columns]

if missing:

    raise ValueError(
        f"Missing required columns: {missing}"
    )

# ============================================================
# Clean Gene Column
# ============================================================

df["Gene"] = (
    df["Gene"]
    .astype(str)
    .str.strip()
)

df = df[
    (~df["Gene"].isin(["", "nan", "None"]))
]

if df.empty:

    raise ValueError(
        "No valid genes found in Docking_Library.csv"
    )

# ============================================================
# Count Compounds per Protein
# ============================================================

protein_counts = (
    df.groupby("Gene")
      .size()
      .reset_index(name="Compound_Count")
      .sort_values(
          "Compound_Count",
          ascending=False
      )
)

# ============================================================
# Select ALL Proteins
# ============================================================

selected_proteins = protein_counts

# ============================================================
# Create Protein Folders
# ============================================================

PROTEIN_DIR = os.path.join(
    PROJECT_DIR,
    "04_proteins"
)

os.makedirs(PROTEIN_DIR, exist_ok=True)

for gene in selected_proteins["Gene"]:

    folder = os.path.join(
        PROTEIN_DIR,
        gene
    )

    os.makedirs(
        folder,
        exist_ok=True
    )

# ============================================================
# Save Safely (Atomic Write)
# ============================================================

fd, temp_file = tempfile.mkstemp(
    suffix=".csv",
    dir=PROJECT_DIR
)

os.close(fd)

selected_proteins.to_csv(
    temp_file,
    index=False
)

os.replace(
    temp_file,
    OUTPUT_CSV
)

# ============================================================
# Verify Output
# ============================================================

if (
    not os.path.exists(OUTPUT_CSV)
    or
    os.path.getsize(OUTPUT_CSV) == 0
):

    raise RuntimeError(
        "proteins.csv was not written correctly."
    )

# ============================================================
# Report
# ============================================================

print("=" * 65)
print("All Proteins Selected")
print("=" * 65)

print(selected_proteins)

print()

print(f"Saved to:\n{OUTPUT_CSV}")

print()

print(f"Total unique proteins : {len(protein_counts)}")
print(f"Proteins selected     : {len(selected_proteins)}")