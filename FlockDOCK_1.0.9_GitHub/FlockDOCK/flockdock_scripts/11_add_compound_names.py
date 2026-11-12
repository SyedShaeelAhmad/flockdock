import os
import sys
import pandas as pd

# ============================================================
# Annotate Docking Results
# ============================================================

PROJECT_DIR = os.environ.get("UNIDOCK_PROJECT", os.environ.get("ONEDOCK_PROJECT", os.getcwd()))

RESULT_DIR = os.path.join(
    PROJECT_DIR,
    "07_results"
)

MASTER_COMPOUND_FILE = os.path.join(
    PROJECT_DIR,
    "01_csv",
    "Docking_Library.csv"
)

# ============================================================
# Filter unrealistic docking scores
# ============================================================

MIN_REALISTIC_AFFINITY = -15.0

# ============================================================
# Input
# ============================================================

if len(sys.argv) != 2:

    print("Usage:")
    print("python 11_add_compound_names.py STAT3")
    sys.exit(1)

protein = sys.argv[1]

# ============================================================
# Files
# ============================================================

docking_file = os.path.join(
    RESULT_DIR,
    protein,
    "Docking_Scores.csv"
)

outfile = os.path.join(
    RESULT_DIR,
    protein,
    "Docking_Scores_Annotated.csv"
)

# ============================================================
# Check Files
# ============================================================

if not os.path.exists(docking_file):

    print("\nERROR")
    print("Docking results not found:")
    print(docking_file)
    sys.exit(1)

if not os.path.exists(MASTER_COMPOUND_FILE):

    print("\nERROR")
    print("Master compound file not found:")
    print(MASTER_COMPOUND_FILE)
    sys.exit(1)

# ============================================================
# Read Files
# ============================================================

dock = pd.read_csv(docking_file)

master = pd.read_csv(MASTER_COMPOUND_FILE)

# ============================================================
# Validate Columns
# ============================================================

required_dock = [

    "Ligand",
    "Best_Affinity"

]

# Compound metadata are optional for docking. Only the ligand identifier is
# required for annotation; Gene/target metadata are not used to select ligands.

missing = [

    c for c in required_dock

    if c not in dock.columns

]

if missing:

    print("\nERROR")
    print("Docking file missing columns:")
    print(missing)
    sys.exit(1)

# Resolve common alternative names from the user's compound CSV.
def resolve_col(df, aliases):
    lower = {str(c).strip().lower(): c for c in df.columns}
    for a in aliases:
        if a.lower() in lower:
            return lower[a.lower()]
    for a in aliases:
        for c in df.columns:
            if a.lower() in str(c).lower():
                return c
    return None

id_col = resolve_col(master, ["LOTUS_ID", "Ligand_LOTUS_ID", "ligand_id", "compound_id", "molecule_id", "id"])
if id_col is None:
    print("\nERROR: Could not identify the compound ID column in the master library.")
    print("Available columns:", list(master.columns))
    sys.exit(1)

# Normalize the master library to the names used by the result table.
rename_map = {id_col: "LOTUS_ID"}
for logical, aliases in {
    "compound_name": ["compound_name", "Compound_Name", "Ligand_name", "name"],
    "Canonical_SMILES": ["Canonical_SMILES", "Ligand_SMILES", "SMILES", "smiles"],
    "Gene": ["Gene", "Target_Gene", "Original_Target_Gene"],
    "target_name": ["target_name", "Target_Name", "LOTUS_target_name"]
}.items():
    col = resolve_col(master, aliases)
    if col is not None and col not in rename_map:
        rename_map[col] = logical
master = master.rename(columns=rename_map)

for col in ["compound_name", "Canonical_SMILES", "Gene", "target_name"]:
    if col not in master.columns:
        master[col] = ""

# ============================================================
# Remove unrealistic affinities
# ============================================================

original_count = len(dock)

dock = dock[
    dock["Best_Affinity"] >= MIN_REALISTIC_AFFINITY
].copy()

removed = original_count - len(dock)

if dock.empty:

    print("\nERROR")
    print("No valid docking results remain after filtering.")
    sys.exit(1)

# ============================================================
# Keep Required Columns
# ============================================================

master = (

    master[

        [

            "LOTUS_ID",
            "compound_name",
            "Canonical_SMILES",
            "Gene",
            "target_name"

        ]

    ]

    .drop_duplicates(
        subset="LOTUS_ID"
    )

)

# ============================================================
# Merge
# ============================================================

annotated = dock.merge(

    master,

    left_on="Ligand",

    right_on="LOTUS_ID",

    how="left"

)

# ============================================================
# Save Missing Annotation
# ============================================================

missing_annotation = annotated[
    annotated["compound_name"].isna()
]

if len(missing_annotation) > 0:

    missing_file = os.path.join(

        RESULT_DIR,

        protein,

        "Missing_Annotation.csv"

    )

    missing_annotation[

        ["Ligand"]

    ].rename(

        columns={

            "Ligand":"LOTUS_ID"

        }

    ).to_csv(

        missing_file,

        index=False

    )

# ============================================================
# Remove duplicate LOTUS_ID
# ============================================================

annotated = annotated.drop(columns=["LOTUS_ID"])

# ============================================================
# Rename Columns
# ============================================================

annotated = annotated.rename(columns={
    "Ligand": "LOTUS_ID",
    "compound_name": "Compound_Name",
    "Gene": "Original_Target_Gene",
    "target_name": "Original_Target_Name",
})

# ============================================================
# Add Docked Protein
# ============================================================

annotated["Docked_To"] = protein

# ============================================================
# Sort Results
# ============================================================

annotated = annotated.sort_values(

    by=[

        "Best_Affinity",

        "LOTUS_ID"

    ],

    ascending=[

        True,

        True

    ]

).reset_index(drop=True)

# ============================================================
# Add Rank
# ============================================================

annotated.insert(

    0,

    "Rank",

    range(

        1,

        len(annotated)+1

    )

)

# ============================================================
# Arrange Columns
# ============================================================

annotated = annotated[

    [

        "Rank",

        "LOTUS_ID",

        "Compound_Name",

        "Canonical_SMILES",

        "Original_Target_Gene",

        "Original_Target_Name",

        "Docked_To",

        "Best_Affinity"

    ]

]

# ============================================================
# Save
# ============================================================

annotated.to_csv(

    outfile,

    index=False

)

# ============================================================
# Annotation Statistics
# ============================================================

annotated_count = (

    annotated["Compound_Name"]

    .notna()

    .sum()

)

missing_count = (

    len(annotated)

    - annotated_count

)

# ============================================================
# Summary
# ============================================================

print()

print("=" * 70)
print("Docking Results Annotated")
print("=" * 70)

print(f"Protein                 : {protein}")
print(f"Original Results        : {original_count}")
print(f"Removed Unrealistic     : {removed}")
print(f"Final Results           : {len(annotated)}")
print(f"Annotated Compounds     : {annotated_count}")
print(f"Missing Annotation      : {missing_count}")

print()

print("Best Hit")

print("-" * 70)

print(

    annotated.iloc[0].to_string()

)

print()

print("Top 20 Hits")

print("-" * 70)

print(

    annotated.head(20).to_string(
        index=False
    )

)

print()

print("=" * 70)
print("Saved")
print("=" * 70)

if os.path.exists(outfile):

    print(outfile)

if missing_count > 0:

    print("\nMissing annotations:")

    print(missing_file)