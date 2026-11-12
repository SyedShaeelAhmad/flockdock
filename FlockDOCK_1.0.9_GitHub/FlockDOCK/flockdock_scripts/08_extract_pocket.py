import argparse
import os
import re
import sys

import pandas as pd

from ligand_filter import (AA3, ligand_instances, parse_pdb_atoms,
                           MIN_LIGAND_ATOMS)

# ============================================================
# Extract the docking pocket.
#
# Priority:
#   1. CO-CRYSTALLISED LIGAND POCKET
#      Step 04 now prefers PDB entries that actually have a ligand bound, so
#      the receptor in use normally carries one. This step re-derives the
#      pocket from that structure rather than trusting whichever HETATM group
#      step 05 happened to pick first:
#         - identify which downloaded structure step 05 built the clean file
#           from, by matching CA coordinates
#         - group its HETATM records into ligand instances
#         - discard crystallisation additives and anything outside
#           MIN/MAX_LIGAND_ATOMS
#         - keep only ligands that actually touch the selected chain
#         - prefer a drug-like ligand over a cofactor, then most contacts,
#           then most heavy atoms
#      The box is the ligand's own extent plus padding. Contact residues are
#      written out for the methods section.
#
#   2. P2RANK best pocket
#      Used only when the structure carries no usable ligand at all, which
#      after the step 04 change means an AlphaFold model or a genuinely
#      ligand-free target.
#
# Nothing else. Writes 06_grid/<P>/Pocket_Config.csv with the original
# 11-column schema, so 09_run_unidock.py needs no change.
#
# Called as before:   python 08_extract_pocket.py PROTEIN
# ============================================================

PROJECT_DIR = os.environ.get(
    "UNIDOCK_PROJECT",
    os.environ.get("ONEDOCK_PROJECT", os.getcwd()))
GRID_DIR = os.path.join(PROJECT_DIR, "06_grid")
PROTEIN_DIR = os.path.join(PROJECT_DIR, "04_proteins")

CONFIG_COLS = [
    "Protein", "Pocket_Source", "Pocket_Rank",
    "Center_X", "Center_Y", "Center_Z",
    "Size_X", "Size_Y", "Size_Z",
    "Score", "Probability",
]

ap = argparse.ArgumentParser(add_help=True)
ap.add_argument("protein")
ap.add_argument("--chain", default=None)
ap.add_argument("--pad", type=float, default=5.0,
                help="Angstrom added to each side of the ligand extent.")
ap.add_argument("--min-box", type=float, default=18.0)
ap.add_argument("--max-box", type=float, default=30.0)
ap.add_argument("--contact-cutoff", type=float, default=4.5)
ap.add_argument("--min-contacts", type=int, default=3,
                help="Protein residues a ligand must touch to count as bound "
                     "in a pocket rather than stuck on the surface.")
ap.add_argument("--contact-tolerance", type=float, default=0.20,
                help="A drug-like ligand outranks a cofactor only if its "
                     "contact count is within this fraction of the best. "
                     "0.20 means it must reach 80%% of the top contact count.")
ap.add_argument("--drugs-only", action="store_true",
                help="Refuse cofactors and substrates such as ATP or G6P.")
ap.add_argument("--custom-center", nargs=3, type=float, default=None,
                metavar=("X", "Y", "Z"),
                help="Use a user-supplied docking-box center instead of automatic pocket selection.")
ap.add_argument("--custom-size", nargs=3, type=float, default=None,
                metavar=("X", "Y", "Z"),
                help="Use a user-supplied docking-box size in Angstrom.")
ap.add_argument("--custom-residues", default=None,
                help="Use a user-supplied residue list to define the pocket, e.g. '151,154,157'.")
ap.add_argument("--custom-chain", default=None,
                help="Protein chain for --custom-residues, e.g. A.")
ap.add_argument("--custom-pad", type=float, default=5.0,
                help="Padding in Angstrom around a custom residue-defined pocket.")
args = ap.parse_args()

protein = args.protein
pdir = os.path.join(PROTEIN_DIR, protein)
gdir = os.path.join(GRID_DIR, protein)
os.makedirs(gdir, exist_ok=True)
outfile = os.path.join(gdir, "Pocket_Config.csv")
report_file = os.path.join(gdir, "Pocket_Report.txt")
evidence_file = os.path.join(gdir, "Pocket_Ligand_Evidence.csv")

TRACE = []


def trace(msg):
    TRACE.append(msg)
    print(f"  {msg}")


def write_config(**kw):
    row = {c: kw.get(c, "") for c in CONFIG_COLS}
    pd.DataFrame([row])[CONFIG_COLS].to_csv(outfile, index=False)


# ============================================================
# Helpers
# ============================================================

def discover_chain(atoms):
    if args.chain:
        return args.chain
    sel = os.path.join(pdir, f"{protein}_selection.txt")
    if os.path.exists(sel):
        m = re.search(r"^Chain\s*:\s*(\S+)", open(sel).read(), re.M)
        if m and m.group(1) not in ("-", "None"):
            return m.group(1)
    chains = {}
    for a in atoms:
        if a["resname"] in AA3:
            chains.setdefault(a["chain"], set()).add(a["resseq"])
    return sorted(chains, key=lambda c: -len(chains[c]))[0] if chains else ""


def ca_fingerprint(atoms, chain=None):
    return {(round(a["x"], 2), round(a["y"], 2), round(a["z"], 2))
            for a in atoms
            if a["name"] == "CA" and a["resname"] in AA3
            and (chain is None or a["chain"] == chain)}


def find_source_structure(clean_atoms, chain):
    """
    Which downloaded structure did step 05 build the clean file from?

    Step 05 does not record it, and the pipeline retries rank1..rank5, so the
    answer is not always rank1. Match on CA coordinates instead of guessing.
    """
    if not os.path.isdir(pdir):
        return None, {}
    target = ca_fingerprint(clean_atoms, chain)
    if not target:
        return None, {}

    candidates = []
    for fn in sorted(os.listdir(pdir)):
        if not fn.endswith(".pdb"):
            continue
        if fn in (f"{protein}_clean.pdb", f"{protein}_ligand.pdb"):
            continue
        candidates.append(fn)
    # rank files first, then the canonical copy, then AlphaFold models
    candidates.sort(key=lambda f: (0 if "_rank" in f else
                                   1 if f == f"{protein}.pdb" else 2, f))

    best, best_frac, atom_cache = None, 0.0, {}
    for fn in candidates:
        path = os.path.join(pdir, fn)
        try:
            atoms = parse_pdb_atoms(path)
        except Exception:
            continue
        if not atoms:
            continue
        atom_cache[fn] = atoms
        overlap = len(target & ca_fingerprint(atoms))
        frac = overlap / len(target)
        if frac > best_frac:
            best, best_frac = fn, frac
        if frac >= 0.98:
            break

    if best is None or best_frac < 0.80:
        trace(f"could not match the clean file to a downloaded structure "
              f"(best overlap {best_frac:.0%})")
        return None, atom_cache
    trace(f"clean file traced to {best} (CA overlap {best_frac:.0%})")
    return best, atom_cache


def contacts(lig_atoms, prot_atoms, cutoff):
    c2 = cutoff * cutoff
    hits = set()
    for p in prot_atoms:
        key = (p["chain"], p["resseq"], p["resname"])
        if key in hits:
            continue
        for l in lig_atoms:
            dx = p["x"] - l["x"]
            if dx * dx > c2:
                continue
            dy = p["y"] - l["y"]
            dz = p["z"] - l["z"]
            if dx * dx + dy * dy + dz * dz <= c2:
                hits.add(key)
                break
    return hits


def box_from_atoms(atoms):
    axes = ([a["x"] for a in atoms], [a["y"] for a in atoms],
            [a["z"] for a in atoms])
    center = [(min(v) + max(v)) / 2.0 for v in axes]
    size = [max(args.min_box, min(args.max_box,
                                  (max(v) - min(v)) + 2 * args.pad))
            for v in axes]
    return center, size


def box_from_atoms_with_pad(atoms, pad):
    axes = ([a["x"] for a in atoms], [a["y"] for a in atoms],
            [a["z"] for a in atoms])
    center = [(min(v) + max(v)) / 2.0 for v in axes]
    size = [max(args.min_box, min(args.max_box,
                                  (max(v) - min(v)) + 2 * pad))
            for v in axes]
    return center, size


def parse_residue_tokens(text):
    """Accept '151,154,157', '151 154 157', or 'A:151,A:154'."""
    tokens = [t.strip() for t in re.split(r"[,;\s]+", text or "") if t.strip()]
    numbers = []
    for token in tokens:
        value = token.split(":", 1)[1] if ":" in token else token
        value = value.strip()
        if not re.fullmatch(r"-?\d+", value):
            raise ValueError(f"Invalid residue '{token}'. Use residue numbers such as 151,154,157.")
        numbers.append(int(value))
    return sorted(set(numbers))


def custom_pocket(clean_atoms):
    """Build a Pocket_Config.csv from the user's explicit GUI choice."""
    if args.custom_center is not None or args.custom_size is not None:
        if args.custom_center is None or args.custom_size is None:
            raise ValueError("Custom coordinates require both center X/Y/Z and box X/Y/Z.")
        if any(v <= 0 for v in args.custom_size):
            raise ValueError("Custom box dimensions must be greater than zero.")
        center = [float(v) for v in args.custom_center]
        size = [float(v) for v in args.custom_size]
        detail = (f"user coordinates; center {center[0]:.3f}, {center[1]:.3f}, "
                  f"{center[2]:.3f}; box {size[0]:.3f}, {size[1]:.3f}, {size[2]:.3f}")
        trace("using user-supplied docking coordinates")
        return dict(center=center, size=size, source="user_coordinates",
                    detail=detail, rank="", score="", probability="")

    if args.custom_residues is not None:
        if not args.custom_chain:
            raise ValueError("A protein chain is required for a custom residue pocket.")
        chain = args.custom_chain.strip()
        residue_numbers = parse_residue_tokens(args.custom_residues)
        if not residue_numbers:
            raise ValueError("No custom residues were supplied.")
        selected = [a for a in clean_atoms
                    if a["resname"] in AA3 and a["chain"] == chain
                    and a["resseq"] in residue_numbers]
        missing = [str(n) for n in residue_numbers
                   if not any(a["resseq"] == n for a in selected)]
        if missing:
            raise ValueError(f"Residue(s) not found in chain {chain}: {', '.join(missing)}")
        center, size = box_from_atoms_with_pad(selected, args.custom_pad)
        labels = [f"{a['resname']}{a['resseq']}" for a in selected]
        labels = list(dict.fromkeys(labels))
        detail = (f"user residues in chain {chain}: {', '.join(labels)}; "
                  f"padding {args.custom_pad:.2f} A")
        trace(f"using user-supplied residues: chain {chain}, {', '.join(labels)}")
        trace(f"custom residue pocket contains {len(labels)} residue(s)")
        return dict(center=center, size=size, source="user_residues",
                    detail=detail, rank="", score="", probability="")

    return None


# ============================================================
# Tier 1 -- co-crystallised ligand pocket
# ============================================================

def tier_cocrystal(clean_atoms, chain):
    source, cache = find_source_structure(clean_atoms, chain)

    if source is None:
        # last resort inside this tier: trust step 05's own ligand box, but
        # only if the ligand it chose passes the filter
        return tier_ligand_box()

    atoms = cache[source]
    prot = [a for a in atoms if a["resname"] in AA3 and a["chain"] == chain]
    if not prot:
        prot = [a for a in atoms if a["resname"] in AA3]
        trace(f"chain {chain} not found in {source}; using all protein chains")
    if not prot:
        trace("no protein atoms in the source structure, Tier 1 declined")
        return None

    insts = ligand_instances(atoms, args.drugs_only)
    if not insts:
        trace(f"{source} contains no HETATM ligand groups at all")
        return None

    rows, viable = [], []
    for i in insts:
        n_contacts = 0
        residues = set()
        if i["verdict"] in ("drug", "cofactor"):
            residues = contacts(i["atoms"], prot, args.contact_cutoff)
            n_contacts = len(residues)
        rows.append(dict(
            Protein=protein, Source_Structure=source, Comp_ID=i["comp_id"],
            Verdict=i["verdict"], Reason=i["reason"],
            Heavy_Atoms=i["heavy_atoms"], Ligand_Chain=i["chain"],
            Ligand_Resseq=i["resseq"], Contact_Residues=n_contacts,
            Residues=";".join(f"{r[2]}{r[1]}"
                              for r in sorted(residues, key=lambda r: r[1]))))
        if i["verdict"] in ("drug", "cofactor") and n_contacts >= args.min_contacts:
            viable.append((i, residues))

    pd.DataFrame(rows).to_csv(evidence_file, index=False)

    kept = [r for r in rows if r["Verdict"] != "reject"]
    dropped = [r for r in rows if r["Verdict"] == "reject"]
    if dropped:
        trace("ignored: " + ", ".join(
            f"{r['Comp_ID']}({r['Heavy_Atoms']}, {r['Reason']})"
            for r in dropped[:8]))
    if kept:
        trace("candidates: " + ", ".join(
            f"{r['Comp_ID']}({r['Verdict']}, {r['Heavy_Atoms']} atoms, "
            f"{r['Contact_Residues']} contacts)" for r in kept[:8]))

    if not viable:
        trace(f"no ligand touches at least {args.min_contacts} residue(s) of "
              "the selected chain, Tier 1 declined")
        return None

    # Contact count decides, because the ligand burying the most of the
    # protein is the one sitting in the real pocket. A drug-like verdict only
    # breaks ties within --contact-tolerance of the best, so a small surface
    # fragment cannot outrank a cofactor occupying the catalytic site.
    best_contacts = max(len(r) for _, r in viable)
    threshold = best_contacts * (1.0 - args.contact_tolerance)
    viable.sort(key=lambda t: (
        0 if len(t[1]) >= threshold and t[0]["verdict"] == "drug" else
        1 if len(t[1]) >= threshold else 2,
        -len(t[1]), -t[0]["heavy_atoms"]))
    if len(viable) > 1:
        trace("ranked by contacts (drug breaks ties within "
              f"{args.contact_tolerance:.0%} of the best): "
              + ", ".join(f"{i['comp_id']}/{i['verdict']}/{len(r)}c"
                          for i, r in viable[:6]))
    best, residues = viable[0]
    center, size = box_from_atoms(best["atoms"])

    labels = [f"{r[2]}{r[1]}" for r in sorted(residues, key=lambda r: r[1])]
    trace(f"selected ligand {best['comp_id']} ({best['verdict']}, "
          f"{best['heavy_atoms']} heavy atoms) in chain {best['chain']}")
    trace(f"pocket lined by {len(labels)} residue(s): {', '.join(labels)}")

    return dict(
        center=center, size=size, source="cocrystal_ligand",
        detail=f"{best['comp_id']} ({best['verdict']}, "
               f"{best['heavy_atoms']} heavy atoms, {len(labels)} contacts) "
               f"from {source}",
        rank="", score="", probability="")


def tier_ligand_box():
    """Fallback within Tier 1: step 05's own box, validated."""
    box = os.path.join(pdir, f"{protein}_ligand_box.csv")
    if not os.path.exists(box):
        trace("no ligand box file from step 05 either")
        return None
    lb = pd.read_csv(box)
    lb.columns = lb.columns.str.strip()
    needed = {"Center_X", "Center_Y", "Center_Z", "Size_X", "Size_Y", "Size_Z"}
    if lb.empty or not needed.issubset(lb.columns):
        trace("step 05 ligand box is malformed")
        return None
    r = lb.iloc[0]
    ligand = str(r["Ligand"]).strip() if "Ligand" in lb.columns else ""

    lp = os.path.join(pdir, f"{protein}_ligand.pdb")
    n = len(parse_pdb_atoms(lp)) if os.path.exists(lp) else None
    from ligand_filter import classify
    verdict, reason = classify(ligand, n, args.drugs_only)
    if verdict == "reject":
        trace(f"step 05 chose {ligand or '?'} "
              f"({n if n is not None else '?'} heavy atoms) - {reason}. "
              "Rejected.")
        return None

    trace(f"using step 05 ligand box: {ligand} ({verdict}, {n} heavy atoms)")
    return dict(
        center=[float(r["Center_X"]), float(r["Center_Y"]), float(r["Center_Z"])],
        size=[float(r["Size_X"]), float(r["Size_Y"]), float(r["Size_Z"])],
        source="cocrystal_ligand",
        detail=f"{ligand} ({verdict}, {n} heavy atoms) via step 05 box",
        rank="", score="", probability="")


# ============================================================
# Tier 2 -- P2Rank
# ============================================================

def tier_p2rank():
    pred = os.path.join(gdir, f"{protein}_clean.pdb_predictions.csv")
    if not os.path.exists(pred):
        trace("no current P2Rank prediction file")
        return None
    try:
        df = pd.read_csv(pred)
    except Exception as exc:
        trace(f"could not parse P2Rank CSV: {exc}")
        return None
    if df.empty:
        trace("P2Rank prediction file is empty")
        return None
    df.columns = df.columns.str.strip()
    required = ["rank", "center_x", "center_y", "center_z", "score", "probability"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        trace("P2Rank file missing columns: " + ", ".join(missing))
        return None
    for col in required:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    df = df.replace([float("inf"), float("-inf")], float("nan")).dropna(subset=required)
    df = df[(df["rank"] >= 1) & (df["rank"] % 1 == 0) &
            (df["probability"] >= 0) & (df["probability"] <= 1)]
    for col in ("center_x", "center_y", "center_z"):
        df = df[df[col].abs() <= 100000]
    if df.empty:
        trace("P2Rank CSV contained no finite, plausible pocket rows")
        return None
    df = df.sort_values(by=["rank", "probability"], ascending=[True, False], kind="stable").reset_index(drop=True)
    b = df.iloc[0]
    prob = float(b["probability"])
    trace(f"P2Rank pocket rank {int(b['rank'])}: score {float(b['score']):.3f}, "
          f"probability {prob:.3f} ({prob*100:.1f}% confidence)")
    if prob < 0.30:
        trace("WARNING low-confidence pocket (probability < 0.30) - inspect the box")
    return dict(
        center=[float(b["center_x"]), float(b["center_y"]), float(b["center_z"])],
        size=["", "", ""], source="p2rank",
        detail=f"rank {int(b['rank'])}, probability {prob:.3f} ({prob*100:.1f}% confidence)",
        rank=int(b["rank"]), score=float(b["score"]), probability=prob)


# ============================================================
# Run
# ============================================================

clean = os.path.join(pdir, f"{protein}_clean.pdb")
if not os.path.exists(clean):
    print("\nERROR")
    print("Clean protein not found. Run 05_extract_protein.py first.")
    print(clean)
    sys.exit(1)

clean_atoms = parse_pdb_atoms(clean)
chain = discover_chain(clean_atoms)

print("=" * 60)
print("Pocket Selection")
print("=" * 60)
print(f"Protein      : {protein}")
print(f"Structure    : {clean}")
print(f"Chain        : {chain}")
print(f"Ligand floor : {MIN_LIGAND_ATOMS} heavy atoms")
if args.drugs_only:
    print("Mode         : --drugs-only, cofactors and substrates refused")
print()

# Explicit user pocket choices bypass both cocrystal selection and P2Rank.
try:
    result = custom_pocket(clean_atoms)
except ValueError as exc:
    print("\nERROR")
    print(f"Invalid custom pocket: {exc}")
    sys.exit(1)

if result is None:
    print("--- Tier 1: co-crystallised ligand")
    result = tier_cocrystal(clean_atoms, chain)
    if not result:
        print()
        print("--- Tier 2: P2Rank")
        result = tier_p2rank()
else:
    print("--- User-defined pocket")

if not result:
    print("\nERROR")
    print("No pocket could be determined. Pocket_Config.csv left untouched.")
    with open(report_file, "w") as f:
        f.write("=" * 60 + "\n")
        f.write("Pocket Extraction Report - FAILED\n")
        f.write("=" * 60 + "\n\n")
        f.write(f"Protein      : {protein}\n\n")
        f.write("\n".join(TRACE) + "\n")
    sys.exit(1)

trace(f"ACCEPTED: {result['source']}")

cx, cy, cz = result["center"]
sx, sy, sz = result["size"]

write_config(
    Protein=protein, Pocket_Source=result["source"], Pocket_Rank=result["rank"],
    Center_X=cx, Center_Y=cy, Center_Z=cz,
    Size_X=sx, Size_Y=sy, Size_Z=sz,
    Score=result["score"], Probability=result["probability"],
)

with open(report_file, "w") as f:
    f.write("=" * 60 + "\n")
    f.write("Pocket Extraction Report\n")
    f.write("=" * 60 + "\n\n")
    f.write(f"Protein      : {protein}\n")
    f.write(f"Structure    : {clean}\n")
    f.write(f"Chain        : {chain}\n")
    f.write(f"Pocket source: {result['source']}\n")
    f.write(f"Detail       : {result['detail']}\n\n")
    f.write(f"Center X     : {cx:.4f}\n")
    f.write(f"Center Y     : {cy:.4f}\n")
    f.write(f"Center Z     : {cz:.4f}\n\n")
    if sx != "":
        f.write(f"Size X       : {sx:.4f}\n")
        f.write(f"Size Y       : {sy:.4f}\n")
        f.write(f"Size Z       : {sz:.4f}\n\n")
    else:
        f.write("Size         : blank, 09_run_unidock.py default applies\n\n")
    f.write("Selection trace\n")
    f.write("-" * 60 + "\n")
    f.write("\n".join(TRACE) + "\n")

print()
print("=" * 60)
print("Best Binding Pocket")
print("=" * 60)
print(f"Protein      : {protein}")
print(f"Pocket source: {result['source']}")
print(f"Detail       : {result['detail']}")
if result["probability"] != "":
    print(f"P2Rank Confidence : {float(result['probability'])*100:.1f}%")
if result["rank"] != "":
    print(f"Pocket Rank  : {result['rank']}")
print()
print(f"Center X     : {cx:.4f}")
print(f"Center Y     : {cy:.4f}")
print(f"Center Z     : {cz:.4f}")
print()
if sx != "":
    print(f"Box Size     : {sx:.2f}  {sy:.2f}  {sz:.2f}")
else:
    print("Box Size     : default from 09_run_unidock.py")
if result["probability"] != "" and float(result["probability"]) < 0.30:
    print("WARNING: low-confidence pocket (<30% confidence) - inspect the box.")
print()
print("Uni-Dock Box Center")
print(f"{cx:.4f}  {cy:.4f}  {cz:.4f}")
print()
print("Saved:")
print(outfile)
print(report_file)
if os.path.exists(evidence_file) and result["source"] == "cocrystal_ligand":
    print(evidence_file)