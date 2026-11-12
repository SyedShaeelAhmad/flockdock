"""
ligand_filter.py

Shared vocabulary for deciding whether a bound HETATM group is a real ligand
that marks a druggable pocket, or a crystallisation artefact that does not.

Used by 04_download_protein.py (to rank PDB entries) and by
08_extract_pocket.py (to build the docking box). Keeping it in one file means
the two steps can never disagree about what counts as a ligand.

Place next to the other scripts in 09_scripts/.
"""

MIN_LIGAND_ATOMS = 8      # heavy atoms; rejects SCN (3), SO4 (5), GOL (6)
MAX_LIGAND_ATOMS = 120    # rejects lipids, polymers, bound peptides

# Standard residues, so polymer records never look like ligands.
AA3 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V", "SEC": "U", "PYL": "O", "MSE": "M",
    "CSO": "C", "CSD": "C", "OCS": "C", "CME": "C", "HIC": "H", "MLY": "K",
}
NUCLEOTIDE_RES = {"A", "U", "G", "C", "T", "DA", "DT", "DG", "DC", "DU", "N"}

# Solvent, ions, buffers, cryoprotectants, detergents, lipids.
ADDITIVES = {
    "HOH", "DOD", "UNX", "UNL", "UNK",
    "NA", "K", "LI", "RB", "CS", "MG", "CA", "SR", "BA", "ZN", "MN", "FE",
    "FE2", "FES", "CU", "CU1", "CD", "NI", "CO", "HG", "AU", "AG", "PT", "PB",
    "CL", "BR", "IOD", "F", "FLO", "XE", "KR", "ARF",
    "SO4", "SUL", "PO4", "2HP", "NO3", "NO2", "CO3", "BCT", "AZI", "CN", "SCN",
    "WO4", "MOO", "VO4", "ALF", "BEF", "AF3", "MGF",
    "ACT", "ACY", "ACE", "FMT", "OXL", "MLI", "MLA", "CIT", "FLC", "TLA", "TAR",
    "SIN", "GLV", "PYR", "LAC", "SCN", "THJ",
    "EDO", "GOL", "PGO", "PEG", "PGE", "P6G", "P33", "1PE", "2PE", "7PE",
    "12P", "15P", "MPD", "MRD", "BU3", "BU1", "IPA", "EOH", "MOH", "DMS",
    "DMF", "ACN", "URE", "GAI", "NHE",
    "TRS", "EPE", "MES", "PIN", "HEZ", "IMD", "NH4", "BME", "DTT", "DTU",
    "DTV", "TCE", "MPO", "POP",
    "BOG", "LDA", "C8E", "OCT", "SDS", "LMT", "DDQ", "PCW", "LPP", "PX4",
    "SGM", "MYR", "PLM", "OLA", "OLC", "STE", "D10", "D12", "HP6", "UND",
    "CPS", "CHD", "BNG", "9BQ",
    # N-/O-linked glycans: covalently attached surface sugars, never the
    # pocket we want. Free substrate sugars (GLC, BGC, G6P, F6P) stay in
    # COFACTORS below.
    "NAG", "NDG", "BMA", "MAN", "GAL", "FUC", "XYP", "SIA", "A2G", "GLA",
    "NGA", "RAM", "XYS", "BGL",
}

# Cofactors, nucleotides and substrates. Not drugs, but they sit in the
# catalytic site and therefore define exactly the pocket we want to dock into.
# Excluded when --drugs-only is passed.
COFACTORS = {
    "ATP", "ADP", "AMP", "ANP", "ACP", "AGS", "APC", "ATR", "ADX",
    "GTP", "GDP", "GMP", "GNP", "GSP", "G2P",
    "CTP", "CDP", "UTP", "UDP", "UMP", "TTP", "TMP", "PAP", "A3P",
    "NAD", "NAI", "NAJ", "NAP", "NDP", "NAX",
    "FAD", "FDA", "FMN", "FNR",
    "SAH", "SAM", "SFG", "COA", "ACO", "COO", "CMX",
    "GSH", "GDS", "GTS", "GSN", "GBX",
    "HEM", "HEC", "HEB", "SRM", "PLP", "PMP", "TPP", "TDP", "BTN", "B12",
    "GLC", "BGC", "G6P", "G6Q", "F6P", "FBP", "G1P", "16G", "AGL",
    "PEP", "3PG", "2PG", "OAA", "AKG", "SUC", "FUM", "MAL",
}


def classify(comp_id, heavy_atoms=None, drugs_only=False):
    """
    Return (verdict, reason).

    verdict is one of: 'drug', 'cofactor', 'reject'.
    heavy_atoms may be None when only the component ID is known, e.g. when
    ranking PDB entries before download.
    """
    cid = (comp_id or "").strip().upper()
    if not cid:
        return "reject", "no component id"
    if cid in AA3 or cid in NUCLEOTIDE_RES:
        return "reject", "standard polymer residue"
    if cid in ADDITIVES:
        return "reject", "crystallisation additive"
    if heavy_atoms is not None:
        if heavy_atoms < MIN_LIGAND_ATOMS:
            return "reject", f"{heavy_atoms} heavy atoms, below {MIN_LIGAND_ATOMS}"
        if heavy_atoms > MAX_LIGAND_ATOMS:
            return "reject", f"{heavy_atoms} heavy atoms, above {MAX_LIGAND_ATOMS}"
    if cid in COFACTORS:
        if drugs_only:
            return "reject", "cofactor or substrate, excluded by --drugs-only"
        return "cofactor", "cofactor or substrate in the catalytic site"
    return "drug", "drug-like ligand"


def is_usable(comp_id, heavy_atoms=None, drugs_only=False):
    return classify(comp_id, heavy_atoms, drugs_only)[0] in ("drug", "cofactor")


def parse_pdb_atoms(source, from_path=True):
    """Heavy atoms only. Returns list of dicts."""
    lines = (open(source).readlines() if from_path
             else source.splitlines(keepends=True))
    atoms = []
    for line in lines:
        if not line.startswith(("ATOM", "HETATM")) or len(line) < 54:
            continue
        element = line[76:78].strip().upper() if len(line) >= 78 else ""
        name = line[12:16].strip()
        if element in ("H", "D"):
            continue
        if not element and name[:1] in ("H", "D"):
            continue
        try:
            atoms.append(dict(
                het=line.startswith("HETATM"),
                name=name,
                resname=line[17:20].strip().upper(),
                chain=line[21].strip(),
                resseq=int(line[22:26]),
                x=float(line[30:38]), y=float(line[38:46]),
                z=float(line[46:54]),
            ))
        except ValueError:
            continue
    return atoms


def ligand_instances(atoms, drugs_only=False):
    """
    Group HETATM records into ligand instances and keep the usable ones.

    Returns list of dicts: comp_id, chain, resseq, heavy_atoms, verdict, atoms.
    """
    groups = {}
    for a in atoms:
        if not a["het"]:
            continue
        groups.setdefault((a["chain"], a["resseq"], a["resname"]), []).append(a)

    out = []
    for (chain, resseq, comp), lat in groups.items():
        verdict, reason = classify(comp, len(lat), drugs_only)
        out.append(dict(comp_id=comp, chain=chain, resseq=resseq,
                        heavy_atoms=len(lat), verdict=verdict, reason=reason,
                        atoms=lat))
    out.sort(key=lambda d: (0 if d["verdict"] == "drug" else
                            1 if d["verdict"] == "cofactor" else 2,
                            -d["heavy_atoms"]))
    return out
