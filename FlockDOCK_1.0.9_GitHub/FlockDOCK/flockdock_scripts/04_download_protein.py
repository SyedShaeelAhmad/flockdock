"""
04_download_protein.py

Download the best receptor structure for a gene.

Order of preference:
    1. Experimental human structure from the RCSB PDB that has a real
       co-crystallised ligand bound  <-- this is the new primary criterion.
       A structure with a drug or cofactor in the pocket tells us where our own
       small molecules should be docked, so ligand-bearing entries outrank
       higher-resolution apo entries.
    2. Experimental structure with no usable ligand (apo), best resolution.
    3. AlphaFold predicted model -> ONLY if no experimental structure exists.

Ranking within group 1:
    drug-like ligand before cofactor/substrate, then X-ray before cryo-EM,
    then best resolution, then newest release.

Every downloaded structure is verified after download by counting heavy atoms
of each candidate ligand in the file itself, so a component that RCSB lists but
that is absent, disordered or too small never survives.

The best structure is copied to  04_proteins/<GENE>/<GENE>.pdb  (canonical).
Two side files are written for step 08:
    <GENE>_metadata.csv          ranked candidates with their ligands
    <GENE>_ligand_candidates.csv verified ligand instances, per structure

Usage:
    python 04_download_protein.py GENE [UNIPROT] [--drugs-only]

    GENE     required   e.g. ATP2A2
    UNIPROT  optional   e.g. P16615  (skips the gene -> accession lookup;
                                      pass it, gene-symbol lookup can pick a
                                      pseudogene paralogue)
    --drugs-only        refuse cofactors and substrates such as ATP or G6P,
                        accepting only drug-like ligands. Off by default,
                        because for several enzymes the nucleotide IS the only
                        thing that marks the catalytic site.
"""

import sys
import os
import time
import shutil
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import requests
import pandas as pd

from ligand_filter import classify, ligand_instances, parse_pdb_atoms

# ============================================================
# Settings
# ============================================================

MAX_STRUCTURES = 5
MAX_METADATA_ENTRIES = 60        # cap the per-entry metadata calls

UNIPROT_SEARCH_URL = "https://rest.uniprot.org/uniprotkb/search"
SEARCH_URL         = "https://search.rcsb.org/rcsbsearch/v2/query"
ENTRY_URL          = "https://data.rcsb.org/rest/v1/core/entry/{}"
DOWNLOAD_URL       = "https://files.rcsb.org/download/{}.pdb"

ALPHAFOLD_API      = "https://alphafold.ebi.ac.uk/api/prediction/{acc}"
ALPHAFOLD_FILE     = "https://alphafold.ebi.ac.uk/files/AF-{acc}-F1-model_v{ver}.pdb"
ALPHAFOLD_VERSIONS = [6, 5, 4]

PROJECT = Path(os.environ.get("UNIDOCK_PROJECT", os.environ.get("ONEDOCK_PROJECT", os.getcwd())))

DRUGS_ONLY = False               # set from the command line


# ============================================================
# Robust request
# ============================================================

def safe_request(method, url, **kwargs):
    for attempt in range(5):
        try:
            response = requests.request(method, url, **kwargs)
            response.raise_for_status()
            return response
        except requests.exceptions.RequestException:
            if attempt == 4:
                raise
            print(f"   Connection failed (attempt {attempt + 1}/5). Retrying...")
            time.sleep(5)


def result_set(response):
    """
    Safely pull 'result_set' from an RCSB search response.

    The RCSB Search API returns HTTP 204 with an EMPTY body when a query
    matches zero structures. raise_for_status() does not flag 204, so calling
    .json() on the empty body throws. Guard against that here.
    """
    if response is None:
        return []
    if response.status_code == 204 or not response.text.strip():
        return []
    try:
        return response.json().get("result_set", [])
    except ValueError:
        return []


# ============================================================
# Gene -> UniProt accession
# ============================================================

def resolve_uniprot(gene):
    params = {
        "query": f"gene_exact:{gene} AND organism_id:9606 AND reviewed:true",
        "fields": "accession",
        "format": "json",
        "size": 1,
    }
    try:
        response = safe_request("GET", UNIPROT_SEARCH_URL, params=params,
                                timeout=30)
    except requests.exceptions.RequestException:
        return None
    if response is None or not response.text.strip():
        return None
    try:
        results = response.json().get("results", [])
    except ValueError:
        return None
    return results[0].get("primaryAccession") if results else None


# ============================================================
# PDB search
# ============================================================

def search_pdb_by_uniprot(acc):
    query = {
        "query": {
            "type": "group", "logical_operator": "and",
            "nodes": [
                {"type": "terminal", "service": "text", "parameters": {
                    "attribute": "rcsb_polymer_entity_container_identifiers."
                                 "reference_sequence_identifiers.database_name",
                    "operator": "exact_match", "value": "UniProt"}},
                {"type": "terminal", "service": "text", "parameters": {
                    "attribute": "rcsb_polymer_entity_container_identifiers."
                                 "reference_sequence_identifiers."
                                 "database_accession",
                    "operator": "exact_match", "value": acc}},
            ],
        },
        "return_type": "entry",
        "request_options": {"return_all_hits": True},
    }
    return result_set(safe_request("POST", SEARCH_URL, json=query, timeout=30))


def search_pdb_by_gene(gene):
    query = {
        "query": {
            "type": "group", "logical_operator": "and",
            "nodes": [
                {"type": "terminal", "service": "text", "parameters": {
                    "attribute": "rcsb_entity_source_organism."
                                 "rcsb_gene_name.value",
                    "operator": "exact_match", "value": gene}},
                {"type": "terminal", "service": "text", "parameters": {
                    "attribute": "rcsb_entity_source_organism."
                                 "ncbi_scientific_name",
                    "operator": "exact_match", "value": "Homo sapiens"}},
            ],
        },
        "return_type": "entry",
        "request_options": {"return_all_hits": True},
    }
    return result_set(safe_request("POST", SEARCH_URL, json=query, timeout=30))


# ============================================================
# Metadata + ligand classification + ranking
# ============================================================

def collect_metadata(hits):
    rows = []
    for hit in hits[:MAX_METADATA_ENTRIES]:
        pdb = hit["identifier"]
        try:
            info = safe_request("GET", ENTRY_URL.format(pdb), timeout=30).json()
        except (requests.exceptions.RequestException, ValueError):
            continue

        if info.get("rcsb_accession_info", {}).get("status_code", "REL") != "REL":
            continue

        entry_info = info.get("rcsb_entry_info", {})
        try:
            resolution = entry_info["resolution_combined"][0]
        except (KeyError, IndexError):
            resolution = None

        # Component IDs of everything non-polymer that is actually bound.
        bound = entry_info.get("nonpolymer_bound_components") or []
        drugs, cofactors, rejected = [], [], []
        for comp in bound:
            verdict, _ = classify(comp, None, DRUGS_ONLY)
            if verdict == "drug":
                drugs.append(comp)
            elif verdict == "cofactor":
                cofactors.append(comp)
            else:
                rejected.append(comp)

        organisms = set()
        for src in info.get("rcsb_entity_source_organism", []) or []:
            name = src.get("ncbi_scientific_name")
            if name:
                organisms.add(name)

        rows.append({
            "PDB": pdb,
            "Resolution": resolution,
            "Method": info.get("exptl", [{}])[0].get("method", ""),
            "Release_Date": info.get("rcsb_accession_info", {})
                                .get("initial_release_date", ""),
            "Chains": entry_info.get("polymer_entity_count", None),
            "Bound_Components": ";".join(bound),
            "Drug_Ligands": ";".join(drugs),
            "Cofactor_Ligands": ";".join(cofactors),
            "Ignored_Components": ";".join(rejected),
            "Ligand_Class": "drug" if drugs else ("cofactor" if cofactors
                                                  else "apo"),
            "Organisms": ";".join(sorted(organisms)),
            "Title": info.get("struct", {}).get("title", ""),
        })
    return pd.DataFrame(rows)


LIGAND_PRIORITY = {"drug": 0, "cofactor": 1, "apo": 2}


def rank_experimental(df):
    """
    Ligand first, then X-ray over cryo-EM, then resolution, then newest.

    Ligand-bearing entries deliberately outrank a better-resolution apo entry:
    an apo structure gives no information about where the pocket is, which is
    how a 3-atom thiocyanate ends up defining a docking box.
    """
    if df.empty:
        return df
    df = df[df["Method"].str.contains("X-RAY|ELECTRON", case=False,
                                      na=False)].copy()
    if df.empty:
        return df

    df["_lig"] = df["Ligand_Class"].map(LIGAND_PRIORITY).fillna(2).astype(int)
    df["_method"] = df["Method"].str.contains("X-RAY", case=False, na=False) \
                               .map({True: 0, False: 1})
    df["_res"] = pd.to_numeric(df["Resolution"], errors="coerce").fillna(9999)
    df = (df.sort_values(["_lig", "_method", "_res", "Release_Date"],
                         ascending=[True, True, True, False])
            .drop(columns=["_lig", "_method", "_res"])
            .reset_index(drop=True))
    return df


# ============================================================
# AlphaFold fallback
# ============================================================

def mean_plddt(pdb_text):
    vals = []
    for line in pdb_text.splitlines():
        if line.startswith("ATOM") and line[12:16].strip() == "CA":
            try:
                vals.append(float(line[60:66]))
            except ValueError:
                pass
    return round(sum(vals) / len(vals), 1) if vals else None


def fetch_alphafold(acc, out_dir, gene):
    pdb_url = version = None
    try:
        r = safe_request("GET", ALPHAFOLD_API.format(acc=acc), timeout=30)
        meta = r.json() if (r is not None and r.text.strip()) else []
        if meta:
            pdb_url = meta[0].get("pdbUrl")
            version = meta[0].get("latestVersion")
    except (requests.exceptions.RequestException, ValueError):
        pdb_url = None

    if pdb_url:
        try:
            pr = requests.get(pdb_url, timeout=60)
            if pr.status_code == 200 and len(pr.text) > 1000:
                af = out_dir / f"{gene}_AF_{acc}.pdb"
                af.write_text(pr.text)
                return af, version or "current", mean_plddt(pr.text)
        except requests.exceptions.RequestException:
            pass

    for ver in ALPHAFOLD_VERSIONS:
        try:
            r = requests.get(ALPHAFOLD_FILE.format(acc=acc, ver=ver), timeout=60)
        except requests.exceptions.RequestException:
            continue
        if r.status_code == 200 and len(r.text) > 1000:
            af = out_dir / f"{gene}_AF_{acc}_v{ver}.pdb"
            af.write_text(r.text)
            return af, ver, mean_plddt(r.text)
    return None, None, None


# ============================================================
# Main
# ============================================================

def main():
    global DRUGS_ONLY

    argv = [a for a in sys.argv[1:] if a != "--drugs-only"]
    DRUGS_ONLY = "--drugs-only" in sys.argv

    if len(argv) not in (1, 2):
        print("\nUsage:")
        print("  python 04_download_protein.py GENE [UNIPROT] [--drugs-only]\n")
        sys.exit(1)

    gene = argv[0].upper()
    uniprot = argv[1].upper() if len(argv) == 2 else None

    out_dir = PROJECT / "04_proteins" / gene
    canonical = out_dir / f"{gene}.pdb"
    cache_file = out_dir / f"{gene}_structure_cache.json"
    out_dir.mkdir(parents=True, exist_ok=True)

    # ---- 1. Resolve identity before trusting a cached structure ----
    if uniprot is None:
        print(f"\nResolving UniProt accession for {gene}...")
        uniprot = resolve_uniprot(gene)
        print(f"UniProt : {uniprot if uniprot else 'not found'}")
        print("NOTE: gene-symbol lookup can return a pseudogene paralogue. "
              "Pass the accession explicitly where you know it.")

    def valid_pdb(path):
        try:
            if not path.is_file() or path.stat().st_size < 1000:
                return False
            with path.open("r", encoding="utf-8", errors="replace") as handle:
                atom_count = sum(1 for line in handle if line.startswith(("ATOM  ", "HETATM")))
            return atom_count >= 20
        except OSError:
            return False

    def file_sha256(path):
        digest = hashlib.sha256()
        with path.open("rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(chunk)
        return digest.hexdigest()

    # If UniProt lookup is temporarily unavailable, prefer a previously
    # verified cache identity over falling back to an ambiguous gene-symbol PDB search.
    if uniprot is None and valid_pdb(canonical) and cache_file.is_file():
        try:
            cached_identity = json.loads(cache_file.read_text(encoding="utf-8"))
            cached_selected = out_dir / str(cached_identity.get("selected_structure", ""))
            if (cached_identity.get("cache_format_version") == 1
                    and cached_identity.get("uniprot")
                    and cached_identity.get("canonical_sha256") == file_sha256(canonical)
                    and valid_pdb(cached_selected)
                    and cached_identity.get("selected_structure_sha256") == file_sha256(cached_selected)):
                uniprot = cached_identity.get("uniprot")
                print(f"UniProt lookup unavailable; retaining verified cached accession {uniprot}.")
        except Exception:
            pass

    # A canonical PDB by itself is not proof that the cache matches the current
    # UniProt mapping. Only skip network discovery when an integrity-checked
    # cache manifest records the same accession and current structure digest.
    if valid_pdb(canonical) and cache_file.is_file():
        try:
            cache = json.loads(cache_file.read_text(encoding="utf-8"))
            selected_cached = out_dir / str(cache.get("selected_structure", ""))
            selected_cached_valid = (valid_pdb(selected_cached)
                                     and cache.get("selected_structure_sha256") == file_sha256(selected_cached))
            if (cache.get("cache_format_version") == 1
                    and cache.get("uniprot") == uniprot
                    and cache.get("canonical_sha256") == file_sha256(canonical)
                    and selected_cached_valid):
                print(f"\n{gene} structure cache is valid ({canonical.name}; UniProt {uniprot or 'unresolved'}).")
                print(f"Source: {cache.get('source', 'unknown')} | selected: {cache.get('selected_structure', 'unknown')}")
                sys.exit(0)
        except Exception as exc:
            print(f"\nStructure cache manifest is invalid; refreshing: {exc}")

    if canonical.exists():
        print(f"\nExisting structure for {gene} has no matching verified cache; refreshing candidates.")

    # Stale ranked/AlphaFold candidates can otherwise be selected by step 05
    # after a new structure lookup. Preserve the canonical file until a new
    # candidate succeeds, but remove outdated generated candidate artefacts.
    stale_patterns = [f"{gene}_rank*.pdb", f"{gene}_AF_*.pdb", f"{gene}_ligand.pdb",
                      f"{gene}_ligand_box.csv", f"{gene}_metadata.csv", f"{gene}_ligand_candidates.csv"]
    for pattern in stale_patterns:
        for stale_path in out_dir.glob(pattern):
            try:
                stale_path.unlink()
            except OSError:
                pass

    # ---- 2. experimental search ----
    print(f"\nSearching PDB for {gene}...")
    hits = search_pdb_by_uniprot(uniprot) if uniprot else search_pdb_by_gene(gene)
    print(f"Experimental entries found : {len(hits)}")
    if len(hits) > MAX_METADATA_ENTRIES:
        print(f"Examining the first {MAX_METADATA_ENTRIES} for ligand content.")

    df = rank_experimental(collect_metadata(hits)) if hits else pd.DataFrame()

    best_structure = None
    source = None
    saved = 0
    ligand_rows = []

    # ---- 3. download, preferring ligand-bearing entries ----
    if not df.empty:
        source = "Experimental"
        counts = df["Ligand_Class"].value_counts().to_dict()
        print()
        print("=" * 70)
        print("Ligand content across candidate structures")
        print("=" * 70)
        print(f"with drug-like ligand : {counts.get('drug', 0)}")
        print(f"with cofactor only    : {counts.get('cofactor', 0)}")
        print(f"apo (no usable ligand): {counts.get('apo', 0)}")
        if DRUGS_ONLY:
            print("--drugs-only is set: cofactors and substrates are refused.")

        print()
        print("=" * 70)
        print("Top candidate structures (ligand first, then resolution)")
        print("=" * 70)
        show = ["PDB", "Resolution", "Method", "Ligand_Class", "Drug_Ligands",
                "Cofactor_Ligands", "Title"]
        print(df.head(10)[show].to_string(index=False))
        print()

        out = df.copy()
        out.insert(0, "Source", source)
        out.to_csv(out_dir / f"{gene}_metadata.csv", index=False)

        attempt = 0
        for _, row in df.iterrows():
            if saved >= MAX_STRUCTURES:
                break
            attempt += 1
            pdb = row["PDB"]
            print(f"[Attempt {attempt}] {pdb}  ({row['Ligand_Class']})")
            try:
                download = requests.get(DOWNLOAD_URL.format(pdb), timeout=30)
            except requests.exceptions.RequestException as e:
                print(f"   Download failed ({e}). Skipping.")
                continue
            if download.status_code != 200 or len(download.text) < 1000:
                print("   PDB file unavailable / invalid (may be mmCIF-only). "
                      "Skipping.")
                continue

            # verify the ligand really is there, in the file, at usable size
            atoms = parse_pdb_atoms(download.text, from_path=False)
            insts = ligand_instances(atoms, DRUGS_ONLY)
            usable = [i for i in insts if i["verdict"] in ("drug", "cofactor")]

            outfile = out_dir / f"{gene}_rank{saved + 1}_{pdb}.pdb"
            outfile.write_text(download.text)
            saved += 1

            if usable:
                names = ", ".join(f"{i['comp_id']}({i['heavy_atoms']})"
                                  for i in usable[:6])
                print(f"   Saved {outfile.name} ({row['Resolution']} A)  "
                      f"verified ligands: {names}")
            else:
                dropped = ", ".join(f"{i['comp_id']}({i['heavy_atoms']}, "
                                    f"{i['reason']})" for i in insts[:6])
                print(f"   Saved {outfile.name} ({row['Resolution']} A)  "
                      f"NO usable ligand"
                      + (f"; ignored: {dropped}" if dropped else ""))

            for i in usable:
                ligand_rows.append(dict(
                    Gene=gene, UniProt=uniprot, PDB=pdb, File=outfile.name,
                    Rank=saved, Comp_ID=i["comp_id"], Verdict=i["verdict"],
                    Heavy_Atoms=i["heavy_atoms"], Ligand_Chain=i["chain"],
                    Ligand_Resseq=i["resseq"]))

            # canonical = first structure that carries a verified ligand
            if best_structure is None and usable:
                best_structure = outfile

        # no downloaded structure had a ligand -> fall back to the best apo one
        if best_structure is None and saved > 0:
            first = sorted(out_dir.glob(f"{gene}_rank1_*.pdb"))
            if first:
                best_structure = first[0]
                print("\nNo downloaded structure carried a usable ligand. "
                      "Using the best-ranked apo structure; step 08 will fall "
                      "through to P2Rank.")

    if ligand_rows:
        pd.DataFrame(ligand_rows).to_csv(
            out_dir / f"{gene}_ligand_candidates.csv", index=False)

    # ---- 4. AlphaFold fallback ----
    if best_structure is None:
        print("\nNo usable experimental structure.")
        if uniprot:
            print(f"Fetching AlphaFold model for {uniprot}...")
            af_file, ver, plddt = fetch_alphafold(uniprot, out_dir, gene)
            if af_file is not None:
                best_structure = af_file
                source = "AlphaFold"
                print(f"   Saved {af_file.name}  (model v{ver}, "
                      f"mean CA pLDDT = {plddt})")
                print("   No co-crystal ligand exists for this target, so "
                      "step 08 will use P2Rank.")
                pd.DataFrame([{
                    "Source": "AlphaFold", "Model": af_file.name,
                    "UniProt": uniprot, "Version": f"v{ver}",
                    "Mean_pLDDT_CA": plddt, "Ligand_Class": "apo",
                }]).to_csv(out_dir / f"{gene}_metadata.csv", index=False)
            else:
                print("   AlphaFold model not available.")
        else:
            print("   No UniProt accession - cannot fetch AlphaFold model.")

    # ---- 5. finalise ----
    if best_structure is None:
        print("\n" + "=" * 70)
        print("ERROR")
        print("=" * 70)
        print(f"No experimental or AlphaFold structure available for {gene}.")
        sys.exit(1)

    shutil.copyfile(best_structure, canonical)
    cache_payload = {
        "cache_format_version": 1,
        "gene": gene,
        "uniprot": uniprot,
        "source": source,
        "selected_structure": best_structure.name,
        "selected_structure_sha256": file_sha256(best_structure),
        "canonical_sha256": file_sha256(canonical),
        "created_utc": datetime.now(timezone.utc).isoformat(),
    }
    cache_tmp = cache_file.with_suffix(".json.tmp")
    cache_tmp.write_text(json.dumps(cache_payload, indent=2), encoding="utf-8")
    os.replace(cache_tmp, cache_file)

    print("\n" + "=" * 70)
    print("Download Complete")
    print("=" * 70)
    print(f"Gene            : {gene}")
    print(f"UniProt         : {uniprot if uniprot else 'n/a'}")
    print(f"Source          : {source}")
    if source == "Experimental":
        print(f"Structures Saved: {saved}")
        print(f"Ligand-bearing  : "
              f"{len({r['File'] for r in ligand_rows})} of {saved}")
    print(f"Canonical File  : {canonical}")
    print(f"Best Structure  : {best_structure.name}")
    print("\nDone.")


if __name__ == "__main__":
    main()
