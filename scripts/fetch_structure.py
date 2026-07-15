#!/usr/bin/env python3
"""
Resolve a gene / protein name to a 3D protein structure via UniProt, then download it.

Selection policy:
  1. Search UniProt for the name (default organism = human, taxon 9606; reviewed
     Swiss-Prot entries preferred) -> get the UniProt accession.
  2. Among the entry's PDB cross-references, PREFER experimentally-solved
     structures (X-ray / EM / NMR), ranked by resolution (best = smallest number).
     Pick the highest-resolution experimental structure and download it from RCSB.
  3. If NO experimental PDB exists, fall back to the AlphaFold predicted model.

Usage:
    python scripts/fetch_structure.py --gene EGFR
    python scripts/fetch_structure.py --protein_name "carbonic anhydrase 2"
    python scripts/fetch_structure.py --gene TP53 --organism 10090   # mouse
    python scripts/fetch_structure.py --gene EGFR --out egfr.pdb
"""
import argparse
import json
import os
import sys
import urllib.parse
import urllib.request
import urllib.error

UNIPROT_SEARCH = "https://rest.uniprot.org/uniprotkb/search"
UNIPROT_ENTRY = "https://rest.uniprot.org/uniprotkb/{acc}.json"
RCSB_PDB = "https://files.rcsb.org/download/{pdb}.pdb"
ALPHAFOLD_PDB = "https://alphafold.ebi.ac.uk/files/AF-{acc}-F1-model_v4.pdb"

HUMAN_TAXON = "9606"
# methods considered "experimental" (vs. predicted models)
EXPERIMENTAL_METHODS = {"X-ray", "EM", "NMR", "Neutron", "Fiber", "IR spectroscopy"}


def _get_json(url, retries=3):
    last_err = None
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={"Accept": "application/json"})
            with urllib.request.urlopen(req, timeout=60) as r:
                return json.loads(r.read().decode())
        except Exception as e:  # SSL/network hiccups are common; retry
            last_err = e
            if attempt < retries - 1:
                import time
                time.sleep(2 * (attempt + 1))
    raise last_err


def search_uniprot(gene=None, protein_name=None, organism=HUMAN_TAXON):
    """Return (accession, entry_id, description) for the best-matching entry."""
    if gene:
        q = f"gene:{gene} AND organism_id:{organism} AND reviewed:true"
    else:
        q = f'protein_name:"{protein_name}" AND organism_id:{organism} AND reviewed:true'

    params = urllib.parse.urlencode({
        "query": q,
        "format": "json",
        "fields": "accession,id,protein_name,gene_names,organism_name",
        "size": "5",
    })
    data = _get_json(f"{UNIPROT_SEARCH}?{params}")
    results = data.get("results") or []

    # If nothing under reviewed, retry without the reviewed filter
    if not results:
        q2 = q.replace(" AND reviewed:true", "")
        params = urllib.parse.urlencode({
            "query": q2, "format": "json",
            "fields": "accession,id,protein_name,gene_names,organism_name", "size": "5",
        })
        data = _get_json(f"{UNIPROT_SEARCH}?{params}")
        results = data.get("results") or []

    if not results:
        sys.exit(f"ERROR: UniProt 未找到匹配条目 (gene={gene}, protein_name={protein_name}, organism={organism})")

    r = results[0]
    acc = r["primaryAccession"]
    entry_id = r.get("uniProtkbId", "")
    desc = ""
    try:
        desc = r["proteinDescription"]["recommendedName"]["fullName"]["value"]
    except (KeyError, TypeError):
        desc = entry_id
    org = r.get("organism", {}).get("scientificName", "")
    print(f"[uniprot] {acc} ({entry_id})  {desc}  [{org}]")
    return acc, entry_id, desc


def list_pdb_structures(accession):
    """Return list of dicts: {pdb, method, resolution(float or None), chains}."""
    data = _get_json(UNIPROT_ENTRY.format(acc=accession))
    out = []
    for x in data.get("uniProtKBCrossReferences", []):
        if x.get("database") != "PDB":
            continue
        props = {p["key"]: p["value"] for p in x.get("properties", [])}
        res_raw = props.get("Resolution", "")
        res = None
        if res_raw and res_raw[0].isdigit():
            try:
                res = float(res_raw.split()[0])
            except ValueError:
                res = None
        out.append({
            "pdb": x["id"],
            "method": props.get("Method", ""),
            "resolution": res,
            "chains": props.get("Chains", ""),
        })
    return out


def choose_best(structures):
    """
    Prefer experimental structures with a numeric resolution, best (smallest) first.
    Returns the chosen dict, or None if no experimental structure is usable.
    """
    exp = [s for s in structures if s["method"] in EXPERIMENTAL_METHODS]
    with_res = [s for s in exp if s["resolution"] is not None]
    if with_res:
        with_res.sort(key=lambda s: s["resolution"])
        return with_res[0]
    # experimental but no resolution parsed (e.g. NMR) — still prefer over predicted
    if exp:
        return exp[0]
    return None


def download(url, out_path, retries=3):
    last_err = None
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={"Accept": "text/plain"})
            with urllib.request.urlopen(req, timeout=120) as r:
                content = r.read().decode()
            break
        except urllib.error.HTTPError:
            raise  # a real 404 etc — don't retry
        except Exception as e:
            last_err = e
            if attempt < retries - 1:
                import time
                time.sleep(2 * (attempt + 1))
            else:
                raise last_err
    with open(out_path, "w") as f:
        f.write(content)
    return len(content)


def resolve_and_download(gene=None, protein_name=None, organism=HUMAN_TAXON, out=None):
    """Full pipeline. Returns (out_path, info_dict)."""
    acc, entry_id, desc = search_uniprot(gene, protein_name, organism)
    structures = list_pdb_structures(acc)
    print(f"[pdb] UniProt 记录 {len(structures)} 个 PDB 交叉引用")

    best = choose_best(structures)
    label = gene or protein_name or acc

    if best:
        pdb_id = best["pdb"]
        res = best["resolution"]
        res_str = f"{res} Å" if res is not None else "N/A"
        print(f"[choose] 选中实验结构 {pdb_id}  方法={best['method']}  分辨率={res_str}  链={best['chains']}")
        out_path = out or f"{label}_{pdb_id}.pdb"
        try:
            n = download(RCSB_PDB.format(pdb=pdb_id), out_path)
            print(f"[ok] 已下载 {pdb_id} 到 {out_path} ({n} 字符)")
            return out_path, {"source": "experimental", "pdb": pdb_id,
                              "method": best["method"], "resolution": res,
                              "accession": acc, "description": desc}
        except urllib.error.HTTPError as e:
            print(f"[warn] 下载 {pdb_id} 失败 (HTTP {e.code})，回退 AlphaFold")

    # fallback: AlphaFold predicted model
    print(f"[choose] 无可用实验结构，回退 AlphaFold 预测模型 (AF-{acc})")
    out_path = out or f"{label}_AF_{acc}.pdb"
    try:
        n = download(ALPHAFOLD_PDB.format(acc=acc), out_path)
        print(f"[ok] 已下载 AlphaFold 模型到 {out_path} ({n} 字符)")
        print("[note] 这是预测结构，非实验解析，对接结果请谨慎解读")
        return out_path, {"source": "alphafold", "accession": acc, "description": desc}
    except urllib.error.HTTPError as e:
        sys.exit(f"ERROR: AlphaFold 也无此结构 (HTTP {e.code})，无法获取 {label} 的结构")


def main():
    ap = argparse.ArgumentParser(description="Fetch protein structure by gene/protein name via UniProt")
    ap.add_argument("--gene", help="Gene symbol, e.g. EGFR")
    ap.add_argument("--protein_name", help="Protein name, e.g. 'carbonic anhydrase 2'")
    ap.add_argument("--organism", default=HUMAN_TAXON, help="NCBI taxon id (default 9606 = human)")
    ap.add_argument("--out", help="Output PDB path")
    ap.add_argument("--list", action="store_true", help="Only list available PDB structures, don't download")
    args = ap.parse_args()

    if not args.gene and not args.protein_name:
        sys.exit("ERROR: 需要 --gene 或 --protein_name")

    if args.list:
        acc, _, _ = search_uniprot(args.gene, args.protein_name, args.organism)
        structures = list_pdb_structures(acc)
        exp = [s for s in structures if s["method"] in EXPERIMENTAL_METHODS and s["resolution"] is not None]
        exp.sort(key=lambda s: s["resolution"])
        print(f"\n实验结构(按分辨率排序，前 15):")
        for s in exp[:15]:
            print(f"  {s['pdb']}  {s['method']:6s}  {s['resolution']} Å  {s['chains']}")
        return 0

    out_path, info = resolve_and_download(args.gene, args.protein_name, args.organism, args.out)
    print(f"\n[result] 结构文件: {out_path}")
    print(f"[result] 来源: {info}")
    return 0


if __name__ == "__main__":
    main()
