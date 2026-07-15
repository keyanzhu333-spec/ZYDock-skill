#!/usr/bin/env python3
"""
Single protein-ligand docking via the NVIDIA NIM DiffDock API (cloud, no local GPU).

Usage:
    export NVIDIA_API_KEY="nvapi-xxxxxxxx"

    # protein from PDB file, ligand as SMILES
    python scripts/dock_single.py --protein protein.pdb \
        --ligand "CC(=O)Oc1ccccc1C(=O)O" --out_dir results/single/

    # protein from sequence
    python scripts/dock_single.py --protein_sequence "MSKGEEL..." \
        --ligand ligand.sdf --out_dir results/single/
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from diffdock_api import read_protein, read_ligand, dock, save_poses, confidence_level  # noqa: E402
from fetch_structure import resolve_and_download  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description="Cloud DiffDock single docking (NVIDIA API)")
    ap.add_argument("--protein", help="Path to protein PDB file")
    ap.add_argument("--protein_sequence", help="Protein amino-acid sequence (alternative to --protein)")
    ap.add_argument("--gene", help="Gene symbol (auto-fetch structure from UniProt, e.g. EGFR)")
    ap.add_argument("--protein_name", help="Protein name (auto-fetch structure from UniProt)")
    ap.add_argument("--organism", default="9606", help="NCBI taxon id for auto-fetch (default 9606 = human)")
    ap.add_argument("--ligand", required=True, help="Ligand SMILES string OR path to .sdf/.mol2 file")
    ap.add_argument("--out_dir", default="results/single", help="Output directory")
    ap.add_argument("--num_poses", type=int, default=20)
    ap.add_argument("--steps", type=int, default=18)
    ap.add_argument("--time_divisions", type=int, default=20)
    args = ap.parse_args()

    # If gene/protein_name given, auto-fetch a structure from UniProt first.
    if args.gene or args.protein_name:
        if args.protein or args.protein_sequence:
            sys.exit("ERROR: --gene/--protein_name 不能与 --protein/--protein_sequence 同时使用")
        os.makedirs(args.out_dir, exist_ok=True)
        label = args.gene or args.protein_name
        fetched = os.path.join(args.out_dir, f"{label}_structure.pdb")
        print(f"[fetch] 按名称自动获取结构: {label} (organism={args.organism})")
        args.protein, _info = resolve_and_download(
            gene=args.gene, protein_name=args.protein_name,
            organism=args.organism, out=fetched,
        )
        print()

    protein_text, is_seq = read_protein(args.protein, args.protein_sequence)
    ligand_text, ligand_type = read_ligand(args.ligand)

    src = args.protein if args.protein else "sequence"
    print(f"[info] 蛋白: {src}  |  配体: {args.ligand[:60]} ({ligand_type})")
    print(f"[info] 调用 NVIDIA DiffDock API (num_poses={args.num_poses}) ...")

    try:
        data = dock(
            protein_text, ligand_text, ligand_file_type=ligand_type,
            num_poses=args.num_poses, steps=args.steps,
            time_divisions=args.time_divisions, is_sequence=is_seq,
        )
    except RuntimeError as e:
        print(f"[error] {e}")
        if "401" in str(e) or "403" in str(e):
            print("  → API key 无效或过期，请去 build.nvidia.com 重新生成")
        sys.exit(1)

    if data.get("status") != "success":
        print(f"[warn] status={data.get('status')}  details={data.get('details')}")

    ranking = save_poses(data, args.out_dir)
    print(f"[ok] 结果已保存到 {args.out_dir}/  (共 {len(ranking)} 个位姿)")

    if ranking:
        print("[result] 置信度 top 5:")
        for rank, pose, c in ranking[:5]:
            print(f"    rank{rank}: pose#{pose}  conf={c:.3f}  ({confidence_level(c)})")


if __name__ == "__main__":
    main()
