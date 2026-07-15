#!/usr/bin/env python3
"""
Batch / virtual-screening docking via the NVIDIA NIM DiffDock API (cloud, no local GPU).

Reads a CSV of protein-ligand pairs, docks each one, saves per-complex results,
and writes a summary.csv ranked by best confidence.

CSV columns (see assets/batch_template.csv):
    complex_name,protein_path,ligand_smiles,protein_sequence

Usage:
    export NVIDIA_API_KEY="nvapi-xxxxxxxx"
    python scripts/dock_batch.py --csv my_input.csv --out_dir results/batch/ --num_poses 10
"""
import argparse
import csv
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from diffdock_api import read_protein, read_ligand, dock, save_poses, confidence_level  # noqa: E402
from fetch_structure import resolve_and_download  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description="Cloud DiffDock batch docking (NVIDIA API)")
    ap.add_argument("--csv", required=True, help="Input CSV file")
    ap.add_argument("--out_dir", default="results/batch", help="Output directory")
    ap.add_argument("--num_poses", type=int, default=10)
    ap.add_argument("--steps", type=int, default=18)
    ap.add_argument("--time_divisions", type=int, default=20)
    args = ap.parse_args()

    if not os.path.isfile(args.csv):
        sys.exit(f"ERROR: CSV not found: {args.csv}")

    os.makedirs(args.out_dir, exist_ok=True)

    with open(args.csv, newline="") as f:
        rows = list(csv.DictReader(f))
    if not rows:
        sys.exit("ERROR: CSV has no data rows.")

    print(f"[info] 共 {len(rows)} 个化合物待对接\n")
    summary = []
    gene_cache = {}  # (gene/name, organism) -> fetched pdb path, avoid re-downloading

    for idx, row in enumerate(rows, 1):
        name = (row.get("complex_name") or f"complex_{idx}").strip()
        protein_path = (row.get("protein_path") or "").strip()
        smiles = (row.get("ligand_smiles") or row.get("ligand_description") or "").strip()
        sequence = (row.get("protein_sequence") or "").strip()
        gene = (row.get("gene") or "").strip()
        protein_name = (row.get("protein_name") or "").strip()
        organism = (row.get("organism") or "9606").strip()

        print(f"[{idx}/{len(rows)}] {name} ...", end=" ", flush=True)

        if not smiles:
            print("跳过 (无配体)")
            summary.append((name, None, "no_ligand"))
            continue

        # auto-fetch structure by gene/protein_name if no explicit protein given
        if (gene or protein_name) and not protein_path and not sequence:
            cache_key = (gene or protein_name, organism)
            if cache_key in gene_cache:
                protein_path = gene_cache[cache_key]
            else:
                try:
                    fetched = os.path.join(args.out_dir, f"_struct_{gene or protein_name}.pdb")
                    protein_path, _ = resolve_and_download(
                        gene=gene or None, protein_name=protein_name or None,
                        organism=organism, out=fetched,
                    )
                    gene_cache[cache_key] = protein_path
                except SystemExit as e:
                    print(f"结构获取失败: {str(e)[:80]}")
                    summary.append((name, None, "no_structure"))
                    continue

        try:
            protein_text, is_seq = read_protein(
                protein_path or None, sequence or None
            )
            ligand_text, ligand_type = read_ligand(smiles)
            data = dock(
                protein_text, ligand_text, ligand_file_type=ligand_type,
                num_poses=args.num_poses, steps=args.steps,
                time_divisions=args.time_divisions, is_sequence=is_seq,
            )
        except (RuntimeError, SystemExit) as e:
            print(f"失败: {str(e)[:80]}")
            summary.append((name, None, "error"))
            continue

        sub = os.path.join(args.out_dir, name)
        ranking = save_poses(data, sub)
        best = ranking[0][2] if ranking else None
        if best is not None:
            print(f"完成  best_conf={best:.3f} ({confidence_level(best)})")
        else:
            print("完成 (无位姿)")
        summary.append((name, best, "ok"))

    # write summary ranked by best confidence
    summary_path = os.path.join(args.out_dir, "summary.csv")
    ranked = sorted(
        summary,
        key=lambda x: (x[1] if x[1] is not None else -1e9),
        reverse=True,
    )
    with open(summary_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["rank", "complex_name", "best_confidence", "level", "status"])
        for r, (name, best, st) in enumerate(ranked, 1):
            lvl = confidence_level(best) if best is not None else "N/A"
            bc = f"{best:.4f}" if best is not None else ""
            w.writerow([r, name, bc, lvl, st])

    print(f"\n[ok] 汇总已保存: {summary_path}")
    print("[result] 按最优置信度排名:")
    for r, (name, best, st) in enumerate(ranked[:10], 1):
        bc = f"{best:.3f}" if best is not None else "N/A"
        print(f"    {r}. {name:20s} best_conf={bc}")


if __name__ == "__main__":
    main()
