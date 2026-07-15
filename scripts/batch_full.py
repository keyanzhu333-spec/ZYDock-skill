#!/usr/bin/env python3
"""
End-to-end BATCH pipeline: one CSV -> for each compound run the full workflow
(docking -> PLIP interaction analysis -> two PyMOL figures -> bilingual report),
then aggregate everything into a master summary table.

Each compound gets its own subfolder under --out_dir:
    out_dir/
        <compound1>/
            docking/       (poses, confidence, structure)
            interactions/  (complex, PLIP report, 2 figures)
            report/        (report_CN.md, report_EN.md, interactions_summary.csv)
        <compound2>/
            ...
        batch_summary.csv  (all compounds ranked by best confidence + interaction counts)

CSV columns (see assets/batch_template.csv):
    complex_name, gene, protein_path, ligand_smiles, protein_sequence, organism, accession

Usage:
    export NVIDIA_API_KEY=...   # or use set_key.py first
    python scripts/batch_full.py --csv compounds.csv --out_dir results/screen/ --num_poses 10

    # Skip figures/report for a fast docking-only screen:
    python scripts/batch_full.py --csv compounds.csv --out_dir results/screen/ --dock_only
"""
import argparse
import csv
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from diffdock_api import read_protein, read_ligand, dock, save_poses, confidence_level  # noqa: E402
from fetch_structure import resolve_and_download  # noqa: E402


def _safe_name(name):
    """Make a filesystem-safe folder name."""
    return "".join(c if c.isalnum() or c in "-_" else "_" for c in name)


def run_one(row, idx, args, struct_cache):
    """
    Run the full pipeline for a single compound row.
    Returns a summary dict.
    """
    name = _safe_name((row.get("complex_name") or f"complex_{idx}").strip())
    protein_path = (row.get("protein_path") or "").strip()
    smiles = (row.get("ligand_smiles") or row.get("ligand_description") or "").strip()
    sequence = (row.get("protein_sequence") or "").strip()
    gene = (row.get("gene") or "").strip()
    protein_name = (row.get("protein_name") or "").strip()
    organism = (row.get("organism") or "9606").strip()
    accession = (row.get("accession") or "N/A").strip()

    result = {"name": name, "gene": gene or protein_name, "smiles": smiles,
              "best_conf": None, "n_interactions": None, "status": "ok"}

    if not smiles:
        result["status"] = "no_ligand"
        return result

    comp_dir = os.path.join(args.out_dir, name)
    docking_dir = os.path.join(comp_dir, "docking")
    interactions_dir = os.path.join(comp_dir, "interactions")
    os.makedirs(docking_dir, exist_ok=True)

    # ---- 1. Resolve protein structure (cache by gene to avoid re-download) ----
    label = gene or protein_name
    if (gene or protein_name) and not protein_path and not sequence:
        cache_key = (label, organism)
        if cache_key in struct_cache:
            protein_path = struct_cache[cache_key]
        else:
            try:
                fetched = os.path.join(docking_dir, f"{label}_structure.pdb")
                protein_path, _ = resolve_and_download(
                    gene=gene or None, protein_name=protein_name or None,
                    organism=organism, out=fetched,
                )
                struct_cache[cache_key] = protein_path
            except (SystemExit, Exception) as e:
                print(f"    结构获取失败: {str(e)[:100]}")
                result["status"] = "no_structure"
                return result
    elif gene or protein_name:
        # structure will be fetched; but if protein_path/sequence given, use them
        pass

    # If a gene structure was cached but this compound's docking dir needs the file,
    # ensure the structure file lives in this compound's docking dir for the report.
    if protein_path and gene and os.path.isfile(protein_path):
        local_struct = os.path.join(docking_dir, f"{label}_structure.pdb")
        if os.path.abspath(protein_path) != os.path.abspath(local_struct):
            import shutil
            shutil.copy2(protein_path, local_struct)
            protein_path = local_struct

    # ---- 2. Docking ----
    try:
        protein_text, is_seq = read_protein(protein_path or None, sequence or None)
        ligand_text, ligand_type = read_ligand(smiles)
        data = dock(
            protein_text, ligand_text, ligand_file_type=ligand_type,
            num_poses=args.num_poses, steps=args.steps,
            time_divisions=args.time_divisions, is_sequence=is_seq,
        )
    except (RuntimeError, SystemExit) as e:
        print(f"    对接失败: {str(e)[:80]}")
        result["status"] = "dock_error"
        return result

    ranking = save_poses(data, docking_dir)
    if ranking:
        result["best_conf"] = ranking[0][2]
        print(f"    对接完成 best_conf={ranking[0][2]:.3f} ({confidence_level(ranking[0][2])})")
    else:
        print("    对接完成 (无位姿)")
        result["status"] = "no_pose"
        return result

    if args.dock_only:
        return result

    # ---- 3. Interaction analysis + figures + report (reuse analyze_interactions) ----
    top_pose = None
    for f in sorted(os.listdir(docking_dir)):
        if f.startswith("rank1_"):
            top_pose = os.path.join(docking_dir, f)
            break
    if not top_pose:
        result["status"] = "no_top_pose"
        return result

    try:
        import analyze_interactions as ai
        # Merge + PLIP + PyMOL
        os.makedirs(interactions_dir, exist_ok=True)
        complex_pdb = os.path.join(interactions_dir, "complex.pdb")
        ai.merge_protein_ligand(protein_path, top_pose, complex_pdb)
        plip_ok = ai.run_plip(complex_pdb, interactions_dir)
        if plip_ok:
            if not args.no_pymol:
                ai.run_pymol_visualization(complex_pdb, interactions_dir)
            # clean up PLIP small pngs + intermediates
            keep = {"interaction_bindingsite.png", "complex_overview.png"}
            for f in os.listdir(interactions_dir):
                if f.endswith(".png") and f not in keep:
                    os.remove(os.path.join(interactions_dir, f))
                elif f.startswith("plipfixed.") or f.endswith("_protonated.pdb"):
                    os.remove(os.path.join(interactions_dir, f))
            # count interactions from PLIP report
            plip_txt = None
            for f in os.listdir(interactions_dir):
                if f.endswith("_report.txt"):
                    plip_txt = os.path.join(interactions_dir, f)
                    break
            if plip_txt:
                import generate_report as gr
                sections, _ = gr.parse_plip_report(plip_txt)
                result["n_interactions"] = len(gr.summarize_interactions(sections))
    except Exception as e:
        print(f"    相互作用分析失败: {str(e)[:100]}")
        result["status"] = "plip_error"
        return result

    # ---- 4. Report ----
    if not args.no_report:
        try:
            import generate_report as gr
            report_dir = os.path.join(comp_dir, "report")
            gr.generate(
                docking_dir=docking_dir,
                interactions_dir=interactions_dir,
                out_dir=report_dir,
                gene=label or None,
                accession=accession,
            )
        except Exception as e:
            print(f"    报告生成失败: {str(e)[:100]}")

    return result


def main():
    ap = argparse.ArgumentParser(description="End-to-end batch docking + analysis + report")
    ap.add_argument("--csv", required=True, help="Input CSV file")
    ap.add_argument("--out_dir", default="results/batch_full", help="Output directory")
    ap.add_argument("--num_poses", type=int, default=10)
    ap.add_argument("--steps", type=int, default=18)
    ap.add_argument("--time_divisions", type=int, default=20)
    ap.add_argument("--dock_only", action="store_true", help="Only dock, skip PLIP/figures/report")
    ap.add_argument("--no-pymol", action="store_true", help="Skip PyMOL figures")
    ap.add_argument("--no-report", action="store_true", help="Skip report generation")
    args = ap.parse_args()

    if not os.path.isfile(args.csv):
        sys.exit(f"ERROR: CSV not found: {args.csv}")
    os.makedirs(args.out_dir, exist_ok=True)

    with open(args.csv, newline="") as f:
        rows = list(csv.DictReader(f))
    if not rows:
        sys.exit("ERROR: CSV has no data rows.")

    print(f"[batch] 共 {len(rows)} 个化合物，完整流水线 "
          f"({'仅对接' if args.dock_only else '对接→PLIP→出图→报告'})\n")

    struct_cache = {}
    results = []
    for idx, row in enumerate(rows, 1):
        cname = (row.get("complex_name") or f"complex_{idx}").strip()
        print(f"[{idx}/{len(rows)}] {cname}")
        try:
            res = run_one(row, idx, args, struct_cache)
        except Exception as e:
            print(f"    未预期错误: {str(e)[:120]}")
            res = {"name": _safe_name(cname), "gene": "", "smiles": "",
                   "best_conf": None, "n_interactions": None, "status": "error"}
        results.append(res)
        print()

    # ---- master summary ----
    summary_path = os.path.join(args.out_dir, "batch_summary.csv")
    ranked = sorted(
        results,
        key=lambda r: (r["best_conf"] if r["best_conf"] is not None else -1e9),
        reverse=True,
    )
    with open(summary_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["rank", "compound", "target", "best_confidence", "reliability",
                    "n_interactions", "status", "smiles"])
        for r, res in enumerate(ranked, 1):
            bc = f"{res['best_conf']:.4f}" if res["best_conf"] is not None else ""
            lvl = confidence_level(res["best_conf"]) if res["best_conf"] is not None else "N/A"
            ni = res["n_interactions"] if res["n_interactions"] is not None else ""
            w.writerow([r, res["name"], res["gene"], bc, lvl, ni, res["status"], res["smiles"]])

    print("=" * 64)
    print(f"[batch] 完成! 汇总表: {summary_path}")
    print("[batch] 按最优置信度排名:")
    print(f"  {'#':<3} {'化合物':<20} {'置信度':>8} {'相互作用':>8}  状态")
    for r, res in enumerate(ranked[:15], 1):
        bc = f"{res['best_conf']:.3f}" if res["best_conf"] is not None else "N/A"
        ni = str(res["n_interactions"]) if res["n_interactions"] is not None else "-"
        print(f"  {r:<3} {res['name'][:20]:<20} {bc:>8} {ni:>8}  {res['status']}")
    print("=" * 64)


if __name__ == "__main__":
    main()
