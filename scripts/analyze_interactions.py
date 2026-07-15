#!/usr/bin/env python3
"""
Analyze protein-ligand interactions using PLIP and optionally visualize with PyMOL.

Pipeline:
  1. Merge protein PDB + ligand SDF into a combined complex PDB (required by PLIP)
  2. Run `plip -i com.pdb -ypxt` in the PLIP conda environment
  3. Optionally generate PyMOL visualization (interaction diagram PNG)

Usage:
    # After a docking run, analyze the top-ranked pose:
    python scripts/analyze_interactions.py \
        --protein protein.pdb \
        --ligand results/single/rank1_conf-0.150.sdf \
        --out_dir results/interactions/

    # With gene name (auto-fetch structure) + specific pose:
    python scripts/analyze_interactions.py \
        --gene EGFR \
        --ligand results/single/rank1_conf-0.150.sdf \
        --out_dir results/interactions/

    # Skip PyMOL visualization:
    python scripts/analyze_interactions.py \
        --protein protein.pdb \
        --ligand ligand.sdf \
        --out_dir results/interactions/ \
        --no-pymol
"""
import argparse
import os
import subprocess
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

ENV_NAME = "PLIP"
CONDA_SH = os.path.expanduser("~/anaconda3/etc/profile.d/conda.sh")
# Alternative conda paths
CONDA_SH_CANDIDATES = [
    CONDA_SH,
    os.path.expanduser("~/miniconda3/etc/profile.d/conda.sh"),
    "/opt/conda/etc/profile.d/conda.sh",
]


def find_conda_sh():
    for path in CONDA_SH_CANDIDATES:
        if os.path.isfile(path):
            return path
    return None


def run_in_plip_env(cmd, cwd=None):
    """Run a shell command inside the PLIP conda environment."""
    conda_sh = find_conda_sh()
    if conda_sh:
        full = f"source {conda_sh} && conda activate {ENV_NAME} && {cmd}"
    else:
        # fallback: use conda run
        conda_exe = os.environ.get("CONDA_EXE", "conda")
        full = f"{conda_exe} run -n {ENV_NAME} {cmd}"

    result = subprocess.run(
        full, shell=True, executable="/bin/bash",
        capture_output=True, text=True, cwd=cwd
    )
    return result


def merge_protein_ligand(protein_pdb, ligand_sdf, out_pdb):
    """
    Merge protein PDB + ligand SDF into a single complex PDB that PLIP can parse.
    Uses openbabel (in PLIP env) to convert ligand SDF → PDB with HETATM records,
    then concatenates with protein.
    """
    # Convert ligand SDF to PDB format (as HETATM records)
    ligand_pdb_tmp = out_pdb + ".lig.pdb"
    result = run_in_plip_env(f"obabel '{ligand_sdf}' -O '{ligand_pdb_tmp}' -opdb")
    if result.returncode != 0:
        print(f"[warn] obabel 转换失败: {result.stderr.strip()}")
        # Try python-based approach as fallback
        result2 = run_in_plip_env(
            f"python -c \""
            f"from openbabel import openbabel as ob;"
            f"conv=ob.OBConversion();"
            f"conv.SetInAndOutFormats('sdf','pdb');"
            f"mol=ob.OBMol();"
            f"conv.ReadFile(mol,'{ligand_sdf}');"
            f"conv.WriteFile(mol,'{ligand_pdb_tmp}')\""
        )
        if result2.returncode != 0:
            sys.exit(f"ERROR: 无法将配体 SDF 转为 PDB: {result2.stderr}")

    # Rename ligand residue to UNL (unknown ligand) and chain to Z for clarity
    lig_lines = []
    if os.path.isfile(ligand_pdb_tmp):
        with open(ligand_pdb_tmp) as f:
            for line in f:
                if line.startswith(("ATOM", "HETATM")):
                    # Force HETATM, residue name UNL, chain Z
                    edited = "HETATM" + line[6:17] + "UNL" + " Z" + line[22:]
                    lig_lines.append(edited)
                elif line.startswith("CONECT"):
                    lig_lines.append(line)
        os.remove(ligand_pdb_tmp)

    if not lig_lines:
        sys.exit("ERROR: 配体转换后无原子记录")

    # Read protein PDB — CLEAN: remove all existing HETATM (crystallization additives,
    # co-crystallized ligands, waters, ions) to keep only the pure protein + our docked ligand.
    prot_lines = []
    removed_hetatm = 0
    with open(protein_pdb) as f:
        for line in f:
            if line.startswith(("END", "MASTER")):
                continue
            # Remove existing HETATM records (waters, ions, co-crystallized ligands, etc.)
            if line.startswith("HETATM"):
                removed_hetatm += 1
                continue
            # Also skip CONECT records that reference removed atoms
            if line.startswith("CONECT"):
                continue
            prot_lines.append(line)
    if removed_hetatm:
        print(f"[clean] 已去除受体中 {removed_hetatm} 个原有杂原子记录（结晶添加剂/共结晶配体/水/离子）")

    # Write combined
    with open(out_pdb, "w") as f:
        for line in prot_lines:
            f.write(line)
        f.write("TER\n")
        for line in lig_lines:
            f.write(line)
        f.write("END\n")

    n_prot = sum(1 for l in prot_lines if l.startswith(("ATOM", "HETATM")))
    n_lig = len([l for l in lig_lines if l.startswith("HETATM")])
    print(f"[merge] 合并完成: {n_prot} 蛋白原子 + {n_lig} 配体原子 → {out_pdb}")
    return out_pdb


def run_plip(complex_pdb, out_dir):
    """Run PLIP analysis: plip -f com.pdb -ypxt"""
    os.makedirs(out_dir, exist_ok=True)
    cmd = f"plip -f '{os.path.abspath(complex_pdb)}' -ypxt -o '{os.path.abspath(out_dir)}'"
    print(f"[plip] 运行: {cmd}")
    result = run_in_plip_env(cmd, cwd=out_dir)
    if result.returncode != 0:
        print(f"[error] PLIP 运行失败:")
        print(result.stderr[:2000])
        return False
    if result.stdout:
        # Print key interaction summary lines
        for line in result.stdout.splitlines():
            if any(k in line.lower() for k in ["interaction", "bond", "contact", "pi"]):
                print(f"  {line}")
    print(f"[plip] 分析完成，结果在 {out_dir}/")
    return True


def run_pymol_visualization(complex_pdb, out_dir):
    """
    Generate a publication-quality PyMOL interaction visualization.
    Loads the PLIP .pse file + complex PDB, applies aesthetic settings,
    labels key interacting residues from the PLIP report, and renders at 300 DPI.
    """
    # Find the .pse file that PLIP generated (only UNL ones, our docked ligand)
    pse_files = [f for f in os.listdir(out_dir) if f.endswith(".pse") and "UNL" in f]
    if not pse_files:
        # fallback: any .pse
        pse_files = [f for f in os.listdir(out_dir) if f.endswith(".pse")]
    if not pse_files:
        print("[pymol] 未找到 PLIP 生成的 .pse 文件，跳过 PyMOL 绘图")
        return False

    pse_path = os.path.abspath(os.path.join(out_dir, pse_files[0]))
    complex_path = os.path.abspath(complex_pdb)
    png_path = os.path.abspath(os.path.join(out_dir, "interaction_bindingsite.png"))
    overview_path = os.path.abspath(os.path.join(out_dir, "complex_overview.png"))

    # Parse PLIP report to get interacting residues for labeling
    interacting_residues = _parse_interacting_residues(out_dir)
    label_commands = _build_label_commands(interacting_residues)

    # Use relative paths in pymol script — .pse files sometimes have issues with
    # absolute paths. Script runs with cwd=out_dir so relative names work.
    pse_basename = pse_files[0]
    complex_basename = os.path.basename(complex_pdb)
    # Copy complex to out_dir if not already there
    complex_in_outdir = os.path.join(out_dir, complex_basename)
    if not os.path.isfile(complex_in_outdir):
        import shutil
        shutil.copy2(complex_pdb, complex_in_outdir)

    pymol_script = f"""
import pymol
from pymol import cmd

# Load PLIP session (contains interaction visualizations) and complex
cmd.load('{pse_basename}')
cmd.load('{complex_basename}', 'complex_struct')

# ===== 1. Background and overall style =====
cmd.bg_color('white')
cmd.set('ray_trace_mode', 1)
cmd.set('ray_trace_fog', 0)
cmd.set('depth_cue', 0)

# ===== 2. Anti-aliasing and sharpness =====
cmd.set('antialias', 2)
cmd.set('line_smooth', 1)
cmd.set('hash_max', 300)

# ===== 3. Lighting =====
cmd.set('light_count', 8)
cmd.set('specular', 0.3)
cmd.set('shininess', 50)
cmd.set('ambient', 0.7)
cmd.set('direct', 0.4)
cmd.set('reflect', 0.1)

# ===== 4. Cartoon and surface rendering =====
cmd.set('cartoon_oval_length', 1.2)
cmd.set('cartoon_oval_width', 0.4)
cmd.set('cartoon_rect_length', 1.2)
cmd.set('cartoon_rect_width', 0.5)
cmd.set('cartoon_loop_radius', 0.7)
cmd.set('cartoon_transparency', 0.8)
cmd.set('surface_quality', 1)

# ===== 5. Label key interacting residues =====
cmd.set('label_font_id', 7)
cmd.set('label_size', 18)
cmd.set('label_color', 'black')
cmd.set('label_shadow_mode', 1)

{label_commands}

# ===== 6. IMAGE 1: binding-site close-up (cartoon transparency 0.8, with labels) =====
cmd.set('cartoon_transparency', 0.8)
cmd.zoom('resn UNL', buffer=8)
cmd.ray(2400, 1800)
cmd.png('interaction_bindingsite.png', dpi=300)

# ===== 7. IMAGE 2: overall view — full protein cartoon (transparency 0) + ligand =====
# Hide residue labels for the clean overview
cmd.label('all', '""')
# Solid cartoon for the whole protein
cmd.set('cartoon_transparency', 0.0)
# Make sure the whole protein is shown as cartoon and ligand as sticks
cmd.show('cartoon', 'complex_struct')
cmd.show('sticks', 'resn UNL')
cmd.color('cyan', 'resn UNL and elem C')
# Zoom out to show the entire complex
cmd.orient()
cmd.zoom('all', buffer=3)
cmd.ray(2400, 1800)
cmd.png('complex_overview.png', dpi=300)

cmd.quit()
"""
    script_path = os.path.join(out_dir, "_pymol_render.py")
    with open(script_path, "w") as f:
        f.write(pymol_script)

    # Run pymol with cwd set to out_dir. Use the script basename since cwd is out_dir.
    abs_out_dir = os.path.abspath(out_dir)
    cmd_str = f"cd '{abs_out_dir}' && pymol -cq '_pymol_render.py'"
    print(f"[pymol] 渲染高质量交互图 (300 DPI, 含残基标签): {png_path}")
    if interacting_residues:
        print(f"[pymol] 标注残基: {', '.join(interacting_residues)}")

    # Use system pymol directly (not through PLIP env) — pymol rendering
    # doesn't need PLIP's python packages, and conda activate can break pymol paths.
    result = subprocess.run(
        cmd_str, shell=True, executable="/bin/bash",
        capture_output=True, text=True,
    )

    if not os.path.isfile(png_path):
        # Fallback: try PLIP env pymol
        if result.stderr:
            print(f"[pymol] 系统 pymol stderr: {result.stderr[:300]}")
        print("[pymol] 系统 PyMOL 未生成图片，尝试 PLIP 环境 PyMOL...")
        result2 = run_in_plip_env(cmd_str)
        if not os.path.isfile(png_path) and result2.stderr:
            print(f"[pymol] PLIP env pymol stderr: {result2.stderr[:300]}")

    # Clean up script
    if os.path.isfile(script_path):
        os.remove(script_path)

    made = []
    if os.path.isfile(png_path):
        size_mb = os.path.getsize(png_path) / 1024 / 1024
        print(f"[pymol] ✓ 结合位点特写图: {png_path} ({size_mb:.1f} MB)")
        made.append(png_path)
    if os.path.isfile(overview_path):
        size_mb = os.path.getsize(overview_path) / 1024 / 1024
        print(f"[pymol] ✓ 整体图: {overview_path} ({size_mb:.1f} MB)")
        made.append(overview_path)

    if made:
        return True
    else:
        print(f"[pymol] ✗ 渲染失败，请检查 PyMOL 是否正确安装")
        print(f"[pymol]   可手动打开 .pse 文件查看: {pse_path}")
        return False


def _parse_interacting_residues(out_dir):
    """
    Parse PLIP text report to extract interacting residue names (e.g. 'LEU718', 'MET793').
    Returns a sorted list of unique 'RESTYPE+RESNR' strings.
    """
    report_file = None
    for f in os.listdir(out_dir):
        if f.endswith("_report.txt"):
            report_file = os.path.join(out_dir, f)
            break
    if not report_file:
        return []

    residues = set()
    in_unl_section = False
    with open(report_file) as f:
        for line in f:
            if "UNL:Z:1" in line:
                in_unl_section = True
                continue
            # Stop at the next ligand section (if any)
            if in_unl_section and line.strip() and not line.startswith((" ", "|", "+", "=", "*", "-")):
                if "SMALLMOLECULE" in line or "POLYMER" in line:
                    break
            if in_unl_section and line.startswith("|"):
                parts = [p.strip() for p in line.split("|")]
                # Table rows have: | RESNR | RESTYPE | RESCHAIN | ...
                if len(parts) >= 4:
                    try:
                        resnr = parts[1].strip()
                        restype = parts[2].strip()
                        if resnr.isdigit() and restype.isalpha() and len(restype) <= 4:
                            residues.add(f"{restype}{resnr}")
                    except (IndexError, ValueError):
                        pass

    return sorted(residues)


def _build_label_commands(residues):
    """
    Build PyMOL commands to label each interacting residue's CA atom
    with 'RESTYPE RESNR' (e.g. 'LEU 718').
    """
    if not residues:
        return "# No interacting residues found to label"

    lines = []
    for res in residues:
        # Split into restype and resnr (e.g. 'LEU718' -> 'LEU', '718')
        i = 0
        while i < len(res) and res[i].isalpha():
            i += 1
        restype = res[:i]
        resnr = res[i:]
        if not resnr:
            continue
        sel_name = f"res_{restype}{resnr}"
        lines.append(f"cmd.select('{sel_name}', 'resi {resnr} and name CA and chain A')")
        lines.append(f"cmd.label('{sel_name}', '\"{restype} {resnr}\"')")
    return "\n".join(lines)


def main():
    ap = argparse.ArgumentParser(description="PLIP interaction analysis + PyMOL visualization")
    ap.add_argument("--protein", help="Protein PDB file")
    ap.add_argument("--gene", help="Gene symbol (auto-fetch structure)")
    ap.add_argument("--protein_name", help="Protein name (auto-fetch structure)")
    ap.add_argument("--organism", default="9606", help="Organism taxon (default: human)")
    ap.add_argument("--ligand", required=True, help="Ligand SDF/MOL2 file (e.g. docking pose)")
    ap.add_argument("--out_dir", default="results/interactions", help="Output directory")
    ap.add_argument("--no-pymol", action="store_true", help="Skip PyMOL visualization")
    ap.add_argument("--no-report", action="store_true", help="Skip SCI report generation")
    ap.add_argument("--docking_dir", help="Docking output dir (for report; auto-detected if omitted)")
    ap.add_argument("--gene_label", help="Gene/protein name to use in the report")
    ap.add_argument("--accession", default="N/A", help="UniProt accession for the report")
    args = ap.parse_args()

    # Resolve protein
    if args.gene or args.protein_name:
        from fetch_structure import resolve_and_download
        os.makedirs(args.out_dir, exist_ok=True)
        label = args.gene or args.protein_name
        fetched = os.path.join(args.out_dir, f"{label}_structure.pdb")
        print(f"[fetch] 自动获取结构: {label}")
        args.protein, _ = resolve_and_download(
            gene=args.gene, protein_name=args.protein_name,
            organism=args.organism, out=fetched,
        )
        print()

    if not args.protein:
        sys.exit("ERROR: 需要 --protein 或 --gene/--protein_name")
    if not os.path.isfile(args.protein):
        sys.exit(f"ERROR: 蛋白文件不存在: {args.protein}")
    if not os.path.isfile(args.ligand):
        sys.exit(f"ERROR: 配体文件不存在: {args.ligand}")

    os.makedirs(args.out_dir, exist_ok=True)

    # 1. Check PLIP environment
    print("[check] 验证 PLIP 环境...")
    result = run_in_plip_env("plip -h")
    if result.returncode != 0:
        print("[error] PLIP 环境未就绪，请先运行: python scripts/setup_plip.py")
        sys.exit(1)
    print("[check] ✓ PLIP 可用\n")

    # 2. Merge protein + ligand → complex PDB
    complex_pdb = os.path.join(args.out_dir, "complex.pdb")
    merge_protein_ligand(args.protein, args.ligand, complex_pdb)

    # 3. Run PLIP
    print()
    success = run_plip(complex_pdb, args.out_dir)
    if not success:
        sys.exit(1)

    # 4. PyMOL visualization
    if not args.no_pymol:
        print()
        pymol_ok = run_pymol_visualization(complex_pdb, args.out_dir)
        if not pymol_ok:
            print("[pymol] ⚠ 高清图渲染失败，PLIP 分析结果仍可用(见 .pse 文件)")

    # 5. Clean up PLIP auto-generated small PNGs (we have our own high-res ones)
    #    Keep only our two output images
    keep_pngs = {"interaction_bindingsite.png", "complex_overview.png"}
    for f in os.listdir(args.out_dir):
        if f.endswith(".png") and f not in keep_pngs:
            os.remove(os.path.join(args.out_dir, f))

    # 6. Clean up intermediate files users don't need
    for f in os.listdir(args.out_dir):
        if f.startswith("plipfixed.") or f.endswith("_protonated.pdb"):
            os.remove(os.path.join(args.out_dir, f))

    # Summary
    print("\n" + "=" * 60)
    print("分析完成! 结果文件:")
    print("-" * 60)
    for f in sorted(os.listdir(args.out_dir)):
        fpath = os.path.join(args.out_dir, f)
        size = os.path.getsize(fpath)
        if size > 0 and not f.startswith("_"):
            ext = os.path.splitext(f)[1]
            desc = {
                ".pdb": "复合物结构(纯蛋白+对接配体)",
                ".pse": "PyMOL 会话(双击打开可交互查看)",
                ".xml": "相互作用详情(结构化 XML)",
                ".txt": "相互作用报告(文本,含残基/距离/类型)",
            }.get(ext, "")
            if f == "interaction_bindingsite.png":
                desc = "★ 结合位点特写图(300 DPI, cartoon透明0.8, 含残基标签)"
            elif f == "complex_overview.png":
                desc = "★ 整体图(300 DPI, cartoon透明0, 蛋白+配体全貌)"
            size_h = f"{size/1024/1024:.1f} MB" if size > 1024*1024 else f"{size/1024:.0f} KB"
            print(f"  {f:40s} {size_h:>8s}  {desc}")
    print("=" * 60)

    # 7. Generate SCI-grade bilingual report (report/ sibling of interactions/)
    if not args.no_report:
        print()
        # Auto-detect docking dir: sibling 'docking' folder next to out_dir
        docking_dir = args.docking_dir
        if not docking_dir:
            parent = os.path.dirname(os.path.abspath(args.out_dir))
            candidate = os.path.join(parent, "docking")
            if os.path.isdir(candidate):
                docking_dir = candidate
        report_dir = os.path.join(
            os.path.dirname(os.path.abspath(args.out_dir)), "report"
        )
        if docking_dir and os.path.isdir(docking_dir):
            try:
                from generate_report import generate  # noqa: E402
                generate(
                    docking_dir=docking_dir,
                    interactions_dir=os.path.abspath(args.out_dir),
                    out_dir=report_dir,
                    gene=args.gene_label or args.gene or args.protein_name,
                    accession=args.accession,
                )
            except Exception as e:
                print(f"[report] ⚠ 报告生成失败: {e}")
                print(f"[report]   可手动运行: python scripts/generate_report.py "
                      f"--docking_dir {docking_dir} --interactions_dir {args.out_dir}")
        else:
            print("[report] 未找到 docking 目录，跳过报告生成。")
            print("[report] 可手动指定: python scripts/generate_report.py "
                  f"--docking_dir <路径> --interactions_dir {args.out_dir}")


if __name__ == "__main__":
    main()
