#!/usr/bin/env python3
"""
Generate an SCI-publication-grade bilingual (中文/English) analysis report
from a completed ZYDock docking + PLIP interaction workflow.

Reads computational parameters and results from:
  - docking/raw_response.json      (ligand SMILES, poses, DiffDock params)
  - docking/confidence_scores.txt  (ranked confidence scores)
  - docking/*_structure.pdb        (receptor; header has PDB id, method, resolution)
  - interactions/complex_report.txt (PLIP interactions)
  - interactions/*.png             (figures)

Writes to a `report/` subfolder (sibling of docking/ and interactions/):
  - report_EN.md   (English Materials & Methods + Results)
  - report_CN.md   (Chinese 材料与方法 + 结果分析)
  - interactions_summary.csv (machine-readable interaction table)

Usage:
    python scripts/generate_report.py --run_dir results/demo_gefitinib
    # or point at docking/ and interactions/ explicitly:
    python scripts/generate_report.py --docking_dir A/docking --interactions_dir A/interactions --out_dir A/report
"""
import argparse
import json
import os
import re
import sys


# ---------- parsers ----------

def parse_structure_header(pdb_path):
    """Extract PDB id / title / resolution / method from a PDB header."""
    info = {"pdb_id": None, "title": "", "resolution": None, "method": None}
    if not pdb_path or not os.path.isfile(pdb_path):
        return info
    title_parts = []
    with open(pdb_path) as f:
        for line in f:
            rec = line[:6].strip()
            if rec == "HEADER":
                # last token of HEADER is often the PDB id
                tok = line[62:66].strip()
                if tok:
                    info["pdb_id"] = tok
            elif rec == "TITLE":
                title_parts.append(line[10:].strip())
            elif rec == "EXPDTA":
                exp = line[10:].strip()
                if "X-RAY" in exp.upper():
                    info["method"] = "X-ray diffraction"
                elif "ELECTRON MICROSCOPY" in exp.upper() or "EM" in exp.upper():
                    info["method"] = "cryo-EM"
                elif "NMR" in exp.upper():
                    info["method"] = "NMR"
                else:
                    info["method"] = exp.title()
            elif rec == "REMARK" and "RESOLUTION." in line and "ANGSTROMS" in line:
                m = re.search(r"([0-9]+\.[0-9]+)\s+ANGSTROMS", line)
                if m:
                    info["resolution"] = float(m.group(1))
            elif rec == "ATOM":
                break
    info["title"] = " ".join(title_parts).strip()
    return info


def parse_confidence(conf_path):
    """Return list of (rank, pose, confidence, level)."""
    rows = []
    if not os.path.isfile(conf_path):
        return rows
    with open(conf_path) as f:
        next(f, None)
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 4:
                try:
                    rows.append((parts[0], parts[1], float(parts[2]), parts[3]))
                except ValueError:
                    pass
    return rows


def parse_plip_report(report_path):
    """
    Parse PLIP text report into a dict keyed by interaction type.
    Each value is a list of dicts with residue + geometry fields.
    """
    if not os.path.isfile(report_path):
        return {}, None
    with open(report_path) as f:
        content = f.read()

    # PDB id from header line "Prediction ... for PDB structure XXXX"
    plip_pdb = None
    m = re.search(r"for PDB structure (\w+)", content)
    if m:
        plip_pdb = m.group(1)

    sections = {}
    # Split into interaction blocks by "**Name**"
    blocks = re.split(r"\*\*(.+?)\*\*", content)
    # blocks: [pre, name1, body1, name2, body2, ...]
    for i in range(1, len(blocks) - 1, 2):
        name = blocks[i].strip()
        body = blocks[i + 1]
        rows = _parse_plip_table(body)
        if rows:
            sections[name] = rows
    return sections, plip_pdb


def _parse_plip_table(body):
    """Parse a markdown-ish PLIP table body into list of dict rows."""
    lines = [l for l in body.splitlines() if l.strip().startswith("|")]
    if len(lines) < 2:
        return []
    header = [h.strip() for h in lines[0].strip().strip("|").split("|")]
    rows = []
    for line in lines[1:]:
        # skip separator lines like |=====|
        if set(line.strip()) <= set("|=+- "):
            continue
        cells = [c.strip() for c in line.strip().strip("|").split("|")]
        if len(cells) != len(header):
            continue
        row = dict(zip(header, cells))
        # only keep genuine data rows (RESNR numeric)
        if row.get("RESNR", "").isdigit():
            rows.append(row)
    return rows


def parse_docking_params(raw_json_path):
    """Extract ligand SMILES + pose count from raw API response."""
    info = {"ligand_smiles": None, "num_poses": None}
    if not os.path.isfile(raw_json_path):
        return info
    with open(raw_json_path) as f:
        d = json.load(f)
    info["ligand_smiles"] = d.get("ligand")
    info["num_poses"] = len(d.get("ligand_positions") or [])
    return info


# ---------- interaction summarization ----------

def summarize_interactions(sections):
    """
    Return a flat list of interaction records and counts by type.
    Record: {type, residue, restype, resnr, distance, extra}
    """
    records = []
    for itype, rows in sections.items():
        for r in rows:
            resnr = r.get("RESNR", "")
            restype = r.get("RESTYPE", "")
            residue = f"{restype}{resnr}"
            dist = r.get("DIST") or r.get("DIST_H-A") or r.get("DIST_D-A") or ""
            rec = {
                "type": itype,
                "residue": residue,
                "restype": restype,
                "resnr": resnr,
                "distance": dist,
            }
            records.append(rec)
    return records


TYPE_EN = {
    "Hydrophobic Interactions": "hydrophobic interaction",
    "Hydrogen Bonds": "hydrogen bond",
    "Salt Bridges": "salt bridge",
    "Halogen Bonds": "halogen bond",
    "pi-Stacking": "π-stacking",
    "pi-Cation Interactions": "π-cation interaction",
    "Water Bridges": "water bridge",
}
TYPE_CN = {
    "Hydrophobic Interactions": "疏水相互作用",
    "Hydrogen Bonds": "氢键",
    "Salt Bridges": "盐桥",
    "Halogen Bonds": "卤键",
    "pi-Stacking": "π-π 堆积",
    "pi-Cation Interactions": "π-阳离子相互作用",
    "Water Bridges": "水桥",
}


def _fmt_residue_list(records, itype):
    res = []
    for r in records:
        if r["type"] == itype:
            d = f" ({r['distance']} Å)" if r["distance"] else ""
            res.append(f"{r['residue']}{d}")
    return res


# ---------- report writers ----------

def write_english_report(ctx, out_path):
    s = ctx
    lines = []
    lines.append("# Molecular Docking and Protein–Ligand Interaction Analysis")
    lines.append("")
    lines.append("## Materials and Methods")
    lines.append("")
    lines.append("### Receptor structure preparation")
    res_txt = f"{s['resolution']:.2f} Å" if s.get("resolution") else "N/A"
    method = s.get("method") or "experimental"
    lines.append(
        f"The three-dimensional structure of the target protein "
        f"({s['gene_or_name']}, UniProt {s['accession']}) was retrieved from the "
        f"Protein Data Bank (PDB ID: {s['pdb_id']}), an {method} structure "
        f"determined at {res_txt} resolution. Where multiple experimental "
        f"structures were available, the entry with the highest resolution was "
        f"selected. All heteroatoms (crystallographic water molecules, ions, "
        f"buffer additives and co-crystallized ligands) were removed prior to "
        f"docking, retaining only the protein chain(s)."
    )
    lines.append("")
    lines.append("### Ligand preparation")
    lines.append(
        f"The ligand was provided as a SMILES string "
        f"(`{s['ligand_smiles']}`) and processed into 3D conformers internally "
        f"by the docking engine."
    )
    lines.append("")
    lines.append("### Molecular docking")
    lines.append(
        f"Blind molecular docking was performed with DiffDock, a diffusion "
        f"generative model for protein–ligand pose prediction, accessed through "
        f"the NVIDIA NIM cloud API (model *mit/diffdock*). A total of "
        f"{s['num_poses']} poses were sampled per complex; no binding pocket was "
        f"specified a priori. Predicted poses were ranked by the model's "
        f"confidence score (higher is more reliable). The top-ranked pose "
        f"(confidence = {s['top_conf']:.3f}) was retained for interaction analysis."
    )
    lines.append("")
    lines.append("### Interaction analysis and visualization")
    lines.append(
        "Noncovalent protein–ligand interactions of the top-ranked complex were "
        "profiled with the Protein–Ligand Interaction Profiler (PLIP v3.0.1). "
        "Binding-mode figures were rendered with open-source PyMOL at 300 DPI "
        "using ray-traced rendering: a whole-complex overview (opaque cartoon) "
        "and a binding-site close-up (transparent cartoon with labelled "
        "interacting residues)."
    )
    lines.append("")
    lines.append("## Results")
    lines.append("")
    lines.append("### Docking confidence")
    lines.append(
        f"DiffDock generated {s['num_poses']} candidate binding poses. The "
        f"confidence scores of the top poses ranged from {s['top_conf']:.3f} to "
        f"{s['worst_conf']:.3f} ({s['level_desc']}). The highest-confidence pose "
        f"was selected as the representative binding mode (Table 1)."
    )
    lines.append("")
    lines.append("**Table 1. Confidence scores of the top-ranked docking poses.**")
    lines.append("")
    lines.append("| Rank | Confidence | Reliability |")
    lines.append("|------|-----------|-------------|")
    for rank, pose, conf, level in s["conf_rows"][:5]:
        lines.append(f"| {rank} | {conf:.3f} | {level} |")
    lines.append("")
    lines.append("### Binding mode and key interactions")
    total = len(s["records"])
    type_present = [t for t in s["sections"].keys()]
    type_str = ", ".join(TYPE_EN.get(t, t) for t in type_present)
    lines.append(
        f"PLIP analysis of the top-ranked complex identified {total} noncovalent "
        f"interaction(s), comprising {type_str}. The interacting residues and "
        f"their geometries are summarized in Table 2 and illustrated in "
        f"Figure 1 (overall view) and Figure 2 (binding-site close-up)."
    )
    lines.append("")
    for itype in type_present:
        reslist = _fmt_residue_list(s["records"], itype)
        if reslist:
            lines.append(
                f"- **{TYPE_EN.get(itype, itype).capitalize()}s** were observed "
                f"with {', '.join(reslist)}."
            )
    lines.append("")
    lines.append("**Table 2. Noncovalent interactions detected by PLIP.**")
    lines.append("")
    lines.append("| Interaction type | Residue | Distance (Å) |")
    lines.append("|------------------|---------|--------------|")
    for r in s["records"]:
        lines.append(f"| {TYPE_EN.get(r['type'], r['type'])} | {r['residue']} | {r['distance']} |")
    lines.append("")
    lines.append("### Figures")
    lines.append("")
    lines.append("![Figure 1. Overall view of the protein–ligand complex.](../interactions/complex_overview.png)")
    lines.append("")
    lines.append("**Figure 1.** Overall view of the predicted protein–ligand complex. "
                 "The receptor is shown as an opaque cartoon and the ligand as cyan sticks.")
    lines.append("")
    lines.append("![Figure 2. Binding-site close-up.](../interactions/interaction_bindingsite.png)")
    lines.append("")
    lines.append("**Figure 2.** Close-up of the binding site. The receptor is shown as a "
                 "transparent cartoon; interacting residues are labelled, and noncovalent "
                 "interactions are drawn as dashed lines by PLIP.")
    lines.append("")
    lines.append("---")
    lines.append("*Report auto-generated by the `ZYDock` skill. "
                 "Confidence scores reflect model certainty about the predicted pose, "
                 "not binding affinity; experimental validation is recommended.*")
    with open(out_path, "w") as f:
        f.write("\n".join(lines) + "\n")


def write_chinese_report(ctx, out_path):
    s = ctx
    lines = []
    lines.append("# 分子对接与蛋白–配体相互作用分析报告")
    lines.append("")
    lines.append("## 材料与方法")
    lines.append("")
    lines.append("### 受体结构准备")
    res_txt = f"{s['resolution']:.2f} Å" if s.get("resolution") else "未知"
    method_cn = {"X-ray diffraction": "X 射线衍射", "cryo-EM": "冷冻电镜", "NMR": "核磁共振"}.get(
        s.get("method"), s.get("method") or "实验解析")
    lines.append(
        f"目标蛋白（{s['gene_or_name']}，UniProt 登录号 {s['accession']}）的三维结构"
        f"取自蛋白质数据库（PDB ID：{s['pdb_id']}），为 {method_cn}结构，"
        f"分辨率 {res_txt}。当存在多个实验结构时，优先选取分辨率最高者。"
        f"对接前去除所有杂原子（结晶水、离子、缓冲添加剂及共结晶配体），"
        f"仅保留蛋白链。"
    )
    lines.append("")
    lines.append("### 配体准备")
    lines.append(
        f"配体以 SMILES 字符串形式提供（`{s['ligand_smiles']}`），"
        f"由对接引擎内部自动生成三维构象。"
    )
    lines.append("")
    lines.append("### 分子对接")
    lines.append(
        f"采用 DiffDock（一种基于扩散生成模型的蛋白–配体构象预测方法）进行盲对接，"
        f"通过 NVIDIA NIM 云端 API（模型 *mit/diffdock*）调用。"
        f"每个复合物采样 {s['num_poses']} 个结合构象，未预先指定结合口袋。"
        f"预测构象按模型置信度分数排序（分数越高越可靠），"
        f"选取置信度最高的构象（置信度 = {s['top_conf']:.3f}）用于后续相互作用分析。"
    )
    lines.append("")
    lines.append("### 相互作用分析与可视化")
    lines.append(
        "使用蛋白–配体相互作用分析工具 PLIP（v3.0.1）对最优复合物的非共价相互作用"
        "进行解析。结合模式图由开源 PyMOL 以光线追踪方式在 300 DPI 下渲染，"
        "包括整体视图（不透明 cartoon）和结合位点特写（半透明 cartoon 并标注关键残基）。"
    )
    lines.append("")
    lines.append("## 结果")
    lines.append("")
    lines.append("### 对接置信度")
    lines.append(
        f"DiffDock 共生成 {s['num_poses']} 个候选结合构象，"
        f"排名靠前构象的置信度分数介于 {s['top_conf']:.3f} 至 {s['worst_conf']:.3f} 之间"
        f"（{s['level_desc_cn']}）。选取置信度最高的构象作为代表性结合模式（表 1）。"
    )
    lines.append("")
    lines.append("**表 1. 排名靠前对接构象的置信度分数。**")
    lines.append("")
    lines.append("| 排名 | 置信度 | 可靠性 |")
    lines.append("|------|--------|--------|")
    level_cn = {"High": "高", "Moderate": "中等", "Low": "低"}
    for rank, pose, conf, level in s["conf_rows"][:5]:
        lines.append(f"| {rank} | {conf:.3f} | {level_cn.get(level, level)} |")
    lines.append("")
    lines.append("### 结合模式与关键相互作用")
    total = len(s["records"])
    type_present = [t for t in s["sections"].keys()]
    type_str = "、".join(TYPE_CN.get(t, t) for t in type_present)
    lines.append(
        f"PLIP 分析在最优复合物中共识别出 {total} 个非共价相互作用，"
        f"包括{type_str}。相互作用残基及其几何参数汇总于表 2，"
        f"并在图 1（整体视图）和图 2（结合位点特写）中展示。"
    )
    lines.append("")
    for itype in type_present:
        reslist = _fmt_residue_list(s["records"], itype)
        if reslist:
            lines.append(
                f"- **{TYPE_CN.get(itype, itype)}**：与 {', '.join(reslist)} 形成。"
            )
    lines.append("")
    lines.append("**表 2. PLIP 检测到的非共价相互作用。**")
    lines.append("")
    lines.append("| 相互作用类型 | 残基 | 距离 (Å) |")
    lines.append("|--------------|------|----------|")
    for r in s["records"]:
        lines.append(f"| {TYPE_CN.get(r['type'], r['type'])} | {r['residue']} | {r['distance']} |")
    lines.append("")
    lines.append("### 图示")
    lines.append("")
    lines.append("![图 1. 蛋白–配体复合物整体视图。](../interactions/complex_overview.png)")
    lines.append("")
    lines.append("**图 1.** 预测的蛋白–配体复合物整体视图。受体以不透明 cartoon 表示，"
                 "配体以青色棒状表示。")
    lines.append("")
    lines.append("![图 2. 结合位点特写。](../interactions/interaction_bindingsite.png)")
    lines.append("")
    lines.append("**图 2.** 结合位点特写。受体以半透明 cartoon 表示，标注了相互作用残基，"
                 "PLIP 以虚线绘制各类非共价相互作用。")
    lines.append("")
    lines.append("---")
    lines.append("*本报告由 `ZYDock` 技能自动生成。置信度分数反映模型对预测构象的确定性，"
                 "不代表结合亲和力；建议进行实验验证。*")
    with open(out_path, "w") as f:
        f.write("\n".join(lines) + "\n")


def write_interactions_csv(records, out_path):
    import csv
    with open(out_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["interaction_type", "residue", "restype", "resnr", "distance_A"])
        for r in records:
            w.writerow([r["type"], r["residue"], r["restype"], r["resnr"], r["distance"]])


# ---------- main ----------

def find_file(directory, predicate):
    if not directory or not os.path.isdir(directory):
        return None
    for f in sorted(os.listdir(directory)):
        if predicate(f):
            return os.path.join(directory, f)
    return None


def generate(docking_dir, interactions_dir, out_dir, gene=None, accession="N/A"):
    """
    Core report generation: parse all data and write EN/CN reports + CSV.
    Returns a summary dict.
    """
    os.makedirs(out_dir, exist_ok=True)

    raw_json = os.path.join(docking_dir, "raw_response.json")
    conf_txt = os.path.join(docking_dir, "confidence_scores.txt")
    struct_pdb = find_file(docking_dir, lambda f: f.endswith("_structure.pdb"))
    plip_txt = find_file(interactions_dir, lambda f: f.endswith("_report.txt"))

    dock = parse_docking_params(raw_json)
    conf_rows = parse_confidence(conf_txt)
    struct = parse_structure_header(struct_pdb)
    sections, plip_pdb = parse_plip_report(plip_txt) if plip_txt else ({}, None)
    records = summarize_interactions(sections)

    if not conf_rows:
        print("[report] ⚠ 未找到置信度数据")
    if not sections:
        print("[report] ⚠ 未找到 PLIP 相互作用数据")

    top_conf = conf_rows[0][2] if conf_rows else float("nan")
    worst_conf = conf_rows[-1][2] if conf_rows else float("nan")
    levels = {r[3] for r in conf_rows}
    level_desc = " and ".join(sorted(levels)) + " confidence" if levels else "N/A"
    level_desc_cn = "、".join(
        {"High": "高", "Moderate": "中等", "Low": "低"}.get(l, l) for l in sorted(levels)
    ) + "置信度" if levels else "N/A"

    gene_or_name = gene or (struct.get("title") or "the target protein").split(" IN ")[0].title()

    ctx = {
        "gene_or_name": gene_or_name,
        "accession": accession,
        "pdb_id": struct.get("pdb_id") or plip_pdb or "N/A",
        "resolution": struct.get("resolution"),
        "method": struct.get("method"),
        "ligand_smiles": dock.get("ligand_smiles") or "N/A",
        "num_poses": dock.get("num_poses") or len(conf_rows),
        "top_conf": top_conf,
        "worst_conf": worst_conf,
        "level_desc": level_desc,
        "level_desc_cn": level_desc_cn,
        "conf_rows": conf_rows,
        "sections": sections,
        "records": records,
    }

    en_path = os.path.join(out_dir, "report_EN.md")
    cn_path = os.path.join(out_dir, "report_CN.md")
    csv_path = os.path.join(out_dir, "interactions_summary.csv")
    write_english_report(ctx, en_path)
    write_chinese_report(ctx, cn_path)
    write_interactions_csv(records, csv_path)

    print("[report] ✓ SCI 级分析报告已生成:")
    print(f"[report]   英文: {en_path}")
    print(f"[report]   中文: {cn_path}")
    print(f"[report]   相互作用表: {csv_path}")
    print(f"[report]   摘要: PDB={ctx['pdb_id']}  分辨率={ctx['resolution']}  "
          f"最优置信度={top_conf:.3f}  相互作用数={len(records)}")
    return ctx


def main():
    ap = argparse.ArgumentParser(description="Generate bilingual SCI-grade docking report")
    ap.add_argument("--run_dir", help="Run dir containing docking/ and interactions/")
    ap.add_argument("--docking_dir", help="Docking output dir (override)")
    ap.add_argument("--interactions_dir", help="Interactions output dir (override)")
    ap.add_argument("--out_dir", help="Report output dir (default: <run_dir>/report)")
    ap.add_argument("--gene", help="Gene/protein name for the report (optional)")
    ap.add_argument("--accession", default="N/A", help="UniProt accession (optional)")
    args = ap.parse_args()

    if args.run_dir:
        docking_dir = args.docking_dir or os.path.join(args.run_dir, "docking")
        interactions_dir = args.interactions_dir or os.path.join(args.run_dir, "interactions")
        out_dir = args.out_dir or os.path.join(args.run_dir, "report")
    else:
        docking_dir = args.docking_dir
        interactions_dir = args.interactions_dir
        out_dir = args.out_dir or "report"

    if not docking_dir or not interactions_dir:
        sys.exit("ERROR: 需要 --run_dir，或同时指定 --docking_dir 和 --interactions_dir")

    generate(docking_dir, interactions_dir, out_dir, gene=args.gene, accession=args.accession)


if __name__ == "__main__":
    main()
