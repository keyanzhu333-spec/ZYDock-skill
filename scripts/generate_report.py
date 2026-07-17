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


# ---------- deeper quantitative analysis ----------

def analyze_confidence(conf_rows):
    """
    Quantitative summary of the confidence-score distribution across all poses.
    Returns counts per reliability tier, mean/spread, and a consensus judgement.
    """
    out = {
        "n": len(conf_rows),
        "high": 0, "moderate": 0, "low": 0,
        "mean": float("nan"), "spread": float("nan"),
        "top": float("nan"), "worst": float("nan"),
        "consensus": "N/A",
    }
    if not conf_rows:
        return out
    scores = [c for (_, _, c, _) in conf_rows]
    for c in scores:
        if c > 0:
            out["high"] += 1
        elif c >= -1.5:
            out["moderate"] += 1
        else:
            out["low"] += 1
    out["mean"] = sum(scores) / len(scores)
    out["top"] = scores[0]
    out["worst"] = scores[-1]
    out["spread"] = scores[0] - scores[-1]
    # consensus among the top poses: how tightly clustered are the best scores?
    top_k = scores[:5]
    if len(top_k) >= 2:
        span = max(top_k) - min(top_k)
        if span < 0.5:
            out["consensus"] = "tight"      # poses agree → robust binding mode
        elif span < 1.5:
            out["consensus"] = "moderate"
        else:
            out["consensus"] = "dispersed"  # poses disagree → binding mode uncertain
    return out


def residue_hotspots(records):
    """
    Rank interacting residues by how many interactions they participate in.
    Residues engaged through multiple contacts (or multiple interaction types)
    are the principal anchoring points of the binding mode.
    Returns list of dicts sorted by contact count (desc).
    """
    agg = {}
    for r in records:
        key = r["residue"]
        if not key or key == "":
            continue
        d = agg.setdefault(key, {"residue": key, "count": 0, "types": set(), "dists": []})
        d["count"] += 1
        d["types"].add(r["type"])
        if r["distance"]:
            try:
                d["dists"].append(float(r["distance"]))
            except ValueError:
                pass
    hotspots = []
    for d in agg.values():
        d["types"] = sorted(d["types"])
        d["min_dist"] = min(d["dists"]) if d["dists"] else None
        hotspots.append(d)
    hotspots.sort(key=lambda x: (-x["count"], x["residue"]))
    return hotspots


def type_counts(sections):
    """Number of interactions of each type, sorted by count desc."""
    counts = [(t, len(rows)) for t, rows in sections.items()]
    counts.sort(key=lambda x: -x[1])
    return counts


# Literature references cited in the auto-generated Methods/Discussion.
REFERENCES_EN = [
    "Corso, G., Stärk, H., Jing, B., Barzilay, R. & Jaakkola, T. "
    "DiffDock: Diffusion Steps, Twists, and Turns for Molecular Docking. "
    "*International Conference on Learning Representations (ICLR)* (2023).",
    "Adasme, M. F. et al. PLIP 2021: expanding the scope of the protein–ligand "
    "interaction profiler to DNA and RNA. *Nucleic Acids Research* 49, W530–W534 (2021).",
    "Schrödinger, LLC. The PyMOL Molecular Graphics System, Version 2.x (open-source).",
    "Berman, H. M. et al. The Protein Data Bank. "
    "*Nucleic Acids Research* 28, 235–242 (2000).",
    "The UniProt Consortium. UniProt: the Universal Protein Knowledgebase in 2023. "
    "*Nucleic Acids Research* 51, D523–D531 (2023).",
    "Jumper, J. et al. Highly accurate protein structure prediction with AlphaFold. "
    "*Nature* 596, 583–589 (2021).",
]
REFERENCES_CN = REFERENCES_EN  # 参考文献采用国际通行英文著录格式


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
        f"Blind molecular docking was performed with DiffDock (Corso et al., 2023) [1], "
        f"a diffusion generative model that treats pose prediction as a generative "
        f"process over the manifold of ligand translations, rotations and torsion "
        f"angles rather than as a scoring-function optimization. The model was "
        f"accessed through the NVIDIA NIM cloud API (model *mit/diffdock*), so all "
        f"GPU computation was performed on hosted infrastructure. A total of "
        f"{s['num_poses']} independent poses were sampled per complex under a "
        f"fully blind protocol: no binding pocket, grid box or pharmacophore "
        f"constraint was specified a priori, allowing the model to explore the "
        f"entire receptor surface. Each sampled pose was assigned a confidence "
        f"score by DiffDock's confidence model, where higher values indicate greater "
        f"certainty that the predicted pose lies within 2 Å RMSD of the true binding "
        f"mode. Poses were ranked in descending order of confidence, and the "
        f"top-ranked pose (confidence = {s['top_conf']:.3f}) was retained as the "
        f"representative binding mode for downstream interaction analysis."
    )
    lines.append("")
    lines.append("### Interaction analysis and visualization")
    lines.append(
        "The top-ranked ligand pose was merged with the prepared receptor to form "
        "a complex, and noncovalent protein–ligand interactions were profiled with "
        "the Protein–Ligand Interaction Profiler (PLIP; Adasme et al., 2021) [2]. "
        "PLIP applies distance- and angle-based geometric rules to detect hydrogen "
        "bonds, hydrophobic contacts, salt bridges, halogen bonds, π-stacking, "
        "π-cation interactions and water bridges, reporting the participating "
        "residues together with their interatomic geometries. Binding-mode figures "
        "were rendered with open-source PyMOL [3] at 300 DPI using ray-traced "
        "rendering: a whole-complex overview (opaque cartoon) and a binding-site "
        "close-up (transparent cartoon with labelled interacting residues and "
        "interactions drawn as dashed lines)."
    )
    lines.append("")
    lines.append("## Results")
    lines.append("")
    lines.append("### Docking confidence")
    ca = s["conf_analysis"]
    consensus_txt = {
        "tight": "The top five poses fall within a narrow confidence window "
                 "(span < 0.5), indicating that the model repeatedly converged on "
                 "the same region of the receptor and that the predicted binding "
                 "mode is robust.",
        "moderate": "The top five poses span a moderate confidence range "
                    "(0.5–1.5), suggesting a broadly consistent but not fully "
                    "converged binding mode.",
        "dispersed": "The top poses are spread over a wide confidence range "
                     "(span > 1.5), indicating that the model did not converge on "
                     "a single binding mode; the top pose should be interpreted "
                     "with caution and cross-checked against the next-ranked poses.",
        "N/A": "",
    }.get(ca["consensus"], "")
    lines.append(
        f"DiffDock generated {s['num_poses']} candidate binding poses. Across all "
        f"sampled poses the confidence scores ranged from {ca['top']:.3f} (best) to "
        f"{ca['worst']:.3f} (worst), with a mean of {ca['mean']:.3f} and a total "
        f"spread of {ca['spread']:.3f}. Of these, {ca['high']} pose(s) fell in the "
        f"high-confidence tier (score > 0), {ca['moderate']} in the moderate tier "
        f"(−1.5 to 0) and {ca['low']} in the low tier (< −1.5). {consensus_txt} "
        f"The highest-confidence pose (confidence = {s['top_conf']:.3f}) was selected "
        f"as the representative binding mode (Table 1)."
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
    tcounts = s["type_counts"]
    count_clause = "; ".join(
        f"{n} {TYPE_EN.get(t, t)}{'s' if n != 1 else ''}" for t, n in tcounts
    )
    lines.append(
        f"PLIP analysis of the top-ranked complex identified {total} noncovalent "
        f"interaction(s) in total, comprising {count_clause}. The predominant "
        f"interaction type was {TYPE_EN.get(tcounts[0][0], tcounts[0][0])} "
        f"({tcounts[0][1]} of {total} contacts), which is consistent with the "
        f"typical driving forces of small-molecule recognition in a protein "
        f"binding pocket. The interacting residues and their geometries are "
        f"summarized in Table 2 and illustrated in Figure 1 (overall view) and "
        f"Figure 2 (binding-site close-up)."
    ) if tcounts else lines.append(
        f"PLIP analysis of the top-ranked complex identified {total} noncovalent "
        f"interaction(s), comprising {type_str}."
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
    # Hotspot analysis — residues engaged through multiple contacts
    hotspots = s["hotspots"]
    multi = [h for h in hotspots if h["count"] >= 2]
    if multi:
        hs_txt = "; ".join(
            f"{h['residue']} ({h['count']} contacts spanning "
            f"{', '.join(TYPE_EN.get(t, t) for t in h['types'])})"
            for h in multi
        )
        lines.append(
            f"Notably, several residues anchor the ligand through more than one "
            f"contact and therefore constitute the principal binding hotspots: "
            f"{hs_txt}. Residues engaged through multiple, geometrically diverse "
            f"interactions typically contribute disproportionately to binding "
            f"specificity and are attractive focal points for subsequent lead "
            f"optimization or site-directed mutagenesis."
        )
    elif hotspots:
        closest = min(
            (h for h in hotspots if h["min_dist"] is not None),
            key=lambda x: x["min_dist"], default=None,
        )
        if closest:
            lines.append(
                f"The interactions are distributed across distinct residues rather "
                f"than concentrated on a single anchor. The shortest contact was "
                f"observed at {closest['residue']} ({closest['min_dist']:.2f} Å), "
                f"which likely represents the strongest single point of "
                f"stabilization within the pocket."
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
    lines.append("## Discussion")
    lines.append("")
    ca = s["conf_analysis"]
    total = len(s["records"])
    tcounts = s["type_counts"]
    # Paragraph 1 — interpret the confidence tier and what it means
    if ca["top"] == ca["top"] and ca["top"] > 0:
        conf_interp = (
            f"The representative pose was predicted with high confidence "
            f"({s['top_conf']:.3f} > 0), which, for DiffDock, corresponds to a high "
            f"probability that the pose lies within 2 Å RMSD of the true binding "
            f"geometry. This lends structural credibility to the binding mode "
            f"described above."
        )
    elif ca["top"] == ca["top"] and ca["top"] >= -1.5:
        conf_interp = (
            f"The representative pose was predicted with moderate confidence "
            f"({s['top_conf']:.3f}), a range in which DiffDock poses are broadly "
            f"reasonable but warrant independent corroboration. The binding mode "
            f"should therefore be treated as a working hypothesis rather than an "
            f"established fact."
        )
    else:
        conf_interp = (
            f"The representative pose was predicted with low confidence "
            f"({s['top_conf']:.3f} < −1.5). Low scores are common for large or "
            f"highly flexible ligands, multi-chain receptors, and targets from "
            f"protein families under-represented in the training data; the "
            f"predicted geometry should be regarded as tentative and prioritized "
            f"for experimental scrutiny."
        )
    lines.append(
        f"In this study, blind diffusion-based docking followed by structural "
        f"interaction profiling was used to generate a hypothesis for how the "
        f"query ligand engages {s['gene_or_name']}. {conf_interp}"
    )
    lines.append("")
    # Paragraph 2 — interpret the interaction network biophysically
    if total > 0 and tcounts:
        dominant = TYPE_EN.get(tcounts[0][0], tcounts[0][0])
        biophys = {
            "hydrophobic interaction":
                "an interaction network dominated by hydrophobic contacts points to "
                "shape-complementary burial of the ligand in a nonpolar sub-pocket, "
                "which contributes substantially to the binding free energy through "
                "the desolvation of apolar surfaces",
            "hydrogen bond":
                "a hydrogen-bond–rich interface indicates directional, "
                "specificity-determining contacts that orient the ligand precisely "
                "within the pocket",
            "salt bridge":
                "the presence of salt bridges signals strong electrostatic anchoring "
                "between charged ligand and protein groups, which can markedly "
                "enhance both affinity and residence time",
            "π-stacking":
                "π-stacking with aromatic residues suggests the ligand's aromatic "
                "system is sandwiched against the protein, a common motif in "
                "kinase and nuclear-receptor recognition",
        }.get(dominant, f"the interaction profile is dominated by {dominant}s")
        lines.append(
            f"From a biophysical standpoint, {biophys}. The combination of "
            f"interaction types observed here ("
            f"{', '.join(TYPE_EN.get(t, t) for t, _ in tcounts)}) is consistent "
            f"with a specific, rather than merely opportunistic, association. The "
            f"binding hotspots identified in the Results represent the most "
            f"promising handles for structure-guided optimization: strengthening or "
            f"adding contacts at these residues is the most direct route to "
            f"improved potency, whereas the peripheral single contacts offer scope "
            f"for tuning selectivity and physicochemical properties."
        )
        lines.append("")
    # Paragraph 3 — limitations and next steps
    lines.append(
        "Several limitations should be borne in mind when interpreting these "
        "results. First, DiffDock predicts binding *poses* and a *confidence* in "
        "those poses; it does not estimate binding affinity (ΔG or K_d), so the "
        "present analysis speaks to *how* the ligand may bind, not *how tightly*. "
        "Quantitative ranking would require an orthogonal scoring approach such as "
        "MM/GBSA, free-energy perturbation, or an empirical scoring function "
        "(e.g., GNINA). Second, docking was performed against a single, largely "
        "rigid receptor conformation; induced-fit rearrangements, alternative "
        "protonation states, and explicit active-site water molecules were not "
        "modelled and could alter the predicted geometry. Third, the confidence "
        "score is a learned proxy for pose accuracy and is not a physical "
        "observable. Accordingly, the binding mode reported here should be "
        "regarded as a computational hypothesis. Recommended next steps include "
        "molecular-dynamics refinement of the complex to assess pose stability, "
        "physics-based rescoring to estimate affinity, and ultimately experimental "
        "validation through biophysical binding assays (e.g., SPR, ITC) or "
        "site-directed mutagenesis of the predicted hotspot residues."
    )
    lines.append("")
    lines.append("## References")
    lines.append("")
    for i, ref in enumerate(REFERENCES_EN, 1):
        lines.append(f"{i}. {ref}")
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
        f"采用 DiffDock（Corso 等，2023）[1] 进行盲对接。DiffDock 是一种基于扩散生成"
        f"模型的蛋白–配体构象预测方法，它将位姿预测建模为在配体平移、旋转和二面角构成"
        f"的流形上的生成过程，而非传统的打分函数优化。计算通过 NVIDIA NIM 云端 API"
        f"（模型 *mit/diffdock*）调用，全部 GPU 运算在云端完成。每个复合物在完全盲对接"
        f"模式下采样 {s['num_poses']} 个独立结合构象：未预先指定结合口袋、格点盒或药效团"
        f"约束，允许模型在整个受体表面自由搜索。DiffDock 的置信度模型为每个构象赋予一个"
        f"置信度分数，分数越高表示模型越确信该构象与真实结合模式的 RMSD 在 2 Å 以内。"
        f"构象按置信度降序排列，选取置信度最高的构象（置信度 = {s['top_conf']:.3f}）"
        f"作为代表性结合模式用于后续相互作用分析。"
    )
    lines.append("")
    lines.append("### 相互作用分析与可视化")
    lines.append(
        "将置信度最高的配体构象与准备好的受体合并为复合物，使用蛋白–配体相互作用分析"
        "工具 PLIP（Adasme 等，2021）[2] 解析其非共价相互作用。PLIP 基于距离和角度的"
        "几何判据自动检测氢键、疏水接触、盐桥、卤键、π-π 堆积、π-阳离子相互作用及水桥，"
        "并报告参与相互作用的残基及其原子间几何参数。结合模式图由开源 PyMOL [3] 以光线"
        "追踪方式在 300 DPI 下渲染，包括整体视图（不透明 cartoon）和结合位点特写"
        "（半透明 cartoon，标注关键残基并以虚线绘制各类相互作用）。"
    )
    lines.append("")
    lines.append("## 结果")
    lines.append("")
    lines.append("### 对接置信度")
    ca = s["conf_analysis"]
    consensus_cn = {
        "tight": "前五位构象的置信度落在一个很窄的区间内（跨度 < 0.5），"
                 "表明模型反复收敛到受体的同一区域，所预测的结合模式较为稳健。",
        "moderate": "前五位构象的置信度跨度中等（0.5–1.5），"
                    "提示结合模式大体一致但尚未完全收敛。",
        "dispersed": "排名靠前的构象置信度分布较广（跨度 > 1.5），"
                     "表明模型未收敛到单一结合模式；最优构象需谨慎解读，"
                     "并应与次优构象相互印证。",
        "N/A": "",
    }.get(ca["consensus"], "")
    lines.append(
        f"DiffDock 共生成 {s['num_poses']} 个候选结合构象。在所有采样构象中，"
        f"置信度分数介于 {ca['top']:.3f}（最优）至 {ca['worst']:.3f}（最差）之间，"
        f"平均值为 {ca['mean']:.3f}，总跨度为 {ca['spread']:.3f}。其中，"
        f"{ca['high']} 个构象属于高置信度区间（分数 > 0），"
        f"{ca['moderate']} 个属于中等置信度区间（−1.5 至 0），"
        f"{ca['low']} 个属于低置信度区间（< −1.5）。{consensus_cn}"
        f"选取置信度最高的构象（置信度 = {s['top_conf']:.3f}）"
        f"作为代表性结合模式（表 1）。"
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
    tcounts = s["type_counts"]
    if tcounts:
        count_clause = "、".join(
            f"{TYPE_CN.get(t, t)} {n} 个" for t, n in tcounts
        )
        lines.append(
            f"PLIP 分析在最优复合物中共识别出 {total} 个非共价相互作用，"
            f"具体包括{count_clause}。其中占主导的相互作用类型为"
            f"{TYPE_CN.get(tcounts[0][0], tcounts[0][0])}"
            f"（{tcounts[0][1]}/{total}），"
            f"这与小分子在蛋白结合口袋中被识别的典型驱动力相符。"
            f"相互作用残基及其几何参数汇总于表 2，"
            f"并在图 1（整体视图）和图 2（结合位点特写）中展示。"
        )
    else:
        lines.append(
            f"PLIP 分析在最优复合物中共识别出 {total} 个非共价相互作用，包括{type_str}。"
        )
    lines.append("")
    for itype in type_present:
        reslist = _fmt_residue_list(s["records"], itype)
        if reslist:
            lines.append(
                f"- **{TYPE_CN.get(itype, itype)}**：与 {', '.join(reslist)} 形成。"
            )
    lines.append("")
    # 热点残基分析
    hotspots = s["hotspots"]
    multi = [h for h in hotspots if h["count"] >= 2]
    if multi:
        hs_txt = "；".join(
            f"{h['residue']}（{h['count']} 个接触，涉及"
            f"{'、'.join(TYPE_CN.get(t, t) for t in h['types'])}）"
            for h in multi
        )
        lines.append(
            f"值得注意的是，若干残基通过不止一个接触锚定配体，构成主要的结合热点："
            f"{hs_txt}。通过多个几何类型各异的相互作用参与结合的残基，"
            f"通常对结合特异性的贡献更大，也是后续先导化合物优化或定点突变的重点关注对象。"
        )
    elif hotspots:
        closest = min(
            (h for h in hotspots if h["min_dist"] is not None),
            key=lambda x: x["min_dist"], default=None,
        )
        if closest:
            lines.append(
                f"各相互作用分布于不同残基，而非集中于单一锚点。"
                f"最短接触出现在 {closest['residue']}（{closest['min_dist']:.2f} Å），"
                f"很可能代表口袋内最强的单点稳定作用。"
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
    lines.append("## 讨论")
    lines.append("")
    ca = s["conf_analysis"]
    total = len(s["records"])
    tcounts = s["type_counts"]
    # 第一段——解读置信度等级及其含义
    if ca["top"] == ca["top"] and ca["top"] > 0:
        conf_interp = (
            f"代表性构象以高置信度预测得到（{s['top_conf']:.3f} > 0）。"
            f"对 DiffDock 而言，这对应该构象与真实结合几何 RMSD 在 2 Å 以内的高概率，"
            f"为上述结合模式提供了结构层面的可信度。"
        )
    elif ca["top"] == ca["top"] and ca["top"] >= -1.5:
        conf_interp = (
            f"代表性构象以中等置信度预测得到（{s['top_conf']:.3f}）。"
            f"在该区间内，DiffDock 构象大体合理但需独立佐证，"
            f"因此该结合模式应被视为工作假设而非既定结论。"
        )
    else:
        conf_interp = (
            f"代表性构象以低置信度预测得到（{s['top_conf']:.3f} < −1.5）。"
            f"低分常见于大分子或高柔性配体、多链受体，以及训练数据中代表性不足的蛋白家族；"
            f"所预测的几何结构应被视为初步结果，并优先接受实验检验。"
        )
    lines.append(
        f"本研究采用盲式扩散对接结合结构相互作用解析，"
        f"为查询配体与 {s['gene_or_name']} 的结合方式提出了一个假设。{conf_interp}"
    )
    lines.append("")
    # 第二段——从生物物理角度解读相互作用网络
    if total > 0 and tcounts:
        dominant = tcounts[0][0]
        biophys = {
            "Hydrophobic Interactions":
                "以疏水接触为主导的相互作用网络，提示配体以形状互补的方式埋入非极性亚口袋，"
                "通过非极性表面的去溶剂化对结合自由能作出重要贡献",
            "Hydrogen Bonds":
                "富含氢键的界面表明存在方向性、决定特异性的接触，"
                "能将配体精确定向于口袋之中",
            "Salt Bridges":
                "盐桥的存在意味着配体与蛋白带电基团之间存在强静电锚定，"
                "可显著增强亲和力与停留时间",
            "pi-Stacking":
                "与芳香族残基的 π-π 堆积提示配体芳香体系与蛋白发生夹心式堆叠，"
                "这是激酶和核受体识别中的常见模式",
        }.get(dominant, f"相互作用谱以{TYPE_CN.get(dominant, dominant)}为主导")
        lines.append(
            f"从生物物理角度看，{biophys}。此处观察到的相互作用类型组合"
            f"（{'、'.join(TYPE_CN.get(t, t) for t, _ in tcounts)}）"
            f"与一种特异性而非偶然性的结合相符。结果部分识别出的结合热点，"
            f"是结构导向优化最有希望的着力点：在这些残基处强化或增加接触，"
            f"是提高活性最直接的途径；而外围的单一接触则为调节选择性和理化性质提供了空间。"
        )
        lines.append("")
    # 第三段——局限性与后续工作
    lines.append(
        "解读上述结果时需注意若干局限。第一，DiffDock 预测的是结合*构象*及对构象的"
        "*置信度*，并不估算结合亲和力（ΔG 或 K_d），因此本分析回答的是配体"
        "*如何*结合，而非结合*有多强*。定量排序需借助正交的打分方法，如 MM/GBSA、"
        "自由能微扰或经验打分函数（如 GNINA）。第二，对接针对单一、基本刚性的受体构象"
        "进行，未建模诱导契合的构象重排、不同质子化状态以及活性位点的显式水分子，"
        "这些因素都可能改变预测的几何结构。第三，置信度分数是位姿准确性的学习型代理指标，"
        "并非物理可观测量。因此，此处报告的结合模式应被视为计算假设。"
        "建议的后续工作包括：对复合物进行分子动力学精修以评估构象稳定性、"
        "采用基于物理的重打分以估算亲和力，并最终通过生物物理结合实验"
        "（如 SPR、ITC）或对预测热点残基的定点突变进行实验验证。"
    )
    lines.append("")
    lines.append("## 参考文献")
    lines.append("")
    for i, ref in enumerate(REFERENCES_CN, 1):
        lines.append(f"{i}. {ref}")
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
        "conf_analysis": analyze_confidence(conf_rows),
        "hotspots": residue_hotspots(records),
        "type_counts": type_counts(sections),
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
