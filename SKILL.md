---
name: ZYDock
description: Cloud-based molecular docking via NVIDIA NIM DiffDock API. Predict protein-ligand binding poses from PDB/SMILES with confidence scores and virtual screening — runs entirely on NVIDIA's hosted GPUs, NO local GPU or DiffDock install needed. Requires an NVIDIA API key. Use when the user wants to dock ligands, predict binding poses, run virtual screening, or find binding sites without local compute. Not for binding affinity prediction.
license: MIT
metadata:
    skill-author: zy
---

# DiffDock_ZY: 云端分子对接（NVIDIA API 版）

## 概述 Overview

本技能通过 **NVIDIA NIM 托管的 DiffDock API** 做蛋白-配体分子对接。所有计算都在 NVIDIA 云端 GPU 上完成，**本地无需 GPU、无需安装 DiffDock**，只要有网络和一个 NVIDIA API key 即可。

**核心能力：**
- 预测小分子配体与蛋白靶点的 3D 结合位姿（盲对接，无需指定口袋）
- 为每个位姿输出 confidence 置信度分数
- 支持单个对接和批量虚拟筛选
- 蛋白输入支持 PDB 文件或氨基酸序列；配体输入支持 SMILES 或 SDF

**关键区分：** DiffDock 预测的是**结合位姿**（3D 结构）和**置信度**，**不是结合亲和力**（ΔG、Kd）。亲和力评估需另配打分函数（GNINA、MM/GBSA）。

## 前置要求：NVIDIA API Key（首次使用时提供一次）

⚠️ **本技能全部计算在 NVIDIA 云端完成，需要一个 API key。**

获取方式：登录 [https://build.nvidia.com](https://build.nvidia.com) → 找到 DiffDock 模型 → 生成 API key（形如 `nvapi-xxxxxxxx`）。

### Key 如何工作（重要）

**首次运行**任何对接脚本时，如果还没提供过 key，脚本会**交互式提示用户粘贴 key**，并把它保存到 skill 目录**之外**的用户配置文件：

```
~/.config/ZYDock/credentials   （权限 600，仅本人可读）
```

**之后所有运行都自动读取这个保存的 key，无需再次输入。**

Key 的解析优先级：
1. 环境变量 `NVIDIA_API_KEY`（临时覆盖用）
2. 保存的凭证文件 `~/.config/ZYDock/credentials`
3. 首次运行的交互式提示（自动保存）

### 管理 key

```bash
python scripts/set_key.py nvapi-xxxxxxxx   # 主动保存/更换 key
python scripts/set_key.py --show           # 查看当前 key 状态（打码显示）
python scripts/set_key.py --clear          # 删除已保存的 key
```

### 🔒 关键安全设计

- **这个 skill 目录本身绝不包含任何 key** —— key 只存在用户配置目录 `~/.config/ZYDock/`，与 skill 分离。因此 skill 可以安全地分享给学生/同事，每个人用自己的 key。
- 凭证文件权限设为 600（仅本人可读）。
- 不要把 key 明文贴进聊天或提交到 git。如泄露，去 build.nvidia.com 吊销重建。

> **给 AI 的提示：** 若用户尚未提供过 key（`set_key.py --show` 显示未保存、且环境变量为空），在跑对接前先向用户索取 NVIDIA API key，拿到后用 `python scripts/set_key.py <key>` 保存一次即可。

## 何时使用本技能 When to Use

当用户表达以下需求时使用：
- "把这个配体对接到蛋白上" / "预测结合位姿"
- "跑分子对接" / "蛋白-配体对接"
- "虚拟筛选" / "筛一批化合物"
- "这个分子结合在哪里" / "预测结合位点"
- 强调"不想用本地 GPU""用在线/云端算力""用英伟达 API"对接
- 手上有 PDB 文件/蛋白序列 + SMILES/配体结构

## 对接需要提供什么

| 输入 | 格式 | 说明 |
|---|---|---|
| **① NVIDIA API key** | `nvapi-...` | 必需，首次使用时提供一次 |
| **② 蛋白（受体）** | **基因/蛋白名** 或 PDB 文件 或 氨基酸序列 | 三选一 |
| **③ 小分子（配体）** | SMILES 字符串 或 SDF 文件 | 二选一 |

**蛋白输入的三种方式：**
- **只给名字（推荐，最省事）**：给基因名（如 `EGFR`）或蛋白名，技能自动从 UniProt 检索并下载结构（见下方"自动获取结构"）
- 给现成 PDB 文件
- 给氨基酸序列（API 内部用 ESMFold 折叠）

**配体从哪来：** 从 PubChem 按分子名查 Canonical SMILES；或用 RDKit 从 SMILES 生成 SDF。

**无需指定结合口袋**——DiffDock 是盲对接，自动寻找结合位点。

## 自动获取蛋白结构（按基因/蛋白名）

只要给基因名或蛋白名，技能会自动完成"检索 → 选结构 → 下载"，无需手动找 PDB。

**结构选择原则：**
1. 在 UniProt 检索该名称（**默认人源** organism_id=9606，可通过 `--organism` 指定其他物种），优先 reviewed（Swiss-Prot）条目，得到 UniProt accession
2. 在该条目的 PDB 交叉引用中，**优先已报道的实验解析结构**（X-ray / EM / NMR），并**按分辨率从高到低排序**（分辨率数值最小者最优），选出最佳结构从 RCSB 下载
3. 若无可用实验结构，**回退到 AlphaFold 预测模型**（会标注为预测结构，提示谨慎解读）

**单独获取结构（不对接）：**
```bash
# 列出某基因的实验结构（按分辨率排序）
python scripts/fetch_structure.py --gene EGFR --list

# 下载最优结构
python scripts/fetch_structure.py --gene EGFR --out egfr.pdb

# 指定物种（如小鼠 10090）
python scripts/fetch_structure.py --gene Trp53 --organism 10090
```

**直接给名字对接（一步到位）：**
```bash
python scripts/dock_single.py --gene EGFR \
  --ligand "Cc1ccc(cc1Nc2nccc(n2)c3cccnc3)NC(=O)c4ccc(cc4)CN5CCN(CC5)C" \
  --out_dir results/egfr_imatinib/
```

## API 技术细节

- **端点：** `POST https://health.api.nvidia.com/v1/biology/mit/diffdock`
- **认证头：** `Authorization: Bearer $NVIDIA_API_KEY`
- **关键点：** 请求体里 `protein` 字段传的是 **PDB 文件的文本内容本身**（不是路径）
- **蛋白自动清理：** 对接前脚本会自动去除 PDB 中所有 HETATM 记录（水分子、离子、结晶添加剂、共结晶配体），只保留纯蛋白发送给 DiffDock。这避免了杂原子干扰结合位点预测。清理对用户透明，无需手动操作。
- **请求体字段：**
  - `protein`：PDB 文本内容
  - `ligand`：SMILES 串 或 SDF 文本
  - `ligand_file_type`：`"txt"`（SMILES）或 `"sdf"`
  - `num_poses`：生成位姿数（默认 20）
  - `steps`、`time_divisions`：扩散采样参数
- **响应字段：** `status`、`ligand_positions`（每个位姿的 SDF 文本）、`position_confidence`（对应置信度）、`trajectory`、`protein`、`ligand`

## 核心工作流 Core Workflows

### 工作流 0：环境自检（推荐先跑）

```bash
python scripts/check_api.py
```
检查 `NVIDIA_API_KEY` 是否设置、网络是否可达、API 是否响应正常。

### 工作流 1：单个配体对接

**方式 A（推荐）：只给基因/蛋白名，自动获取结构**
```bash
python scripts/dock_single.py \
  --gene EGFR \
  --ligand "CC(=O)Oc1ccccc1C(=O)O" \
  --out_dir results/single/ \
  --num_poses 20
```

**方式 B：用已有 PDB 文件**
```bash
python scripts/dock_single.py \
  --protein protein.pdb \
  --ligand "CC(=O)Oc1ccccc1C(=O)O" \
  --out_dir results/single/ \
  --num_poses 20
```

若没有 PDB 只有序列：
```bash
python scripts/dock_single.py \
  --protein_sequence "MSKGEELFTGVVPILVEL..." \
  --ligand "CC(=O)Oc1ccccc1C(=O)O" \
  --out_dir results/single/
```

**输出：** `results/single/` 下按置信度排好序的 `rank1_conf*.sdf` … + `confidence_scores.txt` + 原始 `raw_response.json`。

### 工作流 2：批量虚拟筛选

准备 CSV（见 `assets/batch_template.csv`）。蛋白可用 `gene` 列自动获取，也可用 `protein_path` 指定文件：
```csv
complex_name,gene,protein_path,ligand_smiles,protein_sequence,organism
egfr_aspirin,EGFR,,CC(=O)Oc1ccccc1C(=O)O,,
egfr_caffeine,EGFR,,Cn1cnc2c1c(=O)n(C)c(=O)n2C,,
byfile_cmpd,,protein.pdb,COc1ccc(C#N)cc1,,
```
> 同一个 `gene` 只会下载一次结构并复用；`gene` 为空时用 `protein_path` 或 `protein_sequence`。`organism` 留空默认人源。

运行批量对接：
```bash
python scripts/dock_batch.py --csv my_input.csv --out_dir results/batch/ --num_poses 10
```
逐个调用 API，每个化合物结果存到独立子目录，并生成一张汇总表 `summary.csv`（含每个化合物的最优 confidence）。

### 工作流 3：分析与排序结果

```bash
python scripts/analyze_results.py results/batch/
python scripts/analyze_results.py results/batch/ --export ranking.csv
```
解析所有 confidence 分数，按 High(>0)/Moderate(-1.5~0)/Low(<-1.5) 分级，跨化合物排名。

### 工作流 4：相互作用分析与绘图（PLIP + PyMOL）

对接得到位姿后，用 PLIP 分析蛋白-配体的相互作用（氢键、疏水、π-堆积、盐桥等），并用 PyMOL 出图。

**依赖环境（首次自动准备）：** 本工作流在专用 conda 环境 `PLIP` 中运行。

```bash
# 检测/创建 PLIP 环境(已存在则跳过创建)
python scripts/setup_plip.py
```
`setup_plip.py` 会检测名为 `PLIP` 的 conda 环境是否存在：存在就直接用；不存在则自动创建（Python 3.10 + imagemagick + pymol-open-source + plip）。

**分析相互作用（自动合并复合物 → 跑 PLIP → PyMOL 出图）：**
```bash
# 分析对接的最优位姿(蛋白用文件)
python scripts/analyze_interactions.py \
  --protein protein.pdb \
  --ligand results/single/rank1_conf-0.150.sdf \
  --out_dir results/interactions/

# 蛋白也可用基因名自动获取
python scripts/analyze_interactions.py \
  --gene EGFR \
  --ligand results/single/rank1_conf-0.150.sdf \
  --out_dir results/interactions/

# 只做 PLIP 分析、不出 PyMOL 图
python scripts/analyze_interactions.py --protein protein.pdb \
  --ligand pose.sdf --out_dir results/interactions/ --no-pymol
```

**脚本自动完成五步：**
1. **清理受体 + 合并复合物**——先去除受体 PDB 中原有的所有杂原子（HETATM：结晶添加剂、共结晶配体、水、离子），只保留纯蛋白；再把配体 SDF 转为 HETATM（残基名 UNL、链 Z）合并进来，得到 `complex.pdb`。合并用 openbabel 自动完成。
2. **PLIP 分析**——运行 `plip -f complex.pdb -ypxt`，输出 XML、文本报告、PyMOL 会话 `.pse`。
3. **PyMOL 绘图（两张 300 DPI 图）**——加载 PLIP 的 `.pse` 和复合物，光追渲染：
   - `complex_overview.png` —— 整体图（不透明 cartoon，看配体在蛋白全貌中的位置）
   - `interaction_bindingsite.png` —— 结合位点特写（半透明 cartoon + 相互作用线 + 关键残基标签）
4. **清理中间文件**——删除 PLIP 自带小图和中间产物。
5. **自动生成 SCI 级报告**——见"工作流 5"，输出到并列的 `report/` 文件夹。

> **为什么要清理受体：** 实验解析的 PDB 结构常带有结晶缓冲剂、共结晶配体、水和离子。若不清理，PLIP 会把它们也当作配体逐个分析。清理后只保留纯蛋白 + 我们对接的配体（UNL:Z）。

**输出文件（interactions/ 文件夹）：**
- `complex.pdb` —— 合并的纯蛋白-配体复合物
- `complex_report.txt` —— 文本相互作用报告（氢键/疏水/盐桥/卤键，含残基和距离）
- `complex_report.xml` —— 结构化 XML
- `*.pse` —— PyMOL 会话
- `complex_overview.png` —— 整体图（300 DPI）
- `interaction_bindingsite.png` —— 结合位点特写图（300 DPI，含残基标签）

### 工作流 5：生成 SCI 级中英文分析报告

`analyze_interactions.py` 跑完会**自动**生成报告（除非加 `--no-report`），读取整个流程的计算参数（结构来源/分辨率/对接参数/置信度/PLIP 相互作用），生成可发表级别的**材料方法 + 结果分析**，中英文双语，输出到与 `docking/`、`interactions/` **并列的 `report/` 文件夹**。

也可单独运行（对已完成目录补生成）：
```bash
python scripts/generate_report.py --run_dir results/demo --gene EGFR --accession P00533
```

**report/ 输出：**
- `report_CN.md` —— 中文报告（材料与方法 + 结果 + 讨论 + 参考文献）
- `report_EN.md` —— 英文报告（Materials and Methods + Results + Discussion + References）
- `interactions_summary.csv` —— 相互作用汇总表（机器可读）

**报告包含的深度分析：**
- **材料与方法**：受体来源/分辨率、配体处理、DiffDock 盲对接原理与参数、PLIP + PyMOL 流程，含文献引用 [1–3]
- **结果**：置信度分布定量分析（各等级构象数、均值、跨度、top 构象一致性判断）、相互作用类型计数与主导作用力、**结合热点残基分析**（多接触残基自动识别）、相互作用表、图示
- **讨论**：置信度分级解读、主导相互作用的生物物理意义、优化着力点、局限性与后续实验建议（MD 精修 / MM-GBSA 重打分 / SPR·ITC 验证）
- **参考文献**：DiffDock、PLIP、PyMOL、PDB、UniProt、AlphaFold（国际通行著录格式）

## 置信度分数解读 Confidence Interpretation

| 分数范围 | 置信等级 | 含义 |
|---|---|---|
| **> 0** | High 高 | 强预测，很可能准确 |
| **-1.5 ~ 0** | Moderate 中 | 合理预测，需谨慎验证 |
| **< -1.5** | Low 低 | 不确定，需要验证 |

**注意：** 置信度 ≠ 亲和力。高置信度表示模型对结构的确定性高，不代表结合强。大配体（>500 Da）、多链蛋白、新颖蛋白家族通常置信度偏低。建议看 top 3-5 位姿找共识。

## 参数调整 Parameters

- `--num_poses`：位姿数，难例可增到 30-40
- `--steps`：采样步数（默认 18-20），增大可能提精度但更慢
- `--time_divisions`：时间划分（默认 20）

详见 `references/api_reference.md`。

## 局限性 Limitations

**适用：** 小分子配体（100-1000 Da）、类药有机物、小肽（<20 残基）、单/多链蛋白。

**不适用：**
- 蛋白-蛋白对接 → 用 DiffDock-PP 或 AlphaFold-Multimer
- 大肽（>20 残基）
- 共价对接
- 结合亲和力预测 → 需配打分函数
- 膜蛋白（未专门训练，谨慎使用）

## 故障排查 Troubleshooting

- **报错 "set NVIDIA_API_KEY"**：没设环境变量，先 `export NVIDIA_API_KEY=nvapi-...`
- **HTTP 401/403**：key 无效或已过期，去 build.nvidia.com 重新生成
- **HTTP 422**：请求字段格式错，检查 `ligand_file_type` 与实际配体是否匹配
- **网络超时**：检查能否访问 health.api.nvidia.com；大蛋白上传较慢，脚本已设 300s 超时
- **全部低置信度**：配体过大/结合位点不清/蛋白结构差，尝试增大 `--num_poses`、清理蛋白结构

## 资源清单 Resources

### scripts/
- `set_key.py` —— 保存/查看/删除 NVIDIA API key（存到 skill 之外的用户配置目录）
- `check_api.py` —— 检查 API key、网络、端点连通性
- `fetch_structure.py` —— 按基因/蛋白名从 UniProt 检索、选最优实验结构（按分辨率）并下载，无实验结构则回退 AlphaFold
- `dock_single.py` —— 单个蛋白-配体对接，结果落盘为排序 SDF
- `dock_batch.py` —— 批量/虚拟筛选，读 CSV 逐个对接并汇总
- `analyze_results.py` —— 解析置信度、分级、跨化合物排名、导出 CSV
- `setup_plip.py` —— 检测/创建 PLIP conda 环境（含 plip + pymol-open-source + openbabel）
- `analyze_interactions.py` —— 合并蛋白+配体 → PLIP 分析 → PyMOL 两张图 → 自动生成报告
- `generate_report.py` —— 读取全流程参数，生成 SCI 级中英文材料方法+结果分析报告

### references/
- `api_reference.md` —— API 端点、完整请求/响应字段、参数详解
- `getting_structures.md` —— 如何获取蛋白结构（RCSB/AlphaFold/序列）和配体 SMILES（PubChem/RDKit）

### assets/
- `batch_template.csv` —— 批量对接输入模板

## 参考 Citations

- DiffDock-L: Stärk et al. (2024), arXiv:2402.18396
- DiffDock: Corso et al. (2023), ICLR 2023, arXiv:2210.01776
- NVIDIA NIM DiffDock: https://build.nvidia.com/mit/diffdock
