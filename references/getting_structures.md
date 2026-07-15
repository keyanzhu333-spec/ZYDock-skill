# 如何获取蛋白结构和配体结构

对接需要两样东西：**蛋白（受体）** 和 **小分子（配体）**。

## 零、最省事：让技能按名字自动获取蛋白结构（推荐）

如果你只知道基因名或蛋白名，直接用 `fetch_structure.py`（或在对接脚本里用 `--gene`），
技能会自动完成"UniProt 检索 → 选最优结构 → 下载"：

```bash
python scripts/fetch_structure.py --gene EGFR --list    # 先看有哪些实验结构
python scripts/fetch_structure.py --gene EGFR --out egfr.pdb   # 下载最优
```

**结构选择原则：**
1. UniProt 检索该名称（默认人源 9606，`--organism` 可改物种），优先 reviewed 条目
2. 优先实验解析结构（X-ray/EM/NMR），按分辨率从高到低排序，选分辨率数值最小者
3. 无实验结构则回退 AlphaFold 预测模型（标注为预测，谨慎解读）

下面是各来源的手动获取途径（供理解原理或特殊情况使用）。

## 一、蛋白结构从哪来

### 途径 A：从 RCSB PDB 下载实验结构（最常用）

知道 PDB ID（4 位，如 1IEP）即可直接下载：

```bash
curl -s https://files.rcsb.org/download/1IEP.pdb -o protein.pdb
```

去 https://www.rcsb.org 按蛋白名/基因名搜索，找到 4 位 ID。

### 途径 B：AlphaFold 预测结构（没有实验结构时）

去 https://alphafold.ebi.ac.uk ，按 UniProt ID 下载预测好的 PDB。

### 途径 C：只有序列

DiffDock API 支持直接喂氨基酸序列，内部用 ESMFold 现折叠。
用脚本的 `--protein_sequence` 参数，或 CSV 里填 `protein_sequence` 列。

### ⚠️ 蛋白准备建议

- 去掉多余水分子、无关配体/离子
- 补齐缺失残基
- 结构不干净会导致对接出不合理位姿

## 二、配体（小分子）结构从哪来

### 途径 A：直接用 SMILES 字符串（最省事）

SMILES 是一串文本就能表示分子，DiffDock 直接吃。

获取 SMILES：
- **PubChem**（https://pubchem.ncbi.nlm.nih.gov）：搜分子名 → 复制 Canonical SMILES
- **ChEMBL / DrugBank**：药物类分子
- ChemDraw 等结构编辑器画完导出

常见示例：
| 分子 | SMILES |
|---|---|
| 阿司匹林 | `CC(=O)Oc1ccccc1C(=O)O` |
| 咖啡因 | `Cn1cnc2c1c(=O)n(C)c(=O)n2C` |
| 伊马替尼 | `Cc1ccc(cc1Nc2nccc(n2)c3cccnc3)NC(=O)c4ccc(cc4)CN5CCN(CC5)C` |

### 途径 B：用 SDF/MOL2 结构文件

已有 3D 结构文件直接传（脚本会自动识别扩展名设 `ligand_file_type="sdf"`）。

用 RDKit 从 SMILES 生成 SDF：
```python
from rdkit import Chem
from rdkit.Chem import AllChem
m = Chem.AddHs(Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O"))
AllChem.EmbedMolecule(m)
Chem.MolToMolFile(m, "ligand.sdf")
```

## 三、DiffDock 不需要你提供的东西

和传统对接（AutoDock Vina）不同：
- **不需要指定结合口袋 / grid box 坐标** —— DiffDock 是盲对接，自动找结合位点。

## 四、最短流程

```
蛋白：知道 PDB ID → RCSB 下 .pdb ；没有 → AlphaFold 或直接给序列
配体：知道分子名 → PubChem 抄 SMILES ；直接把 SMILES 交给脚本
→ 设好 NVIDIA_API_KEY → 调 dock_single.py / dock_batch.py
```
