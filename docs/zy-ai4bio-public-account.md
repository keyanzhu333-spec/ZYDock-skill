# ZY-AI4Bio：面向生命科学与药物发现的一站式 AI 工作流

GitHub：<https://github.com/keyanzhu333-spec/ZYDock-skill/tree/codec>

在药物发现和生物信息学研究中，研究者往往需要在蛋白质数据库、化合物数据库、分子对接、虚拟筛选、结构分析和文献平台之间反复切换。ZY-AI4Bio 将这些任务组织成一套可以用自然语言调用的 AI for Biology 工作流。

## 结构生物学与分子对接

这是 ZY-AI4Bio 的核心能力之一，覆盖从靶点结构准备到候选分子分析的连续流程。

### 1. 靶点与结构准备

可以从 UniProt、RCSB PDB 或 AlphaFold 获取蛋白信息，下载序列和结构，并检查 PDB、mmCIF、SDF、MOL2、PDBQT 等常用格式。对于受体结构，还可以完成基础清理、链选择和对接前准备。

### 2. DiffDock 蛋白-小分子盲对接

通过 NVIDIA DiffDock API 预测蛋白与小分子的潜在结合姿态，适合没有明确口袋、需要快速探索结合模式的任务。输出包括候选姿态和模型置信度，可用于单个配体或少量配体的初步研究。置信度是姿态预测信号，不等同于实验亲和力。

### 3. UniDock-Pro 大规模虚拟筛选

对于数千到数十万个已经准备好的 PDBQT 小分子，调用本地 GPU 上的 UniDock-Pro 进行批量筛选。当前本地库已经整理了 FDA/DrugCentral 药物、天然产物、海洋天然产物、微生物天然产物、中药和药用植物相关成分、食品化学、脂质、化学探针、药物再利用、Tox21 和风味分子等专题资源。

筛选流程保留原始结构、标准化 SDF/SMILES、PDBQT、统一化合物主表、来源映射、许可证信息和失败日志；每个活动库提供 `ligand_index.txt`，可以直接接入 UniDock-Pro。

### 4. GNINA/GVINA 风格重新打分

初筛后可用 GNINA 的 CNN 模型进行重新打分，分别保留传统对接分数、CNNscore 和 CNNaffinity，避免将不同模型的结果混成一个不可解释的分数。

### 5. PyMOL 与 PLIP 相互作用分析

对排名靠前的复合物，可以继续查看氢键、疏水作用、盐桥、卤键、π-π 相互作用和关键残基接触，并生成结构图、相互作用报告和候选排名表。

### 6. 对接后的动态验证

如果需要进一步评价复合物稳定性，可以衔接 GROMACS，完成小分子参数化、能量最小化、平衡、生产模拟、RMSD/RMSF、氢键和 MM/PBSA 等分析。蛋白-蛋白体系则可以使用 LightDock；抗体、Fab 和纳米抗体可先用 IgFold 建模。

## 其他能力

除结构生物学与分子对接外，ZY-AI4Bio 还可以协调：

- 蛋白、序列和结构查看，包括 UniProt、PDB、AlphaFold、FASTA/FASTQ 和多序列比对；
- PubChem、ChEMBL、BindingDB、ChEBI、DrugCentral、Open Targets 等化合物、药物和靶点数据库；
- Ensembl、ClinVar、gnomAD、GWAS、GTEx、HPA、cBioPortal、STRING、Reactome 和 QuickGO 等基因、变异、表达、网络和通路资源；
- PRIDE、ProteomeXchange、MetaboLights、MGnify、BioStudies/ArrayExpress 和 NCBI 等组学数据集；
- PubMed、PMC、bioRxiv、medRxiv 与 RSS/Atom 文献检索、追踪和总结。

## 一句话总结

ZY-AI4Bio 不是单一的软件，而是一个面向生命科学研究的工作流入口：从靶点和结构查询，到小分子库准备、虚拟筛选、对接复核、相互作用分析和分子动力学，再到数据库与文献证据整理，帮助研究者把分散的工具连接成可复用的研究流程。

> 所有对接分数、模型置信度和动力学结果都属于计算预测，应与数据库事实、文献证据和实验结果明确区分。