# ZY-AI4Bio 虚拟筛选库

这是与 ZY-AI4Bio 配套的本地虚拟筛选库说明。完整 SDF、SMILES、PDBQT、SQLite 和日志文件保存在本地数据目录，不随技能仓库提交大体积原始数据。

本地目录：`/media/keyanzhu/data/others/ZY-AI4Bio_VScreen_Library`

## 当前规模

- 统一化合物记录：598,214 条
- 已登记筛选/标准化库：38 个
- 活动筛选目录：26 个
- 可直接用于 UniDock-Pro 的 PDBQT：319,630 个
- COCONUT 多样性子库：8,822 个 PDBQT

## 主要专题

FDA/DrugCentral 药物、Broad Drug Repurposing Hub、NPASS、NPAtlas、CMNPD 海洋天然产物、COCONUT、Dr. Duke 药用植物、FooDB 食品化学、FlavorDB、LIPID MAPS、ChEBI、Probes & Drugs、NIH Tox21，以及微生物、脂质、糖脂和特色天然产物子库。

## 文件和质量控制

每个库尽量同时保留原始下载文件、规范化 SDF/SMILES、筛选级 PDBQT、化合物主表、来源映射、许可证/版本信息和失败结构日志。COCONUT 子库采用全量去重、最大有机片段提取、元素/电荷和分子量过滤、MW/TPSA 分层、Morgan 指纹 MaxMin 多样性选择及固定随机种子，避免简单截取前若干条记录。

UniDock-Pro 运行项目仍使用三个顶层目录：`inputs/`、`docking/`、`report/`。筛选前应检查受体 PDBQT、配体目录和 `ligand_index.txt`，并将对接评分视为候选排序信号，而不是实验亲和力。

## 本地库与技能仓库

- 技能仓库：[ZYDock-skill / codec](https://github.com/keyanzhu333-spec/ZYDock-skill/tree/codec)
- 根技能：[zy-ai4bio/SKILL.md](../zy-ai4bio/SKILL.md)
- 本地库目录：`/media/keyanzhu/data/others/ZY-AI4Bio_VScreen_Library`