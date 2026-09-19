# ZYDock

**Cloud-based molecular docking + interaction analysis for Codex, all from a gene name and a SMILES string.**

ZYDock is a Codex skill that runs the full structure-based drug discovery pipeline. No local GPU or DiffDock installation is required. All docking computation runs on NVIDIA's hosted GPUs via the [NVIDIA NIM DiffDock API](https://build.nvidia.com/mit/diffdock).

## What it does

Give it a gene name and a ligand SMILES, and it automatically runs:

```
Gene name (e.g. EGFR)
    ↓  fetch structure from UniProt (prefers highest-resolution experimental structure)
Protein structure (auto-downloaded)
    ↓  clean heteroatoms, dock via NVIDIA cloud API
Binding poses + confidence scores
    ↓  merge complex, run PLIP interaction analysis
Interactions + two 300-DPI PyMOL figures
    ↓  parse all parameters
Bilingual (EN/CN) SCI-grade analysis report
```

## ZY-AI4Bio extension

The repository also includes the `zy-ai4bio` workflow router for combining docking,
UniDock-Pro virtual screening, GNINA/PLIP analysis, molecular dynamics, structure
and sequence viewers, public life-science databases, and literature workflows.

- [ZY-AI4Bio public-account article](docs/zy-ai4bio-public-account.md)
- [Virtual-screening library status](docs/virtual-screening-library.md)
- [ZY-AI4Bio skill](zy-ai4bio/SKILL.md)

## Key features

- **Codex-ready** — includes `SKILL.md` and `agents/openai.yaml` so Codex can discover and invoke the skill
- **No local compute** — docking runs on NVIDIA's hosted GPUs; you only need an API key
- **Auto structure retrieval** — give a gene/protein name; it searches UniProt and picks the best experimentally-solved structure (highest resolution), falling back to AlphaFold if none exists (defaults to human, configurable species)
- **Automatic receptor cleaning** — strips crystallographic waters, ions, buffer additives and co-crystallized ligands before docking
- **Interaction profiling** — PLIP detects hydrogen bonds, hydrophobic contacts, salt bridges, halogen bonds, π-stacking, etc.
- **Publication-quality figures** — two 300-DPI ray-traced PyMOL images (whole-complex overview + labelled binding-site close-up)
- **Bilingual SCI reports** — auto-generated Materials & Methods + Results in both English and Chinese
- **Batch virtual screening** — one CSV → full pipeline for every compound + a ranked master summary
- **Secure by design** — the API key is stored *outside* the skill directory (`~/.config/ZYDock/`, chmod 600); the repo contains no keys and is safe to share

## Requirements

- An [NVIDIA API key](https://build.nvidia.com) (free to obtain; `nvapi-...`)
- Python 3
- For interaction analysis/figures: a conda environment named `PLIP` (auto-created by `scripts/setup_plip.py`) with `plip`, `pymol-open-source`, and `openbabel`

## Quick start

```bash
# 1. Save your NVIDIA API key (stored outside the skill, only needed once)
python scripts/set_key.py nvapi-xxxxxxxx

# 2. Dock a ligand to a target by gene name
python scripts/dock_single.py --gene EGFR \
  --ligand "COc1cc2ncnc(Nc3ccc(F)c(Cl)c3)c2cc1OCCCN1CCOCC1" \
  --out_dir results/docking --num_poses 10

# 3. Analyze interactions + render figures + generate report
python scripts/analyze_interactions.py \
  --protein results/docking/EGFR_structure.pdb \
  --ligand results/docking/rank1_conf*.sdf \
  --out_dir results/interactions --gene_label EGFR --accession P00533
```

### Batch virtual screening (end-to-end)

```bash
python scripts/batch_full.py --csv compounds.csv --out_dir results/screen/ --num_poses 10
```

CSV columns: `complex_name, gene, protein_path, ligand_smiles, protein_sequence, organism, accession`
(see [`assets/batch_template.csv`](assets/batch_template.csv)). The ligand only needs a **SMILES** string.

## Scripts

| Script | Purpose |
|--------|---------|
| `set_key.py` | Save/show/clear the NVIDIA API key (stored outside the skill) |
| `check_api.py` | Verify API connectivity |
| `fetch_structure.py` | Gene/protein name → UniProt → best experimental structure → download |
| `dock_single.py` | Single protein–ligand docking |
| `dock_batch.py` | Batch docking (docking only) |
| `analyze_results.py` | Confidence-score analysis and ranking |
| `setup_plip.py` | Detect/create the `PLIP` conda environment |
| `analyze_interactions.py` | Merge complex → PLIP → PyMOL figures → report |
| `generate_report.py` | Generate the bilingual SCI-grade report |
| `batch_full.py` | End-to-end batch pipeline for every compound |

## Important notes

- DiffDock predicts **binding poses** and **confidence**, **not binding affinity** (ΔG/Kd). For affinity, combine with a scoring function (GNINA, MM/GBSA).
- Confidence scores reflect model certainty about the pose, not binding strength.
- Predicted poses should be experimentally validated.

## Security

Never commit API keys. This repository contains only placeholders (`nvapi-xxxxxxxx`). Your real key lives in `~/.config/ZYDock/credentials` (outside the repo, chmod 600). If a key is ever exposed, revoke and regenerate it at [build.nvidia.com](https://build.nvidia.com).

## License

MIT