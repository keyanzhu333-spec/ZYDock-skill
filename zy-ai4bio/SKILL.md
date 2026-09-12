---
name: zy-ai4bio
description: AI-for-biology workflow router that coordinates installed Codex skills for molecular docking, GNINA/GVINA-style CNN rescoring, molecular dynamics, antibody structure prediction, protein and sequence visualization, life-science databases, and biomedical literature. Use when the user asks for a broad bioinformatics, computational biology, structural biology, drug discovery, docking-to-MD, antibody modeling, database/literature lookup, or multi-step AI4Bio workflow and wants Codex to choose the right skill(s).
---

# ZY-AI4bio

Use this as the top-level coordinator for biology and drug-discovery tasks. Route to the most specific installed skill first; use this skill to combine them when a task spans multiple domains.

## Routing

Read `references/skill-map.md` when the task is broad, ambiguous, or involves more than one biology workflow.

- For protein-small molecule blind docking with NVIDIA DiffDock API, use `diffdock`.
- For prepared receptor-ligand docking with GNINA CNN rescoring, PyMOL visualization, and PLIP interaction analysis, use `gnina-dock`.
- For protein-protein docking, use `lightdock`.
- For docking follow-up MD, trajectory analysis, or MM/PBSA, use `gromacs-md`.
- For antibody, Fab, or nanobody structure prediction, use `igfold`.
- For local molecular structure inspection or rendering, use `structure-viewer:structure-viewer`.
- For FASTA, FASTQ, GenBank, alignments, variants, or chromatograms, use `sequence-viewer:biological-sequence-viewer`.
- For public biological database lookups, use the relevant `life-sciences-databases:*` skill.
- For RSS/Atom literature monitoring, daily research alerts, and topic-filtered paper updates, use `literature-skill`.
- For biomedical literature and PubMed/PMC/preprint lookups, use the relevant `life-sciences-literature:*` skill.

## Workflow Pattern

1. Identify the biological object: gene/protein, ligand, sequence, structure, antibody chain, variant, sample set, study, or dataset.
2. Choose the narrowest child skill that can complete the immediate task.
3. When the request spans steps, chain skills explicitly and keep artifacts organized by stage.
4. Preserve each child skill's own safety, credential, citation, and validation instructions.
5. Report which child skill(s) were used and where outputs were saved.

## Common Chains

- Target plus compound: `life-sciences-databases:uniprot-skill` or `life-sciences-databases:pubchem-pug-skill` for identifiers, then `diffdock` for blind DiffDock docking or `gnina-dock` when prepared structures and a binding box/reference ligand are available, then `gromacs-md` for MD if requested.
- Protein-protein mechanism: database lookup, structure retrieval via `life-sciences-databases:rcsb-pdb-skill` or local files, then `lightdock`.
- Antibody workflow: `igfold` for VH/VL modeling, `structure-viewer:structure-viewer` for inspection, then a docking skill if the user asks for interaction modeling.
- Evidence workflow: database skills for structured facts, literature skills for papers, `literature-skill` for RSS monitoring or daily alerts, then summarize evidence levels and limitations without overclaiming.

## Output Discipline

- Do not merge unrelated outputs into one folder. Prefer stage folders such as `docking/`, `interactions/`, `report/`, `md/`, `figures/`, `tables/`, and `literature/`.
- Do not expose or store API keys in skill folders, output folders, command text, reports, or Git history.
- For scientific claims, distinguish prediction, annotation, database fact, literature evidence, and experimental validation.
