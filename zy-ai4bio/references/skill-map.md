# ZY-AI4bio Skill Map

Use this map to choose the right installed child skill. Prefer the narrowest skill that satisfies the request. Chain multiple skills only when the user asks for a multi-step workflow or when one stage naturally produces the next stage's input.

## Structural Biology And Drug Discovery

| Need | Use skill | Notes |
|---|---|---|
| Protein-small molecule docking, blind binding pose prediction, DiffDock via NVIDIA NIM | `diffdock` | Requires NVIDIA API key. Outputs docking poses, confidence, PLIP interactions, PyMOL images, and bilingual reports. |
| Prepared receptor-ligand docking with GNINA/GVINA-style workflow, CNN rescoring, pose ranking, PyMOL images, PLIP reports | `gnina-dock` | Use when receptor/ligand structures and a binding box or reference ligand are available. Keeps GNINA scores, CNNscore, and CNNaffinity separate. |
| Protein-protein docking, PPI complex prediction, receptor.pdb plus ligand.pdb | `lightdock` | Local LightDock workflow; CPU-based and can take time. |
| Protein-ligand MD, ligand parameterization, GROMACS workflow, trajectory analysis, MM/PBSA | `gromacs-md` | Use after docking when the user wants stability, dynamics, or binding free-energy post-processing. |
| Antibody/Fab/nanobody structure from VH/VL sequences | `igfold` | Uses local igfold-gpu environment when available. |
| Inspect, style, compare, animate, or render molecular structures | `structure-viewer:structure-viewer` | Use for PDB, CIF/mmCIF, SDF/MOL2, PDBQT, GRO, XYZ, topology and trajectory-linked views. |

## Sequence

| Need | Use skill | Notes |
|---|---|---|
| FASTA/FASTQ, GenBank, alignments, chromatograms, sequence edits or exports | `sequence-viewer:biological-sequence-viewer` | Use for local sequence files and visual sequence workbench tasks. |

## Public Databases

| Need | Use skill |
|---|---|
| UniProt accessions, protein function, FASTA | `life-sciences-databases:uniprot-skill` |
| PDB structure metadata, Search API, FASTA downloads | `life-sciences-databases:rcsb-pdb-skill` |
| AlphaFold DB prediction and metadata summaries | `life-sciences-databases:alphafold-skill` |
| PubChem compound properties, descriptions, assays | `life-sciences-databases:pubchem-pug-skill` |
| ChEMBL molecules, targets, activities, mechanisms | `life-sciences-databases:chembl-skill` |
| BindingDB ligand-target binding by PDB, UniProt, or similarity | `life-sciences-databases:bindingdb-skill` |
| ChEBI compound ontology and structure metadata | `life-sciences-databases:chebi-skill` |
| Open Targets target, disease, drug, variant evidence | `life-sciences-databases:opentargets-skill` |
| STRING networks, partners, enrichment | `life-sciences-databases:string-skill` |
| Reactome pathways and events | `life-sciences-databases:reactome-skill` |
| QuickGO terms and annotations | `life-sciences-databases:quickgo-skill` |
| Ensembl lookup, overlaps, xrefs, variations | `life-sciences-databases:ensembl-skill` |
| ClinVar/Variation variant summaries | `life-sciences-databases:clinvar-variation-skill` |
| gnomAD frequency and constraint | `life-sciences-databases:gnomad-graphql-skill` |
| GWAS Catalog studies, associations, SNPs, traits | `life-sciences-databases:gwas-catalog-skill` |
| GTEx eQTL | `life-sciences-databases:gtex-eqtl-skill` |
| eQTL Catalogue | `life-sciences-databases:eqtl-catalogue-skill` |
| FinnGen, UKB-TOPMed, BBJ, TPMI PheWAS | `life-sciences-databases:finngen-phewas-skill`, `life-sciences-databases:ukb-topmed-phewas-skill`, `life-sciences-databases:biobankjapan-phewas-skill`, `life-sciences-databases:tpmi-phewas-skill` |
| Human Protein Atlas | `life-sciences-databases:human-protein-atlas-skill` |
| Bgee expression | `life-sciences-databases:bgee-skill` |
| cBioPortal cancer genomics | `life-sciences-databases:cbioportal-skill` |
| CIViC cancer variant evidence | `life-sciences-databases:civic-skill` |
| PharmGKB pharmacogenomics | `life-sciences-databases:pharmgkb-skill` |
| ClinicalTrials.gov | `life-sciences-databases:clinicaltrials-skill` |
| ENCODE | `life-sciences-databases:encode-skill` |
| CELLxGENE Discover | `life-sciences-databases:cellxgene-skill` |
| BioStudies and ArrayExpress | `life-sciences-databases:biostudies-arrayexpress-skill` |
| MGnify microbiome metadata | `life-sciences-databases:mgnify-skill` |
| MetaboLights | `life-sciences-databases:metabolights-skill` |
| PRIDE or ProteomeXchange proteomics | `life-sciences-databases:pride-skill`, `life-sciences-databases:proteomexchange-skill` |
| NCBI Datasets, Entrez, Clinical Tables gene autocomplete | `life-sciences-databases:ncbi-datasets-skill`, `life-sciences-databases:ncbi-entrez-skill`, `life-sciences-databases:ncbi-clinicaltables-skill` |
| EVA archived variants, IPD HLA, EFO, Rhea, RNAcentral, EpiGraphDB, Genebass | `life-sciences-databases:eva-skill`, `life-sciences-databases:ipd-skill`, `life-sciences-databases:efo-ontology-skill`, `life-sciences-databases:rhea-skill`, `life-sciences-databases:rnacentral-skill`, `life-sciences-databases:epigraphdb-skill`, `life-sciences-databases:genebass-gene-burden-skill` |

## Literature

| Need | Use skill |
|---|---|
| RSS/Atom literature monitoring, daily paper alerts, journal feed tracking, topic-filtered research updates | `literature-skill` |
| PubMed publication search, summaries, fetches, links | `life-sciences-literature:ncbi-entrez-skill` |
| PMC open-access article metadata, license, retraction status, files | `life-sciences-literature:ncbi-pmc-skill` |
| bioRxiv/medRxiv preprint metadata and DOI lookups | `life-sciences-literature:biorxiv-skill` |

## Decision Rules

- Use database and literature skills for facts that may change or need citations.
- Use `literature-skill` when the source is a feed URL or when the user wants recurring/daily topic monitoring.
- Use local workflow skills for computation from files on disk.
- Preserve generated results in clearly named project directories.
- When a workflow needs credentials, ask the user to configure them through the child skill's own secure mechanism.
- Never paste, save, or commit real API tokens.
