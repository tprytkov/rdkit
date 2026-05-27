# Drug Design Data Analysis Scripts

This repository contains Python scripts for cheminformatics and molecular modeling workflows used in computational drug discovery. The scripts support analysis of AI-generated molecules, docking results, RDKit-based property filtering, scaffold diversity, protein-ligand interaction fingerprints, and preparation of molecular property batches for SwissADME.

## Project Motivation

The goal of this repository is to organize reusable scripts for prioritizing drug-like compounds after structure-based or AI-guided molecular generation. The workflow combines cheminformatics, docking post-processing, molecular property filtering, scaffold analysis, and interaction-pattern analysis.

## Main Capabilities

- RDKit-based calculation of QED and synthetic accessibility scores
- Filtering and exporting prioritized molecules from SDF files
- Merging docking results with molecular property tables
- Extracting and summarizing AutoDock Vina docking poses
- Analyzing protein-ligand interaction fingerprints from generated molecules
- Visualizing interaction-pattern space using dimensionality reduction
- Comparing scaffold diversity across generated compound sets
- Preparing compound batches for SwissADME analysis

## Repository Contents

| Script | Purpose |
|---|---|
| `qed_sa_script.py` | Calculates QED and synthetic accessibility scores for molecules. |
| `qed_sa_script_filter_and_export.py` | Filters molecules by drug-likeness and exports selected compounds. |
| `qed_sa_script_filter_and_export_input.py` | Input-controlled version of QED/SA filtering and export. |
| `qed_sa_script_no_sanit.py` | Handles molecules that may fail RDKit sanitization. |
| `sascorer.py` | Helper script for synthetic accessibility score calculation. |
| `batch_docking_Vina.py` | Automates batch docking workflows using AutoDock Vina. |
| `vina_first_pose_summary.py` | Extracts docking scores and first-pose summaries from Vina outputs. |
| `merge_dock_sdf.py` | Merges docking results with molecular structure files. |
| `merge_docking_qedsa.py` | Combines docking scores with QED/SA filtering results. |
| `analyze_interaction_fingerprints.py` | Analyzes protein-ligand interaction fingerprint data. |
| `bind_interaction_pkl_summary.py` | Summarizes interaction data stored in BInD-generated `.pkl` files. |
| `bind_interaction_space.py` | Analyzes interaction-pattern space for generated molecules. |
| `bind_interaction_space_tsne_random.py` | Performs t-SNE-based visualization of interaction fingerprints. |
| `compare_scaffold_diversity.py` | Compares scaffold diversity across molecular datasets. |
| `make_swissadme_batches.py` | Splits molecules into batches for SwissADME submission. |

## Example Workflow

1. Generate or collect candidate molecules in SDF format.
2. Calculate QED and synthetic accessibility scores using RDKit.
3. Dock compounds into the target binding pocket using AutoDock Vina.
4. Merge docking results with property-filtering results.
5. Analyze scaffold diversity and interaction fingerprints.
6. Prioritize compounds for further computational or experimental evaluation.

## Tools and Libraries

- Python
- RDKit
- AutoDock Vina
- pandas
- NumPy
- scikit-learn
- matplotlib
- SwissADME-compatible output formatting

## Example Use Cases

These scripts may be useful for:

- Computational drug discovery
- Cheminformatics workflows
- Structure-based virtual screening
- AI-generated molecule prioritization
- Molecular docking result analysis
- Scaffold diversity analysis
- Drug-likeness filtering

## Notes

This repository is intended to demonstrate reproducible scripting workflows for computational chemistry and drug discovery data analysis. Input files and project-specific compound structures should be reviewed before public sharing.
