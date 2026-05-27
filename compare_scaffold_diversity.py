import os
import glob
import json
import argparse

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem
from rdkit.Chem.Scaffolds import MurckoScaffold


def natural_key(path: str):
    import re
    base = os.path.basename(path)
    return [int(x) if x.isdigit() else x.lower() for x in re.split(r"(\d+)", base)]


def load_molecules_from_sdf_dir(input_dir: str, pattern: str = "*.sdf"):
    sdf_files = sorted(glob.glob(os.path.join(input_dir, pattern)), key=natural_key)
    records = []
    for sdf_path in sdf_files:
        suppl = Chem.SDMolSupplier(sdf_path, sanitize=True, removeHs=False)
        for i, mol in enumerate(suppl):
            if mol is None:
                continue
            name = mol.GetProp("_Name") if mol.HasProp("_Name") else os.path.splitext(os.path.basename(sdf_path))[0]
            smiles = Chem.MolToSmiles(mol, isomericSmiles=True)
            scaffold = MurckoScaffold.MurckoScaffoldSmiles(mol=mol)
            records.append({
                "molecule_name": name,
                "source_sdf": sdf_path,
                "smiles": smiles,
                "scaffold": scaffold,
                "mol": mol,
            })
    return records


def morgan_fp(mol, radius=2, n_bits=2048):
    return AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=n_bits)


def nearest_neighbor_similarity(fps):
    n = len(fps)
    if n <= 1:
        return [np.nan] * n
    sims = []
    for i, fp in enumerate(fps):
        row = DataStructs.BulkTanimotoSimilarity(fp, fps)
        row[i] = -1.0
        sims.append(max(row))
    return sims


def summarize_run(records, run_name):
    if not records:
        return {
            "run_name": run_name,
            "n_molecules": 0,
            "n_unique_smiles": 0,
            "n_unique_scaffolds": 0,
            "top_scaffold_fraction": np.nan,
            "median_nearest_neighbor_tanimoto": np.nan,
            "mean_nearest_neighbor_tanimoto": np.nan,
        }, pd.DataFrame(), pd.DataFrame()

    df = pd.DataFrame([{k: v for k, v in r.items() if k != "mol"} for r in records])

    unique_smiles = df["smiles"].nunique()
    scaffold_counts = df["scaffold"].fillna("").value_counts(dropna=False).reset_index()
    scaffold_counts.columns = ["scaffold", "count"]
    scaffold_counts.insert(0, "run_name", run_name)

    fps = [morgan_fp(r["mol"]) for r in records]
    nn = nearest_neighbor_similarity(fps)
    df["nearest_neighbor_tanimoto"] = nn
    df.insert(0, "run_name", run_name)

    top_scaffold_fraction = scaffold_counts["count"].iloc[0] / len(df) if len(df) else np.nan

    summary = {
        "run_name": run_name,
        "n_molecules": int(len(df)),
        "n_unique_smiles": int(unique_smiles),
        "n_unique_scaffolds": int(df["scaffold"].nunique()),
        "top_scaffold_fraction": float(top_scaffold_fraction),
        "median_nearest_neighbor_tanimoto": float(np.nanmedian(df["nearest_neighbor_tanimoto"])),
        "mean_nearest_neighbor_tanimoto": float(np.nanmean(df["nearest_neighbor_tanimoto"])),
    }
    return summary, df, scaffold_counts


def make_bar_plot(summary_df, out_png):
    metrics = ["n_unique_smiles", "n_unique_scaffolds", "median_nearest_neighbor_tanimoto"]
    labels = {
        "n_unique_smiles": "Unique molecules",
        "n_unique_scaffolds": "Unique scaffolds",
        "median_nearest_neighbor_tanimoto": "Median NN Tanimoto",
    }

    fig, axes = plt.subplots(1, 3, figsize=(14, 4))
    for ax, metric in zip(axes, metrics):
        ax.bar(summary_df["run_name"], summary_df[metric])
        ax.set_title(labels[metric])
        ax.tick_params(axis="x", rotation=30)
    plt.tight_layout()
    plt.savefig(out_png, dpi=300)
    plt.close()


def make_top_scaffold_plot(scaffold_df, out_png, top_n=10):
    runs = scaffold_df["run_name"].unique().tolist()
    fig, axes = plt.subplots(len(runs), 1, figsize=(10, 4 * max(1, len(runs))), squeeze=False)

    for ax, run in zip(axes[:, 0], runs):
        sub = scaffold_df[scaffold_df["run_name"] == run].head(top_n)
        labels = [s if s else "[empty scaffold]" for s in sub["scaffold"]]
        ax.bar(range(len(sub)), sub["count"])
        ax.set_xticks(range(len(sub)))
        ax.set_xticklabels(labels, rotation=90, fontsize=8)
        ax.set_title(f"Top {min(top_n, len(sub))} scaffolds: {run}")
        ax.set_ylabel("Count")

    plt.tight_layout()
    plt.savefig(out_png, dpi=300)
    plt.close()


def main():
    parser = argparse.ArgumentParser(
        description="Compare scaffold diversity across one or more folders of generated SDF files."
    )
    parser.add_argument(
        "--runs",
        nargs="+",
        required=True,
        help='Run definitions in the form name=folder, e.g. exact=C:\\path\\gen none=C:\\path\\gen'
    )
    parser.add_argument(
        "--output_dir",
        required=True,
        help="Directory for output CSV/JSON/PNG files"
    )
    parser.add_argument(
        "--sdf_pattern",
        default="*.sdf",
        help='Glob pattern for SDF files (default: "*.sdf")'
    )
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    summaries = []
    all_manifest = []
    all_scaffold_counts = []

    for run_def in args.runs:
        if "=" not in run_def:
            raise ValueError(f'Invalid --runs entry "{run_def}". Use name=folder')
        run_name, run_dir = run_def.split("=", 1)
        if not os.path.isdir(run_dir):
            raise FileNotFoundError(f"Run directory not found: {run_dir}")

        records = load_molecules_from_sdf_dir(run_dir, pattern=args.sdf_pattern)
        summary, manifest_df, scaffold_df = summarize_run(records, run_name)

        summaries.append(summary)
        if len(manifest_df):
            all_manifest.append(manifest_df)
        if len(scaffold_df):
            all_scaffold_counts.append(scaffold_df)

    summary_df = pd.DataFrame(summaries).sort_values("run_name")
    summary_csv = os.path.join(args.output_dir, "scaffold_diversity_summary.csv")
    summary_df.to_csv(summary_csv, index=False)

    if all_manifest:
        manifest_df = pd.concat(all_manifest, ignore_index=True)
        manifest_csv = os.path.join(args.output_dir, "scaffold_diversity_manifest.csv")
        manifest_df.to_csv(manifest_csv, index=False)
    else:
        manifest_df = pd.DataFrame()

    if all_scaffold_counts:
        scaffold_counts_df = pd.concat(all_scaffold_counts, ignore_index=True)
        scaffold_counts_csv = os.path.join(args.output_dir, "scaffold_counts_by_run.csv")
        scaffold_counts_df.to_csv(scaffold_counts_csv, index=False)
    else:
        scaffold_counts_df = pd.DataFrame()

    summary_json = os.path.join(args.output_dir, "scaffold_diversity_summary.json")
    with open(summary_json, "w", encoding="utf-8") as f:
        json.dump(summaries, f, indent=2)

    bar_png = os.path.join(args.output_dir, "scaffold_diversity_comparison.png")
    make_bar_plot(summary_df, bar_png)

    if not scaffold_counts_df.empty:
        top_scaffold_png = os.path.join(args.output_dir, "top_scaffolds_by_run.png")
        make_top_scaffold_plot(scaffold_counts_df, top_scaffold_png, top_n=10)

    print("Done.")
    print(f"Summary CSV: {summary_csv}")
    if all_manifest:
        print(f"Manifest CSV: {manifest_csv}")
    if all_scaffold_counts:
        print(f"Scaffold counts CSV: {scaffold_counts_csv}")
        print(f"Top scaffolds PNG: {top_scaffold_png}")
    print(f"Comparison PNG: {bar_png}")
    print(f"Summary JSON: {summary_json}")


if __name__ == "__main__":
    main()
