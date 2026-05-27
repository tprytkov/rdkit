import os
import re
import glob
import pickle
import argparse

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from sklearn.manifold import TSNE
import umap


def natural_key(path: str):
    base = os.path.basename(path)
    return [int(x) if x.isdigit() else x.lower() for x in re.split(r"(\d+)", base)]


def load_pkl(pkl_path: str):
    with open(pkl_path, "rb") as f:
        return pickle.load(f)


def normalize_pairs(value):
    """
    Convert interaction entry into a list of (rec_idx, lig_idx) pairs.
    Handles:
      - [[r1,r2,...],[l1,l2,...]]
      - [[r1,l1],[r2,l2],...]
      - empty values
    """
    if value is None:
        return []

    arr = np.asarray(value, dtype=object)
    if arr.size == 0:
        return []

    if arr.ndim == 2 and arr.shape[0] == 2:
        return [(int(arr[0, i]), int(arr[1, i])) for i in range(arr.shape[1])]

    if arr.ndim == 2 and arr.shape[1] == 2:
        return [(int(arr[i, 0]), int(arr[i, 1])) for i in range(arr.shape[0])]

    return []


def build_feature_space(all_dicts):
    """
    Build a global list of binary interaction features:
    INTERACTIONTYPE__recIDX__ligIDX
    """
    features = set()
    for _, data in all_dicts:
        for inter_type, value in data.items():
            pairs = normalize_pairs(value)
            for rec_idx, lig_idx in pairs:
                features.add(f"{inter_type}__{rec_idx}__{lig_idx}")
    return sorted(features)


def build_fingerprint(data, feature_to_idx):
    vec = np.zeros(len(feature_to_idx), dtype=np.uint8)
    for inter_type, value in data.items():
        pairs = normalize_pairs(value)
        for rec_idx, lig_idx in pairs:
            feat = f"{inter_type}__{rec_idx}__{lig_idx}"
            idx = feature_to_idx.get(feat)
            if idx is not None:
                vec[idx] = 1
    return vec


def main():
    parser = argparse.ArgumentParser(
        description="Convert BInD interaction PKL files into binary fingerprints and visualize them."
    )
    parser.add_argument("--input_dir", required=True, help="Folder containing gen_*.pkl files")
    parser.add_argument("--output_dir", required=True, help="Directory for CSV and plot outputs")
    parser.add_argument("--glob_pattern", default="gen_*.pkl", help='Pattern for PKL files')
    parser.add_argument("--make_umap", action="store_true", help="Generate UMAP coordinates and plot")
    parser.add_argument("--make_tsne", action="store_true", help="Generate t-SNE coordinates and plot")
    parser.add_argument("--umap_n_neighbors", type=int, default=15)
    parser.add_argument("--umap_min_dist", type=float, default=0.1)
    parser.add_argument("--tsne_perplexity", type=float, default=10.0,
                        help="Requested t-SNE perplexity (default: 10.0)")
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    pkl_files = sorted(
        glob.glob(os.path.join(args.input_dir, args.glob_pattern)),
        key=natural_key
    )

    if not pkl_files:
        raise FileNotFoundError(f'No files matching "{args.glob_pattern}" found in {args.input_dir}')

    loaded = []
    for p in pkl_files:
        data = load_pkl(p)
        if not isinstance(data, dict):
            raise ValueError(f"{p} does not contain a dictionary")
        loaded.append((p, data))

    features = build_feature_space(loaded)
    feature_to_idx = {f: i for i, f in enumerate(features)}

    if not features:
        raise ValueError("No interaction features were found in the PKL files.")

    rows = []
    X = []

    for pkl_path, data in loaded:
        mol_name = os.path.splitext(os.path.basename(pkl_path))[0]
        vec = build_fingerprint(data, feature_to_idx)
        X.append(vec)

        row = {
            "molecule_name": mol_name,
            "pkl_file": pkl_path,
        }
        for feat, val in zip(features, vec):
            row[feat] = int(val)
        rows.append(row)

    X = np.vstack(X)
    fp_df = pd.DataFrame(rows)

    fp_csv = os.path.join(args.output_dir, "interaction_fingerprints.csv")
    fp_df.to_csv(fp_csv, index=False)

    print("Fingerprint CSV written:", fp_csv)
    print("Number of molecules:", X.shape[0])
    print("Number of interaction features:", X.shape[1])

    base_df = pd.DataFrame({
        "molecule_name": [os.path.splitext(os.path.basename(p))[0] for p, _ in loaded]
    })

    if args.make_umap:
        reducer = umap.UMAP(
            n_neighbors=args.umap_n_neighbors,
            min_dist=args.umap_min_dist,
            metric="jaccard",
            random_state=42
        )
        coords = reducer.fit_transform(X)

        umap_df = base_df.copy()
        umap_df["umap1"] = coords[:, 0]
        umap_df["umap2"] = coords[:, 1]

        umap_csv = os.path.join(args.output_dir, "interaction_umap_coordinates.csv")
        umap_df.to_csv(umap_csv, index=False)

        plt.figure(figsize=(8, 6))
        plt.scatter(umap_df["umap1"], umap_df["umap2"], s=30)
        plt.xlabel("UMAP-1")
        plt.ylabel("UMAP-2")
        plt.title("Interaction-Space UMAP")
        plt.tight_layout()
        umap_png = os.path.join(args.output_dir, "interaction_umap.png")
        plt.savefig(umap_png, dpi=300)
        plt.close()

        print("UMAP CSV written:", umap_csv)
        print("UMAP PNG written:", umap_png)

    if args.make_tsne:
        perplexity = min(float(args.tsne_perplexity), max(5.0, X.shape[0] // 20))
        print(f"Running t-SNE with perplexity={perplexity} and init='random' ...")

        tsne = TSNE(
            n_components=2,
            perplexity=perplexity,
            init="random",
            random_state=42
        )
        coords = tsne.fit_transform(X)

        tsne_df = base_df.copy()
        tsne_df["tsne1"] = coords[:, 0]
        tsne_df["tsne2"] = coords[:, 1]

        tsne_csv = os.path.join(args.output_dir, "interaction_tsne_coordinates.csv")
        tsne_df.to_csv(tsne_csv, index=False)

        plt.figure(figsize=(8, 6))
        plt.scatter(tsne_df["tsne1"], tsne_df["tsne2"], s=30)
        plt.xlabel("t-SNE 1")
        plt.ylabel("t-SNE 2")
        plt.title("Interaction-Space t-SNE")
        plt.tight_layout()
        tsne_png = os.path.join(args.output_dir, "interaction_tsne.png")
        plt.savefig(tsne_png, dpi=300)
        plt.close()

        print("t-SNE CSV written:", tsne_csv)
        print("t-SNE PNG written:", tsne_png)

    print("Done.")


if __name__ == "__main__":
    main()
