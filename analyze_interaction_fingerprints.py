import os
import argparse
import pandas as pd
import matplotlib.pyplot as plt


def main():
    parser = argparse.ArgumentParser(
        description="Analyze BInD interaction fingerprint CSV: unique patterns, repeats, and feature frequencies."
    )
    parser.add_argument(
        "--input_csv",
        required=True,
        help="Path to interaction_fingerprints.csv"
    )
    parser.add_argument(
        "--output_dir",
        required=True,
        help="Directory for summary CSVs and plots"
    )
    parser.add_argument(
        "--top_n",
        type=int,
        default=20,
        help="Number of top repeated fingerprints / features to save and plot (default: 20)"
    )
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    df = pd.read_csv(args.input_csv)

    id_cols = [c for c in ["molecule_name", "pkl_file"] if c in df.columns]
    feature_cols = [c for c in df.columns if c not in id_cols]

    if not feature_cols:
        raise ValueError("No interaction feature columns found.")

    fp = df[feature_cols].copy()

    n_molecules = len(fp)
    n_features = len(feature_cols)
    n_unique = fp.drop_duplicates().shape[0]

    # Pattern frequency table
    pattern_series = fp.astype(str).agg("|".join, axis=1)
    pattern_counts = pattern_series.value_counts().reset_index()
    pattern_counts.columns = ["pattern_signature", "count"]

    # Attach one representative molecule name per pattern
    rep_map = {}
    for i, sig in enumerate(pattern_series):
        if sig not in rep_map:
            rep_map[sig] = df.iloc[i]["molecule_name"] if "molecule_name" in df.columns else f"row_{i}"
    pattern_counts["representative_molecule"] = pattern_counts["pattern_signature"].map(rep_map)

    pattern_counts_csv = os.path.join(args.output_dir, "interaction_pattern_counts.csv")
    pattern_counts.to_csv(pattern_counts_csv, index=False)

    # Feature frequency table
    feature_counts = fp.sum(axis=0).sort_values(ascending=False).reset_index()
    feature_counts.columns = ["interaction_feature", "count"]

    feature_counts_csv = os.path.join(args.output_dir, "interaction_feature_counts.csv")
    feature_counts.to_csv(feature_counts_csv, index=False)

    # Summary text / CSV
    summary_df = pd.DataFrame([{
        "n_molecules": n_molecules,
        "n_interaction_features": n_features,
        "n_unique_interaction_fingerprints": n_unique,
        "fraction_unique": n_unique / n_molecules if n_molecules else 0.0
    }])
    summary_csv = os.path.join(args.output_dir, "interaction_fingerprint_summary.csv")
    summary_df.to_csv(summary_csv, index=False)

    # Plot top repeated interaction patterns
    top_patterns = pattern_counts.head(args.top_n).copy()
    plt.figure(figsize=(10, 6))
    plt.bar(range(len(top_patterns)), top_patterns["count"])
    plt.xticks(
        range(len(top_patterns)),
        top_patterns["representative_molecule"],
        rotation=90
    )
    plt.ylabel("Count")
    plt.xlabel("Representative molecule of pattern")
    plt.title(f"Top {len(top_patterns)} repeated interaction fingerprints")
    plt.tight_layout()
    pattern_png = os.path.join(args.output_dir, "top_repeated_interaction_patterns.png")
    plt.savefig(pattern_png, dpi=300)
    plt.close()

    # Plot top interaction features
    top_features = feature_counts.head(args.top_n).copy()
    plt.figure(figsize=(10, 6))
    plt.bar(range(len(top_features)), top_features["count"])
    plt.xticks(
        range(len(top_features)),
        top_features["interaction_feature"],
        rotation=90
    )
    plt.ylabel("Count")
    plt.xlabel("Interaction feature")
    plt.title(f"Top {len(top_features)} most frequent interaction features")
    plt.tight_layout()
    feature_png = os.path.join(args.output_dir, "top_interaction_features.png")
    plt.savefig(feature_png, dpi=300)
    plt.close()

    print("Done.")
    print(f"Input CSV:                         {args.input_csv}")
    print(f"Molecules:                         {n_molecules}")
    print(f"Interaction features:              {n_features}")
    print(f"Unique interaction fingerprints:   {n_unique}")
    print(f"Fraction unique:                   {n_unique / n_molecules if n_molecules else 0.0:.3f}")
    print(f"Summary CSV:                       {summary_csv}")
    print(f"Pattern counts CSV:                {pattern_counts_csv}")
    print(f"Feature counts CSV:                {feature_counts_csv}")
    print(f"Pattern plot PNG:                  {pattern_png}")
    print(f"Feature plot PNG:                  {feature_png}")


if __name__ == "__main__":
    main()
