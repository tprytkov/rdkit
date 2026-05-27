import argparse
import pandas as pd


def main():
    parser = argparse.ArgumentParser(
        description="Merge docking summary with QED/SA results using molecule name."
    )
    parser.add_argument(
        "--docking_file",
        required=True,
        help="Path to docking_summary.csv"
    )
    parser.add_argument(
        "--qed_file",
        required=True,
        help="Path to qed_sa_results_filtered.csv"
    )
    parser.add_argument(
        "--output_file",
        required=True,
        help="Path to output merged CSV"
    )
    args = parser.parse_args()

    dock = pd.read_csv(args.docking_file)
    qed = pd.read_csv(args.qed_file)

    # Normalize the join key from the QED/SA file: gen_3.sdf -> gen_3
    if "file" not in qed.columns:
        raise KeyError('QED/SA file must contain a "file" column.')
    qed["molecule_name"] = qed["file"].astype(str).str.replace(".sdf", "", regex=False)

    # Detect the molecule name column in the docking file
    if "molecule_name" in dock.columns:
        docking_name_col = "molecule_name"
    elif "molecule" in dock.columns:
        docking_name_col = "molecule"
    else:
        raise KeyError('Docking file must contain either "molecule_name" or "molecule".')

    dock = dock.rename(columns={docking_name_col: "molecule_name"})

    # Detect affinity column
    affinity_candidates = [
        "first_pose_affinity_kcal_per_mol",
        "binding_affinity_kcal_per_mol",
        "affinity",
    ]
    affinity_col = None
    for col in affinity_candidates:
        if col in dock.columns:
            affinity_col = col
            break
    if affinity_col is None:
        raise KeyError(
            "Docking file must contain one of these affinity columns: "
            + ", ".join(affinity_candidates)
        )

    dock = dock[["molecule_name", affinity_col]].copy()
    dock = dock.rename(columns={affinity_col: "first_pose_affinity_kcal_per_mol"})

    # Keep useful columns from the QED/SA file when present
    wanted_qed_cols = ["molecule_name", "smiles", "qed", "sa", "pass", "status"]
    present_qed_cols = [c for c in wanted_qed_cols if c in qed.columns]
    qed = qed[present_qed_cols].copy()

    merged = dock.merge(qed, on="molecule_name", how="left")

    # Best binders first (most negative affinity at top)
    merged = merged.sort_values("first_pose_affinity_kcal_per_mol", ascending=True)

    merged.to_csv(args.output_file, index=False)

    print("Done.")
    print(f"Docking file: {args.docking_file}")
    print(f"QED/SA file:  {args.qed_file}")
    print(f"Output file:  {args.output_file}")
    print(f"Rows in docking file: {len(dock)}")
    print(f"Rows in merged file:  {len(merged)}")
    if "qed" in merged.columns:
        print(f"Matched rows:         {merged['qed'].notna().sum()}")


if __name__ == "__main__":
    main()