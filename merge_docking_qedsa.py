import argparse
from pathlib import Path
import pandas as pd


def normalize_name(value):
    """
    Convert file names like:
    gen_578.sdf -> gen_578
    gen_578 -> gen_578
    """
    if pd.isna(value):
        return None
    return Path(str(value).strip()).stem


def main():
    parser = argparse.ArgumentParser(
        description="Merge docking data with QED/SA/SMILES data"
    )
    parser.add_argument("--docking", required=True, help="Path to docking CSV")
    parser.add_argument("--qedsa", required=True, help="Path to QED/SA CSV")
    parser.add_argument("--output", required=True, help="Path to output merged CSV")
    args = parser.parse_args()

    docking_df = pd.read_csv(args.docking)
    qedsa_df = pd.read_csv(args.qedsa)

    # Create normalized merge keys
    docking_df["compound_id"] = docking_df["molecule_name"].astype(str).str.strip()
    qedsa_df["compound_id"] = qedsa_df["file"].apply(normalize_name)

    # Merge: keep all molecules from QED/SA file
    merged_df = pd.merge(
        qedsa_df,
        docking_df,
        on="compound_id",
        how="left"
    )

    # Add flag showing whether docking result exists
    merged_df["has_docking"] = merged_df["molecule_name"].notna()

    # Reorder useful columns first
    preferred_cols = [
        "compound_id",
        "file",
        "smiles",
        "qed",
        "sa",
        "pass",
        "status_x",
        "molecule_name",
        "first_pose_affinity_kcal_per_mol",
        "status_y",
        "log_file",
        "has_docking",
    ]

    existing_preferred = [c for c in preferred_cols if c in merged_df.columns]
    remaining_cols = [c for c in merged_df.columns if c not in existing_preferred]
    merged_df = merged_df[existing_preferred + remaining_cols]

    # Sort: docked first, then by better affinity, higher QED, lower SA
    sort_cols = []
    ascending = []

    if "has_docking" in merged_df.columns:
        sort_cols.append("has_docking")
        ascending.append(False)
    if "first_pose_affinity_kcal_per_mol" in merged_df.columns:
        sort_cols.append("first_pose_affinity_kcal_per_mol")
        ascending.append(True)   # more negative is better
    if "qed" in merged_df.columns:
        sort_cols.append("qed")
        ascending.append(False)  # higher is better
    if "sa" in merged_df.columns:
        sort_cols.append("sa")
        ascending.append(True)   # lower is easier synthesis

    if sort_cols:
        merged_df = merged_df.sort_values(by=sort_cols, ascending=ascending, na_position="last")

    merged_df.to_csv(args.output, index=False)

    print(f"Merged file saved to: {args.output}")
    print(f"Rows in QED/SA file: {len(qedsa_df)}")
    print(f"Rows in docking file: {len(docking_df)}")
    print(f"Rows in merged file: {len(merged_df)}")
    print(f"Molecules with docking data: {merged_df['has_docking'].sum()}")


if __name__ == "__main__":
    main()