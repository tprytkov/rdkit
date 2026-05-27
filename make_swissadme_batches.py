import argparse
import math
from pathlib import Path
import pandas as pd


def normalize_name(value):
    """Convert names like gen_578.sdf -> gen_578"""
    if pd.isna(value):
        return None
    return Path(str(value).strip()).stem.replace(" ", "_")


def choose_name_column(df):
    for col in ["compound_id", "file", "molecule_name"]:
        if col in df.columns:
            return col
    return None


def main():
    parser = argparse.ArgumentParser(
        description="Create SwissADME batch input files from PASS=True molecules only"
    )
    parser.add_argument("--input", required=True, help="Input CSV file")
    parser.add_argument("--output_dir", required=True, help="Output folder")
    parser.add_argument("--smiles_col", default="smiles", help="SMILES column name")
    parser.add_argument("--pass_col", default="pass", help="PASS column name")
    parser.add_argument("--name_col", default=None, help="Optional name column")
    parser.add_argument("--batch_size", type=int, default=15, help="Molecules per batch")
    parser.add_argument("--prefix", default="swissadme_batch_pass", help="Output file prefix")
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(args.input)

    if args.smiles_col not in df.columns:
        raise ValueError(f"Missing SMILES column: {args.smiles_col}")
    if args.pass_col not in df.columns:
        raise ValueError(f"Missing PASS column: {args.pass_col}")

    name_col = args.name_col or choose_name_column(df)
    if not name_col or name_col not in df.columns:
        raise ValueError("Could not find a valid name column")

    # Keep only pass == True
    if df[args.pass_col].dtype == bool:
        df = df[df[args.pass_col]]
    else:
        df = df[df[args.pass_col].astype(str).str.strip().str.lower() == "true"]

    df = df[[args.smiles_col, name_col]].copy()
    df.columns = ["smiles", "raw_name"]

    df = df.dropna(subset=["smiles"])
    df["smiles"] = df["smiles"].astype(str).str.strip()
    df = df[df["smiles"] != ""]

    df["compound_name"] = df["raw_name"].apply(normalize_name)
    df["compound_name"] = df["compound_name"].fillna("")

    missing = df["compound_name"] == ""
    if missing.any():
        df.loc[missing, "compound_name"] = [f"compound_{i+1}" for i in range(missing.sum())]

    df = df.drop_duplicates(subset=["smiles", "compound_name"]).reset_index(drop=True)

    total = len(df)
    if total == 0:
        raise ValueError("No PASS=True molecules with valid SMILES were found.")

    n_batches = math.ceil(total / args.batch_size)

    for i in range(n_batches):
        start = i * args.batch_size
        end = min((i + 1) * args.batch_size, total)
        batch_df = df.iloc[start:end]

        out_file = output_dir / f"{args.prefix}_{i+1:03d}.txt"

        with open(out_file, "w", encoding="utf-8", newline="\n") as f:
            for _, row in batch_df.iterrows():
                f.write(f"{row['smiles']} {row['compound_name']}\n")

        print(f"Wrote {out_file} with {len(batch_df)} molecules")

    print("\nDone.")
    print(f"Total PASS=True molecules written: {total}")
    print(f"Batch size: {args.batch_size}")
    print(f"Number of batch files: {n_batches}")
    print(f"Output folder: {output_dir}")


if __name__ == "__main__":
    main()