import os
import re
import csv
import glob
import argparse
from typing import Optional, List, Dict

AFFINITY_LINE_RE = re.compile(r"^\s*1\s+(-?\d+(?:\.\d+)?)\s+")
DEFAULT_LOG_PATTERN = "*_vina.log"


def extract_first_pose_affinity(log_path: str) -> Optional[float]:
    """
    Parse an AutoDock Vina log file and return the affinity (kcal/mol)
    for mode 1 (the first/best pose), or None if not found.
    """
    with open(log_path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            m = AFFINITY_LINE_RE.match(line)
            if m:
                return float(m.group(1))
    return None


def infer_molecule_name_from_log(log_path: str) -> str:
    """
    Convert e.g. gen_3_vina.log -> gen_3
    """
    base = os.path.basename(log_path)
    if base.endswith("_vina.log"):
        return base[:-9]
    if base.endswith(".log"):
        return base[:-4]
    return os.path.splitext(base)[0]


def collect_logs(input_dir: str, recursive: bool = True, pattern: str = DEFAULT_LOG_PATTERN) -> List[str]:
    if recursive:
        search_pattern = os.path.join(input_dir, "**", pattern)
        return sorted(glob.glob(search_pattern, recursive=True))
    search_pattern = os.path.join(input_dir, pattern)
    return sorted(glob.glob(search_pattern))


def main():
    parser = argparse.ArgumentParser(
        description="Create a summary CSV of AutoDock Vina first-pose affinities from docking log files."
    )
    parser.add_argument(
        "--input_dir",
        required=True,
        help="Root folder containing docking run subfolders and *_vina.log files."
    )
    parser.add_argument(
        "--output_csv",
        required=True,
        help="Path to the output summary CSV file."
    )
    parser.add_argument(
        "--recursive",
        action="store_true",
        help="Search subfolders recursively (recommended)."
    )
    parser.add_argument(
        "--log_pattern",
        default=DEFAULT_LOG_PATTERN,
        help='Filename pattern for Vina logs (default: "*_vina.log").'
    )
    parser.add_argument(
        "--sort_by_affinity",
        action="store_true",
        help="Sort output by affinity (best binders first, most negative first)."
    )

    args = parser.parse_args()

    input_dir = args.input_dir
    output_csv = args.output_csv

    if not os.path.isdir(input_dir):
        raise FileNotFoundError(f"Input directory not found: {input_dir}")

    log_files = collect_logs(input_dir, recursive=args.recursive, pattern=args.log_pattern)

    if not log_files:
        raise FileNotFoundError(
            f'No log files matching "{args.log_pattern}" were found in: {input_dir}'
        )

    rows: List[Dict[str, object]] = []

    for log_path in log_files:
        affinity = extract_first_pose_affinity(log_path)
        mol_name = infer_molecule_name_from_log(log_path)

        rows.append({
            "molecule_name": mol_name,
            "first_pose_affinity_kcal_per_mol": affinity,
            "log_file": log_path,
        })

    if args.sort_by_affinity:
        rows.sort(
            key=lambda r: (
                r["first_pose_affinity_kcal_per_mol"] is None,
                r["first_pose_affinity_kcal_per_mol"] if r["first_pose_affinity_kcal_per_mol"] is not None else 9999.0
            )
        )

    os.makedirs(os.path.dirname(output_csv) or ".", exist_ok=True)

    with open(output_csv, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "molecule_name",
                "first_pose_affinity_kcal_per_mol",
                "log_file",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)

    n_found = sum(r["first_pose_affinity_kcal_per_mol"] is not None for r in rows)
    n_missing = len(rows) - n_found

    print("Done.")
    print(f"Input folder: {input_dir}")
    print(f"Log files found: {len(log_files)}")
    print(f"Affinities parsed: {n_found}")
    print(f"Missing affinities: {n_missing}")
    print(f"CSV written: {output_csv}")


if __name__ == "__main__":
    main()