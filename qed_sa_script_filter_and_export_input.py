# save as: qed_sa_script_filter_and_export.py
# Run from Anaconda Prompt:
#   cd C:\rdkit
#   python qed_sa_script_filter_and_export.py --input "C:\...\gen" --output "C:\...\out"

import os
import glob
import csv
import argparse
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import QED, AllChem
from rdkit import RDLogger

# Quiet RDKit warnings/errors (optional). Comment out if you want full logs.
RDLogger.DisableLog("rdApp.error")
RDLogger.DisableLog("rdApp.warning")

# --- SA scorer import ---
# Put these two files next to THIS script:
#   sascorer.py
#   fpscores.pkl.gz
import sascorer  # type: ignore


def try_load_first_mol(sdf_path: str):
    """Load first molecule from SDF, robustly. Returns (mol, status)."""
    try:
        suppl = Chem.SDMolSupplier(sdf_path, sanitize=False, removeHs=False)
    except Exception as e:
        return None, f"SDMolSupplier_error:{type(e).__name__}"

    mol = None
    for m in suppl:
        if m is not None:
            mol = m
            break
    if mol is None:
        return None, "no_molecule_parsed"

    # Try full sanitize, then fallback without kekulization (common for generated aromatics)
    try:
        Chem.SanitizeMol(mol)
        return mol, "ok"
    except Exception as e_full:
        try:
            Chem.SanitizeMol(
                mol,
                sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_KEKULIZE
            )
            return mol, f"ok_no_kekulize:{type(e_full).__name__}"
        except Exception as e_light:
            return None, f"sanitize_failed:{type(e_full).__name__}/{type(e_light).__name__}"


def ensure_3d(mol: Chem.Mol) -> Chem.Mol:
    """
    Ensure molecule has 3D coordinates for PDB export.
    If coords exist in SDF, keep them.
    If not, generate 3D with ETKDG + (optional) UFF optimization.
    """
    mol = Chem.AddHs(mol, addCoords=True)
    confs = mol.GetConformers()
    has_conf = len(confs) > 0 and confs[0].Is3D()

    if not has_conf:
        params = AllChem.ETKDGv3()
        params.randomSeed = 0xC0FFEE
        ok = AllChem.EmbedMolecule(mol, params)
        if ok != 0:
            # If embedding fails, keep as-is (PDB may still fail)
            return mol
        try:
            AllChem.UFFOptimizeMolecule(mol, maxIters=200)
        except Exception:
            pass
    return mol


def parse_args():
    p = argparse.ArgumentParser(
        description="Filter gen_*.sdf files by QED and SA, export filtered SDF + PDB + CSV."
    )
    p.add_argument(
        "--input", "-i", required=True,
        help="Input folder containing gen_*.sdf files"
    )
    p.add_argument(
        "--output", "-o", required=True,
        help="Output folder where filtered files and CSV will be written"
    )
    p.add_argument(
        "--qed", type=float, default=0.6,
        help="Minimum QED threshold (default: 0.6)"
    )
    p.add_argument(
        "--sa", type=float, default=4.5,
        help="Maximum SA threshold (default: 4.5)"
    )
    p.add_argument(
        "--pattern", default="gen_*.sdf",
        help='Glob pattern for input SDFs (default: "gen_*.sdf")'
    )
    return p.parse_args()


def main():
    args = parse_args()

    input_dir = os.path.abspath(args.input)
    output_dir = os.path.abspath(args.output)

    qed_min = float(args.qed)
    sa_max = float(args.sa)

    Path(output_dir).mkdir(parents=True, exist_ok=True)

    # Where filtered molecules go (kept separate inside output folder)
    filtered_dir = os.path.join(output_dir, "filtered")
    Path(filtered_dir).mkdir(parents=True, exist_ok=True)

    # CSV in output folder
    csv_out = os.path.join(output_dir, "qed_sa_results_filtered.csv")

    sdf_files = glob.glob(os.path.join(input_dir, args.pattern))

    def sort_key(pth: str):
        base = os.path.basename(pth)
        # try to sort gen_123.sdf numerically if possible
        try:
            if base.lower().startswith("gen_") and base.lower().endswith(".sdf"):
                return int(base[4:-4])
        except Exception:
            pass
        return base

    sdf_files.sort(key=sort_key)

    if not sdf_files:
        print(f"ERROR: No files found matching {args.pattern} in: {input_dir}")
        return

    rows = []
    kept = 0

    for sdf_path in sdf_files:
        fname = os.path.basename(sdf_path)

        mol, status = try_load_first_mol(sdf_path)
        if mol is None:
            rows.append({"file": fname, "smiles": "", "qed": "", "sa": "", "pass": False, "status": status})
            continue

        # SMILES for logging
        try:
            smiles = Chem.MolToSmiles(mol, canonical=True)
        except Exception:
            smiles = ""

        # Compute QED
        try:
            qed_val = float(QED.qed(mol))
        except Exception as e:
            rows.append({"file": fname, "smiles": smiles, "qed": "", "sa": "", "pass": False,
                         "status": f"qed_failed:{type(e).__name__}"})
            continue

        # Compute SA
        try:
            sa_val = float(sascorer.calculateScore(mol))
        except Exception as e:
            rows.append({"file": fname, "smiles": smiles, "qed": qed_val, "sa": "", "pass": False,
                         "status": f"sa_failed:{type(e).__name__}"})
            continue

        passed = (qed_val >= qed_min) and (sa_val <= sa_max)

        rows.append({
            "file": fname,
            "smiles": smiles,
            "qed": qed_val,
            "sa": sa_val,
            "pass": passed,
            "status": status
        })

        if not passed:
            continue

        kept += 1

        # ---- Save filtered SDF ----
        out_sdf = os.path.join(filtered_dir, fname)
        w = Chem.SDWriter(out_sdf)
        mol.SetProp("QED", f"{qed_val:.6f}")
        mol.SetProp("SA", f"{sa_val:.6f}")
        w.write(mol)
        w.close()

        # ---- Save filtered PDB ----
        out_pdb = os.path.join(filtered_dir, fname.replace(".sdf", ".pdb"))
        try:
            mol3d = ensure_3d(mol)
            Chem.MolToPDBFile(mol3d, out_pdb)
        except Exception as e:
            # Keep going; record PDB failure in status
            rows[-1]["status"] = rows[-1]["status"] + f"|pdb_failed:{type(e).__name__}"

    # ---- Write CSV ----
    with open(csv_out, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=["file", "smiles", "qed", "sa", "pass", "status"])
        writer.writeheader()
        writer.writerows(rows)

    total = len(rows)
    print("Done.")
    print(f"Input folder:   {input_dir}")
    print(f"Output folder:  {output_dir}")
    print(f"Filtered dir:   {filtered_dir}")
    print(f"CSV written:    {csv_out}")
    print(f"Passed filter:  {kept} / {total} (QED >= {qed_min}, SA <= {sa_max})")


if __name__ == "__main__":
    main()