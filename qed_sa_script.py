# save as: calc_qed_sa_from_sdf_folder.py
import os
import glob
import csv
import sys
from typing import Optional, Tuple

from rdkit import Chem
from rdkit.Chem import QED

def load_sascorer() -> Optional[object]:
    """
    SA score is usually provided by RDKit's contrib script 'sascorer.py'.
    Common ways to make it importable:
      1) Put sascorer.py in the SAME folder as this script, OR
      2) Put it somewhere on PYTHONPATH, OR
      3) Add its folder to sys.path below.
    """
    # If you know where sascorer.py is, you can hardcode it like:
    # sys.path.append(r"C:\path\to\sascorer_folder")

    try:
        import sascorer  # type: ignore
        return sascorer
    except Exception:
        return None

def read_first_mol_from_sdf(sdf_path: str) -> Tuple[Optional[Chem.Mol], str]:
    """
    Reads the first valid molecule from an SDF file.
    Returns (mol, message). If mol is None, message explains why.
    """
    try:
        suppl = Chem.SDMolSupplier(sdf_path, removeHs=False)
    except Exception as e:
        return None, f"SDMolSupplier error: {e}"

    if suppl is None:
        return None, "SDMolSupplier returned None"

    for mol in suppl:
        if mol is None:
            continue
        return mol, "ok"

    return None, "no valid molecule found in SDF"

def main():
    input_dir = r"C:\BiND_nicotinic_receptor_calculations\sample_0\gen"
    output_csv = r"C:\BiND_nicotinic_receptor_calculations\sample_0\gen\qed_sa_results.csv"

    # Collect gen_1.sdf ... gen_100.sdf (or any gen_*.sdf in folder)
    sdf_files = glob.glob(os.path.join(input_dir, "gen_*.sdf"))

    # Sort naturally by the number after gen_
    def sort_key(p: str):
        base = os.path.basename(p)
        # gen_12.sdf -> 12
        try:
            n = int(base.replace("gen_", "").replace(".sdf", ""))
        except Exception:
            n = 10**9
        return n

    sdf_files.sort(key=sort_key)

    if not sdf_files:
        print(f"ERROR: No files matching gen_*.sdf found in: {input_dir}")
        sys.exit(1)

    sascorer = load_sascorer()
    if sascorer is None:
        print(
            "ERROR: Could not import 'sascorer'.\n"
            "SA score is not part of core RDKit.\n\n"
            "Fix options:\n"
            "  1) Download RDKit contrib sascorer.py and place it next to this script, OR\n"
            "  2) Add the folder containing sascorer.py to PYTHONPATH, OR\n"
            "  3) Edit this script and add: sys.path.append(r\"C:\\\\path\\\\to\\\\sascorer_folder\")\n\n"
            "After that, re-run."
        )
        sys.exit(2)

    rows = []
    for sdf_path in sdf_files:
        fname = os.path.basename(sdf_path)
        mol, status = read_first_mol_from_sdf(sdf_path)
        if mol is None:
            rows.append({
                "file": fname,
                "smiles": "",
                "qed": "",
                "sa": "",
                "status": f"FAILED: {status}",
            })
            continue

        # Create canonical SMILES for logging
        try:
            smiles = Chem.MolToSmiles(mol, canonical=True)
        except Exception:
            smiles = ""

        # QED
        try:
            qed = float(QED.qed(mol))
        except Exception as e:
            qed = ""
            status = f"FAILED_QED: {e}"

        # SA score
        try:
            sa = float(sascorer.calculateScore(mol))
        except Exception as e:
            sa = ""
            status = f"FAILED_SA: {e}"

        rows.append({
            "file": fname,
            "smiles": smiles,
            "qed": qed,
            "sa": sa,
            "status": "ok" if status == "ok" else status,
        })

    # Write CSV
    os.makedirs(os.path.dirname(output_csv), exist_ok=True)
    fieldnames = ["file", "smiles", "qed", "sa", "status"]
    with open(output_csv, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(rows)

    print(f"Done. Wrote {len(rows)} rows to:\n  {output_csv}")

if __name__ == "__main__":
    main()