# save as: qed_sa_script.py
import os
import glob
import csv
import sys
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import QED
from rdkit import RDLogger

# Optional: reduce RDKit error spam to keep output readable
RDLogger.DisableLog("rdApp.error")
RDLogger.DisableLog("rdApp.warning")

# --- SA scorer import ---
# Put sascorer.py and fpscores.pkl.gz in the SAME folder as this script.
try:
    import sascorer  # type: ignore
except Exception as e:
    print("ERROR: Could not import sascorer.")
    print("Put sascorer.py and fpscores.pkl.gz next to this script.")
    raise

def try_load_and_sanitize_first_mol(sdf_path: str):
    """
    Load first molecule from SDF with sanitize=False, then try sanitization.
    Returns (mol, status_message). mol=None if cannot be used.
    """
    try:
        suppl = Chem.SDMolSupplier(sdf_path, sanitize=False, removeHs=False)
    except Exception as e:
        return None, f"SDMolSupplier_error: {e}"

    mol = None
    for m in suppl:
        if m is not None:
            mol = m
            break
    if mol is None:
        return None, "no_molecule_parsed"

    # Try full sanitization first
    try:
        Chem.SanitizeMol(mol)
        return mol, "ok"
    except Exception as e_full:
        # Try a "lighter" sanitization that avoids kekulization step
        # (kekulize is often what fails for generated aromatics)
        try:
            Chem.SanitizeMol(
                mol,
                sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_KEKULIZE
            )
            return mol, f"ok_no_kekulize ({type(e_full).__name__})"
        except Exception as e_light:
            return None, f"sanitize_failed: {type(e_full).__name__} / {type(e_light).__name__}"

def main():
    input_dir = r"C:\BiND_nicotinic_receptor_calculations\sample_0\gen"
    out_csv = os.path.join(input_dir, "qed_sa_results.csv")
    outdir = Path(input_dir)

    sdf_files = glob.glob(os.path.join(input_dir, "gen_*.sdf"))

    def sort_key(p: str):
        base = os.path.basename(p)
        try:
            return int(base.replace("gen_", "").replace(".sdf", ""))
        except Exception:
            return 10**9

    sdf_files.sort(key=sort_key)

    if not sdf_files:
        print(f"ERROR: No gen_*.sdf files found in {input_dir}")
        sys.exit(1)

    rows = []
    for sdf_path in sdf_files:
        fname = os.path.basename(sdf_path)

        mol, status = try_load_and_sanitize_first_mol(sdf_path)
        if mol is None:
            rows.append({
                "file": fname,
                "smiles": "",
                "qed": "",
                "sa": "",
                "status": status
            })
            continue

        # SMILES (may still fail sometimes)
        try:
            smiles = Chem.MolToSmiles(mol, canonical=True)
        except Exception:
            smiles = ""

        # QED
        try:
            qed = float(QED.qed(mol))
        except Exception as e:
            rows.append({
                "file": fname,
                "smiles": smiles,
                "qed": "",
                "sa": "",
                "status": f"qed_failed: {type(e).__name__}"
            })
            continue

        # SA
        try:
            sa = float(sascorer.calculateScore(mol))
        except Exception as e:
            rows.append({
                "file": fname,
                "smiles": smiles,
                "qed": qed,
                "sa": "",
                "status": f"sa_failed: {type(e).__name__}"
            })
            continue

        rows.append({
            "file": fname,
            "smiles": smiles,
            "qed": qed,
            "sa": sa,
            "status": status
        })

    with open(out_csv, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["file", "smiles", "qed", "sa", "status"])
        w.writeheader()
        w.writerows(rows)

    n_ok = sum(1 for r in rows if str(r["qed"]) != "" and str(r["sa"]) != "")
    n_fail = len(rows) - n_ok
    print(f"Done. Wrote: {out_csv}")
    print(f"Success: {n_ok}   Failed/partial: {n_fail}")

if __name__ == "__main__":
    main()