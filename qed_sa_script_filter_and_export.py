# save as: qed_sa_filter_and_export.py
# Run from Anaconda Prompt:
#   cd C:\rdkit
#   python qed_sa_filter_and_export.py

import os
import glob
import csv
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

# ====== SETTINGS ======
INPUT_DIR = r"C:\BiND_nicotinic_receptor_calculations\sample_0\gen"
QED_MIN = 0.6
SA_MAX = 4.5

# Output names/folders (same directory as input)
FILTERED_DIR = os.path.join(INPUT_DIR, "filtered")
CSV_OUT = os.path.join(INPUT_DIR, "qed_sa_results_filtered.csv")  # has _filtered in name
# ======================


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


def main():
    Path(FILTERED_DIR).mkdir(parents=True, exist_ok=True)

    sdf_files = glob.glob(os.path.join(INPUT_DIR, "gen_*.sdf"))

    def sort_key(p: str):
        base = os.path.basename(p)
        try:
            return int(base.replace("gen_", "").replace(".sdf", ""))
        except Exception:
            return 10**9

    sdf_files.sort(key=sort_key)

    if not sdf_files:
        print(f"ERROR: No gen_*.sdf files found in: {INPUT_DIR}")
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

        passed = (qed_val >= QED_MIN) and (sa_val <= SA_MAX)

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

        # ---- Save filtered SDF (copy molecule into new SDF) ----
        out_sdf = os.path.join(FILTERED_DIR, fname)
        w = Chem.SDWriter(out_sdf)
        # Optionally annotate properties in SDF
        mol.SetProp("QED", f"{qed_val:.6f}")
        mol.SetProp("SA", f"{sa_val:.6f}")
        w.write(mol)
        w.close()

        # ---- Save filtered PDB ----
        # RDKit writes PDB; needs coordinates. Use SDF coords if present, else generate.
        out_pdb = os.path.join(FILTERED_DIR, fname.replace(".sdf", ".pdb"))
        try:
            mol3d = ensure_3d(mol)
            Chem.MolToPDBFile(mol3d, out_pdb)
        except Exception as e:
            # Keep going; record PDB failure in status
            # (We dont want to lose the SDF export if PDB fails.)
            rows[-1]["status"] = rows[-1]["status"] + f"|pdb_failed:{type(e).__name__}"

    # ---- Write filtered CSV ----
    with open(CSV_OUT, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=["file", "smiles", "qed", "sa", "pass", "status"])
        writer.writeheader()
        writer.writerows(rows)

    total = len(rows)
    print(f"Done.")
    print(f"Input folder:  {INPUT_DIR}")
    print(f"Filtered dir:  {FILTERED_DIR}")
    print(f"CSV written:   {CSV_OUT}")
    print(f"Passed filter: {kept} / {total} (QED >= {QED_MIN}, SA <= {SA_MAX})")


if __name__ == "__main__":
    main()