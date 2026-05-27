# -*- coding: utf-8 -*-
"""
Batch dock all gen_*.pdb ligands with AutoDock Vina (Windows) using Open Babel conversion.

What it does:
- Takes INPUT_DIR as a command-line parameter
- Finds ligand files like gen_1.pdb, gen_2.pdb, ...
- Converts each ligand PDB -> PDBQT using Open Babel
- Reads docking box/settings from a template config file (conf_nonivamide_1.txt)
- Runs vina_1.2.7_win.exe for each ligand
- Saves each molecule's outputs in a separate folder
- Writes docking_summary.csv

USAGE (Anaconda Prompt):
    conda activate docking_vina
    cd C:\\rdkit
    python batch_docking_Vina.py "C:\\BiND_nicotinic_receptor_calculations\\sample_7ekt\\gen\\filered\\filtered"

Assumptions:
- In INPUT_DIR you have:
    vina_1.2.7_win.exe
    7ekt_no_lig.pdbqt        (required by Vina)
    conf_nonivamide_1.txt
    gen_*.pdb ligand files
- Open Babel is installed here:
    C:\\Users\\tpryt\\miniconda3\\envs\\docking_vina\\Library\\bin\\obabel.exe
"""

import sys
import re
import glob
import shutil
import subprocess
from pathlib import Path


# ===== Fixed Open Babel path (your confirmed path) =====
OBABEL_EXE = Path(r"C:\Users\tpryt\miniconda3\envs\docking_vina\Library\bin\obabel.exe")

# ===== Expected filenames inside input folder =====
VINA_EXE_NAME = "vina_1.2.7_win.exe"
RECEPTOR_PDBQT_NAME = "7ekt_no_lig.pdbqt"       # required
RECEPTOR_PDB_NAME = "7ekt_no_lig.pdb"           # optional/reference
TEMPLATE_CONF_NAME = "conf_nonivamide_1.txt"
LIGAND_PATTERN = "gen_*.pdb"

# ===== Behavior =====
DOCKING_ROOT_NAME = "docking_runs"
SKIP_IF_DONE = True


def parse_vina_template_conf(conf_path: Path) -> dict:
    """Parse simple Vina config lines of the form key = value."""
    params = {}
    with conf_path.open("r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            if "=" not in s:
                continue
            k, v = s.split("=", 1)
            params[k.strip()] = v.strip()
    return params


def natural_sort_key(path_obj: Path) -> int:
    """Sort gen_1.pdb ... gen_10.pdb numerically."""
    m = re.search(r"gen_(\d+)\.pdb$", path_obj.name, re.IGNORECASE)
    return int(m.group(1)) if m else 10**9


def require_exists(path_obj: Path, label: str):
    if not path_obj.exists():
        raise FileNotFoundError(f"{label} not found: {path_obj}")


def convert_ligand_pdb_to_pdbqt(pdb_file: Path, pdbqt_file: Path):
    """
    Convert ligand PDB -> PDBQT using Open Babel.
    Adds hydrogens and Gasteiger charges.
    """
    cmd = [
        str(OBABEL_EXE),
        str(pdb_file),
        "-O", str(pdbqt_file),
        "-xh",
        "--partialcharge", "gasteiger",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(
            f"Open Babel conversion failed for {pdb_file.name}\n"
            f"CMD: {' '.join(cmd)}\n"
            f"STDOUT:\n{result.stdout}\n"
            f"STDERR:\n{result.stderr}"
        )
    if not pdbqt_file.exists():
        raise RuntimeError(f"Expected PDBQT not created: {pdbqt_file}")


def write_vina_conf(out_conf: Path, receptor_pdbqt: Path, ligand_pdbqt: Path, out_pose_pdbqt: Path, params: dict):
    """
    Write per-ligand Vina config using box/settings from template params.
    """
    required = ["center_x", "center_y", "center_z", "size_x", "size_y", "size_z"]
    missing = [k for k in required if k not in params]
    if missing:
        raise ValueError(f"Template config missing keys: {missing}")

    lines = [
        f"receptor = {receptor_pdbqt}",
        f"ligand = {ligand_pdbqt}",
        "",
        f"center_x = {params['center_x']}",
        f"center_y = {params['center_y']}",
        f"center_z = {params['center_z']}",
        "",
        f"size_x = {params['size_x']}",
        f"size_y = {params['size_y']}",
        f"size_z = {params['size_z']}",
        "",
    ]

    # Carry over common optional settings from template if present
    for key in ["energy_range", "exhaustiveness", "num_modes", "cpu", "seed"]:
        if key in params:
            lines.append(f"{key} = {params[key]}")

    lines += [
        "",
        f"out = {out_pose_pdbqt}",
    ]

    out_conf.write_text("\n".join(lines) + "\n", encoding="utf-8")


def run_vina(vina_exe: Path, conf_file: Path, log_file: Path) -> int:
    """
    Run Vina and write stdout/stderr to log_file.
    Returns process return code.
    """
    cmd = [str(vina_exe), "--config", str(conf_file)]
    result = subprocess.run(cmd, capture_output=True, text=True)

    with log_file.open("w", encoding="utf-8", errors="ignore") as f:
        f.write("COMMAND:\n")
        f.write(" ".join(cmd) + "\n\n")
        f.write("STDOUT:\n")
        f.write(result.stdout or "")
        f.write("\n\nSTDERR:\n")
        f.write(result.stderr or "")
        f.write(f"\n\nRETURN_CODE: {result.returncode}\n")

    return result.returncode


def main():
    if len(sys.argv) < 2:
        print("Usage:")
        print(r'  python dock_all_gen_pdb_vina.py "C:\path\to\folder_with_gen_pdb_and_vina_files"')
        sys.exit(1)

    input_dir = Path(sys.argv[1])
    docking_root = input_dir / DOCKING_ROOT_NAME

    # Resolve main files expected in input folder
    vina_exe = input_dir / VINA_EXE_NAME
    receptor_pdbqt = input_dir / RECEPTOR_PDBQT_NAME
    receptor_pdb = input_dir / RECEPTOR_PDB_NAME   # optional
    template_conf = input_dir / TEMPLATE_CONF_NAME

    # Basic checks
    require_exists(input_dir, "Input folder")
    require_exists(OBABEL_EXE, "Open Babel executable")
    require_exists(vina_exe, "Vina executable")
    require_exists(receptor_pdbqt, "Receptor PDBQT")
    require_exists(template_conf, "Template config file")
    if not receptor_pdb.exists():
        print(f"[WARN] Optional receptor PDB not found: {receptor_pdb}")

    print(f"Input folder:    {input_dir}")
    print(f"Using Open Babel:{OBABEL_EXE}")
    print(f"Using Vina:      {vina_exe}")
    print(f"Receptor PDBQT:  {receptor_pdbqt}")
    print(f"Template config: {template_conf}")

    # Read template box/settings
    params = parse_vina_template_conf(template_conf)

    # Find ligands
    ligand_files = [Path(p) for p in glob.glob(str(input_dir / LIGAND_PATTERN))]
    ligand_files.sort(key=natural_sort_key)

    if not ligand_files:
        print(f"[ERROR] No ligands found matching {LIGAND_PATTERN} in {input_dir}")
        sys.exit(2)

    print(f"Found {len(ligand_files)} ligand PDB files.")
    docking_root.mkdir(parents=True, exist_ok=True)

    summary = []  # list of tuples: (molecule, status)

    for lig_pdb in ligand_files:
        stem = lig_pdb.stem  # gen_1
        mol_dir = docking_root / stem
        mol_dir.mkdir(parents=True, exist_ok=True)

        lig_copy = mol_dir / lig_pdb.name
        lig_pdbqt = mol_dir / f"{stem}.pdbqt"
        conf_out = mol_dir / f"conf_{stem}.txt"
        pose_out = mol_dir / f"{stem}_out.pdbqt"
        log_out = mol_dir / f"{stem}_vina.log"

        # Keep original ligand copy in the molecule folder
        if not lig_copy.exists():
            shutil.copy2(lig_pdb, lig_copy)

        if SKIP_IF_DONE and pose_out.exists() and log_out.exists():
            print(f"[SKIP] {stem} (already docked)")
            summary.append((stem, "skipped_existing"))
            continue

        try:
            # 1) Convert ligand
            convert_ligand_pdb_to_pdbqt(lig_pdb, lig_pdbqt)

            # 2) Write per-ligand Vina config
            write_vina_conf(
                out_conf=conf_out,
                receptor_pdbqt=receptor_pdbqt.resolve(),
                ligand_pdbqt=lig_pdbqt.resolve(),
                out_pose_pdbqt=pose_out.resolve(),
                params=params,
            )

            # 3) Run Vina
            rc = run_vina(vina_exe, conf_out, log_out)

            if rc == 0 and pose_out.exists():
                print(f"[OK]   {stem}")
                summary.append((stem, "ok"))
            else:
                print(f"[FAIL] {stem} (return code {rc})")
                summary.append((stem, f"vina_failed_rc_{rc}"))

        except Exception as e:
            print(f"[ERR]  {stem}: {e}")
            # Append Python exception to per-ligand log
            try:
                with log_out.open("a", encoding="utf-8", errors="ignore") as f:
                    f.write("\n\nPYTHON_EXCEPTION:\n")
                    f.write(str(e) + "\n")
            except Exception:
                pass
            summary.append((stem, f"python_error:{type(e).__name__}"))

    # Write batch summary CSV
    summary_csv = docking_root / "docking_summary.csv"
    with summary_csv.open("w", encoding="utf-8", newline="") as f:
        f.write("molecule,status\n")
        for mol, status in summary:
            f.write(f"{mol},{status}\n")

    n_ok = sum(1 for _, s in summary if s == "ok")
    n_total = len(summary)

    print("\nDone.")
    print(f"Docking results root: {docking_root}")
    print(f"Summary CSV:          {summary_csv}")
    print(f"Successful dockings:  {n_ok}/{n_total}")


if __name__ == "__main__":
    main()