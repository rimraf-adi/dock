"""
Prepare ligands: SDF -> Energy-minimized 3D -> PDBQT.
Uses RDKit for chemistry & minimization, and Meeko / OpenBabel for PDBQT conversion.
Handles metallodrugs (e.g. Auranofin with Au) gracefully.
Saves to ligands/prepared/{compound_id}.pdbqt
"""
import os
import pandas as pd
import logging
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
from meeko import MoleculePreparation, PDBQTWriterLegacy
from openbabel import openbabel as ob

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] %(message)s',
    handlers=[
        logging.StreamHandler(),
        logging.FileHandler("logs/execution.log", mode='a')
    ]
)

def prepare_ligand_with_openbabel(sdf_path, pdbqt_path):
    """Fallback ligand preparation using OpenBabel Python API."""
    conv = ob.OBConversion()
    conv.SetInFormat("sdf")
    conv.SetOutFormat("pdbqt")
    conv.AddOption("p", ob.OBConversion.OUTOPTIONS, "7.4")  # protonation at pH 7.4
    conv.AddOption("h", ob.OBConversion.OUTOPTIONS)        # add hydrogens

    mol = ob.OBMol()
    if conv.ReadFile(mol, sdf_path):
        mol.AddHydrogens()
        charge_model = ob.OBChargeModel.FindType("gasteiger")
        if charge_model:
            charge_model.ComputeCharges(mol)
        if conv.WriteFile(mol, pdbqt_path):
            return True
    return False

def prepare_ligand(sdf_path, pdbqt_path, compound_id):
    """Standard pipeline using RDKit and Meeko, falling back to OpenBabel."""
    supplier = Chem.SDMolSupplier(sdf_path, removeHs=False, sanitize=False)
    mol = supplier[0] if len(supplier) > 0 else None
    if mol is None:
        raise ValueError(f"Could not parse {sdf_path}")

    # Sanitize if possible
    try:
        Chem.SanitizeMol(mol)
    except Exception as e:
        logging.warning(f"{compound_id}: Sanitize warning ({e}), continuing...")

    # Check for metals
    has_metal = any(atom.GetSymbol() in ['Au', 'Pt', 'Fe', 'Cu', 'Zn', 'Ru'] for atom in mol.GetAtoms())

    if not has_metal:
        mol = Chem.AddHs(mol, addCoords=True)
        # Attempt MMFF minimization
        try:
            res = AllChem.MMFFOptimizeMolecule(mol, maxIters=500)
            if res != 0:
                AllChem.UFFOptimizeMolecule(mol, maxIters=500)
        except Exception:
            AllChem.UFFOptimizeMolecule(mol, maxIters=500)

        # Attempt Meeko preparation
        try:
            preparator = MoleculePreparation()
            mol_setups = preparator.prepare(mol)
            for setup in mol_setups:
                pdbqt_string, is_ok, error_msg = PDBQTWriterLegacy.write_string(setup)
                if is_ok and len(pdbqt_string.strip()) > 0:
                    with open(pdbqt_path, 'w') as f:
                        f.write(pdbqt_string)
                    return True, "Meeko"
                else:
                    logging.warning(f"{compound_id}: Meeko warning ({error_msg})")
        except Exception as e:
            logging.warning(f"{compound_id}: Meeko exception ({e}), falling back to OpenBabel")

    # Fallback to OpenBabel (ideal for metallodrugs like Auranofin or complex macrocycles)
    success = prepare_ligand_with_openbabel(sdf_path, pdbqt_path)
    if success:
        return True, "OpenBabel"

    return False, "Failed"

def main():
    compounds = pd.read_csv("configs/compounds.csv")
    os.makedirs("ligands/prepared", exist_ok=True)
    os.makedirs("logs", exist_ok=True)

    report = []
    for _, row in compounds.iterrows():
        cid = row['compound_id']
        name = row['name']
        sdf = f"ligands/raw/{cid}.sdf"
        pdbqt = f"ligands/prepared/{cid}.pdbqt"

        if not os.path.exists(sdf):
            logging.error(f"Missing raw SDF: {sdf}")
            continue

        try:
            success, method = prepare_ligand(sdf, pdbqt, cid)
            size_bytes = os.path.getsize(pdbqt) if (success and os.path.exists(pdbqt)) else 0
            
            # Verify PDBQT content
            n_torsions = 0
            if os.path.exists(pdbqt):
                with open(pdbqt) as f:
                    content = f.read()
                    n_torsions = content.count("BRANCH")

            logging.info(f"{'✓' if success else '✗'} {name} ({cid}) -> {method} (size: {size_bytes}B, torsions: {n_torsions})")
            report.append({
                'compound_id': cid,
                'name': name,
                'status': 'OK' if success else 'FAIL',
                'method': method,
                'torsions': n_torsions,
                'file_size_bytes': size_bytes
            })
        except Exception as e:
            logging.error(f"Error preparing {name} ({cid}): {e}")
            report.append({
                'compound_id': cid,
                'name': name,
                'status': f'ERROR: {e}',
                'method': 'None',
                'torsions': 0,
                'file_size_bytes': 0
            })

    df = pd.DataFrame(report)
    df.to_csv("logs/ligand_prep_report.csv", index=False)
    print("\n--- Ligand Preparation Report ---")
    print(df.to_string())

if __name__ == "__main__":
    main()
