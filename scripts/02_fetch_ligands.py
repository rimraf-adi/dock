"""
Fetch 3D SDF structures from PubChem for all compounds.
Validates 3D coordinates using RDKit; if only 2D coordinates exist,
generates 3D conformer with ETKDG and MMFF energy minimization.
Saves to ligands/raw/{compound_id}.sdf
"""
import requests
import pandas as pd
import os
import time
import logging
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] %(message)s',
    handlers=[
        logging.StreamHandler(),
        logging.FileHandler("logs/execution.log", mode='a')
    ]
)

def fetch_pubchem_sdf(cid, outpath):
    """Download 3D or 2D SDF from PubChem and ensure 3D coordinates."""
    # 1. Try 3D conformer
    url_3d = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/SDF?record_type=3d"
    res = requests.get(url_3d, timeout=30)
    has_3d = False
    if res.status_code == 200 and len(res.text.strip()) > 0:
        with open(outpath, 'w') as f:
            f.write(res.text)
        supplier = Chem.SDMolSupplier(outpath, removeHs=False)
        mol = supplier[0] if len(supplier) > 0 else None
        if mol is not None and mol.GetNumConformers() > 0:
            coords = mol.GetConformer().GetPositions()
            if any(c[2] != 0.0 for c in coords):
                has_3d = True
                logging.info(f"CID {cid}: downloaded verified 3D conformer")
                return True

    # 2. Fallback to 2D and generate 3D conformer with RDKit
    logging.warning(f"CID {cid}: 3D conformer not directly available. Fetching 2D...")
    url_2d = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/SDF"
    res2 = requests.get(url_2d, timeout=30)
    if res2.status_code == 200 and len(res2.text.strip()) > 0:
        with open(outpath, 'w') as f:
            f.write(res2.text)
        supplier = Chem.SDMolSupplier(outpath, removeHs=False)
        mol = supplier[0] if len(supplier) > 0 else None
        if mol is not None:
            mol = Chem.AddHs(mol)
            embed_res = AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
            if embed_res == 0:
                AllChem.MMFFOptimizeMolecule(mol, maxIters=500)
                writer = Chem.SDWriter(outpath)
                writer.write(mol)
                writer.close()
                logging.info(f"CID {cid}: generated 3D conformer via RDKit ETKDGv3")
                return True
            else:
                logging.error(f"CID {cid}: failed to embed 3D coordinates")
                return False

    logging.error(f"CID {cid}: failed to fetch SDF from PubChem")
    return False

def main():
    compounds = pd.read_csv("configs/compounds.csv")
    os.makedirs("ligands/raw", exist_ok=True)
    os.makedirs("logs", exist_ok=True)

    summary = []
    for _, row in compounds.iterrows():
        cid = row['pubchem_cid']
        comp_id = row['compound_id']
        name = row['name']
        outpath = f"ligands/raw/{comp_id}.sdf"

        logging.info(f"Processing {name} (CID: {cid})...")
        success = fetch_pubchem_sdf(cid, outpath)
        
        mol_wt = None
        n_atoms = None
        if success and os.path.exists(outpath):
            suppl = Chem.SDMolSupplier(outpath, removeHs=False)
            if len(suppl) > 0 and suppl[0] is not None:
                mol = suppl[0]
                mol_wt = round(Descriptors.MolWt(mol), 2)
                n_atoms = mol.GetNumAtoms()

        summary.append({
            'compound_id': comp_id,
            'name': name,
            'pubchem_cid': cid,
            'success': success,
            'atoms': n_atoms,
            'molecular_weight': mol_wt
        })
        time.sleep(0.3)  # PubChem courtesy delay

    df_summary = pd.DataFrame(summary)
    df_summary.to_csv("logs/ligand_fetch_summary.csv", index=False)
    print("\n--- Ligand Fetch Summary ---")
    print(df_summary.to_string())

if __name__ == "__main__":
    main()
