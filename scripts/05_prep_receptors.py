"""
Prepare receptors: PDB -> clean protein -> PDBQT.
Extract co-crystallized ligands to receptors/prepared/{target_id}_cocrystal_ligand.pdb
Saves cleaned protein PDB to receptors/prepared/{target_id}_clean.pdb
Converts cleaned protein to PDBQT using OpenBabel Python API.
"""
import os
import pandas as pd
import logging
from openbabel import openbabel as ob

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] %(message)s',
    handlers=[
        logging.StreamHandler(),
        logging.FileHandler("logs/execution.log", mode='a')
    ]
)

# Specific chains to retain for relevant functional binding unit
RECEPTOR_CHAINS = {
    'BTUB': ['A', 'B'],   # Alpha-Beta heterodimer for colchicine interface
    'IMPA': ['A'],        # Importin-alpha monomer
    'PAK1': ['1'],        # Kinase domain (chain 1 in 2HY8)
    'GLU1': ['A'],        # GLUT1 monomer
    'STA3': ['A'],        # STAT3 core monomer
    'MDM2': ['A'],        # MDM2 monomer
    'WFZD': ['A'],        # FZD8 CRD monomer
    'TRXR': ['A', 'B'],   # Functional homodimer
}

# Special reference ligands (e.g. bound peptide in 4WV6 chain B)
PEPTIDE_LIGAND_CHAINS = {
    'IMPA': 'B',  # TAF8 NLS peptide binding in the armadillo groove
}

def clean_pdb_structure(raw_pdb, clean_pdb, ligand_pdb, keep_chains, target_id):
    """
    Extracts protein ATOM records for selected chains,
    isolates co-crystallized ligand/inhibitors, removes waters & buffers.
    """
    protein_lines = []
    ligand_lines = []

    # Common non-ligand buffer/crystallization molecules to exclude
    EXCLUDE_RESNAMES = {
        'HOH', 'WAT', 'DOD', 'SO4', 'PO4', 'ACT', 'GOL', 'EDO', 'PEG', 'PG4', 'MPD',
        'DMS', 'TRS', 'FMT', 'CL', 'NA', 'MG', 'CA', 'ZN', 'K', 'IOD', 'BMA', 'NAG'
    }

    peptide_chain = PEPTIDE_LIGAND_CHAINS.get(target_id)

    with open(raw_pdb) as f:
        for line in f:
            record = line[:6].strip()
            if record in ('ATOM', 'TER'):
                chain = line[21]
                if peptide_chain and chain == peptide_chain:
                    ligand_lines.append(line)
                elif keep_chains and chain in keep_chains:
                    protein_lines.append(line)
            elif record == 'HETATM':
                resname = line[17:20].strip()
                chain = line[21]
                if resname in EXCLUDE_RESNAMES:
                    continue
                # Keep selenomethionine as part of protein
                if resname in ('MSE', 'CSE'):
                    if keep_chains and chain in keep_chains:
                        protein_lines.append(line)
                else:
                    if keep_chains and chain in keep_chains:
                        ligand_lines.append(line)

    with open(clean_pdb, 'w') as f:
        f.writelines(protein_lines)
        f.write('END\n')

    if ligand_lines:
        with open(ligand_pdb, 'w') as f:
            f.writelines(ligand_lines)
            f.write('END\n')
        logging.info(f"{target_id}: saved {len(ligand_lines)} ligand lines to {ligand_pdb}")
    else:
        logging.warning(f"{target_id}: no co-crystal ligand found in selected chains")

    return len(protein_lines), len(ligand_lines)

def convert_pdb_to_pdbqt(clean_pdb, output_pdbqt):
    """Converts cleaned protein PDB into rigid receptor PDBQT using OpenBabel."""
    conv = ob.OBConversion()
    conv.SetInFormat("pdb")
    conv.SetOutFormat("pdbqt")
    conv.AddOption("r", ob.OBConversion.OUTOPTIONS)  # rigid receptor
    conv.AddOption("h", ob.OBConversion.OUTOPTIONS)  # add polar hydrogens only

    mol = ob.OBMol()
    if not conv.ReadFile(mol, clean_pdb):
        logging.error(f"OpenBabel could not read {clean_pdb}")
        return False

    mol.AddHydrogens(True)  # polar hydrogens
    charge_model = ob.OBChargeModel.FindType("gasteiger")
    if charge_model:
        charge_model.ComputeCharges(mol)

    if not conv.WriteFile(mol, output_pdbqt):
        logging.error(f"OpenBabel could not write {output_pdbqt}")
        return False

    return True

def main():
    targets = pd.read_csv("configs/targets.csv")
    os.makedirs("receptors/prepared", exist_ok=True)
    os.makedirs("logs", exist_ok=True)

    summary = []
    for _, row in targets.iterrows():
        tid = row['target_id']
        pdb_id = row['primary_pdb']
        raw_pdb = f"receptors/raw/{tid}_{pdb_id}.pdb"
        clean_pdb = f"receptors/prepared/{tid}_clean.pdb"
        ligand_pdb = f"receptors/prepared/{tid}_cocrystal_ligand.pdb"
        output_pdbqt = f"receptors/prepared/{tid}.pdbqt"

        if not os.path.exists(raw_pdb):
            logging.error(f"Missing raw PDB: {raw_pdb}")
            continue

        keep_chains = RECEPTOR_CHAINS.get(tid, ['A'])
        n_prot_atoms, n_lig_atoms = clean_pdb_structure(raw_pdb, clean_pdb, ligand_pdb, keep_chains, tid)

        success = convert_pdb_to_pdbqt(clean_pdb, output_pdbqt)
        pdbqt_size = os.path.getsize(output_pdbqt) if (success and os.path.exists(output_pdbqt)) else 0

        logging.info(f"{'✓' if success else '✗'} {tid} ({pdb_id}): {n_prot_atoms} protein lines, {n_lig_atoms} ligand lines -> PDBQT ({round(pdbqt_size/1024, 1)} KB)")
        summary.append({
            'target_id': tid,
            'pdb_id': pdb_id,
            'chains': "".join(keep_chains),
            'protein_records': n_prot_atoms,
            'ligand_records': n_lig_atoms,
            'pdbqt_size_kb': round(pdbqt_size / 1024, 1),
            'status': 'OK' if success else 'FAIL'
        })

    df = pd.DataFrame(summary)
    df.to_csv("logs/receptor_prep_report.csv", index=False)
    print("\n--- Receptor Preparation Summary ---")
    print(df.to_string())

if __name__ == "__main__":
    main()
