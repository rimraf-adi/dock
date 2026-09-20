"""
Compute ADMET and Drug-Likeness Profiles for the 13 repurposed candidates using RDKit.
Calculates:
- Molecular Weight (MW)
- LogP (Lipophilicity)
- Hydrogen Bond Donors (HBD)
- Hydrogen Bond Acceptors (HBA)
- Rotatable Bonds (RotB)
- Topological Polar Surface Area (TPSA)
- Lipinski Rule of 5 Violations
- Veber Violations
- Quantitative Estimate of Drug-likeness (QED)
Outputs:
- results/admet_properties.csv
"""

import os
import glob
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, QED

os.makedirs("results", exist_ok=True)

def analyze_admet():
    records = []
    comp_df = pd.read_csv("configs/compounds.csv")
    name_map = dict(zip(comp_df['compound_id'], comp_df['name']))
    class_map = dict(zip(comp_df['compound_id'], comp_df['drug_class']))

    for sdf_path in sorted(glob.glob("ligands/raw/*.sdf")):
        cid = os.path.basename(sdf_path).replace(".sdf", "")
        suppl = Chem.SDMolSupplier(sdf_path)
        if not suppl or len(suppl) == 0 or suppl[0] is None:
            continue
        mol = suppl[0]

        mw = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        hbd = Descriptors.NumHDonors(mol)
        hba = Descriptors.NumHAcceptors(mol)
        rotb = Descriptors.NumRotatableBonds(mol)
        tpsa = Descriptors.TPSA(mol)
        ha = mol.GetNumHeavyAtoms()
        qed_score = QED.qed(mol)

        # Lipinski Rule of 5 violations: MW <= 500, LogP <= 5, HBD <= 5, HBA <= 10
        ro5_violations = 0
        if mw > 500: ro5_violations += 1
        if logp > 5: ro5_violations += 1
        if hbd > 5: ro5_violations += 1
        if hba > 10: ro5_violations += 1

        # Veber rules: RotB <= 10, TPSA <= 140
        veber_violations = 0
        if rotb > 10: veber_violations += 1
        if tpsa > 140: veber_violations += 1

        # Blood-brain barrier (BBB) rule of thumb (Clark 1999 / Kelder 1999): TPSA < 90 Å² and MW < 400
        bbb_likely = (tpsa < 90.0) and (mw < 450.0) and (logp > 1.0) and (logp < 4.0)

        records.append({
            'compound_id': cid,
            'name': name_map.get(cid, cid),
            'drug_class': class_map.get(cid, 'Other'),
            'mw': round(mw, 2),
            'heavy_atoms': ha,
            'clogp': round(logp, 2),
            'hbd': hbd,
            'hba': hba,
            'rotb': rotb,
            'tpsa': round(tpsa, 2),
            'ro5_violations': ro5_violations,
            'veber_violations': veber_violations,
            'qed': round(qed_score, 3),
            'bbb_permeant_predicted': bbb_likely
        })

    df = pd.DataFrame(records)
    df.sort_values(by='qed', ascending=False, inplace=True)
    df.to_csv("results/admet_properties.csv", index=False)
    print("✓ Saved results/admet_properties.csv")
    print(df[['name', 'drug_class', 'mw', 'clogp', 'tpsa', 'qed', 'ro5_violations', 'bbb_permeant_predicted']].to_string(index=False))

if __name__ == "__main__":
    analyze_admet()
