"""
Analyze protein-ligand residue contacts and hydrogen bonding networks for top docked hits.
Calculates:
1. Interatomic contact distances (< 4.0 Å) between ligand and receptor atoms.
2. Putative hydrogen bonds (< 3.5 Å, heteroatom-heteroatom: N/O/S).
3. Hydrophobic contacts (< 4.0 Å, carbon-carbon).
Outputs:
- results/residue_contacts.csv
- figures/residue_contacts_top_hits.png
"""

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from Bio.PDB import PDBParser

plt.rcParams.update({
    'font.size': 10,
    'font.family': 'sans-serif',
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
})

os.makedirs("results", exist_ok=True)
os.makedirs("figures", exist_ok=True)

TOP_PAIRS = [
    ('PRA', 'BTUB', 'Praziquantel', 'β-Tubulin', 'Colchicine pocket'),
    ('MEB', 'BTUB', 'Mebendazole', 'β-Tubulin', 'Colchicine pocket'),
    ('FEN', 'BTUB', 'Fenbendazole', 'β-Tubulin', 'Colchicine pocket'),
    ('IVE', 'GLU1', 'Ivermectin', 'GLUT1', 'Glucose translocation cavity'),
    ('KET', 'TRXR', 'Ketoconazole', 'TXNRD1', 'NADPH/FAD interfacial domain'),
    ('NIC', 'STA3', 'Niclosamide', 'STAT3', 'SH2 dimerization domain'),
    ('MEB', 'MDM2', 'Mebendazole', 'MDM2', 'p53-transactivation pocket'),
    ('IVE', 'IMPA', 'Ivermectin', 'Importin-α', 'Major NLS cargo pocket')
]

def parse_pdbqt_model1(pdbqt_file):
    """Extract 3D coordinates and elements of MODEL 1 from a docked PDBQT file."""
    atoms = []
    if not os.path.exists(pdbqt_file):
        return atoms
    with open(pdbqt_file, 'r') as f:
        in_model1 = False
        for line in f:
            if line.startswith('MODEL 1'):
                in_model1 = True
                continue
            if line.startswith('ENDMDL') and in_model1:
                break
            if in_model1 and (line.startswith('ATOM') or line.startswith('HETATM')):
                atom_name = line[12:16].strip()
                element = atom_name[0]
                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    atoms.append({
                        'name': atom_name,
                        'element': element,
                        'coord': np.array([x, y, z])
                    })
                except ValueError:
                    continue
    return atoms

def analyze_pair(comp_id, target_id, comp_name, target_name, pocket_desc, parser):
    pdbqt_path = f"results/vina/{comp_id}_{target_id}_poses.pdbqt"
    rec_path = f"receptors/prepared/{target_id}_clean.pdb"

    if not os.path.exists(pdbqt_path) or not os.path.exists(rec_path):
        return []

    lig_atoms = parse_pdbqt_model1(pdbqt_path)
    if not lig_atoms:
        return []

    rec_struct = parser.get_structure(target_id, rec_path)
    contacts = []

    for res in rec_struct.get_residues():
        res_name = res.get_resname()
        res_num = res.id[1]
        chain = res.get_parent().id
        for r_atom in res.get_atoms():
            r_coord = r_atom.get_coord()
            r_elem = r_atom.element
            for l_atom in lig_atoms:
                dist = np.linalg.norm(r_coord - l_atom['coord'])
                if dist <= 4.0:
                    is_hbond = (dist <= 3.5 and r_elem in ['N', 'O', 'S'] and l_atom['element'] in ['N', 'O', 'S'])
                    is_hydrophobic = (r_elem == 'C' and l_atom['element'] == 'C')
                    contacts.append({
                        'compound_id': comp_id,
                        'compound_name': comp_name,
                        'target_id': target_id,
                        'target_name': target_name,
                        'pocket': pocket_desc,
                        'chain': chain,
                        'res_num': res_num,
                        'res_name': res_name,
                        'dist': dist,
                        'is_hbond': is_hbond,
                        'is_hydrophobic': is_hydrophobic
                    })
    return contacts

def plot_residue_fingerprints(summary_df):
    """Plot multi-panel bar chart showing top contacted residues and hydrogen bonds for key complexes."""
    fig, axes = plt.subplots(4, 2, figsize=(16, 18))
    axes = axes.flatten()

    for idx, (comp_id, target_id, comp_name, target_name, _) in enumerate(TOP_PAIRS):
        ax = axes[idx]
        sub = summary_df[(summary_df['compound_id'] == comp_id) & (summary_df['target_id'] == target_id)]
        if sub.empty:
            ax.set_visible(False)
            continue

        # Top 10 interacting residues sorted by closest distance
        top_res = sub.sort_values(by='min_distance').head(10).copy()
        top_res['res_label'] = top_res['chain'] + ':' + top_res['res_name'] + top_res['res_num'].astype(str)

        # Bar heights: 4.0 - min_distance (taller bar = closer / stronger contact)
        proximity = 4.0 - top_res['min_distance']
        colors = ['#d9534f' if hb > 0 else '#337ab7' for hb in top_res['hbond_count']]

        y_pos = np.arange(len(top_res))
        bars = ax.barh(y_pos, proximity, color=colors, edgecolor='black', linewidth=0.6, height=0.65)
        ax.set_yticks(y_pos)
        ax.set_yticklabels(top_res['res_label'], fontweight='bold', fontsize=9)
        ax.invert_yaxis()

        # Labels
        ax.set_xlabel('Contact Intensity (4.0 Å − Distance)', fontsize=9)
        ax.set_title(f"{comp_name} → {target_name} ({comp_id}–{target_id})", fontweight='bold', fontsize=11, pad=8)
        ax.set_xlim(0, 1.8)

        # Annotate exact distance
        for i, (p, d, hb) in enumerate(zip(proximity, top_res['min_distance'], top_res['hbond_count'])):
            annot = f"{d:.2f} Å" + (" [H-bond]" if hb > 0 else "")
            ax.text(p + 0.04, i, annot, va='center', fontsize=8, color='darkred' if hb > 0 else 'black')

    # Global legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='#d9534f', edgecolor='black', label='Hydrogen Bond (≤ 3.5 Å)'),
        Patch(facecolor='#337ab7', edgecolor='black', label='Van der Waals / Hydrophobic Contact (≤ 4.0 Å)')
    ]
    fig.legend(handles=legend_elements, loc='upper center', bbox_to_anchor=(0.5, 0.99), ncol=2, frameon=True, fontsize=11)

    plt.tight_layout(rect=[0, 0, 1, 0.97])
    plt.savefig("figures/residue_contacts_top_hits.png")
    plt.close()
    print("✓ Generated figures/residue_contacts_top_hits.png")

def main():
    parser = PDBParser(QUIET=True)
    all_contacts = []

    print("Analyzing protein-ligand residue contacts for top hits...")
    for comp_id, target_id, comp_name, target_name, pocket in TOP_PAIRS:
        c_list = analyze_pair(comp_id, target_id, comp_name, target_name, pocket, parser)
        all_contacts.extend(c_list)
        print(f"  Processed {comp_id}-{target_id} ({len(c_list)} atom-atom interactions)")

    if not all_contacts:
        print("No contacts identified.")
        return

    df_raw = pd.DataFrame(all_contacts)
    # Aggregate per residue
    summary = df_raw.groupby([
        'compound_id', 'compound_name', 'target_id', 'target_name', 'pocket',
        'chain', 'res_num', 'res_name'
    ]).agg(
        min_distance=('dist', 'min'),
        hbond_count=('is_hbond', 'sum'),
        hydrophobic_count=('is_hydrophobic', 'sum'),
        total_atom_contacts=('dist', 'count')
    ).reset_index()

    summary.to_csv("results/residue_contacts.csv", index=False)
    print("✓ Saved results/residue_contacts.csv")

    plot_residue_fingerprints(summary)
    print("Residue contact analysis complete!")

if __name__ == "__main__":
    main()
