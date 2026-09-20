"""
Generate advanced publication-grade figures and metrics:
1. figures/polypharmacology_radar.png - Spider/Radar chart comparing top multi-target repurposing hits
2. figures/ligand_efficiency.png - Ligand Efficiency (LE) and Lipophilic Efficiency (LipE) vs Molecular Weight
3. figures/target_selectivity_zscores.png - Target-normalized Z-score selectivity heatmap
4. results/ligand_efficiency_metrics.csv - Detailed per-pair medicinal chemistry parameters
"""

import os
import glob
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from rdkit import Chem
from rdkit.Chem import Descriptors

# Styling
plt.rcParams.update({
    'font.size': 11,
    'font.family': 'sans-serif',
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'axes.labelsize': 12,
    'axes.titlesize': 13,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10
})

os.makedirs("figures", exist_ok=True)
os.makedirs("results", exist_ok=True)

def compute_ligand_descriptors():
    """Extract 2D/3D physicochemical descriptors for each compound using RDKit."""
    data = []
    for sdf_path in sorted(glob.glob("ligands/raw/*.sdf")):
        cid = os.path.basename(sdf_path).replace(".sdf", "")
        suppl = Chem.SDMolSupplier(sdf_path)
        if not suppl or len(suppl) == 0 or suppl[0] is None:
            continue
        mol = suppl[0]
        mw = Descriptors.MolWt(mol)
        ha = mol.GetNumHeavyAtoms()
        logp = Descriptors.MolLogP(mol)
        hbd = Descriptors.NumHDonors(mol)
        hba = Descriptors.NumHAcceptors(mol)
        rotb = Descriptors.NumRotatableBonds(mol)
        tpsa = Descriptors.TPSA(mol)

        data.append({
            'compound_id': cid,
            'mw': mw,
            'heavy_atoms': ha,
            'clogp': logp,
            'hbd': hbd,
            'hba': hba,
            'rotb': rotb,
            'tpsa': tpsa
        })
    return pd.DataFrame(data)

def generate_polypharmacology_radar(df):
    """
    Radar/Spider chart showing target-binding profiles of top candidates.
    Axes represent the 8 cancer targets.
    Radial distance represents calculated binding affinity (-ΔG, kcal/mol).
    """
    targets = sorted(df['target_id'].unique())
    target_labels = [
        f"{t}\n({df[df['target_id']==t]['target_name'].iloc[0]})"
        for t in targets
    ]
    num_vars = len(targets)
    angles = np.linspace(0, 2 * np.pi, num_vars, endpoint=False).tolist()
    angles += angles[:1] # complete the loop

    # Key candidates to highlight
    key_compounds = [
        ('PRA', 'Praziquantel', '#d9534f', '-'),
        ('MEB', 'Mebendazole', '#0275d8', '-'),
        ('FEN', 'Fenbendazole', '#f0ad4e', '--'),
        ('KET', 'Ketoconazole', '#5cb85c', '-.'),
        ('NIC', 'Niclosamide', '#6f42c1', ':')
    ]

    fig, ax = plt.subplots(figsize=(9, 9), subplot_kw=dict(polar=True))
    ax.set_theta_offset(np.pi / 2)
    ax.set_theta_direction(-1)

    plt.xticks(angles[:-1], target_labels, size=10, fontweight='bold')
    ax.set_rlabel_position(0)
    plt.yticks([4, 6, 8, 10, 12], ["4", "6", "8", "10", "12 kcal/mol"], color="grey", size=9)
    plt.ylim(3, 13)

    # Reference circle for hit threshold (-7.0 kcal/mol)
    hit_threshold = [7.0] * (num_vars + 1)
    ax.plot(angles, hit_threshold, color='red', linestyle='--', linewidth=1.2, label='Hit Threshold (7.0 kcal/mol)', alpha=0.8)

    for cid, name, color, style in key_compounds:
        sub = df[df['compound_id'] == cid]
        if sub.empty:
            continue
        val_dict = dict(zip(sub['target_id'], -sub['vina_affinity_kcal']))
        values = [val_dict.get(t, 0.0) for t in targets]
        values += values[:1] # complete loop

        ax.plot(angles, values, color=color, linewidth=2.0, linestyle=style, label=name)
        ax.fill(angles, values, color=color, alpha=0.08)

    ax.set_title("Multi-Target Polypharmacology Radar Profile\n(-ΔG Binding Affinity, Higher = Stronger)", pad=25, fontweight='bold', fontsize=13)
    plt.legend(loc='upper right', bbox_to_anchor=(1.25, 1.1), framealpha=0.9, fontsize=10)
    plt.tight_layout()
    plt.savefig("figures/polypharmacology_radar.png")
    plt.close()
    print("✓ Generated figures/polypharmacology_radar.png")

def generate_ligand_efficiency_plot(df, desc_df):
    """
    Two-panel figure:
    Panel A: Binding Energy (-ΔG) vs Molecular Weight (MW).
    Panel B: Ligand Efficiency (LE = -ΔG / Heavy Atoms) per compound across all targets.
    """
    merged = pd.merge(df, desc_df, on='compound_id', how='inner')
    merged['neg_affinity'] = -merged['vina_affinity_kcal']
    # Ligand efficiency (LE): kcal/mol per heavy atom
    merged['LE'] = merged['neg_affinity'] / merged['heavy_atoms']
    # Lipophilic efficiency: pKd - cLogP (pKd approx neg_affinity / 1.363 at 298.15 K)
    merged['pKd'] = merged['neg_affinity'] / 1.363
    merged['LipE'] = merged['pKd'] - merged['clogp']

    # Save metrics table
    merged.to_csv("results/ligand_efficiency_metrics.csv", index=False)
    print("✓ Saved results/ligand_efficiency_metrics.csv")

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 7))

    # Panel A: -ΔG vs MW
    scatter_colors = {
        'Anthelmintic': '#d9534f',
        'Antiparasitic': '#f0ad4e',
        'Antifungal': '#5cb85c',
        'Antidiabetic': '#0275d8',
        'Antimalarial': '#6f42c1',
        'Antirheumatic': '#e83e8c',
        'Other': '#6c757d'
    }

    for drug_cls, grp in merged.groupby('drug_class'):
        c = scatter_colors.get(drug_cls, '#6c757d')
        ax1.scatter(grp['mw'], grp['neg_affinity'], label=drug_cls, color=c, alpha=0.75, s=65, edgecolors='black', linewidth=0.5)

    # Highlight top hits in Panel A
    annot_candidates = ['PRA', 'MEB', 'FEN', 'KET', 'IVE', 'NIC']
    for cid in annot_candidates:
        c_sub = merged[merged['compound_id'] == cid]
        if not c_sub.empty:
            best_pair = c_sub.sort_values(by='neg_affinity', ascending=False).iloc[0]
            ax1.annotate(
                f"{best_pair['compound_id']}-{best_pair['target_id']}",
                (best_pair['mw'], best_pair['neg_affinity']),
                xytext=(8, 4), textcoords='offset points',
                fontsize=9, fontweight='bold',
                arrowprops=dict(arrowstyle='->', lw=0.8, color='black')
            )

    ax1.axhline(7.0, color='red', linestyle='--', linewidth=1.2, label='Hit Cutoff (7.0 kcal/mol)')
    ax1.axhline(8.5, color='darkred', linestyle=':', linewidth=1.2, label='Top-Tier Cutoff (8.5 kcal/mol)')
    ax1.axvline(500, color='gray', linestyle='-.', linewidth=1.0, label='Lipinski MW Limit (500 Da)')

    ax1.set_xlabel('Molecular Weight (Da)', fontweight='bold')
    ax1.set_ylabel('Calculated Affinity -ΔG (kcal/mol)', fontweight='bold')
    ax1.set_title('A: Raw Binding Affinity vs Molecular Weight', fontweight='bold', pad=12)
    ax1.legend(loc='lower right', fontsize=9, framealpha=0.9)
    ax1.grid(True, linestyle=':', alpha=0.5)

    # Panel B: Ligand Efficiency Distribution per Compound
    le_order = merged.groupby('compound_name')['LE'].median().sort_values(ascending=False).index

    sns.boxplot(
        data=merged,
        x='compound_name',
        y='LE',
        order=le_order,
        palette='Spectral',
        ax=ax2,
        linewidth=1.2,
        fliersize=3
    )
    sns.stripplot(
        data=merged,
        x='compound_name',
        y='LE',
        order=le_order,
        color='black',
        alpha=0.5,
        jitter=0.2,
        size=4,
        ax=ax2
    )

    ax2.axhline(0.30, color='forestgreen', linestyle='--', linewidth=1.5, label='Medicinal Chemistry Threshold (LE ≥ 0.30)')
    ax2.set_xlabel('Compound Name', fontweight='bold')
    ax2.set_ylabel('Ligand Efficiency (kcal/mol per Heavy Atom)', fontweight='bold')
    ax2.set_title('B: Ligand Efficiency Across All 8 Targets', fontweight='bold', pad=12)
    ax2.tick_params(axis='x', rotation=40)
    ax2.legend(loc='upper right', fontsize=9, framealpha=0.9)
    ax2.grid(True, linestyle=':', alpha=0.5)

    plt.tight_layout()
    plt.savefig("figures/ligand_efficiency.png")
    plt.close()
    print("✓ Generated figures/ligand_efficiency.png")

def generate_selectivity_zscores(df):
    """
    Target-normalized Z-score heatmap:
    Z = (affinity - mean_target) / std_target
    Negative Z-scores indicate higher-than-average affinity for that target.
    """
    pivot = df.pivot_table(index='compound_name', columns='target_name', values='vina_affinity_kcal')
    # Target column Z-scores
    zscore_df = (pivot - pivot.mean(axis=0)) / pivot.std(axis=0)

    fig, ax = plt.subplots(figsize=(11, 8))
    sns.heatmap(
        zscore_df,
        annot=True,
        fmt='.2f',
        cmap='vlag_r', # Blue = preferential binding (Z < 0), Red = weak/disfavored
        center=0.0,
        cbar_kws={'label': 'Selectivity Z-Score (Lower = More Selective Affinity)'},
        linewidths=0.75,
        linecolor='white',
        ax=ax
    )
    ax.set_title('Target-Normalized Selectivity Profile (Z-Score Matrix)\nIdentifies true target preferences independent of binding cavity depth', pad=15, fontweight='bold')
    ax.set_xlabel('Cancer-Relevant Target Protein', fontweight='bold', labelpad=10)
    ax.set_ylabel('Repurposed Candidate Compound', fontweight='bold', labelpad=10)
    plt.xticks(rotation=30, ha='right')
    plt.tight_layout()
    plt.savefig("figures/target_selectivity_zscores.png")
    plt.close()
    print("✓ Generated figures/target_selectivity_zscores.png")

def main():
    if not os.path.exists("results/consensus_scores.csv"):
        print("results/consensus_scores.csv missing!")
        return

    df = pd.read_csv("results/consensus_scores.csv")
    desc_df = compute_ligand_descriptors()

    generate_polypharmacology_radar(df)
    generate_ligand_efficiency_plot(df, desc_df)
    generate_selectivity_zscores(df)
    print("Advanced figures completed successfully!")

if __name__ == "__main__":
    main()
