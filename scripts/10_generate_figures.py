"""
Generate publication-quality figures for the docking study.
Outputs:
1. figures/heatmap_vina.png - Compound x Target binding affinity heatmap
2. figures/top_hits_barplot.png - Top 20 compound-target pairs colored by drug class
3. figures/polypharmacology_profile.png - Multi-target hit breadth per compound
4. figures/drug_class_boxplots.png - Binding energy distributions across drug classes
5. figures/hierarchical_clustermap.png - Unsupervised clustering of compounds and targets
"""
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import os

# Set global publication styling
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

CLASS_COLORS = {
    'Anthelmintic': '#d9534f',
    'Antiparasitic': '#f0ad4e',
    'Antifungal': '#5cb85c',
    'Antidiabetic': '#0275d8',
    'Antimalarial': '#6f42c1',
    'Antirheumatic': '#e83e8c',
    'Other': '#6c757d'
}

def generate_heatmap(consensus_df):
    """Figure 1: Full Compound x Target binding affinity matrix."""
    pivot = consensus_df.pivot_table(
        index='compound_name',
        columns='target_name',
        values='vina_affinity_kcal'
    )

    fig, ax = plt.subplots(figsize=(11, 8))
    # Colormap: darker green = stronger affinity (more negative kcal/mol)
    sns.heatmap(
        pivot,
        annot=True,
        fmt='.1f',
        cmap='YlGnBu_r',
        cbar_kws={'label': 'Calculated Binding Affinity (kcal/mol)'},
        linewidths=0.75,
        linecolor='white',
        ax=ax
    )
    ax.set_title('AutoDock Vina Binding Affinity Screen Matrix\n(More negative = stronger binding)', pad=15, fontweight='bold')
    ax.set_xlabel('Cancer-Relevant Target Protein', fontweight='bold', labelpad=10)
    ax.set_ylabel('Repurposed Candidate Compound', fontweight='bold', labelpad=10)
    plt.xticks(rotation=30, ha='right')
    plt.tight_layout()
    plt.savefig("figures/heatmap_vina.png")
    plt.close()
    print("✓ Generated figures/heatmap_vina.png")

def generate_top_hits_barplot(consensus_df):
    """Figure 2: Top 20 highest-affinity compound-target pairs."""
    top20 = consensus_df.sort_values(by='vina_affinity_kcal').head(20).copy()
    top20['pair_label'] = top20['compound_name'] + ' → ' + top20['target_name']

    fig, ax = plt.subplots(figsize=(10, 8))
    bar_colors = [CLASS_COLORS.get(cls, '#6c757d') for cls in top20['drug_class']]

    y_pos = np.arange(len(top20))
    bars = ax.barh(y_pos, top20['vina_affinity_kcal'], color=bar_colors, edgecolor='black', linewidth=0.6, height=0.7)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(top20['pair_label'])
    ax.invert_yaxis()
    ax.set_xlabel('Binding Energy (kcal/mol)', fontweight='bold')
    ax.set_title('Top 20 Predicted Compound–Target Interactions', pad=15, fontweight='bold')

    # Strong binding cutoff line
    ax.axvline(-7.0, color='red', linestyle='--', linewidth=1.2, label='Strong Binder Threshold (-7.0 kcal/mol)')
    ax.axvline(-8.5, color='darkred', linestyle=':', linewidth=1.2, label='Top-Tier Threshold (-8.5 kcal/mol)')

    # Legend
    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor=color, edgecolor='black', label=cls) for cls, color in CLASS_COLORS.items() if cls in top20['drug_class'].values]
    legend_elements.extend([
        plt.Line2D([0], [0], color='red', linestyle='--', label='Strong (≤ -7.0)'),
        plt.Line2D([0], [0], color='darkred', linestyle=':', label='Top-Tier (≤ -8.5)')
    ])
    ax.legend(handles=legend_elements, loc='lower right', framealpha=0.9, fontsize=9)
    plt.tight_layout()
    plt.savefig("figures/top_hits_barplot.png")
    plt.close()
    print("✓ Generated figures/top_hits_barplot.png")

def generate_polypharmacology_profile(consensus_df):
    """Figure 3: Multi-target hit breadth per compound."""
    hits = consensus_df[consensus_df['vina_affinity_kcal'] <= -7.0]
    hit_counts = hits.groupby(['compound_name', 'drug_class'])['target_id'].count().reset_index()
    hit_counts = hit_counts.rename(columns={'target_id': 'hit_count'}).sort_values(by='hit_count', ascending=True)

    fig, ax = plt.subplots(figsize=(9, 6))
    colors = [CLASS_COLORS.get(cls, '#6c757d') for cls in hit_counts['drug_class']]
    ax.barh(hit_counts['compound_name'], hit_counts['hit_count'], color=colors, edgecolor='black', linewidth=0.6)
    ax.set_xlabel('Number of Oncological Targets Bound (≤ -7.0 kcal/mol)', fontweight='bold')
    ax.set_title('Polypharmacology Profile: Multi-Target Hit Capacity', pad=15, fontweight='bold')
    ax.set_xlim(0, 8)
    ax.xaxis.set_major_locator(plt.MaxNLocator(integer=True))

    for i, count in enumerate(hit_counts['hit_count']):
        ax.text(count + 0.15, i, f"{count}/8", va='center', fontweight='bold', fontsize=10)

    plt.tight_layout()
    plt.savefig("figures/polypharmacology_profile.png")
    plt.close()
    print("✓ Generated figures/polypharmacology_profile.png")

def generate_class_boxplots(consensus_df):
    """Figure 4: Binding affinity distributions across drug classes."""
    fig, ax = plt.subplots(figsize=(10, 6))
    order = consensus_df.groupby('drug_class')['vina_affinity_kcal'].median().sort_values().index

    sns.boxplot(
        data=consensus_df,
        x='drug_class',
        y='vina_affinity_kcal',
        order=order,
        palette=CLASS_COLORS,
        ax=ax,
        fliersize=3,
        linewidth=1.2
    )
    sns.stripplot(
        data=consensus_df,
        x='drug_class',
        y='vina_affinity_kcal',
        order=order,
        color='black',
        alpha=0.4,
        jitter=0.2,
        size=4,
        ax=ax
    )
    ax.axhline(-7.0, color='red', linestyle='--', linewidth=1.2, label='Hit Threshold (-7.0 kcal/mol)')
    ax.set_title('Binding Energy Distribution by Pharmacological Drug Class', pad=15, fontweight='bold')
    ax.set_xlabel('Drug Repurposing Class', fontweight='bold')
    ax.set_ylabel('Calculated Binding Energy (kcal/mol)', fontweight='bold')
    plt.xticks(rotation=20, ha='right')
    ax.legend(loc='lower right')
    plt.tight_layout()
    plt.savefig("figures/drug_class_boxplots.png")
    plt.close()
    print("✓ Generated figures/drug_class_boxplots.png")

def generate_clustermap(consensus_df):
    """Figure 5: Hierarchical clustering of compounds and targets."""
    pivot = consensus_df.pivot_table(
        index='compound_name',
        columns='target_name',
        values='vina_affinity_kcal'
    )
    # Fill any missing pairs with the compound's mean or overall neutral value for robust clustering
    pivot_filled = pivot.apply(lambda row: row.fillna(row.mean()), axis=1).fillna(0.0)

    g = sns.clustermap(
        pivot_filled,
        cmap='YlGnBu_r',
        figsize=(10, 9),
        linewidths=0.5,
        annot=True,
        fmt='.1f',
        cbar_kws={'label': 'Binding Affinity (kcal/mol)'}
    )
    g.fig.suptitle('Hierarchical Clustering of Compounds and Oncological Targets', y=1.02, fontweight='bold', fontsize=13)
    plt.savefig("figures/hierarchical_clustermap.png")
    plt.close()
    print("✓ Generated figures/hierarchical_clustermap.png")

def main():
    if not os.path.exists("results/consensus_scores.csv"):
        print("results/consensus_scores.csv not ready yet.")
        return

    df = pd.read_csv("results/consensus_scores.csv")
    generate_heatmap(df)
    generate_top_hits_barplot(df)
    generate_polypharmacology_profile(df)
    generate_class_boxplots(df)
    generate_clustermap(df)
    print("All 5 publication figures generated successfully in figures/!")

if __name__ == "__main__":
    main()
