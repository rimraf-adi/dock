"""
Benchmark Standard Chemotherapy Drugs and Clinical Controls Against Repurposed Hits.
Fetches canonical benchmark drugs, prepares them, docks them against their respective targets,
and generates direct head-to-head comparison figures and tables.

Benchmark Controls:
1. Colchicine (CID 6167) -> β-Tubulin (BTUB)
2. Nocodazole (CID 4122) -> β-Tubulin (BTUB)
3. Nutlin-3a (CID 11433190) -> MDM2 (MDM2)
4. Stattic (CID 2779853) -> STAT3 (STA3)
5. BAY-876 (CID 73292410) -> GLUT1 (GLU1)
6. Staurosporine (CID 44259) -> PAK1 (PAK1)
"""

import os
import time
import json
import subprocess
import requests
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
from meeko import MoleculePreparation, PDBQTWriterLegacy

os.makedirs("ligands/baselines", exist_ok=True)
os.makedirs("results/baselines", exist_ok=True)
os.makedirs("figures", exist_ok=True)

BENCHMARKS = [
    {
        "id": "COL",
        "name": "Colchicine",
        "cid": 6167,
        "target_id": "BTUB",
        "target_name": "β-Tubulin",
        "role": "Standard Tubulin Inhibitor / Gout & Cancer Chemo",
        "repurposed_comparator": "Praziquantel & Mebendazole"
    },
    {
        "id": "NOC",
        "name": "Nocodazole",
        "cid": 4122,
        "target_id": "BTUB",
        "target_name": "β-Tubulin",
        "role": "Canonical Microtubule Depolymerizer",
        "repurposed_comparator": "Fenbendazole & Albendazole"
    },
    {
        "id": "NUT",
        "name": "Nutlin-3a",
        "cid": 11433190,
        "target_id": "MDM2",
        "target_name": "MDM2 (p53 pathway)",
        "role": "Gold-Standard Clinical MDM2 Antagonist",
        "repurposed_comparator": "Praziquantel & Mebendazole"
    },
    {
        "id": "STT",
        "name": "Stattic",
        "cid": 2779853,
        "target_id": "STA3",
        "target_name": "STAT3 SH2 Domain",
        "role": "Canonical Small-Molecule STAT3 Inhibitor",
        "repurposed_comparator": "Niclosamide"
    },
    {
        "id": "BAY",
        "name": "BAY-876",
        "cid": 73292410,
        "target_id": "GLU1",
        "target_name": "GLUT1 transporter",
        "role": "High-Affinity Nanomolar GLUT1 Inhibitor",
        "repurposed_comparator": "Praziquantel & Ivermectin"
    },
    {
        "id": "STU",
        "name": "Staurosporine",
        "cid": 44259,
        "target_id": "PAK1",
        "target_name": "PAK1 kinase",
        "role": "Broad-Spectrum Kinase Benchmark",
        "repurposed_comparator": "Ketoconazole & Fenbendazole"
    }
]

def fetch_and_prep_baseline(bench):
    cid = bench["cid"]
    bid = bench["id"]
    sdf_path = f"ligands/baselines/{bid}.sdf"
    pdbqt_path = f"ligands/baselines/{bid}.pdbqt"

    if not os.path.exists(sdf_path):
        print(f"Fetching {bench['name']} (CID {cid}) from PubChem...")
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/SDF?record_type=3d"
        r = requests.get(url, timeout=30)
        if r.status_code == 200 and len(r.text.strip()) > 0:
            with open(sdf_path, 'w') as f:
                f.write(r.text)
        else:
            # Fallback 2D
            url_2d = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/SDF"
            r2 = requests.get(url_2d, timeout=30)
            if r2.status_code == 200:
                with open(sdf_path, 'w') as f:
                    f.write(r2.text)
                suppl = Chem.SDMolSupplier(sdf_path, removeHs=False)
                if suppl and len(suppl) > 0 and suppl[0] is not None:
                    mol = Chem.AddHs(suppl[0])
                    AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
                    AllChem.MMFFOptimizeMolecule(mol)
                    writer = Chem.SDWriter(sdf_path)
                    writer.write(mol)
                    writer.close()

    # Prep PDBQT via Meeko
    if not os.path.exists(pdbqt_path) and os.path.exists(sdf_path):
        suppl = Chem.SDMolSupplier(sdf_path, removeHs=False)
        if suppl and len(suppl) > 0 and suppl[0] is not None:
            mol = suppl[0]
            mol = Chem.AddHs(mol, addCoords=True)
            prep = MoleculePreparation(rigid_macrocycles=True)
            mol_setups = prep.prepare(mol)
            for setup in mol_setups:
                pdbqt_str, is_ok, error_msg = PDBQTWriterLegacy.write_string(setup)
                if is_ok and len(pdbqt_str.strip()) > 0:
                    with open(pdbqt_path, 'w') as f:
                        f.write(pdbqt_str)
                    print(f"  Prepared {pdbqt_path}")
                    break

    # Compute RDKit properties
    if os.path.exists(sdf_path):
        suppl = Chem.SDMolSupplier(sdf_path, removeHs=False)
        if suppl and len(suppl) > 0 and suppl[0] is not None:
            m = suppl[0]
            bench['mw'] = round(Descriptors.MolWt(m), 2)
            bench['ha'] = m.GetNumHeavyAtoms()
            bench['clogp'] = round(Descriptors.MolLogP(m), 2)

def dock_baseline(bench, binding_sites):
    bid = bench["id"]
    target = bench["target_id"]
    lig_pdbqt = f"ligands/baselines/{bid}.pdbqt"
    rec_pdbqt = f"receptors/prepared/{target}.pdbqt"
    out_pdbqt = f"results/baselines/{bid}_{target}_poses.pdbqt"
    log_path = f"results/baselines/{bid}_{target}.log"

    if not os.path.exists(lig_pdbqt) or not os.path.exists(rec_pdbqt):
        print(f"Missing input for {bid} -> {target}")
        return None

    site = binding_sites.get(target)
    if not site:
        return None

    cx, cy, cz = site["center_x"], site["center_y"], site["center_z"]
    sx, sy, sz = site["size_x"], site["size_y"], site["size_z"]

    vina_bin = "bin/vina"
    cmd = [
        "nice", "-n", "15",
        vina_bin,
        "--receptor", rec_pdbqt,
        "--ligand", lig_pdbqt,
        "--center_x", str(cx),
        "--center_y", str(cy),
        "--center_z", str(cz),
        "--size_x", str(sx),
        "--size_y", str(sy),
        "--size_z", str(sz),
        "--exhaustiveness", "16",
        "--cpu", "2",
        "--out", out_pdbqt
    ]

    print(f"Docking benchmark {bench['name']} into {target}...")
    with open(log_path, 'w') as lf:
        proc = subprocess.run(cmd, stdout=lf, stderr=subprocess.STDOUT)

    time.sleep(1.0) # Thermal throttle delay

    # Parse best affinity
    best_affinity = None
    if os.path.exists(log_path):
        with open(log_path, 'r') as f:
            for line in f:
                parts = line.split()
                if len(parts) >= 4 and parts[0] == '1':
                    try:
                        best_affinity = float(parts[1])
                        break
                    except ValueError:
                        pass
    return best_affinity

def main():
    with open("configs/binding_sites.json", "r") as f:
        binding_sites = json.load(f)

    results = []
    for b in BENCHMARKS:
        fetch_and_prep_baseline(b)
        affinity = dock_baseline(b, binding_sites)
        if affinity is not None:
            ha = b.get('ha', 20)
            clogp = b.get('clogp', 2.5)
            le = -affinity / ha if ha > 0 else 0.0
            pKd = -affinity / 1.363
            lipE = pKd - clogp

            results.append({
                'baseline_id': b['id'],
                'drug_name': b['name'],
                'role': b['role'],
                'target_id': b['target_id'],
                'target_name': b['target_name'],
                'mw': b.get('mw', 0),
                'heavy_atoms': ha,
                'clogp': clogp,
                'vina_affinity_kcal': affinity,
                'LE': round(le, 3),
                'LipE': round(lipE, 2),
                'repurposed_comparator': b['repurposed_comparator']
            })

    df_base = pd.DataFrame(results)
    df_base.to_csv("results/baselines/chemo_baselines_docking.csv", index=False)
    print("\n✓ Saved results/baselines/chemo_baselines_docking.csv")
    print(df_base[['drug_name', 'target_name', 'vina_affinity_kcal', 'LE', 'LipE']].to_string(index=False))

    # Load consensus repurposed hits to create comparison figure
    df_rep = pd.read_csv("results/ligand_efficiency_metrics.csv")

    comparisons = [
        ("BTUB", "β-Tubulin", [
            ("Colchicine (Chemo Baseline)", -df_base[df_base['baseline_id']=='COL']['vina_affinity_kcal'].iloc[0], df_base[df_base['baseline_id']=='COL']['LE'].iloc[0], '#2c3e50'),
            ("Nocodazole (Chemo Baseline)", -df_base[df_base['baseline_id']=='NOC']['vina_affinity_kcal'].iloc[0], df_base[df_base['baseline_id']=='NOC']['LE'].iloc[0], '#34495e'),
            ("Praziquantel (Repurposed)", -df_rep[(df_rep['compound_id']=='PRA')&(df_rep['target_id']=='BTUB')]['vina_affinity_kcal'].iloc[0], df_rep[(df_rep['compound_id']=='PRA')&(df_rep['target_id']=='BTUB')]['LE'].iloc[0], '#e74c3c'),
            ("Mebendazole (Repurposed)", -df_rep[(df_rep['compound_id']=='MEB')&(df_rep['target_id']=='BTUB')]['vina_affinity_kcal'].iloc[0], df_rep[(df_rep['compound_id']=='MEB')&(df_rep['target_id']=='BTUB')]['LE'].iloc[0], '#3498db'),
            ("Fenbendazole (Repurposed)", -df_rep[(df_rep['compound_id']=='FEN')&(df_rep['target_id']=='BTUB')]['vina_affinity_kcal'].iloc[0], df_rep[(df_rep['compound_id']=='FEN')&(df_rep['target_id']=='BTUB')]['LE'].iloc[0], '#f39c12')
        ]),
        ("MDM2", "MDM2 (p53 pathway)", [
            ("Nutlin-3a (Clinical Baseline)", -df_base[df_base['baseline_id']=='NUT']['vina_affinity_kcal'].iloc[0], df_base[df_base['baseline_id']=='NUT']['LE'].iloc[0], '#2c3e50'),
            ("Praziquantel (Repurposed)", -df_rep[(df_rep['compound_id']=='PRA')&(df_rep['target_id']=='MDM2')]['vina_affinity_kcal'].iloc[0], df_rep[(df_rep['compound_id']=='PRA')&(df_rep['target_id']=='MDM2')]['LE'].iloc[0], '#e74c3c'),
            ("Mebendazole (Repurposed)", -df_rep[(df_rep['compound_id']=='MEB')&(df_rep['target_id']=='MDM2')]['vina_affinity_kcal'].iloc[0], df_rep[(df_rep['compound_id']=='MEB')&(df_rep['target_id']=='MDM2')]['LE'].iloc[0], '#3498db')
        ]),
        ("STA3", "STAT3 SH2 Domain", [
            ("Stattic (Benchmark Baseline)", -df_base[df_base['baseline_id']=='STT']['vina_affinity_kcal'].iloc[0], df_base[df_base['baseline_id']=='STT']['LE'].iloc[0], '#2c3e50'),
            ("Niclosamide (Repurposed)", -df_rep[(df_rep['compound_id']=='NIC')&(df_rep['target_id']=='STA3')]['vina_affinity_kcal'].iloc[0], df_rep[(df_rep['compound_id']=='NIC')&(df_rep['target_id']=='STA3')]['LE'].iloc[0], '#9b59b6')
        ]),
        ("GLU1", "GLUT1 Transporter", [
            ("BAY-876 (Nanomolar Baseline)", -df_base[df_base['baseline_id']=='BAY']['vina_affinity_kcal'].iloc[0], df_base[df_base['baseline_id']=='BAY']['LE'].iloc[0], '#2c3e50'),
            ("Praziquantel (Repurposed)", -df_rep[(df_rep['compound_id']=='PRA')&(df_rep['target_id']=='GLU1')]['vina_affinity_kcal'].iloc[0], df_rep[(df_rep['compound_id']=='PRA')&(df_rep['target_id']=='GLU1')]['LE'].iloc[0], '#e74c3c'),
            ("Mebendazole (Repurposed)", -df_rep[(df_rep['compound_id']=='MEB')&(df_rep['target_id']=='GLU1')]['vina_affinity_kcal'].iloc[0], df_rep[(df_rep['compound_id']=='MEB')&(df_rep['target_id']=='GLU1')]['LE'].iloc[0], '#3498db')
        ])
    ]

    # Generate Comparison Figure: Raw Affinity & Ligand Efficiency
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    axes = axes.flatten()

    for idx, (target_code, target_title, comp_list) in enumerate(comparisons):
        ax = axes[idx]
        labels = [c[0] for c in comp_list]
        affinities = [c[1] for c in comp_list]
        les = [c[2] for c in comp_list]
        colors = [c[3] for c in comp_list]

        y_pos = np.arange(len(labels))
        bars = ax.barh(y_pos, affinities, color=colors, edgecolor='black', linewidth=0.7, height=0.6)
        ax.set_yticks(y_pos)
        ax.set_yticklabels(labels, fontweight='bold', fontsize=10)
        ax.invert_yaxis()

        ax.set_xlabel('Calculated Affinity -ΔG (kcal/mol)', fontweight='bold')
        ax.set_title(f"{target_title} ({target_code}): Repurposed vs Chemotherapy Baselines", pad=12, fontweight='bold', fontsize=11)
        ax.grid(True, linestyle=':', alpha=0.5, axis='x')

        for i, (aff, le) in enumerate(zip(affinities, les)):
            ax.text(aff + 0.15, i, f"-{aff:.2f} kcal/mol  (LE: {le:.3f})", va='center', fontsize=9, fontweight='bold')

    plt.suptitle("Benchmark Validation: Repurposed Candidate Drugs vs Standard Chemotherapy & Clinical Baselines", y=1.01, fontweight='bold', fontsize=14)
    plt.tight_layout()
    plt.savefig("figures/repurposed_vs_chemo_baselines.png")
    plt.close()
    print("✓ Generated figures/repurposed_vs_chemo_baselines.png")

if __name__ == "__main__":
    main()
