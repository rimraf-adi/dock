"""
Define docking grid boxes for each of the 8 targets.
Computes centroids and bounding boxes directly from co-crystallized
ligand/inhibitor coordinates extracted from experimental PDB structures.
Saves to configs/binding_sites.json
"""
import json
import numpy as np
import os
import pandas as pd
import logging

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] %(message)s',
    handlers=[
        logging.StreamHandler(),
        logging.FileHandler("logs/execution.log", mode='a')
    ]
)

def compute_ligand_centroid_and_box(ligand_pdb_path):
    """
    Parses ligand PDB and computes:
    - Centroid (center_x, center_y, center_z)
    - Span in x, y, z
    - Recommended grid box size (max span + padding, min 22.0 A)
    """
    coords = []
    with open(ligand_pdb_path) as f:
        for line in f:
            if line.startswith(('ATOM', 'HETATM')):
                try:
                    x = float(line[30:38].strip())
                    y = float(line[38:46].strip())
                    z = float(line[46:54].strip())
                    coords.append([x, y, z])
                except ValueError:
                    continue

    if not coords:
        return None

    coords = np.array(coords)
    centroid = np.mean(coords, axis=0)
    mins = np.min(coords, axis=0)
    maxs = np.max(coords, axis=0)
    spans = maxs - mins

    # For general docking of diverse compounds (including macrocycles like ivermectin),
    # box size should be at least max(span + 14.0, 24.0) Å
    box_size = np.clip(spans + 14.0, a_min=24.0, a_max=32.0)

    return {
        'center_x': round(float(centroid[0]), 3),
        'center_y': round(float(centroid[1]), 3),
        'center_z': round(float(centroid[2]), 3),
        'size_x': round(float(box_size[0]), 1),
        'size_y': round(float(box_size[1]), 1),
        'size_z': round(float(box_size[2]), 1),
        'span_x': round(float(spans[0]), 1),
        'span_y': round(float(spans[1]), 1),
        'span_z': round(float(spans[2]), 1),
        'num_ligand_atoms': len(coords)
    }

def main():
    targets = pd.read_csv("configs/targets.csv")
    binding_sites = {}

    print("\n--- Binding Site Grid Box Definitions ---")
    for _, row in targets.iterrows():
        tid = row['target_id']
        name = row['name']
        ligand_file = f"receptors/prepared/{tid}_cocrystal_ligand.pdb"

        if os.path.exists(ligand_file):
            box_info = compute_ligand_centroid_and_box(ligand_file)
            if box_info:
                box_info['target_name'] = name
                box_info['pdb_id'] = row['primary_pdb']
                box_info['source'] = f"Co-crystal ligand in {row['primary_pdb']} ({row.get('known_ligand_resname', '')})"
                box_info['site_description'] = row.get('binding_site_description', '')
                binding_sites[tid] = box_info
                
                logging.info(f"✓ {tid} ({name}): Center = ({box_info['center_x']}, {box_info['center_y']}, {box_info['center_z']}) | Box = {box_info['size_x']}x{box_info['size_y']}x{box_info['size_z']} Å")
                continue

        raise RuntimeError(f"Could not define binding site for {tid} ({name})! No ligand file found.")

    with open("configs/binding_sites.json", 'w') as f:
        json.dump(binding_sites, f, indent=2)

    df_sites = pd.DataFrame([
        {
            'target_id': tid,
            'name': v['target_name'],
            'pdb_id': v['pdb_id'],
            'center': f"({v['center_x']}, {v['center_y']}, {v['center_z']})",
            'box_size': f"{v['size_x']} x {v['size_y']} x {v['size_z']}",
            'source': v['source']
        }
        for tid, v in binding_sites.items()
    ])
    print(df_sites.to_string())
    df_sites.to_csv("logs/binding_sites_summary.csv", index=False)
    logging.info("All 8 binding sites computed from experimental co-crystal ligands and saved to configs/binding_sites.json")

if __name__ == "__main__":
    main()
