"""
Fetch PDB structures from RCSB for all verified targets.
Checks resolution, experimental method, and download integrity.
Saves to receptors/raw/{target_id}_{pdb_id}.pdb
"""
import requests
import pandas as pd
import os
import logging
import time

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] %(message)s',
    handlers=[
        logging.StreamHandler(),
        logging.FileHandler("logs/execution.log", mode='a')
    ]
)

def fetch_pdb_file(pdb_id, output_path):
    """Download PDB file directly from RCSB."""
    url = f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb"
    res = requests.get(url, timeout=30)
    if res.status_code == 200 and len(res.text) > 100:
        with open(output_path, 'w') as f:
            f.write(res.text)
        return True
    return False

def get_pdb_rcsb_metadata(pdb_id):
    """Fetch official metadata from RCSB REST API."""
    url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id.upper()}"
    try:
        res = requests.get(url, timeout=15)
        if res.status_code == 200:
            data = res.json()
            entry_info = data.get('rcsb_entry_info', {})
            resolution = entry_info.get('resolution_combined', [None])[0]
            method = entry_info.get('experimental_method', 'Unknown')
            deposition_date = entry_info.get('deposition_date', 'Unknown')
            struct_title = data.get('struct', {}).get('title', '')
            return {
                'resolution': resolution,
                'method': method,
                'date': deposition_date,
                'title': struct_title
            }
    except Exception as e:
        logging.warning(f"Metadata fetch failed for {pdb_id}: {e}")
    return {'resolution': None, 'method': 'Unknown', 'date': 'Unknown', 'title': ''}

def main():
    targets = pd.read_csv("configs/targets.csv")
    os.makedirs("receptors/raw", exist_ok=True)
    os.makedirs("logs", exist_ok=True)

    summary = []
    for _, row in targets.iterrows():
        tid = row['target_id']
        name = row['name']
        pdb_id = row['primary_pdb']
        outpath = f"receptors/raw/{tid}_{pdb_id}.pdb"

        logging.info(f"Fetching {name} ({tid}) PDB {pdb_id}...")
        success = fetch_pdb_file(pdb_id, outpath)
        meta = get_pdb_rcsb_metadata(pdb_id)

        file_size_kb = round(os.path.getsize(outpath) / 1024, 1) if os.path.exists(outpath) else 0

        logging.info(f"{'✓' if success else '✗'} {pdb_id}: {meta['method']} | Res: {meta['resolution']} Å | Size: {file_size_kb} KB")
        logging.info(f"   Title: {meta['title'][:80]}...")

        summary.append({
            'target_id': tid,
            'name': name,
            'pdb_id': pdb_id,
            'success': success,
            'resolution_A': meta['resolution'],
            'method': meta['method'],
            'file_size_kb': file_size_kb,
            'title': meta['title']
        })
        time.sleep(0.3)

    df = pd.DataFrame(summary)
    df.to_csv("logs/receptor_fetch_report.csv", index=False)
    print("\n--- Receptor Fetch Summary ---")
    print(df[['target_id', 'name', 'pdb_id', 'resolution_A', 'method', 'file_size_kb']].to_string())

if __name__ == "__main__":
    main()
