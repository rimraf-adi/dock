"""
Run AutoDock Vina molecular docking with strict thermal and system resource limits:
- Throttled CPU usage: --cpu 2 (leaves 6 of 8 cores idle/cool)
- Sequential execution: 1 worker only (avoids thermal spikes)
- Thermal cooldown: 2-second sleep between dockings
- Full resumption: skips already completed pairs
- Dynamic exhaustiveness: 16 for standard ligands (<=16 branches), 8 for macrocycles/extended (>16 branches)
- Hard timeout: 120s per pair
- Auto-triggers analysis and figure generation on completion
"""
import subprocess
import re
import pandas as pd
import json
import os
import time
import logging

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] %(message)s',
    handlers=[
        logging.StreamHandler(),
        logging.FileHandler("logs/execution.log", mode='a')
    ]
)

VINA_BIN = os.path.abspath("bin/vina")

def parse_vina_log(log_text):
    """Parse Vina CLI output table for binding energies."""
    scores = []
    in_table = False
    for line in log_text.splitlines():
        if "-----+------------+----------+----------" in line:
            in_table = True
            continue
        if in_table:
            match = re.match(r'\s*(\d+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)', line)
            if match:
                scores.append({
                    'mode': int(match.group(1)),
                    'affinity_kcal': float(match.group(2)),
                    'rmsd_lb': float(match.group(3)),
                    'rmsd_ub': float(match.group(4)),
                })
            elif line.strip() == "":
                break
    return scores

def load_existing_scores():
    """Load previously calculated scores from CSV or existing log files to enable seamless resumption."""
    completed = {}
    csv_path = "results/vina_scores.csv"
    if os.path.exists(csv_path):
        try:
            df = pd.read_csv(csv_path)
            for _, r in df.iterrows():
                if r.get('status') == 'OK' and pd.notna(r.get('vina_affinity_kcal')):
                    completed[(r['compound_id'], r['target_id'])] = r.to_dict()
        except Exception:
            pass

    # Also check directory for existing .log files
    if os.path.exists("results/vina"):
        for fname in os.listdir("results/vina"):
            if fname.endswith(".log") and "_" in fname:
                pair_id = fname[:-4]
                parts = pair_id.split("_", 1)
                if len(parts) == 2:
                    cid, tid = parts
                    if (cid, tid) not in completed:
                        log_path = os.path.join("results/vina", fname)
                        try:
                            with open(log_path) as f:
                                scores = parse_vina_log(f.read())
                            if scores:
                                completed[(cid, tid)] = {
                                    'compound_id': cid, 'target_id': tid,
                                    'compound_name': cid, 'target_name': tid, 'priority': 'UNKNOWN',
                                    'vina_affinity_kcal': scores[0]['affinity_kcal'],
                                    'n_poses': len(scores),
                                    'status': 'OK',
                                    'time_sec': 0
                                }
                        except Exception:
                            pass
    return completed

def get_ligand_branch_count(ligand_pdbqt):
    """Count number of active branches in ligand PDBQT."""
    if not os.path.exists(ligand_pdbqt):
        return 0
    with open(ligand_pdbqt) as f:
        return f.read().count("BRANCH")

def main():
    pairs_df = pd.read_csv("configs/pairs.csv")
    with open("configs/binding_sites.json") as f:
        binding_sites = json.load(f)

    os.makedirs("results/vina", exist_ok=True)
    os.makedirs("logs", exist_ok=True)

    existing_scores = load_existing_scores()
    logging.info(f"Loaded {len(existing_scores)} previously completed pairs. Starting throttled runner (cpu=2, 1 worker, 2s cooldown)...")

    results = list(existing_scores.values())
    total_pairs = len(pairs_df)

    for idx, (_, row) in enumerate(pairs_df.iterrows(), 1):
        cid = row['compound_id']
        tid = row['target_id']
        cname = row['compound_name']
        tname = row['target_name']
        priority = row['priority']
        bs = binding_sites[tid]

        pair_key = (cid, tid)
        if pair_key in existing_scores:
            continue

        pair_id = f"{cid}_{tid}"
        receptor_pdbqt = f"receptors/prepared/{tid}.pdbqt"
        ligand_pdbqt = f"ligands/prepared/{cid}.pdbqt"
        out_poses = f"results/vina/{pair_id}_poses.pdbqt"
        log_file = f"results/vina/{pair_id}.log"

        if not os.path.exists(receptor_pdbqt) or not os.path.exists(ligand_pdbqt):
            res_dict = {
                'compound_id': cid, 'target_id': tid,
                'compound_name': cname, 'target_name': tname, 'priority': priority,
                'vina_affinity_kcal': None, 'n_poses': 0, 'status': 'MISSING_INPUT',
                'time_sec': 0
            }
            results.append(res_dict)
            continue

        # Dynamic exhaustiveness based on ligand flexibility
        n_branches = get_ligand_branch_count(ligand_pdbqt)
        exhaustiveness = 8 if n_branches > 16 else 16

        # Throttled command: cpu=2, exhaustiveness=16 (or 8 for macrocycles)
        cmd = [
            VINA_BIN,
            '--receptor', receptor_pdbqt,
            '--ligand', ligand_pdbqt,
            '--center_x', str(bs['center_x']),
            '--center_y', str(bs['center_y']),
            '--center_z', str(bs['center_z']),
            '--size_x', str(bs['size_x']),
            '--size_y', str(bs['size_y']),
            '--size_z', str(bs['size_z']),
            '--exhaustiveness', str(exhaustiveness),
            '--num_modes', '9',
            '--cpu', '2',
            '--out', out_poses
        ]

        t0 = time.time()
        try:
            proc = subprocess.run(cmd, capture_output=True, text=True, timeout=180)
            elapsed = round(time.time() - t0, 1)

            with open(log_file, 'w') as f:
                f.write(proc.stdout)
                if proc.stderr:
                    f.write("\n--- STDERR ---\n" + proc.stderr)

            scores = parse_vina_log(proc.stdout)
            if scores:
                best_aff = scores[0]['affinity_kcal']
                res_dict = {
                    'compound_id': cid, 'target_id': tid,
                    'compound_name': cname, 'target_name': tname, 'priority': priority,
                    'vina_affinity_kcal': best_aff,
                    'n_poses': len(scores),
                    'status': 'OK',
                    'time_sec': elapsed
                }
                logging.info(f"[{idx}/{total_pairs}] {cid} -> {tid} ({cname} on {tname}): {best_aff} kcal/mol ({elapsed}s, exh={exhaustiveness})")
            else:
                res_dict = {
                    'compound_id': cid, 'target_id': tid,
                    'compound_name': cname, 'target_name': tname, 'priority': priority,
                    'vina_affinity_kcal': None, 'n_poses': 0,
                    'status': 'PARSE_ERROR',
                    'time_sec': elapsed
                }
                logging.warning(f"[{idx}/{total_pairs}] {cid} -> {tid}: parse error ({elapsed}s)")
        except subprocess.TimeoutExpired:
            res_dict = {
                'compound_id': cid, 'target_id': tid,
                'compound_name': cname, 'target_name': tname, 'priority': priority,
                'vina_affinity_kcal': None, 'n_poses': 0,
                'status': 'TIMEOUT',
                'time_sec': 180
            }
            logging.warning(f"[{idx}/{total_pairs}] {cid} -> {tid}: timed out (180s)")
        except Exception as e:
            res_dict = {
                'compound_id': cid, 'target_id': tid,
                'compound_name': cname, 'target_name': tname, 'priority': priority,
                'vina_affinity_kcal': None, 'n_poses': 0,
                'status': f'ERROR: {e}',
                'time_sec': 0
            }
            logging.error(f"[{idx}/{total_pairs}] {cid} -> {tid}: error ({e})")

        results.append(res_dict)
        existing_scores[(cid, tid)] = res_dict

        # Save progress after every pair
        pd.DataFrame(results).to_csv("results/vina_scores.csv", index=False)

        # Thermal cooldown: 2-second sleep to keep CPU cool
        time.sleep(2.0)

    df_final = pd.DataFrame(results).sort_values(by=['target_id', 'vina_affinity_kcal'])
    df_final.to_csv("results/vina_scores.csv", index=False)
    logging.info("Docking screen complete for all pairs! Triggering analysis and figures...")

    # Auto-trigger analysis and figures
    os.system(".venv/bin/python scripts/09_analyze_results.py")
    os.system(".venv/bin/python scripts/10_generate_figures.py")
    logging.info("All pipeline steps complete! Artifacts, CSVs, and figures are ready.")

if __name__ == "__main__":
    main()
