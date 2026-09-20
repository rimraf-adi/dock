# Molecular Docking Study — Full Execution Plan

## Objective

Perform a systematic, dual-engine molecular docking screen of **13 repurposed non-oncology drugs** against **8 validated cancer-relevant protein targets** using **AutoDock-Vina** and **Gnina**, then generate a consensus binding affinity matrix, ranked hit list, and publication-ready figures.

**All outputs are saved in the current working directory (CWD):**
```
/Users/adityakinjawadekar/Documents/100xcode/docking/
```

---

## System Context

- **OS:** macOS (Apple Silicon / arm64)
- **Python:** 3.13.0
- **Package manager:** `uv` 0.8.13 (sole tool — no conda/mamba/pip)
- **Homebrew:** available at `/opt/homebrew/bin/brew` (for non-Python binaries only)
- **Already available via system Python:** numpy, pandas, scipy, matplotlib, seaborn, requests
- **Must be installed via `uv`:** RDKit, Open Babel, Meeko, Biopython
- **Must be installed via Homebrew/binary:** AutoDock Vina CLI, Gnina (Docker)

---

## Directory Structure (create before starting)

```
docking/                          ← CWD (project root)
├── PLAN.md                       ← this file
├── pyproject.toml                ← uv project spec (dependencies)
├── scripts/                      ← all Python scripts
│   ├── 01_setup_environment.sh
│   ├── 02_fetch_ligands.py
│   ├── 03_fetch_receptors.py
│   ├── 04_prep_ligands.py
│   ├── 05_prep_receptors.py
│   ├── 06_define_binding_sites.py
│   ├── 07_run_autodock_vina.py
│   ├── 08_run_gnina.py
│   ├── 09_analyze_results.py
│   └── 10_generate_figures.py
├── ligands/
│   ├── raw/                      ← SDF files from PubChem
│   └── prepared/                 ← PDBQT files ready for docking
├── receptors/
│   ├── raw/                      ← PDB files from RCSB
│   └── prepared/                 ← cleaned PDBQT files
├── configs/
│   ├── compounds.csv             ← compound registry
│   ├── targets.csv               ← target registry
│   ├── pairs.csv                 ← compound-target pairing matrix
│   └── binding_sites.json        ← grid box definitions
├── results/
│   ├── vina/                     ← AutoDock Vina output poses + logs
│   ├── gnina/                    ← Gnina output poses + logs
│   ├── vina_scores.csv           ← parsed Vina binding energies
│   ├── gnina_scores.csv          ← parsed Gnina scores (affinity + CNN)
│   ├── consensus_scores.csv      ← merged + ranked results
│   └── summary_report.md         ← human-readable summary
├── figures/
│   ├── heatmap_vina.png
│   ├── heatmap_gnina.png
│   ├── heatmap_consensus.png
│   ├── scatter_vina_vs_gnina.png
│   ├── top_hits_barplot.png
│   └── polypharmacology_radar.png
└── logs/
    └── execution.log             ← timestamped log of all steps
```

---

## Phase 0 — Environment Setup

### Step 0.1: Create `pyproject.toml`

Create `pyproject.toml` in project root:

```toml
[project]
name = "docking-study"
version = "0.1.0"
description = "Multi-compound multi-target molecular docking screen"
requires-python = ">=3.11"
dependencies = [
    "rdkit",
    "openbabel-wheel",
    "meeko",
    "numpy",
    "pandas",
    "scipy",
    "matplotlib",
    "seaborn",
    "requests",
    "biopython",
]
```

### Step 0.2: Initialize uv project and install dependencies

```bash
# Create virtual environment with uv (uses system Python or downloads one)
uv venv .venv --python 3.13

# Activate the virtual environment
source .venv/bin/activate

# Install all dependencies from pyproject.toml
uv pip install -e .
```

> **ALL subsequent `python` commands in this plan MUST be run with the venv active** (`source .venv/bin/activate`) or prefixed with `uv run`.
> The agent can use either pattern — just be consistent. Recommended: `uv run python scripts/XX_script.py`

### Step 0.3: Install AutoDock Vina CLI binary

The `vina` Python binding does not build on arm64 macOS. Use the **standalone CLI binary** instead.

```bash
# Download the Vina 1.2.5 macOS binary
curl -Lo /tmp/vina_1.2.5_mac_arm64 https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v1.2.5/vina_1.2.5_mac_arm64
chmod +x /tmp/vina_1.2.5_mac_arm64
# Place it somewhere on PATH
sudo mv /tmp/vina_1.2.5_mac_arm64 /usr/local/bin/vina
```

> **AGENT: If the above URL is unavailable or the binary doesn't work on this macOS version:**
> 1. Check https://github.com/ccsb-scripps/AutoDock-Vina/releases for the latest arm64 macOS release
> 2. If no arm64 build exists, try the x86_64 build under Rosetta: `arch -x86_64 ./vina_1.2.5_mac_x86_64`
> 3. As last resort, install via Homebrew: `brew install autodock-vina`
> 4. Document which approach worked in `logs/execution.log`

### Step 0.4: Install Gnina

> **IMPORTANT:** Gnina does not have a native macOS ARM build. Two options:

**Option A — Use Gnina via Docker (recommended):**
```bash
brew install --cask docker
# Then pull the gnina image:
docker pull gnina/gnina:latest
```

**Option B — Use Gnina on a Linux machine/server and copy results back.**

**Option C — Skip Gnina, use AutoDock Vina only** (weakens the dual-engine novelty but is still publishable if discussed as a limitation).

> **AGENT DECISION POINT:** If Gnina cannot be installed, proceed with Vina-only and note this. If Docker is available, use the Docker approach. Document the decision in `logs/execution.log`.

### Step 0.5: Verify installations

```bash
# Python packages (run from project root)
uv run python -c "from rdkit import Chem; print(f'RDKit OK: {Chem.rdBase.rdkitVersion}')"
uv run python -c "from openbabel import openbabel; print(f'OpenBabel OK: {openbabel.OBReleaseVersion()}')"
uv run python -c "from meeko import MoleculePreparation; print('Meeko OK')"
uv run python -c "import numpy; import pandas; import scipy; print('Core libs OK')"

# Vina CLI
vina --version

# Gnina (if using Docker)
docker run --rm gnina/gnina gnina --help 2>&1 | head -5
```

**Log all version strings to `logs/execution.log`.**

### Step 0.6: Create directory structure

```bash
mkdir -p ligands/{raw,prepared} receptors/{raw,prepared} configs results/{vina,gnina} figures logs scripts
```

---

## Phase 1 — Data Registry (configs/)

### Step 1.1: Create `configs/compounds.csv`

```csv
compound_id,name,pubchem_cid,drug_class,rationale
FEN,Fenbendazole,3334,Anthelmintic,β-tubulin inhibition; p53 activation
IVE,Ivermectin,6321424,Anthelmintic,PAK1/importin-β1/P-gp inhibition
MEB,Mebendazole,4030,Anthelmintic,Tubulin binding; well-studied in cancer
ALB,Albendazole,2082,Anthelmintic,Tubulin binding; GLUT inhibition
NIC,Niclosamide,4477,Anthelmintic,WNT/STAT3/mTOR pathway inhibition
NIT,Nitazoxanide,41684,Antiparasitic,STAT3/c-Myc inhibition
PRA,Praziquantel,4891,Anthelmintic,Under-studied — gap to fill
ITR,Itraconazole,55283,Antifungal,Hedgehog pathway; angiogenesis inhibition
KET,Ketoconazole,456201,Antifungal,Steroidogenesis inhibition
MET,Metformin,4091,Antidiabetic,AMPK activation; mTOR inhibition
DIS,Disulfiram,3117,Other,Proteasome inhibition (copper-dependent)
CLQ,Chloroquine,2719,Antimalarial,Autophagy inhibition
AUR,Auranofin,6333901,Antirheumatic,Thioredoxin reductase inhibition
```

### Step 1.2: Create `configs/targets.csv`

```csv
target_id,name,gene,relevance,primary_pdb,fallback_pdb,known_ligand_resname
BTUB,β-Tubulin,TUBB,Fenbendazole/mebendazole primary target,4TV8,1TUB,COL
IMPB,Importin-β1,KPNB1,Ivermectin nuclear transport target,2BKU,3ND2,
PAK1,PAK1 kinase,PAK1,Ivermectin kinase target,2HY8,,IPA
GLU1,GLUT1 transporter,SLC2A1,Glucose transport inhibition,4PYP,,BGC
STA3,STAT3,STAT3,Niclosamide/nitazoxanide target,6NJS,1BG1,
MDM2,MDM2 (p53 pathway),MDM2,p53 re-activation,1YCR,,
WFZD,Frizzled (WNT pathway),FZD4,Niclosamide/ivermectin target,6BD4,,
TRXR,Thioredoxin reductase,TXNRD1,Auranofin target,2ZZB,,AUF
```

> **NOTE on PDB IDs:** The agent MUST verify each PDB ID is valid by querying RCSB before downloading. If a PDB ID is unavailable or the structure quality is poor (resolution > 3.0 Å), search RCSB for a better alternative using the gene name. Prefer structures:
> 1. With a co-crystallized ligand (for binding site definition)
> 2. With resolution ≤ 2.5 Å
> 3. Human species origin
> 4. Monomeric or relevant biological assembly

### Step 1.3: Create `configs/pairs.csv`

Generate ALL 13 × 8 = 104 compound-target pairs:

```csv
compound_id,target_id,priority,rationale
FEN,BTUB,HIGH,Known primary target
FEN,GLU1,HIGH,Published GLUT inhibition
FEN,MDM2,MEDIUM,p53 pathway activation
IVE,IMPB,HIGH,Known primary target
IVE,PAK1,HIGH,Known kinase target
IVE,WFZD,MEDIUM,WNT pathway
MEB,BTUB,HIGH,Known primary target
...
```

**AGENT: Generate all 104 pairs. Mark as HIGH if there is direct literature evidence for that compound-target interaction, MEDIUM if the pathway is implicated, LOW for exploratory cross-screening.**

---

## Phase 2 — Fetch Ligand Structures

### Step 2.1: Script `scripts/02_fetch_ligands.py`

For each compound in `configs/compounds.csv`:

```python
"""
Fetch 3D SDF structures from PubChem for all compounds.
Save to ligands/raw/{compound_id}.sdf
"""
import requests
import pandas as pd
import os
import time
import logging

def fetch_pubchem_3d_sdf(cid, output_path):
    """Download 3D conformer SDF from PubChem."""
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/SDF?record_type=3d"
    response = requests.get(url, timeout=30)
    if response.status_code == 200:
        with open(output_path, 'w') as f:
            f.write(response.text)
        return True
    else:
        # Fallback: get 2D and generate 3D later with RDKit
        url_2d = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/SDF"
        response = requests.get(url_2d, timeout=30)
        if response.status_code == 200:
            with open(output_path, 'w') as f:
                f.write(response.text)
            logging.warning(f"CID {cid}: only 2D available, will need 3D generation")
            return True
    return False

# Main execution
compounds = pd.read_csv("configs/compounds.csv")
os.makedirs("ligands/raw", exist_ok=True)

for _, row in compounds.iterrows():
    outpath = f"ligands/raw/{row['compound_id']}.sdf"
    if os.path.exists(outpath):
        logging.info(f"Skipping {row['name']} — already downloaded")
        continue
    success = fetch_pubchem_3d_sdf(row['pubchem_cid'], outpath)
    logging.info(f"{'✓' if success else '✗'} {row['name']} (CID {row['pubchem_cid']})")
    time.sleep(0.5)  # PubChem rate limit courtesy
```

### Step 2.2: Validation

After fetching, verify each SDF file:
- File exists and is > 0 bytes
- Can be parsed by RDKit: `Chem.SDMolSupplier(path)[0] is not None`
- Has 3D coordinates: check z-coordinates are not all zero
- Log compound name, CID, atom count, molecular weight, SMILES to `logs/execution.log`

**If any compound fails, retry once. If still failing, try fetching SMILES from PubChem and generating 3D with RDKit `AllChem.EmbedMolecule()`.**

---

## Phase 3 — Fetch Receptor Structures

### Step 3.1: Script `scripts/03_fetch_receptors.py`

For each target in `configs/targets.csv`:

```python
"""
Fetch PDB structures from RCSB for all targets.
Save to receptors/raw/{target_id}_{pdb_id}.pdb
"""
import requests
import pandas as pd
import os
import logging

def fetch_pdb(pdb_id, output_path):
    """Download PDB file from RCSB."""
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    response = requests.get(url, timeout=30)
    if response.status_code == 200:
        with open(output_path, 'w') as f:
            f.write(response.text)
        return True
    return False

def get_pdb_info(pdb_id):
    """Get resolution, species, method from RCSB API."""
    url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
    response = requests.get(url, timeout=30)
    if response.status_code == 200:
        data = response.json()
        resolution = data.get('rcsb_entry_info', {}).get('resolution_combined', [None])[0]
        method = data.get('rcsb_entry_info', {}).get('experimental_method', 'Unknown')
        return {'resolution': resolution, 'method': method}
    return None

def search_better_pdb(gene_name):
    """Search RCSB for best available structure by gene name."""
    query = {
        "query": {
            "type": "group",
            "logical_operator": "and",
            "nodes": [
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "rcsb_entity_source_organism.rcsb_gene_name.value",
                        "operator": "exact_match",
                        "value": gene_name
                    }
                },
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "rcsb_entity_source_organism.ncbi_scientific_name",
                        "operator": "exact_match",
                        "value": "Homo sapiens"
                    }
                },
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "exptl.method",
                        "operator": "exact_match",
                        "value": "X-RAY DIFFRACTION"
                    }
                }
            ]
        },
        "return_type": "entry",
        "request_options": {
            "sort": [{"sort_by": "rcsb_entry_info.resolution_combined", "direction": "asc"}],
            "paginate": {"start": 0, "rows": 5}
        }
    }
    url = "https://search.rcsb.org/rcsbsearch/v2/query"
    response = requests.post(url, json=query, timeout=30)
    if response.status_code == 200:
        data = response.json()
        return [hit['identifier'] for hit in data.get('result_set', [])]
    return []

# Main execution
targets = pd.read_csv("configs/targets.csv")
os.makedirs("receptors/raw", exist_ok=True)

for _, row in targets.iterrows():
    pdb_id = row['primary_pdb']

    # Validate PDB exists and check resolution
    info = get_pdb_info(pdb_id)
    if info and info['resolution'] and info['resolution'] > 3.0:
        logging.warning(f"{pdb_id} resolution {info['resolution']}Å > 3.0 — searching better")
        alternatives = search_better_pdb(row['gene'])
        if alternatives:
            pdb_id = alternatives[0]
            logging.info(f"Using {pdb_id} instead")

    outpath = f"receptors/raw/{row['target_id']}_{pdb_id}.pdb"
    success = fetch_pdb(pdb_id, outpath)
    logging.info(f"{'✓' if success else '✗'} {row['name']} → {pdb_id}")
```

### Step 3.2: Validation

For each downloaded PDB:
- Log: PDB ID, resolution, method, organism, chain count, residue count
- Verify the file is valid PDB format (starts with HEADER or ATOM lines)
- Warn if resolution > 2.5 Å
- Warn if no co-crystallized ligand found (affects binding site definition)

---

## Phase 4 — Ligand Preparation

### Step 4.1: Script `scripts/04_prep_ligands.py`

Convert each SDF → PDBQT for docking.

```python
"""
Prepare ligands: SDF → energy-minimized → PDBQT
Uses RDKit for chemistry and Meeko for PDBQT conversion.
"""
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
from meeko import MoleculePreparation, PDBQTWriterLegacy
import pandas as pd
import os
import logging

def prepare_ligand(sdf_path, pdbqt_path, compound_id):
    """Full ligand preparation pipeline."""
    # Load molecule
    supplier = Chem.SDMolSupplier(sdf_path, removeHs=False)
    mol = supplier[0]
    if mol is None:
        raise ValueError(f"Failed to parse {sdf_path}")

    # Add hydrogens if missing
    mol = Chem.AddHs(mol)

    # Check/generate 3D coordinates
    conf = mol.GetConformer()
    coords = conf.GetPositions()
    if all(c[2] == 0 for c in coords):
        logging.info(f"{compound_id}: 2D only — generating 3D coordinates")
        AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())

    # Energy minimize
    AllChem.MMFFOptimizeMolecule(mol, maxIters=500)

    # Prepare PDBQT with Meeko
    preparator = MoleculePreparation()
    mol_setups = preparator.prepare(mol)
    for setup in mol_setups:
        pdbqt_string, is_ok, error_msg = PDBQTWriterLegacy.write_string(setup)
        if is_ok:
            with open(pdbqt_path, 'w') as f:
                f.write(pdbqt_string)
            return True
        else:
            logging.error(f"{compound_id}: Meeko error — {error_msg}")

    return False

# Main
compounds = pd.read_csv("configs/compounds.csv")
os.makedirs("ligands/prepared", exist_ok=True)

results = []
for _, row in compounds.iterrows():
    cid = row['compound_id']
    sdf = f"ligands/raw/{cid}.sdf"
    pdbqt = f"ligands/prepared/{cid}.pdbqt"
    try:
        success = prepare_ligand(sdf, pdbqt, cid)
        mol = Chem.SDMolSupplier(sdf)[0]
        mw = Descriptors.MolWt(mol) if mol else None
        results.append({'compound_id': cid, 'status': 'OK' if success else 'FAIL', 'mw': mw})
    except Exception as e:
        logging.error(f"{cid}: {e}")
        results.append({'compound_id': cid, 'status': f'ERROR: {e}', 'mw': None})

pd.DataFrame(results).to_csv("logs/ligand_prep_report.csv", index=False)
```

### Step 4.2: Validation
- Every compound must have a `.pdbqt` in `ligands/prepared/`
- Verify PDBQT has torsion tree (ROOT/BRANCH/ENDROOT keywords)
- Log atom count, torsions, molecular weight per compound

---

## Phase 5 — Receptor Preparation

### Step 5.1: Script `scripts/05_prep_receptors.py`

```python
"""
Prepare receptors: PDB → clean → PDBQT
- Remove water molecules (HOH)
- Remove co-crystallized ligands (save separately for binding site ref)
- Remove non-standard residues / ions
- Add polar hydrogens
- Convert to PDBQT using OpenBabel Python API (openbabel-wheel)
"""
import os
import pandas as pd
import logging

def clean_pdb(input_pdb, output_pdb, ligand_pdb=None, keep_chain=None):
    """Clean PDB: remove waters, extract ligand, keep protein only."""
    protein_lines = []
    ligand_lines = []

    with open(input_pdb) as f:
        for line in f:
            record = line[:6].strip()
            if record in ('ATOM', 'TER'):
                chain = line[21]
                if keep_chain and chain != keep_chain:
                    continue
                protein_lines.append(line)
            elif record == 'HETATM':
                resname = line[17:20].strip()
                if resname == 'HOH':
                    continue  # skip water
                if resname in ('MSE', 'CSE'):
                    protein_lines.append(line)  # selenomethionine etc
                else:
                    ligand_lines.append(line)  # ligands, ions, cofactors

    with open(output_pdb, 'w') as f:
        f.writelines(protein_lines)
        f.write('END\n')

    if ligand_pdb and ligand_lines:
        with open(ligand_pdb, 'w') as f:
            f.writelines(ligand_lines)

    return len(protein_lines)

def pdb_to_pdbqt_receptor(input_pdb, output_pdbqt):
    """Convert cleaned PDB to PDBQT using OpenBabel Python API.

    Uses openbabel-wheel (PyPI) — no CLI binary needed.
    """
    from openbabel import openbabel as ob

    conv = ob.OBConversion()
    conv.SetInFormat("pdb")
    conv.SetOutFormat("pdbqt")
    conv.AddOption("r", ob.OBConversion.OUTOPTIONS)  # rigid receptor
    conv.AddOption("h", ob.OBConversion.OUTOPTIONS)  # add polar H

    mol = ob.OBMol()
    if not conv.ReadFile(mol, input_pdb):
        logging.error(f"OpenBabel failed to read {input_pdb}")
        return False

    # Add hydrogens and compute Gasteiger charges
    mol.AddHydrogens(True)  # polar only
    charge_model = ob.OBChargeModel.FindType("gasteiger")
    if charge_model:
        charge_model.ComputeCharges(mol)

    if not conv.WriteFile(mol, output_pdbqt):
        logging.error(f"OpenBabel failed to write {output_pdbqt}")
        return False
    return True

# Main
targets = pd.read_csv("configs/targets.csv")
os.makedirs("receptors/prepared", exist_ok=True)

for _, row in targets.iterrows():
    tid = row['target_id']
    pdb_id = row['primary_pdb']
    raw = f"receptors/raw/{tid}_{pdb_id}.pdb"
    clean = f"receptors/prepared/{tid}_clean.pdb"
    ligand = f"receptors/prepared/{tid}_cocrystal_ligand.pdb"
    pdbqt = f"receptors/prepared/{tid}.pdbqt"

    if not os.path.exists(raw):
        logging.error(f"Missing: {raw}")
        continue

    # Clean
    n_atoms = clean_pdb(raw, clean, ligand)
    logging.info(f"{tid}: {n_atoms} protein atoms retained")

    # Convert to PDBQT
    success = pdb_to_pdbqt_receptor(clean, pdbqt)
    logging.info(f"{'✓' if success else '✗'} {tid} → PDBQT")
```

---

## Phase 6 — Define Binding Sites

### Step 6.1: Script `scripts/06_define_binding_sites.py`

```python
"""
Define docking grid boxes for each target.

Strategy:
1. If co-crystallized ligand exists → center box on ligand centroid
2. If no ligand → use known active site residues from literature
3. Box size: 25×25×25 Å (standard for blind/semi-blind screening)
"""
import json
import numpy as np
import os
import logging

def get_ligand_centroid(ligand_pdb):
    """Calculate centroid of co-crystallized ligand."""
    coords = []
    with open(ligand_pdb) as f:
        for line in f:
            if line.startswith(('ATOM', 'HETATM')):
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                coords.append([x, y, z])
    if not coords:
        return None
    return np.mean(coords, axis=0).tolist()

# Known active site centers (fallback, from literature)
KNOWN_SITES = {
    'BTUB': {'center': [17.0, 12.0, 33.0], 'size': [25, 25, 25], 'source': 'colchicine site'},
    'IMPB': {'center': [0.0, 0.0, 0.0], 'size': [30, 30, 30], 'source': 'importin groove'},
    'PAK1': {'center': [15.0, 45.0, 25.0], 'size': [25, 25, 25], 'source': 'ATP binding site'},
    'GLU1': {'center': [0.0, 0.0, 0.0], 'size': [25, 25, 25], 'source': 'central cavity'},
    'STA3': {'center': [0.0, 0.0, 0.0], 'size': [25, 25, 25], 'source': 'SH2 domain'},
    'MDM2': {'center': [0.0, 0.0, 0.0], 'size': [25, 25, 25], 'source': 'p53 binding cleft'},
    'WFZD': {'center': [0.0, 0.0, 0.0], 'size': [25, 25, 25], 'source': 'CRD domain'},
    'TRXR': {'center': [0.0, 0.0, 0.0], 'size': [25, 25, 25], 'source': 'C-terminal active site'},
}

binding_sites = {}

for tid, fallback in KNOWN_SITES.items():
    ligand_file = f"receptors/prepared/{tid}_cocrystal_ligand.pdb"

    if os.path.exists(ligand_file) and os.path.getsize(ligand_file) > 0:
        centroid = get_ligand_centroid(ligand_file)
        if centroid:
            binding_sites[tid] = {
                'center_x': round(centroid[0], 2),
                'center_y': round(centroid[1], 2),
                'center_z': round(centroid[2], 2),
                'size_x': 25,
                'size_y': 25,
                'size_z': 25,
                'source': 'co-crystallized ligand centroid'
            }
            logging.info(f"{tid}: binding site from co-crystal ligand at {centroid}")
            continue

    # Fallback to literature values
    # AGENT: If fallback center is [0,0,0], this is a PLACEHOLDER.
    # You MUST look up the actual active site residues and calculate the center
    # from the PDB file, or use a cavity detection approach.
    binding_sites[tid] = {
        'center_x': fallback['center'][0],
        'center_y': fallback['center'][1],
        'center_z': fallback['center'][2],
        'size_x': fallback['size'][0],
        'size_y': fallback['size'][1],
        'size_z': fallback['size'][2],
        'source': fallback['source']
    }
    if fallback['center'] == [0.0, 0.0, 0.0]:
        logging.warning(f"{tid}: PLACEHOLDER center — must be determined from structure!")

with open("configs/binding_sites.json", 'w') as f:
    json.dump(binding_sites, f, indent=2)
```

> **CRITICAL AGENT INSTRUCTION:** Any binding site with center `[0, 0, 0]` is a placeholder. The agent MUST compute the actual center by:
> 1. Finding known active-site residues from literature/UniProt annotations
> 2. Extracting those residue coordinates from the PDB
> 3. Computing their centroid
> **Do not run docking with placeholder coordinates.**

---

## Phase 7 — Run AutoDock Vina

### Step 7.1: Script `scripts/07_run_autodock_vina.py`

```python
"""
Run AutoDock Vina for all compound-target pairs.
Uses the Vina CLI binary via subprocess (the python-vina binding
does not build on arm64 macOS).
Saves poses and scores.
"""
import subprocess
import re
import pandas as pd
import json
import os
import logging

logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s %(levelname)s %(message)s',
                    handlers=[logging.StreamHandler(),
                              logging.FileHandler("logs/execution.log", mode='a')])

# Load configs
pairs = pd.read_csv("configs/pairs.csv")
binding_sites = json.load(open("configs/binding_sites.json"))

os.makedirs("results/vina", exist_ok=True)

def parse_vina_log(log_text):
    """Parse Vina CLI stdout/log for binding affinities.

    Vina outputs a table like:
    -----+------------+----------+----------
     mode |   affinity | dist from best mode
          | (kcal/mol) | rmsd l.b.| rmsd u.b.
    -----+------------+----------+----------
       1       -8.3         0.0       0.0
       2       -7.9         1.234     2.456
    """
    scores = []
    for line in log_text.splitlines():
        match = re.match(r'\s+(\d+)\s+([-\d.]+)\s+([\d.]+)\s+([\d.]+)', line)
        if match:
            scores.append({
                'mode': int(match.group(1)),
                'affinity': float(match.group(2)),
                'rmsd_lb': float(match.group(3)),
                'rmsd_ub': float(match.group(4)),
            })
    return scores

all_scores = []

for _, pair in pairs.iterrows():
    cid = pair['compound_id']
    tid = pair['target_id']
    pair_id = f"{cid}_{tid}"

    receptor_pdbqt = f"receptors/prepared/{tid}.pdbqt"
    ligand_pdbqt = f"ligands/prepared/{cid}.pdbqt"
    output_pdbqt = f"results/vina/{pair_id}_poses.pdbqt"
    log_file = f"results/vina/{pair_id}.log"

    if not os.path.exists(receptor_pdbqt) or not os.path.exists(ligand_pdbqt):
        logging.error(f"Missing files for {pair_id}")
        continue

    bs = binding_sites.get(tid)
    if not bs:
        logging.error(f"No binding site for {tid}")
        continue

    cmd = [
        'vina',
        '--receptor', receptor_pdbqt,
        '--ligand', ligand_pdbqt,
        '--center_x', str(bs['center_x']),
        '--center_y', str(bs['center_y']),
        '--center_z', str(bs['center_z']),
        '--size_x', str(bs['size_x']),
        '--size_y', str(bs['size_y']),
        '--size_z', str(bs['size_z']),
        '--exhaustiveness', '32',
        '--num_modes', '10',
        '--out', output_pdbqt,
    ]

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)

        # Save raw log
        with open(log_file, 'w') as f:
            f.write(result.stdout)
            if result.stderr:
                f.write("\n--- STDERR ---\n")
                f.write(result.stderr)

        if result.returncode != 0:
            raise RuntimeError(f"Vina exited with code {result.returncode}: {result.stderr[:200]}")

        # Parse scores from stdout
        scores = parse_vina_log(result.stdout)
        if scores:
            best_score = scores[0]['affinity']
            all_scores.append({
                'compound_id': cid,
                'target_id': tid,
                'vina_score_kcal': round(best_score, 2),
                'n_poses': len(scores),
                'exhaustiveness': 32,
                'status': 'OK'
            })
            logging.info(f"✓ {pair_id}: {best_score:.2f} kcal/mol ({len(scores)} poses)")
        else:
            raise RuntimeError("No scores parsed from Vina output")

    except subprocess.TimeoutExpired:
        logging.error(f"✗ {pair_id}: timeout (600s)")
        all_scores.append({
            'compound_id': cid, 'target_id': tid,
            'vina_score_kcal': None, 'n_poses': 0,
            'exhaustiveness': 32, 'status': 'TIMEOUT'
        })
    except Exception as e:
        logging.error(f"✗ {pair_id}: {e}")
        all_scores.append({
            'compound_id': cid,
            'target_id': tid,
            'vina_score_kcal': None,
            'n_poses': 0,
            'exhaustiveness': 32,
            'status': f'ERROR: {e}'
        })

# Save all scores
scores_df = pd.DataFrame(all_scores)
scores_df.to_csv("results/vina_scores.csv", index=False)
logging.info(f"Vina docking complete: {len(scores_df)} pairs processed")
```

**Key parameters:**
- `exhaustiveness=32` — higher than default (8) for publication quality
- `n_poses=10` — keep top 10 poses per pair
- Scoring: more negative = stronger predicted binding

---

## Phase 8 — Run Gnina

### Step 8.1: Script `scripts/08_run_gnina.py`

```python
"""
Run Gnina for all compound-target pairs.
Gnina provides both Vina-like scores AND CNN-based scores.

If running via Docker:
    docker run --rm -v $(pwd):/data gnina/gnina gnina \
        -r /data/receptors/prepared/{tid}.pdbqt \
        -l /data/ligands/prepared/{cid}.pdbqt \
        --center_x X --center_y Y --center_z Z \
        --size_x 25 --size_y 25 --size_z 25 \
        --exhaustiveness 32 \
        --num_modes 10 \
        -o /data/results/gnina/{pair_id}_poses.sdf \
        --log /data/results/gnina/{pair_id}.log
"""
import subprocess
import pandas as pd
import json
import os
import re
import logging

pairs = pd.read_csv("configs/pairs.csv")
binding_sites = json.load(open("configs/binding_sites.json"))
os.makedirs("results/gnina", exist_ok=True)

all_scores = []
project_root = os.getcwd()

for _, pair in pairs.iterrows():
    cid = pair['compound_id']
    tid = pair['target_id']
    pair_id = f"{cid}_{tid}"
    bs = binding_sites.get(tid, {})

    receptor = f"receptors/prepared/{tid}.pdbqt"
    ligand = f"ligands/prepared/{cid}.pdbqt"
    output = f"results/gnina/{pair_id}_poses.sdf"
    log_path = f"results/gnina/{pair_id}.log"

    if not os.path.exists(receptor) or not os.path.exists(ligand):
        logging.error(f"Missing files for {pair_id}")
        continue

    # Build Docker command
    cmd = [
        'docker', 'run', '--rm',
        '-v', f'{project_root}:/data',
        'gnina/gnina', 'gnina',
        '-r', f'/data/{receptor}',
        '-l', f'/data/{ligand}',
        '--center_x', str(bs['center_x']),
        '--center_y', str(bs['center_y']),
        '--center_z', str(bs['center_z']),
        '--size_x', str(bs['size_x']),
        '--size_y', str(bs['size_y']),
        '--size_z', str(bs['size_z']),
        '--exhaustiveness', '32',
        '--num_modes', '10',
        '-o', f'/data/{output}',
        '--log', f'/data/{log_path}'
    ]

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)

        # Parse log for scores
        if os.path.exists(log_path):
            with open(log_path) as f:
                log_text = f.read()
            # Extract best CNN score and affinity
            # Gnina outputs: mode | affinity | CNN_score | CNN_affinity
            lines = [l for l in log_text.split('\n') if re.match(r'\s+\d+\s+', l)]
            if lines:
                parts = lines[0].split()
                affinity = float(parts[1])
                cnn_score = float(parts[2]) if len(parts) > 2 else None
                cnn_affinity = float(parts[3]) if len(parts) > 3 else None
            else:
                affinity = cnn_score = cnn_affinity = None
        else:
            affinity = cnn_score = cnn_affinity = None

        all_scores.append({
            'compound_id': cid,
            'target_id': tid,
            'gnina_affinity_kcal': affinity,
            'gnina_cnn_score': cnn_score,
            'gnina_cnn_affinity': cnn_affinity,
            'status': 'OK' if affinity else 'PARSE_ERROR'
        })

        logging.info(f"✓ {pair_id}: affinity={affinity}, CNN={cnn_score}")

    except subprocess.TimeoutExpired:
        logging.error(f"✗ {pair_id}: timeout (600s)")
        all_scores.append({
            'compound_id': cid, 'target_id': tid,
            'gnina_affinity_kcal': None, 'gnina_cnn_score': None,
            'gnina_cnn_affinity': None, 'status': 'TIMEOUT'
        })
    except Exception as e:
        logging.error(f"✗ {pair_id}: {e}")
        all_scores.append({
            'compound_id': cid, 'target_id': tid,
            'gnina_affinity_kcal': None, 'gnina_cnn_score': None,
            'gnina_cnn_affinity': None, 'status': f'ERROR: {e}'
        })

scores_df = pd.DataFrame(all_scores)
scores_df.to_csv("results/gnina_scores.csv", index=False)
```

> **AGENT: If Docker/Gnina is not available**, skip this phase entirely. Update the analysis script to work with Vina-only results and note this as a limitation.

---

## Phase 9 — Analyze Results

### Step 9.1: Script `scripts/09_analyze_results.py`

```python
"""
Merge Vina and Gnina scores, compute consensus rankings,
identify top hits, and generate summary report.
"""
import pandas as pd
import numpy as np
import json
import os

# Load scores
vina = pd.read_csv("results/vina_scores.csv")
compounds = pd.read_csv("configs/compounds.csv")
targets = pd.read_csv("configs/targets.csv")

# Load Gnina if available
gnina_path = "results/gnina_scores.csv"
has_gnina = os.path.exists(gnina_path)
if has_gnina:
    gnina = pd.read_csv(gnina_path)
    merged = vina.merge(gnina, on=['compound_id', 'target_id'], suffixes=('_vina', '_gnina'))
else:
    merged = vina.copy()

# ---- Scoring Thresholds ----
VINA_HIT_THRESHOLD = -7.0          # kcal/mol (strong binder)
VINA_MODERATE_THRESHOLD = -6.0     # kcal/mol (moderate)
CNN_SCORE_THRESHOLD = 0.5          # Gnina CNN confidence (0-1)

# ---- Rank compounds per target ----
merged['vina_rank'] = merged.groupby('target_id')['vina_score_kcal'].rank()

if has_gnina:
    merged['gnina_rank'] = merged.groupby('target_id')['gnina_affinity_kcal'].rank()
    merged['consensus_rank'] = (merged['vina_rank'] + merged['gnina_rank']) / 2
    merged['is_consensus_hit'] = (
        (merged['vina_score_kcal'] <= VINA_HIT_THRESHOLD) &
        (merged['gnina_cnn_score'] >= CNN_SCORE_THRESHOLD)
    )
else:
    merged['consensus_rank'] = merged['vina_rank']
    merged['is_consensus_hit'] = merged['vina_score_kcal'] <= VINA_HIT_THRESHOLD

# ---- Add compound/target metadata ----
merged = merged.merge(compounds[['compound_id', 'name', 'drug_class']], on='compound_id')
merged = merged.merge(targets[['target_id', 'name', 'gene']], on='target_id',
                       suffixes=('_compound', '_target'))

# ---- Sort by consensus rank ----
merged = merged.sort_values('consensus_rank')

# ---- Save full results ----
merged.to_csv("results/consensus_scores.csv", index=False)

# ---- Pivot table: compound × target heatmap data ----
heatmap_vina = merged.pivot_table(
    index='name_compound', columns='name_target',
    values='vina_score_kcal'
)
heatmap_vina.to_csv("results/heatmap_vina_data.csv")

if has_gnina:
    heatmap_gnina = merged.pivot_table(
        index='name_compound', columns='name_target',
        values='gnina_cnn_score'
    )
    heatmap_gnina.to_csv("results/heatmap_gnina_data.csv")

# ---- Polypharmacology: compounds hitting multiple targets ----
hit_counts = merged[merged['is_consensus_hit']].groupby('name_compound')['target_id'].count()
hit_counts = hit_counts.sort_values(ascending=False)
hit_counts.to_csv("results/polypharmacology_hits.csv")

# ---- Generate summary report ----
report = []
report.append("# Molecular Docking Results Summary\n")
report.append(f"**Date:** {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M')}\n")
report.append(f"**Total pairs evaluated:** {len(merged)}\n")
report.append(f"**Docking engines:** AutoDock Vina" +
              (" + Gnina (consensus)\n" if has_gnina else " (single engine)\n"))
report.append(f"**Consensus hits (Vina ≤ {VINA_HIT_THRESHOLD} kcal/mol" +
              (f", CNN ≥ {CNN_SCORE_THRESHOLD}): " if has_gnina else "): ") +
              f"{merged['is_consensus_hit'].sum()}\n")

report.append("\n## Top 10 Compound-Target Pairs\n")
report.append("| Rank | Compound | Target | Vina (kcal/mol) |" +
              (" CNN Score | CNN Affinity |" if has_gnina else "") + " Class |\n")
report.append("|------|----------|--------|-----------------|" +
              ("-----------|--------------|" if has_gnina else "") + "-------|\n")

for i, (_, row) in enumerate(merged.head(10).iterrows(), 1):
    line = f"| {i} | {row['name_compound']} | {row['name_target']} | {row['vina_score_kcal']} |"
    if has_gnina:
        line += f" {row.get('gnina_cnn_score', 'N/A')} | {row.get('gnina_cnn_affinity', 'N/A')} |"
    line += f" {row['drug_class']} |"
    report.append(line + "\n")

report.append("\n## Polypharmacology Candidates\n")
report.append("Compounds hitting multiple targets:\n\n")
for compound, count in hit_counts.items():
    report.append(f"- **{compound}**: {count} targets\n")

report.append("\n## Hit Classification\n")
report.append(f"- Strong binders (≤ {VINA_HIT_THRESHOLD} kcal/mol): "
              f"{(merged['vina_score_kcal'] <= VINA_HIT_THRESHOLD).sum()}\n")
report.append(f"- Moderate binders ({VINA_HIT_THRESHOLD} to {VINA_MODERATE_THRESHOLD} kcal/mol): "
              f"{ ((merged['vina_score_kcal'] > VINA_HIT_THRESHOLD) & (merged['vina_score_kcal'] <= VINA_MODERATE_THRESHOLD)).sum()}\n")
report.append(f"- Weak/no binding (> {VINA_MODERATE_THRESHOLD} kcal/mol): "
              f"{(merged['vina_score_kcal'] > VINA_MODERATE_THRESHOLD).sum()}\n")

with open("results/summary_report.md", 'w') as f:
    f.writelines(report)

print("Analysis complete. See results/summary_report.md")
```

---

## Phase 10 — Generate Figures

### Step 10.1: Script `scripts/10_generate_figures.py`

```python
"""
Generate publication-quality figures for the docking study.
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import seaborn as sns
import os
import json

# Style
plt.rcParams.update({
    'font.size': 12,
    'font.family': 'sans-serif',
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight'
})

os.makedirs("figures", exist_ok=True)
merged = pd.read_csv("results/consensus_scores.csv")
has_gnina = 'gnina_cnn_score' in merged.columns

# ---- Figure 1: Vina Binding Affinity Heatmap ----
fig, ax = plt.subplots(figsize=(12, 8))
heatmap_data = merged.pivot_table(index='name_compound', columns='name_target',
                                   values='vina_score_kcal')
sns.heatmap(heatmap_data, annot=True, fmt='.1f', cmap='RdYlGn',
            center=-6.0, linewidths=0.5, ax=ax,
            cbar_kws={'label': 'Binding Affinity (kcal/mol)'})
ax.set_title('AutoDock Vina Binding Affinities\n(more negative = stronger binding)')
ax.set_xlabel('Target Protein')
ax.set_ylabel('Compound')
plt.tight_layout()
plt.savefig("figures/heatmap_vina.png")
plt.close()

# ---- Figure 2: Gnina CNN Score Heatmap (if available) ----
if has_gnina:
    fig, ax = plt.subplots(figsize=(12, 8))
    heatmap_cnn = merged.pivot_table(index='name_compound', columns='name_target',
                                      values='gnina_cnn_score')
    sns.heatmap(heatmap_cnn, annot=True, fmt='.2f', cmap='YlOrRd',
                vmin=0, vmax=1, linewidths=0.5, ax=ax,
                cbar_kws={'label': 'CNN Pose Score (0-1)'})
    ax.set_title('Gnina CNN Binding Scores\n(higher = more drug-like pose)')
    ax.set_xlabel('Target Protein')
    ax.set_ylabel('Compound')
    plt.tight_layout()
    plt.savefig("figures/heatmap_gnina.png")
    plt.close()

# ---- Figure 3: Vina vs Gnina Scatter (if available) ----
if has_gnina:
    fig, ax = plt.subplots(figsize=(8, 8))
    valid = merged.dropna(subset=['vina_score_kcal', 'gnina_affinity_kcal'])
    scatter = ax.scatter(valid['vina_score_kcal'], valid['gnina_affinity_kcal'],
                         c=valid['gnina_cnn_score'], cmap='coolwarm',
                         s=60, edgecolors='black', linewidth=0.5, alpha=0.8)
    plt.colorbar(scatter, label='CNN Score')

    # Add diagonal
    lims = [min(ax.get_xlim()[0], ax.get_ylim()[0]),
            max(ax.get_xlim()[1], ax.get_ylim()[1])]
    ax.plot(lims, lims, 'k--', alpha=0.3, label='y=x')

    # Annotate top hits
    for _, row in valid.nsmallest(5, 'vina_score_kcal').iterrows():
        ax.annotate(f"{row['name_compound']}\n→{row['name_target']}",
                    (row['vina_score_kcal'], row['gnina_affinity_kcal']),
                    fontsize=7, ha='center')

    ax.set_xlabel('AutoDock Vina Affinity (kcal/mol)')
    ax.set_ylabel('Gnina Affinity (kcal/mol)')
    ax.set_title('Vina vs Gnina: Cross-Validation')
    ax.legend()
    plt.tight_layout()
    plt.savefig("figures/scatter_vina_vs_gnina.png")
    plt.close()

# ---- Figure 4: Top Hits Bar Plot ----
fig, ax = plt.subplots(figsize=(14, 6))
top20 = merged.nsmallest(20, 'vina_score_kcal')
top20['pair_label'] = top20['name_compound'] + ' → ' + top20['name_target']
colors = {'Anthelmintic': '#e74c3c', 'Antiparasitic': '#e67e22',
          'Antifungal': '#2ecc71', 'Antidiabetic': '#3498db',
          'Antimalarial': '#9b59b6', 'Antirheumatic': '#f39c12',
          'Other': '#95a5a6'}
bar_colors = [colors.get(cls, '#95a5a6') for cls in top20['drug_class']]

bars = ax.barh(range(len(top20)), top20['vina_score_kcal'], color=bar_colors,
               edgecolor='black', linewidth=0.5)
ax.set_yticks(range(len(top20)))
ax.set_yticklabels(top20['pair_label'], fontsize=9)
ax.set_xlabel('Binding Affinity (kcal/mol)')
ax.set_title('Top 20 Compound-Target Pairs by Vina Score')
ax.axvline(x=-7.0, color='red', linestyle='--', alpha=0.5, label='Strong binding threshold')
ax.invert_yaxis()

# Legend
from matplotlib.patches import Patch
legend_elements = [Patch(facecolor=c, label=l) for l, c in colors.items()]
ax.legend(handles=legend_elements, loc='lower right', fontsize=8)
plt.tight_layout()
plt.savefig("figures/top_hits_barplot.png")
plt.close()

# ---- Figure 5: Polypharmacology Radar ----
hit_matrix = merged[merged['is_consensus_hit']].pivot_table(
    index='name_compound', columns='name_target',
    values='vina_score_kcal', aggfunc='count', fill_value=0
)
if len(hit_matrix) > 0:
    top_compounds = hit_matrix.sum(axis=1).nlargest(6).index.tolist()
    targets_list = hit_matrix.columns.tolist()
    N = len(targets_list)
    angles = np.linspace(0, 2 * np.pi, N, endpoint=False).tolist()
    angles += angles[:1]  # close the polygon

    fig, ax = plt.subplots(figsize=(8, 8), subplot_kw=dict(polar=True))
    for compound in top_compounds:
        values = hit_matrix.loc[compound].tolist()
        values += values[:1]
        ax.plot(angles, values, 'o-', label=compound, linewidth=1.5)
        ax.fill(angles, values, alpha=0.1)

    ax.set_xticks(angles[:-1])
    ax.set_xticklabels(targets_list, fontsize=9)
    ax.set_title('Polypharmacology Profiles\n(targets hit per compound)', pad=20)
    ax.legend(loc='upper right', bbox_to_anchor=(1.3, 1.1), fontsize=8)
    plt.tight_layout()
    plt.savefig("figures/polypharmacology_radar.png")
    plt.close()

# ---- Figure 6: Consensus Heatmap (combined metric) ----
if has_gnina:
    # Normalize both scores to 0-1 scale and average
    from sklearn.preprocessing import MinMaxScaler
    # If sklearn not available, do manual normalization
    try:
        from sklearn.preprocessing import MinMaxScaler
        scaler = MinMaxScaler()
    except ImportError:
        pass

    # Simple normalization without sklearn
    vina_norm = (merged['vina_score_kcal'] - merged['vina_score_kcal'].max()) / \
                (merged['vina_score_kcal'].min() - merged['vina_score_kcal'].max())
    cnn_norm = merged['gnina_cnn_score']  # already 0-1
    merged['consensus_metric'] = (vina_norm + cnn_norm) / 2

    fig, ax = plt.subplots(figsize=(12, 8))
    consensus_heatmap = merged.pivot_table(index='name_compound', columns='name_target',
                                            values='consensus_metric')
    sns.heatmap(consensus_heatmap, annot=True, fmt='.2f', cmap='YlOrRd',
                linewidths=0.5, ax=ax,
                cbar_kws={'label': 'Consensus Score (0-1, higher = better)'})
    ax.set_title('Consensus Binding Score\n(Vina affinity + Gnina CNN, normalized)')
    plt.tight_layout()
    plt.savefig("figures/heatmap_consensus.png")
    plt.close()

print("All figures saved to figures/")
```

---

## Execution Order Summary

```
PHASE 0  ▸ Environment Setup (uv only — no conda/pip)
             ├─ 0.1  Create pyproject.toml
             ├─ 0.2  uv venv + uv pip install (rdkit, openbabel-wheel, meeko)
             ├─ 0.3  Install Vina CLI binary (GitHub release)
             ├─ 0.4  Install Gnina (Docker) — or decide to skip
             ├─ 0.5  Verify all tools (`uv run python -c ...`, `vina --version`)
             └─ 0.6  Create directory structure

PHASE 1  ▸ Data Registry
             ├─ 1.1  Write configs/compounds.csv
             ├─ 1.2  Write configs/targets.csv (validate PDB IDs via RCSB API)
             └─ 1.3  Generate configs/pairs.csv (all 104 pairs)

PHASE 2  ▸ Fetch Ligands
             ├─ 2.1  Download SDF from PubChem
             └─ 2.2  Validate all SDFs (RDKit parseable, 3D coords)

PHASE 3  ▸ Fetch Receptors
             ├─ 3.1  Download PDB from RCSB (with quality checks)
             └─ 3.2  Validate all PDBs (resolution, format, organism)

PHASE 4  ▸ Prepare Ligands
             ├─ 4.1  SDF → add H → minimize → PDBQT (RDKit + Meeko)
             └─ 4.2  Validate all PDBQTs

PHASE 5  ▸ Prepare Receptors
             ├─ 5.1  PDB → remove water/ligands → add H → PDBQT (OpenBabel)
             └─ 5.2  Save co-crystal ligands separately

PHASE 6  ▸ Define Binding Sites
             ├─ 6.1  Compute grid boxes from co-crystal ligands or literature
             └─ 6.2  Verify NO placeholder [0,0,0] coordinates remain

PHASE 7  ▸ AutoDock Vina Docking
             ├─ 7.1  Run all 104 pairs (exhaustiveness=32, 10 poses each)
             └─ 7.2  Parse scores → results/vina_scores.csv

PHASE 8  ▸ Gnina Docking (if available)
             ├─ 8.1  Run all 104 pairs via Docker
             └─ 8.2  Parse scores → results/gnina_scores.csv

PHASE 9  ▸ Analysis
             ├─ 9.1  Merge scores, compute consensus rankings
             ├─ 9.2  Identify hits (Vina ≤ -7.0, CNN ≥ 0.5)
             ├─ 9.3  Polypharmacology profiling
             └─ 9.4  Generate results/summary_report.md

PHASE 10 ▸ Figures
             ├─ 10.1  Vina heatmap (compound × target)
             ├─ 10.2  Gnina CNN heatmap
             ├─ 10.3  Vina vs Gnina scatter
             ├─ 10.4  Top 20 hits bar plot
             ├─ 10.5  Polypharmacology radar chart
             └─ 10.6  Consensus heatmap
```

---

## Error Handling Rules for the Agent

1. **If a PubChem CID returns no 3D SDF:** Fetch 2D SDF, then generate 3D with RDKit `AllChem.EmbedMolecule()` + MMFF minimization.

2. **If a PDB ID is invalid or unavailable:** Search RCSB by gene name for the best alternative (human, X-ray, ≤ 2.5 Å, with ligand).

3. **If Meeko fails on a compound:** Fall back to OpenBabel Python API:
   ```python
   from openbabel import openbabel as ob
   conv = ob.OBConversion()
   conv.SetInFormat("sdf"); conv.SetOutFormat("pdbqt")
   mol = ob.OBMol(); conv.ReadFile(mol, "input.sdf")
   mol.AddHydrogens(); conv.WriteFile(mol, "output.pdbqt")
   ```

4. **If Vina CLI fails on a pair:** Log the error, skip, and continue. Do NOT stop the entire run.

5. **If Gnina/Docker is unavailable:** Proceed with Vina-only analysis. This is still publishable — note as a limitation.

6. **If a binding site is [0,0,0]:** STOP. Compute the real center from active site residues before docking.

7. **All errors must be logged** to `logs/execution.log` with timestamps.

8. **At completion:** Run `uv run python scripts/09_analyze_results.py` and `uv run python scripts/10_generate_figures.py` even if some pairs failed. The scripts handle NaN values.

9. **All scripts MUST be run via `uv run`** to ensure the correct virtualenv is used:
   ```bash
   uv run python scripts/02_fetch_ligands.py
   uv run python scripts/03_fetch_receptors.py
   # ... etc
   ```

---

## Expected Final Outputs (in CWD)

| File | Description |
|---|---|
| `results/vina_scores.csv` | All Vina binding energies (104 rows) |
| `results/gnina_scores.csv` | All Gnina scores (if run) |
| `results/consensus_scores.csv` | Merged, ranked results with metadata |
| `results/summary_report.md` | Human-readable summary with top hits |
| `results/heatmap_vina_data.csv` | Pivot table for Vina heatmap |
| `results/polypharmacology_hits.csv` | Multi-target hit compounds |
| `figures/heatmap_vina.png` | Compound × target heatmap (Vina) |
| `figures/heatmap_gnina.png` | Compound × target heatmap (Gnina CNN) |
| `figures/heatmap_consensus.png` | Combined consensus heatmap |
| `figures/scatter_vina_vs_gnina.png` | Cross-validation scatter |
| `figures/top_hits_barplot.png` | Top 20 pairs ranked |
| `figures/polypharmacology_radar.png` | Multi-target profiles |
| `logs/execution.log` | Full execution log |
| `logs/ligand_prep_report.csv` | Ligand preparation status |

---

## Paper-Ready Metadata

After execution, the agent should also generate a `results/methods_snippet.md` containing:
- Exact software versions used
- Exact parameters (exhaustiveness, box size, scoring function)
- Number of pairs attempted vs successful
- Hit criteria used
- Any deviations from this plan

This text can be directly adapted into the Methods section of the manuscript.
