"""
Generate a standalone interactive 3D WebGL molecular viewer using 3Dmol.js.
Produces results/interactive_viewer.html which can be opened in any web browser
to inspect docked poses in 3D, compare against co-crystallized reference inhibitors,
and visualize key residue interactions.
"""

import os
import json
import base64

PAIRS = [
    {
        "id": "PRA_BTUB",
        "compound": "Praziquantel",
        "target": "β-Tubulin",
        "pdb_id": "4O2B",
        "affinity": "-10.08 kcal/mol",
        "le": "0.438 kcal/mol/HA",
        "pocket": "Colchicine pocket at α/β-tubulin interface",
        "rec_file": "receptors/prepared/BTUB_clean.pdb",
        "lig_file": "results/vina/PRA_BTUB_poses.pdbqt",
        "ref_file": "receptors/prepared/BTUB_cocrystal_ligand.pdb",
        "key_residues": ["ASN 258", "LEU 248", "LYS 352", "VAL 238", "LEU 255", "ILE 318", "THR 179"]
    },
    {
        "id": "MEB_BTUB",
        "compound": "Mebendazole",
        "target": "β-Tubulin",
        "pdb_id": "4O2B",
        "affinity": "-9.56 kcal/mol",
        "le": "0.435 kcal/mol/HA",
        "pocket": "Colchicine pocket at α/β-tubulin interface",
        "rec_file": "receptors/prepared/BTUB_clean.pdb",
        "lig_file": "results/vina/MEB_BTUB_poses.pdbqt",
        "ref_file": "receptors/prepared/BTUB_cocrystal_ligand.pdb",
        "key_residues": ["LYS 254", "THR 145", "GLN 11", "ASN 101", "LEU 248"]
    },
    {
        "id": "FEN_BTUB",
        "compound": "Fenbendazole",
        "target": "β-Tubulin",
        "pdb_id": "4O2B",
        "affinity": "-9.10 kcal/mol",
        "le": "0.433 kcal/mol/HA",
        "pocket": "Colchicine pocket at α/β-tubulin interface",
        "rec_file": "receptors/prepared/BTUB_clean.pdb",
        "lig_file": "results/vina/FEN_BTUB_poses.pdbqt",
        "ref_file": "receptors/prepared/BTUB_cocrystal_ligand.pdb",
        "key_residues": ["ASN 258", "LEU 248", "ALA 250", "LEU 255"]
    },
    {
        "id": "NIC_STA3",
        "compound": "Niclosamide",
        "target": "STAT3",
        "pdb_id": "6NJS",
        "affinity": "-8.44 kcal/mol",
        "le": "0.402 kcal/mol/HA",
        "pocket": "SH2 domain phosphotyrosine binding pocket",
        "rec_file": "receptors/prepared/STA3_clean.pdb",
        "lig_file": "results/vina/NIC_STA3_poses.pdbqt",
        "ref_file": "receptors/prepared/STA3_cocrystal_ligand.pdb",
        "key_residues": ["ARG 609", "SER 611", "SER 613", "SER 636", "GLU 612"]
    },
    {
        "id": "KET_TRXR",
        "compound": "Ketoconazole",
        "target": "TXNRD1",
        "pdb_id": "2ZZB",
        "affinity": "-10.17 kcal/mol",
        "le": "0.282 kcal/mol/HA",
        "pocket": "NADPH/FAD interfacial domain",
        "rec_file": "receptors/prepared/TRXR_clean.pdb",
        "lig_file": "results/vina/KET_TRXR_poses.pdbqt",
        "ref_file": "receptors/prepared/TRXR_cocrystal_ligand.pdb",
        "key_residues": ["HIS 108", "TYR 116", "GLU 477", "TRP 407"]
    },
    {
        "id": "MEB_MDM2",
        "compound": "Mebendazole",
        "target": "MDM2",
        "pdb_id": "4HG7",
        "affinity": "-8.24 kcal/mol",
        "le": "0.375 kcal/mol/HA",
        "pocket": "p53 transactivation hydrophobic cleft",
        "rec_file": "receptors/prepared/MDM2_clean.pdb",
        "lig_file": "results/vina/MEB_MDM2_poses.pdbqt",
        "ref_file": "receptors/prepared/MDM2_cocrystal_ligand.pdb",
        "key_residues": ["LEU 54", "ILE 61", "VAL 93", "HIS 96", "TYR 100"]
    }
]

def extract_model1_pdb(pdbqt_path):
    """Convert MODEL 1 of PDBQT to standard PDB format string."""
    lines = []
    if not os.path.exists(pdbqt_path):
        return ""
    with open(pdbqt_path, 'r') as f:
        in_m1 = False
        for line in f:
            if line.startswith('MODEL 1'):
                in_m1 = True
                continue
            if line.startswith('ENDMDL') and in_m1:
                break
            if in_m1 and (line.startswith('ATOM') or line.startswith('HETATM')):
                # Convert to clean PDB line format
                lines.append(line[:66] + "\n")
    return "".join(lines)

def build_html():
    complexes_data = {}
    for p in PAIRS:
        rec_str = ""
        if os.path.exists(p["rec_file"]):
            with open(p["rec_file"], 'r') as f:
                rec_str = f.read()

        ref_str = ""
        if os.path.exists(p["ref_file"]):
            with open(p["ref_file"], 'r') as f:
                ref_str = f.read()

        lig_str = extract_model1_pdb(p["lig_file"])

        complexes_data[p["id"]] = {
            "meta": p,
            "receptor_pdb": rec_str,
            "ref_ligand_pdb": ref_str,
            "docked_ligand_pdb": lig_str
        }

    html_content = f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <title>Repurposed Oncological Drug Docking — 3D Interactive Viewer</title>
  <script src="https://cdnjs.cloudflare.com/ajax/libs/3Dmol/2.0.4/3Dmol-min.js"></script>
  <style>
    :root {{
      --bg-dark: #121820;
      --panel-bg: #1c2430;
      --border-color: #2c384a;
      --accent: #00d2d3;
      --accent-hover: #01a3a4;
      --text-main: #f5f6fa;
      --text-muted: #a4b0be;
    }}
    * {{ box-sizing: border-box; margin: 0; padding: 0; }}
    body {{
      font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica, Arial, sans-serif;
      background: var(--bg-dark);
      color: var(--text-main);
      display: flex;
      flex-direction: column;
      height: 100vh;
      overflow: hidden;
    }}
    header {{
      background: var(--panel-bg);
      padding: 12px 24px;
      border-bottom: 1px solid var(--border-color);
      display: flex;
      justify-content: space-between;
      align-items: center;
    }}
    header h1 {{
      font-size: 1.25rem;
      font-weight: 600;
      color: var(--accent);
      display: flex;
      align-items: center;
      gap: 10px;
    }}
    header .subtitle {{
      font-size: 0.85rem;
      color: var(--text-muted);
    }}
    .container {{
      display: flex;
      flex: 1;
      overflow: hidden;
    }}
    #viewer-container {{
      flex: 1;
      position: relative;
      background: #080c10;
    }}
    #viewer {{
      width: 100%;
      height: 100%;
    }}
    .sidebar {{
      width: 380px;
      background: var(--panel-bg);
      border-left: 1px solid var(--border-color);
      display: flex;
      flex-direction: column;
      padding: 20px;
      gap: 16px;
      overflow-y: auto;
    }}
    .card {{
      background: #141b24;
      border-radius: 8px;
      padding: 14px;
      border: 1px solid var(--border-color);
    }}
    .card h3 {{
      font-size: 0.95rem;
      color: var(--accent);
      margin-bottom: 8px;
      text-transform: uppercase;
      letter-spacing: 0.5px;
    }}
    label {{
      display: block;
      font-size: 0.85rem;
      color: var(--text-muted);
      margin-bottom: 6px;
    }}
    select, button {{
      width: 100%;
      padding: 10px 12px;
      border-radius: 6px;
      border: 1px solid var(--border-color);
      background: #1e2836;
      color: var(--text-main);
      font-size: 0.9rem;
      cursor: pointer;
      outline: none;
      transition: border 0.2s;
    }}
    select:focus, button:hover {{
      border-color: var(--accent);
    }}
    .stat-grid {{
      display: grid;
      grid-template-columns: 1fr 1fr;
      gap: 8px;
      margin-top: 8px;
    }}
    .stat-box {{
      background: #0e141c;
      padding: 8px;
      border-radius: 4px;
      text-align: center;
    }}
    .stat-val {{
      font-size: 1.1rem;
      font-weight: 700;
      color: #ff6b6b;
    }}
    .stat-val.le {{
      color: #1dd1a1;
    }}
    .stat-lbl {{
      font-size: 0.75rem;
      color: var(--text-muted);
      text-transform: uppercase;
    }}
    .controls {{
      display: flex;
      flex-direction: column;
      gap: 10px;
    }}
    .btn-row {{
      display: flex;
      gap: 8px;
    }}
    .btn-row button {{
      flex: 1;
    }}
    .legend-item {{
      display: flex;
      align-items: center;
      gap: 10px;
      font-size: 0.85rem;
      margin-bottom: 6px;
    }}
    .legend-color {{
      width: 14px;
      height: 14px;
      border-radius: 3px;
    }}
    .chip {{
      display: inline-block;
      padding: 2px 6px;
      background: #253344;
      border-radius: 4px;
      font-size: 0.75rem;
      margin: 2px;
      font-family: monospace;
    }}
  </style>
</head>
<body>
  <header>
    <div>
      <h1>🔬 Molecular Docking Polypharmacology 3D Explorer</h1>
      <span class="subtitle">Systematic Multi-Target In Silico Screen of Repurposed Non-Oncology Drugs</span>
    </div>
    <div style="font-size: 0.8rem; color: var(--text-muted);">AutoDock Vina Engine • RCSB High-Resolution Crystal Structures</div>
  </header>

  <div class="container">
    <div id="viewer-container">
      <div id="viewer"></div>
    </div>

    <div class="sidebar">
      <div class="card">
        <h3>Select Complex</h3>
        <select id="complexSelect" onchange="loadSelectedComplex()">
          <option value="PRA_BTUB">Praziquantel → β-Tubulin (-10.08 kcal/mol)</option>
          <option value="MEB_BTUB">Mebendazole → β-Tubulin (-9.56 kcal/mol)</option>
          <option value="FEN_BTUB">Fenbendazole → β-Tubulin (-9.10 kcal/mol)</option>
          <option value="NIC_STA3">Niclosamide → STAT3 (-8.44 kcal/mol)</option>
          <option value="KET_TRXR">Ketoconazole → TXNRD1 (-10.17 kcal/mol)</option>
          <option value="MEB_MDM2">Mebendazole → MDM2 (-8.24 kcal/mol)</option>
        </select>
      </div>

      <div class="card">
        <h3>Interaction Metrics</h3>
        <div id="metaInfo">
          <div style="font-weight: 600; font-size: 1.05rem;" id="infoTitle">Praziquantel → β-Tubulin</div>
          <div style="font-size: 0.82rem; color: var(--text-muted); margin-top: 2px;" id="infoPocket">Colchicine pocket at α/β-tubulin interface (PDB: 4O2B)</div>
          <div class="stat-grid">
            <div class="stat-box">
              <div class="stat-val" id="infoAffinity">-10.08 kcal/mol</div>
              <div class="stat-lbl">Binding Affinity</div>
            </div>
            <div class="stat-box">
              <div class="stat-val le" id="infoLE">0.438</div>
              <div class="stat-lbl">Ligand Efficiency</div>
            </div>
          </div>
        </div>
      </div>

      <div class="card">
        <h3>Display Controls</h3>
        <div class="controls">
          <div class="btn-row">
            <button onclick="toggleSurface()">Toggle Pocket Surface</button>
            <button onclick="resetView()">Reset Camera</button>
          </div>
          <div class="btn-row">
            <button onclick="toggleSpin()">Toggle Spin</button>
            <button onclick="zoomToPocket()">Zoom to Pocket</button>
          </div>
        </div>
      </div>

      <div class="card">
        <h3>Legend</h3>
        <div class="legend-item">
          <div class="legend-color" style="background: #00d2d3;"></div>
          <span>Docked Candidate Pose (Cyan Sticks)</span>
        </div>
        <div class="legend-item">
          <div class="legend-color" style="background: #feca57;"></div>
          <span>Co-Crystallized Reference Inhibitor (Yellow Sticks)</span>
        </div>
        <div class="legend-item">
          <div class="legend-color" style="background: #1dd1a1;"></div>
          <span>Pocket Contact Residues &lt; 4.0 Å (Green Sticks)</span>
        </div>
        <div class="legend-item">
          <div class="legend-color" style="background: #576574;"></div>
          <span>Receptor Backbone Ribbon (Gray/Blue)</span>
        </div>
      </div>

      <div class="card">
        <h3>Key Pocket Residues</h3>
        <div id="keyResiduesChips"></div>
      </div>
    </div>
  </div>

  <script>
    const DATA = {json.dumps(complexes_data)};
    let glviewer = null;
    let surfaceObj = null;
    let isSpinning = false;
    let surfaceVisible = false;

    window.onload = function() {{
      let element = document.getElementById('viewer');
      let config = {{ backgroundColor: '#080c10' }};
      glviewer = $3Dmol.createViewer(element, config);
      loadSelectedComplex();
    }};

    function loadSelectedComplex() {{
      let key = document.getElementById('complexSelect').value;
      let item = DATA[key];
      if (!item) return;

      // Update sidebar
      document.getElementById('infoTitle').innerText = item.meta.compound + " → " + item.meta.target;
      document.getElementById('infoPocket').innerText = item.meta.pocket + " (PDB: " + item.meta.pdb_id + ")";
      document.getElementById('infoAffinity').innerText = item.meta.affinity;
      document.getElementById('infoLE').innerText = item.meta.le.split(' ')[0];

      let chips = document.getElementById('keyResiduesChips');
      chips.innerHTML = '';
      item.meta.key_residues.forEach(r => {{
        let span = document.createElement('span');
        span.className = 'chip';
        span.innerText = r;
        chips.appendChild(span);
      }});

      // Render 3D scene
      glviewer.clear();
      surfaceObj = null;
      surfaceVisible = false;

      // Add receptor
      glviewer.addModel(item.receptor_pdb, "pdb");
      glviewer.setStyle({{}}, {{ cartoon: {{ color: '#3d566e', opacity: 0.85 }} }});

      // Highlight contact residues
      let resNums = item.meta.key_residues.map(r => parseInt(r.split(' ')[1])).filter(n => !isNaN(n));
      glviewer.setStyle({{ resi: resNums }}, {{
        cartoon: {{ color: '#1dd1a1' }},
        stick: {{ colorscheme: 'greenCarbon', radius: 0.15 }}
      }});

      // Add reference ligand if present
      if (item.ref_ligand_pdb && item.ref_ligand_pdb.trim().length > 0) {{
        glviewer.addModel(item.ref_ligand_pdb, "pdb");
        glviewer.setStyle({{ model: 1 }}, {{
          stick: {{ colorscheme: 'yellowCarbon', radius: 0.25 }}
        }});
      }}

      // Add docked ligand pose
      let dockedModelIdx = (item.ref_ligand_pdb && item.ref_ligand_pdb.trim().length > 0) ? 2 : 1;
      if (item.docked_ligand_pdb && item.docked_ligand_pdb.trim().length > 0) {{
        glviewer.addModel(item.docked_ligand_pdb, "pdb");
        glviewer.setStyle({{ model: dockedModelIdx }}, {{
          stick: {{ colorscheme: 'cyanCarbon', radius: 0.28 }}
        }});
      }}

      // Center and zoom on docked ligand
      glviewer.zoomTo({{ model: dockedModelIdx }});
      glviewer.render();
    }}

    function toggleSurface() {{
      if (!glviewer) return;
      if (surfaceVisible) {{
        if (surfaceObj) glviewer.removeSurface(surfaceObj);
        surfaceObj = null;
        surfaceVisible = false;
      }} else {{
        surfaceObj = glviewer.addSurface($3Dmol.SurfaceType.VDW, {{
          opacity: 0.35,
          color: '#00d2d3'
        }}, {{ model: 0 }});
        surfaceVisible = true;
      }}
      glviewer.render();
    }}

    function resetView() {{
      if (!glviewer) return;
      let key = document.getElementById('complexSelect').value;
      let item = DATA[key];
      let dockedModelIdx = (item && item.ref_ligand_pdb && item.ref_ligand_pdb.trim().length > 0) ? 2 : 1;
      glviewer.zoomTo({{ model: dockedModelIdx }});
      glviewer.render();
    }}

    function zoomToPocket() {{
      resetView();
    }}

    function toggleSpin() {{
      if (!glviewer) return;
      isSpinning = !isSpinning;
      glviewer.spin(isSpinning);
    }}
  </script>
</body>
</html>
"""
    with open("results/interactive_viewer.html", "w") as f:
        f.write(html_content)
    print("✓ Successfully generated results/interactive_viewer.html")

if __name__ == "__main__":
    build_html()
