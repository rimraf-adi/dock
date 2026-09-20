# Systematic *In Silico* Screening & Polypharmacology Profiling of Repurposed Oncology Therapeutics

[![Preprint](https://img.shields.io/badge/Preprint-PDF-red.svg)](manuscript_preprint.pdf)
[![Python 3.11+](https://img.shields.io/badge/Python-3.11%2B-blue.svg)](https://www.python.org/)
[![Package Manager](https://img.shields.io/badge/uv-Fast%20Python%20Tooling-purple.svg)](https://docs.astral.sh/uv/)
[![Engine](https://img.shields.io/badge/AutoDock%20Vina-v1.2.7-green.svg)](https://vina.scripps.edu/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

A systematic structure-based virtual screen and polypharmacology characterization of **13 clinically approved non-oncology therapeutics** across **8 hallmark cancer target proteins**, grounded in high-resolution crystallographic co-complexes and benchmarked directly against canonical chemotherapy controls.

---

## 🌟 Key Discoveries

1. **Praziquantel (PRA) Emerges as an Unprecedented Microtubule & GLUT1 Disruptor:**
   - **$\beta$-Tubulin:** Achieved an affinity of **$-10.08$ kcal/mol**, surpassing standard clinical chemotherapy controls **Colchicine ($-7.06$ kcal/mol)** and **Nocodazole ($-8.77$ kcal/mol)**.
   - **GLUT1 Transporter:** Bound the central cavity with **$-9.63$ kcal/mol**, surpassing the nanomolar inhibitor **BAY-876 ($-9.21$ kcal/mol)**.
   - **Exceptional Ligand Efficiency:** $LE = 0.438$ kcal/mol/HA and $\text{LipE} = 4.86$, indicating potency driven by spatial complementarity rather than molecular grease.
   - **CNS Ready:** Low polar surface area ($\text{TPSA} = 40.62$ Å²) and proven blood-brain barrier permeability in humans.

2. **Mebendazole (MEB) Displays True Pan-Cancer Polypharmacology:**
   - Bound **8 out of 8 targets** with $\le -7.0$ kcal/mol (mean affinity $-8.09$ kcal/mol).
   - High-affinity engagement on $\beta$-tubulin ($-9.56$ kcal/mol) and MDM2 ($-8.24$ kcal/mol), providing a structural mechanism for p53 reactivation.

3. **Niclosamide (NIC) Recapitulates STAT3 Dimerization Blockade:**
   - Bound the STAT3 SH2 domain with **$-8.44$ kcal/mol**, establishing a dense 4-hydrogen-bond network with the critical Arg609/Ser611/Ser613/Ser636 triad (substantially outperforming tool compound **Stattic at $-5.62$ kcal/mol**).

4. **Internal Negative Control Validation:**
   - **Ivermectin** bound GLUT1 ($-12.17$ kcal/mol) and Importin-$\alpha$ ($-9.47$ kcal/mol), but was sterically excluded from $\beta$-tubulin ($+1.22$ kcal/mol), proving pocket parameterization fidelity.

---

## 🔬 Benchmark Comparison: Repurposed Hits vs Chemotherapy Controls

| Target Protein | Standard Chemotherapy / Benchmark | Benchmark $\Delta G$ | Benchmark $LE$ | Top Repurposed Hit | Repurposed $\Delta G$ | Repurposed $LE$ | Outcome |
|---|---|:---:|:---:|---|:---:|:---:|---|
| **$\beta$-Tubulin** | **Colchicine** | $-7.06$ kcal/mol | 0.243 | **Praziquantel** | **$-10.08$ kcal/mol** | **0.438** | $+3.02$ kcal/mol stronger; $1.8\times$ higher $LE$ |
| **$\beta$-Tubulin** | **Nocodazole** | $-8.77$ kcal/mol | 0.417 | **Mebendazole** | **$-9.56$ kcal/mol** | **0.435** | Outperforms benzimidazole chemo |
| **MDM2** | **Nutlin-3a** | $-8.39$ kcal/mol | 0.210 | **Praziquantel** | **$-8.82$ kcal/mol** | **0.383** | Matches affinity with nearly double $LE$ |
| **STAT3 SH2** | **Stattic** | $-5.62$ kcal/mol | 0.401 | **Niclosamide** | **$-8.44$ kcal/mol** | **0.402** | $+2.82$ kcal/mol stronger binding |
| **GLUT1** | **BAY-876** | $-9.21$ kcal/mol | 0.279 | **Praziquantel** | **$-9.63$ kcal/mol** | **0.419** | Surpasses nanomolar metabolic lead |

---

## 📊 Publication Figures Gallery (300 DPI)

| Figure | Description | File Link |
|---|---|---|
| **Figure 1** | Calculated Binding Affinity Matrix Heatmap | [`figures/heatmap_vina.png`](figures/heatmap_vina.png) |
| **Figure 2** | Multi-Target Polypharmacology Radar Chart | [`figures/polypharmacology_radar.png`](figures/polypharmacology_radar.png) |
| **Figure 3** | Ligand Efficiency ($LE$) and MW Decoupling | [`figures/ligand_efficiency.png`](figures/ligand_efficiency.png) |
| **Figure 4** | Chemotherapy Baseline Benchmark Comparison | [`figures/repurposed_vs_chemo_baselines.png`](figures/repurposed_vs_chemo_baselines.png) |
| **Figure 5** | Unsupervised Hierarchical Biclustering | [`figures/hierarchical_clustermap.png`](figures/hierarchical_clustermap.png) |
| **Figure 6** | Target Selectivity Z-Score Matrix | [`figures/target_selectivity_zscores.png`](figures/target_selectivity_zscores.png) |
| **Figure 7** | Atom-Level Residue Contact Networks & H-Bonds | [`figures/residue_contacts_top_hits.png`](figures/residue_contacts_top_hits.png) |
| **Figure 8** | Pharmacological Class Distributions & Jitter | [`figures/drug_class_boxplots.png`](figures/drug_class_boxplots.png) |

---

## 🌐 Interactive 3D WebGL Molecular Viewer

We have built a single-file, zero-dependency 3D molecular viewer using **3Dmol.js**:
- **File:** [`results/interactive_viewer.html`](results/interactive_viewer.html)
- **Usage:** Simply open the file in Safari, Chrome, or Firefox.
- **Features:**
  - Select between top complexes (`PRA_BTUB`, `MEB_BTUB`, `FEN_BTUB`, `NIC_STA3`, `KET_TRXR`, `MEB_MDM2`).
  - Superimpose docked drug poses (cyan) over co-crystallized reference inhibitors (yellow).
  - Inspect binding pocket surfaces, residue contact sticks ($< 4.0$ Å), and hydrogen bonds.

---

## 📑 Preprint Manuscript & Future Roadmaps

- **Compiled PDF Preprint (9 Pages):** [`manuscript_preprint.pdf`](manuscript_preprint.pdf)
- **LaTeX Source Code:** [`manuscript_preprint.tex`](manuscript_preprint.tex)
- **Deep Learning Technical Roadmap:** [`DEEP_LEARNING_ROADMAP.md`](DEEP_LEARNING_ROADMAP.md) (Covers GNINA 3D-CNN rescoring, ESM-2 zero-shot DTI, DiffSBDD 3D generative scaffold morphing, and DepMap transcriptomic synergism).

---

## 🚀 Quickstart & Reproducibility

This project uses [`uv`](https://github.com/astral-sh/uv) for fast, reproducible dependency management.

```bash
# 1. Clone repository
git clone https://github.com/rimraf-adi/dock.git
cd dock

# 2. Synchronize environment with uv
uv sync

# 3. Re-run screening pipeline
uv run python scripts/07_run_autodock_vina.py
uv run python scripts/09_analyze_results.py
uv run python scripts/10_generate_figures.py
uv run python scripts/11_advanced_figures.py
uv run python scripts/12_residue_contacts.py
uv run python scripts/15_benchmark_chemo_controls.py

# 4. Compile LaTeX manuscript
pdflatex manuscript_preprint.tex
```

---

## ⚖️ License
Distributed under the MIT License. See `LICENSE` for details.
