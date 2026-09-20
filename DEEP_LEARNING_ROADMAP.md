# Deep Learning & Generative AI Roadmap for Computational Oncology & Drug Repurposing

**Document Status:** Standalone Research Proposal & Technical Architecture  
**Accompanying Paper:** *Systematic In Silico Screening and Polypharmacology Profiling of Repurposed Therapeutics Against Eight Hallmark Cancer Targets*  
**Scope:** Advanced Deep Learning architectures, pre-trained biophysical foundation models, and generative chemistry pipelines designed to build upon and extend empirical AutoDock Vina findings.

---

## Executive Summary

While empirical scoring functions (e.g., AutoDock Vina) provide rapid, literature-grounded conformational search, they rely on simplified additive energy terms (Coulombic, Van der Waals, hydrophobic contacts) that cannot model:
1. Multi-body electronic polarization and quantum charge transfer.
2. Solvent desolvation entropy and non-local conformational rearrangement.
3. De novo structural generation beyond existing chemical libraries.

This roadmap outlines **four state-of-the-art Deep Learning (DL) tracks** directly tailored to the findings of our 13-compound $\times$ 8-target screen, specifically focusing on **Praziquantel (PRA)** as a novel microtubule/GLUT1 lead and **Mebendazole (MEB)** as a pan-cancer polypharmacology scaffold.

```
                  ┌─────────────────────────────────────────────────────────┐
                  │    Empirical Docking Screen Findings (AutoDock Vina)     │
                  │   • Praziquantel: β-Tubulin (-10.08 kcal/mol, LE: 0.438)│
                  │   • Mebendazole: 8/8 Targets Bound (Mean -8.09 kcal/mol)│
                  └────────────────────────────┬────────────────────────────┘
                                               │
             ┌───────────────────┬─────────────┴───────┬───────────────────┐
             ▼                   ▼                     ▼                   ▼
    ┌─────────────────┐ ┌─────────────────┐   ┌─────────────────┐ ┌─────────────────┐
    │     Track 1     │ │     Track 2     │   │     Track 3     │ │     Track 4     │
    │  SE(3) Neural   │ │  PLM Zero-Shot  │   │  3D Generative  │ │ Transcriptomic  │
    │    Rescoring    │ │   Cross-Attn    │   │  Scaffold Diff. │ │ GNN Synergism   │
    │  (GNINA/DiffD)  │ │ (ESM-2 + GNN)   │   │  (DiffSBDD)     │ │ (DepMap / TGSA) │
    └─────────────────┘ └─────────────────┘   └─────────────────┘ └─────────────────┘
```

---

## Track 1: SE(3)-Equivariant Neural Rescoring & Generative Docking

### 1.1 The Scientific Rationale
Empirical force fields frequently yield false positives in large, hydrophobic pockets (e.g., the GLUT1 translocation cavity and TXNRD1 active site) because raw contact counting scales with molecular volume. Deep learning models decouple molecular grease from true spatial complementarity.

### 1.2 Model Architectures
1. **3D-Convolutional Neural Networks (GNINA):**
   - **Mechanism:** Voxelizes the protein-ligand binding pocket into a 3D grid ($0.5$ Å resolution) with multi-channel atomic features (aliphatic carbons, donor nitrogens, acceptor oxygens, charges).
   - **Metrics Extracted:**
     - $\text{CNNscore} \in [0, 1]$: Neural probability that the binding pose is within $2.0$ Å RMSD of the experimental ground truth.
     - $\text{CNNaffinity}$: Direct non-linear predicted binding affinity ($pK_d / pIC_{50}$).
   - **Target Pairs to Evaluate:** `PRA_BTUB`, `MEB_BTUB`, `IVE_GLU1`, `NIC_STA3`, `MEB_MDM2`.
2. **SE(3)-Equivariant Diffusion Docking (DiffDock):**
   - **Mechanism:** Formulates docking as a generative reverse diffusion process over the continuous manifold of rigid-body translations, rotations ($SO(3)$), and internal torsion angles ($SO(2)^m$).
   - **Validation Goal:** Test whether an unbiased generative diffusion prior independently converges onto the exact colchicine pocket binding mode identified by AutoDock Vina for Praziquantel.

### 1.3 Execution Blueprint
```bash
# Example GNINA rescoring command for Praziquantel in β-Tubulin
gnina -r receptors/prepared/BTUB_clean.pdb \
      -l results/vina/PRA_BTUB_poses.pdbqt \
      --score_only \
      --cnn_scoring rescore \
      --cnn crossdock_default2018
```

---

## Track 2: Zero-Shot Multi-Modal Drug-Target Interaction (PLM + GNN)

### 2.1 The Scientific Rationale
Can deep learning predict our multi-target polypharmacology results directly from the 1D protein amino acid sequence without requiring 3D crystal structures or manual grid-box definition?

### 2.2 Model Architecture: Cross-Attention Transformer
- **Protein Encoder:** Meta's **ESM-2** (`esm2_t33_650M_UR50D` or `esm2_t36_3B_UR50D`) extracts residue-level context embeddings $\mathbf{H}_P \in \mathbb{R}^{L \times d_p}$ directly from primary FASTA sequences (`TUBB`, `KPNA2`, `PAK1`, `SLC2A1`, `STAT3`, `MDM2`, `FZD4`, `TXNRD1`).
- **Ligand Encoder:** Graph Isomorphism Network (GIN) or ChemBERTa Transformer extracts atom/substructure embeddings $\mathbf{H}_L \in \mathbb{R}^{M \times d_l}$ from SMILES.
- **Cross-Attention Interaction Module (MolTrans / DeepDTI):**
  $$\mathbf{A}_{i,j} = \text{Softmax}\left(\frac{(\mathbf{H}_P \mathbf{W}_Q)(\mathbf{H}_L \mathbf{W}_K)^T}{\sqrt{d_k}}\right)$$
  Computes explicit attention weights between individual amino acids ($i$) and chemical pharmacophores ($j$).

### 2.3 Scientific Deliverable
- A comparative correlation matrix between **1D Zero-Shot PLM Affinity** vs **3D Vina Calculated Free Energy**.
- Attention heatmaps revealing whether the model pays attention to $\beta$-tubulin residues 240–260 (colchicine pocket) when exposed to Praziquantel and Mebendazole.

---

## Track 3: 3D Structure-Conditioned De Novo Generative Optimization (Scaffold Morphing)

### 3.1 The Scientific Rationale
Praziquantel (PRA) is a veterinary/human anthelmintic optimized for parasitic flatworms. It has an outstanding **Ligand Efficiency ($LE = 0.438$ kcal/mol/HA)** on human $\beta$-tubulin, meaning every single atom contributes maximal binding free energy. However, it leaves unoccupied subpockets in the colchicine cleft.

Instead of manual medicinal chemistry iterations, **3D Equivariant Denoising Diffusion** can grow optimized substituents directly inside the crystal cavity.

### 3.2 Model Architectures
- **DiffSBDD (Diffusion for Structure-Based Drug Design)** & **TargetDiff:**
  - Treats atom types and 3D Euclidean coordinates as continuous/categorical diffusion variables conditioned on receptor atoms within 10 Å of the colchicine binding centroid.
  - **Constrained Scaffold Inpainting:**
    - Fix the core tetrahydroisoquinoline and cyclohexyl core of Praziquantel that forms hydrogen bonds with $\beta$-Asn258 and $\beta$-Lys352.
    - Diffuse novel functional groups into the empty pocket volume towards $\alpha$-Thr179, $\alpha$-Val181, and $\beta$-Val318.
- **Output:** 1,000 generated de novo analogs filtered by:
  1. Synthetic accessibility score ($\text{SA} \le 3.5$).
  2. QED drug-likeness ($\ge 0.70$).
  3. Predicted nanomolar $\beta$-tubulin binding affinity ($\le -11.5$ kcal/mol).

---

## Track 4: Transcriptome-Grounded GNNs for Tumor Lineage Sensitivity & Combination Synergism

### 4.1 The Scientific Rationale
Does multi-target engagement translate into in vitro cancer cell lethality, and which cancer types are vulnerable?

### 4.2 Architecture: TGSA (Target-Gene-Set Attention Network)
- **Data Integration:**
  - **Input 1:** Polypharmacological Target Profile: 8-dimensional $-\Delta G$ vector per compound ($[\Delta G_{\text{BTUB}}, \dots, \Delta G_{\text{TRXR}}]$).
  - **Input 2:** DepMap / CCLE (Cancer Cell Line Encyclopedia) baseline RNA-seq and CRISPR dependency scores across 1,000+ patient-derived cell lines.
  - **Input 3:** STRING Protein-Protein Interaction (PPI) knowledge graph.
- **Predictive Objectives:**
  1. **Lineage Specificity:** Predict $IC_{50}$ distributions across cancer lineages (e.g., glioblastoma, pancreatic ductal adenocarcinoma, triple-negative breast cancer).
  2. **Combination Synergy:** Predict Loewe / Bliss synergy scores for combinations:
     - **Praziquantel + Temozolomide** (Glioblastoma).
     - **Mebendazole + Paclitaxel** (Taxane-resistant ovarian/breast cancer).
     - **Niclosamide + Cisplatin** (STAT3-driven solid tumors).

---

## Recommended Software & Environment Setup

To execute this Deep Learning roadmap, the following stack is recommended:

```bash
# Create dedicated GPU-accelerated environment (or MPS for Apple Silicon)
uv venv .venv-dl --python 3.11
source .venv-dl/bin/activate

# Core Deep Learning & GNN frameworks
uv pip install torch torchvision torchaudio
uv pip install torch-geometric torch-scatter torch-cluster
uv pip install transformers fair-esm dgl

# Molecular Generative Models
uv pip install rdkit openbabel-wheel meeko
git clone https://github.com/arneschneuing/DiffSBDD.git
git clone https://github.com/gnina/gnina.git
```

---

## Conclusion & Integration with Wet-Lab Assays

By pairing our completed biophysical docking screen with this Deep Learning framework:
1. **Track 1 (GNINA/DiffDock)** computationally validates poses without human bias.
2. **Track 2 (ESM-2/MolTrans)** scales the screen to the entire human druggable proteome ($>20,000$ proteins).
3. **Track 3 (DiffSBDD)** generates novel, patentable praziquantel analogs.
4. **Track 4 (TGSA/DepMap)** identifies the exact patient-derived cell lines to order for in vitro validation.
