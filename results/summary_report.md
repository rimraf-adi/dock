# Comprehensive Molecular Docking Screening Report
**Screening Date:** 2026-09-20 07:42
**Total Pairs Evaluated:** 96 / 104 (92.3% success)
**Docking Engine:** AutoDock Vina v1.2.7 (Exhaustiveness = 32)
**Strong Hit Threshold:** Binding energy $\le -7.0$ kcal/mol
**Total Identified Strong Hits:** 50 pairs (52.1% hit rate)

## 1. Top 15 Highest Affinity Hits Overall

| Rank | Compound | Drug Class | Target | PDB ID | Binding Affinity (kcal/mol) | Tier | Priority |
|:---:|:---|:---|:---|:---:|:---:|:---:|:---:|
| 1 | **Ivermectin** | Anthelmintic | GLUT1 transporter | `5EQI` | **-12.17** | Very Strong (<= -8.5) | LOW |
| 2 | **Ketoconazole** | Antifungal | β-Tubulin | `4O2B` | **-11.32** | Very Strong (<= -8.5) | LOW |
| 3 | **Itraconazole** | Antifungal | GLUT1 transporter | `5EQI` | **-10.76** | Very Strong (<= -8.5) | LOW |
| 4 | **Ketoconazole** | Antifungal | GLUT1 transporter | `5EQI` | **-10.43** | Very Strong (<= -8.5) | LOW |
| 5 | **Ketoconazole** | Antifungal | Thioredoxin reductase 1 | `2ZZB` | **-10.17** | Very Strong (<= -8.5) | MEDIUM |
| 6 | **Praziquantel** | Anthelmintic | β-Tubulin | `4O2B` | **-10.08** | Very Strong (<= -8.5) | MEDIUM |
| 7 | **Ivermectin** | Anthelmintic | Thioredoxin reductase 1 | `2ZZB` | **-9.779** | Very Strong (<= -8.5) | LOW |
| 8 | **Praziquantel** | Anthelmintic | GLUT1 transporter | `5EQI` | **-9.629** | Very Strong (<= -8.5) | LOW |
| 9 | **Mebendazole** | Anthelmintic | β-Tubulin | `4O2B` | **-9.563** | Very Strong (<= -8.5) | HIGH |
| 10 | **Ivermectin** | Anthelmintic | Importin-α (NLS) | `4WV6` | **-9.473** | Very Strong (<= -8.5) | LOW |
| 11 | **Fenbendazole** | Anthelmintic | β-Tubulin | `4O2B` | **-9.099** | Very Strong (<= -8.5) | UNKNOWN |
| 12 | **Praziquantel** | Anthelmintic | Thioredoxin reductase 1 | `2ZZB` | **-8.985** | Very Strong (<= -8.5) | LOW |
| 13 | **Niclosamide** | Anthelmintic | β-Tubulin | `4O2B` | **-8.981** | Very Strong (<= -8.5) | LOW |
| 14 | **Itraconazole** | Antifungal | Importin-α (NLS) | `4WV6` | **-8.835** | Very Strong (<= -8.5) | LOW |
| 15 | **Praziquantel** | Anthelmintic | MDM2 (p53 pathway) | `4HG7` | **-8.819** | Very Strong (<= -8.5) | LOW |

## 2. Top Binder for Each Target Protein

| Target Protein | Gene | PDB | Co-crystallized Inhibitor Reference | Top Compound | Affinity (kcal/mol) |
|:---|:---:|:---:|:---|:---|:---:|
| **β-Tubulin** | `TUBB` | `4O2B` | Colchicine pocket at alpha-beta tubulin interface | **Ketoconazole** (Antifungal) | **-11.32** |
| **GLUT1 transporter** | `SLC2A1` | `5EQI` | Inward-facing glucose/cytochalasin B binding cavity | **Ivermectin** (Anthelmintic) | **-12.17** |
| **Importin-α (NLS)** | `KPNA2` | `4WV6` | Armadillo-repeat NLS recognition groove | **Ivermectin** (Anthelmintic) | **-9.473** |
| **MDM2 (p53 pathway)** | `MDM2` | `4HG7` | Nutlin-3a small-molecule binding pocket | **Praziquantel** (Anthelmintic) | **-8.819** |
| **PAK1 kinase** | `PAK1` | `2HY8` | ATP-binding cleft between N- and C-terminal lobes | **Ketoconazole** (Antifungal) | **-8.787** |
| **STAT3** | `STAT3` | `6NJS` | SH2 domain phosphotyrosine binding pocket | **Ivermectin** (Anthelmintic) | **-8.648** |
| **Thioredoxin reductase 1** | `TXNRD1` | `2ZZB` | C-terminal selenocysteine/redox active site | **Ketoconazole** (Antifungal) | **-10.17** |
| **Frizzled CRD (WNT)** | `FZD8` | `6TFB` | Druggable lipid groove of Frizzled CRD | **Praziquantel** (Anthelmintic) | **-7.465** |

## 3. Polypharmacology & Multi-Target Repurposing Profile

Compounds demonstrating favorable binding ($\le -7.0$ kcal/mol) across multiple oncological targets:

| Compound | Drug Class | Targets Hit ($\le -7.0$ kcal/mol) | Mean Affinity | Best Affinity |
|:---|:---|:---:|:---:|:---:|
| **Mebendazole** | Anthelmintic | **8** / 8 targets | -8.09 kcal/mol | **-9.56 kcal/mol** |
| **Ivermectin** | Anthelmintic | **7** / 8 targets | -9.16 kcal/mol | **-12.17 kcal/mol** |
| **Praziquantel** | Anthelmintic | **7** / 8 targets | -8.69 kcal/mol | **-10.08 kcal/mol** |
| **Ketoconazole** | Antifungal | **6** / 8 targets | -9.38 kcal/mol | **-11.32 kcal/mol** |
| **Itraconazole** | Antifungal | **5** / 8 targets | -8.71 kcal/mol | **-10.76 kcal/mol** |
| **Fenbendazole** | Anthelmintic | **5** / 8 targets | -8.21 kcal/mol | **-9.1 kcal/mol** |
| **Niclosamide** | Anthelmintic | **5** / 8 targets | -8.02 kcal/mol | **-8.98 kcal/mol** |
| **Nitazoxanide** | Antiparasitic | **4** / 8 targets | -7.79 kcal/mol | **-8.53 kcal/mol** |
| **Albendazole** | Anthelmintic | **2** / 8 targets | -7.16 kcal/mol | **-7.22 kcal/mol** |
| **Chloroquine** | Antimalarial | **1** / 8 targets | -7.66 kcal/mol | **-7.66 kcal/mol** |

## 4. Cross-Class Comparative Pharmacology

| Drug Class | Evaluated Pairs | Mean Affinity (kcal/mol) | Std Dev | Peak Affinity (kcal/mol) |
|:---|:---:|:---:|:---:|:---:|
| **Antifungal** | 16 | -8.22 | 1.76 | **-11.32** |
| **Anthelmintic** | 48 | -7.62 | 1.86 | **-12.17** |
| **Antiparasitic** | 8 | -7.0 | 0.94 | **-8.53** |
| **Antimalarial** | 8 | -6.39 | 0.83 | **-7.66** |
| **Antidiabetic** | 8 | -4.68 | 0.5 | **-5.46** |
| **Other** | 8 | -4.6 | 0.59 | **-5.21** |

## 5. Methodological & Biological Implications

- **Anthelmintics Performance:** Benzimidazoles (fenbendazole, mebendazole, albendazole) and salicylanilides (niclosamide) exhibited high-affinity binding across multiple targets beyond their primary tubulin/STAT3 axes, supporting polypharmacology hypotheses.
- **Antifungals (Azoles):** Itraconazole and ketoconazole display bulky hydrophobic structures that dock favorably into deep hydrophobic cavities including MDM2 and GLUT1.
- **Experimental Validation Priorities:** The prioritized top consensus candidates warrant in vitro validation via tubulin polymerization assays, STAT3 phosphorylation Western blots, and glucose uptake assays.
