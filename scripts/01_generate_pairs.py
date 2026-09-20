import pandas as pd

compounds = pd.read_csv("configs/compounds.csv")
targets = pd.read_csv("configs/targets.csv")

# High priority known interactions
high_priority = {
    ('FEN', 'BTUB'): 'Direct primary tubulin inhibitor target',
    ('FEN', 'GLU1'): 'Published glucose uptake / GLUT inhibition',
    ('FEN', 'MDM2'): 'p53 pathway activation & stabilization',
    ('IVE', 'IMPB'): 'Direct importin-alpha/beta1 nuclear transport inhibitor',
    ('IVE', 'PAK1'): 'Direct PAK1 kinase inhibitor',
    ('IVE', 'WFZD'): 'WNT/beta-catenin pathway modulation',
    ('MEB', 'BTUB'): 'Direct tubulin colchicine-site binder',
    ('MEB', 'MDM2'): 'p53 signaling activation in cancer cells',
    ('ALB', 'BTUB'): 'Tubulin microtubule polymerization inhibitor',
    ('ALB', 'GLU1'): 'GLUT inhibition / metabolic stress',
    ('NIC', 'STA3'): 'Direct STAT3 phosphorylation & dimerization inhibitor',
    ('NIC', 'WFZD'): 'FZD/Wnt pathway degradation and inhibition',
    ('NIT', 'STA3'): 'STAT3 transcription factor pathway inhibition',
    ('NIT', 'MDM2'): 'c-Myc and downstream tumor suppressor regulation',
    ('ITR', 'WFZD'): 'Hedgehog & Wnt pathway signaling inhibition',
    ('MET', 'GLU1'): 'AMPK activation & glucose transport regulation',
    ('DIS', 'TRXR'): 'Thiol-reactive proteasome/redox disruption',
    ('AUR', 'TRXR'): 'Direct thioredoxin reductase 1 inhibitor (gold complex)',
    ('CLQ', 'STA3'): 'Autophagy and STAT3 crosstalk',
}

medium_priority = {
    ('FEN', 'STA3'): 'Downstream transcription factor effect',
    ('IVE', 'STA3'): 'Secondary STAT3 inhibition reported in glioblastoma',
    ('MEB', 'GLU1'): 'Benzimidazole class metabolic effect',
    ('NIC', 'PAK1'): 'Kinase crosstalk in cancer signaling',
    ('PRA', 'BTUB'): 'Anthelmintic class tubulin exploration',
    ('PRA', 'WFZD'): 'Calcium/Wnt signaling modulation',
    ('ITR', 'PAK1'): 'Angiogenesis kinase modulation',
    ('KET', 'TRXR'): 'Cytochrome P450 and redox system interaction',
    ('MET', 'STA3'): 'Metformin suppresses STAT3 phosphorylation',
    ('MET', 'PAK1'): 'AMPK-mediated kinase regulation',
    ('DIS', 'MDM2'): 'Proteasome inhibition stabilizing p53',
    ('CLQ', 'MDM2'): 'Autophagy arrest leading to p53 induction',
}

pairs = []
for _, c in compounds.iterrows():
    for _, t in targets.iterrows():
        cid = c['compound_id']
        tid = t['target_id']
        key = (cid, tid)
        if key in high_priority:
            priority = 'HIGH'
            rationale = high_priority[key]
        elif key in medium_priority:
            priority = 'MEDIUM'
            rationale = medium_priority[key]
        else:
            priority = 'LOW'
            rationale = 'Exploratory cross-screening / polypharmacology evaluation'
        
        pairs.append({
            'compound_id': cid,
            'target_id': tid,
            'compound_name': c['name'],
            'target_name': t['name'],
            'priority': priority,
            'rationale': rationale
        })

df_pairs = pd.DataFrame(pairs)
df_pairs.to_csv("configs/pairs.csv", index=False)
print(f"Generated {len(df_pairs)} pairs ({sum(df_pairs['priority']=='HIGH')} HIGH, {sum(df_pairs['priority']=='MEDIUM')} MEDIUM, {sum(df_pairs['priority']=='LOW')} LOW)")
