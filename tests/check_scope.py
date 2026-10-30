"""Cheap archive audit; no fitting, downloads, or patient-data access."""
from pathlib import Path
import csv,json
P=Path(__file__).resolve().parents[1]
sim=list(csv.DictReader((P/'reference/simulation_summary.csv').open()))
clin=list(csv.DictReader((P/'reference/clinical_summary.csv').open()))
assert len(sim)==900 and len(clin)==51
assert not any(r['method']=='DST' for r in sim)
assert len({r['method'] for r in clin})==17
for bad in ['GG-SrcTP','Unified','NoRes-cal__MinCV','SingleWeight-TauCV-large']:
 assert not any(r['method']==bad for r in sim+clin)
t=json.loads((P/'reference/displayed_tables.json').read_text())
assert sum(x['document']=='main.tex'for x in t)==2
assert sum(x['document']=='supplement.tex'for x in t)==16
assert set(r['scenario']for r in csv.DictReader((P/'reference/ablation_contrasts.csv').open()))=={'level','shape','reverse','D'}
print('Displayed scope: PASS')
