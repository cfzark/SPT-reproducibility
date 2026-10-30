"""Restore ONLY the 18 currently displayed table bodies; no experimental fitting."""
from pathlib import Path
import json
root=Path(__file__).resolve().parent;out=root/'output';out.mkdir(exist_ok=True)
tables=json.loads((root/'reference/displayed_tables.json').read_text())
for doc in ['main.tex','supplement.tex']:
 t=[x for x in tables if x['document']==doc]
 (out/(doc.replace('.tex','_tables.tex'))).write_text('% Exact current table bodies, including original means/rounding/bold flags.\n'+'\n\n'.join(x['tex']for x in t)+'\n')
print('Restored',len(tables),'table bodies. These are frozen published displays, not summaries of new fits.')
