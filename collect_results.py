"""Summarize explicitly listed metric files, without traversing a filesystem."""
from pathlib import Path
import sys,csv,collections,statistics,math
if len(sys.argv)<3:raise SystemExit('Usage: python collect_results.py OUTPUT.csv metrics1.csv [metrics2.csv ...]')
rows=[]
for name in sys.argv[2:]:rows.extend(csv.DictReader(open(name)))
keys=['design','setting','distribution','n_target','p','quantile','horizon','method']
# Simulation horizons vary by replication; group by quantile, not realized time.
groups=collections.defaultdict(list);seen=set()
for r in rows:
 key=tuple(r.get(k,'') if k!='horizon' or r['design']=='clinical' else '' for k in keys)
 ident=key+(r['replication'],)
 if ident in seen:raise ValueError('Duplicate method/replication: '+str(ident))
 seen.add(ident);groups[key].append(r)
out=[]
for key,g in groups.items():
 z=dict(zip(keys,key));z['R']=len(g)
 for m in ['brier','uno','auc','true_brier','true_auc','kl_risk']:
  if m not in g[0]:continue
  x=[float(r[m])for r in g];assert all(math.isfinite(v)for v in x)
  z[m]=statistics.mean(x);sd=statistics.stdev(x)if len(x)>1 else float('nan')
  z[m+('_partition_sd'if key[0]=='clinical'else'_mcse')]=sd if key[0]=='clinical'else sd/math.sqrt(len(x))
 out.append(z)
fields=list(dict.fromkeys(k for r in out for k in r))
with open(sys.argv[1],'w')as f:
 w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(out)
