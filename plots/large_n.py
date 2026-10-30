"""Typeset the approved 60-point plot at a readable manuscript width."""
from pathlib import Path
import hashlib,json
import pandas as pd,numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=Path(__file__).resolve().parents[1];O=R/'output';O.mkdir(exist_ok=True)
source=R/'reference/large_n_summary.csv';s=pd.read_csv(source)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})
fig,axs=plt.subplots(1,3,figsize=(9.6,4.15));plotted=[]
styles={'DST-tau0':('#cc6138','-','o',r'DST, $\tau=0$'),'DST-CV-half':('#cc6138','--','s','DST, CV'),
 'SPT-tau0':('#167c80','-','o',r'SPT, $\tau=0$'),'SPT-CV-half':('#167c80','--','s','SPT, CV')}
for ax,metric,title in zip(axs,['kl_risk','brier','auc'],['KL risk','Brier score','Time-dependent AUC']):
 for method,(color,ls,mk,label) in styles.items():
  z=s[s.method.eq(method)].sort_values('target_n')
  line=ax.errorbar(z.target_n,z['scaled_'+metric],yerr=1.96*z['scaled_'+metric+'_mcse'],color=color,ls=ls,marker=mk,lw=1.55,ms=4,capsize=2.5,label=label)
  assert np.array_equal(line.lines[0].get_ydata(),z['scaled_'+metric].to_numpy())
  plotted.extend(dict(metric=metric,method=method,n=int(n),mean=float(v),mcse=float(se)) for n,v,se in zip(z.target_n,z['scaled_'+metric],z['scaled_'+metric+'_mcse']))
 ax.axhline(0,color='#666666',lw=.8);ax.set_xscale('log');ax.minorticks_off()
 ax.set_xticks([500,1000,2000,4000,10000]);ax.set_xticklabels(['500','1,000','2,000','4,000','10,000'],rotation=35,fontsize=9)
 ax.set_xlabel(r'Target sample size $n_T$',fontsize=10);ax.set_ylabel(r'$n_T\times$(method $-$ T-GAM)',fontsize=10)
 ax.set_title(title,fontsize=11);ax.grid(axis='y',alpha=.17);ax.tick_params(axis='y',labelsize=9)
fig.legend(*axs[0].get_legend_handles_labels(),loc='upper center',ncol=4,frameon=False,bbox_to_anchor=(.52,1.015),fontsize=10)
fig.tight_layout(rect=[0,0,1,.90])
fig.savefig(O/'sample_size_risk_differences.pdf',bbox_inches='tight')
fig.savefig(O/'sample_size_risk_differences.png',bbox_inches='tight',dpi=180)
pd.DataFrame(plotted).to_csv(O/'plotted_points.csv',index=False)
(O/'plot_audit.json').write_text(json.dumps({'points':len(plotted),'metrics':['kl_risk','brier','auc'],'source':str(source.relative_to(R)),
 'sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'data_unchanged':True,'errorbars':'1.96 paired MCSE','format':'1x3, readable manuscript fonts'},indent=2)+'\n')
