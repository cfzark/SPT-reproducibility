from pathlib import Path
import csv,collections as co,statistics as st,json,hashlib,os
os.environ.setdefault('MPLCONFIGDIR','/tmp/spt_3x4_mpl_20261027')
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from matplotlib.lines import Line2D
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'output';OUT.mkdir(exist_ok=True)
methods=[('NoRes-cal__OneSE-half','SPT'),('Teacher-cal__OneSE-half','PR-GAM'),('Rho1-cal__OneSE-half',r'Fixed $\rho=1$'),('NoRes-raw__OneSE-half','Uncalibrated'),('NoRes-cal__fixed0','No soft labels'),('SingleWeight-TauCV-half','DST')]
scenarios=['level','shape','reverse','D'];metrics=['brier_relative','uno_pp','auc_pp']
rows=list(csv.DictReader((ROOT/'reference/ablation_contrasts.csv').open()))
values={(r['scenario'],r['method'],r['metric'],int(r['rep'])):float(r['contrast'])for r in rows}
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':13,'axes.titlesize':13.5,'axes.labelsize':13,'xtick.labelsize':12,'ytick.labelsize':14,'pdf.fonttype':42,'ps.fonttype':42,'axes.unicode_minus':False})
fig,axes=plt.subplots(3,4,figsize=(10.6,8.2),sharey=True)
fig.subplots_adjust(left=.15,right=.99,bottom=.13,top=.93,wspace=.24,hspace=.46)
limits={}
oldlims={'brier_relative': [-7.910802898371238, 16.984198958749218], 'uno_pp': [-7.871798177122378, 3.4694234704885965], 'auc_pp': [-9.269350983124182, 3.9980834853181495]}
for i,metric in enumerate(metrics):
 for j,sc in enumerate(scenarios):
  ax=axes[i,j];arrays=[np.array([values[sc,m,metric,rep] for rep in range(1,101)]) for m,_ in methods]
  bp=ax.boxplot(arrays,orientation='horizontal',positions=np.arange(1,7),widths=.53,patch_artist=True,whis=1.5,showfliers=True,medianprops=dict(color='#25313b',linewidth=1.15),whiskerprops=dict(color='#737e87',linewidth=.85),capprops=dict(color='#737e87',linewidth=.85),flierprops=dict(marker='.',markersize=2.3,markerfacecolor='#87929b',markeredgecolor='none',alpha=.65))
  for k,b in enumerate(bp['boxes']):b.set(facecolor='#78c4bb' if k==0 else '#e0e5e8',edgecolor='#146b62' if k==0 else '#737e87',linewidth=1.1 if k==0 else .85)
  ax.scatter([v.mean() for v in arrays],np.arange(1,7),marker='D',s=14,c=['#095d53']+['#35434d']*5,zorder=4)
  ax.axvline(0,color='#53616c',ls='--',lw=.85,zorder=0)
  lo=min(0,min(v.min() for v in arrays));hi=max(0,max(v.max() for v in arrays));pad=(hi-lo)*.06
  lim=(lo-pad,hi+pad)
  if sc=='reverse':
   assert oldlims[metric][0]<=lo and oldlims[metric][1]>=hi
   lim=oldlims[metric]
  limits[sc+'__'+metric]=list(lim);ax.set_xlim(lim);ax.set_ylim(6.65,.35)
  ax.set_yticks(np.arange(1,7),[name for _,name in methods]);ax.tick_params(axis='y',length=0,labelleft=j==0)
  for tick in ax.get_yticklabels():
   if tick.get_text()=='SPT':tick.set_fontweight('bold');tick.set_color('#095d53')
  if j==2:ax.set_xticks({'brier_relative':[-5,0,5,10,15],'uno_pp':[-6,-3,0,3],'auc_pp':[-9,-6,-3,0,3]}[metric])
  else:ax.xaxis.set_major_locator(MaxNLocator(nbins=4))
  ax.grid(axis='x',color='#e8ecee',lw=.5,zorder=0);ax.spines[['top','right','left']].set_visible(False);ax.spines['bottom'].set_color('#bdc6cc')
  if i==0:ax.set_title(['Level shift','Shape change','Direction reversal','Structural mismatch\n(setting D)'][j],pad=15)
 y=axes[i,0].get_position().y0-.045
 fig.text(.57,y,['Relative Brier (%)',r'$\Delta C_U(h)$ (pp)','AUC change (pp)'][i],ha='center',va='top',fontsize=13)
fig.legend(handles=[Line2D([0],[0],marker='D',linestyle='none',color='#35434d',markersize=5,label='Mean'),Line2D([0],[0],color='#25313b',linewidth=1.4,label='Median'),Line2D([0],[0],color='#53616c',linestyle='--',linewidth=.85,label='No change from T-GAM')],loc='lower center',bbox_to_anchor=(.57,.017),ncol=3,frameon=False,fontsize=12)
fig.savefig(OUT/'source_mismatch_boxplots.pdf',bbox_inches='tight')
fig.savefig(OUT/'source_mismatch_boxplots.png',dpi=170,bbox_inches='tight');plt.close(fig)
