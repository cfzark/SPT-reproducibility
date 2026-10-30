from pathlib import Path
import os,csv,json,statistics as st
os.environ.setdefault('MPLCONFIGDIR','/tmp/spt_source_clinical_20261027')
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator,FormatStrFormatter
from matplotlib.lines import Line2D
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'output';OUT.mkdir(exist_ok=True)
rows=list(csv.DictReader((ROOT/'reference/clinical_plot_replications.csv').open()));data={(int(r['horizon']),r['method'],int(r['repeat_id'])):{m:float(r[m]) for m in ['brier','uno','auc']} for r in rows}
methods=['T-GAM','Source GAM','Stratified GAM','SPT-GAM','SPT-Cox','SPT-NN','Score recalibration','CoxTL','CoxKL']
colors={'SPT-GAM':'#78c4bb','SPT-Cox':'#a3d5cf','SPT-NN':'#c4e4df'}
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11.5,'axes.labelsize':11.5,'axes.titlesize':13,'xtick.labelsize':10,'ytick.labelsize':11.5,'pdf.fonttype':42,'ps.fonttype':42,'axes.unicode_minus':False})
fig,axs=plt.subplots(3,3,figsize=(9,9.5),sharey=True)
fig.subplots_adjust(left=.24,right=.985,bottom=.115,top=.94,wspace=.22,hspace=.31)
rowlimits={};limits={};summary=[]
for metric in ['uno','auc']:
 v=[x[metric] for x in data.values()];lo,hi=min(v),max(v);pad=(hi-lo)*.05;rowlimits[metric]=(lo-pad,hi+pad)
for i,metric in enumerate(['brier','uno','auc']):
 for j,h in enumerate([24,30,36]):
  ax=axs[i,j];values=[np.array([data[h,m,rep][metric] for rep in range(1,101)]) for m in methods]
  bp=ax.boxplot(values,orientation='horizontal',positions=np.arange(1,10),widths=.54,patch_artist=True,whis=1.5,showfliers=True,medianprops=dict(color='#25313b',linewidth=1.15),whiskerprops=dict(color='#737e87',linewidth=.85),capprops=dict(color='#737e87',linewidth=.85),flierprops=dict(marker='.',markersize=2.3,markerfacecolor='#87929b',markeredgecolor='none',alpha=.65))
  for m,box in zip(methods,bp['boxes']):box.set(facecolor=colors.get(m,'#e0e5e8'),edgecolor='#146b62' if m in colors else '#737e87',linewidth=1 if m in colors else .85)
  ax.scatter([v.mean() for v in values],np.arange(1,10),marker='D',s=14,c=['#095d53' if m in colors else '#35434d' for m in methods],zorder=4)
  ax.set_ylim(9.65,.35);ax.set_yticks(np.arange(1,10),methods);ax.tick_params(axis='y',length=0,labelleft=j==0)
  for tick in ax.get_yticklabels():
   if tick.get_text() in colors:tick.set_color('#095d53');tick.set_fontweight('bold')
  if metric in rowlimits:lim=rowlimits[metric];ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
  else:
   v=np.concatenate(values);pad=(max(v)-min(v))*.09;lim=(min(v)-pad,max(v)+pad);ax.xaxis.set_major_formatter(FormatStrFormatter('%.4f'))
  ax.set_xlim(lim);limits[f'{metric}_{h}']=lim;ax.xaxis.set_major_locator(MaxNLocator(nbins=3));ax.grid(axis='x',color='#e8ecee',lw=.5);ax.set_axisbelow(True)
  ax.spines[['top','right','left']].set_visible(False);ax.spines['bottom'].set_color('#bdc6cc')
  if i==0:ax.set_title(f'{h} months',pad=15)
  if j==1:ax.set_xlabel(['IPCW Brier score',r'$C_U(h)$','AUC'][i],labelpad=8)
  for m,v in zip(methods,values):summary.append(dict(horizon=h,method=m,metric=metric,n=100,mean=float(v.mean()),sd=st.stdev(v)))
fig.legend(handles=[Line2D([0],[0],marker='D',linestyle='none',color='#35434d',markersize=5,label='Mean'),Line2D([0],[0],color='#25313b',linewidth=1.4,label='Median')],loc='lower center',bbox_to_anchor=(.61,.025),ncol=2,frameon=False,fontsize=10)
fig.savefig(OUT/'clinical_boxplots.pdf',bbox_inches='tight');fig.savefig(OUT/'clinical_boxplots.png',dpi=170,bbox_inches='tight');plt.close(fig)
with (OUT/'plot_summary.csv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=list(summary[0]));w.writeheader();w.writerows(summary)
