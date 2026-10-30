from pathlib import Path
import pandas as pd,numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=Path(__file__).resolve().parents[1];O=R/'output';O.mkdir(exist_ok=True)
s=pd.read_csv(R/'reference/large_n_diagnostics.csv');rates=pd.read_csv(R/'reference/large_n_tau_zero.csv')
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'savefig.bbox':'tight'})
ns=[500,1000,2000,4000,10000];colors=['#167c80','#cc6138','#685ba3'];figs=[]
def setup(ax):
    ax.set_xscale('log');ax.set_xticks(ns);ax.set_xticklabels(['500','1,000','2,000','4,000','10,000'],rotation=25)
    ax.minorticks_off();ax.grid(axis='y',alpha=.17);ax.set_xlabel(r'Target sample size $n_T$')
def save(fig,name):
    fig.savefig(O/(name+'.pdf'));fig.savefig(O/(name+'.png'),dpi=190);figs.append(fig)
fig,axs=plt.subplots(1,3,figsize=(12.6,3.8),sharey=True)
for j in range(3):
    x=s[s.horizon_index.eq(j+1)].sort_values('n');ax=axs[j]
    ax.plot(x.n,x.A_norm_median,color=colors[0],marker='o',lw=1.8)
    ax.set_yscale('log');setup(ax);ax.set_title(['Q25','Q50','Q75'][j])
axs[0].set_ylabel(r'Median $\|S_n\|_2/\sqrt{n_T}$')
fig.suptitle('Penalty magnitude across target sample sizes',y=1.04,fontsize=14);fig.tight_layout();save(fig,'penalty_norm')

fig,axs=plt.subplots(1,3,figsize=(12.6,4),sharey=True)
for j in range(3):
    x=s[s.horizon_index.eq(j+1)].sort_values('n');ax=axs[j]
    ax.errorbar(x.n,x.obs_DST_TG,yerr=1.96*x.obs_DST_TG_mcse,color=colors[0],marker='o',capsize=3,lw=1.7,label=r'Observed, $\tau=0$')
    ax.errorbar(x.n,x.D,yerr=1.96*x.D_mcse,color=colors[1],marker='s',ls='--',capsize=3,lw=1.7,label=r'Plug-in $D(\hat{k})$')
    ax.axhline(0,color='#666666',lw=.8);setup(ax);ax.set_title(['Q25','Q50','Q75'][j])
axs[0].set_ylabel(r'$n_T\{R_{\rm DST}-R_{\rm T-GAM}\}$')
fig.legend(*axs[0].get_legend_handles_labels(),loc='upper center',bbox_to_anchor=(.53,1.03),ncol=2,frameon=False)
fig.tight_layout(rect=[0,0,1,.93]);save(fig,'interpolation_overlay')

fig,axs=plt.subplots(1,3,figsize=(12.6,3.8),sharey=True)
for j in range(3):
    for k,m in enumerate(['SPT-CV-half','DST-CV-half']):
        x=rates[rates.horizon_index.eq(j+1)&rates.method.eq(m)].sort_values('target_n')
        axs[j].plot(x.target_n,x.percent,color=colors[k],marker='o' if k==0 else 's',lw=1.8,label=m.split('-')[0])
    setup(axs[j]);axs[j].set_ylim(-3,103);axs[j].set_title(['Q25','Q50','Q75'][j])
axs[0].set_ylabel(r'CV selecting $\tau=0$ (%)');axs[0].legend(frameon=False,loc='upper left')
fig.suptitle('Automatic rejection of the soft-label channel',y=1.04,fontsize=14);fig.tight_layout();save(fig,'tau_zero_selection')

for fig in figs:plt.close(fig)
