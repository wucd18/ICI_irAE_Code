from common import *
from analysis_functions import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import TwoSlopeNorm,ListedColormap
from matplotlib.patches import Rectangle,Patch
from matplotlib.lines import Line2D
import shutil
cfg=jr(ROOT/'analysis_config.json');families=cfg['families'];flabel={'Stem':'Stem','Transit_Amplifying':'TA','Colonocyte':'Colonocyte','Goblet':'Goblet'}
mlabel={UP:'Disease-up',DEP:'Disease-depleted',IFN:'IFN reference',RES:'IFN-excluded residual'}
fcolors=dict(zip(families,['#4B79A5','#C47729','#398D80','#A46991']))
donors=['DNW14','DNW15','DNW9'];markers=dict(zip(donors,['o','s','^']))
plt.rcParams.update({'font.family':'Arial','font.size':8,'axes.titlesize':10,'axes.labelsize':8.5,'xtick.labelsize':8,'ytick.labelsize':8,'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False,'savefig.facecolor':'white'})
F=ROOT/'figures';P=ROOT/'figure_source_data';(ROOT/'review/rendered').mkdir(exist_ok=True)
g=pd.read_csv(ROOT/'results/all_ligand_programme_responses.tsv',sep='\t',float_precision='round_trip');s=pd.read_csv(ROOT/'results/response_pattern_summary.tsv',sep='\t',float_precision='round_trip');d=pd.read_csv(ROOT/'results/all_ligand_donor_deltas.tsv',sep='\t',float_precision='round_trip');bench=pd.read_csv(ROOT/'results/programme_benchmark.tsv',sep='\t');bp=pd.read_csv(ROOT/'results/benchmark_all_pairs.tsv',sep='\t');ex=pd.read_csv(ROOT/'results/prespecified_display_examples.tsv',sep='\t')
def writepage(fig,pdf,name,page):
    pdf.savefig(fig)
    fig.savefig(F/f'{name}_p{page}.svg')
    fig.savefig(ROOT/'review/rendered'/f'{name}_p{page}.png',dpi=170)
    plt.close(fig)
def panel(ax,letter,title):
    ax.set_title(title,loc='left',pad=12,fontweight='normal');ax.text(-.18,1.08,letter,transform=ax.transAxes,fontweight='bold',fontsize=13,va='bottom')
def zero(ax,both=False):
    ax.axhline(0,color='#8c8c8c',lw=.7,ls='--',zorder=0)
    if both:ax.axvline(0,color='#8c8c8c',lw=.7,ls='--',zorder=0)
def donors_legend():return [Line2D([],[],color='#333333',marker=markers[x],linestyle='',label=x,markersize=5) for x in donors]
def donor_panel(ax,data,col):
    zero(ax)
    for i,f in enumerate(families):
        b=data[data.celltype==f]
        for k,r in enumerate(b.itertuples()):ax.scatter(i+(-.12+donors.index(r.donor)*.12),getattr(r,col)*1000,c=fcolors[f],marker=markers[r.donor],s=35,edgecolors='white',linewidths=.5,zorder=3)
        ax.plot([i-.25,i+.25],[b[col].mean()*1000]*2,color='black',lw=1.4,zorder=4)
    ax.set_xticks(range(4),[flabel[f] for f in families],rotation=20,ha='right');ax.set_xlim(-.5,3.5);ax.set_ylabel('Treatment − control score (×10$^{-3}$)')

# Revised Figure 4 uses the atlas primary FDR, never historical targeted FDR.
name='Figure4_TNFSF12';gg=g[(g.ligand=='TNFSF12')&g.module.isin([UP,DEP,IFN])];dd=d[d.ligand=='TNFSF12']
gg.to_csv(P/(name+'_a.tsv'),sep='\t',index=False);dd.to_csv(P/(name+'_bcd.tsv'),sep='\t',index=False)
fig,axs=plt.subplots(2,2,figsize=(7.5,6.7));fig.subplots_adjust(left=.12,right=.98,bottom=.13,top=.94,hspace=.64,wspace=.43)
ax=axs[0,0];mods=[UP,DEP,IFN];v=gg.pivot(index='module',columns='celltype',values='NES').reindex(index=mods,columns=families).to_numpy();lim=np.ceil(np.nanmax(np.abs(v))*2)/2
im=ax.imshow(v,cmap='RdBu_r',vmin=-lim,vmax=lim,aspect='auto')
for iy,m in enumerate(mods):
    for ix,f in enumerate(families):
        r=gg[(gg.module==m)&(gg.celltype==f)].iloc[0];c='white' if abs(r.NES)>.65*lim else '#111111'
        ax.text(ix,iy,f'{r.NES:.2f}\nq={r.FDR_atlas:.2g}',ha='center',va='center',fontsize=7.5,color=c)
        if r.FDR_atlas<.05:ax.add_patch(Rectangle((ix-.47,iy-.47),.94,.94,fill=False,ec='black',lw=1.15))
ax.set_xticks(range(4),[flabel[f] for f in families],rotation=20,ha='right');ax.set_yticks(range(3),['Disease-up','Disease-\ndepleted','IFN']);ax.tick_params(length=0);panel(ax,'a','Enrichment: atlas FDR')
cb=fig.colorbar(im,ax=ax,orientation='horizontal',fraction=.045,pad=.24,aspect=25);cb.set_label('NES; outlined cells: FDR < 0.05',fontsize=7)
donor_panel(axs[0,1],dd,'delta_up');panel(axs[0,1],'b','Disease-up: observed donor pairs')
donor_panel(axs[1,0],dd,'delta_depleted');panel(axs[1,0],'c','Disease-depleted: observed donor pairs')
ax=axs[1,1];zero(ax,True)
for r in dd.itertuples():ax.scatter(-r.delta_up*1000,r.delta_depleted*1000,c=fcolors[r.celltype],marker=markers[r.donor],s=44,ec='white',lw=.5)
ax.set_xlabel('−Δ disease-up score (×10$^{-3}$)');ax.set_ylabel('Δ disease-depleted score (×10$^{-3}$)');panel(ax,'d','Joint changes in the same donor')
fig.legend(handles=donors_legend(),loc='lower center',bbox_to_anchor=(.5,.035),ncol=3,frameon=False)
fig.text(.5,.012,'Three shared donors; black ticks in b,c show means. Directional concordance is not significance.',ha='center',fontsize=7.3)
with PdfPages(F/(name+'.pdf')) as pdf:writepage(fig,pdf,name,1)

# Complete matrices, alphabetically ordered, including NE and negative findings.
name='Figure_AllLigand_Atlas';a=g[g.condition_class=='ligand'];a.to_csv(P/(name+'.tsv'),sep='\t',index=False)
with PdfPages(F/(name+'.pdf')) as pdf:
    for page,start in enumerate([0,40],1):
        ligs=cfg['ligands'][start:start+40];fig,axs=plt.subplots(1,3,figsize=(8.27,11.69));fig.subplots_adjust(left=.13,right=.97,top=.87,bottom=.21,wspace=.19)
        for ax,m in zip(axs,[UP,DEP,IFN]):
            full=a[a.module==m];lim=float(np.nanmax(np.abs(full.loc[full.eligibility=='ELIGIBLE','delta_mean'])));b=full[full.ligand.isin(ligs)]
            matrix=b.pivot(index='ligand',columns='celltype',values='delta_mean').reindex(index=ligs,columns=families).to_numpy()
            mask=b.pivot(index='ligand',columns='celltype',values='eligibility').reindex(index=ligs,columns=families).to_numpy()!='ELIGIBLE';matrix[mask]=np.nan
            cm=plt.get_cmap('RdBu_r').copy();cm.set_bad('#dfdfdf');im=ax.imshow(matrix,cmap=cm,vmin=-lim,vmax=lim,aspect='auto')
            for iy,l in enumerate(ligs):
                for ix,f in enumerate(families):
                    r=b[(b.ligand==l)&(b.celltype==f)].iloc[0]
                    if r.eligibility!='ELIGIBLE':ax.text(ix,iy,'NE',ha='center',va='center',fontsize=6.7,color='#555555');continue
                    n=r.n_delta_negative if m==UP else r.n_delta_positive;c='white' if abs(r.delta_mean)>.65*lim else '#111111'
                    ax.text(ix,iy,str(int(n)),ha='center',va='center',fontsize=7.2,color=c,fontweight='bold' if n==3 else 'normal')
                    if r.FDR_atlas<.05:ax.add_patch(Rectangle((ix-.46,iy-.44),.92,.88,fill=False,ec='black',lw=.65))
            ax.set_xticks(range(4),[flabel[f] for f in families],rotation=35,ha='left');ax.xaxis.tick_top();ax.set_yticks(range(len(ligs)),ligs if m==UP else []);ax.tick_params(length=0);ax.set_title(mlabel[m],pad=47,fontweight='bold')
            cb=fig.colorbar(im,ax=ax,orientation='horizontal',fraction=.03,pad=.04);cb.ax.tick_params(labelsize=7);cb.set_label('Mean paired raw-rank Δ',fontsize=7.5)
        fig.suptitle(f'All-ligand response atlas | {page}/2',x=.13,ha='left',fontsize=15,y=.977)
        fig.text(.13,.943,'All 79 released ligands; 4 epithelial families; 3 observed donors per eligible comparison.',fontsize=9)
        fig.text(.13,.06,'Colour: mean treatment − matched-control score, with a separate scale for each programme.\nCell numbers: donors with UP Δ < 0, depleted Δ > 0, or IFN Δ > 0 (0–3).\nBlack outline: enrichment FDR < 0.05 across all 948 planned primary tests.\nNE: fewer than 3 matched donors; no zero imputation. Ligands differ in tested dose.\nRaw score magnitudes across programmes are not calibrated biological effect sizes.',fontsize=8.2,linespacing=1.6)
        writepage(fig,pdf,name,page)

name='Figure_Response_Patterns';s.to_csv(P/(name+'.tsv'),sep='\t',index=False)
patterns=['BOTH','UP_ONLY','DEPLETED_ONLY','OTHER_OBSERVED','NE'];pcol=['#30867B','#4C84B8','#C99A41','#EEEEEE','#B0B0B0'];ptxt=['Both directions in all 3','UP only in all 3','Depleted only in all 3','Other observed directions','NE (<3 pairs)']
with PdfPages(F/(name+'.pdf')) as pdf:
    for page,start in enumerate([0,40],1):
        ligs=cfg['ligands'][start:start+40];fig=plt.figure(figsize=(8.27,11.69));ax=fig.add_axes([.15,.12,.49,.76]);side=fig.add_axes([.76,.54,.20,.30]);b=s[s.ligand.isin(ligs)];arr=b.pivot(index='ligand',columns='celltype',values='pattern').reindex(index=ligs,columns=families).map(patterns.index).to_numpy()
        ax.imshow(arr,cmap=ListedColormap(pcol),vmin=-.5,vmax=4.5,aspect='auto');ax.set_xticks(range(4),[flabel[f] for f in families],rotation=20,ha='left');ax.xaxis.tick_top();ax.set_yticks(range(len(ligs)),ligs);ax.tick_params(length=0)
        for iy,l in enumerate(ligs):
            for ix,f in enumerate(families):
                r=b[(b.ligand==l)&(b.celltype==f)].iloc[0];ax.text(ix,iy,'NE' if r.pattern=='NE' else str(r.n_same_donor_both),ha='center',va='center',fontsize=8,color='white' if r.pattern=='BOTH' else '#333333')
        count=s[s.condition_class=='ligand'].pattern.value_counts();side.barh(range(5),[count.get(p,0) for p in patterns],color=pcol,ec='#aaaaaa',lw=.5);side.set_yticks(range(5),['Both','UP only','Depleted only','Other','NE'],fontsize=8);side.invert_yaxis();side.set_xlabel('Ligand–family combinations');side.set_xlim(0,max(count)*1.3)
        for i,p in enumerate(patterns):side.text(count.get(p,0)+3,i,str(count.get(p,0)),va='center')
        fig.legend(handles=[Patch(facecolor=c,label=t,edgecolor='#999999',lw=.5) for c,t in zip(pcol,ptxt)],loc='upper left',bbox_to_anchor=(.68,.46),frameon=False,fontsize=8,labelspacing=1.5)
        fig.suptitle(f'Dual-readout donor patterns | {page}/2',x=.15,ha='left',y=.97,fontsize=15)
        fig.text(.72,.30,'BOTH requires, in each\nof the same 3 donors:\nUP Δ < 0 and depleted Δ > 0.\n\nNumbers inside cells:\ndonors satisfying both\ndirections simultaneously.\n\nDirection categories are\ndescriptive, not P values.',fontsize=8.2,linespacing=1.5,va='top')
        fig.text(.15,.055,'278 eligible / 316 planned comparisons. Shared donors and controls; no independent-grid inference.\n“Other” is not “no response”. Missing comparisons cannot establish family specificity.',fontsize=8.2,linespacing=1.5)
        writepage(fig,pdf,name,page)

name='Figure_Reference_Benchmark';bench.to_csv(P/(name+'_correlations.tsv'),sep='\t',index=False);bp.to_csv(P/(name+'_all_pairs.tsv'),sep='\t',index=False)
refs=['delta_IFN','delta_IBD1','delta_IBD2'];refnames=['IFN','IBD1 membership','IBD2 membership'];diseases=['delta_up','delta_depleted']
with PdfPages(F/(name+'.pdf')) as pdf:
    fig,axs=plt.subplots(2,2,figsize=(8.27,8.4));fig.subplots_adjust(left=.17,right=.97,bottom=.24,top=.87,hspace=.55,wspace=.43)
    for ax,f in zip(axs.flat,families):
        b=bench[(bench.celltype==f)&(bench.stratum=='MEAN_OF_3')];a=np.array([[b[(b.disease_readout==di)&(b.reference==re)].iloc[0].spearman_rho for re in refs] for di in diseases]);im=ax.imshow(a,vmin=-1,vmax=1,cmap='RdBu_r',aspect='auto')
        for iy,di in enumerate(diseases):
            for ix,re in enumerate(refs):
                r=b[(b.disease_readout==di)&(b.reference==re)].iloc[0];ax.text(ix,iy,f'ρ={r.spearman_rho:.2f}\nn={r.n_common_ligands}',ha='center',va='center',color='white' if abs(r.spearman_rho)>.65 else 'black',fontsize=9)
        ax.set_title(flabel[f],loc='left',fontweight='bold');ax.set_xticks(range(3),refnames,rotation=20,ha='right');ax.set_yticks(range(2),['Disease-up','Disease-depleted']);ax.tick_params(length=0)
    cax=fig.add_axes([.35,.14,.4,.018]);fig.colorbar(im,cax=cax,orientation='horizontal',label='Spearman correlation of ligand mean paired changes')
    fig.suptitle('Reference programme comparison',x=.17,ha='left',fontsize=15,y=.97)
    fig.text(.17,.92,'Within-family descriptive correlations; each ligand contributes a 3-donor mean.',fontsize=9)
    fig.text(.17,.035,'Shared donors and controls: no independence-based correlation P values.\nIBD1/2 use published membership sets with the same rank rule, not publisher activity scores.',fontsize=8.2,linespacing=1.5)
    writepage(fig,pdf,name,1)
    for page,di in enumerate(diseases,2):
        fig,axs=plt.subplots(4,3,figsize=(8.27,11.69));fig.subplots_adjust(left=.12,right=.97,bottom=.17,top=.91,hspace=.7,wspace=.45)
        for iy,f in enumerate(families):
            for ix,(re,rn) in enumerate(zip(refs,refnames)):
                ax=axs[iy,ix];b=bp[(bp.celltype==f)&(bp.stratum=='MEAN_OF_3')&(bp.disease_readout==di)&(bp.reference==re)];zero(ax,True);ax.scatter(b.disease_delta*1000,b.reference_delta*1000,c=fcolors[f],s=15,alpha=.72,lw=0)
                t=b[b.ligand=='TNFSF12'];ax.scatter(t.disease_delta*1000,t.reference_delta*1000,facecolors='none',edgecolors='black',s=60,lw=1)
                r=bench[(bench.celltype==f)&(bench.stratum=='MEAN_OF_3')&(bench.disease_readout==di)&(bench.reference==re)].iloc[0]
                ax.set_title(f'{flabel[f]} | ρ={r.spearman_rho:.2f}, n={r.n_common_ligands}',fontsize=9);ax.set_xlabel(('UP' if di=='delta_up' else 'Depleted')+' Δ (×10$^{-3}$)');ax.set_ylabel(rn+' Δ (×10$^{-3}$)');ax.ticklabel_format(style='plain',axis='both',useOffset=False)
        fig.suptitle(('Disease-up' if di=='delta_up' else 'Disease-depleted')+' vs reference programmes',x=.12,ha='left',fontsize=14,y=.97)
        fig.text(.12,.045,'Each point: one ligand mean across 3 paired donors; black ring: TNFSF12.\nNo fitted efficacy ranking. All raw score axes retain their own programme scales.\nDonor-stratified coefficients and every discordant comparison are provided in source tables.',fontsize=8.2,linespacing=1.6)
        writepage(fig,pdf,name,page)

name='Supplementary_Selected_Examples';keys=ex[['ligand','celltype']].drop_duplicates().sort_values(['celltype','ligand']);source=d.merge(keys,validate='many_to_one');source.to_csv(P/(name+'_donors.tsv'),sep='\t',index=False);ex.to_csv(P/(name+'_selection_rules.tsv'),sep='\t',index=False)
with PdfPages(F/(name+'.pdf')) as pdf:
    batches=np.array_split(np.arange(len(keys)),int(np.ceil(len(keys)/9)))
    for page,indices in enumerate(batches,1):
        fig,axs=plt.subplots(3,3,figsize=(8.27,10));fig.subplots_adjust(left=.12,right=.97,top=.83,bottom=.14,wspace=.48,hspace=.8)
        for ax,(_,r) in zip(axs.flat,keys.iloc[indices].iterrows()):
            b=d[(d.ligand==r.ligand)&(d.celltype==r.celltype)];zero(ax,True)
            for x in b.itertuples():ax.scatter(-x.delta_up*1000,x.delta_depleted*1000,color=fcolors[x.celltype],marker=markers[x.donor],s=45,ec='white',lw=.5)
            counts=s[(s.ligand==r.ligand)&(s.celltype==r.celltype)].iloc[0]
            ax.set_title(f'{r.ligand} | {flabel[r.celltype]}\nUP− {counts.n_up_negative}/3 | Dep+ {counts.n_depleted_positive}/3 | IFN+ {counts.n_IFN_positive}/3',loc='left',fontsize=8);ax.set_xlabel('−Δ UP (×10$^{-3}$)');ax.set_ylabel('Δ depleted (×10$^{-3}$)')
        for ax in list(axs.flat)[len(indices):]:ax.axis('off')
        fig.suptitle(f'Prespecified display examples | {page}',x=.12,ha='left',fontsize=14,y=.96)
        fig.text(.12,.91,'Selection: alphabetical first within each family and observed pattern / IFN-discordance class.',fontsize=8.2)
        fig.legend(handles=donors_legend(),loc='lower center',bbox_to_anchor=(.5,.07),ncol=3,frameon=False)
        fig.text(.12,.025,'Every example category and all eligible comparisons are preserved in tables.\nThese examples are not a best-ligand selection; shared donors do not add independent replicates.',fontsize=8.2,linespacing=1.5)
        writepage(fig,pdf,name,page)

name='Supplementary_Residual_and_AssayReference';ss=g[(g.module==RES)|(g.condition_class=='assay_reference')];ss.to_csv(P/(name+'.tsv'),sep='\t',index=False)
with PdfPages(F/(name+'.pdf')) as pdf:
    fig,axs=plt.subplots(1,2,figsize=(8.27,4.8));fig.subplots_adjust(left=.16,right=.96,bottom=.27,top=.8,wspace=.46)
    for ax,l in zip(axs,['TNFSF12','No_Cytomix']):
        b=g[(g.ligand==l)&(g.module==RES)];zero(ax)
        for ix,f in enumerate(families):
            r=b[b.celltype==f].iloc[0];ax.scatter(ix,r.NES,c=fcolors[f] if r.FDR_atlas<.05 else 'white',edgecolors=fcolors[f],s=65,lw=1.3);ax.annotate(f'q={r.FDR_atlas:.2g}',(ix,r.NES),xytext=(0,9),textcoords='offset points',ha='center',fontsize=8)
        ax.set_xticks(range(4),[flabel[f] for f in families],rotation=25,ha='right');ax.set_xlim(-.45,3.45);ax.margins(y=.3);ax.set_ylabel('Residual programme NES');ax.set_title(l.replace('_',' '),loc='left');
    fig.suptitle('IFN-excluded residual sensitivity',x=.16,ha='left',y=.95,fontsize=14)
    fig.text(.16,.055,'TNFSF12: separate 316-test all-ligand residual family. No Cytomix: separate 16-test reference family.\nFilled points: FDR < 0.05. All-ligand residual results and coverage are in the accompanying source table.',fontsize=8,linespacing=1.6)
    writepage(fig,pdf,name,1)
    fig,axs=plt.subplots(1,3,figsize=(8.27,4.6));fig.subplots_adjust(left=.12,right=.97,bottom=.28,top=.78,wspace=.45)
    for ax,col,title in zip(axs,['delta_up','delta_depleted','delta_IFN'],['Disease-up','Disease-depleted','IFN reference']):donor_panel(ax,d[d.ligand=='No_Cytomix'],col);ax.set_title(title,loc='left');ax.set_ylabel('Paired Δ (×10$^{-3}$)')
    fig.suptitle('No-cytomix assay reference',x=.12,ha='left',y=.95,fontsize=14);fig.legend(handles=donors_legend(),loc='lower center',bbox_to_anchor=(.5,.07),ncol=3,frameon=False);fig.text(.12,.025,'Reference condition is not counted among the 79 ligands; all four programme results remain in source tables.',fontsize=8)
    writepage(fig,pdf,name,2)
print('PLOTS_COMPLETE',len(list(F.glob('*.pdf'))),flush=True)
