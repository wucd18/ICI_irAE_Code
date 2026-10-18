"""Pure presentation of frozen result tables. No fitting, testing or score creation."""
from pathlib import Path
import json, sys, os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle, FancyArrowPatch
from matplotlib.colors import ListedColormap, Normalize
from theme_01 import apply_style, signed_cmap, load_style

R=Path(os.environ.get('ICI_PRODUCTION_ROOT',Path(__file__).resolve().parents[1]))
STYLE=load_style(); FONT=apply_style(); P=STYLE['palette']; FC=STYLE['family_colors']; DM=STYLE['donor_markers']
UP='ICI_COLITIS_EPITHELIAL_UP_STRICT'; DEP='ICI_COLITIS_EPITHELIAL_LOSS_STRICT'; IFN='CURATED_IFN_VISIBILITY'; RES='ICI_COLITIS_EPITHELIAL_UP_NON_IFN'; IBD1='IBD1_MEMBERSHIP_RANK'; IBD2='IBD2_MEMBERSHIP_RANK'
MODS=[UP,DEP,IFN,RES,IBD1,IBD2]; MN={UP:'Disease-up',DEP:'Depleted',IFN:'IFN',RES:'Residual',IBD1:'IBD1',IBD2:'IBD2'}
FAM=['Stem','Transit_Amplifying','Colonocyte','Goblet']; FN={f:f.replace('Transit_Amplifying','TA') for f in FAM}
PF=['Absorptive epithelial','Immature epithelial','Mature absorptive epithelial','Secretory epithelial']; PN={p:p.replace(' epithelial','').replace('Mature absorptive','Mature abs.') for p in PF};PN['whole_spot_sample']='Whole tissue'
METHODS=['raw','bg_relative','nonrank']; METH={'raw':'Raw rank','bg_relative':'Background relative','nonrank':'Non-rank'}
UNITS={'raw':'Rank change','bg_relative':'Background-relative rank change','nonrank':'Mean gene log2 change'}
CLASSES=['BOTH','UP_ONLY','DEPLETED_ONLY','OTHER_OBSERVED','NE']; CN={'BOTH':'Both','UP_ONLY':'UP only','DEPLETED_ONLY':'Depleted only','OTHER_OBSERVED':'Other','NE':'NE'}
CM=signed_cmap(); READS={}; PANELS=[]
def read(name,old=False):
    p=Path(os.environ.get('ICI_DISPLAY_INPUT_ROOT',R.parent/'display_revision_v2'))/'source_data'/('old' if old else 'new')/name
    if p not in READS:
        d=pd.read_csv(p,sep=',' if p.suffix=='.csv' else '\t',encoding='utf-8-sig')
        d=d.rename(columns={c:'upstream_'+c for c in ['source_file','source_row'] if c in d})
        d.insert(0,'source_row',np.arange(2,len(d)+2));d.insert(0,'source_file',p.relative_to(Path(os.environ['ICI_DISPLAY_INPUT_ROOT'])).as_posix());READS[p]=d
    return READS[p].copy()
D=read('three_method_donor_effects.tsv'); C=read('classification_stability.tsv'); G=read('original_GSEA_results.tsv'); SCO=read('Source_055_patient_level_frozen_module_scores.csv',True)
def title(fig,letter,text):
    fig.text(.015,.976,letter,fontsize=8,fontweight='bold',va='top')
    fig.text(.085,.976,text,fontsize=8.5,fontweight='bold',va='top')
def fignew(w,h,letter,text):
    f=plt.figure(figsize=(w/25.4,h/25.4));title(f,letter,text);return f
def save(fig,pid,data,desc):
    vector_only(fig)
    for e in ['pdf','svg']:fig.savefig(R/'panels'/f'{pid}.{e}')
    fig.savefig(R/'qa'/f'{pid}.png',dpi=180)
    plt.close(fig)
    data.to_csv(R/'panel_data'/f'{pid}.tsv',sep='\t',index=False,na_rep='NE')
    PANELS.append({'panel':pid,'rows':len(data),'display':desc,'data':f'panel_data/{pid}.tsv'})
def axis(fig,rect):return fig.add_axes(rect)
def vector_only(fig):
    for artist in fig.findobj():
        if hasattr(artist,'set_rasterized'):artist.set_rasterized(False)
    # Adjacent vector colourbar cells overlap by a hairline to avoid PDF-viewer seams.
    for ax in fig.axes:
        cb=getattr(ax,'_colorbar',None)
        if cb is not None and cb.solids is not None:
            cb.solids.set_edgecolor('face');cb.solids.set_linewidth(.3)

def vector_matrix(ax,a,cmap,vmin,vmax):
    a=np.asarray(a,float)
    mesh=ax.pcolormesh(np.arange(a.shape[1]+1)-.5,np.arange(a.shape[0]+1)-.5,a,cmap=cmap,vmin=vmin,vmax=vmax,shading='flat',rasterized=False)
    ax.set_xlim(-.5,a.shape[1]-.5);ax.set_ylim(a.shape[0]-.5,-.5)
    return mesh
def zero(ax):ax.axvline(0,color=P['muted'],lw=.5,ls=(0,(3,3)),zorder=0)
def leg_donor(fig,y=.02):
    hs=[Line2D([],[],marker=m,ls='',color=P['ink'],markersize=3.8,label=d) for d,m in DM.items()]
    fig.legend(handles=hs,loc='lower center',ncol=3,bbox_to_anchor=(.5,y),handletextpad=.35,columnspacing=1)
def patient_points(ax,d,groups,groupcol='condition',val='disease_oriented_score',labels=None):
    for i,g in enumerate(groups):
        dd=d[d[groupcol]==g].sort_values('patient_id' if 'patient_id' in d else 'patient')
        color=P['teal'] if g=='UC' else (P['red'] if g in ['CPI_colitis','irColitis','irAE_colitis','ICI','ICI colitis'] else P['muted'])
        x=np.linspace(-.13,.13,len(dd))+i
        ax.scatter(x,dd[val],s=13,c=color,edgecolor='white',lw=.25,zorder=3)
        if len(dd):ax.plot([i-.22,i+.22],[dd[val].median()]*2,color=P['ink'],lw=1)
    ax.set_xticks(range(len(groups)),labels or groups);ax.set_xlim(-.55,len(groups)-.45)
    ax.tick_params(axis='x',length=0);ax.ticklabel_format(axis='y',style='plain',useOffset=False)
def heat(ax,mat,rows,cols,limit=1,numbers=True,fmt='.2f',cmap=None):
    a=np.asarray(mat,float);im=vector_matrix(ax,a,cmap or CM,-limit,limit)
    ax.set_xticks(range(len(cols)),cols);ax.set_yticks(range(len(rows)),rows)
    ax.tick_params(length=0,pad=3)
    for sp in ax.spines.values():sp.set_visible(False)
    ax.set_xticks(np.arange(-.5,len(cols),1),minor=True);ax.set_yticks(np.arange(-.5,len(rows),1),minor=True);ax.grid(which='minor',color='white',lw=.7);ax.tick_params(which='minor',length=0)
    for i in range(a.shape[0]):
        for j in range(a.shape[1]):
            if np.isnan(a[i,j]):ax.text(j,i,'NE',ha='center',va='center',fontsize=7.5)
            elif numbers:ax.text(j,i,format(a[i,j],fmt),ha='center',va='center',fontsize=7.5,color='white' if abs(a[i,j])>.75*limit else P['ink'])
    return im
def donor_axis(ax,d,groups,mod,method='nonrank',labels=None,show_y=True,limit=None,show_mean=True):
    vals=[]
    for i,(lig,fam) in enumerate(groups):
        s=d[(d.ligand==lig)&(d.celltype==fam)&(d.module==mod)]
        if s.empty:ax.text(.5,i,'NE',transform=ax.get_yaxis_transform(),ha='center',va='center',fontsize=7.5)
        else:
            for k,(don,mark) in enumerate(DM.items()):
                r=s[s.donor==don]
                if r.empty:continue
                v=r[method].iloc[0];vals.append(v)
                ax.scatter(v,i+(k-1)*.15,marker=mark,s=15,c=FC[fam],edgecolor='white',lw=.2,zorder=3)
            if show_mean:
                mean=s[method].mean();ax.plot([mean,mean],[i-.28,i+.28],lw=1.1,color=P['ink'])
    zero(ax);ax.set_ylim(len(groups)-.55,-.55)
    ax.set_yticks(range(len(groups)),labels if show_y else ['']*len(groups));ax.tick_params(axis='y',length=0)
    ax.set_xlabel(UNITS[method]);ax.locator_params(axis='x',nbins=3)
    if limit is not None:ax.set_xlim(-limit,limit)
    elif vals:
        lo=min(0,min(vals));hi=max(0,max(vals));pad=(hi-lo)*.1 or .001;ax.set_xlim(lo-pad,hi+pad)

def f1():
    f=fignew(170,42,'a','Fixed patient programmes linked to ligand responses')
    ax=axis(f,[.04,.03,.93,.76]);ax.axis('off')
    labels=['Patient discovery','External assessment','Organoid perturbations','Response context']
    sub=['GSE253720\nPatient pseudobulks','GSE206300 / GSE210037\nOxford / Xenium','GSE313368\nMatched donor controls','Maturation states\nIFN and IBD references']
    for i in range(4):
        x=.12+i*.25
        ax.text(x,.68,labels[i],ha='center',fontsize=8,fontweight='bold');ax.text(x,.29,sub[i],ha='center',va='center',fontsize=7.5)
        if i<3:ax.annotate('',(x+.18,.53),(x+.09,.53),arrowprops={'arrowstyle':'->','color':P['muted'],'lw':.7})
    ax.set_xlim(0,1);ax.set_ylim(0,1)
    t=read('Table_1_full_cohort_details.csv',True);save(f,'F1a',t[t.context.str.contains('colon|organoid',case=False)],'Design arrows link data stages, not causal effects.')
    f=fignew(96,75,'b','Biological replicate coverage')
    ax=axis(f,[.40,.27,.54,.58]);rows=[]
    for ds,ct in [('GSE253720','Epithelial')]+[('GSE206300',p) for p in PF]+[('GSE210037','whole_spot_sample')]:
        ss=SCO[(SCO.dataset==ds)&(SCO.celltype==ct)&(SCO.module==UP)]
        for case in [False,True]:
            z=ss[ss.condition.isin(['CPI_colitis','irColitis','irAE_colitis'])==case]
            rows.append(dict(dataset=ds,celltype=ct,case=case,n=z.patient_id.nunique(),source_file=z.source_file.iloc[0],source_rows=';'.join(z.source_row.astype(str))))
    b=pd.DataFrame(rows)
    for i,(_,z) in enumerate(b.groupby(['dataset','celltype'],sort=False)):
        for _,r in z.iterrows():
            co=P['red'] if r['case'] else P['muted'];y=i+(.12 if r['case'] else -.12)
            ax.plot([0,r.n],[y,y],color=co,lw=.7);ax.scatter(r.n,y,s=13,c=co);ax.text(r.n+.2,y,str(r.n),va='center',fontsize=7.5,color=co)
    ax.set_yticks(range(6),['Discovery epithelium']+[PN[p] for p in PF]+['Whole tissue']);ax.invert_yaxis();ax.set_xlim(0,b.n.max()+1.3);ax.set_xlabel('Patients or released samples');ax.tick_params(axis='y',length=0)
    f.text(.05,.05,'Red: ICI cases   Grey: controls\nExternal families share patients; do not add counts.',fontsize=7.5)
    save(f,'F1b',b,'Patient counts from saved score rows; organoid coverage in F3a.')
    f=fignew(74,75,'c','Fixed programme definitions');ax=axis(f,[.34,.29,.56,.54]);defs=read('Source_021_frozen_module_definitions.csv',True);d=defs[defs.module.isin([UP,DEP,IFN,RES])]
    sizes=d.groupby('module').gene.nunique().reindex([UP,DEP,IFN,RES]);cols=[P['gold'],P['blue'],P['magenta'],P['purple']]
    ax.barh(range(4),sizes,color=cols,height=.55);ax.set_yticks(range(4),[MN[m] for m in sizes.index]);ax.invert_yaxis();ax.set_xlim(0,max(sizes)*1.22);ax.set_xlabel('Fixed members');ax.tick_params(axis='y',length=0)
    for i,v in enumerate(sizes):ax.text(v+4,i,str(v),va='center',fontsize=7.5)
    ov=len(set(d[d.module==UP].gene)&set(d[d.module==IFN].gene));f.text(.08,.06,f'UP / IFN overlap: {ov} genes\nResidual excludes these IFN members.',fontsize=7.5)
    save(f,'F1c',d,'Distinct frozen gene counts; residual is an IFN-exclusion sensitivity.')
    for pid,mod,let,tt in [('F1d',UP,'d','Discovery disease-up'),('F1e',DEP,'e','Discovery depleted: lower in disease')]:
        f=fignew(85,65,let,tt);ax=axis(f,[.21,.20,.73,.60]);d=SCO[(SCO.dataset=='GSE253720')&(SCO.module==mod)]
        patient_points(ax,d,['HC','CPI_colitis'],val='raw_mean_within_sample_gene_rank',labels=['Healthy','ICI colitis']);ax.set_ylabel('Raw depleted rank\n(lower in disease)' if mod==DEP else 'Raw disease-up rank\n(unitless)')
        f.text(.1,.035,'Definition-set display; points are patients.',fontsize=7.5);save(f,pid,d,'All discovery patient scores, no orientation inversion.')

def f2():
    v=read('Source_022_frozen_module_validation.csv',True);v=v[v.contrast.isin(['irColitis_vs_ICI_control_PD1_only','irAE_colitis_vs_healthy_colon'])&v.module.isin([UP,DEP,IFN])].copy();v['display_NES']=np.where(v.module==DEP,-v.NES,v.NES)
    f=fignew(170,55,'a','External programme transfer');ax=axis(f,[.29,.20,.54,.60]);cts=PF+['whole_spot_sample'];mods=[IFN,UP,DEP]
    norm=Normalize(-3,3)
    for _,r in v.iterrows():
        x=mods.index(r.module);y=cts.index(r.celltype);ax.scatter(x,y,s=55,c=[CM(norm(r.display_NES))],edgecolor=P['ink'],lw=.5)
        if r.padj<.05:ax.text(x,y,'*',ha='center',va='center',fontsize=8)
    ax.set_xticks(range(3),[MN[m] for m in mods]);ax.set_yticks(range(5),[PN[c] for c in cts]);ax.set_xlim(-.5,2.5);ax.set_ylim(4.6,-.6);ax.tick_params(length=0);[s.set_visible(False) for s in ax.spines.values()]
    ca=axis(f,[.9,.29,.025,.32]);cb=f.colorbar(plt.cm.ScalarMappable(norm=norm,cmap=CM),cax=ca,ticks=[-3,0,3]);ca.set_title('NES',fontsize=7.5)
    f.text(.09,.035,'Disease-oriented NES: depleted = -NES; * source FDR < 0.05.',fontsize=7.5);save(f,'F2a',v,'Saved NES; depleted sign reversed only for patient disease orientation.')
    for pid,mod,let,tt in [('F2b',UP,'b','GSE206300 disease-up'),('F2c',DEP,'c','GSE206300 depleted')]:
        f=fignew(85,73,let,tt);d=SCO[(SCO.dataset=='GSE206300')&(SCO.module==mod)]
        ax=axis(f,[.2,.30,.75,.51])
        for i,ct in enumerate(PF):
            z=d[d.celltype==ct]
            for j,g in enumerate(['On ICI therapy','irColitis']):
                zz=z[z.condition==g].sort_values('patient_id');xx=i+(j-.5)*.3+np.linspace(-.06,.06,len(zz));co=P['muted'] if j==0 else P['red']
                ax.scatter(xx,zz.disease_oriented_score,s=11,c=co,edgecolor='white',lw=.2);ax.plot([i+(j-.5)*.3-.12,i+(j-.5)*.3+.12],[zz.disease_oriented_score.median()]*2,c=P['ink'],lw=.8)
        ax.set_xticks(range(4),['Absorptive','Immature','Mature abs.','Secretory'],rotation=30,ha='right');ax.set_ylabel('1 - raw depleted rank\n(higher = disease direction)' if mod==DEP else 'Raw disease-up rank\n(unitless)');ax.tick_params(axis='x',length=0)
        f.text(.08,.065,'Grey: on-ICI controls; red: colitis\n'+('UP = raw rank.' if mod==UP else 'Depleted = 1 - raw rank.'),fontsize=7.5);save(f,pid,d,'All primary PD1-only patients in four released families; patient families not organoid families.')
    f=fignew(85,75,'d','GSE210037 whole-tissue pseudobulk');d=SCO[(SCO.dataset=='GSE210037')&SCO.module.isin([UP,DEP])]
    for i,m in enumerate([UP,DEP]):
        ax=axis(f,[.19+i*.44,.27,.34,.51]);patient_points(ax,d[d.module==m],['healthy_colon','irAE_colitis'],labels=['Healthy','Colitis']);ax.set_title(MN[m],fontsize=8);ax.tick_params(axis='x',labelrotation=30)
        if i==0:ax.set_ylabel('Disease-oriented rank\n(unitless)')
    f.text(.07,.055,'Released samples; no spatial localisation.\nDepleted = 1 - raw rank.',fontsize=7.5);save(f,'F2d',d,'Whole-tissue released samples, not spots or epithelial-only measurements.')
    ox=read('Source_065_new_spatial_patient_scores.csv',True);q=read('Source_066_new_spatial_frozen_transfer.csv',True)
    f=fignew(85,75,'e','Oxford donor-level assessment')
    for i,m in enumerate([UP,DEP,IFN]):
        ax=axis(f,[.17+i*.275,.31,.18,.45]);patient_points(ax,ox[ox.module==m],['Healthy','ICI colitis','UC'],val='oriented_score',labels=['H','ICI','UC']);ax.set_title(MN[m],fontsize=8);ax.locator_params(axis='y',nbins=3)
        ax.text(.5,-.32,'q = '+format(q[q.module==m].FDR_new_four_module_family.iloc[0],'.3f'),transform=ax.transAxes,ha='center',fontsize=7.5)
    f.text(.06,.055,'Unitless ranks; depleted = 1 - raw.\nH: healthy; UC descriptive. Original four-programme FDR.',fontsize=7.5);save(f,'F2e',ox.merge(q[['module','FDR_new_four_module_family']],on='module'),'All Oxford donors; repeated sections aggregated in saved table.')

if __name__=='__main__':
    f1();f2()
    pd.DataFrame(PANELS).to_csv(R/'qa/patient_panel_display.tsv',sep='\t',index=False)
    (R/'qa/font_report.json').write_text(json.dumps(FONT,indent=2))
