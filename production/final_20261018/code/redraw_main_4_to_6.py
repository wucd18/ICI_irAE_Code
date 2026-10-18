"""Redraw integrated Figures 3-6 from saved panel/source tables only.

No DE/GSEA, scoring, eligibility, gene selection or significance is recomputed.
Each chart is a separate Matplotlib figure; assemble_main.py combines PDF/SVG vectors.
Usage: python redraw_display_revision.py --input REVISION_DIR --output REVISION_DIR
"""
from __future__ import annotations
import argparse, json
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.lines import Line2D
from matplotlib.ticker import PercentFormatter, MaxNLocator


UP='ICI_COLITIS_EPITHELIAL_UP_STRICT'; DEP='ICI_COLITIS_EPITHELIAL_LOSS_STRICT'
IFN='CURATED_IFN_VISIBILITY'; IBD1='IBD1_MEMBERSHIP_RANK'; IBD2='IBD2_MEMBERSHIP_RANK'
FAMILIES=['Stem','Transit_Amplifying','Colonocyte','Goblet']
FN={'Stem':'Stem','Transit_Amplifying':'TA','Colonocyte':'Colonocyte','Goblet':'Goblet'}
MN={UP:'Disease-up',DEP:'Disease-depleted',IFN:'IFN',IBD1:'IBD1',IBD2:'IBD2'}
LOG=[]
LAYOUT={}
TEXT_QA=[]

def run(inp:Path,out:Path)->None:
    LOG.clear()
    if not (inp/'panel_data/F3e.tsv').is_file():
        raise FileNotFoundError('Input must contain panel_data/F3e.tsv and source_data/new/.')
    out.mkdir(parents=True,exist_ok=True)
    for k in ['panels','panel_data','previews']: (out/k).mkdir(exist_ok=True)
    s=json.loads((inp/'config/style_01.json').read_text(encoding='utf-8'))
    P=s['palette']; FC=s['family_colors']; DM=s['donor_markers']
    names={f.name for f in font_manager.fontManager.ttflist}
    font=next((x for x in ['Arial','Arimo','Liberation Sans','DejaVu Sans'] if x in names),'DejaVu Sans')
    plt.rcParams.update({'font.family':font,'font.size':7.5,'axes.titlesize':8,'axes.labelsize':7.5,
      'xtick.labelsize':7,'ytick.labelsize':7,'axes.linewidth':.55,'lines.linewidth':.85,
      'axes.spines.top':False,'axes.spines.right':False,'axes.edgecolor':P['muted'],
      'text.color':P['ink'],'axes.labelcolor':P['ink'],'xtick.color':P['ink'],'ytick.color':P['ink'],
      'figure.facecolor':'white','axes.facecolor':'white','pdf.fonttype':42,'ps.fonttype':42,
      'svg.fonttype':'none','savefig.facecolor':'white'})
    def panel(n): return pd.read_csv(inp/f'archive/v1_panel_data/{n}.tsv',sep='\t',na_values=['NE'])
    D=pd.read_csv(inp/'source_data/new/three_method_donor_effects.tsv',sep='\t')
    D.insert(0,'source_row',np.arange(2,len(D)+2));D.insert(0,'source_file','source_data/new/three_method_donor_effects.tsv')
    COR=pd.read_csv(inp/'source_data/new/reference_correlations_by_method.tsv',sep='\t')
    def new(w,h,letter,title,rect):
        f=plt.figure(figsize=(w/25.4,h/25.4));ax=f.add_axes(rect)
        f.text(.025,.975,letter,ha='left',va='top',fontsize=8,fontweight='bold')
        f.text(.09,.975,title,ha='left',va='top',fontsize=8.5,fontweight='bold')
        ax.tick_params(length=2.5,width=.5);return f,ax
    def save(f,pid,data,note):
        for artist in f.findobj():
            if hasattr(artist,'set_rasterized'):artist.set_rasterized(False)
        f.canvas.draw()
        from matplotlib.text import Text
        ren=f.canvas.get_renderer(); fw,fh=f.canvas.get_width_height()
        for t in f.findobj(match=Text):
            if not t.get_visible() or not t.get_text():continue
            bb=t.get_window_extent(ren)
            outside=bb.x0 < -1 or bb.y0 < -1 or bb.x1 > fw+1 or bb.y1 > fh+1
            TEXT_QA.append(dict(panel=pid,text=t.get_text(),font_pt=t.get_fontsize(),outside_canvas=outside))
        for ext in ['pdf','svg']:f.savefig(out/'panels'/f'{pid}.{ext}')
        f.savefig(out/'qa'/f'{pid}.png',dpi=180)
        plt.close(f);data.to_csv(out/'panel_data'/f'{pid}.tsv',sep='\t',index=False)
        LOG.append({'panel':pid,'rows':len(data),'note':note})
    def zero(ax): ax.axvline(0,color=P['muted'],lw=.55,ls=(0,(3,3)),zorder=0)
    def family_legend(f,y=.035,ncol=4):
        hs=[Line2D([],[],ls='',marker='o',markersize=3.5,color=FC[x],label=FN[x]) for x in FAMILIES]
        f.legend(handles=hs,loc='lower center',bbox_to_anchor=(.52,y),ncol=ncol,frameon=False,
                 fontsize=7,handletextpad=.3,columnspacing=.8)
    def donor_legend(f,y=.015):
        hs=[Line2D([],[],ls='',marker=ma,markersize=3.7,color=P['ink'],label=dn) for dn,ma in DM.items()]
        f.legend(handles=hs,loc='lower center',bbox_to_anchor=(.5,y),ncol=3,frameon=False,
                 fontsize=7,handletextpad=.3,columnspacing=.75)
    def assemble(num,w,h,slots):
        LAYOUT[str(num)]={'size':[w,h],'panels':[[path.stem,x,y,ww,hh] for path,x,y,ww,hh in slots]}
    # FIGURE 4: two response blocks, one expression scale, three readouts; no GSEA overlay.
    groups_top=[(l,fam) for l in ['IFNG','TNFA'] for fam in FAMILIES]
    groups_low=[('TNFSF12',fam) for fam in FAMILIES]
    slots=[];mods=[UP,DEP,IFN]; widths=[70,50,50];xs=[0,70,120]
    titles={UP:'Disease-up',DEP:'Disease-depleted',IFN:'IFN reference'}
    limits={}
    for mod in mods:
        v=D[(D.ligand.isin(['IFNG','TNFA','TNFSF12']))&(D.module==mod)].nonrank
        lo=min(0,v.min());hi=max(0,v.max());pad=(hi-lo)*.09;limits[mod]=(lo-pad,hi+pad)
    for row,groups in enumerate([groups_top,groups_low]):
        for col,mod in enumerate(mods):
            let='abcdef'[row*3+col];pid=f'F4{let}';w=widths[col]
            f,ax=new(w,85,let,titles[mod],[.39 if col==0 else .10,.25,.55 if col==0 else .84,.57])
            data=D[(D.ligand.isin(['IFNG','TNFA'] if row==0 else ['TNFSF12']))&(D.module==mod)]
            for i,(lig,fam) in enumerate(groups):
                q=data[(data.ligand==lig)&(data.celltype==fam)]
                if q.empty:ax.text(.55,i,'NE',transform=ax.get_yaxis_transform(),ha='center',va='center',fontsize=7.5,color=P['muted'])
                for k,(dn,marker) in enumerate(DM.items()):
                    rr=q[q.donor==dn]
                    if len(rr):ax.scatter(rr.nonrank,i+(k-1)*.15,s=20,color=FC[fam],marker=marker,edgecolors='white',lw=.25,zorder=3)
            zero(ax);ax.set_xlim(limits[mod]);ax.set_ylim(len(groups)-.65,-.65)
            labels=[(f'{lig} / {FN[fam]}' if row==0 else FN[fam]) for lig,fam in groups]
            ax.set_yticks(range(len(groups)),labels if col==0 else ['']*len(groups));ax.tick_params(axis='y',length=0)
            ax.xaxis.set_major_locator(MaxNLocator(3));ax.set_xlabel('Mean gene log2 change')
            if col==0:f.text(.08,.875,'IFNG and TNFA' if row==0 else 'TNFSF12',fontsize=7.5,fontweight='bold')
            if row==0:
                lab='Decrease opposes disease' if mod==UP else 'Increase opposes disease' if mod==DEP else 'Direction shown; not a benefit score'
                f.text(.08,.09,lab.replace('Direction shown; not a benefit score','IFN reference'),fontsize=7)
            if row==1 and col==0:
                v=D[(D.ligand=='TNFSF12')&(D.celltype=='Transit_Amplifying')&(D.module==UP)&(D.donor=='DNW9')].nonrank.iloc[0]
                f.text(.08,.075,f'TA DNW9: {v:+.6f} (near zero)',fontsize=7)
            if row==1 and col==1:donor_legend(f,.035)
            if row==1 and col==2:f.text(.08,.075,'All comparisons: matched\ncytomix control',fontsize=7,linespacing=1.3)
            save(f,pid,data,'Same saved nonrank donor effects; common x scale within programme across rows. NE retained. GSEA and all representations retained in Tables S4 and S5.')
            slots.append((out/'panels'/f'{pid}.pdf',xs[col],row*85,w,85))
    assemble(4,170,170,slots)
    # FIGURE 5: effect -> broad gene response -> composition -> state-weighted response.
    fa=panel('F5a');f,ax=new(65,82,'a','BMP4 programme changes',[.37,.23,.58,.57])
    groups=[(UP,'BMP4'),(UP,'No_Cytomix'),(DEP,'BMP4'),(DEP,'No_Cytomix')]
    for i,(mod,lig) in enumerate(groups):
        q=fa[(fa.module==mod)&(fa.ligand==lig)]
        for k,(dn,marker) in enumerate(DM.items()):
            v=q[q.donor==dn].nonrank.iloc[0];ax.scatter(v,i+(k-1)*.13,color=FC['Colonocyte'],marker=marker,s=21,edgecolor='white',lw=.25)
    zero(ax);ax.axhline(1.5,c=P['ne'],lw=.75);ax.set_ylim(3.7,-.7);ax.set_xlim(-1.18,1.35)
    ax.set_yticks(range(4),['UP / BMP4','UP / No-cytomix','Depleted / BMP4','Depleted /\nNo-cytomix']);ax.tick_params(axis='y',length=0)
    ax.xaxis.set_major_locator(MaxNLocator(4));ax.set_xlabel('Mean gene log2 change');donor_legend(f,.055)
    save(f,'F5a',fa,'Unchanged nonrank effects. A shared scale is used for the two programmes. This is not a potency comparison.')
    fb=panel('F5b');f,ax=new(105,82,'b','Changes extend across programme members',[.25,.23,.60,.58]);rng=np.random.default_rng(1734)
    summary=[]
    for i,(mod,dn) in enumerate([(m,dn) for m in [UP,DEP] for dn in DM]):
        a=fb[(fb.module==mod)&(fb.donor==dn)].sort_values('gene').gene_log_delta.to_numpy()
        ax.scatter(a,i+rng.uniform(-.19,.19,len(a)),s=3.6,c=FC['Colonocyte'],alpha=.40,linewidths=0)
        med=np.median(a);ax.plot([med,med],[i-.30,i+.30],color=P['ink'],lw=1.05)
        fraction=float(np.mean(a<0) if mod==UP else np.mean(a>0))
        ax.text(1.025,i,f'{fraction:.0%} '+('down' if mod==UP else 'up'),transform=ax.get_yaxis_transform(),ha='left',va='center',fontsize=7)
        summary.append({'module':mod,'donor':dn,'n_members':len(a),'fraction_in_expected_perturbation_direction':fraction,'median':med})
    zero(ax);ax.axhline(2.5,color=P['ne'],lw=.8);ax.set_ylim(5.7,-.65)
    ax.set_yticks(range(6),[('UP' if m==UP else 'Depleted')+' / '+dn for m in [UP,DEP] for dn in DM]);ax.tick_params(axis='y',length=0)
    ax.set_xlabel('Member gene log2 normalized-count change');ax.xaxis.set_major_locator(MaxNLocator(5))
    f.text(.1,.09,f'{fb[fb.module==UP].gene.nunique()} UP and {fb[fb.module==DEP].gene.nunique()}/'+str(pd.read_csv(inp/'source_data/old/Source_021_frozen_module_definitions.csv').query('module == @DEP').gene.nunique())+' depleted members per donor.',fontsize=7)
    f.text(.1,.045,'Each small point is a gene; black marks are medians.',fontsize=7)
    save(f,'F5b',fb,'All 1,032 gene-donor values plotted with deterministic vertical jitter only; member fractions are descriptive.')
    pd.DataFrame(summary).to_csv(out/'panel_data/F5b_member_summary.tsv',sep='\t',index=False)
    fc=panel('F5c');f,ax=new(56.67,94,'c','Recovered mature-state fraction',[.29,.31,.63,.48])
    for r in fc.itertuples():ax.plot([0,1],[r.control_proportion,r.treatment_proportion],color=FC['Colonocyte'],marker=DM[r.donor],ms=4,lw=.85)
    ax.set_ylim(0,1);ax.set_xlim(-.15,1.15);ax.set_xticks([0,1],['Cytomix\ncontrol','BMP4']);ax.yaxis.set_major_formatter(PercentFormatter(1,decimals=0));ax.set_ylabel('Recovered mature-state fraction')
    delta=(fc.treatment_proportion-fc.control_proportion)*100
    f.text(.09,.13,f'Increase: {delta.min():.1f}-{delta.max():.1f} percentage points',fontsize=7)
    # Recovered composition is defined in the adjacent legend.
    save(f,'F5c',fc,'Proportions only converted to percentage-axis labels; source values unchanged.')
    for mod,pid,let,txt in [(UP,'F5d','d','UP decrease persists'),(DEP,'F5e','e','Depleted increase attenuates\nbut remains positive')]:
        dt=panel(pid);f,ax=new(56.67,94,let,txt,[.28,.31,.65,.48])
        for r in dt.itertuples():ax.plot([0,1],[r.unstandardized_common_delta,r.standardized_common_delta],color=FC['Colonocyte'],marker=DM[r.donor],ms=4,lw=.85)
        ax.axhline(0,color=P['muted'],lw=.5);ax.set_xlim(-.16,1.16)
        ax.set_ylim((-.28,.015) if mod==UP else (-.012,.235));ax.yaxis.set_major_locator(MaxNLocator(4))
        ax.set_xticks([0,1],['Observed\nmixture','Control-state\nweights']);ax.set_ylabel('Mean-cell ln-normalized change')
        f.text(.09,.13,'Same cells and substates;\nfixed control-state weights.',fontsize=7)
        # Interpretation is stated in the legend and Discussion.
        save(f,pid,dt,'Exactly the same common-support cells; only saved state-weighting summaries compared, in the cell ln scale.')
    assemble(5,170,176,[(out/'panels/F5a.pdf',0,0,65,82),(out/'panels/F5b.pdf',65,0,105,82),
     (out/'panels/F5c.pdf',0,82,56.67,94),(out/'panels/F5d.pdf',56.67,82,56.67,94),(out/'panels/F5e.pdf',113.34,82,56.66,94)])
    # FIGURE 6: overlap bars, 24 correlation dots, same-family scatter examples, removal contrast.
    ova=panel('F6a');ov=ova[ova.scope=='measured'].copy();pairs=[(a,b) for a in [UP,DEP] for b in [IFN,IBD1,IBD2]]
    f,ax=new(82,76,'a','Direct gene overlap',[.36,.24,.55,.57])
    for i,(a,b) in enumerate(pairs):
        q=ov[(ov.programme_a==a)&(ov.programme_b==b)].iloc[0]
        co=P['gold'] if a==UP else P['blue_point'];ax.barh(i,q.n_intersection,color=co,height=.53)
        ax.text(q.n_intersection+.3,i,f'{int(q.n_intersection)}/{int(q.n_a)}',va='center',fontsize=7)
    ax.set_yticks(range(6),[('UP' if a==UP else 'Depleted')+' / '+MN[b] for a,b in pairs]);ax.set_ylim(5.65,-.65);ax.set_xlim(0,23);ax.set_xticks([0,5,10,15,20]);ax.tick_params(axis='y',length=0)
    ax.set_xlabel('Shared measured genes');f.text(.08,.08,'Labels: intersection / measured ICI members.',fontsize=7)
    save(f,'F6a',ova,'Intersection counts and denominators are saved measured memberships; full fixed/assayed coverage stays in source table.')
    co=panel('F6c');f,ax=new(88,76,'b','Response correlations across families',[.36,.25,.60,.56])
    for i,(a,b) in enumerate(pairs):
        for j,fam in enumerate(FAMILIES):
            q=co[(co.programme_a==a)&(co.programme_b==b)&(co.celltype==fam)].iloc[0]
            ax.scatter(q.Spearman,i+(j-1.5)*.145,s=17,color=FC[fam],edgecolors='white',lw=.2)
    zero(ax);ax.set_xlim(-1.04,1.04);ax.set_xticks([-1,-.5,0,.5,1]);ax.set_yticks(range(6),[('UP' if a==UP else 'Depleted')+' / '+MN[b] for a,b in pairs]);ax.set_ylim(5.6,-.65);ax.tick_params(axis='y',length=0)
    ax.set_xlabel('Spearman correlation (non-rank)');family_legend(f,.015,2)
    save(f,'F6b',co,'All 24 saved full-membership nonrank, MEAN_OF_3 correlations shown without invented error bars or P values.')
    fam='Transit_Amplifying';sm=D[(D.condition_class=='ligand')&(D.celltype==fam)].groupby(['ligand','module']).nonrank.mean().unstack('module')
    for a,pid,let in [(UP,'F6c','c'),(DEP,'F6d','d')]:
        zz=sm[[IBD2,a]].dropna().reset_index();rr=co[(co.celltype==fam)&(co.programme_a==a)&(co.programme_b==IBD2)].iloc[0]
        if len(zz)!=rr.n_common_ligands:raise ValueError('Scatter population differs from the saved correlation population.')
        f,ax=new(85,71,let,'TA: '+MN[a]+' and IBD2',[.19,.27,.74,.55])
        ax.scatter(zz[IBD2],zz[a],s=12,color=FC[fam],alpha=.65,edgecolors='white',lw=.15)
        zero(ax);ax.axhline(0,color=P['muted'],lw=.5);ax.set_xlabel('IBD2 non-rank change');ax.set_ylabel(MN[a]+' non-rank change')
        ax.xaxis.set_major_locator(MaxNLocator(4));ax.yaxis.set_major_locator(MaxNLocator(4))
        f.text(.19,.125,f'n = {int(rr.n_common_ligands)} ligands; Spearman = {rr.Spearman:.3f}',fontsize=7.5)
        f.text(.19,.067,'Each point: mean of the same three donors.',fontsize=7)
        zz['saved_Spearman']=rr.Spearman
        zz['source_file']='source_data/new/three_method_donor_effects.tsv'
        zz['source_rows']=zz.ligand.map(lambda lig:';'.join(D[(D.ligand==lig)&(D.celltype==fam)&(D.module.isin([a,IBD2]))].source_row.astype(str)))
        zz['correlation_source_file']=str(rr.source_file);zz['correlation_source_row']=rr.source_row
        save(f,pid,zz,'Same TA population for both illustrative programme comparisons; all families are represented in panel b and full source tables.')
    rem=panel('F6d');f,ax=new(170,57,'e','Depleted-IBD2 correlation after shared-member removal',[.17,.27,.62,.48]);rrs=[]
    for i,fam in enumerate(FAMILIES):
        full=co[(co.celltype==fam)&(co.programme_a==DEP)&(co.programme_b==IBD2)].iloc[0]
        rm=rem[(rem.celltype==fam)&(rem.method=='nonrank')&(rem.programme_a==DEP)&(rem.programme_b==IBD2)].iloc[0]
        ax.plot([full.Spearman,rm.Spearman],[i,i],color=FC[fam],lw=1)
        ax.scatter(full.Spearman,i,s=27,facecolors='white',edgecolors=FC[fam],lw=1,zorder=3)
        ax.scatter(rm.Spearman,i,s=27,color=FC[fam],edgecolors='white',lw=.35,zorder=3)
        ax.text(1.04,i,f'{full.Spearman:.2f} to {rm.Spearman:.2f}',transform=ax.get_yaxis_transform(),fontsize=7.5,va='center')
        rrs.append({'family':fam,'full_Spearman':full.Spearman,'overlap_removed_Spearman':rm.Spearman,'n_full':full.n_common_ligands,'n_removed':rm.n_common_ligands,'full_source_file':full.source_file,'full_source_row':full.source_row,'removed_source_file':rm.source_file,'removed_source_row':rm.source_row})
    ax.set_yticks(range(4),[FN[x] for x in FAMILIES]);ax.set_ylim(3.5,-.5);ax.set_xlim(-.02,1.02);ax.set_xticks([0,.25,.5,.75,1]);ax.tick_params(axis='y',length=0);ax.set_xlabel('Spearman correlation (non-rank)')
    hs=[Line2D([],[],marker='o',ls='',markerfacecolor='white',markeredgecolor=P['ink'],label='Full memberships',markersize=4),Line2D([],[],marker='o',ls='',color=P['ink'],label='Shared genes removed from both sets',markersize=4)]
    f.legend(handles=hs,loc='upper center',bbox_to_anchor=(.5,.86),ncol=2,frameon=False,fontsize=7.5)
    save(f,'F6e',pd.DataFrame(rrs),'Saved correlations, not new tests. Positive-only display axis explicitly runs 0-1; all-programme panel b uses -1 to +1.')
    assemble(6,170,204,[(out/'panels/F6a.pdf',0,0,82,76),(out/'panels/F6b.pdf',82,0,88,76),
      (out/'panels/F6c.pdf',0,76,85,71),(out/'panels/F6d.pdf',85,76,85,71),(out/'panels/F6e.pdf',0,147,170,57)])
    layout=json.loads((out/'config/layout.json').read_text())
    layout.update(LAYOUT);(out/'config/layout.json').write_text(json.dumps(layout,indent=2))
    pd.DataFrame(TEXT_QA).to_csv(out/'qa/main_text_geometry.tsv',sep='\t',index=False)
    pd.DataFrame(LOG).to_csv(out/'qa/proposal_panel_provenance.tsv',sep='\t',index=False)
    (out/'qa/redraw_environment.json').write_text(json.dumps({'font':font,'font_requested':'Arial','fallback':font!='Arial','matplotlib':matplotlib.__version__,'note':'No font files are included. Composite PDFs remain vector.'},indent=2))
    print('Redrawn Figures 3-6 from saved inputs; vector assembly remains separate. Font:',font)
    print('Output:',out)

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input',type=Path,default=Path(__file__).parent/'source_subset')
    parser.add_argument('--output',type=Path,default=Path(__file__).parent/'display_proposals')
    args=parser.parse_args();run(args.input,args.output)
