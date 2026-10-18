from common import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import Normalize,ListedColormap
from matplotlib.lines import Line2D
plt.rcParams.update({'font.family':'Arial','font.size':8,'axes.titlesize':10,'axes.labelsize':8,'pdf.fonttype':42,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False})
COL=dict(zip(FAMS,['#4B79A5','#C47729','#398D80','#A46991']));SHAPE=dict(zip(DONORS,['o','s','^']))
FL={'Stem':'Stem','Transit_Amplifying':'TA','Colonocyte':'Colonocyte','Goblet':'Goblet'}
METHODS=['raw','bg_relative','nonrank'];ML={'raw':'Raw rank change','bg_relative':'Background-relative rank change','nonrank':'Normalized-expression change'}
CATS=['BOTH','UP_ONLY','DEPLETED_ONLY','OTHER_OBSERVED','NE'];CC=['#178A7B','#467BA4','#C29032','#C9CCD1','#ECECEC'];AB={'BOTH':'B','UP_ONLY':'U','DEPLETED_ONLY':'D','OTHER_OBSERVED':'O','NE':'NE'}
cl=read(ROOT/'results/classification_stability.tsv');ref=read(ROOT/'results/No_Cytomix_classification.tsv');allcl=pd.concat([cl,ref]);effects=read(ROOT/'results/three_method_donor_effects.tsv');means=effects.groupby(['ligand','celltype','module'],as_index=False)[METHODS].mean();original=cl[cl.original_both].copy()
def footer(fig,text):fig.text(.035,.022,text,fontsize=7,color='#444444',va='bottom')
def title(fig,text,sub):fig.suptitle(text,x=.035,ha='left',y=.975,fontsize=15,fontweight='bold');fig.text(.035,.933,sub,fontsize=8,ha='left',color='#444444')
def saveone(fig,name):
 fig.savefig(ROOT/'figures'/(name+'.pdf'));fig.savefig(ROOT/'figures'/(name+'.png'),dpi=180);fig.savefig(ROOT/'figures'/(name+'.svg'));plt.close(fig)
def legends(fig,y=.08):
 handles=[Line2D([],[],color=COL[f],marker='o',lw=0,label=FL[f]) for f in FAMS]+[Line2D([],[],color='#333333',marker=SHAPE[d],lw=0,label=d) for d in DONORS]
 fig.legend(handles=handles,loc='lower center',bbox_to_anchor=(.5,y),ncol=7,frameon=False,fontsize=8)

# Figure 4: direction, matched background, stability and superfamily context.
fig=plt.figure(figsize=(11.69,8.27));gs=fig.add_gridspec(2,3,left=.065,right=.975,bottom=.13,top=.87,hspace=.57,wspace=.4)
ax=fig.add_subplot(gs[0,0]);m=read(ROOT/'results/raw_marginal_direction_counts.tsv');m=m[m.family=='ALL']
for i,(a,label,c) in enumerate([('up_down','Disease-up < 0','#467BA4'),('depleted_up','Disease-depleted > 0','#C29032')]):
 b=m[m.axis==a].sort_values('donor_count');x=np.arange(4)+(i-.5)*.35;bars=ax.bar(x,b.n_combinations,.35,label=label,color=c);ax.bar_label(bars,fontsize=7,padding=2)
ax.set(xticks=np.arange(4),xlabel='Concordant donors / 3',ylabel='Ligand-family combinations',title='A  Raw directional asymmetry');ax.legend(fontsize=7,frameon=False);ax.set_ylim(0,m.n_combinations.max()*1.24)
ax=fig.add_subplot(gs[0,1]);j=read(ROOT/'results/raw_joint_direction_table.tsv');j=j[j.family=='ALL'];mat=j.pivot(index='depleted_up_donor_count',columns='up_down_donor_count',values='n_combinations').reindex(index=range(4),columns=range(4)).to_numpy();ax.imshow(mat,origin='lower',cmap='Blues',aspect='auto')
for y in range(4):
 for x in range(4):ax.text(x,y,str(mat[y,x]),ha='center',va='center',color='white' if mat[y,x]>mat.max()*.55 else '#222222')
ax.set(xticks=range(4),yticks=range(4),xlabel='Disease-up: donors < 0',ylabel='Disease-depleted: donors > 0',title='B  Same-combination joint counts')
ax=fig.add_subplot(gs[0,2]);b=effects[effects.condition_class=='ligand'];groups=[]
for m0 in [UP,DEP]:
 q=b[b.module==m0];groups += [q.raw.values,q.B.values]
bp=ax.boxplot(groups,positions=[1,2,4,5],widths=.6,showfliers=False,patch_artist=True)
for i,box in enumerate(bp['boxes']):box.set_facecolor('#467BA4' if i%2==0 else '#B9BDC5')
ax.axhline(0,color='#777777',lw=.7);ax.set(xticks=[1.5,4.5],xticklabels=['Disease-up','Disease-depleted'],ylabel='Paired rank delta',title='C  Raw change and matched background');ax.legend(handles=[Line2D([],[],color='#467BA4',lw=6,label='Raw'),Line2D([],[],color='#B9BDC5',lw=6,label='Background median B')],fontsize=7,frameon=False);ax.text(.02,-.25,'Boxes: donor-pair distributions; shared controls.\nGene-background B is not a biological null.',transform=ax.transAxes,fontsize=7)
ax=fig.add_subplot(gs[1,0]);bottom=np.zeros(3)
for cat,c in zip(CATS,CC):
 v=np.array([(cl[me+'_class']==cat).sum() for me in METHODS]);ax.bar(range(3),v,bottom=bottom,color=c,label=cat.replace('_',' '));
 for k,n in enumerate(v):
  if n>=7:ax.text(k,bottom[k]+n/2,str(n),ha='center',va='center',fontsize=7)
 bottom+=v
ax.set(xticks=range(3),xticklabels=['Raw','Background\nrelative','Nonrank'],ylabel='Planned combinations',title='D  Classification across methods');fig.legend(*ax.get_legend_handles_labels(),loc='lower center',bbox_to_anchor=(.36,.045),ncol=5,fontsize=7,frameon=False)
ax=fig.add_subplot(gs[1,1]);a=original.sort_values(['ligand','family']);cm=ListedColormap(CC);v=np.array([[CATS.index(r[me+'_class']) for me in METHODS] for _,r in a.iterrows()]);ax.imshow(v,cmap=cm,vmin=-.5,vmax=4.5,aspect='auto')
for y,row in enumerate(v):
 for x,k in enumerate(row):ax.text(x,y,AB[CATS[k]],ha='center',va='center',fontsize=8)
ax.set(yticks=range(len(a)),yticklabels=[r.ligand+' / '+FL[r.family] for r in a.itertuples()],xticks=range(3),xticklabels=['Raw','Background\nrelative','Nonrank'],title='E  Every original Both tracked');ax.tick_params(axis='y',labelsize=7)
ax=fig.add_subplot(gs[1,2]);fr=read(ROOT/'results/family_response_summary.tsv');q=fr[(fr.family=='ALL')&(fr.classification=='BOTH')]
for k,gr in enumerate(['TGF_beta_superfamily','Other']):
 z=q[q.ligand_group==gr].set_index('method').loc[METHODS];v=z.n_combinations/z.n_eligible_combinations*100;x=np.arange(3)+(k-.5)*.34;bars=ax.bar(x,v,.34,color=['#825B8D','#8A929D'][k],label='TGF-beta group' if k==0 else 'Other ligands')
 for xx,yy,r in zip(x,v,z.itertuples()):ax.text(xx,yy+.65,f'{r.n_combinations}/{r.n_eligible_combinations}',rotation=0,ha='center',fontsize=6.5)
ax.set(xticks=range(3),xticklabels=['Raw','Background\nrelative','Nonrank'],ylabel='Both / eligible combinations (%)',title='F  Exploratory superfamily context');ax.set_ylim(0,max(10,ax.get_ylim()[1]*1.25));ax.legend(fontsize=7,frameon=False,loc='upper left')
title(fig,'Figure 4 | Directional and scoring sensitivity','Fixed members and eligibility; matched cytomix controls. Both requires disease-up < 0 and disease-depleted > 0 in all three donors.')
footer(fig,'B/U/D/O = Both / UP-only / depleted-only / other. NE retained. No new inferential P values; GSEA q values do not test these scores.')
saveone(fig,'Figure4_Directional_Background_Sensitivity')

# Full atlas: every planned combination, all six programmes, all three methods, No_Cytomix reference.
limits=[];lim={}
for me in METHODS:
 for mod in MODS:
  v=means.loc[means.module==mod,me].values;l=float(np.quantile(abs(v),.95));assert l>0;lim[(me,mod)]=l
  limits.append(dict(method=me,module=mod,symmetric_limit=l,raw_min=v.min(),raw_max=v.max(),display_rule='all-family all-condition absolute mean effect 95th percentile; saturated cells marked *'))
save(pd.DataFrame(limits),'results/atlas_display_scales.tsv')
with PdfPages(ROOT/'figures/Figure_Full_Three_Method_Atlas.pdf') as pdf:
 for me in METHODS:
  for fam in FAMS:
   labels=sorted(cl.ligand.unique())+['No_Cytomix'];fig=plt.figure(figsize=(8.27,11.69));title(fig,f'{FL[fam]} | {ML[me]}','Mean of three donors; original signs. Grey = NE. * = display saturation; numeric values remain untruncated.')
   gs=fig.add_gridspec(1,7,left=.16,right=.96,bottom=.125,top=.895,wspace=.08,width_ratios=[1]*6+[.58]);q=means[means.celltype==fam];ct=allcl[allcl.family==fam].set_index('ligand');axs=[]
   for j,mod in enumerate(MODS):
    ax=fig.add_subplot(gs[0,j]);axs.append(ax);v=q[q.module==mod].set_index('ligand')[me].reindex(labels).to_numpy();l=lim[(me,mod)];cmap=plt.get_cmap('RdBu_r').copy();cmap.set_bad('#E7E7E7');im=ax.imshow(v[:,None],cmap=cmap,vmin=-l,vmax=l,aspect='auto')
    for y,val in enumerate(v):
     if np.isfinite(val):ax.text(0,y,format(val,'.2g')+('*' if abs(val)>l else ''),ha='center',va='center',fontsize=5.9,color='white' if abs(val)>.65*l else '#202020')
    ax.axhline(len(labels)-1.5,color='#111111',lw=.8);ax.set_xticks([]);ax.set_title(SHORT[mod].replace(' ','\n'),fontsize=8,pad=7);ax.set_yticks(range(len(labels)));ax.set_yticklabels(labels if j==0 else [],fontsize=6.4);ax.tick_params(length=0)
    cax=fig.add_axes([ax.get_position().x0,.083,ax.get_position().width,.009]);cb=fig.colorbar(im,cax=cax,orientation='horizontal',ticks=[-l,0,l]);cb.ax.tick_params(labelsize=5.5);cb.ax.set_xticklabels([format(-l,'.2g'),'0',format(l,'.2g')]);cb.ax.get_xticklabels()[0].set_ha('left');cb.ax.get_xticklabels()[-1].set_ha('right')
   ax=fig.add_subplot(gs[0,6]);ax.set(xlim=(-.5,.5),ylim=(len(labels)-.5,-.5));ax.axis('off');ax.set_title('Class',fontsize=8,pad=7)
   for y,lig in enumerate(labels):ax.text(0,y,AB[ct.loc[lig,me+'_class']],ha='center',va='center',fontsize=6.5)
   footer(fig,'Separate colour scale for each programme and method, fixed across family pages. B/U/D/O/NE: see Figure 4.\nNo_Cytomix is an assay reference outside the 316 planned ligand-family combinations. Full precision: source TSV tables.')
   pdf.savefig(fig);plt.close(fig)

# All required examples, then explicit joint donor plots.
example_keys=set(zip(original.ligand,original.family))|{(l,f) for l in ['IFNG','TNFA','TNFSF12','No_Cytomix'] for f in FAMS}
examples=allcl[[ (r.ligand,r.family) in example_keys for r in allcl.itertuples()]].sort_values(['family','ligand']);save(examples,'results/prespecified_display_examples.tsv')
with PdfPages(ROOT/'figures/Figure_Donor_Examples_and_Joint_Directions.pdf') as pdf:
 for start in range(0,len(examples),12):
  ex=examples.iloc[start:start+12];fig,axes=plt.subplots(1,6,figsize=(11.69,8.27));fig.subplots_adjust(left=.22,right=.97,bottom=.18,top=.82,wspace=.42)
  for col,(me,mod) in enumerate((me,mo) for me in METHODS for mo in [UP,DEP]):
   ax=axes[col];ax.axvline(0,color='#777777',lw=.7)
   for y,r in enumerate(ex.itertuples()):
    b=effects[(effects.ligand==r.ligand)&(effects.celltype==r.family)&(effects.module==mod)]
    if b.empty:ax.text(0,y,'NE',ha='center',va='center',fontsize=7)
    for d,z in b.set_index('donor').iterrows():ax.scatter(z[me],y+(DONORS.index(d)-1)*.17,color=COL[r.family],marker=SHAPE[d],s=23,edgecolor='white',lw=.25)
   ax.set_ylim(len(ex)-.5,-.5);ax.set_yticks(range(len(ex)));ax.set_yticklabels([r.ligand+' / '+FL[r.family] for r in ex.itertuples()] if col==0 else [],fontsize=8);ax.tick_params(axis='x',labelsize=7);ax.set_title(('Raw' if me=='raw' else 'Background relative' if me=='bg_relative' else 'Nonrank')+'\n'+SHORT[mod],fontsize=9);ax.set_xlabel('Actual delta');ax.grid(axis='y',alpha=.15)
  title(fig,'Prespecified examples | All donor effects',f'Examples {start+1}-{start+len(ex)} of {len(examples)}. All original Both, IFNG/TNFA, TNFSF12 and No_Cytomix; NE retained.')
  legends(fig,.065);footer(fig,'These are three different questions and units. No GSEA significance outlines. No_Cytomix shares the corresponding method scale.');pdf.savefig(fig);plt.close(fig)
 for fam in FAMS:
  e=examples[examples.family==fam];fig,axes=plt.subplots(1,3,figsize=(11.69,8.27));fig.subplots_adjust(left=.075,right=.97,bottom=.30,top=.80,wspace=.3)
  for ax,me in zip(axes,METHODS):
   pts=[];annotations=[]
   for idx,r in enumerate(e.itertuples(),1):
    b=effects[(effects.ligand==r.ligand)&(effects.celltype==fam)].pivot(index='donor',columns='module',values=me)
    if b.empty:continue
    for d,z in b.iterrows():ax.scatter(-z[UP],z[DEP],color=COL[fam],marker=SHAPE[d],s=33,edgecolor='white',lw=.35)
    annotations.append(ax.annotate(str(idx),(-b[UP].mean(),b[DEP].mean()),xytext=(4,4),textcoords='offset points',fontsize=8,fontweight='bold',arrowprops=dict(arrowstyle='-',lw=.35,color='#777777')))
   ax.axhline(0,color='#777777',lw=.7);ax.axvline(0,color='#777777',lw=.7);ax.margins(x=.25,y=.22);ax.set_title(ML[me],fontsize=9);ax.set_xlabel('-delta disease-up');ax.set_ylabel('+delta disease-depleted');ax.tick_params(labelsize=7)
   fig.canvas.draw();renderer=fig.canvas.get_renderer();placed=[]
   for annotation in annotations:
    for dx,dy in [(4,4),(4,-12),(-13,4),(-13,-12),(14,14),(14,-22),(-23,14),(-23,-22),(25,0),(-30,0)]:
     annotation.set_position((dx,dy));box=matplotlib.text.Text.get_window_extent(annotation,renderer).expanded(1.25,1.2)
     if not any(box.overlaps(previous) for previous in placed):break
    placed.append(box)
   for x,y,t in [(1,1,'UP down\nDEP up'),(0,1,'UP up\nDEP up'),(1,0,'UP down\nDEP down'),(0,0,'UP up\nDEP down')]:ax.text(.98 if x else .02,.98 if y else .02,t,transform=ax.transAxes,ha='right' if x else 'left',va='top' if y else 'bottom',fontsize=6.5,color='#555555',bbox=dict(facecolor='white',alpha=.7,edgecolor='none',pad=1))
  title(fig,f'Joint donor directions | {FL[fam]}','Positive x and y denote anti-disease direction only for the two disease programmes. Labels identify ligand means.')
  fig.text(.07,.20,'  |  '.join(f'{i}: {r.ligand}' for i,r in enumerate(e.itertuples(),1)),fontsize=8,wrap=True);legends(fig,.10);footer(fig,'Dots are paired donors; numbers label three-donor means. Shared controls and repeated donors preclude independent-grid inference.');pdf.savefig(fig);plt.close(fig)

# Pool background distributions (never interpreted as additive confound removal).
pool=read(ROOT/'results/pool_background_summary.tsv');pool=pool[pool.ligand!='No_Cytomix'];pools=sorted(pool.pool.unique());reuse=read(ROOT/'results/control_reuse_summary.tsv')
with PdfPages(ROOT/'figures/Figure_Score_Background_by_Pool.pdf') as pdf:
 for fam in FAMS:
  fig,axes=plt.subplots(2,3,figsize=(11.69,8.27));fig.subplots_adjust(left=.07,right=.97,bottom=.14,top=.84,hspace=.48,wspace=.34)
  for ax,mod in zip(axes.flat,MODS):
   vals=[];pos=[];colors=[]
   for j,p in enumerate(pools):
    b=pool[(pool.celltype==fam)&(pool.module==mod)&(pool.pool==p)]
    for k,col in enumerate(['delta_raw','B','matched_background_relative_delta']):vals.append(b[col].to_numpy());pos.append(j+(k-1)*.23);colors.append(['#467BA4','#ADB2BC','#178A7B'][k])
   bp=ax.boxplot(vals,positions=pos,widths=.2,patch_artist=True,showfliers=False)
   for b,c in zip(bp['boxes'],colors):b.set_facecolor(c)
   ax.axhline(0,color='#777777',lw=.6);ax.set(xticks=range(len(pools)),xticklabels=[p.replace('Pool ','') for p in pools],xlabel='Original pool',ylabel='Rank delta',title=SHORT[mod]);ax.tick_params(labelsize=7)
  title(fig,f'Pool sensitivity | {FL[fam]}','Pool-level ranks are a separate sensitivity: a pooled-pseudobulk rank change is not the average of pool-level rank changes.')
  fig.legend(handles=[Line2D([],[],color=c,lw=6,label=l) for c,l in zip(['#467BA4','#ADB2BC','#178A7B'],['Raw delta','Matched-background B','Background-relative delta'])],loc='lower center',bbox_to_anchor=(.5,.065),ncol=3,frameon=False)
  r=reuse[reuse.celltype==fam];footer(fig,f'{len(r)} unique donor-pool controls; ligand reuse per control {r.n_ligands.min()}-{r.n_ligands.max()}. Boxes describe dependent donor-ligand rows, not independent replicates.');pdf.savefig(fig);plt.close(fig)

# Composition and correlation context: full source tables retain all programmes and donors.
sub=read(ROOT/'results/within_family_substate_proportions.tsv');ce=read(ROOT/'results/composition_standardized_expression.tsv');co=read(ROOT/'results/reference_correlations_by_method.tsv')
with PdfPages(ROOT/'figures/Figure_Composition_and_Reference_Context.pdf') as pdf:
 fig,axes=plt.subplots(1,2,figsize=(11.69,8.27));fig.subplots_adjust(left=.14,right=.93,bottom=.18,top=.85,wspace=.30)
 for ax,fam,state in zip(axes,['Colonocyte','Goblet'],['Colonocyte_Mature','Goblet_Mature']):
  w=sub[(sub.celltype==fam)&(sub.substate==state)].pivot(index='ligand',columns='donor',values='proportion_delta').reindex(index=SELECTED,columns=DONORS);cm=plt.get_cmap('RdBu_r').copy();cm.set_bad('#E7E7E7');im=ax.imshow(w,cmap=cm,vmin=-1,vmax=1,aspect='auto')
  for y in range(len(w)):
   for x in range(3):
    v=w.iloc[y,x];ax.text(x,y,'NE' if pd.isna(v) else f'{v:+.2f}',ha='center',va='center',fontsize=8,color='white' if abs(v)>.65 else '#222222')
  ax.set(xticks=range(3),xticklabels=DONORS,yticks=range(len(w)),yticklabels=SELECTED,title=fam+' mature substate');ax.tick_params(labelsize=8)
 fig.colorbar(im,ax=axes,shrink=.7,label='Treatment - matched control recovered proportion');title(fig,'Composition context | Released mature substates','All specified TGF-beta members plus INS, TNFSF12, IFNG, TNFA and No_Cytomix. Fixed eligible family pools.')
 footer(fig,'Recovered proportions are not survival or cell yield. Stem/TA have no finer final annotations. Missing matched family inputs remain NE.');pdf.savefig(fig);plt.close(fig)
 fig,axes=plt.subplots(2,3,figsize=(11.69,8.27));fig.subplots_adjust(left=.075,right=.97,bottom=.15,top=.84,hspace=.5,wspace=.38)
 for ax,mod in zip(axes.flat,MODS):
  for fam in ['Colonocyte','Goblet']:
   b=ce[(ce.celltype==fam)&(ce.programme==mod)&(ce.status=='ESTIMATED')]
   for d in DONORS:
    z=b[b.donor==d];ax.scatter(z.unstandardized_common_delta,z.standardized_common_delta,s=17,c=COL[fam],marker=SHAPE[d],alpha=.75)
  lo=min(ax.get_xlim()[0],ax.get_ylim()[0]);hi=max(ax.get_xlim()[1],ax.get_ylim()[1]);ax.plot([lo,hi],[lo,hi],color='#888888',lw=.7,ls='--');ax.axhline(0,lw=.5,color='#BBBBBB');ax.axvline(0,lw=.5,color='#BBBBBB');ax.set(title=SHORT[mod],xlabel='Unstandardized common-support delta',ylabel='Control-weight standardized delta');ax.tick_params(labelsize=7)
 title(fig,'Composition sensitivity | Before and after fixed weights','Separate cell LogNormalize scale (ln[1 + 10000 counts/library]); identical supported substate weights in both arms.')
 legends(fig,.065);footer(fig,'At least 20 cells per arm/substate and >=80% retained family coverage per arm. A composition-associated change is not automatically an artefact.');pdf.savefig(fig);plt.close(fig)
 for definition in ['fixed_members','symmetric_overlap_removed']:
  fig,axes=plt.subplots(1,3,figsize=(11.69,8.27));fig.subplots_adjust(left=.22,right=.93,bottom=.19,top=.82,wspace=.2);pairs=[(a,b) for a in [UP,DEP] for b in [IFN,IBD1,IBD2]]
  for ax,me in zip(axes,METHODS):
   b=co[(co.membership==definition)&(co.method==me)&(co.donor=='MEAN_OF_3')];v=np.array([[b[(b.programme_a==a)&(b.programme_b==bb)&(b.celltype==f)].Spearman.iloc[0] for f in FAMS] for a,bb in pairs]);im=ax.imshow(v,cmap='RdBu_r',vmin=-1,vmax=1,aspect='auto')
   for y in range(6):
    for x in range(4):ax.text(x,y,f'{v[y,x]:.2f}',ha='center',va='center',fontsize=9,color='white' if abs(v[y,x])>.6 else '#222222')
   ax.set(xticks=range(4),xticklabels=[FL[f] for f in FAMS],yticks=range(6),yticklabels=[SHORT[a]+' vs '+SHORT[bb] for a,bb in pairs] if ax==axes[0] else [],title=ML[me]);ax.tick_params(axis='x',rotation=45,labelsize=8)
  fig.colorbar(im,ax=axes,shrink=.6,label='Descriptive Spearman rho');title(fig,'Reference relationships | '+definition.replace('_',' '),'Across eligible ligands within each family; three-donor mean. Common IDs and each donor are reported in the source table.')
  footer(fig,'No correlation P values. Low rho does not imply independence or ICI specificity. Overlap removal also changes member-set size.');pdf.savefig(fig);plt.close(fig)
print('FIVE_PDF_FIGURES_COMPLETE',flush=True)
