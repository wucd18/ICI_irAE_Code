"""Create independent vector panels from saved results, then describe clean pages."""
from plot_main import *
from matplotlib.ticker import MaxNLocator,PercentFormatter
import string

def class_text_color(key):
 from matplotlib.colors import to_rgb
 rgb=to_rgb(STYLE['class_colors'][key])
 linear=[v/12.92 if v<=.04045 else ((v+.055)/1.055)**2.4 for v in rgb]
 lum=sum(a*b for a,b in zip(linear,[.2126,.7152,.0722]))
 return 'white' if 1.05/(lum+.05) >= (lum+.05)/.05 else 'black'

pages=[]; counts={}; records=[]
TITLE={1:'Cross-cohort assessment of fixed programmes',2:'IFN-adjusted disease associations and observed support',3:'Oxford and Xenium programme assessments',4:'Complete continuous responses across the ligand atlas',5:'Direction, matched-background and pool sensitivities',6:'Representation-dependent concordance and donor responses',7:'Recovered composition and reference-response context'}
def panel(s,letter,w,h,titletext,rect):
 f=plt.figure(figsize=(w/25.4,h/25.4));ax=f.add_axes(rect)
 f.text(.02,.985,letter,fontsize=8,fontweight='bold',va='top');f.text(.09,.985,titletext,fontsize=8,fontweight='bold',va='top')
 ax.tick_params(labelsize=7,length=2,width=.5);return f,ax
def put(f,s,letter,w,h,data):
 pid=f'S{s}{letter}';vector_only(f)
 f.savefig(R/'panels'/f'{pid}.pdf');f.savefig(R/'panels'/f'{pid}.svg');plt.close(f)
 data.to_csv(R/'panel_data'/f'{pid}.tsv',sep='\t',index=False,na_rep='NA');records.append(dict(panel=pid,rows=len(data)))
 return [pid,w,h]
def add(s,slots,legend,land=False):
 counts[s]=counts.get(s,0)+1
 pages.append(dict(figure=f'S{s}',part=counts[s],title=TITLE[s],landscape=land,panels=slots,legend=legend))
def slot(p,x,y):return [p[0],x,y,p[1],p[2]]
def tidy(ax):
 ax.xaxis.set_major_locator(MaxNLocator(3));ax.yaxis.set_major_locator(MaxNLocator(4))
def patient_legend(f,y=.025):
 hs=[Line2D([],[],ls='',marker='o',color=co,label=lab,markersize=3) for co,lab in [(P['muted'],'Control'),(P['red'],'ICI colitis')]]
 f.legend(handles=hs,loc='lower center',bbox_to_anchor=(.52,y),ncol=3,fontsize=7,frameon=False)

# S1: complete enrichment contrasts, saved AUCs and nonduplicated IFN/residual scores.
v=read('Source_022_frozen_module_validation.csv',True)
settings=[('a','External PD1-only',v[v.contrast=='irColitis_vs_ICI_control_PD1_only']),('b','External all-regimen sensitivity',v[v.contrast=='irColitis_vs_ICI_control_all_regimens']),('c','Whole-tissue assessment',v[v.validation_dataset=='GSE210037']),('d','Discovery self-check',v[v.validation_dataset=='GSE253720'])]
slots=[]
for i,(letter,titletext,z) in enumerate(settings):
 assert len(z),(titletext,v.contrast.unique())
 mods=[UP,DEP,IFN,RES,'CURATED_EPITHELIAL_BARRIER','CURATED_COLON_METABOLIC_FUNCTION'];cts=list(z.celltype.unique());z=z.copy();z['display']=np.where(z.expected_direction_in_irAE=='down',-z.NES,z.NES)
 a=z.pivot(index='celltype',columns='module',values='display').reindex(index=cts,columns=mods)
 f,ax=panel(1,letter,185,46,titletext,[.24,.24,.64,.52]);im=heat(ax,a,[PN.get(c,c) for c in cts],['UP','Depleted','IFN','Residual','Barrier','Metabolic'],limit=3,fmt='.2f')
 for ri,ct in enumerate(cts):
  for ci,m in enumerate(mods):
   q=z[(z.celltype==ct)&(z.module==m)]
   if len(q) and q.padj.iloc[0]<.05:ax.text(ci+.37,ri-.17,'*',fontsize=7,ha='center')
 ca=f.add_axes([.93,.26,.012,.47]);f.colorbar(im,cax=ca,ticks=[-3,0,3],extend='both');ca.tick_params(labelsize=7);ca.set_title('NES',fontsize=7)
 slots.append(slot(put(f,1,letter,185,46,z),12.5,30+i*46))
add(1,slots,'a-d, Complete enrichment contrasts for the six fixed patient programmes, separating primary PD1-only assessment, all-regimen sensitivity, whole-tissue samples and the discovery self-check. Colours and numbers show disease-oriented NES: expected-down programmes use minus NES for display. Asterisks retain source-analysis q < 0.05. NE denotes unavailable enrichment, including the original set-size restriction. Metabolic and barrier discordance is retained. Cohort-specific denominators and untransformed NES, P and q are in Table S2.')
auc=read('Source_056_patient_level_frozen_module_classification_performance.csv',True);slots=[]
for i,m in enumerate([UP,DEP,IFN,RES]):
 letter=string.ascii_lowercase[4+i];z=auc[auc.module==m];f,ax=panel(1,letter,92.5,83,MN[m]+' descriptive AUC',[.45,.22,.50,.60]);y=np.arange(len(z))
 ax.hlines(y,z.bootstrap_auc_ci_low,z.bootstrap_auc_ci_high,color=P['blue'],lw=.8);ax.scatter(z.auc_disease_higher,y,c=P['blue'],s=15);ax.axvline(.5,c=P['muted'],ls='--',lw=.5)
 labs=[('Discovery' if r.dataset=='GSE253720' else PN.get(r.celltype,r.celltype))+f' ({r.n_case}/{r.n_control})' for r in z.itertuples()]
 ax.set_yticks(y,labs);ax.invert_yaxis();ax.set_xlim(0,1.03);ax.set_xticks([0,.5,1]);ax.set_xlabel('AUC (95% interval)')
 slots.append(slot(put(f,1,letter,92.5,83,z),12.5+(i%2)*92.5,30+(i//2)*83))
add(1,slots,'e-h, Descriptive case-control AUCs and the 95% percentile intervals for disease-up, depleted, IFN and residual programmes. Labels show cases/controls; GSE206300 families use the primary PD1-only contrast. Whole tissue denotes GSE210037. Discovery denotes the programme-definition cohort. All intervals, including degenerate intervals, retain the existing 5,000 stratified bootstrap result. Depleted patient scores are disease-oriented as 1 minus raw rank. These estimates describe separation in the available samples. Table S2 contains all values.')
slots=[]
for i,m in enumerate([IFN,RES]):
 letter='ij'[i];z=SCO[SCO.module==m];f,ax=panel(1,letter,185,75,MN[m]+' patient scores',[.12,.30,.84,.50]);contexts=list(z[['dataset','celltype']].drop_duplicates().itertuples(index=False,name=None))
 for j,(ds,ct) in enumerate(contexts):
  zz=z[(z.dataset==ds)&(z.celltype==ct)]
  for k,case in enumerate([False,True]):
   s=zz[zz.condition.isin(['CPI_colitis','irColitis','irAE_colitis'])==case].sort_values('patient_id');xx=j+(k-.5)*.30+np.linspace(-.075,.075,len(s));ax.scatter(xx,s.disease_oriented_score,c=P['red'] if case else P['muted'],s=11,edgecolor='white',lw=.2);ax.plot([j+(k-.5)*.3-.12,j+(k-.5)*.3+.12],[s.disease_oriented_score.median()]*2,c=P['ink'],lw=.8)
 ax.set_xticks(range(len(contexts)),[('Discovery' if ds=='GSE253720' else PN.get(ct,ct)) for ds,ct in contexts],rotation=20,ha='right');ax.set_ylabel('Mean expression rank\n(unitless)');patient_legend(f)
 slots.append(slot(put(f,1,letter,185,75,z),12.5,30+i*75))
add(1,slots,'i,j, Every patient or released-sample IFN and residual score. Black segments mark group medians; points retain all within-group variation. External family columns are from GSE206300; whole tissue is GSE210037. Discovery is GSE253720. These two programmes are expected-up, so raw and disease-oriented ranks coincide. Disease-up and depleted individual scores appear in Figures 1 and 2. Patients shared between families are not additive biological replicates.')

# S2: five columns, residual c-g above depleted h-l, sorted explicitly by context.
models=read('Source_072_IFN_conditional_models.csv',True); lp=read('Source_073_IFN_leave_patient_out.csv',True)
f,ax=panel(2,'a',185,100,'Before and after IFN adjustment',[.39,.16,.56,.66]);y=np.arange(len(models));ax.hlines(y,models.ci_low,models.ci_high,color=P['blue'],lw=.8);ax.scatter(models.adjusted_beta,y,c=P['blue'],s=16);ax.scatter(models.unadjusted_beta,y,c=P['muted'],s=14,marker='s');zero(ax)
ax.set_yticks(y,[('Discovery' if r.dataset=='GSE253720' else PN.get(r.celltype,r.celltype))+' / '+MN[r.module] for r in models.itertuples()]);ax.invert_yaxis();ax.set_xlabel('Disease coefficient (rank-score units)');a=put(f,2,'a',185,100,models)
f,ax=panel(2,'b',185,64,'Observed IFN-score support',[.39,.24,.56,.53]);ax.scatter(models.IFN_range_overlap,y,c=P['purple'],s=16);ax.set_yticks(y,[('Discovery' if r.dataset=='GSE253720' else PN.get(r.celltype,r.celltype))+' / '+MN[r.module] for r in models.itertuples()]);ax.invert_yaxis();ax.set_xlim(-.03,1.03);ax.set_xlabel('IFN range overlap');b=put(f,2,'b',185,64,models)
add(2,[slot(a,12.5,30),slot(b,12.5,130)],'a, Unadjusted disease coefficients (grey squares) and IFN-adjusted coefficients (blue circles) with pointwise HC3 95% intervals. Rows distinguish external GSE206300 families from the GSE253720 discovery self-check. b, Observed IFN-range overlap for the same rows. Source-model P/q, case/control denominators and omit-one sensitivity ranges are retained in Table S9. Limited observed support is visible, including zero-overlap rows.')
z=lp[(lp.record_type=='overlap_leverage_diagnostic')&lp.IFN_score.notna()&lp.oriented_score.notna()]
contexts=[('GSE206300',ct) for ct in PF]+[('GSE253720','Epithelial')];slots=[];s2order=[]
for rr,m in enumerate([RES,DEP]):
 for cc,(ds,ct) in enumerate(contexts):
  letter=string.ascii_lowercase[2+rr*5+cc];d=z[(z.dataset==ds)&(z.celltype==ct)&(z.module==m)];assert len(d),(ds,ct,m)
  titletext=('Discovery epithelium' if ds=='GSE253720' else PN[ct])+'\n'+('Residual' if m==RES else 'Depleted')
  f,ax=panel(2,letter,54.6,61,titletext,[.27,.28,.66,.45]);co=[P['red'] if c in ['CPI_colitis','irColitis'] else P['muted'] for c in d.condition];ax.scatter(d.IFN_score,d.oriented_score,c=co,s=11,edgecolor='white',lw=.2);tidy(ax);ax.set_xlabel('IFN rank');ax.set_ylabel('Oriented rank' if cc==0 else '')
  slots.append(slot(put(f,2,letter,54.6,61,d),12+cc*54.6,28+rr*61));s2order.append(dict(letter=letter,row=rr+1,column=cc+1,cohort=ds,family=ct,programme=m,n=len(d)))
add(2,slots,'c-g, Residual disease-oriented rank versus IFN rank in the four external epithelial families and discovery epithelium, respectively. h-l, The corresponding disease-depleted scores in the identical column order. Each point is a patient; red denotes ICI colitis and grey controls. Depleted orientation is 1 minus raw rank. These are observed support displays of the observed patients. All data points and the five context columns are preserved.',True)
pd.DataFrame(s2order).to_csv(R/'qa/S2_panel_order.tsv',sep='\t',index=False)

# S3: Oxford all four programmes; Xenium epithelial measurements only.
ox=read('Source_065_new_spatial_patient_scores.csv',True);oq=read('Source_066_new_spatial_frozen_transfer.csv',True);oc=read('Source_067_new_spatial_module_coverage.csv',True);x=read('Source_069_xenium_patient_programme_scores.csv',True);xc=read('Source_071_xenium_targeted_panel_coverage.csv',True)
for source in ['Oxford','Xenium']:
 slots=[]
 for i,m in enumerate([UP,DEP,IFN,RES]):
  letter=string.ascii_lowercase[i+(0 if source=='Oxford' else 4)];d=ox[ox.module==m] if source=='Oxford' else x[(x.module==m)&(x.compartment=='IEC')];cov=(oc if source=='Oxford' else xc);cov=cov[cov.module==m].iloc[0]
  f,ax=panel(3,letter,92.5,83,source+' / '+MN[m],[.22,.29,.71,.49]);groups=['Healthy','ICI colitis','UC'] if source=='Oxford' else list(d.condition.unique())
  patient_points(ax,d,groups,val='oriented_score' if source=='Oxford' else 'oriented_rank',labels=['Healthy','ICI','UC'] if source=='Oxford' else groups);ax.set_ylabel('Disease-oriented rank');ax.locator_params(axis='y',nbins=4)
  f.text(.1,.14,f'Measured / fixed members: {int(cov.available)}/{int(cov.requested)}',fontsize=7)
  if source=='Oxford':f.text(.1,.07,f'Four-programme q = {oq[oq.module==m].FDR_new_four_module_family.iloc[0]:.3f}',fontsize=7)
  else:f.text(.1,.07,f'{d.patient.nunique()} donor labels; descriptive measurement',fontsize=7)
  slots.append(slot(put(f,3,letter,92.5,83,d),12.5+(i%2)*92.5,30+(i//2)*83))
 add(3,slots,('a-d, All four Oxford programmes in all 16 donors: 7 ICI colitis, 3 healthy and 6 UC. Repeated sections are pooled by donor. None of the original four tests passes its FDR threshold; UC is descriptive. ' if source=='Oxford' else 'e-h, All four programme scores in the released Xenium IEC (intestinal epithelial cell) compartment. Donors overlap the discovery study, and targeted-panel membership is incomplete. No Xenium P value is calculated. ')+'Points are donors, black segments group medians; red is ICI colitis, grey controls, teal UC. Depleted scores are 1 minus raw rank; the other programmes retain raw ranks. Exact coverage and results are in Table S3.')

# S4: one full atlas representation, all ligands and all NE, one family per page.
means=D.groupby(['ligand','celltype','module'],sort=False).nonrank.mean().reset_index();sc=read('atlas_display_scales.tsv');conditions=sorted(C.ligand.unique())+['No_Cytomix']
for fi,fa in enumerate(FAM):
 slots=[]
 for half in range(2):
  ligs=conditions[half*40:(half+1)*40];letter=string.ascii_lowercase[fi*2+half];f,ax=panel(4,letter,92.5,198,FN[fa]+f' / {half*40+1}-{half*40+len(ligs)}',[.26,.16,.70,.70]);ax.remove();dd=[]
  for j,m in enumerate(MODS):
   a=f.add_axes([.28+j*.116,.16,.11,.70]);z=means[(means.celltype==fa)&(means.module==m)].set_index('ligand').reindex(ligs);lim=sc[(sc.method=='nonrank')&(sc.module==m)].symmetric_limit.iloc[0]
   im=heat(a,z[['nonrank']], [x.replace('No_Cytomix','No-cytomix') for x in ligs] if j==0 else ['']*len(ligs),[MN[m]],limit=lim,numbers=False);a.xaxis.tick_top();a.set_xticklabels([MN[m]],rotation=55,ha='left',fontsize=7);a.tick_params(axis='y',labelsize=7)
   ca=f.add_axes([.29+j*.116,.095,.085,.008]);f.colorbar(im,cax=ca,orientation='horizontal',ticks=[0],extend='both');ca.tick_params(labelsize=7,length=1);ca.set_xlabel(f'+/-{lim:.2g}',fontsize=7,labelpad=1)
   d=z.reset_index();d['module']=m;d['family']=fa;d['display_limit']=lim;dd.append(d)
  slots.append(slot(put(f,4,letter,92.5,198,pd.concat(dd)),12.5+half*92.5,30))
 add(4,slots,f'{slots[0][0][-1]},{slots[1][0][-1]}, All {len(conditions)-1} ligands and the separate no-cytomix reference in {FN[fa]}. Colours show arithmetic means of three donor changes in log2(size-factor-normalized count + 1), separately for the six fixed programmes. Grey NE cells denote unavailable three-donor comparisons. Per-programme colour limits are fixed across families; colourbar extensions indicate saturation, while Table S5 preserves every value without clipping. Disease-up decreases and depleted increases oppose their patient directions. Raw and background-relative representations are in Table S5; all GSEA NES and current q values, with separate correction scopes, are in Table S4.')

# S5: joint directions, background imbalance and shared-control sensitivity.
jd=read('raw_joint_direction_table.tsv');ref=read('random_reference_summary.tsv');bal=read('matching_balance.tsv');reuse=read('control_reuse_summary.tsv');lo=read('leave_one_pool_out_classification.tsv')
f,ax=panel(5,'a',85,79,'Raw joint donor-direction counts',[.24,.25,.68,.56]);z=jd[jd.family=='ALL'];mat=z.pivot(index='depleted_up_donor_count',columns='up_down_donor_count',values='n_combinations').reindex(index=range(4),columns=range(4));heat(ax,mat,range(4),range(4),limit=mat.max().max(),fmt='.0f');ax.set_xlabel('Donors with UP decrease');ax.set_ylabel('Donors with depleted increase');a=put(f,5,'a',85,79,z)
f,ax=panel(5,'b',100,79,'Negative changes and matched backgrounds',[.19,.32,.77,.49]);z=ref[ref.stratum=='module'];ii=np.arange(len(z))
for off,col,color,lab in [(-.24,'raw_negative_fraction',P['blue'],'Raw'),(0,'random_set_negative_fraction_mean',P['muted'],'Matched sets'),(.24,'background_relative_negative_fraction',P['teal'],'Relative')]:ax.bar(ii+off,z[col],width=.23,color=color,label=lab)
ax.set_xticks(ii,[MN[m] for m in z.module],rotation=30,ha='right');ax.set_ylim(0,1);ax.set_ylabel('Fraction negative');f.legend(loc='lower center',ncol=3,bbox_to_anchor=(.55,.01),fontsize=7,frameon=False);b=put(f,5,'b',100,79,z)
f,ax=panel(5,'c',185,80,'Remaining continuous matching imbalance',[.30,.29,.63,.53]);dims=list(bal.dimension.unique())
for i,dim in enumerate(dims):
 z=bal[bal.dimension==dim];ax.scatter(z.standardized_mean_difference,np.linspace(i-.18,i+.18,len(z)),c=[FC[x] for x in z.celltype],s=12,alpha=.75)
zero(ax);ax.set_yticks(range(len(dims)),['log2(control CPM + 1)' if d.startswith('log2') else 'Cell detection fraction' for d in dims]);ax.set_xlabel('Standardized mean difference');c=put(f,5,'c',185,80,bal)
add(5,[slot(a,12.5,30),slot(b,97.5,30),slot(c,12.5,109)],'a, Counts across all eligible ligand-family combinations by the number of donors with UP decrease and depleted increase. This joint count table differs from requiring both directions in the same donor. b, Fractions of negative raw changes, expression/detection-matched random-set changes and background-relative changes. c, All continuous matching-balance values, coloured by family; matching does not eliminate all imbalance. These summaries are descriptive and share donors and controls. Matched backgrounds provide context; cytomix remains the experimental control. Detailed effects, balance and direction summaries appear in Table S6.')
f,ax=panel(5,'d',185,91,'Shared matched controls across experimental pools',[.22,.20,.72,.64]);z=reuse.copy();z['unit']=z.donor+' / '+z.celltype.map(FN);mat=z.pivot(index='unit',columns='pool',values='n_ligands');heat(ax,mat,mat.index,[x.replace('Pool ','') for x in mat.columns],limit=mat.max().max(),fmt='.0f');ax.set_xlabel('Experimental pool');d=put(f,5,'d',185,91,z)
orig=C[C.original_both][['ligand','family']].rename(columns={'family':'celltype'});z=lo.merge(orig,on=['ligand','celltype']);f,ax=panel(5,'e',185,112,'Original eight combinations after each pool deletion',[.43,.15,.50,.72]);a0=np.array([[CLASSES.index(r[m+'_class']) if r[m+'_class'] in CLASSES else 4 for m in METHODS] for _,r in z.iterrows()]);vector_matrix(ax,a0,ListedColormap([STYLE['class_colors'][c] for c in CLASSES]),-.5,4.5);ax.set_yticks(range(len(z)),[r.ligand+' / '+FN[r.celltype]+' / omit '+r.dropped_pool.replace('Pool ','') for _,r in z.iterrows()]);ax.set_xticks(range(3),['Raw rank','Relative rank','Non-rank']);ax.tick_params(length=0)
for i in range(len(z)):
 for j in range(3):ax.text(j,i,CN[CLASSES[a0[i,j]]],ha='center',va='center',fontsize=7,color=class_text_color(CLASSES[a0[i,j]]))
e=put(f,5,'e',185,112,z)
add(5,[slot(d,12.5,30),slot(e,12.5,121)],'d, Number of ligand comparisons sharing each donor-family-pool cytomix control. e, Every pool-deletion scenario for all original eight raw-rank Both combinations. Both requires both disease-opposing directions in each of three donors. NE is retained when the remaining donor requirement is unmet, including BMP6-Goblet. All independent-pool and deletion results, including additional and non-Both combinations, are in Table S6. No pool median replaces the matched cytomix control.')

# S6: complete original/additional Both set plus explanatory near-zero and reference responses.
z=C[C.original_both|C.new_bg_both|C.new_nonrank_both];f,ax=panel(6,'a',185,100,'All original and additional Both combinations',[.34,.18,.59,.65]);aa=np.array([[CLASSES.index(r[m+'_class']) for m in METHODS] for _,r in z.iterrows()]);vector_matrix(ax,aa,ListedColormap([STYLE['class_colors'][c] for c in CLASSES]),-.5,4.5);ax.set_yticks(range(len(z)),[r.ligand+' / '+FN[r.family]+(' (original)' if r.original_both else ' (additional)') for _,r in z.iterrows()]);ax.set_xticks(range(3),['Raw rank','Relative rank','Non-rank']);ax.tick_params(length=0)
for i in range(len(z)):
 for j in range(3):ax.text(j,i,CN[CLASSES[aa[i,j]]],ha='center',va='center',fontsize=7,color=class_text_color(CLASSES[aa[i,j]]))
a=put(f,6,'a',185,100,z)
z=C[C.original_both];f,ax=panel(6,'b',185,78,'UP-decreasing donor counts for the original eight',[.31,.24,.42,.57])
for i,r in enumerate(z.itertuples()):
 for j,(col,ma) in enumerate([('raw_up_concordance','o'),('bg_up_concordance','s'),('nonrank_up_concordance','^')]):ax.scatter(getattr(r,col),i+(j-1)*.16,marker=ma,s=18,color=FC[r.family],edgecolor='white',lw=.2)
ax.set_yticks(range(len(z)),[r.ligand+' / '+FN[r.family] for r in z.itertuples()]);ax.invert_yaxis();ax.set_xlim(-.1,3.2);ax.set_xticks([0,1,2,3]);ax.set_xlabel('Donors with disease-up decrease');ax.tick_params(axis='y',length=0)
f.text(.76,.60,'Circle: raw\nSquare: relative\nTriangle: non-rank',fontsize=7,linespacing=1.5);f.text(.76,.22,'Depleted increases:\n3/3 in all eight\nunder all methods.',fontsize=7,linespacing=1.4);b=put(f,6,'b',185,78,z)
add(6,[slot(a,12.5,30),slot(b,12.5,130)],'a, Every combination classified Both under at least one representation, retaining all eight original raw-rank Both combinations and all six additions. b, UP-decreasing donor counts under the three representations for the original eight; vertical offsets separate coincident symbols. TGFB2-Stem has one UP-decreasing donor under relative ranks and two under non-rank scores; TNFSF12-TA has two under both. Both have three under raw ranks. All eight retain depleted increases in three donors in each representation. Table S5 contains every planned comparison, including Other, discordant and NE combinations. These classes describe directions in the same three donors.')
groups=[('IFNG','Stem'),('TNFA','Stem'),('TNFSF12','Transit_Amplifying'),('TNFSF12','Stem')]+[('No_Cytomix',fa) for fa in FAM];slots=[]
for rr,me in enumerate(['raw','bg_relative']):
 for cc,m in enumerate([UP,DEP,IFN]):
  letter=string.ascii_lowercase[2+rr*3+cc];w=[85,50,50][cc];xpos=[12.5,97.5,147.5][cc];f,ax=panel(6,letter,w,83,METH[me]+'\n'+MN[m],[.45 if cc==0 else .15,.23,.51 if cc==0 else .80,.55]);donor_axis(ax,D,groups,m,me,labels=[l.replace('No_Cytomix','No-cytomix')+' / '+FN[fa] for l,fa in groups],show_y=cc==0);ax.set_xlabel('Rank change' if me=='raw' else 'Relative rank change');ax.tick_params(labelsize=7)
  if cc==1:leg_donor(f,.005)
  d=D.merge(pd.DataFrame(groups,columns=['ligand','celltype']),on=['ligand','celltype']);d=d[d.module==m];slots.append(slot(put(f,6,letter,w,83,d),xpos,30+rr*83))
add(6,slots,'c-e, Raw-rank UP, depleted and IFN changes for inflammatory references, TNFSF12 contexts with discordance, and all four no-cytomix families. f-h, Corresponding matched-background-relative rank changes. All three donors are displayed with the same circle/square/triangle identities as Figure 4 (DNW14/DNW15/DNW9); black ticks denote donor means. Changes compare each treatment with matched cytomix controls. Near-zero sign differences are retained with no tolerance. The non-rank counterparts appear in Figures 4 and 5 and Table S5. No GSEA significance is attached to donor-score points.')

# S7: all finer-state combinations; no identity reweighting or cNMF process plots.
props=read('within_family_substate_proportions.tsv');comp=read('composition_standardized_expression.tsv');slots=[]
for i,fa in enumerate(['Colonocyte','Goblet']):
 letter='ab'[i];z=props[(props.celltype==fa)&props.substate.str.endswith('_Mature')];ligs=list(z.ligand.unique());f,ax=panel(7,letter,92.5,175,FN[fa]+' mature fraction',[.30,.18,.64,.69])
 for j,lig in enumerate(ligs):
  for k,(dn,ma) in enumerate(DM.items()):
   q=z[(z.ligand==lig)&(z.donor==dn)]
   if q.empty:continue
   r=q.iloc[0];y=j+(k-1)*.17;ax.plot([r.control_proportion,r.treatment_proportion],[y,y],c=FC[fa],lw=.55,alpha=.65);ax.scatter(r.control_proportion,y,marker=ma,s=14,facecolors='white',edgecolors=FC[fa],lw=.6);ax.scatter(r.treatment_proportion,y,marker=ma,s=14,c=FC[fa],edgecolor='white',lw=.2)
 ax.set_yticks(range(len(ligs)),[x.replace('No_Cytomix','No-cytomix') for x in ligs]);ax.invert_yaxis();ax.set_xlim(-.025,1.025);ax.set_xticks([0,.5,1]);ax.xaxis.set_major_formatter(PercentFormatter(1,decimals=0));ax.set_xlabel('Recovered mature-state fraction');ax.tick_params(axis='y',length=0);f.text(.1,.09,'Open: cytomix control\nFilled: treatment',fontsize=7);leg_donor(f,.025);slots.append(slot(put(f,7,letter,92.5,175,z),12.5+i*92.5,30))
add(7,slots,'a,b, Every available mature-substate proportion for Colonocyte and Goblet contexts, including BMP4 and all other assessed ligands. A line connects matched control and treatment for one donor; donor identities use the same shapes as Figure 4. All recovered proportions remain visible, including comparisons with insufficient common support for expression standardization. These proportions measure recovered composition, not cell yield. Table S7 records the substate proportions, common support and all NE reasons.')
slots=[]
for i,(fa,m) in enumerate([(fa,m) for fa in ['Colonocyte','Goblet'] for m in [UP,DEP]]):
 letter='cdef'[i];z=comp[(comp.celltype==fa)&(comp.programme==m)];ok=z[z.status=='ESTIMATED'];f,ax=panel(7,letter,92.5,82,FN[fa]+' / '+MN[m],[.22,.28,.70,.51]);vals=pd.concat([ok.unstandardized_common_delta,ok.standardized_common_delta]);low=min(0,vals.min());high=max(0,vals.max());pad=(high-low)*.08 or .01
 ax.plot([low-pad,high+pad],[low-pad,high+pad],c=P['muted'],lw=.5)
 for dn,ma in DM.items():
  d=ok[ok.donor==dn];ax.scatter(d.unstandardized_common_delta,d.standardized_common_delta,c=FC[fa],marker=ma,s=13,alpha=.65)
 ax.set_xlim(low-pad,high+pad);ax.set_ylim(low-pad,high+pad);tidy(ax);ax.set_xlabel('Observed common-support weights');ax.set_ylabel('Fixed control-state weights');f.text(.1,.09,f'{len(ok)} estimated donor-ligand effects; {len(z)-len(ok)} NE',fontsize=7);slots.append(slot(put(f,7,letter,92.5,82,z),12.5+(i%2)*92.5,30+(i//2)*82))
add(7,slots,'c-f, UP and depleted expression changes before versus after fixing control-state weights in Colonocyte and Goblet contexts. Each point is a donor-ligand result on exactly the same common-support cells; the diagonal denotes equality. Both axes use mean-cell ln(1 + 10,000 counts/library) change, distinct from pseudobulk log2 changes. Unsupported donor-ligand results remain NE in Table S7 and are counted here. Stem and TA lack finer substate annotations, and are not included in this finer-state comparison. All other programme responses are retained in Table S7.')
cor=read('reference_correlations_by_method.tsv');slots=[]
for i,(fa,m) in enumerate([(fa,m) for fa in ['Stem','Colonocyte','Goblet'] for m in [UP,DEP]]):
 letter='ghijkl'[i];d=D[(D.celltype==fa)&(D.condition_class=='ligand')].groupby(['ligand','module']).nonrank.mean().unstack('module');d=d[[m,IBD2]].dropna().reset_index();q=cor[(cor.celltype==fa)&(cor.programme_a==m)&(cor.programme_b==IBD2)&(cor.method=='nonrank')&(cor.donor=='MEAN_OF_3')&(cor.membership=='fixed_members')].iloc[0];assert len(d)==q.n_common_ligands
 f,ax=panel(7,letter,92.5,61,FN[fa]+' / '+MN[m]+' versus IBD2',[.22,.32,.72,.47]);ax.scatter(d[IBD2],d[m],s=10,c=FC[fa],alpha=.65);zero(ax);ax.axhline(0,c=P['muted'],lw=.5);tidy(ax);ax.set_xlabel('IBD2 non-rank change');ax.set_ylabel('UP change' if m==UP else 'Depleted change');f.text(.10,.06,f'n = {len(d)} ligands; Spearman = {q.Spearman:.3f}',fontsize=7);d['saved_rho']=q.Spearman;slots.append(slot(put(f,7,letter,92.5,61,d),12.5+(i%2)*92.5,30+(i//2)*61))
add(7,slots,'g-l, Complete non-rank mean-response scatterplots for Stem, Colonocyte and Goblet, pairing the ICI UP or depleted programme with IBD2 membership. TA is shown in Figure 6c,d. Each point is a common ligand represented by the mean of three paired donors; n is the number of ligands, not patients. Spearman coefficients describe the corresponding ligand responses. All three representations, individual-donor levels and overlap-removed pairs appear in Table S8.')
(R/'inputs/supplement_pages.json').write_text(json.dumps(pages,indent=2),encoding='utf8');pd.DataFrame(records).to_csv(R/'qa/supplement_panel_sources.tsv',sep='\t',index=False)
print('Supplementary pages:',len(pages),'panels:',len(records))
