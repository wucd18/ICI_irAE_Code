"""Format saved results as nine publication tables; no model or score is run."""
from pathlib import Path
import json,hashlib,math,shutil,csv,os
import pandas as pd
import numpy as np
P=Path(os.environ['ICI_PRODUCTION_ROOT']); B=P.parent; R=Path(os.environ['ICI_DISPLAY_INPUT_ROOT'])
UP='ICI_COLITIS_EPITHELIAL_UP_STRICT'; DEP='ICI_COLITIS_EPITHELIAL_LOSS_STRICT'; IFN='CURATED_IFN_VISIBILITY'; RES='ICI_COLITIS_EPITHELIAL_UP_NON_IFN'
IBD1='IBD1_MEMBERSHIP_RANK';IBD2='IBD2_MEMBERSHIP_RANK'
MODS=[UP,DEP,IFN,RES,IBD1,IBD2]
MN={UP:'Disease-up',DEP:'Disease-depleted',IFN:'IFN reference',RES:'Residual (IFN-excluded)',IBD1:'IBD1 membership',IBD2:'IBD2 membership','CURATED_COLON_METABOLIC_FUNCTION':'Metabolic reference','CURATED_EPITHELIAL_BARRIER':'Barrier reference'}
FN={'Transit_Amplifying':'TA','whole_spot_sample':'Whole tissue','Epithelial':'Epithelium'}
ME={'raw':'Raw rank','bg_relative':'Background-relative rank','nonrank':'Non-rank'}
UN={'raw':'rank change','bg_relative':'rank change relative to matched genes','nonrank':'mean member log2(normalized count + 1) change'}
CLASS={'BOTH':'Both','UP_ONLY':'UP only','DEPLETED_ONLY':'Depleted only','OTHER_OBSERVED':'Other observed','NE':'NE'}
DON=['DNW14','DNW15','DNW9']; used={}; tables=[]
titles=['Fixed epithelial and reference gene programmes','Cross-cohort assessment of fixed epithelial programmes','Oxford and Xenium assessments of epithelial programmes','Programme enrichment across the organoid ligand atlas','Donor-level programme responses and directional concordance','Sensitivity to matched backgrounds and experimental pools','Epithelial-substate composition and standardized responses','Overlap and response correlations with reference programmes','IFN-adjusted disease associations and observed support']
def read(n,old=False):
 p=R/'source_data'/('old' if old else 'new')/n
 used[str(p)]={'path':str(p),'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'role':'Saved result for publication display','source':'GSE313368' if not old else 'Patient datasets; accession in result rows'}
 return pd.read_csv(p,sep=',' if p.suffix=='.csv' else '\t',float_precision='round_trip')
def clean(v):
 if v is None or (not isinstance(v,(list,dict)) and pd.isna(v)):return 'NA'
 if isinstance(v,(np.bool_,bool)):return 'Yes' if v else 'No'
 if isinstance(v,np.integer):return int(v)
 if isinstance(v,np.floating):return float(v)
 return v
def section(title,d,cols,notes):
 missing=[k for k in cols if k not in d];assert not missing,(title,missing)
 out=d[list(cols)].rename(columns=cols).copy()
 for col in out:
  if col in ['Programme','Programme A','Programme B']:out[col]=out[col].map(lambda v:MN.get(v,v))
  if col in ['Family','Compartment','Analysis level']:out[col]=out[col].map(lambda v:FN.get(v,v))
  if col in ['Joint class','Raw-rank class','Background-relative class','Non-rank class']:out[col]=out[col].map(lambda v:CLASS.get(v,v))
  if col=='Score representation':out[col]=out[col].map(lambda v:ME.get(v,v))
 rows=[[clean(v) for v in row] for row in out.itertuples(index=False,name=None)]
 return dict(title=title,headers=list(out),rows=rows,notes=notes)
def tab(i,scope,sections):
 tables.append(dict(name=f'Table S{i}',title=f'Table S{i}. {titles[i-1]}',scope=scope,sections=sections))
fixed=read('Source_021_frozen_module_definitions.csv',True)
ibdp=Path(os.environ['ICI_IBD_MEMBERSHIP'])
sha=hashlib.sha256(ibdp.read_bytes()).hexdigest();assert sha=='73b408e3383ca81d6ff8cc0cf04c5295f07f0a854c7a61045def9600e73683e4'
used[str(ibdp)]={'path':str(ibdp),'bytes':ibdp.stat().st_size,'sha256':sha,'role':'Complete pinned IBD membership definitions','source':'Capeling et al.; doi:10.1038/s41467-025-68247-6'}
ibd=pd.read_csv(ibdp,sep='\t');ibd['derivation_source']='Capeling et al., supplementary IBD memberships (Ref. 14)';ibd['expected_direction_in_irAE']='Reference; no fixed ICI direction';ibd['selection_rule']=ibd.definition
members=pd.concat([fixed,ibd],ignore_index=True);assert not members.duplicated(['module','gene']).any()
cov=read('fixed_programme_coverage.tsv'); cvmap=dict(zip(['Disease-up','Disease-depleted','IFN','Residual','IBD1 membership','IBD2 membership'],MODS));cov['module']=cov.programme.map(cvmap)
missing={r.module:set(str(r.missing_members).split(';')) if pd.notna(r.missing_members) else set() for r in cov.itertuples()}
members['measured']=members.apply(lambda r:'Yes' if r.gene not in missing[r.module] else 'No' if r.module in missing else 'NA',axis=1) if False else members.apply(lambda r:('No' if r.gene in missing[r.module] else 'Yes') if r.module in missing else 'NA',axis=1)
tab(1,'All fixed members used for the epithelial, IFN, residual, barrier, metabolic and IBD membership readouts.',[
 section('A. Complete fixed membership',members,{'module':'Programme','gene':'Gene symbol','expected_direction_in_irAE':'Expected disease direction','derivation_source':'Definition source','measured':'Measured in organoid atlas','selection_rule':'Fixed definition'},['Gene symbols are text. Memberships were not reselected. NA in organoid coverage means this reference was not used for the organoid scores. Disease-up decreases and disease-depleted increases oppose the patient disease direction. IBD readouts use published membership, not the publisher activity score.']),
 section('B. Organoid fixed and measured coverage',cov,{'module':'Programme','frozen_members':'Fixed members (n)','measured_members':'Measured members (n)','missing_members':'Unmeasured fixed members'},['Missing genes are not imputed. Coverage applies to the saved organoid readouts; patient-assessment coverage appears in Tables S2 and S3.'])])
v=read('Source_022_frozen_module_validation.csv',True); auc=read('Source_056_patient_level_frozen_module_classification_performance.csv',True)
v['unit']=v.feature_level.map(lambda x:'Released sample pseudobulk' if x=='whole_spot_sample' else 'Patient-family pseudobulk'); v.loc[v.validation_dataset=='GSE210037','unit']='Released sample pseudobulk'
vc={'validation_dataset':'Cohort','cohort_role':'Assessment role','celltype':'Family','contrast':'Comparison','module':'Programme','unit':'Analysis unit','n_case':'Cases (n)','n_control':'Controls (n)','expected_direction_in_irAE':'Expected direction','observed_direction':'Observed direction','direction_concordant':'Direction concordant','size':'Ranked members (n)','NES':'NES','pval':'GSEA P','padj':'GSEA q'}
ac={'dataset':'Cohort','cohort_role':'Assessment role','celltype':'Family','module':'Programme','case':'Case group','control':'Control group','n_case':'Cases (n)','n_control':'Controls (n)','median_case_oriented_score':'Case median (rank)','median_control_oriented_score':'Control median (rank)','cliffs_delta':"Cliff's delta",'auc_disease_higher':'Descriptive AUC','bootstrap_auc_ci_low':'AUC 95% CI lower','bootstrap_auc_ci_high':'AUC 95% CI upper','mannwhitney_p':'Mann-Whitney P','fdr_across_module_celltype_tests':'Score-test q'}
sec=[]
for kind,mask in [('External assessment and prespecified regimen sensitivity',v.validation_dataset!='GSE253720'),('Discovery self-check',v.validation_dataset=='GSE253720')]:sec.append(section(kind+' - enrichment',v[mask],vc,['NES is not direction-reversed in this table. Source GSEA q values and their original contrast-specific families are preserved. Metabolic and barrier references are included even when discordant.']))
for kind,mask in [('External assessment',auc.dataset!='GSE253720'),('Discovery self-check',auc.dataset=='GSE253720')]:sec.append(section(kind+' - score effects and descriptive AUC',auc[mask],ac,['AUC and 95% percentile intervals use the existing 5,000 stratified bootstrap samples. The disease-oriented depleted score equals 1 minus raw rank. All source P/q and degenerate intervals are retained. Discovery self-checks use the definition cohort.']))
tab(2,'Fixed programmes in discovery and external patient or released-sample pseudobulks; all saved contrasts and sensitivities are retained.',sec)
ox=read('Source_066_new_spatial_frozen_transfer.csv',True);oxcov=read('Source_067_new_spatial_module_coverage.csv',True);ox=ox.merge(oxcov,on='module',validate='one_to_one')
oxd=read('Source_065_new_spatial_patient_scores.csv',True);oxd=oxd.merge(oxcov[['module','requested']],on='module',validate='many_to_one')
x=read('Source_069_xenium_patient_programme_scores.csv',True);xcov=read('Source_071_xenium_targeted_panel_coverage.csv',True)
xe=x[x.compartment=='IEC'].copy();assert len(xe)>0
xe=xe.merge(xcov[['module','requested']],on='module',validate='many_to_one')
tab(3,'Oxford donor-level four-programme tests and all corresponding donor scores; Xenium epithelial scores are descriptive measurements.',[
 section('A. Oxford ICI colitis versus healthy donors',ox,{'module':'Programme','n_case':'ICI donors (n)','n_control':'Healthy donors (n)','requested':'Fixed members (n)','available':'Measured members (n)','median_case':'ICI median (oriented rank)','median_control':'Healthy median (oriented rank)','median_difference':'Median difference (rank)', 'cliffs_delta':"Cliff's delta",'U':'Mann-Whitney U','P':'P','FDR_new_four_module_family':'Four-programme q'},['The original BH family contains four tests; none passes q < 0.05. UC donors are descriptive and not included in the test.']),
 section('B. Oxford donor measurements',oxd,{'patient':'Donor ID','condition':'Group','n_sections':'Pooled sections (n)','n_spots':'Spots (n)','module':'Programme','requested':'Fixed members (n)','available_genes':'Measured members (n)','raw_score':'Raw rank score','oriented_score':'Disease-oriented rank score'},['All 16 donors are shown (7 ICI colitis, 3 healthy, 6 UC). Repeated sections are aggregated within donor; spots are not independent replicates.']),
 section('C. Xenium epithelial donor measurements',xe,{'patient':'Donor ID','condition':'Group','compartment':'Compartment','n_cells':'Measured cells (n)','module':'Programme','requested':'Fixed members (n)','available_genes':'Measured members (n)','raw_rank':'Raw rank score','oriented_rank':'Disease-oriented rank score'},[f'All {xe.patient.nunique()} epithelial donor labels and all four programmes are retained. This targeted assay has partial membership coverage and overlapping discovery donors. No Xenium group P or q value was calculated.'])])
g=read('original_GSEA_results.tsv'); fam=read('ligand_family_membership.tsv');g=g.merge(fam[['ligand','TGF_beta_superfamily','subgroup']],on='ligand',how='left',validate='many_to_one')
gc={'ligand':'Ligand / reference','celltype':'Family','module':'Programme','TGF_beta_superfamily':'TGF-beta superfamily','subgroup':'Ligand subgroup','dose':'Dose','n_donors':'Matched donors (n)','eligibility':'Eligibility','requested':'Fixed members (n)','size':'Ranked members (n)','NES':'NES','pval':'P','FDR_atlas':'Atlas q','status':'Estimate status'}
gcnotes=['NES and P/q are copied from the latest saved atlas results. NE_FEWER_THAN_3_MATCHED_DONORS means the unchanged three-paired-donor requirement was not met. NA numerical fields are unavailable estimates, not zero. Historical targeted q values are excluded. Dose and ligand mapping retain source labels, including the documented IL22/IL11 metadata discrepancy.']
secs=[]
for key,title in [('all_ligand_primary_948','A. Primary ligand tests - disease-up, depleted and IFN (948 planned)'),('all_ligand_residual_316','B. Residual programme sensitivity (316 planned)'),('no_cytomix_reference_16','C. No-cytomix assay reference (16 planned)')]:
 z=g[g.FDR_scope==key]
 if z.empty:
  if key.startswith('no_cytomix'):z=g[g.condition_class!='ligand']
  else:raise ValueError((key,g.FDR_scope.unique()))
 secs.append(section(title,z,gc,gcnotes+[f'This section preserves its original correction family, including unavailable planned estimates. Saved FDR scope: {z.FDR_scope.iloc[0]}.']))
assert sum(len(s['rows']) for s in secs)==len(g)
tab(4,'All planned comparisons in the 79-ligand, four-family atlas, with the separate no-cytomix reference and unchanged FDR scopes.',secs)
D=read('three_method_donor_effects.tsv');C=read('classification_stability.tsv');NC=read('No_Cytomix_classification.tsv');raw=read('raw_donor_effects.tsv');allc=pd.concat([C,NC],ignore_index=True)
rows=[]
for cr in allc.itertuples():
 for me in ME:
  for m in MODS:
   d=D[(D.ligand==cr.ligand)&(D.celltype==cr.family)&(D.module==m)]
   if d.empty and me=='raw':
    d=raw[(raw.ligand==cr.ligand)&(raw.celltype==cr.family)&(raw.module==m)].rename(columns={'delta':'raw'})
   vals={dn:(float(d.loc[d.donor==dn,me].iloc[0]) if len(d.loc[d.donor==dn]) else np.nan) for dn in DON}
   avail=[x for x in vals.values() if pd.notna(x)]
   cl=getattr(cr,me+'_class')
   rows.append(dict(ligand=cr.ligand,family=cr.family,programme=m,method=ME[me],unit=UN[me],**vals,n=len(avail),mean=np.mean(avail) if len(avail)==3 else np.nan,negative=sum(v<0 for v in avail),positive=sum(v>0 for v in avail),zero=sum(v==0 for v in avail),joint=CLASS[cl],status='Eligible' if cr.eligible else 'NE: fewer than 3 paired donors',original=cr.original_both,additional=cr.new_bg_both or cr.new_nonrank_both,condition=cr.condition_class))
eff=pd.DataFrame(rows);ec={'ligand':'Ligand / reference','family':'Family','programme':'Programme','DNW14':'DNW14 change','DNW15':'DNW15 change','DNW9':'DNW9 change','n':'Available donors (n)','mean':'Three-donor mean change','negative':'Negative donors (n)','positive':'Positive donors (n)','zero':'Zero donors (n)','joint':'Joint UP/depleted class','status':'Eligibility','original':'Original raw Both','additional':'Additional Both'}
secs=[]
for me in ME:
 for cc,lab in [('ligand','Ligand atlas'),('assay_reference','No-cytomix reference')]:
  z=eff[(eff.method==ME[me]) & ((eff.condition=='ligand') if cc=='ligand' else (eff.condition!='ligand'))]
  secs.append(section(ME[me]+' - '+lab+' ('+UN[me]+')',z,ec,['All planned combinations remain visible. Means require three saved paired donors. Partial raw observations for NE comparisons are retained; no new background or non-rank scores were created for them. Joint class uses only disease-up and depleted: Both requires UP < 0 and depleted > 0 in every donor. UP only and depleted only meet only the respective three-donor condition; Other observed meets neither. IFN, residual and IBD rows share that combination-level label. Counts are descriptive and have no GSEA significance.']))
tab(5,'Treatment minus the matched cytomix control in three distinct representations. All ligands, all families and the separate assay reference are retained.',secs)
bg=read('matched_background_effects.tsv');bal=read('matching_balance.tsv');ref=read('random_reference_summary.tsv');reuse=read('control_reuse_summary.tsv');ip=read('independent_pool_classification.tsv');lo=read('leave_one_pool_out_classification.tsv');fs=read('family_response_summary.tsv')
lo['remaining_donor_count']=lo.remaining_donors.map(lambda s:len(str(s).split(';')) if pd.notna(s) and str(s) else 0)
display_cases=pd.concat([C[C.original_both|C.new_bg_both|C.new_nonrank_both|C.ligand.isin(['IFNG','TNFA','TNFSF12'])][['ligand','family']],NC[['ligand','family']]]).drop_duplicates().rename(columns={'family':'celltype'})
bg=bg.merge(display_cases,on=['ligand','celltype'],validate='many_to_one')
bal['library']=bal.library_id.map(lambda x:x.split('_')[-1]) # display-independent library number below
bal['library_number']=bal.groupby(['celltype','module'],sort=False).cumcount()//2+1
sec=[section('A. Matched gene-background effects (rank change)',bg,{'ligand':'Ligand / reference','celltype':'Family','donor':'Donor','module':'Programme','n_members':'Members (n)','B':'Matched-background median B','random_q025':'Random-set 2.5th percentile','random_q975':'Random-set 97.5th percentile','raw_reference_percentile':'Raw effect percentile in random sets'},['B is the saved median effect of expression/detection-matched gene sets. Random-set quantiles describe reference distributions, not confidence intervals or P values. Raw and background-relative effects are in Table S5.']),
section('B. Remaining expression and detection imbalance',bal,{'celltype':'Family','module':'Programme','library_number':'Matched library','dimension':'Matching dimension','n_target':'Target members (n)','target_mean':'Target mean','background_mean':'Background mean','mean_difference':'Background minus target','standardized_mean_difference':'Standardized mean difference','background_mean_q025':'Background mean 2.5th percentile','background_mean_q975':'Background mean 97.5th percentile'},['Expression dimension: log2(control CPM + 1); detection dimension: fraction of cells detecting a gene. Standardized differences are unitless. Multiple saved matched libraries retain separate rows. Continuous imbalance remains after bin matching.']),
section('C. Direction frequencies in matched backgrounds',ref[ref.stratum.isin(['module','family_module','celltype_module','module_celltype','donor_module','module_donor'])],{'module':'Programme','stratum':'Summary level','celltype':'Family','donor':'Donor','n':'Paired effects (n)','raw_negative_fraction':'Raw-negative fraction','random_set_negative_fraction_mean':'Mean random-set-negative fraction','background_relative_negative_fraction':'Relative-negative fraction'},['Fractions are descriptive summaries from the saved result table. Paired effects share donors and controls.']),
section('D. Shared cytomix controls by donor, family and pool',reuse,{'donor':'Donor','celltype':'Family','pool':'Experimental pool','n_ligands':'Ligands sharing control (n)','control_cells':'Control cells (n)','has_No_Cytomix':'No-cytomix reference present'},['The same control can contribute to multiple ligand comparisons. Counts of ligands are not counts of independent donors.']),
section('E. All saved independent-pool classifications',ip,{'ligand':'Ligand','celltype':'Family','pool':'Experimental pool','method':'Score representation','n_matched_donors':'Paired donors (n)','classification':'Joint class','reason':'Status / NE reason'},['Each estimate uses its matched cytomix control within the pool. Three paired donors are required; NE values and reasons are retained.']),
section('F. All saved leave-one-pool-out classifications',lo,{'ligand':'Ligand','celltype':'Family','dropped_pool':'Omitted pool','remaining_donor_count':'Remaining donors (n)','remaining_pools':'Retained pools','status':'Status / NE reason','raw_class':'Raw-rank class','bg_relative_class':'Background-relative class','nonrank_class':'Non-rank class'},['All available deletion scenarios, including every original raw Both combination and all additional Both combinations, are retained. Pool medians do not replace matched cytomix controls.']),
section('G. Descriptive TGF-beta-family response context',fs,{'family':'Family','ligand_group':'Structural ligand group','method':'Score representation','classification':'Joint class','n_combinations':'Combinations in class (n)','n_planned_combinations':'Planned combinations (n)','n_eligible_combinations':'Eligible combinations (n)','n_ligands_with_category':'Ligands in class (n)','n_ligands_planned':'Planned ligands (n)','n_ligands_with_eligible_family':'Ligands with eligible family (n)'},['Groups are structural memberships and preserve the saved denominators. Combinations share donors, controls and pools. No family-level significance test is supplied.']),
section('H. Raw-rank joint distribution of donor-direction counts',read('raw_joint_direction_table.tsv'),{'family':'Family','up_down_donor_count':'UP-decreasing donors (n)','depleted_up_donor_count':'Depleted-increasing donors (n)','n_combinations':'Ligand-family combinations (n)'},['All 16 count pairs are retained overall and for every family. This distribution is distinct from the same-donor joint directions below.']),
section('I. Raw-rank same-donor joint directions',read('raw_same_donor_joint_counts.tsv'),{'family':'Family','donor':'Donor','n_matched_pairs':'Matched pairs (n)','same_donor_both_decline':'Both scores decrease (n)','same_donor_anti_disease':'UP decreases and depleted increases (n)','same_donor_disease_direction':'UP increases and depleted decreases (n)'},['Counts describe paired score directions in the same donor. Shared controls and donors are not treated as independent combination-level replicates.'])]
sec[0]['notes'].append('This section contains every original or additional Both combination, IFNG/TNFA/TNFSF12 in all families with saved effects, and no-cytomix references. Full atlas responses are in Table S5; other background distributions remain in the project archive.')
assert len(sec[2]['rows'])>0
tab(6,'Saved matched-gene and experimental-pool sensitivities; unchanged eligibility and three-donor rules.',sec)
props=read('within_family_substate_proportions.tsv');comp=read('composition_standardized_expression.tsv');multi=comp[comp.annotation_resolution!='NO_FINER_RELEASED_SUBSTATE'];mp=props[props.finer_substates_available]
mem=pd.read_csv(R/'source_data/derived/BMP4_Colonocyte_available_fixed_members.tsv',sep='\t',float_precision='round_trip');p=R/'source_data/derived/BMP4_Colonocyte_available_fixed_members.tsv';used[str(p)]={'path':str(p),'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'role':'Saved member changes supporting Figure 5b','source':'GSE313368'}
wide=mem.pivot(index=['ligand','celltype','module','gene'],columns='donor',values='gene_log_delta').reset_index()
no=comp[comp.annotation_resolution=='NO_FINER_RELEASED_SUBSTATE'][['ligand','celltype','donor','annotation_resolution']].drop_duplicates()
tab(7,'All saved finer-substate contexts, including BMP4 and other assessed ligands. Composition effects retain their cell ln-normalized scale.',[
section('A. Within-family substate proportions',mp,{'ligand':'Ligand / reference','celltype':'Family','donor':'Donor','substate':'Released substate','treatment_proportion':'Treatment proportion','control_proportion':'Control proportion','proportion_delta':'Proportion difference','common_supported':'Common support','reference_weight':'Fixed control-state weight','treatment_retained_fraction':'Treatment retained fraction','control_retained_fraction':'Control retained fraction','status':'Status / NE reason'},['Proportions describe recovered cells, not cell yield. Common support requires the saved minimum 20 cells per arm/substate and 80% retained in each arm.']),
section('B. Responses before and after standardization on the same common support',multi,{'ligand':'Ligand / reference','celltype':'Family','donor':'Donor','programme':'Programme','unstandardized_full_delta':'Full-population change','unstandardized_common_delta':'Common-support change','standardized_common_delta':'Fixed-weight common-support change','n_common_substates':'Common substates (n)','common_substates':'Common substate names','treatment_retained_fraction':'Treatment retained fraction','control_retained_fraction':'Control retained fraction','status':'Status / NE reason'},['All effects in this section use mean-cell ln(1 + 10,000 counts/library) expression. Compare the two common-support columns on the same cells; do not subtract these effects from pseudobulk log2 changes. NE retains insufficient common support.']),
section('C. Contexts without finer released substate annotation',no,{'ligand':'Ligand / reference','celltype':'Family','donor':'Donor','annotation_resolution':'Reason finer-state sensitivity is unavailable'},['Stem and TA have only family-level labels. Their identity reweighting values are archived; they do not measure sensitivity to finer substate composition.']),
section('D. Complete available BMP4-Colonocyte member changes supporting Figure 5b',wide,{'module':'Programme','gene':'Gene symbol','DNW14':'DNW14 log2 change','DNW15':'DNW15 log2 change','DNW9':'DNW9 log2 change'},['Each value is a saved log2(size-factor-normalized count + 1) change. All 84 disease-up and 260 available disease-depleted members are retained for each of three donors. The fixed depleted definition has 275 members; missing members are not imputed. These gene changes are distinct from the cell ln-normalized composition scale.'])])
ov=read('gene_set_overlap.tsv');cor=read('reference_correlations_by_method.tsv')
tab(8,'Complete fixed and measured membership overlaps and all saved response correlations across families, donors and score representations.',[
section('A. Pairwise membership overlap',ov,{'programme_a':'Programme A','programme_b':'Programme B','scope':'Membership scope','n_a':'A members (n)','n_b':'B members (n)','n_intersection':'Intersection (n)','intersection_fraction_a':'Intersection / A','intersection_fraction_b':'Intersection / B','Jaccard':'Jaccard index','common_genes':'Shared members'},['Both denominators are shown. Fixed and measured memberships are distinct. Shared-member names are retained as text.']),
section('B. Full-membership and symmetric overlap-removed response correlations',cor,{'membership':'Membership comparison','celltype':'Family','donor':'Donor / mean','method':'Score representation','programme_a':'Programme A','programme_b':'Programme B','n_common_ligands':'Common ligands (n)','n_a_remaining':'Remaining A members (n)','n_b_remaining':'Remaining B members (n)','Pearson':'Pearson r','Spearman':'Spearman rho','status':'Estimate status'},['MEAN_OF_3 denotes three-donor mean responses. Sample size is the number of common ligands, not patients. raw = raw rank; bg_relative = matched-background-relative rank; nonrank = mean member log2 change. Symmetric removal excludes shared members from both sets. No new fits, intervals or P values are supplied.'])])
models=read('Source_072_IFN_conditional_models.csv',True)
mc={'dataset':'Cohort','role':'Assessment role','celltype':'Family','module':'Programme','n_case':'Cases (n)','n_control':'Controls (n)','unadjusted_beta':'Unadjusted disease coefficient','unadjusted_se_HC3':'Unadjusted HC3 SE','adjusted_beta':'IFN-adjusted disease coefficient','adjusted_se_HC3':'Adjusted HC3 SE','ci_low':'Adjusted 95% CI lower','ci_high':'Adjusted 95% CI upper','p':'Adjusted-model P','fdr':'Adjusted-model q','IFN_range_overlap':'IFN range overlap','LOPO_beta_min':'Omit-one coefficient minimum','LOPO_beta_max':'Omit-one coefficient maximum','LOPO_all_expected_direction':'All omit-one expected direction'}
tab(9,'Disease associations before and after IFN adjustment in source rank-score units, with observed support and saved sensitivity ranges.',[
section('A. External assessment',models[models.dataset!='GSE253720'],mc,['HC3 intervals are pointwise, and q retains the original eight-test external BH family. IFN range overlap is the saved overlap measure, not an extrapolated support estimate. Omit-one ranges summarize existing fits; individual iterations and model objects are archived.']),
section('B. Discovery self-check',models[models.dataset=='GSE253720'],mc,['Discovery checks use the programme-definition cohort and remain separate from external tests. NA denotes an unavailable or inapplicable source statistic. No omitted interval or P value has been reconstructed.'])])
for t in tables:
 for s in t['sections']:
  assert len(s['rows']),s['title']
  s['notes'].append('NA = unavailable or not applicable; NE = non-estimable. Zero is retained only when the source reports zero. Numerical cells retain source precision; display rounding does not change stored values.')
  for row in s['rows']:
   assert len(row)==len(s['headers'])
(P/'inputs/publication_tables.json').write_text(json.dumps(tables,ensure_ascii=False,allow_nan=False),encoding='utf8')
pd.DataFrame(used.values()).to_csv(P/'qa/table_input_manifest.tsv',sep='\t',index=False)
pd.DataFrame([{'table':t['name'],'section':s['title'],'rows':len(s['rows']),'columns':len(s['headers'])} for t in tables for s in t['sections']]).to_csv(P/'qa/table_sections.tsv',sep='\t',index=False)
eff.to_csv(P/'inputs/donor_display_table.tsv',sep='\t',index=False)
print(json.dumps({t['name']:sum(len(s['rows']) for s in t['sections']) for t in tables}))
