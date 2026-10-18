from common import *
from analysis_functions import *
from itertools import combinations
from scipy.stats import rankdata
cfg=jr(ROOT/'analysis_config.json');inv=pd.read_csv(ROOT/'inputs/comparison_inventory.tsv',sep='\t')
dose=pd.read_csv(ROOT/'inputs/condition_metadata.tsv',sep='\t');assert not dose.ligand.duplicated().any()
scores=pd.concat([pd.read_csv(ROOT/'results/all_ligand_sample_scores.tsv',sep='\t',float_precision='round_trip'),pd.read_csv(ROOT/'results/IBD_reference_sample_scores.tsv',sep='\t',float_precision='round_trip')],ignore_index=True)
delta=paired_deltas(scores);save(delta,'results/all_readout_donor_deltas_long.tsv')
wide=delta.pivot(index=['ligand','celltype','donor'],columns='module',values='delta').reset_index()
names={UP:'delta_up',DEP:'delta_depleted',IFN:'delta_IFN',RES:'delta_residual',IBD1:'delta_IBD1',IBD2:'delta_IBD2'}
wide=wide.rename(columns=names);donors=sorted(delta.donor.unique());assert len(donors)==3
planned=pd.MultiIndex.from_product([cfg['ligands']+['No_Cytomix'],cfg['families'],donors],names=['ligand','celltype','donor']).to_frame(index=False)
wide=planned.merge(wide,how='left',validate='one_to_one').merge(inv[['ligand','celltype','condition_class','eligibility','action']],validate='many_to_one')
wide=wide.merge(delta[delta.module==UP][['ligand','celltype','donor','pools','treatment_cells','control_cells','value_origin']],how='left',validate='one_to_one')
wide=wide.merge(dose[['ligand','dose','published_protein_name','mapping_status']],how='left',validate='many_to_one')
wide.loc[wide.ligand.eq('IL22'),'mapping_status']='BLOCKED_DOSE_IL11_IL22_SOURCE_CONFLICT'
wide.loc[wide.ligand.eq('No_Cytomix'),'mapping_status']='ASSAY_REFERENCE_NOT_A_LIGAND'
wide['paired_score_observed']=wide.delta_up.notna()&wide.delta_depleted.notna();wide['up_negative']=np.where(wide.delta_up.notna(),wide.delta_up<0,np.nan);wide['depleted_positive']=np.where(wide.delta_depleted.notna(),wide.delta_depleted>0,np.nan)
wide['same_donor_both']=np.where(wide.paired_score_observed,(wide.delta_up<0)&(wide.delta_depleted>0),np.nan)
save(wide,'results/all_ligand_donor_deltas.tsv')
records=[]
for i in inv.itertuples():
    b=wide[(wide.ligand==i.ligand)&(wide.celltype==i.celltype)];observed=b[b.paired_score_observed]
    rec=dict(ligand=i.ligand,celltype=i.celltype,condition_class=i.condition_class,n_donors=len(observed),planned_donors=3,eligibility=i.eligibility,value_origin=i.action,pattern=direction_pattern(observed.delta_up,observed.delta_depleted),n_up_negative=int((observed.delta_up<0).sum()),n_depleted_positive=int((observed.delta_depleted>0).sum()),n_same_donor_both=int(((observed.delta_up<0)&(observed.delta_depleted>0)).sum()),n_IFN_positive=int((observed.delta_IFN>0).sum()),n_IFN_negative=int((observed.delta_IFN<0).sum()),n_IFN_zero=int((observed.delta_IFN==0).sum()),IFN_direction=reference_direction(observed.delta_IFN),up_direction=reference_direction(observed.delta_up),depleted_direction=reference_direction(observed.delta_depleted))
    for d in names.values():
        v=b[d].dropna();rec[d+'_mean']=v.mean();rec[d+'_median']=v.median();rec[d+'_min']=v.min();rec[d+'_max']=v.max();rec[d+'_n']=len(v)
    rec['dose']=b.dose.iloc[0];rec['dose_metadata_status']=b.mapping_status.iloc[0];rec['classification_basis']='3 observed donor directions, not a test or efficacy claim'
    records.append(rec)
summary=pd.DataFrame(records);assert summary.n_donors.tolist()==inv.n_donors.tolist()
save(summary,'results/response_pattern_summary.tsv')
lig=summary[summary.condition_class=='ligand'];eligible=lig[lig.pattern!='NE'];countrows=[]
patterns=['BOTH','UP_ONLY','DEPLETED_ONLY','OTHER_OBSERVED','NE']
for f in ['ALL']+cfg['families']:
    b=lig if f=='ALL' else lig[lig.celltype==f];den=(b.pattern!='NE').sum()
    for pat in patterns:
        k=b[b.pattern==pat];countrows.append(dict(celltype=f,pattern=pat,n=len(k),planned_comparisons=len(b),eligible_comparisons=den,denominator_for_fraction=len(b) if pat=='NE' else den,fraction=len(k)/(len(b) if pat=='NE' else den),ligands=';'.join(k.ligand),comparisons=';'.join(k.ligand+'|'+k.celltype)))
save(pd.DataFrame(countrows),'results/response_pattern_counts.tsv')
# Each pair of families is evaluated only on jointly eligible ligands.
cross=[];detail=[]
for a,b in combinations(cfg['families'],2):
    x=eligible[eligible.celltype==a][['ligand','pattern']].merge(eligible[eligible.celltype==b][['ligand','pattern']],on='ligand',suffixes=('_a','_b'),validate='one_to_one');x['family_a']=a;x['family_b']=b;x['different_pattern']=x.pattern_a!=x.pattern_b;detail.append(x)
    cross.append(dict(family_a=a,family_b=b,n_common_ligands=len(x),n_different_pattern=int(x.different_pattern.sum()),different_ligands=';'.join(x.loc[x.different_pattern,'ligand'])))
save(pd.concat(detail),'results/cross_family_pattern_details.tsv');save(pd.DataFrame(cross),'results/cross_family_pattern_summary.tsv')
scope=[]
for l,b in eligible.groupby('ligand'):
    scope.append(dict(ligand=l,n_eligible_families=len(b),n_observed_patterns=b.pattern.nunique(),different_patterns_in_observed_families=b.pattern.nunique()>1,eligible_families=';'.join(b.celltype),patterns=';'.join(b.pattern)))
save(pd.DataFrame(scope),'results/ligand_family_pattern_inventory.tsv')
bench=[];pairs=[]
for family in cfg['families']:
    for donor in ['MEAN_OF_3']+donors:
        if donor=='MEAN_OF_3':x=eligible[eligible.celltype==family].rename(columns={c:c.removesuffix('_mean') for c in eligible if c.endswith('_mean')})
        else:x=wide[(wide.celltype==family)&(wide.donor==donor)&wide.eligibility.eq('ELIGIBLE')&wide.condition_class.eq('ligand')]
        for disease in ['delta_up','delta_depleted']:
            for ref in ['delta_IFN','delta_IBD1','delta_IBD2']:
                b=x[['ligand',disease,ref]].dropna();a=b[disease].to_numpy();v=b[ref].to_numpy()
                r=np.corrcoef(a,v)[0,1] if len(b)>2 and np.std(a)>0 and np.std(v)>0 else np.nan
                rho=np.corrcoef(rankdata(a),rankdata(v))[0,1] if len(b)>2 and np.std(a)>0 and np.std(v)>0 else np.nan
                rec=dict(celltype=family,stratum=donor,disease_readout=disease,reference=ref,n_common_ligands=len(b),pearson_r=r,spearman_rho=rho,n_same_nonzero_sign=int((a*v>0).sum()),n_opposite_sign=int((a*v<0).sum()),n_zero_either=int(((a==0)|(v==0)).sum()),P_value=np.nan,inference='Descriptive; common donors/controls, no independent-grid P value',reference_method='Fixed IFN membership raw-rank score' if ref=='delta_IFN' else 'Published IBD membership raw-rank score, not original publisher activity score')
                bench.append(rec)
                for _,row in b.iterrows():pairs.append(dict(ligand=row.ligand,celltype=family,stratum=donor,disease_readout=disease,reference=ref,disease_delta=row[disease],reference_delta=row[ref],same_nonzero_sign=row[disease]*row[ref]>0,opposite_sign=row[disease]*row[ref]<0))
save(pd.DataFrame(bench),'results/programme_benchmark.tsv');save(pd.DataFrame(pairs),'results/benchmark_all_pairs.tsv')
# All direction conflicts are preserved, then examples selected using the locked rule.
discord=[]
for r in eligible.itertuples():
    for disease,dr in [('UP',r.up_direction),('DEPLETED',r.depleted_direction)]:
        if dr!=r.IFN_direction:discord.append(dict(ligand=r.ligand,celltype=r.celltype,disease_readout=disease,disease_direction=dr,IFN_direction=r.IFN_direction,category=disease+'_'+dr+'_IFN_'+r.IFN_direction,n_donors=r.n_donors))
dis=pd.DataFrame(discord);save(dis,'results/disease_IFN_all_direction_discordances.tsv')
selected=[]
for (f,p),b in eligible.groupby(['celltype','pattern']):selected.append(dict(ligand=b.sort_values('ligand').iloc[0].ligand,celltype=f,selection_stratum='PATTERN_'+p,rule='Lexicographically first eligible ligand within family and pattern, locked before analysis'))
for (f,p),b in dis.groupby(['celltype','category']):selected.append(dict(ligand=b.sort_values('ligand').iloc[0].ligand,celltype=f,selection_stratum=p,rule='Lexicographically first eligible ligand within family and disease/IFN direction-discordance category, locked before analysis'))
save(pd.DataFrame(selected),'results/prespecified_display_examples.tsv')
save(lig[lig.pattern=='NE'][['ligand','celltype','n_donors','eligibility','dose','dose_metadata_status']],'results/non_estimable_contrasts.tsv')
print('PATTERNS',eligible.pattern.value_counts().to_dict(),'ELIGIBLE',len(eligible),'NE',len(lig)-len(eligible),flush=True)
print('TNFSF12',summary[summary.ligand=='TNFSF12'][['celltype','pattern','delta_up_mean','delta_depleted_mean','n_same_donor_both']].to_dict('records'),flush=True)
