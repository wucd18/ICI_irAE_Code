from common import *
cl=read(ROOT/'results/classification_stability.tsv');raw=read(ROOT/'results/raw_donor_effects.tsv');audit=[]
for i,r in cl.iterrows():
 b=raw[(raw.ligand==r.ligand)&(raw.celltype==r.family)]
 cl.loc[i,'raw_observed_donors_in_saved_source']=b.donor.nunique()
 if r.eligible:continue
 data=json.loads(r.continuous_effects_by_donor);n=0
 for v in b.itertuples():
  if v.module in MODS and np.isfinite(v.delta):data['raw'][v.donor][v.module]=float(v.delta);n+=1
 cl.loc[i,'continuous_effects_by_donor']=json.dumps(data,separators=(',',':'),allow_nan=False)
 cl.loc[i,'control_ids']=';'.join(sorted({s for ss in b.control_ids.dropna().unique() for s in ss.split(';')}))
 audit.append(dict(ligand=r.ligand,family=r.family,raw_observed_donors=b.donor.nunique(),raw_observed_programme_donor_values=n,classification='NE',reason=r.exclusion_reason,new_scores='Not estimated without original complete three-donor comparison'))
save(cl,'results/classification_stability.tsv');save(pd.DataFrame(audit),'results/NE_partial_raw_observation_inventory.tsv')
save(cl[cl.original_both],'results/original_eight_Both_tracking.tsv');save(cl[cl.new_bg_both|cl.new_nonrank_both],'results/additional_Both_combinations.tsv')
print('PRESERVED_PARTIAL_RAW_VALUES',sum(r['raw_observed_programme_donor_values'] for r in audit),'while all original NE classifications remain NE')
