from common import *
counts=read(ROOT/'results/recovered_substate_counts.tsv');po=read(BASE/'inputs/all_ligand_pool_eligibility.tsv');po=po[po.eligible]
logs=read(ROOT/'cache/substate_log_expression.tsv');gps=read(ROOT/'cache/substate_saved_programmes.tsv')
props=[];within=[];effects=[];markers=[]
# Recovered major-family proportions: matched original pools, including families outside the four score families in denominator.
for ligand in SELECTED:
 for d in DONORS:
  t=counts[(counts.assignment==d)&(counts.arm==ligand)];c=counts[(counts.assignment==d)&(counts.arm=='Control')]
  ps=sorted(set(t.pool)&set(c.pool))
  for scope,pools in [('POOLED_MATCHED_POOLS',ps)]+[(p,[p]) for p in ps]:
   tb=t[t.pool.isin(pools)].groupby('celltype').n_cells.sum();cb=c[c.pool.isin(pools)].groupby('celltype').n_cells.sum()
   for f in sorted(set(tb.index)|set(cb.index)):
    nt=int(tb.get(f,0));nc=int(cb.get(f,0));props.append(dict(ligand=ligand,donor=d,pool_scope=scope,pools=';'.join(pools),celltype=f,treatment_recovered_cells=nt,control_recovered_cells=nc,treatment_total=int(tb.sum()),control_total=int(cb.sum()),treatment_proportion=nt/tb.sum(),control_proportion=nc/cb.sum(),proportion_delta=nt/tb.sum()-nc/cb.sum()))
 for fam in FAMS:
  for d in DONORS:
   pp=po[(po.ligand==ligand)&(po.celltype==fam)&(po.donor==d)];ps=sorted(pp.pool.unique())
   if not ps:continue
   b=counts[(counts.assignment==d)&(counts.celltype==fam)&counts.pool.isin(ps)&counts.arm.isin([ligand,'Control'])]
   st=b.groupby(['arm','final_annotations']).n_cells.sum().unstack(0,fill_value=0).reindex(columns=[ligand,'Control'],fill_value=0)
   nt=st[ligand].sum();nc=st.Control.sum();support=st[(st[ligand]>=20)&(st.Control>=20)];tcov=support[ligand].sum()/nt;ccov=support.Control.sum()/nc;ok=len(support)>0 and min(tcov,ccov)>=.8
   reason='ESTIMATED' if ok else 'NE_COMMON_SUBSTATE_SUPPORT_BELOW_80_PERCENT_OR_20_CELLS'
   finer=len(st)>1
   for state,r in st.iterrows():
    within.append(dict(ligand=ligand,donor=d,celltype=fam,substate=state,pools=';'.join(ps),treatment_cells=r[ligand],control_cells=r.Control,treatment_proportion=r[ligand]/nt,control_proportion=r.Control/nc,proportion_delta=r[ligand]/nt-r.Control/nc,common_supported=state in support.index,reference_weight=r.Control/support.Control.sum() if state in support.index else np.nan,treatment_retained_fraction=tcov,control_retained_fraction=ccov,status=reason,finer_substates_available=finer))
   for data,value,program,is_marker in [(logs,'mean_cell_log_normalized','module',False),(gps,'saved_cNMF_usage','programme',True)]:
    selected=data[(data.donor==d)&(data.celltype==fam)&data.pool.isin(ps)&data.arm.isin([ligand,'Control'])]
    for m,mb in selected.groupby(program):
     ag=mb.assign(weighted=mb[value]*mb.n_cells).groupby(['arm','substate']).agg(total=('weighted','sum'),n=('n_cells','sum'))
     z=(ag.total/ag.n).unstack(0).reindex(index=st.index,columns=[ligand,'Control']);assert z.notna().all().all()
     unstd=(z[ligand]*st[ligand]/nt).sum()-(z.Control*st.Control/nc).sum()
     if ok:
      zs=z.loc[support.index];w=support.Control/support.Control.sum();standard=(zs[ligand]*w).sum()-(zs.Control*w).sum();common=(zs[ligand]*support[ligand]/support[ligand].sum()).sum()-(zs.Control*w).sum()
     else:standard=common=np.nan
     row=dict(ligand=ligand,donor=d,celltype=fam,programme=m,pools=';'.join(ps),unstandardized_full_delta=unstd,unstandardized_common_delta=common,standardized_common_delta=standard,n_common_substates=len(support),common_substates=';'.join(support.index),treatment_retained_fraction=tcov,control_retained_fraction=ccov,status=reason,annotation_resolution='MULTI_SUBSTATE' if finer else 'NO_FINER_RELEASED_SUBSTATE',units='saved_publisher_cNMF_usage' if is_marker else 'mean_cell_ln1p_10000_library_normalized_expression')
     (markers if is_marker else effects).append(row)
save(pd.DataFrame(props),'results/recovered_family_proportions.tsv');save(pd.DataFrame(within),'results/within_family_substate_proportions.tsv');save(pd.DataFrame(effects),'results/composition_standardized_expression.tsv');save(pd.DataFrame(markers),'results/composition_marker_summary.tsv')
print('COMPOSITION_COMPLETE',len(effects),'normalized-expression summaries;',len(markers),'published programme summaries',flush=True)
