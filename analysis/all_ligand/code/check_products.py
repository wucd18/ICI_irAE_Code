from common import *
from analysis_functions import *
from scipy.stats import false_discovery_control
import py_compile,collections
cfg=jr(ROOT/'analysis_config.json');g=pd.read_csv(ROOT/'results/all_ligand_programme_responses.tsv',sep='\t',float_precision='round_trip');s=pd.read_csv(ROOT/'results/response_pattern_summary.tsv',sep='\t',float_precision='round_trip');d=pd.read_csv(ROOT/'results/all_ligand_donor_deltas.tsv',sep='\t',float_precision='round_trip');bench=pd.read_csv(ROOT/'results/programme_benchmark.tsv',sep='\t')
for p in (ROOT/'code').glob('*.py'):py_compile.compile(str(p),doraise=True)
for scope,b in g.groupby('FDR_scope'):
    pp=b.pval.to_numpy();known=np.isfinite(pp);oracle=false_discovery_control(np.where(known,pp,1),method='bh');assert np.allclose(oracle[known],b.FDR_atlas.to_numpy()[known],atol=2e-15)
    assert b.loc[~known,'FDR_atlas'].isna().all()
assert not g.duplicated(['ligand','celltype','module']).any() and not d.duplicated(['ligand','celltype','donor']).any()
assert bench.P_value.isna().all()
for r in s.itertuples():
    b=d[(d.ligand==r.ligand)&(d.celltype==r.celltype)&d.paired_score_observed]
    assert direction_pattern(b.delta_up,b.delta_depleted)==r.pattern
    assert len(b)==r.n_donors
    for n in ['up','depleted','IFN','residual']:
        assert np.isclose(b['delta_'+n].mean(),getattr(r,'delta_'+n+'_mean'),rtol=0,atol=1e-16,equal_nan=True)
    gg=g[(g.ligand==r.ligand)&(g.celltype==r.celltype)]
    for col,m in [('delta_up',UP),('delta_depleted',DEP),('delta_IFN',IFN),('delta_residual',RES)]:
        z=gg[gg.module==m].iloc[0];assert np.isclose(z.delta_mean,b[col].mean(),atol=1e-16,rtol=0,equal_nan=True)
    assert not ((b.delta_up==0)&(b.up_negative==1)).any()
assert len(g[g.condition_class=='ligand'])==len(cfg['ligands'])*len(cfg['families'])*len(cfg['modules_in_original_fgsea_order'])
# Data used by Figure 4 are exact copies of the final machine-readable tables.
figa=pd.read_csv(ROOT/'figure_source_data/Figure4_TNFSF12_a.tsv',sep='\t',float_precision='round_trip');fdat=g[(g.ligand=='TNFSF12')&g.module.isin(cfg['primary_modules'])]
for c in ['NES','pval','FDR_atlas','old_targeted_FDR']:assert np.array_equal(figa[c],fdat[c],equal_nan=True)
assert len(list((ROOT/'figures').glob('*.pdf')))==6
nr=pd.read_csv(ROOT/'review/new_model_checks_and_provenance.tsv',sep='\t');assert len(nr)==len(pd.read_csv(ROOT/'inputs/new_model_jobs.tsv',sep='\t')) and nr.status.eq('PASS').all()
# Validate source immutability without reproducing prior analyses.
sourcecheck=[]
for r in pd.read_csv(ROOT/'review/input_manifest.tsv',sep='\t').itertuples():
    actual=sha(r.path);sourcecheck.append(dict(path=r.path,unchanged=actual==r.sha256,sha256=actual));assert actual==r.sha256,r.path
save(pd.DataFrame(sourcecheck),'review/input_immutability_check.tsv')
status=dict(status='PASS',scope='New outputs only; independent BH check, all donor pattern/mean checks, reused estimate checks, model/schema checks, figure source consistency, unchanged input hashes, Python static compilation',primary_estimable=int(g[(g.FDR_scope=='all_ligand_primary_948')].pval.notna().sum()),new_models_checked=len(nr),old_models_refitted=0,original_files_unchanged=True,created_utc=now())
jw(ROOT/'review/product_numerical_checks.json',status);print('PRODUCT_CHECKS_PASS',flush=True)
