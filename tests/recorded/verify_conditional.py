from pathlib import Path
import hashlib,json
import numpy as np,pandas as pd
from scipy.stats import t
O=Path(__file__).resolve().parents[1]
d=pd.read_csv(O/'source_data/conditional_model_input.tsv',sep='\t')
r=pd.read_csv(O/'source_data/conditional_results.tsv',sep='\t')
l=pd.read_csv(O/'source_data/conditional_leave_patient_out.tsv',sep='\t')
checks=[]
def fit(y,X):
 beta=np.linalg.lstsq(X,y,rcond=None)[0];inv=np.linalg.inv(X.T@X);h=np.sum((X@inv)*X,axis=1);u=y-X@beta
 v=inv@(X.T@((u/(1-h))[:,None]**2*X))@inv
 return beta,np.sqrt(np.diag(v))
for z in r.itertuples():
 b=d[(d.dataset==z.dataset)&(d.celltype==z.celltype)];case='irColitis' if z.dataset=='GSE206300' else 'CPI_colitis'
 x=b.CURATED_IFN_VISIBILITY.to_numpy();X=np.column_stack([np.ones(len(b)),b.condition.eq(case).astype(float),x-x.mean()]);y=b[z.module].to_numpy()
 beta,se=fit(y,X);p=2*t.sf(abs(beta[1]/se[1]),len(b)-3)
 for key,val,original in [('beta',beta[1],z.adjusted_beta),('se',se[1],z.adjusted_se_HC3),('p',p,z.p)]:
  assert np.isclose(val,original,rtol=1e-9,atol=1e-12),(z.dataset,z.celltype,z.module,key,val,original)
  checks.append(dict(dataset=z.dataset,celltype=z.celltype,module=z.module,quantity=key,python=val,R=original,status='PASS'))
 for i,patient in enumerate(b.patient_id):
  bb,ss=fit(np.delete(y,i),np.delete(X,i,axis=0));v=l[(l.dataset==z.dataset)&(l.celltype==z.celltype)&(l.module==z.module)&(l.omitted_patient==patient)].adjusted_beta
  assert len(v)==1 and np.isclose(bb[1],v.iloc[0],rtol=1e-9,atol=1e-12)
for ds in r.dataset.unique():
 ix=r.dataset.eq(ds);p=r.loc[ix,'p'].to_numpy();order=np.argsort(p);q=np.minimum.accumulate((p[order]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1];expected=np.empty_like(q);expected[order]=np.minimum(q,1)
 assert np.allclose(r.loc[ix,'fdr'],expected,rtol=1e-12,atol=1e-12)
pd.DataFrame(checks).to_csv(O/'review/conditional_independent_numerical_checks.tsv',sep='\t',index=False)
summary=dict(status='PASS',R_estimates_checked=len(checks),leave_patient_out_coefficients_checked=len(l),primary_tests=int(r.dataset.eq('GSE206300').sum()),primary_FDR05=int(((r.dataset=='GSE206300')&(r.fdr<.05)).sum()),primary_nonoverlap_families=int(r[(r.dataset=='GSE206300')&(r.IFN_range_overlap==0)].celltype.nunique()))
(O/'review/conditional_verification.json').write_text(json.dumps(summary,indent=2),encoding='utf-8');print(summary)
