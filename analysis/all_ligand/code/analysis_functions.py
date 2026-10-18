import numpy as np
import pandas as pd
UP='ICI_COLITIS_EPITHELIAL_UP_STRICT'
DEP='ICI_COLITIS_EPITHELIAL_LOSS_STRICT'
IFN='CURATED_IFN_VISIBILITY'
RES='ICI_COLITIS_EPITHELIAL_UP_NON_IFN'
IBD1='IBD1_MEMBERSHIP_RANK';IBD2='IBD2_MEMBERSHIP_RANK'
def bh_planned(values,n):
    p=np.asarray(values,dtype=float);valid=np.isfinite(p);v=p[valid]
    assert len(p)<=n and ((v>=0)&(v<=1)).all()
    out=np.full(p.shape,np.nan)
    if len(v):
        order=np.argsort(v,kind='stable');s=v[order];q=np.minimum.accumulate((s*n/np.arange(1,len(s)+1))[::-1])[::-1].clip(0,1)
        q0=np.empty(len(v));q0[order]=q;out[valid]=q0
    return out
def direction_pattern(up,dep,required=3):
    up=np.asarray(up);dep=np.asarray(dep)
    assert up.shape==dep.shape
    if len(up)!=required or not (np.isfinite(up).all() and np.isfinite(dep).all()):return 'NE'
    a=(up<0).all();b=(dep>0).all()
    return 'BOTH' if a and b else 'UP_ONLY' if a else 'DEPLETED_ONLY' if b else 'OTHER_OBSERVED'
def reference_direction(deltas):
    d=np.asarray(deltas)
    if len(d)!=3 or not np.isfinite(d).all():return 'NE'
    return 'ALL_NEGATIVE' if (d<0).all() else 'ALL_POSITIVE' if (d>0).all() else 'ALL_ZERO' if (d==0).all() else 'MIXED_OR_ZERO'
def paired_deltas(scores):
    keys=['ligand','celltype','donor','module']
    assert not scores.duplicated(keys+['treatment']).any()
    t=scores[scores.treatment=='Treatment'].set_index(keys);c=scores[scores.treatment=='Control'].set_index(keys)
    assert set(t.index)==set(c.index)
    out=t[['raw_mean_rank','pools','n_cells','value_origin']].rename(columns={'raw_mean_rank':'treatment_score','n_cells':'treatment_cells'}).join(c[['raw_mean_rank','pools','n_cells']].rename(columns={'raw_mean_rank':'control_score','pools':'control_pools','n_cells':'control_cells'}))
    assert (out.pools==out.control_pools).all()
    out['delta']=out.treatment_score-out.control_score
    return out.reset_index()
