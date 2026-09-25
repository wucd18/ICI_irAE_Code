from pathlib import Path
import io,struct,json,tempfile
import numpy as np,pandas as pd
from scipy.sparse import csc_matrix
from stream_r_object import Reader,Node
from export_r_capture import unpack
O=Path(__file__).resolve().parents[1];D=O/'review/serialization_fixture_v3'
a,obs,genes,shape=unpack(D/'parsed')
m=csc_matrix((a['x'],a['i'].astype(int),a['p'].astype(int)),shape=shape)
exp=pd.read_csv(D/'fixture_expected_counts.tsv',sep='\t',index_col=0)
assert np.array_equal(m.toarray(),exp.to_numpy())
assert genes==list(exp.index) and list(obs.index)==list(exp.columns)
em=pd.read_csv(D/'fixture_expected_metadata.tsv',sep='\t',index_col=0)
pd.testing.assert_frame_equal(obs,em,check_dtype=False,check_index_type=False,check_names=False,rtol=1e-13)
# Exercise the production streaming threshold with exact big-endian doubles.
v=np.arange(600001,dtype='>f8')
b=struct.pack('>ii',14,len(v))+v.tobytes()
r=Reader(io.BytesIO(b),D/'large_vector_test')
n=r.item('/assays/0/counts/x')
assert np.array_equal(np.memmap(n.v['file'],dtype=n.v['dtype'],mode='r'),v)
print('PASS: every fixture count, gene, cell, metadata value, factor and large numeric payload')
(D/'fixture_verification.json').write_text(json.dumps(dict(status='PASS',counts=int(m.shape[0]*m.shape[1]),metadata=int(obs.size),large_vector=len(v)),indent=2),encoding='utf-8')
