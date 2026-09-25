"""Acquire primary metadata only, retaining status, content type and hashes."""
from pathlib import Path
import requests,json,hashlib
from datetime import datetime,timezone
O=Path(__file__).resolve().parents[1]
jobs=[('zenodo16948925.json','https://zenodo.org/api/records/16948925'),('jci202488.html','https://www.jci.org/articles/view/202488'),('UNICIT_tree.json','https://api.github.com/repos/mickvaneijs/UNICIT/git/trees/main?recursive=1')]
out=[]
for name,url in jobs:
 rec={'name':name,'url':url,'start_utc':datetime.now(timezone.utc).isoformat()}
 try:
  q=requests.get(url,timeout=(20,90));q.raise_for_status();p=O/'new_data'/name
  if p.exists():assert p.read_bytes()==q.content,'Existing immutable source differs'
  else:p.write_bytes(q.content)
  rec.update(status='PASS',bytes=len(q.content),sha256=hashlib.sha256(q.content).hexdigest(),content_type=q.headers.get('Content-Type'))
 except Exception as e:rec.update(status='FAIL',error=repr(e))
 out.append(rec);print(name,rec,flush=True)
(O/'review/public_cohort_probe.json').write_text(json.dumps(out,indent=2))
raise SystemExit(0 if all(r['status']=='PASS' for r in out) else 1)
