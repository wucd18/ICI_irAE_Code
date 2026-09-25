from pathlib import Path
import sys
sys.stdout.reconfigure(encoding='utf-8')
import requests, json, re, hashlib, time
from datetime import datetime,timezone
O=Path(__file__).resolve().parents[1]
urls={
 'GSE313368_soft.txt':'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE313368&targ=self&form=text&view=full',
 'GSE313368_supplement_listing.html':'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE313nnn/GSE313368/suppl/',
 'secretome_repository_tree.json':'https://api.github.com/repos/Genentech/secretome_dictionary/git/trees/3c80f66552d9c3ea344ba02304eb588b717e9ed8?recursive=1',
 'Capeling_publisher.html':'https://www.nature.com/articles/s41467-025-68247-6'
}
rows=[]
for name,url in urls.items():
 try:
  r=requests.get(url,timeout=45);r.raise_for_status();p=O/'provenance'/name;p.write_bytes(r.content)
  rows.append(dict(file=name,url=url,status=r.status_code,size=len(r.content),sha256=hashlib.sha256(r.content).hexdigest(),utc=datetime.now(timezone.utc).isoformat()))
  if name.endswith('soft.txt'):print('\n'.join(x for x in r.text.splitlines() if 'supplement' in x.lower() or 'Series_' in x))
  elif name.endswith('listing.html'):print(r.text)
  elif name.endswith('tree.json'):print('\n'.join(x['path'] for x in r.json()['tree']))
  else:print(name,len(r.content))
 except Exception as e:rows.append(dict(file=name,url=url,error=repr(e)));print(name,repr(e))
(O/'provenance/acquisition_metadata_manifest.json').write_text(json.dumps(rows,indent=2),encoding='utf-8')
