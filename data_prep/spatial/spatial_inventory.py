from pathlib import Path
import re,json,csv,hashlib
O=Path(__file__).resolve().parents[1];p=O/'new_data/GSE189184_family.soft'
samples=[];d=None
for line in p.read_text(encoding='utf-8').splitlines():
 if line.startswith('^SAMPLE = '):
  if d:samples.append(d)
  d={'gsm':line.split(' = ',1)[1],'files':[],'characteristics':{}}
 elif d is not None:
  if line.startswith('!Sample_title = '):d['title']=line.split(' = ',1)[1]
  if line.startswith('!Sample_description = '):d['section']=line.split(' = ',1)[1]
  if line.startswith('!Sample_characteristics_ch1 = '):
   s=line.split(' = ',1)[1];k,_,v=s.partition(': ');d['characteristics'][k]=v
  if line.startswith('!Sample_supplementary_file') and ' = ' in line:d['files'].append(line.split(' = ',1)[1].replace('ftp://','https://'))
if d:samples.append(d)
assert len(samples)==len({s['gsm'] for s in samples})
print(json.dumps(samples,indent=2))
(O/'review/spatial_inventory.json').write_text(json.dumps(samples,indent=2),encoding='utf-8')
rows=[{'gsm':d['gsm'],'title':d.get('title'),'section':d.get('section'),**d['characteristics'],'files':';'.join(d['files'])} for d in samples]
with (O/'review/spatial_inventory.tsv').open('w',newline='',encoding='utf-8') as f:
 w=csv.DictWriter(f,fieldnames=list(dict.fromkeys(k for r in rows for k in r)),delimiter='\t');w.writeheader();w.writerows(rows)
