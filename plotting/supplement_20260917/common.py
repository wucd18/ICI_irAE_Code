import os
from pathlib import Path
import csv,io,json,hashlib,zipfile,re,sys,datetime,shutil
ROOT=Path(__file__).resolve().parents[1]
C=ROOT/'03_Editorial_Code'
I=ROOT/'02_Internal_Review'
U=ROOT/'01_Journal_Upload'
SOURCE=Path(os.environ['ICI_SUPPLEMENT_SOURCE_ROOT']).resolve()
if not SOURCE.is_dir(): raise ValueError('BLOCKED_INPUT: ICI_SUPPLEMENT_SOURCE_ROOT')
PYTHON=Path(sys.executable)
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(p,b):
 p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
 p.write_bytes(b) if isinstance(b,bytes) else p.write_text(b,encoding='utf8')
def jw(p,x):write(p,json.dumps(x,ensure_ascii=False,indent=2))
def jr(p):return json.loads(Path(p).read_text(encoding='utf8'))
def cw(p,rows):
 rows=list(rows);p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
 with p.open('w',encoding='utf-8-sig',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def cr(p):
 with Path(p).open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
def src(n):
 matches=list((SOURCE/'full_source_tables').glob(f'Source_{n:03d}_*.csv'))
 if len(matches)!=1:raise ValueError(f'Ambiguous/missing Source_{n:03d}')
 return matches[0]
