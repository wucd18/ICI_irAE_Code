import os
from pathlib import Path
import csv,io,json,hashlib,zipfile,re,sys,datetime,shutil
ROOT=Path(__file__).resolve().parents[1]
C=ROOT/'05_Rebuild';A=ROOT/'04_Audit';U=ROOT/'01_Journal_Upload';R=ROOT/'02_Review';F=ROOT/'03_Figshare_Update'
INP=A/'inputs_readonly';ZIN=Path(os.environ['ICI_FINAL_WORKBOOK_INPUT_ROOT']).resolve();SIN=INP/'source_zip';SOUT=C/'source_package_working'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def bh(b):return hashlib.sha256(b).hexdigest()
def put(p,s):
 p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
 p.write_bytes(s) if isinstance(s,bytes) else p.write_text(s,encoding='utf8')
def jw(p,x):put(p,json.dumps(x,ensure_ascii=False,indent=2)+'\n')
def jr(p):return json.loads(Path(p).read_text(encoding='utf-8-sig'))
def cw(p,rows):
 rows=list(rows);p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
 with p.open('w',encoding='utf-8-sig',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def cr(p):
 with Path(p).open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
