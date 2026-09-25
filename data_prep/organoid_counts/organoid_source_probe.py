from pathlib import Path
import requests,json,sys,ctypes
sys.stdout.reconfigure(encoding='utf-8')
O=Path(__file__).resolve().parents[1]
base='https://ftp.ncbi.nlm.nih.gov/geo/series/GSE313nnn/GSE313368/suppl/'
for f in ['filelist.txt']:
 r=requests.get(base+f,timeout=30);r.raise_for_status();(O/'provenance'/f).write_bytes(r.content);print(r.text)
for f in ['Cytomix_Secretome_Processing.Rmd','README.md']:
 u='https://raw.githubusercontent.com/Genentech/secretome_dictionary/3c80f66552d9c3ea344ba02304eb588b717e9ed8/'+f
 r=requests.get(u,timeout=30);r.raise_for_status();(O/'provenance'/f).write_bytes(r.content)
 print(f,len(r.content))
for f in ['GSE313368_RAW.tar','GSE313368_Combined_Pool_Merge_Harmonized_Donor_Only.Robj.gz']:
 r=requests.get(base+f,headers={'Range':'bytes=0-2047','Accept-Encoding':'identity'},timeout=45);print(f,r.status_code,dict(r.headers));(O/'provenance'/(f+'.head.bin')).write_bytes(r.content[:2048]);print(repr(r.content[:512]))
class Mem(ctypes.Structure):
 _fields_=[('length',ctypes.c_ulong),('load',ctypes.c_ulong)]+[(x,ctypes.c_ulonglong) for x in ['total','avail','page_total','page_avail','virtual_total','virtual_avail','extended_avail']]
m=Mem();m.length=ctypes.sizeof(m);ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(m));print('RAM',m.total,m.avail)
print('SOURCE_METHOD_EXCERPTS')
from bs4 import BeautifulSoup
soup=BeautifulSoup((O/'provenance/Capeling_publisher.html').read_text(encoding='utf-8'),'html.parser')
for p in soup.find_all(['p']):
 t=p.get_text(' ',strip=True)
 if any(x in t.lower() for x in ['pseudobulk','deseq','donor-pooled','control wells','raw counts']):print(t)
