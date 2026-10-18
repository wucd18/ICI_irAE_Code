from common import *
import zipfile,csv,io,ast
assert json.loads((ROOT/'review/targeted_validation.json').read_text())['status']=='PASS'
doc=json.loads((ROOT/'review/document_and_figure_QA.json').read_text());assert doc['status']=='PASS'
cl=read(ROOT/'results/classification_stability.tsv');limitations=read(ROOT/'review/module_status_and_limitations.tsv')
jw(ROOT/'review/FINAL_STATUS.json',dict(status='COMPLETE_WITH_DISCLOSED_LIMITATIONS',numerical_validation='PASS',document_figure_QA='PASS',planned_combinations=len(cl),eligible=int(cl.eligible.sum()),NE=int((~cl.eligible).sum()),Both_counts={me:int((cl[me+'_class']=='BOTH').sum()) for me in ['raw','bg_relative','nonrank']},original_Both_retained_in_all_three=int(cl.raw_bg_nonrank_both.sum()),total_PDF_pages=doc['total_pdf_pages'],limitations=limitations.to_dict('records'),original_manuscript_changed=False,original_DESeq2_or_GSEA_rerun=False,new_inferential_P_values=False,public_upload=False,package_verification='See sibling ZIP.verification.json after archive creation',finished_utc=now()))
# Complete new-run deliverable, including newly derived caches, without copying original raw source blobs.
files=sorted([p for p in ROOT.rglob('*') if p.is_file() and '__pycache__' not in p.parts and p.name!='SHA256_manifest.tsv'],key=lambda p:p.relative_to(ROOT).as_posix())
manifest=[]
for p in files:manifest.append(dict(relative_path=p.relative_to(ROOT).as_posix(),bytes=p.stat().st_size,sha256=sha(p)))
save(pd.DataFrame(manifest),'SHA256_manifest.tsv');files.append(ROOT/'SHA256_manifest.tsv')
dest=ROOT.parents[1]/(ROOT.name+'_DELIVERY.zip');assert not dest.exists(),dest
print('BUILDING_ZIP',str(dest),'files',len(files),flush=True)
with zipfile.ZipFile(dest,'x',compression=zipfile.ZIP_DEFLATED,compresslevel=4,allowZip64=True) as z:
 for i,p in enumerate(files):
  z.write(p,ROOT.name+'/'+p.relative_to(ROOT).as_posix())
  if i%50==0:print('ARCHIVED',i,'/',len(files),flush=True)
print('VERIFYING_ARCHIVE',flush=True)
with zipfile.ZipFile(dest) as z:
 assert z.testzip() is None
 inside=pd.read_csv(io.BytesIO(z.read(ROOT.name+'/SHA256_manifest.tsv')),sep='\t')
 for i,r in enumerate(inside.itertuples()):
  h=hashlib.sha256();n=0
  with z.open(ROOT.name+'/'+r.relative_path) as f:
   while b:=f.read(8*1024*1024):h.update(b);n+=len(b)
  assert h.hexdigest()==r.sha256 and n==r.bytes,r.relative_path
 assert len(z.namelist())==len(inside)+1
receipt=dict(status='PASS',zip=str(dest),bytes=dest.stat().st_size,sha256=sha(dest),members=len(files),manifest_verified_members=len(manifest),CRC='PASS',all_member_SHA256='PASS',completed=now(),scientific_release_status='COMPLETE_WITH_DISCLOSED_LIMITATIONS')
jw(Path(str(dest)+'.verification.json'),receipt);Path(str(dest)+'.sha256').write_text(receipt['sha256']+'  '+dest.name+'\n',encoding='ascii');print(json.dumps(receipt),flush=True)
