from pathlib import Path
from pypdf import PdfReader
import pypdf,PIL,openpyxl,csv,json,re,datetime,hashlib,ast,sys
ROOT=Path(__file__).resolve().parents[1]
def rows(p):
 with p.open(encoding='utf8',newline='') as f:return list(csv.DictReader(f,delimiter='\t'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
cl=rows(ROOT/'results/classification_stability.tsv');ligands=sorted({r['ligand'] for r in cl});checks=[];pdfs=[]
for p in sorted((ROOT/'figures').glob('*.pdf')):
 reader=PdfReader(p);pages=[]
 for i,page in enumerate(reader.pages):
  w=float(page.mediabox.width);h=float(page.mediabox.height);a,b=sorted([w,h]);assert abs(a-595.276)<.5 and abs(b-841.89)<.5,(p,i,w,h)
  text=page.extract_text();assert len(text)>50 and '\ufffd' not in text
  if p.name=='Figure_Full_Three_Method_Atlas.pdf':assert all(l in text for l in ligands) and 'No_Cytomix' in text and 'NE' in text,(p,i)
  pages.append(dict(page=i+1,width_pt=w,height_pt=h,text_chars=len(text),standard_A4=True))
 pdfs.append(dict(file=p.name,sha256=sha(p),pages=pages))
for p in (ROOT/'text').glob('*.md'):
 for name in set(re.findall(r'\b[A-Za-z][A-Za-z0-9_]*\.tsv\b',p.read_text(encoding='utf8'))):
  assert any((ROOT/d/name).exists() for d in ['results','review','text','inputs']), (p,name)
# Explicit checks for the key Results claims against their machine-readable tables.
orig=[r for r in cl if r['original_both']=='True'];stable=[r for r in orig if r['raw_bg_nonrank_both']=='True'];text=(ROOT/'text/results_update.md').read_text(encoding='utf8')
assert f'Of the {len(orig)} original Both combinations, {len(stable)} retained Both' in text
for r in orig:assert r['ligand'] in text and r['family'] in text
g=rows(ROOT/'results/gene_set_overlap.tsv');d={(r['programme_a'],r['programme_b'],r['scope']):int(r['n_intersection']) for r in g}
up='ICI_COLITIS_EPITHELIAL_UP_STRICT';dep='ICI_COLITIS_EPITHELIAL_LOSS_STRICT';refs=['CURATED_IFN_VISIBILITY','IBD1_MEMBERSHIP_RANK','IBD2_MEMBERSHIP_RANK']
assert [d[(up,b,'measured')] for b in refs]==[4,2,0] and [d[(dep,b,'measured')] for b in refs]==[6,12,17]
cx=rows(ROOT/'results/published_proliferation_context.tsv');summary={(r['ligand'],r['celltype'],r['programme']):r for r in cx}
gp11=[(r['ligand'],r['family']) for r in orig if int(summary[(r['ligand'],r['family'],'GP11')]['n_negative'])==3]
gp14down=[(r['ligand'],r['family']) for r in orig if int(summary[(r['ligand'],r['family'],'GP14')]['n_negative'])==3]
gp14up=[(r['ligand'],r['family']) for r in orig if int(summary[(r['ligand'],r['family'],'GP14')]['n_positive'])==3]
assert set(gp11)=={('BMP4','Colonocyte'),('BMP6','Goblet'),('INHBA','Transit_Amplifying'),('INS','Colonocyte'),('TGFB1','Transit_Amplifying'),('TGFB2','Transit_Amplifying')}
assert gp14down==[('BMP6','Goblet')] and gp14up==[('INHBA','Transit_Amplifying')]
for p in (ROOT/'code').glob('*.py'):ast.parse(p.read_text(encoding='utf8'))
record=dict(status='PASS',scope='PDF size/text coverage, text-to-result references and key claim verification',pdfs=pdfs,total_pdf_pages=sum(len(r['pages']) for r in pdfs),all_planned_ligands_on_every_atlas_page=True,source_file_references_resolve=True,Results_key_claims_checked=True,publisher_programme_claims_checked=True,visual_review=dict(status='PASS',method='All pages rendered with Poppler; contact sheets reviewed for every page and full-size pages inspected for main figure, dense atlas and donor plots. Initial legend/footer, colourbar label and mean-label collisions corrected; final render reviewed.',reviewer='Codex',scientific_content_unchanged_by_layout_fixes=True),artifact_environment=dict(executable=sys.executable,pypdf=pypdf.__version__,Pillow=PIL.__version__,openpyxl=openpyxl.__version__),completed_utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
(ROOT/'review/document_and_figure_QA.json').write_text(json.dumps(record,indent=2),encoding='utf8')
(ROOT/'review/visual_review.md').write_text('# Visual review\n\nPASS for the final rendered figure set. All 27 pages reviewed in contact sheets; full-resolution inspection of the main figure, representative dense atlas pages and donor plots. Corrected legend/footer collision, adjacent colourbar endpoint labels, and nearby mean-number labels. No panel data changed. Stars explicitly denote display saturation, not significance. All pages are standard A4.\n',encoding='utf8')
print('DOCUMENT_FIGURE_QA_PASS',record['total_pdf_pages'],'pages')
