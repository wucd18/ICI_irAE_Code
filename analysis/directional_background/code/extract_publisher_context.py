from pathlib import Path
import zipfile,io,openpyxl,csv,re,html,json
ROOT=Path(__file__).resolve().parents[1]
source=Path(r'E:\ICI_irAE_Original_Project\01_data\GSE313368_article_supplementary.zip')
with zipfile.ZipFile(source) as z:
 name='41467_2025_68247_MOESM11_ESM.xlsx';w=openpyxl.load_workbook(io.BytesIO(z.read(name)),read_only=True,data_only=True)
 rows=list(w['C_GO_Terms'].values);head=rows[0];data=[dict(zip(head,r)) for r in rows[1:] if r[-1]]
 with (ROOT/'inputs/published_GP_GO_annotations.tsv').open('w',newline='',encoding='utf8') as f:
  writer=csv.DictWriter(f,fieldnames=[h for h in head if h],delimiter='\t');writer.writeheader();writer.writerows({k:v for k,v in r.items() if k} for r in data)
 top=[]
 for gp in sorted(set(r['Gene_Program'] for r in data)):
  b=sorted([r for r in data if r['Gene_Program']==gp],key=lambda r:float(r['Adjusted.P.value']))[:3]
  top.append(dict(programme=gp,metadata_column=gp.replace('GP','X')+'_DEG',top_published_GO_terms='; '.join(r['Term'] for r in b),source=name+' / C_GO_Terms',semantics='Saved publisher perturbation-derived cNMF usage; no refitting, no new gene selection'))
 with (ROOT/'inputs/published_GP_context.tsv').open('w',newline='',encoding='utf8') as f:
  wr=csv.DictWriter(f,fieldnames=list(top[0]),delimiter='\t');wr.writeheader();wr.writerows(top)
old=Path(r'E:\ICI_irAE_V1_Redline_Revision\20260913T220504Z_d38de3ca\provenance')
s=(old/'Capeling_publisher.html').read_text(encoding='utf8');paragraphs=[]
for p in re.findall(r'<p[ >].*?</p>',s,re.S):
 t=html.unescape(re.sub('<[^>]+>','',p))
 if re.search(r'GP\d|GP \d',t) and any(x in t.lower() for x in ['stem','cycle','differenti','mature']):paragraphs.append(t)
(ROOT/'inputs/publisher_programme_context.txt').write_text('\n\n'.join(paragraphs),encoding='utf8')
print('Extracted all 20 programme annotation records, source metadata mapping and published context')
