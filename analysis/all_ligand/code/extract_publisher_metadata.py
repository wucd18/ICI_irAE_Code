from pathlib import Path
import zipfile,io,xml.etree.ElementTree as ET,json,csv,hashlib
import openpyxl
ROOT=Path(__file__).resolve().parents[1]
src=Path(r'E:\ICI_irAE_Original_Project\01_data\GSE313368_article_supplementary.zip')
z=zipfile.ZipFile(src)
xname=next(n for n in z.namelist() if n.endswith('MOESM5_ESM.xlsx'))
book=openpyxl.load_workbook(io.BytesIO(z.read(xname)),read_only=True,data_only=True)
sheet=book['B Genes in IBD1 and IBD2']; rows=list(sheet.values)
defs=[];duplicates=[]
for col,name in enumerate(['IBD1_MEMBERSHIP_RANK','IBD2_MEMBERSHIP_RANK']):
    members=[str(r[col]).strip() for r in rows[1:] if r[col] is not None]
    from collections import Counter
    duplicates.extend(dict(module=name,gene=g,occurrences=n,action='Count once for a membership set; no published weights in this sheet') for g,n in Counter(members).items() if n>1)
    members=list(dict.fromkeys(members))
    defs.extend(dict(module=name,gene=g,weight=1,definition='Published IBD gene membership; no weight supplied in this sheet',source=xname,sheet=sheet.title) for g in members)
with (ROOT/'inputs/IBD_membership_definitions.tsv').open('w',encoding='utf8',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(defs[0]),delimiter='\t');w.writeheader();w.writerows(defs)
(ROOT/'inputs/IBD_duplicate_members.json').write_text(json.dumps(duplicates,indent=2),encoding='utf8')
docname=next(n for n in z.namelist() if n.endswith('MOESM16_ESM.docx'))
doc=zipfile.ZipFile(io.BytesIO(z.read(docname)));ns={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main'}
tree=ET.fromstring(doc.read('word/document.xml'));tables=tree.findall('.//w:tbl',ns)
rr=[[''.join(c.itertext()) if False else ''.join(t.text or '' for t in c.findall('.//w:t',ns)) for c in r.findall('./w:tc',ns)] for r in tables[0].findall('./w:tr',ns)]
with (ROOT/'inputs/publisher_protein_doses_verbatim.json').open('w',encoding='utf8') as f:json.dump(rr,f,ensure_ascii=False,indent=2)
# Explicit source-name mappings only; no fuzzy matching or inferred concentrations.
aliases={'IFN-alpha':'IFNA','IFN-epsilon':'IFNE','IFN-gamma':'IFNG','TNF-alpha':'TNFA','TWEAK (TNFSF12)':'TNFSF12','TRAIL (TNFSF10)':'TNFSF10','RANKL (TNFSF11)':'TNFSF11'}
aliases.update({'CT-1':'CTF1','GH':'GH1','PGRN':'GRN','IFNɣ':'IFNG','Il13':'IL13','Activin A':'INHBA','Insulin':'INS','Metastin':'KISS1','SCF':'KITLG','Galectin-9':'LGALS9','PBEF':'NAMPT','NT-3':'NTF3','NT-4':'NTF4','OPN':'SPP1','TGFβ1':'TGFB1','TGFβ2':'TGFB2','TGFβ3':'TGFB3','TNFɑ':'TNFA','IFNɑ':'IFNA','IFNε':'IFNE','Il1ɑ':'IL1A'})
ligands=json.loads((ROOT/'analysis_config.json').read_text())['ligands'];dose=[]
for r in rr[2:]:
    raw=r[0].strip();gene=aliases.get(raw,raw)
    if gene not in ligands:
        matches=[g for g in ligands if '('+g+')' in raw]
        if len(matches)==1:gene=matches[0]
    dose.append(dict(ligand=gene,published_protein_name=raw,vendor=r[1],catalog_number=r[2],dose=r[3],source=docname,mapping='exact' if gene==raw else 'explicit alias',mapping_status='MATCHED' if gene in ligands else 'UNMATCHED'))
with (ROOT/'inputs/condition_metadata.tsv').open('w',encoding='utf8',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(dose[0]),delimiter='\t');w.writeheader();w.writerows(dose)
print('IBD_MEMBERS',[(n,sum(d['module']==n for d in defs)) for n in sorted({d['module'] for d in defs})])
print('DOSE_ROWS',len(dose),'UNMATCHED',[r for r in dose if r['mapping_status']!='MATCHED'])
print('MISSING',sorted(set(ligands)-{r['ligand'] for r in dose}))
