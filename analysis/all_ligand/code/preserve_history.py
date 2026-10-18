from common import *
import zipfile,xml.etree.ElementTree as ET,shutil
src=ROOT/'inputs/Manuscript_CSBJ_revised_clean.docx';z=zipfile.ZipFile(src)
ns={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main','a':'http://schemas.openxmlformats.org/drawingml/2006/main','r':'http://schemas.openxmlformats.org/officeDocument/2006/relationships'}
tree=ET.fromstring(z.read('word/document.xml'));paras=tree.findall('.//w:p',ns)
txt=lambda p:''.join(t.text or '' for t in p.findall('.//w:t',ns))
idx=next(i for i,p in enumerate(paras) if txt(p).startswith('Figure 4. TNFSF12'))
blip=paras[idx-1].find('.//a:blip',ns);rid=blip.attrib['{'+ns['r']+'}embed'];rels=ET.fromstring(z.read('word/_rels/document.xml.rels'));target=next(r.attrib['Target'] for r in rels if r.attrib['Id']==rid);media='word/'+target
dest=ROOT/'inputs/historical_Figure4';dest.mkdir(exist_ok=True)
(dest/('Figure4_embedded_original'+Path(media).suffix)).write_bytes(z.read(media))
(dest/'original_Figure4_legend_and_context.txt').write_text('\n\n'.join(txt(p) for p in paras[idx-9:idx+2]),encoding='utf8')
jw(dest/'media_provenance.json',dict(source=str(src),source_sha256=sha(src),relationship_id=rid,zip_member=media,sha256=hashlib.sha256(z.read(media)).hexdigest(),extraction='Unmodified embedded media from the current CSBJ manuscript, including historical panel c',instruction='Retain selected-set transfer and original selection caveats; proposed supplementary move only'))
shutil.copy2(CODEBASE/'analysis/organoid/analyze_organoid.R',ROOT/'inputs/approved_original_analyze_organoid.R')
shutil.copy2(SAVED/'sessionInfo.txt',ROOT/'inputs/targeted_sessionInfo.txt')
shutil.copy2(OLD/'provenance/Capeling_Figure_4_Code_pinned.Rmd',ROOT/'inputs/publisher_Figure4_code_pinned.Rmd')
print('HISTORICAL_MEDIA_PRESERVED',media,flush=True)
