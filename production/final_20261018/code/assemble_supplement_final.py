from pathlib import Path
import os
import io,json,shutil
from reportlab.pdfgen import canvas
from reportlab.lib.styles import ParagraphStyle
from reportlab.platypus import Paragraph
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from pypdf import PdfReader,PdfWriter,Transformation
P=Path(os.environ['ICI_PRODUCTION_ROOT']);O=Path(os.environ['ICI_SUBMISSION_ROOT']);pt=72/25.4
pdfmetrics.registerFont(TTFont('Arial',os.environ.get('ICI_ARIAL_REGULAR','C:/Windows/Fonts/arial.ttf')));pdfmetrics.registerFont(TTFont('ArialBold',os.environ.get('ICI_ARIAL_BOLD','C:/Windows/Fonts/arialbd.ttf')))
style=ParagraphStyle('legend',fontName='Arial',fontSize=9,leading=12,spaceAfter=5)
titleStyle=ParagraphStyle('title',fontName='ArialBold',fontSize=11,leading=14)
spec=json.loads((P/'inputs/supplement_pages.json').read_text(encoding='utf8'));writer=PdfWriter();records=[]
def text(c,s,x,y,w,sty=style):
 p=Paragraph(s,sty);_,h=p.wrap(w,1000);p.drawOn(c,x,y-h);return h
buf=io.BytesIO();c=canvas.Canvas(buf,pagesize=(210*pt,297*pt));x=17*pt;y=278*pt;w=176*pt
y-=text(c,'Supplementary figures and tables',x,y,w,titleStyle)+14
y-=text(c,'Mapping recurrent epithelial programmes of checkpoint inhibitor colitis to human organoid ligand responses',x,y,w,titleStyle)+18
intro=[
'Figures S1-S7 provide the complementary patient, spatial-assay and organoid results. Tables S1-S9 are supplied as nine editable result worksheets in Supplementary_Tables.xlsx.',
'Patients, released tissue samples and organoid donors are the biological replicates. Cells, genes, spatial spots and repeated sections are nested measurements. Three matched donors support eligible organoid contrasts.',
'Patient disease-oriented depleted scores equal 1 minus the raw rank. Organoid changes are treatment minus the matched cytomix control, with depleted changes retaining their raw direction. Disease-up decreases and depleted increases oppose the patient disease directions.',
'Raw rank, matched-background-relative rank and non-rank expression changes describe different representations of the same experiments. Non-rank changes use mean member log2(size-factor-normalized count + 1); composition standardization uses the distinct mean-cell ln(1 + 10,000 counts/library) scale.',
'Grey NE marks a non-estimable comparison, not zero. All planned ligand-family comparisons and necessary negative, discordant and unavailable results are retained in the result tables. GSEA q values retain the existing test families and are not assigned to donor-score or composition points.',
]
for t in intro:y-=text(c,t,x,y,w)+10
y-=8
for i in range(1,8):
 p=next(q for q in spec if q['figure']==f'S{i}');start=2+next(j for j,q in enumerate(spec) if q['figure']==f'S{i}');n=sum(q['figure']==f'S{i}' for q in spec)
 y-=text(c,f'<b>Figure S{i}.</b> {p["title"]} (pages {start}-{start+n-1}).',x,y,w)+5
c.setFont('Arial',8);c.drawCentredString(105*pt,12*pt,'1');c.showPage();c.save();buf.seek(0);writer.add_page(PdfReader(buf).pages[0])
for i,s in enumerate(spec,2):
 W,H=(297,210) if s['landscape'] else (210,297);buf=io.BytesIO();c=canvas.Canvas(buf,pagesize=(W*pt,H*pt))
 heading='Figure '+s['figure']+(' (continued)' if s['part']>1 else '')
 c.setFont('ArialBold',11);c.drawString(12.5*pt,(H-13)*pt,heading)
 text(c,s['title'],12.5*pt,(H-18)*pt,(W-25)*pt,ParagraphStyle('small',fontName='ArialBold',fontSize=9,leading=11))
 bottom=max(y+h for pid,x,y,w,h in s['panels']);caption_top=(H-bottom-5)*pt
 used=text(c,s['legend'],12.5*pt,caption_top,(W-25)*pt)
 assert caption_top-used>18*pt,(s['figure'],s['part'],caption_top/pt,used/pt)
 c.setFont('Arial',8);c.drawCentredString(W*pt/2,12*pt,str(i));c.showPage();c.save();buf.seek(0)
 page=writer.add_page(PdfReader(buf).pages[0])
 for pid,x,y,w,h in s['panels']:
  src=PdfReader(P/'panels'/f'{pid}.pdf').pages[0]
  assert abs(float(src.mediabox.width)-w*pt)<.1
  page.merge_transformed_page(src,Transformation().translate(x*pt,(H-y-h)*pt))
 records.append({'page':i,**s,'caption_height_mm':used/pt})
writer.write(O/'Supplementary_Figures.pdf')
(P/'qa/supplement_page_map.json').write_text(json.dumps(records,indent=2))
for i in range(1,7):shutil.copy2(P/'main_figures'/f'Figure_{i}.pdf',O/'Figures'/f'Figure_{i}.pdf')
print('Supplementary figures:',len(writer.pages),'clean vector pages; six main PDFs copied.')
