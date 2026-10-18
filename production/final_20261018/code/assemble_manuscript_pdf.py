from pathlib import Path
import os
import io,json,re
from pypdf import PdfReader,PdfWriter,Transformation
from reportlab.pdfgen import canvas
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
P=Path(os.environ['ICI_PRODUCTION_ROOT']);O=Path(os.environ['ICI_SUBMISSION_ROOT']);pt=72/25.4
pdfmetrics.registerFont(TTFont('Arial',os.environ.get('ICI_ARIAL_REGULAR','C:/Windows/Fonts/arial.ttf')))
pdfmetrics.registerFont(TTFont('ArialBold',os.environ.get('ICI_ARIAL_BOLD','C:/Windows/Fonts/arialbd.ttf')))
word=PdfReader(P/'qa/word/Manuscript.pdf');writer=PdfWriter();start=len(word.pages)-6
for p in word.pages[:start]:writer.add_page(p)
records=[]
for i,p in enumerate(word.pages[start:],1):
 text=p.extract_text();assert f'Figure {i}' in text and len(p.images)==1,text
 assert re.sub(r'[\s\d]+','',text)=='Figure',text
 buf=io.BytesIO();c=canvas.Canvas(buf,pagesize=(210*pt,297*pt))
 c.setFont('ArialBold',11);c.drawString(20*pt,279*pt,f'Figure {i}')
 c.setFont('Arial',9);c.drawCentredString(105*pt,12*pt,str(start+i));c.showPage();c.save();buf.seek(0)
 page=writer.add_page(PdfReader(buf).pages[0]);fig=PdfReader(O/'Figures'/f'Figure_{i}.pdf').pages[0]
 assert abs(float(fig.mediabox.width)/pt-170)<.05
 height=float(fig.mediabox.height)
 page.merge_transformed_page(fig,Transformation().translate(20*pt,267*pt-height))
 records.append(dict(page=start+i,figure=i,width_mm=170,height_mm=height/pt,method='Native vector PDF merged onto a clean A4 page'))
writer.add_metadata({'/Title':'Mapping recurrent epithelial programmes of checkpoint inhibitor colitis to human organoid ligand responses','/Subject':'Manuscript with six vector figures'})
writer.write(O/'Manuscript_with_Figures.pdf')
(P/'qa/manuscript_vector_assembly.json').write_text(json.dumps({'pages':len(writer.pages),'word_text_pages':start,'figures':records},indent=2))
print('Manuscript assembled:',len(writer.pages),'pages;',start,'unchanged Word text pages and six native vector figure pages.')
