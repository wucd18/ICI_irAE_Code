from pathlib import Path
import subprocess,json
from PIL import Image,ImageOps,ImageDraw
import pypdf
ROOT=Path(__file__).resolve().parents[1];out=ROOT/'review/pdf_rendered';out.mkdir(exist_ok=True)
exe=Path(r'C:\Users\frede\.cache\codex-runtimes\codex-primary-runtime\dependencies\native\poppler\Library\bin\pdftoppm.exe')
records=[]
for pdf in sorted((ROOT/'figures').glob('*.pdf')):
    r=subprocess.run([str(exe),'-r','150','-png',str(pdf),str(out/pdf.stem)],capture_output=True,text=True)
    if r.returncode:raise RuntimeError(r.stderr)
    reader=pypdf.PdfReader(pdf);pngs=sorted(out.glob(pdf.stem+'-*.png'));assert len(pngs)==len(reader.pages)
    for i,p in enumerate(reader.pages):
        text=p.extract_text();assert len(text)>80
        records.append(dict(file=pdf.name,page=i+1,width_pt=float(p.mediabox.width),height_pt=float(p.mediabox.height),text_chars=len(text),rendered_png=str(pngs[i]),render_exit_code=r.returncode))
        (out/f'{pdf.stem}-{i+1}.txt').write_text(text,encoding='utf8')
    thumbs=[]
    for p in pngs:
        im=Image.open(p).convert('RGB');im.thumbnail((690,950));thumb=Image.new('RGB',(720,1000),'#e8e8e8');thumb.paste(im,((720-im.width)//2,30));ImageDraw.Draw(thumb).text((10,5),p.name,fill='black');thumbs.append(thumb)
    for start in range(0,len(thumbs),4):
        batch=thumbs[start:start+4];cols=min(2,len(batch));rows=(len(batch)+1)//2;sheet=Image.new('RGB',(720*cols,1000*rows),'white')
        for j,im in enumerate(batch):sheet.paste(im,((j%2)*720,(j//2)*1000))
        sheet.save(out/f'{pdf.stem}_contact_{start//4+1}.png')
(ROOT/'review/pdf_render_checks.json').write_text(json.dumps(records,indent=2),encoding='utf8')
print('PDF_RENDERED',len(records),'pages',flush=True)
