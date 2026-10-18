from pathlib import Path
import subprocess,json
from concurrent.futures import ThreadPoolExecutor
from PIL import Image,ImageOps,ImageDraw
from pypdf import PdfReader
ROOT=Path(__file__).resolve().parents[1];POP=Path(r'C:\Users\frede\.cache\codex-runtimes\codex-primary-runtime\dependencies\native\poppler\Library\bin\pdftoppm.exe')
out=ROOT/'review/rendered';out.mkdir(exist_ok=True)
def render(p):
 d=out/p.stem;d.mkdir(exist_ok=True);cmd=[str(POP),'-r','150','-png',str(p),str(d/'page')];r=subprocess.run(cmd,capture_output=True,text=True);assert r.returncode==0,(p,r.stderr)
 images=sorted(d.glob('page-*.png'));reader=PdfReader(p);assert len(images)==len(reader.pages)
 for start in range(0,len(images),4):
  sheet=Image.new('RGB',(1800,1800),'#dddddd');draw=ImageDraw.Draw(sheet)
  for j,f in enumerate(images[start:start+4]):
   im=Image.open(f).convert('RGB');im.thumbnail((880,840));x=(j%2)*900+(900-im.width)//2;y=(j//2)*900+40;sheet.paste(im,(x,y));draw.text(((j%2)*900+20,(j//2)*900+10),p.stem+' / '+f.stem,fill='black')
  sheet.save(out/(p.stem+f'_contact_{start//4+1}.png'))
 return dict(file=str(p),pages=len(images),page_sizes_pt=[list(map(float,page.mediabox[2:])) for page in reader.pages],render_command=cmd,render_exit_code=r.returncode,rendered_pages=[str(x) for x in images])
with ThreadPoolExecutor(max_workers=3) as ex:records=list(ex.map(render,sorted((ROOT/'figures').glob('*.pdf'))))
(ROOT/'review/pdf_render_record.json').write_text(json.dumps(records,indent=2),encoding='utf8');print([(Path(r['file']).name,r['pages']) for r in records])
