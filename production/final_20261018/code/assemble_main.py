import os
from pathlib import Path
import json,xml.etree.ElementTree as ET
from pypdf import PdfReader,PdfWriter,Transformation
R=Path(os.environ['ICI_PRODUCTION_ROOT'])
layout=json.loads((R/'config/layout.json').read_text());pt=72/25.4
for num,spec in layout.items():
 w,h=spec['size'];writer=PdfWriter();page=writer.add_blank_page(w*pt,h*pt)
 root=ET.Element('svg',{'xmlns':'http://www.w3.org/2000/svg','width':f'{w}mm','height':f'{h}mm','viewBox':f'0 0 {w*pt} {h*pt}'})
 for pid,x,y,pw,ph in spec['panels']:
  src=PdfReader(R/'panels'/f'{pid}.pdf').pages[0]
  assert abs(float(src.mediabox.width)-pw*pt)<.1
  page.merge_transformed_page(src,Transformation().translate(x*pt,(h-y-ph)*pt))
  # Nested vector SVG retains native text; prefix IDs to avoid clipping/marker collisions.
  s=(R/'panels'/f'{pid}.svg').read_text(encoding='utf-8');import re
  ids=re.findall(r' id="([^"]+)"',s)
  for iid in ids:s=s.replace('id="'+iid+'"','id="'+pid+'_'+iid+'"').replace('#'+iid,'#'+pid+'_'+iid)
  n=ET.fromstring(s);n.set('x',str(x*pt));n.set('y',str(y*pt));n.set('width',str(pw*pt));n.set('height',str(ph*pt));root.append(n)
 writer.write(R/'main_figures'/f'Figure_{num}.pdf')
 ET.ElementTree(root).write(R/'main_figures'/f'Figure_{num}.svg',encoding='utf-8',xml_declaration=True)
print('Six PDF/SVG figures assembled without rasterisation.')
