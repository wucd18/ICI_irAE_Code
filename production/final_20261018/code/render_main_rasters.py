import os
from pathlib import Path
import subprocess,os,shutil,json
from PIL import Image
R=Path(os.environ['ICI_PRODUCTION_ROOT'])
exe=os.environ.get('PDFTOPPM') or shutil.which('pdftoppm')
if not exe:raise RuntimeError('Set PDFTOPPM to an installed executable.')
records=[]
for i in range(1,7):
    p=R/'main_figures'/f'Figure_{i}.pdf';out=p.with_name(f'Figure_{i}_600dpi')
    cp=subprocess.run([exe,'-r','600','-singlefile','-png',str(p),str(out)],capture_output=True,text=True)
    if cp.returncode:raise RuntimeError(cp.stderr)
    with Image.open(out.with_suffix('.png')) as im:records.append(dict(figure=i,width_px=im.width,height_px=im.height,requested_dpi=600,origin='Rendered directly from source-redrawn vector PDF'))
(R/'qa/raster_600dpi_record.json').write_text(json.dumps(records,indent=2))
print('Six 600 dpi PNGs rendered from vector figures.')
