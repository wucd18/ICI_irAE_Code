"""New path-binding wrapper around recovered, actual publication generators.

This wrapper executes presentation stages only; no scientific analysis entrypoint
is imported. An existing output directory is rejected.
"""
from pathlib import Path
import argparse,os,json,shutil,subprocess,datetime,hashlib
def main():
 ap=argparse.ArgumentParser()
 ap.add_argument('--bindings',type=Path,required=True)
 ap.add_argument('--out',type=Path,required=True)
 ap.add_argument('--stage',choices=['all','figures','tables','documents'],default='all')
 a=ap.parse_args();cfg=json.loads(a.bindings.read_text(encoding='utf-8-sig'))
 out=a.out.resolve()
 if out.exists():raise FileExistsError('Choose a new output directory: '+str(out))
 inp=Path(cfg['ICI_DISPLAY_INPUT_ROOT']).resolve();base=Path(cfg['ICI_SUBMISSION_BASELINE']).resolve()
 if not (inp/'source_data/new/three_method_donor_effects.tsv').is_file():raise FileNotFoundError('Missing approved display source tables')
 if not (base/'Manuscript.docx').is_file():raise FileNotFoundError('Missing the exact reviewed manuscript baseline')
 for k in ['plot_python','artifact_python','node','pdftoppm','powershell']:
  if not Path(cfg[k]).is_file():raise FileNotFoundError(k+': '+cfg[k])
 out.mkdir(parents=True)
 for sub in ['inputs','panels','panel_data','main_figures','qa','config']:(out/sub).mkdir()
 shutil.copytree(base,out/'submission_final')
 here=Path(__file__).resolve().parent;code=here/'code'
 for p in (here/'config').glob('*.json'):shutil.copy2(p,out/'config'/p.name)
 env=os.environ.copy();env.update({k:str(v) for k,v in cfg.items() if k.startswith('ICI_')})
 env.update(ICI_PRODUCTION_ROOT=str(out),ICI_SUBMISSION_ROOT=str(out/'submission_final'),PDFTOPPM=cfg['pdftoppm'],PYTHONUTF8='1')
 records=[]
 def run(exe,name,args=()):
  cmd=[cfg[exe],str(code/name),*map(str,args)]
  if name.endswith('.ps1'):cmd=[cfg[exe],'-NoProfile','-File',str(code/name),'-RunRoot',str(out)]
  begin=datetime.datetime.now(datetime.timezone.utc).isoformat()
  q=subprocess.run(cmd,env=env,capture_output=True)
  (out/'qa'/(name+'.stdout.log')).write_bytes(q.stdout);(out/'qa'/(name+'.stderr.log')).write_bytes(q.stderr)
  records.append(dict(command=cmd,started=begin,exit_code=q.returncode,script_sha256=hashlib.sha256((code/name).read_bytes()).hexdigest()))
  (out/'qa/executed_commands.json').write_text(json.dumps(records,indent=2))
  if q.returncode:raise RuntimeError(q.stderr.decode('utf8','replace'))
  print(name,'complete',flush=True)
 if a.stage in ['all','figures']:
  for name in ['plot_main.py','redraw_figure3.py']:run('plot_python',name)
  run('plot_python','redraw_main_4_to_6.py',['--input',inp,'--output',out])
  run('plot_python','plot_supplement_final.py')
  run('artifact_python','assemble_main.py')
  run('artifact_python','assemble_supplement_final.py')
  run('artifact_python','render_main_rasters.py')
 if a.stage in ['all','tables']:
  run('artifact_python','prepare_tables.py')
  # Resolve the preinstalled ESM module explicitly; no installation or upgrade.
  module=Path(cfg['artifact_tool_module']).resolve()
  if not module.is_file():raise FileNotFoundError(module)
  script=(code/'build_workbook.mjs').read_text(encoding='utf8')
  script=script.replace("from '@oai/artifact-tool'","from "+json.dumps(module.as_uri()))
  temp=out/'qa/build_workbook_bound.mjs';temp.write_text(script,encoding='utf8')
  q=subprocess.run([cfg['node'],str(temp),str(out)],env=env,capture_output=True)
  (out/'qa/workbook.stdout.log').write_bytes(q.stdout);(out/'qa/workbook.stderr.log').write_bytes(q.stderr)
  if q.returncode:raise RuntimeError(q.stderr.decode('utf8','replace'))
  for p in (out/'submission_final').glob('*.inspect.ndjson'):shutil.move(p,out/'qa'/p.name)
 if a.stage in ['all','documents']:
  if a.stage=='documents':
   for i in range(1,7):shutil.copy2(base/'Figures'/f'Figure_{i}.pdf',out/'main_figures'/f'Figure_{i}.pdf')
   run('artifact_python','render_main_rasters.py')
  run('artifact_python','revise_manuscript.py')
  run('powershell','export_documents.ps1')
  run('artifact_python','assemble_manuscript_pdf.py')
 print('Presentation stages completed. Visually inspect the new outputs; no analysis model was run.')
if __name__=='__main__':main()
