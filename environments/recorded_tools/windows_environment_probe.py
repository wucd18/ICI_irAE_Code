from pathlib import Path
import os,subprocess,json,ctypes
from datetime import datetime,timezone
O=Path(__file__).resolve().parents[1]
buf=ctypes.create_unicode_buffer(32768)
assert ctypes.windll.kernel32.GetWindowsDirectoryW(buf,len(buf))>0
win=Path(buf.value);assert (win/'System32/cmd.exe').exists()
def folder(csidl):
 b=ctypes.create_unicode_buffer(32768)
 hr=ctypes.windll.shell32.SHGetFolderPathW(None,csidl,None,0,b)
 if hr!=0:
  fallback={26:Path(os.environ['USERPROFILE'])/'AppData/Roaming',28:Path(os.environ['USERPROFILE'])/'AppData/Local',35:Path(win.anchor)/'ProgramData'}[csidl]
  assert fallback.is_dir(),(hr,fallback)
  return str(fallback)
 return b.value
standard={'SystemRoot':str(win),'windir':str(win),'ComSpec':str(win/'System32/cmd.exe'),'APPDATA':folder(26),'LOCALAPPDATA':folder(28),'PROGRAMDATA':folder(35),'TEMP':str(O/'temp'),'TMP':str(O/'temp'),'TMPDIR':str(O/'temp'),'R_USER':str(O),'LANG':'C','LC_ALL':'C'}
write=lambda name,obj:(O/'review'/name).write_text(json.dumps(obj,ensure_ascii=False,indent=2),encoding='utf-8')
write('child_environment_restore.json',{'utc':datetime.now(timezone.utc).isoformat(),'scope':'process only; no machine/user environment mutation','original_present':{k:bool(os.environ.get(k)) for k in standard},'restored':standard})
R=os.environ['ICI_RSCRIPT_TARGETED']
records=[]
for label,extra in [('baseline',{}),('systemroot_only',{'SystemRoot':str(win),'windir':str(win)}),('standard',standard)]:
 env={**os.environ,**extra};env.update(TEMP=str(O/'temp'),TMP=str(O/'temp'))
 cmd=[R,'--vanilla','-e','library(rlang);cat("LOADED_OK\\n")']
 p=subprocess.run(cmd,env=env,cwd=O,capture_output=True,timeout=30)
 (O/'logs'/f'{label}_R.stdout.txt').write_bytes(p.stdout);(O/'logs'/f'{label}_R.stderr.txt').write_bytes(p.stderr)
 row={'probe':label,'command':cmd,'exit_code':p.returncode};records.append(row);print(row,flush=True)
cmd=[str(win/'System32/curl.exe'),'--fail','--location','--connect-timeout','10','--max-time','30','--output',str(O/'new_data/GSE189184_metadata.soft'),'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE189184&targ=self&view=full&form=text']
p=subprocess.run(cmd,env={**os.environ,**standard},cwd=O,capture_output=True,timeout=35)
(O/'logs/curl_restored.stdout.txt').write_bytes(p.stdout);(O/'logs/curl_restored.stderr.txt').write_bytes(p.stderr)
records.append({'probe':'curl_restored','command':cmd,'exit_code':p.returncode});print('curl',p.returncode,flush=True)
write('environment_AB_test.json',records)
