from common import *
import subprocess,time,platform
step=sys.argv[1];cmd=[sys.executable,'-X','utf8',str(ROOT/'code'/step)]+sys.argv[2:]
started=now();t=time.time();stem=Path(step).stem
if (ROOT/'logs'/f'{stem}.execution.json').exists():stem+='_'+datetime.datetime.now(datetime.timezone.utc).strftime('%H%M%S')
with (ROOT/'logs'/f'{stem}.stdout.log').open('w',encoding='utf8') as out,(ROOT/'logs'/f'{stem}.stderr.log').open('w',encoding='utf8') as err:
 p=subprocess.run(cmd,stdout=out,stderr=err,cwd=ROOT)
record=dict(command=cmd,started=started,finished=now(),seconds=time.time()-t,exit_code=p.returncode,python=platform.python_version(),numpy=np.__version__,pandas=pd.__version__,script_sha256=sha(ROOT/'code'/step))
jw(ROOT/'logs'/f'{stem}.execution.json',record);print(json.dumps(record),flush=True)
print((ROOT/'logs'/f'{stem}.stdout.log').read_text(encoding='utf8')[-3000:]);print((ROOT/'logs'/f'{stem}.stderr.log').read_text(encoding='utf8')[-3000:]);sys.exit(p.returncode)
