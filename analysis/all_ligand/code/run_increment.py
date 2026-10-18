from common import *
import subprocess,time,concurrent.futures
exe=r'D:\Software\Bioinformatics\R\R-4.4.1\bin\x64\Rscript.exe'
env=os.environ.copy();env.update(SystemRoot=r'C:\Windows',windir=r'C:\Windows',LC_ALL='C',LANG='C',TEMP=str(ROOT/'temp'),TMP=str(ROOT/'temp'),OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
def run(shard):
 command=[exe,'--vanilla',str(ROOT/'code/increment_models.R'),str(ROOT),str(shard)];start=now()
 with (ROOT/'logs'/f'R_{shard}_stdout.log').open('w',encoding='utf8') as out,(ROOT/'logs'/f'R_{shard}_stderr.log').open('w',encoding='utf8') as err:
  p=subprocess.run(command,cwd=ROOT,env=env,stdout=out,stderr=err)
 completion=ROOT/('review/R_fixture.json' if shard=='fixture' else f'logs/shard_{shard}_complete.json')
 record=dict(command=command,started=start,finished=now(),exit_code=p.returncode,completed_records_exist=completion.is_file(),status='PASS' if p.returncode==0 and completion.is_file() else 'ARTIFACTS_COMPLETE_PROCESS_EXIT_FAILED' if completion.is_file() else 'FAILED_INCOMPLETE')
 jw(ROOT/'logs'/f'execution_{shard}.json',record)
 if not completion.is_file():raise RuntimeError(record)
 return record
if __name__=='__main__':
 fixture=run('fixture');print('FIXTURE',fixture,flush=True)
 assert jr(ROOT/'review/R_fixture.json')['status']=='PASS'
 jobs=pd.read_csv(ROOT/'inputs/new_model_jobs.tsv',sep='\t');workers=4
 for i in range(workers):save(jobs.iloc[i::workers],f'inputs/jobs_{i}.tsv')
 futures=[]
 with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
  futures=[pool.submit(run,i) for i in range(workers)]
  for f in concurrent.futures.as_completed(futures):print('WORKER_COMPLETE',f.result(),flush=True)
 print('ALL_NEW_MODELS_COMPLETE',len(jobs),flush=True)
