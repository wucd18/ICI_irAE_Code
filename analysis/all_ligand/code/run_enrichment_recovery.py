from common import *
import subprocess,concurrent.futures
exe=r'D:\Software\Bioinformatics\R\R-4.4.1\bin\x64\Rscript.exe'
df=pd.read_csv(ROOT/'inputs/frozen_module_definitions.csv');jw(ROOT/'inputs/expected_module_member_counts.json',df.groupby('module').gene.nunique().to_dict())
env=os.environ.copy();env.update(SystemRoot=r'C:\Windows',windir=r'C:\Windows',LC_ALL='C',LANG='C',TEMP=str(ROOT/'temp'),TMP=str(ROOT/'temp'),OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
def run(shard):
    command=[exe,'--vanilla',str(ROOT/'code/recover_enrichment.R'),str(ROOT),str(shard)];start=now()
    with (ROOT/'logs'/f'enrichment_{shard}_stdout.log').open('w',encoding='utf8') as out,(ROOT/'logs'/f'enrichment_{shard}_stderr.log').open('w',encoding='utf8') as err:
        p=subprocess.run(command,cwd=ROOT,env=env,stdout=out,stderr=err)
    complete=(ROOT/'logs'/f'enrichment_{shard}_complete.json').exists()
    record=dict(command=command,started=start,finished=now(),exit_code=p.returncode,completed_records_exist=complete,status='PASS' if p.returncode==0 and complete else 'ARTIFACTS_COMPLETE_PROCESS_EXIT_FAILED' if complete else 'FAILED_INCOMPLETE')
    jw(ROOT/'logs'/f'execution_enrichment_{shard}.json',record)
    if not complete:raise RuntimeError(record)
    return record
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
    for record in pool.map(run,range(4)):print(record,flush=True)
print('ENRICHMENT_RECOVERY_COMPLETE',flush=True)
