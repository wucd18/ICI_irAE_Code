from common import *
import subprocess
exe=r'D:\Software\Bioinformatics\R\R-4.4.1\bin\x64\Rscript.exe'
env=os.environ.copy();env.update(SystemRoot=r'C:\Windows',windir=r'C:\Windows',LC_ALL='C',LANG='C',TEMP=str(ROOT/'temp'),TMP=str(ROOT/'temp'),OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
tests={'base':'cat("BODY_COMPLETE\\n")','fgsea_load':'suppressPackageStartupMessages(library(fgsea));cat("BODY_COMPLETE\\n")','fgsea_unload':'suppressPackageStartupMessages(library(fgsea));unloadNamespace("fgsea");cat("BODY_COMPLETE\\n")','fgsea_tinythread':'Sys.setenv(RCPP_PARALLEL_BACKEND="tinythread",RCPP_PARALLEL_NUM_THREADS="1");suppressPackageStartupMessages(library(fgsea));cat("BODY_COMPLETE\\n")'}
rows=[]
for name,expr in tests.items():
    cmd=[exe,'--vanilla','-e',expr];t=now();p=subprocess.run(cmd,cwd=ROOT,env=env,capture_output=True,text=True,errors='replace');row=dict(test=name,command=cmd,started=t,finished=now(),exit_code=p.returncode,stdout=p.stdout,stderr=p.stderr);rows.append(row);print(name,p.returncode,flush=True)
jw(ROOT/'review/runtime_exit_diagnostics.json',rows)
