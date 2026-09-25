#!/usr/bin/env python3
"""Safe static repository checks; no imports/execution of scientific scripts.

This is a structure/source-integrity check, NOT a model reproducibility test.
"""
from __future__ import annotations
import argparse
import ast
import csv
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BAN_EXT = {'.zip','.gz','.tar','.7z','.rar','.pdf','.png','.jpg','.jpeg','.tif','.tiff','.eps','.svg',
           '.xlsx','.xls','.docx','.pptx','.rds','.rdata','.robj','.h5ad','.h5','.hdf5','.mtx','.loom','.parquet','.feather','.npy','.npz'}
BAN_ROOT = {'results','result','outputs','output','data','input','inputs','source_data','figures','tables','manuscript','logs',
            'visual_input_snapshot','input_snapshot','organoid_donor_extension','spatial_extension','xenium_extension','01_data','02_results'}
SKIP_DIR = {'.git','__pycache__','.venv','venv','node_modules'}
TOKEN = re.compile(r'ghp_[A-Za-z0-9]{30,}|github_pat_[A-Za-z0-9_]{35,}|-----BEGIN (?:RSA |EC )?PRIVATE KEY-----|AKIA[A-Z0-9]{16}')

def main() -> int:
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--r-parse', action='store_true', help='Also run R parse() only if Rscript is available')
    parser.add_argument('--rscript', help='Explicit Rscript executable for parse-only checks')
    args=parser.parse_args()
    failures=[]; filelist=[]
    for p in sorted(ROOT.rglob('*')):
        rel=p.relative_to(ROOT)
        if any(part in SKIP_DIR for part in rel.parts):
            continue
        if p.is_symlink():
            failures.append({'file':rel.as_posix(),'problem':'Symlink not allowed in public source collection'})
            continue
        if p.is_file():
            if p.name=='paths.local.json' or p.name.startswith('.env'):
                # Local ignored config may exist, but must never be in a tracked set.
                if (ROOT/'.git').exists():
                    c=subprocess.run(['git','ls-files','--error-unmatch',rel.as_posix()],cwd=ROOT,capture_output=True)
                    if c.returncode==0:failures.append({'file':rel.as_posix(),'problem':'Local config tracked by Git'})
                continue
            filelist.append(p)
            if p.suffix.lower() in BAN_EXT or rel.parts[0] in BAN_ROOT:
                failures.append({'file':rel.as_posix(),'problem':'Data/result/document/archive artifact in code-only tree'})
            if p.stat().st_size > 5 * 1024 * 1024:
                failures.append({'file':rel.as_posix(),'problem':'File exceeds repository policy size (5 MiB)'})
    with (ROOT/'docs/script_index.csv').open(encoding='utf-8',newline='') as f:
        rows=list(csv.DictReader(f))
    for row in rows:
        p=ROOT/row['repository_path']
        if not p.is_file() or hashlib.sha256(p.read_bytes()).hexdigest()!=row.get('adapted_sha256',row['sha256']):
            failures.append({'file':row['repository_path'],'problem':'Missing or changed recorded source'})
    parsed=0
    for p in filelist:
        if p.suffix=='.py':
            try:ast.parse(p.read_text(encoding='utf-8-sig'),filename=str(p));parsed+=1
            except (SyntaxError,UnicodeError) as e:failures.append({'file':str(p.relative_to(ROOT)),'problem':'Python syntax: '+str(e)})
        try:text=p.read_text(encoding='utf-8-sig')
        except (UnicodeError,OSError):continue
        if TOKEN.search(text):
            failures.append({'file':str(p.relative_to(ROOT)),'problem':'Possible credential pattern; value not logged'})
    # Staged local imports must resolve within their original script folder.
    recorded_modules={Path(r['workspace_relative_path']).stem for r in rows if r['repository_path'].endswith('.py')}
    staged={(r['execution_family'],r['workspace_relative_path']) for r in rows}
    for r in rows:
        if not r['repository_path'].endswith('.py'):continue
        tree=ast.parse((ROOT/r['repository_path']).read_text(encoding='utf-8-sig'))
        for n in ast.walk(tree):
            if isinstance(n,ast.ImportFrom) and n.module and n.module.split('.')[0] in recorded_modules:
                module=n.module.split('.')[0]
                target=(Path(r['workspace_relative_path']).parent/(module+'.py')).as_posix()
                if (r['execution_family'],target) not in staged:
                    failures.append({'file':r['repository_path'],'problem':'Staged local import unresolved: '+target})
    r_status='NOT_RUN';r_count=0
    if args.r_parse:
        exe=args.rscript or os.environ.get('ICI_RSCRIPT')
        if exe is None:r_status='NOT_AVAILABLE'
        else:
            r_status='PASS'
            for p in filelist:
                if p.suffix=='.R':
                    expr='parse(file=commandArgs(trailingOnly=TRUE)[1]);cat("PARSE_ONLY_OK\\n")'
                    c=subprocess.run([exe,'--vanilla','-e',expr,str(p)],capture_output=True,text=True,timeout=30)
                    r_count+=1
                    if c.returncode:
                        failures.append({'file':str(p.relative_to(ROOT)),'problem':'R parse failed','message':c.stderr[:500]});r_status='FAIL'
    report={'status':'STRUCTURE_SYNTAX_PASS' if not failures else 'FAIL','scope':'Static code-only/source-integrity validation; no biological analyses',
            'retained_original_scripts':len(rows),'python_files_parsed':parsed,'R_parse_status':r_status,'R_files_parsed':r_count,
            'files_checked':len(filelist),'failures':failures,
            'not_verified':['data availability and file schemas','absolute path portability','all dynamic code dependencies',
                            'historical model versions','final CII figure-version identity','end-to-end numerical reproduction'],
            'note':'Credential-pattern scan is limited and is not a complete secret-security certification.'}
    print(json.dumps(report,indent=2))
    return 0 if not failures else 1

if __name__=='__main__':
    raise SystemExit(main())
