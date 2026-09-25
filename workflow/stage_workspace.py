#!/usr/bin/env python3
"""Plan/stage code only. Never execute science, copy data or download.
--help and --plan are read-only. --write creates a fresh external run.
--verify-inputs checks bytes/schema only; it does not establish biological validity.
"""
from __future__ import annotations
import argparse, json, re, sys, uuid
from datetime import datetime,timezone
from pathlib import Path
sys.dont_write_bytecode = True
from runtime import CODE_ROOT,load_config,read_index,safe_relative,sha256,within,verify_manifest
def plan(config,family='all',run_id=None):
    cfg,roots=load_config(config)
    run_id=run_id or datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')+'_'+uuid.uuid4().hex[:8]
    if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_-]{0,79}',run_id): raise ValueError('Unsafe run ID')
    run=roots['WORK_ROOT']/'runs'/run_id
    if run.exists(): raise ValueError('Refusing to overwrite existing run: '+str(run))
    if not within(run.resolve(),roots['WORK_ROOT']): raise ValueError('Run escapes WORK_ROOT')
    rows=[];seen=set()
    for r in read_index():
        if family!='all' and r['execution_family']!=family: continue
        src=(CODE_ROOT/safe_relative(r['repository_path'])).resolve()
        rel=safe_relative(r['execution_family'])/safe_relative(r['workspace_relative_path'])
        if not within(src,CODE_ROOT) or not src.is_file() or rel.as_posix().casefold() in seen: raise ValueError('Missing/ambiguous source')
        digest=sha256(src)
        if digest!=r.get('adapted_sha256',r['sha256']): raise ValueError('Source checksum changed: '+r['repository_path'])
        seen.add(rel.as_posix().casefold())
        rows.append(dict(source=str(src),target=str(run/rel),sha256=digest,family=r['execution_family']))
    if not rows: raise ValueError('No code selected')
    env={'ICI_INPUT_ROOT':str(roots['INPUT_ROOT']),'ICI_ARCHIVE_ROOT':str(roots['ARCHIVE_ROOT']),
         'ICI_WORK_ROOT':str(run),'WCD_PROJECT_ROOT':str(run/'original_v1'),
         'ICI_ORIGINAL_CODE_ROOT':str(run/'original_v1'),'PYTHONDONTWRITEBYTECODE':'1'}
    for key,value in cfg.get('executables',{}).items():
        if value:
            if not Path(value).is_absolute() or not Path(value).is_file(): raise ValueError('Invalid executable: '+key)
            env['ICI_'+key.upper()]=value
    return dict(run_id=run_id,run_root=str(run),files=rows,environment=env,scope='CODE_STAGING_ONLY',
                analysis_executed=False,data_copied=False,scientific_pipeline_status='PORTABILITY_PENDING',
                warning='Read docs/PORTABILITY.md. Legacy data read/write contracts remain incomplete.')
def write_plan(result):
    run=Path(result['run_root']);run.mkdir(parents=True,exist_ok=False)
    for r in result['files']:
        p=Path(r['target']);p.parent.mkdir(parents=True,exist_ok=True)
        with p.open('xb') as f: f.write(Path(r['source']).read_bytes())
        if sha256(p)!=r['sha256']: raise ValueError('Staged checksum mismatch')
    with (run/'staging_manifest.json').open('x',encoding='utf-8') as f: json.dump(result,f,indent=2)
def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--config',type=Path,required=True);p.add_argument('--family',default='all')
    p.add_argument('--run-id');mode=p.add_mutually_exclusive_group()
    mode.add_argument('--plan',action='store_true');mode.add_argument('--write',action='store_true')
    p.add_argument('--verify-inputs',type=Path);a=p.parse_args()
    try:
        result=plan(a.config,a.family,a.run_id)
        if a.verify_inputs:
            _,roots=load_config(a.config);result['verified_inputs']=verify_manifest(a.verify_inputs,roots)
        if a.write:write_plan(result)
        print(json.dumps(result,indent=2));return 0
    except (ValueError,OSError,KeyError) as e:
        print(str(e),file=sys.stderr);return 2
if __name__=='__main__': raise SystemExit(main())
