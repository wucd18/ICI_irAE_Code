"""Engineering contracts only: no scientific imports or downloads."""
from __future__ import annotations
import csv, hashlib, json, re
from pathlib import Path
CODE_ROOT = Path(__file__).resolve().parents[1]
def within(path, root):
    return path == root or root in path.parents
def sha256(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(4*1024*1024),b''): h.update(b)
    return h.hexdigest()
def safe_relative(raw):
    p=Path(raw)
    if not raw or p.is_absolute() or '..' in p.parts or ':' in raw or '\\' in raw:
        raise ValueError('Unsafe relative path: '+raw)
    return p
def load_config(path):
    cfg=json.loads(Path(path).read_text(encoding='utf-8-sig')); roots={}
    for key in ('CODE_ROOT','INPUT_ROOT','ARCHIVE_ROOT','WORK_ROOT'):
        value=cfg.get(key)
        if not isinstance(value,str) or not value.strip() or not Path(value).is_absolute():
            raise ValueError(key+' must be an explicit absolute path')
        p=Path(value).resolve()
        if p==Path(p.anchor) or p==Path.home().resolve(): raise ValueError('Root/home not allowed: '+key)
        roots[key]=p
    if roots['CODE_ROOT']!=CODE_ROOT.resolve(): raise ValueError('Incorrect CODE_ROOT')
    for a,b in [('WORK_ROOT','CODE_ROOT'),('WORK_ROOT','INPUT_ROOT'),('WORK_ROOT','ARCHIVE_ROOT'),('CODE_ROOT','INPUT_ROOT'),('CODE_ROOT','ARCHIVE_ROOT')]:
        if within(roots[a],roots[b]) or within(roots[b],roots[a]):
            raise ValueError('Overlapping roots: '+a+'/'+b)
    for k in ('INPUT_ROOT','ARCHIVE_ROOT'):
        if not roots[k].is_dir(): raise ValueError('BLOCKED_INPUT: missing '+k)
    if cfg.get('OUTPUT_ROOT'):
        p=Path(cfg['OUTPUT_ROOT'])
        if not p.is_absolute() or not within(p.resolve(),roots['WORK_ROOT']) or p.resolve()==roots['WORK_ROOT']:
            raise ValueError('OUTPUT_ROOT must be strictly below WORK_ROOT')
    return cfg,roots
def read_index(root=CODE_ROOT):
    with (root/'docs/script_index.csv').open(encoding='utf-8-sig',newline='') as f: rows=list(csv.DictReader(f))
    if len({r['repository_path'] for r in rows})!=len(rows): raise ValueError('Duplicate source path')
    return rows
def verify_manifest(path,roots):
    obj=json.loads(Path(path).read_text(encoding='utf-8-sig')); rows=obj.get('files')
    if not obj.get('version') or not isinstance(rows,list) or not rows:
        raise ValueError('BLOCKED_INPUT: nonempty versioned input manifest required')
    ids=set();paths=set();report=[]
    for row in rows:
        required=('id','root','path','size','sha256','source','role','stage')
        if not all(k in row for k in required): raise ValueError('Incomplete input specification')
        if row['id'] in ids: raise ValueError('Duplicate input id')
        ids.add(row['id'])
        if row['root'] not in ('INPUT_ROOT','ARCHIVE_ROOT'): raise ValueError('Input requires read-only root')
        base=roots[row['root']];p=(base/safe_relative(row['path'])).resolve()
        if not within(p,base) or p in paths: raise ValueError('Escaping/duplicate input path')
        paths.add(p)
        if not row['source'] or not row['role'] or not isinstance(row['size'],int) or row['size']<1 or not re.fullmatch('[0-9a-f]{64}',row['sha256']):
            raise ValueError('BLOCKED_INPUT: unresolved provenance/hash/size: '+row['id'])
        if not p.is_file(): raise ValueError('BLOCKED_INPUT: missing '+row['path'])
        if p.stat().st_size!=row['size'] or sha256(p)!=row['sha256']:
            raise ValueError('BLOCKED_INPUT: changed '+row['path'])
        if row.get('columns') or row.get('unique_keys'):
            with p.open(encoding='utf-8-sig',newline='') as f:
                reader=csv.DictReader(f,delimiter=row.get('delimiter','\t'))
                if not reader.fieldnames or len(reader.fieldnames)!=len(set(reader.fieldnames)):
                    raise ValueError('Missing/duplicate columns: '+row['id'])
                if not (set(row.get('columns',[]))|set(row.get('unique_keys',[]))).issubset(reader.fieldnames):
                    raise ValueError('Missing required columns: '+row['id'])
                keys=set()
                for item in reader:
                    if None in item or any(v is None for v in item.values()): raise ValueError('Malformed row')
                    if row.get('unique_keys'):
                        key=tuple(item[k] for k in row['unique_keys'])
                        if not all(key) or key in keys: raise ValueError('Empty/duplicate critical key')
                        keys.add(key)
        report.append({k:row[k] for k in required})
    return report
