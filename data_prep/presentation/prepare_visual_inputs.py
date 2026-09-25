"""Extracted code-only subset of recorded prepare_visual.py: verified input copies
and F3C member-file alias. No selection-rule rewrite or publication-tree copying.
Original-member identity and extraction diff are in docs/script_index.csv / local review.
"""
from __future__ import annotations
import argparse, hashlib, json
from pathlib import Path

def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def copy_plan(source, destination, manifest):
    source=source.resolve();destination=destination.resolve()
    code=Path(__file__).resolve().parents[2]
    if destination==code or code in destination.parents or destination in code.parents or destination in (Path.home().resolve(),Path(destination.anchor)):
        raise ValueError('Cannot write to code tree, its ancestors, home or drive root')
    if source==destination or source in destination.parents or destination in source.parents:
        raise ValueError('Input/output roots must be separate')
    if destination.exists():raise ValueError('Output must be a new directory')
    rows=json.loads(manifest.read_text(encoding='utf-8-sig'))
    if not isinstance(rows,list) or not rows:raise ValueError('Nonempty selected-input manifest required')
    seen=set();selected=[]
    for row in rows:
        raw=row['path'];rel=Path(raw)
        if rel.is_absolute() or '..' in rel.parts or ':' in raw or '\\' in raw or raw in seen:
            raise ValueError('Unsafe/duplicate snapshot member')
        seen.add(raw);p=(source/rel).resolve()
        if source not in p.parents or not p.is_file():raise ValueError('Missing/escaping snapshot member')
        if p.stat().st_size!=row['bytes'] or digest(p)!=row['sha256']:raise ValueError('Snapshot hash mismatch')
        if raw.startswith('source_data/') and rel.suffix=='.csv':
            selected.append((p,destination/rel,row['sha256']))
    if not selected:raise ValueError('No source CSVs selected')
    return selected

def prepare(source,destination,manifest,write=False):
    rows=copy_plan(source,destination,manifest)
    # Explicit alias from the original preparer; no gene selection/recalculation.
    member=next((r for r in rows if r[1].relative_to(destination).as_posix()=='source_data/panel_sources/F3C_members.csv'),None)
    alias=destination/'source_data/panel_sources/F3C.csv'
    current=next((r for r in rows if r[1]==alias),None)
    if member and current and member[2]!=current[2]:raise ValueError('Conflicting F3C alias')
    if member and not current:rows.append((member[0],alias,member[2]))
    if write:
        destination.mkdir(parents=True,exist_ok=False)
        for src,dst,sha in rows:
            dst.parent.mkdir(parents=True,exist_ok=True)
            with dst.open('xb') as f:f.write(src.read_bytes())
            if digest(dst)!=sha:raise ValueError('Output copy checksum mismatch')
    return [{'source':str(a),'destination':str(b),'sha256':c} for a,b,c in rows]

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for n in ('source','destination','manifest'):p.add_argument('--'+n,type=Path,required=True)
    mode=p.add_mutually_exclusive_group();mode.add_argument('--plan',action='store_true');mode.add_argument('--write',action='store_true')
    a=p.parse_args()
    try:print(json.dumps(prepare(a.source,a.destination,a.manifest,a.write),indent=2));return 0
    except (ValueError,OSError,KeyError) as e:p.exit(2,str(e)+'\n')
if __name__=='__main__':raise SystemExit(main())
