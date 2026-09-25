"""Prepare display records only. No inferential tests, enrichment, or model fitting."""
from common import *
from collections import Counter
from decimal import Decimal
from statistics import median

MODULES=[('CURATED_IFN_VISIBILITY','IFN reference'),('ICI_COLITIS_EPITHELIAL_UP_STRICT','Disease-up'),('ICI_COLITIS_EPITHELIAL_UP_NON_IFN','IFN-gene-excluded residual'),('ICI_COLITIS_EPITHELIAL_LOSS_STRICT','Disease-depleted')]
LIGANDS=['TNFSF12','BMP2','BMP4','INHBA','TGFB1','TGFB2','TGFB3','No_Cytomix']
FAMILIES=[('Stem','Stem'),('Transit_Amplifying','TA'),('Colonocyte','Colonocyte'),('Goblet','Goblet')]
CONDITIONS=['Healthy','ICI colitis','UC']
CONDITION_NAMES={'HC':'Healthy','ICI':'ICI colitis','UC':'UC','Healthy':'Healthy','ICI colitis':'ICI colitis'}
COMPARTMENTS=[('IEC','Epithelium'),('Immune','Immune'),('Stromal','Stromal'),('Endothelial','Endothelium'),('Enteric_glia','Enteric glia')]
GENES=['SORL1','TNFSF12','TNFRSF12A']
tables={n:cr(src(n)) for n in range(64,79)}
coords=cr(SOURCE/'current_panel_data/display/S5B_all_donor_plot_coordinates.csv')
def records(n):
 return [dict(r,source_file=src(n).relative_to(SOURCE).as_posix(),source_row=i+2) for i,r in enumerate(tables[n])]
def unique(rows,keys):
 d={tuple(r[k] for k in keys):r for r in rows}
 assert len(d)==len(rows),keys
 return d
G=unique(records(76),['ligand','celltype','module'])
S=unique(records(77),['ligand','celltype','module','donor','treatment'])
P=unique(records(78),['ligand','celltype','module'])
models=unique(records(75),['ligand','celltype'])
assert {r['ligand'] for r in G.values()}==set(LIGANDS)
assert {r['module'] for r in G.values()}=={m for m,t in MODULES}
assert {r['celltype'] for r in G.values()}=={m for m,t in FAMILIES}
unique(coords,['ligand','celltype','module','donor'])
spec=dict(modules=MODULES,ligands=LIGANDS,families=FAMILIES,conditions=CONDITIONS,condition_names=CONDITION_NAMES,compartments=COMPARTMENTS,genes=GENES,S5_heatmaps=[],S5_donors=[],S5_reference=[],S7_Oxford=[],S7_genes=[],S7_epithelium=[],S7_medians=[],S7_annotations=[])
for i,(mod,label) in enumerate(MODULES):
 for ligand in LIGANDS:
  for cell,display in FAMILIES:
   r=G[ligand,cell,mod].copy();r.update(panel=chr(ord('a')+i),display_row=LIGANDS.index(ligand),display_col=[f[0] for f in FAMILIES].index(cell))
   r['display_text']='NE' if r['status']!='ESTIMATED' else f"{float(r['NES']):.2f}"+('*' if Decimal(r['FDR_predeclared_family'])<Decimal('0.05') else '')
   spec['S5_heatmaps'].append(r)
   if ligand=='No_Cytomix':spec['S5_reference'].append(dict(r,panel='i',display_row=i))
for i,r0 in enumerate(coords):
 r=r0.copy();r.update(source_file='current_panel_data/display/S5B_all_donor_plot_coordinates.csv',source_row=i+2,panel=chr(ord('e')+[m for m,t in MODULES].index(r['module'])))
 for arm in ['Control','Treatment']:
  raw=S[r['ligand'],r['celltype'],r['module'],r['donor'],arm]
  assert r[arm]==raw['raw_mean_rank']
 assert abs(Decimal(r['change'])-(Decimal(r['Treatment'])-Decimal(r['Control'])))<Decimal('2e-15')
 assert (r['ligand'],r['celltype'],r['module']) in P
 assert (r['ligand'],r['celltype']) in models
 spec['S5_donors'].append(r)
assert len(spec['S5_donors'])==len(S)//2
ox=records(65);xf=records(69);xg=records(70);oxstats=unique(records(66),['module']);xcov=unique(records(71),['module'])
unique(ox,['patient','module']);unique(xf,['patient','compartment','module']);unique(xg,['patient','compartment','gene'])
for i,(mod,label) in enumerate(MODULES):
 for r in ox:
  if r['module']==mod:spec['S7_Oxford'].append(dict(r,panel=chr(ord('a')+i),display_condition=CONDITION_NAMES[r['condition']],value=r['oriented_score']))
 for r in xf:
  if r['module']==mod and r['compartment']=='IEC':spec['S7_epithelium'].append(dict(r,panel=chr(ord('h')+i),display_condition=CONDITION_NAMES[r['condition']],value=r['oriented_rank']))
 st=oxstats[mod,];cv=xcov[mod,]
 spec['S7_annotations'].append(dict(st,panel=chr(ord('a')+i),annotation=f"ICI colitis vs Healthy: P={float(st['P']):.4f}; FDR={float(st['FDR_new_four_module_family']):.4f}"))
 spec['S7_annotations'].append(dict(cv,panel=chr(ord('h')+i),annotation=f"Targeted coverage: {cv['available']}/{cv['requested']} fixed genes"))
for i,gene in enumerate(GENES):
 for r in xg:
  if r['gene']==gene:spec['S7_genes'].append(dict(r,panel=chr(ord('e')+i),display_condition=CONDITION_NAMES[r['condition']],value=r['mean_molecules_per_cell']))
for key in ['S7_Oxford','S7_epithelium']:
 for panel in sorted({r['panel'] for r in spec[key]}):
  for cond in CONDITIONS:
   group=sorted([r for r in spec[key] if r['panel']==panel and r['display_condition']==cond],key=lambda r:r['patient'])
   mid=median(Decimal(r['value']) for r in group)
   # Medians are display summaries of existing donor measurements, not new tests.
   spec['S7_medians'].append(dict(panel=panel,condition=cond,n=len(group),value=str(mid),source_file=group[0]['source_file'],source_rows=';'.join(str(r['source_row']) for r in group),operation='median of unchanged donor scores',patients=';'.join(r['patient'] for r in group)))
   if key=='S7_Oxford' and cond!='UC':
    old=Decimal(oxstats[group[0]['module'],]['median_control' if cond=='Healthy' else 'median_case'])
    assert abs(mid-old)<Decimal('2e-15')
spec['donors']=sorted({r['donor'] for r in coords})
spec['NES_limit']=str(__import__('math').ceil(max(abs(float(r['NES'])) for r in G.values() if r['status']=='ESTIMATED')))
spec['source_counts']=dict(S5_GSEA=len(G),S5_estimable=sum(r['status']=='ESTIMATED' for r in G.values()),S5_NE=sum(r['status']!='ESTIMATED' for r in G.values()),S5_paired_donor_points=len(coords),S7_Oxford_points=len(ox),S7_gene_points=len(xg),S7_epithelial_points=len(spec['S7_epithelium']),Xenium_donors=len({r['patient'] for r in xf}))
jw(I/'plot_data/plot_spec.json',spec)
captions=(I/'inputs_readonly/Supplementary_figure_legends.md').read_text(encoding='utf8')
for figure in ['S5','S7']:
 text=re.search(r'## Supplementary Figure '+figure+r'\s+([\s\S]*?)(?=\n## |\Z)',captions)[1].strip()
 if figure=='S7':
  new=text.replace('3 healthy, 7 ICI-colitis','3 Healthy, 7 ICI colitis').replace('with healthy controls','with Healthy controls').replace('colours distinguish healthy controls','colours distinguish Healthy controls')
 else:new=text
 write(I/'plot_data'/f'{figure}_caption.txt',new)
 captions=captions.replace(text,new)
write(I/'Supplementary_figure_legends.md',captions)
print(json.dumps(spec['source_counts'],indent=2))
