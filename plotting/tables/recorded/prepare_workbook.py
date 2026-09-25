from pathlib import Path
import json,pandas as pd
O=Path(__file__).resolve().parents[1];I=O/'input_snapshot'
defs=json.loads((I/'supplementary/reader_tables/table_definitions.json').read_text(encoding='utf-8'))
items=[]
for d in defs:
 p=I/Path(d['file'].replace('\\','/'));f=pd.read_csv(p)
 for col in f.select_dtypes('object').columns:
  f[col]=f[col].map(lambda s:s.replace('Platform samples','GSE210037 sample pseudobulk').replace('whole_spot_sample','GSE210037 sample pseudobulk') if isinstance(s,str) else s)
 items.append({'sheet':d['table'],'title':d['title'],'columns':list(f),'rows':f.astype(object).where(pd.notnull(f),None).values.tolist(),'note':d['note'],'source':'reader_tables/'+p.name+'; full sources: '+', '.join(d['full_source_ids'])})
(O/'supplement/workbook_input.json').write_text(json.dumps(items,ensure_ascii=False,allow_nan=False),encoding='utf-8')
