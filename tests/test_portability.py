"""Synthetic engineering fixtures only. Never import scientific scripts."""
import ast,csv,importlib.util,json,os,subprocess,sys,tempfile,unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'workflow'))
from runtime import load_config,sha256,verify_manifest
from stage_workspace import plan,write_plan
spec=importlib.util.spec_from_file_location('visual_input_copy',ROOT/'data_prep/presentation/prepare_visual_inputs.py')
visual=importlib.util.module_from_spec(spec);spec.loader.exec_module(visual)

class EngineeringTests(unittest.TestCase):
 def setUp(self):
  self.tmp=tempfile.TemporaryDirectory(prefix='ici_engineering_');self.base=Path(self.tmp.name)
  self.raw=self.base/'raw';self.raw.mkdir();self.archive=self.base/'archive';self.archive.mkdir()
  self.work=self.base/'work';self.cfg=self.base/'paths.json'
  self.values={'CODE_ROOT':str(ROOT),'INPUT_ROOT':str(self.raw),'ARCHIVE_ROOT':str(self.archive),'WORK_ROOT':str(self.work)}
  self.save()
 def tearDown(self):self.tmp.cleanup()
 def save(self):self.cfg.write_text(json.dumps(self.values),encoding='utf-8')
 def manifest(self,text='sample\tvalue\na\t1\nb\t2\n'):
  p=self.raw/'synthetic.tsv';p.write_text(text,encoding='utf-8')
  row={'id':'synthetic_engineering_only','root':'INPUT_ROOT','path':p.name,'size':p.stat().st_size,'sha256':sha256(p),'source':'synthetic engineering fixture','role':'schema test','stage':'test','columns':['sample','value'],'unique_keys':['sample']}
  m=self.base/'manifest.json';m.write_text(json.dumps({'version':'fixture','files':[row]}),encoding='utf-8');return m
 def test_plan_has_no_writes(self):
  before=set(self.base.rglob('*'));result=plan(self.cfg);self.assertTrue(result['files']);self.assertEqual(before,set(self.base.rglob('*')))
 def test_new_run_and_collision(self):
  a=plan(self.cfg);b=plan(self.cfg);self.assertNotEqual(a['run_id'],b['run_id']);write_plan(a)
  self.assertRaises(ValueError,plan,self.cfg,'all',a['run_id'])
  for row in a['files']:self.assertEqual(sha256(row['target']),row['sha256'])
 def test_overlap_roots(self):
  for root in [ROOT,self.raw,self.archive,self.base]:
   self.values['WORK_ROOT']=str(root);self.save();self.assertRaises(ValueError,load_config,self.cfg)
 def test_relative_root_rejected(self):
  self.values['WORK_ROOT']='relative';self.save();self.assertRaises(ValueError,load_config,self.cfg)
 def test_output_must_be_inside_work(self):
  self.values['OUTPUT_ROOT']=str(self.raw/'out');self.save();self.assertRaises(ValueError,load_config,self.cfg)
 def test_missing_root(self):
  self.values['INPUT_ROOT']=str(self.base/'missing');self.save();self.assertRaises(ValueError,load_config,self.cfg)
 def test_run_traversal(self):self.assertRaises(ValueError,plan,self.cfg,'all','../escape')
 def test_manifest_success(self):
  m=self.manifest();_,r=load_config(self.cfg);self.assertEqual(len(verify_manifest(m,r)),1)
 def test_changed_input(self):
  m=self.manifest();(self.raw/'synthetic.tsv').write_text('changed');_,r=load_config(self.cfg);self.assertRaises(ValueError,verify_manifest,m,r)
 def test_missing_input(self):
  m=self.manifest();(self.raw/'synthetic.tsv').rename(self.raw/'other.tsv');_,r=load_config(self.cfg);self.assertRaises(ValueError,verify_manifest,m,r)
 def test_duplicate_metadata(self):
  m=self.manifest('sample\tvalue\na\t1\na\t2\n');_,r=load_config(self.cfg);self.assertRaises(ValueError,verify_manifest,m,r)
 def test_bad_schema(self):
  m=self.manifest('wrong\tvalue\na\t1\n');_,r=load_config(self.cfg);self.assertRaises(ValueError,verify_manifest,m,r)
 def test_empty_manifest(self):
  m=self.base/'empty.json';m.write_text('{"version":"x","files":[]}');_,r=load_config(self.cfg);self.assertRaises(ValueError,verify_manifest,m,r)
 def test_duplicate_manifest(self):
  m=self.manifest();o=json.loads(m.read_text());o['files']*=2;m.write_text(json.dumps(o));_,r=load_config(self.cfg);self.assertRaises(ValueError,verify_manifest,m,r)
 def test_cli_help_plan_no_io(self):
  env={**os.environ,'PYTHONDONTWRITEBYTECODE':'1'};before=set(self.base.rglob('*'))
  for args in [['--help'],['--config',str(self.cfg),'--plan']]:
   c=subprocess.run([sys.executable,'-B',str(ROOT/'workflow/stage_workspace.py'),*args],cwd=self.base,capture_output=True,text=True,env=env)
   self.assertEqual(c.returncode,0,c.stderr)
  self.assertEqual(before,set(self.base.rglob('*')))
 def test_visual_plan_copy_and_alias(self):
  p=self.archive/'source_data/panel_sources/F3C_members.csv';p.parent.mkdir(parents=True);p.write_text('label\nsynthetic\n')
  m=self.base/'selected.json';m.write_text(json.dumps([{'path':p.relative_to(self.archive).as_posix(),'bytes':p.stat().st_size,'sha256':sha256(p)}]))
  dst=self.work/'display';rows=visual.prepare(self.archive,dst,m);self.assertEqual(len(rows),2);self.assertFalse(dst.exists())
  visual.prepare(self.archive,dst,m,True);self.assertEqual(sha256(dst/'source_data/panel_sources/F3C.csv'),sha256(p))
  self.assertRaises(ValueError,visual.prepare,self.archive,dst,m,True)
 def test_scientific_config_tail_unchanged(self):
  # Exact scientific tail is tracked by the local diff review; constants are not re-invented here.
  text=(ROOT/'config/recorded/project_config.R').read_text();self.assertNotIn('F:/WCD',text);self.assertIn('Sys.getenv(key',text)
 def annotation_loader(self):
  # Execute only the inspected pure JSON loader, never the scientific module/imports/main.
  tree=ast.parse((ROOT/'plotting/tables/historical/build_final_synthesis_tables.py').read_text())
  fn=next(x for x in tree.body if isinstance(x,ast.FunctionDef) and x.name=='load_annotations')
  ns={'json':json,'Path':Path};exec(compile(ast.Module(body=[fn],type_ignores=[]),'pure_annotation_loader','exec'),ns)
  return ns['load_annotations']
 def test_external_annotations_preserve_synthetic_records(self):
  obj={'ligands':['synthetic_a'],'gate_rows':[['g','label','synthetic_status','text','source']], 'claims':[['c','text','synthetic_status','source','boundary']], 'candidate_decisions':{'synthetic_a':'synthetic_value','SORL1':'synthetic_value'}}
  p=self.base/'annotations.json';p.write_text(json.dumps(obj));self.assertEqual(self.annotation_loader()(p),obj)
 def test_external_annotations_reject_duplicate_keys(self):
  p=self.base/'annotations.json';p.write_text('{"ligands":[],"ligands":[]}');self.assertRaises(ValueError,self.annotation_loader(),p)
 def test_external_annotations_require_records(self):
  p=self.base/'annotations.json';p.write_text('{}');self.assertRaises(ValueError,self.annotation_loader(),p)
 def test_visual_output_rejects_code_tree(self):
  self.assertRaises(ValueError,visual.copy_plan,self.archive,ROOT/'uncreated_fixture',self.base/'missing.json')

if __name__=='__main__':unittest.main(verbosity=2)
