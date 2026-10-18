"""Presentation-wrapper safety checks only; no scientific module is imported."""
from pathlib import Path
import unittest, tempfile, subprocess, sys, json, hashlib
ROOT=Path(__file__).resolve().parents[1]
WRAPPER=ROOT/'production/final_20261018/rebuild_publication.py'
class ProductionBindings(unittest.TestCase):
 def test_existing_output_is_rejected_before_inputs(self):
  with tempfile.TemporaryDirectory() as s:
   p=Path(s);cfg=p/'bindings.json';cfg.write_text('{}')
   q=subprocess.run([sys.executable,str(WRAPPER),'--bindings',str(cfg),'--out',str(p)],capture_output=True)
   self.assertNotEqual(q.returncode,0);self.assertIn(b'Choose a new output directory',q.stderr)
 def test_missing_source_table_fails_without_creating_output(self):
  with tempfile.TemporaryDirectory() as s:
   p=Path(s);out=p/'new';cfg=p/'bindings.json';cfg.write_text(json.dumps({'ICI_DISPLAY_INPUT_ROOT':str(p),'ICI_SUBMISSION_BASELINE':str(p)}))
   q=subprocess.run([sys.executable,str(WRAPPER),'--bindings',str(cfg),'--out',str(out)],capture_output=True)
   self.assertNotEqual(q.returncode,0);self.assertIn(b'Missing approved display source tables',q.stderr);self.assertFalse(out.exists())
 def test_recovered_code_matches_packaged_manifest(self):
  for r in json.loads((ROOT/'docs/source_code_manifest.json').read_text()):
   p=ROOT/r['file']
   self.assertTrue(p.is_file(),r['file'])
   self.assertEqual(hashlib.sha256(p.read_bytes()).hexdigest(),r['packaged_sha256'],r['file'])
if __name__=='__main__':unittest.main()

