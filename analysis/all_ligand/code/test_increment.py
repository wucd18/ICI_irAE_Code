from common import *
from analysis_functions import *
import unittest
class IncrementTests(unittest.TestCase):
    def test_bh_planned_ne_denominator(self):
        np.testing.assert_allclose(bh_planned([.01,.04,np.nan,.03],5),[.05,1/15,np.nan,1/15],equal_nan=True)
        np.testing.assert_allclose(bh_planned([0,.5,.5,1],4),[0,2/3,2/3,1])
    def test_direction_same_donors(self):
        self.assertEqual(direction_pattern([-1,-2,-3],[1,2,3]),'BOTH')
        self.assertEqual(direction_pattern([-1,-2,-3],[1,0,3]),'UP_ONLY')
        self.assertEqual(direction_pattern([-1,0,-3],[1,2,3]),'DEPLETED_ONLY')
        self.assertEqual(direction_pattern([-1,1,-3],[1,-1,3]),'OTHER_OBSERVED')
        self.assertEqual(direction_pattern([-1,-2],[1,2]),'NE')
        self.assertEqual(direction_pattern([-1,-2,-3],[1,2,np.nan]),'NE')
    def test_tiny_values_not_rounded(self):
        self.assertEqual(direction_pattern([-1e-12]*3,[1e-12]*3),'BOTH')
    def test_depleted_not_inverted(self):
        a=pd.DataFrame([dict(ligand='L',celltype='C',donor='d',module=DEP,treatment=x,raw_mean_rank=y,pools='p',n_cells=20,value_origin='fixture') for x,y in [('Treatment',.7),('Control',.6)]])
        self.assertAlmostEqual(paired_deltas(a).delta.iloc[0],.1)
    def test_invalid_duplicate(self):
        a=pd.DataFrame([dict(ligand='L',celltype='C',donor='d',module=DEP,treatment='Treatment',raw_mean_rank=.7,pools='p',n_cells=20,value_origin='fixture')]*2)
        with self.assertRaises(AssertionError):paired_deltas(a)
if __name__=='__main__':
    result=unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(IncrementTests))
    jw(ROOT/'review/python_fixture.json',dict(status='PASS' if result.wasSuccessful() else 'FAIL',tests=result.testsRun,failures=len(result.failures),errors=len(result.errors)))
    sys.exit(0 if result.wasSuccessful() else 1)
