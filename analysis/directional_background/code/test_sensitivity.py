from common import *
from scipy.stats import rankdata
import unittest
class Fixtures(unittest.TestCase):
 def test_strict_same_donor(self):
  self.assertEqual(pattern([-1,-1,-1],[1,1,1]),'BOTH')
  self.assertEqual(pattern([-1,0,-1],[1,1,1]),'DEPLETED_ONLY')
  self.assertEqual(pattern([-1,-1,-1],[1,np.nan,1]),'NE')
  self.assertEqual(pattern([-1,1,-1],[1,1,-1]),'OTHER_OBSERVED')
 def test_ratio_factors(self):
  c=np.array([10,20,40,80])[:,None]*np.array([1,2,4,1,2,4])[None,:]
  sf=size_factors(c);np.testing.assert_allclose(sf,[.5,1,2,.5,1,2])
  np.testing.assert_allclose((c/sf).std(axis=1),0,atol=1e-12)
 def test_rank_ties(self):np.testing.assert_allclose(rankdata([0,0,2,4],method='average')/4,[.375,.375,.75,1])
 def test_matching_exclusion_and_repeatability(self):
  x=np.r_[np.zeros(90),np.arange(10)];d=(x>0).astype(float);ix=np.array([95,96,97,98,99])
  a,s,m,b=bg_library(x,d,ix,'fixture',20);a2,_,_,_=bg_library(x,d,ix,'fixture',20)
  np.testing.assert_array_equal(a,a2);self.assertFalse(np.isin(a,ix).any());self.assertTrue(all(len(set(v))==5 for v in a))
 def test_less_decline_not_recovery(self):
  v=np.array([-.01,-.01,-.004]);lib=np.array([[0],[1]])
  r=bg_stats(v,lib,-.004);self.assertGreater(r['matched_background_relative_delta'],0);self.assertLess(-.004,0)
 def test_fixed_composition_weights(self):
  t=np.array([2.,8.]);c=np.array([2.,8.]);w=np.array([.8,.2]);tw=np.array([.2,.8])
  self.assertNotEqual(float(t@tw-c@w),0);self.assertEqual(float(t@w-c@w),0)
if __name__=='__main__':unittest.main()
