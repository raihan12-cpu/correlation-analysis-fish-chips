"""Known-answer and invalid-input tests; no private workbooks required."""
import unittest
import numpy as np
import openpyxl
from correlation_core import bh, ranks, correlations, keyed, worksheet_records
from clustered_heatmap import average_linkage

class CorrelationTests(unittest.TestCase):
    def test_average_tied_ranks(self):
        np.testing.assert_allclose(ranks([9,1,1,5]),[4,1.5,1.5,3])

    def test_bh_known_answer_and_ties(self):
        np.testing.assert_allclose(bh([.01,.04,.03,.002]),[.02,.04,.04,.008])
        np.testing.assert_allclose(bh([.01,.01,.9]),[.015,.015,.9])

    def test_exact_zero_and_perfect(self):
        matrix=np.array([[-1,1,1],[0,-2,0],[1,1,-1]])
        _,rows=correlations(matrix,'Pearson',['a','b','c'])
        self.assertAlmostEqual(rows[0]['r'],0)
        self.assertAlmostEqual(rows[0]['p'],1)
        self.assertAlmostEqual(rows[1]['r'],-1)
        self.assertAlmostEqual(rows[1]['p'],0)
        self.assertEqual(rows[1]['n'],3)

    def test_known_nontrivial_probability(self):
        # r=0.5 with n=4: the Pearson null distribution is uniform on [-1,1].
        x=np.array([-1,-1,1,1]);z=np.array([-1,1,-1,1])
        _,rows=correlations(np.column_stack([x,.5*x+np.sqrt(.75)*z]),'Pearson',['x','y'])
        self.assertAlmostEqual(rows[0]['r'],.5)
        self.assertAlmostEqual(rows[0]['p'],.5)

    def test_spearman_monotonic(self):
        _,rows=correlations([[1,1],[2,4],[3,9],[4,16]],'Spearman',['x','y'])
        self.assertAlmostEqual(rows[0]['r'],1)

    def test_rejects_invalid_matrix(self):
        for a in [[[1,2],[1,3],[1,4]],[[1,2],[2,np.nan],[3,4]],[[1,2],[2,3]]]:
            with self.assertRaises(ValueError):correlations(a,'Pearson',['x','y'])

    def test_duplicate_and_named_headers(self):
        with self.assertRaises(ValueError):keyed([{'id':'A'},{'id':'A'}],['id'])
        w=openpyxl.Workbook();s=w.active;s.append(['value','id']);s.append([3,'A'])
        self.assertEqual(worksheet_records(s,['id','value']),[{'id':'A','value':3}])
        s.append([4,None])
        with self.assertRaises(ValueError):worksheet_records(s,['id','value'])

    def test_clustering_known_answer(self):
        links,order=average_linkage(np.array([[0,.2,.8],[.2,0,.6],[.8,.6,0]]))
        np.testing.assert_allclose([r[2] for r in links],[.2,.7])
        self.assertEqual(sorted(order),[0,1,2])

if __name__=='__main__':unittest.main()
