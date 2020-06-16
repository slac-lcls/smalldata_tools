import unittest
import smalldata_tools.utilities_FitCenter as fc
import numpy as np

# To run, run following command from top level directory
# $ python -m unittest discover test

# Constants
IMAGE = np.load('test/test_data/AvImg_xppo5616_Run716.npy')
MASK = np.load('test/test_data/AvImg_xppo5616_Run716_mask.npy')
SHAPE = (1692, 1691)
MEAN = 0.0022
RADII = [326, 336, 171, 161]
XCEN = 846.8528
YCEN = 846.8122
RADII_RANSAC = [
    327.7875, 
    336.5052, 
    173.3988, 
    164.8353
]
XCEN_RANSAC = 847.0735
YCEN_RANSAC = 846.9079

class FitCenterTest(unittest.TestCase):

    def test_find_edges(self):
        """Verify canny returns expected results with default params"""
        edges, sparse_edges = fc.find_edges(IMAGE, MASK)
        self.assertEqual(edges.shape, SHAPE)
        self.assertEqual(sparse_edges.shape, SHAPE)
        self.assertEqual(round(edges.mean(), 4), MEAN)
        self.assertEqual(round(sparse_edges.mean(), 4), MEAN)

    def test_iterate_center(self):
        """Verify we get expected results from iterate center"""
        _, sparse_edges = fc.find_edges(IMAGE, MASK)
        r_vals, x, y = fc.iterate_center(sparse_edges)
        self.assertEqual(r_vals, RADII)
        self.assertEqual(round(x, 4), XCEN)
        self.assertEqual(round(y, 4), YCEN)

    def test_ransac_result(self):
        _, sparse_edges = fc.find_edges(IMAGE, MASK)
        res, _, _ = fc.ransac_result(sparse_edges, [XCEN, YCEN], RADII)
        self.assertEqual([round(r, 4) for r in res['R']], RADII_RANSAC)
        self.assertEqual(round(res['x_cen'], 4), XCEN_RANSAC)
        self.assertEqual(round(res['y_cen'], 4), YCEN_RANSAC)

if __name__ == '__main__':
    unittest.main()