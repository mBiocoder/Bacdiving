import unittest
import numpy as np

try:
    from  BacDiving.treeplots_maker import newick_to_linkage
except ImportError:
    raise ImportWarning('Bacdiving not installed.')

class TestUtil(unittest.TestCase):
    def test_newick_to_linkage(self):
        newick = '((D,F)E,(B,H)B);' #newick format 8 should be used
        linkage_matrix, labels = newick_to_linkage(newick=newick)

        resulting_matrix = np.array([[0., 3., 2.82842712, 2.], [1., 2., 2.82842712, 2.], [4., 5., 6.32455532, 4.]])

        # For testing purposes decimals of only up to 4 should match
        linkage_matrix = np.around(linkage_matrix, 4)
        resulting_matrix = np.around(resulting_matrix, 4)

        self.assertEqual(linkage_matrix.tolist(), resulting_matrix.tolist())
        self.assertEqual(labels, ['B', 'D', 'F', 'H'])


if __name__ == '__main__':
    unittest.main()
