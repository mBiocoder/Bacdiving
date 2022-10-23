import unittest
import pandas as pd

try:
    from  BacDiving.treeplots_maker import df2newick
except ImportError:
    raise ImportWarning('Bacdiving not installed.')

class TestUtil(unittest.TestCase):
    def test_df2newick(self):
        data = {
            "x": ['A', 'A', 'B', 'B', 'B'],
            "y": ['Ab', 'Ac', 'Ba', 'Ba', 'Bd'],
            "z": ['Abb', 'Acc', 'Bad', 'Bae', 'Bdd']
        }
        df = pd.DataFrame(data)

        self.assertEqual(df2newick(df, "z"), "(Abb,Acc,Bad,Bae,Bdd);")


if __name__ == '__main__':
    unittest.main()
