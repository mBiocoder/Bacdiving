import unittest

try:
    from  BacDiving.visualizations_maker import symbol_remover
except ImportError:
    raise ImportWarning('Bacdiving not installed.')

class TestUtil(unittest.TestCase):
    def test_symbol_remover(self):
        values = ['30%', '100 %', '0.5%']
        self.assertEqual(symbol_remover(values, "%"), ['30', '100 ', '0.5'])

        values = ['30__', '100 __', '0.5__']
        self.assertEqual(symbol_remover(values, "__"), ['30', '100 ', '0.5'])


if __name__ == '__main__':
    unittest.main()
