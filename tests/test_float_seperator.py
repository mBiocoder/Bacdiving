import unittest

try:
    from  BacDiving.visualizations_maker import float_seperator
except ImportError:
    raise ImportWarning('Bacdiving not installed.')

class TestUtil(unittest.TestCase):
    def test_float_seperator(self):
        values = ['30-35', '29', '100-101', '10.5 - 12.5']
        self.assertEqual(float_seperator(values), [32.5, 29.0, 100.5, 11.5])


if __name__ == '__main__':
    unittest.main()
