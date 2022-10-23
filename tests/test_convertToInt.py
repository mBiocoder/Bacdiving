import unittest

try:
    from  BacDiving.visualizations_maker import convertToInt
except ImportError:
    raise ImportWarning('Bacdiving not installed.')

class TestUtil(unittest.TestCase):
    def test_convertToInt(self):
        values = ['30-35', '29', '100-101', '10.5 - 12.5']
        self.assertEqual(convertToInt(values), [30, 31, 32, 33, 34, 35, 29, 100, 101, 10, 11, 12])


if __name__ == '__main__':
    unittest.main()
