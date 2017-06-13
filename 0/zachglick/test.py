import unittest
from molecule import Molecule as mol
print("hey")

class Test(unittest.TestCase):
    """Tests for `primes.py`."""
    
    def test_something_important(self):
        """Is five successfully determined to be prime?"""
        mol1 = mol.Molecule("todo: things")

if __name__ == '__main__':
    unittest.main()
