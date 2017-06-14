import unittest
from molecule import Molecule

class Test(unittest.TestCase):
    """ Tests for Molecule class in 'molecule.py' """
    
    def test_self_variables(self):
        """ Are the initial self variables (units, natom, masses, charges, geom)
            initialized correctly? """
        
        with open('input1.xyz', 'r') as file1:
            mol1 = Molecule(file1.read())
        self.assertEqual(mol1.units, 'Angstrom', 'checking units')
        self.assertEqual(mol1.natom, 3, 'checking natom')
        self.assertEqual(mol1.labels, ['O', 'H', 'H'], 'checking labels')
        self.assertEqual(mol1.masses, [15.99491461956, 1.00782503207, 1.00782503207], 'checking masses')
        self.assertEqual(mol1.charges, [8, 1, 1], 'checking charges')
        self.assertEqual(mol1.geom,
                         [[0.0, 0.0, -0.0711762954],
                          [0.0, -0.8916195680, 0.5648097613],
                          [0.0, 0.8916195680, 0.5648097613]],
                         'checking geometry')

    def test_copy_function(self):
        """ Does the copy function initialize all variables of the new molecule?
            
            Also, is the new molecule truly a different object? (ie: changing properties of
            the original does not affect the copy)"""
                
        with open('input1.xyz', 'r') as file1:
            mol1 = Molecule(file1.read())
        mol2 = mol1.copy()
        self.assertEqual(mol1.units, mol2.units, 'checking units')
        self.assertEqual(mol1.natom, mol2.natom, 'checking natom')
        self.assertEqual(mol1.labels, mol2.labels, 'checking labels')
        self.assertEqual(mol1.masses, mol2.masses, 'checking masses')
        self.assertEqual(mol1.charges, mol2.charges, 'checking charges')
        self.assertEqual(mol1.geom, mol2.geom, 'checking geometry')

        mol2.to_bohr()
        self.assertNotEqual(mol1.units, mol2.units)
        self.assertNotEqual(mol1.geom, mol2.geom)




if __name__ == '__main__':
    unittest.main()
