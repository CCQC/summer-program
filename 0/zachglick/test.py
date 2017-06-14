import unittest
import numpy
from molecule import Molecule

class Test(unittest.TestCase):
    """ 
    Tests for Molecule class in 'molecule.py'
    """
    
    def test_self_variables(self):
        """ 
        Are the self variables (units, natom, masses, charges, geom) from the xyz file initialized correctly?
        """
        
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
        """ 
        Does the copy function correctly initialize all variables of a new molecule? Also, is the new molecule
        truly a different object? (ie: changing properties of the original does not affect the copy and vice versa)?
        """
                
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

    def test_unit_conversion_accuracy(self):
        """ 
        1.0 Angstrom is approximately 1.889725989 Bohr. Is this conversion (and its reverse) carried out correctly?
        """

        with open('input1.xyz', 'r') as file1:
            mol1 = Molecule(file1.read())
        mol2 = mol1.copy()
        mol2.to_bohr()
        for i in range(mol1.natom) :
            self.assertAlmostEqual(mol1.geom[i][0] * 1.889725989, mol2.geom[i][0])
            self.assertAlmostEqual(mol1.geom[i][1] * 1.889725989, mol2.geom[i][1])
            self.assertAlmostEqual(mol1.geom[i][2] * 1.889725989, mol2.geom[i][2])

    def test_unit_conversion_symmetry(self):
        """ 
        Does converting back and forth between bohrs and angstroms introduce and compound rounding errors? Each 
        operation should exactly reverse its counterpart.
        """
                
        with open('input1.xyz', 'r') as file1:
            mol1 = Molecule(file1.read())
        mol2 = mol1.copy()
        for count in range(10000) :
            mol2.to_bohr()
            mol2.to_angstrom()
        self.assertEqual(mol1.geom, mol2.geom)

    def test_iteration(self):
        """ 
        Can the atoms of the molecule be iterated over? Each atom should be a tuple containing the atom's label 
        and a numpy array of the 3 cartesian coordinates.
        """

        with open('input1.xyz', 'r') as file1:
            mol1 = Molecule(file1.read())
        
        newlabels = []
        for atom in mol1:
            newlabels += atom[0]

        self.assertEqual(mol1.labels, newlabels, 'checking labels')
        #TODO: assert that the coordinates are also correct

if __name__ == '__main__':
    unittest.main()
