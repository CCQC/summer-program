import numpy as np
import sys
sys.path.insert(0,'../../extra-files')
import masses

class Molecule(object):
    """
    Python class used to store the geometry, masses, and nuclear charges of a molecule
    """

    def __init__(self, geom):
        """
        Creates a molecule with a geometry and corresponding units
        """

        lines = geom.strip().split("\n")

        units = lines[1]
        
        natom = int(lines[0])

        atoms = []
        geom = []

        for line in lines[2:]:
            atom, x, y, z = line.split()
            atoms.append(atom)
            geom.append([float(x), float(y), float(z)])

        self.units = units
        self.natom = natom
        self.atoms = atoms
        self.geom = np.array(geom)

        mass = []
        charges =[]

        for atom in atoms:
            mass.append(float(masses.get_mass(atom)))
            charges.append(int(masses.get_charge(atom))) #might want to include symbols

        self.mass = mass
        self.charges = charges

    def to_bohr(self):
        """
        Converts geometry to Bohr
        """

        if self.units == "Angstrom":
            self.units = "Bohr"
            self.geom *= 1.889725989
    
    def to_angstrom(self):
        """
        Converts geometry to Angstroms
        """

        if self.units == "Bohr":
            self.units = "Angstrom"
            self.geom /= 1.889725989

    def xyz_string(self):
        """
        Returns a string representing the molecule in xyz format
        """
        formatting = "{:3s} {: >15.10f} {: >15.10f} {: >15.10f} | {: >15.10f}\t{:d}\n"
        out = "\nNumber of atoms: {:d}\n\nUnits: {:s}\n\nAtom\tPostition x\tPosition y\tPosition z  | \t\tMass\tCharge\n".format(self.natom, self.units)
        for atom, (x, y, z), mass, charge in zip(self.atoms, self.geom, self.mass, self.charges):
            out += formatting.format(atom, x, y, z, mass, charge)
        
        return out

if __name__ == "__main__":
    geom = open("../../extra-files/molecule.xyz").read()
    mol = Molecule(geom)
    print(mol.xyz_string())
