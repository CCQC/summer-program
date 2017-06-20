import sys
sys.path.append('../../extra-files')
from numpy import array

import masses as m

class Molecule :
    """
    Represenation of a molecule
    """
    # Default Methods

    def __init__(self, data) :
        """
        data -- contents of an .xyz file as a str
        """
        data = data.split('\n')
        self.natom = int(data[0])
        self.units = data[1]
        
        self.labels = []
        self.masses = []
        self.charges = []
        self.geom = []
        
        for atom in data[2:] :
            
            atom = atom.split()
            
            if len(atom) != 4 :
                continue # Throw Error?
            self.labels.append(atom[0])
            self.masses.append(m.get_mass(atom[0]))
            self.charges.append(m.get_charge(atom[0]))
            self.geom.append(atom[1:])
        
        self.geom = [[float(coord) for coord in atom] for atom in self.geom]
        
    def to_bohr(self) :
        """
        Converts geometry to units of Bohrs
        """
        if self.units == 'Angstrom' :
            self.units = 'Bohr'
            self.geom =[[coord * 1.889725989 for coord in atom] for atom in self.geom]

    def to_angstrom(self) :
        """
        Converts geometry to units of Angstroms
        """
        if self.units == 'Bohr' :
            self.units = 'Angstrom'
            self.geom =[[coord / 1.889725989 for coord in atom] for atom in self.geom]

    def xyz_string(self) :
        """
        Returns (str) a representation of the molecular geometry in .xyz format
        """
        coordinates = ["\t".join(str(coord) for coord in atom) for atom in self.geom]
        atoms = [zipped[0] + "\t" + zipped[1]
                 for zipped in zip(self.labels, coordinates)]
    
        return "\n".join([str(self.natom), self.units, "\n".join(atoms)])


    def copy(self) :
        """
        Returns (Molecule) a new molecule identical to this one
        """
        return Molecule(self.xyz_string())

    # Advanced Methods

    def __len__(self) :
        """
        Returns (int) the number of atoms in the molecule
        """
        return self.natom

    def __str__(self) :
        """
        Returns (str) a nice represenation of the molecule
        """
        coordinates = ["(" +
                       ", ".join(str(coord) for coord in atom) +
                       ")"
                       for atom in self.geom]
        atoms = [zipped[0] + " : " + zipped[1]
                for zipped in zip(self.labels, coordinates)]
        summary = "\nMolecule with %s atoms given in units of %s:\n" % (self.natom, self.units)
        return summary + "\n" + "\n".join(atoms)
    
    def __repr__(self) :
        """
        Returns (str) a representation of the molecular geometry in .xyz format
        """
        return self.xyz_string()

    def __iter__(self) :
        """
        Returns (iter) an iterator over the atoms of molecule
        Each atom is formatted as a tuple: (label, array of coordinates)
        """
        self.iteration_index = 0
        return self

    def __next__(self) :
        """
        Returns (tuple) an atom and increments the iteration over atoms of the molecule.
        Each atom is formatted as a tuple: (label, array of coordinates)
        """
        if(self.iteration_index < self.natom) :
            ans = (self.labels[self.iteration_index], array(self.geom[self.iteration_index]))
            self.iteration_index += 1
            return ans
        else :
            raise StopIteration

    def __add__(self, mol) :
        """
        mol -- another object of the class molecule
        
        Adds mol to this object, increasing the number of atoms
        """
        
        if self.units == 'Angstroms' :
            mol.to_bohr()
        else :
            mol.to_angstrom()
        newmol = self.copy()
        newmol.natom += mol.natom
        newmol.labels += mol.labels
        newmol.masses += mol.masses
        newmol.charges += mol.charges
        newmol.geom += mol.geom

        return newmol

if __name__ == '__main__':
    print("This is a supporting class and not meant to be run on its own")
