import masses as m
import copy

class Molecule :

    # Default Methods

    def __init__(self, data) :

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
        if self.units == 'Angstrom' :
            self.units = 'Bohr'
            self.geom =[[coord * 1.889725989 for coord in atom] for atom in self.geom]

    def to_angstrom(self) :
        if self.units == 'Bohr' :
            self.units = 'Angstrom'
            self.geom =[[coord / 1.889725989 for coord in atom] for atom in self.geom]

    def xyz_string(self) :
        ans = "Molecule: " + " ".join(self.labels)
        return ans
    
    def copy(self) :
        return copy.deepcopy(self)

    # Advanced Methods

    def __len__(self) :
        return self.natom

    def __str__(self) :
        ans = "Molecule: " + " ".join(self.labels)
        return ans
    
    def __repr__(self) :
        return "\n".join(self.labels)

    def __iter__(self) :
        pass

    def __add__(self, mol) :
        if self.units == 'Angstroms' :
            mol.to_bohr()
        else :
            mol.to_angstrom()
        newmol = self.copy()
        newmol.natom += mol.natom
        newmol.geom += mol.geom
        newmol.labels += mol.labels
        newmol.masses += mol.masses
        return newmol

if __name__ == '__main__':
    print("This is a supporting class and not meant to be run on its own")
