import masses
import numpy as np

class Molecule(object):
    
    def __init__(self, geom_str):
        lines=geom_str.splitlines()
        self.natom = int(lines[0])
        self.units = lines[1]
        self.labels = []
        self.mass = []
        #self.charge = []
        self.geom = []
        for line in lines[2:]:
            atom,x,y,z = line.split()
            self.labels.append(str(atom))
            self.mass.append(float(masses.get_mass(atom)))
           # self.charge.append(int(masses.get_charge(atom)))
            self.geom.append([float(x),float(y),float(z)])
        self.geom = np.array(self.geom)
    
    def to_bohr(self):
        if self.units == "Bohr":
            print("Already in Bohr")
        else:
            self.geom *= 1.8897 
    
    def to_angstrom(self):
        if self.units == "Angstrom":
            print("Already in Angstrom")
        else:
            self.geom *= 0.5292
    
    def __str__(self): 
        out = "{}\n\n".format(self.natom)
        for atom,line in zip(self.labels, self.geom):
            x,y,z=line
            out += "{:2s}{:>10.5f}{:>10.5f}{:>10.5f}\n".format(atom,x,y,z)
        return out
    
    def copy(self):
        return Molecule(self)
 
    @staticmethod 
    def from_file(infile = "geom.xyz"):
        with open(infile) as f:
            return Molecule(f.read())

       
