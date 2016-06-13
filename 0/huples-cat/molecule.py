import numpy as np
import sys
sys.path.insert(0, '../../extra_files')
import masses as atomlib                   ##contains get_mass('x') and get_charge('x') methods

Class Molecule(object):
    def __init__(self):
        f = open('molecule.xyz', r)        ##file object for molecule.xyz is f
        for line in f:                     ##iterate over lines
            input = line.split()           ##returns array of strings
            if(len(input) == 1):           ##if input array has only one element then is either:
                try:
                    self.natom = int(input[0])                                ##an int natom
                except ValueError:
                    if(input[0] == 'Angstrom' or input[0] == 'Bohr'):         ##or string units
                        self.units = input[0]
                    else:
                        print('Format error')                                 ##or neither? then error and exit
                        exit(0)
            elif(len(input) == 4):         ##if input array has four elements then process xyz:
                try:
                    self.charges[] += atomlib.get_charge(input[0])            ##is atomic symbol recognized in masses.py?
                    self.masses[] += atomlib.get_mass(input[0])               ##then assign charge and mass
                except KeyError:
                    print('Unrecognized atomic symbol')
                    exit(0)
                self.labels[] += input[0]                                     ##assign atomic symbol
                for i in range(1, 3):
                    try:
                        self.xyz[].append(float(input[i]))                    ##check and assign coords
                    except ValueError:
                        print('Format error in coordinate for %s'), input[0]
            else:                          ##if input array has unrecognized amount of elements, then error and exit
                print('Format error')
                exit(0)
        self.geom = np.matrix(self.xyz)                            ##make numpy array
        self.geom = np.reshape(self.geom, (3, self.natom))         ##reformat geom matrix to be 3 by number of atoms
        f.close()
