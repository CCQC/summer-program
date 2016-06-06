      #!/usr/bin/python 

import sys 
import masses as M
import numpy as np

class Molecule(object):
## Initialises the class and its member variables
   
   def __init__(self, units, natom, labels, masses, charges, geom):
      self.units = units #Done
      self.natom = natom #Done
      self.labels = labels #Done
      self.masses = masses #Done
      self.charges = charges #Done
      self.geom = geom #Done
   
           file = open ('molecule.xyz', 'r')

           natom = int(file.readline()) #These 2 lines populuate natom and units
           units = file.readline()

           if 'Bohr' in units:
              pass 
           elif 'Angstrom' in units:
              pass
           else:
              sys.exit('Your file specifies %s. Only Angstrom and Bohr are acceptable units of distance.' % units.split('\n')[0])
        #   elif units ==' ':
        #      units = 'Bohr'
        #   elif units == '':
        #      units = 'Bohr'
           
           labels = []             #These 4 lines create empty lists for the rest of the variables
           masses = []
           charges = []
           geom = []

           for i in range(natom):           #These lines populate the empty lists in 33-36
              temp = (file.readline().split())   
              labels.append(temp[0])
              geom.append(temp[1:4])
              charges.append(M.get_charge(labels[i]))
              masses.append(M.get_mass(labels[i]))
           print(self.geom)
           file.close()  

           print (masses, charges)

           def to_bohr(geom):
              if units == 'Bohr':
                 print ('Bohr units already employed.')
              elif units == 'Angstrom':
                 units = 'Bohr'
                 for i in range(len(geom)):           #This executes over the number of lists in geom
                    for j in range(3):                #This executes over the number of cartesian coordinations (3)
                       geom[i][j] = float(geom[i][j])/0.529177208     #This substitutes the jth index of the ith 
                       print ('Changed to Bohrs')                     #list of geom for its value in Bohrs

           def to_angstrom(geom):
              if units == 'Angstrom':
                 print ('Angstrom units already employed.')
              elif units == 'Bohr':
                 for i in range(len(geom)): #See lines 20-22 for explanation
                    for j in range(3):
                       geom[i][j] = float(geom[i][j])*0.529177208
                 print ('Changed to Angstroms')
                       
if __name__ == "__main__":
   mol = Molecule('molecule.xyz')
