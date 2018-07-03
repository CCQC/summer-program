import sys
sys.path.append('../../extra-files')
sys.path.append('../../1/mdav2')
sys.path.append('../../0/mdav2')
import os
import re
from subprocess import call
import numpy as np
import frequencies
import molecule

std_energy_prefix = '@DF-RHF Final Energy:'
mol = molecule.Molecule('../../extra-files/molecule.xyz')

def generate_inputs(mol, template, disp_size = 0.005, directory = 'DISPS'):
    #make directory and populate with all unique geometries
    if not os.path.exists(directory):
        os.makedirs(directory)
    with open(template,'r') as rawtemplate:
        temp = ''.join(rawtemplate.readlines())
    #create list of entries necessary
    entries = [(i,z) for i in range(0,3*mol.natom) for z in range(i,3*mol.natom)]
    with open(directory+'/XX','w') as origgeom:
        out = '\n'.join(mol.xyz_string().split('\n')[2:])
        origgeom.write(temp.format(out)) 
    for i in entries:
        tp = mol.copy()
        tm = mol.copy()
        posgeom = np.reshape(tp.geom,(9,1))
        neggeom = np.reshape(tm.geom,(9,1))
                
        if i[0] != i[1]:
            posgeom[i[0]]  += disp_size
            posgeom[i[1]]  += disp_size
            neggeom[i[0]]  += -1*disp_size
            neggeom[i[1]]  += -1*disp_size
        
        else:
            posgeom[i[0]] += disp_size
            neggeom[i[0]] += -1*disp_size
        
        tp.geom = posgeom.reshape(3,3)
        tm.geom = neggeom.reshape(3,3)
        
        with open(directory+'/'+str(i[0])+str(i[1])+'p.xyz','w') as positive:
            out ='\n'.join(tp.xyz_string().split('\n')[2:])
            positive.write(temp.format(out))
        with open(directory+'/'+str(i[0])+str(i[1])+'n.xyz','w') as negative:
            out = '\n'.join(tm.xyz_string().split('\n')[2:])
            negative.write(temp.format(out))
        
def run_jobs(mol,command='psi4',directory = 'DISPS'):
    #walk through directory and execute command for each file
    #works, but is there a way to do this with the importable psi4 module perhaps?
    files = os.listdir(directory)
    for f in files:
       
        call(["psi4",directory+'/'+f]) 

def build_hessian(mol, energy_prefix, disp_size = 0.005, directory = 'DISPS'):
    files = os.listdir(directory)
    energies = {}
    entries = [(i,z) for i in range(0,3*mol.natom) for z in range(i,3*mol.natom)]
    hessian = np.zeros((mol.natom * 3, mol.natom * 3))
    disp_bohr = disp_size/0.529177
    for f in files:
        if f.endswith('.dat'):
            with open(directory+'/'+f,'r') as data:
                datalines = data.readlines()
                for line in datalines:
                    if energy_prefix in line:
                        energies[f.strip('.xyz.dat')] = float(re.compile(r'[^\d.]+').sub('', line))
    
    for entry in entries:
        key = str(entry[0]) + str(entry[1])
        key0 = str(entry[0]) + str(entry[0])
        key1 = str(entry[1]) + str(entry[1])
        
        if entry[0] == entry[1]:
            e = (energies[key+'p']+energies[key+'n']-2*energies['XX'])/(disp_bohr ** 2)
            hessian[entry[0]][entry[1]] = e
        else:
            e = energies[key+'p']+energies[key+'n']
            e += -1*energies[key0+'p']-energies[key0+'n']
            e += -1*energies[key1+'p']-energies[key1+'n']
            e += 2*energies['XX']
            e = e/(2*(disp_bohr**2))
            hessian[entry[0]][entry[1]] = e 
            hessian[entry[1]][entry[0]] = e 
    return(hessian)

if __name__ == '__main__':
    templ = '../../extra-files/template.dat' 
    #generate_inputs(mol,templ)
    #run_jobs(mol)
    hessian = build_hessian(mol, std_energy_prefix)
    modes, freqs = frequencies.frequencies(mol, hessian)
    print(frequencies.printfrequencies(mol, freqs,modes))

