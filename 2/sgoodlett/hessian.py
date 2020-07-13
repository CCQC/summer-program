from geom import molecule
from frequencies import frequencies
import numpy as np
import os
import subprocess as sp
import re

class hessian():

    def __init__(self, mol, disp_size=0.005, directory='DISPS', command='psi4', energy_prefix = '@DF-RHF Final Energy:'):
        self.mol = mol
        self.disp_size = disp_size
        self.directory = directory
        self.command = command
        self.energy_prefix = energy_prefix
        self.cwd = os.getcwd()
        self.hes_directory = '{}/{}'.format(self.cwd,self.directory)
        if self.mol.units != 'Angstrom':
            mol.to_angstrom()
        with open('template.dat','r') as fn:
            self.template = fn.read()

    def write_file(self,mol,filename):
        with open(filename,'w') as fn:
            fn.write(self.template.format(mol.xyz_string(option=1)))

    def generate_inputs(self):
        try:
            os.mkdir(self.directory)
        except:
            pass
        os.chdir(self.hes_directory)
        try:
            os.mkdir('initial_geom')
        except:
            pass
        os.chdir('{}/initial_geom'.format(self.hes_directory))
        self.write_file(self.mol,'input.dat')
        os.chdir(self.hes_directory)
        #Single displacements
        for i in range(len(self.mol.labels)*3):
            try:
                os.mkdir(str(i)+'_p')
            except:
                pass
            os.chdir('{}/{}_p'.format(self.hes_directory,i))
            new_mol = self.mol.copy()
            flat = new_mol.geom.flatten()
            flat[i] = flat[i] + self.disp_size
            new_mol.geom = flat.reshape(len(new_mol.labels),3)
            self.write_file(new_mol,'input.dat')
            os.chdir(self.hes_directory)

            try:
                os.mkdir(str(i) + '_n')
            except:
                pass
            os.chdir('{}/{}_n'.format(self.hes_directory, i))
            new_mol = self.mol.copy()
            flat = new_mol.geom.flatten()
            flat[i] = flat[i] - self.disp_size
            new_mol.geom = flat.reshape(len(new_mol.labels), 3)
            self.write_file(new_mol, 'input.dat')
            os.chdir(self.hes_directory)

        #Double displacements
        for i in range(len(self.mol.labels)*3):
            for j in range(len(self.mol.labels)*3):
                if i < j:
                    try:
                        os.mkdir('{}_{}_pp'.format(i,j))
                    except:
                        pass
                    os.chdir('{}/{}_{}_pp/'.format(self.hes_directory, i, j))
                    new_mol = self.mol.copy()
                    flat = new_mol.geom.flatten()
                    flat[i] = flat[i] + self.disp_size
                    flat[j] = flat[j] + self.disp_size
                    new_mol.geom = flat.reshape(len(new_mol.labels), 3)
                    self.write_file(new_mol, 'input.dat')
                    os.chdir(self.hes_directory)

                    try:
                        os.mkdir('{}_{}_nn'.format(i,j))
                    except:
                        pass
                    os.chdir('{}/{}_{}_nn/'.format(self.hes_directory, i,j))
                    new_mol = self.mol.copy()
                    flat = new_mol.geom.flatten()
                    flat[i] = flat[i] - self.disp_size
                    flat[j] = flat[j] - self.disp_size
                    new_mol.geom = flat.reshape(len(new_mol.labels), 3)
                    self.write_file(new_mol, 'input.dat')
                    os.chdir(self.hes_directory)
        self.dirs = os.listdir()
        os.chdir(self.cwd)



    def run_jobs(self):
        os.chdir(self.hes_directory)
        for i in self.dirs:
            os.chdir('{}/{}'.format(self.hes_directory,i))
            #run execution command
            sp.run([self.command])
            os.chdir(self.hes_directory)
        os.chdir(self.cwd)

    def grab_e(self,e_dir):
        os.chdir(self.hes_directory)
        os.chdir(e_dir)
        with open('output.dat','r') as fn:
            for line in fn:
                if re.search(self.energy_prefix,line):
                    return float(line.split()[3])
        os.chdir(self.hes_directory)

    def build_hessian(self):
        a = len(self.mol.labels)*3
        #grab energies and parse to arrays
        initial_energy = self.grab_e('initial_geom')
        os.chdir(self.hes_directory)
        p_displacements = np.zeros((a,a))
        n_displacements = np.zeros((a,a))
        for i in range(a):
            for j in range(a):
                if i == j:
                    pd = self.grab_e('{}_p'.format(i))
                    nd = self.grab_e('{}_n'.format(i))
                    p_displacements[i][j] = pd
                    n_displacements[i][j] = nd
                elif i < j:
                    pp = self.grab_e('{}_{}_pp'.format(i,j))
                    nn = self.grab_e('{}_{}_nn'.format(i,j))
                    p_displacements[i][j] = pp
                    n_displacements[i][j] = nn
                else:
                    pass
        #Calc derivatives and parse to Hessian
        #Hessian is symmetric, so do not repeat comps but utilize symmetry accordingly
        hess = np.zeros((a,a))
        for i in range(a):
            for j in range(a):
                if i == j:
                    hess[i][j] = (p_displacements[i][j]+n_displacements[i][j]-(2*initial_energy))/(self.disp_size**2)
                elif i < j:
                    hess[i][j] = (p_displacements[i][j]+n_displacements[i][j]-p_displacements[i][i]-p_displacements[j][j]
                                  -n_displacements[i][i]-n_displacements[j][j]+(2*initial_energy))/(2*(self.disp_size**2))
                    hess[j][i] = hess[i][j]
                else:
                    pass
        os.chdir(self.cwd)
        np.savetxt('hessian_out.dat',hess,delimiter=' ')
        return hess


if __name__ == '__main__':
    #sp.run(['v','load','psi4@master'])
    mol = molecule(units='Bohr')
    hes = hessian(mol)
    hes.generate_inputs()
    hes.run_jobs()
    hes.build_hessian()
    f = frequencies('hessian_out.dat',units='Bohr')

