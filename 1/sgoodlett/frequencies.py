from geom import molecule
import numpy as np
from masses import get_mass

def frequencies(filename,units='Angstrom'):
    hess_list = []
    with open(filename,'r') as fn:
        for line in fn:
            hess_list.append(line.split())
    hess = np.array(hess_list).astype(float)

    mol = molecule(filename='molecule.xyz',units=units)
    mol.to_angstrom()
    mw_hess = hess.copy()
    shape = mw_hess.shape


    for i in range(len(mol.labels)):
        for j in range(len(mol.labels)):
            a1 = mol.labels[i]
            a2 = mol.labels[j]
            m1 = get_mass(a1)
            m2 = get_mass(a2)
            for k in range(3):
                for l in range(3):
                    mw_hess[(3*i+k)][(3*j+l)] = hess[(3*i+k)][(3*j+l)] / ((m1 * m2)**0.5)
    eig = np.linalg.eigh(mw_hess)

    mass_matrix = np.zeros(shape)
    mass_list = []

    for i in mol.labels:
        for j in range(3):
            mass_list.append(get_mass(i))

    for i in range(len(mass_list)):
        mass_matrix[i][i] = (mass_list[i])**(-0.5)

    q = np.dot(mass_matrix,eig[1])
    force_constants = (np.lib.scimath.sqrt(eig[0]*(4.35974e-18)*((1.89e10)**2)*(6.022e26)) / (2*np.pi)) / (3e10)
    fc_strings = []
    for i in range(len(force_constants)):
        if np.real(force_constants[i]) == 0.0:
            fc_strings.append('{}i'.format(np.imag(force_constants[i])))
        else:
            fc_strings.append('{}'.format(np.real(force_constants[i])))

    fout = open('normal_modes.xyz','w')
    for i in range(shape[0]):
        fout.write(str(mol.natom)+'\n'+'{}. Frequency: {} cm^-1'.format(i,fc_strings[i])+'\n')
        for j in range(len(mol.labels)):
            string = '{} {:10f} {:10f} {:10f} {:10f} {:10f} {:10f}'.format(mol.labels[j],mol.geom[j][0],mol.geom[j][1],mol.geom[j][2],
                                                   q[3*j,i],q[(3*j+1),i],q[(3*j+2),i])
            fout.write(string+'\n')
        fout.write('\n')
    fout.close()

if __name__ == "__main__":
    frequencies('hessian.dat',units='Bohr')