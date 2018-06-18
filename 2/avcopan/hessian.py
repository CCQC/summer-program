#!/usr/bin/python
import os
import numpy as np
import re

def generate_inputs(mol, template, disp_size = 0.005, directory = "DISPS"):

  mol = mol.copy()
  mol.to_bohr()

  def write_input(coord1, coord2, k1, k2, disp_mol):
    path = directory + "/x{:d}x{:d}_{:d}{:d}".format(coord1, coord2, k1, k2)
    if not os.path.isdir(path): os.makedirs(path)
    inpf = open(path + "/input.dat", 'w+')
    inpf.write(template.format(str(disp_mol)))
    inpf.close()

  # origin
  write_input(0, 0, 0, 0, mol)

  # displacements
  for coord1 in range(3*mol.natom):
    for k in [-1, +1]:
      # single displacement
      disp1          = np.zeros(3*mol.natom)
      disp1[coord1] += k * disp_size
      disp1_mol      = mol.copy()
      disp1_mol.displace(disp1.reshape(mol.natom,3))
      write_input(coord1, 0, k, 0, disp1_mol)

      for coord2 in range(coord1):
        # double displacement
        disp2          = np.zeros(3*mol.natom)
        disp2[coord2] += k * disp_size
        disp2_mol      = disp1_mol.copy()
        disp2_mol.displace(disp2.reshape(mol.natom,3))
        write_input(coord1, coord2, k, k, disp2_mol)


def run_jobs(mol, command = "psi4", directory = "DISPS"):
  natom = mol.natom
  wd    = os.getcwd()

  def run_input(coord1, coord2, k1, k2):
    os.chdir(directory + "/x{:d}x{:d}_{:d}{:d}".format(coord1, coord2, k1, k2))
    print("running ..." + (2 * " coord {:d} displaced by {:+d}*disp_size").format(coord1, k1, coord2, k2))
    os.system(command)
    os.chdir(wd)

  run_input(0, 0, 0, 0)
  for coord1 in range(3*mol.natom):
    for k in [-1, +1]:
      run_input(coord1, 0, k, 0)
      for coord2 in range(coord1):
        run_input(coord1, coord2, k, k)


def build_hessian(mol, energy_prefix, disp_size = 0.005, directory = "DISPS"):
  ncart = 3 * mol.natom
  h = disp_size

  def E(coord1, coord2, k1, k2):
    coord1 = coord1 if not k1 is 0 else 0
    coord2 = coord2 if not k2 is 0 else 0
    if coord2 > coord1 or k1 is 0: coord1, coord2, k1, k2 = coord2, coord1, k2, k1
    path = directory + "/x{:d}x{:d}_{:d}{:d}".format(coord1, coord2, k1, k2)
    try:
      outf = open(path + "/output.dat").read()
      for match in re.finditer(energy_prefix + '\s+(-?\d+\.\d+)', outf): pass
      energy = float(match.group(1))
      return energy
    except:
      raise Exception("Failed to parse energy in {:s}".format(path))

  H = np.zeros((ncart, ncart))

  for A in range(ncart):
    for B in range(ncart):
      if A == B:
        H[A,A] = (E(A,A,+1, 0) + E(A,A,-1, 0) - E(A,A, 0, 0) * 2) / (h**2)
      else:
        H[A,B] = (E(A,B,+1,+1) + E(A,B,-1,-1) - E(A,B,+1, 0) - E(A,B,-1, 0) - E(A,B, 0,+1) - E(A,B, 0,-1) + 2 * E(A,B, 0, 0)) / (2 * h**2)

  np.savetxt("hessian.dat", H)

  return H
  


if __name__ == "__main__" :
  import sys
  sys.path.insert(0, "../../0/avcopan/")
  sys.path.insert(0, "../../1/avcopan/")

  from molecule_advanced import Molecule
  from frequencies       import frequencies

  # 1. read in molecule object
  mol = Molecule(open("../../extra-files/molecule.xyz").read())

  # 2. read in template file
  template = open("../../extra-files/template.dat").read()

  # call input file generator
  generate_inputs(mol, template)

  # execute jobs
  run_jobs(mol)

  # build hessian
  H = build_hessian(mol, '@DF-RHF Final Energy:')

  v, Q = frequencies(mol, H)
