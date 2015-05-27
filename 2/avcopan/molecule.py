filepath       = '../../extra-files/molecule.xyz'
bohr2angstroms = 0.52917721092

def get_N():
  return int(open(filepath).readline())

def get_labels():
  return [ l.split()[0] for l in open(filepath).readlines()[2:] ]

def get_xyz():
  xyz =  []
  for l in open(filepath).readlines()[2:]:
    xyz += [float(x)/bohr2angstroms for x in l.split()[1:4]] 
  return xyz

