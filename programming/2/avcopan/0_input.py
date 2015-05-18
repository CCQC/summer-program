#!/usr/bin/python
import os
from molecule import get_N, get_labels, get_xyz

template =  """
molecule h2o {
units bohr
%s
}

set basis sto-3g
set e_convergence 11
set d_convergence 9
set df_scf_guess false
set scf_type pk
energy('scf')
"""

def make_input(dirname,labels,coords):
  coordstring = ""
  for a in range(len(labels)):
    coordstring += ("%3s %10.7f %10.7f %10.7f\n" % (labels[a],coords[3*a+0],coords[3*a+1],coords[3*a+2]))
  if not os.path.exists(dirname):
    os.makedirs(dirname)
  f = open(dirname+'/input.dat','w+')
  f.write(template % coordstring)
  print dirname + template % coordstring


h = 0.005
N = get_N()
make_input("x0x0_00",get_labels(),get_xyz())
for i in range(3*N):
  disp_forward  = get_xyz()
  disp_backward = get_xyz()
  disp_forward[i]  += h
  disp_backward[i] -= h
  make_input("x%dx0_10"%i,  get_labels(), disp_forward)
  make_input("x%dx0_-10"%i, get_labels(), disp_backward)

for i in range(3*N):
  for j in range(0,i):
    disp_forward  = get_xyz()
    disp_backward = get_xyz()
    disp_forward[i]  += h
    disp_backward[i] -= h
    disp_forward[j]  += h
    disp_backward[j] -= h
    make_input("x%dx%d_11"  %(i,j), get_labels(), disp_forward)
    make_input("x%dx%d_-1-1"%(i,j), get_labels(), disp_backward)

