#!/usr/bin/python
import os
from molecule import get_N, get_labels, get_xyz

def E(i,j,hi,hj):
  path = "x%dx%d_%d%d/output.dat" % (i,j,hi,hj)
  if os.path.isfile(path):
    for line in reversed( open(path).readlines() ):
      if line[:18] == "    Total Energy =":
        return float(line[19:])

###########################
# BEGIN COMPUTING HESSIAN #
###########################
h     = 0.005   # displacement size
N     = get_N() # number of atoms
H     = [[0.0 for i in range(3*N)] for j in range(3*N)]  # Hessian matrix -- initialize with zeros
# diagonal elements
for i in range(3*N):
  H[i][i] = ( E(i,0,+1,0)+E(i,0,-1,0)-2*E(0,0,0,0) ) / (h*h)

# off-diagonal elements
for i in range(3*N):
  for j in range(i):
    H[i][j] = ( E(i,j,+1,+1)+E(i,j,-1,-1)-E(i,0,+1,0)-E(i,0,-1,0)-E(j,0,+1,0)-E(j,0,-1,0)+2*E(0,0,0,0) ) / (2*h*h)
    H[j][i] = H[i][j]
#########################
# END COMPUTING HESSIAN #
#########################

outfile = open('hessian.txt','w')
for i in range(3*N):
  for j in range(3*N):
    outfile.write(" %12.7f" % H[i][j])
  outfile.write("\n")
