#!/usr/bin/python
import os
from molecule import get_N, get_labels, get_xyz

def run_input(dirname):
  os.chdir(dirname)
  os.system("psi4")
  print "running in ... " + os.getcwd()
  os.chdir("../")

run_input("x0x0_00")
for i in range(3*get_N()):
  run_input("x%dx0_10" % i)
  run_input("x%dx0_-10" % i)
for i in range(3*get_N()):
  for j in range(i):
    run_input("x%dx%d_11" %(i,j))
    run_input("x%dx%d_-1-1" %(i,j))
