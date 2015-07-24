#################################################################
# Numerical Hessian - written by quantumbutterfly 2015-07-24
# Last updated: 2015-07-24
#################################################################

import sys
sys.path.append("../../0/quantumbutterfly/")
import molecule from molecule
import numpy as np

#Displacements fixed to 0.005 Bohr
#1 directory for the reference config
#2*3N directories for single displacements
#
