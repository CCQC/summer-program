import numpy as np
import sys
sys.path.insert(0, '../../extra_files')
import masses as atomlib                   ##contains get_mass('x') and get_charge('x') methods

Class Molecule(object):
    def __init__:
        f = open('molecule.xyz', r)        ##file object molecule.xyz is f
	f.seek(-1)
        file_length = f.tell()
	f.seek(0)
	while(f.tell() < file_length):
	    input = ''
            while(f.read(1) != '\n'):      ##cat byte by byte to input string until '\n'
	        f.seek(-1, 1)
                input += f.read(1)
	    if (input.length == 0):        ##if blank move on
	    elif (input.length == 1):      ##if one byte try for natom
	        try:
	            int(input)
		except ValueError:
