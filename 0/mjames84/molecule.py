import numpy as np

class Molecule(object):
    """
    A simple molecule
    """

    def __init__(self, geom_str, units = "Angstrom"):
        self.read(geom_string)
        self.units = units

    def read(self, geom_string):
        lines = geom_str.split("\n")
        num = int(lines[0])
        atoms = []
        coords = []
        for line in lines: 
     
