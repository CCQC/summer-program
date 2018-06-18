#!/usr/bin/env python
import os, re, numpy as np
from molecule import Molecule
from processzmat import ProcessZmat

class Hessian(object):

	def __init__(self,mol,templateString,disp=0.005):

		self.mol = mol
		self.N   = len(self.mol)
		self.h   = disp

		zmat          = ProcessZmat(templateString)
		self.template = zmat.template

		self.dimH  = 3*self.N-6
		if self.N == 2:
			self.dimH = 1
	
	def makeInput(self,dirname,coords):

		Template  = self.template.format(*coords) 
		temp      = Molecule(Template,lengthUnits="Bohr",angleUnits="Radian")
		temp.toAngstrom()
		temp.toDegree()

		os.mkdir("{:s}".format(dirname))
		f = open("{:s}/input.dat".format(dirname),"w")
		f.write( self.template.format(*temp.coords) )
		f.close()


	def runInput(self,dirname):
		
		os.chdir(dirname)
		os.system("psi4")
		os.chdir("..")


	def find_E(self,i,j,hi,hj):
		"""
		Pull energy from output.dat and throw an error if not found
		:params i,j: indices of atoms 1,2
		:params hi,hj: displacements of atoms 1,2 (-1, 0, or 1, corresponds to -h, 0, or h)
		"""
		dirname = "Q%dQ%d_%d%d" % (i,j,hi,hj)
		out_str = open("HESS/%s/output.dat" % dirname, "r").read()
		match = re.findall("Total Energy\s=\s+-\d+.\d+",out_str)
		if match == []:
			out = "Cannot find energy!"
		else:
			out = float(match[0].split()[-1])
		return out


	def runDisps(self):
		
		h,N = self.h, self.N
		os.mkdir("HESS")
		os.chdir("HESS")

		##  run reference configuration  ##
		self.makeInput("Q0Q0_00",self.mol.coords)
		self.runInput("Q0Q0_00")

		##  run single displacements  ##
		for i in range(self.dimH):

			forward   = "Q%dQ0_10" % i
			reverse   = "Q%dQ0_-10" % i
			coordCopy = self.mol.copy().coords
			coordCopy[i] += h
	
			self.makeInput(forward,coordCopy)
			self.runInput(forward)

			coordCopy[i] -=2*h
			self.makeInput(reverse,coordCopy)
			self.runInput(reverse)
		
		##  run double displacements  ##	

		for i in range(self.dimH):
			for j in range(i):
				forward    = "Q{:d}Q{:d}_11".format(i,j)
				reverse    = "Q{:d}Q{:d}_-1-1".format(i,j)
				coordCopy2 = self.mol.copy().coords
				
				coordCopy2[i] += h
				coordCopy2[j] += h

				self.makeInput(forward,coordCopy2)
				self.runInput(forward)

				coordCopy2[i] -= 2*h
				coordCopy2[j] -= 2*h

				self.makeInput(reverse,coordCopy2)
				self.runInput(reverse)

		os.chdir("..")

	def makeHessian(self):

		self.runDisps()

		h, N = self.h, self.N
		E0 = self.find_E(0,0,0,0)
		self.H = np.zeros((self.dimH,self.dimH))

		for i in range(self.dimH):
			self.H[i,i]= (self.find_E(i,0,1,0)+self.find_E(i,0,-1,0)-2*E0)/(h**2)
			for j in range(0,i):
				self.H[i,j] = (self.find_E(i,j,1,1)+self.find_E(i,j,-1,-1)-self.find_E(i,0,1,0)-self.find_E(j,0,1,0)-self.find_E(j,0,-1,0)-self.find_E(i,0,-1,0)+2*E0)
				self.H[i,j] /= 2*h**2
				self.H[j,i] = self.H[i,j]


	def writeHessian(self):
		"""
		write Hessian matrix to hessian.dat file
		"""
		self.makeHessian()
		np.savetxt("hessian.dat",self.H,"%15.7f"," ","\n")
		
if __name__ == '__main__':
	f    = open("template.dat","r").read()
	mol  = Molecule(f)
	mol.toBohr()
	mol.toRadian()
	test = Hessian(mol,f)
	test.writeHessian()

