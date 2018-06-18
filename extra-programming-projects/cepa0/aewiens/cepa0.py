#!/usr/bin/env python3
import psi4, configparser, sys, numpy as np
sys.path.insert(0,"../../../5/aewiens/")
sys.path.insert(0,"../../../6/aewiens/")
from uhf import UHF
from mp2 import MP2


class CEPA0:

	def __init__(self,options):

		self.conv    = 10**( -int( options['CEPA0']['conv'] )) 
		self.tConv   = 10**( -int( options['CEPA0']['t_conv'] )) 
		self.maxiter = int( options['CEPA0']['max_iter'] )
		self.options = options

		uhf = UHF(options)
		mp2 = MP2(options)

		self.E0    = uhf.computeEnergy()
		self.uhf   = uhf
		self.mp2   = mp2

		self.t = np.zeros(( uhf.nelec, uhf.nelec, uhf.nvirt, uhf.nvirt ))
		self.Ec = 0.0 


	def computeEnergy(self):

		uhf = self.uhf

		t = self.t
		g = self.mp2.transformTEI( uhf.G, uhf.C )
		o = slice( 0, uhf.nelec )
		v = slice( uhf.nelec, uhf.nelec + uhf.nvirt )
		x = np.newaxis
		Ep = 1.0 / ( uhf.e[o,x,x,x] + uhf.e[x,o,x,x] - uhf.e[x,x,v,x] - uhf.e[x,x,x,v] )

		print('       Iter         Energy                 ΔE                ‖ΔT‖')
		print('----------------------------------------------------------------------')

		for i in range(self.maxiter):

			t1 = g[o,o,v,v]
			t2 = np.einsum("abcd,ijcd->ijab",g[v,v,v,v],t)
			t3 = np.einsum("klij,klab->ijab",g[o,o,o,o],t)
			t4 = np.einsum("akic,jkbc->ijab",g[v,o,o,v],t)

			##  permute t4
			t4P = t4 - t4.transpose((1,0,2,3)) - t4.transpose((0,1,3,2)) + t4.transpose((1,0,3,2))
			t   = t1 + 0.5*t2 + 0.5*t3 + t4PA
			t  *= Ep

			Ec  = 0.25 * np.sum(g[o,o,v,v] * t )
			dE  = np.fabs(Ec - self.Ec)
			dT  = np.fabs( np.linalg.norm(t) - np.linalg.norm(self.t) )

			self.Ec = Ec
			self.t  = t
	
			print("@CEPA {:3d} {:20.12f}{:20.12f}{:20.12f}".format(i, self.E0 + self.Ec, dE, dT) )

			if dE < self.conv and dT < self.tConv: 
				print('----------------------------------------------------------------------')
				print("\nCEPA equations have converged.")
				print("Final CEPA0 energy: {:21.13f}".format( self.E0 + self. Ec ) )
				print("Compare to psi4:    {:21.13f}".format( self.psiCEPA0( self.options)) )
				break

		return self.E0 + self.Ec
	

	def psiCEPA0(self,options):

		psi4.core.set_output_file("output.dat",False)

		psi4.set_options({'basis': self.uhf.basisName,
						  'scf_type': 'pk',
						  'mp2_type': 'conv',
						  'reference': 'uhf',
						  'e_convergence': 12,
						  'r_convergence': 12,
						  'puream': 0 })

		return psi4.energy('lccd')


if __name__ == '__main__':

	config = configparser.ConfigParser()
	config.read('Options.ini')
	cepa = CEPA0(config)
	cepa.computeEnergy()
