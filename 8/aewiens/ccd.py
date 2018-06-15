#!/usr/bin/env python3
import psi4, configparser, sys, numpy as np
sys.path.insert(0,"../../5/aewiens")
sys.path.insert(0,"../../6/aewiens")
from uhf import UHF
from mp2 import MP2

class CCD:

	def __init__(self,options):

		self.conv    = 10**( -int( options['CCD']['conv'] )) 
		self.tConv   = 10**( -int( options['CCD']['t_conv'] )) 
		self.maxiter = int( options['CCD']['max_iter'] )
		self.options = options

		uhf = UHF(options)
		mp2 = MP2(options)

		self.E0    = uhf.computeEnergy()
		self.e     = uhf.e
		self.Ec    = 0.0 
		self.nocc  = uhf.nelec
		self.nvirt = uhf.norb - uhf.nelec
		self.GMO   = mp2.transformTEI( uhf.G, uhf.C )
		self.uhf   = uhf

		##  initialize t-amplitudes  ##
		self.t = np.zeros((self.nocc,self.nocc,self.nvirt,self.nvirt))
    

	def transformTEI(self,g,C):
		return np.einsum('Pp,Pqrs->pqrs', C, \
			   np.einsum('Qq,PQrs->Pqrs', C,  \
			   np.einsum('Rr,PQRs->PQrs', C, \
			   np.einsum('Ss,PQRS->PQRs', C, g))))


	def getAmplitudes(self):
		g = self.GMO
		t  = self.t
		e  = self.e
		o  = slice(None,self.nocc)
		v  = slice(self.nocc,None)
		x  = np.newaxis
		Ep = 1.0/ (e[o,x,x,x] + e[x,o,x,x] - e[x,x,v,x] - e[x,x,x,v] )


		for k in range(self.maxiter):
			##  terms  ## t1 = g[o,o,v,v]
			t2 = np.einsum("abcd,ijcd->ijab",g[v,v,v,v],t)
			t3 = np.einsum("klij,klab->ijab",g[o,o,o,o],t)
			t4 = np.einsum("akic,jkbc->ijab",g[v,o,o,v],t)
			t5 = np.einsum("klcd,ijac,klbd->ijab",g[o,o,v,v],t,t)
			t6 = np.einsum("klcd,ikab,jlcd->ijab",g[o,o,v,v],t,t)
			t7 = np.einsum("klcd,ijcd,klab->ijab",g[o,o,v,v],t,t)
			t8 = np.einsum("klcd,ikac,jlbd->ijab",g[o,o,v,v],t,t)


			##   permutations  ## 
			t4P = t4 - t4.transpose((1,0,2,3)) - t4.transpose((0,1,3,2)) + t4.transpose((1,0,3,2))
			t5P = t5 - t5.transpose((0,1,3,2))
			t6P = t6 - t6.transpose((1,0,2,3))
			t8P = t8 - t8.transpose((1,0,2,3))


			##   update t2 amplitudes  ##
			t = t1 + 0.5*t2 + 0.5*t3 + t4P - 0.5*t5P - 0.5*t6P + 0.25*t7 + t8P
			t *= Ep 

		return t
    
	def computeEnergy(self):

		g = self.GMO
		nocc, t, e, = self.nocc, self.t, self.e

		o = slice(None,nocc)
		v = slice(nocc,None)
		x = np.newaxis
		Ep = 1.0/ (e[o,x,x,x] + e[x,o,x,x] - e[x,x,v,x] - e[x,x,x,v] )

		print('           Iter         Energy                 ΔE                   ‖ΔT‖')
		print('---------------------------------------------------------------------------')

		for k in range(self.maxiter):
			# terms
			t1 = g[o,o,v,v]
			t2 = np.einsum("abcd,ijcd->ijab",g[v,v,v,v],t)
			t3 = np.einsum("klij,klab->ijab",g[o,o,o,o],t)
			t4 = np.einsum("akic,jkbc->ijab",g[v,o,o,v],t)
			t5 = np.einsum("klcd,ijac,klbd->ijab",g[o,o,v,v],t,t)
			t6 = np.einsum("klcd,ikab,jlcd->ijab",g[o,o,v,v],t,t)
			t7 = np.einsum("klcd,ijcd,klab->ijab",g[o,o,v,v],t,t)
			t8 = np.einsum("klcd,ikac,jlbd->ijab",g[o,o,v,v],t,t)

			# permutations
			t4P = t4 - t4.transpose((1,0,2,3)) - t4.transpose((0,1,3,2)) + t4.transpose((1,0,3,2))
			t5P = t5 - t5.transpose((0,1,3,2))
			t6P = t6 - t6.transpose((1,0,2,3))
			t8P = t8 - t8.transpose((1,0,2,3))

			# update t2 amplitudes
			t = t1 + 0.5*t2 + 0.5*t3 + t4P - 0.5*t5P - 0.5*t6P + 0.25*t7 + t8P
			t *= Ep 

			# evaluate Ecorr
			Ec = 0.25 * np.sum(g[o,o,v,v] * t )
			dE = np.fabs(Ec - self.Ec)
			dT = np.fabs( np.linalg.norm(t) - np.linalg.norm(self.t) )

			self.Ec = Ec
			self.t = t

			print       ("@CCD-iter {:3d} {:20.12f} {:20.12f} {:20.12f}"  .format(k, self.E0 + Ec, dE, dT) )

			if dE < self.conv and dT < self.tConv:
				print('----------------------------------------------------------------------')
				print("\nCCD equations have converged.")
				print("Final CCD energy: {:20.12f}".format( self.E0 + self. Ec ) )
				print('----------------------------------------------------------------------')
				break

		return self.E0 + self.Ec

	"""
	def psiCCD(self):

		psi4.core.set_output_file("output.dat",False)

		psi4.set_options({'basis': self.uhf.basisName,
						'scf_type': 'pk',
						'reference': 'uhf',
						'e_convergence': 12,
						'r_convergence': 12,
						'puream': 0 })

		return psi4.energy('ccd')
	"""


if __name__ == '__main__':

	config = configparser.ConfigParser()
	config.read('Options.ini')
	ccd = CCD(config)
	ccd.computeEnergy()
