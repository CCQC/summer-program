#__author__ = "mbowman"

import STO3Gbasis as sto
import molecule
import masses
import math
import numpy as np
from scipy import special, misc
import collections

np.set_printoptions(threshold=np.inf, precision=3, linewidth=200, suppress=True)

class Integral(object):
	
	def __init__(self, mol_str):
		"""
		""" 
		self.mol = molecule.Molecule(mol_str)
		self.mol.to_bohr()
		self.atoms = self.mol.labels
		self.coord = np.array(self.mol.geom)
		self.bd = self.createBasisDict() 
		self.K = len(self.bd) #dimension of molecular integrals		
		self.VNuc = self.nuclearRepulsionEnergy()
		print self.VNuc
		self.S = self.overlapIntegral()
		self.T = self.kineticEnergyIntegral()
		self.V = self.nuclearAttractionIntegral()
		self.G = self.electronRepulsionIntegral()
		
	def subshells(self, z):
		s = ['1s']
		if z > 2:
			s.extend(['2s','2px', '2py', '2pz'])
			if z > 10:
				s.extend(['3s','3px', '3py', '3pz']) 
		return s

	def createBasisDict(self):
		bd = {} 
		duplicates = [item for item, count in collections.Counter(self.atoms).items() if count > 1]
		duplicateCounter = [1] * len(duplicates)
		for i, atom in enumerate(self.atoms):
			tempCharge = masses.get_charge(atom) 	#atomic number of atom, determines subshells
			if atom in duplicates:
				tas = atom +  str(duplicateCounter[duplicates.index(atom)]) 	#creates a string to denote duplicate atoms e.g. H_2
				duplicateCounter[duplicates.index(atom)] += 1
			else:
				tas = atom							#creates a string to denote unique atoms e.g. O
			tempShells = self.subshells(tempCharge) 	#creates temp array of subshells from atom
			for subshell in tempShells:
				shellInfo = [c for c in subshell] 	#splits shell into three sections based on three quantum numbers
				sss = tas + '_' + subshell
				bd[sss] = {} 		#initializes dictionary entry for subshell of atom
				bd[sss]['a'] = sto.expFact(tempCharge, int(shellInfo[0])) #assigns expansion factors based on atomic number and shell
				bd[sss]['d'] = sto.contCoeff(int(shellInfo[0]),(2 if shellInfo[1] == 'p' else 1))
				bd[sss]['z'] = tempCharge		#atomic mass of atom
				bd[sss]['c'] = 0 if not atom in duplicates else duplicateCounter[duplicates.index(atom)] - 1 #number of atom ie H1 vs H2
				bd[sss]['s'] = shellInfo[0]		#orbital shell
				if shellInfo[1] == 'p':	
					bd[sss]['ss'] = 2		#orbital subshell
					bd[sss]['l'] = 1 if shellInfo[2] == 'x' else 0 #orbital angular momentum
					bd[sss]['m'] = 1 if shellInfo[2] == 'y' else 0
					bd[sss]['n'] = 1 if shellInfo[2] == 'z' else 0
				else:
					bd[sss]['ss'] = 1
					bd[sss]['l'], bd[sss]['m'], bd[sss]['n'] = 0, 0, 0
				bd[sss]['R'] = [self.coord[i][c] for c in range(3)] 
		
		return collections.OrderedDict(sorted(bd.items(), key=lambda t: (-t[1]['z'], t[1]['c'], t[1]['s'], t[1]['ss'], -t[1]['l'],-t[1]['m'],-t[1]['n'] )))

	def nuclearRepulsionEnergy(self):
		ts = 0
		for i, atom1 in enumerate(self.atoms):
			for j, atom2 in enumerate(self.atoms[i+1:]):
				if(i != j):
					ts += masses.get_charge(atom1) * masses.get_charge(atom2) / np.linalg.norm(self.coord[j] - self.coord[i])
		return ts

	def overlapIntegral(self):			 
		S = np.zeros((self.K, self.K))
		print "Calculating overlap integral (S):"
		for i, orb1 in enumerate(self.bd):
			for j, orb2 in enumerate(self.bd):
				if i > j:
					S[i][j] = S[j][i] #Turnover rule guarantees these matrices are symmetric same with T, V
				else:
					ts = 0
					for alphaA, da in zip(self.bd[orb1]['a'], self.bd[orb1]['d']):
						for alphaB, db in zip(self.bd[orb2]['a'], self.bd[orb2]['d']):
							ts += da* db * self.normalization(alphaA, self.bd[orb1]['l'], self.bd[orb1]['m'],self.bd[orb1]['n'] ) * self.normalization(alphaB, self.bd[orb2]['l'], self.bd[orb2]['m'], self.bd[orb2]['n']) * self.f3(orb1, orb2, alphaA, alphaB, da, db, 0, 0, 0) 
					S[i][j] = ts
		print "\tDone"
		return S	
	
	def kineticEnergyIntegral(self):
		T = np.zeros((self.K, self.K))
		print "Calculating kinetic energy integral (T):"
		for i, orb1 in enumerate(self.bd):
			for j, orb2 in enumerate(self.bd):
				if i > j:
					T[i][j] = T[j][i] #Turnover rule	
				else:	
					ts = 0
					for alphaA, da in zip(self.bd[orb1]['a'], self.bd[orb1]['d']):
						for alphaB, db in zip(self.bd[orb2]['a'], self.bd[orb2]['d']):
							ss  = alphaB * (2 * (self.bd[orb2]['l'] + self.bd[orb2]['m'] + self.bd[orb2]['n']) + 3) * self.f3(orb1, orb2, alphaA, alphaB, da, db, 0, 0, 0)
							ss += -2.0 * (alphaB ** 2) * self.f3(orb1, orb2, alphaA, alphaB, da, db, 2, 0, 0) 
							ss += -2.0 * (alphaB ** 2) * self.f3(orb1, orb2, alphaA, alphaB, da, db, 0, 2, 0)
							ss += -2.0 * (alphaB ** 2) * self.f3(orb1, orb2, alphaA, alphaB, da, db, 0, 0, 2)
							ss += -0.5 * self.bd[orb2]['l'] * ( self.bd[orb2]['l'] - 1 ) * self.f3(orb1, orb2, alphaA, alphaB, da, db, -2, 0, 0)
							ss += -0.5 * self.bd[orb2]['m'] * ( self.bd[orb2]['m'] - 1 ) * self.f3(orb1, orb2, alphaA, alphaB, da, db, 0, -2, 0)
							ss += -0.5 * self.bd[orb2]['n'] * ( self.bd[orb2]['n'] - 1 ) * self.f3(orb1, orb2, alphaA, alphaB, da, db, 0, 0, -2)
							ss *= da* db * self.normalization(alphaA, self.bd[orb1]['l'], self.bd[orb1]['m'],self.bd[orb1]['n'] ) * self.normalization(alphaB, self.bd[orb2]['l'], self.bd[orb2]['m'], self.bd[orb2]['n']) 
							ts += ss
					T[i][j] = ts
		print "\tDone"
		return T
	
	def nuclearAttractionIntegral(self):
		V = np.zeros((self.K, self.K))
		print "Calculating nuclear attraction integral (V):"
		for atom, c in zip(self.atoms, self.coord):
			tempV = np.zeros((self.K, self.K))
			for a, orb1 in enumerate(self.bd):
				for b, orb2 in enumerate(self.bd):
					if a > b:
						tempV[a][b] = tempV[b][a] #Turnover rule
					else:
						ts = 0.0 #total sum
						r1 = np.array(self.bd[orb1]['R']) #coordinates of atom of orb1
                                                r2 = np.array(self.bd[orb2]['R']) #coordinates of atom of orb2
						norm = (r1[0] - r2[0])**2  + (r1[1] - r2[1])**2 + (r1[2] - r2[2])**2 #norm square of vector AB
						la, ma, na, lb, mb, nb = self.bd[orb1]['l'], self.bd[orb1]['m'], self.bd[orb1]['n'], self.bd[orb2]['l'], self.bd[orb2]['m'], self.bd[orb2]['n']
						for alphaA, da in zip(self.bd[orb1]['a'], self.bd[orb1]['d']):
							for alphaB, db in zip(self.bd[orb2]['a'], self.bd[orb2]['d']):
								g = alphaA + alphaB
								p = self.p(r1, r2, alphaA, alphaB) 
								cnorm = (p[0] - c[0])**2 + (p[1] - c[1])**2 + (p[2] - c[2])**2 #norm square of vector PC 
								ss = da* db * self.normalization(alphaA, la, ma, na ) * self.normalization(alphaB, lb, mb, nb) #second sum
								ss *= -masses.get_charge(atom)*2*math.pi*math.exp(-1*alphaA*alphaB*norm/g)/g
								ss *= (math.fsum([math.fsum([math.fsum([self.nu(l,r,i,la,lb,r1[0],r2[0],c[0],p[0],g)*
								math.fsum([math.fsum([math.fsum([self.nu(m,s,j,ma,mb,r1[1],r2[1],c[1],p[1],g)*
								math.fsum([math.fsum([math.fsum([self.nu(n,t,k,na,nb,r1[2],r2[2],c[2],p[2],g)*
								self.Boys(l+m+n-2*(r+s+t)-(i+j+k),g*cnorm) 
								for k in range(int((n-2*t)/2)+1)]) for t in range(int(n/2)+1)]) for n in range(na+nb+1)]) 
								for j in range(int((m-2*s)/2)+1)]) for s in range(int(m/2)+1)]) for m in range(ma+mb+1)]) 
								for i in range(int((l-2*r)/2)+1)]) for r in range(int(l/2)+1)]) for l in range(la+lb+1)])) 
								ts += ss
						tempV[a][b] = ts
			V += tempV
		print "\tDone"
		return V
	
	def electronRepulsionIntegral(self):
		G = np.zeros((self.K,self.K,self.K,self.K))
		print "Calculating electron repulsion integral (G):"
		for a, orb1 in enumerate(self.bd):
			for b, orb2 in enumerate(self.bd):
				for c, orb3 in enumerate(self.bd):
					for d, orb4 in enumerate(self.bd):
						if c > d:
							pass #These values are already calculated :)	
						elif a > b:
							pass 
						elif a > c and b > d:
							pass 
						elif a > d and b > c:
							pass 
						else:
							ts = 0
							r1, r2, r3, r4 = np.array(self.bd[orb1]['R']), np.array(self.bd[orb2]['R']), np.array(self.bd[orb3]['R']), np.array(self.bd[orb4]['R'])
							norm1, norm2 = ((r1[0]-r2[0])**2 + (r1[1]-r2[1])**2 + (r1[2]-r2[2])**2), ((r3[0]-r4[0])**2 + (r3[1]-r4[1])**2 + (r3[2]-r4[2])**2)
							la, lb, lc, ld = self.bd[orb1]['l'], self.bd[orb2]['l'], self.bd[orb3]['l'], self.bd[orb4]['l']
                                                        ma, mb, mc, md = self.bd[orb1]['m'], self.bd[orb2]['m'], self.bd[orb3]['m'], self.bd[orb4]['m']
                                                        na, nb, nc, nd = self.bd[orb1]['n'], self.bd[orb2]['n'], self.bd[orb3]['n'], self.bd[orb4]['n']
							for alphaA, da in zip(self.bd[orb1]['a'], self.bd[orb1]['d']):
								for alphaB, db in zip(self.bd[orb2]['a'], self.bd[orb2]['d']):
									for alphaC, dc in zip(self.bd[orb3]['a'], self.bd[orb3]['d']):
										for alphaD, dd in zip(self.bd[orb4]['a'], self.bd[orb4]['d']):
											g1, g2 = (alphaA + alphaB), (alphaC + alphaD) #equation 20
											delta = (g1 + g2)/(4*g1*g2) #equation 21
											p, q = self.p(r1,r2,alphaA,alphaB), self.p(r3,r4,alphaC,alphaD)
											pqnorm = (p[0] - q[0])**2 + (p[1] - q[1])**2 + (p[2] -q[2])**2
											ss = (da*db*dc*dd*self.normalization(alphaA,la,ma,na)*self.normalization(alphaB,lb,mb,nb)*
											self.normalization(alphaC,lc,mc,nc)*self.normalization(alphaD,ld,md,nd))
											ss *= 2 * math.pi**2 *math.sqrt(math.pi/(g1+g2))/(g1*g2)
											ss *= math.exp(-1*alphaA*alphaB*norm1/g1)*math.exp(-1*alphaC*alphaD*norm2/g2)
											ss *= (math.fsum([self.f4(l,lPrime,r,rPrime,i,delta,la,lb,r1[0],r2[0],p[0],g1,lc,ld,r3[0],r4[0],q[0],g2)*
											math.fsum([self.f4(m,mPrime,s,sPrime,j,delta,ma,mb,r1[1],r2[1],p[1],g1,mc,md,r3[1],r4[1],q[1],g2)*
											math.fsum([self.f4(n,nPrime,t,tPrime,k,delta,na,nb,r1[2],r2[2],p[2],g1,nc,nd,r3[2],r4[2],q[2],g2)*
											self.Boys(l+lPrime+m+mPrime+n+nPrime-2*(r+rPrime+s+sPrime+t+tPrime)-(i+j+k),pqnorm/(4*delta))
											for n in range(na+nb+1) for t in range(int(n/2)+1) for nPrime in range(nc+nd+1) for tPrime in range(int(nPrime/2)+1) for k in range(int((n+nPrime-2*(t+tPrime))/2)+1)])
											for m in range(ma+mb+1) for s in range(int(m/2)+1) for mPrime in range(mc+md+1) for sPrime in range(int(mPrime/2)+1) for j in range(int((m+mPrime-2*(s+sPrime))/2)+1)])
											for l in range(la+lb+1) for r in range(int(l/2)+1) for lPrime in range(lc+ld+1) for rPrime in range(int(lPrime/2)+1) for i in range(int((l+lPrime-2*(r+rPrime))/2)+1)]))
											ts += ss
							G[a][b][c][d]=G[c][d][a][b]=G[b][a][d][c]=G[d][c][b][a]=G[b][a][c][d]=G[d][c][a][b]=G[a][b][d][c]=G[c][d][b][a] = ts	
							# These permutations are all equal if the orbitals are real
			print "\t{:3.1f}% Complete".format(100*(a*len(self.bd)+b)/(len(self.bd)**2))
		print "\tDone"			
		return G	
		 
	def f1(self, j, l, m, a, b):
		"""
		corresponds to equation 8
		"""
		s = 0
		for k in range(max(0,(j-m)), min(j, l) + 1):
			s += special.comb(l,k)*special.comb(m,j-k)*math.pow(a,l - k)*math.pow(b,m + k - j)
			#print "      k: " + str(s)	
		return s

	def p(self, a, b, alphaA, alphaB):
		"""
		corresponds to equation 5
		"""
		return ( np.multiply(alphaA, a) + np.multiply(alphaB, b) )/ (alphaA + alphaB)

	def f2(self, gamma, l, m, a, b, p):
		"""
		corresponds to equation 3
		""" 
		s = self.f1(0, l, m, (p-a) ,(p-b))   
		#print "   s: " + str(s)
		for x in range(1, int( (l + m) / 2 ) +1 ):
			s += (self.f1(2*x, l, m, (p-a) ,(p-b) ) * misc.factorial2((2*x-1),exact=True) / ((2* gamma) ** x)) 
			#print "   s: " + str(s)
		s = s * math.sqrt(( math.pi / gamma) )
		return s	
	
	def f3(self, orb1, orb2, alphaA, alphaB, da, db, ls, ms, ns):
		"""
		corresponds to equation 1 and 2
		"""
		if any((i < 0 for i in {(self.bd[orb2]['l'] + ls) + (self.bd[orb2]['m'] + ms) + (self.bd[orb2]['n'] + ns)})):
			return 0

		gamma = alphaA + alphaB
                p = self.p(self.bd[orb1]['R'], self.bd[orb2]['R'], alphaA, alphaB)
		norm = (self.bd[orb1]['R'][0] - self.bd[orb2]['R'][0])**2  + (self.bd[orb1]['R'][1] - self.bd[orb2]['R'][1])**2 + (self.bd[orb1]['R'][2] - self.bd[orb2]['R'][2])**2
		t = math.exp( -1*alphaA*alphaB*norm / gamma  )
		t *= self.f2(gamma, self.bd[orb1]['l'], self.bd[orb2]['l'] + ls, self.bd[orb1]['R'][0], self.bd[orb2]['R'][0], p[0])
		t *= self.f2(gamma, self.bd[orb1]['m'], self.bd[orb2]['m'] + ms, self.bd[orb1]['R'][1], self.bd[orb2]['R'][1], p[1])
		t *= self.f2(gamma, self.bd[orb1]['n'], self.bd[orb2]['n'] + ns, self.bd[orb1]['R'][2], self.bd[orb2]['R'][2], p[2])
		#print ""
		return t

	def nu(self, l, r, i, la, lb, a, b, c, p, gamma):
		"""
		corresponds to triply nested sums in equation 13/ 14
		"""
		ts = ((-1) ** (l+ i)) * self.f1(l, la, lb, (p - a), (p - b)) * misc.factorial(l, exact=True) * ( (p - c) ** (l - 2*r - 2*i))
		ts /= misc.factorial(r, exact=True) * misc.factorial(i, exact=True) * misc.factorial(l - 2*r - 2*i, exact=True) * ( (4*gamma) ** (r + i)) 
		return ts	
		
	def Boys(self, nu, x):
		"""
		corresponds to equation 16 and 17
		"""
		if x < 0.000001:
			return ( 1.0 / (2 * nu + 1)) - ( x / (2 * nu + 3))
		else:
			return 0.5 * x**(-(nu + 0.5)) * special.gamma(nu + 0.5) * special.gammainc(nu + 0.5, x)  
	
	def theta(self, l, la, lb, a, b, r, g):
		return self.f1(l,la,lb,a,b) * misc.factorial(l,exact=True) * g ** (r - l) / (misc.factorial(r,exact=True)*misc.factorial(l-2*r,exact=True))

	def f4(self,l,lp,r,rp,i,delta,la,lb,a,b,p,g1,lc,ld,c,d,q,g2):
		ts = (-1)**l * self.theta(l,la,lb,p-a,p-b,r,g1) * self.theta(lp,lc,ld,q-c,q-d,rp,g2)
		ts *= (-1)**i * (2*delta)**(2*(r+rp)) * misc.factorial(l+lp-2*(r+rp),exact=True) * delta**i * (p-q)**(l+lp-2*(r+rp+i))
		ts /= (4*delta)**(l+lp) * misc.factorial(i,exact=True) * misc.factorial(l+lp-2*(r+rp+i),exact=True)
		return ts	
	
	def normalization(self, alpha, l, m, n):
		"""
		corresponds to equation 9
		"""
		s = (4 * alpha) ** (l + m + n)
		s /= ( misc.factorial2( (2*l -1), exact=True) * misc.factorial2( (2*m -1), exact=True) * misc.factorial2( (2*n -1), exact=True) )
		s *= math.pow( (2*alpha/math.pi), 1.5) 
		return math.sqrt(s)			

	def printIntegral(self, I, name):
		"""
		returns Overlap integral in neat tabular form
		"""
		l = "+-------+"
		for i in range(len(self.bd)):
			l += "-------+"
		l += "\n"
		s = l
		s += "|   " + name + "   |"
		for orb in self.bd:
			s += orb.ljust(7) + "|"
		s += "\n" + l  
		for i, orb1 in enumerate(self.bd):
			s += "|" + orb1.ljust(7) + "|"
			for j, orb2 in enumerate(self.bd):
				s += "  1.0  |" if np.absolute(I[i][j]-1) < 0.0000001 else "  0.0  |" if np.absolute(I[i][j]) < 0.0000001 else "{:>7.3f}|".format(I[i][j])
		 	s += "\n" + l

		return s
				
