import numpy     as np
import itertools as it

class Numerov:

  def __init__(self, Vfunction):
    self.options   = {"mu":918.92865, "xmin":-1.0, "xmax":1.0, "pts":101, "maxiter": 100, "tol": 1e-7}
    self.Vfunction = Vfunction

  def init_state(self, state):
    xmax, xmin = self.options["xmax"], self.options["xmin"]
    self.pts   = self.options["pts"]
    self.h     = h = (xmax - xmin) / self.pts
    self.grid  = np.array([xmin + h * i for i in range(self.pts)])
    self.V     = np.array([self.Vfunction(x) for x in self.grid])
    self.Y     = np.ones(self.pts)
    self.Y[0]  = 0.0000
    self.Y[1]  = 0.0001
    minrange = (1.0 * E for E in it.count(0))
    maxrange = (1.0 * E for E in it.count(1))
    self.Emin, self.Emax = next(it.dropwhile(lambda (Emin, Emax): self.integrate_Y(Emax) < state + 1, it.izip(minrange, maxrange)))

  def integrate_Y(self, E):
    Y, h, m = self.Y, self.h, self.options["mu"]
    G = 2 * m * (self.V - E)
    nodecount = 0
    for i in range(2, self.pts):
      Y[i] =  ( 2*Y[i-1] - Y[i-2] + 5*G[i-1]*Y[i-1]*h**2/6. + G[i-2]*Y[i-2]*h**2/12. ) / (1 - G[i]*h**2/12.)
      if Y[i] * Y[i-1] < 0: nodecount += 1
    self.Y = Y
    self.nodecount = nodecount
    return nodecount

  def compute_state(self, state):
    self.init_state(state)
    Emax, Emin = self.Emax, self.Emin
    for i in range(self.options["maxiter"]):
      E = (Emax + Emin) / 2.
      nodecount = self.integrate_Y(E)
      if nodecount == state and abs(self.Y[-1]) < self.options["tol"]: break
      if   nodecount >  state:  Emax = E
      elif nodecount <= state:  Emin = E
    self.E = Emin
    return self.E




