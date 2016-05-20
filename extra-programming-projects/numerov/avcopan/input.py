from numerov import Numerov

def V(x): return 0.5 * 0.399720901267 * x**2

num = Numerov(V)
num.options["tol"] = 1e-8
num.options["xmin"] = -1.0
num.options["xmax"] = +1.0
num.compute_state(0)

for i in range(80):
  print num.compute_state(i)
  print num.nodecount
