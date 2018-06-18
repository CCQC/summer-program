import itertools as it

class P(object):

  def __init__(self, *axis_tuples, **kwargs):
    self.weight = 1.0 if not 'weight' in kwargs else kwargs['weight']
    self.axis_tuples = axis_tuples

  def iter_permutations(self, ndim):
    ref   = tuple(range(ndim))
    pools = []
    for ax in range(ndim):
      pool = (ax,)
      for tup in self.axis_tuples:
        if ax in tup: pool = tup
      pools.append(pool)
    for prod in it.product(*pools):
      if len(prod) == len(set(prod)):
        yield parity(ref, prod), prod

  def __mul__(self, other):
    if hasattr(other, "transpose") and hasattr(other, "ndim"):
      return self.weight * sum(sgn * other.transpose(per) for sgn, per in self.iter_permutations(other.ndim))
    elif isinstance(other, (float, int)):
      return P(*self.axis_tuples, weight = self.weight * other)
    else:
      raise Exception("Cannot left-multiply P permutation object with {:s}".format(type(other).__name__))

  def __rmul__(self, other):
    if isinstance(other, (float, int)):
      return P(*self.axis_tuples, weight = self.weight * other)
    else:
      raise Exception("Cannot right-multiply P permutation object with {:s}".format(type(other).__name__))
    

def parity(ref, per):
  sgn, per = +1, list(per)
  for elem in ref:
    i, j = ref.index(elem), per.index(elem)
    sgn *= -1 if not i is j else +1
    per[i], per[j] = per[j], per[i]
  return sgn
