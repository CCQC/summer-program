import numpy as np

class tensor(np.ndarray):

  def __new__(cls, shape):
    return np.ndarray.__new__(cls, shape)

  def __mod__ (self, other):
    return self.dot(other)

  def __rmod__(self, other):
    return other.dot(self)
