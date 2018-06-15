import numpy     as np
import itertools as it

class P(object):

  def __init__(self, axis_strings_with_bars, **kwargs):
    """
    Permutation object constructor
    :param axis_string: a string of the form "0,1/2,3|4,5,6/7|..." defining a permutation of axes
           according to Bartlett and Shavitt
    """
    self.axis_strings_with_bars = axis_strings_with_bars
    self.weight = 1.0 if not 'weight' in kwargs else kwargs['weight']
    axis_strings = axis_strings_with_bars.split('|')
    element_sets, compositions = zip(*(P.process_axis_string(axis_string) for axis_string in axis_strings))
    self.block_permutations = [BlockPermutations(element_set, composition) for element_set, composition in zip(element_sets, compositions)]
    self.elements = sum(element_sets, ())

  def iter_permutations_with_signature(self, ndim):
    axes = tuple(range(ndim))
    pools = tuple(block_permutation.iter_permutations_with_signature() for block_permutation in self.block_permutations)
    for prod in it.product(*pools):
      signs, permuted_subsets = zip(*prod)
      sign = np.product(signs)
      permuted_elements = sum(permuted_subsets, ())
      permuted_axes = self.permute(permuted_elements)(axes)
      yield sign, permuted_axes

  def __mul__(self, other):
    if hasattr(other, "transpose") and hasattr(other, "ndim"):
      return self.weight * sum(sgn * other.transpose(per) for sgn, per in self.iter_permutations_with_signature(other.ndim))
    elif isinstance(other, (float, int)):
      return P(self.axis_strings_with_bars, weight = self.weight * other)
    else:
      raise Exception("Cannot left-multiply P permutation object with {:s}".format(type(other).__name__))

  def __rmul__(self, other):
    if isinstance(other, (float, int)):
      return P(self.axis_strings_with_bars, weight = self.weight * other)
    else:
      raise Exception("Cannot right-multiply P permutation object with {:s}".format(type(other).__name__))

  def permute(self, permuted_elements):
    return lambda x: tuple(permuted_elements[self.elements.index(item)] if item in self.elements else item for item in x)

  @staticmethod
  def process_axis_string(axis_string):
    try:
      axis_equivalent_sets = [[int(axis) for axis in substring.split(',')] for substring in axis_string.split('/')]
      elements = tuple(sum(axis_equivalent_sets, []))
      composition = tuple(len(axis_equivalent_set) for axis_equivalent_set in axis_equivalent_sets)
      return elements, composition
    except:
      raise Exception("Invalid string {:s} passed as constructor argument.".format(axis_string))



class BlockPermutations(object):
  """
  Given a list of elements, iterate over all unique permutations, treating subsets of the elements as equivalent.
  This is done by mapping equivalent elements into a single element, applying all unique permutations to that,
  and applying the same changes to the original set.
  Uses Algorithm L from Knuth 'The Art of Computer Programming: Volume 4A: Pre-Fascicle 2B: Draft of
  Section 7.1.1.2 - Generating All Permutations' which can be found at http://www.cs.utsa.edu/~wagner/knuth/
  """

  def __init__(self, elements, composition):
    """
    :param elements: an iterable; for example ('a','b','c','d')
    :param composition: an integer composition of the number of elements, indicating which subsets of elements
           should be treated as equivalent; with the above example, (2,2) would indicate that ('a','b') and
           ('c','d') should be treated as equivalent.
    """
    self.elements = tuple(elements)
    self.nelem = len(elements)
    self.indices = sum((block_size * [ block_index ] for block_index, block_size in enumerate(composition)), [])

  def iter_permutations_with_signature(self):
    elements = list(self.elements)
    indices  = list(self.indices)
    sgn = 1
    p0 = 0
    while p0 >= 0:
      yield sgn, tuple(elements)
      p0 = next((p for p in reversed(range(self.nelem-1)) if indices[p ] < indices[p+1]), -1)
      p1 = next((p for p in reversed(range(self.nelem))   if indices[p0] < indices[p  ]),  0)
      indices [p0], indices [p1] = indices [p1], indices [p0]
      elements[p0], elements[p1] = elements[p1], elements[p0]
      sgn *= -1
      indices [p0+1:] = indices [p0+1:][::-1]
      elements[p0+1:] = elements[p0+1:][::-1]
      sgn *= (-1) ** ( (self.nelem-p0-1)/2 )


