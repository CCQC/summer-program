class P(object):

  def __init__(self, ax1, ax2):
    self.ax1 = ax1
    self.ax2 = ax2

  def __mul__(self, array):
    return array - array.swapaxes(self.ax1,self.ax2)






















