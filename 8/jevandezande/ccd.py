import numpy as np
from cc import CC


class CCD(CC):
    """
    Coupled Cluster Doubles
    """
    
    def __init__(self, uhf):
        super().__init__(uhf)
        self.method_name = 'CCD'

    def update_amplitudes(self, t):
        """
        Update the CC amplitudes (differs for each theory level)
        """
        g = self.gmo
        o  = slice(None, self.nocc)
        v  = slice(self.nocc, None)
        w0 = g[o, o, v, v]
        w1 = 1/2 * np.einsum('abcd,ijcd->ijab', g[v, v, v, v], t)
        w2 = 1/2 * np.einsum('klij,klab->ijab', g[o, o, o, o], t)
        w3 = np.einsum('akic,jkbc->ijab', g[v, o, o, v], t)
        w3 += -w3.transpose((0, 1, 3, 2)) - w3.transpose(1, 0, 2, 3) - w3.transpose((1, 0, 3, 2))
        w4 = np.einsum('klcd,ijac,klbd->ijab', g[o, o, v, v], t, t)
        w4 += -w4.transpose((0, 1, 3, 2))
        w5 = 1/2 * np.einsum('klcd,ikab,jlcd->ijab', g[o, o, v, v], t, t)
        w5 += w5.transpose((1, 0, 2, 3))
        w6 = 1/4 * np.einsum('klcd,ijcd,klab->ijab', g[o, o, v, v], t, t)
        w7 = 1/4 * np.einsum('klcd,ikac,jlbd->ijab', g[o, o, v, v], t, t)
        w7 += w7.transpose((1, 0, 2, 3))

        return w0 + w1 + w2 + w3 + w4 + w5 + w6 + w7


if __name__ == "__main__":
    import sys
    sys.path.insert(0, '../../5/jevandezande')
    from uhf import UHF
    uhf = UHF('../../3/jevandezande/Options.ini')
    uhf.energy()
    ccd = CCD(uhf)
    ccd.energy()
    ccd.plot_convergence()
