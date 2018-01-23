import numpy as np

def tensor_basis_transform(int_tensor, C_bra, C_ket):
    """Transform a tensor of bra-ket products from one basis to another.
    tensordot works by transforming the first row and moving that to the back.
    e.g., mnrs -> nrsP

    Normally, C_bra = C_ket.conjugate(), but not necessarily."""
    n_dim = len(int_tensor.shape)
    if n_dim % 2:
        raise Exception("Cannot transform basis of an odd-dimensional tensor.")
    for dim in range(n_dim // 2):
        int_tensor = np.tensordot(int_tensor, C_bra, axes=(0,0))
    for dim in range(n_dim // 2):
        int_tensor = np.tensordot(int_tensor, C_ket, axes=(0,0))
    return int_tensor