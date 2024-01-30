import pywt
import numpy as np
import scipy as sp

def create_sparse_dict_from_atoms(atoms: list, n: int, steps: int = 1, **kwargs) -> sp.sparse.csc_matrix:
    """ 
        atoms: list of np.ndarray
        n: amount of samples
        step: subsampling
    """
    sizes = np.array([a.size for a in atoms], dtype=np.uint16)
    if np.any(n < sizes):
        raise ValueError(f'n value should be greater than the larger atom {sizes}.')

    if isinstance(steps, int):
        steps = len(atoms)*[steps]

    if len(steps) != len(atoms):
        raise ValueError(f'The amount of steps and atoms should be equal.')


    blocks = []
    for atom, step in zip(atoms, steps):
        cols = n + 1 - atom.size
        blocks.append(sparse_convolution_matrix(atom, cols, step=step, **kwargs))

    return sp.sparse.hstack(blocks, format='csc')

def sparse_convolution_matrix(a: np.ndarray, n: int, step: int = 1, normalize: bool = True) -> sp.sparse.csc_matrix:
    """
    Implementation of sparse convolution matrix in full mode only [1].
    
    a: 1D np.ndarray represent an atom
    n: amount of columns. Is the size of v in A@v

    [1] https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.convolution_matrix.html#scipy.linalg.convolution_matrix
    """
    m = a.size
    k = m + n - 1

    if normalize:
        a /= np.linalg.norm(a)

    data = np.broadcast_to(a[:,np.newaxis], shape=(m, n))
    offsets = np.arange(start=0, stop=-m, step=-1)

    D = sp.sparse.dia_matrix((data, offsets), shape=(k, n))

    return D.tocsc(copy=False)[:,::step]

def get_swt_atoms(wname: str = 'db4', max_level: int = 1) -> list:
    atoms = []
    wobj = pywt.Wavelet(wname)
    atoms.append(get_discrete_scale(wobj, max_level))
    for level in range(max_level, 0, -1):
        atoms.append(get_discrete_wavelet(wobj, level))
    return atoms

def get_wp_atoms(wname: str = 'db4', max_level: int = 1) -> list:
    atoms = []
    wobj = pywt.Wavelet(wname)
    for level in range(max_level, 0, -1):
        atoms.append(get_discrete_scale(wobj, level))
        atoms.append(get_discrete_wavelet(wobj, level))
    return atoms

def get_cwt_atoms(wname: str = 'morl', max_level: int = 1) -> list:
    atoms = []
    wobj = pywt.ContinuousWavelet(wname)
    for level in range(max_level, 0, -1):
        atoms.append(get_continuous_wavelet(wobj, level))
    return atoms

def get_discrete_scale(w: pywt.Wavelet, level: int):
    if w.orthogonal:
        scale, _, _ = w.wavefun(level)
    else:
        _, _, scale, _, _ = w.wavefun(level)
    return scale

def get_discrete_wavelet(w: pywt.Wavelet, level: int):
    if w.orthogonal:
        _, wavelet, _ = w.wavefun(level)
    else:
        _, _, _, wavelet, _ = w.wavefun(level)
    return wavelet

def get_continuous_wavelet(w: pywt.ContinuousWavelet, level: int):
    wavelet, _ = w.wavefun(level)
    return wavelet

