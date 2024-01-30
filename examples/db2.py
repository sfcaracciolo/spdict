from src.spdict  import * 
from geometric_plotter import Plotter 
import numpy as np 
import pathlib
filename = pathlib.Path(__file__).stem
p = Plotter(_2d=True, figsize=(8,8))

atoms = get_swt_atoms('db2', 4) # get atoms
D = create_sparse_dict_from_atoms(
    atoms,
    n = 512,
    steps = 40,
    normalize = False
)

offsets = np.arange(D.shape[1])[np.newaxis, :]
p.axs.plot(D.toarray(order='F') + offsets, '-k', linewidth=.5)
p.axs.set_xlabel('n')
p.axs.set_ylabel('number of atom')
Plotter.show()
p.save(folder='figs/', name=filename, format = 'png' )