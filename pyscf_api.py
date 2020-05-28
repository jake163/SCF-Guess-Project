from pyscf import gto
from pyscf import scf

def run_scf(atoms, coords, basis, density = None):
    assert(len(atoms) == len(coords))

    atomic_coords = [[atom, coord] for atom, coord in zip(atoms, coords)]
    mol = gto.Mole()
    mol.atom = atomic_coords
    mol.basis = basis
    mol.build()

    mf = scf.HF(mol)
    mf.kernel()
    
    pass