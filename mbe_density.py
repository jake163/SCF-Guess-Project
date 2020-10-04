import os
import itertools
import pyscf
import numpy as np
from wrap_pyscf import run_scf

def load_molecule(file_path):
    mol = []
    with open(file_path, 'r') as f:
        for line in f:
            split_line = line.split()
            if len(split_line) == 4:
                mol += [' '.join(split_line)]
    return mol

def make_nmers(mono2atom, mol, n):
    nmers      = {}
    n_monomers = len(mono2atom)
    for nmer in itertools.combinations(range(n_monomers), n):
        nmers[tuple(nmer)] = []
        for mono in nmer:
            for atom in mono2atom[mono]:
                nmers[(nmer)] += [mol[atom]]
    return nmers

def compute_densities(nmers, bs):
    data = {}
    for tag, nmer in nmers.items():
        mol       = pyscf.gto.Mole()
        mol.atom  = nmer
        mol.basis = bs
        mol.build()
        n_cycles = 0
        def get_cycles(locl):
            nonlocal n_cycles
            if 'cycle' in locl: n_cycles = locl['cycle'] + 1
        mf = pyscf.scf.HF(mol)
        mf.callback = get_cycles
        mf.kernel()
        data[tuple(tag)] = (n_cycles, mf.make_rdm1())
    return data

def atom2ao(mol, bs):
    ss = pyscf.gto.Mole()
    ss.atom = mol
    ss.basis = bs
    ss.build()
    return ss.offset_ao_by_atom()

def make_guess(mono_densities, dimer_densities, mol, mono2atom,  bs):
    a2ao = atom2ao(mol, bs)
    nbf = a2ao[-1][-1]
    rho = np.zeros([nbf, nbf])

    factor = -(len(mono_densities) - 2) if dimer_densities else 1
    for tag, data in mono_densities.items():
        atoms = []
        for mono in tag:
            atoms += mono2atom[mono]

        guess_map = []
        for atom in atoms:
            atom_bs_start = a2ao[atom][2]
            atom_bs_end   = a2ao[atom][3]
            guess_map += range(atom_bs_start, atom_bs_end)

        for tag_i, rho_i in enumerate(guess_map):
            for tag_j, rho_j in enumerate(guess_map):
                rho[rho_i, rho_j] += factor * data[1][tag_i, tag_j]

    for tag, data in dimer_densities.items():
        atoms = []
        for mono in tag:
            atoms += mono2atom[mono]

        guess_map = []
        for atom in atoms:
            atom_bs_start = a2ao[atom][2]
            atom_bs_end   = a2ao[atom][3]
            guess_map += range(atom_bs_start, atom_bs_end)

        for tag_i, rho_i in enumerate(guess_map):
            for tag_j, rho_j in enumerate(guess_map):
                rho[rho_i, rho_j] += data[1][tag_i, tag_j]

    return rho


if __name__ == '__main__':
    root_dir  = os.getcwd()
    geom_dir  = os.path.join(root_dir, 'geometries')
    file_path = os.path.join(geom_dir, 'Trimer.xyz')
    mol = load_molecule(file_path)
    mono2atom = [[0, 1, 2], [3, 4, 5], [6, 7, 8]]
    bs = 'sto-3g'

    monomers = make_nmers(mono2atom, mol, 1)
    dimers   = make_nmers(mono2atom, mol, 2)

    rhos = {}
    rhos[1] = compute_densities(monomers, bs)
    rhos[2] = compute_densities(dimers, bs)

    guess = make_guess(rhos[1], rhos[2], mol, mono2atom, bs)

    run_scf(mol, bs, guess)
    run_scf(mol, bs)
