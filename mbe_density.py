import os
import itertools
import pyscf

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
    monomers = [i for i in range(n_monomers)]
    for nmer in itertools.combinations(monomers, n):
        nmers[tuple(nmer)] = []
        for mono in nmer:
            for atom in mono2atom[mono]:
                nmers[(mono,)] += [mol[atom]]
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
        mf.callback= get_cycles
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
    nbf  = a2ao[-1][-1]
    rho = numpy.zeros([nbf, nbf])
    for tag, data in dimer_densities:
        for mono_i in tag:
            for atom_i in mono2atom[mono_i]:
                ss_start_i = a2ao[atom_i][2]
                ss_end_i   = a2ao[atom_i][3]
                nbf_i      = ss_end_i - ss_start_i
                for offset_i in range(nbf_i):
                    for mono_j in tag:
                        for atom_j in mono2atom[mono_j]:
                            ss_start_j = a2ao[atom_j][2]
                            ss_end_j   = a2ao[atom_j][3]
                            nbf_j      = ss_end_j - ss_start_j 
                            for offset_j in range(nbf_j):
                                rho[ss_start_i + offset_i, ss_start_j + offset_j] += data[1][offset_i, offset_j]


if __name__ == '__main__':
    root_dir  = os.getcwd()
    geom_dir  = os.path.join(root_dir, 'geometries')
    file_path = os.path.join(geom_dir, 'Trimer.xyz')
    mol = load_molecule(file_path)
    mono2atom = [[0, 1, 2], [3, 4, 5], [6, 7, 8]]
    bs = 'sto-3g'

    monomers  = make_nmers(mono2atom, mol, 1)
    dimers    = make_nmers(mono2atom, mol, 2)
    trimer    = make_nmers(mono2atom, mol, 3)
    
    rhos = {}
    rhos[1]  = compute_densities(monomers, bs)
    rhos[2]  = compute_densities(dimers, bs)
    rhos[3]  = compute_densities(trimer, bs)
    
    guess = make_guess(rhos[1], rhos[2], mol, bs)
    run_scf(mol, bs, guess)

