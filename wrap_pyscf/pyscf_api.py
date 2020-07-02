import pyscf
import os
from scipy.linalg import block_diag
import time

def read_geometry(file_name):
    cwd = os.getcwd()
    rv = []
    with open(os.path.join(cwd, "geometries", file_name)) as f:
        counter = 0
        natoms = 0
        for line in f:
            if counter == 3:
                natoms = int(line.split()[0])
            elif counter == 4 + natoms:
                break
            elif counter > 3:
                atom = line.split()
                at_sym = atom[3]
                x,y,z = float(atom[0]), float(atom[1]), float(atom[2])
                rv.append([at_sym, x, y, z])
            counter += 1
        return rv

def fragment_molecule(mol, atoms):
    return [[mol[atom] for atom in frag] for frag in atoms]

def get_guess(frags, unpaired_es, basis):
    if len(frags) != len(unpaired_es):
        raise Exception("Must provide the number of unpaired electrons for each"
                        " fragment")

    rhos = []
    for frag, ne in zip(frags, unpaired_es):
        rho = density(frag, basis, unpaired_electrons = ne)
        if len(rho.shape) == 3:
            rhos.append(rho[0] + rho[1])
        elif len(rho.shape) == 2:
            rhos.append(rho)
    return block_diag(*rhos)

def make_molecule(atoms, basis, unpaired_electrons = 0):
    mol = pyscf.gto.Mole()
    mol.atom = atoms
    mol.basis = basis
    mol.spin = unpaired_electrons
    mol.build()
    return mol

def density(molecule, basis, verbose = 0, unpaired_electrons = 0):
    mol = make_molecule(molecule, basis, unpaired_electrons)
    mf = pyscf.scf.HF(mol)
    mf.verbose = verbose
    mf.kernel()
    return mf.make_rdm1()

def run_scf(molecule, basis, density = None, unpaired_electrons = 0,
            verbose = 4):
    n_cycles = 0
    def get_cycles(locl):
        nonlocal n_cycles
        if "cycle" in locl: n_cycles = locl["cycle"] + 1

    mol = make_molecule(molecule, basis, unpaired_electrons)
    mf = pyscf.scf.HF(mol)
    mf.verbose = verbose
    mf.callback = get_cycles
    mf.kernel(density)
    return n_cycles

def guess_driver(file_name, frags, unpaired_es, basis = "aug-cc-pvdz"):
    mol = read_geometry(file_name)

    # Reorder molecule into same order as fragments
    new_mol = []
    for frag in frags:
        for atom in frag:
            new_mol.append(mol[atom])
    if len(new_mol) != len(mol):
        raise Exception("Was each atom in the molecule assigned to one (and "
                        "only one) fragment?")

    fragments = fragment_molecule(mol, frags)

    # Assemble the guess and time how long it takes
    start = time.time()
    rho_guess = get_guess(fragments, unpaired_es, basis)
    end = time.time()
    print("Time to assemble guess: {} s".format(end - start))

    # Run the full calculation for a point of comparison
    start = time.time()
    orig_n_cycles = run_scf(new_mol, basis)
    end = time.time()
    print("Time for full calculation: {} s".format(end - start))

    # Run it with our guess
    start = time.time()
    new_n_cycles  = run_scf(new_mol, basis, density = rho_guess)
    end = time.time()
    print("Time for better guess: {} s".format(end - start))

    return orig_n_cycles, new_n_cycles


if __name__ == "__main__":
    mol = [
        ("H", 0.0, 0.0, 0.89),
        ("H", 0.0, 0.0, -0.89),
        ("O", 0.0, 0.0, 0.0),
        ("H", 1.0, 0.0, 0.89),
        ("H", 1.0, 0.0, -0.89),
        ("O", 1.0, 0.0, 0.0)
    ]

    rho_total  = density(mol, "sto-3g", 0)
    n_itr_none = run_scf(mol, "sto-3g")
    n_itr_rho  = run_scf(mol, "sto-3g", rho_total)
    print(n_itr_none)
    print(n_itr_rho)

