import pyscf

atom2z = {'h' : 1, 'he' : 2,
          'li' : 3, ' be' : 4, 'n' : 5, 'c': 6, 'n':7, ' o' : 8, 'f': 9, 'ne':10,
          }

def make_molecule(atoms, basis):
    nelec = 0
    for a in atoms:
        at_sym = a[0].lower()
        if at_sym not in atom2z:
            raise "Atomic symbol {} is not in atom2z. Please add it".format(at_sym)
        nelec += atom2z[at_sym]
    mol = pyscf.gto.Mole()
    mol.atom = atoms
    mol.basis = basis
    mol.spin = 0 if nelec % 2 == 0 else 1
    mol.build()
    return mol

def density(molecule, basis, verbose = 0):
    mol = make_molecule(molecule, basis)
    mf = pyscf.scf.HF(mol)
    mf.verbose = verbose
    mf.kernel()
    return mf.make_rdm1()

def run_scf(molecule, basis, density = None, verbose = 4):   
    n_cycles = 0
    def get_cycles(locl):
        nonlocal n_cycles
        if "cycle" in locl: n_cycles = locl["cycle"] + 1

    mol = make_molecule(molecule, basis)
    mf = pyscf.scf.HF(mol)
    mf.verbose = verbose
    mf.callback = get_cycles
    mf.kernel(density)
    return n_cycles


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

