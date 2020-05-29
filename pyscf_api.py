import pyscf

def density(molecule, basis, verbose = 0):
    mol = pyscf.M(atom = molecule, basis = basis)
    mf = pyscf.scf.HF(mol)
    mf.verbose = verbose
    mf.kernel()
    return mf.make_rdm1()

def run_scf(molecule, basis, density = None, verbose = 4):   
    n_cycles = 0
    def get_cycles(locl):
        nonlocal n_cycles
        if "cycle" in locl: n_cycles = locl["cycle"] + 1

    mol = pyscf.M(atom = molecule, basis = basis)
    mf = pyscf.scf.HF(mol)
    mf.verbose = verbose
    mf.callback = get_cycles
    mf.kernel(density)
    print(n_cycles)


if __name__ == "__main__":
    mol = [
        ("H", 0.0, 0.0, 0.89),
        ("H", 0.0, 0.0, -0.89),
	    ("O", 0.0, 0.0, 0.0),
	    ("H", 1.0, 0.0, 0.89),
        ("H", 1.0, 0.0, -0.89),
	    ("O", 1.0, 0.0, 0.0)
    ]

    run_scf(mol, "sto-3g")