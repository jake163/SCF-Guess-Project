import time
from mbe_density import *

bs = "aug-cc-pvdz"
root_dir = os.getcwd()
geom_dir = os.path.join(root_dir, 'geometries')
n = list(range(30))

def time_monomer_densities(monomers, bs):
    start = time.time()
    monomer_dens = compute_densities(monomers, bs)
    end = time.time()
    elapsed = end-start
    return elapsed

def time_dimer_densities(dimers, bs):
    start = time.time()
    dimer_dens = compute_densities(dimers, bs)
    end = time.time()
    elapsed = end-start
    return elapsed

def atomic_scf_time(mol, bs):
    start = time.time()
    nitr = run_scf(mol, bs)
    end = time.time()
    elapsed = end-start
    return nitr, elapsed

def mbed_scf_time(mono2atom, mol, bs, n):
    start = time.time()
    monomers = make_nmers(mono2atom, mol, n=1)
    dimers = make_nmers(mono2atom, mol, n=2)
    guess = make_guess(compute_densities(monomers, bs), compute_densities(dimers, bs), mol, mono2atom, bs)
    nitr = run_scf(mol, bs, guess)
    end = time.time()
    elapsed = end-start
    return nitr, elapsed