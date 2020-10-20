import os
import time
import pickle
from mbe_density import *

wtr_mbe_d_scf_times = {}
filename = "wtr_mbe_d_scf_times.py"
bs = "aug-cc-pvdz"

root_dir  = os.getcwd()
geom_dir  = os.path.join(root_dir, 'geometries')

n = list(range(30))

def atomic_scf_time(mol, bs):
    start = time.time()
    nitr = run_scf(mol, bs)
    end = time.time()
    elapsed = end-start
    return elapsed

def mbe_d_scf_time(mono2atom, mol, bs, n):
    start = time.time()
    monomers = make_nmers(mono2atom, mol, n=1)
    dimers = make_nmers(mono2atom, mol, n=2)
    guess = make_guess(compute_densities(monomers, bs), compute_densities(dimers, bs), mol, mono2atom, bs)
    nitr = run_scf(mol, bs, guess)
    end = time.time()
    elapsed = end-start
    return elapsed
    
if "Trimer" not in wtr_mbe_d_scf_times:
    file_path = os.path.join(geom_dir, "Trimer.xyz")
    mol = load_molecule(file_path)
    mono2atom = [n[0:3], n[3:6], n[6:9]]
    wtr_mbe_d_scf_times["Trimer-atom_g_scf_time"] = atomic_scf_time(mol, bs)
    wtr_mbe_d_scf_times["Trimer-mbe_d_scf_time"] = mbe_d_scf_time(mono2atom, mol, bs, n)
    
if "Tetramer" not in wtr_mbe_d_scf_times:
    file_path = os.path.join(geom_dir, "Tetramer.xyz")
    mol = load_molecule(file_path)
    mono2atom = [n[0:3], n[3:6], n[6:9], n[9:12]]
    wtr_mbe_d_scf_times["Tetramer-atom_g_scf_time"] = atomic_scf_time(mol, bs)
    wtr_mbe_d_scf_times["Tetramer-mbe_d_scf_time"] = mbe_d_scf_time(mono2atom, mol, bs, n)
    
if "Pentamer" not in wtr_mbe_d_scf_times:
    file_path = os.path.join(geom_dir, "Pentamer.xyz")
    mol = load_molecule(file_path)
    mono2atom = [n[0:3], n[3:6], n[6:9], n[9:12], n[12:15]]
    wtr_mbe_d_scf_times["Pentamer-atom_g_scf_time"] = atomic_scf_time(mol, bs)
    wtr_mbe_d_scf_times["Pentamer-mbe_d_scf_time"] = mbe_d_scf_time(mono2atom, mol, bs, n)
    
#if "Hexamer" not in wtr_mbe_d_scf_times:
 #   file_path = os.path.join(geom_dir, "Hexamer_Prism.xyz")
  #  mol = load_molecule(file_path)
   # mono2atom = [n[0:3], n[3:6], n[6:9], n[9:12], n[12:15], n[15:18]]
   # wtr_mbe_d_scf_times["Hexamer_Prism-atom_g_scf_time"] = atomic_scf_time(mol, bs)
   # wtr_mbe_d_scf_times["Hexamer_Prism-mbe_d_scf_time"] = mbe_d_scf_time(mono2atom, mol, bs, n)
    
if "Heptamer" not in wtr_mbe_d_scf_times:
    file_path = os.path.join(geom_dir, "Heptamer.xyz")
    mol = load_molecule(file_path)
    mono2atom = [n[0:3], n[3:6], n[6:9], n[9:12], n[12:15], n[15:18], n[18:21]]
    wtr_mbe_d_scf_times["Heptamer-atom_g_scf_time"] = atomic_scf_time(mol, bs)
    wtr_mbe_d_scf_times["Heptamer-mbe_d_scf_time"] = mbe_d_scf_time(mono2atom, mol, bs, n)
    
if "Octamer" not in wtr_mbe_d_scf_times:
    file_path = os.path.join(geom_dir, "Octamer.xyz")
    mol = load_molecule(file_path)
    mono2atom = [n[0:3], n[3:6], n[6:9], n[9:12], n[12:15], n[15:18], n[18:21], n[21:24]]
    wtr_mbe_d_scf_times["Octamer-atom_g_scf_time"] = atomic_scf_time(mol, bs)
    wtr_mbe_d_scf_times["Octamer-mbe_d_scf_time"] = mbe_d_scf_time(mono2atom, mol, bs, n)  
                
if "Nonamer" not in wtr_mbe_d_scf_times:
    file_path = os.path.join(geom_dir, "Nonamer.xyz")
    mol = load_molecule(file_path)
    mono2atom = [n[0:3], n[3:6], n[6:9], n[9:12], n[12:15], n[15:18], n[18:21], n[21:24], n[24:27]]
    wtr_mbe_d_scf_times["Nonamer-atom_g_scf_time"] = atomic_scf_time(mol, bs)
    wtr_mbe_d_scf_times["Nonamer-mbe_d_scf_time"] = mbe_d_scf_time(mono2atom, mol, bs, n)
                
if "Decamer" not in wtr_mbe_d_scf_times:
    file_path = os.path.join(geom_dir, "Decamer.xyz")
    mol = load_molecule(file_path)
    mono2atom = [n[0:3], n[3:6], n[6:9], n[9:12], n[12:15], n[15:18], n[18:21], n[21:24], n[24:27], n[27:30]]
    wtr_mbe_d_scf_times["Decamer-atom_g_scf_time"] = atomic_scf_time(mol, bs)
    wtr_mbe_d_scf_times["Decamer-mbe_d_scf_time"] = mbe_d_scf_time(mono2atom, mol, bs, n)
                
with open(filename, 'wb') as f:
    pickle.dump(wtr_mbe_d_scf_times, f)