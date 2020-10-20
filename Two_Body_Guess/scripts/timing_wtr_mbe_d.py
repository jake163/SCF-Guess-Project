import os
import time
import pickle
from mbe_density import *

wtr_mbe_d_times = {}
filename = "wtr_mbe_d_times.py"
bs = "aug-cc-pvdz"

root_dir  = os.getcwd()
geom_dir  = os.path.join(root_dir, 'geometries')
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
    
if "Trimer" not in wtr_mbe_d_times:
    file_path = os.path.join(geom_dir, "Trimer.xyz")
    mol = load_molecule(file_path)
    mono2atom = [n[0:3], n[3:6], n[6:9]]
    monomers = make_nmers(mono2atom, mol, 1)
    dimers = make_nmers(mono2atom, mol, 2)
    wtr_mbe_d_times["Trimer-time_monomers"] = time_monomer_densities(monomers, bs)
    wtr_mbe_d_times["Trimer-time_dimers"] = time_dimer_densities(dimers, bs)
    
if "Tetramer" not in wtr_mbe_d_times:
    file_path = os.path.join(geom_dir, "Tetramer.xyz")
    mol = load_molecule(file_path)
    mono2atom = [n[0:3], n[3:6], n[6:9], n[9:12]]
    monomers = make_nmers(mono2atom, mol, 1)
    dimers = make_nmers(mono2atom, mol, 2)
    wtr_mbe_d_times["Tetramer-time_monomers"] = time_monomer_densities(monomers, bs)
    wtr_mbe_d_times["Tetramer-time_dimers"] = time_dimer_densities(dimers, bs)
    
if "Pentamer" not in wtr_mbe_d_times:
    file_path = os.path.join(geom_dir, "Pentamer.xyz")
    mol = load_molecule(file_path)
    mono2atom = [n[0:3], n[3:6], n[6:9], n[9:12], n[12:15]]
    monomers = make_nmers(mono2atom, mol, 1)
    dimers = make_nmers(mono2atom, mol, 2)
    wtr_mbe_d_times["Pentamer-time_monomers"] = time_monomer_densities(monomers, bs)
    wtr_mbe_d_times["Pentamer-time_dimers"] = time_dimer_densities(dimers, bs)
    
if "Hexamer" not in wtr_mbe_d_times:
    file_path = os.path.join(geom_dir, "Hexamer_Prism.xyz")
    mol = load_molecule(file_path)
    mono2atom = [n[0:3], n[3:6], n[6:9], n[9:12], n[12:15], n[15:18]]
    monomers = make_nmers(mono2atom, mol, 1)
    dimers = make_nmers(mono2atom, mol, 2)
    wtr_mbe_d_times["Hexamer_Prism-time_monomers"] = time_monomer_densities(monomers, bs)
    wtr_mbe_d_times["Hexamer_Prism-time_dimers"] = time_dimer_densities(dimers, bs)
    
if "Heptamer" not in wtr_mbe_d_times:
    file_path = os.path.join(geom_dir, "Heptamer.xyz")
    mol = load_molecule(file_path)
    mono2atom = [n[0:3], n[3:6], n[6:9], n[9:12], n[12:15], n[15:18], n[18:21]]
    monomers = make_nmers(mono2atom, mol, 1)
    dimers = make_nmers(mono2atom, mol, 2)
    wtr_mbe_d_times["Heptamer-time_monomers"] = time_monomer_densities(monomers, bs)
    wtr_mbe_d_times["Heptamer-time_dimers"] = time_dimer_densities(dimers, bs)
    
if "Octamer" not in wtr_mbe_d_times:
    file_path = os.path.join(geom_dir, "Octamer.xyz")
    mol = load_molecule(file_path)
    mono2atom = [n[0:3], n[3:6], n[6:9], n[9:12], n[12:15], n[15:18], n[18:21], n[21:24]]
    monomers = make_nmers(mono2atom, mol, 1)
    dimers = make_nmers(mono2atom, mol, 2)
    wtr_mbe_d_times["Octamer-time_monomers"] = time_monomer_densities(monomers, bs)
    wtr_mbe_d_times["Octamer-time_dimers"] = time_dimer_densities(dimers, bs)  
                
if "Nonamer" not in wtr_mbe_d_times:
    file_path = os.path.join(geom_dir, "Nonamer.xyz")
    mol = load_molecule(file_path)
    mono2atom = [n[0:3], n[3:6], n[6:9], n[9:12], n[12:15], n[15:18], n[18:21], n[21:24], n[24:27]]
    monomers = make_nmers(mono2atom, mol, 1)
    dimers = make_nmers(mono2atom, mol, 2)
    wtr_mbe_d_times["Nonamer-time_monomers"] = time_monomer_densities(monomers, bs)
    wtr_mbe_d_times["Nonamer-time_dimers"] = time_dimer_densities(dimers, bs)
                
if "Decamer" not in wtr_mbe_d_times:
    file_path = os.path.join(geom_dir, "Decamer.xyz")
    mol = load_molecule(file_path)
    mono2atom = [n[0:3], n[3:6], n[6:9], n[9:12], n[12:15], n[15:18], n[18:21], n[21:24], n[24:27], n[27:30]]
    monomers = make_nmers(mono2atom, mol, 1)
    dimers = make_nmers(mono2atom, mol, 2)
    wtr_mbe_d_times["Decamer-time_monomers"] = time_monomer_densities(monomers, bs)
    wtr_mbe_d_times["Decamer-time_dimers"] = time_dimer_densities(dimers, bs)
                
with open(filename, 'wb') as f:
    pickle.dump(wtr_mbe_d_times, f)