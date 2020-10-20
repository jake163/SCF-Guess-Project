import time
import pickle
from mbe_density import *

wtr_mbe_d_scf = {}
filename = "wtr_mbe_d_scf.py"
bs = "aug-cc-pvdz"

root_dir  = os.getcwd()
geom_dir  = os.path.join(root_dir, 'geometries')
n = list(range(30))

if "Trimer" not in wtr_mbe_d_scf:
    file_path = os.path.join(geom_dir, "Trimer.xyz")
    mol = load_molecule(file_path)
    mono2atom = [n[0:3], n[3:6], n[6:9]]
    monomers = make_nmers(mono2atom, mol, 1)
    dimers = make_nmers(mono2atom, mol, 2)
    guess = make_guess(compute_densities(monomers, bs), compute_densities(dimers, bs), mol, mono2atom, bs)
    wtr_mbe_d_scf["Trimer-mbe_d_g"] = run_scf(mol, bs, guess)
    wtr_mbe_d_scf["Trimer-atom_g"] = run_scf(mol, bs)
    
if "Tetramer" not in wtr_mbe_d_scf:
    file_path = os.path.join(geom_dir, "Tetramer.xyz")
    mol = load_molecule(file_path)
    mono2atom = [n[0:3], n[3:6], n[6:9], n[9:12]]
    monomers = make_nmers(mono2atom, mol, 1)
    dimers = make_nmers(mono2atom, mol, 2)
    guess = make_guess(compute_densities(monomers, bs), compute_densities(dimers, bs), mol, mono2atom, bs)
    wtr_mbe_d_scf["Tetramer-mbe_d_g"] = run_scf(mol, bs, guess)
    wtr_mbe_d_scf["Tetramer-atom_g"] = run_scf(mol, bs)
    
if "Pentamer" not in wtr_mbe_d_scf:
    file_path = os.path.join(geom_dir, "Pentamer.xyz")
    mol = load_molecule(file_path)
    mono2atom = [n[0:3], n[3:6], n[6:9], n[9:12], n[12:15]]
    monomers = make_nmers(mono2atom, mol, 1)
    dimers = make_nmers(mono2atom, mol, 2)
    guess = make_guess(compute_densities(monomers, bs), compute_densities(dimers, bs), mol, mono2atom, bs)
    wtr_mbe_d_scf["Pentamer-mbe_d_g"] = run_scf(mol, bs, guess)
    wtr_mbe_d_scf["Pentamer-atom_g"] = run_scf(mol, bs)
    
#if "Hexamer_Prism" not in wtr_mbe_d_scf:
 #   file_path = os.path.join(geom_dir, "Hexamer_Prism.xyz")
  #  mol = load_molecule(file_path)
   # mono2atom = [n[0:3], n[3:6], n[6:9], n[9:12], n[12:15], n[15:18]]
   # monomers = make_nmers(mono2atom, mol, 1)
   # dimers = make_nmers(mono2atom, mol, 2)
   # guess = make_guess(compute_densities(monomers, bs), compute_densities(dimers, bs), mol, mono2atom, bs)
   # wtr_mbe_d_scf["Hexamer_Prism-mbe_d_g"] = run_scf(mol, bs, guess)
   # wtr_mbe_d_scf["Hexamer_Prism-atom_g"] = run_scf(mol, bs)
    
if "Heptamer" not in wtr_mbe_d_scf:
    file_path = os.path.join(geom_dir, "Heptamer.xyz")
    mol = load_molecule(file_path)
    mono2atom = [n[0:3], n[3:6], n[6:9], n[9:12], n[12:15], n[15:18], n[18:21]]
    monomers = make_nmers(mono2atom, mol, 1)
    dimers = make_nmers(mono2atom, mol, 2)
    guess = make_guess(compute_densities(monomers, bs), compute_densities(dimers, bs), mol, mono2atom, bs)
    wtr_mbe_d_scf["Heptamer-mbe_d_g"] = run_scf(mol, bs, guess)
    wtr_mbe_d_scf["Heptamer-atom_g"] = run_scf(mol, bs)
    
if "Octamer" not in wtr_mbe_d_scf:
    file_path = os.path.join(geom_dir, "Octamer.xyz")
    mol = load_molecule(file_path)
    mono2atom = [n[0:3], n[3:6], n[6:9], n[9:12], n[12:15], n[15:18], n[18:21], n[21:24]]
    monomers = make_nmers(mono2atom, mol, 1)
    dimers = make_nmers(mono2atom, mol, 2)
    guess = make_guess(compute_densities(monomers, bs), compute_densities(dimers, bs), mol, mono2atom, bs)
    wtr_mbe_d_scf["Octamer-mbe_d_g"] = run_scf(mol, bs, guess)
    wtr_mbe_d_scf["Octamer-atom_g"] = run_scf(mol, bs)  
                
if "Nonamer" not in wtr_mbe_d_scf:
    file_path = os.path.join(geom_dir, "Nonamer.xyz")
    mol = load_molecule(file_path)
    mono2atom = [n[0:3], n[3:6], n[6:9], n[9:12], n[12:15], n[15:18], n[18:21], n[21:24], n[24:27]]
    monomers = make_nmers(mono2atom, mol, 1)
    dimers = make_nmers(mono2atom, mol, 2)
    guess = make_guess(compute_densities(monomers, bs), compute_densities(dimers, bs), mol, mono2atom, bs)
    wtr_mbe_d_scf["Nonamer-mbe_d_g"] = run_scf(mol, bs, guess)
    wtr_mbe_d_scf["Nonamer-atom_g"] = run_scf(mol, bs)
                
if "Decamer" not in wtr_mbe_d_scf:
    file_path = os.path.join(geom_dir, "Decamer.xyz")
    mol = load_molecule(file_path)
    mono2atom = [n[0:3], n[3:6], n[6:9], n[9:12], n[12:15], n[15:18], n[18:21], n[21:24], n[24:27], n[27:30]]
    monomers = make_nmers(mono2atom, mol, 1)
    dimers = make_nmers(mono2atom, mol, 2)
    guess = make_guess(compute_densities(monomers, bs), compute_densities(dimers, bs), mol, mono2atom, bs)
    wtr_mbe_d_scf["Decamer-mbe_d_g"] = run_scf(mol, bs, guess)
    wtr_mbe_d_scf["Decamer-atom_g"] = run_scf(mol, bs)
                
with open(filename, 'wb') as f:
    pickle.dump(wtr_mbe_d_scf, f)