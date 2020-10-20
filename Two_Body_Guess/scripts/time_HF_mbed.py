import pickle
from time_mbed import *

HF_mbed_times = {}
filename = "HF_mbed_times.py"
class_dir  = os.path.join(geom_dir, 'HF_clusters')

if "Dimer" not in HF_mbed_times:
    file_path = os.path.join(class_dir, "2.xyz")
    mol = load_molecule(file_path)
    mono2atom = [n[0:2], n[2:4]]
    monomers = make_nmers(mono2atom, mol, 1)
    dimers = make_nmers(mono2atom, mol, 2)
    HF_mbed_times["Dimer-time_monomers"] = time_monomer_densities(monomers, bs)
    HF_mbed_times["Dimer-time_dimers"] = time_dimer_densities(dimers, bs)
    
if "Tetramer" not in HF_mbed_times:
    file_path = os.path.join(class_dir, "4.xyz")
    mol = load_molecule(file_path)
    mono2atom = [n[0:2], n[2:4], n[4:6], n[6:8]]
    monomers = make_nmers(mono2atom, mol, 1)
    dimers = make_nmers(mono2atom, mol, 2)
    HF_mbed_times["Tetramer-time_monomers"] = time_monomer_densities(monomers, bs)
    HF_mbed_times["Tetramer-time_dimers"] = time_dimer_densities(dimers, bs)
    
if "Pentamer" not in HF_mbed_times:
    file_path = os.path.join(class_dir, "5.xyz")
    mol = load_molecule(file_path)
    mono2atom = [n[0:2], n[2:4], n[4:6], n[6:8], n[8:10]]
    monomers = make_nmers(mono2atom, mol, 1)
    dimers = make_nmers(mono2atom, mol, 2)
    HF_mbed_times["Pentamer-time_monomers"] = time_monomer_densities(monomers, bs)
    HF_mbed_times["Pentamer-time_dimers"] = time_dimer_densities(dimers, bs)

if "Hexamer" not in HF_mbed_times:
    file_path = os.path.join(class_dir, "6.xyz")
    mol = load_molecule(file_path)
    mono2atom = [n[0:2], n[2:4], n[4:6], n[6:8], n[8:10], n[10:12]]
    monomers = make_nmers(mono2atom, mol, 1)
    dimers = make_nmers(mono2atom, mol, 2)
    HF_mbed_times["Hexamer-time_monomers"] = time_monomer_densities(monomers, bs)
    HF_mbed_times["Hexamer-time_dimers"] = time_dimer_densities(dimers, bs)
    
if "Heptamer" not in HF_mbed_times:
    file_path = os.path.join(class_dir, "7.xyz")
    mol = load_molecule(file_path)
    mono2atom = [n[0:2], n[2:4], n[4:6], n[6:8], n[8:10], n[10:12], n[12:14]]
    monomers = make_nmers(mono2atom, mol, 1)
    dimers = make_nmers(mono2atom, mol, 2)
    HF_mbed_times["Heptamer-time_monomers"] = time_monomer_densities(monomers, bs)
    HF_mbed_times["Heptamer-time_dimers"] = time_dimer_densities(dimers, bs)
    
if "Octamer" not in HF_mbed_times:
    file_path = os.path.join(class_dir, "8.xyz")
    mol = load_molecule(file_path)
    mono2atom = [n[0:2], n[2:4], n[4:6], n[6:8], n[8:10], n[10:12], n[12:14], n[14:16]]
    monomers = make_nmers(mono2atom, mol, 1)
    dimers = make_nmers(mono2atom, mol, 2)
    HF_mbed_times["Octamer-time_monomers"] = time_monomer_densities(monomers, bs)
    HF_mbed_times["Octamer-time_dimers"] = time_dimer_densities(dimers, bs)
    
other_dir = os.path.join(geom_dir, 'H2O_clusters')
    
if "Wtr_Hexamer_Prism" not in HF_mbed_times:
    file_path = os.path.join(other_dir, "Hexamer_Prism.xyz")
    mol = load_molecule(file_path)
    mono2atom = [n[0:3], n[3:6], n[6:9], n[9:12], n[12:15], n[15:18]]
    monomers = make_nmers(mono2atom, mol, 1)
    dimers = make_nmers(mono2atom, mol, 2)
    HF_mbed_times["Wtr_Hexamer_Prism-time_monomers"] = time_monomer_densities(monomers, bs)
    HF_mbed_times["Wtr_Hexamer_Prism-time_dimers"] = time_dimer_densities(dimers, bs)
    
with open(filename, 'wb') as f:
    pickle.dump(HF_mbed_times, f)