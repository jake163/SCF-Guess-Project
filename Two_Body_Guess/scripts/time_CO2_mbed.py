import pickle
from time_mbed import *

CO2_mbed_times = {}
filename = "CO2_mbed_times.py"
class_dir  = os.path.join(geom_dir, 'CO2_clusters')

if "Tetramer" not in CO2_mbed_times:
    file_path = os.path.join(class_dir, "4.xyz")
    mol = load_molecule(file_path)
    mono2atom = [n[0:3], n[3:6], n[6:9], n[9:12]]
    monomers = make_nmers(mono2atom, mol, 1)
    dimers = make_nmers(mono2atom, mol, 2)
    CO2_mbed_times["Tetramer-time_monomers"] = time_monomer_densities(monomers, bs)
    CO2_mbed_times["Tetramer-time_dimers"] = time_dimer_densities(dimers, bs)
    
if "Pentamer" not in CO2_mbed_times:
    file_path = os.path.join(class_dir, "5.xyz")
    mol = load_molecule(file_path)
    mono2atom = [n[0:3], n[3:6], n[6:9], n[9:12], n[12:15]]
    monomers = make_nmers(mono2atom, mol, 1)
    dimers = make_nmers(mono2atom, mol, 2)
    CO2_mbed_times["Pentamer-time_monomers"] = time_monomer_densities(monomers, bs)
    CO2_mbed_times["Pentamer-time_dimers"] = time_dimer_densities(dimers, bs)
    
if "Hexamer" not in CO2_mbed_times:
    file_path = os.path.join(class_dir, "6.xyz")
    mol = load_molecule(file_path)
    mono2atom = [n[0:3], n[3:6], n[6:9], n[9:12], n[12:15], n[15:18]]
    monomers = make_nmers(mono2atom, mol, 1)
    dimers = make_nmers(mono2atom, mol, 2)
    CO2_mbed_times["Hexamer-time_monomers"] = time_monomer_densities(monomers, bs)
    CO2_mbed_times["Hexamer-time_dimers"] = time_dimer_densities(dimers, bs)
    
if "Heptamer" not in CO2_mbed_times:
    file_path = os.path.join(class_dir, "7.xyz")
    mol = load_molecule(file_path)
    mono2atom = [n[0:3], n[3:6], n[6:9], n[9:12], n[12:15], n[15:18], n[18:21]]
    monomers = make_nmers(mono2atom, mol, 1)
    dimers = make_nmers(mono2atom, mol, 2)
    CO2_mbed_times["Heptamer-time_monomers"] = time_monomer_densities(monomers, bs)
    CO2_mbed_times["Heptamer-time_dimers"] = time_dimer_densities(dimers, bs)
    
if "Octamer" not in CO2_mbed_times:
    file_path = os.path.join(class_dir, "8.xyz")
    mol = load_molecule(file_path)
    mono2atom = [n[0:3], n[3:6], n[6:9], n[9:12], n[12:15], n[15:18], n[18:21], n[21:24]]
    monomers = make_nmers(mono2atom, mol, 1)
    dimers = make_nmers(mono2atom, mol, 2)
    CO2_mbed_times["Octamer-time_monomers"] = time_monomer_densities(monomers, bs)
    CO2_mbed_times["Octamer-time_dimers"] = time_dimer_densities(dimers, bs)
    
if "Nonamer" not in CO2_mbed_times:
    file_path = os.path.join(class_dir, "9.xyz")
    mol = load_molecule(file_path)
    mono2atom = [n[0:3], n[3:6], n[6:9], n[9:12], n[12:15], n[15:18], n[18:21], n[21:24], n[24:27]]
    monomers = make_nmers(mono2atom, mol, 1)
    dimers = make_nmers(mono2atom, mol, 2)
    CO2_mbed_times["Nonamer-time_monomers"] = time_monomer_densities(monomers, bs)
    CO2_mbed_times["Nonamer-time_dimers"] = time_dimer_densities(dimers, bs)
    
if "Decamer" not in CO2_mbed_times:
    file_path = os.path.join(class_dir, "10.xyz")
    mol = load_molecule(file_path)
    mono2atom = [n[0:3], n[3:6], n[6:9], n[9:12], n[12:15], n[15:18], n[18:21], n[21:24], n[24:27], n[27:30]]
    monomers = make_nmers(mono2atom, mol, 1)
    dimers = make_nmers(mono2atom, mol, 2)
    CO2_mbed_times["Decamer-time_monomers"] = time_monomer_densities(monomers, bs)
    CO2_mbed_times["Decamer-time_dimers"] = time_dimer_densities(dimers, bs)
    
with open(filename, 'wb') as f:
    pickle.dump(CO2_mbed_times, f)