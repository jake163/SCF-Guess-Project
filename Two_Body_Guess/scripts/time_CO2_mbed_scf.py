import pickle
from time_mbed import *

CO2_mbed_scf_times = {}
filename = "CO2_mbed_scf_times.py"
class_dir  = os.path.join(geom_dir, 'CO2_clusters')

if "Tetramer" not in CO2_mbed_scf_times:
    file_path = os.path.join(class_dir, "4.xyz")
    mol = load_molecule(file_path)
    mono2atom = [n[0:3], n[3:6], n[6:9], n[9:12]]
    CO2_mbed_scf_times["Tetramer-atom_g_scf_time"] = atomic_scf_time(mol, bs)
    CO2_mbed_scf_times["Tetramer-mbed_scf_time"] = mbed_scf_time(mono2atom, mol, bs, n)

if "Pentamer" not in CO2_mbed_scf_times:
    file_path = os.path.join(class_dir, "5.xyz")
    mol = load_molecule(file_path)
    mono2atom = [n[0:3], n[3:6], n[6:9], n[9:12], n[12:15]]
    CO2_mbed_scf_times["Pentamer-atom_g_scf_time"] = atomic_scf_time(mol, bs)
    CO2_mbed_scf_times["Pentamer-mbed_scf_time"] = mbed_scf_time(mono2atom, mol, bs, n)
    
if "Hexamer" not in CO2_mbed_scf_times:
    file_path = os.path.join(class_dir, "6.xyz")
    mol = load_molecule(file_path)
    mono2atom = [n[0:3], n[3:6], n[6:9], n[9:12], n[12:15], n[15:18]]
    CO2_mbed_scf_times["Hexamer-atom_g_scf_time"] = atomic_scf_time(mol, bs)
    CO2_mbed_scf_times["Hexamer-mbed_scf_time"] = mbed_scf_time(mono2atom, mol, bs, n)
                 
if "Heptamer" not in CO2_mbed_scf_times:
    file_path = os.path.join(class_dir, "7.xyz")
    mol = load_molecule(file_path)
    mono2atom = [n[0:3], n[3:6], n[6:9], n[9:12], n[12:15], n[15:18], n[18:21]]
    CO2_mbed_scf_times["Heptamer-atom_g_scf_time"] = atomic_scf_time(mol, bs)
    CO2_mbed_scf_times["Heptamer-mbed_scf_time"] = mbed_scf_time(mono2atom, mol, bs, n)
                 
if "Octamer" not in CO2_mbed_scf_times:
    file_path = os.path.join(class_dir, "8.xyz")
    mol = load_molecule(file_path)
    mono2atom = [n[0:3], n[3:6], n[6:9], n[9:12], n[12:15], n[15:18], n[18:21], n[21:24]]
    CO2_mbed_scf_times["Octamer-atom_g_scf_time"] = atomic_scf_time(mol, bs)
    CO2_mbed_scf_times["Octamer-mbed_scf_time"] = mbed_scf_time(mono2atom, mol, bs, n)
                
if "Nonamer" not in CO2_mbed_scf_times:
    file_path = os.path.join(class_dir, "9.xyz")
    mol = load_molecule(file_path)
    mono2atom = [n[0:3], n[3:6], n[6:9], n[9:12], n[12:15], n[15:18], n[18:21], n[21:24], n[24:27]]
    CO2_mbed_scf_times["Nonamer-atom_g_scf_time"] = atomic_scf_time(mol, bs)
    CO2_mbed_scf_times["Nonamer-mbed_scf_time"] = mbed_scf_time(mono2atom, mol, bs, n)
                 
if "Decamer" not in CO2_mbed_scf_times:
    file_path = os.path.join(class_dir, "10.xyz")
    mol = load_molecule(file_path)
    mono2atom = [n[0:3], n[3:6], n[6:9], n[9:12], n[12:15], n[15:18], n[18:21], n[21:24], n[24:27], n[27:30]]
    CO2_mbed_scf_times["Decamer-atom_g_scf_time"] = atomic_scf_time(mol, bs)
    CO2_mbed_scf_times["Decamer-mbed_scf_time"] = mbed_scf_time(mono2atom, mol, bs, n) 
                 
with open(filename, 'wb') as f:
    pickle.dump(CO2_mbed_scf_times, f)