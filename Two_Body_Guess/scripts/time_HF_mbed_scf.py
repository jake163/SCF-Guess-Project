import pickle
from time_mbed import *

HF_mbed_scf_times = {}
filename = "HF_mbed_scf_times.py"
class_dir  = os.path.join(geom_dir, 'HF_clusters')

if "Dimer" not in HF_mbed_scf_times:
    file_path = os.path.join(class_dir, "2.xyz")
    mol = load_molecule(file_path)
    mono2atom = [n[0:2], n[2:4]]
    HF_mbed_scf_times["Dimer-atom_g_scf_time"] = atomic_scf_time(mol, bs)
    HF_mbed_scf_times["Dimer-mbed_scf_time"] = mbed_scf_time(mono2atom, mol, bs, n)

if "Tetramer" not in HF_mbed_scf_times:
    file_path = os.path.join(class_dir, "4.xyz")
    mol = load_molecule(file_path)
    mono2atom = [n[0:2], n[2:4], n[4:6], n[6:8]]
    HF_mbed_scf_times["Tetramer-atom_g_scf_time"] = atomic_scf_time(mol, bs)
    HF_mbed_scf_times["Tetramer-mbed_scf_time"] = mbed_scf_time(mono2atom, mol, bs, n)
    
if "Pentamer" not in HF_mbed_scf_times:
    file_path = os.path.join(class_dir, "5.xyz")
    mol = load_molecule(file_path)
    mono2atom = [n[0:2], n[2:4], n[4:6], n[6:8], n[8:10]]
    HF_mbed_scf_times["Pentamer-atom_g_scf_time"] = atomic_scf_time(mol, bs)
    HF_mbed_scf_times["Pentamer-mbed_scf_time"] = mbed_scf_time(mono2atom, mol, bs, n)
                 
if "Hexamer" not in HF_mbed_scf_times:
    file_path = os.path.join(class_dir, "6.xyz")
    mol = load_molecule(file_path)
    mono2atom = [n[0:2], n[2:4], n[4:6], n[6:8], n[8:10], n[10:12]]
    HF_mbed_scf_times["Hexamer-atom_g_scf_time"] = atomic_scf_time(mol, bs)
    HF_mbed_scf_times["Hexamer-mbed_scf_time"] = mbed_scf_time(mono2atom, mol, bs, n)
                 
#if "Heptamer" not in HF_mbed_scf_times:
 #   file_path = os.path.join(class_dir, "7.xyz")
  #  mol = load_molecule(file_path)
   # mono2atom = [n[0:2], n[2:4], n[4:6], n[6:8], n[8:10], n[10:12], n[12:14]]
   # HF_mbed_scf_times["Heptamer-atom_g_scf_time"] = atomic_scf_time(mol, bs)
   # HF_mbed_scf_times["Heptamer-mbed_scf_time"] = mbed_scf_time(mono2atom, mol, bs, n)
                
if "Octamer" not in HF_mbed_scf_times:
    file_path = os.path.join(class_dir, "8.xyz")
    mol = load_molecule(file_path)
    mono2atom = [n[0:2], n[2:4], n[4:6], n[6:8], n[8:10], n[10:12], n[12:14], n[14:16]]
    HF_mbed_scf_times["Octamer-atom_g_scf_time"] = atomic_scf_time(mol, bs)
    HF_mbed_scf_times["Octamer-mbed_scf_time"] = mbed_scf_time(mono2atom, mol, bs, n)
                 
#other_dir = os.path.join(geom_dir, 'H2O_clusters')                 
                 
#if "Wtr_Hexamer_Prism" not in HF_mbed_scf_times:
 #   file_path = os.path.join(other_dir, "Hexamer_Prism_2.xyz")
  #  mol = load_molecule(file_path)
   # mono2atom = [n[0:3], n[3:6], n[6:9], n[9:12], n[12:15], n[15:18]]
    #HF_mbed_scf_times["Wtr_Hexamer_Prism-atom_g_scf_time"] = atomic_scf_time(mol, bs)
    #HF_mbed_scf_times["Wtr_Hexamer_Prism-mbed_scf_time"] = mbed_scf_time(mono2atom, mol, bs, n) 
                 
with open(filename, 'wb') as f:
    pickle.dump(HF_mbed_scf_times, f)