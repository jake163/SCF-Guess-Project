from wrap_pyscf import *

import pickle

scf_guess = {}

filename = "Dimer_scf-vs-radius-record.pickle"

funct_grps = [[1,0,2], [3,4,5]]

if "Dimer equill" not in scf_guess:
    scf_guess["Dimer equill"] = guess_driver("Dimer.xyz", funct_grps, [0,0])
    
if "Dimer ~4ang" not in scf_guess:
    scf_guess["Dimer ~4ang"] = guess_driver("Dimer_1.xyz", funct_grps, [0,0])
    
if "Dimer ~6ang" not in scf_guess:
    scf_guess["Dimer ~6ang"] = guess_driver("Dimer_2.xyz", funct_grps, [0,0])
    
if "Dimer ~8ang" not in scf_guess:
    scf_guess["Dimer ~8ang"] = guess_driver("Dimer_3.xyz", funct_grps, [0,0])
    
if "Dimer ~10ang" not in scf_guess:
    scf_guess["Dimer ~10ang"] = guess_driver("Dimer_4.xyz", funct_grps, [0,0])
    
if "Dimer ~12ang" not in scf_guess:
    scf_guess["Dimer ~12ang"] = guess_driver("Dimer_5.xyz", funct_grps, [0,0])
    
if "Dimer ~14ang" not in scf_guess:
    scf_guess["Dimer ~14ang"] = guess_driver("Dimer_6.xyz", funct_grps, [0,0])
    
if "Dimer ~15ang" not in scf_guess:
    scf_guess["Dimer ~15ang"] = guess_driver("Dimer_7.xyz", funct_grps, [0,0])
        
with open(filename, 'wb') as infile:
    pickle.dump(scf_guess, infile)
