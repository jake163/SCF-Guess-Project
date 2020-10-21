from wrap_pyscf import *

import pickle

scf_guess = {}

filename = "water_cluster_record.pickle"

if "Dimer 1wtrx2" not in scf_guess:
    funct_grps = [[0,1,2], [3,4,5]]
    scf_guess["Dimer 1wtrx2"] = guess_driver("Dimer.xyz", funct_grps, [0,0])

if "Trimer 1wtrx3" not in scf_guess:
    funct_grps = [[0,1,2], [3,4,5], [6,7,8]]
    scf_guess["Trimer 1wtrx3"] = guess_driver("Trimer.xyz", funct_grps, [0,0,0])

if "Tetramer 1wtrx4" not in scf_guess:
    funct_grps = [[0,1,2], [3,4,5], [6,7,8], [9,10,11]]
    scf_guess["Tetramer 1wtrx4"] = guess_driver("Tetramer.xyz", funct_grps, [0,0,0,0])

if "Tetramer 2wtrx2" not in scf_guess:
    funct_grps = [[0,1,2,3,4,5], [6,7,8,9,10,11]]
    scf_guess["Tetramer 2wtrx2"] = guess_driver("Tetramer.xyz", funct_grps, [0,0])

if "Pentamer 1wtrx5" not in scf_guess:
    funct_grps = [[4,9,14], [3,13,8], [2,12,7], [1,11,6], [0,5,10]]
    scf_guess["Pentamer 1wtrx5"] = guess_driver("Pentamer.xyz", funct_grps, [0,0,0,0,0])

if "Hexamer_Prism 1wtrx6" not in scf_guess:
    funct_grps = [[0,1,2], [3,4,5], [6,7,8], [9,10,11], [12,13,14], [15,16,17]]
    scf_guess["Hexamer_Prism 1wtrx6"] = guess_driver("Hexamer_Prism.xyz", funct_grps, [0,0,0,0,0,0])

if "Hexamer_Prism 2wtrx3" not in scf_guess:
    funct_grps = [[0,1,2,3,4,5], [6,7,8,9,10,11], [12,13,14,15,16,17]]
    scf_guess["Hexamer_Prism 2wtrx3"] = guess_driver("Hexamer_Prism.xyz", funct_grps, [0,0,0])

if "Hexamer_Prism 3wtrx2" not in scf_guess:
    funct_grps = [[0,1,2,3,4,5,6,7,8], [9,10,11,12,13,14,15,16,17]]
    scf_guess["Hexamer_Prism 3wtrx2"] = guess_driver("Hexamer_Prism.xyz", funct_grps, [0,0])

with open(filename, 'wb') as infile:
    pickle.dump(scf_guess, infile)
    
print(scf_guess)
