import pickle

from wrap_pyscf import *

scf_guess = {}

filename = "scf_record"

#with open(filename, 'rb') as outfile:
    #scf_guess = pickle.load(outfile)

if "EthaneA 1Cx2" not in scf_guess:
    funct_grps = [[0,2,3,4], [1,5,6,7]]
    scf_guess["EthaneA 1Cx2"] = guess_driver("EthaneA.sdf", funct_grps, [1,1])
        
if "PropaneA 1Cx3" not in scf_guess:
    funct_grps = [[2,0,1,3], [4,5,6], [7,8,9,10]]
    scf_guess["PropaneA 1Cx3"] = guess_driver("PropaneA.sdf", funct_grps, [1,2,1])
        
if "PropaneA cc-c" not in scf_guess:
    funct_grps = [[2,0,1,3,4,5,6], [7,8,9,10]]
    scf_guess["PropaneA cc-c"] = guess_driver("PropaneA.sdf", funct_grps, [1,1])
        
if "ButaneA 1Cx4" not in scf_guess:
    funct_grps = [[2,6,5,7], [0,3,4], [1,9,10], [8,11,12,13]]
    scf_guess["ButaneA 1Cx4"] = guess_driver("ButaneA.sdf", funct_grps, [1,2,2,1])
        
if "ButaneA 2Cx2" not in scf_guess:
    funct_grps = [[2,6,5,7,0,3,4], [1,9,10,8,11,12,13]]
    scf_guess["ButaneA 2Cx2"] = guess_driver("ButaneA.sdf", funct_grps, [1,1])
        
if "ButaneA ccc-c" not in scf_guess:
    funct_grps = [[2,6,5,7,0,3,4,1,9,10], [8,11,12,13]]
    scf_guess["ButaneA ccc-c"] = guess_driver("ButaneA.sdf", funct_grps, [1,1])
    
if "HexaneA 1Cx6" not in scf_guess:
    funct_grps = [[0,12,15,4], [2,18,19], [14,17,13], [11,3,16], [1,6,7], [5,8,9,10]]
    scf_guess["HexaneA 1Cx6"] = guess_driver("HexaneA.sdf", funct_grps, [1,2,2,2,2,1])

if "HexaneA 2Cx3" not in scf_guess:
    funct_grps = [[0,12,15,4,2,18,19], [14,17,13,11,3,16], [1,6,7,5,8,9,10]]
    scf_guess["Hexane 2Cx3"] = guess_driver("HexaneA.sdf", funct_grps, [1,2,1])

if "HexaneA 3Cx2" not in scf_guess:
    funct_grps = [[0,12,15,4,2,18,19,14,17,13], [11,3,16,1,6,7,5,8,9,10]]
    scf_guess["Hexane 3Cx3"] = guess_driver("HexaneA.sdf", funct_grps, [1,1])

with open(filename, 'wb') as infile:
    pickle.dump(scf_guess, infile)
