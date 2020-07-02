from wrap_pyscf import *

methyl_1 = [0,2,3,4]
methyl_2 = [1,5,6,7]

funct_grps = [methyl_1, methyl_2]

guess_driver("Ethane.sdf", funct_grps, [1,1])
