# SCF Guess Project

## Contents

- [Prerequisites](#Prerequisites)
- [Installing](#Installing)
- [Usage](#Usage)
- [API](#API)
- [Acknowledgements](#Acknowledgements)

## Prerequisites

To use this repository you will need to have the following installed:

- git
- Python 3

Both of these dependencies are very common dependencies and are likely available
in your operating system's package repository (*e.g.*, Aptitude on Ubuntu)

## Installing

To obtain and install this repo in a Linux-like environment run the following
commands in a terminal:

```bash
# Obtain the source code and enter the directory
git clone https://github.com/rmrresearch/SCF-Guess-Project.git scf_guess
cd scf_guess

# Create a virtual Python environment to keep our dependencies isolated from the
# system (also ensures reproducibility)
python3 -m venv venv
source venv/bin/activate

# Install this repository's dependencies
pip install -r requirements.txt
```

The above assumes that `git` and `python3` are installed in a location included
in your path. If they are not you will need to specify the full path to the
`git` and `python3` executables.

## Usage

Once installed this repo is designed to be usable as a Python module by another
Python script or Jupyter notebook. For example, put the following inside a
script `h2_aug_cc_pvdz.py`

```python3
from wrap_pyscf import *

# Define the molecule in Cartesian coordinates
mol = [("H", 0.0, 0.0, 0.0),
       ("H", 0.0, 0.0, 0.89)]

# The name of the atomic basis set
bs  = "aug-cc-pvdz"

# Compute the density of the molecule at the SCF/aug-cc-pVDZ level of theory
rho = density(mol, bs)

# Test that the density is converged (return the number of iterations needed to)
# converge the SCF using rho as a guess)
nitr = run_scf(mol, bs, rho)

print(nitr) # This should print out 1 (1 cycle is needed to test convergence)
```

This script is then run as:
```bash
python3 h2_aug_cc_pvdz.py
```
and the output should be something along the lines of (all but the last line are
output from PySCF's SCF run):

```bash
******** <class 'pyscf.scf.hf.RHF'> ********
method = RHF
initial guess = minao
damping factor = 0
level_shift factor = 0
DIIS = <class 'pyscf.scf.diis.CDIIS'>
diis_start_cycle = 1
diis_space = 8
SCF conv_tol = 1e-09
SCF conv_tol_grad = None
SCF max_cycles = 50
direct_scf = True
direct_scf_tol = 1e-13
chkfile to save SCF result = /home/ryan/CLionProjects/wrap_pyscf/tmpmc9v4foe
max_memory 4000 MB (current use 115 MB)
Set gradient conv threshold to 3.16228e-05
init E= -1.11798956774344
  HOMO = -0.551284046862827  LUMO = 0.0586573242914685
cycle= 1 E= -1.1179895677435  delta_E= -6.31e-14  |g|= 7.54e-08  |ddm|= 3.68e-07
  HOMO = -0.551284016599658  LUMO = 0.0586573291817753
Extra cycle  E= -1.1179895677435  delta_E= -4.44e-15  |g|= 2.03e-08  |ddm|= 1.02e-07
converged SCF energy = -1.1179895677435
1
```

## API

The `wrap_pyscf` module defines two functions `density` and `run_scf` whose APIs
are detailed here.

### density

`density` is used to compute the atomic density of the provided molecular system
using the provided atomic basis set.

```python3
the_density density(molecule, basis_set, verbose=0)
```

where:

- `the_density` is a `numpy.ndarray` instance containing the converged SCF
  density matrix in the AO basis set.
- `molecule` is a `[(str, float, float, float)]` (*i.e.* a list of four element
   tuples) such that the `i`-th tuple is the atomic symbol and `x`, `y`, and `z`
   Cartesian coordinates (in Angstroms) for the `i`-th atom
- `basis_set` is the name of the atomic basis set to use
- `verbose` is an optional integer controlling how much printing is output by
  PySCF. Default is 0, *i.e.*, no printing.

### run_scf

`run_scf` is used to compute the SCF energy of the provided molecular system
using the provided atomic basis set and optionally the provided atomic density.

```python3
nitr run_scf(molecule, basis_set, the_density = None, verbose = 4)
```

where:

- `nitr` is the number of iterations that the SCF took to converge
- `molecule` is a `[(str, float, float, float)]` (*i.e.* a list of four element
   tuples) such that the `i`-th tuple is the atomic symbol and `x`, `y`, and `z`
   Cartesian coordinates (in Angstroms) for the `i`-th atom
- `basis_set` is the name of the atomic basis set to use
- `the_density` is a `numpy.ndarray` instance containing the density to use as a
   guess for the SCF density matrix. Passing a value of `None`, the default,
   indicates that the program should use the default guess.
- `verbose` is an optional integer controlling how much printing is output by
  PySCF. Default is 4, which will print out basic SCF information.

## Acknowledgements

- Jonathan Waldrop - Created the simplified API for calling PySCF
