"An Extended Hückel Theory. I. Hydrocarbons"
R. Hoffmann, J. Chem. Phys. 39, 1397 (1963); https://doi.org/10.1063/1.1734456

This seems to be the seminal Extended Huckel Theory paper. A brief description of the guess is that it is Huckel theory with an
extended basis set comprised of 2s and 2p Carbon, and 1s Hydrogen orbitals, with overlap and interaction considerations. LCAO's
are calculated for MO's. The EHT Hamiltonian is simply a summation of one-electron effective Hamiltonians, and the matrix
elements are systematically determined. For Hii, they use reported values of valence state ionization potentials of Carbon sp3
valence states C2p, C2s, and H1s. Hij is approximated with eqn 3: Hij = 0.5K(Hii+Hjj)Sij, where K = 1.75 in this paper. Sij is
internally computed, and the energy minimization Huckel equation (eqn 2) is solved with two matrix diagonalizations. There is no
emphasis on SCF efficiency.

Energy barriers to numerous hydrocarbons, e.g. internal rotation, ring conformations, and geometric isomerization, are
calculated. They highlight their successful charge distribution calculations for aliphatics, and sigma and pi system orbital
energy calculations in aromatics. They acknowledge failures in the overemphasis of steric factors, which lead to inacurrate
isomerization energies, and the failure to predict strain energies.  



"Approximate Self‐Consistent Molecular Orbital Theory. III. CNDO Results for AB2 and AB3 Systems"
J.A. Pople, J. Chem. Phys. 44, 3289 (1966); https://doi.org/10.1063/1.1727227

Here the authors examine an approximate method self-consistent MO's for all valence electrons in molecules that contain no atoms
heavier than fluorine. The key aspect of this method is that the overlap distribution of any two atomic orbitals is neglected in
all electron repulsion integrals: a "complete neglect of differential overlap" (CNDO). Some parameters are chosen
semiempirically to best fit the data. This paper proposes a modification to CNDO they call "CNDO/2." Their method fixes two
failures of CNDO (and maintains other similarities); bond lengths were too short and binding energies were too large because of
interatomic electron shell penetration and subsequent attractions unaccounted for, and the way in which the local core matrix
element was determined with the atomic ionization potentials. Their respective fixes involve adding "penetration terms" to the
expression for the Hartree-Fock matrix F, and using the average of atomic electron affinities and atomic ionization potentials.
They also extended their method to open shell configurations. There is no emphasis on SCF efficiency.

They tested symmetrical triatomic molecules (AB2) of D∞h and bent C2v point groups, as well as certain AB3 molecules. They
subsequently calculated equillibrium bond anges, force constants, dipole moments, and atomic charges. When B is Oxygen or
Fluorine, AB2 equillibrium bond angles are in fairly good agreement with experimental values. Dipole moments and bending force
constants were in moderate agreement with experimental values, with exception to BeF2 in which the force constant was too low.
AB3 bond angles were also in good agreement with experiments done, and certain geometries were correctly predicted.



"Molecular orbital set determined by a localization procedure" 
H. Weinstein, R. Pauncz, Symp. Faraday Soc., 1968, 2, 23-31; https://doi.org/10.1039/SF9680200023

Localized molecular orbitals (LMOs) are examined for a group of hydride molecules LiH, BH, and BH3. Since LMOs describe
electron densities localized to well defined regions, they are more in tune with known chemcial properties. These orbitals are
described as a set of molecular maximum overlap orbitals that are localized using a local density maximation as a criterion for
transformation. The external localization procedure they use are similar to previously reported "localization functions",
however their method does not require two electron integral calculations to obtain good starting wavefunctions for SCF
calculations. They define orbitals as any electron pair, so their set of orbitals contain molecular functions that describe
bonds between atoms, molecular orbitals for inner shells, and molecular orbitals for lone pairs. The sum of these orbitals are
half the total number of electrons a given molecule. The starting point of their localization calculation begins with a set of
eigenvectors from the diagonalization of the atomic overlap matrix (they call them maximum overlap orbitals). They then
separated  functions with large localization numerical values into separate sets, followed by orthogonal transformations between
pairs of orbitals in the separate sets. 

They used an STO atomic basis set for the diatomic molecules, consisting of 1s, 2s, and 2p(x,y,z) on the heavy atom, and 1s for
the hydrogen, with exponents optimized for energy. For the polyatomic, the same type of basis set is used except that teh
exponents used were reported. Energies and overlap of orbitals obtained from their localization procedure are more accurate
compared to their starting set (maximum overlap orbitals). They observed a small improvement for BH because of its small atomic
basis. But for their other molecules, the large localization sets and smaller localization sets were larger, and they therefore
were able to show improved MO sets suitable for beginning SCF calculations.



"A 'Level–Shifting' method for converging closed shell Hartree–Fock wave functions"
V.R. Saunders, I.H. Hillier, Int. J. Quantum Chem. 1973, 7, 4, 699-705; https://doi.org/10.1002/qua.560070407

Once again my Iowa State affiliation doesn't allow me to climb the paywall here. 



"Convergence acceleration of iterative sequences. the case of scf iteration"
P. Pulay, Chem. Phys. Lett. 1980, 73, 2, 393-398; https://doi.org/10.1016/0009-2614(80)80396-4

In SCF iteration Pulay uses a quasi-Newton-Raphson procedure for energy minimzation. A trial set of parameters of the energy
functional is iteratively improved by acting on the gradient vector by an approximate inverse Hessian matrix. This method is
called simple relaxation (SR). Unless the approximate Hessian is a good approximate, convergence is slow. To fix this method,
they define the Direct Inversion of Iterative Subspace Method (DIIS). The procedure is this: choose a starting vector trial
parameter (these vectors are chosen such that residuum vector approaches the zero vector in a mean-square way) and approximate
inverse Hessian, iterate by SR until the quadratic region is reached -> The parameter vector is stored after each SR iteration,
and solve a linear system of equations presented paper's method section -> if convergence below predetermined limit has not
occurred, build a new parameter vector and check with SR for convergence. In this paper, elements of the Fock matrix 
transformed to an orthonormal basis are worked with in the algorithm, i.e. S(^-1/2)FS(^-1/2). 

Results are presented as mean-square deviations of the SR and DIIS polynomials against the zero functions, with data points
for n number of terms. As n goes from 2 to 10, SR has deviations that are 1.1 to 3000 times larger than DIIS on a (-0.5, +0.5)
interval. The (-0.3, +0.7) interval for SR showed a 1.1 to 7900 larger size, n from 2 to 10. They also obtained results for
CNDO/2 and ab initio SCF calculations for butadiene and CN+ ion. CNDO/2 calculations for butadiene showed basic SCF convergence
is fair, and that DIIS gains were moderate. Although slowly converging, DIIS showed dramatic improvement for CNDO/2 SCF
convergence on CN+. PF3 ab initio calculations showed equally promising results, with the 7-term DIIS procedure yielding a 3
orders of magnitude improvement in the convergence parameter, which is about 10 SCF cycles.



"Improved SCF convergence acceleration" 
P. Pulay, J. Comput. Chem. 1982, 3, 4, 556-560; https://doi.org/10.1002/jcc.540030413

This paper proposes a method that improves upon the above method, DIIS. They say that this improved method can be applied to
what they consider to be divergent or slowly converging cases. (Old) DIIS is briefly redescribed: an error vector (related to
the gradient of electronic energy from w.r.t. SCF parameters) is constructed in each SCF calculation, where the vanishing of the
vector below an established threshold constitutes convergence. Three errors are declared for their original algorithm: A Fock
matrix (for no other purpose) had to be made for the determination of the error vector, their error vector (calculated as a
difference between two consecutive Fock matrices) was not accurate enough, and DIIS extrapolation was performed periodically
after an arbitrary number of cycles (5-12) thereby restricting it's efficacy for slowly converging or divergent systems. Their
new algorithm, what I call new DIIS, addresses these issues. They derive different equations for the calculation of the error
matrices (eqn 4 is the result), evaluate scalar products found with the integrated square of the error vectors, and lastly find
and diagonalize the Fock matrix for an arbitrary SCF cycle where a new density matrix and new Fock matrix can be determined.

Their new method has shown to be superior to old DIIS. Generally the results show that the maximum element of the error vector
decreases an order of magnitude per iteration step, where ordinary SCF shows a decrease of a factor of two per step.
Subsequently, SCF cycles are reduced to half (8 cycles) for acetylene (compared to ordinary SCF), and reduced to more than a
third (8 cycles) for water. They also compare their method to the Quadratically Convergent Hartree-Fock method (QC-SCF). For the
CO molecule and the H2C2N+ ion, via lowering of the error parameter, new DIIS outperforms QC-SCF by more than halving the number
of cycles for CO, and reducing SCF cycles for H2C2N+ by 58%. 

 
 
"Convergence properties of Hartree–Fock SCF molecular calculations" 
M.A. Natiello, G.E. Scuseria, Int. J. Quantum Chem. 1984, 26, 6, 1039-1049; https://doi.org/10.1002/qua.560260608

Paywall.



"Extension of the PS‐GVB electronic structure code to transition metal complexes"
D.T. Mainz, J.J. Klicic, R.A. Friesner, J.-M. Langlois, J.K. Perry, J. Comput. Chem. 1997, 18, 15, 1863-1874; https://doi.or
/10.1002/(SICI)1096-987X(19971130)18:15<1863::AID-JCC3>3.0.CO;2-M

Via the PS-GVB program, a parameterization enabling ab initio electronic structure calculations were done for systems containing
transition metal complexes using two standard effective core potential basis sets. Among improved computational efficiency for
single point energies and geometry optimization, their initial guess strategy for HF and DFT SCF in PS-GVB showed more reliable
convergence to ground state. Their Pseudospectral methods (PS) with Generalized Valence Bond techniques (GVB) approach involves
the standard use of contracted Gaussian basis functions to represent MO's, but they use a numerical grid to construct Coulomb
and exchange matrix elements. The key to the efficiency of their method, they claim, is to define local least-squares operators
for each atom, with systematic optimizing for each of the local operators. For accuracy, a small fraction of two-electron
integrals are done analytically. For any given system, numerical grids (radial) and fitting basis functions (normal Gaussian
basis functions with some FT dealiasing functions) must be defined, which they did for 25 different transition metal complexes.
With respect to their SCF guess, their approach attempts to find the best set of linear combination of valence orbitals with
maximum amount of bonding character, because the d orbitals in transition metals make semiempirical methods impractical. They do
HF-SCF calculations for each atomic basis independently. They then diagonalize the matrix of overlaps between atomic valence
orbitals. The eigenvectors with the largest eigenvalues have the greatest bonding character.

All comparisons were made with the Gaussian-92 ab initio electronic structure package, both programs were run in direct mode.
Their SCF guess is used for both programs. At the HF level of theory with LACVP** and LAV3P** basis sets, PS-GVB tends to save a
significant amount of time and iterations. In virtually all molecules for both basis', time is reduced by a factor of one to
three, with some complexes having user time reduced by a factor of 4, 5, 7, and 11. The total iterations
are reduced by a few iterations to half of iterations needed for convergence by Gaussian-92. The gain in time is most evident in
their larger systems (300+ basis functions), they observed a two to three factor speedup. For the total energies, converged
total energies mostly show agreement to better than 0.25 Kcal/mol. They did similar calculations for BLYP and B3LYP DFT
calculations, and the trends generally hold.



"A black-box self-consistent field convergence algorithm: One step closer"
K.N. Kudin, G.E. Scuseria, J. Chem. Phys. 116, 8255 (2002); https://doi.org/10.1063/1.1470195

Here a successor to the DIIS method is proposed. A vert brief summary of the DIIS method: An acceleration technique based on the
commutator of the Fock and density matrices, where SCF convergence is accomplished with a reasonable initial guess. These are
readily accomplished in most molecular systems, but convergence difficulties arise with systems that contain unusual bonding
patterns (transition states) and for transition metal/lantinide/actinide species with many eigenstates clustered around the
Fermi level. The relaxed constrained algorithm (RCA), which has been mathematically proven to converge to a solution to HF
equations regardless of initial guess, is extended to include techniques from the DIIS method. This method minimizes a function
that is based on energy, hence its name energy-DIIS or EDIIS. This paper implements the RCA equations (Section 2, eqns 1 and 2)
using the optimal damping algorithm (ODA). The ODA procedure is described in Section 2A, eqns 1-4. They briefly describe the
algorithm as a optimal step steepest descent algorithm, where a Fock matrix that is a function of a density satisfying relaxed 
constraints (pseudodensity) is calculated and a density matrix (specific pseudodensity matrix) is assembled using the aufbau
principle, a Fock matrix is assembled and E is computed, E is minimized by solving a line search problem, and finally obtain a
new Fock matrix with pseudodensity input. This makes EDIIS an interpolative procedure, where DIIS is an extrapolative procedure.
They note that occupation numbers can very freely (given the sum of occupation numbers in the set is equal to N electron pairs)
from range [0,1] ; at convergence they are equal to 0 or 1. 

They examine a few systems, where the performance of DIIS and EDIIS were examined. Their data consists of plots of log(En-Ec)
vs. Number of SCF cycles. For acetaldehyde with a RHF/6-31G(d) basis set using an initial guess obtained by diagonalizing the
core Hamiltonian, DIIS converged after 15 cycles versus EDIIS converging at 24. The second system examined was UF4 with a
RB3LYP/LANL2DZ basis set. Using a "Projected New-EHT Guess" (reference 14), DIIS failed to converge after 800 cycles, where
EDIIS needed 600. They note that the final energy is quickly reached, and then much cycles are spending converging to the
minimum. Lastly, CrC is examined with a RBLYP/6-31G(d) basis. They contend that there is likely no solution to KS equations
without violating aufbau principle, DIIS fails to converge (they did use level-shifting with DIIS to get integer occupations
with aufbau violations). EDIIS also failes to converge, arriving at a solution with fractional occupations (occupations jump
from cycle to cycle).



"The Field-Adapted ADMA Approach:  Introducing Point Charges"
T.E. Exner, P.G. Mezey, J. Phys. Chem. A 2004, 108, 19, 4301–4309; https://doi.org/10.1021/jp037447p

The method presented in this paper is an augmentation of the ADMA (Adjustable Density Matrix Assembler) guess density method,
the paper of which I included in the previous scf literature document. A recapitulation of ADMA: a macromolecule is divided into
fragments with fuzzy boundaries, in which SCF calculations can be performed on "parent" fragments (whose size is less than or
equal to the macromolecule size), that contain the smaller fragments ("families" of nuclei). The guess density matrix is a sum
of weighted parent densities (see eqns 1-3); a scalar factor of 1 is applied if atomic orbitals i and j are both centered on the
parent molecule's nuclei family, a scalar of 0.5 is applied if only 1 is, and 0 is applied if neither are centered on the parent
molecule's family fragment. This density is multiplied by the total number of electrons divided by the integrated electron
density to address error. Local interactions can be accounted for, and for a given size of parent molecule, computation time
scales linear with macromolecule size (the O(Nb^3) high power scaling bottleneck is avoided). With this approach (FA-ADMA), for
the parent molecules, they use point charges to simulate distant portions of the macromolecule; more accurate energies can be
obtained with smaller parent molecules. To accomplish this, they obtain partial charges of all atoms in macromolecule from an
ADMA calculation. These are then incorporated in new quantum chemical calculations of parent molecules, which gives a new ADMA
guess, which gives more accurate point charges, and so on until point charges converge below a threshold. That is then the
initial guess for HF SCF. 

In cases where hydrogens were placed to fixed dangling bonds on parent molecules, and effectively replace the bonded atom, they
looked at the best ways to include partial charges of the replaced atom (which were centered on the hydrogen in the
calculations). They also looked at parent molecule distances of 3, 4, and 5, angstroms. In all cases a partial scaling factor
(intermediate between the full partial charge and no partial charge of replaced atom) of 60-80% were found to give lowest
energies using an STO-3G basis for testing the hexapeptide and pentapeptide proteins, and 5 angstrom distances gave the lowest
energies. They used 16 protein structures from the protein data bank (PDB) with 6-31G** basis sets, and an 80% border atom
partial charge factor. It was found that FA-ADMA yields smaller and more accurate energies in proteins, especially proteins with
a relatively high degree of polarity and proteins that are formally charged. Compared to a direct Hartree-Fock calculation, the
mean absolute error for the 16 proteins was 15.08 mhartree for ADMA and 3.99-3.76 mhartree when Mulliken and Lowden charges are
used respectively.



"Accelerating self-consistent field convergence with the augmented Roothaan–Hall energy function" 
X. Hu, W. Yang, J. Chem. Phys. 132, 054109 (2010); https://doi.org/10.1063/1.3304922

Here a new method of SCF convergence acceleration is based off of the DIIS method. It takes the Augmented Roothaan-Hall energy
function and treats it as the object of minimization to obtain linear coefficients of Fock matrices for DIIS, as opposed to
minimizing the commutator of the Fock and density matrices in the original DIIS method. Pulay's DIIS method obtains optimized
linear coefficients for each density matrix by minimizing the orbital rotation gradient based on the commutator matrix of the
Fock and density matrices. However, this minimization does not always produce lower energies if the SCF is not close to
convergence, which EDIIS fixes. However their main critique is that EDIIS (and EDIIS+DIIS) is only accurate for HF SCF. The ARH
DIIS (ADIIS) method is based on the second order taylor expansion of total energy with respect to the density matrix. The linear
coefficients of the density are precise for HF because the HF energy is quadratic in the density matrix, but only an approximate
quadratic expression for KS-DFT is obtained because the exchange-correlation term in KS-DFT is nonlinear in the density matrix.
ADIIS is interpolative, like EDIIS and unlike DIIS.

All molecules had their SCF convergence examined with DIIS, EDIIS+DIIS, and ADIIS+DIIS (acetaldehyde also had just EDIIS and
ADIIS used as well). They examined several molecules, acetaldehyde, a 51-monomer water cluster, a Cd(^2+)-imidazole complex, a
polyanaline peptide, [Ru(tpy)(bpm)(OH2)](^2.5+) and Ru4(CO). Interestingly, EDIIS+DIIS is outperformed moderately or
significantly by ADIIS+DIIS in all systems. Core Hamiltonians of the acetaldehyde and cadmium complexes were used to generate
initial density guesses, whereas atomic density matrices were used as initial guess densities for the other molecules.
Acetaldehyde had a RHF/6-31G* basis set. EDIIS and ADIIS were slower to reach SCF convergence than DIIS, as the energy functions
are not as sensitive to small variations in density matrices (particularly when convergence is close) as DIIS is. DIIS was much
faster to converge, and convergence was slightly faster when ADIIS was used with DIIS, particularly when DIIS was switched on
within 10 SCF iterations. A B3LYP(VWN5)/6-31G* basis set (1275 Gaussian functions) was used for the water cluster. DIIS
moderately outformed AD+D moderately, as the atomic guess was a good initial guess. For Cadmium, DIIS didn't converge for RHF
and B3LYP(VWN5)/3–21G, and AD+D moderately and significatly outperformed ED+D with HF and DFT respectively. The same basis sets
were used for the peptide system. In HF, all molecules converged and DIIS was slightly fastest. In B3LYP, ED+D didn't converge,
and at two different energy thresholds (0.01 and 2.0 AU) AD+D significantly outperformed DIIS (35-40 vs 57 iterations). Finally,
the larger Ru complex has an open-shell (over 120 alpha and beta electrons) fractional electron occupation, with hundreds of
thousands of TIP3P water molecules. Only AD+D yielded a converged energy. The smaller Ru complex didn't converge with ED+D, and
AD+D moderately outperformed DIIS. 
