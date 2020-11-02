"Starting SCF Calculations by Superposition of Atomic Densities"
J. H. Van Lenthe, R. Zwaans, H. J. J. Van Dam, M. F. Guest, J. Comput. Chem. 2006, 27, 8, 926-932; https://doi.org/10.100
/jcc.20393

   Start SCF calculations by addition of atomic electron densities. Orbitals are generated without first needing to complete
closed shell HF or DFT. They have a nice flowchart of the procedure: Determine atomic configuration using optional user input
(charge, orbital occupations, spin states, etc) -> atomic scf -> compute density, with optional user input, and expand from
atomic to molecular basis -> make molecular block diagonal density by adding atomic densities -> calculate molecular Fock matrix
with molecular density -> diagonalize the Fock matrix to obtain MO's -> order orbitals by energy and optional user input to
generate start orbitals for molecular (scf) calculations. 
   
   They tested a G2 Set of 200 molecules, and 23 transition metal complexes with a 6-31G* basis and also a Ahlrichs TZVP Basis.
Their data is comprised of (E1 - Econverged) and 'SCF iteration characteristics'. 
In the G2 set the SAD guess outperforms HCORE method (obtain start orbitals by diagonalizing one-electron Hamiltonian) with
regard to accuracy of initial energy and number of cycles that lead to convergence. 
HCORE largely is terrible (needed an average 36% more cycles for HF convergence and 31% more cycles for DFT convergence). 
MINGUESS method (make minimal basis SCF, and using the standard minimal basis start, and projecting the resulting orbitals onto
the actual basis) and SAD are comparable, but SAD is slightly better in terms of number of cycles and intial cycle deviation
from converged (MINGUESS needed on average 5% more cycles to converge for HF SCF.). 
For the two transitional metal basis', HCORE failed to converge in most cases, MINGUESS needed 7% more iterations and 16% more
iterations respectively to converge for the 6-31G* and the larger Ahlrichs basis, relative to SAD.  

G2 Test set ref: J. A. Chem. Phys. Lett. 1997, 270, 419



"Assessment of Initial Guesses for Self-Consistent Field Calculations. Superposition of Atomic Potentials: Simple yet Efficient"
S. Lehtola, J. Chem. Theory Comput. 2019, 15, 3, 1593–1604; https://doi.org/10.1021/acs.jctc.8b01089

The idea behind the SAP guess is that they are trying to make an atomic guess from a potential that can give correct atomic
electron density. This guess accounts for interatomic interactions where SAD does not. They make two points on implementing this
guess, first that atomic potentials are obtained from calculations near the limit of the basis set, and second that forming the
SAP potential at an arbitrary point involves a summation over tabulated atomic potentials, which can be narrowed to a finite
range. They then develop their theory with an extensive mathematical argument. 

Molecular Data Set:
"In the present work, we study the 183 nonmultireference molecules from the high-level W4−17 test set of first- and second-row
molecules (^87), which we furthermore augment with a data set composed of 50 transition metal complexes from refs 10 and 88
(referred to as TMC), as well as 28 complexes containing third or fourth period elements from the MOR41 database of single
reference systems (^89)"

Generation of Atomic Potentials:
Atomic potentials generated by KS-DFT calculations, which is a part of "GPAW program package (^100, 101) and which uses the
LIBXC library (^102) to evaluate the exchange-correlation function." Kohn-Sham equations used for self-consistent solutions in
atomic program. 

Results:
They assess their inital guesses by calculating the projection of the guess orbitals onto the converged SCF wavefunction, and
they then calculate a fraction quantity with that. Check out equation 14 and 15. The SAP guess is the best compared to numerous
other guess methods, used with HF/aug-pcseg-2 and HF/pcseg-0 basis sets; it yielded excellent guess wavefunctions in conjuction
with exchange functionals such as Chachiyo (106), CAP (105), and LDA (103, 104). The Extended Huckel guess (which is really an
approximation of the SAP, in that the atomic Huckel Hamiltonian is like an atomic Fock that has an orthonormal single-center
atomic basis) had comparably good results. The SAP guess showed more scatter than their SADNO guess (Take the SAD guess and
diagonalize it's density) and Extended Huckel.
    


"Advanced initial-guess algorithm for self-consistent-field calculations on organometallic systems"
G. Vacek, J. K. Perry, J.-M. Langlois, Chem. Phys. Lett. 1999, 310, 1-2, 189-194; https://doi.org/10.1016/S0009-2614(99)00722-8

This paper uses an algorithm developed using Ligand Field Theory to produce inital guess wavefunctions for organometallic
systems. The guess separates the molecule into metal and ligand fragments where the user can assign spins and formal charges,
estimate ligand-field splitting of the non-bonding d orbitals of the T.M.'s, include effects of unoccupied metal s orbitals, and
estimate d-d repulsion of metal orbitals if necessary. This all gets incorporated into their effective Hamiltonian, which is
diagonalized for new a wavefunction, which gives new overlap matrix. 

They look at HF and DFT SCF for 88 organometallic molecules using exclusively a LACVP** basis  set. The initial guess
wavefunctions produced by their algorithm led to ground state convergence 88% of the time for HF and 92% for DFT. They compared
with Huckel guess, with ground state convergence occurring 14% and 9% of the time, respectively. Traditional AO overlap guess
converged 39% and 42% respectively.
    


"Use of promolecular ASA density functions as a general algorithm to obtain starting MO in SCF calculations"
L. Amat, R. Carbo-Dorca, Int. J. Quantum Chem. 2001, 87, 59-67; https://doi.org/10.1002/qua.10068

Their guess is called the Atomic Shell Approximation, which is described as a way to fit first order density functions to a
linear combination of spherical functions. Positive definite expansion coefficients are used to fit features of a real
probablility distribution. Atoms in a molecule are represented as neutral, and of spherical shape. An initial guess molecular
density can be obtained from a summation of these isolated, spherically symmetric atomic densities. On a theoretical basis, they
calculate their initial Fock matrix by adding their Promolecular ASA density to a one-electron core Hamiltonian. 

They test a fairly good range of organic, organometallic, and inorganic geometrically optimitized molecules (at the H-F level),
and they use a 6-311G basis set implemented in Gaussian. The reduction in the number of cycles isn't impressive compared to use
of core Hamiltonian only, and the use of CNDO, INDO, and extended Huckel Hamiltonians. Since this method isn't impressive, they
conclude that this initial MO method is good on the sole basis that it is a new method of making an initial MO guess for SCF
calculations.
    


"Use of exchange maximization to generate starting vectors for self-consistent field calculations on metal cluster/adsorbate
systems"
R. J. Buenker, J. L. Whitten, E. I. Izgorodina, H.‐P. Liebermann, D. B. Kokh, J. Comput. Chem. 2002, 23, 10, 943-949;
https://doi.org/10.1002/jcc.10094

They derive Localized Molecular Orbitals (LMOs) by exchange maximization in regard to all atom centered basis functions to
generate a good SCF inital guess on metal cluster systems, where there aren't any defined chemical bonds. They also report that
they can obtain better starting vectors. They diagonalize a Fock matrix, computed with null electronic field, to remix a subset
of LMOs with largest exchange eigenvalues. They use these LMOs in inital cycles, and then extend the space of the
diagonalizations in further cycles.

They look at a handful of systems; they have results for a 20 Ni atom cluster, along with a platinum, two silver, and
vinylidene/Ni adsorbate systems. They don't give comparative results, only specific results. For the Ni cluster, they provide
data on SCF energies and exchange energies for two cases, the first with a 'pure' exchange operator, and the other with an
exchange operator augmented with a fraction of kinetic
energy. They report effective SCF convergence. For Pt(97)CO, they provide iteration and respective SCF energies using exchange
maximization technique. Same with Methyl Nitrate adsorbed on Ag(43). They then report improved convergence, for formaldehyde
adosrbed on silver and vinylidene/Ni systems, by expanding the number of LMOs in the basis, and reporting the number of
iterations needed for convergence. Converged SCF energy improved each time they did this, without linearly increasing the number
of iterations needed for convergence.
    


"Diagonalization-free initial guess to SCF calculations for large molecules"
Z. Szekeres, P. G. Mezey, P. R. Surján, Chem. Phys. Lett. 2006, 424, 4-6, 420-424; https://doi.org/10.1016/j.cplett.2006.04.089

Their HF and DFT SCF guess is the adjustable density matrix assembler (ADMA) method. Its an approximation of molecular density
by a subdivision of the molecule into fragments. Neighboring fragments are then grouped into 'parent molecules' that overlap.
The cutoff size of these 'parent molecules' can be adjusted, and fragments can be chosen via chemical intuition (functional
groups, amino acid residues in proteins, etc). The dangling bonds, connecting the parent molecule to the molecule itself are
closed by attaching some hydrogen atoms, or any other atoms. Compute the density of the parent fragments by computing the
density of the constituent fragments (exclude elements that close up the dangling bonds), and then sum the parent fragments. In
the summation there is a weighting to the fragments which limits overcounting the overlapping parent fragments (equations 3 and
4).

They showed faster convergence than Extended Huckel Theory in traditional SCF (HF and DFT), and while looking at P iteration.
Generally, ADMA showed a lower energy error per iteration step. In an hexagonal ice cluster containing 216 atoms, they saved one
iteration using Pople's 6-31G* basis set (5 and 6 iterations respectively). In an Alanine Oligomer containing 413 atoms, they
converged in 4 steps while EHT converged in 16. They show 4 (ice) and 9 (alanine) cycle savings respectively when converging to
idempotent P matrices. They use GAMESS and MUNGAUSS packages. 



"Improvement of initial guess via grid‐cutting for efficient grid‐based density functional calculations"
J. Lim, S. Choi, S. Kang, J. Kim, K. Hong, W. Y. Kim, Int. J. Quantum Chem. 2016, 116, 19, 1397-1403; https://doi.org/10.100
/qua.25193

Their SCF "grid cutting" intial guess method they say is specialized for grid-based DFT calculations. Initial density, as well
as an inital set of eigenvectors are generated through DFT in a simulation box that is smaller than the full box. The idea is
that electron density decays exponentially as the radius to the nucleus increases. More specifically, they do make a (DFT) pre
SCF guess by calculating a SAD density in their "inner box" simulation. They mention that their GC guess takes a longer time to
generate an intial density and set of orbitals than other methods, but that this is compensated for by noted reduction in cycles
needed for convergence. 

They carried these calculations out on 117 molecules (only in singlet state) in the G2-1 test set. They look at two methods of
mixing, Pulay and Broyden. With respect to the a SAD initial guess, GC method needed a 20% smaller number of steps for SCF
convergence, as compared to the extended Huckel guess needing 7-15% more steps when the maximum number of iterations for
diagonalization at each SCF step (Imax) was set at 100. When it was set to be 10000, Broyden overall was the most effective
mixing method, and GC was 35% faster than SAD to convergence. LCAO inital guess method was also considerably faster than SAD,
but slower than GC. Extended Huckel was slower than SAD. They looked at virtual orbitals as well, and GC was around 5% faster
than sad at calculating 5, 10, and 20 excited state orbitals for 40 molecules in the G2-1 set. For large molecules such as
Fullerene and Porphyrin, they saw nearly a 50% reduction in the number of cycles, and a less than 15% reduction in time. They
found an optimal inner box size of 70%, but other sizes were tested. 
 
 

"Improving self-consistent field convergence by varying occupation numbers"
A. D. Rabuck, G. E. Scuseria, J. Chem. Phys. 110, 695 (1999); https://doi.org/10.1063/1.478177

This paper looks at two methods that fractionally occupy orbitals around Fermi energy during SCF cycles. Occupations at
convergence are forced to be zeros and ones (they use a temperature factor that cools to 0 K at SCF convergence). They report
improved convergence at no significant overhead. Specifically, they describe a Fractional Occupation Number (FON) and a pseudo
FON that can help with other convergence methods such as DIIS (build an extrapolated Fock matrix based on error matrices formed
during SCF iterations) and level shifting (increase HOMO-LUMO gap by shifting virtual orbitals to a higher energy). 

They use four different .../3-21G basis sets and four different .../6-31G** basis sets for CrC, Cr2, TiC2, NiC2, PN, and H2O. In
most data points, FON and pFON helped either slighly or not at all. Some cases showed a significant reduction in the number of
cyles, a lot of cases showed a minor reduction or no reduction, in other cases the number of cycles increased. As mentioned
before, they test each basis set using a no FON, FON, and pFON when used in conjuction with level shifting and DIIS methods.
They claim a radical reduction in the number of cycles overall, but the spread of the data suggests more of a mixed bag to me. 



"Extended Hückel and Slater’s rule initial guess for real space grid-based density functional theory" (CANNOT ACCESS)
M. Lee, K. Leiter, C. Eisner, J. Crone, J. Knap, Comput. Theor. Chem. 2015, 1062, 24-29; https://doi.org/10.1016/j.comptc.2015.03.011

   My Iowa State login can only provide me this paper's abstract. It's contents over Extended Huckel and Slater's Rules intial guesses for SCF calculations would be interesting to look at.