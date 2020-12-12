"A quadratically convergent Hartree-Fock (QC-SCF) method. Application to closed shell systems"
G.B. Bacskay, Chem. Phys. 1981, 61, 3, 384-404; https://doi.org/10.1016/0301-0104(81)85156-7

This paper presents a procedure to present an orbital optimization procedure that is quadratically convergent for single
determinantal Hartree-Fock wavefunctions. They compare their findings to the Roothaan SCF method with respect to computational
efficiency. Essentially, this method accomplishes a refinement of the Hartree-Fock wavefunction through a modified CI
calculation that is a single excitation. Here, the Hamiltonian matrix includes the interaction between the reference state and
its double excitations. The method, the author states, is directly applicable to UHF and closed-shell RHF cases. For a
spin-orbital basis, the author derives an energy expression (eqn 6), using single determinantal wavefunctions from 'improved'
spin orbitals which are calculated with summations over weighted (by a coefficient C) virtual orbitals. The energy minimization
equations (7a & 7b) includes a two electron integral calculated with 'hole' and 'particle' states (neglecting this term from eqn
6 gives a standard single excitation CI calculation). Equation 11 is the final matrix equation that is iteratively treated in
the QC-SCF procedure. Equation 13 is used to generate a new (improved) set of occupied orbitals to be normalized and
orthogonalized. The energy expression of equation 6 is correct in the first order, convergence will then be quadratic. For a
closed shell RHF case, the wavefunction is defined with spin-adapted configurations. From this, an eigenvalue equation (#20) is
derived. 

The paper outlines another method that complete avoids a four-index transformation step, which changes the order of summations
in the calculations to generate vectors from an arbitrary set of basis vectors. The energy eigenvalue equations are rewritten
(eqns 32 and 33). Rayleigh-Schrodinger Perturbation Theory is then used to generate a set of basis vectors in making a
variational solution to equation 20. The relationship of this method to coupled-perturbed Hartree-Fock theory (CPHF) is
discussed, where a Newton-Raphson form of QC-SCF is derived with CPHF (for QC-SCf, minimizing the expanded energy excitation
value gives the Newton-Raphson equations).

The author first examines H2O (before the latter method was introduced), where the Newton-Raphson equations were solved with
varied level-shift parameters. QC-SCF is compared to standard SCF and SX-CI (single excitement CI) via errors in the values of
the expectation values of the total energy and electron repulsion energy. A 4-21G basis was used, here and in subsequent data
(equillibrium geometry) that shows level-shifted hessians in Newton-Raphson calculations. The diagonalization approach worked
well once the energy was converged to 0.1-0.2 hartree. However results also show that a well damped QC-SCF calculation converges
in four iterations to under 10(^-10) hartree. Using a different approach that avoids the four-index transformation step, closed
shell RHF calculations were done for H2O, CO, and C2H2N^+. A universal level shift was used to avoid divergence in initial
iterations. Results show that QC-SCF are inferior for initial SCF convergence to conventional HF; convergence can be improved
with level shifts and solving Newton-Raphson. This is due to starting orbitals needing to have the correct symmetry for QC-SCF.
DIIS was incorporated into the SCF algorithm. QC-SCF for H2O at 4-31G converged faster (even slightly better than DIIS), given
that the starting orbitals were partly converged by traditional SCF (5 iterations). The CO molecule showed oscillations with
conventional HF calculations. QC-SCF converged readily with the partially converged starting orbitals. For convergence to
10(^-8) au, QC-SCF needed two-thirds of microiterations compared to conventional SCF (comparable to DIIS). Lastly, C2H2N^+
converged in 4 iterations with partially converged starting orbitals (iterations = 8) compared to conventional SCF needing 6,
with QC-SCF needing 18 fewer microiterations. The polarizability of N2 was calculated with QC-SCF (44 orbital Gaussian basis),
and comparable to other results (three iterations needed). 



"Direct inversion in the iterative subspace (DIIS) optimization of open‐shell, excited‐state, and small multiconfiguration SCF
wave functions"
T.P. Hamilton, P. Pulay, J. Chem. Phys. 84, 5728 (1986); https://doi.org/10.1063/1.449880

DIIS, in this work, is applied to more general wavefunctions, such as a high-spin open-shell SCF case, two-configuration ground
state SCF and excited-state singlets, represented by a three determinant expansion. DIIS is capable of converging to excited 
states, due to it being a gradient-norm minimizing method. The original DIIS is recapitulated, however I have done this other 
times summarizing various other papers, so I will skip doing it here. It can be the case that during the course of method, the 
error vectors become nearly linear dependent. C coefficients can become quizzically large, and errors accumulate in the 
interpolated Fock matrix. In the old procedure, the error vectors would have been discarded, and the DIIS procedure would have 
been restarted. However this slows down convergence. They propose a different method, using first order information, where only 
the early iterates are discarded until the condition of the normalized system is improved satisfactorily. The minimum condition 
is modified with equation 5, which corresponds to multiplying the linear coefficients in equation 3 with (1+d), where d is a 
damping factor of ~0.02. Three configuration 2 x 2 CAS wavefunctions are used to describe excited states. This form of the CAS 
wavefunction is maintained because the identification of the excited state is possible and convergence to a wrong state does 
not occur, as opposed to rewriting them (such as in the S1 case) as a two-configuration equation. Their code stores the two 
electron integrals in triplets, which improves the efficiency of the construction of the closed-shell Fock matrix.

They look at water, with one O-H bond stretched, two-configuration GVB with a 6-31G** basis set. Two-configuration GVB were used 
again for other results. They point out that their convergence was much better than a contemporary calculation that used a 
straight first order Fock method (for convergence such that the maximum error matrix element was under 10(^-8) was 14-23 
iterations depending on the state). Ethylene at the 6-31G* level of theory, where the twist angle was varied from 0 to 90 
degrees needed 15-23 iterations (same convergence criterion). At 6-31G**, H2O2 needed 25 iterations to converge (with SCF 
convergence slowing down as final convergence is approached), while 250 iterations were needed when DIIS was not used. The two-
configuration calculation of F2 at 6-311G* needed 14-23 iterations. To my knowledge, they do not mention the kind of initial SCF 
guess they used. Cases involving high-spin-restricted open-shell Hartree-Fock were examined with NO2 and an ethyl radical. The 
NO2 calculations are comparable to contemporary calculations, yet they included a distorted ^2B1 state to show convergence to an 
excited state that is not the lowest of its symmetry. For the Ethyl radical DIIS worked well, except for the 3px and 3pz Rydberg 
states, to which they state that second-order methods may be needed. They then looked at excited singlet wave functions, where 
they compared their vertical transition energies of 20 formaldehyde low-lying states (Rydberg basis) to accurate CI results from 
Harding and Goddard. They achieved 14-37 iterations, however the pi to pi* transition because the wavefunction always assumed 
Rydberg character. They remark that most of the errors in the excitation energy arise from correlation energy in the ground 
state. Lastly, geometries and force constants for the lowest excited singlet and triplet states of formaldehyde were examined. 



"The C^2‐DIIS convergence acceleration algorithm"
H. Sellers, Int. J. Quantum Chem. 1993, 45, 1, https://doi.org/10.1002/qua.560450106

Paywall.



"Linear scaling computation of the Fock matrix"
M. Challacombe, E. Schwegler, J. Chem. Phys. 106, 5526 (1997); https://doi.org/10.1063/1.473575

This is a new method concerning the computation of the Fock matrix, and when used in conjuction with a method concerning the 
computation of the Hartree-Fock exchange matrix of insulators (K matrix computed in O(N) CPU time) [J. Chem. Phys. 105, 2726 
(1996)], produces a linear  scaling algorithm for calculating the Fock matrix. They examine water clusters and proteins with 
RHF/3-21G basis sets. They develop a novel version of Quantum Chemical Tree Code (QCTC), which implements their new method. It 
essentially invovles constructed multipole contraction libraries and multiple translation libraries that have various 
admissibility criteria. For their method, they describe a new method of computing J, the Coulomb matrix. To start, it takes an 
existing method of decomposing the density into families of Hermite Gaussian-Type Functions (HGTF). With this scheme, 
insignificant contributions to the Coulomb matrix and density are done away with, which leaves an O(N) complexes for both of 
those calculations. Then by structuring the density, thresholding and evaluation of primitive ERI contributions to J is 
evaluated by a competitive algorithm. In the calculation of each element contribution to J to within an absolute error 
(thresholding), they use statistical error estimates that are size dependent. This dependence prevents build up of small errors 
that can collectively become significant; these can accumulate with fixed tolerances. They also structure the density by group 
family HGTFs by primitive exponent and sorting on the basis of significance, so that the computation of zero primitive ERIs are 
avoided.   

There results for n water clusters (n = 50, 70, 90, and so on up to 150) were computed with MONDO. Using multiresolution 
improved timings (fast J build, direct J build, tree traversal, and near-field integrals) by 5%. Single point MONDO calculations 
were performed for endothelin, charybdotoxin, and tetramerization monomer of P53 (698 atoms, 3836 basis functions), and the K 
and J matrix results were shown (hours of calculation time). Electrostatic potentials of charybdotoxin and P53 were computed, 
with respective iso-surfaces displayed. They conclude that further research should be applied to improving efficacy/rigor of K 
builds.



"What is the best alternative to diagonalization of the Hamiltonian in large scale semiempirical calculations?"
A.D. Daniels, G.E. Scuseria, J. Chem. Phys. 110, 1321 (1999); https://doi.org/10.1063/1.478008 

"Fragment molecular orbital method: an approximate computational method for large molecules"
K. Kitaura, E. Ikeo, T. Asada, T. Nakano, M. Uebayasi, Chem. Phys. Lett. 1999, 313, 3-4, 701-706; https://doi.org/10.101
/S0009-2614(99)00874-X

"The trust-region self-consistent field method: towards a black-box optimization in Hartree-Fock and Kohn-Sham theories"
L. Thøgersen, J. Olsen, D. Yeager, P. Jørgensen, J. Chem. Phys. 121, 16 (2004); https://doi.org/10.1063/1.1755673

"Converging Self-Consistent Field Equations in Quantum Chemistry - Recent Achievements and Remaining Challenges"
K.N. Kudin, G.E. Scuseria, ESAIM:M2AN 2007, 41, 2, 281-296; https://doi.org/10.1051/m2an:2007022

"Generalized Energy-Based Fragmentation Approach for Computing the Ground-State Energies and Properties of Large Molecules"
W. Li, S. Li, Y. Jiang, J. Phys. Chem. A 2007, 111, 11, 2193–2199; https://doi.org/10.1021/jp067721q

"The augmented Roothaan-Hall method for optimizing Hartree-Fock and Kohn-Sham density matrices"
S. Høst, J. Olsen, B. Jansík, L. Thøgersen, P. Jørgensen, T. Helgaker, J. Chem. Phys. 129, 124106 (2008); https://doi.org/10.106
/1.2974099 

"Perturbative total energy evaluation in self-consistent field iterations: tests on molecular systems"
Y.A. Zhang, Y.A. Wang, J. Chem. Phys. 130, 144116 (2009); https://doi.org/10.1063/1.3104662

"LISTb: a Better Direct Approach to LIST"
Y.K. Chen, Y.A. Wang, J. Chem. Theory Comput. 2011, 7, 10, 3045-3048; https://doi.org/10.1021/ct2004512

"An analysis for the DIIS acceleration method used in quantum chemistry calculations"
T. Rohwedder, R. Schneider, J. Math. Chem. 2011, 49, 1889-1914; https://doi.org/DOI 10.1007/s10910-011-9863-y

