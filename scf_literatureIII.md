"A dynamical damping scheme for converging molecular SCF calculations"
M.Z. Zerner, M. Hehenberger, Chem. Phys. Lett. 1979, 62, 3, 550-554; https://doi.org/10.1016/0009-2614(79)80761-7

A method of accelerating SCF convergence is presented, where methods of extrapolation and damping are combined; it is an
extrapolation technique that uses Mulliken gross atomic populations in the calculation of a damping factor. In each cycle the
factor is adjusted so that obtained densities correspond to P's that are extrapolated. This allows them to avoid independent
extrapolation on each of highly interrelated elements of D. The idea behind the calculation of their damping factor is that for
any iteratively calculated quantity, there will be different input and output values. With input as x-coordinate and output as
y-coordinate, a graph would show a 45 degree line at convergence (input would equal output). They derive an expression for the 
extrapolated quantity in terms of a damping factor, alpha, where alpha = m/(m-1). This method (slightly modified) was
implemented in an ab initio UHF program. The Mulliken gross populations were extrapolated, and the alpha values were weighted
based on the number of basis functions for a given atom and averaged. They were then used to damp the densities. 

By forming weighted averages of damping factors based on atomic characters of basis functions and by treating different
irreducible representations separately, this method has been generally shown to adjust automatically according to the inherent
convergence difficulties in each system. The excitation of CuCl4(^2-) from ground state was extensively studied. This method
lead to convergence after 13 cycles, where other contemporary methods had failed. One result showed that after 4 cycles, energy
is 7 x 10(^-4) a.u. from the converged result. The calculation of CO ground state was examined with a triple zeta basis set, and
convergence was reached after 18 cycles where contemporary methods were much slower. Their damping factor was reduced from 0.5
to 0.4, as the former damping factor slowed down convergence.



"A general relation between the intrinsic properties of the SCF Hartree-Fock calculations and the stability conditions of their
solutions"
J.C. Facelli, R.H. Contreras, J. Chem. Phys. 79, 3421 (1983); https://doi.org/10.1063/1.446190

For this paper, the authors are attempting to identify the relationship between a convergence requirement for iterative SCF for
even the most general one determinantal HF wavefunction, which accords with aufbau principle and the attending stability
conditions. At the outset, they state that previous research (Cizek and J. paldus, J. Chem. Phys. 47, 3976 (1967)) has indicated
all converged SCF calculations to a minimum of the energy hypersurface of closed-shell configurations, although they could be a
local or absolute minima. They then do some extensive derivations of HF wavefunctions and SCF equations that lead to their
subsequent conclusions, which I describe.

The results of this paper confirms that when an SCF calculation converges, the solution is stable in the given Fock subspace,
and a local minimal of the energy hypersurface is reached. Given convergence in a certain Fock subspace, they stress that
nothing can be ascertained about solutions obtained in other Fock subspaces. An example is that a closed-shell SCF convergence
solution will guarantee there are no imaginary instabilities, there is no information relating to nonsinglet instabilites. Given
the relationship between the stability of the solution and the SCF convergence, enforcing convergence in an intrinsically
nonconvergent result could be erroneous. The authors attribute slow convergence to the need of configuration interations as
opposed to one-determinantal representations. 



"Simplified method for calculating the energy of weakly interacting fragments"
J. Harris, Phys. Rev. B 1985, 31, 4, 1770-1779; https://doi.org/10.1103/PhysRevB.31.1770

Though not specifically derived as a method to calculate an intial guess density for SCF calculations, the method has been
utilized for this purpose in subsequent quantum chemistry literature. The method was proposed to tackle systems that require
calculations of that the type that, for example, involve interactions of a molecule and a surface. The author derives a
expression that expands the Kohn-Sham energy about a reference density, typically a sum of overlapped frozen fragment densities,
with a reference potential being generated from this reference density. There has to be a variational calculation of an
eigenvalue sum of an independent-electron Hamiltonian whose potential energy term is the reference potential. From these
operations equation 2.17, Harris' energy function, results. This is an approximate calculation of energy that does not require
determination of molecular density and coulomb potential. To make an expression that determines interaction energy (2.18),
energies E1 and E2 of the overlapped frozen fragments F1 and F2 must be subtracted from 2.17, assuming F1 and F2 are stable at
an internuclear radius R = infinity and that there is a continuous change of energy to the summed isolate fragment energies. For
systems where this isn't the case, he introduces an expression with an approximation of frozen core wavefunctions, which can
allow for the elimination of considerations of core electron eigenvalue shifts from 2.18 (eqn 2.21). He derived equation 2.17
(and subsequent) for a system of two fragments, but it can be generalized to N fragments. This can't be said of 2.18 or others. 

He examines homonuclear dimers, He2, in comparison to the work of Gordon and Kim's work (R.G. Gordon, Y.S. Kim, J. Chem. Phys.
56, 3122 (1972)). From eqn 2.18, he introduces a sum of second order eigenvalue shifts for a value that, as far as I can tell,
free particle kinetic energy to yield eqn 3.6 and applies the eigenvalue sum with an LCAO basis for He2 (3.13). They present an
energy curve (Potential (V) versus Radius (R)) which displays energies from different theories. The Harris energy is the most
attractive (similar result to full self-consistent KS-LDA calculation), more than the Gordon-Kim energy, and a CI energy of Liu
and McLean (B. Liu and A.D. McLean, J. Chem. Phys. 59, 4557 (1973)). He notes that LDA was a relatively poor approximation for
exchange and correlation energy/free-particle kinetic energy alone. More diatomic systems are examined (singlet dimer states
only) with 2.21, in terms of the spectroscopic parameters of well depth, equillibrium separation, and vibration frequency. They
are compared to experimental values to highlight accuracy. Be2 was excellent, and was attributed to small density differences
due to closed shell structure, C2 (doubly pi bonded ground state) was least accurate due to conjugation (sum of overlapped
atomic densities was a poor approximation). Absolute binding energy was too large for C2, F2, and Cu2, while too small for N2,
which has been attributed to the sign of the difference density at the bond center. Ground state energy curves of N2 and Be2 are
shown and analyzed, calculated with eqn. 4.1. The results tract with the accuracy of binding energy, Be2 had a high and accurate
binding energy, while N2 had a low and less accurate one due to improper inclusion of nuclear screening in the calculations. To
conclude, errors are negligible in closed shell dimers such as Be2, and errors increase in covalently bonded dimers. With these
dimers, he notes that the Harris functional is still accurate enough to be useful as an alternative to a complete Kohn-Sham
solution, if the latter is not feasible.



"Ab Initio Quality Electron Densities for Proteins: A MEDLA Approach"
P.D. Walker, P.G. Mezey, J. Am. Chem. Soc. 1994, 116, 26, 12022-12032; https://doi.org/10.1021/ja00105a050

This method of assembling molecular electron densities from fragments is an earlier version of the ADMA method, developed by the
same author. Ab initio quality electron density distributions for large molecules with the Molecular Electron Density Lego
Assembler from superimposing the (fuzzy) densities of smaller fragments. There is mutual interpenetration of the fragment
densities, hence the analogy of the joining of two legos. They have shown that the electron distribution is quantitatively very
similar to an ab initio density obtained from a 6-31G** basis. Despite requiring a small fraction of computational time, the
density provided by this method is more accurate than direct ab initio results performed with a 3-21G basis. If a database of
densities is constructed at a higher level of theory, they claim, there will be a complimentary increase in the accuracy of the
MEDLA density. Similar to the ADMA approach, a molecule's density can be represented as a summation of scaled fragment
densities. Fragments are scaled by 1 when both orbitals are centered on that fragment, 0.5 when only 1 is, and 0 if neither are.
They have MIDCOs (molecular isodensity contour) for each of their molecules, noting that the contours that have higher density
thresholds show a tight skeleton of the molecule while inaccurately portraying the electron density (which is actually fuzzy of
course). MEDLA computational time increases linearly with the number of fragments. A contemporary direct ab initio calculation
of the gene 5 protein would have taken a century of CPU time on a CRAY supercomputer, but their lab computer did the MEDLA
calculation for the protein in 21 minutes.

They used a database of 21 fragments (some fragments have multiple versions in the database, so the accuracy of the MEDLA
density can be as accurate as possible) to compute the densities of polypeptides and proteins that are composed of the residues
of the 20 amino acids. A numerical similarity of MEDLA and 3-21G densities with 6-31G** (Sxy) was obtained, and a point-by-point
grid comparison was calculated (L(a,a',X,Y) where L is 1 minus the average relative difference when X is compared to Y at a and
a' thresholds). Where L is the most sensitive (a and a' at 0.01 and 0.001 au), MEDLA was comparatively more accurate. For the
dipeptide system (glycinal alanine), a MEDLA MIDCO was a highly accurate representation of the ab initio 6-31G** MIDCO, better
than an ab initio 3-21G MIDCO. Hydrogen bonding systems were tested, a helical peptide that had a hydrogen bond between the
first and fourth amino acid (both ends 'tied off' with methyl groups), and the MEDLA MIDCO was more accurate to 6-31G**. This
behavior was replicated in metenkephalin (which had nonbonding interactions between -SCH3 and -Ph fragments). It should be
pointed out that the 3-21G has a slightly lower S value than MEDLA, and the L values are fair at low sensitivities. They show a
MEDLA density for an alpha helix polypeptide with 13 amino acids ((Gly)2-(Ser)10-Gly), as well as bioactive peptides (two
enkaphalins, leuenkephalin and metenkephalin) and low energy conformers of proteins (bradykinin, crambin, and gene 5 protein).
The helical nature of the first polypeptide as accurately captured by MEDLA for various thresholds in a <5 minute calculation.
The enkaphalins are 5 amino acids each, and they differ only in the last (C terminus) amino acid. Not all chemical features were
captured in the MEDLA MIDCOs. Bradkinin had 59 MEDLA fragments, and the low density MIDCOs were shown, and the authors speculate
that they are more accurate than fused sphere Van der Waals surfaces. They needed 11 minutes to get a MEDLA density for crambin
(46 amino acids). The MIDCOs for this and gene 5 protein (87 amino acids) are both detailed, with bonding and nonbonding
interactions captured for g5p. 



"Can we outperform the DIIS approach for electronic structure calculations?"
E. Cancès, C.L. Bris, Int. J. Quantum Chem. 2000, 79, 2, 82-90; https://doi.org/10.1002/1097-461X(2000)79:2<82::AID
QUA3>3.0.CO;2-I

Several SCF algorithms are explored, such as the Roothaan, level-shifting, and direct inversion of iterative subspace (DIIS)
algorithms. New SCF algorithms are presented in this paper, and they have been mathematically proven to converge to a critical
point of the energy functional, and that calculations that have been performed with RHF methods show comparable and sometimes
superior efficiency to DIIS w.r.t. CPU time and memory usage. The author calls them relaxed constraint algorithms (RCA), and
they are applied to closed-shell HF problems exclusively in this paper. One of the problems with the Roothann algorithm is that
it converges to a stationary point in the HF energy (in good cases), or the final convergence behavior is an oscillation between
two states, neither being solutions of HF equations. The response to this behavior is to instead minimize E by relaxation, where
minimizations are done alternatively with two different densities (D and D-tilde), which is proven (under some fulfilled
assumptions) to allow for convergence to a critical point. Oscillation (doesn't converge to stationary point) can still occur,
which can be remedied by adding a penalization term to off-diagonal pairs to induce convergence to a point where D-tilde = D.
Incidentally, this is what the level-shifting (in the self-named algorithm) parameter b does, though those authors stated that
it prevents virtual and occupied orbitals from mixing together. This paper's authors improved the algorithm and proved there is
a nonnegative integer level-shift parameter that decreases HF energy to a converged stationary value, for any inital guess. 

Now to introduce an RCA called the optimal damping algorithm (ODA). Principally, the constraints DSD = D are relaxed (noted by D = D-tilde), and PN-tilde is convex. Convergence arises from HF energy decerasing at each step, and DSD = D constraint returns at
convergence. The algorithm is initialized by choosing an initial guess, assemble an inital fock matrix with guess, and compute
an inital E. To iterate, diagonalize Fk-tilde and construct Dk+1 via aufbau. If Dk+1 - Dk is "small enough", stop iterating and
compute E. Then assemble F(Dk+1) and compute Ek+1. Compute s, c, E(^2e)k+1-tilde, and finally set k = k+1 to repeat. They
explore an iterative subspace method as well, where an updated density matrix is assembled with all other stored density
matrices.

They compare the ODA and DIIS SCF convergence accelerators for calculating RHF ground states of various systems. They use two
guess methodologies, a semiempirical method from (Frisch, A.; Frisch, M. J. Gaussian 98 User’s Reference;
Gaussian Inc.: Pittsburgh, PA, 1999), and the core Hamiltonian diagonalization guess. The calculations were carried out in
Gaussian98, and efficacy is calculated by computing the log of the difference of the energy for the given iteration's density
and the converged RHF ground state energy (completed for each iteration). They examined Acetaldehyde, Cr2, and the (E)-isomer of
N-methyl-2-nitrovinylamine (CH3—NH—CH=CH—NO2) with RHF/6-31G(d), RHF/6-31G, and RHF/6-31G(d) basis sets respectively. For the
first molecule, ODA an DIIS had similar convergence properties, with ODA being more efficient with initial iterations and DIIS
being more efficient with later iterations that are closer to convergence (DIIS has consistently shown this property in the
papers I have seen). DIIS was faster to converge with the semiempirical guess (by 0.4s) than the diagonalized core Hamiltonian
guess (by 0.1s). With respect to Cr2, both ODA and DIIS found aufbau solutions, but ODA converged to a lower energy minimum than
DIIS (DIIS energy is 0.253 a.u. higher). DIIS converged faster than ODA, and the convergence patterns didn't change for either
methods with either guess (DIIS was about 1.2s faster). For the last molecule, DIIS converged faster than ODA with the
semiempircal guess (by 1s, but the graph is cut off so I extrapolated). ODA converged with the cruder core Hamiltonian guess
(after 2.6s), while DIIS failed to converge. In the conclusion, application to open-shell systems and KS-DFT methods are
discussed.



"Electrostatically Embedded Many-Body Expansion for Large Systems, with Applications to Water Clusters"
E.E. Dahlke, D.G. Truhlar, J. Chem. Theory. Comput. 2007, 3, 1, 46-53; https://doi.org/10.1021/ct600253j

This paper implements a many body expansion with predefined point charges to represent the electrostatic field created by the
other molecules. They use their method to calculate binding energies for water trimer up to pentamer. The total energy of a
system can be represented as a sum of many-body energies. They derive an approximate energy expression (eqn 6) for their three
body method, and another for a pairwise addition method (eqn 5). All equations are linear in the energies. The added point
charges (the use is denoted by EE - 'Electrostatically Embedded') are put in place to represent the other N-n particles, where n
is an n-mer. They use AM1 Mulliken charges in two ways. The first way is to determine a charge representation for the whole
cluster, and represent the other N-1, N-2, or N-3 water molecules with corresponding point charges from the charge
representation calculation. The second way is to determine the gas phase charges for a monomer, and represent the other N-nmer
with gas phase point charges (abbreviated AM1M). The monomers of this method are represented with three choices of monomers:
B3LYP/6-31G** Mulliken Charges (abbreviated B3LYPM), B3LYP/6-31G* Class IV charges (CM4M), and TIP3P molecular mechanics charges
(TIP3P). AM1 calculations were carried out in Gaussian03, other charge calculations were carried out in MNGSM v.6.0. They give a
summary of their method: first determine the partial atomic charges for all N molecules (one of two ways), calculate the energy
of each n-body cluster in the expansion by represent the N-n molecules with point charges centered on the nuclei in the method
previously determined. Lastly, use a two/three body expansion equation (5 or 6) to calculate the total energy of the system. All
of their partial charges correspond to equillibrium geometry of water monomer.

Their 15 water cluster structures come from Monte Carlo and molecular mechanics simulations of bulk water and ice, and also
optimizations in the gas phase. The efficacy of their method's test results are how accurately total energy from full
calculations can be reproduced, at a given level of theory. The trimers and tetramers were bounded by optimized monomers, the
pentamer was taken from a simulation of water. They look at different levels of theory: BLYP, B3LYP, PBE, PBE1W, and MP2. For
each density functional of water, the MG3S basis set is used. Their results are presented in terms of mean errors, with EE-PA
B3LYPM and EE-PA-TIP3P showing a 10 fold error reduction. Their worst result (an EE-PA with AM1M charges, the average mean
unsigned error is reduced by a factor of 2). Their best results use gas phase monomer charges, and they try to optimize these
charges. When these optimized charges' errors to the other EE-PA methods using monomer charges. They did not get significant
improvement over TIP3P charges, and didn't see large differences in error for the other EE-PA monomers with various gas phase
point charges. They then proceeded with using Mulliken charges (B3LYP/6-31G*). They then look at the mean unsigned errors of the
three body energies and the electrostatically embedded three body energies. The EE-3B method, in all cases, gives a mean
unsigned error (MUE) of less than 0.1 kCal/mol, an MUE of 0.01 kCal/mol was observed in the MP2/aug′-cc-pVTZ. Once again the
minimum error decrease was 2-fold when point charges are introduced, and a maximum of 10-fold for MP2. For the nine methods,
average MUE was 0.05 kCal/mol, compared to 0.31 kCal/mol for EE-PA for tetramers and pentamer. Lastly, they looked at large
water clusters (n=21 monomers) at the MP2/aug′-cc-pVTZ basis. EE-PA gave a 2.97 kCal/mol error relative to true energy, while
EE-3B gives a 0.38 kCal/mol error. Their EE-3B method computationally scales at N^3, while EE-PA scales at N^2. 



"Communication: Linear-expansion shooting techniques for accelerating self-consistent field convergence"
Y.A. Wang, C.Y. Yam, Y.K. Chen, G.H. Chen, J. Chem. Phys. 134, 241103 (2011); https://doi.org/10.1063/1.3609242

The method of SCF convergence acceleration presented here are two linear expansion shooting techniques (LIST), direct and
indirect (LISTd and LISTi). DIIS is the predominant convergence acceleration technique, but occassional failure arises from the
coefficients that minimize the norm of the error vector (obtained as the commutator of F and D) are not optimal for the linear
expansion of the Fock matrices, or for the convergence of the totoal energy. In situations where DIIS fails, the authors state
that AD+D is generally superior to ED+D in virtually every case. For their LIST methods, the authors derive their own energy
expression from the corrected Hohenberg-Kohn-Sham funcitonal, and introduce a linear expansion of iterated densities. The
expression is compact, but linear in expansion coefficients, so they have to minimize via a shooting technique in numerical
analysis. Their direct LIST method constitutes solving a matrix equation, AC=O, where A has elements of an energy expression
built with iterated output densities, C is the expansion coefficient matrix, O is a zero matrix. LISTd will start to fail close
to convergence as the density matrix residual error becomes small. LISTi leads to another matrix equation, but instead of doing
a direct shooting with each linear expansion separately, they do a disting shooting scheme between two linear expansions, an
indirect approach. 

They implemented the methods in NWChem code, the input Fock matrix per iteration was a linear combination of four output Fock
matrices. The cHKS energy functional was evaluated at each iteration for total energy. They looked at six molecular systems:
water monomer (LDA/6-31G), benzene (LDA/6-31G), silane molecule with extended Si-H bond (LDA-6-31G*), cadmium(II)-imidazole
complex (B3LYP/3-21G), tetrahedral uranium fluoride (B3LYP/LanL2DZ), and Ru4(CO) (B3LYP/LanL2DZ). The cadmium SCF calculations
utilized the core Hamiltonian initial guess, whereas all others used the inital atomic guess for SCF calculations. Generally,
they show that LISTi outperforms other methods, and converges in all systems. For water and benzene, LISTi and DIIS had
comparable performances in terms of the number of iterations needed to converge (LISTi being slightly better in both cases).
DIIS failed to converge in other systems, and LISTd could converge for the cadmium system (worse than their LISTi by 12
iterations). LISTd requires the most iterations in the water and benzene systems, and didn't converge in other systems besides
the cadmium complex.



"Comparison of self-consistent field convergence acceleration techniques"
A.J. Garza, G.E. Scuseria, J. Chem. Phys. 137, 054110 (2012); https://doi.org/10.1063/1.4740249

In this interesting paper, ADIIS and LIST methods for accerlerating SCF convergence is compared to EDIIS+DIIS. The results in
this paper do not support the conclusion that ADIIS and LIST methods are not superior to ED+D. In all the tested cases of this
paper, ED+D surpasses LIST. They mathematically prove that the energy functions of ADIIS and EDIIS are identical in the case of
HF wavefunctions. Since EDIIS results actually match the ADIIS results of J. Chem. Phys. 132, 054109 (2010), they actually
conclude that EDIIS was not correctly implemented in that paper. W.r.t. DIIS, an inherent flaw is that, while [F, D] = 0
(commutator) is a good condition for SCF convergence, minimizing ||[F, D]|| (frobenius norm) does not force convergence. In 
EDIIS, the HF energy functional is minimized, and the coefficents of linear expansion are subject to c ∈ [0, 1], thus allowing
D-tilde ∈ ̃Pn-tilde. The fact that EDIIS can only interpolate between matrices in linear expansion, chances of finding a
minimum decrease as convergence is approached. Therefore, it is best to use EDIIS to start from a (poor) guess (as the DIIS
error is greater than 10(^-1) a.u.), and use DIIS as convergence is approached (when DIIS error goes below 10(^-4) a.u.). The
LIST methods are derived from the corrected Hohenberg-Kohn-Sham energy functional, and lead to a system of linear equations that
are similar to DIIS (eqns 8 and 9). 

The calculations in this paper were carried out in Gaussian-09, and the Harris functional
diagonalization was used for an initial guess. In all methods, the Fock matrix is extrapolated and diagonalized to yield an
idempotent density. This is approximate for KS-DFT because Fock matrix is not linear w.r.t. density. They note that since the
Fock matrix is not linearly dependent on density in KS-DFT, ADIIS does not exactly minimze the ARH energy functional (it isn't
even clear that AD+D is better than ED+D in DFT, which they show it really isn't in their results).

High symmetry systems that contain transition metal and actinide species were studied where competition for certain bonds lead
to spontaneous structure changes. Specifically, cadmium-imidazole ([Cd(Im)]^+2), Ru4(CO), SiH4, CrC, and UF4. All calculations
were closed shell. On the whole, ED+D was the most effective convergence accelerator. For the cadmium species at RHF/3-21G
theory level, all four methods display similar convergence rate. But ED+D was faster for the cadmium species, CrC, and SiH4 at
B3LYP/3-21G, B3LYP/6-31G, and LDA/6-31G* theory levels respectively. List and ED+D were compared with UF4 and Ru4(CO)
calculations. With UF4 there were three energetic solutions, with there being two broken symmetry solutions, lowest energy at
ms=1 (triplet) and another at ms=0 (singlet) that is lower than any restricted singlet solutions. The unrestricted UF4
calculation with the triplet initial guess showed ED+D being far faster at convergence. When the ED+D subspace is expanded
(number of matrices in linear expansion), convergence improves. With a 20 matrix subspace, ED+D is the only method converges; it
goes to a restricted solution for Ru species, with Harris guess B3LYP/LANL2DZ basis. They conclude that ED+D can get trapped on
high energy solutions when most other methods fail to converge. For HF/6-31G*, ADIIS and EDIIS results were virtually identical.
For KS-DFT, ED+D was faster for Cd species at B3LYP/3-21G, but ED+D could only be faster than AD+D for Ru at B3LYP/LANL2DZ with
a 20 matrix subspace.



"Optimization of convergence criteria for fragmentation methods"
Z. Sun, T. Zhu, X. Wang, Y. Mei, J.Z.H. Zhang, Chem. Phys. Lett. 2017, 687, 163-170; https://doi.org/10.101
/j.cplett.2017.08.059

The big idea of this paper is that the relaxing of convergence criteria for ab initio calculations of molecule fragments
(implemented with a molecular electron density fragmentation method) can accelerate convergence without sabotaging the accuracy
of the total energy. This work takes a fragment ab initio molecular dynamic (AIMD) approach to study protein dynamics, where a
protein is broken up into fragments embedded in a field of point charge, and atomic forces are calculated with EE-GMFCC
(electrostatically embedded generalized molecular fractionation with conjugate caps) to yield protein potential energy. Included
is the QM interaction energies from interactions of two short-range and non-neighboring fragments. The associated error with
this method, Basis Set Superposition Error (BSSE), they claim is small, and is therefore not calculated or corrected for. EE
GMFCC (MFCC for short) is combined with Gaussian-09 in this paper. They discuss that Gaussian 09 requires for SCF convergence
that the RMS density matrix change be < 10(^-N), where N ranges from 4 to 8. For MFCC, the convergence criterion needs to be
under 10^(-8). Even when convergence criteria is infinitely small, fragmentation estimates will still be inaccurate relative to
the full system counterparts, which gives rise to the so-called fragmentation error. Convergence errors in MFCC, by ignoring
cross correlations of properties in fragments, can be calculated by summing the convergence errors of all fragments and caps.
Since the convergence criterion can give a ball-park estimate of the convergence error (though these quantities are not the
same), and because of how they calculate total MFCC convergence error, a 10(^-N) precision could be accomplished in a small
system with a fragment precision of, say, 10(^-(N-1)). They therefore posit that a loose fragment convergence criterion (from
their data, they put forward 10(^-5)) is precise enough to calculate molecule total energy via a fragmentation method. 

They perform calculations (which are mentioned as being single configuration) are done with M062X/6-31G* basis sets for protein
systems such as 1LE1, 2I9 M, 2OMM, 3DGJ, polypeptide 4R0W, and a 2OMM fragment which respectively contain 218, 246, 107, 106,
96, and 36 atoms. Due to energy differences between fragment and full-system instantiations of the fragment being less than 1
kCal/mol, no potential energy surface calculations are performed. They use elapsed time to measure efficacy of SCF calculations
with varied convergence criteria. The Harris Functional initial guess is used for all calculations, and larger basis sets
6-311G** and B3LYP are also tested. For all systems, the MFCC fragmentation method shows a linear increase in computation time
as the convergence criteria is tightened from N=5 to N=8 (as in 10(^-N)), where a full system ones show a much lower
correlation, and this fluctuation increases as the density matrices become bigger. Specifically, when N decreases from 8 to 5,
the B3LYP calculation becomes about 1.5 times faster, and the M062X becomes about 1.8 times faster.  The speed up is
'functional' dependent (energy functional(?), bigger systems have greater decreases in time needed for convergence as the
convergence criterion is relaxed) and doesn't seem to depend on basis set used. The fragmentation errors in MFCC estimates lead
to relatively small SCF iterational errors, and these errors are smaller than the fragmentation errors themselves. Energy
gradients for the loosened convergence criteria relative to the N=8 criterion of MFCC and full system are mostly comparable, as
are the size of the gradients for each N. Gradients increase as system size increases, and as N decreases, for both methods.
Notably, MFCC still has small convergence gradients, even when N <= 5, for 2I9 M. The full system calculation has far larger
gradients at these relaxed criteria. 