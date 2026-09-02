# OES Novelty / Prior-Art Source Ledger

Date: 2026-09-02
Scope: source ledger for the claim-by-claim novelty audit of Orbital-Eclipse-Spectroscopy (OES).
Status semantics: literature provenance only. This file does not itself promote any OES claim to established novelty.

## 1. Automated active-space selection

Christopher J. Stein, Markus Reiher, "Automated Selection of Active Orbital Spaces," Journal of Chemical Theory and Computation 12(4), 1760-1771 (2016).
DOI: 10.1021/acs.jctc.6b00156
Primary URL: https://doi.org/10.1021/acs.jctc.6b00156
Supporting public record: https://pubmed.ncbi.nlm.nih.gov/26959891/

Relevance to OES:
- prior art for automated / systematic active-orbital selection,
- demonstrates that poor orbital choice can qualitatively damage multiconfigurational calculations,
- constrains novelty claims around OES orbital-selection logic.

Audit tag: KNOWN_PRIOR_ART__ACTIVE_SPACE_SELECTION

## 2. Fixed / reduced quantum resources for active-space chemistry

Tyler Takeshita, Nicholas C. Rubin, Zhang Jiang, Eunseok Lee, Ryan Babbush, Jarrod R. McClean, "Increasing the Representation Accuracy of Quantum Simulations of Chemistry without Extra Quantum Resources," Physical Review X 10, 011004 (2020).
DOI: 10.1103/PhysRevX.10.011004
Primary URL: https://doi.org/10.1103/PhysRevX.10.011004
Google Research record: https://research.google/pubs/pub47848
Preprint: https://arxiv.org/abs/1902.10679

Relevance to OES:
- explicit prior art for improving active-space quantum-chemistry accuracy without increasing quantum resources,
- directly constrains generic claims of novelty based only on a fixed qubit budget.

Audit tag: KNOWN_PRIOR_ART__FIXED_QUANTUM_RESOURCE_ACTIVE_SPACE

## 3. Quantum-computing spectroscopy with small active spaces

Shih-Kai Chou, Jyh-Pin Chou, Alice Hu, Yuan-Chung Cheng, Hsi-Sheng Goan, "Accurate harmonic vibrational frequencies for diatomic molecules via quantum computing," Physical Review Research 5, 043216 (2023).
DOI: 10.1103/PhysRevResearch.5.043216
Primary URL: https://doi.org/10.1103/PhysRevResearch.5.043216
Journal URL: https://journals.aps.org/prresearch/abstract/10.1103/PhysRevResearch.5.043216

Relevance to OES:
- very close prior art for active-space quantum chemistry -> PES -> harmonic vibrational spectroscopy,
- constrains novelty claims around using small qubit registers to recover diatomic vibrational observables.

Audit tag: KNOWN_PRIOR_ART__QUANTUM_ACTIVE_SPACE_SPECTROSCOPY

## 4. Density-matrix embedding and algebraic bath construction

Gerald Knizia, Garnet Kin-Lic Chan, "Density Matrix Embedding: A Simple Alternative to Dynamical Mean-Field Theory," Physical Review Letters 109, 186404 (2012).
DOI: 10.1103/PhysRevLett.109.186404
Primary URL: https://doi.org/10.1103/PhysRevLett.109.186404
Journal URL: https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.109.186404

Related chemical-Hamiltonian extension:
Gerald Knizia, Garnet Kin-Lic Chan, "Density matrix embedding: A strong-coupling quantum embedding theory," arXiv:1212.2679 (2012).
URL: https://arxiv.org/abs/1212.2679

Relevance to OES:
- prior art for mapping a large interacting system to a smaller embedded subsystem plus algebraically constructible bath,
- constrains generic novelty claims around 'active core + external bath'.

Audit tag: KNOWN_PRIOR_ART__EMBEDDING_BATH

## 5. Bath natural orbitals / systematically improvable embedding

Max Nusspickel, George H. Booth, "Systematic Improvability in Quantum Embedding for Real Materials," Physical Review X 12, 011046 (2022).
DOI: 10.1103/PhysRevX.12.011046
Primary URL: https://doi.org/10.1103/PhysRevX.12.011046
Journal URL: https://journals.aps.org/prx/abstract/10.1103/PhysRevX.12.011046
Preprint: https://arxiv.org/abs/2107.04916

Additional expectation-value / BNO implementation context:
"Effective Reconstruction of Expectation Values from Ab Initio Quantum Embedding," Journal of Chemical Theory and Computation (2023).
DOI: 10.1021/acs.jctc.2c01063
URL: https://doi.org/10.1021/acs.jctc.2c01063

Relevance to OES:
- close prior art for systematic bath enlargement through correlated bath natural orbitals,
- close comparator for OES external natural-orbital bath,
- the OES novelty question therefore concerns the exact response-state construction, complete-state classes, complete degenerate eigenspaces and explicit Q-basis rotation gauge, not the general use of bath natural orbitals.

Audit tag: CLOSE_PRIOR_ART__BATH_NATURAL_ORBITALS

## 6. State-averaged embedding for excited states

Zhe-Bin Guan, Hong Jiang, "State-Averaged Density Matrix Embedding Theory for Local Excitations," arXiv:2607.08178 (2026).
URL: https://arxiv.org/abs/2607.08178

Relevance to OES:
- very recent close prior art addressing ground-state bias in embedding baths by using state-averaged information,
- important comparator to the OES complete-class response construction for ground / triplet / dark / bright classes,
- as of this ledger date this is a preprint and must be identified as such.

Audit tag: CLOSE_PRIOR_ART__STATE_AVERAGED_EXCITED_STATE_EMBEDDING

## 7. Active-space downfolding / effective Hamiltonians

Nicholas P. Bauman, Eric J. Bylaska, Sriram Krishnamoorthy, Guang Hao Low, Nathan Wiebe, et al., "Downfolding of many-body Hamiltonians using active-space models: Extension of the sub-system embedding sub-algebras approach to unitary coupled cluster formalisms," Journal of Chemical Physics 151, 014107 (2019).
DOI: 10.1063/1.5094643
Primary URL: https://doi.org/10.1063/1.5094643
Preprint: https://arxiv.org/abs/1902.01553

Nicholas P. Bauman, Guang Hao Low, Karol Kowalski, "Quantum simulations of excited states with active-space downfolded Hamiltonians," Journal of Chemical Physics 151, 234114 (2019).
DOI: 10.1063/1.5128103
Primary URL: https://doi.org/10.1063/1.5128103
Preprint: https://arxiv.org/abs/1909.06404

Relevance to OES:
- strong prior art for integrating external correlation into a reduced active-space effective Hamiltonian,
- directly constrains generic novelty claims around OES Feshbach/Schur/downfolded Hamiltonians,
- OES-specific novelty, if established, must lie in how P and Q are constructed / gauged / validated and in the complete protocol rather than in downfolding itself.

Audit tag: KNOWN_PRIOR_ART__ACTIVE_SPACE_DOWNFOLDING

## 8. Full configuration interaction as exact active-space reference

Jeppe Olsen, Poul Jorgensen, Henrik Koch, Anna Balkova, Rodney J. Bartlett, "Full Configuration-Interaction and state of the art correlation calculations on water in a valence double-zeta basis with polarization functions," Journal of Chemical Physics 104, 8007 (1996).
DOI: 10.1063/1.471518
Primary URL: https://doi.org/10.1063/1.471518
Institutional record: https://pure.au.dk/portal/en/publications/full-configuration--interaction-and-state-of-the-art-correlation-calculations-on-water-in-a-valence-doublezeta-basis-with-polarization-functions

Relevance to OES:
- historical reference confirming determinant-FCI as established methodology,
- constrains novelty claims around exact determinant-sector diagonalization and Slater-Condon construction.

Audit tag: KNOWN_PRIOR_ART__FCI

## 9. Hamiltonian-generated external response / first-order interacting space

Celestino Angeli, Renzo Cimiraglia, Stefano Evangelisti, Thierry Leininger, Jean-Paul Malrieu, "Introduction of n-electron valence states for multireference perturbation theory," Journal of Chemical Physics 114, 10252-10264 (2001).
DOI: 10.1063/1.1361246
Primary URL: https://doi.org/10.1063/1.1361246

Celestino Angeli, Renzo Cimiraglia, Jean-Paul Malrieu, "N-electron valence state perturbation theory: a fast implementation of the strongly contracted variant," Chemical Physics Letters 350, 297-305 (2001).
DOI: 10.1016/S0009-2614(01)01303-3
Primary URL: https://doi.org/10.1016/S0009-2614(01)01303-3

Celestino Angeli, Renzo Cimiraglia, Jean-Paul Malrieu, "N-electron valence state perturbation theory: A spinless formulation and an efficient implementation of the strongly contracted and of the partially contracted variants," Journal of Chemical Physics 117, 9138-9153 (2002).
DOI: 10.1063/1.1515317
Primary URL: https://doi.org/10.1063/1.1515317

Relevance to OES:
- decisive prior art for the primitive `projected Hamiltonian response` itself,
- strongly contracted NEVPT2 constructs perturber functions by projecting the Hamiltonian acting on the multireference state into external excitation subspaces, schematically `|Psi_l^(k)> = P_l^(k) H |Psi0>`,
- therefore the OES primitive `|chi_c> = Q H |Psi_c^P>` MUST NOT be claimed as novel by itself,
- NEVPT2 also establishes important orbital-rotation invariance properties in active/inactive/virtual subspaces depending on contraction scheme.

Audit tag: KNOWN_PRIOR_ART__PROJECTED_HAMILTONIAN_RESPONSE

## 10. State-averaged natural-orbital spaces for complete excited-state manifolds

Chong Peng, Marjory C. Clement, Edward F. Valeev, "State-Averaged Pair Natural Orbitals for Excited States: A Route toward Efficient Equation of Motion Coupled-Cluster," Journal of Chemical Theory and Computation 14(11), 5597-5607 (2018).
DOI: 10.1021/acs.jctc.8b00171
Primary URL: https://doi.org/10.1021/acs.jctc.8b00171
Supporting record: https://pubmed.ncbi.nlm.nih.gov/30252467/

Relevance to OES:
- strong prior art for constructing natural-orbital truncation spaces from densities averaged over a target excited-state manifold,
- explicitly motivated by robustness to root flipping and state degeneracies,
- directly constrains novelty claims based solely on `average densities over a complete manifold -> natural orbitals`,
- OES remains structurally different because its densities are built in Q from Hamiltonian coupling-response states and are then used as a bath coupled back to a fixed P-space.

Audit tag: CLOSE_PRIOR_ART__STATE_AVERAGED_MANIFOLD_NATURAL_ORBITALS

## 11. Unitary-rotation invariance as an established multistate design requirement

Toru Shiozaki, Werner Gyorffy, Paolo Celani, Hans-Joachim Werner, "Communication: Extended multi-state complete active space second-order perturbation theory: energy and nuclear gradients," Journal of Chemical Physics 135, 081106 (2011).
DOI: 10.1063/1.3633329
Primary URL: https://doi.org/10.1063/1.3633329
Supporting record: https://pubmed.ncbi.nlm.nih.gov/21895152/

Relevant earlier multistate formulation:
Alexander A. Granovsky, "Extended multi-configuration quasi-degenerate perturbation theory: The new approach to multi-state multi-reference perturbation theory," Journal of Chemical Physics 134, 214113 (2011).
DOI: 10.1063/1.3596699
Primary URL: https://doi.org/10.1063/1.3596699

Relevance to OES:
- strong prior art showing that invariance under unitary rotations of a model/reference subspace is a known correctness requirement, not an OES invention,
- complements the orbital-rotation invariance properties of NEVPT2,
- therefore an abstract claim of `rotation invariance / covariance` is not novel,
- OES's explicit randomized Q-basis rotation gauge remains potentially useful as a reproducible validation protocol for the selected bath projector, but should be framed as a validation implementation unless deeper search finds a genuinely new mathematical property.

Audit tag: KNOWN_PRIOR_ART__UNITARY_ROTATION_INVARIANCE_REQUIREMENT

## 12. Fixed register across multiple chemical compositions

Panagiotis Kl. Barkoutsos, Fotios Gkritsis, Pauline J. Ollitrault, Igor O. Sokolov, Stefan Woerner, Ivano Tavernelli, "Quantum algorithm for alchemical optimization in material design," Chemical Science 12, 4345-4352 (2021).
DOI: 10.1039/D0SC05718E
Primary URL: https://doi.org/10.1039/D0SC05718E
Open record: https://pmc.ncbi.nlm.nih.gov/articles/PMC8179438/
Preprint: https://arxiv.org/abs/2008.06449

Relevance to OES:
- important prior art for deliberately keeping the number of qubits / Hamiltonian matrix-element structure constant across a family of different molecular compositions by extending active spaces to a common size,
- directly weakens any generic claim that `same qubit count across different molecules` is novel,
- OES differs in purpose and mechanism: fixed 20-spin-orbital spectroscopy/model-chemistry ladder, varying frozen-core electron count, symmetry-complete active blocks, blind benchmark separation and common property/validation contract rather than alchemical superposition/optimization.

Audit tag: CLOSE_PRIOR_ART__FIXED_REGISTER_ACROSS_CHEMICAL_SPACE

## 13. Additional fixed-resource / optimized-orbital quantum chemistry comparators

Luca A. A. et al., "Complete Active Space Methods for NISQ Devices: The Importance of Canonical Orbital Optimization for Accuracy and Noise Resilience," Journal of Chemical Theory and Computation 19, 2863-2872 (2023).
DOI: 10.1021/acs.jctc.3c00123
Primary URL: https://doi.org/10.1021/acs.jctc.3c00123

"Improving the Accuracy of Variational Quantum Eigensolvers with Fewer Qubits Using Orbital Optimization" (2023).
Supporting record: https://pubmed.ncbi.nlm.nih.gov/36696487/

Relevance to OES:
- additional prior art for optimizing orbital representations under strict qubit constraints,
- reinforces that fixed-qubit accuracy improvement is an established research direction.

Audit tag: KNOWN_PRIOR_ART__FIXED_QUBIT_ORBITAL_OPTIMIZATION

## Updated OES novelty verdict after deeper prior-art pass

The following primitives are now explicitly DOWNGRADED from possible novelty:

- `Q H |Psi>` / projected Hamiltonian response by itself -> KNOWN (SC-NEVPT2 / first-order interacting space),
- state-averaged natural orbitals over a target excited-state manifold -> KNOWN / CLOSE PRIOR ART (state-averaged PNO),
- unitary rotation invariance as a correctness principle -> KNOWN (NEVPT2, XMS-CASPT2 and related multistate methods),
- fixed qubit/register size across multiple molecules -> CLOSE PRIOR ART (including alchemical fixed-width constructions),
- bath natural orbitals, embedding and effective-Hamiltonian downfolding -> KNOWN.

The OES candidates that remain potentially differentiable are therefore narrower combinations / protocols:

1. `COMPLETE_CLASS_Q_RESPONSE_TO_BATH_PIPELINE`
   - start from complete physical classes (ground singlet, complete triplet manifold, dark singlet, complete bright manifold),
   - form external Hamiltonian responses in Q,
   - construct spin-summed external 1-RDM for each complete class,
   - trace-normalize class densities and combine them without using benchmark energies,
   - diagonalize to external natural orbitals,
   - admit numerically degenerate occupation eigenspaces only as complete blocks,
   - couple that bath back to the fixed symmetry-complete P core.
   Status: POTENTIALLY_NOVEL_AS_EXACT_COMBINATION; individual ingredients have prior art.

2. `EXPLICIT_RANDOM_Q_BASIS_PROJECTOR_GAUGE`
   - deliberate arbitrary orthogonal rotation of the complete external orbital basis,
   - rebuild the response bath from transformed integrals,
   - compare selected subspaces by principal cosines / projector distance and require recovery of the same bath dimension and projector.
   Status: POTENTIALLY_DISTINCT_VALIDATION_PROTOCOL; rotation invariance itself is KNOWN.

3. `SYMMETRY_COMPLETE_FIXED_20Q_MODEL_CHEMISTRY`
   - exactly one 20-spin-orbital interface across atoms and molecules,
   - variable frozen-core partition and active electron count,
   - fixed-N / fixed-Ms exact reference sectors,
   - symmetry/degeneracy blocks are not split merely to satisfy width,
   - same implementation and predeclared validation contract across systems.
   Status: POTENTIALLY_NOVEL_FRAMEWORK; fixed-width chemistry alone has prior art.

4. `FLAVOR_CORE_PLUS_COLOR_DRESSING_DIAGNOSTIC`
   - empirical observation that selected P-space preserves relative manifold structure / transition strength while Q-space mainly repairs a common absolute energy offset in the tested He construction.
   Status: POTENTIALLY_NOVEL_DIAGNOSTIC_FORMULATION; requires broader-system validation before any general claim.

5. `END_TO_END_FIXED_20Q_BLIND_SPECTROSCOPY_PROTOCOL`
   - experimental benchmark values excluded from active-space/bath selection,
   - one fixed-width computational contract across atoms and molecules,
   - common outputs spanning energies, complete manifolds, oscillator strength, PES, geometry, vibration, rotation and density-derived observables,
   - explicit PASS/FAIL receipts and regression gates.
   Status: POTENTIALLY_NOVEL_BENCHMARK_PROTOCOL / MODEL-CHEMISTRY PACKAGE; individual computations are known.

Current strongest defensible novelty target after this pass:

`THE EXACT OES COMPOSITION AND VALIDATION CONTRACT`, not any single primitive.

No priority claim is promoted by this ledger. A dedicated patent-style and citation-network search remains required before using FIRST / NOVEL / UNPRECEDENTED language in a manuscript.