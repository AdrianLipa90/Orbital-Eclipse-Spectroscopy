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

## OES-specific novelty candidates retained after this literature pass

The literature above strongly covers:
- frozen core / CAS,
- automated active-space selection,
- fixed-resource quantum active spaces,
- sparse / determinant FCI,
- bath embedding,
- bath natural orbitals,
- state-averaged embedding,
- effective-Hamiltonian downfolding,
- active-space vibrational spectroscopy.

The strongest OES-specific candidates requiring deeper dedicated search are therefore:

1. `COMPLETE_CLASS_HAMILTONIAN_RESPONSE_BATH`
   - response states `|chi_c> = Q H |Psi_c^P>` built for complete physical state classes,
   - class-wise normalized external 1-RDMs,
   - complete degenerate natural-orbital eigenspaces.

2. `EXPLICIT_Q_BASIS_ROTATION_GAUGE`
   - active verification that arbitrary external-orbital rotations recover the same selected bath projector / principal subspace.

3. `SYMMETRY_COMPLETE_FIXED_QUBIT_MODEL_CHEMISTRY`
   - one fixed 20-spin-orbital register across systems with different electron counts and frozen-core partitions,
   - complete symmetry / degeneracy blocks retained under the fixed width.

4. `FLAVOR_CORE_PLUS_COLOR_DRESSING_DIAGNOSTIC`
   - empirical decomposition where the fixed P-space retains relative manifold geometry / transition structure and external Q primarily repairs the common absolute-energy offset.

5. `END_TO_END_FIXED_20Q_BLIND_SPECTROSCOPY_PROTOCOL`
   - common register and solver contract across atoms and molecules,
   - active-space selection independent of experimental benchmark outputs,
   - common validation of energies, manifolds, transition strength, PES, geometry, vibration, rotation, density-derived observables and external-bath recovery.

These are marked POTENTIALLY_NOVEL_PENDING_DEEPER_PRIOR_ART_SEARCH. Absence from this source ledger is not proof of novelty.
