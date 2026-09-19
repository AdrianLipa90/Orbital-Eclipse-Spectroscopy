# OES Novelty Audit — Patent Prior Art

Date: 2026-09-02
Scope: preliminary patent-style search for claims overlapping active-space quantum chemistry, fixed/reduced quantum resources, RDM-based hybrid reconstruction, downfolding and spectroscopy/property recovery.
Status: preliminary landscape only; not a legal freedom-to-operate or patentability opinion.

## 1. Hybrid quantum-classical simulation of chemical systems

US20230377693A1 — "Hybrid quantum-classical computing simulation of chemical systems"
Google Patents: https://patents.google.com/patent/US20230377693A1/en
Related international family: WO2023227470A1
Google Patents: https://patents.google.com/patent/WO2023227470A1/en

Relevant disclosed concepts include:
- total and active-space electronic Hamiltonians;
- translation of the active-space Hamiltonian to a qubit basis;
- active-space wavefunctions and RDMs represented/measured on a quantum component;
- classical use of active-space RDM information to approximate additional RDMs / characteristics of the total system;
- structural, interaction and response characteristics as downstream outputs;
- inactive orbitals including core and virtual orbitals;
- NEVPT2 / related classical correction pathways in family disclosures.

OES relevance:
- strong patent prior art against a broad system claim of `active-space quantum core + classical reconstruction of full-system properties/response`;
- OES differentiation must be narrower than generic hybrid active-space/RDM architecture.

Audit tag: PATENT_CLOSE_PRIOR_ART__HYBRID_ACTIVE_SPACE_RDM_PROPERTIES

## 2. Increasing representation accuracy without extra quantum resources

CA3253674A1 — "Increasing representation accuracy of quantum simulations without additional quantum resources"
Google Patents: https://patents.google.com/patent/CA3253674A1/en

The disclosure explicitly frames quantum chemistry as solving the hard part inside a reduced active space to lower qubit requirements while recovering effects associated with a larger representation without simply adding qubits.

OES relevance:
- directly constrains generic claims around recovering external-space accuracy at a fixed quantum-resource budget;
- reinforces the literature prior art from Takeshita et al. (PRX 2020).

Audit tag: PATENT_CLOSE_PRIOR_ART__FIXED_RESOURCE_REPRESENTATION_RECOVERY

## 3. Downfolded electronic Hamiltonians in hybrid quantum-classical architectures

US20240202561A1 — "Method and system for downfolding electronic hamiltonians using a hybrid quantum-classical architecture"
Google Patents: https://patents.google.com/patent/US20240202561A1/en

The disclosure explicitly describes an effective Hamiltonian as a lower-dimensional representation of a parent electronic Hamiltonian and references Feshbach-Löwdin-Fano, Schrieffer-Wolff, perturbative and other downfolding traditions.

OES relevance:
- strong patent prior art against broad claims on `large Hamiltonian -> reduced active Hamiltonian for quantum execution`;
- OES Feshbach/Schur machinery must be treated as known machinery, with any differentiation sought in P/Q construction, validation contract, symmetry/class balancing or end-to-end workflow.

Audit tag: PATENT_CLOSE_PRIOR_ART__HYBRID_HAMILTONIAN_DOWNFOLDING

## Preliminary patent-landscape verdict

Broad OES claims that should NOT be treated as proprietary novelty without substantially narrower language:
- active-space Hamiltonian mapped to qubits;
- fixed/reduced qubit resources for chemistry;
- quantum active-space RDM measurement followed by classical property reconstruction;
- classical external-correlation correction of active quantum calculations;
- effective-Hamiltonian/downfolding reduction for quantum execution;
- generic structural/response property computation from the hybrid calculation.

Potentially differentiable OES areas remain narrow implementation/protocol combinations:
1. complete physical state-class balancing of Q-response densities;
2. trace-normalize-per-class then equal-class weighting before external NO bath construction;
3. complete-degenerate-block admission under one fixed-width model-chemistry contract;
4. explicit randomized Q-basis projector recovery as a required reproducibility gate;
5. one blind, predeclared, append-only validation ladder spanning atomic and molecular spectroscopy under a fixed 20-spin-orbital interface.

These remain SEARCH TARGETS, not promoted novelty claims.

No `first`, `unprecedented`, patentability, freedom-to-operate or priority conclusion is established by this file.