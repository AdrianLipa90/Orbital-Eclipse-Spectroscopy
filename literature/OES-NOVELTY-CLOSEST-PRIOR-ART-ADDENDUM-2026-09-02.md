# OES Novelty Audit — Closest Prior-Art Addendum

Date: 2026-09-02
Scope: deep-search addendum focused specifically on the remaining OES response-bath / natural-orbital / projector-gauge novelty candidates.
Status: literature provenance and claim narrowing only.

## A. PNO-NEVPT2: multireference perturbation plus external natural-orbital compression

Yang Guo, Kantharuban Sivalingam, Edward F. Valeev, Frank Neese,
"SparseMaps—A systematic infrastructure for reduced-scaling electronic structure methods. III. Linear-scaling multireference domain-based pair natural orbital N-electron valence perturbation theory," Journal of Chemical Physics 144, 094111 (2016).
DOI: 10.1063/1.4942769
Primary URL: https://doi.org/10.1063/1.4942769
Supporting record: https://pubmed.ncbi.nlm.nih.gov/26957161/

Relevance:
- combines NEVPT2 with pair-natural-orbital compression of external spaces,
- is substantially closer to OES than generic DMET/BNO prior art because it joins a multireference first-order / perturbative construction with natural-orbital reduction of the external orbital domain,
- strongly constrains any broad claim that OES uniquely combines multireference Hamiltonian response with external natural-orbital compression.

Audit tag: VERY_CLOSE_PRIOR_ART__PNO_NEVPT2

## B. Local PNO-NEVPT2 refinements

"Local N-electron valence state perturbation theory using pair-natural orbitals based on localized virtual molecular orbitals," Journal of Chemical Physics (2023).
DOI: 10.1063/5.0143793
Primary URL: https://doi.org/10.1063/5.0143793
Supporting record: https://pubmed.ncbi.nlm.nih.gov/37094010/

Relevance:
- further develops PNO compression inside NEVPT2 external/semi-internal spaces,
- confirms that natural-orbital compression of multireference perturbative external sectors is an active, established methodology rather than an isolated precedent.

Audit tag: VERY_CLOSE_PRIOR_ART__LOCAL_PNO_NEVPT2

## C. PNO-CASPT2 comparator

Felipe Menezes, Daniel Kats, Hans-Joachim Werner,
"Local complete active space second-order perturbation theory using pair natural orbitals (PNO-CASPT2)," Journal of Chemical Physics 145, 124115 (2016).
DOI: 10.1063/1.4963019
Primary URL: https://doi.org/10.1063/1.4963019

Relevance:
- independent multireference perturbation framework combining CAS reference information and PNO reduction,
- reinforces that `multireference active state -> compressed external natural-orbital domain` is known prior art.

Audit tag: VERY_CLOSE_PRIOR_ART__PNO_CASPT2

## D. NEVPT2 first-order densities and natural orbitals

Current ORCA NEVPT2 documentation records FIC-NEVPT2 unrelaxed densities for a state

`gamma(p,q) = <Psi_I|E^p_q|Psi_I>`

with `Psi_I = |0> + |1>`, including a first-order density option and generation of natural orbitals from the resulting density. The implementation also includes state selection, state-average canonicalization options, quasi-degenerate SC-NEVPT2 and revised natural orbitals.

Documentation URL:
https://orca-manual.mpi-muelheim.mpg.de/contents/modelchemistries/NEVPT2.html

Relevance:
- strongly constrains novelty claims based merely on `first-order/response wavefunction -> 1-RDM -> natural orbitals`,
- demonstrates that response/perturbative densities and natural-orbital generation are already operational in mainstream multireference quantum-chemistry software.

Audit tag: KNOWN_IMPLEMENTATION__NEVPT2_RESPONSE_DENSITY_TO_NATURAL_ORBITALS

## E. State-averaged excited-state PNO manifold

Chong Peng, Marjory C. Clement, Edward F. Valeev,
"State-Averaged Pair Natural Orbitals for Excited States: A Route toward Efficient Equation of Motion Coupled-Cluster," Journal of Chemical Theory and Computation 14, 5597-5607 (2018).
DOI: 10.1021/acs.jctc.8b00171
Primary URL: https://doi.org/10.1021/acs.jctc.8b00171
Supporting record: https://pubmed.ncbi.nlm.nih.gov/30252467/

Relevance:
- state-averaged pair densities are used to construct natural orbitals for a target manifold of excited states,
- method is explicitly motivated by robustness to state degeneracies and root flipping,
- strongly constrains any broad OES claim based on using complete excited-state manifolds to build a shared natural-orbital representation.

Audit tag: VERY_CLOSE_PRIOR_ART__STATE_AVERAGED_MANIFOLD_PNO

## F. Projected-Hamiltonian response primitive

Celestino Angeli et al.,
"Introduction of n-electron valence states for multireference perturbation theory," Journal of Chemical Physics 114, 10252-10264 (2001).
DOI: 10.1063/1.1361246

Celestino Angeli, Renzo Cimiraglia, Jean-Paul Malrieu,
"N-electron valence state perturbation theory: a fast implementation of the strongly contracted variant," Chemical Physics Letters 350, 297-305 (2001).
DOI: 10.1016/S0009-2614(01)01303-3

Relevant SC-NEVPT2 schematic:

`|Psi_l^(k)> = P_l^(k) H |Psi0>`

Relevance:
- the primitive `project Hamiltonian action into an external class` is established prior art,
- OES `|chi_c> = Q H |Psi_c^P>` is therefore not independently novel.

Audit tag: KNOWN_PRIOR_ART__PROJECTED_H_RESPONSE

## Updated narrow novelty assessment

After this addendum, the following sequence is NOT defensibly novel at a broad level:

`multireference reference -> H-generated external response -> response/perturbative density -> natural orbitals -> compressed external space`.

Literature contains all of these ingredients, including combinations that are materially close (especially PNO-NEVPT2, PNO-CASPT2 and state-averaged PNO-EOM-CC).

The remaining OES-specific differentiators are therefore narrower:

### 1. Exact class construction
OES uses physically complete classes spanning ground singlet, full triplet manifold, dark singlet and full bright p manifold rather than conventional excitation-pair classes or merely a set of numerical roots.

Status: POTENTIALLY_DISTINCT, but not yet established as novel.

### 2. Class normalization / equalization rule
OES trace-normalizes each complete-class external response density before equal class weighting, preventing a class with larger raw coupling norm from automatically dominating bath selection.

Status: POTENTIALLY_DISTINCT. Dedicated search required for identical weighting rules in state-averaged PNO / embedding literature.

### 3. Complete-degenerate-eigenspace admission under a hard fixed-width budget
OES treats numerical occupation degeneracies as indivisible blocks when constructing the bath / active representation rather than truncating inside the degenerate eigenspace.

Status: MATHEMATICALLY_NATURAL / POSSIBLY_DISTINCT_AS_PROTOCOL, not yet a novelty claim.

### 4. Explicit randomized Q-basis projector recovery test
OES deliberately rotates the complete external orbital basis, reconstructs transformed integrals, rebuilds the selected bath and checks principal cosines / projector distance / retained dimension.

Status: POTENTIALLY_DISTINCT_VALIDATION_GATE. Orbital invariance itself is known; explicit end-to-end randomized gauge testing may be a reproducibility contribution rather than new theory.

### 5. Fixed-20Q blind model-chemistry validation ladder
OES applies one fixed 20-spin-orbital interface across atoms and molecules, with varying active-electron/frozen-core partitions, predeclared blind spectroscopy gates, manifold/intensity/property validation and append-only receipts.

Status: currently the strongest candidate for a publishable FRAMEWORK / BENCHMARK contribution rather than a claim of a wholly new electronic-structure primitive.

## Current strongest defensible framing

The evidence now favors framing OES as a novel or potentially novel `composition + validation architecture` only after further search, rather than claiming a novel primitive response-bath algorithm.

Recommended manuscript-level claim wording at this stage:

`We introduce a fixed-width, symmetry-complete spectroscopy workflow that combines established active-space, response-space and natural-orbital ideas in a specific class-balanced, gauge-audited construction, and evaluate it under a common blind validation contract across chemically distinct systems.`

Do not use `first`, `unprecedented`, or `new response-bath theory` without a further patent-style/citation-network audit.