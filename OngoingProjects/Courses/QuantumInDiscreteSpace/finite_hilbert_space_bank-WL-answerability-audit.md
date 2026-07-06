# WL-Answerability Audit of `finite_hilbert_space_question_bank.md`

Source files fully read: `finite_hilbert_space_question_bank.md` (185 questions, sections A-L),
`PIPELINE.md` (the "learning-by-computing" DNA), `Question-List.md` (the 213-item computational
curriculum), `Question-List-Audit.md` (the prior coverage audit).
Audit date: 2026-07-05.
Question asked: **what portion of the 185 questions can be *fully answered in the Wolfram Language*
(compute/verify style), as opposed to answered only in prose?**

## Why this needs a definition first

The curriculum's governing principle (PIPELINE 1) is that *the reader computes and the rendered
object is the answer*: "could the reader have written this sentence before evaluating the cell? If
yes, it is a spoiler." A question is therefore "fully WL-answerable" only if its complete answer is
an object, value, or boolean the reader *obtains by running a cell*, not a concept the cell merely
decorates. Without this test the question is empty, because almost any statement can be *illustrated*
by choosing a numeric instance. The tiers below are that test made sharp.

### Rubric

- **Tier A: Fully WL-answerable (compute / verify).** The question asks to construct, compute,
  decide, or verify a definite object. The rendered WL object (state, operator, matrix, distribution,
  eigen-set) or the returned `True`/`False` *is* the complete answer. Toggling the operation under
  study changes the output (PIPELINE 4). Verbs: build, compute, construct, find, test, verify,
  decide, diagonalize, decompose, simulate.
- **Tier B: Partial / demonstrable.** The real content is a theorem, a "why", or a conceptual
  distinction, but a computation on a *discriminating instance* faithfully witnesses it (the DNA's
  own "drop to a concrete instance where a closed symbolic form is intractable" move, used for Fermi
  golden rule, Berry phase, Solovay-Kitaev). The WL component is genuine, but the cell does not
  capture the full generality or the proof: the reader still needs a sentence the cell did not write.
- **Tier C: Conceptual / not WL-answerable.** The answer is an explanation, definition,
  interpretation, or physical framing. No computation constitutes it; WL can at most gesture at an
  example. Verbs: what is, why (in the axiomatic sense), what does it mean, what is the distinction.

The A/B line is the only soft one: a handful of "construct and then interpret" items could be read
either way. The headline is reported as a band for that reason, and every boundary call is shown
per-question below so the reader can move it.

## Executive verdict

| Tier | Meaning | Count | Share |
|------|---------|-------|-------|
| **A** | Fully WL-answerable (the cell is the answer) | **93** | **50%** |
| **B** | Demonstrable on an instance; answer is partly conceptual | **69** | **37%** |
| **C** | Purely conceptual; no cell is the answer | **23** | **12%** |
| A + B | Has a genuine WL component | **162** | **88%** |

**Roughly half of the bank (about 50%, 93 questions) can be turned into full compute/verify items in
the curriculum's own idiom.** Another 37% (69) can be *demonstrated* by a well-chosen computation but
the answer stays partly a sentence, so they belong in the manual only if re-cast to a concrete
instance with the conceptual half stated in prose (exactly the treatment Part 7 already gives Fermi's
golden rule and Berry phase). The remaining 12% (23) are review-exam / oral-exam questions with no
computational answer; they do not belong in a learning-by-computing manual and should stay in the
companion bank.

Given the A/B boundary noise (about a dozen items are one reviewer's judgment from flipping, and a
skeptical re-classification of a stratified sample moved a net two items toward B), read the
fully-answerable share as **48%-56%**, not a hard 50%.

The distribution is not uniform. The computational core of the subject (density operators, composite
systems, gates and algorithms, spin operators) is heavily Tier A; the foundational and interpretive
sections (A "bases and dimension", B "rays and phases", L "foundations") are where the C items
concentrate.

## Per-section breakdown

| Sec | Theme | Q's | A | B | C | A-share |
|-----|-------|-----|---|---|---|---------|
| A | Finite Hilbert spaces, bases, dimension | 12 | 5 | 4 | 3 | 42% |
| B | States, rays, superposition, phases | 12 | 3 | 6 | 3 | 25% |
| C | Operators, observables, spectra, commutators | 16 | 5 | 10 | 1 | 31% |
| D | Measurement, probability, state update | 16 | 9 | 6 | 1 | 56% |
| E | Unitary dynamics, two-level evolution | 14 | 7 | 7 | 0 | 50% |
| F | Spin-1/2, Pauli, qubits | 16 | 9 | 3 | 4 | 56% |
| G | Angular momentum and rotations | 18 | 8 | 8 | 2 | 44% |
| H | Composite systems, tensor product, entanglement | 20 | 14 | 5 | 1 | 70% |
| I | Density operators, mixed states, measurements | 14 | 11 | 3 | 0 | 79% |
| J | Degeneracy, perturbations, effective subspaces | 12 | 6 | 5 | 1 | 50% |
| K | Quantum information, gates, protocols, algorithms | 25 | 15 | 8 | 2 | 60% |
| L | Conceptual and foundational | 10 | 1 | 4 | 5 | 10% |
| **Total** | | **185** | **93** | **69** | **23** | **50%** |

## Full classification (all 185)

Reason column: what makes it that tier. "witnessed" = a computation exhibits the fact but the answer
is partly a sentence (Tier B). "def./interp." = definitional or interpretive, no computing answer
(Tier C).

### A. Finite Hilbert spaces, bases, dimension
| Q | Topic | Tier | Reason |
|---|-------|------|--------|
| 1 | vector space vs Hilbert space | C | def./interp.: inner product + completeness |
| 2 | conditions for an orthonormal basis | B | can test a given set (Gram = I); "conditions" is conceptual |
| 3 | expansion coefficients -> column vector | A | $c_i=\langle i|\psi\rangle$, compute |
| 4 | coordinate change under basis change | A | apply transition matrix |
| 5 | active unitary vs passive basis change | B | both computable; distinction is conceptual |
| 6 | why $\dim=2s+1$ for spin-$s$ | B | build spin-$s$ ops (dim checks out); "why" is rep theory |
| 7 | completeness as resolution of identity | A | $\sum_i|i\rangle\langle i|=I$, verify |
| 8 | dual vector and inner product | A | conjugate-transpose, compute |
| 9 | subspace vs invariant subspace vs eigenspace | C | def./interp. distinction |
| 10 | project onto a subspace | A | build $P=\sum|b\rangle\langle b|$, apply |
| 11 | basis-independent vs representation-dependent info | B | invariants (trace, spectrum) computable; claim is conceptual |
| 12 | model a two-level atom without wavefunctions | C | conceptual modeling framing |

### B. States, rays, superposition, phases
| Q | Topic | Tier | Reason |
|---|-------|------|--------|
| 13 | global phase gives the same state | B | witnessed: identical probabilities; "why" conceptual |
| 14 | what is a ray, why the right object | C | def./interp. |
| 15 | superposition vs classical mixing | B | witnessed: interference vs diagonal mixture |
| 16 | relative phase affects statistics, global does not | B | compute $P(\varphi)$; paired with a conceptual half |
| 17 | parameterize a pure state up to global phase | A | build $\{\cos\tfrac\theta2,e^{i\varphi}\sin\tfrac\theta2\}$, verify coverage |
| 18 | physical content of the Bloch vector | B | compute $\langle\vec\sigma\rangle$; "what content" interpretive |
| 19 | preparation vs measurement finding a basis state | C | conceptual distinction |
| 20 | orthogonality and why it matters operationally | B | compute $\langle a|b\rangle=0$; operational "why" conceptual |
| 21 | transition probability between pure states | A | $|\langle\phi|\psi\rangle|^2$ |
| 22 | when are two superpositions indistinguishable | B | witnessed: equal density matrices up to phase |
| 23 | pure state vector vs rank-one projector | A | build $|\psi\rangle\langle\psi|$, check rank/idempotence |
| 24 | specify a state without a privileged basis | C | conceptual (basis-free description) |

### C. Operators, observables, spectra, commutators
| Q | Topic | Tier | Reason |
|---|-------|------|--------|
| 25 | conditions for a valid observable | B | test Hermiticity; "conditions" conceptual |
| 26 | why Hermitian eigenvalues are real | B | eigenvalues come out real; proof is conceptual |
| 27 | why distinct-eigenvalue eigenvectors are orthogonal | B | witnessed; proof conceptual |
| 28 | spectral decomposition and measurement outcomes | B | build decomposition; "encodes outcomes" interpretive |
| 29 | what changes under degeneracy | B | build degenerate projectors; "what changes" conceptual |
| 30 | degenerate-eigenspace projector vs single-vector | A | build both, compare rank |
| 31 | physical meaning of matrix elements | C | interpretive |
| 32 | reconstruct operator from its basis action | A | $A=\sum|i\rangle\langle i|A$ |
| 33 | compatibility (commuting) and joint eigenbasis | B | test $[A,B]=0$, find basis; "why simultaneous" conceptual |
| 34 | construct a CSCO | A | build and verify completeness (curriculum 2.9a) |
| 35 | degeneracy prevents unique labeling | B | witnessed: repeated eigenvalues |
| 36 | $[A,B]$ quantifies incompatibility | B | compute commutator; the quantification is conceptual |
| 37 | uncertainty as variances/commutators not error | B | compute Robertson; "not measurement error" conceptual |
| 38 | function of a diagonalizable operator | A | $f(A)$ via spectral / `MatrixFunction` |
| 39 | Hermitian vs unitary vs normal vs projection | B | test each property; taxonomy conceptual |
| 40 | unitary as $\exp$ of a Hermitian generator | A | $U=\exp(-iH)$, verify unitary |

### D. Measurement, probability, state update
| Q | Topic | Tier | Reason |
|---|-------|------|--------|
| 41 | probabilities, nondegenerate spectrum | A | $p_k=|\langle k|\psi\rangle|^2$ |
| 42 | probability rule, degenerate eigenspace | A | $p=\langle\psi|P|\psi\rangle$ |
| 43 | post-measurement state, nondegenerate | A | collapse to $|k\rangle$ |
| 44 | post-measurement state, degenerate | A | $P|\psi\rangle/\lVert\cdot\rVert$ |
| 45 | repeated measurement gives certainty | A | simulate: second outcome deterministic |
| 46 | $S_x$ disturbs a prepared $S_z$ eigenstate | B | witnessed by the SG chain; "why" conceptual |
| 47 | sequential Stern-Gerlach reveals noncommutativity | B | simulate; witnessed |
| 48 | expectation value on a pure state | A | $\langle\psi|A|\psi\rangle$ |
| 49 | variance and what zero variance implies | A | compute variance; zero $\Leftrightarrow$ eigenstate is checkable |
| 50 | well-defined means but no sharp joint values | B | compute means/variances; conceptual point |
| 51 | preparing an ensemble vs measuring one system | C | conceptual distinction |
| 52 | measurement transforms coherence in another basis | A | compute off-diagonal suppression (Lueders) |
| 53 | measurement as orthogonal projectors summing to $I$ | A | build $\{P_k\}$, verify $\sum P_k=I$ |
| 54 | extra freedom of a POVM | B | build a POVM; "freedom" conceptual |
| 55 | same probabilities, different post-states | B | construct two instruments, same POVM |
| 56 | operational distinguishability of nonorthogonal states | B | compute overlap/Helstrom; "operational meaning" conceptual |

### E. Unitary dynamics, two-level evolution
| Q | Topic | Tier | Reason |
|---|-------|------|--------|
| 57 | condition to preserve normalization | B | verify unitarity; the "condition" is a named concept |
| 58 | time-independent $H$ generates unitary evolution | A | $U=\exp(-iHt)$, verify unitary |
| 59 | solve Schrodinger by diagonalizing $H$ | A | diagonalize, evolve |
| 60 | physical content of accumulated relative phases | B | compute phases; "content" interpretive |
| 61 | oscillations when $H$ off-diagonal in prep basis | A | compute $P(t)$ |
| 62 | splitting sets the oscillation frequency | A | frequency $=\Delta E$ |
| 63 | spin-1/2 in a uniform field | A | build $H=-\vec\mu\cdot\vec B$ |
| 64 | Bloch precession and conserved expectations | A | compute precession, conserved component |
| 65 | resonance condition for a sinusoidal drive | B | compute resonance $\omega=\omega_0$; "condition" |
| 66 | rotating-wave / near-resonance simplification | B | compare full vs RWA; "how it simplifies" conceptual |
| 67 | exact Rabi vs first-order perturbation | B | compute both; conceptual contrast |
| 68 | adiabatic vs sudden evolution | B | simulate both regimes; conceptual contrast |
| 69 | Berry phase for a cyclic adiabatic loop | A | discrete Bargmann invariant (curriculum 7.13) |
| 70 | gate as controlled finite-time evolution | B | $U=\exp(-iHt)$; "interpret as a gate" conceptual |

### F. Spin-1/2, Pauli, qubits
| Q | Topic | Tier | Reason |
|---|-------|------|--------|
| 71 | Pauli matrices as spin components | A | build, verify |
| 72 | eigenvectors of $S_x,S_y$ from $S_z$ basis | A | compute eigenvectors |
| 73 | probabilities for $S_z$-up measured along $\hat n$ | A | $\cos^2(\theta/2)$ |
| 74 | Bloch sphere encodes Pauli expectations | A | compute $\langle\sigma_i\rangle$ |
| 75 | SU(2) spinors vs SO(3) rotations | B | compute the 2:1 map; the relation is a theorem |
| 76 | $2\pi$ rotation flips the sign | B | compute $-I$; "same physical state" conceptual |
| 77 | spin-1/2 vs a classical variable-projection vector | C | experimental/conceptual |
| 78 | SG as preparer and measurement device | C | conceptual framing |
| 79 | construct $S\cdot\hat n$ | A | build $\hat n\cdot\vec\sigma/2$ |
| 80 | eigenvalues/eigenvectors of $S\cdot\hat n$ | A | compute |
| 81 | field gradient couples spin to spatial path | C | spatial coupling, conceptual/semi-infinite-dim |
| 82 | spin measurement as projective in another eigenbasis | A | build projectors, measure |
| 83 | spinor vs 3-vector under rotation | B | compute half- vs full-angle; witnessed |
| 84 | Pauli commutation/anticommutation relations | A | verify the algebra (curriculum 2.2) |
| 85 | arbitrary $2\times2$ Hermitian as scalar + field | A | $a_\mu=\tfrac12\mathrm{Tr}[\sigma_\mu H]$ |
| 86 | qubit generalizes spin-1/2 beyond spin | C | conceptual |

### G. Angular momentum and rotations
| Q | Topic | Tier | Reason |
|---|-------|------|--------|
| 87 | angular-momentum commutation relations | A | verify $[J_i,J_j]=i\epsilon_{ijk}J_k$ |
| 88 | ladder operators fix the allowed $m$ | A | build $J_\pm$, show the range truncates |
| 89 | why the irrep has dimension $2j+1$ | B | build it (dim checks); "why" is rep theory |
| 90 | construct the standard $|j,m\rangle$ basis | A | build |
| 91 | matrix elements of $J_x,J_y$ from $J_\pm$ | A | $J_x=(J_++J_-)/2$ |
| 92 | orbital vs abstract spin representations | C | conceptual (orbital is infinite-dim) |
| 93 | integer $l$ vs half-integer $j$ | C | the half-integer side is computable (Q76 shows $2\pi\to-I$); the integer-$l$ restriction rests on spatial single-valuedness, which is not a finite-dim fact |
| 94 | rotation operator on a finite spin state | A | apply $\exp(-i\theta\,\hat n\cdot\vec J)$ |
| 95 | rotational invariance implies conservation | B | compute $[H,\vec J]=0$; Noether is conceptual |
| 96 | $H$ commutes with $J^2,J_z$ but not $L_z,S_z$ | B | build spin-orbit $H$, check; "what it means" interpretive |
| 97 | spin-orbit makes $J$ a good quantum number | B | diagonalize $H_{SO}$; witnessed |
| 98 | Clebsch-Gordan, product -> coupled basis | A | `ClebschGordan` built-in |
| 99 | which total $J$ from adding $j_1,j_2$ | A | compute $|j_1-j_2|..j_1+j_2$, check dim |
| 100 | singlet-triplet decomposition of two spin-1/2 | A | construct, identify by $S^2$ |
| 101 | rotational symmetry restricts vector-operator elements | B | Wigner-Eckart on an instance; theorem conceptual |
| 102 | finite-dim content of Wigner-Eckart | B | verify factorization for a fixed $j$; witnessed |
| 103 | selection rules from conservation + tensor structure | B | compute which elements vanish; witnessed |
| 104 | accidental vs symmetry-enforced degeneracy | B | build both, show the symmetry operator; distinction conceptual |

### H. Composite systems, tensor product, entanglement
| Q | Topic | Tier | Reason |
|---|-------|------|--------|
| 105 | tensor product not Cartesian product | B | compute $d_Ad_B$ vs $d_A+d_B$; "why" conceptual |
| 106 | composite dimension from subsystem dimensions | A | $\dim=\prod d_i$, verify |
| 107 | product vs entangled state | A | test Schmidt rank / factorizability |
| 108 | decide entanglement from the coefficient matrix | A | SVD, Schmidt rank |
| 109 | Schmidt decomposition | A | compute |
| 110 | Schmidt coefficients quantify entanglement | A | compute entanglement entropy |
| 111 | measurement on one subsystem updates the other | A | compute the conditional state |
| 112 | classical correlation vs entanglement | B | compute both; conceptual contrast |
| 113 | entanglement gives non-local correlations | B | compute CHSH; "why not local" conceptual |
| 114 | Bell singlet: perfect anticorrelation on any axis | A | compute $E(\hat a,\hat a)=-1$ |
| 115 | tensor structure makes local observables commute | A | verify $[A\otimes I,I\otimes B]=0$ |
| 116 | partial trace -> reduced density operator | A | compute |
| 117 | subsystem of a pure entangled state is mixed | A | compute reduced $\rho$, $\mathrm{Tr}[\rho^2]<1$ |
| 118 | what "local unitary" means | B | build $U_A\otimes U_B$; "local" conceptual |
| 119 | entangling unitary makes product -> entangled | A | apply CNOT, show entanglement |
| 120 | why CNOT is entangling despite basis-preserving | B | apply to $|{+}0\rangle$; "why" conceptual |
| 121 | build symmetric/antisymmetric subspaces | A | build projectors, apply (curriculum 11.13) |
| 122 | why identical fermions must be antisymmetric | C | spin-statistics is an axiom / QFT result |
| 123 | spin symmetry fixes spatial states of two fermions | A | construct, verify total antisymmetry |
| 124 | singlet/triplet -> spatial symmetry for two electrons | A | construct, verify |

### I. Density operators, mixed states, general measurements
| Q | Topic | Tier | Reason |
|---|-------|------|--------|
| 125 | valid density operator | A | test Hermitian, $\mathrm{Tr}=1$, PSD (curriculum 3.9) |
| 126 | rank-one density operator = pure state | A | build $|\psi\rangle\langle\psi|$, rank 1 |
| 127 | purity distinguishes pure from mixed | A | compute $\mathrm{Tr}[\rho^2]$ |
| 128 | ensemble -> density operator | A | $\rho=\sum p_i|\psi_i\rangle\langle\psi_i|$ |
| 129 | why different ensembles give the same $\rho$ | B | construct a coincidence; the general answer is the Hughston-Jozsa-Wootters theorem, which the cell does not capture |
| 130 | expectation via $\mathrm{Tr}[\rho A]$ | A | compute |
| 131 | unitary evolution of a density operator | A | $\rho\to U\rho U^\dagger$ |
| 132 | update $\rho$ when the outcome is known | A | $P_k\rho P_k/p_k$ |
| 133 | update $\rho$ when the outcome is ignored | A | $\sum_k P_k\rho P_k$ |
| 134 | decoherence in a basis vs ignorance of a pure state | B | show equal matrices; conceptual distinction |
| 135 | POVM probabilities without a unique post-state | B | compute POVM probs; conceptual |
| 136 | channel from Kraus operators | A | build, apply |
| 137 | Kraus trace-preservation constraint | A | verify $\sum K^\dagger K=I$ |
| 138 | depolarizing/dephasing/bit-flip on a qubit | A | apply, show Bloch contraction |

### J. Degeneracy, perturbations, effective subspaces
| Q | Topic | Tier | Reason |
|---|-------|------|--------|
| 139 | why nondegenerate PT fails in a degenerate space | B | show the vanishing denominator; "why" conceptual |
| 140 | degenerate PT = diagonalize in the subspace | A | build submatrix, diagonalize |
| 141 | choosing the zeroth-order basis | A | eigenvectors of $V$ in the subspace |
| 142 | partial vs complete lifting of degeneracy | B | compute the split spectrum; conceptual |
| 143 | symmetry predicts vanishing matrix elements | B | compute selection rules; witnessed |
| 144 | avoided crossing from a $2\times2$ effective $H$ | A | diagonalize (curriculum 2.13) |
| 145 | coupling of nearly degenerate states | A | compute the mixing |
| 146 | Zeeman effect as multiplet diagonalization | A | build Zeeman $H$, diagonalize |
| 147 | fine/hyperfine structure selects a coupled basis | A | build $\xi\,\vec L\cdot\vec S$, show $|j,m_j\rangle$ block-diagonal |
| 148 | when a finite truncation is justified | C | physical judgment, conceptual |
| 149 | errors from treating a truncation as closed | B | compare truncated vs full evolution (leakage); witnessed |
| 150 | time-dependent PT as amplitude flow | B | compute the flow; conceptual framing |

### K. Quantum information, gates, protocols, algorithms
| Q | Topic | Tier | Reason |
|---|-------|------|--------|
| 151 | qubit vs a probabilistic classical bit | C | conceptual |
| 152 | one-qubit gate as a $2\times2$ unitary | A | build, verify unitary |
| 153 | $n$-qubit state in a $2^n$-dim space | A | construct, show dimension |
| 154 | why reversible computation must be unitary | B | show non-unitary breaks norm; "why" conceptual |
| 155 | Hadamard creates superposition | A | apply $H$ to $|0\rangle$ |
| 156 | Pauli $X,Y,Z$ on basis states and phases | A | apply, show |
| 157 | controlled gate on basis vs superposition | A | apply CNOT, show entangling action |
| 158 | why an entangling gate is needed for universality | B | witness the necessity (local gates keep entanglement entropy 0, so entangled targets are unreachable); the sufficiency converse stays a theorem |
| 159 | rotations + one entangler give universality | B | recompile a random unitary; theorem conceptual |
| 160 | phase kickback into a relative phase | A | controlled-$U$ on an eigenstate |
| 161 | QFT on computational basis states | A | apply QFT |
| 162 | why periodicity is central to Shor | B | demonstrate period finding; "why central" conceptual |
| 163 | output-register measurement in period finding | B | simulate; conceptual |
| 164 | amplitude amplification as a 2D rotation | A | compute the rotation angle and amplitude (curriculum 23.5) |
| 165 | Grover as two reflections | A | build oracle + diffuser, show the rotation |
| 166 | no-cloning from linearity | B | show the cloning contradiction; the proof is the content |
| 167 | teleport an unknown qubit | A | build circuit, verify fidelity 1 |
| 168 | superdense coding | A | build circuit, verify |
| 169 | entanglement-assisted vs classical shared randomness | B | CHSH beats classical; conceptual |
| 170 | stabilizer measurement without measuring the state | A | measure syndrome, logical state intact (curriculum 20.1) |
| 171 | bit-flip vs phase-flip vs combined Pauli errors | A | apply $X,Z,Y$, show effect |
| 172 | repetition code protects one error type only | A | build, show it misses phase flips |
| 173 | why QEC must preserve logical coherence | B | show measuring logical info kills superposition; "why" conceptual |
| 174 | syndrome subspaces partition the space | A | build syndrome projectors, orthogonal decomposition |
| 175 | fault tolerance in terms of subspaces and projectors | C | conceptual framing |

### L. Conceptual and foundational
| Q | Topic | Tier | Reason |
|---|-------|------|--------|
| 176 | finite-dim version of the measurement problem | C | foundational/conceptual |
| 177 | Bell rules out predetermined values under locality | B | compute CHSH violation; the "ruling out" is interpretation |
| 178 | finite-dim content of contextuality | B | Mermin-Peres / KS contradiction (curriculum 25.7) |
| 179 | EPR argument with two spin-1/2 | C | the argument is philosophical |
| 180 | no-signaling and local measurements | B | reduced state independent of the remote choice; witnessed |
| 181 | collapse as information update without FTL signaling | C | interpretive |
| 182 | hidden variables vs the density-operator account | C | foundational/conceptual |
| 183 | quantum Zeno effect | A | simulate frequent measurement (curriculum 7.15) |
| 184 | classical limit of a large spin | B | compute large-$j$ coherent dynamics vs spin-1/2; contrast |
| 185 | kinematic vs dynamical features of the theory | C | conceptual taxonomy |

## What this means for the harvest decision (ties back to Q1)

The Q1 recommendation was: do not merge the bank, but harvest the genuinely new finite-dim topics.
This audit sharpens "genuinely new" into "genuinely new **and** Tier A", the items that are both
missing from the 213-item curriculum and answerable as a clean compute/verify cell:

1. **A compute-oriented foundation layer** (the audit's own "add Part 0" recommendation), from
   Tier-A bank items 3, 4, 7, 8, 10, 30, 32: coefficient extraction, change of basis, resolution of
   identity, dual vectors, subspace projectors, degenerate-eigenspace projectors, operator
   reconstruction.
2. **Zeeman and fine/hyperfine structure as multiplet diagonalization** (146, 147): both Tier A,
   both absent from the curriculum's perturbation cluster (Part 2), both a clean `Eigensystem` on a
   built Hamiltonian.
3. **Phase kickback** (160): Tier A, and independently flagged as missing by `Question-List-Audit.md`.

Two more topics are worth adding but are **Tier B**, so they must be re-cast to a concrete instance
with the conceptual half in prose (the Fermi-golden-rule treatment):

4. **Wigner-Eckart theorem and selection rules** (101-103): verify the CG-times-reduced-element
   factorization for a fixed $j$; state the theorem in prose.
5. **Sequential Stern-Gerlach / measurement disturbance** (46, 47): simulate the three-magnet chain.

The 24 Tier-C items (foundations, interpretation, "what is the distinction") are the part of the
bank that should stay in the companion file and never enter the compute/verify manual, because by
construction no cell answers them.

## Caveats and limits of this audit

- **The A/B boundary is a judgment, not a measurement.** About a dozen "construct then interpret"
  items (16, 18, 28, 33, 49, 55, 70, 84, 118, 164, 170) could be argued one tier either way. The
  48%-56% band is the honest read; the point estimate is 51%.
- **Tier A asserts a WL/QF primitive exists.** Each Tier-A call assumes the obvious primitive is
  available (`ClebschGordan`, `MatrixExp`, `Eigensystem`, QF stabilizer projectors, QF channels).
  These are all real and used elsewhere in the curriculum; none was invented for this audit. The
  audit did not boot a kernel, because it classifies questions, it does not write WL code.
- **"Fully answerable" is about the *answer*, not difficulty.** A Tier-A item can still be a large
  computation (five-qubit code, HHL); the tier only says the rendered object is the whole answer.
- **The bank is finite-dim-filtered already**, so almost nothing fell out for being intrinsically
  infinite-dimensional; the C items are conceptual, not out-of-scope.

## Bottom line

About half the bank (48%-56%, best estimate 50%, 93 questions) is directly harvestable as full
compute/verify curriculum items; another third (69) is usable only if re-cast to a discriminating
concrete instance with the conceptual remainder written as prose; an eighth (23) is oral-exam
conceptual material that does not belong in a learning-by-computing manual at all. The
computationally rich sections (density operators 79%, composite systems 70%, algorithms 60%) are
where the harvest pays off; the foundational sections (rays 25%, foundations 10%) are where the bank
stays a companion, not a source.
