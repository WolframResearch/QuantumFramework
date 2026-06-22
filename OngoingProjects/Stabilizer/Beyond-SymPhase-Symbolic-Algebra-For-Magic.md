---
Title: 'Beyond SymPhase: can symbolic algebra compress magic (T-gate) branching?'
Author: investigation in the Wolfram QuantumFramework
Description: 'A self-contained, computation-first study of whether the stabilizer-simulation idea of carrying randomness as algebra extends from measurement branching to non-Clifford (T-gate) branching. States the complexity wall precisely, maps the literature of algebraic representations of Clifford+T (path sums, F_2 Gauss sums, LIMDD, exact rings, ZX, stabilizer rank/extent), and builds a working symbolic path-sum simulator in Wolfram Language whose rewriting folds Clifford branching to a closed form and isolates the irreducible magic. Every equation is re-derived and every claim verified against exact amplitudes.'
---

# Beyond SymPhase: symbolic algebra for magic branching

## The question

A stabilizer circuit that measures undetermined qubits forks into a tree of $2^r$ outcome histories. The SymPhase technique (Fang and Ying, arXiv:2311.03906) refuses to build that tree: it carries each undetermined outcome as a **fresh formal bit** $s_k \in \mathbb{F}_2$ written into a stabilizer sign, traverses the circuit once, and recovers any branch by substituting values for the bits. The randomness has become algebra, a degree-$1$ polynomial over $\mathbb{F}_2$ in the sign column, and one pass serves all $2^r$ shots.

A non-Clifford $T$ gate forks the same way, but the branches are not classical outcomes: $T\lvert\psi\rangle$ is a coherent superposition of stabilizer states (the source of "magic"), and $k$ gates give up to $2^k$ terms, the stabilizer rank. The question is whether the SymPhase move transfers: can magic branching also be carried as algebra and folded, instead of enumerated?

The answer depends on **which question you ask the classical simulator**, and the three are not equally hard. A *sample* from the output distribution (weak simulation) is easier than an *exact amplitude or probability* (strong simulation), which is easier than nothing only for structured inputs; a *local expectation* is easier still. This report is about the hardest one, the exact amplitude, because that is what a symbolic Gauss sum computes and where the algebra either folds or does not.

> **The worst case is an unbreakable wall, but the reachable territory is far larger than "diagonal circuits."** The right generalization of SymPhase is the **symbolic path sum**: an amplitude written as a Gauss sum $\sum_v \omega^{\varphi(v)}$ over Hadamard branch variables $v$, with $\omega = e^{i\pi/4}$ and $\varphi$ a low-degree polynomial over $\mathbb{F}_2$ valued in $\mathbb{Z}_8$, carried symbolically and *reduced by rewriting* rather than enumerated. SymPhase is its degree-$1$, sign-only restriction. A path-sum simulator built here in Wolfram Language reproduces exact Clifford+$T$ amplitudes, and its Hadamard-elimination rule **folds all Clifford branching to a closed form** (a $9$-variable sum collapses to $1$, a $10$-variable sum to $2$) while leaving **exactly the magic** as an irreducible residue (a magic $9$-variable sum does not reduce). The count of surviving variables is the algebraic shadow of the stabilizer rank.

All numerical claims are verified against `Method -> "Schrodinger"` dense amplitudes on paclet `Wolfram/QuantumFramework 2.0.0`; every gate rule below is independently re-derived from the gate matrix. Driver scripts are in `beyond-symphase-scripts/` next to this file.

---

## 1. The wall that cannot be moved

Two facts bound everything below, so the positive results are not mistaken for a free lunch.

- **Exact strong simulation of Clifford+$T$ is $\#\mathrm{P}$-hard.** Computing a single output amplitude exactly is as hard as the permanent: the permanent, the canonical $\#\mathrm{P}$-complete quantity, is *itself* a quantum transition amplitude of a linear-optical network (Aaronson, arXiv:1109.1674), and the qubit Clifford+$T$ case is no easier. Equivalently (Montanaro, arXiv:1607.08473), the amplitude of a $\{H, Z, CZ, CCZ\}$ circuit is, up to a power-of-$\sqrt2$ normalization, the Gauss sum $\mathrm{gap}(f) = \sum_x (-1)^{f(x)}$ of a polynomial $f$ of degree $\le 3$ over $\mathbb{F}_2$. At degree $2$ the sum is computable in $O(n^3)$, which is Gottesman-Knill in disguise; at degree $3$ it is $\#\mathrm{P}$-hard.
- **Clifford+$T$ is universal.** A polynomial-time symbolic method that folded generic $T$ branching exactly would put $\mathrm{BQP}$ in $\mathrm{P}$ and collapse the polynomial hierarchy.

No representation in this report beats that. The entire game is **detecting and exploiting structure**, and the value of symbolic algebra is that it makes the structure of the Gauss sum explicit and exact rather than hidden in a $2^n$ vector.

---

## 2. The candidate landscape

Two literature surveys were run, papers fetched and read. Ranked by promise as a symbolic-algebra route that compresses *structured* magic beyond eager $2^t$ enumeration:

| Rank | Representation | arXiv | Precise claim | Compression regime |
|---|---|---|---|---|
| **1** | **Symbolic path sums + confluent rewriting** | **1805.06908** (Amy), **2408.02778** (Amy-Stinchcombe) | A path sum $s\sum_V (-1)^P\lvert O\rangle$ with $P,O$ multilinear over $\mathbb{F}_2[V]/(v^2-v)$; a confluent, degree-bounded rewrite system reduces it *without expanding*. Drives hidden-shift circuits with $>1000$ $T$-gates to the answer (zero residual variables) in a linear number of steps: genuine **polynomial time** on a family that belongs to no prior tractable class. | Structured magic (bounded-degree bent oracles); degrades to a stabilizer-sum fallback in the worst case. Named future direction: sum over varieties via Gröbner bases / point counting. |
| **2** | **$\mathbb{F}_2$ low-degree polynomial / Gauss sum** | **1607.08473** (Montanaro) | Amplitude $=\mathrm{gap}(f)/2^{h/2+\ell}$, $f$ degree $\le 3$ over $\mathbb{F}_2$; each internal Hadamard adds one summation variable. Tractable when a *hitting set* of the degree-$3$ terms is small, or the number of *essential variables* (under $GL_n(\mathbb{F}_2)$) is small: both strictly generalize "low $T$-count." | Sparse / low-treewidth phase polynomials. The exact theory underneath candidate 1. |
| **3** | **LIMDD (Pauli-Tower decision diagram)** | **2108.00931** | Poly-size LIMDDs *strictly contain* stabilizer $\cup$ QMDD $\cup$ MPS states; every stabilizer state is $O(n)$ nodes via local-Pauli edge labels; provably simulate $W$-states and Hamming-weight-controlled Cliffords in poly time where Clifford+$T$ needs $\Omega(2^n)$. Whether it dominates stabilizer rank is **open** (numerically suggestive). | Circuits whose intermediate states are local-Clifford-related sub-blocks, even at high $T$-count. The data-structure realization of the same local-Pauli factoring. |
| **4** | **Exact ring $\mathbb{Z}[\omega,1/\sqrt2]$** | **1206.5236** (Kliuchnikov-Maslov-Mosca), Giles-Selinger | Every Clifford+$T$ amplitude is an exact element of $\mathbb{Z}[\omega]$, $\omega=e^{i\pi/4}$. | Not a compressor: the exactness substrate that makes 1-3 exact (no floating-point error, exact rank and equality tests, closed-form Gauss sums). Free in Wolfram Language. |
| **5** | **ZX-calculus + stabilizer decomposition** | **2109.01076** (Kissinger-van de Wetering) | Interleave a $\chi = 2^{\alpha t}$ ($\alpha \approx 0.40$) stabilizer decomposition with ZX simplification; cost independent of qubit count. | Large on *structured* circuits ($>6$ orders of magnitude on hidden shift), but only $\sim 10\times$ (constant) on random magic, where the authors recommend T-count reduction then a tableau simulator. |
| - | **Stabilizer rank / extent / robustness** | 1601.07601, 1808.00128, 1609.07488 | The eager baseline ($\chi = O(2^{\beta t})$, $\beta \le 0.47$) plus exact magic *measures*: approximate rank scales as $2^{0.23 t} \approx 1.17^t$ (Bravyi-Gosset), robustness as $1.685^t$ down to a proven $1.207^t$ floor (Howard-Campbell), both convex programs exact in WL. | Shrink the *base* of the exponential; structure-blind. Excellent as measures, not as structured compressors. |
| - | **MPS / tensor networks** | 2404.18768 | Magic (stabilizer Rényi entropy) is **independent of bond dimension**: modest-$\chi$ MPS are already Haar-magical. | Compresses *entanglement*, an axis orthogonal to magic. Not a magic-branching compressor. |

The leader for a symbolic-algebra attack is candidate 1/2: the **symbolic path sum**, which is SymPhase generalized from degree-$1$ sign polynomials to the full Gauss-sum object, with compression set by the algebraic structure of the phase polynomial. The rest of the report demonstrates it.

---

## 3. What the framework already gives under symbolic algebra

Four properties the framework's own objects deliver once the inputs are symbolic, none of which the eager `StabilizerFrame` exploits. (Script: `exp-symbolic-angles-and-rank.wls`.)

### 3a. The frame carries a symbolic angle: one traversal, every angle setting

`StabilizerFrame` propagates a symbolic phase, so a parametrized phase gate keeps its angle and the state vector comes out as a closed form in it:

```wl
fSym = PauliStabilizer[1][{"H", "P"[\[Theta]]}];
fSym["Coefficients"]                              (* {(1 + E^(I/2 \[Theta]))/2, (1 - E^(I/2 \[Theta]))/2} *)
Simplify[Normal[fSym["StateVector"]]]             (* {1/Sqrt[2], E^(I \[Theta]/2)/Sqrt[2]} : closed form in \[Theta] *)
```

An expectation as a function of the angle then follows from that one traversal, the magic analogue of SymPhase's "one pass serves every shot," now over a continuum of gate settings (the object a variational circuit needs):

```wl
zTheta = FullSimplify[Conjugate[#] . (PauliMatrix[3] . #) &[
    Normal[PauliStabilizer[1][{"H", "P"[\[Theta]], "H"}]["StateVector"]]], Assumptions -> \[Theta] \[Element] Reals];
(* zTheta = Cos[\[Theta]/2] : evaluate at any angle with no re-traversal *)
```

The half-angle is the framework's convention, $\texttt{"P"}[\theta] = \mathrm{diag}(1, e^{i\theta/2})$.

### 3b. Structured magic collapses under `Simplify`

Several phase gates on one qubit eagerly produce $2^k$ components, but the materialized state simplifies to **rank $1$**, the magic absorbed into a single phase, and commuting gates on distinct qubits factor exactly as a product, $O(k)$ rather than $2^k$:

```wl
Simplify[Normal[PauliStabilizer[1][{"H", "P"[a], "P"[b], "P"[c]}]["StateVector"]]]
(* {1/Sqrt[2], E^(I (a + b + c)/2)/Sqrt[2]} : 8 eager components, one phase *)
```

### 3c. The eager frame over-counts the rank; exact linear algebra reveals it

The component list has length $2^{\#T}$, but the components are often linearly dependent, and any $n$-qubit state has stabilizer rank at most $2^n$ (the computational basis is stabilizer states) regardless of $T$-count. An interleaved one-qubit circuit keeps $8$ components spanning a $2$-dimensional space; a four-$T$ two-qubit circuit keeps $16$ spanning a $4$-dimensional space:

```wl
compMatrix[f_] := Chop[N[Normal[#["State"]["StateVector"]] & /@ f["Stabilizers"]]];
fint = PauliStabilizer[1][{"H", "T", "H", "T", "H", "T"}];
{fint["Length"], MatrixRank[compMatrix[fint]]}     (* {8, 2} : 6 of 8 components are redundant *)
```

Whenever $\#T > n$ the frame is provably wasteful; exact linear algebra over $\mathbb{Z}[\omega]$ caps it at the span dimension. Reaching the *optimal* rank ($\approx 2^{0.4 t}$) is the NP-hard $\ell_0$ problem of Bravyi-Smith-Smolin; capping at the span is the free part.

### 3d. Exact cyclotomic arithmetic and closed-form Gauss sums

Every amplitude is an exact element of $\mathbb{Z}[\omega, 1/\sqrt2]$, and a structured Gauss sum evaluates to a closed form rather than a float:

```wl
RootReduce[Normal[PauliStabilizer[1][{"H", "T"}]["StateVector"]]]   (* {1/Sqrt[2], 1/2 + I/2} : exact *)
Table[{n, Simplify[Sum[Exp[I Pi k^2/n], {k, 0, n - 1}]]}, {n, 2, 6}]
(* {{2,1+I},{3,1},{4,2 E^(I \[Pi]/4)},{5,1},{6,(1+I) Sqrt[3]}} : quadratic Gauss sums, closed form *)
```

---

## 4. A symbolic path-sum simulator in Wolfram Language

The leading candidate, built and verified. (Script: `pathsum-simulator-and-elimination.wls`.)

### 4a. The representation

Each qubit's value is an $\mathbb{F}_2$ multilinear polynomial in Hadamard branch variables $x_i$; the global phase is a $\mathbb{Z}_8$-valued pseudo-Boolean polynomial. Each gate acts locally, and each rule is the gate matrix read off on basis states (all four re-derived to zero):
- $H$ on qubit $q$: introduce a fresh $v$, add $4\, e_q v$ to the phase, set $e_q := v$. This is $H\lvert e\rangle = 2^{-1/2}\sum_v \omega^{4 e v}\lvert v\rangle$.
- $T$: add $e_q$, since $T = \mathrm{diag}(\omega^0, \omega^1)$.
- $CZ$: add $4\, e_a e_b$, since $CZ = (-1)^{e_a e_b} = \omega^{4 e_a e_b}$. $CCZ$: add $4\, e_a e_b e_c$.
- $CNOT$: the linear update $e_b := e_a \oplus e_b$.

The amplitude is the Gauss sum

$$
\langle y \lvert C \rvert 0\cdots0\rangle = 2^{-h/2}\sum_{v\,:\,E(v)=y}\omega^{\varphi(v)},
$$

with $h$ the number of Hadamards and $E(v)$ the vector of final qubit values. For $H\,T\,H$ on $\lvert 0\rangle$ the phase is $\varphi = x_1 + 4 x_1 x_2$ and the amplitudes are $(1 \pm \omega)/2$, matching the dense result. The builder was checked against `Method -> "Schrodinger"` on a one-qubit interleaved circuit, a two-qubit `CNOT` circuit, a diagonal circuit, and a `CCZ` gadget: **all amplitudes match**, and an independent reviewer reproduced the match on three further circuits to one part in $10^{17}$.

### 4b. Elimination: Clifford branching folds to a closed form, magic does not

The Hadamard-elimination rule (the [HH] rule of Amy, arXiv:1805.06908) is the place where branching becomes algebra. Write the phase as $\varphi = v\, C_v + \varphi_0$, splitting off the cofactor $C_v = \varphi\rvert_{v=1} - \varphi\rvert_{v=0}$ of an *internal* variable $v$ (one absent from the output). The sum over $v$ is

$$
\sum_{v \in \{0,1\}} \omega^{v C_v + \varphi_0} = \omega^{\varphi_0}\bigl(1 + \omega^{C_v}\bigr).
$$

When $C_v = 4Q$ is a pure sign with $Q$ a linear Boolean form, this is $2\,\omega^{\varphi_0}[Q = 0]$: the variable is summed away in closed form, leaving the linear constraint $Q = 0$, with **no enumeration**. A $T$ gate makes the cofactor $\mathbb{Z}_8$-valued instead of a sign, so $1 + \omega^{C_v}$ does not collapse and the variable **survives**. Reducing to a fixed point and checking the residual amplitude against the dense result:

| circuit | branch variables before | after elimination | matches exact? |
|---|---|---|---|
| Clifford $H(ZH)^8$, 1 qubit | 9 | **1** | yes |
| Clifford 2-qubit deep, $(CZ\,H_1H_2)^4$ | 10 | **2** | yes |
| Magic $H(TH)^8$, 1 qubit | 9 | **9** | yes |
| Mixed, Clifford with one $T$ | 5 | **3** | yes |

This is the computational statement of "branching becomes algebra, generalized past SymPhase." Clifford structure of arbitrary depth folds entirely into a closed form ($9$ variables to $1$, $10$ to $2$, the residue being just the output qubits); the $T$ gates are the only thing the algebra cannot fold, and the count of surviving variables is the algebraic shadow of the stabilizer rank. SymPhase is the corner of this picture with no $T$ gates, where the only symbols are measurement signs.

The rule implemented here handles the sign-cofactor ($\mathbb{Z}_2$) case. The complete confluent system (Amy-Stinchcombe, arXiv:2408.02778) adds the $\mathbb{Z}_4$ case ($S$ gates), quadratic completion, and multi-variable constraints, and is what drives structured-magic families all the way to the answer. Porting it to Wolfram Language, using native $\mathbb{F}_2$ polynomial reduction and `GroebnerBasis` for the variety-valued generalization, is the natural build.

---

## 5. Verdict

**The wall is real; the ceiling is high.** Exact single-amplitude extraction for generic Clifford+$T$ is $\#\mathrm{P}$-hard, and no symbolic method beats that. But the tractable side is not "diagonal circuits only": path-sum rewriting folds *all* Clifford structure, of any depth, into a closed form, and isolates the magic as a residue whose size is the real cost. The framework leaves this on the table: it already carries symbolic angles (closed-form parametric expectations), structured magic already collapses under `Simplify`, the eager component count already over-counts the rank (capped exactly by linear algebra over $\mathbb{Z}[\omega]$), and exact Gauss sums already evaluate to closed forms, yet the `StabilizerFrame` doubles its component list per $T$ regardless of structure.

**The object to build is the symbolic path sum** (Amy; Montanaro's $\mathbb{F}_2$ Gauss-sum theory underneath), with LIMDD as the data-structure alternative and exact $\mathbb{Z}[\omega]$ arithmetic as the substrate. Its compression is real and structure-dependent: polynomial for low-hitting-set, few-essential-variable, or structured-bent-oracle circuits, and degrading gracefully to the stabilizer-sum fallback in the worst case.

**Staged build plan.**
1. **Exact rank cap on `StabilizerFrame`**: reduce the $2^{\#T}$ eager components to the span dimension by exact linear algebra over $\mathbb{Z}[\omega]$ (demonstrated here). Bounds the frame by $2^n$.
2. **Phase-polynomial compression** for the diagonal (CNOT-dihedral) fragment: store one degree-$\le 3$ $\mathbb{Z}_8$ polynomial instead of $2^k$ components (rank $1$, exact).
3. **Full symbolic path-sum simulator** with the confluent rewrite system: the prototype here proves the core works; the complete system adds the $\mathbb{Z}_4$, quadratic, and multi-constraint rules and the variety-valued (`GroebnerBasis`) generalization, targeting structured-magic families in polynomial time.
4. **Exact magic measures**: stabilizer extent and robustness as exact convex programs over the stabilizer polytope, a resource-theory contribution distinct from the simulator.

In one line: the worst case is a $\#\mathrm{P}$-hard Gauss sum no algebra can fold, but wherever the phase polynomial has structure, the symbolic path sum rewrites its Clifford part to a closed form and leaves only the magic, and exact polynomial-and-cyclotomic algebra is an unusually good home for that, one the current framework does not yet use.

---

## References

- M. Amy, *Towards Large-scale Functional Verification of Universal Quantum Circuits*, QPL 2018, arXiv:1805.06908. Path sums; rewriting hidden-shift benchmarks to their answer.
- M. Amy, M. Stinchcombe, *Polynomial-Time Classical Simulation of Hidden Shift Circuits via Confluent Rewriting of Symbolic Sums*, arXiv:2408.02778. The confluent, degree-bounded rewrite system; polynomial time on structured magic; the variety-valued future direction.
- A. Montanaro, *Quantum circuits and low-degree polynomials over $\mathbb{F}_2$*, arXiv:1607.08473. Amplitude as $\mathrm{gap}(f)$ of a degree-$\le 3$ $\mathbb{F}_2$ polynomial; hitting-set and essential-variable tractability; $\#\mathrm{P}$-hardness.
- L. Vinkhuijzen, T. Coopmans, D. Elkouss, V. Dunjko, A. Laarman, *LIMDD: A Decision Diagram for Simulation of Quantum Computing Including Stabilizer States*, arXiv:2108.00931. Poly-size LIMDDs strictly contain stabilizer $\cup$ QMDD $\cup$ MPS.
- V. Kliuchnikov, D. Maslov, M. Mosca, arXiv:1206.5236; B. Giles, P. Selinger, arXiv:1212.0506. The exact Clifford+$T$ number ring $\mathbb{Z}[\omega,1/\sqrt2]$.
- A. Kissinger, J. van de Wetering, *Simulating quantum circuits with ZX-calculus reduced stabiliser decompositions*, arXiv:2109.01076.
- S. Bravyi, D. Gosset, arXiv:1601.07601; S. Bravyi et al., arXiv:1808.00128 (stabilizer extent); M. Howard, E. Campbell, arXiv:1609.07488 (robustness of magic). The eager rank baseline and exact magic measures.
- S. Aaronson, *A linear-optical proof that the permanent is $\#\mathrm{P}$-hard*, arXiv:1109.1674. The permanent as a quantum amplitude; the hardness wall.
- M. Frau, P. S. Tarabunga, et al., *Non-stabilizerness versus entanglement in matrix product states*, arXiv:2404.18768. Magic is independent of bond dimension.
- Y. Fang, M. Ying, *SymPhase: Phase Symbolization for Fast Simulation of Stabilizer Circuits*, arXiv:2311.03906. The degree-$1$ measurement-branching case this report generalizes.

Scripts: `beyond-symphase-scripts/{exp-symbolic-angles-and-rank.wls, pathsum-simulator-and-elimination.wls}`, verified against `Wolfram/QuantumFramework 2.0.0`. A companion note, `SymPhase-vs-StabilizerRank-Algebraic-Branching.md`, gives the worst-case analysis (why the sign trick alone reaches only the diagonal fragment) that this report extends.
