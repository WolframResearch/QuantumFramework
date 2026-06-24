# A Wolfram Stabilizer Package — Synthesis from 28 Papers

**Goal:** Distill what a Wolfram Language stabilizer package should *contain*, *compute*, and *expose symbolically*, based on a complete reading of all 28 TeX papers in `audit/Stabilizer/tex`.

**Bias of this report.** Wolfram Language is a symbolic algebra system on top of a numeric kernel; existing simulators (`Stim`, `CHP`, `GraphSim`, `Qiskit`) are bitwise-optimized C/C++/Rust. A Wolfram package that tries to win on raw qubit count will lose. A Wolfram package wins on (a) *symbolic* phases / parameters / outcomes, (b) interconversion of representations (state vector ↔ tableau ↔ graph ↔ quadratic form), (c) exact closed-form inner products, distances, and entanglement measures, (d) integration with `QuantumState`, `QuantumOperator`, `QuantumCircuitOperator` already in QuantumFramework. Every recommendation below is graded against that bias.

---

## Paper inventory

| Tag | Paper | Role |
|---|---|---|
| Got97 | Gottesman, *Stabilizer Codes and Quantum Error Correction* (PhD thesis, quant-ph/9705052) | Foundational. Stabilizer formalism, $[[n,k,d]]$ codes, normalizer, GF(4) connection. |
| Got98 | Gottesman, *The Heisenberg Representation of Quantum Computers* (quant-ph/9807006) | Heisenberg picture; Clifford rules; Knill's theorem. |
| Got00 | Gottesman, *An Introduction to Quantum Error Correction* (quant-ph/0004072) | Pedagogical intro; CSS construction; quantum Hamming bound. |
| Got09 | Gottesman, *An Introduction to Quantum Error Correction and Fault-Tolerant Quantum Computation* (0904.2557) | Extended QEC + fault-tolerance review. |
| AarGot04 | Aaronson & Gottesman, *Improved Simulation of Stabilizer Circuits* (quant-ph/0406196) | **Tableau algorithm with destabilizers**; canonical form H-C-P-C-P-C-H-P-C-P-C; $\oplus L$-completeness; CHP. |
| AndBri05 | Anders & Briegel, *Fast Simulation … Graph State Representation* (quant-ph/0504117) | **Graph-state + VOP** representation, $\mathcal O(N\log N)$ space. |
| Beaudrap11 | de Beaudrap, *A linearized stabilizer formalism for systems of finite dimension* (1102.3354) | Qudits via Weyl operators; eliminates quadratic phases. |
| Biswas24 | Biswas, *On Classical Simulation of Quantum Circuits Composed of Clifford Gates* (2405.13590) | Pedagogical walkthrough with explicit conjugation rules. |
| DehMoo03 | Dehaene & De Moor, *The Clifford group, stabilizer states, and linear and quadratic operations over GF(2)* (quant-ph/0304125) | Clifford group as $(2n+1)\times(2n+1)$ symplectic matrix + quadratic form $h$. |
| EllEasCav08 | Elliott, Eastin & Caves, *Graphical description of Pauli measurements on stabilizer states* (0806.2651) | Graph-rules for $X/Y/Z$ Pauli-product measurements. |
| FangYing23 | Fang & Ying, *SymPhase: Phase Symbolization for Fast Simulation of Stabilizer Circuits* (2311.03906) | **Symbolic phases** (Pauli faults as bit-symbols); decoupling sampling from circuit traversal. |
| GarMar15 | García & Markov, *Simulation of Quantum Circuits via Stabilizer Frames* (Quipu, 1712.03554) | **Stabilizer frames** = superpositions of stabilizer states; non-Clifford support. |
| GarMarCro12 | García, Markov & Cross, *Efficient Inner-product Algorithm for Stabilizer States* (1210.6646) | Closed-form inner product, canonical $H$-$C$-$CZ$-$P$-$H$ form, nearest neighbors. |
| Gid21 | Gidney, *Stim: a fast stabilizer circuit simulator* (2103.02202) | Modern reference simulator; SIMD layout; noise + sampling APIs. |
| HosDehMoo04 | Hostens, Dehaene & De Moor, *Stabilizer states and Clifford operations for systems of arbitrary dimensions, and modular arithmetic* (quant-ph/0408190) | Qudit Clifford via $\mathbb{Z}_d$, $\mathbb{Z}_{2d}$ matrices; Smith normal form for stabilizer expansion. |
| KocHuaLov17 | Kocia, Huang & Love, *Discrete Wigner Function Derivation of the Aaronson-Gottesman Tableau Algorithm* (1703.04630) | Wigner function ↔ tableau; symplectic phase-space view. |
| KoeSmo14 | Koenig & Smolin, *How to efficiently select an arbitrary Clifford group element* (1406.2170) | $O(n^3)$ subgroup algorithm + symplectic transvections; bijection $i\leftrightarrow C_n$. |
| MonPar23 | Mondal & Parhi, *Quantum Circuits for Stabilizer Error Correcting Codes* (2309.11793) | Encoder/decoder synthesis from $H_s$ standard form. |
| Mueller26 | Müller et al., *PauliEngine: High-Performant Symbolic Arithmetic for Quantum Operations* (2601.02233) | Binary symplectic + symbolic coefficients; commutator via popcount; DLA / Lie-closure benchmarks. |
| Paler14 | Paler et al., *Software Pauli Tracking for Quantum Computation* (1401.5872) | Pauli-frame tracking algorithm with state $\in\{I,X,Z,XZ\}$ per qubit. |
| PatGuh26 | Patil & Guha, *Clifford Manipulations of Stabilizer States* (2312.02377) | Karnaugh-map derivation of CHP-rules; cluster-state graph rule book. |
| PayWin24a | Paykin & Winnick, *Qudit Quantum Programming with Projective Cliffords* (2407.16801) | Type system / lambda calculus $\lambda^{\mathcal Pc}$ for projective Cliffords. |
| Reid24 | Berg, Cross, Wood et al., *Tableaux Manipulation* (2404.19408) | Tableau-by-pivot synthesis; comparison of empirical 2Q gate ratios. |
| RuhDev25 | Ruh & Devitt, *Quantum Circuit Optimisation and MBQC Scheduling with a Pauli Tracking Library* (2405.03970) | $O(n^2)$ Pauli tracking; MBQC measurement scheduling. |
| WinPay24b | Winnick & Paykin, *Condensed Encodings of Projective Clifford Operations in Arbitrary Dimension* (2407.16861) | Even/odd $d$ unified encoding; phase-correction cocycle; condensed product $\star$. |
| Winderl23 | Winderl et al., *Architecture-Aware Synthesis of Stabilizer Circuits from Clifford Tableaus* (2309.08972) | Steiner-tree-aware synthesis preserving connectivity. |
| Yashin25 | Yashin, *A streamlined demonstration that stabilizer circuits simulation reduces to Boolean linear algebra* (2504.14101) | Tableau formalism extended to **arbitrary stabilizer (Clifford) channels** via Choi tableaux. |
| deSilSalYin23 | de Silva, Salmon & Yin, *Fast algorithms for classical specifications of stabiliser states and Clifford gates* (2311.10357) | 10 fast interconversion algorithms across {amplitudes, quadratic-form, check-matrix} × {unitary, tableau}. |

---

## 1 — Core data structures the package must expose

### 1.1 Pauli string

The atomic object. Three coexisting representations, all interconvertible:

1. **Symbolic tensor**, `KroneckerProduct[X, Z, I, Y]` or `QuantumOperator[{"X","Z","I","Y"}]` — friendly, slow, used for printing.
2. **Symplectic bit-vector pair** `(x, z) ∈ 𝔽₂ⁿ × 𝔽₂ⁿ`, where `(0,0)=I, (1,0)=X, (0,1)=Z, (1,1)=Y` (Got97 §4; AarGot04 §2; Mueller26; Yashin25). Multiplication is XOR plus a phase correction.
3. **Phase-tracked symplectic** `(ε, δ, x, z)` where total operator is $i^\delta(-1)^\epsilon \tau_{(x,z)}$. The $\delta$-versus-$\epsilon$ split (DehMoo03 Lemma 1) is *not cosmetic*: it makes the binary phase calculus work mod 2 even when $d$ is even.

**Symbolic value-add (Wolfram-specific):** allow each of $\epsilon, \delta$ to be a **symbolic variable** (FangYing23 §3, Mueller26 §2). For a noise channel $\mathcal E(\rho) = (1-p)\rho + pX\rho X$ the corresponding Pauli is $X^{s}$ with $s$ a `Symbol[]` having sampling distribution `BernoulliDistribution[p]`. Sampling outcomes then becomes substitution (`ReplaceAll`) plus matrix-multiply mod 2 — pure symbolic algebra (FangYing23 Eq 5).

**Phase-correction cocycle.** Multiplication $P(u|c)\cdot P(u'|c') = i^{\beta(u,u')} P(u\oplus u'|c\oplus c')$ with
$\beta(u,u') = \langle z\oplus z', x\oplus x'\rangle - \langle z,x\rangle - \langle z',x'\rangle - 2\langle z',x\rangle$ (Yashin25 Eq 4; WinPay24b §3). This $\beta$ is the load-bearing piece for *all* qudit-aware code — write it once, reuse everywhere.

### 1.2 Tableau (extended / improved)

Aaronson–Gottesman tableau (AarGot04 §3): a $2n\times(2n+1)$ binary matrix
$$\begin{pmatrix}
\bar X & \bar Z & \bar r \\
X & Z & r
\end{pmatrix}$$
Top half = **destabilizers** (added to halve measurement cost from $O(n^3)$ to $O(n^2)$). Bottom = stabilizer generators. The single bit $r_i \in\{0,1\}$ is a sign; for symbolic packages it should be a `Symbol[]` holding a polynomial over $\mathbb{F}_2$ (FangYing23 Eq 4).

**Invariants** (AarGot04 Prop 2, replicate in `VerifyTableau[]`):
- $R_{n+1},\ldots,R_{2n}$ generate the stabilizer; $R_1,\ldots,R_{2n}$ generate $\mathcal P_n$.
- $R_1,\ldots,R_n$ commute pairwise.
- $R_h$ anticommutes with $R_{h+n}$.
- $R_i$ commutes with $R_{h+n}$ for $i\ne h$.

These four conditions are equivalently $\mathsf U \mathsf J \mathsf U^T = \mathsf J$ — symplectic over $\mathbb F_2$ (Yashin25 §2.3).

### 1.3 Graph state representation (Anders–Briegel)

Every stabilizer state is local-Clifford-equivalent to a graph state (AndBri05 §2; cite van den Nest, Dehaene, De Moor 2004). Stored as
1. an adjacency list of the graph $G$, and
2. a list of $n$ vertex operators (VOPs), each one of the **24 single-qubit Clifford operators** indexed `0..23` (AndBri05 §2 footnote, And05).

Memory is $\mathcal O(N\overline d)$ with $\overline d=$ avg degree, which is $\mathcal O(N\log N)$ for QEC and entanglement-purification circuits — beating the $\mathcal O(N^2)$ tableau by orders of magnitude (AndBri05 §6).

In Wolfram, the graph naturally fits into `Graph[]` with `VertexLabels -> { i -> VOPs[[i]] }`, which gives the user `GraphPlot` for free.

### 1.4 Quadratic-form triple (Dehaene–De Moor)

$$|s\rangle \propto \sum_{\vec z\in V} (-1)^{Q(\vec z)} i^{\ell(\vec z)} |\vec z + \vec z_0\rangle$$

with $V\subset\mathbb F_2^n$ vector subspace, $\vec z_0$ a shift, $Q$ quadratic, $\ell$ linear (DehMoo03 §4 Theorem 3; deSilSalYin23 Theorem 2.5). The package should support this representation explicitly. It is *the* form most natural for symbolic computation: amplitudes are exact, arithmetic is over $\mathbb F_2$, and $Q,\ell$ remain symbolic under Clifford action.

### 1.5 Stabilizer frame (Quipu)

A list of stabilizer-state generators sharing a global phase (GarMar15 §3). Used to represent **superpositions of stabilizer states** — i.e. arbitrary states with bounded stabilizer rank. Critical for handling non-Clifford gates (Toffoli, $T$) symbolically.

```
StabilizerFrame[{ {phase_1, generators_1}, {phase_2, generators_2}, ... }]
```

### 1.6 Clifford channel (Yashin)

Choi-state stabilizer tableau $[\mathsf U_A | \mathsf U_B | \mathsf c]$ describing arbitrary stabilizer operations — including measurements, dephasing, qubit discarding, mixed-state preparations (Yashin25 §2.3). Composition becomes "find a basis of the intersection of two vector subspaces" — Gaussian elimination, no special cases. A diagram-of-tableaux contraction (Yashin25 §3.2) reproduces the AarGot04 algorithm as one possible contraction strategy.

This is the most *general* representation in any of the papers and the natural target for a Wolfram package that wants to subsume mixed states, classical control, post-selection, and noise without a forest of special cases.

---

## 2 — What the package must compute

### 2.1 Tableau-update rules for Clifford gates (AarGot04 §3, Biswas24 §3, PatGuh26 §3.3)

For gate $U \in \{H, S, \text{CNOT}\}$ applied to qubit $a$ (or $a\to b$ for CNOT), update each tableau row $i\in\{1,\ldots,2n\}$ as:

| Gate | Update |
|---|---|
| **Hadamard $H_a$** | $r_i \mathrel\oplus= x_{ia} z_{ia}$; swap $x_{ia}\leftrightarrow z_{ia}$ |
| **Phase $S_a$** | $r_i \mathrel\oplus= x_{ia} z_{ia}$; $z_{ia}\mathrel\oplus= x_{ia}$ |
| **CNOT $a\to b$** | $r_i \mathrel\oplus= x_{ia} z_{ib}(x_{ib}\oplus z_{ia}\oplus 1)$; $x_{ib}\mathrel\oplus= x_{ia}$; $z_{ia}\mathrel\oplus= z_{ib}$ |
| **Pauli $X_a, Z_a$** | flip $r_i$ when stabilizer has corresponding anticommuting literal (FangYing23 Fact 1) |

These are *the* essential primitives. Wolfram-natural realization: keep these as pattern-rewrite rules over a `SparseArray`-backed binary matrix; vectorize across rows with `BitXor[..., ...]` (one `BitXor` per gate, $\mathcal O(n)$ per gate, no loops).

PatGuh26 §3.3 derives all of these from Karnaugh maps — useful documentation for the package's notebooks but not for the runtime.

### 2.2 The 24-element local Clifford group $\mathcal C_1$

For graph-state representation: every single-qubit Clifford acts on Pauli operators as a $3!\times \pm = 24$-fold action on $\{X,Y,Z\}$ (AndBri05 §2). Encode as a $24\times 24$ multiplication table (lookup) plus a 24-entry decomposition table (each $C_a = \sqrt{\pm iX} \cdot \sqrt{\pm iZ} \cdots$ at most 5 factors, AndBri05 Eq 7).

The package needs:
- `LocalCliffordGroup[]` — list of 24 elements as 2×2 matrices over $\mathbb Z[i]$.
- `LocalCliffordIndex[m]` — find the index 0..23.
- `LocalCliffordCompose[a, b]` — table lookup.
- `LocalCliffordDecompose[a]` — return list of $\sqrt{\pm iX}, \sqrt{\pm iZ}$ factors (AndBri05 Listing 1).

### 2.3 Local complementation

Given a vertex $a$ in graph $G$, `LocalComplement[G, a]` complements all edges among $a$'s neighbors (AndBri05 Def 1; PatGuh26 §3.2). Theorem (AndBri05 Thm 1): the resulting state differs from $|G\rangle$ by a known local unitary $U \propto \sqrt{K_G^{(a)}}$, which is captured by VOP updates per Eq (8) of AndBri05.

In Wolfram, both the graph and the VOP correction are first-class outputs. `Graph[VertexComplement[g, a]]` plus a list-rewrite of VOPs.

### 2.4 Measurement (random and deterministic outcomes)

Z-basis measurement of qubit $a$ on a tableau state (AarGot04 §3, ten lines of code):

**Random** (some $p\in\{n+1,\ldots,2n\}$ has $x_{pa}=1$):
1. `rowsum(i, p)` for every $i\ne p$ with $x_{ia}=1$.
2. Copy row $p$ into row $p-n$ (replace destabilizer).
3. Zero row $p$, set $z_{pa}=1$, set $r_p$ uniform random (or a fresh `Symbol[]`).
4. Return $r_p$.

**Deterministic** (no such $p$):
1. Zero the scratch row $2n+1$.
2. `rowsum(2n+1, i+n)` for every $i$ with $x_{ia}=1$.
3. Return $r_{2n+1}$.

The `rowsum(h, i)` subroutine multiplies stabilizer row $i$ into row $h$ and tracks the phase via the function $g(x_1,z_1,x_2,z_2)$ (AarGot04 explicit formula). Implement once.

`X` and `Y` measurements: precondition with $H_a$ or $SH_a$, measure $Z_a$, postcondition (PatGuh26 §3.4; EllEasCav08 §3.1). The package should expose `Measure[state, basis_, qubit_]` polymorphically.

**Symbolic outcomes (Wolfram-specific):** for random measurements, instead of a coin flip, allocate a fresh `Symbol[]` and let later `MeasurementOutcomes[]` invoke `RandomVariate` on the joint distribution. This is exactly the FangYing23 SymPhase trick: it amortizes circuit traversal across many samples.

### 2.5 Inner products

`StabilizerInnerProduct[ψ, φ]`: zero if the stabilizer groups have a Pauli with opposite signs; otherwise $2^{-s/2}$ where $s$ is the minimal symmetric difference of generators (AarGot04 §3 last paragraph; GarMarCro12 §3, $\mathcal O(n^3)$ algorithm).

deSilSalYin23 §3 gives more refined algorithms: amplitudes ↔ check matrix runs in $\mathcal O(Nn)$ (their "Sone→Sthree" Theorem) — exponential improvement over naive $\mathcal O(N^4/Nn^2)$. For a Wolfram package this is what makes "give me the dot product of two state vectors I built symbolically" feasible at $n=10\text{--}15$.

### 2.6 Distance / nearest neighbors

For an $n$-qubit stabilizer state, there are exactly $4(2^n-1)$ nearest-neighbor stabilizer states with $|\langle \psi|\varphi\rangle|=2^{-1/2}$ (GarMarCro12 §5). Useful as a `NearestNeighbors[]`-style method; a sanity check; and a building block for stabilizer rank algorithms.

### 2.7 Counting and enumeration

- $|\mathcal P_n| = 4^{n+1}$.
- Number of $n$-qubit stabilizer states: $N(n) = 2^n \prod_{k=0}^{n-1}(2^{n-k}+1)$ (AarGot04 Prop 1; equivalently $N(n)=2(2^n+1)N(n-1)$ with $N(1)=6$).
- $|\mathcal C_n / \mathcal P_n| = |\text{Sp}(2n, \mathbb F_2)| = 2^{n^2}\prod_{j=1}^n(4^j-1)$ (KoeSmo14 §1.1).
- $|\mathcal C_n| = 2^{n^2+2n}\prod_{j=1}^n(4^j-1)$ (KoeSmo14 Eq 2).

These are useful as `StabilizerStateCount[n]`, `CliffordGroupOrder[n]`, etc. — exact symbolic results that demonstrate the framework's "no approximation" stance.

### 2.8 Random Clifford / random stabilizer state

`RandomClifford[n]` and `RandomStabilizerState[n]`:
- KoeSmo14 §3.2 gives an $\mathcal O(n^3)$ algorithm via symplectic transvections that maps `Range[CliffordGroupOrder[n]]` → $\text{Sp}(2n,\mathbb F_2)$ bijectively. It returns the Clifford as a product of $\le 4n$ transvections, which immediately yields a circuit.
- The Pauli factor is a fresh random $4n$-bit vector ($r,s$ in KoeSmo14 §1.1).

This is a 2-design: for many use cases (randomized benchmarking, data hiding, average-case stabilizer rank decompositions, scrambling estimates) you don't need Haar-random unitaries; uniform-Clifford suffices.

### 2.9 Canonical forms

**Aaronson–Gottesman canonical form** (AarGot04 §5 Theorem 8): every unitary stabilizer circuit is equivalent to a circuit of 11 rounds H-C-P-C-P-C-H-P-C-P-C, with $\mathcal O(n^2/\log n)$ gates total. The `Cleve-Gottesman` and `Bravyi-Maslov` variants give shorter forms (5- and 4-round) at cost of less generality (deSilSalYin23 §3.4 references).

**Garcia–Markov–Cross canonical form** (GarMarCro12 §3): H-C-CZ-P-H form, derived from a stabilizer state by symplectic Gaussian elimination; useful for inner-product computation.

**Standard form $H_s$ for QEC encoders** (Got97; MonPar23 §IV): $H_s = [I_1, A_1, A_2 | B, C_1, C_2; 0, 0, 0 | D, I_2, E]$. Drives encoder-circuit synthesis directly (MonPar23 Algorithm 1).

The package should expose `CanonicalForm[circuit, "AaronsonGottesman"|"GarciaMarkov"|"BravyiMaslov"|"StandardForm"]`.

### 2.10 Synthesis from a tableau

Given a Clifford tableau $C$, generate a circuit:
- **Connectivity-free:** Gottesman's $\mathcal O(n^3)$ Gaussian-elimination synthesis (Got97 §6.4); produces $\mathcal O(n^2)$ gates.
- **Hardware-aware** (Reid24, Winderl23): *iteratively pivot* a qubit, sanitize destabilizer then stabilizer rows by single-qubit Cliffords, then remove inter-qubit interactions using **CNOT gates restricted to a Steiner tree** over the connectivity graph (Winderl23 §IV.A, Algorithm 1). Worst case $\mathcal O(q^4)$ time, $\mathcal O(q^2)$ gates, but typically much fewer than the unconstrained canonical form once routing is included.

For Wolfram, expose `CliffordSynthesize[tableau, ConnectivityGraph -> g, GateSet -> {"H","S","CNOT"}]`.

### 2.11 Pauli tracking

Given a `Clifford ∘ Measurement ∘ Clifford ∘ ...` circuit, propagate Pauli corrections without applying them on hardware (Paler14, RuhDev25). Per qubit, the Pauli "frame" is one of $\{I, X, Z, XZ\}$. The update rules are in Paler14 Tables 1–3:

| Frame at input | After CNOT (control,target) |
|---|---|
| $(s_c, s_t)$ | $(s_c \oplus_Z [s_t\text{ has Z}], s_t \oplus_X [s_c\text{ has X}])$ |

For arbitrary rotations $R_z(\pi/4), R_z(\pi/8), R_x(\pi/4)$, post-measurement frame updates are tabulated (Paler14 Table 2).

The Wolfram-natural object is `PauliFrame[{"I","X","Z","XZ", ...}]` carried alongside the circuit and updated by `PauliTrack[]`. RuhDev25 demonstrates that for MBQC, this frame *also* defines a strict partial order over measurements which determines minimum qubit count and parallelizability.

### 2.12 QEC code extraction

Given a stabilizer subgroup, return:
- $[[n,k,d]]$ parameters: `n` = number of qubits, `k = n - rank(generators)`, `d` = min weight of $N(S)\setminus S$ (Got97 §3; Got00 §4).
- Logical $\bar X, \bar Z$ operators (Got97; MonPar23 §IV.B).
- Encoder circuit (MonPar23 Algorithm 1, automatic).
- Syndrome decoder (lookup-table by error syndrome).
- Pretty-printed stabilizer table.

`QuantumErrorCorrectingCode[5]`, `QuantumErrorCorrectingCode[7,"Steane"]`, `QuantumErrorCorrectingCode[9,"Shor"]` would be natural built-ins, all derivable from their stabilizer tables (Got97 §3.3, §3.4).

---

## 3 — What is symbolic in this package

This is the section where Wolfram has its biggest comparative advantage. The other simulators in the literature push numerical performance; *none* of them — except FangYing23, Mueller26 (in C++ via SymEngine), and the embedded Mathematica subroutine in AndBri05 — even attempts symbolic computation.

### 3.1 Symbolic Pauli arithmetic

`PauliMultiply[$X_1 Z_2$, $Y_1 Z_2$]` should return $i Z_1$ symbolically, with the phase tracked exactly. The implementation is the binary symplectic XOR plus the Hamming-weight–based phase formula (Mueller26 Eq 6):
$$c = i^{(\#F_+ - \#F_-) \bmod 4}$$
with $F_+, F_-$ the Boolean truth-table outputs of single-qubit phase contributions.

Commutators: $[P, Q] = 0$ iff the symplectic inner product $\sum_i (x_i z'_i + x'_i z_i) \bmod 2 = 0$ (Got98 §2; Mueller26 §3). One XOR + popcount, $\mathcal O(n)$.

This lets `Commutator[]` and `AntiCommutator[]` work on Pauli strings of $n=10^4$ qubits in milliseconds, *symbolically*, while preserving phase — a thing no Mathematica package currently does well.

### 3.2 Symbolic phases & symbolic measurement outcomes

FangYing23 §3 (SymPhase): represent the sign vector $\vec r$ in the tableau as **bit-vectors over $\mathbb F_2^{n_s+1}$**, where $n_s$ = number of fresh symbols introduced so far. Each row's sign $r_i$ is then a polynomial $\bigoplus_j c_{ij} s_j$. Initial value $s_0 \equiv 1$, fresh variables $s_1, s_2, \ldots$ for each random measurement outcome and each Pauli error.

When a measurement is requested, *no traversal of the circuit happens*. Instead, the package collects the linear combinations $\vec m_k = (m_{k,0}, m_{k,1}, \ldots) \in \mathbb F_2^{n_s+1}$ for each measurement $k$, then realizes a sample by drawing $\vec b\in\mathbb F_2^{n_s+1}$ from the appropriate joint Bernoulli distribution and computing $\vec m_k \cdot \vec b^T$ as a single matrix multiply. For $S$ samples, this is one $n_m\times (n_s+1)$ × $(n_s+1)\times S$ matrix product.

In Wolfram, this realizes as:

```wolfram
TableauSymbolicMeasure[tableau, qubit_, "RandomSymbol" :> $rand[++$rand]]
SampleOutcomes[tableau, distribution, n_]
```

This is the natural Mathematica analog of Stim's measurement-record system, *with the distinction* that the joint distribution is stored exactly as a polynomial in $\mathbb F_2[s_1,\ldots,s_{n_s}]$ until materialized. This buys the user `CountDistinct`, `Probability`, `Expectation`, etc., over discrete random variables — Wolfram's strong suit.

### 3.3 Symbolic Clifford parameters

For variational quantum algorithms and ansatz design, the Clifford gates depend on parameters $\theta$ that are symbolic until fixed. Mueller26 PauliEngine demonstrates this for Pauli rotations $e^{-i\theta P/2}$. The Wolfram package should support:

```wolfram
ParametricPauliRotation[p_, theta_]
ParametricCliffordCircuit[gates_, parameters_]
```

with `D[expectation, theta]` returning the gradient via the parameter-shift rule (Schuld 2019), which for Pauli generators reduces to $\frac14 \langle [H,G]\rangle$ (Mueller26 §3 paragraph "Generator Gradients") — a single commutator, computed symbolically in $\mathcal O(n)$.

### 3.4 Symbolic dynamical Lie algebras

Mueller26 §3 builds the DLA $\mathfrak g = \langle i\mathcal G\rangle_\text{Lie}$ for a parameterized circuit by computing nested commutators of Pauli strings until closure. Each commutator is a single XOR + popcount; closure is detected by hash-set lookup. This scales to dimensions in the millions where dense-matrix methods would need $2^n\times 2^n$ matrices.

```wolfram
DynamicalLieAlgebra[generators_]
LieClosure[generators_, "MaxDimension" -> ∞]
StructureConstants[generators_]
```

This is the right interface for analyzing barren plateaus, Goh-style efficient simulation of poly-DLA ansätze, etc. — research a Wolfram user is plausibly doing.

### 3.5 Symbolic interconversion between representations

deSilSalYin23 gives 10 algorithms for going between {amplitudes (S₁), quadratic form triple (S₂), check matrix (S₃)} × {unitary matrix (C₁), tableau (C₂)}. Their key insight: these conversions don't require running a circuit; they exploit the algebraic structure directly.

For Wolfram, this means:
- `StabilizerStateFromAmplitudes[vec]` → `(check matrix, quadratic form triple)`.
- `StabilizerStateFromCheckMatrix[H]` → amplitude vector (closed-form sum over the affine subspace, Theorem 2.5).
- `CliffordTableauFromMatrix[U]` → stabilizer tableau, in $\mathcal O(n2^n)$ time without reading every matrix entry (deSilSalYin23 §4.3) — exponentially better than the brute-force $N^3$ method.
- `MatrixFromCliffordTableau[T]` → unitary matrix.

### 3.6 Symbolic qudits

HosDehMoo04 and Beaudrap11 generalize everything to $d$-dimensional qudits using $\mathbb Z_d$ (or $\mathbb Z_{2d}$ when $d$ even). WinPay24b Condensed encodings unifies the parity treatment via the "fake $\frac12$" $\mathbb Z_d$-module $\frac12 \mathbb Z_{d'}$.

For the Wolfram package, this means the same code path handles qubits ($d=2$) and qutrits ($d=3$) and ququarts ($d=4$) etc., parameterized on `d`. Implementation hints:
- The Pauli rotation uses $\omega = e^{2\pi i/d}$ symbolically (`Exp[2 Pi I / d]`) and $\tau = e^{i\pi(d^2+1)/d}$ (Beaudrap11 Definition).
- Smith normal form (HosDehMoo04 Appendix) is built into Wolfram via `SmithDecomposition[mat, Modulus -> d]`.
- The condensed product $\star$ (WinPay24b §3) preserves elements of order $\le d$ in the Pauli group. It is an associative operation with a phase-correction cocycle; one place to implement it.

The Wolfram package would be the only one in the literature handling all $d$ in a single uniform code path — clear positioning.

---

## 4 — What to compute for stabilizer states (the user-facing menu)

| Function | Inputs | Output | Reference |
|---|---|---|---|
| `StabilizerStateQ[ψ]` | state vector or graph | Boolean | deSilSalYin23 §3 (S₁→S₂ in $\mathcal O(N\log N)$) |
| `StabilizerCheckMatrix[ψ]` | state | $n\times(2n+1)$ binary matrix | Got97; deSilSalYin23 |
| `StabilizerGenerators[ψ]` | state | list of $n$ Pauli operators | Got97 |
| `Destabilizers[ψ]` | state | list of $n$ Paulis anticommuting w/ stab. | AarGot04 §3 |
| `StabilizerTableau[ψ]` | state | full $2n\times(2n+1)$ tableau | AarGot04 |
| `StabilizerToGraph[ψ]` | state | `(Graph, VOPs)` | AndBri05; cite Van den Nest |
| `GraphToStabilizer[g, vops]` | graph + ops | tableau | AndBri05 |
| `StabilizerToQuadraticForm[ψ]` | state | $(\mathcal A, Q, \ell)$ triple | DehMoo03 §4; deSilSalYin23 |
| `QuadraticFormToStabilizer[a, Q, ℓ]` | triple | state vector / tableau | DehMoo03; deSilSalYin23 |
| `StabilizerInnerProduct[ψ, φ]` | two states | exact value $\in \mathbb Z[i]/\sqrt 2^n$ | GarMarCro12 |
| `StabilizerDistance[ψ, φ]` | two states | symbolic | Got00 §3 (quantum Singleton) |
| `StabilizerEntanglement[ψ, partition]` | state + bipartition | exact entropy | use cluster→stabilizer + bipartite stabilizer entropy formulas |
| `StabilizerNearestNeighbors[ψ]` | state | list of $4(2^n-1)$ states | GarMarCro12 §5 |
| `LocalCliffordEquivalent[ψ, φ]` | two states | Boolean + transformation | AndBri05; cluster-state LC literature (Van den Nest) |
| `LocalComplement[g, vertex]` | graph | new graph + VOP correction | AndBri05; PatGuh26 |

---

## 5 — What to compute for Clifford operations / channels

| Function | Inputs | Output | Reference |
|---|---|---|---|
| `CliffordOperatorQ[U]` | unitary matrix | Boolean | deSilSalYin23 §4 |
| `CliffordTableau[U]` | matrix or symbolic gate sequence | $2n\times(2n+1)$ tableau | AarGot04; deSilSalYin23 |
| `CliffordMatrix[T]` | tableau | unitary matrix in $\mathcal O(N\,n)$ | deSilSalYin23 §4 |
| `CliffordCompose[T1, T2]` | two tableaux | tableau | AarGot04, DehMoo03 §2 |
| `CliffordInverse[T]` | tableau | tableau, $\mathcal O(n^2)$ | DehMoo03 §2 Theorem 3 |
| `CliffordSymplectic[T]` | tableau | symplectic matrix $\in\text{Sp}(2n,\mathbb F_2)$ | DehMoo03; KoeSmo14 |
| `CliffordCircuit[T, gateSet, connectivity]` | tableau | circuit | Reid24; Winderl23; canonical form via AarGot04 |
| `CliffordCanonicalForm[T, "AaronsonGottesman"]` | tableau | 11-round H-C-P-C-P-C-H-P-C-P-C circuit | AarGot04 §5 |
| `RandomClifford[n]` | int | tableau | KoeSmo14 |
| `IndexClifford[T]` | tableau | int $0\le i < |\mathcal C_n|$ | KoeSmo14 §3.3 (inverse) |
| `ParametricCliffordRotation[P, θ]` | Pauli + symbol | parametric symbolic operator | Mueller26 |
| `CliffordChannel[circuit]` | circuit with measurements & noises | Choi tableau $[\mathsf U_A | \mathsf U_B|\mathsf c]$ | Yashin25 |
| `ChannelCompose[Φ, Ψ]` | two channels | channel via vector-space intersection | Yashin25 §2.3 |
| `PauliTrack[circuit, frame]` | circuit + initial frame | final frame | Paler14; RuhDev25 |

---

## 6 — Specifically *symbolic* operations

| Function | Reference |
|---|---|
| `SymbolicMeasurementOutcome[tableau, qubit]` | FangYing23 §3 |
| `SubstituteOutcomes[symbolic_state, rules_]` | FangYing23 §4 |
| `SymbolicPauliCommutator[P, Q]` | Mueller26 §3 |
| `DynamicalLieAlgebra[generators]` | Mueller26 §3 |
| `LieClosure[generators]` | Mueller26 |
| `StructureConstants[basis]` | Mueller26 |
| `StabilizerEntropy[ψ, n]` | for nondegenerate codes |
| `MagicMonotone[ρ]` | stabilizer rank, robustness of magic |
| `ParametricCircuit[gates, params]` | Mueller26 |
| `GradientPauli[expectation, parameter]` | Mueller26; parameter-shift rule |
| `StabilizerRankDecomposition[ρ]` | GarMar15; Bravyi 2016 |

---

## 7 — Specific things that surprised me (worth including)

- **Yashin25's contribution is the cleanest**: every stabilizer object — pure state, mixed state, channel, measurement, post-selection, parameter-dependent operation — is a single Boolean matrix $[\mathsf U_A|\mathsf U_B|\mathsf c]$. Composition is intersection of vector spaces (Gaussian elimination). This collapses the AarGot04 / GarMar15 / FangYing23 / RuhDev25 separate special cases into one. **A Wolfram package would be wise to use this as the core type.**

- **AndBri05 implementation cost claims** (§4): "The implementation has approximately 1400 lines"; testing against CHP found no discrepancies after $4\times 10^6$ ops. The Wolfram package can plausibly reuse the same correctness oracle (cross-check against `QuantumOperator[]` direct simulation for small $n$).

- **Patil & Guha cluster-state rule book** (PatGuh26 §3, §4): explicit graph rewrites for $X, Y, Z$ measurements on cluster states. For users doing measurement-based quantum computing or photonic quantum networks, this is irreplaceable. *Implement these as `Graph` rewrites with direct visualization output.*

- **Beaudrap11's Weyl-operator trick** (§3): replacing $X^a Z^b$ with $W_{a,b} = \tau^{-ab}Z^aX^b$ removes all quadratic phase computations. The cost is having to track everything mod $D=2d$ for $d$ even, but the codebase is uniform. *This is the right form to use internally.*

- **deSilSalYin23 §4.3**: tableau-from-unitary in $\mathcal O(n2^n)$ instead of $\mathcal O(2^{3n})$ — *without reading every matrix entry*. Useful when a user constructs a Clifford symbolically and needs its tableau form.

- **MonPar23 §IV** gives a complete encoder-circuit synthesis algorithm from $H_q$ (parity-check matrix) for *any* stabilizer code — not just 5-qubit and Steane. Worth implementing as `EncoderCircuit[code]`.

- **GarMar15's stabilizer frames**: superpositions of stabilizer states. The right object for $T$-gate-rich circuits, magic-state distillation analysis, and stabilizer-rank simulation. The data structure is just a list of `(phase, generators)` pairs with operations defined to keep it compact.

- **KoeSmo14 §3.3** gives the *inverse* map: take a Clifford and produce the integer index. Useful for hashing, deduplication, randomized testing.

- **No paper covers**: efficient *symbolic* commutators between long Pauli strings with parametric coefficients on both. Mueller26 PauliEngine in C++ is closest; the Wolfram analog is straightforward but doesn't exist yet.

---

## 8 — Suggested architecture for the package

```
QuantumFramework/Stabilizer/
├── Pauli/
│   ├── PauliString.wl          (* symplectic encoding, mult, commutator *)
│   ├── PauliGroup.wl           (* enumeration, group ops *)
│   └── SymbolicPauli.wl        (* parametric coefficients via Symbol[] *)
├── Tableau/
│   ├── Tableau.wl              (* AarGot04 with destabilizers *)
│   ├── ExtendedTableau.wl      (* Yashin25 channel tableaux *)
│   ├── SymbolicTableau.wl      (* FangYing23 phase symbolization *)
│   └── TableauUpdate.wl        (* H, S, CNOT, X, Y, Z, measurement *)
├── Graph/
│   ├── GraphState.wl           (* AndBri05 *)
│   ├── LocalClifford.wl        (* 24-element table *)
│   ├── LocalComplement.wl
│   └── ClusterStateRules.wl    (* PatGuh26 *)
├── QuadraticForm/
│   ├── QuadraticFormTriple.wl  (* DehMoo03 *)
│   └── Conversions.wl          (* deSilSalYin23 *)
├── Frame/
│   ├── StabilizerFrame.wl      (* GarMar15 *)
│   └── FrameOperations.wl
├── Synthesis/
│   ├── CanonicalForm.wl        (* AarGot04 11-round *)
│   ├── ConnectivityAware.wl    (* Reid24, Winderl23 *)
│   ├── EncoderSynthesis.wl     (* MonPar23 *)
│   └── PauliTracking.wl        (* Paler14, RuhDev25 *)
├── RandomClifford/
│   ├── KoenigSmolin.wl         (* KoeSmo14 *)
│   └── SymplecticTransvection.wl
├── Codes/
│   ├── FiveQubit.wl
│   ├── Steane.wl
│   ├── Shor.wl
│   ├── CSS.wl                  (* Got00 §4 *)
│   └── Stabilizer.wl
├── Qudit/
│   ├── QuditPauli.wl           (* HosDehMoo04 *)
│   ├── WeylOperator.wl         (* Beaudrap11 *)
│   └── CondensedEncoding.wl    (* WinPay24b *)
└── Symbolic/
    ├── ParametricCircuit.wl    (* Mueller26 *)
    ├── DynamicalLieAlgebra.wl
    └── SymbolicGradient.wl
```

The `Stabilizer` package should integrate with the existing `QuantumState`, `QuantumOperator`, and `QuantumCircuitOperator` types via:

- A `"StabilizerForm"` property on `QuantumState[]` that returns a `StabilizerTableau[]`.
- A `"GraphForm"` property returning `(Graph, VOPs)`.
- A `"QuadraticForm"` property returning `(A, Q, ℓ)`.
- A `"Clifford"` test on `QuantumOperator[]` returning Boolean (using deSilSalYin23 algorithm).
- A `"CliffordTableau"` property on `QuantumOperator[]` (when the test passes).

These are all *deterministic, exact, side-effect-free* operations — Mathematica's natural style.

---

## 9 — What *not* to do

- **Don't try to outrun Stim/Qiskit on raw qubit counts.** Stim does $\sim 10^4$ qubits in C++ with SIMD and a custom column-major layout (Gid21 §3; FangYing23 §4). Wolfram cannot match that. Stop competing on that axis.
- **Don't write loops** when broadcasting / `BitXor` over `SparseArray`s suffices. The whole point of using Wolfram is functional, vectorized code.
- **Don't store global phases of stabilizer states.** GarMar15 §2: they are unobservable, and tracking them creates bugs (Got97 §4). Track them only when *relative* phases matter (i.e. inside `StabilizerFrame[]` superpositions).
- **Don't reimplement the qubit-only case.** Use HosDehMoo04 / Beaudrap11 modular arithmetic from the start, parameterized on `d`. The qubit case is `d=2`, no special code path.
- **Don't expose Karnaugh-map-derived gate rules to the end user.** They belong to internal documentation (PatGuh26 §3.3).
- **Don't try to symbolically diagonalize $2^n\times 2^n$ matrices.** The whole point of stabilizer formalism is to avoid that. Stay in the binary symplectic representation.
- **Don't separate `pure stabilizer state` and `Clifford channel` into different code paths.** Yashin25 unifies them. Use the unified Choi-tableau encoding.

---

## 10 — Direct citations driving each design choice

- **Tableau core**: AarGot04 §3 (the algorithm), with phase-symbolization extension from FangYing23 §3 and the channel generalization from Yashin25 §2.3.
- **Graph state representation**: AndBri05 §2 (graph + 24-VOP encoding), with measurement rules from EllEasCav08 §3 and PatGuh26 §3.
- **Quadratic form triple**: DehMoo03 §4 Theorem 3 (existence), HosDehMoo04 §V Theorem 1 (qudit generalization), deSilSalYin23 §3 (algorithms for the conversions).
- **Inner product**: GarMarCro12 §3 ($\mathcal O(n^3)$, with $H$-$C$-$CZ$-$P$-$H$ canonical form), with deSilSalYin23 §3 improvements.
- **Random Clifford**: KoeSmo14 §3 (subgroup algorithm via symplectic transvections, $\mathcal O(n^3)$).
- **Symbolic Pauli arithmetic**: Mueller26 §2 (phase tracking via popcount), §3 (commutators / Lie closure / DLA / parameter-shift gradients).
- **Hardware-aware synthesis**: Winderl23 §IV (Steiner-tree-based pivot), Reid24 §III (related architecture-aware approach).
- **Pauli tracking**: Paler14 §V (algorithm), RuhDev25 §3 (modern library + MBQC scheduling).
- **Encoder synthesis**: MonPar23 §IV (algorithm + standard form $H_s$), Got97 §6 (foundational).
- **Cluster-state rule book**: PatGuh26 §3, §4 (graph rewrites for measurements + fusions).
- **Stabilizer frames**: GarMar15 §3 (data structure for non-Clifford support).
- **Qudit unification**: Beaudrap11 §2–3 (Weyl-operator linearization), HosDehMoo04 §III–IV (qudit Clifford), WinPay24b (condensed encoding for even $d$).
- **Wigner-function equivalence**: KocHuaLov17 §3 (tableau ⇔ discrete Wigner function for odd $d$).
- **Lambda calculus / type theory**: PayWin24a (projective Cliffords as a typed programming abstraction; useful as a *correctness* model, not for runtime).

---

## 11 — Final priorities for v1

1. **Pauli string + symplectic + symbolic phase** (Mueller26 + FangYing23). Without this, nothing works.
2. **Tableau with destabilizers + AarGot04 update rules** (AarGot04). The workhorse.
3. **Graph-state + 24-VOP** (AndBri05). The compact storage for QEC / entanglement-purification simulations.
4. **Quadratic form triple** (DehMoo03 / deSilSalYin23) and the 6 interconversion algorithms among {amplitudes, quadratic form, check matrix}.
5. **Clifford channels via Choi tableau** (Yashin25). The unifier.
6. **Random Clifford** (KoeSmo14). One-liner once symplectic is in.
7. **Inner product** (GarMarCro12 / deSilSalYin23).
8. **Stabilizer frames** (GarMar15) — required for non-Clifford.
9. **Hardware-aware synthesis** (Winderl23) — for QF integration with hardware backends.
10. **Pauli tracking + MBQC scheduling** (RuhDev25) — for QF measurement-based pipeline.

Everything else — the qudit case (Beaudrap11, HosDehMoo04, WinPay24b), the symbolic gradient / DLA tools (Mueller26 §3), the cluster-state rule book (PatGuh26), the QEC code library (Got00, MonPar23) — sits on top of these primitives without modification.
