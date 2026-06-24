---
Title: 'Symbolic path sums for Clifford+T in the Wolfram Language: a pedagogical implementation, with honest benchmarks and provenance'
Author: investigation in the Wolfram QuantumFramework
Description: 'A computation-first, arXiv-style note. It introduces the sum-over-paths (path-sum) method for classically computing Clifford+T amplitudes, traces its provenance honestly (Feynman; Dawson et al. 2004; Montanaro 2016; Koh-Penney-Spekkens 2017; Amy 2018; Amy-Stinchcombe 2024), explains a Wolfram Language implementation, and benchmarks it. The contribution claimed is an implementation, not a method: an exact, computer-algebra-native, symbolic-parameter, QuantumFramework-integrated realization. Every literature claim was checked against the fetched TeX sources; every amplitude is exact against QuantumFramework dense simulation.'
---

# Symbolic path sums for Clifford+T in the Wolfram Language

> **Revised verdict (see `../WL-Symbolic-Novelty-Reexamination.md`).** A later adversarial re-examination reaffirms this note's central claim, that there is no new method here that breaks the $\#\mathrm{P}$ wall or enlarges the tractable class, but it re-weights two items the contribution table in Section 3 files as "implementation": exact **symbolic-parameter path sums** are a capability that *no competing simulator has* (feynman, quizx, de Colnet, Koh-Penney-Spekkens, and the stabilizer-rank tools are all numeric coefficient-ring engines), and the **Gröbner/variety point count** is a distinct, unbuilt realization of a problem Amy-Stinchcombe explicitly leave open, with its own tractability parameter. See that note's Sections 3(a) and 3(b).

## Abstract

We give a Wolfram Language realization of the *path-sum* (Feynman sum-over-paths) method for exactly computing amplitudes of Clifford+$T$ quantum circuits. An amplitude is carried as a Gauss sum $\sum_v \omega^{\varphi(v)}$ with $\omega = e^{i\pi/4}$, the path variables $v$ ranging over $\mathbb{F}_2$ and $\varphi$ a degree-$\le 3$ phase polynomial valued in $\mathbb{Z}_8$; a rewrite engine folds the Clifford structure and a brute-force sum completes the magic residual. **None of this is new.** The representation, the rewrite rules, the degree bound, and the residual-completion idea are due to Feynman, Dawson et al. (2004), Montanaro (2016), and above all Amy (2018), who also ships an open-source Haskell toolkit (`feynman`); confluence of a related system is Amy-Stinchcombe (2024). What this note contributes is an *implementation and an exposition*: a computer-algebra-native engine with exact cyclotomic arithmetic, symbolic gate angles, and live integration into the QuantumFramework object model, together with an honest account of where carrying the amplitude symbolically helps and where it does not. We benchmark within the Wolfram Language (against dense `Method -> "Schrodinger"` and the `StabilizerFrame`) and name the external state of the art (quizx, the rank-width sum-of-powers simulator of de Colnet et al. 2026, LIMDD, stabilizer rank) that this implementation does not out-run. Every literature claim below was verified against the fetched arXiv TeX sources; every amplitude is exact, equal to QuantumFramework dense simulation under `RootReduce`.

> **What this is, and is not.** This is an *exact, symbolic, structure-exploiting* amplitude calculator built inside a computer-algebra system. It is *not* a performance-competitive general simulator: on unstructured or narrow-deep magic it completes its residual by a brute-force Gauss sum and is slower than dedicated tools. Its niche is exactness in $\mathbb{Z}[\omega, 1/\sqrt2]$, amplitudes as closed-form *functions of gate angles*, and the path sum as one representation among many in QuantumFramework. The reusable code is `pathsum.wl`, `harness.wl`, `elimination.wl`, `diagonal.wl`, `rankcap.wl`; benchmarks in `benchmarks/`.

---

## 1. The idea, and where it comes from

The question this note answers is narrow and was sharpened by an adversarial review: *given that symbolic algebra can carry randomness (SymPhase) and magic (a phase polynomial) instead of branching, is there a usable, exact, computer-algebra realization of that for Clifford+$T$, and what does it actually buy?* Answering it honestly requires first stating, without inflation, how much of the idea already exists. We fetched and read the TeX sources of the anchor papers (folder `pathsum-literature/`); the provenance below is anchored to them, not recalled.

### 1.1 The path sum

A quantum amplitude is a Feynman sum over computational paths. For a circuit $C$ on $n$ qubits built from $H$ and diagonal gates, inserting a resolution of the identity at each Hadamard gives

$$
\langle y \lvert C \rvert 0\cdots 0\rangle \;=\; \frac{1}{\sqrt{2^{h}}}\sum_{v\,:\,E(v)=y}\omega^{\varphi(v)},
\qquad \omega = e^{i\pi/4},
$$

where $h$ is the number of Hadamards, each contributing one path variable $v_i \in \mathbb{F}_2$ and one factor $1/\sqrt2$; $E(v)$ is the (multilinear, $\mathbb{F}_2$) vector of output-qubit values; and $\varphi$ is a phase polynomial valued in $\mathbb{Z}_8$. For Clifford+$T$ on $n$ qubits, $\varphi$ has degree at most $3$ (a $T$ contributes a linear term, $CZ$ a quadratic term, $CCZ$ a cubic term). Evaluating the sum exactly is $\#\mathrm{P}$-hard in general; the quadratic (Clifford) part is polynomial.

### 1.2 Provenance (verified against the sources)

- **Feynman** introduced the sum-over-paths formulation that underlies all of this.
- **Dawson, Haselgrove, Hines, Mortimer, Nielsen, Osborne (2004; QIC 2005), arXiv:quant-ph/0408129** established the core simulation idea: a circuit amplitude equals a balanced-gate Feynman sum-over-paths, which is a signed count of solutions of a degree-$\le 3$ polynomial system over $\mathbb{Z}_2$ (their Eq. for $\langle b\lvert U\rvert a\rangle$, Sec. III). They already write the $H,T,\mathrm{CNOT}$ case as $\sum_y e^{i\pi\varphi(y)/4}$ with $\omega=e^{i\pi/4}$, $\varphi$ "of order at most two" in mixed $\mathbb{Z}_2/\mathbb{Z}_8$ arithmetic (their Sec. VI), and they sketch the stabilizer version with $T$ as the residual non-Clifford generator. They state the general case is $\#\mathrm{P}$-hard.
- **Montanaro (2016), arXiv:1607.08473** cleaned the qubit case to a *single* degree-$\le 3$ polynomial $f_C$ over $\mathbb{F}_2$ with amplitude $\propto \mathrm{gap}(f_C)=\sum_x(-1)^{f_C(x)}$, and pinned the dichotomy: degree-2 gap is computable in $O(n^3)$ ("Gottesman-Knill in polynomial form"), degree-3 is $\#\mathrm{P}$-hard. He notes this is a reframing of Dawson et al.
- **Koh, Penney, Spekkens (2017), arXiv:1702.03316** closed the quadratic case for *quopit* ($p$ odd prime) Clifford circuits by evaluating the resulting Weil/Gauss sums in closed form after diagonalizing the quadratic form over $\mathbb{F}_p$, with an explicit polynomial-time algorithm. They note the construction *fails at characteristic 2* (the qubit $\mathbb{Z}_8$ case) and defer it.
- **Amy (2018), arXiv:1805.06908** is the direct antecedent of the implementation in this note. His "path-sum" object (Def. 2.1) is exactly ours: amplitude $\frac{1}{\sqrt{2^m}}\sum_{y}e^{2\pi i P(x,y)}\lvert f(x,y)\rangle$, with $P$ multilinear over $\mathbb{F}_2$ in the dyadic ring, the Clifford+$T$ phase of degree $\le 3$ over $\mathbb{Z}_8$, one path variable and one $1/\sqrt2$ per Hadamard. He defines the rewrite rules we reuse: **[Elim]** (drop a free variable, factor 2), **[$\omega$]** (the $\mathbb{Z}_8$/$S$ case), **[HH]** (Hadamard elimination with a linear substitution), and notes that a complete reducer would brute-force-expand the residual variables. He provides an open-source Haskell toolkit, `feynman` (github.com/meamy/feynman), used there to *verify* (and effectively simulate) circuits up to 96 qubits and tens of thousands of gates, and to handle the hidden-shift algorithm in seconds. His system is terminating but, he states, *not confluent*.
- **Amy, Stinchcombe (2024), arXiv:2408.02778** make a related rewrite system **confluent** (modulo a degree-preserving equivalence) by restricting the [HH] substitution to a single variable, and prove polynomial-time simulation of the hidden-shift family ($3n$ rule applications). Crucially, their confluence result is for the **$\mathbb{F}_2$ sign fragment** (Toffoli-Hadamard, phase $(-1)^P$), not the $\mathbb{Z}_8$ Clifford+$T$ phase. They list the "sum over varieties via Gröbner bases / point counting" generalization as explicit future work.

The honest conclusion: the path-sum object, the degree bound, the [HH]/[$\omega$]/[Elim] rewriting, the residual-completion idea, and the hidden-shift application are all prior art, principally Amy's. A Wolfram Language version re-derives none of this mathematics.

### 1.3 Three neighbors that frame the method

- **Gottesman-Knill, exactly.** The Clifford fragment is not something the path sum must brute-force. Guan and Regan (2019, arXiv:1904.00101) show the Clifford amplitude is a $\mathbb{Z}_4$ Gauss sum whose magnitude is $2^{-r/2}$, with $r$ the $\mathbb{F}_2$ rank of a quadratic form, computable in matrix-multiplication time; Aaronson-Gottesman (2004, arXiv:quant-ph/0406196) is the tableau baseline our "fold the Clifford structure" recovers. Adding $T$ pushes the phase to $\mathbb{Z}_8$ and breaks the closed form, leaving the magic as the residual.
- **SymPhase is the sign analogue.** Fang and Ying (2023, arXiv:2311.03906) symbolize the *sign* column of the stabilizer tableau as a degree-1 $\mathbb{F}_2$ polynomial to carry **measurement and Pauli-fault** branching (their Facts 1 and 2: Paulis touch only the sign, and the control flow ignores the sign). There is no $T$ gate anywhere in SymPhase. The path sum is its magic analogue: a degree-$\le 3$ $\mathbb{Z}_8$ *phase* polynomial in place of a degree-1 $\mathbb{Z}_2$ *sign* polynomial.
- **The magic residual has a known cost.** Bravyi-Gosset (2016, arXiv:1601.07601) and Bravyi-Browne-Calpin-Campbell-Gosset-Howard (2019, arXiv:1808.00128) price the magic as a stabilizer rank $\chi \sim 2^{\beta t}$ ($\beta \le 0.47$ exact, $\approx 0.23$ approximate; in general $\prod_j \xi(V_j)$ via the stabilizer extent). The path sum carries the same exponential-in-$t$ residual as one brute-forced phase polynomial rather than a linear combination of stabilizer states.

---

## 2. The Wolfram Language package

The package is small and self-loads QuantumFramework. The amplitude carrier is one association.

### 2.1 Representation and gates (`pathsum.wl`)

Each qubit value is an $\mathbb{F}_2$ multilinear polynomial in path variables `x[i]`; the phase is a $\mathbb{Z}_8$ pseudo-Boolean polynomial; a scale field carries the exact $\mathbb{Z}[\omega,1/\sqrt2]$ prefactor. Idempotency ($x_i^2=x_i$) and the modular reductions are native one-liners:

```wl
PathSum`idem[p_] := Expand[p] /. Power[x[i_], _] :> x[i];
PathSum`f2[p_]   := PolynomialMod[PathSum`idem[p], 2];     (* F_2 qubit values *)
PathSum`z8[p_]   := PolynomialMod[PathSum`idem[p], 8];     (* Z_8 phase *)
```

Each gate is the gate matrix read off on basis states: a Hadamard introduces a fresh `x[i]` and adds $4\,e_q v$ to the phase; $T$ adds $e_q$; $CZ$ adds $4\,e_a e_b$; $CCZ$ adds $4\,e_a e_b e_c$; $CNOT$ does $e_b := e_a \oplus e_b$. The amplitude is the brute-force Gauss sum, kept permanently as the exact oracle and fallback:

```wl
PathSum`amp[ps_, n_, y_] := (ps["scale"]/Sqrt[2]^ps["h"]) Total[
    (PathSum`omega^Mod[ps["phase"] /. Thread[ps["vars"] -> #], 8]) Boole[
        And @@ Table[Mod[ps["val"][q] /. Thread[ps["vars"] -> #], 2] == y[[q]], {q, n}]] & /@
    Tuples[{0, 1}, Length[ps["vars"]]]];
```

### 2.2 The rewrite engine (`elimination.wl`)

The engine sums out a path variable $v$ via the exact identity $\sum_{v\in\{0,1\}}\omega^{\varphi_0 + v C} = \omega^{\varphi_0}(1+\omega^{C})$, with $C = \mathrm{Coefficient}[\varphi, v]$ the cofactor. It factors cleanly in three degree-checked cases, each a proven identity: a constant cofactor ($1+\omega^c$); the $\mathbb{Z}_2$ **[HH]** sign ($C=4Q$, $Q$ a linear Boolean form, giving $2[Q=0]$, eliminating $v$ and imposing $Q=0$); and the $\mathbb{Z}_4$/$S$ case ($C\equiv 2 \bmod 4$, giving $\sqrt2\,\omega^{1-2L}$). Anything else (an odd nonconstant cofactor, i.e. genuine magic; a nonlinear constraint; a degree-violating step) is left as a residual that `amp` completes exactly. Correctness therefore never depends on the rewriter's completeness, and termination is guaranteed because every applied rule removes a variable.

A subtle implementation point that the verification caught: when imposing a linear constraint $u = \mathrm{rhs}$, the substitution into the *phase* must use the integer value `boolToInt[rhs]`, not the bare $\mathbb{F}_2$ form, because the phase treats variables as integer $0/1$ values; this is invisible for pure-$\mathbb{Z}_2$ phases and only bites once $S$ gates inject $\mathbb{Z}_4$ terms.

### 2.3 Three extras

- **`diagonal.wl`** compresses the diagonal (CNOT-dihedral) fragment on $\lvert+\cdots+\rangle$: with no internal Hadamard the output map $E$ is a linear $\mathbb{F}_2$ bijection, so $E(v)=y$ has a unique solution (found by `LinearSolve[M, y, Modulus -> 2]`) and the amplitude is one closed-form term $2^{-n/2}\omega^{\psi(y)}$, a single degree-$\le 3$ phase polynomial, with no enumeration.
- **`rankcap.wl`** caps the eager `StabilizerFrame` over-count: it materializes the component state vectors exactly and uses `MatrixRank` over $\mathbb{Z}[\omega]$ to collapse $2^{\#T}$ components to the span dimension ($\le 2^n$). This is a *span-dimension* witness, not the NP-hard optimal stabilizer rank.
- An **`Inactive[Sum]`** carrier renders the residual as an inert sum on a `HoldAll`, `NonThreadable` head that `Activate` collapses to the amplitude, the computer-algebra realization of "carry the sum, do not enumerate it."

The native Wolfram primitives each step uses (`PolynomialMod`, `Coefficient`, `Collect`, `LinearSolve`/`MatrixRank` with `Modulus -> 2`, `RootReduce`, `Cyclotomic`, `Inactive`/`Activate`, the `Hold*` family, `FiniteFieldElementTrace` for the $GF(2^k)$ character view) are catalogued in `../wl-symbolic-toolkit-docs/WL-Symbolic-Toolkit-for-PathSums.md`.

---

## 3. What is, and is not, the contribution

| Ingredient | Status | Source |
|---|---|---|
| Amplitude as a Gauss sum $\sum_v\omega^{\varphi(v)}$, $v$ over $\mathbb{F}_2$, $\varphi$ degree-$\le3$ over $\mathbb{Z}_8$ | prior art | Dawson 2004; Amy 2018 |
| [HH] / [$\omega$] / [Elim] rewrite rules, degree bound, poly-time reducer | prior art | Amy 2018 |
| Brute-force completion of the residual | prior art | Amy 2018 |
| Confluent rewriting (for the $\mathbb{F}_2$ sign fragment) | prior art | Amy-Stinchcombe 2024 |
| Hidden-shift handled in path-sum form | prior art | Amy 2018, 2024 |
| Clifford amplitude closed-form via $\mathbb{F}_2$ rank | prior art | Guan-Regan 2019 |
| Magic-residual cost $\sim 2^{\beta t}$ | prior art | Bravyi-Gosset 2016; Bravyi et al. 2019 |
| **Exact $\mathbb{Z}[\omega,1/\sqrt2]$ arithmetic, no integer-overflow** | implementation | this note (Amy reports 31-bit overflow) |
| **Symbolic gate angles: amplitudes as closed-form functions of $\theta$** | capability (no competing simulator has it) | this note; re-exam Section 3(b) |
| **Gröbner/variety point count over $\mathbb{F}_2$** (distinct tractability axis: ideal regularity, not Amy's degree bound or de Colnet's rank-width) | open realization of a named problem | re-exam Section 3(a); idea Dawson 2004, future-work Amy-Stinchcombe 2024 |
| **Native integration with QuantumFramework** (frame rank-cap, dense oracle, one representation among many) | implementation | this note |
| **Pedagogical exposition + a verified differential-test harness** | implementation | this note |

The single algorithmic gap worth noting, not claiming as a result: confluence has been *proven* only for the $\mathbb{F}_2$ Toffoli-Hadamard fragment (Amy-Stinchcombe), whereas the engine here rewrites the $\mathbb{Z}_8$ Clifford+$T$ phase, for which we observe confluence empirically (the completed amplitude is invariant under shuffled rule order) but do not prove it.

---

## 4. Benchmarks

All timings are medians of repeated runs with `ClearSystemCache[]`, on an Apple M-series laptop (12 cores), Wolfram Language 15. Amplitudes are exact (`RootReduce`), never floating point. The Wolfram-Language baselines are dense `Method -> "Schrodinger"` and the `StabilizerFrame`. **Caveat for fairness:** the dense single-amplitude time includes QuantumFramework's full circuit-construction pipeline; the meaningful signal is the *scaling*, not the absolute constant, and we report it as such. External state-of-the-art tools are named in Section 5 as reference points, not run here.

### 4.1 A genuine win: the diagonal fragment (capability + scaling)

For the CNOT-dihedral fragment, the closed-form `diagAmp` (one `LinearSolve` plus one phase evaluation) computes a single amplitude in milliseconds, flat in $n$, where dense grows:

| $n$ | gates | `diagAmp` (s) | QF dense, one amplitude (s) |
|---|---|---|---|
| 4 | 12 | 0.0022 | 0.59 |
| 6 | 18 | 0.0033 | 0.93 |
| 8 | 24 | 0.0047 | 1.34 |
| 10 | 30 | 0.0059 | 1.89 |
| 12 | 36 | 0.0074 | 3.00 |
| 14 | 42 | 0.0088 | 7.44 |

This is the $\mathbb{Z}_8$ extension of the Clifford closed forms of Koh-Penney-Spekkens and Guan-Regan; it is a real win on its fragment, and it is exact.

### 4.2 Compression: Clifford structure telescopes to a closed form

On a Clifford $H/S/CZ$ ladder, the engine eliminates every internal path variable (internal residual $\to 0$) in a number of steps linear in depth, leaving only the output variables:

| depth | vars before | vars after | internal residual | steps | reduce (s) |
|---|---|---|---|---|---|
| 1 | 12 | 6 | 0 | 4 | 0.004 |
| 2 | 18 | 6 | 0 | 8 | 0.009 |
| 3 | 24 | 6 | 0 | 11 | 0.015 |
| 4 | 30 | 6 | 0 | 15 | 0.025 |
| 6 | 42 | 6 | 0 | 22 | 0.054 |
| 8 | 54 | 0 | 0 | 32 | 0.100 |

### 4.3 Capability QF cannot match: symbolic gate angles

Carrying a symbolic angle, the engine produces a closed-form expectation once, $\langle Z\rangle(\theta) = \cos(\theta/2)$ for $H\,P(\theta)\,H$, in $\sim 0.009$ s, where re-running dense for 17 angle values costs $\sim 1.95$ s and yields no closed form. This is a capability difference, not merely a speed ratio: the dense path has no symbolic-$\theta$ amplitude at all.

### 4.4 The representational over-count

The eager `StabilizerFrame` stores $2^{\#T}$ components for a diagonal Clifford+$T$ circuit that the compressor represents as one rank-1 phase polynomial:

| $\#T$ | `StabilizerFrame` components | compressor rank |
|---|---|---|
| 1 | 2 | 1 |
| 2 | 4 | 1 |
| 3 | 8 | 1 |
| 4 | 16 | 1 |
| 5 | 32 | 1 |
| 6 | 64 | 1 |

### 4.5 The honest loss: narrow-but-deep magic

The path sum costs $2^{\text{residual}}$, the number of *paths* it must brute-force, which can exceed $2^{n}$, the number of basis states. For $H(TH)^k$ on one qubit the residual is $k$ magic variables, so the engine brute-forces $2^k$ where dense handles a single qubit trivially. It wins for small $k$ and loses decisively past the crossover near $k\approx 16$:

| $k=\#T$ | internal residual | `engineAmp` (s) | QF dense (s) |
|---|---|---|---|
| 4 | 4 | 0.002 | 0.37 |
| 8 | 8 | 0.013 | 0.54 |
| 12 | 12 | 0.25 | 0.78 |
| 16 | 16 | 5.55 | 1.09 |
| 18 | 18 | 26.7 | 1.19 |

This is the central honest statement: when the number of magic paths exceeds the dimension of the state space, the path sum is the wrong tool, and a dense or stabilizer-rank method wins.

---

## 5. Where this stands relative to the state of the art

This implementation does **not** beat the performance-oriented exact simulators; they out-scale a brute-force-completed path sum along different axes, and an honest note must cite them as the reference frame.

- **The same object, evaluated faster: de Colnet, Geerts, Hai, Laarman, Lee, Perez (2026), arXiv:2605.29944.** Their "quadratic sums-of-powers" simulator evaluates the *identical* Clifford+$T$ path sum (citing Dawson) by a fixed-parameter dynamic program in the **rank-width** $k$ of the path-variable interaction graph: $2^{O(k)}\mathrm{poly}(n)$, where the engine here pays $2^{\lvert V\rvert}$. The $T$-count does not enter their exponent. This is the most direct competitor and the central reason not to overclaim.
- **Aggressive $T$-reduction: Kissinger and van de Wetering (2021), arXiv:2109.01076.** ZX-simplified stabilizer decompositions, $2^{\alpha t}$ with $\alpha \approx 0.468$ and 6 to 15 orders of magnitude empirical term reduction, exact in $\mathbb{D}[e^{i\pi/4}]$, reaching 50-100 qubits and $>1000$ $T$ gates in the Rust tool `quizx`.
- **Succinct representation: LIMDD, Vinkhuijzen, Coopmans, Elkouss, Dunjko, Laarman (2023), arXiv:2108.00931.** Poly-size for stabilizer, cluster, coset, and $W$-state families, with proven exponential separations from QMDD, MPS, and Clifford+$T$ stabilizer decompositions.
- **The residual yardstick: Bravyi-Gosset and Bravyi et al.** Stabilizer-rank simulation at $2^{0.47 t}$ (exact) and $2^{0.23 t}$ (approximate), a rigorous randomized estimate where the path sum here is exact but brute-forces the residual.

The honest niche of the Wolfram Language version is therefore not speed. It is (1) exactness with no floating point and no integer overflow, (2) genuinely symbolic gate parameters, and (3) being one representation inside a computer-algebra and QuantumFramework environment, with a `Modulus -> 2` Gröbner/variety path (the generalization Amy-Stinchcombe leave open) sitting one function call away.

---

## 6. Verification

Every package file and benchmark was checked by a fresh-context adversarial reviewer that reproduces each claim in a kernel (no `Quiet`, no `Print`). Across the build: the engine amplitude equals the un-reduced brute-force amplitude and the QuantumFramework dense amplitude exactly (`RootReduce`) on the worked battery, on 30+ random circuits, and on circuits the reviewer chose; the degree-$\le 3$ invariant holds after every rule; the completed amplitude is invariant under shuffled rule order; and the diagonal compressor matches dense on randomized circuits up to $n=6$. The literature claims in Sections 1, 3, and 5 were each checked against the fetched arXiv TeX (folder `pathsum-literature/`).

---

## 7. Conclusion and open directions

We have implemented, exactly and symbolically, a method that is not ours: the Clifford+$T$ path sum of Dawson, Montanaro, Koh-Penney-Spekkens, and especially Amy. The value of doing so in the Wolfram Language is exactness, symbolic parameters, and integration, not performance; the benchmarks show real wins on the diagonal fragment, on Clifford compression, and on parametric capability, a representational win over the eager `StabilizerFrame`, and a clean loss on narrow-deep magic where the path count exceeds the state-space dimension. The open directions are all known: prove (or refute) confluence for the $\mathbb{Z}_8$ Clifford+$T$ rewrite system (Amy-Stinchcombe proved it only for the $\mathbb{F}_2$ fragment); evaluate the Clifford residual in closed form via the Guan-Regan rank rather than brute force; realize the variety/Gröbner-basis generalization with `GroebnerBasis[..., Modulus -> 2]`; and, to be competitive, adopt the rank-width fixed-parameter contraction of de Colnet et al. in place of the brute-force completion. A later adversarial re-examination (`../WL-Symbolic-Novelty-Reexamination.md`) sharpens two of these: the symbolic-parameter capability is not merely an implementation detail but a capability *no cited competitor has*, and the variety/Gröbner-basis route is a *distinct tractability axis*, governed by the Gröbner-basis complexity of the Boolean ideal and separate from Amy's degree bound and de Colnet's rank-width, whose concrete realization remains unbuilt, is bounded to the $\mathbb{F}_2$ real-amplitude fragment, and by the originators' own words is no asymptotic free lunch.

## References

- R. P. Feynman, space-time approach to non-relativistic quantum mechanics, Rev. Mod. Phys. 20, 367 (1948).
- C. M. Dawson, H. L. Haselgrove, A. P. Hines, D. Mortimer, M. A. Nielsen, T. J. Osborne, *Quantum computing and polynomial equations over the finite field $\mathbb{Z}_2$*, QIC 5(2):102-112 (2005), arXiv:quant-ph/0408129.
- A. Montanaro, *Quantum circuits and low-degree polynomials over $\mathbb{F}_2$*, arXiv:1607.08473 (2016).
- D. E. Koh, M. D. Penney, R. W. Spekkens, *Computing quopit Clifford circuit amplitudes by the sum-over-paths technique*, QIC 17(13&14):1081-1095 (2017), arXiv:1702.03316.
- M. Amy, *Towards Large-scale Functional Verification of Universal Quantum Circuits*, QPL 2018, arXiv:1805.06908. Software: `feynman`, github.com/meamy/feynman.
- M. Amy, L. S. Stinchcombe, *Polynomial-Time Classical Simulation of Hidden Shift Circuits via Confluent Rewriting of Symbolic Sums*, Quantum (2025), arXiv:2408.02778.
- S. Aaronson, D. Gottesman, *Improved Simulation of Stabilizer Circuits*, Phys. Rev. A 70, 052328 (2004), arXiv:quant-ph/0406196.
- C. Guan, K. W. Regan, *Stabilizer Circuits, Quadratic Forms, and Computing Matrix Rank*, arXiv:1904.00101 (2019).
- W. Fang, M. Ying, *SymPhase: Phase Symbolization for Fast Simulation of Stabilizer Circuits*, arXiv:2311.03906 (2023).
- S. Bravyi, D. Gosset, *Improved classical simulation of quantum circuits dominated by Clifford gates*, Phys. Rev. Lett. 116, 250501 (2016), arXiv:1601.07601.
- S. Bravyi, D. Browne, P. Calpin, E. Campbell, D. Gosset, M. Howard, *Simulation of quantum circuits by low-rank stabilizer decompositions*, Quantum 3, 181 (2019), arXiv:1808.00128.
- A. Kissinger, J. van de Wetering, *Simulating quantum circuits with ZX-calculus reduced stabiliser decompositions*, arXiv:2109.01076 (2021). Software: `quizx`.
- L. Vinkhuijzen, T. Coopmans, D. Elkouss, V. Dunjko, A. Laarman, *LIMDD: A Decision Diagram for Simulation of Quantum Computing Including Stabilizer States*, Quantum (2023), arXiv:2108.00931.
- A. de Colnet, F. Geerts, R. Hai, A. Laarman, J. H. Lee, G. A. Perez, *Quadratic Sums-of-Powers for Fixed-Parameter Tractable Quantum-Circuit Simulation*, arXiv:2605.29944 (2026).

*Fetched TeX sources for all anchor papers are in `../pathsum-literature/`. Package and benchmark scripts in this folder and `benchmarks/`, verified against `Wolfram/QuantumFramework 2.0.0`.*
