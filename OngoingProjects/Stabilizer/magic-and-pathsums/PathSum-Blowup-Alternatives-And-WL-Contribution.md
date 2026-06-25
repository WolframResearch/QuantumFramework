---
Title: 'The T-gate blow-up: alternatives to the symbolic path sum, and the single most original contribution the Wolfram Language can prove'
Author: investigation in the Wolfram QuantumFramework
Description: 'A computation-first synthesis that answers two questions about the path-sum line of work. (1) Is the symbolic path sum the only solution to the T-gate blow-up (stabilizer rank growing as 2^t)? Answer: no. It is one structure-exploiting representation among several (eager frames; stabilizer rank / extent / robustness decompositions; LIMDD decision diagrams; ZX-reduced stabilizer decompositions; tensor networks on the wrong axis; de Colnet rank-width sums-of-powers), each exploiting a different axis, none removing the exponent in magic, which universality forbids. (2) Given the path sum, what is the single most original contribution the Wolfram Language can PROVE and VERIFY? Answer: exact path sums with FREE-SYMBOL gate parameters, producing closed-form amplitude, expectation, and gradient FUNCTIONS over a whole parametrized circuit family in one pass, a capability no competing simulator has. Demonstrated and kernel-verified, with the Groebner/variety point count as a bounded companion. Every literature claim is anchored to the fetched arXiv TeX; every amplitude is exact against QuantumFramework dense simulation.'
---

# The T-gate blow-up: alternatives to the path sum, and the most original thing WL can prove

## 0. What this note settles, up front

Two questions, answered with computation, building on (not re-deriving) the prior verdict in [`WL-Symbolic-Novelty-Reexamination.md`](WL-Symbolic-Novelty-Reexamination.md) and [`pathsum-build/PathSum-WL-Note.md`](pathsum-build/PathSum-WL-Note.md).

**Question 1: is the symbolic path sum the only solution to the $T$-gate blow-up?** No. The blow-up is the cost of *magic*, and every exact method pays an exponential in the magic somewhere, because $\mathrm{Clifford}+T$ is universal. What differs is the *axis* along which each method finds structure to exploit. The path sum exploits the algebraic structure of a low-degree phase polynomial. The stabilizer-rank / extent / robustness decompositions shrink the *base* of the exponential along a magic monotone. LIMDD decision diagrams quotient out local-Clifford-redundant magic. ZX-reduced stabilizer decompositions cancel non-Clifford generators graphically. Tensor networks compress *entanglement*, which is the wrong axis (magic is independent of bond dimension). The de Colnet sums-of-powers simulator evaluates the *identical* path-sum object by a fixed-parameter program in the *rank-width* of the path-variable interaction graph, a different tractability parameter that strictly out-scales a brute-force-completed path sum. The path sum is one representation among these, with a clear win/lose regime, not a unique or dominant solution.

**Question 2: what is the single most original contribution the Wolfram Language can prove and verify?** Carry the gate angles as *free symbols* through the path sum and return the amplitude, expectation, and even the gradient as exact closed-form *functions* of those angles, computed in one symbolic pass over a whole parametrized circuit family. No competing simulator does this: they are all numeric coefficient-ring engines. This is the most original-yet-defensible item because it is a *capability difference* (not a speed ratio), it is an exact algebraic identity (provable), and it is machine-checkable against dense simulation at every angle. We demonstrate it below and verify it in a clean kernel: $\langle Z\rangle(\theta)=\cos\theta$ and its gradient $-\sin\theta$ in one pass, a two-angle two-qubit $\langle Z_1 Z_2\rangle(\alpha,\beta)=\sin\alpha\,\sin\beta$ exact on a grid, and a free angle composed with a fixed $T$. The Groebner/variety point count is a genuine but *bounded* companion (a distinct tractability axis, the first concrete realization of an Amy-Stinchcombe open problem, but limited to the $\mathbb{F}_2$ real-amplitude fragment and, by the originators' own words, no asymptotic free lunch).

All experiments run against the working-tree paclet,

```wl
PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"];
```

paclet `Wolfram/QuantumFramework 2.0.0`. Every amplitude below is exact, equal to dense `Method -> "Schrodinger"` under `RootReduce` in $\mathbb{Z}[\omega, 1/\sqrt2]$, $\omega = e^{i\pi/4}$. The verification script is [`pathsum-blowup-alternatives-demos.wls`](pathsum-blowup-alternatives-demos.wls) in this folder; fetched arXiv TeX sources are in [`../literature/pathsum-literature/`](../literature/pathsum-literature/).

---

## 1. The blow-up, stated precisely, and why it is the cost of magic

A $T$ gate acting on a stabilizer state is a coherent superposition of two stabilizer states,
$$
T_q\lvert s\rangle \;=\; \frac{1+e^{i\pi/4}}{2}\,\lvert s\rangle \;+\; \frac{1-e^{i\pi/4}}{2}\,Z_q\lvert s\rangle ,
$$
so each non-Clifford gate at most *doubles* the number of stabilizer terms needed to write the state. After $t$ such gates a naive decomposition carries $2^t$ terms. The framework's `StabilizerFrame` realizes exactly this naive doubling, with no merging: the component count is $2^k$ for $k$ non-Clifford gates, confirmed bit-for-bit out to $k=18$ ($262144$ components, $0.6$ GB) in [`StabilizerFrame-NonClifford-Blowup.md`](StabilizerFrame-NonClifford-Blowup.md).

```wl
stabRank[x_] := If[MatchQ[x, _StabilizerFrame], x["Length"], 1];
interleaved[k_] := PauliStabilizer[1][Prepend[Flatten[Table[{"T", "H"}, k]], "H"]];
Table[{k, stabRank[interleaved[k]]}, {k, 0, 8}]   (* {1,2,4,8,16,32,64,128,256} = 2^k *)
```

The number that *controls* the cost is the **stabilizer rank** $\chi$, the fewest stabilizer states whose linear combination is the target. The frame crudely upper-bounds it by $2^t$; the true rank is far smaller, $\chi \lesssim 2^{\beta t}$ with $\beta \le 0.47$ exact and $\approx 0.23$ for approximate (additive-error) decompositions (Bravyi-Gosset, [arXiv:1601.07601](../literature/pathsum-literature/1601.07601/); Bravyi et al., [arXiv:1808.00128](../literature/pathsum-literature/1808.00128/)). The exponent is what matters: it is the same in $2^t$ and in $2^{0.47 t}$, only the base shrinks.

The exponent cannot be removed, and this is not an implementation weakness. Exact strong simulation of $\mathrm{Clifford}+T$ is $\#\mathrm{P}$-hard: the amplitude of a $\{H, Z, CZ, CCZ\}$ circuit is, up to a power-of-$\sqrt2$ normalization, the Gauss sum $\mathrm{gap}(f)=\sum_x(-1)^{f(x)}$ of a degree-$\le 3$ polynomial $f$ over $\mathbb{F}_2$; at degree $2$ this is computable in $O(n^3)$ (Gottesman-Knill in polynomial form), and at degree $3$ it is $\#\mathrm{P}$-hard (Montanaro, [arXiv:1607.08473](../literature/pathsum-literature/1607.08473/)). Because $\mathrm{Clifford}+T$ is universal, any polynomial-time method that folded *generic* $T$ branching exactly would put $\mathrm{BQP}$ in $\mathrm{P}$ and collapse the polynomial hierarchy. So the exponential-in-magic is mandatory; the entire game is to detect and exploit *structure* so that the exponent is governed by something smaller than the raw $T$-count. The sibling note [`SymPhase-vs-StabilizerRank-Algebraic-Branching.md`](SymPhase-vs-StabilizerRank-Algebraic-Branching.md) traces precisely where the magic relocates: a diagonal $\mathrm{Clifford}+T$ circuit absorbs it into one rank-$1$ degree-$\le 3$ $\mathbb{Z}_8$ phase polynomial, but the first Hadamard interleaved between two $T$ gates turns it back into a $\#\mathrm{P}$-hard Gauss sum.

---

## 2. The alternatives: one axis per method (answering Question 1)

The object every entry tries to beat is the eager frame, $\chi = 2^t$. The improvements are *not* the same kind. The substance is the distinction: some shrink the base of an exponential that stays exponential, one compresses along an axis unrelated to magic, one is a graphical heuristic, one is a fixed-parameter program in a graph width, and the path sum exploits the algebra of the phase polynomial. Below, each with its exploited axis and its win/lose regime, anchored to the fetched papers.

### 2.1 Eager stabilizer frames (the baseline to beat)

**Axis: none.** The frame stores one $\{c_i, \lvert s_i\rangle\}$ component per branch and doubles per $T$, with no deduplication or rank reduction. It is the thing the others improve on. **Win:** small $T$-count, where it is exact, faithful, and now (after the kernel fix at `61bc0e39`) a correct dense-readout object for one state. **Lose:** $T$-rich circuits, and even rank-$1$ states: $k$ $T$ gates on $\lvert 0\rangle$ produce $2^k$ identical components for the single state $\lvert 0\rangle$ ([`StabilizerFrame-NonClifford-Blowup.md`](StabilizerFrame-NonClifford-Blowup.md), Section 3). The companion [`StabilizerFrame-Observables-Design.md`](StabilizerFrame-Observables-Design.md) shows that observables can at least be *read off* a frame of whatever rank it has in $\chi^2\,\mathrm{poly}(n)$ time without densifying, but that does not compress $\chi$.

### 2.2 Stabilizer rank, extent, and robustness decompositions

**Axis: a convex magic monotone.** The stabilizer rank $\chi$ is the $\ell_0$ quantity; the stabilizer extent (Bravyi et al., [arXiv:1808.00128](../literature/pathsum-literature/1808.00128/)), multiplicative across few-qubit factors, bounds the approximate rank; the robustness of magic (Howard-Campbell, [arXiv:1609.07488](../literature/pathsum-literature/1609.07488/)) prices the sampling cost as $\mathcal{R}^2$, numerically $1.685^t$ down to a proven floor of $1.207^t$. **Win:** generic magic (these are the right tools when there is no exploitable structure), and weak simulation / sampling, where the approximate exponent $2^{0.23 t}$ is a genuine improvement. **Lose:** they never remove the exponent. They trade $2^t$ for $b^t$ with a fixed base $b>1$. They set the yardstick; they are not what beats it. The framework keeps $2^t$ and does not yet implement any of this compression (the Garcia-Markov frame-merging it cites is not implemented), so relative to an optimized stabilizer-rank simulator the QF frame is exponentially worse in term count ([`StabilizerFrame-NonClifford-Blowup.md`](StabilizerFrame-NonClifford-Blowup.md), Section 5).

### 2.3 LIMDD decision diagrams

**Axis: local-Clifford redundancy of the magic, captured by Pauli-labeled diagram edges.** A LIMDD (Vinkhuijzen et al., [arXiv:2108.00931](../literature/pathsum-literature/2108.00931/)) is a decision diagram whose edges carry local Pauli operators; as a class it strictly contains the stabilizer states together with the states reachable by ordinary decision diagrams (QMDD) and by matrix product states. **Win:** families such as $W$-states and Hamming-weight-controlled Cliffords are poly-size where a $\mathrm{Clifford}+T$ stabilizer decomposition costs $\Omega(2^n)$, because the magic there is local-Clifford-redundant and the Pauli edge labels quotient it out. **Lose:** magic that is *not* local-Clifford-redundant; the diagram then has no compact structure to exploit. Whether the poly-size LIMDD class dominates bounded stabilizer rank outright is open (it would need an explicit provably-exponential-rank family, itself a standing open problem).

### 2.4 ZX-calculus reduced stabilizer decompositions

**Axis: graphical simplification structure.** Interleaving a stabilizer decomposition with ZX-calculus rewriting (Kissinger-van de Wetering, [arXiv:2109.01076](../literature/pathsum-literature/2109.01076/)) cancels non-Clifford generators *between* decomposition steps, pushing the exponent from the raw $2^{0.47 t}$ down to roughly $2^{0.40 t}$ and reaching $50$ to $100$ qubits with more than $1000$ $T$ gates in the Rust tool `quizx`. **Win:** structured circuits, hidden shift especially, where the rewriting buys six or more orders of magnitude in term count. **Lose:** generic random magic, where the authors' own honest reading is that the rewriting buys only a constant factor of order ten over plain $T$-count reduction followed by a tableau simulator. It is real engineering, not an asymptotic gain on unstructured input.

### 2.5 Matrix product states and tensor networks (the wrong axis)

**Axis: entanglement (bond dimension), which is orthogonal to magic.** A $T$-rich but lightly entangled circuit *looks* like it should be cheap for an MPS, and it is not: magic and entanglement are independent quantities, and an MPS of modest bond dimension is already as non-stabilizer as a Haar-random state (Frau et al., [arXiv:2404.18768]). **Win:** genuinely low-entanglement dynamics, including the light-cone economy for *local* observables, where one keeps only the $T$ gates inside the observable's backward cone and the relevant rank is $2^{(\text{cone }T\text{-count})}$, often far below $2^t$ ([`SymPhase-vs-StabilizerRank-Algebraic-Branching.md`](SymPhase-vs-StabilizerRank-Algebraic-Branching.md), Section 5; the same light-cone win the QF tensor-network engine exploits, see the sibling [`../../Why QF/QF-Showcase-TN-LightCone-Draft.md`](../../Why%20QF/QF-Showcase-TN-LightCone-Draft.md)). But that win comes from *entanglement* being local, not from compressing magic; the cone still pays $2^{(\text{cone }T)}$. **Lose:** the magic axis itself. A tensor network neither needs nor exploits the structure at issue here, so it is the wrong tool for taming the $T$-blow-up per se.

### 2.6 The de Colnet rank-width sums-of-powers simulator (the direct competitor)

**Axis: rank-width of the path-variable interaction graph.** de Colnet, Geerts, Hai, Laarman, Lee, Perez ([arXiv:2605.29944](../literature/pathsum-literature/2605.29944/)) evaluate the *identical* $\mathrm{Clifford}+T$ path-sum object (citing Dawson, not Amy) by a fixed-parameter dynamic program, equivalently graphical-model bucket elimination dispatched onto BDD and weighted-model-counting engines: time $2^{O(k)}\mathrm{poly}(n)$ in the rank-width $k$ of the graph of path-variable interactions. **The $T$-count does not enter their exponent.** **Win:** low-rank-width circuits, where it out-scales a brute-force-completed path sum decisively. **Lose:** high-rank-width circuits. This is the most direct competitor to the path-sum representation and the central reason an honest account of the WL path sum does not overclaim performance: where the WL build pays $2^{\lvert V\rvert}$ to brute-force its residual, de Colnet pays $2^{O(k)}$ in a graph width that can be far smaller. The rank-width axis is *distinct* from both Amy's degree bound and the Groebner-ideal regularity of Section 4.2.

### 2.7 The symbolic path sum itself

**Axis: the algebraic structure of the degree-$\le 3$ phase polynomial.** The amplitude is carried as a Gauss sum $2^{-h/2}\sum_v \omega^{\varphi(v)}$ and reduced by the $[\mathrm{HH}]/[\omega]/[\mathrm{Elim}]$ rewrite rules (Amy, [arXiv:1805.06908](../literature/pathsum-literature/1805.06908/); confluence for the $\mathbb{F}_2$ sign fragment, Amy-Stinchcombe, [arXiv:2408.02778](../literature/pathsum-literature/2408.02778/)). The cost is set not by the gate count but by the structure of $\varphi$: Montanaro makes this exact, with the sum computable in time $2^{\lvert S\rvert}\mathrm{poly}(n)$ for a hitting set $S$ of the cubic terms, and $2^v\mathrm{poly}(n)$ for the number $v$ of essential variables surviving an optimal $GL_n(\mathbb{F}_2)$ change of basis, both of which can fall far below the $T$-count. **Win:** structured magic (hidden shift in genuine polynomial time, zero residual), Clifford telescoping (a $54$-variable sum collapses to $0$ internal residual), and the diagonal fragment (rank $1$, one closed-form term). **Lose:** narrow-but-deep magic, where the residual is $k$ magic variables and the engine brute-forces $2^k$ paths; the WL build crosses over near $k \approx 16$ and loses to dense past it ([`pathsum-build/PathSum-WL-Note.md`](pathsum-build/PathSum-WL-Note.md), Section 4.5). On its own axis the path sum is real and exact; it is *not* uniquely best, and de Colnet out-scales its brute-force completion.

### 2.8 Summary table

| Method | Exploited axis | Wins on | Loses on |
|---|---|---|---|
| Eager stabilizer frame | none | small $T$-count, exact small states | $T$-rich, even rank-$1$ states ($2^k$ copies) |
| Stabilizer rank / extent / robustness | a convex magic monotone (shrinks base $b^t$) | generic magic, sampling ($2^{0.23 t}$) | never removes the exponent |
| LIMDD decision diagram | local-Clifford redundancy of the magic | $W$-states, Hamming-weight-controlled Cliffords | non-redundant magic |
| ZX-reduced stabilizer decomposition | graphical (ZX) simplification | structured (hidden shift, $>10^6\times$) | generic random magic (constant factor) |
| MPS / tensor network | entanglement (bond dimension) | low-entanglement, light-cone local observables | the magic axis itself (wrong axis) |
| de Colnet sums-of-powers | rank-width of path-variable graph | low-rank-width ($2^{O(k)}\mathrm{poly}(n)$) | high-rank-width |
| Symbolic path sum | algebra of the phase polynomial (hitting set / essential variables) | structured magic, Clifford telescoping, diagonal fragment | narrow-deep magic ($2^k$ paths $> 2^n$) |

**The answer to Question 1 is therefore: no, the symbolic path sum is not the only solution, and it is not even the best on its own object.** It is one structure-exploiting representation among several, each keyed to a different axis. None removes the exponent in magic, because universality forbids that. The right way to read the path sum is as the *algebraic* member of this family: the one that makes the structure of the Gauss sum explicit and exact rather than hidden in a $2^n$ vector. That is what sets up Question 2.

---

## 3. The most original contribution WL can prove: symbolic-parameter path sums (Question 2)

### 3.1 The claim, and why it is the right altitude

Every competing simulator named in Section 2 is built on a *fixed numeric coefficient ring*: `feynman` works in the dyadic ring $\mathbb{D}[\omega]$ with concrete gate phases (and in fact overflows its machine-integer dyadic arithmetic at $31$ bits, [arXiv:1805.06908](../literature/pathsum-literature/1805.06908/)); `quizx` in $\mathbb{D}[e^{i\pi/4}]$; de Colnet over $\mathbb{Z}_r$ with numeric weights; Koh-Penney-Spekkens ([arXiv:1702.03316](../literature/pathsum-literature/1702.03316/)) evaluate numeric Weil sums; SymPhase symbolizes only the *measurement-outcome sign bits*, never the gate angles. A grep of all thirteen fetched papers for symbolic / free / parametric *gate* parameters returns nothing: the path-sum literature is "symbolic" only in the path variables, never in the gate angles ([`WL-Symbolic-Novelty-Reexamination.md`](WL-Symbolic-Novelty-Reexamination.md), C1).

Hosting the path sum in a computer-algebra system removes that restriction. Carry a gate angle $\theta$ as a *free symbol* through the Gauss sum, and the amplitude comes out as a finite sum of complex exponentials $e^{i\Phi(v)}$ with $\Phi$ linear in each path variable, which `FullSimplify` collapses to a closed form *as a function of $\theta$*. From that single object one reads the amplitude, the Born probability, the expectation, and, by symbolic differentiation, the exact gradient, over the *entire parametrized circuit family*, in one pass. This is the object a variational circuit (VQE, QAOA) needs: an exact $\langle H\rangle(\vec\theta)$ with exact gradients and exact landscape geometry, computed once rather than sampled on a grid of numerical values.

This is the single most original *and* most defensible contribution for three reasons, all of which the competing items fail:

1. **It is a capability no competitor has** (Section 2 ring-by-ring), not a speed ratio. The dense path produces no symbolic-$\theta$ amplitude at all.
2. **It is provable**: it is an exact algebraic identity in the ring $\mathbb{Z}[\omega, 1/\sqrt2]$ extended by $e^{i\theta}$, not a heuristic or an empirical speedup.
3. **It is machine-verifiable**: the closed form can be checked against dense `Method -> "Schrodinger"` at *every* angle (exactly, via `RootReduce`, at special angles where the state lies in $\mathbb{Z}[\omega,1/\sqrt2]$), and the symbolic gradient against a dense finite difference.

I prefer this over the Groebner/variety point count (Section 4.2) because the latter, while genuinely a distinct algorithm with its own tractability axis and the first concrete realization of an Amy-Stinchcombe open problem, is *bounded* to the $\mathbb{F}_2$ real-amplitude fragment, does not cleanly cover the $\mathbb{Z}_8$ $T$-rich case, and is by the originators' own statement no asymptotic free lunch. The symbolic-parameter capability has no such boundary: it applies to any parametrized $\mathrm{Clifford}+T$ family and is a clean exact identity throughout.

### 3.2 The worked demonstration, kernel-verified

All of the following is reproduced in a clean kernel by [`pathsum-blowup-alternatives-demos.wls`](pathsum-blowup-alternatives-demos.wls), no `Quiet`, no `Print`. The script defines a standalone symbolic-angle path-sum builder: each gate adds its phase to a symbolic real polynomial $\Phi$ with amplitude factor $e^{i\Phi}$ (a Hadamard on qubit $q$ introduces a fresh path variable $v$, the factor $1/\sqrt2$, and the sign $(-1)^{e_q v}=e^{i\pi e_q v}$, then sets $e_q := v$; a phase gate $P(\theta)=\mathrm{diag}(1,e^{i\theta})$ adds $\theta\,e_q$; $T$ adds $\tfrac{\pi}{4}e_q$; $CZ$ adds $\pi\,e_a e_b$), where $e_q$ is the integer $0/1$ value of the $\mathbb{F}_2$ wire polynomial.

**One free angle.** For $U(\theta) = H\,P(\theta)\,H$ on $\lvert 0\rangle$ the builder returns, in one pass,
$$
\langle x \lvert H\,P(\theta)\,H \rvert 0\rangle \;=\; \Bigl\{\tfrac{1+e^{i\theta}}{2},\ \tfrac{1-e^{i\theta}}{2}\Bigr\},
$$
exactly the `StabilizerFrame` component coefficients, and the local expectation as a closed-form trig function,
$$
\langle Z\rangle(\theta) \;=\; \cos\theta ,
$$
both verified `=== True` symbolically. The closed form equals dense `Method -> "Schrodinger"` *exactly* (`RootReduce`, no floating point) at the eight special angles $\{0, \tfrac{\pi}{4}, \tfrac{\pi}{3}, \tfrac{\pi}{2}, \tfrac{2\pi}{3}, \pi, \tfrac{\pi}{5}, \tfrac{\pi}{7}\}$.

**The gradient (the variational payoff).** Differentiating the one closed form gives the exact analytic gradient
$$
\frac{d}{d\theta}\langle Z\rangle(\theta) \;=\; -\sin\theta ,
$$
which matches a central finite difference of dense simulation at four sampled angles to better than $10^{-5}$. No competing simulator yields this object; they would each need a re-run grid and a numerical derivative.

**Two free angles, two qubits.** For the layer $U(\alpha,\beta) = (H\!\otimes\! H)\,P(\alpha)_1 P(\beta)_2\,CZ_{12}\,(H\!\otimes\! H)$ on $\lvert 00\rangle$, one symbolic pass produces the closed-form two-qubit correlator
$$
\langle Z_1 Z_2\rangle(\alpha,\beta) \;=\; \sin\alpha\,\sin\beta ,
$$
which matches dense *exactly* on a $25$-point grid of special angles $\{0,\tfrac{\pi}{4},\tfrac{\pi}{3},\tfrac{\pi}{2},\pi\}^2$ and at a generic numeric point $(0.7, 1.3)$ to $< 10^{-12}$. This is precisely the energy-landscape object a two-parameter variational ansatz needs, delivered as a closed form rather than a grid of samples.

**A free angle composed with fixed magic.** For $U(\theta) = H\,T\,P(\theta)\,H$ on $\lvert 0\rangle$, the closed form carries *both* the fixed eighth-root phase $e^{i\pi/4}$ and the free $\theta$ in the same pass, and matches dense exactly at $\{0,\tfrac{\pi}{5},\tfrac{\pi}{2},\pi\}$. The symbolic capability survives composition with genuine $T$-magic.

The point is not the particular trig identities; it is that they are produced *once*, *exactly*, *as functions*, over a whole family, by an engine hosted in a computer-algebra system, and that this is verifiable against dense simulation at every angle. That is the cleanest defensible sense in which the Wolfram Language host enables something no other listed $\mathrm{Clifford}+T$ simulator can.

---

## 4. Two honest companions

### 4.1 The symbolic capability inherits the wall off the structured fragment

The closed form is only as compact as the structure allows. Off the Clifford-foldable / low-residual fragment, the amplitude as a function of $\theta$ has $2^{\text{residual}}$ exponential terms in $\theta$ and inherits the $\#\mathrm{P}$ wall: "amplitude as a function of $\theta$" does not enlarge the tractable class. Within the structured fragment it is a capability the competitors lack; outside it, the same blow-up of Section 1 reappears, now as a sum of exponentials in $\theta$ instead of a sum of numbers. The capability is genuine and unique; it is not a loophole in the complexity of the sum.

### 4.2 The Groebner / variety point count (a bounded, distinct second contribution)

For the Toffoli-Hadamard fragment (gate set $\{H, Z, CZ, CCZ\}$, a real phase $(-1)^f$), Montanaro's correctness proof gives the amplitude as a *signed point count* of the $\mathbb{F}_2$ variety $\{f=0\}$:
$$
\langle 0\lvert C\rvert 0\rangle = \frac{\mathrm{gap}(f_C)}{2^{h/2+\ell}}, \qquad \mathrm{gap}(f) = \#\{x: f(x)=0\} - \#\{x: f(x)=1\} = 2\,\#\{f=0\} - 2^n .
$$
The conceptual content (amplitude $=$ variety point count) is the original Dawson 2004 paper ([arXiv:quant-ph/0408129](../literature/pathsum-literature/quant-ph_0408129/), in so many words: "counting the number of points in an algebraic variety"); the technique of absorbing the *whole* phase polynomial into an $\mathbb{F}_2$ ideal and solving by Groebner machinery is what Amy-Stinchcombe name as explicit future work, with the candid caveat that "in one sense this offloads the difficulty of simulation to the similarly difficult problem of computing Groebner bases and counting points in varieties." Its tractability axis is the *regularity of the ideal*, distinct from both Amy's degree bound and de Colnet's rank-width (Section 2.6), and the realization is genuinely natural in the Wolfram Language and un-built in the cited literature.

Kernel-verified in [`pathsum-blowup-alternatives-demos.wls`](pathsum-blowup-alternatives-demos.wls): the Nullstellensatz amplitude-zero test is one call (an inconsistent $\mathbb{F}_2$ system gives `GroebnerBasis[..., Modulus -> 2]` $= \{1\}$, hence amplitude $0$); the nonlinear point count via the Groebner quotient dimension (adjoin the field equations $x_i^2 - x_i$, take the basis, count standard monomials) matches brute force exactly for a quadratic term $x_1 x_2 + x_3$ (count $4$, $\mathrm{gap}=0$) and a cubic $CCZ$-type term $x_1 x_2 x_3 + x_4$.

The boundary is sharp. This is clean only for the $\mathbb{F}_2$ real-amplitude fragment, where "absorb the phase into an $\mathbb{F}_2$ ideal" is literally possible. The full $\mathrm{Clifford}+T$ amplitude carries a $\mathbb{Z}_8$ phase, which is not a function into $\mathbb{F}_2$; recovering $\sum_v \omega^{\varphi(v)}$ would require stratifying by the eight values of $\varphi \bmod 8$ and counting each $\mathbb{Z}_8$ level set, a congruence rather than an $\mathbb{F}_2$ polynomial equation, which blows the system up. That is why this is a *bounded* companion, not the headline.

---

## 5. What this does NOT claim (the boundary)

State plainly, so the positive results are not mistaken for a free lunch:

- **No $\#\mathrm{P}$ break.** Exact strong simulation of generic $\mathrm{Clifford}+T$ remains $\#\mathrm{P}$-hard. Nothing here changes which amplitudes are hard to compute.
- **No enlargement of the efficiently simulable class.** The symbolic-parameter capability inherits the wall off the structured fragment (Section 4.1); the Groebner/variety route does not break the wall and does not cover the $\mathbb{Z}_8$ $T$-rich case (Section 4.2); and the path sum itself is out-scaled on its own object by de Colnet's rank-width program (Section 2.6).
- **No new method.** The path-sum object, the degree bound, the rewrite rules, the residual completion, the hidden-shift application, the variety-point-count *idea* (Dawson 2004), and the rank-width fixed-parameter evaluation (de Colnet 2026) are all prior art. The Wolfram Language contribution is a *capability* (free gate parameters) and a *realization* (the Groebner point count), not a complexity-theoretic method.
- **Symbolic-parameter equality is heuristic.** `RootReduce` canonicalizes only parameter-free algebraic numbers (its `Details` restrict it to "integers and `Root` and `AlgebraicNumber` objects"), so two symbolic-$\theta$ amplitudes cannot in general be machine-decided equal by `RootReduce`; one falls back to `FullSimplify` with assumptions. The exactness is real; the *decision* of symbolic-parameter equality is not free.
- **The characteristic-2 wall blocks the alternative algebraic routes.** A primitive eighth-root phase cannot live in a Pauli sign: $(e^{i\pi/4}X)^2 = iI \ne I$ (verified `=== True`), so it stabilizes no state and cannot be a $\mathbb{Z}_8$ tableau sign; and $\omega$'s cyclotomic structure collapses entirely modulo $2$, $1 + x^4 \equiv (1+x)^4 \pmod 2$ (verified). The $GF(2^k)$ trace / character-sum reformulation and the $\mathbb{Z}_8$-ideal route are dead ends for exactly this reason ([`WL-Symbolic-Novelty-Reexamination.md`](WL-Symbolic-Novelty-Reexamination.md), candidate (c)).

In one line: the Wolfram Language host buys one thing no competing $\mathrm{Clifford}+T$ simulator buys, the freedom to carry the gate parameters as free symbols and return the amplitude, expectation, and gradient as exact functions of them over a whole circuit family in one pass; everything else, including the path sum itself, is a faithful and unusually comfortable re-expression of methods that already exist, and none of it moves the $\#\mathrm{P}$ wall.

---

## References

Fetched arXiv TeX sources are in [`../literature/pathsum-literature/`](../literature/pathsum-literature/) (the bracketed IDs link to the source folders).

- C. M. Dawson, H. L. Haselgrove, A. P. Hines, D. Mortimer, M. A. Nielsen, T. J. Osborne, *Quantum computing and polynomial equations over the finite field $\mathbb{Z}_2$*, QIC 5(2):102-112 (2005), [arXiv:quant-ph/0408129](../literature/pathsum-literature/quant-ph_0408129/). Amplitude as a point count of an $\mathbb{F}_2$ variety.
- A. Montanaro, *Quantum circuits and low-degree polynomials over $\mathbb{F}_2$*, [arXiv:1607.08473](../literature/pathsum-literature/1607.08473/). Amplitude $= \mathrm{gap}(f)$ of a degree-$\le 3$ $\mathbb{F}_2$ polynomial; degree-2 / degree-3 dichotomy; hitting-set and essential-variable tractability.
- D. E. Koh, M. D. Penney, R. W. Spekkens, *Computing quopit Clifford circuit amplitudes by the sum-over-paths technique*, QIC 17(13&14):1081-1095 (2017), [arXiv:1702.03316](../literature/pathsum-literature/1702.03316/). Closed-form Weil sums for odd characteristic; characteristic-2 exclusion.
- M. Amy, *Towards Large-scale Functional Verification of Universal Quantum Circuits*, QPL 2018, [arXiv:1805.06908](../literature/pathsum-literature/1805.06908/). The path-sum object and rewrite rules; 31-bit dyadic overflow.
- M. Amy, L. S. Stinchcombe, *Polynomial-Time Classical Simulation of Hidden Shift Circuits via Confluent Rewriting of Symbolic Sums*, Quantum (2025), [arXiv:2408.02778](../literature/pathsum-literature/2408.02778/). Confluent rewriting for the $\mathbb{F}_2$ fragment; the variety/Groebner future-work passage.
- S. Bravyi, D. Gosset, *Improved classical simulation of quantum circuits dominated by Clifford gates*, PRL 116, 250501 (2016), [arXiv:1601.07601](../literature/pathsum-literature/1601.07601/). Stabilizer rank $2^{\beta t}$, runtime $2^{0.23 t}$.
- S. Bravyi, D. Browne, P. Calpin, E. Campbell, D. Gosset, M. Howard, *Simulation of quantum circuits by low-rank stabilizer decompositions*, Quantum 3, 181 (2019), [arXiv:1808.00128](../literature/pathsum-literature/1808.00128/). Stabilizer extent; exact and approximate rank bounds.
- M. Howard, E. Campbell, *Application of a resource theory for magic states to fault-tolerant quantum computing*, PRL 118, 090501 (2017), [arXiv:1609.07488](../literature/pathsum-literature/1609.07488/). Robustness of magic; sampling cost $\mathcal{R}^2$.
- L. Vinkhuijzen, T. Coopmans, D. Elkouss, V. Dunjko, A. Laarman, *LIMDD: A Decision Diagram for Simulation of Quantum Computing Including Stabilizer States*, Quantum (2023), [arXiv:2108.00931](../literature/pathsum-literature/2108.00931/).
- A. Kissinger, J. van de Wetering, *Simulating quantum circuits with ZX-calculus reduced stabiliser decompositions*, [arXiv:2109.01076](../literature/pathsum-literature/2109.01076/). Software: `quizx`.
- A. de Colnet, F. Geerts, R. Hai, A. Laarman, J. H. Lee, G. A. Perez, *Quadratic Sums-of-Powers for Fixed-Parameter Tractable Quantum-Circuit Simulation*, [arXiv:2605.29944](../literature/pathsum-literature/2605.29944/). Rank-width fixed-parameter evaluation of the same path-sum object.
- C. Guan, K. W. Regan, *Stabilizer Circuits, Quadratic Forms, and Computing Matrix Rank*, [arXiv:1904.00101](../literature/pathsum-literature/1904.00101/). Closed-form Clifford ($\mathbb{Z}_4$) Gauss sum via $\mathbb{F}_2$ rank.
- W. Fang, M. Ying, *SymPhase: Phase Symbolization for Fast Simulation of Stabilizer Circuits*, [arXiv:2311.03906](../literature/pathsum-literature/2311.03906/). The degree-$1$ measurement-branching case the path sum generalizes.
- M. Frau, P. S. Tarabunga, et al., *Non-stabilizerness versus entanglement in matrix product states*, arXiv:2404.18768. Magic is independent of bond dimension.

### Sibling reports (relative paths)

- [`pathsum-build/PathSum-WL-Note.md`](pathsum-build/PathSum-WL-Note.md): the main arXiv-style note, honest provenance and benchmarks.
- [`Beyond-SymPhase-Symbolic-Algebra-For-Magic.md`](Beyond-SymPhase-Symbolic-Algebra-For-Magic.md): the candidate landscape and the symbolic path-sum prototype.
- [`WL-Symbolic-Novelty-Reexamination.md`](WL-Symbolic-Novelty-Reexamination.md): the skeptical re-examination whose verdict this note builds on.
- [`wl-symbolic-toolkit-docs/WL-Symbolic-Toolkit-for-PathSums.md`](wl-symbolic-toolkit-docs/WL-Symbolic-Toolkit-for-PathSums.md): the 91-function toolkit-to-pipeline mapping.
- [`SymPhase-vs-StabilizerRank-Algebraic-Branching.md`](SymPhase-vs-StabilizerRank-Algebraic-Branching.md): when "branching becomes algebra" transfers, and the diagonal-fragment rank-$1$ phase polynomial.
- [`StabilizerFrame-NonClifford-Blowup.md`](StabilizerFrame-NonClifford-Blowup.md): the $2^k$ frame blow-up, measured and mechanized.
- [`StabilizerFrame-Observables-Design.md`](StabilizerFrame-Observables-Design.md): reading observables off a frame in $\chi^2\,\mathrm{poly}(n)$ without densifying.
- [`tutorial/Stabilizer-Formalism-By-Computation.md`](tutorial/Stabilizer-Formalism-By-Computation.md): the "Magic and stabilizer rank" section.

*Verification script: [`pathsum-blowup-alternatives-demos.wls`](pathsum-blowup-alternatives-demos.wls), run against `Wolfram/QuantumFramework 2.0.0`; all assertions return `True` and all amplitudes are exact against dense `Method -> "Schrodinger"`.*
