---
Title: 'Does universality really forbid removing the magic exponent? An adversarial audit of one load-bearing sentence'
Author: investigation in the Wolfram QuantumFramework
Description: 'A first-principles, adversarial audit of the sentence "None removes the exponent in magic, because universality forbids that" (and its sibling "Clifford+T is universal, so any polynomial-time method that folded generic T branching exactly would put BQP in P and collapse the polynomial hierarchy") from PathSum-Blowup-Alternatives-And-WL-Contribution.md. The sentence is reduced to its minimal mathematical content (the amplitude is the gap of a degree-<=3 F_2 polynomial), theorem is separated from rhetoric, and the impossibility is stress-tested rather than restated. Verdict: the slogan is true under one precise reading (no universal poly-time exact strong simulator exists) but overclaimed and imprecise as written, because the report''s own body relocates the exponent off the T-count onto structural parameters that can be far smaller; the BQP/PH-collapse chain also has a broken middle link (the trigger is PostBQP=PP, not BQP in P). Every computational claim is verified in a live kernel.'
---

# Does universality really forbid removing the magic exponent?

## 0. Verdict, up front

The audited sentence is **true under one precise reading and overclaimed (imprecise) as written.**

- The defensible, airtight statement is: *no method makes the worst case polynomial.* Computing an arbitrary $\mathrm{Clifford}+T$ amplitude exactly is $\#\mathrm{P}$-hard, so no algorithm runs in $\mathrm{poly}(n, t)$ on all inputs. That much is a theorem and nothing in the survey threatens it.
- The phrase actually used, "removes **the exponent in magic**," says something stronger and different: that the exponent's *argument* is magic (the $T$-count, or a magic monotone). The report's **own body refutes that** in two places (Section 2.6 de Colnet, Section 2.7 Montanaro): several listed methods keep an exponential but move its argument off the $T$-count onto a **structural** parameter (rank-width, hitting-set size, essential-variable count, local-Clifford redundancy) that is bounded above by the magic yet can be arbitrarily smaller. For the class of circuits where that structural parameter is bounded, the cost stops scaling with magic at all. So "none removes the exponent in magic" is false for those classes; what is true is "none removes the worst-case exponential."
- The companion sentence's complexity chain has a **broken middle link**. "Would put $\mathrm{BQP}$ in $\mathrm{P}$ and collapse the polynomial hierarchy" is a non-sequitur: $\mathrm{BQP}\subseteq\mathrm{P}$ is not known to collapse $\mathrm{PH}$ (it is open whether $\mathrm{BQP}\subseteq\mathrm{PH}$ at all). The collapse comes from the *stronger* consequence of exact strong simulation, $\mathrm{PostBQP}\subseteq\mathrm{P}$, i.e. $\mathrm{P}=\mathrm{PP}$, via Toda's theorem.
- Finally, "because **universality** forbids that" misattributes the cause. The worst-case impossibility follows directly from the $\#\mathrm{P}$-hardness of a *combinatorial* counting problem (the gap of a cubic $\mathbb{F}_2$ polynomial), which needs no universality. Universality is what upgrades the impossibility to the dramatic $\mathrm{PH}$-collapse, not what creates it. Hardness is in fact *more robust* than universality: non-universal families (IQP, commuting circuits) are already hard.

**Proposed corrected wording** is in Section 5. In one line: replace "none removes the exponent in magic, because universality forbids that" with "**none makes the worst case polynomial (exact strong simulation is $\#\mathrm{P}$-hard); what several of them do is move the exponent off the $T$-count onto a structural parameter that can be far smaller.**"

Every concrete computational claim below is verified in a clean kernel by [`magic-exponent-audit-demos.wls`](magic-exponent-audit-demos.wls) (Section 9), no `Quiet`, no `Print`.

---

## 1. The disputed sentences, verbatim

From [`PathSum-Blowup-Alternatives-And-WL-Contribution.md`](PathSum-Blowup-Alternatives-And-WL-Contribution.md), three near-identical forms:

> (abstract) "...each exploiting a different axis, **none removing the exponent in magic, which universality forbids**."

> (Section 1) "Because $\mathrm{Clifford}+T$ is universal, **any polynomial-time method that folded generic $T$ branching exactly would put $\mathrm{BQP}$ in $\mathrm{P}$ and collapse the polynomial hierarchy. So the exponential-in-magic is mandatory**; the entire game is to detect and exploit *structure* so that the exponent is governed by something smaller than the raw $T$-count."

> (Section 2.8) "**None removes the exponent in magic, because universality forbids that.**"

A first observation worth making immediately: the Section 1 sentence already contains, in its own second half, the *correct* statement ("the entire game is to detect and exploit structure so that the exponent is governed by something smaller than the raw $T$-count"). The slogan and the precise statement sit one clause apart and contradict each other. The slogan says the exponent *is* in magic; the next clause says the exponent is governed by *something smaller than the $T$-count*. The audit below is, in large part, the claim that the report should keep the second clause and retire the slogan.

The sibling notes inherit the same slogan: [`SymPhase-vs-StabilizerRank-Algebraic-Branching.md`](SymPhase-vs-StabilizerRank-Algebraic-Branching.md) ("the exponential is **relocated, not removed**, which is exactly what $\mathrm{Clifford}+T$ universality demands") and [`Beyond-SymPhase-Symbolic-Algebra-For-Magic.md`](Beyond-SymPhase-Symbolic-Algebra-For-Magic.md) ("a polynomial-time symbolic method that folded generic $T$ branching exactly would put $\mathrm{BQP}$ in $\mathrm{P}$ and collapse the polynomial hierarchy"). The "relocated, not removed" phrasing is actually the better one and survives this audit largely intact; it is the "in magic" and "$\mathrm{BQP}$ in $\mathrm{P}$" phrasings that do not.

---

## 2. The minimal representation

Strip the problem to its smallest mathematical core, exactly as the reports themselves derive it.

A $\mathrm{Clifford}+T$ amplitude is, up to a power-of-$\sqrt2$ normalization, a **Gauss sum** over the Hadamard branch variables $v\in\mathbb{F}_2^h$,
$$
\langle y\lvert C\rvert 0\cdots0\rangle \;=\; 2^{-h/2}\sum_{v\,:\,E(v)=y}\omega^{\varphi(v)},\qquad \omega=e^{i\pi/4},
$$
where $\varphi$ is a pseudo-Boolean phase polynomial valued in $\mathbb{Z}_8$ and of **degree at most $3$** (the CNOT-dihedral normal form: $T$ a linear term, $CZ$ a quadratic term $4x_ix_j$, $CCZ$ a cubic term $4x_ix_jx_k$). For the real-amplitude (Toffoli-Hadamard) fragment, gate set $\{H,Z,CZ,CCZ\}$, this collapses to a signed point count of an $\mathbb{F}_2$ variety,
$$
\langle 0\lvert C\rvert 0\rangle \;\propto\; \mathrm{gap}(f) \;=\; \sum_{x\in\mathbb{F}_2^n}(-1)^{f(x)} \;=\; \#\{f=0\}-\#\{f=1\} \;=\; 2\,\#\{f=0\}-2^n ,
$$
with $f$ of degree $\le 3$. The gap identity is verified on a random cubic in Section 9 ([demo (1)](magic-exponent-audit-demos.wls)).

So the entire dispute reduces to one clean question:

> **Minimal question.** Is "remove the magic exponent" equivalent to "compute $\mathrm{gap}(f)$ of an arbitrary degree-$3$ polynomial $f$ over $\mathbb{F}_2$ in time $\mathrm{poly}(n)$"?

If yes, the impossibility is exactly the hardness of that counting problem, and the slogan stands or falls with it. The next section shows the equivalence holds **only for the universal-poly reading**, and that the survey's own methods break the equivalence by parameterizing the cost differently.

Two structural facts about the cost, both standard and both re-derivable from the sum:

- **Degree $2$ is easy.** A quadratic $f$ is an $\mathbb{F}_2$ quadratic form; its gap is a classical quadratic Gauss sum fixed by the rank and Arf invariant of the form, computable in $O(n^3)$. This is Gottesman-Knill in polynomial-counting clothes. Verified on the smallest case $\mathrm{gap}(x_1x_2)=2$ ([demo (1)](magic-exponent-audit-demos.wls)).
- **Degree $3$ is the wall.** Counting zeros of a general cubic over $\mathbb{F}_2$ is $\#\mathrm{P}$-hard. This is a purely combinatorial fact about $(\mathrm{AND},\mathrm{XOR})$-counting; it does *not* mention quantum circuits, universality, or $\mathrm{BQP}$. The single $CCZ$ term has $\mathrm{gap}(x_1x_2x_3)=6$ ([demo (1)](magic-exponent-audit-demos.wls)); the hardness is in *combining many* cubic terms with no exploitable structure.

---

## 3. Theorem versus rhetoric

This section pins down precisely what the mathematics forbids and what it does not. There are two distinct impossibility theorems doing the work, plus one rhetorical chain that is broken. Keeping them apart is the whole audit.

### 3.1 The airtight core (no universality required)

**Theorem A (combinatorial).** Computing $\mathrm{gap}(f)$ for an arbitrary degree-$3$ polynomial $f$ over $\mathbb{F}_2$ is $\#\mathrm{P}$-hard. Hence there is no algorithm that computes *every* $\mathrm{Clifford}+T$ amplitude exactly in time $\mathrm{poly}(n,t)$ unless $\mathrm{FP}=\#\mathrm{P}$ (which itself collapses $\mathrm{PH}$ by Toda's theorem, $\mathrm{PH}\subseteq\mathrm{P}^{\#\mathrm{P}}$).

This is the sound core of the slogan. Notice what it does **not** use: it never mentions universality, $\mathrm{BQP}$, postselection, or the polynomial hierarchy as a *premise*. It is a statement about a counting problem. The minimal question of Section 2 answers "yes" precisely at this reading: "remove the magic exponent" $=$ "compute the cubic gap in poly time" $=$ "collapse a $\#\mathrm{P}$-hardness," and that is forbidden (modulo the universally-believed $\mathrm{P}\ne\#\mathrm{P}$).

The phrase "**the exponential-in-magic is mandatory**" is correct **only if** "mandatory" is read as "no single algorithm handles all instances in polynomial time." Read instance-by-instance it is false, because structured instances are easy (Section 3.3).

### 3.2 The $\mathrm{PH}$-collapse argument, and the broken middle link

**Theorem B (postselection).** Exact (or constant-multiplicative-error) strong simulation of a *universal* family in polynomial time implies $\mathrm{PostBQP}\subseteq\mathrm{P}$. Since $\mathrm{PostBQP}=\mathrm{PP}$, this gives $\mathrm{P}=\mathrm{PP}$, and then $\mathrm{PH}\subseteq\mathrm{P}^{\mathrm{PP}}=\mathrm{P}$ collapses by Toda.

This is the theorem the report is *reaching for*. Here universality genuinely enters: the identity $\mathrm{PostBQP}=\mathrm{PP}$ is proved for a universal gate set. The trigger is precise and worth stating exactly, because the report blurs it:

- **Exact strong** simulation (compute the amplitude/probability exactly): triggers the collapse.
- **Multiplicative-error strong** simulation (each probability within a constant relative factor): also triggers it.
- **Exact weak** simulation (sample exactly from the output distribution): also triggers it.
- **Additive-error weak** simulation (sample from a distribution close in total variation): does **not** trigger it on its own. This is the regime of Bravyi-Gosset $2^{0.23 t}$ and of the quantum-supremacy program, where hardness needs extra conjectures (anti-concentration plus average-case hardness).

The report's word "exactly" ("folded generic $T$ branching **exactly**") correctly selects the exact-strong trigger. The conclusion ("collapse the polynomial hierarchy") is correct. But the **middle term is wrong**:

> "would put $\mathrm{BQP}$ in $\mathrm{P}$ and collapse the polynomial hierarchy"

$\mathrm{BQP}\subseteq\mathrm{P}$ does not collapse $\mathrm{PH}$. It is not even known that $\mathrm{BQP}\subseteq\mathrm{PH}$; there is oracle evidence (Raz-Tal) that $\mathrm{BQP}\not\subseteq\mathrm{PH}$ relative to an oracle. So one cannot route "collapse $\mathrm{PH}$" through "$\mathrm{BQP}$ in $\mathrm{P}$." Exact strong simulation gives *more* than $\mathrm{BQP}\subseteq\mathrm{P}$: it gives $\mathrm{PostBQP}\subseteq\mathrm{P}$, and *that* is what collapses $\mathrm{PH}$. The fix is to replace "$\mathrm{BQP}$ in $\mathrm{P}$" with "$\mathrm{PostBQP}\subseteq\mathrm{P}$, i.e. $\mathrm{P}=\mathrm{PP}$."

This is a rhetorical defect, not a fatal one: the conclusion survives, and in any case Theorem A already delivers the worst-case impossibility without the postselection detour. But since the user asked for the chain to be pinned down exactly, the chain as written does not connect.

### 3.3 What is NOT forbidden

Neither theorem forbids any of the following, and the survey itself relies on several of them.

- **Structured instances.** Theorems A and B are worst-case. A circuit with enormous $T$-count can have a poly-time exact amplitude if its phase polynomial has structure. Two such families are exhibited and verified in Section 9; they are the heart of the correction.
- **The diagonal / phase-polynomial fragment.** Diagonal $\mathrm{Clifford}+T$ (no interior Hadamard) is rank $1$: the amplitude is a single closed-form phase regardless of $t$. The QF kernel does not exploit this (it still pays $2^t$, [demo (4)](magic-exponent-audit-demos.wls)), but the *intrinsic* cost is $O(1)$.
- **Weak/additive sampling at a shrunk base.** Bravyi-Gosset's $2^{0.23 t}$ is an additive-error weak simulator. It does not trigger the collapse and it does not make the worst case poly; it shrinks the *base* of the exponential.
- **Moving the exponent's argument.** Nothing forbids an algorithm whose cost is $2^{O(k)}\mathrm{poly}(n)$ with $k$ a parameter that is *upper-bounded by the magic but can be far smaller*, and indeed bounded on whole families with unbounded magic. The de Colnet rank-width program is exactly this, and it is in the report's own Section 2.6.

The last bullet is where the slogan breaks, so it deserves its own section.

---

## 4. The three readings of "removes the exponent in magic"

The slogan is ambiguous across three readings. It is true under one and false under two, and the report's body lives in the readings where it is false.

**Reading R1: "makes all $\mathrm{Clifford}+T$ poly-time" (a universal polynomial algorithm).**
Under R1 the slogan is **true**: such an algorithm would compute the cubic gap in poly time (Theorem A) and collapse $\mathrm{PH}$. No listed method claims this. Verdict under R1: correct.

**Reading R2: "the exponent's argument is the $T$-count."**
Under R2 the slogan is **false**. Several listed methods keep an exponential whose argument is *not* $t$:

- *Montanaro hitting set / essential variables* (Section 2.7): cost $2^{\lvert S\rvert}\mathrm{poly}(n)$ for a hitting set $S$ of the cubic terms, and $2^v\mathrm{poly}(n)$ for the essential-variable count $v$ after an optimal $GL_n(\mathbb{F}_2)$ change of basis. Both $\lvert S\rvert$ and $v$ are $\le t$ but can be $O(1)$ while $t\to\infty$.
- *de Colnet rank-width* (Section 2.6): cost $2^{O(k)}\mathrm{poly}(n)$ in the rank-width $k$ of the path-variable interaction graph. The report states plainly "**the $T$-count does not enter their exponent.**"

The **star witness** (Section 9, [demo (2)](magic-exponent-audit-demos.wls)) makes R2-falseness concrete and exact. Take
$$
f_m(x_0,\dots,x_{2m}) \;=\; \sum_{i=1}^m x_0\,x_{2i-1}\,x_{2i}\quad\text{over }\mathbb{F}_2,
$$
a star of $m$ cubic ($CCZ$) terms sharing the hub $x_0$. The $T$-count proxy (number of cubic terms) is $m$, unbounded; the magic grows with $m$. But the hitting set of the cubic monomials is $\{x_0\}$, size $1$. Conditioning on $x_0$ splits the Gauss sum into two **quadratic** sums:
$$
\mathrm{gap}(f_m) \;=\; \underbrace{2^{2m}}_{x_0=0,\ f\equiv 0} \;+\; \underbrace{2^{m}}_{x_0=1,\ f=\sum_i x_{2i-1}x_{2i}} ,
$$
computed in $O(m)$ arithmetic. Brute force confirms $\mathrm{gap}(f_m)=2^{2m}+2^m$ exactly for $m=1,\dots,5$ ($6,20,72,272,1056$). The exponent here is in the hitting set ($=1$), not in the magic ($=m$).

**Reading R3: "the exponent is in a magic monotone (stabilizer rank, extent, robustness)."**
This is the most charitable reading, and under it the slogan is **partly true and partly false**. The stabilizer-rank/extent/robustness methods (Section 2.2) do keep their exponent in a magic monotone: they trade $2^t$ for $b^t$ with a fixed base $b>1$, exactly as the report says. But de Colnet's rank-width is **not a magic monotone at all**. It is a combinatorial property of the circuit's interaction graph, independent of how much non-stabilizerness the state carries.

The **matching witness** (Section 9, [demo (3)](magic-exponent-audit-demos.wls)) shows rank-width and magic are genuinely independent. Take
$$
g_m(x_1,\dots,x_{3m}) \;=\; \sum_{i=1}^m x_{3i-2}\,x_{3i-1}\,x_{3i}\quad\text{over }\mathbb{F}_2,
$$
a matching of $m$ disjoint $CCZ$ terms on $3m$ variables. The hitting set is now $m$ (no shared variable), so Montanaro's hitting-set axis does *not* see this as easy. But the interaction graph is $m$ disjoint triangles: bounded tree-width and bounded rank-width. The gap factorizes over the components,
$$
\mathrm{gap}(g_m) \;=\; \prod_{i=1}^m \mathrm{gap}(x_{3i-2}x_{3i-1}x_{3i}) \;=\; 6^m ,
$$
computed in $O(m)$ arithmetic; brute force confirms $6,36,216$ for $m=1,2,3$. Here the magic is unbounded ($m$ $CCZ$ gates, robustness growing with $m$), the rank-width is bounded, and the cost is poly. This is a bounded-width, unbounded-magic family: de Colnet's exponent is bounded on it while every magic monotone diverges. So "the exponent is in magic" fails even under R3 for at least one listed method.

**Summary of the readings.**

| Reading of "removes the exponent in magic" | Slogan is | Why |
|---|---|---|
| R1: a universal poly-time algorithm | true | would compute the cubic gap in poly time; $\#\mathrm{P}$-hard (Theorem A) |
| R2: exponent's argument is the $T$-count | false | Montanaro ($2^{\lvert S\rvert}$, $2^v$) and de Colnet ($2^{O(k)}$) put a non-$t$ parameter in the exponent (star, matching witnesses) |
| R3: exponent is in a magic monotone | mixed | true for stabilizer rank/extent/robustness (base-shrinking); false for de Colnet (rank-width is not a magic measure) |

The report's prose uses "magic" and "$T$-count" almost interchangeably, which is the R2/R3 register, which is exactly where the slogan is false.

---

## 5. Verdict and corrected wording

**Verdict: overclaimed and imprecise, with a sound core.** The intended content (no universal polynomial exact strong simulator) is correct and important. The sentence as written asserts more (that the exponent's argument is magic), and that stronger claim is contradicted by the report's own Sections 2.6 and 2.7. Two secondary defects: the $\mathrm{BQP}$/$\mathrm{PH}$ chain has a broken middle link, and "universality forbids" misattributes a hardness that is really combinatorial ($\#\mathrm{P}$) and more robust than universality.

The fix is small and keeps every true thing the report says. Concretely:

**(a) The abstract clause.** Replace
> "...each exploiting a different axis, none removing the exponent in magic, which universality forbids."

with
> "...each exploiting a different axis. None makes the worst case polynomial (exact strong simulation of $\mathrm{Clifford}+T$ is $\#\mathrm{P}$-hard); what several of them do is move the exponent off the $T$-count onto a structural parameter that can be far smaller."

**(b) The Section 2.8 slogan.** Replace
> "None removes the exponent in magic, because universality forbids that."

with
> "None removes the worst-case exponential: exact strong simulation of generic $\mathrm{Clifford}+T$ is $\#\mathrm{P}$-hard. What the structure-exploiting methods do is relocate the exponent off the raw $T$-count onto a parameter (rank-width, hitting set, essential-variable count, local-Clifford redundancy) that is bounded above by the magic but can be much smaller, and is only forced to be large in the worst case. The exponent is relocated, not deleted."

**(c) The Section 1 complexity chain.** Replace
> "any polynomial-time method that folded generic $T$ branching exactly would put $\mathrm{BQP}$ in $\mathrm{P}$ and collapse the polynomial hierarchy."

with
> "computing an arbitrary $\mathrm{Clifford}+T$ amplitude exactly is the gap of an arbitrary cubic $\mathbb{F}_2$ polynomial, a $\#\mathrm{P}$-hard counting problem, so no method does it for all inputs in polynomial time. (Equivalently, exact strong simulation of a universal family would give $\mathrm{PostBQP}\subseteq\mathrm{P}$, i.e. $\mathrm{P}=\mathrm{PP}$, collapsing the polynomial hierarchy by Toda's theorem.)"

The sibling note's existing phrasing "the exponential is **relocated, not removed**" is already correct and should become the canonical slogan across the folder; the "in magic / $\mathrm{BQP}$ in $\mathrm{P}$" forms should be retired.

---

## 6. Genuine exploration of escape routes

The mandate is to look for a crack, not to restate the wall. Below, every candidate that could "tame the cost," with what it buys and where it dies. The recurring pattern: each one **enlarges the tractable class or relocates the exponent**, and each one **dies on the worst case** by the same $\#\mathrm{P}$ argument. The honest distinction is between the moving *frontier* of tractable structure and the fixed *wall*.

**(1) Structural width parameters (rank-width, tree-width, hitting set, essential variables).**
*Buys:* poly-time exact amplitudes on whole families with unbounded $T$-count, whenever the parameter is bounded. The star and matching witnesses are the proof of concept; de Colnet's $2^{O(k)}\mathrm{poly}(n)$ is the general algorithm. This is the strongest honest sense in which "the magic exponent is removed" for a meaningful class.
*Dies:* the parameter is $\Theta(n)$ on expander-like or generic instances. A cubic whose interaction graph is an expander has rank-width $\Theta(n)$; a generic cubic has $\Omega(n)$ essential variables and no small hitting set. So none of these makes the worst case poly; each carves out a tractable sub-class.
*Is there a strictly more general width nobody has named?* Possibly, and this is a live research frontier (rank-width already strictly generalizes tree-width and clique-width for this purpose). But any such parameter is bound by the same constraint: it must blow up on the $\#\mathrm{P}$-hard instances, or it would collapse $\mathrm{PH}$. So a better parameter can *enlarge* the tractable class; it cannot make it everything.

**(2) Non-stabilizer poly-describable representations (LIMDD and successors).**
*Buys:* poly-size representation of states whose magic is *local-Clifford-redundant* ($W$-states, Hamming-weight-controlled Cliffords), where a $\mathrm{Clifford}+T$ stabilizer decomposition costs $\Omega(2^n)$ (Section 2.3). The Pauli edge labels quotient out the redundancy.
*Dies:* a counting argument caps any single poly-size representation class. There are doubly-exponentially many inequivalent $n$-qubit states but only exponentially many poly-size diagrams, so almost all states need exponential diagram size; and any specific $\#\mathrm{P}$-hard amplitude resists a poly-size diagram unless the collapse follows. New data structures keep arriving and each enlarges the class; none can be universal-and-efficient.

**(3) Average-case / smoothed analysis (is worst-case even the right bar?).**
*Buys:* in principle, "practical magic" might be generic, and generic instances might be easy, making worst-case $\#\mathrm{P}$-hardness irrelevant to practice.
*Dies hard, for exact strong simulation.* Worst-case-to-average-case reductions (polynomial interpolation over the amplitude as a low-degree polynomial in the gate entries) show that computing output probabilities of random circuits exactly, or to high enough precision, is $\#\mathrm{P}$-hard *on average*. So "generic instances are easy" is false for the exact-strong task. For *weak/additive* sampling the average-case story is only conjectural (it is the basis of the quantum-supremacy arguments, requiring anti-concentration plus a permanent-of-Gaussians-style conjecture). Net: average-case opens no crack for the exact strong task the slogan is about.

**(4) Approximate strong / weak simulation (Bravyi-Gosset $2^{0.23 t}$).**
*Buys:* a genuinely smaller base. $2^{0.23 t}\approx 1.17^t$ is an additive-error weak (sampling) simulator, a real improvement over $2^t$ and over the exact rank $2^{0.47 t}$.
*Dies as a counterexample to the slogan:* its exponent's argument is still $t$. It shrinks the base; it does not relocate or remove the $t$-exponent. So it is *consistent* with the slogan under R2/R3 and is not a crack. The report states this correctly in Section 2.2; the audit agrees. (It is worth flagging only because "approximate" sometimes gets miscounted as "removing" the exponent: here it shrinks the base, full stop.)

**(5) A loophole in "universal $\Rightarrow$ hard"?**
*Buys:* if some universal family were classically simulable, the wall would have a hole.
*Dies:* universality plus efficient exact strong simulation lets you compute the hard amplitudes directly (Solovay-Kitaev approximation has only poly overhead), so no universal family is exact-strong simulable unless the collapse follows. But the more interesting point cuts the *other* way: hardness does **not** require universality. IQP circuits, commuting circuits, and boson sampling are non-universal yet classically hard (to sample, under standard conjectures). So "because universality forbids that" understates the situation. The hardness of the cubic gap is combinatorial and gate-set-local; universality is sufficient for the dramatic $\mathrm{PH}$-collapse but is neither necessary for hardness nor the cleanest route to the worst-case impossibility.

**(6) Encoded / measurement-based / fault-tolerant structure.**
*Buys:* nothing for classical simulability. Measurement-based computation is polynomially equivalent to the circuit model, so it inherits the same hardness. Fault-tolerant structure (transversal gates, magic-state distillation) is about suppressing noise, not about making simulation easier; the logical computation it protects is exactly the hard one.
*Dies:* immediately, by equivalence of models.

**(7) The diagonal / phase-polynomial fragment (a real win the framework leaves on the table).**
*Buys:* exact, rank-$1$ description of diagonal $\mathrm{Clifford}+T$ regardless of $T$-count; a single closed-form phase per amplitude. This is the cleanest "structure beats $t$" case and it is *not even exploited by QF*: [demo (4)](magic-exponent-audit-demos.wls) shows $T^t\lvert+\rangle$ stored as $2^t$ frame components (fidelity $1$) for a state whose true amplitude $\{1, e^{i t\pi/4}\}/\sqrt2$ is computed in $O(1)$. So here the exponent QF *pays* ($2^t$, nominally "in magic") is demonstrably not the intrinsic cost (which is $O(1)$). This is a concrete, in-kernel instance of the slogan's failure: the implementation conflates "frame component count" with "intrinsic hardness," and the two diverge maximally.
*Dies:* the first Hadamard interleaved between two $T$ gates leaves the diagonal fragment, and the cost reappears as a $\#\mathrm{P}$-hard Gauss sum (the sibling note's Section 4). The fragment is not universal and not even Clifford-closed.

---

## 7. The steelman, run to ground

Construct the strongest case that *some* approach removes the magic exponent for a meaningful class, then test it.

> **Steelman.** "de Colnet's program removes the magic exponent for all bounded-rank-width circuits. On that class the runtime is $2^{O(k)}\mathrm{poly}(n)$ with $k$ bounded and independent of the magic, so a family with unbounded $T$-count, unbounded stabilizer rank, unbounded robustness is simulated in polynomial time. The cost does not scale with magic at all. Therefore the magic exponent is, for this class, genuinely gone."

**This steelman survives.** The matching witness is a concrete, verified member of the class: unbounded magic, bounded rank-width, exact amplitude in $O(m)$. There is no sense in which the cost of that family scales with its magic. The slogan "none removes the exponent in magic" is simply false on this class, and the report's own Section 2.6 endorses the mechanism.

What the steelman **cannot** be stretched to is R1: bounded rank-width is a *promise* on the input. Drop the promise and rank-width is $\Theta(n)$ on generic instances, and de Colnet's $2^{O(k)}$ becomes $2^{\Theta(n)}$, i.e. back to the wall. So the steelman removes the magic exponent *on a class*, not *in general*. That is exactly the corrected statement of Section 5: the worst-case exponential is mandatory; its argument is movable.

The precise step where universality (really, $\#\mathrm{P}$-hardness) bites the steelman: the moment you ask for the *unpromised* worst case, the structural parameter you were exploiting is itself forced to be $\Omega(n)$ on the hard instances, because if it were bounded there you would have a poly-time exact strong simulator and collapse $\mathrm{PH}$. The wall is not "you cannot find structure"; it is "the structural parameter must diverge precisely on the instances that encode $\#\mathrm{P}$-hard counting."

---

## 8. Open questions: what would actually settle this

The *wall* (Theorem A, Theorem B) is settled, conditional only on the universally-believed $\mathrm{P}\ne\#\mathrm{P}$ and non-collapse of $\mathrm{PH}$. What is genuinely open is the *frontier*, i.e. which structures are tractable, and the report should frame its honesty there rather than at the wall.

1. **Is there a structural parameter strictly more general than rank-width for the $\mathrm{Clifford}+T$ Gauss sum?** de Colnet's rank-width is the current frontier; whether a more general width (or a non-graph parameter) enlarges the tractable class is open. Any such parameter must be $\Omega(n)$ on $\#\mathrm{P}$-hard instances.
2. **Does the poly-size LIMDD class strictly dominate bounded stabilizer rank?** This needs an explicit family with provably exponential stabilizer rank, itself a long-standing open problem (Section 2.3).
3. **Where exactly do the three tractability axes (degree/hitting-set, rank-width, ideal regularity) overlap and where does each strictly win?** The report asserts they are distinct (Section 2.6, Section 4.2). The star (hitting set $1$, but its interaction graph is a star, also low rank-width) and matching (hitting set $m$, but bounded rank-width) witnesses show rank-width can win where hitting set does not; a family that is easy for hitting set but hard for rank-width, or vice versa in the other direction, would sharpen the separation. This is a concrete, computable next step.
4. **Engineering, not complexity: would QF implementing (a) phase-polynomial compression for the diagonal fragment and (b) a rank-width/bucket-elimination path change which families QF can actually do?** Yes, by the witnesses; the diagonal fragment alone turns $2^t$ into $O(1)$ on [demo (4)](magic-exponent-audit-demos.wls). This is the actionable item the audit surfaces for the framework.

What would *not* settle anything: producing yet another structured family that is easy. The space of "easy because structured" families is unbounded and already well-populated. The only thing that would overturn the wall is a universal poly-time exact strong simulator, which would collapse $\mathrm{PH}$, which is why nobody expects it.

---

## 9. Verification

All claims marked with a demo reference are reproduced in a clean kernel by [`magic-exponent-audit-demos.wls`](magic-exponent-audit-demos.wls), run with `wolframscript`, no `Quiet`, no `Print`. Brute force over $\mathbb{F}_2^n$ is the ground truth; the live QF cross-check loads the working-tree paclet (`PacletDirectoryLoad`).

| Demo | Claim | Result |
|---|---|---|
| (1) | $\mathrm{gap}(f)=2\#\{f=0\}-2^n$ on a random cubic | `True` |
| (1) | $\mathrm{gap}(x_1x_2)=2$ (degree $2$), $\mathrm{gap}(x_1x_2x_3)=6$ (one $CCZ$) | `{2,2}`, `{6,6}` |
| (2) | star $f_m$: $\mathrm{gap}=2^{2m}+2^m$ (brute vs $O(m)$ hitting-set-$1$ split), $m=1..5$ | $6,20,72,272,1056$; `True` |
| (3) | matching $g_m$: $\mathrm{gap}=6^m$ (brute vs factorization), $m=1..3$ | $6,36,216$; `True` |
| (4) | $T^t\lvert+\rangle$: QF frame components $=2^t$ yet fidelity to dense $=1$, $t=1..6$ | $2,4,8,16,32,64$; all fidelity $1$; `True` |
| (4) | true $O(1)$ closed form $\{1,e^{i t\pi/4}\}/\sqrt2$ equals dense up to phase | `True` |
| (5) | $(e^{i\pi/4}X)^2=iI\ne I$ (no $+1$ eigenspace, stabilizes nothing) | `True` |
| (5) | $1+x^4\equiv(1+x)^4\pmod 2$ ($\omega$'s cyclotomic structure collapses) | `(1+x)^4` |

The star and matching demos are the load-bearing ones: each exhibits a family with **unbounded $T$-count whose exact amplitude is computed in polynomial time**, which is precisely the configuration the slogan "none removes the exponent in magic" forbids and the report's own Sections 2.6 and 2.7 permit. The QF demo shows the framework paying the $T$-count exponent where the intrinsic cost is $O(1)$, the in-kernel face of the same conflation.

---

## A note for a physicist

The quantity in dispute is a Feynman sum-over-paths amplitude written as a finite Gauss sum $\sum_v \omega^{\varphi(v)}$, $\omega=e^{i\pi/4}$, with $\varphi$ a cubic phase over the two-element field. For real (Toffoli-Hadamard) circuits it is a signed count of the zeros of a cubic Boolean polynomial, an object number theorists call the point count of a variety. The report claims that "magic" (the $T$-gate content) sits irremovably in the exponent of the simulation cost, and that universality is what makes this so. The audit's finding is that this is two true statements fused into one misleading one. The true, unconditional statement is that the *worst-case* count is $\#\mathrm{P}$-hard, a combinatorial fact that has nothing to do with universality and that no clever representation can dodge: there is no single algorithm fast on all instances. The false statement is that the cost necessarily scales with the *magic* of a particular circuit. It does not. Just as the partition function of a graphical model is cheap when the interaction graph is a tree and expensive when it is an expander, the cost of this Gauss sum is governed by the *connectivity of the phase polynomial*, not by how much non-stabilizerness the state carries. The two verified witnesses make this physical: a "star" of cubic couplings sharing one hub, and a "matching" of decoupled cubic triples, both carry unbounded magic and both have amplitudes computable in time linear in the number of gates, because in one case a single conditioning variable cuts all the couplings and in the other the system factorizes into independent pieces. The deep wall is real and worth stating: it is the statement that you cannot make *every* circuit cheap, on pain of collapsing the polynomial hierarchy. The shallow slogan, that magic itself is the exponent, is what the report's own better sentence already corrects when it says the cost should be "governed by something smaller than the raw $T$-count."

The deliverable is one Markdown file, `Magic-Exponent-Impossibility-Audit.md`, plus the verification script `magic-exponent-audit-demos.wls`; no kernel code was changed.
