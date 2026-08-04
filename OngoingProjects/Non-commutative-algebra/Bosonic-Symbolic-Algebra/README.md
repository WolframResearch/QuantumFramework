# Symbolic Bosonic Operator Algebra — Project Plan

**Owner:** Bruno Tenorio (lead) · **Editor / co-author:** Mads Bahrami
**Two deliverables:** (A) extend `Wolfram`QuantumFramework`SecondQuantization`` paclet, (B) ship arXiv preprint by **2026-07-31** + *J. Phys. A* submission by **2026-08-15**.
**Working envelope:** ~1 hr/day for 14 weeks (Apr 27 – Aug 7, 2026).

---

## 1. Goal

Make the Wolfram Quantum Framework the only environment where you can:

- Derive a bosonic operator identity symbolically (BCH / Zassenhaus / normal ordering),
- Lift the resulting truncated power series into a **closed-form** generating function,
- Compose with displacement / squeeze / Bogoliubov rules,
- Solve the resulting ODEs in the same kernel.

Then publish that workflow as a citable foundational paper.

---

## 2. What the external symbolic landscape contains (5 facts)

Verified by reading the source papers, repos, and documentation directly.

1. **No Python or Julia package implements Magnus expansion.** Symdyn ([arXiv:2503.22061](https://arxiv.org/abs/2503.22061), Sec. III) explicitly considered Magnus and chose Wei-Norman instead. SymPy, pyBoLaNO, SymPT, OpenFermion, QuTiP, and scqubits do not list it.
2. **No Python or Julia package recovers closed-form generating functions from truncated BCH or Zassenhaus power series.** SymPy quantum and Symdyn produce the symbolic polynomial coefficients but neither offers a sequence-to-generating-function lift; the round-trip ends at the polynomial.
3. **Symbolic composition of `D(α) D(β)` and `S(ξ₁) S(ξ₂)` from BCH primitives is absent.** [SymPy quantum](https://docs.sympy.org/latest/modules/physics/quantum/) can produce `D(α) D(β)` only with manual user intervention; pyBoLaNO, Symdyn, SymPT, and OpenFermion do not expose displacement / squeeze composition rules at all.
4. **Symdyn, pyBoLaNO, and SymPT each stop at symbolic derivation.** Symdyn ([arXiv:2503.22061](https://arxiv.org/abs/2503.22061)) derives Wei-Norman ODEs but does not solve them; pyBoLaNO ([arXiv:2501.01603](https://arxiv.org/abs/2501.01603)) derives Lindblad EOMs `d⟨A⟩/dt` but does not solve them; SymPT ([arXiv:2412.10240](https://arxiv.org/abs/2412.10240)) derives effective Hamiltonians but does not simulate dynamics. Each requires hand-off to a separate numerical environment (typically SciPy or QuTiP).
5. **`qutip-symbolic` is inactive.** The 82-page QuTiP 5 paper ([arXiv:2412.04705](https://arxiv.org/abs/2412.04705)) does not mention it once; QuTiP remains a purely numerical framework. v6 plans for "symbolic representation" are not yet shipped.

---

## 3. Existing material

| Source | URL / path | Reuse target |
|---|---|---|
| Staff Pick #6 — non-commutative algebra, BCH, Zassenhaus, bosonic algebra | https://community.wolfram.com/groups/-/m/t/3533713 | Paper §A, §B (lead), §C displacement subsection |
| Antinormal.md — `SqueezeOperator2[ξ, Ordering -> Normal | Weak | Antinormal]` | local file in this folder | Paper §C squeeze subsection + paclet `SqueezeOperator` extension |
| test-2nd-quantization.md — numerical verification of the three squeeze orderings | local file in this folder | Paper §C verification table |
| WL non-commutative algebra docs (v14.3) | https://reference.wolfram.com/language/guide/NoncommutativeAlgebra.html | Paper §A reference |
| `BakerCampbellHausdorffTerms` resource function | https://resources.wolframcloud.com/FunctionRepository/resources/BakerCampbellHausdorffTerms/ | Paper §B + paclet wrapper |
| Van Brunt & Visser, *J. Phys. A* 48, 225207 (2015) — anchor paper | https://arxiv.org/abs/1501.02506 | Paper §B reproduces Eq. 8 + Eq. 28 mechanically |

---

## 4. Competitor packages (run these on the §D benchmark)

| Tool | arXiv | Repo | Does | Does **not** |
|---|---|---|---|---|
| **Symdyn** | [2503.22061](https://arxiv.org/abs/2503.22061) | https://gitlab.com/VolodyaCO/dynamics-of-sun-systems | Symbolic Wei-Norman (SU(N), su(1,1), arbitrary closed Lie algebras up to L=15); exact BCH | **Does not solve** ODEs; no Magnus; no normal ordering; no D/S composition |
| **pyBoLaNO** | [2501.01603](https://arxiv.org/abs/2501.01603) | https://github.com/hendry24/pyBoLaNO | Bosonic normal ordering via Blasiak — **10–1000× faster than SymPy**; symbolic Lindblad EOM `d⟨A⟩/dt` | bosonic only; no BCH; no Wei-Norman; no Magnus; **does not solve** EOMs |
| **SymPT** | [2412.10240](https://arxiv.org/abs/2412.10240) | https://github.com/qcode-uni-a/SymPT | Symbolic Schrieffer-Wolff to arbitrary order (cQED dispersive, cross-Kerr) | no BCH; no Wei-Norman; no Magnus; no normal ordering; no dynamics |
| **OpenFermion** | [1710.07629](https://arxiv.org/abs/1710.07629) | https://github.com/quantumlib/OpenFermion | Numerical BCH; fermionic normal ordering | no symbolic BCH / Wei-Norman / Magnus; no closed-form lifting |
| **SymPy quantum** | — | https://docs.sympy.org/latest/modules/physics/quantum/ | basic commutator algebra; basic normal ordering | no BCH; no Wei-Norman; no Magnus; no closed-form lifting |
| **QuTiP 5** | [2412.04705](https://arxiv.org/abs/2412.04705) | https://github.com/qutip/qutip | numerical only | nothing symbolic — `qutip-symbolic` is dead |
| **scqubits** | [2107.08552](https://arxiv.org/abs/2107.08552) | https://github.com/scqubits/scqubits | numerical transmon/fluxonium | nothing symbolic |
| **Bosonic Qiskit (C2QA)** | — | https://github.com/C2QA/bosonic-qiskit | qumode + qubit hybrid circuits, Wigner | no operator algebra |

**Honest concessions for the paper:** pyBoLaNO wins on raw normal-ordering speed for large moments; Symdyn covers larger Lie algebras for Wei-Norman *derivation*. WL's edge is composability + closed-form lifting + `DSolve` in one kernel.

---

## 5. Track A — `SecondQuantization` paclet additions

Seven functions. Naming convention matches the existing paclet (`AnnihilationOperator`, `FockState`, `CoherentState`, `SqueezeOperator`, `WignerRepresentation`).

| # | Function | Signature | Source material | Effort | Target |
|---|---|---|---|---|---|
| 1 | **`SqueezeOperator` extension** | `SqueezeOperator[ξ, n, Ordering -> Normal | Weak | Antinormal]` | Antinormal.md (`SqueezeOperator2` already drafted) | 1 day | W19 |
| 2 | **`DisplacementOperator`** (symbolic) | `DisplacementOperator[α, n, Ordering -> Normal | Antinormal]` returning `e^{-|α|²/2} e^{αa†} e^{-α*a}` etc. | Staff Pick #6 lines 449–462 | 1 day | W19 |
| 3 | **`BosonicNormalOrder`** | `BosonicNormalOrder[expr, {a, a†}]` — wrapper around `NonCommutativePolynomialReduce` with `commRel = {Commutator[a, a†] - 1}` baked in | Staff Pick #6 §"Quantum Optics: extra commutators" | 2 days | W20 |
| 4 | **`BosonicBCH` / `BosonicZassenhaus`** | `BosonicBCH[X, Y, n]`, returns the n-th term reduced to normal-ordered bosonic monomials | Staff Pick #6 BCH/Zassenhaus blocks | 2 days | W20 |
| 5 | **`ClosedFormBCH`** ★ lead capability | `ClosedFormBCH[X, Y, relation, OrderTo -> 12]` — applies `BosonicBCH`, extracts coefficient sequence, runs `FindGeneratingFunction` / `FindSequenceFunction`, returns closed form or `$Failed` | Staff Pick #6 lines 219–325 (Van Brunt & Visser Eq. 8 + Eq. 28 already reproduced) | 4 days | W21 |
| 6 | **`BosonicMatrixElement`** | `BosonicMatrixElement[{m, n}, DisplacementOperator[α]]`, `BosonicMatrixElement[{m, n}, SqueezeOperator[ξ]]` returning closed forms via Laguerre polynomials | Staff Pick #7 kernel + Book Ch 7 | 3 days | W22 |
| 7 | **`BogoliubovDecomposition`** | `BogoliubovDecomposition[ξ₁, ξ₂]` returning `{ξ_c, θ_c}` such that `S(ξ₁) S(ξ₂) = S(ξ_c) R(θ_c)` | new — derive symbolically in W21 | 3 days (or downgrade) | W22 |

**Documentation deliverables for the paclet:**
- One reference page per new function (no `XXXX` placeholders).
- One example notebook covering all seven, demonstrated end-to-end on the canonical Heisenberg-Weyl algebra and on `[X, Y] = uX + vY + cI`.
- Cross-link from the existing `SqueezeOperator` and `WignerRepresentation` pages.

---

## 6. Track B — arXiv paper

### 6.1 What this paper is (and is not)

**This is not a top-end original-research paper.** It is a **clean arXiv preprint plus a journal-of-record submission** whose purpose is to put the `SecondQuantization` paclet on the citation map and give users something to point at when they use the new functions from Track A.

**The model:** the QuTiP 5 paper ([arXiv:2412.04705](https://arxiv.org/abs/2412.04705)), the scqubits paper ([arXiv:2107.08552](https://arxiv.org/abs/2107.08552)), the OpenFermion paper ([arXiv:1710.07629](https://arxiv.org/abs/1710.07629)), and the Symdyn paper ([arXiv:2503.22061](https://arxiv.org/abs/2503.22061)) — each is a software / methods paper that anchored a tool's scientific presence without claiming a major theoretical advance. Each is now cited heavily; none of them existed for the QuantumFramework before.

**What this paper does:**

1. **Establishes a citable reference for the `SecondQuantization` paclet.** The framework currently has 7 Google Scholar mentions; this preprint is the first foundational artifact.
2. **Showcases the seven new symbolic functions from Track A** through worked, reproducible examples.
3. **Maps the WL workflow against the four named Python competitors honestly**, including where those competitors do something better.
4. **Reproduces** Van Brunt & Visser 2015 ([arXiv:1501.02506](https://arxiv.org/abs/1501.02506)) Eq. 8 + Eq. 28 mechanically as the demonstration of closed-form lifting — the result is theirs; the *automation* is the paper's contribution.

**What this paper does *not* claim:**

- Original new theorems in operator algebra. The §B "theorem-lite" remark is a sharper statement of what is observable in the Van Brunt & Visser examples, not a new general result.
- A faster bosonic normal-ordering algorithm than [pyBoLaNO](https://arxiv.org/abs/2501.01603). pyBoLaNO is faster on raw normal-ordering of large moments by 10–1000× — the paper concedes this in §D.
- A broader Lie-algebraic engine than [Symdyn](https://arxiv.org/abs/2503.22061) for the Wei-Norman derivation step. Symdyn handles arbitrary closed Lie algebras up to L=15; the WL workflow's edge is integrating the derivation with `DSolve` in one kernel, not algebraic breadth.

**What success looks like:**

- arXiv preprint live by **2026-07-31**.
- Cited by the next users of `SqueezeOperator`, `DisplacementOperator`, `ClosedFormBCH`, etc., not by QFT theorists writing about new operator-algebra results.
- One follow-on paper per quarter for the next 12 months (Magnus, Wei-Norman end-to-end, Schrieffer-Wolff) cites this paper as the framework reference.
- Becomes the "the paper to cite" when anyone publishes WL-based work using the SecondQuantization functions.

### 6.2 Title, length, sections, venue

**Title (working):** *Symbolic bosonic operator algebra in the Wolfram Language: BCH, Zassenhaus, and normal-ordered displacement and squeeze identities for the QuantumFramework `SecondQuantization` paclet.*

**Length:** 18–22 pp.

**Sections:** §1 Intro · §A Setup (2 pp) · §B Closed-form BCH recovery — *main demonstration* (5 pp) · §C Normal-ordered identities (5 pp) · §D Benchmark vs Python tools (2.5 pp) · §6 Conclusion + roadmap (1.5 pp).

**Venue strategy** — software / methods journals are the right target for a presence paper, *not* a top-tier physics journal that expects original results.

| Priority | Venue | Why | Notes |
|---|---|---|---|
| 1 | **SciPost Physics Codebases** (https://scipost.org/SciPostPhysCodeb) | Purpose-built for code-anchored papers; open access; transparent peer review | Submit Aug 2026 |
| 2 | **Computer Physics Communications** | Long-established software journal; historical home for Mathematica methods papers; long methods sections allowed | Fallback if SciPost stalls |
| 3 | **Quantum** (https://quantum-journal.org) | Higher prestige; open access; accepts software/tutorial papers | Fallback 2 |
| — | **arXiv** (cross-list `quant-ph` + `math-ph`) | The deliverable that actually moves the visibility needle | **Hard target: 2026-07-31** |

*J. Phys. A: Math. Theor.* is **not** the right primary target for this paper. J. Phys. A expects original mathematical results; this is a software/methods paper. Reserve J. Phys. A for the future Magnus-expansion paper.

### 6.3 14-week schedule

| Week | Dates | Deliverable | Action items |
|---|---|---|---|
| **W18** | Apr 27 – May 3 | 2-page outline locked | Read [Van Brunt & Visser 2015](https://arxiv.org/abs/1501.02506) cover-to-cover. Skim [Symdyn](https://arxiv.org/abs/2503.22061), [pyBoLaNO](https://arxiv.org/abs/2501.01603), [SymPT](https://arxiv.org/abs/2412.10240) intros + conclusions. Map every section to source material. |
| **W19** | May 4 – May 10 | Section maps + figure list | Identify reused figures from Staff Pick #6. Mads provides paper template. **Paclet items 1–2 land.** |
| **W20** | May 11 – May 17 | §A (2 pp) + §C displacement subsection (2 pp) | Convert Staff Pick #6 lines 11–143 (§A) and lines 354–505 (§C displacement) into journal voice. **Paclet items 3–4 land.** |
| **W21** | May 18 – May 24 | §C squeeze subsection (3 pp) + Bogoliubov derivation | Convert Antinormal.md squeeze table. Run `BogoliubovDecomposition` symbolically. **If not closed in 2 days, defer to Paper 2 and mark numerical-only.** **Paclet item 5 lands.** |
| **W22** | May 25 – May 31 | §B examples 1 + 2 (3 pp) | Convert Staff Pick #6 lines 219–325. Density plot becomes Figure 1. **Paclet items 6–7 land.** |
| **W23** | Jun 1 – Jun 7 | §B example 3 + theorem-lite half-page (2 pp) | Add `[X, Y] = cI` degenerate limit. Cite Bourbaki Ch I, Humphreys, Casas-Murua-Bossa. **§B complete.** |
| **W24** | Jun 8 – Jun 14 | **CUT WEEK** — re-read draft, cut 20% | No new content. Tighten captions and prose. Buffer if W21 ran over. |
| **W25** | Jun 15 – Jun 21 | Python competitors installed; problems 1–3 run | `pip install pybolano sympy openfermion`, clone Symdyn, set up SymPT. Run benchmark problems 1 (`(a+a†)^6`), 2 (`[a†ⁿ¹aᵐ¹, a†ⁿ²aᵐ²]`), 3 (BCH closed-form). **Budget 2 days for env setup alone.** |
| **W26** | Jun 22 – Jun 28 | Problems 4–6 run + §D table written | Run `D(α)D(β)` (problem 4), `S(ξ₁)S(ξ₂)` (problem 5), Wei-Norman su(1,1) (problem 6). Write honest concessions on pyBoLaNO speed and Symdyn breadth. §D drafted (2.5 pp). |
| **W27** | Jun 29 – Jul 5 | Intro + conclusion + bibliography | §1 the "every Python tool stops at step 1" pitch. §6 explicit roadmap reserving Magnus / Wei-Norman end-to-end / Schrieffer-Wolff for follow-up papers. **Full draft, 18–22 pp.** |
| **W28** | Jul 6 – Jul 12 | Mads internal review #1 | Bruno works on book / QEC Staff Picks this week. |
| **W29** | Jul 13 – Jul 19 | Revisions absorbed | Re-run benchmark cells if Mads requests. |
| **W30** | Jul 20 – Jul 26 | Second internal reviewer (Nikolay) + final polish | Submission-ready PDF. Final figure resolution. |
| **W31** | Jul 27 – Aug 7 | **arXiv 2026-07-31. SciPost Physics Codebases submission 2026-08-15.** | Cross-list quant-ph + math-ph. Deposit code (Wolfram Cloud notebook + GitHub mirror). Same week: Wolfram Community announce + LinkedIn. |

**Total: ~85 hours = ~6 hr/week = ~50 min/day.**

---

## 7. The §D benchmark — the six problems Bruno runs on every tool

| # | Problem | SymPy | pyBoLaNO | Symdyn | SymPT | WL |
|---|---|:-:|:-:|:-:|:-:|:-:|
| 1 | Normal-order `(a + a†)^6` | ✓ | ✓ (10–1000× faster) | — | — | ✓ |
| 2 | `[a†ⁿ¹aᵐ¹, a†ⁿ²aᵐ²]` normal-ordered, n_i, m_i ≤ 3 | ✓ | partial | — | — | ✓ |
| 3 | BCH of `[-α*a, αa†]` to 10 orders, recover `e^{-½|α|²}` | poly only | — | — | — | **✓ + closed form** |
| 4 | `D(α) D(β)` composition | manual | — | — | — | **✓ automatic** |
| 5 | `S(ξ₁) S(ξ₂)` Bogoliubov decomposition | — | — | — | — | **✓ WL only** |
| 6 | Wei-Norman for su(1,1) parametric H — derive **and** solve | — | — | derives only | — | **derives + solves** |

---

## 8. Risks (4 that matter)

| Risk | Likelihood | Mitigation |
|---|---|---|
| `BogoliubovDecomposition` (W21) doesn't close symbolically | Medium | Pre-decide: if not closed in 2 days, defer to Paper 2; paper has 5 other §C results. |
| Python tool install hell (W25) | High | Fresh `conda` env per package. 2 days budget. If SymPT refuses, drop from runnable benchmark — keep textual citation only. |
| Reviewer rejects on novelty grounds | Low–Medium | W18 literature review: search every paper citing [Van Brunt & Visser 2015](https://arxiv.org/abs/1501.02506). Confirm no prior `FindGeneratingFunction` lift. |
| Book deadline slips, paper has to pause | Medium | Book wins absolutely. Schedule has ~2 mo buffer to *J. Phys. A* review window. |

---

## 9. Success criteria (must-have)

- [ ] arXiv preprint live by **2026-07-31** (cross-listed `quant-ph` + `math-ph`)
- [ ] Software/methods journal submission by **2026-08-15** (SciPost Physics Codebases primary; CPC or Quantum as fallback)
- [ ] Six-problem benchmark runs honestly — including where WL loses
- [ ] All four Python competitors cited; ≥3 actually run
- [ ] Van Brunt & Visser 2015 ([arXiv:1501.02506](https://arxiv.org/abs/1501.02506)) Eq. 8 + Eq. 28 reproduced mechanically; code public
- [ ] Paclet items 1–6 shipped; item 7 (`BogoliubovDecomposition`) shipped or explicitly deferred
- [ ] Every new paclet function has a reference page (no `XXXX`)
- [ ] Paper title and abstract explicitly identify the artifact as a paclet/methods paper, not an original-results paper
- [ ] Book deadline (Wolfram Media, Summer 2026) not compromised

## 10. Should-have (drop without regret if needed)

- [ ] §B theorem-lite paragraph on closed-form-recoverable commutator classes
- [ ] `BogoliubovDecomposition` symbolic (vs. numerical only)
- [ ] Two internal review rounds (Mads + Nikolay)
- [ ] Companion Wolfram Community post on the day arXiv goes live

---

## 11. What this paper deliberately defers (protects against scope creep)

Reviewer-proof: §6 of the paper explicitly lists these as future work, so reviewers cannot demand them in revision.

- Magnus expansion for bosonic Hamiltonians → **Paper 2**
- Wei-Norman end-to-end for su(1,1) symbolically → **Paper 2**
- Schrieffer-Wolff automation → **Paper 3** (would compete with [SymPT](https://arxiv.org/abs/2412.10240))
- Floquet / stroboscopic Hamiltonian synthesis → future
- Tie-in to quantum Fisher information / metrology → **separate two-mode phase-space paper**

---

## 12. One-line summary

**Ship 7 new symbolic-bosonic functions into the `SecondQuantization` paclet between May and June 2026, then post a clean arXiv presence-paper by July 31 (and a SciPost Physics Codebases / CPC submission by August 15) that showcases those functions, benchmarks honestly against the four named Python competitors, and gives the QuantumFramework its first citable foundational artifact in the symbolic-operator-algebra space.**
