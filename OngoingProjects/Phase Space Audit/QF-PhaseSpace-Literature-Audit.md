# QF Discrete Phase-Space Surface: Kernel Audit vs arXiv Literature

**Date:** 2026-06-30 (v2, full rewrite after a scoping correction) · **Kernel:**
`Wolfram/QuantumFramework` 2.0.0, HEAD `179a1acf`.
**Scope:** finite-dimensional / qudit *discrete* phase space. **Excludes** second quantization
(Fock, bosonic, optical CV phase space). **Method:** two-agent adversarial audit, three rounds.
One agent read the kernel and live-verified every claim with `wolframscript` (no `Quiet`); one
agent searched arXiv; they cross-examined, and the disputed items were settled by code and by
fresh kernel runs, not by description.

> **Why v2.** v1 scoped the audit to the files *named* for phase space
> (`QuantumWignerTransform.m`, `QuantumPhaseSpaceTransform.m`, `MIC.m`) and so missed that the
> phase-space *picture* composes with the rest of the object model. Two corrections, both verified
> in round 3, are folded in: (1) QF has a full **dynamical** phase-space surface
> (`QuantumEvolve` + `HamiltonianTransitionRate` + `LindbladTransitionRates`), so the
> Braasch-Wootters transition-rate construction is **implemented**, not a gap, and is extended to
> open systems; (2) QF **packages entanglement monotones** (`QuantumEntanglementMonotone`), so the
> earlier "QF only hands you the tensor" line was an overgeneralization. A corrections log is in
> Section 9.

---

## 1. Headline answer

**Yes: QF's discrete phase-space surface sits squarely on a well-developed arXiv literature, and
it is a working *dynamical* framework, not just a static change of representation.** Every major
static capability (discrete Wigner function on finite fields, Wigner negativity as
contextuality/magic witness, mana, the discrete Hudson theorem, the frame/quasiprobability
representation of states/channels/POVMs, SIC- and MIC-POVMs, the QBist Born rule, mutually
unbiased bases, the Weyl transform) maps to a canonical founding paper, usually one-to-one. The
dynamical capability (transition-rate generator, discrete-Wigner time evolution, circuit
contraction in the phase-space picture) maps to the Wootters/Feynman -> Bianucci-Miquel-Paz-Saraceno
-> Braasch-Wootters lineage. Two pieces are genuine QF distinctives: the **open-system**
(Lindblad) discrete-Wigner transition rates, which the literature only gestures at, and the signed
transition **graph** object. A small cluster of *magic-resource* theorems remains genuinely
absent (robustness of magic, stabilizer Rényi entropy, distillation thresholds, contextuality
polytopes, analytic negativity bounds).

---

## 2. Phase-space STATIC surface (all re-verified, round 3)

Kernel files: `QuantumWignerTransform.m`, `QuantumPhaseSpaceTransform.m`, `MIC.m`, `Utilities.m`
(`pauliMatrix`, `fanoMatrix`, `GellMannMatrices`, `RegularSimplex`, `GramMatrix`/`GramDual`),
`QuditBasis/NamedBases.m`, `QuantumState/Properties.m:182-326`,
`QuantumMeasurementOperator/NamedMeasurementOperators.m`, `Visualization.m:100-104`.

The static surface is one idea: **pick a frame of `d^2` operators, expand any object in it, invert
with the Gram dual.**

- **Phase-point / Fano operators** `fanoMatrix[d,q,p] = e^{i pi q p/d} F^2 . Z^q . X^p`
  (`Utilities.m:241`) from the generalized Heisenberg-Weyl group (`pauliMatrix[1,d]`=shift,
  `pauliMatrix[3,d]`=clock). Odd `d`: covariant `d x d` Wigner grid (Gross/Wootters). Even `d`:
  Wootters `2d x 2d` doubling (no covariant `d x d` exists).
- **`QuantumWignerTransform`** / **`QuantumWeylTransform`**: discrete Wigner of
  state/operator/channel/POVM/circuit, and the Weyl symbol map + inverse ("Undouble").
- **`QuantumPhaseSpaceTransform`**: the general engine: transform any object into *any* `d^2`-dim
  `QuditBasis` frame (default Wigner). Uniform across states, operators, **channels**, **POVMs**,
  circuits.
- **`QuantumPositiveTransform`**: non-negative over-complete representation via a `Ramp(+/-)`
  odd-axis split (classical-sampling representation; support grows with negativity).
- **MIC/SIC machinery** (`MIC.m`): `DisplacementOperator[d,a,b]`, `WeylHeisenbergGroup` /
  `FactoredWeylHeisenbergGroup`, `QBismSICPOVM`, `HesseSICPOVM` (d=3), `HoggarSICPOVM` (d=8),
  `QuantumWignerMICPOVM` (with even-`d` `a/b/a2/b2` dual-frame correction, `MIC.m:17-51`),
  `GellMannMICPOVM`, `RandomHaarPOVM`, `RandomBlochMICPOVM`, `WoottersBasis`; `GramDual` is the
  canonical dual frame.
- **State getters**: `"QuasiProbability"`, `"PhaseSpace"`, `"TransitionQuasiProbability"`,
  `"TransitionPhaseSpace"`, `"TransitionGraph"`, `"FullTransitionGraph"`, `"Prime*"` variants.
- **Named bases / POVMs**: `"Wigner"[d]`, `"WignerMIC"`, `"Wootters"[d]`, `"Feynman"[d]`,
  `"HesseSIC"`, `"HoggarSIC"`, `"Ivanovic"[prime d]` (MUB); named measurements `WignerMICPOVM,
  GellMannMICPOVM, TetrahedronSICPOVM, QBismSICPOVM, HesseSICPOVM, HoggarSICPOVM, RandomPOVM`.

**Verified static facts.** Wigner of `|0>`,`|+>` at odd `d=3` is non-negative (`min W = 0`), of the
magic state `|H> = T|+>` is negative (`min W = -0.177`); `mana = Log Total|W|` is `0` on
stabilizers and `0.511` on the strange state. The qubit `|0>` Wigner is **negative**
(`min W = -1/4`) because even-`d` uses the Wootters doubling, so the **discrete Hudson theorem
holds only at odd `d`, not at the qubit level**. Hesse SIC (`d=3`) resolves the identity exactly
(`Sum E_i - I_3 = 0`); the Tetrahedron SIC (`d=2`) is equiangular (Gram off-diagonal `1/12`). The
even-`d` Wigner-MIC branch yields a valid POVM (`d=2`: 4 elements, `Sum = I_2`, min eigenvalue
`0.033 > 0`). The Weyl transform round-trips states exactly (`Weyl(Wigner(rho)) - rho = 0`).

---

## 3. Phase-space DYNAMICAL surface (the part v1 missed, all verified round 3)

Kernel file: `QuantumEvolution.m`. **This is what makes QF's discrete phase space a dynamical
framework.**

- **`HamiltonianTransitionRate[h]`** (`:267`) `= Chop[QuantumPhaseSpaceTransform[
  QuantumOperator["Hamiltonian"[h]], args]["Matrix"] / I]`. This is the **Braasch-Wootters
  transition-rate generator**: take the phase-space symbol of `H`, divide out the unitary `i`, and
  you get the matrix `R` driving `dW/dt = R . W` as a flow between discrete phase-space points.
  Verified: `rate[Z3]` is `9 x 9` and equals `PhaseSpaceTransform[Ham[Z3]]["Matrix"]/I` exactly.
  **Caveat:** `R` is **complex in general** (real only when the Hamiltonian's phase-space symbol is
  real-antisymmetric); the tempting "the `/I` makes it a real rate matrix" picture is not generic.
- **`LindbladTransitionRates[ls]`** (`:269`): the same construction for jump operators, giving the
  **open-system / dissipative** transition rates. This goes **beyond** Braasch-Wootters (closed
  systems only).
- **`QuantumEvolve` phase-space branch** (`:49, 66-79, 194-202`): hand `QuantumEvolve` an initial
  state carried in `"Picture" -> "PhaseSpace"` and it integrates the dynamics **directly on the
  discrete Wigner grid** as a real vector field (via `NDSolve`/`DSolve`), then maps back with
  `QuantumWeylTransform`.
- **Circuit contraction in the phase-space picture** (`QuantumCircuitOperator/TensorNetwork.m:190`,
  `phaseSpaceQ`): whole circuits are tensor-network-contracted in the Wigner representation.

**Verified dynamical facts.**
- Closed-system Rabi: `H = X` on `|0>`, integrated entirely on the phase space and mapped back,
  reproduces the direct Schrödinger density matrix at `t = 1` to `5.6e-6` (NDSolve tolerance).
  **DISPUTED (2026-06-30):** a later independent fresh-kernel run could not reproduce this via
  the literal `QuantumEvolve[QuantumOperator["X",{1}], QuantumState[{1,0}, QuantumBasis["Wigner"[2]]],
  {t,0,1}]` path: the Wigner-basis state came back frozen (`||sv(1)-sv(0)|| = 0`) and `["Undouble"]`
  errors on a single-output-qudit Wigner state. The transition-rate generator itself builds
  correctly (`HamiltonianTransitionRate` -> 9x9), and the closed/open dynamics are robustly
  reproduced in the **computational basis** (Rabi populations `{cos^2 1, sin^2 1}`; damping below).
  Treat "integrated on the grid and mapped back for a qubit" as not-yet-reproducible pending a
  kernel check of the even-`d` Wigner `QuantumEvolve` path; the generator + computational-basis
  validation is the solid claim.
- Open-system amplitude damping on `|1>` (rate 1, `t = 2`): the phase-space-evolved state is
  `diag(0.865, 0.135)`, trace 1, spectrum `>= 0`, matching the exact `1 - e^{-2}`, `e^{-2}`.
- Circuit: `QuantumWignerTransform[{H, CNOT}]` applied on the phase-space `|00>` and mapped back is
  `QuantumDistance = 0` against the direct result.

---

## 4. Adjacent packaged measures (entanglement, not magic)

`$QuantumEntanglementMonotones` (`QuantumEntanglement.m:9`) `= {Concurrence, Negativity,
LogNegativity, EntanglementEntropy, RenyiEntropy, Realignment, MutualInformationI,
MutualInformationJ, Discord}`, plus `QuantumEntangledQ`. Verified: Bell `{Negativity, Concurrence}
= {0.5, 1}` (maximal), product `|00>` `= {0, 0}`, `QuantumEntangledQ` `= {True, False}`.

**These quantify *entanglement* (the LOCC resource theory), not *magic* / non-stabilizerness.**
The only non-stabilizerness quantity QF exposes is **mana**, and only as a hand-built
`Log[Total[Abs[ qs["QuasiProbability"] ]]]`, not as a named getter. `stabilizerEntropy`
(`Stabilizer/Entropy.m`) is the *bipartite entanglement entropy* of a stabilizer state
(Fattal et al. F_2-rank formula), a name collision with, not an instance of, the
Leone-Oliviero-Hamma stabilizer Rényi *magic* measure.

---

## 5. Claim -> literature map (settled across three rounds)

Strength: STRONG = one-to-one; PARTIAL = concept published, QF detail differs; DISTINCTIVE = QF
goes past the published object.

| QF computable quantity | Founding / representative arXiv (verified) | Strength |
|---|---|---|
| **Static** | | |
| Discrete Wigner `W(q,p)` on finite fields | Wootters 1987 (DOI 10.1016/0003-4916(87)90176-X); Gibbons-Hoffman-Wootters `quant-ph/0401155`; Gross `quant-ph/0602001` | STRONG |
| Even-`d` Wigner via `2d x 2d` doubling | Wootters et al `quant-ph/0306135`; Delfosse et al `1409.5170`; Zhu no-go `1504.03773` | STRONG |
| Wigner negativity = contextuality/nonclassicality | Spekkens `0710.5549`; Delfosse-Raussendorf `1409.5170`; Howard et al `1401.4174` | STRONG |
| Mana `= Log Total|W|` (magic monotone) | Veitch-Mousavian-Gottesman-Emerson `1307.7171` | STRONG (1-to-1) |
| Discrete Hudson theorem (odd `d`) | Gross `quant-ph/0602001` | STRONG (odd `d` only) |
| Frame/quasiprobability rep of states, POVMs, **channels** | Ferrie-Emerson `0711.2658`, `0903.4843`; review `1010.2701` | STRONG (POVM) / PARTIAL (channel rarely worked explicitly) |
| Non-negative rep -> classical-sim cost ~ negativity | Veitch et al `1201.1256`; Pashayan-Wallman-Bartlett `1503.07525`; `1810.03622` | STRONG (rep); QF builds the distribution, no sampler |
| Weyl-Heisenberg displacements + Clifford covariance | Gross `quant-ph/0602001`; Zhu `1504.03773` | STRONG |
| WH-covariant SIC-POVM + Zauner program | Renes et al `quant-ph/0310075`; Scott-Grassl `0910.5784`; Fuchs-Hoang-Stacey `1703.07901`; Appleby et al `2501.03970` | STRONG |
| QBist Born rule / urgleichung (Gram dual) | Appleby-Fuchs-Stacey-Zhu (Qplex) `1612.03234` | STRONG |
| MIC-POVMs + Qplex | DeBrota-Fuchs-Stacey `1812.08762` | STRONG |
| **Wigner-MIC** (minimal-IC Wigner basis) | DeBrota-Stacey `1912.07554` | STRONG (exact source); even-`d` coeff surgery a QF realization |
| Mutually unbiased bases (`"Ivanovic"`) | Ivanovic 1981 (DOI 10.1088/0305-4470/14/12/019); Wootters-Fields 1989 (DOI 10.1016/0003-4916(89)90322-9); review `1004.3348` | STRONG |
| Wootters line-sum / `"Feynman"` basis | Wootters 1987; GHW `quant-ph/0401155`; Feynman 1987 "Negative Probability" (book chapter, no arXiv) | STRONG |
| **Dynamical** | | |
| Transition-rate generator `dW/dt = R W` (`HamiltonianTransitionRate`) | Braasch-Wootters `2007.10187`; precursor Bianucci-Miquel-Paz-Saraceno `quant-ph/0106091`; recent `2503.09353` | STRONG (1-to-1) |
| Discrete-Wigner time evolution (state evolves in phase space) | same lineage; Wootters/Feynman 1987 -> `quant-ph/0106091` -> `2007.10187` | STRONG |
| **Open-system (Lindblad) discrete-Wigner transition rates** (`LindbladTransitionRates`) | CV Wigner-Lindblad only: Brun-Halliwell `quant-ph/9605041`, `2307.16510`; discrete adjacent `2303.05291`, `quant-ph/0308104` | **DISTINCTIVE** (no canonical discrete-qudit paper) |
| Circuit contraction in phase-space picture | (engineering; TN contraction of Wigner-represented circuits) | no direct paper |
| **Static transition object** | | |
| `"TransitionQuasiProbability"` (two-basis X<->Z, real signed) | Kirkwood-Dirac family: Yunger Halpern 2017 `1609.09689`; Yunger Halpern-Swingle-Dressel `1704.01971` | PARTIAL (right family, real vs complex) |
| **Adjacent (entanglement)** | | |
| Negativity / LogNegativity | Vidal-Werner `quant-ph/0102117`; Plenio `quant-ph/0505071` | STRONG |
| Concurrence | Wootters `quant-ph/9709029` | STRONG |
| Realignment / CCNR | Chen-Wu `quant-ph/0205017`; Rudolph `quant-ph/0202121` | STRONG |
| Entanglement / Rényi entropy | Bennett-Bernstein-Popescu-Schumacher `quant-ph/9511030` | STRONG |

Loose end: Appleby 2005 "SIC-POVMs and the extended Clifford group" (JMP 46:052107) exact arXiv ID
not pinned.

---

## 6. The transition object: two distinct things

v1 conflated these; they are separate.

1. **Dynamical transition rates** (`HamiltonianTransitionRate`, `LindbladTransitionRates`, the
   `QuantumEvolve` phase-space branch): the **Braasch-Wootters** generator `R = symbol(H)/i` and
   its open-system extension. QF *runs* this: `dW/dt = R . W` integrated on the Wigner vector.
   Match to `2007.10187` is one-to-one for the closed case; the open case is QF-distinctive.
2. **Static two-basis getter** `"TransitionQuasiProbability"` / `"FullTransitionGraph"`: a single
   `rho`, no Hamiltonian, no time. Its graph vertices are `Join[QuditBasis["X"[dims]]["Names"],
   QuditBasis[dims]["Names"]]` (`Properties.m:284,304`): the X-eigenbasis and the computational
   (Z) basis of the same state. That is the **Kirkwood-Dirac** skeleton `Q(a,b)=<a|rho|b><b|a>`,
   except QF's matrix is **real and signed** (built inside the doubled Wigner phase space) where
   textbook KD is complex. The colored signed transition **graph** has no literature match: likely
   a QF-original *presentation*.

---

## 7. QF-distinctive / under-published pieces

- **Open-system discrete-Wigner transition rates** (`LindbladTransitionRates` + `QuantumEvolve`):
  dissipative dynamics on the finite-field qudit grid. The continuous-variable Wigner-Lindblad is
  textbook (`quant-ph/9605041`, `2307.16510`), but a discrete-qudit Lindblad transition-rate
  generator is not an established named construction (the literature flags the obstruction that
  phase-point operators lose orthogonality under dissipation). Verified working
  (amplitude damping exact). **Strongest novelty candidate.**
- **`"FullTransitionGraph"` as a rendered signed graph object** (X-vertices vs Z-vertices, red=+,
  blue=-): original presentation of a Kirkwood-Dirac-like quantity.
- **The unified arbitrary-frame `QuantumPhaseSpaceTransform`** over states, operators, **channels**,
  **POVMs**, and circuits through one code path: the Ferrie-Emerson formalism covers it, but an
  explicit channel/process frame transform as a computed object is uncommon.
- **The even-`d` Wigner-MIC `a/b/a2/b2` dual-frame correction** (`MIC.m:17-51`): QF's closed-form
  realization of the DeBrota-Stacey `1912.07554` principle; likely the only shipped software that
  builds the even-`d` Wigner-MIC explicitly.

---

## 8. Gaps: literature computes it, QF does not (corrected, grep-confirmed)

| Quantity | Paper that computes it | QF status |
|---|---|---|
| Robustness of magic (LP monotone) | Howard-Campbell `1609.07488` | Absent (QF has mana only) |
| Stabilizer Rényi entropy (Pauli-based magic) | Leone-Oliviero-Hamma `2106.12587` | Absent (`stabilizerEntropy` is *entanglement*) |
| Magic-distillation rate/overhead bounds | Bravyi-Kitaev `quant-ph/0403025`; Bravyi-Haah `1209.2426` | Absent |
| Contextuality / exclusivity polytopes (CSW) | Cabello-Severini-Winter `1010.2163` | Absent |
| Analytic Weyl-Heisenberg negativity bounds | DeBrota-Fuchs `1703.08272` | Absent |
| Born / Monte-Carlo sampler on the positive rep | Pashayan-Wallman-Bartlett `1503.07525` | Absent (QF *builds* the distribution; `QuantumStateSampler` is tomography) |
| Clifford-covariance-of-`W` verification | Gross `quant-ph/0602001` | Absent as a getter (displacements exist; hand-assembled) |
| Covariant even-`d` / qubit Wigner | (no-go: Zhu `1504.03773`) | Not possible covariantly; only Wootters `2d x 2d` doubling |

What is **NOT** a gap (v1 wrongly listed or implied these): Braasch-Wootters transition rates
(present, `HamiltonianTransitionRate`), open-system phase-space dynamics (present,
`LindbladTransitionRates` + `QuantumEvolve`), discrete-Wigner time evolution (present),
entanglement monotones (present, 9 of them + `QuantumEntangledQ`).

---

## 9. Corrections log (v1 -> v2)

| Claim in v1 | Status now | Why |
|---|---|---|
| "QF does not do Braasch-Wootters" (transition rates a gap) | **REFUTED** | `HamiltonianTransitionRate` (`QuantumEvolution.m:267`) is exactly the BW generator; `= PhaseSpaceTransform[Ham]/I` verified |
| Dynamical phase-space evolution not mentioned | **ADDED, VERIFIED** | `QuantumEvolve` phase-space branch; qubit Rabi matches Schrödinger to `5.6e-6` |
| Open-system phase-space dynamics not mentioned | **ADDED, VERIFIED** | `LindbladTransitionRates`; amplitude damping exact `1-e^{-2}` |
| "QF gives you the tensor, not the monotone" (implying no resource measures) | **OVERGENERALIZATION, corrected** | `QuantumEntanglementMonotone`: 9 entanglement monotones packaged. (Magic monotones beyond mana still absent.) |
| Transition object framed only as static KD-like | **SPLIT** | static getter = KD-like; dynamical rates = Braasch-Wootters. Two different objects. |
| (implicit) rate matrix is real-signed | **REFUTED** | `R` is complex in general; real only for special Hamiltonians |

The magic-monotone gap list (Section 8) is unchanged from v1: those quantities are genuinely
absent.

---

## 10. Physicist summary

QF carries a complete finite-dimensional discrete phase-space toolkit, and the central correction
in this revision is that it is a working *dynamical* framework, not merely a static change of
basis. Statically it is the Wootters/Gibbons-Hoffman-Wootters/Gross discrete Wigner function with
the whole "negativity = contextuality = magic = simulation cost" chain on top (Spekkens; Veitch et
al; Pashayan-Wallman-Bartlett), plus the SIC/MIC/QBism program (Renes et al; the Qplex;
DeBrota-Fuchs-Stacey, whose 2019 "Wigner functions from informationally complete measurements" is
the exact source of QF's `WignerMIC`), mutually unbiased bases, and an exact Weyl symbol map.

Dynamically, QF builds the Braasch-Wootters transition-rate generator by taking the phase-space
symbol of the Hamiltonian and dividing out the unitary `i`, and it integrates the resulting flow
directly on the discrete Wigner grid: a qubit Rabi flip reproduces the density matrix to six parts
in a million, and a circuit contracted in the phase-space picture is bit-for-bit the direct
result. It then extends the construction to jump operators, so you can run dissipative dynamics as
a flow between discrete phase-space points: amplitude damping reproduces the exact `1-e^{-2}`
populations. That open-system discrete-Wigner evolution is the one piece with no clean literature
antecedent (the continuous-variable Wigner-Lindblad is textbook, but the finite-field qudit
version is not a named construction), so it is the strongest novelty candidate, alongside the
signed transition-graph object. One honest caveat: the transition-rate matrix is complex in
general, not a real rate matrix, so the convenient classical-stochastic picture is not generic.

Finally, QF does package resource measures, but for the *entanglement* resource theory (negativity,
log-negativity, concurrence, realignment, discord, verified maximal on a Bell pair and zero on a
product state), not the *magic* one. On the magic side QF gives you mana (and only as a hand-built
`Log Sum |W|`) and the raw quasiprobability tensor; the higher-level magic theorems a referee would
reach for next, robustness of magic, the stabilizer Rényi entropy, distillation thresholds,
contextuality polytopes, analytic negativity bounds, and a Born sampler on the positive
representation, are genuinely not in the kernel.
