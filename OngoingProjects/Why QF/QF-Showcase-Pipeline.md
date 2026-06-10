# QF Showcase: Example-Design Pipeline and Decision Log

A living tracking doc for *how* the showcase examples are generated, why, and how to make them better. The deliverable is [QF-Showcase.md](QF-Showcase.md); this file is the process behind it.

---

## 1. The pipeline (reusable checklist per claim)

Run every candidate example through these gates, in order. An example that fails an early gate is rejected before any kernel time is spent.

1. **Name the unique mechanism.** What can QF do here that Qiskit / Cirq / PennyLane cannot, or can only do clumsily? Write it in one sentence. If you cannot, the claim has no example.
2. **Demonstrate a property of the object model, not a single convenience function.** A good example shows that the *whole algebra* behaves a certain way (states, operators, circuits, channels), not that one function has a nice option. This is the gate the first draft of Claim 1 failed.
3. **Contrast test.** Would a competent Qiskit user, seeing this, say "I cannot easily do that"? If the result is hand-derivable or one `import` away elsewhere, it is not a differentiator, only a feature.
4. **Audience fit, explicit.** Map the example to the four audiences (student, educator, researcher, industry). State the "so what" for each. If it lands for fewer than two, either redesign or label it as deliberately narrow.
5. **Non-triviality vs legibility.** The example must be too involved to do by hand, yet still readable in under ten lines. Trivial-but-clean is a trap: it reads as "toy."
6. **Verify in kernel.** Run it. Capture the real output. Ship only what passed. No `Quiet`, no `Print` (use natural output / `Get` in `-code`).
7. **Annotate honestly.** One line per audience, and an explicit flag for any audience the example serves weakly.

### Anti-patterns this pipeline is meant to catch
- **Function-as-claim.** Picking the one function whose uniqueness is most obvious (here, `QuantumEvolve` because `DSolve` is visible) and letting it stand in for a property that is actually pervasive.
- **Trivial instance.** Choosing the smallest case (H = X, computational initial state, no free parameter in the model) so that nothing in the example is actually parametric except an incidental variable.
- **Audience drift.** Shipping an example that only a student could love and tagging all four audiences anyway.

---

## 2. Assumptions log

Each row is an assumption baked into a first-draft example, and whether it survived scrutiny.

| # | Claim | Assumption I made | Held? | Correction |
|---|---|---|---|---|
| A1 | 1 | "Symbolic-first" is best shown by symbolic *dynamics*, so use `QuantumEvolve`. | No | Symbolic-ness is a property of the whole object algebra (operators, states, circuits, observables). Evolution is one rung, not the headline. |
| A2 | 1 | The smallest instance (evolve $\lvert0\rangle$ under $X$) is "clean and therefore good." | No | It is hand-derivable and has no free parameter in the model. Reads as a toy; hooks no researcher or industry user. |
| A3 | 1 | One example can serve all four audiences. | Partly | True only if the example *ladders* from recognizable (gate matrices) to consequential (closed-form VQE cost). A single trivial instance cannot. |
| A4 | 2 | `QuantumChannel["Operators"]` returns the Kraus list. | No | Undefined for this head; dropped that line, kept the verified density-matrix result. |
| A5 | 3 | Global `VonNeumannEntropy` shows entanglement. | No | It is $0$ for any pure state. Entanglement lives in the reduced state; used concurrence as the witness. |
| A6 | 4 | The Qiskit / IBM / AWS bridges can be verified in this environment. | No | `QuantumCircuitOperatorToQiskit` returns unevaluated (no Python/qiskit). Claim 4 anchored on verified QASM export; the hardware-submission line is shown but explicitly labeled NOT executed. |
| A7 | 2 | `G2Coherence` of a symbolic-`\[Alpha]` coherent state gives a clean closed form. | No | It expands the Fock-truncation series into a huge expression. Used numeric photon statistics instead (Fock $0$, coherent $1$, thermal $2$). |
| A8 | all | One laddered example per claim serves every audience. | Superseded | User chose **one example per audience** (4 per claim, 20 total). Applied. |
| A9 | 1 | The per-audience rebuild of Claim 1 covered all three object types the user named. | No (caught on review) | First per-audience rebuild used `QuantumState["0"]` everywhere: no symbolic `QuantumState`. Fixed: Student is now a symbolic state (Bloch), Industry combines a symbolic Hamiltonian with a symbolic initial state. |
| A10 | 1 | The audience axis (student/educator/...) is the right organizing structure for every claim. | No | For Claim 1 the natural axis is the **object type**, because the claim *is* "this property pervades the object model." Audience labels were noise there. Claim 1 re-cut as an object-type ladder (state -> operator -> evolution -> circuit -> composition), labels dropped. Lesson in section 6. |
| A11 | 1 | `QuantumEvolve` was wrong to use at all (per A1). | Refined | A1's error was making `QuantumEvolve` *the headline* of a pervasive property. Using it as **one rung** (operator -> its symbolic dynamics) is correct and valuable: it shows $U(t)=e^{-iHt}$ in closed form and ties the evolution frequency $\sqrt{J^2+h^2}$ back to the spectrum. Reintroduced in that role. |
| A12 | 2-5 | The Claim 1 lessons transfer cleanly to the other claims by just dropping audience labels. | Mostly, plus surgery | Claims 3 and 5 needed only re-titling to their natural axis (WL operation / view type). Claim 4 needed *consolidation* (3 near-identical QASM rungs collapsed to 2 + an honest bridge rung) and a sharpened Claim-4-vs-5 boundary. Claim 2 needed real surgery (next row). |
| A13 | 2 | Claim 2's examples proved "full theory." | No (cross-contamination) | Three of four leaned on *symbolic* results ($1-2p+2p^2$, $\sqrt{1-\gamma}/2$, $\sqrt{1-p/2}$), which is Claim 1's punch, not Claim 2's. Re-cut on **object category** with **numeric** punches: channel purity $1\to0.875$, SIC-POVM 4 outcomes, $g^{(2)}=0/1/2$, Wigner $W(0,0)=-0.318$. Each foregrounds an object a gate-only framework lacks. |
| A14 | 2 | `TracePreservingQ` is safe to showcase for a CPTP channel. | No | Returns `False` for a genuinely trace-preserving amplitude-damping channel at numeric $\gamma$ (QF numerical-tolerance quirk); `KrausOperators`/`Operators` are undefined. Dropped both; used the purity drop instead. Logged so we do not reach for them again. |
| A15 | 2-5 | Strict claim orthogonality (A13: keep symbolic out of Claim 2-5) is the final word. | Overridden by user | User asked to add a symbolic example to every other claim. Resolution that avoids re-creating A13's problem: each added rung keeps *its own* claim's identity and adds symbolic as a **secondary** dimension, explicitly labelled "and this is symbolic too". Added: Claim 2 symbolic channel ($\rho(\gamma)$), Claim 3 symbolic dispatch (Werner purity/entropy in $p$), Claim 5 symbolic views (parametric diagram + matrix). **Claim 4 cannot**: QASM/hardware need numeric angles (`["QASM"]` errors on a symbolic angle), stated as the honest boundary. |
| A16 | 1 | The manual braket $\langle\psi\|O\|\psi\rangle$ was the only way to get an expectation. | No (user corrected) | A `QuantumCircuitOperator` accepts **states as elements**: `{psi, op -> 1, psi["Dagger"]}` compiles to the scalar $\langle\psi\|O\|\psi\rangle$. More idiomatic and more compelling for Claim 1 (the *object itself* is the braket, all symbols). Adopted in the final rung. Lesson: ask whether the framework already expresses the thing as an object before hand-rolling it. |
| A17 | 1 | A ladder of *small* symbolic examples is enough to prove "symbolic-first." | No (too trivial) | Every rung was single-qubit / $2\times2$: fine as principle, unconvincing to a serious reader. Added a **capstone**: a 2-qubit Heisenberg quench ($H=J(XX{+}YY{+}ZZ)$), symbolic evolution $\|\psi(t)\rangle = e^{iJt}[\cos 2Jt\,\|01\rangle - i\sin 2Jt\,\|10\rangle]$, and the entanglement entropy $S(t) = -\cos^2 2Jt\,\log_2\cos^2 2Jt - \sin^2 2Jt\,\log_2\sin^2 2Jt$ as a closed-form function of time. **Lesson: a uniqueness claim needs at least one example that is hard or impossible on the competitor**, not only clean small ones. The contrast-test win here is concrete: numeric tools give $S$ at one $t$; QF gives the whole $S(t)$. Verified `matches_QF -> True`; new image `claim1_entropy.png`. |

---

## 3. Claim 1 case study: critique and redesign

### 3.1 The questions raised (and the honest answers)

**Who is the target audience for the `QuantumEvolve` example?**
At best a student or educator who recognizes a Rabi oscillation in $\cos t\,\lvert0\rangle - i\sin t\,\lvert1\rangle$. It does **not** hook researchers or industry: the result is textbook, hand-derivable, and contains no parameter in the *model* (only time appears, and only because the solver introduces it). By gate 4 of the pipeline it serves at most two audiences, weakly. The critique is correct: it is the most trivial example in the set.

**Why `QuantumEvolve` and not a symbolic Hamiltonian (`QuantumOperator`) and symbolic initial state (`QuantumState`)?**
This was assumption A1. I reached for the function whose uniqueness is most *visible* (`DSolve` under the hood) and let it represent the claim. That narrows a pervasive property to one function. The deeper, truer demonstration is that QF objects are symbolic *everywhere*: you can build $H(J,h) = J\,Z + h\,X$ as a `QuantumOperator` with free parameters, prepare a symbolic `QuantumState`, and read off exact properties, all without evolving anything.

**Why not a parametric / symbolic `QuantumCircuitOperator`?**
This is the strongest missing piece, and the one that actually reaches researchers and industry. A parametric circuit (the shape of every VQE / QAOA ansatz) compiles in QF to a **closed-form unitary in the angles**, and its observable expectation is a **closed-form cost landscape**. In Qiskit you bind parameters and run a statevector simulation per point; in QF you get the function. That is a genuine contrast-test win, and the first draft omitted it entirely.

### 3.2 The redesign: a symbolic ladder (all rungs kernel-verified)

Instead of one trivial evolution, ladder from recognizable to consequential. Every output below was captured from a live kernel (QF 2.0.0, WL 15.0).

**Rung 1, symbolic operator** (students/educators recognize the gate):
```wolfram
QuantumOperator["RX"[\[Theta]]]["Matrix"] // Normal
```
$$
\begin{pmatrix}\cos\frac{\theta}{2} & -i\sin\frac{\theta}{2}\\[1mm] -i\sin\frac{\theta}{2} & \cos\frac{\theta}{2}\end{pmatrix}
$$

**Rung 2, symbolic Hamiltonian with free parameters** (the model itself is symbolic):
```wolfram
H = J QuantumOperator["Z"] + h QuantumOperator["X"];
H["Matrix"] // Normal
```
$$
\begin{pmatrix}J & h\\ h & -J\end{pmatrix}
$$

**Rung 3, symbolic state property** (one symbolic state, exact geometry):
```wolfram
QuantumState[{Cos[\[Alpha]], Sin[\[Alpha]] Exp[I \[Beta]]}]["BlochVector"]
```
$$
\bigl(\sin 2\alpha\,\cos\beta,\ \ \sin 2\alpha\,\sin\beta,\ \ \cos 2\alpha\bigr)
$$

**Rung 4, parametric circuit to closed-form cost** (the researcher/industry hook):
```wolfram
ansatz = QuantumCircuitOperator[{"RY"[\[Theta]] -> 1, "RZ"[\[Phi]] -> 1}];
psi = ansatz[QuantumState["0"]];
v = Normal[psi["StateVector"]];
FullSimplify[Conjugate[v] . Normal[QuantumOperator["Z"]["Matrix"]] . v,
  Assumptions -> {\[Theta] \[Element] Reals, \[Phi] \[Element] Reals}]
```
$$
\langle Z \rangle(\theta,\phi) = \cos\theta
$$
The full variational cost landscape, in closed form, with zero sampling. (Optionally, rung 5 is `QuantumEvolve` for symbolic dynamics, now demoted to "and time evolution is symbolic too," not the headline.)

### 3.3 Why the ladder serves all four audiences

| Audience | What hooks them | Rung |
|---|---|---|
| Student | Recognizes $R_X(\theta)$ and the Bloch parametrization as exact algebra, not floats | 1, 3 |
| Educator | Can derive each rung on the board and reproduce it identically | 1, 2, 3 |
| Researcher | Symbolic Hamiltonian parameters and a closed-form VQE cost landscape | 2, 4 |
| Industry | Parametric ansatz to analytic cost: optimize the function, do not sample it | 4 |

### 3.4 Verification record
- Probes: `/tmp/qfshow/c1alt.wls`, `c1circ.wls`, `c1exp.wls`, `c1final.wls`.
- `QuantumState["ExpectationValue", ...]` and `QuantumOperator["Expectation", ...]` are **not** valid; expectation was computed as $\langle\psi|Z|\psi\rangle$ directly. (Open item: find the QF-native expectation idiom, if one exists, rather than the manual braket.)
- All symbolic simplifications used explicit real-variable assumptions; without them QF correctly keeps `Conjugate[...]` terms, which is a feature, not a bug.

---

## 4. Status / open decisions

- [x] Apply the Claim 1 redesign to [QF-Showcase.md](QF-Showcase.md). **Done**, restructured as one example per audience.
- [x] Re-audit Claims 2 to 5 and rebuild as one verified example per audience (20 total). **Done.** All but one line (Claim 4 industry, hardware submission) executed in-kernel.
- [x] Switch structure to one example per audience. **Done** (supersedes A8).
- [x] Find a QF-native expectation-value idiom. **Found (user-supplied):** a `QuantumCircuitOperator` accepts **states as elements**. `QuantumCircuitOperator[{psi, op -> 1, psi["Dagger"]}]` compiles to a 0-qudit operator whose scalar value is $\langle\psi|op|\psi\rangle$ (ket first, bra `["Dagger"]` last). Verified: gives $a\cos\alpha+b\sin\alpha$ for the VQE Hamiltonian, and $\cos\frac\theta2\cos\frac\phi2 - i\cos(\alpha+\frac\theta2)\sin\frac\phi2$ for a unitary. Replaces the manual braket in Claim 1's last rung. (`["ExpectationValue"]`/`["Expectation"]` are still not valid; this is the idiomatic route.)
- [ ] Optional: install Python + qiskit so the Claim 4 industry hardware line can be verified (or run it on a machine that has a token).
- [ ] Optional: the Claim 5 multiway graph (GHZ(2)) is sparse; a richer circuit would make the "branchial view" more striking.

## 5. Verified example map (current showcase)

**All five claims are now cut by their own natural axis (audience labels dropped everywhere, A10 + A12).**

| Claim | Axis | Rungs (each kernel-verified unless noted) |
|---|---|---|
| 1 Symbolic-first | **object type**, building to a hard capstone | symbolic `QuantumState` (Bloch) -> symbolic `QuantumOperator` $H$ (spectrum $\pm\sqrt{J^2+h^2}$) -> `QuantumEvolve` ($U(t)$, freq $=\sqrt{J^2+h^2}$ = gap) -> parametric `QuantumCircuitOperator` ($\langle Z\rangle=\cos\theta$) -> bra-ket object ($\langle\psi\|H\|\psi\rangle$, VQE min $-\sqrt{a^2+b^2}$) -> **capstone: 2-qubit Heisenberg quench, closed-form $S(t)$** (impossible on numeric platforms) |
| 2 Full theory | **object category** (numeric punches) | open-system channel (purity $1\to0.875$) -> SIC-POVM (4 outcomes, dim 2) -> optics $g^{(2)}$ ($0/1/2$) -> Wigner negativity ($W(0,0)=-0.318$) |
| 3 First-class WL | **WL operation** | property dispatch -> into `Plot` -> into `Reduce` (Werner threshold $p<\tfrac12$) -> into `Dataset` |
| 4 Interop | **export target** | QASM export -> QASM faithful for QFT -> external bridges (Qiskit/pyzx/QuEST/IBM/AWS, **not executed**) |
| 5 Many views | **view type** | Bloch sphere -> diagram + matrix -> tensor network + multiway -> measurement histogram |

Boundary kept explicit: **Claim 4 = leaves QF for another tool**; **Claim 5 = stays in QF, rendered differently.**

## 6. Process lessons (the reusable essence)

These generalize beyond this task.

1. **Pick the organizing axis from the claim's mechanism, not a fixed template.** A claim about a property that *pervades the object model* (Claim 1, "symbolic-first") is best cut by **object type**; a claim about *who benefits* may be cut by **audience**; a claim about *coverage* may be cut by **physical regime**. Forcing an audience grid onto Claim 1 added noise (A10). Decide the axis first, then fill it.

2. **Keep claims on orthogonal axes, and watch for cross-contamination.** Two claims should differentiate along *different* dimensions (see the Claim 1 vs Claim 2 note below). If a Claim-2 example's wow-factor is actually its symbolic closed form, it is silently doing Claim 1's job. Make each example carry *its own* claim's punch.

3. **A rejected tool can return in a smaller role.** `QuantumEvolve` was wrong as the *headline* of "symbolic-first" (A1) but correct as *one rung* showing operator -> dynamics (A11). The fix to a bad example is often re-scoping, not deletion.

4. **Tie rungs together for a through-line.** The evolution frequency $\sqrt{J^2+h^2}$ equals the spectrum gap from the previous rung; the VQE minimum equals that same eigenvalue. Internal echoes make a showcase read as one argument, not a list.

5. **Self-reported "done" is unreliable; re-read the artifact against the literal ask.** The symbolic-state gap (A9) survived a "done" claim and was only caught on a forced re-read. Verify the deliverable, not your memory of it.

### Claim 1 vs Claim 2: the axes they differentiate on
- **Claim 1 ("symbolic-first") differentiates on the *exactness / computation-mode* axis:** symbolic closed forms vs floating-point. The contrast-test win is *the form of the answer* ($\cos(t\sqrt{J^2+h^2})$, not `0.7071...`). It can apply to *any* object type.
- **Claim 2 ("full theory") differentiates on the *scope / ontology* axis:** the existence of object categories (channels, POVMs, open-system/Lindblad, second-quantized optics, tomography) that a gate-only framework does not have at all. The contrast-test win is *that the object exists*, independent of whether its output is symbolic.
- **Risk:** they overlap when a Claim-2 example leans on a symbolic result (e.g. the amplitude-damping coherence $\sqrt{1-\gamma}/2$ is *also* a Claim-1 flex). To keep them distinct, Claim 2 examples should foreground the *object category* and ideally include a **numeric** one (the $g^{(2)}=0/1/2$ photon statistics is good precisely because it is numeric, isolating "scope" from "symbolic").
