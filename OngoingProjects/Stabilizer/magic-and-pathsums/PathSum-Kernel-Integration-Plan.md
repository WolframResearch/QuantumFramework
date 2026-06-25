# Path-sum to QF kernel: integration plan (minimal master-function design)

**Status.** Design only. No `Kernel/*.m` file is modified and nothing is committed. Every kernel citation below was
re-verified against the live working-tree paclet at HEAD `d0e71690` (a fresh `wl-verifier` reproduced all behavioral
claims and `file:line` cites, OPEN ISSUES: 0). The packaging mechanism was audited from `PacletInfo.wl`,
`Kernel/QuantumFrameworkMain.m`, and `Kernel/QuantumFramework.m`. Behavioral seams were confirmed by an out-of-tree
scratch script (`PacletDirectoryLoad`, no `Quiet`, no `Print`).

**Design principle (the spine of this plan).** Add the path-sum capabilities with **zero new public symbols**. The
existing `StabilizerFrame` head already *is* QF's object for a state past the Clifford boundary, and its dense
readout already branches on an internal backend key. So every capability folds onto that one head as
**backend-aware property dispatch**: `StabilizerFrame` becomes the single master magic-state object carrying two
internal encodings (a component list, and a compressed phase polynomial) behind a unified property API. The result
is a handful of new property strings on one head, no `PacletInfo` change, and no new documentation page.

Sources read in full first: the path-sum build (`pathsum.wl`, `elimination.wl`, `diagonal.wl`, `rankcap.wl`,
`harness.wl`) and its audit; the sibling notes
[`StabilizerFrame-NonClifford-Blowup.md`](StabilizerFrame-NonClifford-Blowup.md),
[`StabilizerFrame-Observables-Design.md`](StabilizerFrame-Observables-Design.md),
[`SymPhase-vs-StabilizerRank-Algebraic-Branching.md`](SymPhase-vs-StabilizerRank-Algebraic-Branching.md); the
symbolic-angle builder [`pathsum-blowup-alternatives-demos.wls`](pathsum-blowup-alternatives-demos.wls); and the QF
audit active set.

---

## 1. New surface added (complete list)

**No new public symbols. No `PacletInfo.wl` change. No new doc page.** All additions are property strings on the
existing `StabilizerFrame` head, one internal backend, and edits to existing rules.

### New master methods on `StabilizerFrame` (new property strings, backend-aware)
- `f_StabilizerFrame["Expectation", P_String]` : exact $\langle\psi\lvert P\rvert\psi\rangle$, no densification.
- `f_StabilizerFrame["Amplitude", y_List]` : the computational amplitude $\langle y\lvert\psi\rangle$.
- `f_StabilizerFrame["Compress"]` : the master post-hoc compressor (span-dimension rank cap).

### Extended existing methods (same names, new backend / polynomial path)
- `f_StabilizerFrame["StateVector"]` / `["State"]` : add the phase-polynomial branch.
- `stabilizerInnerProduct[fA_StabilizerFrame, fB_StabilizerFrame]` : shared-reference polynomial path and
  phase-polynomial backend; the existing densification stays as the fallback.

### New internal backend + widened predicate (no public symbol)
- `StabilizerFrame[<|"PhasePolynomial" -> \[CurlyPhi], "Variables" -> {...}, "LinearMap" -> M,
  "InverseMap" -> Minv, "Scale" -> 2^{-n/2}, "Reference" -> ps0, "Qubits" -> n|>]` : the compressed backend.
- `StabilizerFrameQ` widened to accept it (a second definition; existing component frames are untouched).

### New build-time entry (no public symbol)
- `Method -> {"Stabilizer", "Compress" -> "PhasePolynomial"}` on the circuit-apply path: builds the
  phase-polynomial backend directly from a diagonal circuit (never materializing the $2^t$ frame); falls back to the
  ordinary stabilizer path if an interior Hadamard or non-diagonal gate appears.

### Edited rules (behavior change, no new symbol)
- `ps_PauliStabilizer["P"[phase_], j_]` and `ps["T", j]` : full-angle convention ($T$ unchanged).
- `f_StabilizerFrame["P"[phase_], q_]` and `f["T", q]` : the same on the frame side.
- `_StabilizerFrame["Properties"]` : append `"Expectation"`, `"Amplitude"`, `"Compress"`.

### Internal `PackageScope` helpers (plumbing, not user-facing functions)
- Symplectic-Pauli kernel (one reusable module): `framePauliMultiply`, `framePauliDagger`, `framePauliString`,
  `frameSandwich`, `sfGramMatrix`. Serves `["Expectation"]`, `["Amplitude"]`, `["InnerProduct"]`, and `["Compress"]`.
- Phase-polynomial carrier: `sfPhasePolyFold` (diagonal gate fold), `sfPhasePolyFromCircuit` (builder).
- Compression: `sfCompress`.

### Symbolic angles (no new function)
- After the convention fix, a free $\theta$ flows through `["Expectation"]` / `["Amplitude"]` as exact symbolic
  arithmetic; the gradient is `D[f["Expectation", P], \[Theta]]`.

---

## 2. QF package-style compliance (audited at HEAD `d0e71690`)

The conventions every addition follows, verified from the loader and the existing files:

1. **Auto-scan loader, no manifest.** `QuantumFrameworkMain.m:23` calls
   `PacletManager`Package`loadWolframLanguageCode[..., "AutoloadSymbols" -> {}]` (the WL PackageFramework): it
   discovers and loads every `.m` under `Kernel/` whose first line is `Package["Wolfram`QuantumFramework`"]`. There
   is **no** explicit `Get`/`Needs` list of Stabilizer files. A new file is loaded just by existing with that first
   line; this plan adds at most an optional `PackageScope`-only file and needs no loader edit.
2. **A public symbol needs two registrations**, confirmed on `GraphState`/`LocalComplement`: `PackageExport[Sym]` at
   the top of its file **and** the full name in `PacletInfo.wl` `Symbols` (`:50-55`). **This design adds no public
   symbol, so `PacletInfo.wl` is untouched.** A method on an existing head is a `DownValue`, not a symbol.
3. **No new sub-context.** Every `Kernel/Stabilizer/*.m` file lives in the primary `Wolfram`QuantumFramework`
   context (the `PacletInfo.wl` `Context` list `:17-28` is unchanged); the additions stay in that context.
4. **Cross-file `PackageScope` visibility is load-order-independent** (the framework collects declarations before
   evaluating bodies), so the new handlers may call `stabilizerExpectation`, `sfApplyPauli`, `encodeStabilizerGates`,
   `psConcreteFastQ`, and `QuantumShortcut` from sibling files freely.
5. **Canonical `StabilizerFrame.m` section order** (existing): `Package[...]` line, exports, predicate,
   `_StabilizerFrame["Properties"]` contract, direct lookup, constructors, relating-Pauli helpers, property/method
   handlers, gate updates, `Plus`/`Times`/`Equal` upvalues, `MakeBoxes`. New handlers slot into the matching
   sections.
6. **No `Usage.m` edit** (it is auto-generated from doc pages); the `StabilizerFrame` md2nb doc page gains the new
   properties.
7. **House style** (CLAUDE.md / `wl-style-guide`): `Block`/`With` over `Module`; `Enclose`/`Confirm` for new error
   handling; `SparseArray` and `Modulus -> 2` linear algebra; no `Quiet`, no `Print`.

**Net:** edits to existing files only (`StabilizerFrame.m`, `InnerProduct.m`, `GateUpdates.m`, `PauliStabilizer.m`),
plus an optional self-contained `PackageScope`-only `Kernel/Stabilizer/PhasePolynomial.m` for the carrier fold if
`StabilizerFrame.m` is kept lean. No new public surface either way.

---

## 3. Additions at a glance

Legend: **DV** `DownValue`, **PS** `PackageScope`. Every row is on the existing `StabilizerFrame` head or an existing
file; there is no `PackageExport` row and no `PacletInfo` row by design.

| # | Stage | Surface | Exact signature | Kind | File : seam | Backend behaviour / mechanic | Reuses |
|---|---|---|---|---|---|---|---|
| 1 | 1 | Master expectation | `f_StabilizerFrame["Expectation", P_String]` | DV | `StabilizerFrame.m` (beside `:55`) | component: `frameSandwich`; phase-poly: closed form over the support | `stabilizerExpectation` `InnerProduct.m:156` |
| 2 | 1 | Master amplitude | `f_StabilizerFrame["Amplitude", y_List]` | DV | `StabilizerFrame.m` | component: relating-Pauli sum; phase-poly: $2^{-n/2}\omega^{\varphi(M^{-1}y)}$ | `sfApplyPauli` `:147`; `LinearSolve[_,_,Modulus->2]` |
| 3 | 2 | Master compress | `f_StabilizerFrame["Compress"]` | DV | `StabilizerFrame.m` | component: span-dimension rank cap; phase-poly: identity | `sfGramMatrix` (row 9) |
| 4 | 1 | Extend materialize | edit `f_StabilizerFrame["StateVector"]` | DV (edit) | `StabilizerFrame.m:174` | add the phase-poly branch (table of row-2 amplitudes) | row 2 |
| 5 | 1 | Extend inner product | edit `stabilizerInnerProduct[fA_StabilizerFrame, fB_StabilizerFrame]` | DV (edit) | `InnerProduct.m:111` | shared-ref polynomial via `frameSandwich`; phase-poly supported; densify fallback | `stabilizerExpectation` |
| 6 | 2 | Widen predicate | second def of `StabilizerFrameQ` for the `"PhasePolynomial"` key | DV | `StabilizerFrame.m:39` | accept the compressed backend; existing frames unaffected | - |
| 7 | 2 | Build-time entry | edit `PauliStabilizerApply`: `Method -> {"Stabilizer", "Compress" -> "PhasePolynomial"}` | DV (edit) | `PauliStabilizer.m:176` | diagonal-circuit guard (no interior H) builds the phase-poly frame; else fallback | `QuantumShortcut`, `sfPhasePolyFromCircuit` |
| 8 | 1/2 | Contract | append `"Expectation","Amplitude","Compress"` to `_StabilizerFrame["Properties"]` | edit | `StabilizerFrame.m:47` | catalog | - |
| 9 | 1 | Symplectic-Pauli kernel | `framePauliMultiply`, `framePauliDagger`, `framePauliString`, `frameSandwich`, `sfGramMatrix` | PS | `StabilizerFrame.m` | the one reusable algebra for rows 1, 2, 3, 5 | `sfApplyPauli` encoding `:147`; `stabilizerExpectation` |
| 10 | 2 | Phase-poly carrier | `sfPhasePolyFold`, `sfPhasePolyFromCircuit` | PS | `StabilizerFrame.m` (or optional PS-only `PhasePolynomial.m`) | diagonal gate fold + compression | ported diagonal subset of `pathsum.wl`/`diagonal.wl`; `Inverse[_,Modulus->2]` |
| 11 | 3 | Convention fix (state) | edit `ps["P"[phase_], j_]` $\to$ `c = Exp[I phase]`; `ps["T",j] := ps["P"[Pi/4], j]` | DV (edit) | `GateUpdates.m:182, :196` | full-angle; $T$ unchanged | - |
| 12 | 3 | Convention fix (frame) | edit `f["P"[phase_], q_]` $\to$ `c = Exp[I phase]`; `f["T",q] := f["P"[Pi/4], q]` | DV (edit) | `StabilizerFrame.m:231, :248` | full-angle, frame side | - |
| 13 | 3 | Symbolic angle + gradient | (no new function) free $\theta$ through rows 1/2; `D[f["Expectation", P], \[Theta]]` | - | rows 1, 2 | exact symbolic arithmetic | - |

---

## 4. How each capability surfaces (the design decisions)

### 4.1 The diagonal-fragment amplitude and the phase-polynomial object (path-sum capabilities 1 + 2)

**The physics.** A circuit over $\{T, T^\dagger, S, S^\dagger, Z, CZ, CCZ, CNOT\}$ (no Hadamard) acting on
$|+\rangle^{\otimes n}$ keeps each qubit value $E_q(x)$ a linear $\mathbb{F}_2$ form, so $E:\mathbb{F}_2^n \to
\mathbb{F}_2^n$ is a linear bijection (CNOT generates $\mathrm{GL}(n,\mathbb{F}_2)$) and the phase composes into a
single degree-$\le 3$ polynomial $\varphi:\mathbb{F}_2^n \to \mathbb{Z}_8$. The amplitude is one closed-form term,
$$
\langle y \,|\, U \,|+\rangle^{\otimes n} \;=\; 2^{-n/2}\,\omega^{\varphi(E^{-1}(y))}, \qquad \omega = e^{i\pi/4},
$$
rank one as a phase-polynomial object: the $\mathbb{Z}_8$ analogue of SymPhase, the genuine speed win measured in the
path-sum audit (flat in $n$, milliseconds where dense grows to seconds).

**Decision: a compressed backend of `StabilizerFrame`, not a new head.** `StabilizerFrame` is already "the object for
a state past the Clifford boundary," and its `["StateVector"]` handler already branches on an internal backend key
(`"Paulis"`, `StabilizerFrame.m:174`). A phase polynomial is a third, compressed encoding of such a state, so it
folds in as a `"PhasePolynomial"` backend with the master reads (`["StateVector"]`, `["Expectation"]`,
`["Amplitude"]`, `["InnerProduct"]`, `["Compress"]`) dispatching on which key is present. This adds **zero public
symbols**. The honest cost: backend-specific accessors (`"Coefficients"`/`"Stabilizers"` for the component backend;
`"PhasePolynomial"`/`"LinearMap"` for the compressed one) apply to only their own backend.

**Build-time entry.** The phase polynomial must never materialize the $2^t$ frame, so it is built straight from the
circuit at apply time: `Method -> {"Stabilizer", "Compress" -> "PhasePolynomial"}` on the existing apply path. The
guard scans the gate heads (`Replace[specs, (g_ -> _) :> g, {1}]`) and fires only when no head is `"H"` and every
head is in the diagonal+CNOT set (verified live: a diagonal list gives `False`, an `{"H","T","H"}` list gives
`True`); otherwise it returns `$Failed` and the apply path proceeds on the ordinary trichotomy. The default
`Method -> "Stabilizer"` is unchanged.

### 4.2 The rank cap (path-sum capability 3)

**Decision: one master method `f["Compress"]`.** A gate-built frame keeps $2^t$ components with no merging (blow-up
report §3); `["Compress"]` keeps a maximal linearly-independent subset and re-expresses the state exactly. The
recommended in-kernel form is polynomial in $n$: because all components of a gate-built frame share one tableau and
differ only by a relating Pauli, the Gram matrix $G_{ij} = \langle u_i|u_j\rangle = \langle\mathrm{ref}|P_i^\dagger
P_j|\mathrm{ref}\rangle$ is one `stabilizerExpectation` per entry (`InnerProduct.m:156`), exact and phase-correct,
$\chi^2\,\mathrm{poly}(n)$, no densification. The span dimension is the exact rank of $G$ and the reduced
coefficients are an exact projection. This removes the prototype's numeric-`MatrixRank` float caveat. A densifying
fallback (port of `rankcap.wl` via `sfApplyPauli`, $O(\chi\,2^n)$, small $n$) is available where no Gram is wanted.

### 4.3 The general path-sum engine (path-sum capability 4)

**Decision: keep out of the kernel.** Out-scaled by the existing dense `Method -> "Schrodinger"` past $k \approx 14$
$T$ gates, with a brute-force residual that is the opposite of a kernel fast path. Reasons in the do-not-add list.

### 4.4 The symbolic-angle capability (path-sum capability 5)

**Decision: no new function; it rides `["Expectation"]`/`["Amplitude"]` with a symbolic angle, plus the convention
fix.** The frame already carries symbolic coefficients $(1 \pm e^{i\theta/2})/2$ for a free `"P"[\theta]`. Once
`["Expectation"]` computes $\langle P\rangle$ by exact symbolic arithmetic, feeding a free $\theta$ returns a
closed-form $\langle P\rangle(\theta)$, and the gradient is `D[..., \theta]`. The prerequisite is the half-angle
fix: verified live, the frame's `"P"[\theta]` realizes $\mathrm{diag}(1, e^{i\theta/2})$ (fidelity $1.0$ against the
dense half-angle gate, $0.97$ against the full-angle gate), while `QuantumOperator["P"[\theta]]` and the path-sum
builder use the full angle. The fix is one token per rule with $T$ preserved (§5.5).

---

## 5. How each implementation works, fully

### 5.1 The unified `StabilizerFrame` backends

```
(* existing component backend *)
StabilizerFrame[<|"Components" -> {{c_i, ps_i}, ...}, "Paulis" -> {{x_i, z_i, coeff_i}, ...}|>]
(* new phase-polynomial backend *)
StabilizerFrame[<|"PhasePolynomial" -> phi, "Variables" -> {x[1],...,x[n]}, "Qubits" -> n,
                  "LinearMap" -> M, "InverseMap" -> Minv, "Scale" -> 1/Sqrt[2]^n, "Reference" -> ps0|>]
```
`StabilizerFrameQ` gets a second definition accepting the `"PhasePolynomial"` key; the existing component-frame
definition (`StabilizerFrame.m:39`, requiring `"Components" -> {{_, _PauliStabilizer}..}`) is unchanged, so no
existing frame is affected and the compressed backend is produced only by the explicit build-time / `["Compress"]`
paths. The master reads dispatch on which key is present, exactly as `["StateVector"]` already branches on
`"Paulis"` (`StabilizerFrame.m:174`).

### 5.2 The symplectic-Pauli kernel (one module; serves rows 1, 2, 3, 5)

Relating Paulis are encoded `{x, z, s}` for the operator $s\,X^{x}Z^{z}$ (a relating Pauli `{x, z, coeff}` enters as
$s = \mathrm{coeff}\cdot i^{x\cdot z}$).

```
framePauliMultiply[{x,z,s},{x',z',s'}] := {BitXor[x,x'], BitXor[z,z'], s s' (-1)^(z.x')}
framePauliDagger[{x,z,s}]              := {x, z, Conjugate[s] (-1)^(x.z)}
framePauliString[{x,z,s}]              := (* {x,z,s} -> signed "IXYZ" string for stabilizerExpectation *)
frameSandwich[fA_, pg:{_,_,_}, fB_]    := Module[{a, pa, b, pb, ref},
  a = fA["Coefficients"]; pa = First[fA]["Paulis"]; b = fB["Coefficients"]; pb = First[fB]["Paulis"];
  ref = fA["Components"][[1, 2]];                          (* shared reference, all components share its tableau *)
  Sum[Conjugate[a[[i]]] b[[j]] *
      With[{q = framePauliMultiply[framePauliDagger[pa[[i]]], framePauliMultiply[pg, pb[[j]]]]},
           q[[3]] stabilizerExpectation[ref, framePauliString[q]]],     (* underbrace in {0,+1,-1} *)
      {i, fA["Length"]}, {j, fB["Length"]}]]
```
The mathematics: $\langle\psi|P|\psi\rangle = \sum_{ij}\overline{c_i}c_j\langle\mathrm{ref}|P_i^\dagger P
P_j|\mathrm{ref}\rangle$, every factor a Pauli, so each sandwich is a scalar phase times a single
`stabilizerExpectation` on the reference. No densification, $\chi^2\,\mathrm{poly}(n)$.

### 5.3 Rows 1, 2, 4: `["Expectation"]`, `["Amplitude"]`, `["StateVector"]`

```
f_StabilizerFrame["Expectation", P_String] := If[KeyExistsQ[First[f], "PhasePolynomial"],
  phasePolyExpectation[f, P],                              (* closed form for diagonal P; else via ["StateVector"] *)
  frameSandwich[f, framePauliFromString[P], f]]

f_StabilizerFrame["Amplitude", y_List] := With[{a = First[f]},
  If[KeyExistsQ[a, "PhasePolynomial"],
     a["Scale"] Exp[I Pi/4]^Mod[a["PhasePolynomial"] /. Thread[a["Variables"] ->
        LinearSolve[a["LinearMap"], y, Modulus -> 2]], 8],     (* symbolic phi: a["Scale"] Exp[I (phi /. ...)] *)
     (* component backend: relating-Pauli sum over reference amplitudes (observables-design section 6) *)
     With[{ref = a["Components"][[1, 2]]},
       Sum[f["Coefficients"][[i]] referenceAmplitudeThroughPauli[a["Paulis"][[i]], y, ref], {i, f["Length"]}]]]]

f_StabilizerFrame["StateVector"] := With[{a = First[f]},
  Which[
    KeyExistsQ[a, "PhasePolynomial"],
       SparseArray[Table[f["Amplitude", IntegerDigits[k, 2, a["Qubits"]]], {k, 0, 2^a["Qubits"] - 1}]],
    KeyExistsQ[a, "Paulis"],  (* existing coherent branch, unchanged *) ...,
    True,                     (* existing independent-sum branch, unchanged *) ...]]
```
Born probabilities are the documented one-liner $\Pr(\pm) = (\langle\psi|\psi\rangle \pm \langle P\rangle)/2
\langle\psi|\psi\rangle$ with $\langle\psi|\psi\rangle =$ `f["Expectation", "I..I"]`, not a separate method. The
phase-poly amplitude is `diagonal.wl`'s `diagAmp`; the $\mathbb{F}_2$ `LinearSolve`/`Inverse` are the primitives at
`Constructors.m:55` and `InnerProduct.m:69, :178` (verified live: a CNOT-type map is invertible over $\mathbb{F}_2$
and `LinearSolve[M, y, Modulus -> 2]` returns the unique preimage).

### 5.4 Row 3: `["Compress"]` (master post-hoc compressor)

```
sfGramMatrix[f_] := With[{p = First[f]["Paulis"], ref = f["Components"][[1, 2]]},
  Table[With[{q = framePauliMultiply[framePauliDagger[p[[i]]], p[[j]]]},
             q[[3]] stabilizerExpectation[ref, framePauliString[q]]], {i, f["Length"]}, {j, f["Length"]}]]

f_StabilizerFrame["Compress"] := If[KeyExistsQ[First[f], "PhasePolynomial"], f,   (* already minimal *)
  Module[{G = sfGramMatrix[f], c = f["Coefficients"], basis, d},
    basis = (* greedy maximal independent rows of G, EXACT rank, reference index 1 first *);
    d = LinearSolve[G[[basis, basis]], G[[basis, All]] . c];                       (* exact projection *)
    StabilizerFrame[<|"Components" -> Thread[{d, f["Stabilizers"][[basis]]}],
                      "Paulis" -> First[f]["Paulis"][[basis]]|>]]]
```
Span-dimension cap, $\chi^2\,\mathrm{poly}(n)$, exact (no numeric `MatrixRank`). Requires `"Paulis"`; returns `f`
unchanged otherwise.

### 5.5 Row 7: build-time phase-poly entry (no symbol)

```
sfPhasePolyFold[ps_, gate_] := (* fold one diagonal gate; e_q := boolToInt[val[q]]:
   "T":phase+=e_q;  "S":+=2 e_q;  "Z":+=4 e_q;  Tdg:+=7 e_q;
   "CZ":+=4 e_a e_b;  "CCZ":+=4 e_a e_b e_c;  "CNOT": val[k] = F2[val[j] + val[k]] *)

sfPhasePolyFromCircuit[qco_] := Module[{specs = QuantumShortcut[qco], heads, n = qco["Max"], ps, M},
  heads = Replace[specs, (g_ -> _) :> g, {1}];
  If[MemberQ[heads, "H"] || ! SubsetQ[{"T", SuperDagger["T"], "S", "Z", "CZ", "CCZ", "CNOT", "CX"}, heads],
     $Failed,                                                          (* signal: caller falls back *)
     ps = Fold[sfPhasePolyFold, init[n], Join[hLayer[n], specs]];      (* prepend the |+...+> H-layer *)
     M = phasePolyLinearMatrix[ps, n];
     StabilizerFrame[<|"PhasePolynomial" -> z8[ps["phase"]], "Variables" -> Array[x, n], "Qubits" -> n,
        "LinearMap" -> M, "InverseMap" -> Inverse[M, Modulus -> 2], "Scale" -> 1/Sqrt[2]^n,
        "Reference" -> PauliStabilizer[n]|>]]]
```
`PauliStabilizerApply` (`PauliStabilizer.m:176`), handed `"Compress" -> "PhasePolynomial"`, calls this and uses the
result unless it is `$Failed`. The phase-poly backend is read-oriented (StateVector / Expectation / Amplitude /
InnerProduct / Compress); applying a further gate to it is a clean future extension (a diagonal gate updates
$\varphi$/$M$ in place; an $H$ materializes), out of scope now. `QuantumShortcut[qco]` was verified live to emit
`g -> {q}` specs, and `encodeStabilizerGates` to return `$Failed` on a $T$ (the natural fall-through marker).

### 5.6 Rows 11-12 convention fix; row 13 symbolic angles

One token per rule: `c = Exp[I phase/2]` becomes `c = Exp[I phase]`, and `"T"` retargets from `"P"[Pi/2]` to
`"P"[Pi/4]` so $c = e^{i\pi/4}$ for $T$ is preserved (verified: $T\lvert+\rangle = \{1, e^{i\pi/4}\}/\sqrt2$
unchanged). After the fix a frame built with a free `"P"[θ]` carries $\theta$ symbolically through `frameSandwich`
(row 1) and `["Amplitude"]` (row 2), so `f["Expectation", "Z"]` on `{"H", "P"[θ], "H"}` is $\cos\theta$ and
`D[..., θ]` is $-\sin\theta$, with no new function.

---

## 6. Staged plan (lowest-risk / highest-value first)

| Stage | What | One-line risk | One-line test gate |
|---|---|---|---|
| 1 | Master reads on the component backend: `["Expectation"]`, `["Amplitude"]`, polynomial `["InnerProduct"]`, the symplectic-Pauli kernel | thin reuse of `stabilizerExpectation`; distinct-reference inner-product phase stays the `InnerProduct.m:41` TODO (the shared-reference path avoids it) | $\langle P\rangle$ vs dense over 72 random Clifford+$T$ frames, worst error $< 10^{-15}$; stabilizer suite stays green |
| 2 | `["Compress"]` rank cap; the phase-polynomial backend + `sfPhasePolyFromCircuit` + the `"Compress" -> "PhasePolynomial"` apply option + widened `StabilizerFrameQ` | widening the predicate and backend-branching the read methods (precedent: `["StateVector"]` already branches) | `["Compress"]` exact reproduction; phase-poly `["Amplitude"]` `RootReduce`-equal to dense on diagonal circuits; interior-$H$ circuit falls back |
| 3 | Convention fix + symbolic angles + gradient | a back-compat change to frame `"P"[θ]` (a documented trap, so a bugfix; `CHANGELOG` note) | frame `"P"[θ]` $\equiv$ dense full-angle `"P"[θ]` (fail-before/pass-after); $\langle Z\rangle(\theta) = \cos\theta$, $\partial_\theta = -\sin\theta$ |

Stage 1 and Stage 2 are independent in their component-backend parts; Stage 2's phase-poly backend builds on Stage
1's symplectic-Pauli kernel. Stage 3 depends on Stages 1-2.

---

## 7. Do not add (with reasons)

- **A separate public head for the phase polynomial.** Folded into `StabilizerFrame` as a backend instead (the
  user-confirmed minimal-surface decision); a new head would cost a `PackageExport` symbol, a `PacletInfo` line, a
  doc page, and formatting, for no capability gain.
- **`["M"]` / `["Probability"]` as methods.** Derived one-liners from `["Expectation"]` and the norm; not new
  capability.
- **The brute-force Gauss-sum residual (`PathSum`amp`) as any default.** The path-sum audit measures it losing to
  the existing dense `Method -> "Schrodinger"` at $k \approx 14$ $T$ gates and by $5\times$ at $k = 16$; the kernel
  already computes exact amplitudes. Research oracle only.
- **The general Hadamard-elimination rewrite engine (`elimination.wl`) as a `Method`.** Out-scaled by dense and by
  de Colnet's rank-width program; reach is variable-order-dependent (audit caveat 1); $\mathbb{Z}_8$ confluence is
  empirical only. Research tool.
- **The `Inactive[Sum]` formal carrier.** A presentation device with no kernel consumer.
- **The Groebner / variety point-count path.** $\#\mathrm{P}$ territory with no speed win over enumeration on the
  tractable fragments.
- **The GF$(2^k)$ trace route.** Dead at the characteristic-2 wall ($\Phi_8(x) = (1+x)^4$ over $\mathbb{F}_2$;
  $(e^{i\pi/4}X)^2 = iI \ne I$ stabilizes nothing).
- **Auto-switching the default apply path.** The `"Compress" -> "PhasePolynomial"` route is opt-in; the default
  `Method -> "Stabilizer"` trichotomy (`PauliStabilizer.m:176`) is unchanged.
- **NP-hard rank compression** (Bravyi-Gosset sparsification / extent). `["Compress"]` is a span-dimension witness
  only, not the optimal stabilizer rank.
- **Numeric `MatrixRank` in `["Compress"]`.** Use the exact Gram-matrix rank (removes the prototype's float caveat).
- **Gate application onto the phase-poly backend.** Read-oriented for now; in-place gate updates are a clean future
  extension.

A related small fix, noted but not part of this plan: routing `"RZ"[θ]`/`"RY"[θ]` (currently `$Failed` +
`PauliStabilizer::nonclifford`) to the `"P"` frame handler with the reconciled convention would make the diagonal
rotation representable; orthogonal, its own pass.

---

## 8. Files and tests

| File | Edit / new | Adds | `PacletInfo` change |
|---|---|---|---|
| `Kernel/Stabilizer/StabilizerFrame.m` | edit | `["Expectation"]`, `["Amplitude"]`, `["Compress"]`; phase-poly backend + widened `StabilizerFrameQ`; symplectic-Pauli kernel + Gram; `["StateVector"]` phase-poly branch; convention fix; `["Properties"]` append | none |
| `Kernel/Stabilizer/InnerProduct.m` | edit | polynomial shared-reference + phase-poly inner product (`:111`) | none |
| `Kernel/Stabilizer/GateUpdates.m` | edit | convention fix (`:182, :196`) | none |
| `Kernel/Stabilizer/PauliStabilizer.m` | edit | `"Compress" -> "PhasePolynomial"` apply branch (`:176`) | none |
| `Kernel/Stabilizer/PhasePolynomial.m` | optional new (PackageScope only, auto-loaded) | the carrier fold + builder, if `StabilizerFrame.m` is kept lean; **no public symbol** | none |
| `Tests/Stabilizer/Formulas.wlt` | edit | `["Expectation"]`/`["Amplitude"]` vs dense; `["Compress"]` exact reproduction; phase-poly `["Amplitude"]` `RootReduce`-equal to dense on diagonal circuits + H-guard fallback; convention fail-before/pass-after; symbolic $\langle Z\rangle(\theta)$ + gradient | none |
| `StabilizerFrame` doc page (md2nb) | edit | document the new properties + the phase-poly backend | none |

**No new file is required; no `PacletInfo.wl` change; no new public symbol.**

---

## 9. Precise kernel cites (re-verified at HEAD `d0e71690`)

All line numbers are from a direct read of the working-tree files; the only stabilizer commit past the audit anchor
`b4fcdb31` is `b81f3332` (`SymbolicMeasure.m`), which touches none of these.

| Seam | `file:symbol` (line) | Role in this plan |
|---|---|---|
| Frame predicate / contract | `Stabilizer/StabilizerFrame.m:39` `StabilizerFrameQ`, `:47` `_StabilizerFrame["Properties"]` | widened for the phase-poly backend; new properties appended |
| Frame materialization | `Stabilizer/StabilizerFrame.m:174` `f["StateVector"]` (branches on `"Paulis"`) | the precedent for backend-branching; gains the phase-poly branch |
| Relating-Pauli helpers | `Stabilizer/StabilizerFrame.m:147` `sfApplyPauli`, `:134/:137` `sfConjugatePauli`, `:142` `sfTimesZPauli` | encoding reused by the symplectic-Pauli kernel and `["Amplitude"]`/`["Compress"]` |
| Frame inner product | `Stabilizer/InnerProduct.m:111` `stabilizerInnerProduct[fA,fB]` (densifies) | extended with the polynomial shared-reference and phase-poly paths |
| Single-state expectation | `Stabilizer/InnerProduct.m:156` `stabilizerExpectation[ps, P]` | the exact, phase-correct primitive reused by `["Expectation"]`/`["Amplitude"]`/Gram (verified live: GHZ$_3$ $\langle ZZI\rangle,\langle XXX\rangle,\langle ZII\rangle = \{1,1,0\}$) |
| Closed-form overlap (magnitude) | `Stabilizer/InnerProduct.m:50` `stabilizerInnerProductClosedForm`, phase TODO `:41` | the only missing piece for distinct-reference phase; not needed by the shared-reference paths |
| Method dispatch (state) | `Stabilizer/Properties.m:78` `ps["InnerProduct",...]`, `:80` `ps["Expectation", pauli]` | the dispatch shape the new frame methods mirror |
| Non-Clifford boundary | `Stabilizer/GateUpdates.m:182` `ps["P"[phase], j]` (`c = Exp[I phase/2]`), `ps["T", j] :196` | the half-angle origin; the convention fix edits these two lines |
| Frame doubling | `Stabilizer/StabilizerFrame.m:231` `f["P"[phase], q]`, `f["T", q] :248` | the frame-side half-angle; the convention fix edits these |
| Circuit dispatcher | `Stabilizer/PauliStabilizer.m:176` `PauliStabilizerApply` | the opt-in `"Compress"` branch hooks here; the default trichotomy is unchanged |
| Compiled fast path | `Stabilizer/Compiled.m:196` `encodeStabilizerGates`, `:41` `$stabilizerGateCodes`, `:240` `ps["ApplyCircuit", ...]` | `encodeStabilizerGates` returns `$Failed` on a $T$ (verified live), the natural fall-through for the diagonal guard |
| $\mathbb{F}_2$ linear algebra | `Stabilizer/Constructors.m:55` `agExtendToSymplecticBasis`, `InnerProduct.m:69, :178` (`LinearSolve`/`NullSpace`, `Modulus -> 2`) | the primitives the phase-poly amplitude reuses for $E^{-1}$ |
| Canonical vector | `Stabilizer/Conversions.m:69` `ps["State"]` | the reference materialization the readouts call once |
| Packed-path guard | `Stabilizer/Packed.m:87` `psConcreteFastQ` (verified `True` on a register) | the concreteness gate the fast paths respect |
| Loader / registration | `QuantumFrameworkMain.m:23` `loadWolframLanguageCode` (auto-scan), `PacletInfo.wl:50-55` (stabilizer `Symbols`) | confirms a new file auto-loads and that this design needs no `PacletInfo` edit |

---

## 10. Physicist-audience summary

**The concept.** A Clifford+$T$ circuit leaves the stabilizer world, and the cost of following it exactly shows up
either as a number of stabilizer branches or as the size of an interference (Gauss) sum. The path-sum build records
that magic as algebra: a Boolean polynomial for each qubit value and a degree-$\le 3$ phase polynomial valued in the
eighth roots of unity. This plan promotes the two pieces where the algebra genuinely helps and refuses the piece
where it cannot, and it does so with **no new public object**: the existing `StabilizerFrame`, which is already QF's
representation of a magic state, becomes the single master object carrying a second, compressed encoding behind the
same property interface.

**What is added.** For the **diagonal (no-Hadamard) fragment** the state is rank one, an affine support plus one
low-degree phase polynomial, and any amplitude is one closed-form term $2^{-n/2}\omega^{\varphi(E^{-1}(y))}$, flat in
$n$ where the dense vector grows as $2^n$. That becomes a compressed backend of `StabilizerFrame`, built straight
from the circuit so it never pays the $2^t$ blow-up, read through `["Amplitude"]` and `["StateVector"]`. For any
frame, `["Expectation"]` reads $\langle P\rangle$ without building the $2^n$ vector (a sum of single-reference Pauli
expectations, polynomial in the qubit count for fixed magic), and `["Compress"]` collapses the profligate $2^t$
components down to the true span dimension. With a free angle in place of a number, the same exact arithmetic returns
a closed-form $\langle P\rangle(\theta)$ and its analytic gradient, the natural object for variational algorithms,
once a one-token convention fix aligns the frame's internal phase gate with the gate-level full-angle convention
($T$ unchanged). The general engine that handles interleaved Hadamards by a brute-force, $\#\mathrm{P}$-hard Gauss
sum stays out: the build's own benchmarks show it losing to QF's dense path past about fourteen $T$ gates, and
Clifford+$T$ universality forbids any polynomial shortcut.

**The code.** This pass writes no kernel code and commits nothing. The plan is
[PathSum-Kernel-Integration-Plan.md](PathSum-Kernel-Integration-Plan.md). When implemented, every addition is an
edit to an existing `Kernel/Stabilizer/*.m` file (optionally one self-contained `PackageScope`-only file for the
carrier): a few new property strings on `StabilizerFrame`, one internal backend, a reusable symplectic-Pauli helper
module, and a one-token convention fix. No new public symbol, no `PacletInfo` change, no new documentation page.
