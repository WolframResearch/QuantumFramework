# QuantumInterferometer: symbolic behaviour

Every output below was produced by the installed kernel (WL 15.0) on QuantumFramework `main` at a28fd507 (2026-10-27), plus the fixes to `Interferometer.m` committed with this page (lossy probabilities without `Conjugate`, `Assumptions` in `"Probabilities"` and rejection of contradictory ones from `Assuming`, coherent slots in the all-outputs form). That kernel includes the physical-domain contract (c1901374), the compact squeezed, coherent and mixed-state results (23ac0cf1), lossy and coherent sampling (c91ca9ad, 24f5679f) and the photon-number allowlist (46f9ebf6). Association keys print as `QuditName[{...}]`; the `"Dual" -> False` option is omitted here. The usage message of `QuantumInterferometer` now lists only the call forms, so this page is the reference for symbolic behaviour.

**Physical domain.** Photon numbers are nonnegative integers, beam-splitter angles θ, φ and phase angles α are real, loss transmissivities η are real with 0 ≤ η ≤ 1, coherent amplitudes and squeezing parameters are complex. Symbolic results hold on this domain and carry no condition that restates it.

```wl
bs = QuantumInterferometer[{"BS"[th, ph]}];
```

## 1. Symbolic element parameters

Transfer matrix of a beam splitter; an omitted phase is 0.
```wl
bs["Unitary"]
(* {{Cos[th], -Sin[th]/E^(I ph)}, {E^(I ph) Sin[th], Cos[th]}} *)
QuantumInterferometer[{"BS"[th]}]["Unitary"]
(* {{Cos[th], -Sin[th]}, {Sin[th], Cos[th]}} *)
```

Phase shifter.
```wl
QuantumInterferometer[{"PS"[al]}]["Unitary"]
(* {{E^(I al)}} *)
```

Long element names normalize to the short ones; a beam splitter also takes its angles as a list, and its two modes in either order.
```wl
QuantumInterferometer[{"PhaseShift"[al], "Squeeze"[xi]}]["Elements"]
(* {"PS"[al] -> {1}, "S"[xi] -> {1}} *)
QuantumInterferometer[{"BeamSplitter"[{th, ph}]}]["Elements"]
(* {"BS"[th, ph] -> {1, 2}} *)
QuantumInterferometer[{"BS"[th, ph] -> {2, 1}}]["Unitary"]
(* {{Cos[th], E^(I ph) Sin[th]}, {-Sin[th]/E^(I ph), Cos[th]}} *)
```

Elements compose in time order.
```wl
QuantumInterferometer[{"BS"[th, ph], "PS"[al] -> 2}]["Unitary"]
(* {{Cos[th], -Sin[th]/E^(I ph)}, {E^(I al + I ph) Sin[th], E^(I al) Cos[th]}} *)
```

Composition `qi2[qi1]` is qi1 followed by qi2.
```wl
QuantumInterferometer[{"PS"[al] -> 2}, 2][bs]["Unitary"]
(* same matrix as above *)
```

Inverse.
```wl
bs["Dagger"]["Unitary"]
(* {{Cos[th], Sin[th]/E^(I ph)}, {-E^(I ph) Sin[th], Cos[th]}} *)
```

A passive network's transfer matrix is its unitary.
```wl
bs["TransferMatrix"] === bs["Unitary"]
(* True *)
```

Structure does not depend on the parameter values.
```wl
qp = QuantumInterferometer[{"BS"[th, ph], "PS"[al] -> 2}, "Parameters" -> {th, ph, al}];
{qp["Modes"], qp["Depth"], qp["BeamSplitterCount"], qp["PhaseShifterCount"], qp["Parameters"], qp["Layers"], qp["Method"]}
(* {2, 1, 1, 1, {th, ph, al}, {1, 3/2}, None} *)
Normal[qp["ElementTable"]]
(* {<|"Element" -> "BS", "Parameters" -> {th, ph}, "Modes" -> {1, 2}|>, <|"Element" -> "PS", "Parameters" -> {al}, "Modes" -> {2}|>} *)
```

Declared parameters bind positionally or by name, to numbers or to other symbols.
```wl
qp = QuantumInterferometer[{"BS"[th, ph]}, "Parameters" -> {th, ph}];
qp[Pi/4, 0]["Unitary"]
(* {{1/Sqrt[2], -1/Sqrt[2]}, {1/Sqrt[2], 1/Sqrt[2]}} *)
qp[<|th -> Pi/2|>]["Elements"]
(* {"BS"[Pi/2, ph] -> {1, 2}} *)
qp[th0, ph0]["Elements"]
(* {"BS"[th0, ph0] -> {1, 2}} *)
```

Fock-space circuit with d levels per mode.
```wl
Simplify[Normal[bs["CircuitOperator", 2]["Matrix"]]]
(* {{1, 0, 0, 0}, {0, Cos[th], E^(I ph) Sin[th], 0}, {0, -Sin[th]/E^(I ph), Cos[th], 0}, {0, 0, 0, Cos[2 th]}} *)
```

The diagram labels symbolic parameters.
```wl
Head[bs["Diagram", "ShowParameters" -> True]]
(* Graphics *)
```

Every property name.
```wl
bs["Properties"]
(* {"Unitary", "TransferMatrix", "Modes", "Method", "Elements", "Layers", "Parameters", "Depth", "BeamSplitterCount",
    "PhaseShifterCount", "SqueezerCount", "ElementTable", "Diagram", "MatrixPlot", "CircuitOperator", "Channel", "LossCount",
    "Amplitude", "Probability", "Probabilities", "Samples", "SampleCounts", "Dagger", "BlochMessiah", "GBSForm"} *)
```

## 2. Concrete photon numbers, symbolic parameters

```wl
bs["Amplitude", {1, 1} -> {2, 0}]
(* -Sqrt[2] Cos[th] Sin[th]/E^(I ph) *)

bs["Probability", {1, 1} -> {1, 1}]
(* (Cos[th]^2 - Sin[th]^2)^2 *)
```

Two identical beam splitters double the angle.
```wl
QuantumInterferometer[{"BS"[th, ph], "BS"[th, ph]}]["Amplitude", {1, 0} -> {0, 1}]
(* 2 E^(I ph) Cos[th] Sin[th] *)
```

All output amplitudes of an input.
```wl
bs[{1, 0}]
(* <|QuditName[{1, 0}] -> Cos[th], QuditName[{0, 1}] -> E^(I ph) Sin[th]|> *)
```

All output probabilities.
```wl
bs["Probabilities", {1, 1}]
(* <|QuditName[{2, 0}] -> 2 Cos[th]^2 Sin[th]^2, QuditName[{1, 1}] -> (Cos[th]^2 - Sin[th]^2)^2,
     QuditName[{0, 2}] -> 2 Cos[th]^2 Sin[th]^2|> *)
```

A pure QuantumState goes through the permanents; a mixed one through the Fock-space circuit.
```wl
Normal[bs[FockState[{1, 0}, 2]]["StateVector"]]
(* {0, E^(I ph) Sin[th], Cos[th], 0} *)
Normal[bs[QuantumState[DiagonalMatrix[{0, 1/2, 1/2, 0}], {2, 2}]]["DensityMatrix"]]
(* {{0, 0, 0, 0}, {0, (Cos[th]^2 + Sin[th]^2)/2, 0, 0}, {0, 0, (Cos[th]^2 + Sin[th]^2)/2, 0}, {0, 0, 0, 0}} *)
```

A state with its own parameters keeps them, so the output can be bound later.
```wl
Normal[bs[QuantumState[{0, Cos[g], Sin[g], 0}, {2, 2}, "Parameters" -> {g}]]["StateVector"]]
(* {0, Cos[g] Cos[th] + E^(I ph) Sin[g] Sin[th], Cos[th] Sin[g] - Cos[g] Sin[th]/E^(I ph), 0} *)
bs[QuantumState[{0, Cos[g], Sin[g], 0}, {2, 2}, "Parameters" -> {g}]]["Parameters"]
(* {g} *)
```

## 3. Symbolic photon numbers (passive)

One beam splitter is a Wigner D-function; conservation is a Kronecker delta.
```wl
bs["Amplitude", {n1, n2} -> {t1, t2}]
(* KroneckerDelta[n1 + n2, t1 + t2] WignerD[{(n1 + n2)/2, (n1 - n2)/2, (-n1 - n2 + 2 t1)/2}, ph, 2 th, -ph] *)

bs["Probability", {n1, n2} -> {t1, t2}]
(* Abs[WignerD[{(n1 + n2)/2, (n1 - n2)/2, (-n1 - n2 + 2 t1)/2}, ph, 2 th, -ph]]^2 KroneckerDelta[n1 + n2, t1 + t2] *)
```

Writing the output to conserve photons drops the delta.
```wl
bs["Amplitude", {n1, n2} -> {t1, n1 + n2 - t1}]
(* WignerD[{(n1 + n2)/2, (n1 - n2)/2, (-n1 - n2 + 2 t1)/2}, ph, 2 th, -ph] *)

bs["Probability", {n1, n2} -> {t1, n1 + n2 - t1}]
(* Abs[WignerD[{(n1 + n2)/2, (n1 - n2)/2, (-n1 - n2 + 2 t1)/2}, ph, 2 th, -ph]]^2 *)
```

Hong-Ou-Mandel, any n.
```wl
QuantumInterferometer[{"BS"[Pi/4, 0]}]["Amplitude", {n, n} -> {k, 2 n - k}]
(* WignerD[{n, 0, k - n}, 0, Pi/2, 0] *)
```

Phase shifter; its probability is the selection rule.
```wl
QuantumInterferometer[{"PS"[al]}]["Amplitude", {n} -> {t}]
(* E^(I al n) KroneckerDelta[n, t] *)

QuantumInterferometer[{"PS"[al]}]["Probability", {n} -> {t}]
(* KroneckerDelta[n, t] *)
```

A phase before the beam splitter is a factor in the photon number it acts on.
```wl
QuantumInterferometer[{"PS"[al] -> 2, "BS"[th, ph]}]["Amplitude", {n1, n2} -> {t1, n1 + n2 - t1}]
(* E^(I al n2) WignerD[{(n1 + n2)/2, (n1 - n2)/2, (-n1 - n2 + 2 t1)/2}, ph, 2 th, -ph] *)
```

Two beam splitters on the same modes: an inactive sum over the internal photon number K[1].
```wl
chain = QuantumInterferometer[{"BS"[th1, ph1], "BS"[th2, ph2]}];
chain["Amplitude", {n, 0} -> {t, n - t}]
(* Inactive[Sum][WignerD[{n/2, n/2, (-n + 2 K[1])/2}, ph1, 2 th1, -ph1] *
     WignerD[{n/2, (-n + 2 K[1])/2, (-n + 2 t)/2}, ph2, 2 th2, -ph2], {K[1], 0, n}] *)
```

Substitute photon numbers, then Activate; it agrees with the concrete query.
```wl
FullSimplify[Activate[chain["Amplitude", {n, 0} -> {t, n - t}] /. {n -> 2, t -> 1}] - chain["Amplitude", {2, 0} -> {1, 1}]]
(* 0 *)
```

Three modes: sum bounds follow the photon flow.
```wl
QuantumInterferometer[{"BS"[a1, b1] -> {1, 2}, "BS"[a2, b2] -> {2, 3}, "BS"[a3, b3] -> {1, 2}}, 3]["Amplitude", {n1, n2, 0} -> {t1, t2, t3}]
(* KroneckerDelta[n1 + n2, t1 + t2 + t3] Inactive[Sum][WignerD[{(n1 + n2)/2, (n1 - n2)/2, (-n1 - n2 + 2 K[1])/2}, b1, 2 a1, -b1] *
     WignerD[{(t1 + t2)/2, (-t1 - t2 + 2 K[1])/2, (t1 - t2)/2}, b3, 2 a3, -b3] *
     WignerD[{(n1 + n2 - K[1])/2, (n1 + n2 - K[1])/2, (-n1 - n2 + 2 t1 + 2 t2 - K[1])/2}, b2, 2 a2, -b2],
     {K[1], 0, Min[n1 + n2, t1 + t2]}] *)
```

Photons in mode 3 reach mode 1 only through mode 2, after the first beam splitter; this one condition is not implied by the photon numbers, so it stays.
```wl
QuantumInterferometer[{"BS"[th1, ph1] -> {1, 2}, "BS"[th2, ph2] -> {2, 3}}, 3]["Amplitude", {n1, n2, n3} -> {t1, t2, t3}]
(* KroneckerDelta[n1 + n2 + n3, t1 + t2 + t3] Piecewise[{{WignerD[{(n1 + n2)/2, (n1 - n2)/2, (-n1 - n2 + 2 t1)/2}, ph1, 2 th1, -ph1] *
     WignerD[{(n1 + n2 + n3 - t1)/2, (n1 + n2 - n3 - t1)/2, (-n1 - n2 - n3 + t1 + 2 t2)/2}, ph2, 2 th2, -ph2], n1 + n2 >= t1}}, 0] *)
```

Mixed concrete and symbolic photon numbers; an untouched mode keeps its photons.
```wl
QuantumInterferometer[{"BS"[th, ph] -> {2, 3}, "PS"[al] -> 1}, 3]["Amplitude", {n, 1, 0} -> {t1, t2, t3}]
(* E^(I al n) KroneckerDelta[n, t1] KroneckerDelta[1 + n, t1 + t2 + t3] WignerD[{1/2, 1/2, (-1 + 2 t2)/2}, ph, 2 th, -ph] *)
```

Inexact angles: the selection rules sit in a held condition, so off-shell patterns give 0, not Indeterminate.
```wl
QuantumInterferometer[{"BS"[0.3, 0.7]}]["Amplitude", {n1, n2} -> {t1, t2}]
(* Piecewise[{{WignerD[{(n1 + n2)/2, (n1 - n2)/2, (-n1 - n2 + 2 t1)/2}, 0.7, 0.6, -0.7], n1 + n2 == t1 + t2}}, 0] *)
% /. {n1 -> 1, n2 -> 0, t1 -> 2, t2 -> 0}
(* 0 *)
QuantumInterferometer[{"BS"[0.3, 0.7]}]["Probability", {n1, n2} -> {t1, t2}]
(* Abs[Piecewise[{{WignerD[...], n1 + n2 == t1 + t2}}, 0]]^2 *)
```

Sum indices skip any K[j] already in the query, the parameters or the assumptions.
```wl
chain["Amplitude", {K[1], 0} -> {t, K[1] - t}]
(* Inactive[Sum][WignerD[{K[1]/2, K[1]/2, (-K[1] + 2 K[2])/2}, ph1, 2 th1, -ph1] *
     WignerD[{K[1]/2, (-K[1] + 2 K[2])/2, (2 t - K[1])/2}, ph2, 2 th2, -ph2], {K[2], 0, K[1]}] *)
```

Patterns may be lists, Ket or QuditName.
```wl
bs["Amplitude", Ket[{n1, n2}] -> QuditName[{t1, n1 + n2 - t1}]]
(* WignerD[{(n1 + n2)/2, (n1 - n2)/2, (-n1 - n2 + 2 t1)/2}, ph, 2 th, -ph] *)
```

A photon number may be any polynomial with rational coefficients in photon-number variables: plain symbols, indexed symbols such as `n[1]`, subscripted ones, formal letters and one-letter system symbols such as `K` or `N`. Each such entry is taken to be a nonnegative integer.
```wl
bs["Amplitude", {n + 1, n} -> {t, 2 n + 1 - t}]
(* WignerD[{(1 + 2 n)/2, 1/2, (-1 - 2 n + 2 t)/2}, ph, 2 th, -ph] *)
bs["Amplitude", {2 n, 0} -> {n, n}]
(* WignerD[{n, n, 0}, ph, 2 th, -ph] *)
bs["Amplitude", {n[1], n[2]} -> {t[1], n[1] + n[2] - t[1]}]
(* WignerD[{(n[1] + n[2])/2, (n[1] - n[2])/2, (-n[1] - n[2] + 2 t[1])/2}, ph, 2 th, -ph] *)
bs["Amplitude", {Subscript[n, 1], 0} -> {k, Subscript[n, 1] - k}]
(* WignerD[{Subscript[n, 1]/2, Subscript[n, 1]/2, (2 k - Subscript[n, 1])/2}, ph, 2 th, -ph] *)
chain["Amplitude", {n + 1, n} -> {t, 2 n + 1 - t}]
(* Inactive[Sum][WignerD[{(1 + 2 n)/2, 1/2, (-1 - 2 n + 2 K[1])/2}, ph1, 2 th1, -ph1] *
     WignerD[{(1 + 2 n)/2, (-1 - 2 n + 2 K[1])/2, (-1 - 2 n + 2 t)/2}, ph2, 2 th2, -ph2], {K[1], 0, 1 + 2 n}] *)
```

With no beam splitter, every mode keeps its photons and phases multiply.
```wl
QuantumInterferometer[{}, 2]["Amplitude", {n1, n2} -> {t1, t2}]
(* KroneckerDelta[n1, t1] KroneckerDelta[n1 + n2, t1 + t2] *)
QuantumInterferometer[{"PS"[al] -> 1, "PS"[be] -> 2}, 2]["Amplitude", {n1, n2} -> {n1, n2}]
(* E^(I (al n1 + be n2)) *)
```

Four modes: one sum for the internal photon number of the two beam splitters on modes 1 and 2, and one path condition, since photons in modes 3 and 4 reach the output of mode 4 only through the splitter on modes 2 and 3.
```wl
QuantumInterferometer[{"BS"[a1, b1] -> {1, 2}, "BS"[a2, b2] -> {3, 4}, "BS"[a3, b3] -> {2, 3}, "BS"[a4, b4] -> {1, 2}}, 4]["Amplitude", {n1, n2, n3, n4} -> {t1, t2, t3, t4}]
(* KroneckerDelta[n1 + n2 + n3 + n4, t1 + t2 + t3 + t4] Piecewise[{{Inactive[Sum][
     WignerD[{(n1 + n2)/2, (n1 - n2)/2, (-n1 - n2 + 2 K[1])/2}, b1, 2 a1, -b1] *
     WignerD[{(n3 + n4)/2, (n3 - n4)/2, (-2 n1 - 2 n2 - n3 - n4 + 2 t1 + 2 t2 + 2 t3)/2}, b2, 2 a2, -b2] *
     WignerD[{(t1 + t2)/2, (-t1 - t2 + 2 K[1])/2, (t1 - t2)/2}, b4, 2 a4, -b4] *
     WignerD[{(t1 + t2 + t3 - K[1])/2, (2 n1 + 2 n2 - t1 - t2 - t3 - K[1])/2, (t1 + t2 - t3 - K[1])/2}, b3, 2 a3, -b3],
     {K[1], 0, Min[n1 + n2, t1 + t2]}], n3 + n4 >= t4}}, 0] *)
```

## 4. Assumptions

Optional; they add relations beyond the physical domain. Option, list and Assuming forms agree.
```wl
QuantumInterferometer[{"PS"[al]}]["Amplitude", {n} -> {n}, Assumptions -> al == Pi]
(* (-1)^n *)
Assuming[al == Pi, QuantumInterferometer[{"PS"[al]}]["Amplitude", {n} -> {n}]]
(* (-1)^n *)
QuantumInterferometer[{"PS"[al]}]["Amplitude", {n} -> {n}, Assumptions -> {al == Pi}]
(* (-1)^n *)
```

An assumed photon-flow relation removes the path condition.
```wl
QuantumInterferometer[{"BS"[th1, ph1] -> {1, 2}, "BS"[th2, ph2] -> {2, 3}}, 3]["Amplitude", {n1, n2, n3} -> {t1, t2, t3}, Assumptions -> n1 + n2 >= t1]
(* KroneckerDelta[n1 + n2 + n3, t1 + t2 + t3] WignerD[...] WignerD[...] *)
```

Restating the physical domain changes nothing.
```wl
bs["Amplitude", {n1, n2} -> {t1, t2}, Assumptions -> Element[{n1, n2, t1, t2}, Integers] && Element[{th, ph}, Reals]] === bs["Amplitude", {n1, n2} -> {t1, t2}]
(* True *)
```

Assumptions that contradict the domain fail.
```wl
bs["Amplitude", {n, 0} -> {t, n - t}, Assumptions -> n == 1/2]
(* $Failed, Refine::fas *)
```

Assumptions on squeezing parameters and coherent amplitudes, which the domain leaves complex, simplify them.
```wl
QuantumInterferometer[{"S"[xi]}]["Amplitude", {0} -> {2}, Assumptions -> xi > 0]
(* -Sech[xi]^(3/2) Sinh[xi]/Sqrt[2] *)
bs["Amplitude", "Coherent"[{a, 0}] -> {1, 0}, Assumptions -> a > 0]
(* a E^(-a^2 Cos[th]^2/2 - a^2 Sin[th]^2/2) Cos[th] *)
```

Whole distributions take assumptions the same three ways.
```wl
bs["Probabilities", {1, 1}, Assumptions -> th == Pi/4]
(* <|QuditName[{2, 0}] -> 1/2, QuditName[{1, 1}] -> 0, QuditName[{0, 2}] -> 1/2|> *)
Assuming[th == Pi/4, bs["Probabilities", {1, 1}]] === bs["Probabilities", {1, 1}, Assumptions -> {th == Pi/4}]
(* True *)
QuantumInterferometer[{"S"[xi]}]["Probabilities", {0}, "Cutoff" -> 2, Assumptions -> xi > 0]
(* <|QuditName[0] -> Sech[xi], QuditName[2] -> Sech[xi] Tanh[xi]^2/2|> *)
bs["Probabilities", {1, 1}, Assumptions -> th == I]
(* $Failed, Refine::fas *)
```

## 5. Coherent light (passive)

A coherent state stays coherent: |a> -> |W.a>.
```wl
bs["Coherent"[{a, b}]]
(* "Coherent"[{a Cos[th] - b Sin[th]/E^(I ph), b Cos[th] + a E^(I ph) Sin[th]}] *)
```

Coherent input, Fock output.
```wl
bs["Amplitude", "Coherent"[{a, 0}] -> {1, 0}]
(* a E^(-Abs[a Cos[th]]^2/2 - Abs[a Sin[th]]^2/2) Cos[th] *)
```

Fock input, coherent output.
```wl
bs["Amplitude", {1, 0} -> "Coherent"[{c, d}]]
(* E^((-Abs[c]^2 - Abs[d]^2)/2) (Conjugate[c] Cos[th] + E^(I ph) Conjugate[d] Sin[th]) *)
```

Coherent output from coherent input: the overlap probability is a Gaussian in the distance between the amplitudes.
```wl
bs["Probability", "Coherent"[{a, b}] -> "Coherent"[{c, d}]]
(* E^(-Abs[c - a Cos[th] + b Sin[th]/E^(I ph)]^2 - Abs[d - b Cos[th] - a E^(I ph) Sin[th]]^2) *)
Simplify[bs["Probability", "Coherent"[{a, 0}] -> "Coherent"[{a Cos[th], a E^(I ph) Sin[th]}]], Element[{th, ph}, Reals]]
(* 1 *)
```

Photon distribution up to a cutoff.
```wl
bs["Probabilities", "Coherent"[{a, 0}], "Cutoff" -> 1]
(* <|QuditName[{0, 0}] -> E^(-Abs[a Cos[th]]^2 - Abs[a Sin[th]]^2),
     QuditName[{1, 0}] -> E^(-Abs[a Cos[th]]^2 - Abs[a Sin[th]]^2) Abs[a Cos[th]]^2,
     QuditName[{0, 1}] -> Abs[a E^(I ph - Abs[a Cos[th]]^2/2 - Abs[a Sin[th]]^2/2) Sin[th]]^2|> *)
```

A coherent input may also be written mode by mode.
```wl
bs["Amplitude", {"Coherent"[a], "Coherent"[b]} -> {1, 0}] === bs["Amplitude", "Coherent"[{a, b}] -> {1, 0}]
(* True *)
```

Mixed Fock and coherent input.
```wl
bs["Amplitude", {1, "Coherent"[a]} -> {1, 0}]
(* E^(-Abs[a Cos[th]]^2/2 - Abs[a Sin[th]]^2/2) (1 - Abs[a Sin[th]]^2) Cos[th] +
   a E^(-Abs[a Cos[th]]^2/2 - Abs[a Sin[th]]^2/2) Conjugate[a] Cos[th] Sin[th]^2 *)
```

Its distribution starts at the Fock photons, since the coherent part only adds photons.
```wl
Keys[bs["Probabilities", {1, "Coherent"[a]}, "Cutoff" -> 2]]
(* {QuditName[{1, 0}], QuditName[{0, 1}], QuditName[{2, 0}], QuditName[{1, 1}], QuditName[{0, 2}]} *)
```

The all-outputs form takes coherent slots and a cutoff, and gives the amplitudes whose squares are that distribution.
```wl
bs[{"Coherent"[a], "Coherent"[b]}, "Cutoff" -> 1]
(* <|QuditName[{0, 0}] -> E^(-Abs[a Cos[th] - b Sin[th]/E^(I ph)]^2/2 - Abs[b Cos[th] + a E^(I ph) Sin[th]]^2/2),
     QuditName[{1, 0}] -> E^(-Abs[a Cos[th] - b Sin[th]/E^(I ph)]^2/2 - Abs[b Cos[th] + a E^(I ph) Sin[th]]^2/2) (a Cos[th] - b Sin[th]/E^(I ph)),
     QuditName[{0, 1}] -> E^(-Abs[a Cos[th] - b Sin[th]/E^(I ph)]^2/2 - Abs[b Cos[th] + a E^(I ph) Sin[th]]^2/2) (b Cos[th] + a E^(I ph) Sin[th])|> *)
Union[Values[FullSimplify[(Abs[#]^2 & /@ bs[{1, "Coherent"[a]}, "Cutoff" -> 2]) - bs["Probabilities", {1, "Coherent"[a]}, "Cutoff" -> 2], Element[{th, ph}, Reals]]]]
(* {0} *)
```

## 6. Squeezing

```wl
sq = QuantumInterferometer[{"S"[xi]}];
sq["Amplitude", {0} -> {2}]
(* -E^(I Arg[xi]) Sech[Abs[xi]]^(3/2) Sinh[Abs[xi]]/Sqrt[2] *)

sq["Probability", {0} -> {2}]
(* Sech[Abs[xi]] Tanh[Abs[xi]]^2/2 *)

sq["Amplitude", {1} -> {3}]
(* -Sqrt[3/2] E^(I Arg[xi]) Sech[Abs[xi]]^(5/2) Sinh[Abs[xi]] *)
```

All amplitudes or probabilities up to a cutoff on the total photon number.
```wl
sq[{0}, "Cutoff" -> 2]
(* <|QuditName[0] -> Sqrt[Sech[Abs[xi]]], QuditName[2] -> -E^(I Arg[xi]) Sech[Abs[xi]]^(3/2) Sinh[Abs[xi]]/Sqrt[2]|> *)

sq["Probabilities", {0}, "Cutoff" -> 4]
(* <|QuditName[0] -> Sech[Abs[xi]], QuditName[2] -> Sech[Abs[xi]] Tanh[Abs[xi]]^2/2,
     QuditName[4] -> 3 Sech[Abs[xi]] Tanh[Abs[xi]]^4/8|> *)
```

A QuantumState goes through the Fock-space circuit of SqueezeOperator.
```wl
Normal[sq[FockState[{0}, 3]]["StateVector"]]
(* {1/Sqrt[Cosh[Abs[xi]]], 0, -E^(I Arg[xi]) Sinh[Abs[xi]]/(Sqrt[2] Cosh[Abs[xi]]^(3/2))} *)
Normal[sq[QuantumState[DiagonalMatrix[{1/2, 1/2, 0}], 3]]["DensityMatrix"]]
(* {{Sech[Abs[xi]]/2, 0, -Sech[Abs[xi]] Tanh[Abs[xi]]/(2 Sqrt[2] E^(I Arg[xi]))}, {0, Sech[Abs[xi]]^3/2, 0},
    {-E^(I Arg[xi]) Sech[Abs[xi]] Tanh[Abs[xi]]/(2 Sqrt[2]), 0, Sech[Abs[xi]] Tanh[Abs[xi]]^2/4}} *)
```

Squeezer then a balanced beam splitter.
```wl
Simplify[QuantumInterferometer[{"S"[r], "BS"[Pi/4, 0]}]["Amplitude", {0, 0} -> {1, 1}], r > 0]
(* -Sech[r]^(3/2) Sinh[r]/2 *)
```

Two squeezers interfering: equal squeezing leaves a product of single-mode squeezed vacua, so the coincidence vanishes; opposite squeezing, x2 = -x1, makes a two-mode squeezed vacuum and the coincidence is largest.
```wl
QuantumInterferometer[{"S"[x1] -> 1, "S"[x2] -> 2, "BS"[Pi/4, 0]}]["Amplitude", {0, 0} -> {1, 1}]
(* Sqrt[Sech[Abs[x1]] Sech[Abs[x2]]] (-E^(I Arg[x1]) Tanh[Abs[x1]]/2 + E^(I Arg[x2]) Tanh[Abs[x2]]/2) *)
```

A squeezer before a general beam splitter; the denominator is Cosh[Abs[xi]] once Cos^2 + Sin^2 = 1 is applied.
```wl
QuantumInterferometer[{"S"[xi] -> 1, "BS"[th, ph]}]["Amplitude", {0, 0} -> {1, 1}]
(* -E^(I ph + I Arg[xi]) Cos[th] Sqrt[Sech[Abs[xi]]] Sin[th] Sinh[Abs[xi]]/(Cos[th]^2 Cosh[Abs[xi]] + Cosh[Abs[xi]] Sin[th]^2) *)
```

Counts and inverse; the inverse undoes the squeezer.
```wl
{sq["SqueezerCount"], sq["Dagger"]["Elements"]}
(* {1, {"S"[-xi] -> {1}}} *)
FullSimplify[{sq["Dagger"][sq]["Amplitude", {0} -> {0}], sq["Dagger"][sq]["Amplitude", {1} -> {3}]}]
(* {1, 0} *)
```

The Fock-space circuit with 3 levels; its corner entry is the truncated SqueezeOperator, not an amplitude of the network.
```wl
Simplify[Normal[sq["CircuitOperator", 3]["Matrix"]]]
(* {{1/Sqrt[Cosh[Abs[xi]]], 0, Sech[Abs[xi]]^(3/2) Sinh[Abs[xi]]/(Sqrt[2] E^(I Arg[xi]))}, {0, Cosh[Abs[xi]]^(-3/2), 0},
    {-E^(I Arg[xi]) Sech[Abs[xi]]^(3/2) Sinh[Abs[xi]]/Sqrt[2], 0, -1/4 (-5 + Cosh[2 Abs[xi]])/Cosh[Abs[xi]]^(5/2)}} *)
```

A declared squeezing parameter binds like any other.
```wl
QuantumInterferometer[{"S"[r]}, "Parameters" -> {r}][s0]["Elements"]
(* {"S"[s0] -> {1}} *)
```

## 7. Pure loss

Retained transfer matrix.
```wl
l = QuantumInterferometer[{"Loss"[eta]}];
l["TransferMatrix"]
(* {{Sqrt[eta]}} *)
bs[QuantumInterferometer[{"Loss"[eta] -> 2}, 2]]["TransferMatrix"]
(* {{Cos[th], -Sqrt[eta] Sin[th]/E^(I ph)}, {E^(I ph) Sin[th], Sqrt[eta] Cos[th]}} *)
```

Count law: a bare binomial. Phases do not change it.
```wl
l["Probability", {n} -> {k}]
(* (1 - eta)^(-k + n) eta^k Binomial[n, k] *)
QuantumInterferometer[{"PS"[al], "Loss"[eta]}]["Probability", {n} -> {k}]
(* (1 - eta)^(-k + n) eta^k Binomial[n, k] *)
l["Probability", {3} -> {k}]
(* (1 - eta)^(3 - k) eta^k Binomial[3, k] *)
l["Probability", {1} -> {0}]
(* 1 - eta *)
l["Probabilities", {2}]
(* <|QuditName[0] -> (1 - eta)^2, QuditName[1] -> 2 (1 - eta) eta, QuditName[2] -> eta^2|> *)
```

Independent attenuators: rates on one mode multiply; an untouched mode conserves its photons.
```wl
QuantumInterferometer[{"Loss"[eta] -> 1, "Loss"[mu] -> 1, "Loss"[nu] -> 2}, 3]["Probability", {n1, n2, n3} -> {k1, k2, k3}]
(* (eta mu)^k1 (1 - eta mu)^(-k1 + n1) (1 - nu)^(-k2 + n2) nu^k2 Binomial[n1, k1] Binomial[n2, k2] KroneckerDelta[k3, n3] *)
```

Endpoints: a numeric rate or an assumption gives the endpoint law; a symbolic formula reaches it as a limit.
```wl
{QuantumInterferometer[{"Loss"[0]}]["Probability", {n} -> {k}], l["Probability", {n} -> {k}, Assumptions -> eta == 1]}
(* {KroneckerDelta[0, k], KroneckerDelta[k, n]} *)
Assuming[eta == 1, l["Probability", {n} -> {k}]]
(* KroneckerDelta[k, n] *)
Table[Limit[l["Probability", {2} -> {k}], eta -> 1], {k, 0, 2}]
(* {0, 0, 1} *)
```

Loss with mode mixing, concrete photon numbers.
```wl
lb = QuantumInterferometer[{"BS"[th], "Loss"[eta] -> 1}];
lb["Probability", {1, 0} -> {0, 0}]
(* (1 - eta) Cos[th]^2 *)
lb["Probabilities", {1, 0}]
(* <|QuditName[{0, 0}] -> Cos[th]^2 - eta Cos[th]^2, QuditName[{1, 0}] -> eta Cos[th]^2, QuditName[{0, 1}] -> Sin[th]^2|> *)
lb["TransferMatrix"]
(* {{Sqrt[eta] Cos[th], -Sqrt[eta] Sin[th]}, {Sin[th], Cos[th]}} *)
FullSimplify[Total[lb["Probabilities", {2, 1}]]]
(* 1 *)
```

Loss between two beam splitters: in each probability the two phase terms are complex conjugates of each other, so it is visibly real, and the probabilities total 1 on the physical domain.
```wl
l2 = QuantumInterferometer[{"BS"[th, ph], "Loss"[eta] -> 1, "BS"[th2, ph2]}];
l2["Probabilities", {1, 0}]
(* <|QuditName[{0, 0}] -> Cos[th]^2 - eta Cos[th]^2,
     QuditName[{1, 0}] -> eta Cos[th]^2 Cos[th2]^2 - Sqrt[eta] Cos[th] Cos[th2] Sin[th] Sin[th2]/E^(I (ph - ph2)) - E^(I ph - I ph2) Sqrt[eta] Cos[th] Cos[th2] Sin[th] Sin[th2] + Sin[th]^2 Sin[th2]^2,
     QuditName[{0, 1}] -> Cos[th2]^2 Sin[th]^2 + Sqrt[eta] Cos[th] Cos[th2] Sin[th] Sin[th2]/E^(I (ph - ph2)) + E^(I ph - I ph2) Sqrt[eta] Cos[th] Cos[th2] Sin[th] Sin[th2] + eta Cos[th]^2 Sin[th2]^2|> *)
FullSimplify[Total[l2["Probabilities", {1, 0}]], Element[{th, ph, th2, ph2}, Reals] && 0 <= eta <= 1]
(* 1 *)
```

Coherent light stays pure through loss: the amplitude shrinks by Sqrt[eta].
```wl
l["Coherent"[{a}]]
(* "Coherent"[{a Sqrt[eta]}] *)
lb["Coherent"[{a, b}]]
(* "Coherent"[{a Sqrt[eta] Cos[th] - b Sqrt[eta] Sin[th], b Cos[th] + a Sin[th]}] *)
l["Probabilities", "Coherent"[{a}], "Cutoff" -> 1]
(* <|QuditName[0] -> E^(-eta Abs[a]^2), QuditName[1] -> eta Abs[a]^2/E^(eta Abs[a]^2)|> *)
l["Probability", "Coherent"[{a}] -> {2}]
(* eta^2 Abs[a]^4/(2 E^(eta Abs[a]^2)) *)
lb["Probability", "Coherent"[{a, 0}] -> {1, 0}]
(* E^(-eta Abs[a Cos[th]]^2 - Abs[a Sin[th]]^2) eta Abs[a Cos[th]]^2 *)
l["Probability", "Coherent"[{a}] -> "Coherent"[{b}]]
(* E^(-Abs[b - a Sqrt[eta]]^2) *)
```

Because the output is pure, a coherent input through loss also has amplitudes; the mode-by-mode form works here too.
```wl
l["Amplitude", "Coherent"[{a}] -> {1}]
(* a Sqrt[eta]/E^(eta Abs[a]^2/2) *)
l["Probability", {"Coherent"[a]} -> {2}]
(* eta^2 Abs[a]^4/(2 E^(eta Abs[a]^2)) *)
```

Finite states: coherences scale by Sqrt[eta].
```wl
Normal[l[QuantumState[{1, I}/Sqrt[2], 2]]["DensityMatrix"]]
(* {{1/2 + (1 - eta)/2, -I Sqrt[eta]/2}, {I Sqrt[eta]/2, eta/2}} *)
Simplify[Normal[lb[FockState[{1, 0}, 2]]["DensityMatrix"]]]
(* {{-((-1 + eta) Cos[th]^2), 0, 0, 0}, {0, Sin[th]^2, Sqrt[eta] Cos[th] Sin[th], 0}, {0, Sqrt[eta] Cos[th] Sin[th], eta Cos[th]^2, 0}, {0, 0, 0, 0}} *)
```

The channel is a native QuantumChannel; QF applies it, so its output is not simplified by QI.
```wl
Normal[l["Channel", 2][QuantumState[{0, 1}, 2]]["DensityMatrix"]]
(* {{Sqrt[1 - eta] Conjugate[Sqrt[1 - eta]], 0}, {0, Sqrt[eta] Conjugate[Sqrt[eta]]}} *)
```

Binding eta -> 1 restores a lossless, unitary network; binding a number gives exact counts.
```wl
lq = QuantumInterferometer[{"BS"[th], "Loss"[eta] -> 2}, "Parameters" -> {th, eta}, "DropIdentities" -> False];
{lq["LossCount"], lq[Pi/4, 1]["Unitary"]}
(* {1, {{1/Sqrt[2], -1/Sqrt[2]}, {1/Sqrt[2], 1/Sqrt[2]}}} *)
QuantumInterferometer[{"Loss"[eta]}, "Parameters" -> {eta}][1/2]["Probability", {2} -> {1}]
(* 1/2 *)
```

The diagram labels symbolic loss.
```wl
Head[l["Diagram", "ShowParameters" -> True]]
(* Graphics *)
```

Assumptions outside the domain, or contradicting each other, fail.
```wl
l["Probability", {n} -> {k}, Assumptions -> eta > 1]
(* $Failed, Refine::fas *)
l["Probability", {n} -> {k}, Assumptions -> {eta == 0, eta == 1}]
(* $Failed, Refine::cas *)
```

## 8. Not supported symbolically

Symbolic photon numbers work only in a passive "Amplitude" or "Probability" query, and in a loss "Probability" without mode mixing. Entries that are not polynomials with rational coefficients in photon-number variables are rejected.
```wl
bs["Amplitude", {Sqrt[n], 0} -> {t, 0}]                            (* $Failed, ::occ *)
bs["Amplitude", {n + I, 0} -> {t, 0}]                              (* $Failed, ::occ *)
bs["Amplitude", {n + Pi, 0} -> {t, 0}]                             (* $Failed, ::occ *)
bs["Amplitude", {Abs[n], 0} -> {t, 0}]                             (* $Failed, ::occ *)
bs["Amplitude", {n, 1/2} -> {t, 0}]                                (* $Failed, ::occ *)
bs[{n1, n2}]                                                       (* $Failed, ::occ *)
bs["Probabilities", {n1, n2}]                                      (* $Failed, ::occ *)
bs["Amplitude", "Coherent"[{a, 0}] -> {t1, t2}]                     (* $Failed, ::occ *)
bs["Amplitude", {n, 0} -> "Coherent"[{c, d}]]                       (* $Failed, ::occ *)
bs["Amplitude", {nx, Infinity} -> {tx, 0}]                          (* $Failed, ::occ *)
sq["Amplitude", {n} -> {n}]                                         (* $Failed, ::squeezed *)
l["Probabilities", {n}]                                             (* $Failed, ::occ *)
l["Probability", "Coherent"[{a}] -> {k}]                            (* $Failed, ::occ *)
lb["Probability", {n, 0} -> {k, 0}]                                 (* $Failed, ::lossoccupation *)
QuantumInterferometer[{"BS"[th] -> {1, 2}, "Loss"[eta] -> 3}, 3]["Probability", {1, 0, n} -> {0, 1, k}]
                                                                    (* $Failed, ::lossoccupation *)
```

Automatic cutoffs need numbers.
```wl
bs["Probabilities", "Coherent"[{a, 0}]]                            (* $Failed, ::symcutoff *)
bs[{"Coherent"[a], "Coherent"[b]}]                                 (* $Failed, ::symcutoff *)
sq["Probabilities", {0}]                                           (* $Failed, ::symcutoff *)
```

Active and lossy networks have no unitary; some properties need numbers.
```wl
sq["Unitary"]                                                      (* $Failed, ::squeezed *)
sq["TransferMatrix"]                                               (* $Failed, ::squeezed *)
sq["Amplitude", "Coherent"[{a}] -> {0}]                            (* $Failed, ::squeezed *)
sq["BlochMessiah"]                                                 (* $Failed, ::notnumeric *)
sq["GBSForm"]                                                      (* $Failed, ::notnumeric *)
l["Unitary"]                                                       (* $Failed, ::lossy *)
l["Dagger"]                                                        (* $Failed, ::lossy *)
l["Amplitude", {n} -> {k}]                                         (* $Failed, ::lossy *)
QuantumInterferometer[{"S"[xi], "Loss"[eta]}]                      (* $Failed, ::lossactive *)
QuantumInterferometer[{{Cos[th], -Sin[th]}, {Sin[th], Cos[th]}}]   (* $Failed, ::symbolic *)
bs["MatrixPlot"]                                                   (* $Failed, ::notnumeric *)
QuantumInterferometer[bs]                                          (* $Failed, ::symbolic *)
QuantumInterferometer[sq]                                          (* $Failed, ::squeezed *)
sq["Coherent"[{a}]]                                                (* $Failed, ::squeezed *)
QuantumInterferometer[{"Loss"[eta], "BS"[th]}]["BlochMessiah"]     (* $Failed, ::lossy *)
```

Sampling needs numbers, lossy or not.
```wl
bs["Samples", {1, 0}, 1]                                           (* $Failed, ::notnumeric *)
l["Samples", {1}, 2]                                               (* $Failed, ::notnumeric *)
bs["Samples", "Coherent"[{a, 0}], 2]                               (* $Failed, ::coherentnumeric *)
```

## Gaps seen while running this list

Substitute photon numbers before machine values for symbolic angles: a numeric WignerD outside its index range is Indeterminate. Inexact angles given at construction are safe.
```wl
bs["Amplitude", {n1, n2} -> {t1, t2}] /. {th -> 0.3, ph -> 0.7} /. {n1 -> 1, n2 -> 0, t1 -> 2, t2 -> 0}
(* Indeterminate, Infinity::indet *)
bs["Amplitude", {n1, n2} -> {t1, t2}] /. {n1 -> 1, n2 -> 0, t1 -> 2, t2 -> 0} /. {th -> 0.3, ph -> 0.7}
(* 0 *)
QuantumInterferometer[{"BS"[0.3, 0.7]}]["Amplitude", {n1, n2} -> {t1, t2}] /. {n1 -> 1, n2 -> 0, t1 -> 2, t2 -> 0}
(* 0 *)
```

The all-outputs form `qi[pattern]`, with Fock photons or a coherent slot, takes no `Assumptions` and ignores an enclosing `Assuming`; `"Amplitude"` and `"Probabilities"` honour both.

Simplification the kernel leaves to the caller:
- `qi["Channel", d]` and the `"CircuitOperator"` of a lossy network apply QF's QuantumChannel, so their density matrices keep `Conjugate`, as in `Sqrt[eta] Conjugate[Sqrt[eta]]`. `qi[state]` gives the simplified form for every network and state.
- `Refine` does not apply trigonometric identities, so a mixed state through a beam splitter keeps `(Cos[th]^2 + Sin[th]^2)/2`, and a squeezer before a general beam splitter keeps `Cos[th]^2 Cosh[Abs[xi]] + Cosh[Abs[xi]] Sin[th]^2`.
- Coherent photon distributions keep `Abs[a Cos[th]]^2` rather than `Abs[a]^2 Cos[th]^2`; `Refine` does not factor the modulus.
- A squeezing parameter written in polar form keeps `Arg[E^(I f)]`, since `Arg` reduces f modulo 2 Pi:
```wl
QuantumInterferometer[{"S"[r Exp[I f]]}]["Amplitude", {0} -> {2}, Assumptions -> r > 0 && Element[f, Reals]]
(* -E^(I Arg[E^(I f)]) Sech[r]^(3/2) Sinh[r]/Sqrt[2] *)
```
