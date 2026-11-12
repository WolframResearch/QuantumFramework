# QuantumInterferometer: symbolic behaviour

Every statement below runs on the installed kernel (WL 15.0) with QuantumFramework `main` at c74d5911 (2026-10-28). That kernel includes the physical-domain contract (c1901374), the compact squeezed, coherent and mixed-state results (23ac0cf1), lossy and coherent sampling (c91ca9ad, 24f5679f), the photon-number allowlist (46f9ebf6), and real lossy probabilities, `Assumptions` in `"Probabilities"` and coherent slots in the all-outputs form (ed916377), and squeezed-vacuum sampling (c74d5911). The usage message of `QuantumInterferometer` now lists only the call forms, so this page is the reference for symbolic behaviour.

**Physical domain.** Photon numbers are nonnegative integers, beam-splitter angles θ, φ and phase angles α are real, loss transmissivities η are real with 0 ≤ η ≤ 1, coherent amplitudes and squeezing parameters are complex. Symbolic results hold on this domain and carry no condition that restates it.

```wl
bs = QuantumInterferometer[{"BS"[th, ph]}];
```

## 1. Symbolic element parameters

Transfer matrix of a beam splitter; an omitted phase is 0.
```wl
bs["Unitary"]
QuantumInterferometer[{"BS"[th]}]["Unitary"]
```

Phase shifter.
```wl
QuantumInterferometer[{"PS"[al]}]["Unitary"]
```

Long element names normalize to the short ones; a beam splitter also takes its angles as a list, and its two modes in either order.
```wl
QuantumInterferometer[{"PhaseShift"[al], "Squeeze"[xi]}]["Elements"]
QuantumInterferometer[{"BeamSplitter"[{th, ph}]}]["Elements"]
QuantumInterferometer[{"BS"[th, ph] -> {2, 1}}]["Unitary"]
```

Elements that are exactly the identity are dropped, symbolic ones never; `"DropIdentities" -> False` keeps them.
```wl
QuantumInterferometer[{"BS"[0, ph], "PS"[2 Pi] -> 2, "PS"[al] -> 1}]["Elements"]
QuantumInterferometer[{"BS"[0, ph], "PS"[al] -> 1}, "DropIdentities" -> False]["Elements"]
```

Elements compose in time order.
```wl
QuantumInterferometer[{"BS"[th, ph], "PS"[al] -> 2}]["Unitary"]
```

Composition `qi2[qi1]` is qi1 followed by qi2.
```wl
QuantumInterferometer[{"PS"[al] -> 2}, 2][bs]["Unitary"]
```

Inverse.
```wl
bs["Dagger"]["Unitary"]
```

A passive network's transfer matrix is its unitary.
```wl
bs["TransferMatrix"] === bs["Unitary"]
```

Structure does not depend on the parameter values.
```wl
qp = QuantumInterferometer[{"BS"[th, ph], "PS"[al] -> 2}, "Parameters" -> {th, ph, al}];
{qp["Modes"], qp["Depth"], qp["BeamSplitterCount"], qp["PhaseShifterCount"], qp["Parameters"], qp["Layers"], qp["Method"]}
Normal[qp["ElementTable"]]
```

Declared parameters bind positionally or by name, to numbers or to other symbols.
```wl
qp = QuantumInterferometer[{"BS"[th, ph]}, "Parameters" -> {th, ph}];
qp[Pi/4, 0]["Unitary"]
qp[<|th -> Pi/2|>]["Elements"]
qp[th0, ph0]["Elements"]
qp[Pi/4]["Elements"]
qp[0]["Elements"]
qp["Dagger"]["Parameters"]
```

`N` leaves a symbolic network as it is, and a passive network is its own Bloch-Messiah and GBS form.
```wl
N[bs] === bs
bs["BlochMessiah"] === bs
bs["GBSForm"] === bs
```

An exact unitary decomposes into exact angles, and an exact network re-meshes exactly; the cosine-sine mesh and `"Fourier"[m]` are numeric.
```wl
QuantumInterferometer[{{Cos[Pi/7], -Sin[Pi/7]}, {Sin[Pi/7], Cos[Pi/7]}}]["Elements"]
QuantumInterferometer[FourierMatrix[2]]["Elements"]
QuantumInterferometer[FourierMatrix[2], Method -> "Reck"]["Elements"]
qf = QuantumInterferometer[FourierMatrix[3]];
AllTrue[Flatten[QuantumInterferometer[qf["Elements"], 3]["Unitary"] - FourierMatrix[3]], PossibleZeroQ]
qf["Amplitude", {1, 1, 1} -> {3, 0, 0}]
QuantumInterferometer[QuantumInterferometer[{"BS"[Pi/5, Pi/7]}]]["Elements"]
Precision[QuantumInterferometer[FourierMatrix[2], Method -> "CosineSine"]["Elements"]]
Precision[QuantumInterferometer["Fourier"[2]]["Unitary"]]
```

Fock-space circuit with d levels per mode.
```wl
Simplify[Normal[bs["CircuitOperator", 2]["Matrix"]]]
```

The diagram labels symbolic parameters.
```wl
Head[bs["Diagram", "ShowParameters" -> True]]
Head[QuantumInterferometer[{"S"[xi]}]["Diagram", "ShowParameters" -> True]]
```

Every property name.
```wl
bs["Properties"]
```

## 2. Concrete photon numbers, symbolic parameters

```wl
bs["Amplitude", {1, 1} -> {2, 0}]
bs["Probability", {1, 1} -> {1, 1}]
```

Two identical beam splitters double the angle.
```wl
QuantumInterferometer[{"BS"[th, ph], "BS"[th, ph]}]["Amplitude", {1, 0} -> {0, 1}]
```

All output amplitudes of an input, given as a list, Ket or QuditName.
```wl
bs[{1, 0}]
bs[Ket[{1, 0}]] === bs[QuditName[{1, 0}]] === bs[{1, 0}]
```

All output probabilities.
```wl
bs["Probabilities", {1, 1}]
```

A pure QuantumState goes through the permanents; a mixed one through the Fock-space circuit.
```wl
Normal[bs[FockState[{1, 0}, 2]]["StateVector"]]
Normal[bs[QuantumState[DiagonalMatrix[{0, 1/2, 1/2, 0}], {2, 2}]]["DensityMatrix"]]
```

Too few Fock levels drop the bunched outputs, with a warning; the norm left is the coincidence amplitude.
```wl
bs[FockState[{1, 1}, 2]]["Norm"]
```

A state with its own parameters keeps them, so the output can be bound later.
```wl
Normal[bs[QuantumState[{0, Cos[g], Sin[g], 0}, {2, 2}, "Parameters" -> {g}]]["StateVector"]]
bs[QuantumState[{0, Cos[g], Sin[g], 0}, {2, 2}, "Parameters" -> {g}]]["Parameters"]
```

## 3. Symbolic photon numbers (passive)

One beam splitter is a Wigner D-function; conservation is a Kronecker delta.
```wl
bs["Amplitude", {n1, n2} -> {t1, t2}]
bs["Probability", {n1, n2} -> {t1, t2}]
```

Writing the output to conserve photons drops the delta.
```wl
bs["Amplitude", {n1, n2} -> {t1, n1 + n2 - t1}]
bs["Probability", {n1, n2} -> {t1, n1 + n2 - t1}]
```

Hong-Ou-Mandel, any n.
```wl
QuantumInterferometer[{"BS"[Pi/4, 0]}]["Amplitude", {n, n} -> {k, 2 n - k}]
```

Phase shifter; its probability is the selection rule.
```wl
QuantumInterferometer[{"PS"[al]}]["Amplitude", {n} -> {t}]
QuantumInterferometer[{"PS"[al]}]["Probability", {n} -> {t}]
```

A phase before the beam splitter is a factor in the photon number it acts on.
```wl
QuantumInterferometer[{"PS"[al] -> 2, "BS"[th, ph]}]["Amplitude", {n1, n2} -> {t1, n1 + n2 - t1}]
```

Two beam splitters on the same modes: an inactive sum over the internal photon number K[1].
```wl
chain = QuantumInterferometer[{"BS"[th1, ph1], "BS"[th2, ph2]}];
chain["Amplitude", {n, 0} -> {t, n - t}]
```

Substitute photon numbers, then Activate; it agrees with the concrete query.
```wl
FullSimplify[Activate[chain["Amplitude", {n, 0} -> {t, n - t}] /. {n -> 2, t -> 1}] - chain["Amplitude", {2, 0} -> {1, 1}]]
```

Three modes: sum bounds follow the photon flow.
```wl
QuantumInterferometer[{"BS"[a1, b1] -> {1, 2}, "BS"[a2, b2] -> {2, 3}, "BS"[a3, b3] -> {1, 2}}, 3]["Amplitude", {n1, n2, 0} -> {t1, t2, t3}]
```

Photons in mode 3 reach mode 1 only through mode 2, after the first beam splitter; this one condition is not implied by the photon numbers, so it stays.
```wl
QuantumInterferometer[{"BS"[th1, ph1] -> {1, 2}, "BS"[th2, ph2] -> {2, 3}}, 3]["Amplitude", {n1, n2, n3} -> {t1, t2, t3}]
```

Mixed concrete and symbolic photon numbers; an untouched mode keeps its photons.
```wl
QuantumInterferometer[{"BS"[th, ph] -> {2, 3}, "PS"[al] -> 1}, 3]["Amplitude", {n, 1, 0} -> {t1, t2, t3}]
```

Inexact angles: the selection rules sit in a held condition, so off-shell patterns give 0, not Indeterminate.
```wl
QuantumInterferometer[{"BS"[0.3, 0.7]}]["Amplitude", {n1, n2} -> {t1, t2}]
% /. {n1 -> 1, n2 -> 0, t1 -> 2, t2 -> 0}
QuantumInterferometer[{"BS"[0.3, 0.7]}]["Probability", {n1, n2} -> {t1, t2}]
QuantumInterferometer[{"BS"[N[3/10, 30], N[7/11, 30]]}]["Amplitude", {n1, n2} -> {t1, n1 + n2 - t1}]
```

Sum indices skip any K[j] already in the query, the parameters or the assumptions.
```wl
chain["Amplitude", {K[1], 0} -> {t, K[1] - t}]
QuantumInterferometer[{"BS"[K[5], 0], "BS"[th2, ph2]}]["Amplitude", {n, 0} -> {t, n - t}]
```

Patterns may be lists, Ket or QuditName.
```wl
bs["Amplitude", Ket[{n1, n2}] -> QuditName[{t1, n1 + n2 - t1}]]
```

A photon number may be any polynomial with rational coefficients in photon-number variables: plain symbols, indexed symbols such as `n[1]`, subscripted ones, formal letters and one-letter system symbols such as `K` or `N`. Each such entry is taken to be a nonnegative integer.
```wl
bs["Amplitude", {n + 1, n} -> {t, 2 n + 1 - t}]
bs["Amplitude", {2 n, 0} -> {n, n}]
bs["Amplitude", {n[1], n[2]} -> {t[1], n[1] + n[2] - t[1]}]
bs["Amplitude", {Subscript[n, 1], 0} -> {k, Subscript[n, 1] - k}]
chain["Amplitude", {n + 1, n} -> {t, 2 n + 1 - t}]
```

With no beam splitter, every mode keeps its photons and phases multiply.
```wl
QuantumInterferometer[{}, 2]["Amplitude", {n1, n2} -> {t1, t2}]
QuantumInterferometer[{"PS"[al] -> 1, "PS"[be] -> 2}, 2]["Amplitude", {n1, n2} -> {n1, n2}]
```

Four modes: one sum for the internal photon number of the two beam splitters on modes 1 and 2, and one path condition, since photons in modes 3 and 4 reach the output of mode 4 only through the splitter on modes 2 and 3.
```wl
QuantumInterferometer[{"BS"[a1, b1] -> {1, 2}, "BS"[a2, b2] -> {3, 4}, "BS"[a3, b3] -> {2, 3}, "BS"[a4, b4] -> {1, 2}}, 4]["Amplitude", {n1, n2, n3, n4} -> {t1, t2, t3, t4}]
```

## 4. Assumptions

Optional; they add relations beyond the physical domain. Option, list and Assuming forms agree.
```wl
QuantumInterferometer[{"PS"[al]}]["Amplitude", {n} -> {n}, Assumptions -> al == Pi]
Assuming[al == Pi, QuantumInterferometer[{"PS"[al]}]["Amplitude", {n} -> {n}]]
QuantumInterferometer[{"PS"[al]}]["Amplitude", {n} -> {n}, Assumptions -> {al == Pi}]
```

An assumed photon-flow relation removes the path condition.
```wl
QuantumInterferometer[{"BS"[th1, ph1] -> {1, 2}, "BS"[th2, ph2] -> {2, 3}}, 3]["Amplitude", {n1, n2, n3} -> {t1, t2, t3}, Assumptions -> n1 + n2 >= t1]
```

Restating the physical domain changes nothing.
```wl
bs["Amplitude", {n1, n2} -> {t1, t2}, Assumptions -> Element[{n1, n2, t1, t2}, Integers] && Element[{th, ph}, Reals]] === bs["Amplitude", {n1, n2} -> {t1, t2}]
```

Assumptions that contradict the domain fail.
```wl
bs["Amplitude", {n, 0} -> {t, n - t}, Assumptions -> n == 1/2]
```

Assumptions on squeezing parameters and coherent amplitudes, which the domain leaves complex, simplify them.
```wl
QuantumInterferometer[{"S"[xi]}]["Amplitude", {0} -> {2}, Assumptions -> xi > 0]
bs["Amplitude", "Coherent"[{a, 0}] -> {1, 0}, Assumptions -> a > 0]
```

Whole distributions take assumptions the same three ways.
```wl
bs["Probabilities", {1, 1}, Assumptions -> th == Pi/4]
Assuming[th == Pi/4, bs["Probabilities", {1, 1}]] === bs["Probabilities", {1, 1}, Assumptions -> {th == Pi/4}]
QuantumInterferometer[{"S"[xi]}]["Probabilities", {0}, "Cutoff" -> 2, Assumptions -> xi > 0]
bs["Probabilities", {1, 1}, Assumptions -> th == I]
```

## 5. Coherent light (passive)

A coherent state stays coherent: |a> -> |W.a>.
```wl
bs["Coherent"[{a, b}]]
```

Coherent input, Fock output.
```wl
bs["Amplitude", "Coherent"[{a, 0}] -> {1, 0}]
```

Fock input, coherent output.
```wl
bs["Amplitude", {1, 0} -> "Coherent"[{c, d}]]
```

Coherent output from coherent input: the overlap probability is a Gaussian in the distance between the amplitudes.
```wl
bs["Probability", "Coherent"[{a, b}] -> "Coherent"[{c, d}]]
Simplify[bs["Probability", "Coherent"[{a, 0}] -> "Coherent"[{a Cos[th], a E^(I ph) Sin[th]}]], Element[{th, ph}, Reals]]
```

Photon distribution up to a cutoff.
```wl
bs["Probabilities", "Coherent"[{a, 0}], "Cutoff" -> 1]
```

A coherent input may also be written mode by mode.
```wl
bs["Amplitude", {"Coherent"[a], "Coherent"[b]} -> {1, 0}] === bs["Amplitude", "Coherent"[{a, b}] -> {1, 0}]
```

Mixed Fock and coherent input, to a Fock or a coherent output.
```wl
bs["Amplitude", {1, "Coherent"[a]} -> {1, 0}]
bs["Amplitude", {1, "Coherent"[a]} -> "Coherent"[{c, d}]]
```

Its distribution starts at the Fock photons, since the coherent part only adds photons.
```wl
Keys[bs["Probabilities", {1, "Coherent"[a]}, "Cutoff" -> 2]]
```

The all-outputs form takes coherent slots and a cutoff, and gives the amplitudes whose squares are that distribution.
```wl
bs[{"Coherent"[a], "Coherent"[b]}, "Cutoff" -> 1]
Union[Values[FullSimplify[(Abs[#]^2 & /@ bs[{1, "Coherent"[a]}, "Cutoff" -> 2]) - bs["Probabilities", {1, "Coherent"[a]}, "Cutoff" -> 2], Element[{th, ph}, Reals]]]]
```

## 6. Squeezing

```wl
sq = QuantumInterferometer[{"S"[xi]}];
sq["Amplitude", {0} -> {2}]
sq["Probability", {0} -> {2}]
sq["Amplitude", {1} -> {3}]
```

All amplitudes or probabilities up to a cutoff on the total photon number. Squeezers change the photon number in pairs, so an odd input keeps only the odd shells.
```wl
sq[{0}, "Cutoff" -> 2]
sq["Probabilities", {0}, "Cutoff" -> 4]
sq["Probabilities", {1}, "Cutoff" -> 3]
QuantumInterferometer[{"S"[x1] -> 1, "S"[x2] -> 2, "BS"[Pi/4, 0]}]["Probabilities", {0, 0}, "Cutoff" -> 2]
```

A QuantumState goes through the Fock-space circuit of SqueezeOperator.
```wl
Normal[sq[FockState[{0}, 3]]["StateVector"]]
Normal[sq[QuantumState[DiagonalMatrix[{1/2, 1/2, 0}], 3]]["DensityMatrix"]]
```

Squeezer then a balanced beam splitter.
```wl
Simplify[QuantumInterferometer[{"S"[r], "BS"[Pi/4, 0]}]["Amplitude", {0, 0} -> {1, 1}], r > 0]
```

Two squeezers interfering: equal squeezing leaves a product of single-mode squeezed vacua, so the coincidence vanishes; opposite squeezing, x2 = -x1, makes a two-mode squeezed vacuum and the coincidence is largest.
```wl
QuantumInterferometer[{"S"[x1] -> 1, "S"[x2] -> 2, "BS"[Pi/4, 0]}]["Amplitude", {0, 0} -> {1, 1}]
```

A squeezer before a general beam splitter; the denominator is Cosh[Abs[xi]] once Cos^2 + Sin^2 = 1 is applied.
```wl
QuantumInterferometer[{"S"[xi] -> 1, "BS"[th, ph]}]["Amplitude", {0, 0} -> {1, 1}]
```

Counts and inverse; the inverse undoes the squeezer.
```wl
{sq["SqueezerCount"], sq["Dagger"]["Elements"]}
FullSimplify[{sq["Dagger"][sq]["Amplitude", {0} -> {0}], sq["Dagger"][sq]["Amplitude", {1} -> {3}]}]
```

The Fock-space circuit with 3 levels; its corner entry is the truncated SqueezeOperator, not an amplitude of the network.
```wl
Simplify[Normal[sq["CircuitOperator", 3]["Matrix"]]]
```

A declared squeezing parameter binds like any other.
```wl
QuantumInterferometer[{"S"[r]}, "Parameters" -> {r}][s0]["Elements"]
```

## 7. Pure loss

Retained transfer matrix.
```wl
l = QuantumInterferometer[{"Loss"[eta]}];
l["TransferMatrix"]
bs[QuantumInterferometer[{"Loss"[eta] -> 2}, 2]]["TransferMatrix"]
```

Count law: a bare binomial. Phases do not change it.
```wl
l["Probability", {n} -> {k}]
QuantumInterferometer[{"PS"[al], "Loss"[eta]}]["Probability", {n} -> {k}]
l["Probability", {3} -> {k}]
l["Probability", {1} -> {0}]
l["Probabilities", {2}]
QuantumInterferometer[{"Loss"[Sin[x]^2]}]["Probability", {2} -> {1}]
```

Independent attenuators: rates on one mode multiply; an untouched mode conserves its photons.
```wl
QuantumInterferometer[{"Loss"[eta] -> 1, "Loss"[mu] -> 1, "Loss"[nu] -> 2}, 3]["Probability", {n1, n2, n3} -> {k1, k2, k3}]
```

Endpoints: a numeric rate or an assumption gives the endpoint law; a symbolic formula reaches it as a limit.
```wl
{QuantumInterferometer[{"Loss"[0]}]["Probability", {n} -> {k}], l["Probability", {n} -> {k}, Assumptions -> eta == 1]}
Assuming[eta == 1, l["Probability", {n} -> {k}]]
Table[Limit[l["Probability", {2} -> {k}], eta -> 1], {k, 0, 2}]
```

Loss with mode mixing, concrete photon numbers.
```wl
lb = QuantumInterferometer[{"BS"[th], "Loss"[eta] -> 1}];
lb["Probability", {1, 0} -> {0, 0}]
lb["Probabilities", {1, 0}]
lb["TransferMatrix"]
FullSimplify[Total[lb["Probabilities", {2, 1}]]]
```

Loss between two beam splitters: in each probability the two phase terms are complex conjugates of each other, so it is visibly real, and the probabilities total 1 on the physical domain.
```wl
l2 = QuantumInterferometer[{"BS"[th, ph], "Loss"[eta] -> 1, "BS"[th2, ph2]}];
l2["Probabilities", {1, 0}]
FullSimplify[Total[l2["Probabilities", {1, 0}]], Element[{th, ph, th2, ph2}, Reals] && 0 <= eta <= 1]
```

Coherent light stays pure through loss: the amplitude shrinks by Sqrt[eta].
```wl
l["Coherent"[{a}]]
lb["Coherent"[{a, b}]]
l["Probabilities", "Coherent"[{a}], "Cutoff" -> 1]
l["Probability", "Coherent"[{a}] -> {2}]
lb["Probability", "Coherent"[{a, 0}] -> {1, 0}]
l["Probability", "Coherent"[{a}] -> "Coherent"[{b}]]
```

Because the output is pure, a coherent input through loss also has amplitudes; the mode-by-mode form works here too.
```wl
l["Amplitude", "Coherent"[{a}] -> {1}]
l["Probability", {"Coherent"[a]} -> {2}]
```

Finite states: coherences scale by Sqrt[eta].
```wl
Normal[l[QuantumState[{1, I}/Sqrt[2], 2]]["DensityMatrix"]]
Simplify[Normal[lb[FockState[{1, 0}, 2]]["DensityMatrix"]]]
```

The channel is a native QuantumChannel; QF applies it, so its output is not simplified by QI.
```wl
Normal[l["Channel", 2][QuantumState[{0, 1}, 2]]["DensityMatrix"]]
Normal[l["CircuitOperator", 2][QuantumState[{0, 1}, 2]]["DensityMatrix"]]
```

Binding eta -> 1 restores a lossless, unitary network; binding a number gives exact counts.
```wl
lq = QuantumInterferometer[{"BS"[th], "Loss"[eta] -> 2}, "Parameters" -> {th, eta}, "DropIdentities" -> False];
{lq["LossCount"], lq[Pi/4, 1]["Unitary"]}
QuantumInterferometer[{"Loss"[eta]}, "Parameters" -> {eta}][1/2]["Probability", {2} -> {1}]
```

The diagram labels symbolic loss.
```wl
Head[l["Diagram", "ShowParameters" -> True]]
```

Assumptions outside the domain, or contradicting each other, fail.
```wl
l["Probability", {n} -> {k}, Assumptions -> eta > 1]
l["Probability", {n} -> {k}, Assumptions -> {eta == 0, eta == 1}]
```

## 8. Not supported symbolically

Symbolic photon numbers work only in a passive "Amplitude" or "Probability" query, and in a loss "Probability" without mode mixing. Entries that are not polynomials with rational coefficients in photon-number variables are rejected.
```wl
bs["Amplitude", {Sqrt[n], 0} -> {t, 0}]
bs["Amplitude", {n + I, 0} -> {t, 0}]
bs["Amplitude", {n + Pi, 0} -> {t, 0}]
bs["Amplitude", {Abs[n], 0} -> {t, 0}]
bs["Amplitude", {n, 1/2} -> {t, 0}]
bs[{n1, n2}]
bs["Probabilities", {n1, n2}]
bs["Amplitude", "Coherent"[{a, 0}] -> {t1, t2}]
bs["Amplitude", {n, 0} -> "Coherent"[{c, d}]]
bs["Amplitude", {nx, Infinity} -> {tx, 0}]
sq["Amplitude", {n} -> {n}]
l["Probabilities", {n}]
l["Probability", "Coherent"[{a}] -> {k}]
lb["Probability", {n, 0} -> {k, 0}]
QuantumInterferometer[{"BS"[th] -> {1, 2}, "Loss"[eta] -> 3}, 3]["Probability", {1, 0, n} -> {0, 1, k}]
QuantumInterferometer[{"BS"[th], "Loss"[eta] -> 1}]["Probability", {1, "Coherent"[a]} -> {1, 0}]
l[{1}]
```

Automatic cutoffs need numbers.
```wl
bs["Probabilities", "Coherent"[{a, 0}]]
bs[{"Coherent"[a], "Coherent"[b]}]
sq["Probabilities", {0}]
```

Active and lossy networks have no unitary; some properties need numbers.
```wl
sq["Unitary"]
sq["TransferMatrix"]
sq["Amplitude", "Coherent"[{a}] -> {0}]
sq["BlochMessiah"]
sq["GBSForm"]
l["Unitary"]
l["Dagger"]
l["Amplitude", {n} -> {k}]
QuantumInterferometer[{"S"[xi], "Loss"[eta]}]
QuantumInterferometer[{{Cos[th], -Sin[th]}, {Sin[th], Cos[th]}}]
bs["MatrixPlot"]
QuantumInterferometer[bs]
QuantumInterferometer[sq]
sq["Coherent"[{a}]]
QuantumInterferometer[{"Loss"[eta], "BS"[th]}]["BlochMessiah"]
QuantumInterferometer[{"Loss"[{eta}]}]
sq["MatrixPlot"]
l["MatrixPlot"]
bs["CircuitOperator"]
```

Sampling needs numbers, lossy, squeezed or not; a squeezed network samples vacuum input only.
```wl
bs["Samples", {1, 0}, 1]
l["Samples", {1}, 2]
bs["Samples", "Coherent"[{a, 0}], 2]
sq["Samples", {0}, 2]
QuantumInterferometer[{"S"[0.3]}]["Samples", {1}, 2]
```

## Gaps seen while running this list

Substitute photon numbers before machine values for symbolic angles: a numeric WignerD outside its index range is Indeterminate. Inexact angles given at construction are safe.
```wl
bs["Amplitude", {n1, n2} -> {t1, t2}] /. {th -> 0.3, ph -> 0.7} /. {n1 -> 1, n2 -> 0, t1 -> 2, t2 -> 0}
bs["Amplitude", {n1, n2} -> {t1, t2}] /. {n1 -> 1, n2 -> 0, t1 -> 2, t2 -> 0} /. {th -> 0.3, ph -> 0.7}
QuantumInterferometer[{"BS"[0.3, 0.7]}]["Amplitude", {n1, n2} -> {t1, t2}] /. {n1 -> 1, n2 -> 0, t1 -> 2, t2 -> 0}
```

The all-outputs form `qi[pattern]`, with Fock photons or a coherent slot, takes no `Assumptions` and ignores an enclosing `Assuming`; `"Amplitude"` and `"Probabilities"` honour both.

Simplification the kernel leaves to the caller:
- `qi["Channel", d]` and the `"CircuitOperator"` of a lossy network apply QF's QuantumChannel, which knows nothing of the physical domain, so their density matrices keep `Conjugate`, as in `Sqrt[eta] Conjugate[Sqrt[eta]]`, and for a passive network `Conjugate[Cos[th]]` and `Conjugate[ph]`. `qi[state]` gives the simplified form for every network and state.
```wl
Normal[bs["Channel", 2][QuantumState[{0, 1, 0, 0}, {2, 2}]]["DensityMatrix"]]
```
- An exact unitary decomposes into angles that are exact but left as `Root` objects and `Arg` of them; `FullSimplify` or `ToRadicals[RootReduce[...]]` reduces them, for `FourierMatrix[3]` to `Pi/4`, `ArcTan[Sqrt[2]]` and `Pi/4`.
```wl
qf["Elements"]
FullSimplify[qf["Elements"]]
```
- `Refine` does not apply trigonometric identities, so a mixed state through a beam splitter keeps `(Cos[th]^2 + Sin[th]^2)/2`, and a squeezer before a general beam splitter keeps `Cos[th]^2 Cosh[Abs[xi]] + Cosh[Abs[xi]] Sin[th]^2`.
- Coherent photon distributions keep `Abs[a Cos[th]]^2` rather than `Abs[a]^2 Cos[th]^2`; `Refine` does not factor the modulus.
- A squeezing parameter written in polar form keeps `Arg[E^(I f)]`, since `Arg` reduces f modulo 2 Pi:
```wl
QuantumInterferometer[{"S"[r Exp[I f]]}]["Amplitude", {0} -> {2}, Assumptions -> r > 0 && Element[f, Reals]]
```
