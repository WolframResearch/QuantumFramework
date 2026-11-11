```wl
PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"]; Needs["Wolfram`QuantumFramework`SecondQuantization`"];
```

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### QAOA on MaxCut: the exact landscape and a certified optimum

Color the $n$ vertices of a graph with two colors so that as many edges as possible join different colors; take the triangle, the square and the complete graph $K_4$. One layer of the quantum approximate optimization algorithm (QAOA) prepares a state that depends on two angles $\gamma$ and $\theta$ and scores it by the expected number of cut edges. How close does the best pair of angles come to the true maximum cut?

With $Z_i$ and $X_i$ the Pauli operators on vertex $i$, the cost Hamiltonian, the one-layer state and its landscape are

$$
H_C=\sum_{(i,j)\in\mathcal E}\frac{1-Z_iZ_j}{2},\qquad |\psi(\gamma,\theta)\rangle=\prod_i e^{-i\theta X_i/2}\prod_{(i,j)\in\mathcal E}e^{-i\gamma Z_iZ_j/2}\,|+\rangle^{\otimes n},\qquad E(\gamma,\theta)=\langle\psi|H_C|\psi\rangle,
$$

with $\mathcal E$ the edge set; the largest diagonal entry of $H_C$ is the maximum cut.

QF builds $E(\gamma,\theta)$ in closed form from each graph's edge list, `Maximize` finds its exact maximum (dot), and each title compares that maximum, with its minimal polynomial when it is irrational, with the maximum cut:

```wl
qaoaLandscape[graph_] := With[{n = VertexCount[graph], e = List @@@ EdgeList[graph]},
   With[{hc = Total[(QuantumOperator[StringRepeat["I", n]] - QuantumOperator["ZZ", #])/2 & /@ e],
     ψ = QuantumCircuitOperator[{"H" -> Range[n], Splice[("R"[γ, "ZZ"] -> #) & /@ e], "RX"[θ] -> Range[n]}][]},
    {FullSimplify[ψ["Dagger"][hc[ψ]]["Scalar"], {γ, θ} ∈ Reals], Max[Normal[Diagonal[hc["Matrix"]]]]}]];
qaoaPanel[graph_] := With[{l = qaoaLandscape[graph]}, With[{opt = Maximize[First[l], {γ, θ}]},
    DensityPlot[First[l], {γ, 0, 2 π}, {θ, 0, π}, PlotPoints -> 60, PlotLegends -> Automatic, FrameLabel -> {"γ", "θ"},
     Epilog -> {PointSize[Large], Point[Mod[N[{γ, θ} /. Last[opt]], {2 π, π}]]},
     PlotLabel -> Column[{Row[{Graph[graph, ImageSize -> 40], "   max cut = ", Last[l]}],
        Row[{"max E = ", If[First[opt] ∈ Rationals, First[opt],
           Row[{N[First[opt], 4], ", root of ", TraditionalForm[MinimalPolynomial[First[opt], x] == 0]}]]}]}, Alignment -> Center],
     ImageSize -> 300]]];
Grid[{qaoaPanel /@ {CycleGraph[3], CycleGraph[4], CompleteGraph[4]}}]
```

For the triangle the best angles reach the maximum cut. For the square and for $K_4$ the exact maximum of the landscape falls short of it, so no choice of angles makes the expected cut of one layer equal the maximum cut; because the maximum is exact, not the best point a search happened to visit, the shortfall is certified.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### SU(2) double cover from a controlled rotation

A rotation by $2\pi$ returns every Bloch vector to itself but multiplies the spinor by $-1$, so the rotated qubit alone cannot reveal it. An ancilla in superposition decides whether a target qubit, in an arbitrary pure state, is rotated about $x$ by the angle $\alpha$. How often does the ancilla read 1 as a function of $\alpha$?

With the target $|\psi\rangle=\cos\frac\beta2|0\rangle+e^{i\phi}\sin\frac\beta2|1\rangle$ and $R_X(\alpha)=e^{-i\alpha X/2}$, the ancilla starts in $|0\rangle$ and the circuit is

$$
(H\otimes 1)\left(|0\rangle\langle0|\otimes 1+|1\rangle\langle1|\otimes R_X(\alpha)\right)(H\otimes 1)\;|0\rangle|\psi\rangle,
$$

followed by a measurement of the ancilla.

QF applies the circuit with $\alpha$, $\beta$ and $\phi$ as symbols; the title is the simplified probability of reading 1:

```wl
ancillaProbs = FullSimplify[QuantumCircuitOperator[{"H" -> 1, "C"["RX"[α] -> 2, {1}], "H" -> 1, {1}}][
     QuantumTensorProduct[QuantumState["0"], QuantumState[{Cos[β/2], Exp[I ϕ] Sin[β/2]}]]]["ProbabilitiesList"], {α, β, ϕ} ∈ Reals];
Plot[Last[ancillaProbs], {α, 0, 4 π}, Frame -> True,
 FrameLabel -> {"α", "P(ancilla reads 1)"}, FrameTicks -> {{Automatic, None}, {Range[0, 4 π, π], None}},
 PlotLabel -> Row[{"P(1) = ", TraditionalForm[Last[ancillaProbs]]}], ImageSize -> 450]
```

The title has no trace of $\beta$ or $\phi$: every target state gives the same fringe, because the controlled rotation kicks its phase back onto the ancilla. At $\alpha=2\pi$, where every Bloch vector is back where it started, the $-1$ of the spinor flips the ancilla and it reads 1 with certainty; it returns to 0 only at $4\pi$.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Quantum Fisher information: n for product probes, n² for GHZ

A phase $\phi$ is written on $n$ qubits by $R_Z(\phi)$ on each, and any unbiased estimate of $\phi$ from $\nu$ repetitions has variance at least $1/(\nu F_Q)$. How does the quantum Fisher information $F_Q$ grow with $n$ for the product probe $|+\rangle^{\otimes n}$ and for the GHZ probe $(|0\cdots0\rangle+|1\cdots1\rangle)/\sqrt2$?

For a pure probe the phase family and its Fisher information, four times the Fubini-Study metric, are

$$
|\psi_\phi\rangle=\bigotimes_{k=1}^n e^{-i\phi Z_k/2}\,|\psi_0\rangle,\qquad F_Q=4\left(\langle\partial_\phi\psi_\phi|\partial_\phi\psi_\phi\rangle-\left|\langle\psi_\phi|\partial_\phi\psi_\phi\rangle\right|^2\right).
$$

QF's `FubiniStudyMetricTensor` differentiates the state with $\phi$ kept as a symbol, for $n=1$ to $6$:

```wl
Needs["Wolfram`QuantumFramework`QuantumOptimization`"];
qfiPhase[probe_] := 4 FubiniStudyMetricTensor[
     QuantumCircuitOperator[{"RZ"[ϕ] -> Range[probe["Qudits"]]}][probe], "Matrix", "Parameters" -> {ϕ}][[1, 1]];
ListPlot[Transpose[Table[{{m, qfiPhase[QuantumState[StringRepeat["+", m]]]}, {m, qfiPhase[QuantumState["GHZ"[m]]]}}, {m, 6}]],
 PlotLegends -> {"product", "GHZ"}, PlotMarkers -> Automatic, Joined -> True, Frame -> True, FrameLabel -> {"n", "\!\(\*SubscriptBox[\(F\), \(Q\)]\)"},
 PlotRange -> {0, 40}, ImageSize -> 400]
```

The product probe gives $F_Q=n$, the standard quantum limit, and the GHZ probe $n^2$, the Heisenberg limit; at $n=1$ the two probes are the same state.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### A coherent error on the three-qubit code, undone by measurement records

Encode a qubit in three, rotate the middle one by the angle $2\theta$ about $x$ (a coherent error, not a Pauli flip), measure the two parities, and flip back the qubit the parities point to. QF keeps each measurement outcome as a record qubit, so each correction is a gate controlled by the records. Does the cycle restore the encoded state for every $\theta$?

With the input $|\psi\rangle=\cos\frac\beta2|0\rangle+e^{i\phi}\sin\frac\beta2|1\rangle$, the code word, the error and the measured parity bits are

$$
|\bar\psi\rangle=\cos\tfrac\beta2|000\rangle+e^{i\phi}\sin\tfrac\beta2|111\rangle,\qquad e^{-i\theta X_2},\qquad m_{12}=\frac{1-Z_1Z_2}{2},\quad m_{23}=\frac{1-Z_2Z_3}{2},
$$

each measured by its two eigenspace projectors, which leave superpositions inside an eigenspace intact. The records $(m_{12},m_{23})=(1,0),(1,1),(0,1)$ trigger $X_1$, $X_2$, $X_3$.

QF runs the cycle with $\theta$, $\beta$ and $\phi$ as symbols; the plot shows the record probabilities and the fidelity of the three code qubits with $|\bar\psi\rangle$:

```wl
parityBit = (QuantumOperator["II"] - QuantumOperator["ZZ"])/2;
qecOut = QuantumCircuitOperator[{"CNOT" -> {1, 2}, "CNOT" -> {1, 3}, "RX"[2 θ] -> 2,
     "Lueders"[parityBit] -> {1, 2}, "Lueders"[parityBit] -> {2, 3},
     "C"["X" -> 1, {0}, {-1}], "C"["X" -> 2, {0, -1}, {}], "C"["X" -> 3, {-1}, {0}]}][
   QuantumTensorProduct[QuantumState[{Cos[β/2], Exp[I ϕ] Sin[β/2]}], QuantumState["00"]]];
recordProbs = FullSimplify[Diagonal[Normal[QuantumPartialTrace[qecOut, {3, 4, 5}]["DensityMatrix"]]], {θ, β, ϕ} ∈ Reals];
codeFidelity = FullSimplify[QuantumSimilarity[QuantumPartialTrace[qecOut, {1, 2}],
    QuantumState[{Cos[β/2], 0, 0, 0, 0, 0, 0, Exp[I ϕ] Sin[β/2]}], "Fidelity"], {θ, β, ϕ} ∈ Reals];
With[{curves = {recordProbs[[1]], recordProbs[[2]] + recordProbs[[3]], recordProbs[[4]]}},
 Plot[Evaluate[Append[curves, codeFidelity]], {θ, 0, π}, PlotRange -> {-0.05, 1.05},
  PlotStyle -> {Automatic, Automatic, Automatic, Dashed}, Frame -> True, FrameLabel -> {"θ", "probability, fidelity"},
  FrameTicks -> {{Automatic, None}, {Range[0, π, π/4], None}},
  PlotLegends -> MapThread[Row[{#1, ": ", TraditionalForm[#2]}] &, {{"records (0,0)", "one odd parity", "records (1,1)", "fidelity"}, Append[curves, codeFidelity]}],
  PlotLabel -> Row[{"fidelity of code qubits = ", codeFidelity, " for all β, ϕ"}], ImageSize -> 450]]
```

The records read $(0,0)$ with probability $\cos^2\theta$ and $(1,1)$ with probability $\sin^2\theta$, and never a single odd parity: the measurement turns the rotation into either no error or a full $X_2$ flip. The record-controlled gates undo both branches, so the code qubits are back in $|\bar\psi\rangle$ with fidelity 1 for every $\theta$ and every input; the error is left only in the records.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Teleportation with the classical channel as record wires

Alice holds an unknown qubit and one half of a Bell pair; Bob holds the other half. She measures her two qubits in the Bell basis, and Bob corrects his qubit according to her two outcomes. QF keeps the outcomes as record qubits, so Bob's corrections are gates controlled by the records. Does Bob end up with Alice's state, for every input?

With the input $|\psi\rangle=\cos\frac\beta2|0\rangle+e^{i\phi}\sin\frac\beta2|1\rangle$ on qubit 1, whose Bloch vector is $(\sin\beta\cos\phi,\sin\beta\sin\phi,\cos\beta)$, and $|00\rangle$ on qubits 2 and 3, the protocol is

$$
\underbrace{\mathrm{CNOT}_{23}\,H_2}_{\text{Bell pair}},\qquad \underbrace{H_1\,\mathrm{CNOT}_{12}}_{\text{Bell basis}},\qquad \text{measure qubits 1, 2}\to m_1,m_2,\qquad \text{Bob applies } Z^{m_1}X^{m_2}.
$$

QF runs it with $\beta$ and $\phi$ as symbols, with both record-controlled corrections, with one of them, and with none. The titles give Bob's Bloch vector in closed form; the points compare its $x$, $y$ and $z$ components with the input's for six inputs:

```wl
bobBloch[corrections_] := FullSimplify[QuantumPartialTrace[Replace[
      QuantumCircuitOperator[Join[{"H" -> 2, "CNOT" -> {2, 3}, "CNOT" -> {1, 2}, "H" -> 1, {1}, {2}}, corrections]][
       QuantumTensorProduct[QuantumState[{Cos[β/2], Exp[I ϕ] Sin[β/2]}], QuantumState["00"]]], m_QuantumMeasurement :> m["State"]],
     {1, 2, 3, 4}]["BlochVector"], {β, ϕ} ∈ Reals];
teleportInputs = Thread[{β, ϕ} -> #] & /@ {{π/6, 0}, {π/3, π/2}, {π/2, π}, {2 π/3, 3 π/2}, {5 π/6, π/4}, {π/2, π/4}};
teleportPanel[label_, bob_] := ListPlot[Transpose[{{Sin[β] Cos[ϕ], Sin[β] Sin[ϕ], Cos[β]}, bob}] /. teleportInputs // Transpose,
   PlotLegends -> {"x", "y", "z"}, PlotMarkers -> Automatic,
   PlotRange -> {{-1.1, 1.1}, {-1.1, 1.1}}, AspectRatio -> 1, Frame -> True, FrameLabel -> {"input Bloch component", "Bob's Bloch component"},
   PlotLabel -> Column[{label, TraditionalForm[bob]}, Alignment -> Center], ImageSize -> 300];
teleportPanels = MapThread[teleportPanel[#1, bobBloch[#2]] &, {{"X and Z", "X only", "Z only", "no correction"},
    {{"CNOT" -> {-1, 3}, "CZ" -> {0, 3}}, {"CNOT" -> {-1, 3}}, {"CZ" -> {0, 3}}, {}}}];
Legended[Grid[Partition[First /@ teleportPanels, 2]], Last[First[teleportPanels]]]
```

With both corrections Bob's Bloch vector equals the input's for every $\beta$ and $\phi$, so every point lies on the diagonal. One correction alone keeps only the component its Pauli leaves unflipped, and with neither Bob holds the maximally mixed state. Each record tells Bob which Pauli his qubit picked up; without a record, the Bloch components that its Pauli flips average to zero.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Grover search: the success probability, iteration by iteration

A Boolean oracle accepts $M$ of the $N=2^n$ bit strings. Grover's search starts from the uniform superposition and repeats one iteration, a phase flip of the accepted strings followed by a reflection about the uniform state. How does the probability of finding an accepted string change with the number of iterations $k$?

With $f$ the oracle and $|s\rangle=|+\rangle^{\otimes n}$, one iteration and the state after $k$ of them are

$$
G=\left(2|s\rangle\langle s|-1\right)O,\qquad O|x\rangle=(-1)^{f(x)}|x\rangle,\qquad |\psi_k\rangle=G^k|s\rangle,\qquad P_k=\sum_{f(x)=1}|\langle x|\psi_k\rangle|^2 .
$$

QF builds $G$ from the Boolean function and applies it to the state again and again, up to about $\pi\sqrt N$ iterations; the accepted strings are a few fixed integers below 16:

```wl
groverSuccess[n_, accepted_] := With[{g = QuantumCircuitOperator["GroverPhase"[BooleanFunction[Total[2^accepted], n]]]},
   Total[#["ProbabilitiesList"][[accepted + 1]]] & /@ NestList[g, QuantumState[StringRepeat["+", n]], Ceiling[π Sqrt[2^n]]]];
groverAccepted = {{5}, {5, 9}, {5, 9, 12, 14}};
groverPanels = Table[ListLinePlot[groverSuccess[n, #] & /@ groverAccepted, DataRange -> {0, Ceiling[π Sqrt[2^n]]}, Mesh -> All,
    PlotRange -> {0, 1.02}, Frame -> True, FrameLabel -> {"iterations k", "success probability"}, PlotLabel -> Row[{"N = ", 2^n}],
    PlotLegends -> (Row[{"M = ", Length[#]}] & /@ groverAccepted), ImageSize -> 300], {n, {4, 6, 8}}];
Legended[Grid[{First /@ groverPanels}], Last[First[groverPanels]]]
```

The success probability oscillates with $k$: it climbs close to 1 and falls again as the iterations rotate the state past the accepted strings. The first peak comes later for a larger search space and earlier for more accepted strings; for $N=16$ and $M=4$ a single iteration finds an accepted string with certainty.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Shor's discrete logarithm on a toy elliptic curve

The points of $y^2=x^3+x+1$ over the integers mod 23, with the point at infinity $\mathcal O$, form a cyclic group of order $r$ under point addition. A public key is $Q=\ell G$ for a generator $G$ and a secret integer $\ell$. Can the interference pattern of Shor's algorithm reveal $\ell$ from $G$ and $Q$ alone?

Two registers of $t=7$ qubits hold $a$ and $b$ in uniform superposition, a third register receives $f(a,b)=aG+bQ$ and is measured, and a Fourier transform on each input register gives the outcome $(c,d)$ with probability

$$
P(c,d)=\sum_{w}\Big|\frac{1}{2^{2t}}\sum_{f(a,b)=w}e^{-2\pi i(ac+bd)/2^t}\Big|^2 .
$$

Since $f$ is constant along $a+\ell b\equiv\text{const} \pmod r$, the peaks sit near $(c,d)\approx(j,\,j\ell)\,2^t/r$, and every shot with $j$ invertible mod $r$ gives a candidate $\ell=\mathrm{round}(rd/2^t)\,j^{-1} \bmod r$, kept only if $\ell G=Q$.

WL adds points by the chord-and-tangent rule, simulates $P(c,d)$ with `Fourier`, and draws 20 shots (dots) for three secret keys; each title gives the candidates that pass the check:

```wl
ecAdd[∞, q_] := q; ecAdd[q_, ∞] := q;
ecAdd[{x1_, y1_}, {x2_, y2_}] /; x1 == x2 && Mod[y1 + y2, 23] == 0 := ∞;
ecAdd[{x1_, y1_}, {x2_, y2_}] := With[{s = Mod[If[x1 == x2, 3 x1^2 + 1, y2 - y1] PowerMod[If[x1 == x2, 2 y1, x2 - x1], -1, 23], 23]},
   With[{x3 = Mod[s^2 - x1 - x2, 23]}, {x3, Mod[s (x1 - x3) - y1, 23]}]];
ecMultiples[pt_] := NestList[ecAdd[#, pt] &, ∞, 127];
ecPoints = Select[Tuples[Range[0, 22], 2], Mod[#[[2]]^2 - #[[1]]^3 - #[[1]] - 1, 23] == 0 &];
ecOrder = Length[ecPoints] + 1;
ecGen = SelectFirst[ecPoints, DuplicateFreeQ[Take[ecMultiples[#], ecOrder]] &];
ecGenMultiples = ecMultiples[ecGen];
shorRecover[pub_, shots_] := Select[DeleteDuplicates[
    With[{j = Mod[Round[ecOrder #[[1]]/128], ecOrder], l = Mod[Round[ecOrder #[[2]]/128], ecOrder]},
       If[CoprimeQ[j, ecOrder], Mod[l PowerMod[j, -1, ecOrder], ecOrder], Nothing]] & /@ shots], ecGenMultiples[[# + 1]] === pub &];
shorPanel[key_] := With[{pub = ecGenMultiples[[key + 1]]},
   With[{fab = Outer[ecAdd, ecGenMultiples, ecMultiples[pub], 1]},
    With[{probs = Total[Function[w, Abs[Fourier[Boole[Map[# === w &, fab, {2}]], FourierParameters -> {0, -1}]]^2] /@
         DeleteDuplicates[Flatten[fab, 1]]]/128^2},
     With[{shots = RandomChoice[Flatten[probs] -> Tuples[Range[0, 127], 2], 20]},
      ListDensityPlot[Transpose[probs], DataRange -> {{0, 127}, {0, 127}}, InterpolationOrder -> 0, PlotRange -> All,
       Epilog -> {PointSize[0.008], Point[shots]}, FrameLabel -> {"c", "d"}, ImageSize -> 400,
       PlotLabel -> Row[{"secret key = ", key, ",   recovered: ", shorRecover[pub, shots]}]]]]]];
SeedRandom[1]; Grid[{shorPanel /@ {3, 7, 20}}]
```

The probability sits on a lattice of isolated peaks whose arrangement changes with $\ell$, and the shots land on those peaks. In every panel at least one of the 20 shots gives a candidate that passes $\ell G=Q$: the peak positions supply $\ell$, and the check only discards wrong candidates.
