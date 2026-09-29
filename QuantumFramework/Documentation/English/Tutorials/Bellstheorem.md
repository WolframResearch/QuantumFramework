---
Template: TechNote
Name: Bellstheorem
Title: Bell's Theorem: CHSH inequality
Context: Wolfram`QuantumFramework`
CellContext: Global`
Paclet: Wolfram/QuantumFramework
URI: Wolfram/QuantumFramework/tutorial/Bellstheorem
RelatedGuides: [WolframQuantumComputationFramework]
---

```wl
Column@EntityValue[
  Entity["Person", "JohnStewartBell::364v4"], {"Image", "BirthDate", 
   "DeathDate", "NotableFacts"}]
```

Bell's inequalities are examined within the Wolfram quantum framework, using the standard language of quantum computation. A straightforward computational approach is employed to comprehend the underlying derivations of those equations, fully derived from standard quantum computation.

## Entangled states, correlations, and measurement

Consider a 2-qubit quantum system that is prepared in the following quantum state

```wl
\[Psi]m = QuantumState["PsiMinus"]
```

Check if it is an entangled state:

```wl
QuantumEntangledQ[\[Psi]m]
```

which is, of course, maximally entangled:

```wl
QuantumEntanglementMonotone[\[Psi]m]
```

An entangled state is a state that cannot be written as product of two states $\rho_{1}\otimes \rho_{2}$. This possibility can be tested using [QuantumEntangledQ](https://resources.wolframcloud.com/PacletRepository/resources/Wolfram/QuantumFramework/ref/QuantumEntangledQ.html) which returns as True or False.

Let’s define a unit vector in 3D (which is parametrized by two angles):

```wl
a[\[Theta]_, \[Phi]_] := {Sin[\[Theta]] Cos[\[Phi]], 
  Sin[\[Theta]] Sin[\[Phi]], Cos[\[Theta]]}
```

Given two unit vectors, let’s find the angle between:

```wl
vecAngle = 
 FullSimplify[VectorAngle[a[\[Theta]1, \[Phi]1], a[\[Theta]2, \[Phi]2]],
   Assumptions -> {\[Theta]1 \[Element] Reals, \[Theta]2 \[Element] 
     Reals, \[Phi]1 \[Element] Reals, \[Phi]2 \[Element] Reals}]
```

Now let us calculate the quantum expectation (i.e. the mean value) of a composite operator as follows: $\langle \psi^{-}|\vec{a}_{1}.\vec{\sigma }\otimes \vec{a}_{2}.\vec{\sigma }|\psi^{-}\rangle $ with $|\psi^{-}\rangle =\frac{1}{\sqrt{2}}(|01\rangle -|10\rangle )$ and $\vec{\sigma }$ the Pauli vector.

Define a Pauli vector:

```wl
\[Sigma][\[Theta]_, \[Phi]_] := 
 QuantumOperator[a[\[Theta], \[Phi]] . Table[PauliMatrix[i], {i, 3}]]
```

Calculate $\langle \psi^{-}|\vec{a}_{1}.\vec{\sigma }\otimes \vec{a}_{2}.\vec{\sigma }|\psi^{-}\rangle $:

```wl
\[Psi]m["Dagger"][
 QuantumTensorProduct[\[Sigma][\[Theta]1, \[Phi]1], \
\[Sigma][\[Theta]2, \[Phi]2]][\[Psi]m]]
```

Note that the result is a scalar, but as an object, we treat it as a quantum state object, with 0 number of qudits (look at the summary box of above result). The actual scalar number can be extracted like this:

```wl
\[CapitalPhi]12 = FullSimplify[%["Number"]]
```

The expectation value is related to the vector angle as follows:

```wl
-Cos[vecAngle] == \[CapitalPhi]12
```

In other words: $\langle \psi^{-}|\vec{a}_{1}.\vec{\sigma }\otimes \vec{a}_{2}.\vec{\sigma }|\psi^{-}\rangle =-cos(\Phi_{12})$ Remember this result. We will use it later within Bell’s inequality. But before jumping to Bell’s inequality, let’s explore how above expectation value can be calculated or obtained experimentally using quantum circuits.

For simplicity, let’s consider the case where qubit-1 is measured in Pauli-X basis, and qubit-2 in another basis which obtained by rotating Pauli-X basis by π/8 around z-axis

Define new basis (Pauli-X which is rotated by π/8 around z-axis) and label it:

```wl
newBasis = 
  QuantumBasis[N@QuantumOperator[{"RZ", \[Pi]/8}][QuantumBasis["X"]], 
   "Label" -> 
    "\!\(\*SubscriptBox[\((\*FractionBox[\(\[Pi]\), \(8\)])\), \
\(xz\)]\)"];
```

Note I added N, just to enforce numeric calculations, which are faster.

Measure in Pauli-X basis on qubit-1, in new basis on qubit-2, when the system is prepared in $\psi^{-}$:

```wl
meas = QuantumMeasurementOperator["X"]@
   QuantumMeasurementOperator[newBasis, {2}]@QuantumState["PsiMinus"];
```

Find the measurement probabilities:

```wl
prob = meas["Probabilities"]
```

Note the first and last results correspond to + eigenvalue (their multiplication) and middle ones -. So for the mean value, we do as follows:

```wl
prob[[1]] + prob[[4]] - (prob[[2]] + prob[[3]])
```

Taking $cos^{-1}$ of above result (note $\langle \psi^{-}|\vec{a}_{1}.\vec{\sigma }\otimes \vec{a}_{2}.\vec{\sigma }|\psi^{-}\rangle =-cos(\Phi_{12})$) one can confirm:

```wl
ArcCos[-%] == \[Pi]/8
```

which is what we expected. Now let’s focus on the circuit version.

The following quantum circuit prepare a 2-qubit quantum system in $\psi^{-}$ first, and then two measurements are performed on 2-qubits (as described previously).

```wl
qc1 = QuantumCircuitOperator[{QuantumCircuitOperator[{"X" -> {1, 2}, 
      "H", "CNOT"}, 
     "Circuit to create \!\(\*SuperscriptBox[\(\[Psi]\), \(-\)]\)"], \
{"RZ", \[Pi]/8} -> 2, "H" -> {1, 2}, {1}, {2}}];
qc1["Diagram"]
```

The above circuit is the one you see in many literature, mostly because that is how you send it to a Quantum Processing Unit (QPU, ie a quantum hardware). However, in our framework, you can define the measurement in any basis, for example:

```wl
qc2 = QuantumCircuitOperator[{QuantumCircuitOperator[{"X" -> {1, 2}, 
      "H", "CNOT"}, 
     "Circuit to create \!\(\*SuperscriptBox[\(\[Psi]\), \(-\)]\)"], 
    QuantumMeasurementOperator["X"], 
    QuantumMeasurementOperator[newBasis, {2}]}];
qc2["Diagram"]
```

Note that those light purple boxes representing measurement are fundamentally different from other boxes representing usual gates. But that is a different story ([refer to this for more details](https://arxiv.org/pdf/0806.0647.pdf)).

The above circuits are equivalent. However, for the sake of communicating with a QPU, some transpilers may have some issues with defining measurements for a customized basis. So let’s focus on first circuit:

```wl
qc1["Diagram"]
```

In above circuit, $\vec{a}_{1} \to \{\theta =\pi /2,\phi =0\}$ (measurement on qubit-1) and $\vec{a}_{2} \to \{\theta =\pi /2,\phi =\pi /8\}$ (measurement on qubit-2). Note the Hadamard (H) boxes before measurements are transforming the measurement to Pauli-X, and $R_{z}$ rotation rotates it by an angle in xy-plane. So the correlation function should be $-cos(\pi /8)$

Find the probabilities from above circuit:

```wl
prob1 = FullSimplify /@ qc1[]["Probabilities"]
```

Calculate the correlation function:

```wl
prob1[
Wolfram`QuantumFramework`QuditName[{0, 0}, "Dual" -> False]] + prob1[
Wolfram`QuantumFramework`QuditName[{1, 1}, 
    "Dual" -> False]] - (prob1[
Wolfram`QuantumFramework`QuditName[{0, 1}, "Dual" -> False]] + prob1[
Wolfram`QuantumFramework`QuditName[{1, 0}, 
      "Dual" -> False]]) // FullSimplify
```

which is, of course, what we expected.

Note that the experimental results that one gets out of a QPU is counters per measurement results; something like this:

```wl
exp = qc1[]
```

Simulate measurement results for 1024 shots:

```wl
counts = Counts[exp["SimulatedMeasurement", 1024]]
```

Calculate the corresponding correlation (using frequency of occurrence for each result):

```wl
(counts[
Wolfram`QuantumFramework`QuditName[{0, 0}, "Dual" -> False]] + 
  counts[
Wolfram`QuantumFramework`QuditName[{1, 1}, "Dual" -> False]] - (counts[
Wolfram`QuantumFramework`QuditName[{0, 1}, "Dual" -> False]] + 
    counts[
Wolfram`QuantumFramework`QuditName[{1, 0}, "Dual" -> False]]))/1024 //
  N
```

which is what we expected (close enough, given number of shots)

```wl
-Cos[\[Pi]/8] // N
```

Now let’s see how those correlations are used in Bell inequality

## Bell’s inequality

There are many good literature on the Bell’s inequality and what it implies. If you are looking for a short paper, with some nice [personal] history, we highly recommend this preprint by [GianCarlo Ghirardi](https://en.wikipedia.org/wiki/Giancarlo_Ghirardi): [John Stuart Bell: recollections of a great scientist and a great man](https://arxiv.org/pdf/1411.1425.pdf).

Let us rephrase the reasoning behind the derivation of Bell’s inequality (following Ghirardi):

“Experimental Perfect Correlations” ∧ “Bell’s Locality” $\Rightarrow$ “Determinism”  “Determinism” ∧ “Bell’s Locality” $\Rightarrow$ “Bell’s Inequality”  “General Quantum Correlations” $\Rightarrow$ ¬ “Bell’s Inequality”  ¬“Bell’s Inequality” $\Rightarrow$ ¬ “Determinism” ∨ ¬“Bell’s Locality”  Summarizing: “Natural Processes” $\Rightarrow$ ¬“Bell’s Locality”, i.e.: *Nature is non-locally causal*.

From locality to determinism. If distant measurements on the singlet state always give opposite results (perfect correlations), and if we insist that what happens at one site cannot depend on what happens at the other (Bell's Locality), then the outcomes must have been fixed in advance, determinism follows as a consequence. If outcomes are predetermined and Bell's Locality holds, then the correlations between measurement results are constrained, they must satisfy Bell's inequality. Quantum mechanics violates the inequality. The correlations predicted by quantum theory (and confirmed by experiment) exceed the bound imposed by Bell's inequality. Something must give. Since the inequality is violated, at least one of its premises fails: either determinism is wrong, or Bell's Locality is wrong (or both). Locality fails regardless. Suppose we try to save things by abandoning determinism. But then we lose the only local explanation for why perfect correlations occur. Without locality, we can't derive determinism, and without determinism, we can't locally account for the perfect correlations we actually observe. So dropping determinism alone doesn't help; locality must fail too.

Well, if anything, let's comment on Bell’s Locality, which implies that measurement probabilities (accessible or not, to observer) at two space-like regions are independent.

Bell’s inequality in the [CHSH form](https://en.wikipedia.org/wiki/CHSH_inequality) can be written as follows: $|E_{\mu }(\vec{a}_{1},\vec{a}_{2})-E_{\mu }(\vec{a}_{1},\vec{a}_{4})+E_{\mu }(\vec{a}_{3},\vec{a}_{2})+E_{\mu }(\vec{a}_{3},\vec{a}_{4})| \le |E_{\mu }(\vec{a}_{1},\vec{a}_{2})-E_{\mu }(\vec{a}_{1},\vec{a}_{4})|+|E_{\mu }(\vec{a}_{3},\vec{a}_{2})+E_{\mu }(\vec{a}_{3},\vec{a}_{4})| \le 2$ with μ the accessible variable of the system that specifies the quantum state (eg the preparation step). Note there is another variable in the derivation of Bell’s inequality, usually denoted b λ, which is inaccessible variable specifying the quantum state (the quantum world/level is obtained by averaging over λ).

Also, $E_{\mu }(\vec{a}_{i},\vec{a}_{j})$ is identified as the quantum expectation value. Assuming that the state is $\psi^{-}$ and using the result for the correlation function ($\langle \psi^{-}|\vec{a}_{1}.\vec{\sigma }\otimes \vec{a}_{2}.\vec{\sigma }|\psi^{-}\rangle =-cos(\Phi_{12})$ with $\Phi_{12}$ the vector angle between $\vec{a}_{1}$ and $\vec{a}_{2}$), the left-side of above inequality can be written as:

```wl
Abs[Cos[VectorAngle[a1, a2]] - Cos[VectorAngle[a1, a4]]] + 
  Abs[Cos[VectorAngle[a3, a2]] + Cos[VectorAngle[a3, a4]]] <= 2
```

For example, let us consider this case: $VectorAngle[a1,a2]=VectorAngle[a3,a2]=VectorAngle[a3,a4]=\phi $ and $VectorAngle[a1,a4]=3\phi $. Plugging this into above equation, and one can plot regions of ϕ such that the Bell’s inequality is violated:

```wl
NumberLinePlot[
 Not[Abs[Cos[\[Phi]] - Cos[3 \[Phi]]] + Abs[2 Cos[\[Phi]]] <= 
   2], {\[Phi], 0, 2 \[Pi]/3}, 
 Ticks -> {{0, \[Pi]/6, \[Pi]/3, \[Pi]/2, 2 \[Pi]/3}}, 
 PlotLabel -> "Regions within the Bell's inequality are violated"]
```

One can obtain those regions (where the Bell’s inequality is violated) explicitly:

```wl
Reduce[Abs[Cos[\[Phi]] - Cos[3 \[Phi]]] + Abs[2 Cos[\[Phi]]] > 2 && 
  0 <= \[Phi] <= 2 \[Pi]/3, \[Phi]]
```

One can also find the maximum violation, which happens at π/4:

```wl
NMaximize[
 Abs[Cos[\[Phi]] - Cos[3 \[Phi]]] + Abs[2 Cos[\[Phi]]], \[Phi]]
```

## Explore Bell’s inequality in a quantum circuit

Given four unit vectors $\vec{a}_{i}$, we will focus on this case: $VectorAngle[a1,a2]=VectorAngle[a3,a2]=VectorAngle[a3,a4]=\phi $ and $VectorAngle[a1,a4]=3\phi $. To find those four correlation functions, we will design four circuit as follows:

```wl
corQC12 = 
  QuantumCircuitOperator[{"X" -> {1, 2}, "H", 
    "CNOT", {"RZ", \[Phi]} -> 2, "H" -> {1, 2}, {1, 2}}];
corQC12["Diagram"]
```

```wl
corQC32 = 
  QuantumCircuitOperator[{"X" -> {1, 2}, "H", 
    "CNOT", {"RZ", \[Phi]}, {"RZ", 2 \[Phi]} -> 2, 
    "H" -> {1, 2}, {1, 2}}];
corQC32["Diagram"]
```

```wl
corQC34 = 
  QuantumCircuitOperator[{"X" -> {1, 2}, "H", 
    "CNOT", {"RZ", 2 \[Phi]} -> 1, {"RZ", 3 \[Phi]} -> 2, 
    "H" -> {1, 2}, {1, 2}}];
corQC34["Diagram"]
```

```wl
corQC14 = 
  QuantumCircuitOperator[{"X" -> {1, 2}, "H", 
    "CNOT", {"RZ", 3 \[Phi]} -> 2, "H" -> {1, 2}, {1, 2}}];
corQC14["Diagram"]
```

Now let’s calculate the corresponding correlations, using above circuits. We will find those correlations analytically, using our quantum framework. Of course, one can find the numerical value of those correlations (theoretically, by some simulation as shown in the previous section, or experimentally using some QPUs).

```wl
$Assumptions = 0 <= \[Phi] <= 2 \[Pi]/3;
```

Measurement probabilities from circuit to find $E_{\mu }(\vec{a}_{1},\vec{a}_{2})$:

```wl
meas12 = FullSimplify /@ corQC12[]["Probabilities"]
```

```wl
cor12 = FullSimplify[(#[
Wolfram`QuantumFramework`QuditName[{0, 0}, "Dual" -> False]] + #[
Wolfram`QuantumFramework`QuditName[{1, 1}, "Dual" -> False]] - #[
Wolfram`QuantumFramework`QuditName[{0, 1}, "Dual" -> False]] - #[
Wolfram`QuantumFramework`QuditName[{1, 0}, "Dual" -> False]]) &@
   meas12]
```

Measurement probabilities from circuit to find $E_{\mu }(\vec{a}_{3},\vec{a}_{2})$:

```wl
meas32 = FullSimplify /@ corQC32[]["Probabilities"]
```

```wl
cor32 = FullSimplify[(#[
Wolfram`QuantumFramework`QuditName[{0, 0}, "Dual" -> False]] + #[
Wolfram`QuantumFramework`QuditName[{1, 1}, "Dual" -> False]] - #[
Wolfram`QuantumFramework`QuditName[{0, 1}, "Dual" -> False]] - #[
Wolfram`QuantumFramework`QuditName[{1, 0}, "Dual" -> False]]) &@
   meas32]
```

Measurement probabilities from circuit to find $E_{\mu }(\vec{a}_{3},\vec{a}_{4})$:

```wl
meas34 = FullSimplify /@ corQC34[]["Probabilities"]
```

```wl
cor34 = FullSimplify[(#[
Wolfram`QuantumFramework`QuditName[{0, 0}, "Dual" -> False]] + #[
Wolfram`QuantumFramework`QuditName[{1, 1}, "Dual" -> False]] - #[
Wolfram`QuantumFramework`QuditName[{0, 1}, "Dual" -> False]] - #[
Wolfram`QuantumFramework`QuditName[{1, 0}, "Dual" -> False]]) &@
   meas34]
```

Measurement probabilities from circuit to find $E_{\mu }(\vec{a}_{1},\vec{a}_{4})$:

```wl
meas14 = FullSimplify /@ corQC14[]["Probabilities"]
```

```wl
cor14 = FullSimplify[(#[
Wolfram`QuantumFramework`QuditName[{0, 0}, "Dual" -> False]] + #[
Wolfram`QuantumFramework`QuditName[{1, 1}, "Dual" -> False]] - #[
Wolfram`QuantumFramework`QuditName[{0, 1}, "Dual" -> False]] - #[
Wolfram`QuantumFramework`QuditName[{1, 0}, "Dual" -> False]]) &@
   meas14]
```

Plug them into Bell’s inequality:

```wl
Abs[cor12 - cor14] + Abs[cor32 + cor34] <= 2
```

Plot the region that Bell’s inequality (above) is violated:

```wl
NumberLinePlot[
 Not[Abs[cor12 - cor14] + Abs[cor32 + cor34] <= 2], {\[Phi], 0, 
  2 \[Pi]/3}, Ticks -> {{0, \[Pi]/6, \[Pi]/3, \[Pi]/2, 2 \[Pi]/3}}, 
 PlotLabel -> "Regions within the Bell's inequality is voilated"]
```

Certainly, the aforementioned circuits may be executed using existing QPUs to obtain empirical correlations. However, such findings hold little value unless we know the details of experimental setting employed, to address problems commonly known as “loopholes” (such as the locality loophole).

## Reference

GianCarlo Ghirardi: *John Stuart Bell: recollections of a great scientist and a great man*, arXiv:1411.1425v1 (2014).

John S. Bell: *The Trieste lecture of John Stewart Bell*, Journal of Physics A: Mathematical and Theoretical 40, 2919-2933 (2007).
