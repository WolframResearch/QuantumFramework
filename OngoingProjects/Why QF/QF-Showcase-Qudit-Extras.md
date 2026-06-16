# QF Showcase: Extended Qudit Examples (held for future use)

These three worked examples were developed and kernel-verified for the main showcase ([QF-Showcase.md](QF-Showcase.md)) and then lifted out of it to keep Claim 2 focused on the *object model* rather than on a tour of qutrit physics. They are preserved here intact, with their original prose and verified code, in case a future version wants a dedicated higher-dimensional section or claim. Each is self-contained: a qutrit ($d=3$) or a single $d=8$ qudit, using `QuantumBasis[{d}]`/`QuantumBasis[3]` throughout.

All cells ran as part of the showcase suite (Wolfram Language 15.0, `Wolfram/QuantumFramework` 2.0.0), document order, one kernel, zero messages. Source circuits with traps and provenance: [Qudit-Showcase-Circuits.md](../Higher-Dimensional%20Circuits/Qudit-Showcase-Circuits.md).

---

## Discrete-Wigner magic and the Clifford/stabilizer boundary

A natural coda to the continuous Wigner function of a cat state: the Wigner function has a discrete sibling, and at odd dimension it turns into a fault-tolerance resource. For a qutrit, QF computes the discrete Wigner function as a phase-space basis change, exposed through the state property `"PhaseSpace"`, a genuine $3\times3$ quasi-probability. Gross's discrete Hudson theorem makes it a magic detector: a pure qutrit state has a nonnegative Wigner function *if and only if* it is a stabilizer state, so negativity certifies *magic*, the resource a fault-tolerant computer must distill. Take the sum-negativity $\sum_{W<0}|W|$ of two stabilizer states, $|0\rangle$ and $|+\rangle = H|0\rangle$, against the "Strange" magic state $(|1\rangle - |2\rangle)/\sqrt2$:
```wolfram
b3 = QuantumBasis[{3}];
zero = QuantumState[{1, 0, 0}, b3];
plus = QuantumOperator["H"[3]][zero];
strange = QuantumState[{0, 1, -1}/Sqrt[2], b3];
sneg[qs_] := Total[Abs @ Select[Flatten[Chop @ N @ qs["PhaseSpace"]], # < 0 &]];
{sneg[zero], sneg[plus], sneg[strange]}
```
The negativities are `{0, 0, 1/3}`: both stabilizer states are strictly nonnegative, while the Strange state carries genuine negativity. Its Wigner matrix is a single dip:
```wolfram
Normal @ Chop @ N @ strange["PhaseSpace"]
```
A lone $-\tfrac13$ sits at the phase-space origin, the rest smeared to $\tfrac16$, the entries still summing to $1$ as any quasi-probability must. The associated magic measure is the mana $\mathcal M = \log(1 + 2\,\mathrm{sn})$:
```wolfram
Log[1 + 2 sneg[strange]]
```
It returns $\log(5/3) \approx 0.511$, a strictly positive resource that lower-bounds the cost of distilling this state. A Clifford gate cannot create that resource; it can only move it. Feed the two-qutrit `SUM` gate (the qutrit analog of a CNOT) a stabilizer input versus the magic state and watch the negativity:
```wolfram
csum = QuantumOperator["SUM"[3]];
{sneg[csum[QuantumTensorProduct[plus, zero]]], sneg[csum[QuantumTensorProduct[plus, strange]]]}
```
The pair is `{0, 1/3}`: the Clifford `SUM` keeps a stabilizer input inside the nonnegative, classically simulable polytope, and merely transports the magic when the input already has it. This is the resource theory of fault tolerance made arithmetic, computed by a discrete-Wigner transform with no qubit counterpart, since at $d=2$ there is no clean Wigner-stabilizer equivalence.

## A qutrit breaks classicality: contextuality

That qutrit negativity is one face of nonclassicality; a sharper one, carrying no quasi-probability bookkeeping at all, is *contextuality*. The Kochen-Specker theorem says no model can assign every observable a context-independent pre-existing value, and the smallest witness lives on a qutrit: the KCBS inequality. It reads five yes/no measurements arranged in a pentagon, each one jointly measurable with its two neighbors; every noncontextual hidden-variable model obeys $\sum_{i=0}^4 \langle\Pi_i\rangle \le 2$, and a qutrit overshoots it.

Each measurement projects onto a ray, and the five rays sit evenly on a cone. Rather than write five direction vectors out by hand, take one seed ray and turn it about the axis by the pentagram step $4\pi/5$, four times over, with the built-in `RotationMatrix`:
```wolfram
th = ArcCos[5^(-1/4)];
rays = Table[RotationMatrix[4 Pi k/5, {0, 0, 1}] . {Sin[th], 0, Cos[th]}, {k, 0, 4}];
Chop[N[rays]] // MatrixForm
```
Five unit vectors share a single height, the cone's $\cos\theta = 5^{-1/4}$ standing in the third column, while the first two coordinates walk around the pentagram: the whole configuration generated from one seed, with no per-ray trigonometry.

The pentagon is well posed only where cyclically adjacent measurements are compatible, which for rank-one projectors means their rays are orthogonal. Pair the list with its own one-step rotation, `RotateLeft`, and read every neighboring overlap at once:
```wolfram
Chop @ N @ MapThread[Dot, {rays, RotateLeft[rays]}]
```
The five inner products come back `{0, 0, 0, 0, 0}`: each ray is orthogonal to its successor, so neighbors are jointly measurable and the pentagon closes.

Now turn each ray into its projective measurement $\Pi_i = |v_i\rangle\langle v_i|$, the outer product of the ray with itself, hand all five to QF as measurement operators, and sum their expectations in the cone-axis state $|2\rangle$:
```wolfram
apex = QuantumState[{0, 0, 1}];
Total[QuantumMeasurementOperator[N @ Outer[Times, #, #], {1}, QuantumBasis[3]][apex]["Mean"] & /@ rays]
```
The sum returns $\sqrt5 \approx 2.236$, above the noncontextual ceiling of $2$ by the gap $\sqrt5 - 2$. This violation is *impossible for a qubit*: every set of projective measurements on $\mathbb{C}^2$ admits a noncontextual model, so KCBS-type contextuality first appears at $d=3$, and the qutrit is its minimal carrier, the same nonclassicality the Wigner negativity above measured in a different currency.

## Search on one qudit, with no entangling gate

Beyond foundations, the extra dimension does computational work. Grover search is usually a multi-qubit story whose power is often pinned on entanglement; encode the same $N = d$ search space in a *single* qudit of dimension $d = 8$ instead, and the algorithm runs with only single-qudit unitaries, no second subsystem, hence no entanglement at any step. Initialize the uniform superposition with the qudit Fourier transform (here `"Fourier"[8]` is an $8\times 8$ DFT on one wire, not an eight-qubit QFT), reflect about the marked item and the uniform state $\lfloor\tfrac\pi4\sqrt8\rfloor = 2$ times, and read the marked probability:
```wolfram
d = 8; w = 5;
b8 = QuantumBasis[{8}];
uni = ConstantArray[1, d]/Sqrt[d];
oracle = QuantumOperator[DiagonalMatrix[ReplacePart[ConstantArray[1, d], (w + 1) -> -1]], {1}, b8];
diffusion = QuantumOperator[2 KroneckerProduct[uni, uni] - IdentityMatrix[d], {1}, b8];
grover = QuantumCircuitOperator[Join[{QuantumOperator["Fourier"[d]]},
   Flatten @ Table[{oracle, diffusion}, Round[Pi/4 Sqrt[d]]]]];
grover[QuantumState[UnitVector[d, 1], b8]]["Probabilities"][[w + 1]]
```
The marked item comes back with probability exactly $121/128 \approx 0.945$, the optimal Grover amplitude for $N = 8$, leaving $1/128$ on each of the other seven. That every gate is single-qudit, with no entangler anywhere, is itself checkable:
```wolfram
Max[#["Arity"] & /@ grover["Operators"]]
```
The maximum arity is `1`: the entire search ran without a single two-body interaction, on one $d = 8$ particle. The circuit is logically identical to three-qubit Grover, yet holds zero entanglement at every step; the work entanglement does in the qubit version is done here by the dimension of one carrier.

---

**Setup.** Evaluate once before the cells:

```wolfram
Needs["Wolfram`QuantumFramework`"]
Needs["Wolfram`QuantumFramework`SecondQuantization`"];
Needs["Wolfram`TensorNetworks`"];
```
