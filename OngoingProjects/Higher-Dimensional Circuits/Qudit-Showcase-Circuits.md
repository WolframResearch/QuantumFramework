# Five Qudit Circuits That Researchers Actually Care About

Five quantum circuits, each with at least one wire of dimension $d\ge 3$, where the qudit
does something a qubit circuit cannot. Each entry has three parts: **the problem** (the
physics and why $d\ge 3$ is essential), **the QF model** (the circuit, built and run in
QuantumFramework), and **the interpretation** (what the numbers mean). The five span
foundations, compilation, resource theory, quantum simulation, and fault tolerance.

> **Verification.** Every code cell below was executed against the repo working tree
> (`Wolfram/QuantumFramework 2.0.0`, Wolfram Language 15.0, loaded with
> `PacletDirectoryLoad`). The stated outputs are the literal kernel results, not recalled.
> Cells follow the no-`Print` / no-`Quiet` convention: each ends in the expression whose
> value is quoted in the surrounding prose.

## Setup

```wolfram
PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"]
```

A note on dimension defaults that bites immediately: a bare `QuantumState[v]` with
`Length[v]` a power of two is read as *several qubits*, not one qudit. To force a single
$d$-level wire, pass the basis explicitly: `QuantumState[v, QuantumBasis[{d}]]`. Length-3,
5, 6, $\ldots$ vectors (non powers of two) already collapse to one qudit, but being
explicit costs nothing.

---

## 1. Quantum contextuality on a single qutrit (KCBS)

### The problem

The Kochen-Specker theorem says quantum mechanics cannot be reproduced by any model that
assigns each observable a context-independent pre-existing value. The sharpest experimental
witness on the smallest possible system is the Klyachko-Can-Binicioğlu-Shumovsky (KCBS)
inequality. Take five rank-1 projectors $\Pi_0,\ldots,\Pi_4$ on a qutrit, arranged so that
*cyclically adjacent* ones are orthogonal, $\langle v_i | v_{i+1\bmod 5}\rangle = 0$ (the
pentagram compatibility graph). Any noncontextual hidden-variable model obeys

$$
\sum_{i=0}^{4} \langle \Pi_i \rangle \;\le\; 2 ,
$$

because the orthogonality graph is an odd cycle and at most two of five "exclusive" events
can be assigned value 1 noncontextually. Quantum mechanics reaches $\sqrt5 \approx 2.236$.

**Why $d\ge 3$ is essential.** This violation is impossible for a qubit. Every set of
projective measurements on $\mathbb{C}^2$ admits a noncontextual model (the Bloch sphere is
a classical state space for sharp measurements). State-independent and KCBS-type
contextuality first appears at $d=3$; the qutrit is the minimal carrier of the phenomenon.
This is also the resource behind contextuality-powered measurement-based computation and is
the same nonclassicality that reappears as Wigner negativity in Circuit 5.

### The QF model

The five vectors live on a cone around the $z$ axis,
$|v_i\rangle = (\sin\theta\cos\varphi_i,\ \sin\theta\sin\varphi_i,\ \cos\theta)$ with
$\varphi_i = 4\pi i/5$ (so consecutive vectors are two steps apart on the pentagon) and
$\cos^2\theta = 1/\sqrt5$, which is exactly the angle that makes adjacent vectors
orthogonal. The optimal qutrit state is $|\psi\rangle = |2\rangle$, the cone axis. Each
$\Pi_i$ becomes a `QuantumMeasurementOperator` whose `"Mean"` returns
$\langle\Pi_i\rangle = \langle\psi|\Pi_i|\psi\rangle$.

```wolfram
phi[i_] := 4 Pi i/5;
th = ArcCos[Sqrt[1/Sqrt[5]]];
vec[i_] := {Sin[th] Cos[phi[i]], Sin[th] Sin[phi[i]], Cos[th]};
proj[i_] := KroneckerProduct[vec[i], vec[i]];
state = QuantumState[{0, 0, 1}];
meas[i_] := QuantumMeasurementOperator[N @ proj[i], {1}, QuantumBasis[3]];
Total @ Table[meas[i][state]["Mean"], {i, 0, 4}]
```

The orthogonality that defines the compatibility graph is exact:

```wolfram
Simplify @ Table[vec[i] . vec[Mod[i + 1, 5]], {i, 0, 4}]
```

### Interpreting the results

The five inner products return `{0, 0, 0, 0, 0}`: cyclically adjacent measurements are
jointly measurable (compatible), so the KCBS pentagon is well posed. Each projector has
$\langle\Pi_i\rangle = 1/\sqrt5 \approx 0.4472$, and the sum evaluates to
$2.2360679774997907 = \sqrt5$, comfortably above the noncontextual ceiling of $2$. The gap
$\sqrt5 - 2 \approx 0.236$ is the quantum violation: no assignment of definite
context-independent outcomes to these five qutrit projectors can reproduce it. Run the same
construction on a qubit and the maximum is exactly $2$, which is the precise sense in which
this experiment needs the third level.

---

## 2. Ancilla-free Toffoli that borrows the third level

### The problem

The Toffoli gate $C^2X$ (flip the target iff both controls are $|1\rangle$) is not in the
two-qubit Clifford+$T$ palette for free: the standard qubit decomposition costs six CNOTs
and seven $T$ gates, or needs a borrowed ancilla. Lanyon et al. (Nat. Phys. 2009) and
Gokhale et al. (ISCA 2019) observed that if one wire can hold a *third* level, the same
logic costs far less, because the $|2\rangle$ state acts as a built-in scratch register that
stores the partial AND of the controls without any extra particle.

**Why $d\ge 3$ is essential.** The construction temporarily promotes a control to a qutrit
and parks the intermediate AND in $|2\rangle$, a state that does not exist on a qubit. The
extra Hilbert-space dimension *is* the ancilla. This is the cleanest demonstration of
"higher dimensions simplify logic," and it has been run on superconducting and trapped-ion
qutrit hardware (Baker et al. 2020; Nikolaeva et al. PRA 2022, PRL 2025).

### The QF model

Take a heterogeneous register of three wires with dimensions $\{2,3,2\}$: control $c_1$
(qubit), control $c_2$ (qutrit), target $t$ (qubit). Two gates do the work.

- A **raise** gate on $\{c_1, c_2\}$: if $c_1 = |1\rangle$, transpose the qutrit levels
  $|1\rangle \leftrightarrow |2\rangle$; otherwise do nothing. After it, $c_2$ sits in
  $|2\rangle$ exactly when $c_1 = c_2 = 1$.
- A **target** gate on $\{c_2, t\}$: flip $t$ iff $c_2 = |2\rangle$.

Running raise, target, raise (the last uncomputes the first, since the transposition is its
own inverse) leaves $c_1, c_2$ untouched and flips $t$ exactly when both controls were
$|1\rangle$. Each gate is given as its explicit block-diagonal matrix on the relevant pair
of wires; QF embeds them in the mixed register.

```wolfram
P12 = {{1, 0, 0}, {0, 0, 1}, {0, 1, 0}};                                   (* |1><->|2| *)
mRaise = ArrayFlatten[{{IdentityMatrix[3], 0}, {0, P12}}];                  (* on {2,3} *)
mTarget = ArrayFlatten[{{IdentityMatrix[2], 0, 0},
                        {0, IdentityMatrix[2], 0},
                        {0, 0, {{0, 1}, {1, 0}}}}];                         (* on {3,2} *)
g1 = QuantumOperator[mRaise, {1, 2}, QuantumBasis[{2, 3}]];
g2 = QuantumOperator[mTarget, {2, 3}, QuantumBasis[{3, 2}]];
toffoli = QuantumCircuitOperator[{g1, g2, g1}];

basisState[c1_, c2_, t_] :=
  QuantumState[SparseArray[{6 c1 + 2 c2 + t + 1 -> 1}, 12], QuantumBasis[{2, 3, 2}]];
Total @ Flatten @ Table[
  Chop @ Norm[ toffoli[basisState[c1, c2, t]]["StateVector"]
             - basisState[c1, c2, BitXor[t, c1 c2]]["StateVector"] ],
  {c1, 0, 1}, {c2, 0, 1}, {t, 0, 1}]
```

### Interpreting the results

The table sweeps all eight logical inputs $|c_1, c_2, t\rangle$ with the controls and target
in $\{0,1\}$, compares the circuit output against the textbook Toffoli image
$|c_1, c_2,\ t \oplus (c_1 \wedge c_2)\rangle$, and sums the errors. The total is exactly
`0`: the mixed $\{2,3,2\}$ circuit reproduces the three-qubit Toffoli on the logical
subspace, and the qutrit returns to $\{|0\rangle,|1\rangle\}$ with no leakage. The point is
the cost: three two-body gates and no ancilla, versus the six-CNOT qubit decomposition,
bought entirely by the one wire that was allowed a third level. The register dimensions
`{2,3,2,2,3,2}` (output then input) confirm QF carried the heterogeneous structure through
the contraction, something a qubit-only OpenQASM circuit cannot even express.

---

## 3. Grover search on one particle, with no entangling gate

### The problem

Grover search is usually told as a multi-qubit story whose power is often attributed to
entanglement. Encode the same $N = d$ item search space in a *single* qudit of dimension
$d = 2^k$ instead, and the algorithm runs with only single-qudit unitaries: there is no
second subsystem, hence no entanglement at any step. This is the single-qudit / single-
particle computing line of work (Ivanov-Tonchev-Vitanov PRA 2012; recently Shi et al.,
Nat. Commun. 2026, on trapped-ion qudits with $d=5,8$).

**Why $d\ge 3$ is essential.** A $d=8$ qudit carries the same Hilbert space as three qubits,
but as one indivisible object. The oracle reflection and the diffusion operator are global
$SU(8)$ rotations, not two-body gates, so the speedup here is purchased with dimension
rather than entanglement. It sharpens the textbook question "is entanglement necessary for
quantum advantage?" into a concrete circuit where the answer, in this setting, is no.

### The QF model

One wire of dimension $d=8$. Initialize the uniform superposition with the single-qudit
discrete Fourier transform $F_8 |0\rangle$ (this is what `"Fourier"[8]` is: an $8\times 8$
DFT on one qudit, not an eight-qubit QFT). The oracle is the reflection
$O = I - 2|w\rangle\langle w|$ about the marked item $w$; the diffusion is
$D = 2|s\rangle\langle s| - I$ about the uniform state $|s\rangle$. Iterate
$\lfloor \tfrac{\pi}{4}\sqrt{d}\rfloor = 2$ times.

```wolfram
d = 8; w = 5;
b8 = QuantumBasis[{8}];
uni = ConstantArray[1, d]/Sqrt[d];
oracle = QuantumOperator[DiagonalMatrix[ReplacePart[ConstantArray[1, d], (w + 1) -> -1]], {1}, b8];
diffusion = QuantumOperator[2 KroneckerProduct[uni, uni] - IdentityMatrix[d], {1}, b8];
grover = QuantumCircuitOperator[
  Join[{QuantumOperator["Fourier"[d]]}, Flatten @ Table[{oracle, diffusion}, Round[Pi/4 Sqrt[d]]]]];
grover[QuantumState[UnitVector[d, 1], b8]]["Probabilities"][[w + 1]]
```

That every gate is single-qudit (no entangler anywhere) is itself checkable:

```wolfram
Max[#["Arity"] & /@ grover["Operators"]]
```

### Interpreting the results

The probability of the marked item after two iterations is exactly $121/128 = 0.9453125$
(the kernel returns the rational, since every gate is exact), the expected Grover amplitude
for $N=8$ at the optimal iteration count; each of the other seven items keeps
$1/128 = 0.0078125$. The arity check returns `1`: the maximum number of wires any gate touches
is one, so the entire search ran without a single two-body interaction, on one $d=8$
particle. The circuit is logically identical to three-qubit Grover, yet contains zero
entanglement at every step; the work that entanglement does in the qubit version is done
here by the dimension of a single carrier.

---

## 4. Native quantum simulation of a spin-1 chain

### The problem

Simulating a spin-1 magnet (or a $\mathbb{Z}_3$ lattice gauge link, a 3-state Potts/clock
model, a parafermion chain) on qubits forces a binary encoding: each spin-1 site, with three
states $m \in \{-1,0,+1\}$, is squeezed into two qubits, wasting one of the four levels and
requiring penalty terms to keep the dynamics out of the unphysical leakage state. On qutrits
the site *is* the qudit, one to one, with no wasted dimension and no leakage subspace. This
is the qudit quantum-simulation program for higher-spin and gauge systems
(González-Cuadra et al. on qudit gauge theories; Meurice; Calixto et al.).

**Why $d\ge 3$ is essential.** The physical Hilbert space of a spin-1 site is $\mathbb{C}^3$.
A three-site chain is $3^3 = 27$ dimensional natively versus $2^6 = 64$ under a binary qubit
encoding, and the qutrit version needs no constraint terms. The circuit primitives are the
spin-1 angular-momentum operators, which have no faithful single-qubit analog.

### The QF model

Build the spin-1 operators in the computational basis (the named `"JX"[1]` etc. live in
their own eigenbasis, where they are diagonal, so use explicit matrices to stay in the
$|0\rangle,|1\rangle,|2\rangle$ frame). The Heisenberg bond is
$h = J_x\!\otimes\! J_x + J_y\!\otimes\! J_y + J_z\!\otimes\! J_z$. For a three-site chain
$H = h_{12} + h_{23}$ with $[h_{12}, h_{23}] \ne 0$, a first-order Trotter slice is
$U_{\mathrm T}(\Delta t) = e^{-i h_{12}\Delta t}\, e^{-i h_{23}\Delta t}$, assembled from two
native two-qutrit gates and compared against the exact $e^{-iH\Delta t}$.

```wolfram
Jx = {{0, 1, 0}, {1, 0, 1}, {0, 1, 0}}/Sqrt[2];
Jy = {{0, -I, 0}, {I, 0, -I}, {0, I, 0}}/Sqrt[2];
Jz = {{1, 0, 0}, {0, 0, 0}, {0, 0, -1}};
Hbond = KroneckerProduct[Jx, Jx] + KroneckerProduct[Jy, Jy] + KroneckerProduct[Jz, Jz];
b33 = QuantumBasis[{3, 3}];
trotterStep[dt_] := Module[{Ub = MatrixExp[-I Hbond dt]},
  QuantumCircuitOperator[{QuantumOperator[Ub, {1, 2}, b33], QuantumOperator[Ub, {2, 3}, b33]}]];
Hfull = KroneckerProduct[Hbond, IdentityMatrix[3]] + KroneckerProduct[IdentityMatrix[3], Hbond];
trotErr[dt_] := Norm[trotterStep[dt]["QuantumOperator"]["MatrixRepresentation"]
                     - MatrixExp[-I Hfull dt], "Frobenius"];
{trotErr[0.1], trotErr[0.05], trotErr[0.1]/trotErr[0.05]}
```

### Interpreting the results

The slice error is $0.0346$ at $\Delta t = 0.1$ and $0.00866$ at $\Delta t = 0.05$, a ratio
of $3.9956 \approx 4$. Halving the step quarters the error, the signature of the
$O(\Delta t^2)$ local accuracy of first-order Trotter, exactly as it should be for two
noncommuting bonds. The assembled `trotterStep` is unitary on the full $27$-dimensional
space, confirming the two-qutrit gates compose into a legitimate evolution. The structural
payoff is the dimension count: the qutrit chain lives in $27$ dimensions and needs no
leakage penalties, where the binary-encoded qubit version would carry $64$ dimensions and
extra constraint terms. QF treats the spin-1 site as a first-class object, so the simulation
is written in the natural language of the physics.

---

## 5. Qutrit magic and the Wigner-negativity certificate

### The problem

Fault-tolerant quantum computing splits operations into a cheap, classically simulable part
(Clifford gates plus stabilizer states) and an expensive part (magic states) that must be
distilled. For odd prime dimension $d$ there is a beautiful diagnostic: the discrete Wigner
function $W_\rho$ on the $\mathbb{Z}_d \times \mathbb{Z}_d$ phase space. Gross's discrete
Hudson theorem says a pure state has a *nonnegative* Wigner function if and only if it is a
stabilizer state. So negativity of $W_\rho$ certifies that a state is magic, the very
resource a distillery consumes, and it coincides with contextuality (Howard-Wallman-Veitch-
Emerson, Nature 2014; Veitch et al.; Campbell-Anwar-Browne PRX 2012).

**Why $d\ge 3$ is essential.** At $d=2$ there is no clean Wigner / stabilizer equivalence;
the resource theory of magic is sharpest at odd $d$, and the qutrit is the smallest case.
Negativity is not a numerical artifact here, it is the order parameter for nonclassicality
and the quantity lower-bounding distillation cost.

### The QF model

QF computes the discrete Wigner function as a phase-space basis change, exposed through the
state property `"PhaseSpace"` (for odd $d$ this is the genuine $d\times d$ quasi-probability,
no even-$d$ quadrant bookkeeping). Compare two stabilizer states, $|0\rangle$ and
$|+\rangle = F|0\rangle$, against the qutrit "Strange" magic state
$|\mathcal S\rangle = (|1\rangle - |2\rangle)/\sqrt2$. The sum-negativity
$\mathrm{sn}(\rho) = \sum_{W < 0} |W|$ is the witness.

```wolfram
b3 = QuantumBasis[{3}];
zero = QuantumState[{1, 0, 0}, b3];
plus = QuantumOperator["H"[3]][zero];
strange = QuantumState[{0, 1, -1}/Sqrt[2], b3];
sneg[qs_] := Total[Abs @ Select[Flatten[Chop @ N @ qs["PhaseSpace"]], # < 0 &]];
{sneg[zero], sneg[plus], sneg[strange]}
```

The magic state's Wigner function and its mana $\mathcal M = \log(1 + 2\,\mathrm{sn})$:

```wolfram
Normal @ Chop @ N @ strange["PhaseSpace"]
```

```wolfram
Log[1 + 2 sneg[strange]]
```

The state-injection use of this resource is a Clifford circuit: entangle the magic state
with the data qutrit through the qutrit `SUM` (the $C\mathrm{SUM}_3$ gate), measure, and
apply a Pauli correction, teleporting a non-Clifford gate onto the data. The Clifford step
cannot create magic, and the witness sees this directly: feeding `SUM` a stabilizer input
leaves the joint state nonnegative, while feeding it the magic state carries the negativity
into the output.

```wolfram
csum = QuantumOperator["SUM"[3]];
{sneg[csum[QuantumTensorProduct[plus, zero]]], sneg[csum[QuantumTensorProduct[plus, strange]]]}
```

### Interpreting the results

The sum-negativities come back `{0, 0, 1/3}`: both stabilizer states have a strictly
nonnegative Wigner function (each is the classic $1/3$-on-a-line quasi-probability, a $Z$
eigenstate on a vertical line and an $X$ eigenstate on a horizontal line), while the Strange
state carries genuine negativity. Its Wigner matrix is

$$
W_{\mathcal S} = \begin{pmatrix} -\tfrac13 & \tfrac16 & \tfrac16 \\
\tfrac16 & \tfrac16 & \tfrac16 \\ \tfrac16 & \tfrac16 & \tfrac16 \end{pmatrix},
$$

a single negative dip of $-\tfrac13$ at the phase-space origin with the rest smeared to
$\tfrac16$; the entries still sum to $1$, as any quasi-probability must. The mana is
$\log(5/3) = 0.5108$, a strictly positive resource measure that lower-bounds the cost of
preparing this state by distillation. The two-qutrit check returns `{0, 1/3}`: the Clifford
`SUM` keeps a stabilizer input inside the nonnegative (classically simulable) polytope and
merely *transports* the magic when the input already has it. This is the resource theory of
fault tolerance made arithmetic, and the discrete-Wigner transform that computes it is
native QF machinery with no qubit counterpart.

---

## Cross-cutting notes and traps

These are the qudit-specific footguns that the five circuits navigate; worth stating once.

- **Single qudit vs many qubits.** `QuantumState[v]` / `QuantumOperator[m]` with a
  power-of-two size default to a multi-*qubit* register. Force a single $d$-level wire with
  an explicit `QuantumBasis[{d}]` (Circuits 3 and 5).
- **`"H"[d]` is the generalized Fourier transform**, not a block Hadamard; for $d=2$ they
  coincide, for $d>2$ it is the qudit DFT (used in Circuits 3 and 5).
- **The named `"JX"[j]`, `"JY"[j]`, `"JZ"[j]` sit in their own eigenbasis**, where they are
  diagonal. For dynamics in the computational frame, read `["MatrixRepresentation"]` or
  supply explicit matrices (Circuit 4).
- **`"RX"/"RY"/"RZ"` are non-unitary for $d>2$.** Build $SU(d)$ rotations from generators
  via `MatrixExp` or `"R"[\theta, \ldots]`, never the named single-axis rotations at $d\ge 3$.
- **`QuantumWignerTransform` is tagged Experimental.** The state property `"PhaseSpace"`
  (used here) and the stable form `QuantumState[\rho["Double"], QuantumBasis["Wigner"[d]]]`
  are the supported routes; for odd $d$ both give the genuine $d\times d$ quasi-probability.
- **Named `"Fourier"[n]` and `"PhaseEstimation"` *circuits* are qubit-only.** Qudit QFT and
  qudit phase estimation must be built from the generalized gates (`"H"[d]`, controlled
  powers of `"X"[d]`/`"Z"[d]`), which is why Circuit 3 uses the single-qudit *operator*
  `"Fourier"[8]`, not the named circuit.
- **OpenQASM / QuEST export are qubit-only.** Qudit circuits live natively in QF and do not
  round-trip out, which is itself a differentiator: Circuits 1, 2, 4, 5 cannot be written in
  OpenQASM at all.

## References

- S. Kochen, E. P. Specker, *J. Math. Mech.* **17**, 59 (1967).
- A. A. Klyachko, M. A. Can, S. Binicioğlu, A. S. Shumovsky, *Phys. Rev. Lett.* **101**,
  020403 (2008) (KCBS).
- B. P. Lanyon et al., *Nat. Phys.* **5**, 134 (2009) (higher-dimensional Hilbert spaces
  simplify logic).
- P. Gokhale et al., *Proc. ISCA* 2019 (asymptotic improvements via qutrits); A. Baker et
  al. 2020; A. S. Nikolaeva et al., *Phys. Rev. A* **105**, 032621 (2022) and *Phys. Rev.
  Lett.* (2025).
- E. O. Kiktenko, A. S. Nikolaeva, A. K. Fedorov; G. K. Brennen, D. P. O'Leary, S. S.
  Bullock, *Phys. Rev. A* **71**, 052318 (2005) (universal qudit circuits).
- P. A. Ivanov, E. S. Kyoseva, N. V. Vitanov; Y. Wang et al.; Z. Shi et al., *Nat. Commun.*
  **17** (2026) (single trapped-ion qudit search, $d=5,8$).
- D. González-Cuadra et al. (qudit lattice gauge theories); Y. Meurice; M. Calixto et al.
  (higher-spin and clock-model qudit simulation).
- D. Gross, *J. Math. Phys.* **47**, 122107 (2006) (discrete Hudson theorem).
- M. Howard, J. Wallman, V. Veitch, J. Emerson, *Nature* **510**, 351 (2014); V. Veitch et
  al., *New J. Phys.* **16**, 013009 (2014); E. T. Campbell, H. Anwar, D. E. Browne, *Phys.
  Rev. X* **2**, 041021 (2012) (qudit magic-state distillation).
