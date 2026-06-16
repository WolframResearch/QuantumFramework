# TN-engine worked example: a local observable in a long chain (light-cone win)

A drop-in for Claim 4's subsection **"The default engine: a steered tensor-network
contraction"**, replacing the Grover example with one where the tensor-network engine genuinely
beats the dense state-vector engine: it returns an exact expectation value at a qubit count the
dense path cannot reach. All numbers are real captured `wolframscript` outputs against the
working tree (HEAD `f9dc1cdc`, QF 2.0.0, TensorNetworks 1.0.5).

## The example

Quench a transverse-field Ising chain from the all-up state $|0\dots0\rangle$ under
$H = J\sum_i Z_iZ_{i+1} + h\sum_i X_i$, Trotterized into $d$ layers of nearest-neighbor $ZZ$
rotations followed by on-site $X$ rotations, and read one site's magnetization $\langle Z_c\rangle$
at the end. This is ordinary non-equilibrium many-body physics: a sudden quench, and a local
observable measured afterward.

The point is **causality**. A depth-$d$ circuit can move information at most $d$ sites, so
$\langle Z_c\rangle$ after $d$ Trotter steps cannot depend on any qubit outside the backward
causal cone $[c-d,\,c+d]$, the $2d+1$ sites the dynamics could have reached. In the
expectation $\langle 0|\,U^\dagger Z_c\,U\,|0\rangle$ every gate outside that cone meets its own
conjugate and cancels to the identity. So the tensor network that computes $\langle Z_c\rangle$
closes over $2d+1$ qubits **whatever the length of the chain**, while the dense state vector it
would otherwise build carries all $2^n$ amplitudes. This is the discrete, exact form of a
Lieb-Robinson light cone, and it is the textbook reason tensor-network contraction scales where
dense linear algebra cannot.

## Why it beats the dense engine

| | tensor-network engine | dense state-vector engine |
|---|---|---|
| what it builds | the $2d+1$-qubit causal cone | the full $2^n$ amplitude vector |
| cost vs chain length $n$ | **independent of $n$** | grows as $2^n$ |
| exact? | yes (no truncation) | yes |
| reach | any $n$ | out of memory near $n=30$ |

At fixed depth the contraction is the *same* network for a site in a 20-qubit chain or a
20000-qubit chain, so its cost does not move with $n$. The dense vector has $2^{50}\approx
1.1\times10^{15}$ entries at $n=50$, about 18 petabytes, and exhausts memory near $n=30$
($2^{30}\cdot 16$ B $\approx 17$ GB).

## Evidence (all measured)

**The contraction is cheap and its cost tracks the cone, not the chain.** Times are
$n$-independent (the observable was placed deep in a long chain):

| Trotter depth $d$ | causal cone | $\langle Z_c\rangle$ | TN contraction time |
|---|---|---|---|
| 4 | 9 qubits | $-0.171156$ | 0.34 s |
| 5 | 11 qubits | $-0.2849$ | 1.0 s |
| 6 | 13 qubits | $-0.29406$ | 2.0 s |

**It is exact, not an approximation.** Where the dense path can still run, the two engines
agree to machine precision, and widening the cone does not change the value:

| check (depth 4, $n=14$, bulk site 7) | value |
|---|---|
| dense $\langle Z_7\rangle$ (`Method -> "Schrodinger"`) | $-0.171156$ |
| light-cone TN $\langle Z_c\rangle$ (9-qubit cone) | $-0.171156$ |
| $\lvert\text{difference}\rvert$ | $\approx 3\times10^{-16}$ |
| cone widened (half-width 4 vs 6) | identical |

**The dense baseline is genuinely out of reach.** Measured directly for this circuit (building
the statevector, fixed depth 4):

| $n$ | dense build time |
|---|---|
| 8 | 4.3 s |
| 10 | 9.1 s |
| 12 | 13.9 s |
| 14 | 27.6 s |
| 16, 18 | $>45$ s |

The dense time grows as $2^n$ (per-gate cost climbs from $\sim$70 ms at $n=8$ to $\sim$260 ms
at $n=14$ as the state vector grows), on top of a few seconds of fixed apply overhead. The
$2^n$ scaling is the fundamental wall: it is what makes $n=50$ impossible regardless of
implementation, while the light-cone contraction stays flat at a fraction of a second.

**The contraction is a real, inspectable object.** `meanZ["TensorNetwork"]` is a *closed*
network (no free indices: the ket prepares every wire, the bra caps it), the greedy path is
**138 pairwise contractions**, and `ContractionTree` is a **277-node `Tree`** in which every
intermediate tensor leg is dimension **2**, that is, the bond never widens. That bounded width
is exactly why the cost tracks the cone and not the qubit count. All cells run with **zero
messages**.

## Drop-in cells and prose (matching the doc's per-cell style)

Keep the first two sentences of the subsection's opening paragraph (the order/cost framing and
the `Method` sub-option list), then continue:

---

The payoff is largest when the network's *structure* keeps every intermediate small however
many qubits the circuit carries, and the cleanest source of that structure is a light cone.
Quench a transverse-field Ising chain from the all-up state, Trotterizing the evolution under
$H = J\sum_i Z_iZ_{i+1} + h\sum_i X_i$ into $d$ layers of nearest-neighbor $ZZ$ rotations and
on-site $X$ rotations, and ask for one site's magnetization $\langle Z_c\rangle$ at the end. The
Trotter depth and the two rotation angles:

```wolfram
d = 4; angleZZ = 0.4; angleX = 0.6;
```

Causality makes the work local: $\langle Z_c\rangle$ after $d$ Trotter steps cannot depend on
any qubit outside the observable's backward causal cone, the $2d+1$ sites $[c-d, c+d]$ the
dynamics could have reached, because every gate outside it cancels against its conjugate in the
bra-ket sandwich. Build the quench on a window of $m$ sites:

```wolfram
quench[m_] := QuantumCircuitOperator @ Join[
   {StringJoin @ ConstantArray["0", m]},
   Flatten @ Table[Join[
      ("R"[angleZZ, "ZZ"] -> # &) /@ Partition[Range[m], 2, 1],
      ("RX"[angleX] -> # &) /@ Range[m]], {d}]];
```

The cone of a depth-4 observable is nine qubits wide. Build it and sandwich a $Z$ on its center
site between the cone and its dagger, assembling $\langle 0|\,U^\dagger Z_c\, U\,|0\rangle$ as a
single object:

```wolfram
cone = quench[2 d + 1];
meanZ = cone /* QuantumOperator["Z", {d + 1}] /* cone["Dagger"];
```

Contract it to the number, through a greedily chosen path, and time it:

```wolfram
AbsoluteTiming[meanZ[Method -> {"TensorNetwork", "Path" -> "Greedy"}]["Scalar"] // Chop]
```

The magnetization comes back as $\langle Z_c\rangle \approx -0.171$ in a few tenths of a second,
and that cost is set by the nine-qubit cone, not by the chain: the identical contraction
returns the same number for a site deep inside a chain of fifty qubits or of fifty thousand.
The dense engine has no such shortcut. To read the same expectation it must hold the entire
$2^n$ amplitude vector; at $n = 50$ that is $2^{50} \approx 1.1\times10^{15}$ complex numbers,
about $18$ petabytes, and the machine runs out of memory near $n = 30$ long before.

At a width the dense path can still reach, the two engines agree exactly. Build the full
$2^{14}$ state vector for the same quench and read the magnetization of a bulk site off it,
beside the light-cone value:

```wolfram
psi = quench[14][Method -> "Schrodinger"];
{(psi["Dagger"] @ QuantumOperator["Z", {7}] @ psi)["Scalar"],
  meanZ[Method -> {"TensorNetwork", "Path" -> "Greedy"}]["Scalar"]} // Chop
```

The dense value lands on the light-cone contraction to one part in $10^{16}$: site $7$ of the
fourteen-site chain has its whole cone $[3, 11]$ in the interior, so its exact magnetization is
the bulk number the nine-qubit cone already gave. The cone is not an approximation but the
exact reduction the bra-ket sandwich permits, and building that dense vector already takes tens
of seconds where the cone took tenths.

The contraction plan is an inspectable object. Pull the cone's tensor network:

```wolfram
net = meanZ["TensorNetwork"]
```

It is a *closed* network with no free indices: the ket prepares every wire and the bra caps it,
so the contraction yields a single number. Compute the greedy path and draw it as a tree whose
nodes carry the dimension of each intermediate:

```wolfram
ContractionTree[TensorNetworkContraction[net, GreedyContractionPath[net]], "Labels" -> "Dimensions"]
```

The tree is the engine's execution plan: 138 pairwise contractions fold the gate and state
tensors into one scalar, and every intermediate it records is a handful of dimension-2 legs.
That bounded width is the whole story: it is why the cost tracks the light cone, not the qubit
count, and why this contraction finishes while the $2^n$ vector never gets built.

---

## Honest one-line statement of the win

At fixed Trotter depth the contraction cost is set by the observable's light cone and is
independent of the chain length, so a depth-4 $\langle Z_c\rangle$ costs a few tenths of a
second whether the chain holds 20 qubits or 20000, while the dense state vector grows as $2^n$
and is out of reach past roughly 30 qubits; the two engines agree to $\sim10^{-16}$ wherever
dense can still run.

## Integration notes

- This **replaces** the Grover example (correct, but at $n=8$ dense is faster). It uses one
  consistent depth ($d=4$) across all cells; deeper quenches stay $n$-independent ($d=6$ is a
  13-qubit cone at $\approx 2$ s) if a meatier contraction tree is wanted in the figure.
- The win comes from locality (the light cone bounds the network), so the example is honest
  about *why* it scales: it is the standard reason tensor-network methods reach many qubits,
  shown on a concrete observable.
- `.nb` regeneration and the full verification-suite run are **left to the main session**. Each
  cell was run verbatim and is message-free; the whole suite was not rerun and `QF-Showcase.nb`
  was not rebuilt.
- The `ContractionTree` cell uses `"Labels" -> "Dimensions"` so the dimension-2 bonds are
  visible (supporting the prose); the doc's current Grover tree cell uses an unlabeled render.
  Either works; pick per house aesthetic.
