# The `"Grover"` named circuit is an iteration operator, not the full algorithm

**Verdict: CORRECT-USAGE-FOUND.** The named Grover circuit in QuantumFramework 2.0.0 is
*not* broken. The reported symptom (the all-zeros state being amplified regardless of which
input the oracle marks) is the expected result of using the operator the wrong way: applying
a single Grover reflection to $|0\dots0\rangle$ instead of running the search loop on the
uniform superposition. A secondary, genuine defect lives in the published example notebook,
which uses dead pre-2.0 list-form syntax; that is documented at the end.

Anchor: HEAD `f9dc1cdc` (matches the audit `last-synced.md`). All claims below were checked
live against the working tree via `wolframscript` (`PacletDirectoryLoad`).

## Kernel-bug assessment: none

**There is no kernel bug.** The Grover implementation in
`QuantumFramework/Kernel/QuantumCircuitOperator/NamedCircuits.m` is mathematically exact. When
the operator is driven correctly (uniform superposition, ancilla $|-\rangle$ for the Boolean
variant, $k$ iterations), the measured success probability reproduces the analytic Grover curve

$$
P(k) \;=\; \sin^2\!\big((2k+1)\,\arcsin\sqrt{M/N}\big)
$$

to machine precision across every configuration tested:

| case | variant | $n$ | $M$ | answer at $k_{\mathrm{opt}}$ | $\max_k \lvert P_{\mathrm{measured}} - P(k)\rvert$ |
|---|---|---|---|---|---|
| marked $101$ | `GroverPhase` | 3 | 1 | $101$ | $3.3\times10^{-16}$ |
| marked $41$ | `GroverPhase` | 6 | 1 | $101001$ | $0$ |
| two solutions $1100,\,0011$ | `GroverPhase` | 4 | 2 | solution subspace | $1.1\times10^{-16}$ |
| marked $101$, ancilla $|-\rangle$ | `Grover` (Boolean) | 3 | 1 | $101$ | $1.1\times10^{-16}$ |
| explicit `varSpec` | `GroverPhase` | 3 | 1 | $101$ | (argmax correct) |

If the kernel had a defect in the oracle, the diffusion sign, the reflection axis, or the
multiply-controlled gate, these curves would diverge from the analytic formula; they do not. The
diffusion-via-ancilla design of the Boolean variant (the same $|-\rangle$ ancilla carries both
the oracle and the diffusion phase kickback) is self-consistent: the index register stays pure
across all iterations and traces out exactly the ideal Grover trajectory.

The two real problems are (1) a **usage error** in the original report (one iteration on
$|0\dots0\rangle$ instead of the search loop on $|+\rangle^{\otimes n}$) and (2) a genuine
defect in the **published example notebook**, which is documentation, not kernel code. Both are
detailed below. No change to `Kernel/` is warranted.

## What `"Grover"[bf]` actually is

`QuantumCircuitOperator["Grover"[formula]]` builds a single **Grover iteration operator**

$$
G \;=\; D \cdot O ,
$$

where $O$ is the oracle and $D$ is the diffusion (Grover amplification) operator. It is one
turn of the Grover rotation, not the whole algorithm. The kernel makes this explicit
(`QuantumFramework/Kernel/QuantumCircuitOperator/NamedCircuits.m`):

- Lines 137-156: `"Grover"[formula]` builds a `BooleanOracle` for `formula`, then wraps it as
  the Grover operator with a `"NOT"` gate on the ancilla.
- Lines 121-134: the operator form `"Grover"[op, gate]` composes `GroverDiffusion[...] @ op`,
  i.e. diffusion after the oracle. One iteration, no state preparation, no repetition.

For the 3-variable example `bf = BooleanFunction[2^5, 3]` (true only on input
$5 = 101$), the resulting circuit has 14 gates: 1 oracle gate plus a 13-gate diffusion block
($H^{\otimes 3}, X^{\otimes 3}$, a multiply-controlled NOT, $X^{\otimes 3}, H^{\otimes 3}$),
on 4 qubits (3 index + 1 ancilla). There is no Hadamard layer creating $|+\rangle^{\otimes n}$
and no ancilla preparation.

## Why $|0\dots0\rangle$ gets "amplified"

The diffusion $D$ is a reflection about the uniform-superposition axis. Apply $O$ then $D$ to
$|0000\rangle$: the bit-flip oracle leaves $|0000\rangle$ untouched (its input $000$ is not the
solution and the ancilla stays $|0\rangle$), and a single diffusion reflection of $|0\dots0\rangle$
piles most of the weight back onto $|0\dots0\rangle$. The $49/64 \approx 0.77$ on $|0000\rangle$
that was observed is an artifact of one reflection, with no connection to the marked item. The
same thing happens for every oracle, which is exactly the "always all-zeros" pattern reported.

## Correct usage (verified minimal working examples)

The algorithm is: prepare the uniform superposition, apply $G$ the optimal number of times,
measure. The optimal count for $N = 2^n$ items with $M$ solutions is

$$
k_{\mathrm{opt}} \;=\; \mathrm{round}\!\left(\frac{\pi}{4\,\arcsin\sqrt{M/N}} - \frac{1}{2}\right),
$$

and the success probability after $k$ iterations is

$$
P(k) \;=\; \sin^2\!\big((2k+1)\,\arcsin\sqrt{M/N}\big).
$$

### Phase-oracle variant `"GroverPhase"` (cleanest: no ancilla)

```wolfram
topBits[st_] := (First @ Keys @ ReverseSort @ st["ProbabilityAssociation"]) /. QuditName[b_, ___] :> b;

bf = BooleanFunction[2^5, 3];                    (* marked input = 101 *)
g  = QuantumCircuitOperator["GroverPhase"[bf]];  (* phase oracle marks by (-1)^f(x), no ancilla *)
psi0 = QuantumState["+++"];                       (* uniform superposition; NOT QuantumState["000"] *)
kopt = 2;                                          (* round(pi/4 * Sqrt[8] - 1/2) for N=8, M=1 *)
final = Nest[g[#] &, psi0, kopt];
topBits[final]                                     (* -> {1, 0, 1} = 101  (the marked input) *)
```

Iterating $k = 0 \ldots 6$ reproduces the textbook oscillation:

| $k$ | 0 | 1 | 2 | 3 | 4 | 5 | 6 |
|---|---|---|---|---|---|---|---|
| measured $P$ | 0.125 | 0.781 | 0.945 | 0.330 | 0.012 | 0.548 | 1.000 |
| formula $P(k)$ | 0.125 | 0.781 | 0.945 | 0.330 | 0.012 | 0.548 | 1.000 |

Agreement is exact to $3.3\times10^{-16}$ (floating-point round-off). The argmax is $101$ at
every $k$ where $P$ is near its peak; over-rotation past $k_{\mathrm{opt}}$ swings it back toward
$|000\rangle$ at $k = 4$, as it must.

### Boolean-oracle variant `"Grover"` (needs the ancilla in $|-\rangle$)

The `"Grover"` (Boolean) variant uses a bit-flip oracle
$|x\rangle|q\rangle \to |x\rangle|q \oplus f(x)\rangle$. For this to act as a phase oracle, the
ancilla (the last qubit, here qubit 4) must be prepared in $|-\rangle$ so the flip becomes a
$(-1)^{f(x)}$ phase kickback:

```wolfram
gb   = QuantumCircuitOperator["Grover"[bf]];
psi0 = QuantumTensorProduct[QuantumState["+++"], QuantumState["-"]];   (* index |+++>, ancilla |-> *)
final = Nest[gb[#] &, psi0, 2];
topBits[QuantumPartialTrace[final, {4}]]                                (* -> {1, 0, 1} = 101 *)
```

Ancilla initialization matters quantitatively:

| ancilla | $P(\text{index} = 101)$ at $k=2$ | index-register purity |
|---|---|---|
| $|-\rangle$ | 0.945 (= ideal Grover) | 1.000 (stays pure, ancilla decouples) |
| $|0\rangle$ | 0.535 (degraded) | 0.508 (ancilla entangles, interference spoiled) |

With $|0\rangle$ the oracle records $f(x)$ into the ancilla, entangling it with the index
register; the marked item is still the argmax for this small case but the success probability
is badly degraded. Use $|-\rangle$.

### A second worked case (4 variables)

```wolfram
bf = a && ! b && c && ! d;                         (* unique solution 1010 *)
g  = QuantumCircuitOperator["GroverPhase"[bf]];     (* 4 qubits *)
final = Nest[g[#] &, QuantumState["++++"], 3];      (* kopt = 3 for N=16, M=1 *)
topBits[final]                                      (* -> {1, 0, 1, 0}; success 0.961 *)
```

## The `MemberQ[$QuantumCircuitOperatorNames, "Grover"]` red herring

`$QuantumCircuitOperatorNames` is a `PackageScope` symbol
(`Wolfram`QuantumFramework`PackageScope`$QuantumCircuitOperatorNames`), not a public export.
Typed bare in `Global`` it resolves to an undefined symbol (`ValueQ` is `False`), so
`MemberQ[undefinedSymbol, "Grover"]` returns `False`. The real catalog does contain `"Grover"`
plus seven other Grover-family names:
`{GroverDiffusion, GroverDiffusion0, GroverPhaseDiffusion, GroverPhaseDiffusion0, Grover,
GroverPhase, Grover0, GroverPhase0}`. This was not evidence of a missing circuit.

## Secondary finding: the published example notebook is broken under 2.0.0

`Notebooks/Examples Landing Page/Published examples/Example-Grovers-search-algorithm-1-0-0-examples.nb`
uses pre-2.0 list-form syntax throughout, which now misparses silently (the general v2.0
list-form regression, audit `mistakes.md` 2.1):

- **State preparation.** `QuantumState[{"UniformSuperposition", 3}]` builds a single-qubit
  garbage state (dimensions `{2}`, amplitudes the symbols `UniformSuperposition` and `3`), not
  the 3-qubit uniform superposition. `QuantumState[{"Register", 3, i}]` builds a 3-element
  garbage vector (dimensions `{3}`). Call forms are required:
  `QuantumState["UniformSuperposition"[3]]` or `QuantumState["+++"]`; `QuantumState["Register"[3, i]]`.
- **Circuits.** `QuantumCircuitOperator[{"GroverPhase", bf}]` is parsed as a two-gate circuit
  (gate 1 = `GroverPhase` with its *default* 3-variable formula; gate 2 = `bf` coerced into an
  operator, with the `BooleanFunction` leaking in as a raw amplitude). It returns the wrong
  answer ($110$ instead of $101$ on the test problem). The same applies to `{"Grover", bf}`,
  `{"BooleanOracle", bf}`, and `{"PhaseOracle", bf}`. Correct forms:
  `QuantumCircuitOperator["GroverPhase"[bf]]`, etc.

The notebook should be migrated to the call-form syntax so its cells run correctly under 2.0.0.
This is independent of the Grover-operator semantics above.

## One-line summary

`"Grover"[bf]` is the Grover iteration operator $G = D\cdot O$; run it as
`Nest[g[#] &, QuantumState["+...+"], kopt]` on the uniform superposition (ancilla $|-\rangle$
for the Boolean variant). Applying it once to $|0\dots0\rangle$ is the usage error behind the
"always all-zeros" report.
