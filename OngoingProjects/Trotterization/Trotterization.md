# `TrotterizationStandalone\`` — A Pure-WL Suzuki–Trotter Decomposition

A standalone Wolfram Language package implementing the Suzuki–Trotter product formula at arbitrary even order, independent of `Wolfram\`QuantumFramework\``. This document explains the math, walks through every public function, demonstrates use on Hamiltonians of one to four qubits, verifies that the package agrees with `QuantumFramework`'s built-in `"Trotterization"` named circuit to machine precision, and tests numerically that the error decays at the rate Suzuki's recursion predicts. Everything below is computation-first: every formula has an evaluable example, and every claim is checked.

We start with the math, build up the API one function at a time, work through four example Hamiltonians, verify against `QuantumFramework`, then close with a convergence study.

## Table of contents

1. [Theory primer — the Suzuki–Trotter formula](#1-theory-primer)
2. [Loading the package](#2-loading-the-package)
3. [`TrotterizationCoefficients`](#3-trotterizationcoefficients)
4. [`TrotterizationOrdering`](#4-trotterizationordering)
5. [`TrotterizationMatrices`](#5-trotterizationmatrices)
6. [Examples on physical Hamiltonians](#6-examples-on-physical-hamiltonians)
7. [Agreement with `QuantumFramework`](#7-agreement-with-quantumframework)
8. [Convergence tests](#8-convergence-tests)
9. [Where this leaves us](#9-where-this-leaves-us)

---

## 1. Theory primer

A time-independent Hamiltonian generates the unitary $U(t) = e^{-i t H}$. Splitting $H$ into a sum of "easy" terms, $H = \sum_{i=1}^{l} A_i$, does not give a simple product: $e^{-i t H} \neq \prod_i e^{-i t A_i}$ unless every pair $[A_i, A_j]$ vanishes. The Suzuki–Trotter product formula approximates $U(t)$ by a product of single-term exponentials with carefully chosen coefficients.

**Order 1 (Lie–Trotter).** The crudest split is

$$
S_1(t) = e^{-i t A_l} \cdots e^{-i t A_2} \, e^{-i t A_1},
\qquad
\bigl\| U(t) - S_1(t) \bigr\| = \mathcal{O}(t^2).
$$

**Order 2 (Strang).** Symmetrising the order-1 product yields a one-step error of order $t^3$:

$$
S_2(t) = e^{-i (t/2) A_1} \cdots e^{-i (t/2) A_l} \, e^{-i (t/2) A_l} \cdots e^{-i (t/2) A_1}.
$$

**Higher orders (Suzuki).** Suzuki's recursion lifts $S_2$ to even order $n$ via a five-fold construction with weight $p = 1/(4 - 4^{1/(n-1)})$:

$$
S_n(t) = S_{n-2}(p\,t)\, S_{n-2}(p\,t)\, S_{n-2}\bigl((1-4p)\,t\bigr)\, S_{n-2}(p\,t)\, S_{n-2}(p\,t).
$$

**Trotter steps.** For $r$ "reps" (Trotter steps), apply $S_n(t/r)$ to itself $r$ times. The error of the assembled approximation $S_n(t/r)^r$ scales as $\mathcal{O}\!\left(t^{n+1}/r^{n}\right)$ (see arXiv:[math-ph/0506007](https://arxiv.org/abs/math-ph/0506007)).

The package exposes Suzuki's recursion as three independent functions: a coefficient list (`TrotterizationCoefficients`), a term-index ordering (`TrotterizationOrdering`), and the assembled gate list (`TrotterizationMatrices`).

---

## 2. Loading the package

The package lives next to this document. Load it directly:

```wolfram
Get["Trotterization.wl"];
```

After loading, the symbols `TrotterizationCoefficients`, `TrotterizationOrdering`, and `TrotterizationMatrices` are available in the public context `TrotterizationStandalone\``. There is no dependency on `Wolfram\`QuantumFramework\`` — the package needs only `MatrixExp`, pattern matching, and `Dot`.

```wolfram
?TrotterizationStandalone`*
(* Output: TrotterizationCoefficients, TrotterizationMatrices, TrotterizationOrdering *)
```

---

## 3. `TrotterizationCoefficients`

`TrotterizationCoefficients[l, order, c]` returns the coefficient list of length $|S|$, where $|S|$ is the number of single-term exponentials in the Suzuki product at the given order for a Hamiltonian with `l` terms.

| Argument | Description | Default |
|---|---|---|
| `l` | Number of Hamiltonian terms (positive integer) | required |
| `order` | Suzuki order (positive integer; odd inputs round up to the next even value) | required |
| `c` | Overall prefactor (any number, symbol, or expression) | `1` |

The total weight is conserved at every order: `FullSimplify @ Total[TrotterizationCoefficients[l, order, c]] == c * l`. (At orders 4 and above, the sum is a non-trivial Suzuki identity hidden behind surds — `FullSimplify` is needed to see it.)

**Order 1, two terms.** A Hamiltonian $A_1 + A_2$ at order 1 produces one exponential per term:

```wolfram
TrotterizationCoefficients[2, 1]
(* {1, 1} *)
```

**Order 2, two terms.** Strang splitting halves each coefficient and gives a palindrome:

```wolfram
TrotterizationCoefficients[2, 2]
(* {1/2, 1/2, 1/2, 1/2} *)
```

**Order 4, two terms.** Suzuki's recursion produces 20 exponentials, weighted by $p$ and $1-4p$:

```wolfram
TrotterizationCoefficients[2, 4] // Length
(* 20 *)

TrotterizationCoefficients[2, 4] // Union // FullSimplify
(* {1/(2 - 4 * 2^(1/3)), 1/(8 - 2 * 2^(2/3))} *)

N @ %
(* {-0.328814, 0.207198} *)
```

The two distinct values are the Suzuki weight $p/2 = 1/(8 - 2\cdot 2^{2/3}) \approx 0.207$ (appearing 16 times in the order-4 product for $l=2$) and $(1-4p)/2 = 1/(2 - 4\cdot 2^{1/3}) \approx -0.329$ (appearing 4 times). The negative coefficient is the well-known feature that makes higher-order Suzuki schemes mix forward and backward steps.

**Odd orders are rounded up.** Odd inputs map to the next even order, so `order = 3` is identical to `order = 4`:

```wolfram
TrotterizationCoefficients[2, 3] === TrotterizationCoefficients[2, 4]
(* True *)
```

**Symbolic prefactors propagate.** The `c` argument is treated symbolically; it can be a fraction like `1/reps` or a generic symbol:

```wolfram
TrotterizationCoefficients[3, 2, 1/r]
(* {1/(2 r), 1/(2 r), 1/(2 r), 1/(2 r), 1/(2 r), 1/(2 r)} *)
```

This is exactly how the package threads the per-step weight `time/reps` through Suzuki's recursion when assembling matrices.

---

## 4. `TrotterizationOrdering`

`TrotterizationOrdering[l, order]` returns the sequence of term indices (or operators, if `l` is a list of operators) that the product formula traverses. It mirrors `TrotterizationCoefficients` element-for-element: the $k$-th gate is $\exp(-i\, c_k\, A_{i_k})$ where $c_k$ comes from the coefficient list and $A_{i_k}$ from the ordering list.

```wolfram
TrotterizationOrdering[{a, b}, 1]
(* {a, b} *)

TrotterizationOrdering[{a, b}, 2]
(* {a, b, b, a} *)
```

The order-2 ordering is the palindrome that gives Strang splitting its symmetry.

For higher orders the recursion catenates five copies of the lower-order ordering:

```wolfram
TrotterizationOrdering[{a, b}, 4] // Length
(* 20 *)

TrotterizationOrdering[{a, b, c}, 4] // Length
(* 30 *)
```

In general the length of the ordering at order $n$ is $l \cdot 2 \cdot 5^{(n-2)/2}$. This is what determines the gate count.

---

## 5. `TrotterizationMatrices`

`TrotterizationMatrices[ops, order, reps, time, opts...]` assembles Suzuki's product into a list of $d \times d$ unitaries (or, optionally, into a single product matrix).

| Argument | Description | Default |
|---|---|---|
| `ops` | List of $d \times d$ Hermitian matrices $A_i$ with $H = \sum_i A_i$ | required |
| `order` | Suzuki order | `1` |
| `reps` | Number of Trotter steps | `1` |
| `time` | Total evolution time (numeric or symbolic) | `1` |
| `"Form" -> ...` | Output shape: `"Flat"`, `"Steps"`, or `"Unitary"` | `"Flat"` |

**Output shapes.**

- `"Flat"` (default) — flat list of `reps * Length[TrotterizationOrdering[ops, order]]` matrices, in circuit-application order. The approximate unitary is `Dot @@ Reverse[gates]` (right-most gate applied first).
- `"Steps"` — list of length `reps`; each entry is the gate list of one Trotter step. Useful when you want to inspect or reuse a single step.
- `"Unitary"` — the assembled approximate unitary as one matrix.

**Sign convention.** The package returns gates of the form $\exp(-i\,c_k\,A_{i_k})$. The total approximate unitary equals $\prod_k \exp(-i\,c_k\,A_{i_k})$ (right-to-left), which approximates $\exp(-i\,t\,H)$. This is the same sign convention as `QuantumFramework`.

**A first concrete call.** Build a 1-qubit gate list for $H = X + Z$ at order 2, two Trotter steps, $t = 1/2$:

```wolfram
matOps = N @ {PauliMatrix[1], PauliMatrix[3]};
gates  = TrotterizationMatrices[matOps, 2, 2, 0.5];
Length[gates]
(* 8 *)

(* same call, assembled into a single matrix *)
u = TrotterizationMatrices[matOps, 2, 2, 0.5, "Form" -> "Unitary"];
Chop[u]
(*
  {{0.762658 - 0.464521 I, 0. - 0.450081 I},
   {0. - 0.450081 I,       0.762658 + 0.464521 I}}
*)

(* compare against the exact unitary *)
exact = MatrixExp[-I 0.5 (PauliMatrix[1] + PauliMatrix[3])];
Norm[u - exact, 2]
(* 0.0108902 *)
```

The order-2 / 2-step error at $t=0.5$ is about $10^{-2}$. We will return to convergence behaviour in §8.

**Form options agree on the same algorithm.** The three forms are different presentations of the same gate list:

```wolfram
matOps2 = N @ {
    KroneckerProduct[PauliMatrix[1], IdentityMatrix[2]],
    KroneckerProduct[IdentityMatrix[2], PauliMatrix[1]],
    KroneckerProduct[PauliMatrix[3], PauliMatrix[3]]};

flat  = TrotterizationMatrices[matOps2, 2, 2, 0.3, "Form" -> "Flat"];
steps = TrotterizationMatrices[matOps2, 2, 2, 0.3, "Form" -> "Steps"];
uMat  = TrotterizationMatrices[matOps2, 2, 2, 0.3, "Form" -> "Unitary"];

Length[flat]                                 (* 12  = reps * 6  *)
{Length[steps], Length[steps[[1]]]}          (* {2, 6}          *)
Dimensions[uMat]                             (* {4, 4}          *)
flat === Catenate[steps]                     (* True            *)
Max[Abs[Flatten[uMat - Dot @@ Reverse[flat]]]] < 10^-12   (* True *)
```

**Symbolic time.** The `time` argument is propagated symbolically. Useful for analytic work, parameter shifts, or inspecting how a gate depends on $t$:

```wolfram
matOps = {PauliMatrix[1], PauliMatrix[3]};
gates  = TrotterizationMatrices[matOps, 1, 1, t];

gates[[1]] // FullSimplify
(*  Exp[-I t X]  =  {{ Cos[t], -I Sin[t]},
                     {-I Sin[t],  Cos[t]}}  *)

gates[[2]] // FullSimplify
(*  Exp[-I t Z]  =  {{Exp[-I t], 0},
                     {0,         Exp[I t]}}  *)
```

For order 2, the coefficient list at symbolic time is exactly half-time per gate:

```wolfram
TrotterizationCoefficients[2, 2, t]
(* {t/2, t/2, t/2, t/2} *)
```

---

## 6. Examples on physical Hamiltonians

We work through four Hamiltonians of increasing size. In each case we build the term list, assemble the Trotter approximation, and read off a one-step error against `MatrixExp`.

### 6.1 One qubit: $H = X + Z$

The simplest non-trivial split: two non-commuting Pauli terms.

```wolfram
matOps = N @ {PauliMatrix[1], PauliMatrix[3]};
exact  = MatrixExp[-I 1.0 Total[matOps]];

Table[
  {order, reps,
   Norm[
     TrotterizationMatrices[matOps, order, reps, 1.0, "Form" -> "Unitary"]
     - exact, 2]},
  {order, {1, 2, 4}}, {reps, {1, 2, 4}}] // Catenate // TableForm
(*
  1  1  0.7992
  1  2  0.3624
  1  4  0.1763
  2  1  0.3137
  2  2  0.0711
  2  4  0.0173
  4  1  0.0127
  4  2  6.88 * 10^-4
  4  4  4.21 * 10^-5
*)
```

Doubling `reps` roughly halves the order-1 error, cuts the order-2 error by ~4×, and cuts the order-4 error by ~16× — Suzuki's $1/r^n$ rate.

### 6.2 Two qubits: transverse-field Ising-like

Two single-qubit X fields plus an XX interaction:

```wolfram
xx = KroneckerProduct[PauliMatrix[1], IdentityMatrix[2]];
yy = KroneckerProduct[IdentityMatrix[2], PauliMatrix[1]];
zz = KroneckerProduct[PauliMatrix[3], PauliMatrix[3]];
matOps = N @ {xx, yy, zz};

(* gate count formula: reps * Length[TrotterizationOrdering[ops, order]] *)
{Length[TrotterizationMatrices[matOps, 1, 1, 0.5]],
 Length[TrotterizationMatrices[matOps, 2, 1, 0.5]],
 Length[TrotterizationMatrices[matOps, 4, 1, 0.5]]}
(* {3, 6, 30} *)

exact  = MatrixExp[-I 0.5 Total[matOps]];
err[order_, reps_] := Norm[
  TrotterizationMatrices[matOps, order, reps, 0.5, "Form" -> "Unitary"]
  - exact, 2];

{err[1, 1], err[2, 1], err[4, 1]}
(* {0.4346, 0.1080, 0.00285} *)
```

### 6.3 Three qubits: Ising-style chain

The Hamiltonian from the original computational essay,

$$H = X_1 + X_2 + X_3 - Z_1 Z_2 - Z_2 Z_3,$$

decomposes into 5 Pauli strings. We construct it with `QuantumFramework` (only to extract Pauli matrices — the Trotterization itself is pure WL):

```wolfram
Needs["Wolfram`QuantumFramework`"];
h = QuantumOperator[
        "X" + ("X" -> {2}) + ("X" -> {3}) - "ZZ" - ("ZZ" -> {2, 3})];
qfOps  = KeyValueMap[QuantumOperator[#2 #1] &, (N @ h)["PauliDecompose"]];
matOps = Normal[#["Matrix"]] & /@ qfOps;
Length[matOps]                  (* 5 *)
Dimensions[matOps[[1]]]         (* {8, 8} *)
```

Run a $t = 1$, order 2, four-step Trotter:

```wolfram
exact  = MatrixExp[-I 1.0 Normal[(N @ h)["Matrix"]]];
sU     = TrotterizationMatrices[matOps, 2, 4, 1.0, "Form" -> "Unitary"];
Norm[sU - exact, 2]
(* 0.0637 *)
```

### 6.4 Four qubits: transverse Ising chain with field $g = 0.7$

A more demanding test:

$$H = -\sum_{i=1}^{3} Z_i Z_{i+1} - g \sum_{i=1}^{4} X_i.$$

```wolfram
g  = 0.7;
hT = QuantumOperator[
        -("ZZ" -> {1, 2}) - ("ZZ" -> {2, 3}) - ("ZZ" -> {3, 4})
        - g ("X" -> {1}) - g ("X" -> {2}) - g ("X" -> {3}) - g ("X" -> {4})];
qfOps  = KeyValueMap[QuantumOperator[#2 #1] &, (N @ hT)["PauliDecompose"]];
matOps = Normal[#["Matrix"]] & /@ qfOps;
Length[matOps]                  (* 7  *)
Dimensions[matOps[[1]]]         (* {16, 16} *)
```

Convergence at $t = 0.4$:

```wolfram
exact  = MatrixExp[-I 0.4 Normal[(N @ hT)["Matrix"]]];
err[order_, reps_] := Norm[
  TrotterizationMatrices[matOps, order, reps, 0.4, "Form" -> "Unitary"]
  - exact, 2];

Table[{order, reps, err[order, reps]},
  {order, {1, 2, 4}}, {reps, {1, 2, 4, 8}}] // Catenate // TableForm
(*
  1  1  0.4609
  1  2  0.2433
  1  4  0.1233
  1  8  0.0619
  2  1  0.0570
  2  2  0.0145
  2  4  0.0036
  2  8  9.10 * 10^-4
  4  1  3.10 * 10^-4
  4  2  1.96 * 10^-5
  4  4  1.23 * 10^-6
  4  8  7.67 * 10^-8
*)
```

The order-4 column drops by ~16× per doubling of `reps`, exactly the predicted $1/r^4$ scaling.

---

## 7. Agreement with `QuantumFramework`

The point of the standalone implementation is to be a *transparent* reference: every gate, in every order, in every step, must equal what `QuantumFramework`'s built-in `"Trotterization"` named circuit produces. This section verifies that.

### 7.1 The QF call signature

`QuantumFramework`'s named-circuit dispatch is

```wolfram
QuantumCircuitOperator["Trotterization"[ops, order, reps, time]]
```

— that is, `"Trotterization"` is applied as a *head* to its arguments. The list form `QuantumCircuitOperator[{"Trotterization", ops, order, reps, time}]` is **not** a valid call (it would be parsed as a list of separate gates, with `"Trotterization"` triggering the default 3-qubit X/Y/Z circuit and the trailing arguments becoming `PhaseShift` and `OverHat` gates). All comparisons below use the curried head form.

### 7.2 Coefficient parity

The internal `trotterCoeffs` of QF and the public `TrotterizationCoefficients` produce identical lists for every `(l, order, c)`:

```wolfram
qfCoeffs = Wolfram`QuantumFramework`NamedCircuits`PackagePrivate`trotterCoeffs;

testCases = {{2, 1, 1}, {2, 2, 1}, {2, 4, 1}, {3, 4, 1/2}, {5, 2, 1/r}};

AllTrue[testCases,
  TrotterizationCoefficients @@ # === qfCoeffs @@ # &]
(* True *)
```

### 7.3 Gate-by-gate equality

For the 3-qubit Ising chain:

```wolfram
Do[
  c       = QuantumCircuitOperator["Trotterization"[qfOps, order, reps, 0.7]];
  qfGates = Normal[#["Matrix"]] & /@ c["Flatten"]["Operators"];
  sGates  = TrotterizationMatrices[matOps, order, reps, 0.7];
  Print[order, " ", reps, ": max|qf - std| = ",
    Max[Abs[Flatten[qfGates - sGates]]]],
  {order, {1, 2, 4}}, {reps, {1, 3}}]
(*
  1 1: 1.11 * 10^-16
  1 3: 1.11 * 10^-16
  2 1: 1.24 * 10^-16
  2 3: 2.22 * 10^-16
  4 1: 2.22 * 10^-16
  4 3: 2.22 * 10^-16
*)
```

Every gate matches to within machine epsilon — the largest gap is $\sim 2 \cdot 10^{-16}$, the size of one floating-point rounding error.

### 7.4 Full unitary equality on the 4-qubit transverse Ising

```wolfram
Do[
  c   = QuantumCircuitOperator["Trotterization"[qfOps, order, reps, 0.4]];
  qfU = Normal @ c["Matrix"];
  sU  = TrotterizationMatrices[matOps, order, reps, 0.4, "Form" -> "Unitary"];
  Print[order, " ", reps, ": max|qf - std| = ",
    Max[Abs[Flatten[qfU - sU]]]],
  {order, {1, 2, 4}}, {reps, {1, 4}}]
(*
  1 1: 1.24 * 10^-15
  1 4: 1.17 * 10^-15
  2 1: 2.60 * 10^-15
  2 4: 3.46 * 10^-15
  4 1: 3.14 * 10^-15
  4 4: 3.52 * 10^-14
*)
```

The growth from $10^{-16}$ (per gate) to $10^{-14}$ (full 16×16 unitary on order 4 / 4 reps, $\approx 280$ matrix multiplies) is the expected accumulation of rounding noise over the product chain.

### 7.5 Symbolic-time equality

The package also matches QF when `time` is left symbolic — a useful property when you need a parametric circuit for analytic work:

```wolfram
matOps   = N @ {PauliMatrix[1], PauliMatrix[3]};
qfOps    = QuantumOperator /@ {"X", "Z"};

cSym     = QuantumCircuitOperator["Trotterization"[qfOps, 2, 2, t]];
qfGSym   = Normal[#["Matrix"]] & /@ cSym["Flatten"]["Operators"];
stdGSym  = TrotterizationMatrices[matOps, 2, 2, t];

(* Substitute several real values and check numeric equality. *)
Table[
  Max[Abs[Flatten[N[qfGSym /. t -> tval] - N[stdGSym /. t -> tval]]]],
  {tval, {0.1, 0.3, 0.7, 1.5}}]
(* {5.7*10^-17, 5.7*10^-17, 8.3*10^-17, 1.2*10^-16} *)
```

---

## 8. Convergence tests

Suzuki's product formula at order $n$ with $r$ Trotter steps approximates $e^{-i t H}$ with error $\mathcal{O}(t^{n+1}/r^n)$. We verify both the slope (order $n$) and the convergence ratio (≥ $2^n$ per doubling of $r$) on the 3-qubit Ising chain.

### 8.1 Error vs. `reps` at fixed time

```wolfram
exact = MatrixExp[-I 1.0 Normal[(N @ h)["Matrix"]]];

Table[
  {order, reps,
   Norm[
     TrotterizationMatrices[matOps, order, reps, 1.0, "Form" -> "Unitary"]
     - exact, 2]},
  {order, {1, 2, 4}}, {reps, {1, 2, 4, 8, 16}}] // Catenate // TableForm
(*
  1   1   1.94
  1   2   1.29
  1   4   0.675
  1   8   0.335
  1  16   0.166
  2   1   0.792
  2   2   0.238
  2   4   0.0637
  2   8   0.0162
  2  16   0.00407
  4   1   0.0411
  4   2   0.00279
  4   4   1.80 * 10^-4
  4   8   1.13 * 10^-5
  4  16   7.09 * 10^-7
*)
```

### 8.2 Reduction factors

Going from `reps = 1` to `reps = 16`:

| Order | err(r=1) | err(r=16) | ratio | predicted ($16^n$) |
|------:|---------:|----------:|------:|-------------------:|
| 1     | 1.94     | 0.166     |    12 |                 16 |
| 2     | 0.792    | 0.0041    |   195 |                256 |
| 4     | 0.0411   | 7.1·10⁻⁷  | 5.8·10⁴ |             65 536 |

Order 1 is below the predicted ratio because the leading-order error is large enough that the next-order correction in $t/r$ is still felt at $r=1$. Order 2 and order 4 are within a factor of 2 of the predicted scaling. The $2^n$-per-doubling slope holds cleanly across the table.

### 8.3 Slope check at fixed `reps`

A complementary check: at fixed `reps = 8`, increasing `order` from 1 to 4 cuts the error by ~5 orders of magnitude:

```wolfram
{Norm[TrotterizationMatrices[matOps, 1, 8, 1.0, "Form" -> "Unitary"] - exact, 2],
 Norm[TrotterizationMatrices[matOps, 2, 8, 1.0, "Form" -> "Unitary"] - exact, 2],
 Norm[TrotterizationMatrices[matOps, 4, 8, 1.0, "Form" -> "Unitary"] - exact, 2]}
(* {0.335, 0.0162, 1.13 * 10^-5} *)
```

### 8.4 Monotone decrease per order

For each order separately, the spectral error is monotone decreasing in `reps`:

```wolfram
errCol[order_] := Table[
  Norm[TrotterizationMatrices[matOps, order, reps, 1.0, "Form" -> "Unitary"]
       - exact, 2],
  {reps, {1, 2, 4, 8, 16}}];

AllTrue[{1, 2, 4}, OrderedQ[Reverse[errCol[#]]] &]
(* True *)
```

The reversed list is monotone non-decreasing, so the original is monotone non-increasing — the Trotter approximation always improves with more steps.

---

## 9. Where this leaves us

We have a small, fully transparent Suzuki–Trotter package. Three public functions cover the full surface area: `TrotterizationCoefficients` for the per-gate weights, `TrotterizationOrdering` for the term sequence, and `TrotterizationMatrices` for the assembled gates. Every function accepts symbolic arguments where it makes sense, so the package is usable both as a numerical engine and as an analytic tool.

Key results computed and verified above:

- Coefficient lists match `QuantumFramework`'s internal `trotterCoeffs` exactly, for both numeric and symbolic prefactors.
- Per-gate matrices match `QuantumFramework`'s `"Trotterization"` named circuit to ~$10^{-16}$ on every example tested (1q to 4q, orders 1, 2, 4, with up to 4 Trotter steps).
- The full assembled unitaries match to ~$10^{-14}$, the expected accumulation of rounding noise.
- Spectral error on the 3-qubit Ising chain follows the predicted $\mathcal{O}(t^{n+1}/r^n)$ scaling: order 4 at $r = 16$ reaches error $\sim 7 \cdot 10^{-7}$.
- Symbolic-time gates ($t$ left as a parameter) agree with QF after numerical substitution.

The package is a drop-in standalone reference — useful when you want Trotterization without the rest of the `QuantumFramework` stack, when you want to inspect or modify the recursion, or when you want to compare results against a transparent baseline.

---

*References:* Suzuki, M. *Decomposition formulas of exponential operators and Lie exponentials.* arXiv:[math-ph/0506007](https://arxiv.org/abs/math-ph/0506007). The `Wolfram\`QuantumFramework\`` implementation lives in [`NamedCircuits.m`](../../QuantumFramework/Kernel/QuantumCircuitOperator/NamedCircuits.m) (lines 575–615 in the current source).
