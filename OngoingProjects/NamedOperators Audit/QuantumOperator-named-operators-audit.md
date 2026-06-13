# Audit: Named Operators in `QuantumOperator`

**Date:** 2026-06-11.
**Anchor:** repo HEAD `3225a899` (working tree; `NamedOperators.m` has no uncommitted changes), paclet `Wolfram/QuantumFramework 2.0.0`.
**Scope:** all 108 entries of `$QuantumOperatorNames` (`QuantumFramework/Kernel/QuantumOperator/NamedOperators.m:10-35`), plus the generalized-Pauli helpers in `Kernel/Utilities.m:199-238` that they call. Out of scope: named circuits (`QuantumCircuitOperator`), named channels, named measurement operators, `QuantumHamiltonianOperator` named operators, and the separate `Gates.m` library.
**Method:** every definition was read in source, and every operator family was instantiated in a live kernel (`wolframscript -file`, `PacletDirectoryLoad` of the working tree) and its **computational-basis matrix** (`op["Computational"]["Matrix"]`) compared numerically (tolerance $10^{-10}$) against the conventional definition. Conventions were cross-checked against literature where there was any doubt (citations inline). No source files were modified.

**Important methodological gotcha discovered during the audit:** `op["Matrix"]` returns the matrix *in the operator's own basis*. For operators constructed in a named eigenbasis (the `J` family), this is a diagonal matrix of eigenvalues and looks "wrong" if read as a computational matrix. All verdicts below are based on `["Computational"]["Matrix"]`.

---

## Executive summary

Of the 108 named operators, **101 match their conventional definitions exactly** (verified numerically). The standard gate set (Paulis, $H$, $S$, $T$, $\sqrt{X}$, $U_1/U_2/U_3$, $R_{X,Y,Z}$, multi-Pauli $R$, controlled family, SWAP family, Toffoli, Fredkin, QFT, SUM, spiders, Cup/Cap/Trace/Reset, WignerD, angular momentum family, Liouvillian, Left/Right superoperators, Switch) is correct.

The flagged issues:

| # | Severity | Item |
|---|---|---|
| F1 | **Bug** | Case-insensitive aliases (`"cnot"`, `"toffoli"`, ...) can never fire; always return `Failure` |
| F2 | **Bug (math)** | `RX`/`RY`/`RZ` with qudit dimension $d>2$ are silently **non-unitary** |
| F3 | Convention | Generalized $Z[d]$ is the **inverse** of the standard clock matrix ($\omega^{-n}$, not $\omega^{+n}$) |
| F4 | Inconsistency | `HeisenbergWeyl` uses the standard clock phase, so `HeisenbergWeyl[d,0,1]` $\neq$ `Z[d]` |
| F5 | Convention | `Deutsch[θ]` is the literature gate $D(\theta/2)$ (half-angle); Toffoli at $\theta=\pi$ instead of $\pi/2$ |
| F6 | API trap | `QuantumOperator["SX"]` is $S\otimes X$ (two qubits), **not** $\sqrt{X}$; only `"SX"[]` is $\sqrt{X}$ |
| F7 | Label | `"Fourier"` on $n$ qudits is the per-qudit DFT $F_d^{\otimes n}$, labeled "QFT"; it is not the joint QFT |
| F8 | Ambiguity | `"0"`/`"1"` are $\mp Z$ sign-flip operators as named operators, but state-prep kets in circuit shorthand |
| F9-F13 | Notes | Naming/convention observations (ladder operators, Braid transpose, qubit-controlled qudit gates, `RandomHermitian` ensemble, unlisted aliases) |

---

## Flagged findings

### F1. BUG: case-insensitive name lookup is dead code; lowercase names always fail

[NamedOperators.m:935-941](../../QuantumFramework/Kernel/QuantumOperator/NamedOperators.m):

```wolfram
$upperCasesOperatorNames := AssociationThread[ToUpperCase @ $QuantumOperatorNames, $QuantumOperatorNames]

QuantumOperator[name_String, opts___] /; ToUpperCase[name] =!= name && KeyExistsQ[$upperCasesOperatorNames, name] :=
    QuantumOperator[$upperCasesOperatorNames[name], opts]
```

The keys of `$upperCasesOperatorNames` are the fully **uppercased** names (`"CNOT"`, `"TOFFOLI"`, ...). The guard requires both `ToUpperCase[name] =!= name` (i.e. `name` contains a lowercase letter) **and** `KeyExistsQ[..., name]` (i.e. `name` is itself an all-uppercase key). These two conditions are mutually exclusive, so the rule never fires for any input. The lookup should be keyed on `ToUpperCase[name]`.

Verified live:

```
In:  QuantumOperator["cnot"]
Out: QuantumOperator::invalidName: cnot is not a recognized QuantumOperator constructor
     Failure["InvalidName", ...]
```

Same for `"sx"`, `"Rx"[0.3]` (the parametrized rule at line 940 has the identical broken guard). This silently removes an entire documented convenience layer; it is a dispatch bug, not a physics bug.

### F2. BUG (math): qudit rotations `RX`/`RY`/`RZ` with $d>2$ are not unitary

[NamedOperators.m:156-178](../../QuantumFramework/Kernel/QuantumOperator/NamedOperators.m) define, for any dimension,

$$
R_X(\theta, d) = \exp\!\left(-\tfrac{i\theta}{2}\, X_d\right), \qquad \text{etc.}
$$

For $d=2$ this is the standard rotation (verified). For $d>2$, the generalized $X_d$, $Y_d$, $Z_d$ (shift and clock matrices) are **unitary but not Hermitian**, so $e^{-i\theta X_d/2}$ is not unitary. Verified numerically:

```
RX[1.3, 3]: U U† ≠ I   (False)
RY[1.3, 3]: U U† ≠ I   (False)
RZ[1.3, 3]: U U† ≠ I   (False)
```

The constructor signature explicitly advertises the dimension argument (`("XRotation" | "RX")[angle_, dimension : _Integer ? Positive : 2]`), and no message or assertion fires. Conventional qudit rotations use Hermitian generators (Gell-Mann matrices, or the angular momentum operators that QF itself provides correctly as `"JX"`, `"JY"`, `"JZ"`; $e^{-i\theta J_k}$ is verified unitary). Anyone constructing `RX[θ, 3]` gets a non-unitary "gate" silently.

### F3. Convention: generalized $Z[d]$ is the inverse of the standard clock matrix

[Utilities.m:229](../../QuantumFramework/Kernel/Utilities.m):

```wolfram
pauliMatrix[3, dimension_] := SparseArray[Table[{n, n} + 1 -> Exp[- 2 Pi I n / dimension], ...]]
```

QF implements $Z_d\,|j\rangle = \omega^{-j}|j\rangle$ with $\omega = e^{2\pi i/d}$. The standard convention (Weyl clock matrix; generalized Pauli / qudit stabilizer literature, e.g. [Wikipedia: Generalizations of Pauli matrices](https://en.wikipedia.org/wiki/Generalizations_of_Pauli_matrices)) is $Z_d\,|j\rangle = \omega^{+j}|j\rangle$, i.e. $\mathrm{diag}(1, \omega, \omega^2, \ldots, \omega^{d-1})$. Verified for $d=4$: QF gives $\mathrm{diag}(1, -i, -1, i)$, the conjugate of the standard $\mathrm{diag}(1, i, -1, -i)$. For $d=2$ the two coincide, so qubit code is unaffected.

The shift $X_d\,|j\rangle = |j{+}1 \bmod d\rangle$ **is** standard (verified). Knock-on effects of the nonstandard clock:

- `"PauliZ"`, `"Z"`, `"ShiftPhase"` aliases all carry it; `"Y"[d]` $= -i\,Z_d X_d$ inherits it.
- QF's own Fourier operator diagonalizes the shift into the **standard** clock, not QF's: verified $F X F^\dagger = Z_{\text{std}} \neq Z_{\text{QF}}$ while $F^\dagger X F = Z_{\text{QF}}$. So the textbook identity $F X F^\dagger = Z$ fails with QF's own pair of named operators.
- The Weyl relation holds in the inverted form: QF satisfies $XZ = \omega^{-1} ZX$ instead of $XZ = \omega\, ZX$ (equivalently it is the standard relation for $Z^\dagger$).

This is internally consistent within the Pauli family itself, but it is the mirror image of the literature convention and of QF's own `HeisenbergWeyl` (next item). Worth either documenting prominently or flipping.

### F4. Internal inconsistency: `HeisenbergWeyl` uses the standard clock, contradicting `Z[d]`

[NamedOperators.m:814-825](../../QuantumFramework/Kernel/QuantumOperator/NamedOperators.m): the Heisenberg-Weyl operator is built with phases $e^{+2\pi i a l/d}$, i.e. $T_{(i,a)} = X^i Z_{\text{std}}^{\,a}$ with the **standard** clock. Verified for $d=4$:

```
HeisenbergWeyl[4, 0, 1] == Z_standard      True
HeisenbergWeyl[4, 0, 1] == Z_QF ("Z"[4])   False
HeisenbergWeyl[4, 1, 0] == X ("X"[4])      True
```

So `QuantumOperator["HeisenbergWeyl"[d, 0, 1]]` and `QuantumOperator["Z"[d]]` disagree (they are mutual inverses). `HeisenbergWeyl` is the one that matches the literature; whichever convention is intended, the two named operators should agree.

### F5. Convention: `Deutsch[θ]` is the literature $D(\theta/2)$

[NamedOperators.m:786-793](../../QuantumFramework/Kernel/QuantumOperator/NamedOperators.m) builds the doubly-controlled block as $i\,R_x(\theta) = i\cos(\theta/2)\,\mathbb{1} + \sin(\theta/2)\,X$ (verified numerically). The conventional Deutsch gate ([Deutsch 1989; Wikipedia: Quantum logic gate](https://en.wikipedia.org/wiki/Quantum_logic_gate)) has the block

$$
D(\theta): \;|a,b,c\rangle \mapsto i\cos\theta\,|a,b,c\rangle + \sin\theta\,|a,b,1{-}c\rangle \quad (a=b=1),
$$

with full angle, and $D(\pi/2) = \text{Toffoli}$. QF's `Deutsch[θ]` $= D(\theta/2)$, and Toffoli is reached at $\theta = \pi$ (verified: `Deutsch[π] == Toffoli` exactly). A user transcribing $\theta$ from a paper gets the wrong gate by a factor of 2. (Defensible as "the $i R_x(\theta)$ convention", but it deviates from the gate's defining reference and should be documented.)

### F6. API trap: bare `"SX"` is $S \otimes X$, not $\sqrt{X}$

`"SX"` is in `$QuantumOperatorNames`, and [NamedOperators.m:277](../../QuantumFramework/Kernel/QuantumOperator/NamedOperators.m) defines `("V" | "SX")[] := √X` (the IBM SX gate, verified $\tfrac{1}{2}\begin{psmallmatrix}1+i & 1-i\\ 1-i & 1+i\end{psmallmatrix}$). But the **Pauli-string chain rule** at line 930 matches any multi-character string over `{I,X,Y,Z,H,S,T,V,P}` first, so:

```
QuantumOperator["SX"]     two-qubit operator S ⊗ X   (InputDimension 4)   <- trap
QuantumOperator["SX"[]]   single-qubit √X                                  correct
QuantumOperator["V"]      single-qubit √X                                  correct
```

Verified live. `"SX"` is the only name in the list whose characters are all in the chain alphabet, so it is the unique collision, and it collides with one of the most common hardware-native gates (IBM basis gates are `{id, rz, sx, x}`). Anyone writing `QuantumOperator["SX"]` or `"SX" -> 3` in a circuit gets a silent two-qubit $S\otimes X$.

### F7. Label/semantics: multi-qudit `"Fourier"` is $F_d^{\otimes n}$ labeled "QFT"

[NamedOperators.m:486-495](../../QuantumFramework/Kernel/QuantumOperator/NamedOperators.m): given an order of length $n$, the operator is `KroneckerProduct` of $n$ identical single-qudit DFTs. Verified:

```
QuantumOperator["Fourier", {1, 2}] == H ⊗ H     True
QuantumOperator["Fourier", {1, 2}] == QFT_4     False
Label: "QFT"
```

The $n$-qubit quantum Fourier transform $\mathrm{QFT}_{2^n}$ (Nielsen-Chuang convention, $F_{jk} = \omega^{jk}/\sqrt{N}$ over the joint register) is what the label suggests; that object lives in `QuantumCircuitOperator["Fourier"[n]]` instead. The single-qudit `"Fourier"[d]` itself is the correct DFT with $\omega = e^{+2\pi i/d}$ (verified for $d=2,3$). The per-qudit tensor-power behavior is reasonable but the `"QFT"` label on the multi-qudit form is misleading.

### F8. Ambiguity: `"0"` and `"1"` mean two different things

As named operators ([NamedOperators.m:564-566](../../QuantumFramework/Kernel/QuantumOperator/NamedOperators.m)): `"0"` $= -Z = \mathrm{diag}(-1,1)$ and `"1"` $= Z = \mathrm{diag}(1,-1)$, i.e. sign-flip (phase-oracle) operators on the indicated basis state, consistent with `FlipSign[{0}]` and `FlipSign[{1}]`.

In circuit/operator shorthand ([NamedOperators.m:78](../../QuantumFramework/Kernel/QuantumOperator/NamedOperators.m)), the same strings construct **state-injection operators**: verified `QuantumCircuitOperator[{"0"}]` produces the ket operator $\begin{psmallmatrix}1\\0\end{psmallmatrix}$.

Both are individually sensible; the collision means `QuantumOperator["0"]` and the `"0"` appearing in a circuit list are unrelated objects. Documentation should state this explicitly.

---

## Lower-severity observations

### F9. Ladder operators `"+"`/`"-"`/`"Up"`/`"Down"`: oscillator amplitudes, spin direction

[NamedOperators.m:859-861](../../QuantumFramework/Kernel/QuantumOperator/NamedOperators.m): `"+"[d]` is the superdiagonal $\sqrt{n}$ matrix, i.e. exactly the truncated harmonic-oscillator **annihilation** operator $a$ ($a|n\rangle = \sqrt{n}\,|n{-}1\rangle$; identical to QuTiP's `destroy(d)`), and `"-"[d]` $= a^\dagger$. At $d=2$, `"+"` $= \sigma_+ = |0\rangle\langle 1|$, which **raises** $J_z$ in QF's convention ($|0\rangle$ is the $m=+j$ state; QF's `JZ[1/2]` $= Z/2$, verified), so the names are coherent with the verified `J+`/`J-` family. The footgun is for $d>2$: a physicist reading `"Up"` as a Fock-space raising operator gets the lowering operator $a$, and the amplitudes are oscillator $\sqrt{n}$, not the spin amplitudes $\sqrt{(j \mp m)(j \pm m + 1)}$ of `"J+"` (compare `"+"[3]` $=$ superdiag$(1, \sqrt{2})$ vs `"J+"[1]` $=$ superdiag$(\sqrt{2}, \sqrt{2})$, both verified). Naming is convention-dependent; docs should pin the reading.

### F10. `Braid["Bell"]` is the transpose of the Kauffman-Lomonaco Bell matrix

[NamedOperators.m:523-527](../../QuantumFramework/Kernel/QuantumOperator/NamedOperators.m). QF's $R = \tfrac{1}{\sqrt 2}\begin{psmallmatrix}1&0&0&1\\0&1&1&0\\0&-1&1&0\\-1&0&0&1\end{psmallmatrix}$; the Bell matrix of Kauffman and Lomonaco (PRA **69**, 032311 (2004), "Braiding operators are universal quantum gates") is its transpose (the inverse braid generator). Both are unitary solutions of the Yang-Baxter equation; QF's `"Bell"` and `"SwapPhase"` braids were both **verified unitary and YBE-satisfying** numerically. Convention-level only.

### F11. Qudit-flavored controlled gates have a qubit control

`"CNOT"[d]`, `"CX"[d]`, `"CY"[d]`, `"CZ"[d]`, `"CPHASE"[α,d]` apply a $d$-dimensional target operation under a **2-dimensional** control (verified: `CNOT[3]` has dimensions $\{2,3\} \to \{2,3\}$); the `"C"` machinery always builds $2^k$-dimensional control registers ([NamedOperators.m:384-419](../../QuantumFramework/Kernel/QuantumOperator/NamedOperators.m)). The conventional all-qudit generalization of CNOT is `"SUM"[d]` ($|i,j\rangle \mapsto |i, i{+}j \bmod d\rangle$, verified correct for $d=2,3$). Behavior is self-consistent; just not what "qudit CNOT" means in most references.

### F12. `"RandomHermitian"` is $U + U^\dagger$ with $U \sim$ CUE, not GUE

[NamedOperators.m:625-632](../../QuantumFramework/Kernel/QuantumOperator/NamedOperators.m). The result is Hermitian (verified) with eigenvalues $2\cos\theta_k$ bounded in $[-2,2]$; it is not the Gaussian Unitary Ensemble that "random Hermitian matrix" usually denotes. Fine for generic test matrices; not fine if a user expects GUE statistics.

### F13. Name-list / dispatch inconsistencies (all verified live)

- `"CP"[α]` works (alias of `"CPHASE"` at line 285) but bare `"CP"` fails the membership check (`"CP"` not in `$QuantumOperatorNames`).
- `"Fredkin0"`, `"AngularMomentumX/Y/Z"`, `"JI+"`, `"JI-"`, `"I+"`, `"I-"`, `"BlockDiagonal"`: accepted in bracketed call forms by their rules, but absent from `$QuantumOperatorNames`, so bare-string forms fail and shorthand membership tests skip them.
- `"YSpider"` is handled by the spider rules at line 674 (and works: verified it produces a 2-in/2-out spider) but is excluded from both `$QuantumOperatorNames` and the `$Spider` pattern (line 744), so support is call-form-dependent.
- `"O"` (alias of `"Zero"`, line 153) works bare and bracketed but is unlisted.
- `"Zero"` is implemented as `Log[QuantumOperator["I"...]]`, which does produce the zero operator (verified).

### F14. Quirky defaults (harmless)

`U3[]` $=U(0,\pi,\pi) = \mathbb{1}$; `RX[]` etc. default to $\theta = \pi/2$; `Phase[]` defaults to $\pi$ (giving $Z$); `Marginal`/`Discard` default to dimension 4. None of these are errors, but the `U3` default label renders as a non-identity-looking $U(0,\pi,\pi)$.

---

## Verified correct (selected proofs)

All checks below passed numerically at tolerance $10^{-10}$ in the computational basis.

**Single-qubit:** $X, Y, Z$ standard Paulis; $H = \tfrac{1}{\sqrt2}\begin{psmallmatrix}1&1\\1&-1\end{psmallmatrix}$; $S = \mathrm{diag}(1,i)$; $T = \mathrm{diag}(1, e^{i\pi/4})$; $V = \sqrt{X}$ (IBM SX); `"0"`-projector-free phase ops as flagged above; `PhaseShift[k]` $= \mathrm{diag}(1, e^{\pm 2\pi i/2^{|k|}})$, the QFT $R_k$ gate, with $k<0$ giving the inverse (verified `PhaseShift[-2]` $= S^\dagger$); `GlobalPhase[α,d]` $= e^{i\alpha}\mathbb{1}_d$.

**IBM/OpenQASM parametrized:** `U3[θ,φ,λ]` $= \begin{psmallmatrix}\cos\frac\theta2 & -e^{i\lambda}\sin\frac\theta2\\ e^{i\varphi}\sin\frac\theta2 & e^{i(\varphi+\lambda)}\cos\frac\theta2\end{psmallmatrix}$, verified to equal the OpenQASM 3 `U` gate **including the global phase**, $U = e^{i(\varphi+\lambda)/2} R_z(\varphi) R_y(\theta) R_z(\lambda)$; `U2[φ,λ]` consistent; `U1`/`Phase`/`P` $= \mathrm{diag}(1, e^{i\alpha})$.

**Rotations (qubit):** `RX/RY/RZ[θ]` $= e^{-i\theta P/2}$, standard half-angle convention. Multi-Pauli `R[θ, spec]`: verified `R[θ,"XX"]` $= e^{-i\theta(X\otimes X)/2}$ and `R[θ,"ZZ"]` $= e^{-i\theta (Z \otimes Z)/2}$.

**Two/three-qubit:** CNOT, CX, CY, CZ, CH, CS, CT, CPHASE all standard (control $|1\rangle$, target block in lower-right); `C0NOT` correctly targets on control $|0\rangle$; SWAP; $\sqrt{\mathrm{SWAP}}$ standard; Toffoli; Fredkin/CSWAP; `Multiplexer[A,B]` $= A \oplus B$ (select/uniformly-controlled operator).

**Qudit:** `X[d]` standard shift; `SUM[d]` standard ($|i,j\rangle \to |i, i{+}j\rangle$); `Fourier[d]` standard DFT $F_{jk} = e^{2\pi i jk/d}/\sqrt d$; `InverseFourier` $= F^\dagger$; `Hadamard[d]` $= F_d$ (the common Chrestenson-gate identification); `NOT[d]` $= X_d$; `RootNOT[d]` $= X_d^{1/2}$ (principal); `Permutation[dims, perm]` correct qudit-wire permutation; `Curry`/`Uncurry` identity reshapes.

**Angular momentum (the family that *looks* wrong in the native basis but is right):** in the computational basis, `JX[j]`, `JY[j]`, `JZ[j]` are the standard spin matrices with $|0\rangle = |m{=}{+}j\rangle$: verified for $j = 1/2$ ($J_k = \sigma_k/2$) and $j=1$, and the algebra $[J_x, J_y] = i J_z$ holds. `J+`/`J-` are the standard raising/lowering operators ($\sigma_\pm$ at $j=1/2$, superdiag/subdiag $\sqrt2$ at $j=1$); `JX+`, `JY-`, `JZ+` etc. are the corresponding ladder operators in the $J_x$/$J_y$/$J_z$ eigenbases (verified `JZ+ == J+`, and `JY-` maps $|y_+\rangle \to |y_-\rangle$ with unit amplitude). Native `["Matrix"]` is $\mathrm{diag}(-j..j)$ in the operator's own eigenbasis for all three; use `["Computational"]["Matrix"]` to see the textbook matrices.

**`WignerD[j, {α,β,γ}]`** $= R_z(\alpha) R_y(\beta) R_z(\gamma)$ exactly (verified $j=1/2$ and $j=1$ against $e^{-i\alpha J_z} e^{-i\beta J_y} e^{-i\gamma J_z}$ with the descending-$m$, $|0\rangle = |m{=}{+}j\rangle$ identification), i.e. the standard active Euler-angle rotation operator (Sakurai convention). Single-angle form $= R_y(\beta)$.

**ZX/categorical family:** `ZSpider[φ]` ($1{\to}2$): $|00\rangle\langle 0| + e^{i\varphi}|11\rangle\langle 1|$, textbook; `XSpider[φ]` verified $= (H \otimes H)\, Z\text{-spider}\, H$, textbook; `WSpider` is the standard W (fan-in) spider; `Cup` $= \sum_i |ii\rangle$, `Cap` $= \sum_i \langle ii|$ (unnormalized maximally-entangled pair, standard cups/caps); `Copy` $|i\rangle \to |ii\rangle$; `Decohere`/`Measure` kill off-diagonals in the doubled (density-matrix) picture; `Encode` is the dual; `Marginal[d]` $= \sqrt d\,\langle\text{uniform}| = (1,\ldots,1)$; `Discard`/`Trace[d]` implement $\rho \mapsto \operatorname{Tr}\rho$ on doubled states; `Reset[s]` $= \operatorname{Tr}(\rho)\,|s\rangle\langle s|$ in the doubled picture. All verified.

**Open systems:** `Liouvillian[H, {L}, {γ}]` verified against the standard row-major-vectorized Lindblad generator
$$
\mathcal{L} = -i\left(H \otimes \mathbb{1} - \mathbb{1} \otimes H^{T}\right) + \sum_k \gamma_k \left( L_k \otimes L_k^{*} - \tfrac12 L_k^\dagger L_k \otimes \mathbb{1} - \tfrac12 \mathbb{1} \otimes (L_k^\dagger L_k)^{T} \right)
$$
for random $H$, $L$ (also covered by `Tests/Liouvillian.wlt`, 24/24 as of `3225a899`). `Hamiltonian[...]` $= i\mathcal{L}$, consistent with `QuantumEvolve`. `Left[A]`/`Right[A]` verified as the left/right multiplication superoperators: $\mathrm{Left}[A]\,\mathrm{vec}(\rho) = \mathrm{vec}(A\rho)$, $\mathrm{Right}[A]\,\mathrm{vec}(\rho) = \mathrm{vec}(\rho A)$.

**`Switch[A,B]`** $= |0\rangle\langle 0| \otimes A B + |1\rangle\langle 1| \otimes B A$ (verified for random CUE unitaries, including the ancilla partial-trace construction). This matches the quantum-switch form $|0\rangle\langle 0| \otimes U_1 U_2 + |1\rangle\langle 1| \otimes U_2 U_1$ used in the literature ([arXiv:2106.00034](https://arxiv.org/abs/2106.00034)). Note that some papers state the opposite prose convention ("control $|0\rangle$: first argument applied first"); QF matches the displayed-formula convention.

**Random:** `RandomUnitary` is CUE (`CircularUnitaryMatrixDistribution`), verified unitary; `RandomHermitian` verified Hermitian (see F12 for the ensemble caveat).

**Misc:** `Identity`/`I[dims]` correct including rectangular forms; `Zero` is the zero operator; `Diagonal` and `FlipSign[digits, d]` are the advertised diagonal/phase-oracle matrices (`FlipSign` places $-1$ at the basis state with the given digits, matching the `"0"`/`"1"` operators); `Braid` both variants unitary and Yang-Baxter (F10); `HeisenbergWeyl[d,i,a]` $= X^i Z_{\text{std}}^a$, the standard Heisenberg-Weyl displacement (F4 is about its mismatch with `"Z"`, not about HW itself); `Double`, `Measurement[...]`, `Channel[...]` are delegations to the corresponding heads (no independent math content).

---

## Full verdict table (all 108 names)

| Name(s) | Verdict |
|---|---|
| `Identity`, `I` | OK |
| `Permutation` | OK |
| `Uncurry`, `Curry` | OK |
| `Zero` (+unlisted alias `O`) | OK (F13: `O` unlisted) |
| `Fourier`, `InverseFourier` | OK single-qudit; **F7** multi-qudit label |
| `XRotation`/`RX`, `YRotation`/`RY`, `ZRotation`/`RZ` | OK for $d=2$; **F2** non-unitary for $d>2$ |
| `R` | OK ($e^{-i\theta P/2}$, multi-Pauli verified) |
| `U`, `U3`, `U2`, `U1`, `Phase`, `P` | OK (exact OpenQASM/IBM, incl. global phase) |
| `Diagonal` | OK |
| `GlobalPhase` | OK |
| `PhaseShift` | OK ($R_k$, negative $k$ inverse) |
| `FlipSign` | OK |
| `SUM` | OK (standard qudit SUM) |
| `RootNOT` | OK |
| `X`/`PauliX`/`Shift` | OK (standard shift) |
| `Y`/`PauliY` | OK $d=2$; $d>2$ inherits **F3** |
| `Z`/`PauliZ`/`ShiftPhase` | OK $d=2$; **F3** inverse clock for $d>2$ |
| `H`, `Hadamard` | OK ($d\neq2$: $=F_d$, Chrestenson) |
| `NOT` | OK |
| `0`, `1` | matrices coherent with `FlipSign`; **F8** dual semantics |
| `SWAP`, `RootSWAP` | OK |
| `CSWAP`/`Fredkin`, `C0SWAP` | OK |
| `Braid` | OK (unitary, YBE); **F10** transpose of K-L convention |
| `C`/`Controlled`, `C0`/`Controlled0` | OK (incl. integer-control form, matrix-operator CPM path) |
| `CX`, `CY`, `CZ`, `CH`, `CT`, `CS`, `CPHASE`, `CNOT`, `C0NOT` | OK for qubits; **F11** qubit control for $d>2$ |
| `S`, `T`, `V` | OK |
| `SX` | `"SX"[]` OK; **F6** bare `"SX"` is $S\otimes X$ |
| `Toffoli` | OK |
| `Deutsch` | **F5** half-angle vs literature |
| `RandomUnitary` | OK (CUE) |
| `RandomHermitian` | OK Hermitian; **F12** not GUE |
| `Spider`, `SimpleSpider`, `ZSpider`, `XSpider`, `WSpider` | OK (ZX/ZW semantics verified) |
| `HeisenbergWeyl` | OK vs literature; **F4** inconsistent with `Z[d]` |
| `Measure`, `Encode`, `Copy`, `Decohere`, `Marginal`, `Discard` | OK (doubled-picture semantics) |
| `Cup`, `Cap`, `Trace`, `Reset` | OK |
| `Channel`, `Measurement` | OK (delegations) |
| `Switch` | OK (matches arXiv:2106.00034 form) |
| `Multiplexer` | OK ($\oplus$ blocks) |
| `WignerD` | OK ($=R_z R_y R_z$, verified $j=\tfrac12, 1$) |
| `JX`, `JY`, `JZ` | OK in computational basis (SU(2) algebra verified); native-basis display gotcha |
| `JX+`, `JY+`, `JZ+`, `JX-`, `JY-`, `JZ-`, `J+`, `J-` | OK (standard ladder amplitudes, per-axis eigenbases) |
| `+`/`Up`, `-`/`Down` | OK matrices ($\sigma_\pm$ at $d=2$; truncated $a$, $a^\dagger$); **F9** naming |
| `Double` | OK (delegation) |
| `Left`, `Right` | OK (multiplication superoperators verified) |
| `Hamiltonian`, `Liouvillian` | OK (standard vectorized Lindblad, verified + test suite) |

Also: case-insensitive aliases for **all** names are broken (**F1**).

---

## Suggested priorities (report only; no source changes made)

1. **F1** (one-line fix: key the lookup on `ToUpperCase[name]`).
2. **F2** (either restrict `RX/RY/RZ` to $d=2$, switch the $d>2$ generator to a Hermitian one, or emit a message).
3. **F6** (reorder the chain rule below the named-operator membership rule, or drop `"SX"` from the chain alphabet collision).
4. **F3/F4** (pick one clock convention; `HeisenbergWeyl` currently matches the literature, `Z[d]` does not).
5. **F5, F7, F8** (documentation, or deliberate convention changes).
