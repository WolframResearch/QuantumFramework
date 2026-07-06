---
Template: Default
---

# Quantum in Finite Dimensions: Answers in WL and in QuantumFramework (Part 8)

A companion answer key to `Question-List.md`. For each question it gives **two** worked answers:

1. **WL** : native Wolfram Language only (plain vectors and matrices, `PauliMatrix`, `MatrixExp`,
   `WignerD`, `ClebschGordan`, `Eigensystem`, and friends). Nothing beyond the built-in language.
2. **QF** : the same task done through the **QuantumFramework** paclet, leaning on its objects and
   property downvalues (`QuantumState[...]["..."]`) rather than rebuilding matrices by hand. A QF
   result is a rich object, a state or an operator, not a bare array: it can act on a target, compose
   with others, report its spectrum, change basis, and more. The answers keep results as objects and
   treat a matrix or amplitude vector as just one of the object's many representations, not the thing
   itself.

Three habits run through the answers. They stay **symbolic** wherever possible (general angles
$\theta,\phi$), so each operation is exercised in full generality rather than on a lucky special case.
They write the computation **explicitly** rather than hiding it in a helper function, because the code
is itself part of the explanation: the spin-$j$ operators are built in place from the ladder relation
$J_\pm|j,m\rangle = \sqrt{j(j+1)-m(m\pm1)}\,|j,m\pm1\rangle$, never pasted in as fixed matrices. And
they prefer **built-in features** (a `WignerD`, a `ClebschGordan`) over hand-rolled constructions.
Each cell does **one** thing and is preceded by a sentence saying what. We use the physicist's
convention $\hbar = 1$, so a spin-$j$ has $J^2 = j(j+1)$ and $J_z$ eigenvalues $m = -j,\ldots,j$. Run
the WL answers in a bare kernel; the QF answers need the paclet loaded once.

## Setup

The QF answers load the paclet once:

```wl
Needs["Wolfram`QuantumFramework`"]
```

---

## Part 8. Spin and angular momentum

### 8.1 [BSc] How do I model a spin-1/2 in a magnetic field and reproduce the Stern-Gerlach spin-projection outcomes?

A spin-1/2 in a field along $z$ has $\hat H = -\gamma B\,S_z$ with $S_z = \tfrac12 Z$. A Stern-Gerlach
apparatus measures $S_z$, which has just two eigenvalues $\pm\tfrac12$: the beam splits in two, never a
continuum. Prepare a spin pointing along a general direction $\hat n(\theta,\varphi)$, the state
$\{\cos\tfrac\theta2, e^{i\varphi}\sin\tfrac\theta2\}$, and find the two deflection outcomes and their
Born weights.

**WL** : the outcomes are the eigenvalues of $S_z$, and the weights are the squared amplitudes in the
$S_z$ eigenbasis.

```wl
With[{sz = PauliMatrix[3]/2, \[Psi] = {Cos[\[Theta]/2], Exp[I \[Phi]] Sin[\[Theta]/2]}}, {Eigenvalues[sz], ComplexExpand[Abs[\[Psi]]^2]}]
```

**QF** : a `QuantumMeasurementOperator` for $S_z$ measures the state and reports its outcomes (the
eigenvalues) with their probabilities.

```wl
{QuantumOperator[PauliMatrix[3]/2]["Eigenvalues"], ComplexExpand[QuantumMeasurementOperator[QuantumOperator[PauliMatrix[3]/2]][QuantumState[{Cos[\[Theta]/2], Exp[I \[Phi]] Sin[\[Theta]/2]}]]["ProbabilitiesList"]]}
```

Both give the two spin projections $\pm\tfrac12$ with probabilities $\{\cos^2\tfrac\theta2,
\sin^2\tfrac\theta2\}$: the apparatus reveals the discreteness of angular momentum (only two beams for a
spin-1/2), and a spin tilted at angle $\theta$ from the field lands in the upper beam with probability
$\cos^2\tfrac\theta2$, the quantitative Stern-Gerlach result.

### 8.2 [BSc] How do I build the spin-$j$ angular-momentum operators and verify $[J_x,J_y]=iJ_z$ and the Casimir?

The angular-momentum operators are defined by their algebra: the raising and lowering operators
$J_\pm|j,m\rangle = \sqrt{j(j+1)-m(m\pm1)}\,|j,m\pm1\rangle$ build $J_x = \tfrac12(J_+ + J_-)$ and
$J_y = \tfrac1{2i}(J_+ - J_-)$, and $J_z$ is diagonal with entries $m$. They satisfy $[J_x,J_y]=iJ_z$
and the Casimir $J^2 = J_x^2+J_y^2+J_z^2 = j(j+1)\,I$. Build them for spin $j=1$ (a qutrit) from the
ladder relation.

**WL** : construct $J_z$ and $J_+$ from the definition, then $J_x, J_y$, and check both identities.

```wl
With[{j = 1}, With[{m = Range[j, -j, -1], n = 2 j + 1}, With[{jz = DiagonalMatrix[Range[j, -j, -1]], jp = SparseArray[{r_, c_} /; r == c - 1 :> Sqrt[j (j + 1) - m[[c]] (m[[c]] + 1)], {n, n}]}, With[{jx = (jp + Transpose[jp])/2, jy = (jp - Transpose[jp])/(2 I)}, {jx . jy - jy . jx == I jz, jx . jx + jy . jy + jz . jz == j (j + 1) IdentityMatrix[n]}]]]]
```

**QF** : keep the three components as `QuantumOperator`s and verify the algebra with `Commutator` and the
operator product `@`, at the operator level.

```wl
With[{j = 1}, With[{m = Range[j, -j, -1], n = 2 j + 1}, With[{jz = DiagonalMatrix[Range[j, -j, -1]], jp = SparseArray[{r_, c_} /; r == c - 1 :> Sqrt[j (j + 1) - m[[c]] (m[[c]] + 1)], {n, n}]}, With[{jx = QuantumOperator[(jp + Transpose[jp])/2], jy = QuantumOperator[(jp - Transpose[jp])/(2 I)], jzo = QuantumOperator[jz]}, {Commutator[jx, jy] == I jzo, jx @ jx + jy @ jy + jzo @ jzo == j (j + 1) QuantumOperator[IdentityMatrix[n]]}]]]]
```

Both return $\{$True, True$\}$: the spin-1 operators close under the $\mathfrak{su}(2)$ commutator and
their Casimir is the scalar $j(j+1) = 2$, the defining features of an angular momentum. The same
construction gives any spin $j$ by changing the one number.

### 8.3 [BSc] How do I rotate a spin state and read the rotation off the expectation values?

A spatial rotation of a spin-$j$ acts by $R_{\hat n}(\beta) = e^{-i\beta\,\hat n\cdot\vec J}$. Rotating
the highest-weight state $|j,j\rangle$ (which points along $+z$, with $\langle J_z\rangle = j$) about the
$y$ axis by $\beta$ tilts the spin into the $x$-$z$ plane: the expectation vector becomes
$\langle\vec J\rangle = j\,(\sin\beta, 0, \cos\beta)$. Show this for spin $j=1$.

**WL** : rotate $|1,1\rangle$ by $e^{-i\beta J_y}$ and read the three expectation values.

```wl
With[{j = 1}, With[{m = Range[j, -j, -1], n = 2 j + 1}, With[{jz = DiagonalMatrix[Range[j, -j, -1]], jp = SparseArray[{r_, c_} /; r == c - 1 :> Sqrt[j (j + 1) - m[[c]] (m[[c]] + 1)], {n, n}]}, With[{jx = (jp + Transpose[jp])/2, jy = (jp - Transpose[jp])/(2 I)}, With[{\[Psi] = MatrixExp[-I \[Beta] jy] . UnitVector[n, 1]}, Simplify[ComplexExpand[{Conjugate[\[Psi]] . jx . \[Psi], Conjugate[\[Psi]] . jy . \[Psi], Conjugate[\[Psi]] . jz . \[Psi]}], \[Beta] \[Element] Reals]]]]]]
```

**QF** : the rotation is a `QuantumOperator` exponential acting on the state; the expectations come from
the operator sandwich on the rotated `QuantumState`.

```wl
With[{j = 1}, With[{m = Range[j, -j, -1], n = 2 j + 1}, With[{jz = DiagonalMatrix[Range[j, -j, -1]], jp = SparseArray[{r_, c_} /; r == c - 1 :> Sqrt[j (j + 1) - m[[c]] (m[[c]] + 1)], {n, n}]}, With[{jx = (jp + Transpose[jp])/2, jy = (jp - Transpose[jp])/(2 I)}, With[{\[Psi] = QuantumState[MatrixExp[-I \[Beta] jy, UnitVector[n, 1]]]}, Simplify[ComplexExpand[Table[(\[Psi]["Dagger"] @ QuantumOperator[op] @ \[Psi])["Scalar"], {op, {jx, jy, jz}}]], \[Beta] \[Element] Reals]]]]]]
```

Both return $\langle\vec J\rangle = (\sin\beta, 0, \cos\beta)$: the quantum expectation vector rotates
exactly as a classical arrow would, so a spin measurement tracks the geometric rotation, the content of
the vector (adjoint) representation of $SU(2)$ on angular momentum.

### 8.4 [MSc] How do I construct the Wigner $D$-matrix element $D^j_{m'm}(\alpha,\beta,\gamma)=\langle j\,m'|e^{-i\alpha J_z}e^{-i\beta J_y}e^{-i\gamma J_z}|j\,m\rangle$ for a spin-$j$ rotation?

The Wigner $D$-matrix is the spin-$j$ representation of a rotation parametrized by Euler angles,
$D^j_{m'm}(\alpha,\beta,\gamma) = \langle j\,m'|e^{-i\alpha J_z}e^{-i\beta J_y}e^{-i\gamma J_z}|j\,m\rangle$.
Build the full rotation for spin $j=1$ from the ladder operators and read off its matrix; compare it to
the built-in `WignerD`.

**WL** : assemble $e^{-i\alpha J_z}e^{-i\beta J_y}e^{-i\gamma J_z}$ from the spin-1 operators, and check
its entries against the built-in `WignerD` (index order $\{j, m', m\}$). The built-in uses the opposite
(active, $e^{+iJ}$) sign convention, so its Euler angles are negated relative to the passive $e^{-iJ}$
formula stated here.

```wl
With[{j = 1}, With[{m = Range[j, -j, -1], n = 2 j + 1}, With[{jz = DiagonalMatrix[Range[j, -j, -1]], jp = SparseArray[{r_, c_} /; r == c - 1 :> Sqrt[j (j + 1) - m[[c]] (m[[c]] + 1)], {n, n}]}, With[{jy = (jp - Transpose[jp])/(2 I)}, With[{d = MatrixExp[-I \[Alpha] jz] . MatrixExp[-I \[Beta] jy] . MatrixExp[-I \[Gamma] jz]}, Simplify[d == Table[WignerD[{j, mp, mm}, -\[Alpha], -\[Beta], -\[Gamma]], {mp, j, -j, -1}, {mm, j, -j, -1}], {\[Alpha], \[Beta], \[Gamma]} \[Element] Reals]]]]]]
```

**QF** : the same rotation as a product of `QuantumOperator` exponentials; its matrix is the $D$-matrix.

```wl
With[{j = 1}, With[{m = Range[j, -j, -1], n = 2 j + 1}, With[{jz = DiagonalMatrix[Range[j, -j, -1]], jp = SparseArray[{r_, c_} /; r == c - 1 :> Sqrt[j (j + 1) - m[[c]] (m[[c]] + 1)], {n, n}]}, With[{jy = (jp - Transpose[jp])/(2 I)}, With[{d = QuantumOperator[MatrixExp[-I \[Alpha] jz]] @ QuantumOperator[MatrixExp[-I \[Beta] jy]] @ QuantumOperator[MatrixExp[-I \[Gamma] jz]]}, Simplify[Normal[d["Matrix"]] == Table[WignerD[{j, mp, mm}, -\[Alpha], -\[Beta], -\[Gamma]], {mp, j, -j, -1}, {mm, j, -j, -1}], {\[Alpha], \[Beta], \[Gamma]} \[Element] Reals]]]]]]
```

Both confirm the constructed rotation equals the built-in `WignerD` table entry by entry: the $D$-matrix
is nothing more than the explicit representation of the abstract rotation on the $(2j+1)$-dimensional spin
space, and building it from the ladder operators reproduces the tabulated special functions exactly.

### 8.5 [MSc] How do I add two angular momenta and compute the Clebsch-Gordan coefficients?

Coupling two spins combines their spaces, $\tfrac12\otimes\tfrac12 = 0\oplus1$: the product basis
$|m_1\rangle|m_2\rangle$ recombines into total-spin eigenstates $|J,M\rangle$, and the overlaps
$\langle m_1 m_2|J,M\rangle$ are the Clebsch-Gordan coefficients. Build the coupled basis for two
spin-1/2 and read the coefficients.

**WL** : the coupled states are assembled from the product basis with the built-in `ClebschGordan`
coefficients; here is $|1,0\rangle$, the symmetric $M=0$ triplet state.

```wl
With[{basis = {{1/2, 1/2}, {1/2, -1/2}, {-1/2, 1/2}, {-1/2, -1/2}}}, Sum[If[Total[basis[[k]]] == 0, ClebschGordan[{1/2, basis[[k, 1]]}, {1/2, basis[[k, 2]]}, {1, 0}] UnitVector[4, k], 0], {k, 4}]]
```

**QF** : the coupled basis is the eigenbasis of the total spin $\vec S = \vec S_1 + \vec S_2$; diagonalize
$S^2$ and its overlaps with the product basis *are* the Clebsch-Gordan coefficients. The $|1,0\rangle$
state is the $M=0$ ($S_z^{\mathrm{tot}}=0$) eigenvector inside the $S^2=2$ (triplet) block.

```wl
With[{stot = Table[KroneckerProduct[PauliMatrix[k]/2, IdentityMatrix[2]] + KroneckerProduct[IdentityMatrix[2], PauliMatrix[k]/2], {k, 3}]}, With[{s2 = Sum[stot[[k]] . stot[[k]], {k, 3}], sz = stot[[3]]}, Select[Normalize /@ (Normal /@ QuantumOperator[s2]["Eigenvectors"]), Simplify[# . s2 . # == 2 && # . sz . # == 0] &]]]
```

Both single out $|1,0\rangle = \tfrac1{\sqrt2}(|{\uparrow\downarrow}\rangle + |{\downarrow\uparrow}\rangle)$,
its Clebsch-Gordan coefficients $\{0, \tfrac1{\sqrt2}, \tfrac1{\sqrt2}, 0\}$: adding two spins is a change
of basis from the product basis to the total-spin basis, and the change-of-basis matrix is the table of
Clebsch-Gordan coefficients, recoverable either from the built-in or by diagonalizing $S^2$.

### 8.6 [MSc] How do I build the singlet and triplet states and identify them by total spin?

Two spin-1/2 combine into one singlet ($S=0$, antisymmetric) and three triplet states ($S=1$, symmetric).
They are told apart by the total-spin Casimir $S^2 = (\vec S_1+\vec S_2)^2$, whose eigenvalue is
$S(S+1)$: $0$ for the singlet, $2$ for the triplet. Build the four states and read their $S^2$.

**WL** : form $S^2$ from the two spins, and evaluate it on the singlet and on a triplet representative.

```wl
With[{singlet = {0, 1, -1, 0}/Sqrt[2], triplet = {0, 1, 1, 0}/Sqrt[2]}, With[{stot = Table[KroneckerProduct[PauliMatrix[k]/2, IdentityMatrix[2]] + KroneckerProduct[IdentityMatrix[2], PauliMatrix[k]/2], {k, 3}]}, With[{s2 = Sum[stot[[k]] . stot[[k]], {k, 3}]}, {Conjugate[singlet] . s2 . singlet, Conjugate[triplet] . s2 . triplet}]]]
```

**QF** : the same $S^2$ as a `QuantumOperator`; its spectrum on the whole space shows the multiplicities
directly, one $S^2=0$ level and three $S^2=2$ levels.

```wl
With[{stot = Table[QuantumTensorProduct[QuantumOperator[PauliMatrix[k]/2], QuantumOperator["I"]] + QuantumTensorProduct[QuantumOperator["I"], QuantumOperator[PauliMatrix[k]/2]], {k, 3}]}, Sort[Sum[stot[[k]] @ stot[[k]], {k, 3}]["Eigenvalues"]]]
```

The singlet returns $S^2=0$ and every triplet state returns $S^2=2=1(1+1)$; the QF spectrum is
$\{0,2,2,2\}$, exactly the $1\oplus3$ multiplicity of $\tfrac12\otimes\tfrac12 = 0\oplus1$. Total spin is
the conserved label that block-diagonalizes any rotationally invariant two-spin interaction (as in Part
7, 7.10).

### 8.7 [MSc] How do I construct a spin coherent state and place it on the generalized Bloch sphere?

A spin coherent state is the most classical state of a spin-$j$: the highest-weight state $|j,j\rangle$
rotated to point along a direction $\hat n(\theta,\varphi)$, $|j,\hat n\rangle = e^{-i\varphi J_z}
e^{-i\theta J_y}|j,j\rangle$. Its expectation vector is fully polarized, $\langle\vec J\rangle = j\,\hat n$,
so it sits on the sphere of radius $j$, the generalized Bloch sphere. Build it for spin $j=1$.

**WL** : rotate the highest-weight state and confirm $\langle\vec J\rangle = j\,\hat n$.

```wl
With[{j = 1}, With[{m = Range[j, -j, -1], n = 2 j + 1}, With[{jz = DiagonalMatrix[Range[j, -j, -1]], jp = SparseArray[{r_, c_} /; r == c - 1 :> Sqrt[j (j + 1) - m[[c]] (m[[c]] + 1)], {n, n}]}, With[{jx = (jp + Transpose[jp])/2, jy = (jp - Transpose[jp])/(2 I)}, With[{coh = MatrixExp[-I \[Phi] jz] . MatrixExp[-I \[Theta] jy] . UnitVector[n, 1]}, Simplify[ComplexExpand[{Conjugate[coh] . jx . coh, Conjugate[coh] . jy . coh, Conjugate[coh] . jz . coh}], {\[Theta], \[Phi]} \[Element] Reals]]]]]]
```

**QF** : the same coherent state as a `QuantumState`; its polarization vector, the three spin
expectations, points along $\hat n$ with length $j$.

```wl
With[{j = 1}, With[{m = Range[j, -j, -1], n = 2 j + 1}, With[{jz = DiagonalMatrix[Range[j, -j, -1]], jp = SparseArray[{r_, c_} /; r == c - 1 :> Sqrt[j (j + 1) - m[[c]] (m[[c]] + 1)], {n, n}]}, With[{jx = (jp + Transpose[jp])/2, jy = (jp - Transpose[jp])/(2 I)}, With[{coh = QuantumState[MatrixExp[-I \[Phi] jz] . MatrixExp[-I \[Theta] jy] . UnitVector[n, 1]]}, Simplify[ComplexExpand[Table[(coh["Dagger"] @ QuantumOperator[op] @ coh)["Scalar"], {op, {jx, jy, jz}}]], {\[Theta], \[Phi]} \[Element] Reals]]]]]]
```

Both give $\langle\vec J\rangle = (\sin\theta\cos\varphi, \sin\theta\sin\varphi, \cos\theta)$, exactly
$j\,\hat n$ with $j=1$: the spin coherent state points as definitely as a spin can in the direction
$\hat n$, saturating the angular-momentum uncertainty and serving as the finite-dimensional stand-in for
the (infinite-dimensional) coherent state of light.
