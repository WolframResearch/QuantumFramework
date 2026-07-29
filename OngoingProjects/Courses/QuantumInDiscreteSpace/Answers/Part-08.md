## Part 8. Spin and angular momentum

Spin is the cleanest finite-dimensional laboratory for angular momentum: every operator is a small matrix
and every state a short vector, so the whole $SU(2)$ story fits in a few lines. Each answer below is
self-contained, building the spin operators and states it needs from scratch.

### 8.1 [BSc] How do I model a spin-1/2 in a magnetic field and reproduce the Stern-Gerlach spin-projection outcomes?

A spin-1/2 in a field along $z$ has $\hat H = -\gamma B\,S_z$ with $S_z = \tfrac12 Z$. A Stern-Gerlach
apparatus measures $S_z$, whose two eigenvalues $\pm\tfrac12$ are the only outcomes: the beam splits in two,
never a continuum. Prepare a spin along $\hat n(\theta,\varphi)$ and read the outcomes, their Born weights,
and the mean.

```wl
sz = PauliMatrix[3]/2; \[Psi] = {Cos[\[Theta]/2], Exp[I \[Phi]] Sin[\[Theta]/2]};
```

**WL** : the outcomes are the diagonal entries of $S_z$ and the weights the squared amplitudes; `Thread`
pairs each outcome with its probability:

```wl
Thread[Diagonal[sz] -> ComplexExpand[Abs[\[Psi]]^2]]
```

The mean is $\langle\psi|S_z|\psi\rangle$:

```wl
ComplexExpand[Conjugate[\[Psi]] . sz . \[Psi]]
```

**QF** : a `QuantumMeasurementOperator` measures the state and returns its outcome-keyed `"Probabilities"`
and its `"Mean"`:

```wl
meas = QuantumMeasurementOperator[QuantumOperator[sz]][QuantumState[\[Psi]]];
```

```wl
Simplify@*ComplexExpand /@ meas["Probabilities"]
```

```wl
meas["Mean"] // ComplexExpand // Simplify
```

Both give outcomes $\pm\tfrac12$ with probabilities $\cos^2\tfrac\theta2, \sin^2\tfrac\theta2$ and mean
$\langle S_z\rangle = \tfrac12\cos\theta$: the apparatus reveals the discreteness of angular momentum, and a
spin tilted at $\theta$ from the field lands in the upper beam with probability $\cos^2\tfrac\theta2$. The
azimuth $\varphi$ drops out, because an $S_z$ measurement cannot see it.

### 8.2 [BSc] How do I build the spin-$j$ angular-momentum operators and verify the full algebra $[J_a,J_b]=i\,\epsilon_{abc}J_c$ and the Casimir?

The angular-momentum operators are defined by the ladder relation $J_\pm|j,m\rangle =
\sqrt{j(j+1)-m(m\pm1)}\,|j,m\pm1\rangle$, from which $J_x = \tfrac12(J_+ + J_-)$, $J_y = \tfrac1{2i}(J_+ -
J_-)$, and $J_z = \mathrm{diag}(m)$. They must close under the full $\mathfrak{su}(2)$ algebra
$[J_a,J_b]=i\sum_c\epsilon_{abc}J_c$ ($\epsilon_{abc}$ the Levi-Civita structure constants) and give the
Casimir $J^2 = j(j+1)\,I$. Build them for spin $j=1$.

**WL** : $J_z$ is diagonal and $J_+$ carries one super-diagonal; collect $\vec J = \{J_x, J_y, J_z\}$:

```wl
j = 1; m = Range[j, -j, -1]; jp = SparseArray[{r_, c_} /; r == c - 1 :> Sqrt[j (j + 1) - m[[c]] (m[[c]] + 1)], {2 j + 1, 2 j + 1}]; jvec = {(jp + Transpose[jp])/2, (jp - Transpose[jp])/(2 I), DiagonalMatrix[m]};
```

All nine commutators close on the structure constants, checked in one tensor contraction against
`LeviCivitaTensor[3]`:

```wl
Table[jvec[[a]] . jvec[[b]] - jvec[[b]] . jvec[[a]], {a, 3}, {b, 3}] == I LeviCivitaTensor[3] . jvec
```

And the Casimir is the scalar $j(j+1)\,I$:

```wl
Total[#.# & /@ jvec] == j (j + 1) IdentityMatrix[2 j + 1]
```

The same operators follow from group theory, because $\vec J$ is a rank-1 *spherical tensor*: its
**spherical components** are $J^{\mathrm{sph}}_{+1} = -\tfrac1{\sqrt2}(J_x+iJ_y)$, $J^{\mathrm{sph}}_0 =
J_z$, $J^{\mathrm{sph}}_{-1} = +\tfrac1{\sqrt2}(J_x-iJ_y)$, labeled by $\mu = +1,0,-1$ and rotating like the
$j{=}1$ kets $|1,\mu\rangle$. The Wigner-Eckart theorem fixes each matrix element as one Clebsch-Gordan coefficient times a single reduced
matrix element, equal to $\sqrt{j(j+1)}$ for $\vec J$ itself:

$$\langle j,m'|\,J^{\mathrm{sph}}_\mu\,|j,m\rangle = \sqrt{j(j+1)}\;\big(\langle j,m|\otimes\langle 1,\mu|\big)\,\big|\,J{=}j,\;M{=}m'\,\big\rangle.$$

The coefficient on the right is built by adding two angular momenta. The first is the physical spin: the
operator $\vec J$, with angular momentum quantum number $j$, acting on the $(2j+1)$-dimensional space
$\mathcal H_j$ whose orthonormal basis is $|j,m\rangle$ for $m=j,j-1,\dots,-j$. The second is a spin-1 angular
momentum $\vec S$ attached to the tensor index, acting on the $3$-dimensional space $\mathcal H_1$ with basis
$|1,\mu\rangle$; it enters not as a second physical particle but as the representation space in which the
operator's three spherical components transform. On the product space $\mathcal H_j\otimes\mathcal H_1$, of
dimension $3(2j+1)$, the total angular momentum is the operator $\vec J_{\mathrm{tot}}=\vec J\otimes I+I\otimes\vec S$.
Its square $\vec J_{\mathrm{tot}}^{\,2}$ has eigenvalue $J(J+1)$ on three orthogonal subspaces labeled
$J=j+1,\,j,\,j-1$ (for $j\ge1$; at $j=\tfrac12$ the $J=j-1$ subspace is absent), so the product space is the
direct sum

$$\mathcal H_j\otimes\mathcal H_1=\mathcal H_{j+1}\oplus\mathcal H_j\oplus\mathcal H_{j-1}.$$

A *coupled state* $|J,M\rangle$ is the simultaneous eigenstate of $\vec J_{\mathrm{tot}}^{\,2}$, with eigenvalue
$J(J+1)$, and of the total projection $(J_{\mathrm{tot}})_z$, with eigenvalue $M$. Expanded in the product basis
it reads

$$|J,M\rangle=\sum_{m+\mu=M}\big(\langle j,m|\otimes\langle 1,\mu|\big)|J,M\rangle\;\;|j,m\rangle\otimes|1,\mu\rangle,$$

and the numbers $\big(\langle j,m|\otimes\langle 1,\mu|\big)|J,M\rangle$ weighting each product state are the
Clebsch-Gordan coefficients, ordinary inner products taken entirely inside $\mathcal H_j\otimes\mathcal H_1$.

The coefficient in the Wigner-Eckart formula is this number for the middle subspace $J=j$, at total projection
$M=m'$. The total quantum number equals $j$, the same value as the physical spin, because $J^{\mathrm{sph}}_\mu$
carries the spin-$j$ multiplet back into itself, so only the $J=j$ subspace contributes and the $J=j\pm1$ ones
drop out. The right-hand ket $\big|\,J{=}j,\;M{=}m'\,\big\rangle$ is therefore a vector in the
$3(2j+1)$-dimensional product space, a different object from the physical single-spin state $|j,m'\rangle$ on
the left of the matrix element, which lives in the $(2j+1)$-dimensional space $\mathcal H_j$; they share the
label $(j,m')$ only through this coincidence of values. The coefficient is nonzero only when $M=m+\mu$, that
is $m'=m+\mu$. In `ClebschGordan` the first two arguments $\{j,m\}$ and $\{1,\mu\}$ are the two states that
are combined, and the third, $\{j,m'\}$, is the coupled state
$\big|\,J{=}j,\;M{=}m'\,\big\rangle$. Reading the three components straight from `ClebschGordan`:

```wl
lsph[\[Mu]_] := Table[If[m[[r]] == m[[c]] + \[Mu], Sqrt[j (j + 1)] ClebschGordan[{j, m[[c]]}, {1, \[Mu]}, {j, m[[r]]}], 0], {r, 2 j + 1}, {c, 2 j + 1}];
```

Unpacking $J_x, J_y, J_z$ reproduces the ladder operators exactly:

```wl
{(lsph[-1] - lsph[1])/Sqrt[2], I (lsph[-1] + lsph[1])/Sqrt[2], lsph[0]} == jvec // Simplify
```

**QF** : the built-in named operators `"JX"[j]`, `"JY"[j]`, `"JZ"[j]` carry the same algebra; `Commutator`
and the product `@` reconcile their eigenbases:

```wl
jq = QuantumOperator[#[j]] & /@ {"JX", "JY", "JZ"};
```

```wl
And @@ Flatten@Table[Commutator[jq[[a]], jq[[b]]] == I Sum[LeviCivitaTensor[3][[a, b, c]] jq[[c]], {c, 3}], {a, 3}, {b, 3}]
```

```wl
jq[[1]] @ jq[[1]] + jq[[2]] @ jq[[2]] + jq[[3]] @ jq[[3]] == j (j + 1) QuantumOperator[IdentityMatrix[2 j + 1]]
```

Every route returns True: the spin-1 operators close under $[J_a,J_b]=i\epsilon_{abc}J_c$ with Casimir
$j(j+1)=2$, whether built from the ladder recursion, from the Clebsch-Gordan matrix elements of the
spherical tensor $\vec J$, or as the framework's named operators. Any spin $j$ follows by changing the one
number.

### 8.3 [BSc] How do I rotate a spin state and read the rotation off the expectation values?

A rotation of a spin-$j$ acts by $R_{\hat n}(\beta) = e^{-i\beta\,\hat n\cdot\vec J}$. Rotating the
highest-weight state $|j,j\rangle$ (which points along $+z$) about $y$ by $\beta$ tilts it into the $x$-$z$
plane, so the expectation vector becomes $\langle\vec J\rangle = j\,(\sin\beta, 0, \cos\beta)$. Show it for
spin $j=1$.

**WL** : build the spin-1 operators:

```wl
j = 1; m = Range[j, -j, -1]; jp = SparseArray[{r_, c_} /; r == c - 1 :> Sqrt[j (j + 1) - m[[c]] (m[[c]] + 1)], {2 j + 1, 2 j + 1}]; {jx, jy, jz} = {(jp + Transpose[jp])/2, (jp - Transpose[jp])/(2 I), DiagonalMatrix[m]};
```

Rotate the top state $|1,1\rangle = \{1,0,0\}$ about $y$ by $\beta$:

```wl
\[Psi] = MatrixExp[-I \[Beta] jy] . UnitVector[2 j + 1, 1];
```

and read the three expectations $\langle\psi|J_k|\psi\rangle$:

```wl
Simplify[ComplexExpand[Conjugate[\[Psi]] . # . \[Psi] & /@ {jx, jy, jz}], \[Beta] \[Element] Reals]
```

**QF** : rotate with the framework's own operator, the exponential of the named `"JY"[j]` acting on the top
state $|0\rangle$, then read the expectations from the named operators:

```wl
\[Psi]q = Exp[-I \[Beta] QuantumOperator["JY"[j]]][QuantumState["0", 2 j + 1]];
```

```wl
Simplify[ComplexExpand[Table[(\[Psi]q["Dagger"] @ QuantumOperator[op[j]] @ \[Psi]q)["Scalar"], {op, {"JX", "JY", "JZ"}}]], \[Beta] \[Element] Reals]
```

Both return $\langle\vec J\rangle = (\sin\beta, 0, \cos\beta)$: the quantum expectation vector rotates
exactly as a classical arrow would, the content of the vector (adjoint) representation of $SU(2)$ on angular
momentum.

### 8.4 [MSc] How do I construct the Wigner $D$-matrix element $D^j_{m'm}(\alpha,\beta,\gamma)=\langle j\,m'|e^{-i\alpha J_z}e^{-i\beta J_y}e^{-i\gamma J_z}|j\,m\rangle$ for a spin-$j$ rotation?

The Wigner $D$-matrix is the spin-$j$ representation of a rotation parametrized by Euler angles: the matrix
of the abstract rotation once we choose the $|j,m\rangle$ basis. Build it for spin $j=1$ from the ladder
operators and match it to the built-in `WignerD`.

**WL** : build $J_y, J_z$:

```wl
j = 1; m = Range[j, -j, -1]; jp = SparseArray[{r_, c_} /; r == c - 1 :> Sqrt[j (j + 1) - m[[c]] (m[[c]] + 1)], {2 j + 1, 2 j + 1}]; {jy, jz} = {(jp - Transpose[jp])/(2 I), DiagonalMatrix[m]};
```

Assemble the Euler rotation $e^{-i\alpha J_z}e^{-i\beta J_y}e^{-i\gamma J_z}$:

```wl
dmat = MatrixExp[-I \[Alpha] jz] . MatrixExp[-I \[Beta] jy] . MatrixExp[-I \[Gamma] jz];
```

Compare it entry by entry to the built-in `WignerD` (index order $\{j, m', m\}$). The built-in uses the
opposite (active, $e^{+iJ}$) convention, so its Euler angles are negated relative to the passive $e^{-iJ}$
formula:

```wl
dmat == Table[WignerD[{j, mp, mm}, -\[Alpha], -\[Beta], -\[Gamma]], {mp, j, -j, -1}, {mm, j, -j, -1}] // Simplify
```

**QF** : assemble the same rotation as a product of `QuantumOperator` exponentials:

```wl
dq = QuantumOperator[MatrixExp[-I \[Alpha] jz]] @ QuantumOperator[MatrixExp[-I \[Beta] jy]] @ QuantumOperator[MatrixExp[-I \[Gamma] jz]];
```

QuantumFramework carries the $D$-matrix as the named operator `"WignerD"`, taking the spin and the three
Euler angles as a list, so the check is a direct operator equality (the framework reconciles the index
orderings itself):

```wl
dq == QuantumOperator["WignerD"[j, {\[Alpha], \[Beta], \[Gamma]}]] // Simplify
```

Both confirm the constructed rotation is the Wigner $D$-matrix: the WL side matches the built-in `WignerD`
table, the QF side equals the framework's named `"WignerD"` operator. The $D$-matrix is nothing more than
the explicit representation of the abstract rotation on the $(2j+1)$-dimensional spin space.

### 8.5 [MSc] How do I add two spin-1/2 and read off the full table of Clebsch-Gordan coefficients (the change of basis to the total-spin basis)?

Coupling two spins combines their spaces, $\tfrac12\otimes\tfrac12 = 0\oplus1$: the product basis
$|j_1\,m_1\rangle|j_2\,m_2\rangle$ recombines into total-spin eigenstates $|J,M\rangle$, tied together by
one closed-form expansion,

$$|J,M\rangle \;=\; \sum_{m_1+m_2=M} \langle j_1\,m_1;\, j_2\,m_2 \,|\, J\,M\rangle\;|j_1\,m_1\rangle|j_2\,m_2\rangle,$$

whose coefficients are the Clebsch-Gordan coefficients. The sum runs only over $m_1+m_2 = M$ (the total
$S_z$ is conserved), and they vanish unless $|j_1-j_2|\le J\le j_1+j_2$ (the triangle rule). The full table
of overlaps $\langle j_1\,m_1;\, j_2\,m_2|J,M\rangle$ *is* the change-of-basis matrix.

**WL** : the product basis is every pair of single-spin projections $m_i \in \{+\tfrac12, -\tfrac12\}$
(descending, matching the $J_z$ convention); the coupled states run over $J$ from $j_1+j_2$ down to
$|j_1-j_2|$ and $M$ from $J$ to $-J$:

```wl
prodBasis = Tuples[Range[1/2, -1/2, -1], 2]; coupled = Flatten[Table[{J, M}, {J, 1, 0, -1}, {M, J, -J, -1}], 1];
```

Each entry is $\langle j_1 m_1; j_2 m_2|J,M\rangle$ from the built-in `ClebschGordan`, nonzero only when
$M = m_1+m_2$:

```wl
cg = Outer[If[#2[[1]] + #2[[2]] == #1[[2]], ClebschGordan[{1/2, #2[[1]]}, {1/2, #2[[2]]}, #1], 0] &, coupled, prodBasis, 1];
```

Read the full table with its $(J,M)$ and $(m_1,m_2)$ labels:

```wl
TableForm[cg, TableHeadings -> {coupled, prodBasis}]
```

**QF** : the same coupling built genuinely as operators. Each total-spin component adds a named Pauli on
qubit 1 and qubit 2 (the order argument places each), scaled by $\tfrac12$,
$\vec S = \tfrac12(\vec\sigma_1 + \vec\sigma_2)$:

```wl
Sq = Table[1/2 (QuantumOperator[p, {1}] + QuantumOperator[p, {2}]), {p, {"X", "Y", "Z"}}];
```

The Casimir $S^2$ is the sum of their squares (`@`):

```wl
S2 = Sq[[1]] @ Sq[[1]] + Sq[[2]] @ Sq[[2]] + Sq[[3]] @ Sq[[3]];
```

Its spectrum is the coupling $\tfrac12\otimes\tfrac12 = 0\oplus1$ directly, one singlet ($S^2{=}0$) and a
threefold triplet ($S^2{=}2$):

```wl
Sort[S2["Eigenvalues"]]
```

The Clebsch-Gordan table is exactly the basis that diagonalizes this total spin. Take it as a
change-of-basis `QuantumOperator` `w` on the two qubits:

```wl
w = QuantumOperator[cg, {1, 2}];
```

Conjugating $S^2$ by `w` gives a diagonal reading $S(S+1)$ for each coupled state:

```wl
MatrixForm[Simplify@Normal[(w @ S2 @ w["Dagger"])["Matrix"]]]
```

The table is the unitary change of basis from the product basis to the total-spin basis, and conjugating
$S^2$ by it labels every coupled state by $S(S+1)$. The Clebsch-Gordan coefficients are the entries of the
matrix that adds the two spins.

### 8.6 [MSc] How do I build the singlet and triplet states and identify them by total spin?

Two spin-1/2 combine into one singlet ($S=0$, antisymmetric) and three triplet states ($S=1$, symmetric),
told apart by the total-spin Casimir $S^2 = S(S+1)$: $0$ for the singlet, $2$ for the triplet. Build all
four from the total-spin operators and identify each by its $(S,M)$.

**WL** : form the total spin $\vec S = \vec S_1 + \vec S_2$, then the triplet is the top state
$|{\uparrow\uparrow}\rangle$ lowered twice by $S_- = S_x - iS_y$, and the singlet is the one state the ladder
misses, the null space of $S^2$:

```wl
sk = Table[KroneckerProduct[PauliMatrix[k]/2, IdentityMatrix[2]] + KroneckerProduct[IdentityMatrix[2], PauliMatrix[k]/2], {k, 3}]; s2 = Total[#.# & /@ sk]; sMinus = sk[[1]] - I sk[[2]];
```

```wl
states = Append[NestList[Normalize[sMinus . #] &, UnitVector[4, 1], 2], Normalize[First[NullSpace[s2]]]]
```

Each state's total-spin quantum numbers are its Casimir and $S_z$ expectations, $(\langle S^2\rangle,
\langle S_z\rangle) = (S(S+1), M)$:

```wl
Simplify[{Conjugate[#] . s2 . #, Conjugate[#] . sk[[3]] . #} & /@ states]
```

**QF** : build the total spin genuinely from named Paulis, its components and Casimir:

```wl
Sq = Table[1/2 (QuantumOperator[p, {1}] + QuantumOperator[p, {2}]), {p, {"X", "Y", "Z"}}];
```

```wl
S2q = Sq[[1]] @ Sq[[1]] + Sq[[2]] @ Sq[[2]] + Sq[[3]] @ Sq[[3]];
```

Wrap the same four states as `QuantumState`s and read each $(\langle S^2\rangle, \langle S_z\rangle)$ from
the operators:

```wl
Simplify[{(#["Dagger"] @ S2q @ #)["Scalar"], (#["Dagger"] @ Sq[[3]] @ #)["Scalar"]} & /@ (QuantumState /@ states)]
```

The three ladder states share $S^2=2$ (so $S=1$) with $M=1,0,-1$, and the null-space state has $S^2=0$ (so
$S=0$) with $M=0$: singlet and triplet emerge from the construction already sorted by total spin, the
$1\oplus3$ multiplicity of $\tfrac12\otimes\tfrac12$. The singlet is fixed only up to an overall sign, as any
state is up to a global phase. Total spin is the conserved label that block-diagonalizes any rotationally
invariant two-spin interaction (Part 7, 7.10).

### 8.7 [MSc] How do I construct a spin coherent state and place it on the generalized Bloch sphere?

A spin coherent state is the most classical state of a spin-$j$: the highest-weight state $|j,j\rangle$
rotated to point along $\hat n(\theta,\varphi)$, $|j,\hat n\rangle = e^{-i\varphi J_z}e^{-i\theta
J_y}|j,j\rangle$. Its expectation vector is fully polarized, $\langle\vec J\rangle = j\,\hat n$, so it sits
on the sphere of radius $j$, the generalized Bloch sphere. Build it for spin $j=1$.

**WL** : build the spin-1 operators:

```wl
j = 1; m = Range[j, -j, -1]; jp = SparseArray[{r_, c_} /; r == c - 1 :> Sqrt[j (j + 1) - m[[c]] (m[[c]] + 1)], {2 j + 1, 2 j + 1}]; {jx, jy, jz} = {(jp + Transpose[jp])/2, (jp - Transpose[jp])/(2 I), DiagonalMatrix[m]};
```

Rotate the highest-weight state to point along $\hat n$:

```wl
coh = MatrixExp[-I \[Phi] jz] . MatrixExp[-I \[Theta] jy] . UnitVector[2 j + 1, 1];
```

and read its polarization, the three spin expectations:

```wl
Simplify[ComplexExpand[Conjugate[coh] . # . coh & /@ {jx, jy, jz}], {\[Theta], \[Phi]} \[Element] Reals]
```

**QF** : build the same coherent state with the framework's own rotations, exponentials of the named
`"JZ"[j]` and `"JY"[j]` on the top state, and read its polarization from the named operators:

```wl
cohq = Exp[-I \[Phi] QuantumOperator["JZ"[j]]][Exp[-I \[Theta] QuantumOperator["JY"[j]]][QuantumState["0", 2 j + 1]]];
```

```wl
Simplify[ComplexExpand[Table[(cohq["Dagger"] @ QuantumOperator[op[j]] @ cohq)["Scalar"], {op, {"JX", "JY", "JZ"}}]], {\[Theta], \[Phi]} \[Element] Reals]
```

Both give $\langle\vec J\rangle = (\sin\theta\cos\varphi, \sin\theta\sin\varphi, \cos\theta) = j\,\hat n$
with $j=1$: the spin coherent state points as definitely as a spin can along $\hat n$, saturating the
angular-momentum uncertainty and serving as the finite-dimensional stand-in for the coherent state of
light.
