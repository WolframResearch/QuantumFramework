# Quantum in Discrete Space: A Comprehensive Exam in the Wolfram Language and QuantumFramework

This manual is a self-study comprehensive exam in finite-dimensional quantum theory, the
quantum mechanics of "discrete space": qubits and qudits, observables and measurement,
entanglement, mixed states and channels, circuits and algorithms, the stabilizer formalism,
open-system dynamics, and a short coda on the bosonic (Fock) limit. It is built around a single
conviction: if you can compute a thing, you understand it. Every question is answered with the
shortest, clearest executable computation, almost always plain Wolfram Language built directly on
matrices and vectors, reaching for the QuantumFramework paclet where it is the right tool.

Work through it with a kernel open. Each item states a physics question, gives a compact concept
gloss, shows the code, and closes with one line on what it computes. We start from the bare qubit
and the Born rule, build up operators, measurement, and entanglement, pass through circuits and
the stabilizer formalism, reach open quantum systems and the standard algorithms, and finish in
the harmonic oscillator. The exam is long by design; treat each lettered part as one sitting.

## How to use this manual

The guiding rule is minimality: each answer is the most direct computation that settles the
question. For a single qubit, its operators, its density matrix, and its open-system dynamics,
that is plain Wolfram Language built on `PauliMatrix`: states are ordinary vectors, operators
ordinary matrices, a bra is `Conjugate[v]`, an inner product `Conjugate[u] . v`, an expectation
value `Conjugate[\[Psi]] . A . \[Psi]` or `Tr[A . \[Rho]]`, and a claim is checked by `==` under
`FullSimplify`. For the genuinely many-body or specialized objects, multi-qubit stabilizer
circuits, qudit clock-and-shift algebra, named quantum algorithms, channels, discrete phase
space, and Fock space, the most direct route is the QuantumFramework paclet, which we load once:

```wolfram
Needs["Wolfram`QuantumFramework`"]
```

A small toolkit, reused throughout, sets up the Pauli vector and the two maps between a qubit and
its Bloch picture (the Bloch vector $\vec r = (\langle\sigma_x\rangle, \langle\sigma_y\rangle,
\langle\sigma_z\rangle)$, and the density matrix $\rho = \tfrac12(I + \vec r\cdot\vec\sigma)$):

```wolfram
\[Sigma] = Table[PauliMatrix[j], {j, 3}];                         (* Pauli vector {X, Y, Z} *)
blochVector[\[Psi]_?VectorQ] := Re@Table[Conjugate[\[Psi]] . s . \[Psi], {s, \[Sigma]}];
blochVector[\[Rho]_?MatrixQ] := Re@Table[Tr[\[Rho] . s], {s, \[Sigma]}];
densityMatrix[r_] := 1/2 (IdentityMatrix[2] + r . \[Sigma]);
```

All code is paste-safe plain text. Greek letters are written in their Wolfram input form
(`\[Sigma]`, `\[Theta]`, `\[Psi]`) so the cells paste straight into a notebook.

---

## Part A. The discrete state space and the qubit

A qubit is a unit vector in $\mathbb{C}^2$, taken up to an overall phase. Everything physical,
the probabilities and the expectation values, is quadratic in that vector, so an overall phase
$e^{i\gamma}$ never shows up in a measurement. This part fixes notation by computing with it.

#### A1. How do I represent a qubit and read its measurement probabilities?

A pure qubit state is $|\psi\rangle = \alpha|0\rangle + \beta|1\rangle$ with
$|\alpha|^2 + |\beta|^2 = 1$. The Born rule says the probability of outcome $i$ in the
computational basis is $p_i = |\langle i|\psi\rangle|^2 = |\,\text{amplitude}_i\,|^2$. Construct a
state from its amplitudes and ask for both the vector and the probabilities.

```wolfram
\[Psi] = {Cos[\[Pi]/6], Sin[\[Pi]/6] Exp[I \[Pi]/4]};
\[Psi]
Abs[\[Psi]]^2
```

The amplitudes are $\tfrac{\sqrt3}{2}$ and $\tfrac12 e^{i\pi/4}$; squaring their moduli gives
$p_0 = \tfrac34$, $p_1 = \tfrac14$, exactly the Born rule, and they sum to one.

#### A2. How do I check normalization, and confirm that a global phase is physically invisible?

A legitimate state vector has unit norm, and two vectors differing only by an overall phase
$e^{i\gamma}$ describe the same physical state: they have identical density matrices and zero
distinguishing distance. Verify both facts.

```wolfram
Norm[\[Psi]]
Normalize[{1, 1}]
FullSimplify[Abs[{\[Alpha], \[Beta]}]^2 == Abs[Exp[I \[CurlyPhi]] {\[Alpha], \[Beta]}]^2,
  \[CurlyPhi] \[Element] Reals]
```

The norm is one; `Normalize` rescales any nonzero vector to unit length; and the squared
amplitudes, hence every Born probability, are identical for $\{\alpha,\beta\}$ and
$e^{i\varphi}\{\alpha,\beta\}$. Relative phases between amplitudes are physical, but the overall
phase is pure gauge.

#### A3. How do I get the Bloch vector of a qubit?

Any qubit density matrix is $\rho = \tfrac12(I + \vec r \cdot \vec\sigma)$ with
$\vec r = (\langle X\rangle, \langle Y\rangle, \langle Z\rangle)$, the Bloch vector. Pure states
sit on the unit sphere $|\vec r| = 1$. Read it off directly.

```wolfram
blochVector[{1, 1}/Sqrt[2]]
blochVector[{1, I}/Sqrt[2]]
```

The state $|+\rangle = (|0\rangle+|1\rangle)/\sqrt2$ points along $+x$ and the state
$|R\rangle = (|0\rangle+i|1\rangle)/\sqrt2$ (right-circular) along $+y$, as the eigenstate
structure of $X$ and $Y$ demands.

#### A4. How do I prepare a qubit from Bloch angles and recover them?

The map from sphere to state is
$|\psi\rangle = \cos\tfrac\theta2|0\rangle + e^{i\varphi}\sin\tfrac\theta2|1\rangle$. Map this
general ket to its Bloch vector and read off the angles.

```wolfram
FullSimplify[blochVector[{Cos[\[Theta]/2], Exp[I \[Phi]] Sin[\[Theta]/2]}],
  {\[Theta], \[Phi]} \[Element] Reals]
```

The Bloch vector comes out as $(\sin\theta\cos\varphi, \sin\theta\sin\varphi, \cos\theta)$: the
ket's parameters $\theta, \varphi$ are exactly the spherical-coordinate angles on the Bloch
sphere, which is why the half-angle $\theta/2$ appears in the amplitude.

#### A5. How do I compute an overlap $\langle\phi|\psi\rangle$ and the transition probability?

The inner product $\langle\phi|\psi\rangle$ is the conjugated dot product
$\mathrm{Conjugate}[\phi]\cdot\psi$. The probability of finding $|\psi\rangle$ to pass a test for
$|\phi\rangle$ is $|\langle\phi|\psi\rangle|^2$.

```wolfram
amp = Conjugate[{1, 0}] . ({1, 1}/Sqrt[2])
Abs[amp]^2
```

The overlap $\langle 0|+\rangle = 1/\sqrt2$ gives a transition probability of $\tfrac12$: a state
prepared as $|+\rangle$ and measured in the computational basis yields $0$ half the time.

#### A6. How do I prepare a qudit and a uniform superposition?

Discrete space is not only qubits. A qudit lives in $\mathbb{C}^d$, so it is just a longer unit
vector. The uniform superposition $\tfrac{1}{\sqrt d}\sum_k |k\rangle$ is the natural "flat"
state, the starting point of most algorithms.

```wolfram
UnitVector[3, 1]                  (* qutrit |0> in C^3 *)
Normalize@ConstantArray[1, 3]     (* flat qutrit *)
Normalize@ConstantArray[1, 4]     (* flat two-qubit register in C^4 *)
```

The first is the qutrit $|0\rangle$ in $\mathbb{C}^3$; the second the flat qutrit
$\tfrac{1}{\sqrt3}\sum_k|k\rangle$; the third the flat two-qubit register in $\mathbb{C}^4$.
Nothing about the formalism cared that $d = 2$.

**Summary of Part A.** We can build a qubit from amplitudes or Bloch angles, read its
probabilities (Born rule) and Bloch vector, confirm that normalization holds and global phase is
invisible, compute overlaps and transition probabilities, and lift all of this to qudits in
$\mathbb{C}^d$. These are the nouns; next we build the verbs.

---

## Part B. Observables and single-qubit operators

Observables are Hermitian operators; their real eigenvalues are the possible measurement
outcomes and their eigenvectors the states that yield those outcomes with certainty. Gates are
unitary operators; they rotate the state without changing its norm. The single qubit is where
both ideas are cleanest, because every $2\times2$ Hermitian operator is a real combination of
$I, X, Y, Z$ and every single-qubit gate is a rotation of the Bloch sphere.

#### B1. How do I build the Pauli operators and verify their algebra?

The whole single-qubit algebra is one product rule,
$\sigma_i\sigma_j = \delta_{ij}I + i\,\epsilon_{ijk}\sigma_k$. Pack the Paulis into a vector,
build the table of products with `Outer`, and check it against the Levi-Civita contraction.

```wolfram
\[Sigma] = Table[PauliMatrix[j], {j, 3}];
(* product rule: \[Sigma]i.\[Sigma]j = \[Delta]ij I + i \[Epsilon]ijk \[Sigma]k *)
Outer[Dot, \[Sigma], \[Sigma], 1] ==
  TensorProduct[IdentityMatrix[3], IdentityMatrix[2]] +
   I TensorContract[TensorProduct[LeviCivitaTensor[3], \[Sigma]], {{3, 4}}]
(* commutator: [\[Sigma]i, \[Sigma]j] = 2 i \[Epsilon]ijk \[Sigma]k *)
Outer[Commutator[#1, #2, Dot] &, \[Sigma], \[Sigma], 1] ==
  2 I TensorContract[TensorProduct[LeviCivitaTensor[3], \[Sigma]], {{3, 4}}]
```

Both return `True`. The product rule's antisymmetric part is the commutator
$[\sigma_i,\sigma_j] = 2i\,\epsilon_{ijk}\sigma_k$ (the $SU(2)$ structure, $X^2=I$ on the diagonal,
$XY=iZ$ off it), its symmetric part the anticommutator $\{\sigma_i,\sigma_j\} = 2\delta_{ij}I$.
The Levi-Civita tensor is the cross-product table, so the qubit algebra is literally 3D geometry.

#### B2. How do I compute an expectation value $\langle\psi|A|\psi\rangle$, two ways?

The expectation value of an observable $A$ in state $|\psi\rangle$ is
$\langle A\rangle = \langle\psi|A|\psi\rangle = \mathrm{Tr}(A\rho)$ with $\rho = |\psi\rangle\langle\psi|$.
Take $|\psi\rangle = \tfrac{\sqrt3}{2}|0\rangle + \tfrac12|1\rangle$, for which $\langle X\rangle$
should be $2\,\mathrm{Re}(\alpha^*\beta) = \tfrac{\sqrt3}{2}$.

```wolfram
\[Psi] = {Sqrt[3]/2, 1/2};
Conjugate[\[Psi]] . PauliMatrix[1] . \[Psi]
Tr[PauliMatrix[1] . TensorProduct[\[Psi], Conjugate@\[Psi]]]
blochVector[\[Psi]]
```

Both forms give $\langle X\rangle = \tfrac{\sqrt3}{2}$: the bra-operator-ket contraction
$\langle\psi|X|\psi\rangle$ and the trace against the density matrix $|\psi\rangle\langle\psi|$.
And it is exactly the first component of the Bloch vector $(\tfrac{\sqrt3}{2}, 0, \tfrac12)$: for a
qubit, $\langle\sigma_j\rangle = r_j$. The Bloch vector is just the triple of Pauli expectations.

#### B3. How do I get the eigenvalues of an observable and the outcome probabilities?

Measuring $A$ projects onto its eigenbasis; the outcomes are the eigenvalues and their
probabilities are the squared overlaps with the eigenvectors. For $Y$ the eigenvalues are
$\pm1$, and on $|0\rangle$ each occurs with probability $\tfrac12$ (since $|0\rangle$ is an equal
superposition of the $Y$ eigenstates).

```wolfram
Eigenvalues[PauliMatrix[2]]
Abs[Conjugate[Normalize[#]] . {1, 0}]^2 & /@ Eigenvectors[PauliMatrix[2]]
```

The eigenvalues are $\{-1, 1\}$ and the outcome distribution on $|0\rangle$, the squared overlaps
$|\langle y_\pm|0\rangle|^2$, is flat: the hallmark of measuring in a basis mutually unbiased to
the one you prepared in.

#### B4. How do I build single-qubit rotations and the closed-form exponential?

The rotation $R_n(\theta) = e^{-i\theta\, \hat n\cdot\vec\sigma/2}$ turns the Bloch sphere about
axis $\hat n$ by angle $\theta$. Because $(\hat n\cdot\vec\sigma)^2 = I$ for a unit vector, the
exponential collapses to a closed cosine-sine form.

```wolfram
MatrixExp[-I \[Theta]/2 PauliMatrix[1]]    (* Rx(\[Theta]) *)
MatrixExp[-I \[Theta]/2 PauliMatrix[3]]    (* Rz(\[Theta]) *)
With[{n = {Cos[\[Phi]] Sin[\[Theta]], Sin[\[Phi]] Sin[\[Theta]], Cos[\[Theta]]}},
 FullSimplify[MatrixExp[-I \[Alpha]/2 n . \[Sigma]] ==
   Cos[\[Alpha]/2] IdentityMatrix[2] - I Sin[\[Alpha]/2] n . \[Sigma],
  (\[Theta] | \[Phi] | \[Alpha]) \[Element] Reals]]
```

$R_x$ mixes the amplitudes with an imaginary off-diagonal; $R_z$ is a relative phase
$e^{\pm i\theta/2}$; and the general axis-angle rotation is exactly
$\cos\tfrac\alpha2\,I - i\sin\tfrac\alpha2\,\hat n\cdot\vec\sigma$. Every single-qubit gate is one
of these, up to a global phase.

#### B5. How do I check that an operator is unitary or Hermitian (and avoid the symbolic trap)?

`UnitaryMatrixQ` applies a numeric test, so on a symbolic gate it can wrongly report `False`
because it cannot see that $\sin^2+\cos^2=1$ under an unstated reality assumption. The reliable
symbolic check is to simplify $U U^\dagger$ to the identity with the angle declared real.

```wolfram
ry = MatrixExp[-I \[Theta]/2 PauliMatrix[2]];
UnitaryMatrixQ[ry]
FullSimplify[ry . ConjugateTranspose[ry] == IdentityMatrix[2], \[Theta] \[Element] Reals]
HermitianMatrixQ[PauliMatrix[1]]
```

The naive `UnitaryMatrixQ` says `False` on the symbolic $R_y(\theta)$, but the proper check, with
$\theta\in\mathbb{R}$, confirms unitarity; and $X$ is correctly Hermitian. The lesson is to supply
reality assumptions when verifying symbolic operators.

#### B6. How do I decompose a single-qubit unitary into Euler ($Z$-$Y$-$Z$) angles?

Any single-qubit unitary is a Bloch-sphere rotation, and any rotation factors into three
fixed-axis turns, the Euler ($Z$-$Y$-$Z$) decomposition $R_z(\alpha)R_y(\beta)R_z(\gamma)$.
`EulerAngles` returns the angles; the $SU(2)$ product of half-angle exponentials reproduces the
unitary. Take the rotation by $\pi/3$ about $\hat n = (\cos\tfrac\pi8, \sin\tfrac\pi8, 0)$.

```wolfram
{\[Alpha], \[Beta], \[Gamma]} = EulerAngles[RotationMatrix[\[Pi]/3, {Cos[\[Pi]/8], Sin[\[Pi]/8], 0}]];
{\[Alpha], \[Beta], \[Gamma]}
FullSimplify[
 MatrixExp[-I \[Alpha]/2 PauliMatrix[3]] . MatrixExp[-I \[Beta]/2 PauliMatrix[2]] .
   MatrixExp[-I \[Gamma]/2 PauliMatrix[3]] ==
  MatrixExp[-I (\[Pi]/3)/2 {Cos[\[Pi]/8], Sin[\[Pi]/8], 0} . \[Sigma]]]
```

The middle angle $\beta = \pi/3$ is the rotation angle, and the $Z$-$Y$-$Z$ product of half-angle
exponentials equals the axis-angle unitary exactly. This is how an arbitrary gate is compiled into
the fixed-axis pulses that hardware can actually perform.

#### B7. How do I compute a commutator and read off an uncertainty relation?

The commutator $[A,B] = AB - BA$ measures incompatibility; for the Paulis $[X,Y] = 2iZ$. The
Robertson uncertainty relation $\Delta A\,\Delta B \ge \tfrac12|\langle[A,B]\rangle|$ then bounds
the spreads of $X$ and $Y$. Check the commutator, then the bound on $|0\rangle$.

```wolfram
Commutator[PauliMatrix[1], PauliMatrix[2], Dot] == 2 I PauliMatrix[3]
\[Psi] = {1, 0};
var[A_] := Conjugate[\[Psi]] . A . A . \[Psi] - (Conjugate[\[Psi]] . A . \[Psi])^2;
{Sqrt[var[PauliMatrix[1]] var[PauliMatrix[2]]], Abs[Conjugate[\[Psi]] . PauliMatrix[3] . \[Psi]]}
```

The commutator is $2iZ$, and on $|0\rangle$ the product of spreads $\Delta X\,\Delta Y = 1$ exactly
meets the bound $\tfrac12|\langle[X,Y]\rangle| = |\langle Z\rangle| = 1$: the uncertainty relation
is saturated. Incompatible observables cannot both be sharp.

**Summary of Part B.** Observables are Hermitian, gates are unitary, and on one qubit both are
combinations of Paulis. We encoded the whole algebra in one Levi-Civita product rule, computed
expectation values two equivalent ways and saw they are the Bloch components, built rotations from
the closed-form exponential and decomposed one by Euler angles, tested unitarity symbolically, and
saturated an uncertainty bound.

---

## Part C. Measurement

Measurement is where the linear, reversible world of states and gates meets probability. A
projective measurement of an observable returns one eigenvalue, with the Born probability, and
collapses the state onto the corresponding eigenspace. A generalized (POVM) measurement relaxes
the orthogonality of the outcomes into positive effects $\{E_k\}$. Every probability is the trace
rule $\mathrm{Tr}(E_k\rho)$, and a projective outcome collapses the state onto its eigenspace.

#### C1. How do I perform a projective measurement and read outcomes and probabilities?

The probability of outcome $k$ is $\mathrm{Tr}(P_k\rho)$ with $P_k$ the eigenprojector, the
outcomes are the eigenvalues, and the mean is $\mathrm{Tr}(Z\rho)$.

```wolfram
\[Rho] = densityMatrix[{1, 0, 0}];                           (* |+> *)
proj = TensorProduct[#, Conjugate@#] & /@ {{1, 0}, {0, 1}};   (* {|0><0|, |1><1|} *)
Tr[# . \[Rho]] & /@ proj                                     (* outcome probabilities *)
Eigenvalues[PauliMatrix[3]]                                  (* the outcomes \[PlusMinus]1 *)
Tr[PauliMatrix[3] . \[Rho]]                                  (* mean = <Z> *)
```

Measuring $Z$ on $|+\rangle$ gives $\pm1$ each with probability $\tfrac12$ and mean zero, exactly
$\langle +|Z|+\rangle = 0$. The probabilities come cleanly as a list; the eigenvalues come from
the observable.

#### C2. How do I measure in a basis other than the computational one?

To measure in the $X$ basis, take the squared overlaps of $|0\rangle$ with the $X$ eigenstates
$|x_\pm\rangle = (|0\rangle\pm|1\rangle)/\sqrt2$. Measuring $X$ on $|0\rangle$ is maximally uncertain.

```wolfram
xpm = {{1, 1}, {1, -1}}/Sqrt[2];        (* |x+>, |x-> as rows *)
Abs[Conjugate[xpm] . {1, 0}]^2
```

Because $|0\rangle$ is an equal superposition of $|+\rangle$ and $|-\rangle$, the $X$-basis
outcomes are equally likely: a flat distribution is the signature of mutually unbiased bases.

#### C3. How do I get the post-measurement (collapsed) states?

After outcome $i$ the state collapses to the normalized projection $P_i|\psi\rangle/\sqrt{p_i}$.
Applying the projectors to $|+\rangle$ gives the (sub-normalized) collapsed branches, whose squared
norms are the outcome probabilities.

```wolfram
\[Psi] = {1, 1}/Sqrt[2];
P0 = TensorProduct[{1, 0}, {1, 0}]; P1 = TensorProduct[{0, 1}, {0, 1}];
{P0 . \[Psi], P1 . \[Psi]}
```

The two branches are $|0\rangle$ and $|1\rangle$, each carrying weight $1/\sqrt2$ so that
$\||{\cdot}\rangle\|^2 = \tfrac12$ reproduces the Born probability. Measuring $Z$ on $|+\rangle$
collapses it to a computational eigenstate, as expected.

#### C4. How do I build and apply a POVM (a SIC measurement)?

A positive operator-valued measure replaces orthogonal projectors with positive operators
$\{E_k\}$ summing to the identity. The qubit SIC-POVM has four effects $E_k = \tfrac14(I + \hat
b_k\cdot\vec\sigma)$ whose Bloch vectors $\hat b_k$ point to the vertices of a tetrahedron. Build
them and check positivity, completeness, and that the probabilities on a state sum to one.

```wolfram
verts = {{1, 1, 1}, {1, -1, -1}, {-1, 1, -1}, {-1, -1, 1}}/Sqrt[3];
effects = (1/4) (IdentityMatrix[2] + # . \[Sigma]) & /@ verts;
And @@ (AllTrue[Eigenvalues[#], # >= 0 &] & /@ effects)   (* each effect PSD *)
Total[effects] == IdentityMatrix[2]                        (* completeness *)
Simplify@Total[Tr[# . densityMatrix[{0, 0, 1}]] & /@ effects]   (* probabilities sum to 1 *)
```

The four tetrahedral effects are positive semidefinite, sum to the identity, and give a normalized
probability distribution on $|0\rangle$. A single such informationally complete measurement, on
many copies, determines the full qubit density matrix, something no projective scheme can do.

#### C5. How do I simulate finite-shot statistics and watch the Born rule emerge?

Real measurements are finite samples from the Born distribution. `RandomChoice` with the Born
weights draws $n$ shots; with $2000$ shots the empirical frequencies should hug the theoretical
$\{\tfrac12,\tfrac12\}$.

```wolfram
SeedRandom[42];
Counts@RandomChoice[{1/2, 1/2} -> {0, 1}, 2000]
N@{1/2, 1/2}
```

The $2000$ shots split close to $1000{:}1000$, within the $\sqrt{N}$ fluctuation of the exact
distribution. This is the law of large numbers turning amplitudes into observed frequencies.

#### C6. How do I measure an entangled pair and see the correlations?

Measuring both qubits of a Bell state $|\Phi^+\rangle = (|00\rangle+|11\rangle)/\sqrt2$ in the
computational basis is just reading the squared amplitudes over $\{00, 01, 10, 11\}$.

```wolfram
bell = {1, 0, 0, 1}/Sqrt[2];
Abs[bell]^2
```

The joint distribution puts weight only on $00$ and $11$, each $\tfrac12$: the two qubits are
random individually but perfectly correlated together. This is the raw material of the EPR
argument and of every entanglement-based protocol to come.

**Summary of Part C.** A measurement is the trace rule $\mathrm{Tr}(E_k\rho)$ against projectors
or POVM effects. We read projective probabilities in arbitrary bases, collapsed a state with its
projectors, built a tetrahedral SIC-POVM, sampled finite-shot statistics that converge to the Born
rule, and exhibited the perfect correlations of a measured Bell pair.

---

## Part D. Composite systems, tensor products, and entanglement

Two systems combine by the tensor product, and the dimension multiplies: two qubits live in
$\mathbb{C}^4$. The states that do not factor into a product are entangled, and entanglement is
the genuinely quantum resource behind teleportation, error correction, and Bell nonlocality. This
part builds composite states and measures their entanglement quantitatively.

A toolkit for bipartite (two-qubit) work: the reduced density matrices (trace out one qubit), the
partial transpose, and the von Neumann entropy in bits.

```wolfram
reduceA[\[Rho]_] := TensorContract[ArrayReshape[\[Rho], {2, 2, 2, 2}], {{2, 4}}];
reduceB[\[Rho]_] := TensorContract[ArrayReshape[\[Rho], {2, 2, 2, 2}], {{1, 3}}];
partialT[\[Rho]_] := ArrayReshape[Transpose[ArrayReshape[\[Rho], {2, 2, 2, 2}], {1, 4, 3, 2}], {4, 4}];
vne[\[Rho]_] := -Total[# Log2[#] & /@ Select[Re@Eigenvalues[\[Rho]], # > 10^-12 &]];
```

#### D1. How do I form the tensor product of states and of operators?

The tensor product is the Kronecker product: `Flatten@KroneckerProduct` for two kets (the
amplitudes multiply), `KroneckerProduct` for two operators (the matrices do).

```wolfram
Flatten@KroneckerProduct[{1, 0}, {0, 1}]          (* |0> \[CircleTimes] |1> *)
KroneckerProduct[PauliMatrix[1], PauliMatrix[3]]  (* X \[CircleTimes] Z *)
```

The joint state is $|01\rangle$, the second basis vector of $\mathbb{C}^4$, and $X \otimes Z$
is the $4\times4$ Kronecker product. The ordering convention is left qubit first.

#### D2. How do I build the standard entangled states (Bell, GHZ, W)?

Three families recur throughout quantum information: the Bell states (two-qubit maximally
entangled), the GHZ states (all-or-nothing $n$-party entanglement), and the W states (robust
single-excitation entanglement). Write them as their amplitude vectors.

```wolfram
bell = {1, 0, 0, 1}/Sqrt[2];                (* |\[CapitalPhi]+> *)
ghz = {1, 0, 0, 0, 0, 0, 0, 1}/Sqrt[2];     (* 3-qubit GHZ *)
w = {0, 1, 1, 0, 1, 0, 0, 0}/Sqrt[3];       (* 3-qubit W *)
{bell, ghz, w}
```

$|\Phi^+\rangle$ superposes $|00\rangle$ and $|11\rangle$; the GHZ state superposes $|000\rangle$
and $|111\rangle$; the W state spreads one excitation symmetrically over $|001\rangle$,
$|010\rangle$, $|100\rangle$. GHZ and W are the two inequivalent classes of tripartite
entanglement.

#### D3. How do I take a partial trace and get a reduced density matrix?

Tracing out a subsystem gives the reduced state seen by the rest. For a maximally entangled pair
the reduction is maximally mixed: locally, each half looks completely random.

```wolfram
\[Rho]bell = KroneckerProduct[bell, Conjugate@bell];
reduceA[\[Rho]bell]
```

The reduction is $\tfrac12 I$: one qubit of a Bell pair (and likewise one qubit of a GHZ triple)
carries no local information. The information lives entirely in the correlations.

#### D4. How do I test whether a state is entangled?

A pure bipartite state is entangled exactly when its Schmidt rank exceeds one, that is, when the
reshaped amplitude matrix has more than one nonzero singular value. (For mixed states, use the
positive-partial-transpose test of E3.)

```wolfram
schmidtRank[\[Psi]_, {dA_, dB_}] := Length@Select[SingularValueList[ArrayReshape[\[Psi], {dA, dB}]], # > 10^-10 &];
schmidtRank[bell, {2, 2}]                                 (* 2 -> entangled *)
schmidtRank[Flatten@KroneckerProduct[{1, 0}, {0, 1}], {2, 2}]   (* 1 -> product *)
```

The Bell state has Schmidt rank two (entangled), the product state $|01\rangle$ rank one
(separable). The rank counts how many product terms the state genuinely needs.

#### D5. How do I quantify entanglement (entropy, concurrence, negativity)?

Three monotones: the entanglement entropy (von Neumann entropy of either reduction), the
concurrence $|\langle\psi|(\sigma_y\otimes\sigma_y)|\psi^*\rangle|$, and the negativity (sum of
the absolute negative eigenvalues of the partial transpose). All vanish on products and are
maximal on Bell states.

```wolfram
vne[reduceA[\[Rho]bell]]                                   (* entropy of entanglement, in bits *)
concurrence[\[Psi]_] := Abs[Conjugate[\[Psi]] . KroneckerProduct[PauliMatrix[2], PauliMatrix[2]] . Conjugate[\[Psi]]];
{concurrence[bell], concurrence[Flatten@KroneckerProduct[{1, 0}, {0, 1}]]}
-Total@Select[Re@Eigenvalues@partialT[\[Rho]bell], # < 0 &]   (* negativity *)
```

The Bell state carries one bit of entanglement entropy, concurrence $1$, negativity $\tfrac12$
(all maximal); the product state has concurrence $0$. These numbers are the currency in which
entanglement is spent and measured.

#### D6. How do I compute the Schmidt decomposition and Schmidt rank?

Any bipartite pure state can be written as $\sum_i \lambda_i |a_i\rangle|b_i\rangle$ with
$\lambda_i \ge 0$ (the Schmidt coefficients) and orthonormal local bases. The number of nonzero
$\lambda_i$ is the Schmidt rank; rank one means product, rank above one means entangled. The
coefficients are the singular values of the reshaped amplitude matrix.

```wolfram
SingularValueList[ArrayReshape[bell, {2, 2}]]
```

Two equal Schmidt coefficients $\tfrac{1}{\sqrt2}$: the Bell state has Schmidt rank two and is
maximally entangled, since the entanglement entropy is the Shannon entropy of the squared
coefficients, $-\sum_i \lambda_i^2 \log_2 \lambda_i^2 = 1$ bit.

#### D7. How do I verify no-signalling through reduced states?

A local operation on one half of an entangled pair cannot change the reduced state of the other
half; otherwise entanglement would transmit information faster than light. Apply a flip to qubit
$1$ and confirm qubit $2$'s reduction is untouched.

```wolfram
XI = KroneckerProduct[PauliMatrix[1], IdentityMatrix[2]];
reduceB[\[Rho]bell] == reduceB[XI . \[Rho]bell . ConjugateTranspose[XI]]
```

Qubit $2$'s reduced state is $\tfrac12 I$ both before and after the flip on qubit $1$: the local
density matrix is invariant under remote operations. Entanglement correlates but does not signal.

**Summary of Part D.** Composite systems multiply dimensions through the tensor product; the
non-factoring states are entangled. We built the Bell, GHZ, and W states, reduced them by partial
trace, tested separability, quantified entanglement by entropy, concurrence, and negativity, read
the Schmidt decomposition, and confirmed that entanglement respects no-signalling.

---

## Part E. Mixed states and density matrices

A density matrix $\rho$ is the most general state: a positive operator of unit trace describing
both quantum superposition and classical ignorance. Pure states satisfy $\rho^2 = \rho$;
everything else is mixed. The density matrix becomes unavoidable the moment a system is entangled
with an environment we cannot see, which is the bridge to open dynamics in Part H.

#### E1. How do I build a density matrix from an ensemble and measure its purity?

A classical mixture of preparations is the weighted sum of their density matrices. The purity
$\mathrm{Tr}(\rho^2)$ runs from $1$ (pure) down to $1/d$ (maximally mixed); it is the simplest
measure of how mixed a state is.

```wolfram
\[Rho] = 1/2 densityMatrix[{0, 0, 1}] + 1/2 densityMatrix[{1, 0, 0}];  (* 1/2|0><0| + 1/2|+><+| *)
\[Rho]
Tr[\[Rho] . \[Rho]]                                    (* purity *)
Tr[densityMatrix[{0, 0, 1}] . densityMatrix[{0, 0, 1}]]   (* pure -> 1 *)
Tr[(IdentityMatrix[2]/2) . (IdentityMatrix[2]/2)]      (* maximally mixed -> 1/2 *)
```

The equal mixture of $|0\rangle$ and $|+\rangle$ has purity $\tfrac34$, between the pure value
$1$ and the maximally mixed qubit value $\tfrac12$. Purity $\mathrm{Tr}(\rho^2) = \tfrac12(1+r^2)$
is a one-number summary of mixedness through the Bloch radius $r$.

#### E2. How do I compute the von Neumann entropy (and handle the large-system limit)?

The von Neumann entropy $S(\rho) = -\mathrm{Tr}(\rho \log_2 \rho)$ is zero for pure states and
maximal ($\log_2 d$) for the maximally mixed state. It is also the entanglement entropy of a pure
bipartite state through its reductions.

```wolfram
vne[densityMatrix[{0, 0, 1}]]      (* pure -> 0 *)
vne[IdentityMatrix[2]/2]           (* maximally mixed -> 1 bit *)
vne[reduceA[\[Rho]bell]]           (* Bell reduction -> 1 bit of entanglement entropy *)
```

A pure state has zero entropy, a maximally mixed qubit one bit, and the Bell reduction one bit of
entanglement entropy. The eigenvalue formula $S = -\sum_k \lambda_k \log_2 \lambda_k$ has no size
limit; for very large systems, reduce to a small subsystem first so the eigensolve stays cheap.

#### E3. How do I build a Werner state and find its entanglement threshold?

A Werner state mixes the singlet with the maximally mixed state,
$\rho_W(q) = q\,|\psi^-\rangle\langle\psi^-| + (1-q)\,I/4$. The partial-transpose eigenvalues
locate the boundary where entanglement disappears.

```wolfram
singlet = {0, 1, -1, 0}/Sqrt[2];
werner[q_] := q KroneckerProduct[singlet, singlet] + (1 - q) IdentityMatrix[4]/4;
Eigenvalues[partialT[werner[q]]] // Simplify       (* one eigenvalue is (1 - 3 q)/4 *)
Reduce[(1 - 3 q)/4 < 0 && 0 <= q <= 1, q]           (* entangled iff q > 1/3 *)
```

The partial transpose has eigenvalues $\{(1-3q)/4, (1+q)/4, (1+q)/4, (1+q)/4\}$; the first turns
negative exactly at $q = 1/3$. So the Werner state is entangled precisely for $q > 1/3$, the
Peres-Horodecki threshold, and separable below it: a sharp boundary between quantum and merely
classical correlations.

#### E4. How do I compare two states by fidelity and trace distance?

Two states are compared by their fidelity $F = |\langle\psi|\phi\rangle|^2$ (near $1$ for close
states) or their trace distance $T = \tfrac12\mathrm{Tr}|\rho-\sigma| = \tfrac12\sum_k
|\lambda_k(\rho-\sigma)|$ (near $0$). For pure states the two are tied by $T = \sqrt{1-F}$.

```wolfram
\[Psi]0 = {1, 0}; \[Psi]p = {1, 1}/Sqrt[2];
overlap = Conjugate[\[Psi]0] . \[Psi]p;
{Abs[overlap], Abs[overlap]^2}                       (* root fidelity, fidelity F *)
\[Rho]0 = TensorProduct[\[Psi]0, \[Psi]0]; \[Rho]p = TensorProduct[\[Psi]p, Conjugate@\[Psi]p];
1/2 Total@Abs@Eigenvalues[\[Rho]0 - \[Rho]p]          (* trace distance *)
```

The root fidelity $|\langle 0|+\rangle|$ is $\tfrac{1}{\sqrt2}$, the fidelity $F$ is $\tfrac12$,
and the trace distance is $\tfrac{1}{\sqrt2} = \sqrt{1-F}$. Fidelity and distance are two readings
of the same closeness.

#### E5. How do I purify a mixed state?

Every mixed state is the reduction of a pure state on a larger space (the Church of the larger
Hilbert space). The purification of $\rho = \sum_i\lambda_i|i\rangle\langle i|$ is
$|\Psi\rangle = \sum_i\sqrt{\lambda_i}\,|i\rangle_A|i\rangle_B$; tracing the ancilla back out
recovers $\rho$.

```wolfram
purify[\[Rho]_] := With[{e = Eigensystem[\[Rho]]},
   Total@MapThread[Sqrt[#1] Flatten@KroneckerProduct[#2, #2] &, e]];
purify[IdentityMatrix[2]/2]
reduceA[KroneckerProduct[purify[IdentityMatrix[2]/2], purify[IdentityMatrix[2]/2]]]
```

The maximally mixed qubit purifies to the Bell state $|\Phi^+\rangle$, whose reduction is again
$\tfrac12 I$. Mixedness of a part is entanglement with the whole: the two descriptions are the
same physics.

**Summary of Part E.** The density matrix unifies superposition and ignorance. We mixed ensembles
and measured purity, computed von Neumann entropy with no size limit, found the entanglement
threshold of a Werner family at $q = 1/3$, tied fidelity to trace distance, and purified a mixed
state into a Bell pair. Mixedness of a part is entanglement with the whole.

---

## Part F. Multi-qubit gates, circuits, and universality

Computation is a unitary on many qubits, built from a small repertoire of one- and two-qubit
gates. The decisive fact is that two-qubit gates plus single-qubit rotations are universal: any
unitary, however large, factors into them. This part builds the standard gates, assembles and
inspects circuits, and exhibits the decomposition that makes universality concrete.

#### F1. How do I build the standard multi-qubit gates?

A controlled gate is "apply the target gate when the controls are $|1\rangle$": the control
projector $P_1 = |1\rangle\langle1|$ tensored with the target gate, plus identity on the rest.
Build CNOT and Toffoli this way (and SWAP from the Heisenberg form $\tfrac12\sum_\mu \sigma_\mu
\otimes \sigma_\mu$), and read their action on basis states.

```wolfram
ket[bits__] := Fold[Flatten@KroneckerProduct[#1, #2] &, UnitVector[2, # + 1] & /@ {bits}];
P0 = {{1, 0}, {0, 0}}; P1 = {{0, 0}, {0, 1}};     (* |0><0|, |1><1| *)
cnot = KroneckerProduct[P0, IdentityMatrix[2]] + KroneckerProduct[P1, PauliMatrix[1]];
swap = 1/2 Sum[KroneckerProduct[s, s], {s, Table[PauliMatrix[j], {j, 0, 3}]}];
toffoli = KroneckerProduct[IdentityMatrix[4] - KroneckerProduct[P1, P1], IdentityMatrix[2]] +
   KroneckerProduct[P1, P1, PauliMatrix[1]];
cnot . ket[1, 0]          (* |10> -> |11> *)
swap . ket[1, 0]          (* |10> -> |01> *)
toffoli . ket[1, 1, 0]    (* |110> -> |111> *)
```

CNOT sends $|10\rangle \to |11\rangle$ (flip the target when the control is one); SWAP sends
$|10\rangle \to |01\rangle$; Toffoli sends $|110\rangle \to |111\rangle$ (flip only when both
controls are one). Toffoli alone is universal for reversible classical computation.

#### F2. How do I compose a circuit and apply it?

A circuit is a product of gate unitaries, applied right to left. The canonical example is
Bell-state preparation, a Hadamard on qubit $1$ then a CNOT.

```wolfram
H = {{1, 1}, {1, -1}}/Sqrt[2];
bellU = cnot . KroneckerProduct[H, IdentityMatrix[2]];   (* H on qubit 1, then CNOT *)
bellU . ket[0, 0]
```

The circuit turns $|00\rangle$ into $|\Phi^+\rangle$: entanglement is generated by a single
two-qubit gate acting after a local superposition. To draw the wire diagram instead of multiplying
matrices, QuantumFramework's `CircuitDraw` renders it.

#### F3. How do I get the depth and parallel layers of a circuit?

Gates touching disjoint qubits run in parallel; the depth is the number of such layers and sets
the runtime on hardware. `["Layers"]` partitions the circuit, `["Depth"]` counts the partitions.

```wolfram
qc = QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}, "H" -> 1, "X" -> 2}];
qc["Depth"]
Length[qc["Layers"]]
```

The four gates compress into three layers because the final $H$ on qubit $1$ and $X$ on qubit $2$
act on different wires and share a layer. Depth, not gate count, is the relevant cost.

#### F4. How do I get the overall unitary of a circuit?

The circuit's overall unitary is just the matrix product of its gates, which we already formed as
`bellU`. Displaying it shows the whole gate sequence collapsed into one $4\times4$ operator (this
costs $2^n\times2^n$ memory, so reserve it for small $n$).

```wolfram
bellU // MatrixForm
```

This is the Bell-preparation unitary $\mathrm{CNOT}\,(H \otimes I)$; its first column is
$|\Phi^+\rangle$, confirming the action on $|00\rangle$.

#### F5. How do I decompose a two-qubit gate into elementary gates (KAK)?

The Cartan (KAK) decomposition factors any two-qubit unitary into single-qubit gates and at most
three CNOTs. `["KAK"]` returns this as a circuit; tallying its gate labels shows the structure.

```wolfram
kak = QuantumOperator["CNOT"]["KAK"];
Head[kak]
Length[kak["Operators"]]
Tally[Replace[#["Label"],
     {"U"[__] -> "U", Subscript["C", "X"][__] -> "CX", "GlobalPhase"[_] -> "GlobalPhase"}] & /@
   kak["Operators"]]
```

The decomposition is a circuit of single-qubit $U$ gates and a global phase wrapped around two
controlled-$X$ gates: even CNOT is re-expressed in the canonical local-plus-entangler form (this
build emits two entanglers, within the KAK budget of three). Every two-qubit gate reduces to this
template of local rotations separated by controlled-$X$s.

#### F6. How do I confirm universality by decomposing an arbitrary unitary?

Universality means an arbitrary unitary can be rebuilt from the elementary set; the proof is the
decomposition itself. Take a random two-qubit unitary, decompose it, and check that recomposing
the pieces reproduces it exactly.

```wolfram
SeedRandom[7];
u = QuantumOperator["RandomUnitary"[{2, 2}]];
dec = u["Decompose"];
Norm[(u["MatrixRepresentation"] // Normal) - (dec["MatrixRepresentation"] // Normal), "Frobenius"] // Chop
```

The recomposed circuit matches the original to machine zero: an arbitrary two-qubit unitary is
exactly the elementary gates QF emitted. This is universality made executable, and it scales to
any number of qubits.

**Summary of Part F.** Circuits are unitaries assembled from one- and two-qubit gates. We built
the standard gates and saw their basis action, composed and measured the depth of circuits,
compiled a circuit to its unitary, and decomposed two-qubit gates by KAK, witnessing universality
by exactly rebuilding a random unitary from elementary pieces.

---

## Part G. Qudits and higher-dimensional discrete systems

Nothing in quantum theory forces the dimension to be two. A qudit lives in $\mathbb{C}^d$, and
the natural operators are the Weyl-Heisenberg clock and shift, the discrete cousins of position
and momentum. Higher dimensions are not a curiosity: they carry more information per carrier,
realize spin-$j$ systems directly, and exhibit contextuality that two-level systems cannot.

#### G1. How do I build the clock and shift (Weyl-Heisenberg) operators?

The shift $X$ cycles the basis states $|k\rangle \to |k+1 \bmod d\rangle$ and the clock $Z$
phases them by powers of $\omega = e^{2\pi i/d}$. Together they generate the qudit Pauli group and
obey the Weyl relation $ZX = \omega\, XZ$.

```wolfram
\[Omega] = Exp[2 Pi I/3];
Xd = RotateRight[IdentityMatrix[3]];        (* shift *)
Zd = DiagonalMatrix[\[Omega]^Range[0, 2]];  (* clock *)
Xd
Zd
Zd . Xd == \[Omega] Xd . Zd                 (* Weyl relation *)
```

The shift is the cyclic permutation matrix; the clock is diagonal in cube roots of unity; and they
fail to commute by exactly the phase $\omega$. This noncommutation is the engine of qudit phase
space and the mutually unbiased bases below.

#### G2. How do I build a qutrit entangling gate (the SUM gate)?

The SUM gate is the qudit generalization of CNOT: it adds the control register into the target
modulo $d$, $|a\rangle|b\rangle \to |a\rangle|a+b \bmod d\rangle$. Build it as
$\sum_a |a\rangle\langle a| \otimes X^a$, the shift applied $a$ times conditioned on the control.

```wolfram
ketd[d_, ks__] := Fold[Flatten@KroneckerProduct[#1, #2] &, UnitVector[d, # + 1] & /@ {ks}];
sum = Sum[KroneckerProduct[KroneckerProduct[UnitVector[3, a + 1], UnitVector[3, a + 1]],
     MatrixPower[Xd, a]], {a, 0, 2}];
sum . ketd[3, 1, 0]    (* |1>|0> -> |1>|1> *)
```

The input $|1\rangle|0\rangle$ maps to $|1\rangle|1\rangle$: the target becomes $1 + 0 = 1$. SUM
plays the role of CNOT in qudit error correction and teleportation.

#### G3. How do I work with spin-$j$ operators?

A spin-$j$ particle is a $(2j+1)$-dimensional discrete system with angular-momentum operators
$J_x, J_y, J_z$ obeying $[J_x, J_y] = i J_z$. Build the spin-one matrices and check the algebra.

```wolfram
Jz = {{1, 0, 0}, {0, 0, 0}, {0, 0, -1}};
Jx = {{0, 1, 0}, {1, 0, 1}, {0, 1, 0}}/Sqrt[2];
Jy = {{0, -I, 0}, {I, 0, -I}, {0, I, 0}}/Sqrt[2];
Eigenvalues[Jz]
Commutator[Jx, Jy, Dot] == I Jz
Jx . Jx + Jy . Jy + Jz . Jz
```

The spin-one $J_z$ spectrum is $\{-1, 0, 1\}$, the algebra $[J_x, J_y] = i J_z$ holds, and the
Casimir $J^2 = j(j+1) I = 2 I$. A discrete three-level system carries the full $SU(2)$ structure.

#### G4. How do I build mutually unbiased bases (the qudit Fourier basis)?

Two bases are mutually unbiased if a state sharp in one is uniform in the other; the computational
and Fourier bases are the canonical pair. The qudit Fourier transform $F_{jk} = \omega^{jk}/\sqrt d$
maps between them, and the overlap of any two cross-basis states has constant modulus $1/\sqrt d$.

```wolfram
w = Exp[2 Pi I/3];
F = (1/Sqrt[3]) Table[w^(j k), {j, 0, 2}, {k, 0, 2}];
Table[Simplify[Abs[F[[1, k]]]^2], {k, 1, 3}]
Simplify[ConjugateTranspose[F] . F] == IdentityMatrix[3]
```

Every computational state has probability exactly $1/3$ in each Fourier outcome: the bases are
mutually unbiased, and $F$ is unitary. Mutual unbiasedness is the resource behind quantum key
distribution and optimal state tomography.

#### G5. How do I exhibit quantum contextuality (the KCBS qutrit inequality)?

Contextuality is the impossibility of assigning outcomes to all measurements independently of which
compatible set they are measured with. The KCBS inequality bounds a sum of five projector
expectations by $2$ for any noncontextual model, but a qutrit reaches $\sqrt5 > 2$. Build the five
pentagon vectors and the optimal state.

```wolfram
ct2 = Cos[Pi/5]/(1 + Cos[Pi/5]); ct = Sqrt[ct2]; st = Sqrt[1 - ct2];
lvec[k_] := {ct, st Sin[4 Pi k/5], st Cos[4 Pi k/5]};
psi = {1, 0, 0};
FullSimplify[lvec[0] . lvec[1]]
FullSimplify[Sum[Abs[Conjugate[psi] . lvec[k]]^2, {k, 0, 4}]]
```

Consecutive pentagon vectors are orthogonal (jointly measurable), yet the sum of the five
expectations is $\sqrt5 \approx 2.236$, above the noncontextual bound $2$. A single qutrit already
defeats any noncontextual hidden-variable account.

**Summary of Part G.** Discrete space extends past the qubit. We built the clock and shift and
their Weyl relation, the qutrit SUM gate, the spin-one angular-momentum algebra with its Casimir,
the mutually unbiased Fourier basis, and a qutrit violation of the KCBS contextuality bound. Higher
dimensions bring genuinely new phenomena, not just bigger matrices.

---

## Part H. Quantum channels and open systems

No quantum system is truly isolated. Coupling to an environment turns unitary evolution into a
quantum channel, a completely positive trace-preserving map that can decohere, dissipate, and
destroy entanglement. The density matrix of Part E is exactly the object channels act on, and the
Lindblad equation is their continuous-time generator.

#### H1. How do I build a channel from Kraus operators and apply it?

A channel acts as $\rho \mapsto \sum_k K_k \rho K_k^\dagger$. Amplitude damping, the decay of an
excited state, has Kraus operators $K_0 = \mathrm{diag}(1, \sqrt{1-\gamma})$ and
$K_1 = \sqrt\gamma\,|0\rangle\langle1|$. Apply it to the excited state.

```wolfram
g = 3/10;
K0 = {{1, 0}, {0, Sqrt[1 - g]}}; K1 = {{0, Sqrt[g]}, {0, 0}};
channel[\[Rho]_, kr_] := Sum[k . \[Rho] . ConjugateTranspose[k], {k, kr}];
channel[{{0, 0}, {0, 1}}, {K0, K1}]    (* apply to |1><1| *)
```

The excited population $\rho_{11}$ drops from $1$ to $1 - \gamma = 0.7$ while $0.3$ leaks into the
ground state: the channel models spontaneous emission with decay probability $\gamma$.

#### H2. How do I get the Choi matrix and check complete positivity?

The Choi-Jamiolkowski isomorphism turns a channel into a state, $J = \sum_k |K_k\rangle\rangle
\langle\langle K_k|$ where $|K\rangle\rangle$ stacks the columns of $K$. The channel is completely
positive if and only if its Choi matrix is positive semidefinite.

```wolfram
vec[K_] := Flatten[K];
choi = Sum[KroneckerProduct[vec[K], Conjugate[vec[K]]], {K, {K0, K1}}];
choi
PositiveSemidefiniteMatrixQ[N[choi]]
```

The Choi matrix is positive semidefinite, so amplitude damping is a legitimate completely positive
map. Complete positivity is the condition that the map stays physical even when applied to half of
an entangled pair.

#### H3. How do I apply the standard named noise channels?

The common noise models are just lists of Kraus operators: bit flip $\{\sqrt{1-p}\,I, \sqrt p\,X\}$
and depolarizing $\{\sqrt{1-p}\,I, \sqrt{p/3}\,X, \sqrt{p/3}\,Y, \sqrt{p/3}\,Z\}$.

```wolfram
bitflip[p_] := {Sqrt[1 - p] IdentityMatrix[2], Sqrt[p] PauliMatrix[1]};
depol[p_] := Prepend[Sqrt[p/3] # & /@ Table[PauliMatrix[j], {j, 3}], Sqrt[1 - p] IdentityMatrix[2]];
channel[{{1, 0}, {0, 0}}, bitflip[1/4]]          (* bit flip on |0> *)
channel[densityMatrix[{1, 0, 0}], depol[1/2]]    (* depolarize |+> *)
```

The bit-flip channel moves weight $\tfrac14$ from $|0\rangle$ to $|1\rangle$; the depolarizing
channel shrinks the Bloch vector of $|+\rangle$ toward the maximally mixed state. These are the
elementary error models for quantum hardware.

#### H4. How do I verify trace preservation and get the Stinespring dilation?

A channel is trace preserving when $\sum_k K_k^\dagger K_k = I$, and by Stinespring's theorem every
channel is a unitary on a larger space followed by tracing out the environment. Stacking the Kraus
operators vertically gives that dilation isometry.

```wolfram
Simplify[ConjugateTranspose[K0] . K0 + ConjugateTranspose[K1] . K1]   (* = I, trace preserving *)
stinespring = ArrayFlatten[{{K0}, {K1}}];          (* stack Kraus -> 4x2 isometry V *)
Simplify[ConjugateTranspose[stinespring] . stinespring == IdentityMatrix[2]]
```

The completeness relation gives the identity (trace preserving), and the stacked Kraus operators
form a $4\times2$ isometry $V$ with $V^\dagger V = I$: one input qubit dilated to a
system-plus-environment, whose partial trace over the environment is the channel. Decoherence is
entanglement with an environment we then discard.

#### H5. How do I solve the Lindblad master equation (T1 relaxation)?

Continuous open dynamics obey the Lindblad equation
$\dot\rho = -i[H,\rho] + \sum_k (L_k \rho L_k^\dagger - \tfrac12\{L_k^\dagger L_k, \rho\})$. With no
Hamiltonian and a single lowering jump $L = |0\rangle\langle1|$, `DSolveValue` integrates the matrix
ODE in closed form.

```wolfram
ClearAll[\[Rho]];            (* \[Rho] held a value above; free it to use as the ODE function *)
L = {{0, 1}, {0, 0}};        (* lowering |0><1| *)
diss[m_] := L . m . ConjugateTranspose[L] -
   1/2 (ConjugateTranspose[L] . L . m + m . ConjugateTranspose[L] . L);
\[Rho]t = DSolveValue[{\[Rho]'[t] == diss[\[Rho][t]], \[Rho][0] == {{0, 0}, {0, 1}}},
   \[Rho][t] \[Element] Matrices[{2, 2}], t]
```

The solution is $\rho(t) = \mathrm{diag}(1 - e^{-t},\, e^{-t})$: the excited population decays as
$e^{-t}$, reaching $e^{-1} = 0.368$ after one lifetime. This is the hallmark exponential $T_1$
decay, with the jump operator $|0\rangle\langle1|$ as spontaneous emission in continuous time.

#### H6. How do I find the steady state of an open system?

A dissipative channel drives the system to a steady state independent of the start. For pure
relaxation, that fixed point is the ground state. Take the long-time limit of the closed-form
solution.

```wolfram
Limit[\[Rho]t, t -> Infinity]
```

The limit is $\mathrm{diag}(1, 0) = |0\rangle\langle0|$: the system relaxes entirely to its ground
state, the unique steady state of amplitude damping. Dissipation, not unitary dynamics, sets where
an open system ends up.

#### H7. How does noise destroy entanglement, and can it die at finite time?

Local noise degrades entanglement, and remarkably the entanglement can reach exactly zero at
finite noise, before the noise itself is complete: entanglement sudden death. Track the negativity
of a balanced Bell pair versus an asymmetric state under amplitude damping on both qubits.

```wolfram
neg[vec_, gg_] := Module[{k0, k1, rho, kr},
   k0 = {{1, 0}, {0, Sqrt[1 - gg]}}; k1 = {{0, Sqrt[gg]}, {0, 0}};
   rho = KroneckerProduct[vec, Conjugate@vec];
   kr = Flatten[Table[KroneckerProduct[a, b], {a, {k0, k1}}, {b, {k0, k1}}], 1];
   -Total@Select[Re@Eigenvalues@partialT[Sum[m . rho . ConjugateTranspose[m], {m, kr}]], # < 0 &]];
Table[neg[{1, 0, 0, 1}/Sqrt[2], gg], {gg, {0, 0.3, 0.5, 0.9}}]
Table[neg[{Sqrt[1/5], 0, 0, Sqrt[4/5]}, gg], {gg, {0, 0.3, 0.5, 0.6}}]
```

The balanced state $|\Phi^+\rangle$ keeps a positive negativity until $\gamma \to 1$. The asymmetric
state $\sqrt{1/5}\,|00\rangle + \sqrt{4/5}\,|11\rangle$ instead loses all entanglement at finite
$\gamma = \tfrac12$ and stays dead: entanglement sudden death, a finite-time end with no analogue
in the decay of any single amplitude.

**Summary of Part H.** Open systems evolve by completely positive channels. We built channels from
Kraus operators, certified complete positivity through the Choi matrix, applied named noise,
exhibited trace preservation and the Stinespring dilation, integrated the Lindblad equation to see
$T_1$ decay and relaxation to a ground steady state, and watched entanglement die at finite noise.

---

## Part I. The stabilizer formalism and quantum error correction

A large class of quantum states and circuits, the stabilizer (Clifford) class, can be simulated
efficiently on a classical computer by tracking not the exponential amplitude vector but a compact
tableau of Pauli operators. This is the formalism in which quantum error correction is written, and
it draws the sharp line between Clifford circuits (classically tractable) and the magic that makes
quantum computing hard to simulate.

#### I1. How do I represent a stabilizer state and read its generators?

A stabilizer state is the simultaneous $+1$ eigenstate of $n$ commuting Pauli operators, its
stabilizer generators. `PauliStabilizer` stores the tableau; `["Stabilizers"]` reports the
generators as Pauli strings.

```wolfram
ps = PauliStabilizer[2];
ps["Stabilizers"]
```

The two-qubit register $|00\rangle$ is stabilized by $Z_1$ and $Z_2$ (the strings $ZI$ and $IZ$):
these two operators fix the state uniquely. An $n$-qubit stabilizer state needs only $n$ such
generators, an exponential compression over the $2^n$ amplitudes.

#### I2. How do I apply Clifford gates and watch the generators transform?

Clifford gates (Hadamard, phase, CNOT) map Pauli strings to Pauli strings, so they act directly on
the generators. Preparing a Bell state turns the computational stabilizers into the Bell
stabilizers. The bracket form applies a whole circuit.

```wolfram
psBell = PauliStabilizer[2][{"H", "CNOT"}];
psBell["Stabilizers"]
psBell["QuantumState"]["StateVector"]
```

The Bell state $|\Phi^+\rangle$ is stabilized by $XX$ and $ZZ$, and materializing the tableau
recovers the familiar amplitudes. The entire evolution happened on two Pauli strings, never on the
state vector.

#### I3. How do I simulate a large Clifford circuit efficiently?

Because the tableau grows linearly, Clifford circuits on hundreds of qubits run in milliseconds,
far beyond the reach of dense simulation. Build a sixty-qubit entangling circuit and confirm it
stays a valid stabilizer state.

```wolfram
n = 60;
circ = Join[Table["H" -> q, {q, n}], Table["CNOT" -> {q, q + 1}, {q, n - 1}]];
psBig = PauliStabilizer[n][circ];
psBig["Qubits"]
StabilizerStateQ[psBig]
```

A sixty-qubit state, whose amplitude vector would have $2^{60} \approx 10^{18}$ entries, is handled
exactly through its tableau. This Gottesman-Knill efficiency is why Clifford circuits alone cannot
provide a quantum advantage.

#### I4. How do I build a graph (cluster) state?

A graph state attaches a qubit in $|+\rangle$ to each vertex and entangles neighbors with a
controlled-$Z$. Its stabilizers are $K_i = X_i \prod_{j \sim i} Z_j$, one per vertex. Build the
cluster state of a three-vertex path.

```wolfram
gcirc = {"H" -> 1, "H" -> 2, "H" -> 3, "CZ" -> {1, 2}, "CZ" -> {2, 3}};
PauliStabilizer[3][gcirc]["Stabilizers"]
```

The generators are $X_1 Z_2$, $Z_1 X_2 Z_3$, $Z_2 X_3$: each vertex carries an $X$ dressed by $Z$
on its neighbors, exactly the graph-state stabilizer rule. Cluster states are the substrate of
measurement-based quantum computation.

#### I5. How do I encode the three-qubit code and detect an error by syndrome?

The three-qubit bit-flip code stores one logical qubit as $|0_L\rangle = |000\rangle$,
$|1_L\rangle = |111\rangle$. The stabilizers $Z_1 Z_2$ and $Z_2 Z_3$ check parities without
disturbing the logical state; a flipped qubit shows up as a syndrome of $-1$ values.

```wolfram
enc = QuantumCircuitOperator[{"CNOT" -> {1, 2}, "CNOT" -> {1, 3}}];
errored = QuantumCircuitOperator[{"X" -> 2}][enc[QuantumState["100"]]];
QuantumMeasurementOperator[QuantumOperator["ZZI"]][errored]["Mean"]
QuantumMeasurementOperator[QuantumOperator["IZZ"]][errored]["Mean"]
```

The encoded $|1_L\rangle = |111\rangle$ suffers an $X$ on qubit $2$, and both parity checks read
$-1$. The syndrome $(-1, -1)$ points uniquely to qubit $2$, which is then corrected, all without
measuring the logical information.

#### I6. How do I work with the five-qubit code stabilizers?

The smallest code correcting an arbitrary single-qubit error is the five-qubit code, with four
stabilizer generators that are cyclic shifts of $XZZXI$. For the code space to exist the generators
must commute.

```wolfram
s1 = QuantumOperator["XZZXI"]; s2 = QuantumOperator["IXZZX"];
s3 = QuantumOperator["XIXZZ"]; s4 = QuantumOperator["ZXIXZ"];
commute[a_, b_] := (a @ b)["MatrixRepresentation"] == (b @ a)["MatrixRepresentation"];
And @@ (commute @@@ Subsets[{s1, s2, s3, s4}, {2}])
```

All four generators pairwise commute, so they define a joint $+1$ eigenspace: a two-dimensional
code holding one logical qubit inside five physical ones. This $[[5,1,3]]$ code is the perfect code,
saturating the quantum Hamming bound.

#### I7. How do I detect magic (the non-stabilizer resource)?

Clifford gates keep a state in the stabilizer class; a single non-Clifford $T$ gate takes it out.
A stabilizer state has an exact tableau, while $T|+\rangle$ does not, and it appears instead as a
two-term stabilizer decomposition, the seed of exponential simulation cost.

```wolfram
MatchQ[PauliStabilizer[QuantumState["+"]], _PauliStabilizer]
MatchQ[PauliStabilizer[QuantumOperator["T"][QuantumState["+"]]], _PauliStabilizer]
Length[(PauliStabilizer[1]["T", 1])["Components"]]
```

The state $|+\rangle$ is a stabilizer state, but $T|+\rangle$ is not, and the $T$ gate splits the
tableau into two stabilizer components. This magic, the departure from the Clifford class, is the
resource that makes quantum computation classically hard and is consumed by magic-state
distillation.

**Summary of Part I.** The stabilizer formalism tracks Pauli generators instead of amplitudes,
making Clifford circuits classically efficient. We read and transformed stabilizer generators,
simulated a sixty-qubit circuit, built a cluster state, ran the three-qubit and five-qubit codes,
and detected the magic a $T$ gate injects. This is the language of fault tolerance.

---

## Part J. Discrete phase space and the Wigner function

A finite-dimensional state has a discrete Wigner function, a real quasi-probability distribution on
a $d \times d$ phase-space grid. For odd dimension it is especially clean. Its decisive feature is
that the stabilizer states are exactly the ones with a non-negative Wigner function; negativity is
a witness of non-classicality and a necessary resource for quantum advantage.

#### J1. How do I compute the discrete Wigner function of a qudit?

For odd $d$ the Wigner function is the $d \times d$ phase-space representation, a genuine
probability distribution (non-negative, summing to one) for stabilizer states. Compute it for a
qutrit computational state.

```wolfram
QuantumState["0", 3]["PhaseSpace"]
Total[Flatten[QuantumState["0", 3]["PhaseSpace"] // Normal]]
Min[QuantumState["0", 3]["PhaseSpace"] // Normal]
```

The Wigner function of $|0\rangle$ is non-negative and normalized: a qutrit stabilizer state is
classical in phase space. Its marginals reproduce the position and momentum probability
distributions.

#### J2. How do I detect Wigner negativity?

A non-stabilizer (magic) state has negative Wigner entries. The qutrit "strange" state
$(|1\rangle - |2\rangle)/\sqrt2$ is the textbook example.

```wolfram
strange = QuantumState[{0, 1, -1}/Sqrt[2], 3];
strange["PhaseSpace"]
Min[strange["PhaseSpace"] // Normal]
```

The Wigner function dips to $-\tfrac13$: the strange state is genuinely non-classical. This
negativity, impossible for any stabilizer state, is what a quantum computer must supply beyond the
efficiently simulable Clifford part.

#### J3. How do I compute the mana (a magic monotone)?

The mana quantifies negativity as $\mathcal{M} = \log \sum_{u} |W_u|$, the logarithm of the
Wigner $1$-norm. It is zero exactly for non-negative (stabilizer) states and positive for magic.

```wolfram
mana[qs_] := Log[Total[Abs[Flatten[qs["PhaseSpace"] // Normal]]]];
FullSimplify[mana[strange]]
FullSimplify[mana[QuantumState["0", 3]]]
```

The strange state has mana $\log\tfrac53 \approx 0.511$, while the stabilizer state has mana zero.
Mana is additive and non-increasing under stabilizer operations, making it a genuine resource
measure of magic.

**Summary of Part J.** The discrete Wigner function maps finite-dimensional states onto a
quasi-probability grid. We computed it for a stabilizer qutrit (non-negative), exhibited the
negativity of a magic state, and quantified that negativity by the mana. Wigner negativity is the
phase-space face of the same magic resource Part I found in the tableau.

---

## Part K. Quantum algorithms

Here the pieces assemble into algorithms: superpositions queried in parallel, interference that
amplifies the right answer, entanglement spent as a communication resource. Each item below states
the task and gives the QuantumFramework realization that solves it, from the textbook oracle
problems through search, transforms, optimization, and the entanglement protocols.

#### K1. Deutsch-Jozsa: how do I tell a constant from a balanced function in one query?

A function $f$ on $n$ bits is promised constant or balanced; classically distinguishing them takes
up to $2^{n-1}+1$ queries, but one quantum query suffices. The circuit is Hadamards, a phase
oracle, Hadamards; the all-zeros outcome appears with certainty if and only if $f$ is constant.
Here the balanced oracle for $f(x) = x_1$ is a single $Z$.

```wolfram
ket[bits__] := Fold[Flatten@KroneckerProduct[#1, #2] &, UnitVector[2, # + 1] & /@ {bits}];
H = {{1, 1}, {1, -1}}/Sqrt[2]; H3 = KroneckerProduct[H, H, H];
dj[oracle_] := H3 . oracle . H3 . ket[0, 0, 0];
Abs[dj[IdentityMatrix[8]][[1]]]^2     (* constant f -> 1 *)
Abs[dj[KroneckerProduct[PauliMatrix[3], IdentityMatrix[2], IdentityMatrix[2]]][[1]]]^2   (* balanced -> 0 *)
```

The constant function returns the all-zeros string with probability $1$; the balanced function
returns it with probability $0$. One query reads a global property of $f$ that no single classical
evaluation can.

#### K2. Bernstein-Vazirani: how do I recover a hidden string in one query?

Given $f(x) = s \cdot x$ for a hidden $s$, the same Hadamard-oracle-Hadamard sandwich returns $s$
itself in one query (classically it takes $n$). The oracle places a $Z$ on each bit where
$s_i = 1$; here $s = 101$.

```wolfram
bvOracle = KroneckerProduct[PauliMatrix[3], IdentityMatrix[2], PauliMatrix[3]];  (* Z on bits 1, 3 *)
H3 . bvOracle . H3 . ket[0, 0, 0] // Chop
```

The output is the basis state $|101\rangle$ (index $6$): the hidden string $s = 101$ is read off
directly. Interference concentrates the entire amplitude on the answer.

#### K3. Simon: how do I find a hidden period?

Simon's problem (find $s$ with $f(x) = f(x \oplus s)$) is the first with an exponential quantum
speedup, and the ancestor of Shor's. Each run measures a string $y$ orthogonal to $s$; a handful
of runs determine $s$ by linear algebra over $\mathbb{F}_2$. For $f(x) = x_1 \oplus x_2$ the hidden
period is $s = 11$, so the input register must yield only $y \in \{00, 11\}$.

```wolfram
simon = QuantumCircuitOperator[{"H" -> 1, "H" -> 2, "CNOT" -> {1, 3}, "CNOT" -> {2, 3}, "H" -> 1, "H" -> 2}];
QuantumPartialTrace[simon[QuantumState["000"]], {3}]["ProbabilitiesList"] // Chop
```

The input register is supported only on $00$ and $11$, exactly the strings $y$ with $y \cdot s = 0$
for $s = 11$. Collecting such $y$ and solving the linear system recovers the period in a polynomial
number of runs.

#### K4. Grover: how do I amplify a marked item?

Grover search finds a marked item among $N$ in $O(\sqrt N)$ steps, a quadratic speedup over
classical $O(N)$. Each iteration is the oracle $O$ (phase-flip the target) followed by the diffuser
$D$ (reflect about the uniform state); iterate $\approx \tfrac\pi4\sqrt N$ times. For one marked
item in eight, two iterations nearly saturate.

```wolfram
NN = 8; tgt = 6;
oracle = IdentityMatrix[NN]; oracle[[tgt, tgt]] = -1;   (* phase-flip the marked |101> *)
s = ConstantArray[1, NN]/Sqrt[NN];
diffuser = 2 Outer[Times, s, s] - IdentityMatrix[NN];   (* reflect about uniform *)
grov = Nest[diffuser . oracle . # &, s, 2];
Abs[grov[[tgt]]]^2
```

After two iterations the marked state $|101\rangle$ carries probability $\tfrac{121}{128} \approx
0.945$, matching the exact amplitude $\sin^2(5\arcsin(1/\sqrt8))$. Amplitude amplification rotates
the state toward the target far faster than any classical scan.

#### K5. Quantum Fourier transform: how do I build and verify it?

The QFT is the change of basis at the heart of phase estimation and Shor's algorithm,
$|k\rangle \mapsto \tfrac1{\sqrt N}\sum_j e^{2\pi i jk/N}|j\rangle$. On $|0\dots0\rangle$ it produces
the uniform superposition.

```wolfram
Abs[FourierMatrix[8] . UnitVector[8, 1]]^2
```

The QFT (`FourierMatrix[8]` on three qubits) spreads $|000\rangle$ uniformly over all eight basis
states, as its $k=0$ column demands. Unlike the classical FFT, the QFT runs in $O(n^2)$ gates on
$n$ qubits.

#### K6. Quantum phase estimation: how do I read out an eigenphase?

Phase estimation extracts the eigenphase $\varphi$ of a unitary, the subroutine behind Shor and
quantum chemistry. Its core is the inverse QFT, which maps a phase-kicked register
$\tfrac1{\sqrt N}\sum_k e^{2\pi i\varphi k}|k\rangle$ to the binary digits of $\varphi$. Encode
$\varphi = 3/8$ and read it back.

```wolfram
phi = 3/8;
kicked = Table[Exp[2 Pi I phi k], {k, 0, 7}]/Sqrt[8];
Abs[ConjugateTranspose[FourierMatrix[8]] . kicked]^2 // Chop   (* inverse QFT *)
```

The inverse QFT concentrates all probability on index $3$, the binary fraction $0.011 = 3/8 =
\varphi$. The controlled-unitary ladder of full phase estimation exists only to prepare this kicked
register, which the inverse QFT then decodes.

#### K7. HHL: how do I apply a matrix to a vector coherently?

Quantum linear-algebra algorithms (such as HHL) hinge on applying a matrix $A$ to an amplitude-
encoded vector. They write $A = \sum_\mu c_\mu \sigma_\mu$ (its Pauli decomposition, with
$c_\mu = \tfrac12\mathrm{Tr}(\sigma_\mu A)$) and apply each Pauli term coherently, weighted by
$c_\mu$ through an ancilla, a block-encoding.

```wolfram
A = {{1, 1/2}, {1/2, 1}}; xvec = {3/5, 4/5};
c = Table[Tr[s . A]/2, {s, Table[PauliMatrix[j], {j, 0, 3}]}];   (* Pauli coefficients *)
c
Sum[c[[j + 1]] PauliMatrix[j] . xvec, {j, 0, 3}]                 (* recombine -> A.xvec *)
```

The decomposition is $A = I + \tfrac12 X$, and recombining the Pauli terms reproduces
$A x = (1, \tfrac{11}{10})$. On hardware an ancilla-controlled multiplexer applies this weighted sum
coherently; QuantumFramework's `QuantumLinearSolve` packages the full variational solver around this
block-encoding.

#### K8. VQE: how do I find a ground-state energy variationally?

The variational quantum eigensolver minimizes $\langle\psi(\theta)|H|\psi(\theta)\rangle$ over a
parametrized state to approximate the ground energy, the workhorse of near-term quantum chemistry.
For $H = Z$ with an $R_y(\theta)$ ansatz the minimum is the lowest eigenvalue $-1$.

```wolfram
ClearAll[ry];                (* ry held the B5 matrix above; free it to define the ansatz *)
ry[\[Theta]_] := MatrixExp[-I \[Theta]/2 PauliMatrix[2]];
cost[\[Theta]_] := Re[Conjugate[#] . PauliMatrix[3] . #] &[ry[\[Theta]] . {1, 0}];
NMinimize[{cost[\[Theta]], -\[Pi] <= \[Theta] <= \[Pi]}, \[Theta]][[1]]
```

The optimizer drives the energy to $-1$, the ground state $|1\rangle$ of $Z$. Use `NMinimize`
rather than a local method: variational landscapes have barren plateaus that defeat gradient
descent.

#### K9. QAOA: how do I solve Max-Cut?

The quantum approximate optimization algorithm encodes a combinatorial problem in a cost
Hamiltonian and alternates cost and mixer layers. For Max-Cut the cost is
$H_C = \sum_{(i,j)\in E} \tfrac12(I - Z_i Z_j)$, whose largest eigenvalue is the size of the maximum
cut. On a triangle that maximum is two.

```wolfram
Zat[positions_, n_] := KroneckerProduct @@ Table[If[MemberQ[positions, k], PauliMatrix[3], IdentityMatrix[2]], {k, n}];
Hc = Total[1/2 (IdentityMatrix[8] - Zat[#, 3]) & /@ {{1, 2}, {2, 3}, {1, 3}}];
Max[Eigenvalues[Hc]]
```

The cost Hamiltonian's top eigenvalue is $2$: any cut of the triangle separates at most two of its
three edges. QAOA prepares a state concentrated on the optimizing bit strings by tuning the
layer angles toward this maximum.

#### K10. Teleportation: how do I teleport a qubit?

Teleportation moves an unknown state across a shared Bell pair using two classical bits. Replacing
the measurement-and-correction with coherent controlled gates (deferred measurement) gives a
unitary circuit whose third qubit ends in the input state.

```wolfram
psiIn = QuantumState[{Cos[Pi/5], Sin[Pi/5] Exp[I Pi/3]}];
init = QuantumTensorProduct[psiIn, QuantumState["PhiPlus"]];
tele = QuantumCircuitOperator[{"CNOT" -> {1, 2}, "H" -> 1, "CNOT" -> {2, 3}, "CZ" -> {1, 3}}];
q3 = QuantumPartialTrace[tele[init], {1, 2}];
N@QuantumSimilarity[q3, psiIn]
N@q3["Purity"]
```

The output qubit matches the input with fidelity $1$ and is pure: the unknown amplitudes have moved
from qubit $1$ to qubit $3$ with no copy left behind, consistent with no-cloning. Entanglement plus
classical communication transmits a quantum state.

#### K11. Superdense coding: how do I send two bits with one qubit?

Superdense coding is teleportation's dual: sharing a Bell pair, one party encodes two classical
bits by applying one of $\{I, X, Z, XZ\}$ to a single qubit, and a Bell measurement decodes both.

```wolfram
encode[g_] := QuantumCircuitOperator[Join[{"H" -> 1, "CNOT" -> {1, 2}}, g, {"CNOT" -> {1, 2}, "H" -> 1}]];
decode[g_] := encode[g][QuantumState["00"]]["ProbabilitiesList"];
{decode[{}], decode[{"X" -> 1}], decode[{"Z" -> 1}], decode[{"Z" -> 1, "X" -> 1}]}
```

The four Pauli encodings decode deterministically to the four messages $00, 01, 10, 11$. One
transmitted qubit carries two bits, because the shared entanglement was prepaid.

#### K12. Entanglement swapping: how do I entangle two qubits that never interacted?

Entanglement swapping creates entanglement between particles with no common past. Start with two
independent Bell pairs $(1,2)$ and $(3,4)$; a Bell measurement on $(2,3)$ leaves $1$ and $4$
entangled.

```wolfram
twoPairs = QuantumTensorProduct[QuantumState["PhiPlus"], QuantumState["PhiPlus"]];
N@QuantumEntanglementMonotone[QuantumPartialTrace[twoPairs, {2, 3}], {1}, "Concurrence"]
bm = QuantumCircuitOperator[{"CNOT" -> {2, 3}, "H" -> 2}][twoPairs];
branch = QuantumMeasurementOperator[{2, 3}][bm]["States"][[1]]["Normalized"];
N@QuantumEntanglementMonotone[QuantumPartialTrace[branch, {2, 3}], {1}, "Concurrence"]
```

Qubits $1$ and $4$ start with concurrence $0$ (independent pairs) and, after the Bell measurement
on $(2,3)$, reach concurrence $1$, maximally entangled. The measurement projected the swap, an
entanglement created without any interaction between $1$ and $4$.

**Summary of Part K.** Algorithms turn superposition, interference, and entanglement into
computational and communication advantage. We ran the Deutsch-Jozsa, Bernstein-Vazirani, and Simon
oracle problems, Grover search, the QFT and phase-estimation readout, the coherent matrix
primitive of HHL, VQE and QAOA optimization, and the teleportation, superdense-coding, and
entanglement-swapping protocols. Each is a few lines once the objects of Parts A through J are in
hand.

---

## Part L. Foundations and the limits of quantum information

Quantum theory is hemmed in by sharp impossibilities and sharp excesses. You cannot copy an
unknown state, you cannot extract more than $\log d$ bits from a $d$-level system, and yet you can
violate every local-realistic bound by a definite margin. These are not engineering limits but
theorems, and each is a short computation.

#### L1. No-cloning: why can no machine copy an unknown state?

A would-be universal cloner acting unitarily would have to satisfy
$\langle\psi|\phi\rangle = \langle\psi|\phi\rangle^2$ for all pairs, since inner products are
preserved; this forces every overlap to be $0$ or $1$. Two genuinely different non-orthogonal
states violate it.

```wolfram
ov = Conjugate[{1, 0}] . ({1, 1}/Sqrt[2]);   (* <0|+> *)
{ov, ov^2, ov == ov^2}
```

The overlap $\tfrac1{\sqrt2}$ does not equal its square $\tfrac12$, so no unitary can clone both
$|0\rangle$ and $|+\rangle$. No-cloning is what makes quantum money and quantum key distribution
secure, and what keeps teleportation honest.

#### L2. CHSH: how large is the quantum violation of local realism?

The CHSH quantity $S = \langle A_0 B_0\rangle + \langle A_0 B_1\rangle + \langle A_1 B_0\rangle -
\langle A_1 B_1\rangle$ is bounded by $2$ for any local hidden-variable theory, but quantum
mechanics on the singlet reaches Tsirelson's bound $2\sqrt2$ with the right measurement angles.

```wolfram
rhoS = QuantumState["PsiMinus"]["DensityMatrix"] // Normal;
pX = {{0, 1}, {1, 0}}; pZ = {{1, 0}, {0, -1}};
Edir[A_, B_] := Re@Tr[KroneckerProduct[A, B] . rhoS];
B0 = (pZ + pX)/Sqrt[2]; B1 = (pZ - pX)/Sqrt[2];
Edir[pZ, B0] + Edir[pZ, B1] + Edir[pX, B0] - Edir[pX, B1] // FullSimplify
```

The magnitude $2\sqrt2 \approx 2.83$ exceeds the classical bound $2$: no local hidden-variable model
can reproduce the singlet's correlations. The violation is the experimental signature that nature
is not locally real.

#### L3. Holevo: how much classical information fits in a qubit?

The Holevo bound caps the accessible information of an ensemble by
$\chi = S(\bar\rho) - \sum_i p_i S(\rho_i)$, which for an ensemble of pure states is just
$S(\bar\rho) \le \log_2 d$. Even cleverly chosen non-orthogonal signals cannot beat one bit per
qubit.

```wolfram
rhoEns = 1/2 densityMatrix[{0, 0, 1}] + 1/2 densityMatrix[{1, 0, 0}];   (* 1/2|0><0| + 1/2|+><+| *)
N@vne[rhoEns]
```

The ensemble $\{|0\rangle, |+\rangle\}$ carries Holevo information $\approx 0.60$ bits, below the
one-bit ceiling because the two states are non-orthogonal and hence not perfectly distinguishable.
A qubit holds one qubit of quantum state but at most one bit of retrievable classical message.

#### L4. The GHZ paradox: how do three qubits refute local realism with certainty?

Where CHSH is statistical, the GHZ state gives a single-shot contradiction. The observables
$XYY$, $YXY$, $YYX$ each have value $-1$ on $|\mathrm{GHZ}\rangle$, so any local assignment forces
their product $XXX = (-1)^3 = -1$; quantum mechanics says $XXX = +1$.

```wolfram
ghz = {1, 0, 0, 0, 0, 0, 0, 1}/Sqrt[2];
X = PauliMatrix[1]; Y = PauliMatrix[2];
ev[ops_] := Re[Conjugate[ghz] . KroneckerProduct @@ ops . ghz];
{ev[{X, X, X}], ev[{X, Y, Y}], ev[{Y, X, Y}], ev[{Y, Y, X}]}
```

The three mixed correlators are $-1$ and would force $XXX = -1$ in any local-realistic account, yet
$XXX = +1$. The contradiction is deterministic, needing no inequality and no statistics: local
realism fails outright.

**Summary of Part L.** Quantum information has hard edges. We proved no-cloning from overlap
preservation, measured the $2\sqrt2$ Tsirelson violation of CHSH, bounded retrievable information by
Holevo, and exhibited the all-or-nothing GHZ refutation of local realism. The impossibilities and
the excesses are two sides of the same non-classical structure.

---

## Part M. Coda: the bosonic limit (Fock space and second quantization)

The discrete state space has a natural continuation: let the dimension grow without bound and you
reach the harmonic oscillator, the quantum of a field mode. Its Fock states $|n\rangle$ are still a
countable discrete basis, so the tools transfer, but new objects appear that have no
finite-dimensional analogue: coherent states, squeezing, and the photon statistics of light. We
work in a Fock space truncated to $n = 12$ levels, where the annihilation operator is just the
super-diagonal of $\sqrt{k}$ values.

```wolfram
n = 12;
a = DiagonalMatrix[Sqrt@Range[n - 1], 1];   (* annihilation: a|k> = Sqrt[k] |k-1> *)
fock[k_] := UnitVector[n, k + 1];
```

#### M1. How do I build Fock states and the ladder operators?

The Fock state $|n\rangle$ has exactly $n$ quanta; the annihilation operator lowers it,
$a|n\rangle = \sqrt n\,|n-1\rangle$, and the number operator $N = a^\dagger a$ counts,
$N|n\rangle = n|n\rangle$.

```wolfram
Chop[a . fock[2] - Sqrt[2] fock[1]]                       (* a|2> = Sqrt[2] |1> *)
Chop[(ConjugateTranspose[a] . a) . fock[3] - 3 fock[3]]   (* N|3> = 3|3> *)
```

The ladder lowers $|2\rangle$ to $\sqrt2\,|1\rangle$ and the number operator returns the eigenvalue
$3$ on $|3\rangle$. These two relations generate the entire oscillator algebra.

#### M2. How do I build a coherent state and check it is an eigenstate of $a$?

A coherent state $|\alpha\rangle = D(\alpha)|0\rangle$ is the displaced vacuum, the most classical
state of light. It is the eigenstate of annihilation, $a|\alpha\rangle = \alpha|\alpha\rangle$, with
mean photon number $\langle N\rangle = |\alpha|^2$.

```wolfram
al = 1.2;
disp[\[Alpha]_] := MatrixExp[\[Alpha] ConjugateTranspose[a] - Conjugate[\[Alpha]] a];   (* D(\[Alpha]) *)
coh = disp[al] . fock[0];
Chop[(a . coh - al coh)][[1 ;; 6]]                          (* eigenstate of a *)
Chop[Re[Conjugate[coh] . (ConjugateTranspose[a] . a) . coh]]   (* mean photon |al|^2 *)
```

Annihilation returns $\alpha$ times the state (zero difference on the low components, up to
truncation), and the mean photon number is $1.44 = (1.2)^2 = |\alpha|^2$. The coherent state stays
coherent under photon loss, which is why a laser beam is coherent.

#### M3. How do I show the coherent-state photon distribution is Poissonian?

The probability of finding $n$ photons in $|\alpha\rangle$ is
$|\langle n|\alpha\rangle|^2 = e^{-|\alpha|^2}|\alpha|^{2n}/n!$, a Poisson distribution of mean
$|\alpha|^2$. Compare the amplitudes against the Poisson formula.

```wolfram
Chop[Abs[coh]^2][[1 ;; 6]]
Table[N@PDF[PoissonDistribution[al^2], k], {k, 0, 5}]
```

The photon-number probabilities match the Poisson distribution to all shown digits. Poissonian
counting statistics, with variance equal to the mean, are the signature of classical-like coherent
light.

#### M4. How does squeezing push a quadrature below the vacuum noise?

The quadratures $X$ and $P$ obey an uncertainty relation, and the vacuum spreads its noise equally
between them. A squeezed state redistributes that noise, pushing one quadrature variance below the
vacuum floor at the expense of the other.

```wolfram
squeeze[\[Xi]_] := MatrixExp[1/2 (Conjugate[\[Xi]] a . a - \[Xi] ConjugateTranspose[a] . ConjugateTranspose[a])];
quadX = (a + ConjugateTranspose[a])/Sqrt[2];
var[\[Psi]_, op_] := Re[Conjugate[\[Psi]] . op . op . \[Psi] - (Conjugate[\[Psi]] . op . \[Psi])^2];
var[fock[0], quadX]                  (* vacuum variance 1/2 *)
var[squeeze[1.0] . fock[0], quadX]   (* squeezed < 1/2 *)
```

The vacuum quadrature variance is $1/2$; squeezing drops it to $0.26$, well below the vacuum level.
Sub-vacuum noise in one quadrature is a genuinely non-classical resource, used to push
interferometers like gravitational-wave detectors past the standard quantum limit.

#### M5. How do I distinguish light by its photon statistics ($g^{(2)}(0)$)?

The second-order coherence $g^{(2)}(0) = \langle a^\dagger a^\dagger a a\rangle /
\langle a^\dagger a\rangle^2$ classifies light: thermal light bunches ($g^{(2)} = 2$), coherent
light is Poissonian ($g^{(2)} = 1$), and a single photon antibunches ($g^{(2)} = 0$), which no
classical field can do.

```wolfram
ClearAll[proj, g2];          (* proj held the C1 projector list above; free it *)
g2[\[Rho]_] := Re[Tr[\[Rho] . ConjugateTranspose[a] . ConjugateTranspose[a] . a . a]]/Re[Tr[\[Rho] . ConjugateTranspose[a] . a]]^2;
proj[\[Psi]_] := KroneckerProduct[\[Psi], Conjugate@\[Psi]];
thermal = With[{x = 0.5}, (1 - x) DiagonalMatrix[x^Range[0, n - 1]]]; thermal = thermal/Tr[thermal];
{g2[thermal], g2[proj[coh]], g2[proj[fock[1]]]}
```

Thermal light gives $g^{(2)} \approx 2$ (photon bunching), the coherent state gives $1$ (Poissonian),
and the single-photon Fock state gives $0$ (perfect antibunching). The value $g^{(2)}(0) < 1$ is the
operational definition of non-classical light and the standard test of a single-photon source.

**Summary of Part M.** Letting the dimension grow yields the harmonic oscillator, a discrete Fock
basis with new structure. We built the ladder operators, the coherent state and its Poisson
statistics, sub-vacuum squeezing, and the $g^{(2)}(0)$ that separates classical from non-classical
light. The discrete toolkit extends continuously into quantum optics.

---

## Where this leaves us

Working top to bottom, we have built every standard object of finite-dimensional quantum theory and
made each one compute: states and the Born rule, observables and measurement, entanglement and its
measures, mixed states and channels, circuits and the stabilizer formalism, open-system dynamics,
the discrete Wigner function, the canonical algorithms, the information-theoretic limits, and the
bosonic continuation. The recurring lesson is one of altitude: a physics question deserves the most
direct computation, and for the single qubit, its operators, its density matrix, and its
open-system dynamics, that is plain Wolfram Language on matrices and vectors. A state is a vector, an
operator a matrix, an expectation value a dot product, a channel a Kraus sum, the Lindblad equation
a one-line `DSolveValue`, and a claim is settled by `==` under `FullSimplify`.

QuantumFramework earned its place exactly where that direct route stops being the shortest one: the
multi-qubit stabilizer tableaux and error-correcting codes of Part I, the discrete Wigner and mana
of Part J, the circuit-structure operations (depth, layers, the KAK decomposition) of Part F, and
the many-wire circuit protocols (Simon, teleportation, superdense coding, entanglement swapping) of
Part K. The skill is not choosing one tool over the other but reaching for whichever expresses the
physics in the fewest, clearest symbols.

From here the natural next steps split three ways. Toward hardware: compile the circuits to a real
device with QuantumFramework's `QuantumQASM` and the IBM or AWS backends. Toward scale: route large
circuits through its tensor-network and stabilizer engines that keep them tractable. Toward depth:
push any single part into its own study, whether fault-tolerant codes, variational chemistry, or
continuous-variable quantum optics. Every one of those continuations is the same move repeated:
state the physics, then compute it.
