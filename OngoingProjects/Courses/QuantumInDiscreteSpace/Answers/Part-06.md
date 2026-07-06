---
Template: Default
---

# Quantum in Finite Dimensions: Answers in WL and in QuantumFramework (Part 6)

A companion answer key to `Question-List.md`. For each question it gives **two** worked answers:

1. **WL** : native Wolfram Language only (plain vectors and matrices, `PauliMatrix`, `Conjugate`,
   `Normalize`, `Eigenvectors`, `MatrixExp`, and friends). Nothing beyond the built-in language.
2. **QF** : the same task done through the **QuantumFramework** paclet, leaning on its objects and
   property downvalues (`QuantumState[...]["..."]`) rather than rebuilding matrices by hand. A QF
   result is a rich object, a state or an operator, not a bare array: it can act on a target, compose
   with others, report its spectrum, change basis, and more. The answers keep results as objects and
   treat a matrix or amplitude vector as just one of the object's many representations, not the thing
   itself.

Three habits run through the answers. They stay **symbolic** wherever possible (general amplitudes
$\alpha,\beta$, general angles $\theta,\phi$), so each operation is exercised in full generality
rather than on a lucky special case. They write the computation **explicitly** rather than hiding it
in a helper function, because the code is itself part of the explanation. And they prefer **built-in
features** (a `Normalize`, a named state) over hand-rolled constructions, so nothing is hardcoded:
every cell computes its result. Each cell does **one** thing and is preceded by a sentence saying
what, so the notebook reads as a sequence of single, captioned computations. Run the WL answers in a
bare kernel; the QF answers need the paclet loaded once.

## Setup

The QF answers load the paclet once:

```wl
Needs["Wolfram`QuantumFramework`"]
```

---

## Part 6. Uncertainty and incompatibility

### 6.1 [MSc] How do I verify the Robertson-Schrodinger uncertainty relation for a given state?

For observables $A, B$ on *any* state $\rho$ (pure or mixed), with standard deviations
$\sigma_A = \sqrt{\langle A^2\rangle - \langle A\rangle^2}$ and $\sigma_B$ likewise, the Robertson-Schrodinger
relation bounds the product of variances,
$$
\sigma_A^2\,\sigma_B^2 \;\ge\; \Big(\tfrac12\langle\{A,B\}\rangle - \langle A\rangle\langle B\rangle\Big)^2
  + \Big(\tfrac1{2i}\langle[A,B]\rangle\Big)^2 ,
$$
the first term a covariance, the second the Robertson (commutator) piece, with $\langle\cdot\rangle =
\mathrm{Tr}[\cdot\,\rho]$. Every qubit observable is $a_0 I + a\,\hat n\cdot\vec\sigma$; the identity shift
drops out of every variance, covariance, and commutator, and the two scales multiply both sides of the
relation equally, so the fully general pair is a pair of unit directions, $A = \hat n_1\cdot\vec\sigma$ and
$B = \hat n_2\cdot\vec\sigma$, each built from its spherical angles, on the general qubit
$\rho = \tfrac12(I + \vec r\cdot\vec\sigma)$, $|\vec r|\le1$. The telling quantity is the slack, how far
the state sits above the bound.

**WL** : build the two directions with `FromSphericalCoordinates`, then, with $\langle\cdot\rangle$ the
trace against $\rho$, compute the left side minus the right side.

```wl
{n1, n2} = FromSphericalCoordinates /@ {{1, \[Theta]1, \[Phi]1}, {1, \[Theta]2, \[Phi]2}};
slack = With[{a = n1 . PauliMatrix[{1, 2, 3}], b = n2 . PauliMatrix[{1, 2, 3}], rho = 1/2 (IdentityMatrix[2] + {rx, ry, rz} . PauliMatrix[{1, 2, 3}])},
  With[{ev = Tr[# . rho] &},
   FullSimplify[(ev[a . a] - ev[a]^2) (ev[b . b] - ev[b]^2) -
     ((ev[(a . b + b . a)/2] - ev[a] ev[b])^2 + (ev[a . b - b . a]/(2 I))^2)]]]
```

The relation itself is the claim that this slack is never negative. So confirm the inequality first: ask
the kernel for any physical Bloch vector ($|\vec r|\le1$) and any pair of directions that violate it.

```wl
Reduce[slack < 0 && rx^2 + ry^2 + rz^2 <= 1, {rx, ry, rz, \[Theta]1, \[Phi]1, \[Theta]2, \[Phi]2}, Reals]
```

The violating set is empty, so the relation holds for every qubit state and every observable pair at once.
To see why, look at the slack itself: one factor is the mixedness $1 - |\vec r|^2$, the other is
trigonometric in the four angles. Since the qubit commutator is geometric,
$[\hat n_1\cdot\vec\sigma, \hat n_2\cdot\vec\sigma] = 2i\,(\hat n_1\times\hat n_2)\cdot\vec\sigma$, test
whether that factor is the squared cross product.

```wl
FullSimplify[slack == Cross[n1, n2] . Cross[n1, n2] (1 - {rx, ry, rz} . {rx, ry, rz})]
```

**QF** : the same slack, the general state as a `QuantumState` and its density matrix, the observables as
`QuantumOperator`s built from the same spherical directions, the commutator via `Commutator`.

```wl
With[{qs = QuantumState[1/2 (IdentityMatrix[2] + {rx, ry, rz} . PauliMatrix[{1, 2, 3}])],
  a = QuantumOperator[FromSphericalCoordinates[{1, \[Theta]1, \[Phi]1}] . PauliMatrix[{1, 2, 3}]],
  b = QuantumOperator[FromSphericalCoordinates[{1, \[Theta]2, \[Phi]2}] . PauliMatrix[{1, 2, 3}]]},
 With[{ev = Tr[Normal[#["Matrix"]] . Normal[qs["DensityMatrix"]]] &},
  FullSimplify[(ev[a @ a] - ev[a]^2) (ev[b @ b] - ev[b]^2) -
    ((ev[(a @ b + b @ a)/2] - ev[a] ev[b])^2 + (ev[Commutator[a, b]]/(2 I))^2)]]]
```

Both compute the identical factored slack; the violation search returns False and the equality test
returns True: the slack is exactly $|\hat n_1\times\hat n_2|^2\,(1 - |\vec r|^2)$, the observables'
incompatibility ($|\hat n_1\times\hat n_2|^2 = \sin^2\gamma$, with $\gamma$ the angle between the two
directions) times the state's mixedness ($1 - |\vec r|^2 = 2(1 - \mathrm{Tr}[\rho^2])$). A product of two
non-negative factors, which is the structural reason no violation exists, and it saturates exactly when a
factor vanishes: a pure state ($|\vec r| = 1$) or a commuting pair ($\hat n_1\times\hat n_2 = 0$, the
directions parallel or antiparallel). The perpendicular choice $A = X$,
$B = Z$ has $|\hat n_1\times\hat n_2| = 1$, so its slack is the mixedness alone. The commutator term alone
is the weaker Robertson relation $\sigma_A\sigma_B \ge \tfrac12|\langle[A,B]\rangle|$; Schrodinger's
covariance term is what tightens it.

### 6.2 [MSc] How do I find a minimum-uncertainty state that saturates the Robertson bound?

The Robertson bound $\sigma_A\sigma_B \ge \tfrac12|\langle[A,B]\rangle|$ is saturated exactly where the
*Robertson gap* $g = \sigma_A^2\sigma_B^2 - \big(\tfrac1{2i}\langle[A,B]\rangle\big)^2$ vanishes. Take the
generic observable pair of 6.1, $A = \hat n_1\cdot\vec\sigma$ and
$B = \hat n_2\cdot\vec\sigma$, compute $g$ over the whole pure-qubit family
$\{\cos\tfrac\theta2, e^{i\phi}\sin\tfrac\theta2\}$, and solve $g = 0$. Purity is no restriction here: by
6.1 the slack of the full relation is $|\hat n_1\times\hat n_2|^2(1 - |\vec r|^2)$, so for a non-commuting
pair a mixed state sits strictly above even the stronger Robertson-Schrodinger bound, and every
Robertson-saturating state must be pure.

**WL** : the gap over the whole pure family for the generic pair; the expectation is the sandwich
$\langle\psi|\cdot|\psi\rangle$, and `ComplexExpand` keeps it explicitly real (all angles are real).

```wl
{n1, n2} = FromSphericalCoordinates /@ {{1, \[Theta]1, \[Phi]1}, {1, \[Theta]2, \[Phi]2}};
gap = With[{psi = {Cos[\[Theta]/2], Exp[I \[Phi]] Sin[\[Theta]/2]}, a = n1 . PauliMatrix[{1, 2, 3}], b = n2 . PauliMatrix[{1, 2, 3}]},
  With[{ev = Simplify[ComplexExpand[Conjugate[psi] . # . psi]] &},
   (ev[a . a] - ev[a]^2) (ev[b . b] - ev[b]^2) - (ev[a . b - b . a]/(2 I))^2]]
```

The raw gap is a tangle of trigonometry in six angles; identify it. For a pure state the
Robertson-Schrodinger slack of 6.1 is zero, so the gap should be exactly the squared covariance: with
$\hat m$ the state's Bloch vector, $\tfrac12\langle\{A,B\}\rangle - \langle A\rangle\langle B\rangle =
\hat n_1\cdot\hat n_2 - (\hat m\cdot\hat n_1)(\hat m\cdot\hat n_2)$. Test the identity: the difference
should expand to exactly zero.

```wl
Expand[TrigExpand[gap - With[{m = FromSphericalCoordinates[{1, \[Theta], \[Phi]}]}, (n1 . n2 - (m . n1) (m . n2))^2]]]
```

So $g = \big(\hat n_1\cdot\hat n_2 - (\hat m\cdot\hat n_1)(\hat m\cdot\hat n_2)\big)^2$, and the
minimum-uncertainty states are the zero-covariance states, the solutions of
$(\hat m\cdot\hat n_1)(\hat m\cdot\hat n_2) = \hat n_1\cdot\hat n_2$. The eigenstates
$\hat m = \pm\hat n_1, \pm\hat n_2$ always solve it (the two signs cancel and one variance vanishes). For
the mutually unbiased pair $A = X$, $B = Z$ ($\hat n_1\cdot\hat n_2 = 0$), solve the condition on the sphere.

```wl
Reduce[(With[{m = FromSphericalCoordinates[{1, \[Theta], \[Phi]}]}, (m . n1) (m . n2) == n1 . n2] /. {\[Theta]1 -> Pi/2, \[Phi]1 -> 0, \[Theta]2 -> 0}) && 0 < \[Theta] < Pi && -Pi < \[Phi] <= Pi, {\[Theta], \[Phi]}]
```

One such state sits where the two solution branches cross, $\theta = \phi = \tfrac\pi2$.

```wl
{Cos[\[Theta]/2], Exp[I \[Phi]] Sin[\[Theta]/2]} /. {\[Theta] -> Pi/2, \[Phi] -> Pi/2} // FullSimplify
```

**QF** : the same gap through the paclet's objects (state, operators, `Commutator`); it is the
`===`-identical expression.

```wl
With[{psi = QuantumState[{Cos[\[Theta]/2], Exp[I \[Phi]] Sin[\[Theta]/2]}],
  a = QuantumOperator[FromSphericalCoordinates[{1, \[Theta]1, \[Phi]1}] . PauliMatrix[{1, 2, 3}]],
  b = QuantumOperator[FromSphericalCoordinates[{1, \[Theta]2, \[Phi]2}] . PauliMatrix[{1, 2, 3}]]},
 With[{ev = Simplify[ComplexExpand[(psi["Dagger"] @ # @ psi)["Scalar"]]] &},
  (ev[a @ a] - ev[a]^2) (ev[b @ b] - ev[b]^2) - (ev[Commutator[a, b]]/(2 I))^2]] === gap
```

For a qubit, then, *every* pure state already saturates the full Robertson-Schrodinger relation, and the
plain Robertson bound is additionally reached exactly on the zero-covariance set. For $X, Z$ the solve
returns that set as the equator $\theta = \tfrac\pi2$ ($\langle Z\rangle = 0$) and the meridians
$\phi = \pm\tfrac\pi2$ ($\langle X\rangle = 0$); at their crossing both means vanish and the state is the
$Y$ eigenstate $\{1, i\}/\sqrt2$, the most symmetric minimum-uncertainty state.

### 6.3 [BSc] How do I compute the classical Shannon entropy $H(p) = -\sum_i p_i\log_2 p_i$ of a measurement distribution?

A measurement yields a probability distribution; its Shannon entropy $H(p) = -\sum_i p_i\log_2 p_i$ (in
bits) measures the outcome uncertainty. Take the $Z$ distribution of a state with populations
$\{\tfrac34, \tfrac14\}$.

**WL** : the entropy is the negative weighted sum of $\log_2$ probabilities.

```wl
-Total[# Log2[#] & /@ {3/4, 1/4}]
```

**QF** : a `QuantumMeasurement` reports its `"Entropy"` directly, as a `Quantity` in bits.

```wl
QuantumMeasurementOperator[QuantumOperator["Z"]][QuantumState[{Sqrt[3]/2, 1/2}]]["Entropy"]
```

Both give $-\tfrac34\log_2\tfrac34 - \tfrac14\log_2\tfrac14 \approx 0.811$ bits: a biased two-outcome
distribution carries less than one bit, reaching exactly one bit only for the uniform $\{\tfrac12,\tfrac12\}$.

### 6.4 [MSc] How do I construct a maximal set of mutually unbiased bases for a qubit?

Two bases are mutually unbiased when every cross-overlap is flat, $|\langle e|f\rangle|^2 = 1/d$: a definite
state in one is maximally uncertain in the other. A set of bases is mutually unbiased when every pair in it
is, and *maximal* when no further basis can be added that stays unbiased to every member. The construction
has three steps: derive the form of a basis unbiased to the computational one, solve for the phase
separation that makes two such families mutually unbiased, then show the set cannot be extended.

**WL** : a second basis $\{u, u_\perp\}$ is unbiased to the computational one when
$|\langle 0|u\rangle|^2 = |\langle 1|u\rangle|^2 = \tfrac12$. Write the candidate amplitudes in polar form,
$u = \{r_a e^{i\delta_a}, r_b e^{i\delta_b}\}$, and solve the two conditions.

```wl
Reduce[Abs[ra Exp[I \[Delta]a]]^2 == 1/2 && Abs[rb Exp[I \[Delta]b]]^2 == 1/2 && ra >= 0 && rb >= 0, {ra, rb, \[Delta]a, \[Delta]b}, Reals]
```

Both moduli are forced to $\tfrac1{\sqrt2}$ and neither phase is constrained. Discarding the physically
irrelevant global phase (Part 1) leaves the one-parameter family $u = \{1, e^{i\alpha}\}/\sqrt2$, and
orthogonality then fixes its partner to $u_\perp = \{1, -e^{i\alpha}\}/\sqrt2$ up to phase. Two such
families $\alpha$ and $\beta$ are themselves mutually unbiased only at a special separation; only the
difference matters, so fix $\alpha = 0$ and solve for $\beta$.

```wl
With[{ov = FullSimplify[Abs[Conjugate[{1, 1}/Sqrt[2]] . ({1, Exp[I \[Beta]]}/Sqrt[2])]^2, \[Beta] \[Element] Reals]},
 Reduce[ov == 1/2 && 0 <= \[Beta] < 2 Pi, \[Beta]]]
```

The two solutions $\beta = \tfrac\pi2$ and $\tfrac{3\pi}2$ are one and the same basis with its two vectors
exchanged ($e^{3i\pi/2} = -e^{i\pi/2}$), so exactly one new basis appears and the set holds three:
computational, $\alpha = 0$, $\alpha = \tfrac\pi2$. Maximality is the claim that no fourth basis can join.
A fourth basis would again be unbiased to the computational one, hence of the same family form at some
phase $\gamma$, and simultaneously unbiased to both the $\alpha = 0$ and $\alpha = \tfrac\pi2$ families;
solve for $\gamma$.

```wl
With[{ov = Abs[Conjugate[{1, Exp[I #1]}/Sqrt[2]] . ({1, Exp[I #2]}/Sqrt[2])]^2 &},
 Reduce[FullSimplify[ov[0, \[Gamma]] == 1/2 && ov[Pi/2, \[Gamma]] == 1/2, \[Gamma] \[Element] Reals] && 0 <= \[Gamma] < 2 Pi, \[Gamma]]]
```

No phase survives both conditions: the three bases cannot be extended, and the set is maximal. Assemble it.

```wl
With[{mub = Function[a, {{1, Exp[I a]}, {1, -Exp[I a]}}/Sqrt[2]]},
 Prepend[mub /@ {0, Pi/2}, IdentityMatrix[2]]]
```

**QF** : those three constructed bases are (up to phase) the eigenbases of $Z, X, Y$; read them off the
operators.

```wl
Function[A, Normal /@ QuantumOperator[A]["Eigenvectors"]] /@ {"Z", "X", "Y"}
```

Solving the unbiasedness conditions delivers three bases and proves no fourth exists: the maximal count for
a qubit is $d+1 = 3$ (the qubit case of the general prime-power-dimension result), and the bases are the
computational ($Z$) basis and the phase families at $\alpha = 0$ and $\tfrac\pi2$, which are the $X$ and $Y$
eigenbases.

### 6.4b [BSc] How do I verify that given bases are pairwise mutually unbiased, i.e. every cross-overlap satisfies $|\langle e_i|f_j\rangle|^2 = 1/d$?

Given candidate bases $\{|e_i\rangle\}$, $\{|f_j\rangle\}, \ldots$, the test is the definition itself: for
every pair of *different* bases, *every* cross-overlap between their members must equal
$|\langle e_i|f_j\rangle|^2 = 1/d$; a single deviating entry disqualifies the set. Run the test on the three
bases constructed in 6.4 in their operator form, the eigenbases of $X, Y, Z$ (here $d = 2$, so the target
value is $\tfrac12$).

**WL** : form the three eigenbases and compute the full overlap matrix for each of the three pairs.

```wl
ovs = With[{bx = Normalize /@ Eigenvectors[PauliMatrix[1]], by = Normalize /@ Eigenvectors[PauliMatrix[2]], bz = Normalize /@ Eigenvectors[PauliMatrix[3]]},
  Outer[Abs[Conjugate[#1] . #2]^2 &, #[[1]], #[[2]], 1] & /@ {{bx, by}, {bx, bz}, {by, bz}}]
```

The verdict collapses all twelve cross-overlaps to their set of distinct values; the test passes exactly
when that set is $\{1/d\}$.

```wl
Union[Flatten[ovs]]
```

**QF** : the same test from the operators' `"Eigenvectors"`, over all pairs at once.

```wl
Union @ Flatten @ With[{b = Function[A, Normal /@ QuantumOperator[A]["Eigenvectors"]] /@ {"X", "Y", "Z"}},
   Outer[Abs[Conjugate[#1] . #2]^2 &, #[[1]], #[[2]], 1] & /@ Subsets[b, {2}]]
```

Every cross-overlap between different bases collapses to the single value $\tfrac12 = 1/d$: the $X, Y, Z$
eigenbases are pairwise mutually unbiased, confirming they are the three MUBs of a qubit. The collapse to a
one-element set is what gives the test teeth: any deviating overlap would survive as a second element.

### 6.5 [MSc] How do I verify the Maassen-Uffink entropic uncertainty relation $H(A)+H(B)\ge-\log_2 c^2$, with $c$ the largest overlap between eigenstates of $A$ and $B$?

Let $H(A)$ denote the Shannon entropy (6.3) of the outcome distribution of measuring the observable $A$,
and let $\{|a_j\rangle\}$ and $\{|b_k\rangle\}$ be the eigenstates of $A$ and of $B$. The Maassen-Uffink
relation bounds the summed entropies for *every* state, $H(A) + H(B) \ge -\log_2 c^2$, where
$c = \max_{j,k}|\langle a_j|b_k\rangle|$ is the largest overlap between the two eigenbases; the bound
involves only the observables, not the state. The pair of observables is itself kept general: as in 6.1,
identity shifts and scales drop out, and a global rotation changes no entropy and no overlap, so a general
qubit pair reduces to $B = Z$ and $A = \hat n_\gamma\cdot\vec\sigma$ with
$\hat n_\gamma = (\sin\gamma, 0, \cos\gamma)$; the one genuine invariant is the angle $\gamma$ between the
two Bloch axes, folded into $[0, \tfrac\pi2]$ since $\hat n \to -\hat n$ only relabels outcomes. Verify the
relation for every pair and every state at once: derive $c(\gamma)$ and the bound, find which $\gamma$
carries the strongest bound, then minimize the entropy sum minus the bound over the state and the pair
together. The minimum over mixed states lies on the pure sphere: scaling a Bloch vector out to the sphere
grows the projection on both axes at once, lowering both entropies.

**WL** : every cross-overlap between the two eigenbases at once, from the rank-one projector identity
$|\langle a_j|b_k\rangle|^2 = \mathrm{Tr}[P^A_j P^B_k]$, with $P_\pm = \tfrac12(I \pm A)$ the spectral
projectors of a $\pm1$ observable.

```wl
FullSimplify[Table[Tr[((IdentityMatrix[2] + s1 {Sin[\[Gamma]], 0, Cos[\[Gamma]]} . PauliMatrix[{1, 2, 3}])/2) . ((IdentityMatrix[2] + s2 PauliMatrix[3])/2)], {s1, {1, -1}}, {s2, {1, -1}}], 0 <= \[Gamma] <= Pi/2]
```

The four overlaps are $\cos^2\tfrac\gamma2$ and $\sin^2\tfrac\gamma2$; on $[0, \tfrac\pi2]$ the largest is
$c^2 = \cos^2\tfrac\gamma2$, so the bound is $-\log_2\cos^2\tfrac\gamma2$. Ask which pair carries the
strongest bound.

```wl
Maximize[{-Log2[Cos[\[Gamma]/2]^2], 0 <= \[Gamma] <= Pi/2}, \[Gamma]]
```

The bound is strongest, one full bit, exactly at $\gamma = \tfrac\pi2$: the mutually unbiased pair of 6.4b,
which is $X, Z$ up to rotation. A commuting pair ($\gamma = 0$) has $c = 1$ and bound $0$, a vacuous
statement. Now build the entropy sum for the generic pair over the whole pure family: the Born weights of a
$\pm1$ observable are $p_\pm = \tfrac12(1 \pm \langle A\rangle)$, `ComplexExpand` renders each mean as a
manifestly real function of the angles, and the `If` keeps the convention $0\log_2 0 = 0$ explicit.

```wl
hsum = With[{psi = {Cos[\[Theta]/2], Exp[I \[Phi]] Sin[\[Theta]/2]}, a = {Sin[\[Gamma]], 0, Cos[\[Gamma]]} . PauliMatrix[{1, 2, 3}], z = PauliMatrix[3], h = -Total[If[# == 0, 0, # Log2[#]] &[#] & /@ #] &},
  With[{ev = Simplify[ComplexExpand[Conjugate[psi] . # . psi]] &},
   h[{(1 + ev[a])/2, (1 - ev[a])/2}] + h[{(1 + ev[z])/2, (1 - ev[z])/2}]]]
```

Minimize the slack, the entropy sum minus the bound, over the state and the pair together.

```wl
NMinimize[{hsum + Log2[Cos[\[Gamma]/2]^2], 0 <= \[Gamma] <= Pi/2}, {\[Theta], \[Phi], \[Gamma]}]
```

The slack never drops below zero (the minimum is $0$ to machine precision, $\sim 3\times10^{-16}$):
Maassen-Uffink holds for every qubit pair and every state at once. The minimizer sits at
$\gamma = \tfrac\pi2$, so park there, the strongest-bound pair, and minimize the entropy sum itself.

```wl
NMinimize[hsum /. \[Gamma] -> Pi/2, {\[Theta], \[Phi]}]
```

**QF** : the same objective through the paclet's objects, each mean from the object sandwich
`(psi["Dagger"] @ op @ psi)["Scalar"]`; it is the `===`-identical expression, so both scans above apply to
it verbatim.

```wl
With[{psi = QuantumState[{Cos[\[Theta]/2], Exp[I \[Phi]] Sin[\[Theta]/2]}],
  a = QuantumOperator[{Sin[\[Gamma]], 0, Cos[\[Gamma]]} . PauliMatrix[{1, 2, 3}]], z = QuantumOperator["Z"],
  h = -Total[If[# == 0, 0, # Log2[#]] &[#] & /@ #] &},
 With[{ev = Simplify[ComplexExpand[(psi["Dagger"] @ # @ psi)["Scalar"]]] &},
  h[{(1 + ev[a])/2, (1 - ev[a])/2}] + h[{(1 + ev[z])/2, (1 - ev[z])/2}]]] === hsum
```

At the maximally incompatible pair the minimum of $H(X) + H(Z)$ is exactly the computed bound, $1$ bit,
attained at $\theta \approx 0$, a $Z$ eigenstate: no qubit state is sharp in both $X$ and $Z$, and the total
uncertainty of a mutually unbiased pair can never drop below one full bit. The $X, Z$ of 6.6 is therefore
not a convenience choice but the extremal pair the `Maximize` step singled out; which states attain its
bound is 6.6's question.

### 6.6 [MSc] How do I find the states that saturate the Maassen-Uffink entropic bound $H(X)+H(Z)\ge 1$ bit for the complementary $X,Z$ (maximal complementarity)?

The bound of 6.5, $H(X) + H(Z) \ge 1$ bit for the mutually unbiased $X, Z$ (the extremal pair 6.5's
`Maximize` singled out), is attained: 6.5's minimization bottomed out exactly there. Find *which* states
attain it. Work on the pure sphere (by 6.5's argument the minimum lives there): locate the candidates from
the definition of the entropy, then check saturation exactly.

**WL** : the summed entropy over the whole pure family, the two Born distributions read from the
eigenbases (6.4b).

```wl
hsum = With[{psi = {Cos[\[Theta]/2], Exp[I \[Phi]] Sin[\[Theta]/2]}, h = -Total[If[# == 0, 0, # Log2[#]] &[#] & /@ #] &},
  h[Simplify[ComplexExpand[Abs[Conjugate[#] . psi]^2]] & /@ (Normalize /@ Eigenvectors[PauliMatrix[1]])] +
   h[Simplify[ComplexExpand[Abs[Conjugate[#] . psi]^2]] & /@ (Normalize /@ Eigenvectors[PauliMatrix[3]])]]
```

The landscape over the sphere shows a floor at one bit, touched only at isolated states.

```wl
Plot3D[hsum, {\[Theta], 0, Pi}, {\[Phi], -Pi, Pi}, PlotRange -> All, AxesLabel -> {"\[Theta]", "\[Phi]", "H(X)+H(Z)"}]
```

Locate those states exactly. A Shannon entropy vanishes only on a point-mass distribution ($-\sum_i p_i
\log_2 p_i = 0$ forces every $p_i \in \{0, 1\}$), so the candidates are the states *sharp* in $Z$ or in
$X$: for each observable, solve for where the product of its two outcome probabilities vanishes.

```wl
With[{psi = {Cos[\[Theta]/2], Exp[I \[Phi]] Sin[\[Theta]/2]}},
 Function[A, Reduce[Simplify[ComplexExpand[Times @@ (Abs[Conjugate[#] . psi]^2 & /@ (Normalize /@ Eigenvectors[PauliMatrix[A]]))]] == 0 &&
     0 <= \[Theta] <= Pi && -Pi < \[Phi] <= Pi, {\[Theta], \[Phi]}]] /@ {3, 1}]
```

The sharp states are the $Z$ eigenstates (the poles $\theta = 0, \pi$, where $\phi$ is a redundant
coordinate) and the $X$ eigenstates (on the equator at $\phi = 0, \pi$). Check saturation exactly at all of
them.

```wl
FullSimplify[{hsum /. \[Theta] -> 0, hsum /. \[Theta] -> Pi, hsum /. {\[Theta] -> Pi/2, \[Phi] -> 0}, hsum /. {\[Theta] -> Pi/2, \[Phi] -> Pi}}, \[Phi] \[Element] Reals]
```

**QF** : the same four states by name; each measurement reports its `"Entropy"` (as in 6.3), and each sum
lands exactly on the bound.

```wl
With[{sum = QuantityMagnitude[QuantumMeasurementOperator[QuantumOperator["X"]][#]["Entropy"]] +
    QuantityMagnitude[QuantumMeasurementOperator[QuantumOperator["Z"]][#]["Entropy"]] &},
 sum /@ (QuantumState /@ {"0", "1", "Plus", "Minus"})]
```

The saturating states are exactly the eigenstates of one of the two observables: sharp in $Z$ at the cost
of total ignorance of $X$ ($H(Z) = 0$, $H(X) = 1$), or the reverse. Mutually unbiased bases realize maximal
complementarity, and the entropic bound is met with equality only there; every other state pays strictly
more total uncertainty, as the landscape shows.
