# Part 4 Plan: The harmonic oscillator

8 questions. Class census, per the `Route-Table.md` class census: C1 = {4.1, 4.5}, both entries
inheriting the C1 verdict's binding traps in their prose; C0 = {4.2, 4.3, 4.4, 4.6, 4.7, 4.8}, WL
machinery named directly per the C0 row.

## Common ground

The part rests on: the oscillator Hamiltonian $\hat H=\tfrac12(\hat p^2+\hat x^2)$ with
$\hat p=-i\,d/dx$ (natural units $\hbar=m=\omega=1$ throughout; only 4.5 restores units); the
ladder pair $a=(\hat x+i\hat p)/\sqrt2$, $a^\dagger=(\hat x-i\hat p)/\sqrt2$ with
$[a,a^\dagger]=1$, the number operator $\hat n=a^\dagger a$, and $\hat H=\hat n+\tfrac12$; the
ladder action $a\vert n\rangle=\sqrt n\,\vert n-1\rangle$,
$a^\dagger\vert n\rangle=\sqrt{n+1}\,\vert n+1\rangle$; the normalized eigenfunctions
$\psi_n(x)=H_n(x)\,e^{-x^2/2}/\sqrt{2^n\,n!\,\sqrt\pi}$ with $E_n=n+\tfrac12$; and the vacuum
quadrature variances $\Delta x^2=\Delta p^2=\tfrac12$, the reference disk against which the
coherent and squeezed states of 4.6 through 4.8 are measured.

### 4.1 [BSc] How do I solve the oscillator's stationary equation and obtain the Hermite-function eigenstates with energies $E_n=n+\tfrac12$?

Bind the Hamiltonian as an operator on an expression, then solve at symbolic energy so the general
solution still carries its constants:

```wl
h[f_] := -D[f, {x, 2}]/2 + x^2 f/2;
generic = DSolveValue[h[u[x]] == en u[x], u, x]
```

Returns a `Function` with body
`C[2] ParabolicCylinderD[(-1-2 en)/2, I Sqrt[2] x] + C[1] ParabolicCylinderD[(-1+2 en)/2, Sqrt[2] x]`:
two constants, every energy still admitted, nothing quantized yet. Do not attach decay conditions at
$\pm\infty$ to force the issue; ask what happens rather than assuming a hang:

```wl
AbsoluteTiming[DSolve[{h[u[x]] == en u[x], u[-Infinity] == 0, u[Infinity] == 0}, u, x]]
```

Returns `{27.8, {}}`, an empty solution set rather than a hang, with `DSolveValue` on the same system
returning unevaluated; the cost is paid once per kernel and `ClearSystemCache[]` restores it. (The
hang recorded in the C1 verdict belongs to the Coulomb radial problem, a different equation.) Decay
has to enter through the eigensolver's domain instead, and quantization arrives there as a return
value:

```wl
{energies, eigenfunctions} = DEigensystem[h[u[x]], u[x], {x, -Infinity, Infinity}, 8]
```

Returns `{1/2, 3/2, 5/2, 7/2, 9/2, 11/2, 13/2, 15/2}` with $H_n(x)e^{-x^2/2}$. With no boundary
condition supplied the eigensolver imposes Neumann-zero on the domain boundary, harmless on the whole
line but not the same as decay, so this phrasing does not transfer to a finite interval without an
explicit condition. Locate the free constants rather than counting them:

```wl
constants = DeleteDuplicates@Cases[generic[x], C[_], Infinity]
```

Returns `{C[2], C[1]}`, traversal order rather than naming order, and the names follow the
`GeneratedParameters -> C` default. Now confront the two families so the reduction is computed rather
than chosen. Note the two evaluation idioms inside one equation: `generic` is a `Function` so it
takes `generic[0]`, while each eigenfunction is an expression so it takes `/. x -> 0`, and the energy
must be substituted before solving:

```wl
matchingRules = Table[
   Solve[{(generic[0] /. en -> energies[[j]]) == (eigenfunctions[[j]] /. x -> 0),
          (D[generic[x], x] /. {x -> 0, en -> energies[[j]]}) ==
            (D[eigenfunctions[[j]], x] /. x -> 0)},
     constants],
   {j, Length[energies]}]
```

Returns one solution per level, `{{C[2] -> 0, C[1] -> 1}}`, `{{C[2] -> 0, C[1] -> Sqrt[2]}}`,
`{{C[2] -> 0, C[1] -> 2}}`, on through `{{C[2] -> 0, C[1] -> 8 Sqrt[2]}}`, so $C_1=2^{(j-1)/2}$ and
`C[2] -> 0` at every level: the branch growing like $e^{+x^2/2}$ dies by admissibility, and "select
the decaying branch" becomes an output instead of an instruction. Two things this does not show, both
worth denying explicitly. The surviving $C_1$ is not a normalization, since the eigensolver does not
normalize; confirm that rather than take it on trust with
`Table[Integrate[eigenfunctions[[j]]^2, {x, -Infinity, Infinity}], {j, 8}]`, which returns
`{Sqrt[Pi], 2 Sqrt[Pi], 8 Sqrt[Pi], 48 Sqrt[Pi], ...}`, that is $\|f_n\|^2=2^n n!\sqrt\pi$, so $C_1$
is matching an arbitrary scale. And because the energies were substituted in before solving, this
cell cannot return discreteness; it returns branch death at energies already quantized one cell
earlier. Then the spectrum as a law:

```wl
FindSequenceFunction[energies, n]
```

Returns `(-1 + 2 n)/2`. State the offset in the same breath, since the Wolfram index is the quantum
number plus one: $E_n=n+\tfrac12$ for $n=0,1,2,\dots$

Refute on a carrier that is not the ground state. For $\psi_3\propto(2x^3-3x)e^{-x^2/2}$ (three
nodes, at $0$ and $\pm\sqrt{3/2}$; the vibrational ladder),
`FullSimplify[h[psi3] - (7/2) psi3]` returns the exact integer `0`, with `FullSimplify` doing real
work since the raw difference is a nonzero three-term expression. Pair it with
`Asymptotic[ParabolicCylinderD[-1, I Sqrt[2] x], x -> Infinity]`, returning
`-I E^(x^2/2)/(Sqrt[2] x)`, which is why admissibility killed `C[2]`; prefer `Asymptotic` to
`Series`, whose raw output hides the growth factor inside `SeriesData` until `Normal` is applied.
Close on the count: level $n$ carries exactly $n$ nodes, and the returned spectrum is uniformly
spaced by 1.

### 4.2 [BSc] How do I build the ladder operators $a=(\hat x+i\hat p)/\sqrt2$ and $a^\dagger$ as differential operators and verify $[a,a^\dagger]=1$?

Bind the pair as operators on an expression, with $\hat p=-i\,d/dx$ making $a=(x+\partial_x)/\sqrt2$
and $a^\dagger=(x-\partial_x)/\sqrt2$:

```wl
a[f_]  := (x f + D[f, x])/Sqrt[2];
ad[f_] := (x f - D[f, x])/Sqrt[2];
```

A commutator is a statement about operators, so verify it on an unassigned function and never on a
state:

```wl
Simplify[a[ad[f[x]]] - ad[a[f[x]]]]
```

Returns `f[x]`. That one return is the entire content of $[a,a^\dagger]=1$, and it refutes a sign
slip in $\hat p$, which would return $-f[x]$ or an $x$-dependent expression instead. The same
generic-$f$ move ties the algebra to the Hamiltonian rather than asserting the connection: with
`h[f_] := -D[f, {x, 2}]/2 + x^2 f/2;`, `Simplify[ad[a[f[x]]] - (h[f[x]] - f[x]/2)]` returns `0`, so
$\hat n=a^\dagger a=\hat H-\tfrac12$ holds as an operator identity, on any function at all.

Now the ladder factors, at arbitrary level rather than one carrier. With
`psi[n_, xx_] := HermiteH[n, xx] Exp[-xx^2/2]/Sqrt[2^n n! Sqrt[Pi]]`:

```wl
FullSimplify[a[psi[n, x]] - Sqrt[n] psi[n - 1, x],
  Assumptions -> n > 0 && Element[n, Integers]]
FullSimplify[ad[psi[n, x]] - Sqrt[n + 1] psi[n + 1, x],
  Assumptions -> n >= 0 && Element[n, Integers]]
```

Both return `0`, and `ad[a[psi[n, x]]] - n psi[n, x]` does too. The $\sqrt n$ down and $\sqrt{n+1}$
up are therefore a closed result at symbolic $n$, so nothing is generalized by hand and no law is
asserted; the ledger's nodal carrier $\psi_2=(2x^2-1)e^{-x^2/2}/\sqrt{2\sqrt\pi}$ is then read off as
the instance $a\psi_2=\sqrt2\,\psi_1$, $a^\dagger\psi_2=\sqrt3\,\psi_3$. Do not reach for
`FindSequenceFunction` on the ratios here: on the surd sequence $\{1,\sqrt2,\sqrt3,2,\dots\}$ it
returns unevaluated, so symbolic $n$ is what replaces inference in this entry.

Two closing returns. The edge that bounds the ladder below: `Simplify[a[psi[0, x]]]` is exactly `0`
while `ad[psi[0, x]] - psi[1, x]` is also `0`, so the ladder terminates downward and not upward, and
that asymmetry is the whole reason a lowest state exists. Then a normalization discriminator that
comes free from the eigensolver: hand `DEigensystem`'s own eigenfunctions to the same operators and
the ratios return `2 Sqrt[2]` and `3 Sqrt[2]`, not $\sqrt2$ and $\sqrt3$, because that output is
unnormalized with $\|f_n\|^2=2^n n!\sqrt\pi$. The $\sqrt n$ law is a statement about normalized
states, and the mismatch is what says so rather than a warning in prose. Keep this a
generic-function argument throughout; 4.4 owns the truncated matrices and the corner defect that no
finite representation can avoid.

### 4.3 [BSc] How do I generate the whole spectrum algebraically from $a|0\rangle=0$ and the number operator $\hat n=a^\dagger a$?

Solve $a\psi_0=0$ as a first-order ODE with DSolve, normalize, and build the whole ladder up to
$n=4$ by $\psi_{n+1}=a^\dagger\psi_n/\sqrt{n+1}$ with FoldList (algebraic spectrum): the whole
family, not one level. Verify with checks that share no machinery with the first-order
construction: the second-order residual $\hat H\psi_4-\tfrac92\psi_4=0$ under Simplify,
$\hat n\,\psi_n=n\,\psi_n$ across the family with Table, the laddered $\psi_4$ equal to the closed
HermiteH form under FullSimplify, and the orthonormality Gram matrix by Integrate with a numeric
spot check. Descend as well: four applications of $a$ return $\sqrt{4!}\,\psi_0$ and a fifth
annihilates, the ladder floor. Close with the algebraic reading: $E_n=n+\tfrac12$ follows from
$[a,a^\dagger]=1$ and the floor alone, the Hermite functions falling out of the algebra rather
than being put in.

### 4.4 [BSc] How do I compute the matrix elements of $\hat x$ and $\hat p$ in the number basis, giving the truncated banded-matrix representation of the oscillator (with the truncation made explicit)?

Build $a$, $a^\dagger$ at truncation $N=8$ (levels $0$ through $8$, $9\times9$) with SparseArray
and Band, form $\hat x,\hat p$ by matrix algebra with entries kept as exact surds,
$\langle m|\hat x|n\rangle=(\sqrt n\,\delta_{m,n-1}+\sqrt{n+1}\,\delta_{m,n+1})/\sqrt2$, and
exhibit the banded structure (truncation diagnostics); cross-check the entries against the
independent quadrature route, Integrate of $\psi_m\,x\,\psi_n$ over HermiteH states. The honest
truncation exhibit is the commutator corner defect
$[\hat x,\hat p]=i\,(\mathbb 1-(N{+}1)\,\vert N\rangle\langle N\vert)$ with corner entry $-iN$:
the trace of any finite commutator vanishes identically, so the corner must carry exactly $-N$, a
refuting check on the whole construction; sweep $N$ with Table to show the pattern. Then
Eigenvalues of the truncated $\tfrac12(\hat x^2+\hat p^2)$: every level below the cutoff is
exactly $n+\tfrac12$ while the top one is $N/2$ instead of $N+\tfrac12$. Close on the theorem the
defect encodes: the canonical commutation relation has no finite-dimensional representation, so a
corner of this exact size is forced, not a bug.

### 4.5 [BSc] How do I show the ground state is a Gaussian, and restore $\omega$ and $m$ to read off the oscillator length $\sqrt{\hbar/m\omega}$?

Keep the Gaussian ground state as the carrier (the question names it, so the named-trivial policy
holds): DSolve the dimensionful ODE with symbolic $E$ and symbolic $\hbar,m,\omega$ (C1: no decay
BCs, hang trap; quantization is the manual termination read-off at $n=0$), giving
$\psi_0(x)=(m\omega/\pi\hbar)^{1/4}e^{-m\omega x^2/2\hbar}$ with $E_0=\hbar\omega/2$ and
oscillator length $\ell=\sqrt{\hbar/m\omega}$; if the three-parameter symbolic solve stalls,
nondimensionalize by $x=\ell\xi$ first and restore units by substitution, since the rescaling is
the physics being asked. Close the dimensionful residual and unit norm with FullSimplify under
$\hbar,m,\omega>0$, anchor the rescaled equation on the exact DEigensystem domain form, and refute
by scaling: recompute $\ell$ from $\Delta x^2=\ell^2/2$ by Integrate and demand
$\sqrt{\hbar/m\omega}$; Limit shows $\ell\to\infty$ as $\omega\to0$, the bound state flattening
away. The discriminating contrast is two real laboratory scales via Quantity and UnitConvert: the
CO stretch ($\omega/2\pi c\approx2143\ \mathrm{cm}^{-1}$, reduced mass $\approx6.86\,\mathrm u$)
gives $\ell\approx5\ \mathrm{pm}$, while a Rb atom in a $2\pi\times100\ \mathrm{kHz}$ optical
tweezer gives $\ell\approx34\ \mathrm{nm}$, nearly four orders of magnitude apart. Close there:
one formula spans molecular vibration to cold atoms.

### 4.6 [MSc] How do I build a coherent state $|\alpha\rangle$ as the eigenstate of $a$, and equivalently by applying $e^{\alpha a^\dagger-\alpha^* a}$ to the ground state, and show it is a minimum-uncertainty state?

Build the coherent state at complex $\alpha=2e^{i\pi/3}$, so
$x_0=\sqrt2\,\operatorname{Re}\alpha=\sqrt2$ and $p_0=\sqrt2\,\operatorname{Im}\alpha=\sqrt6$
(quantum optics), twice. Eigenstate route: DSolve the first-order ODE
$\psi'=(\sqrt2\,\alpha-x)\psi$ and normalize with Integrate plus ComplexExpand. Displaced route:
split $D(\alpha)=e^{\alpha a^\dagger-\alpha^*a}$ by BCH as
$e^{-ix_0p_0/2}\,e^{ip_0\hat x}\,e^{-ix_0\hat p}$, legitimate because $[a,a^\dagger]$ is central
(the MSc ordering point), with $e^{-ix_0\hat p}$ acting as translation on $\psi_0$. The refuting
check is exact equality first: the ratio of the two constructions must FullSimplify to a constant
of modulus 1, since a modulus-only overlap would silently pass a missing $e^{-ix_0p_0/2}$ phase,
and both must satisfy $a\psi=\alpha\psi$; NIntegrate spot-checks the overlap. Then earn minimum
uncertainty with $\alpha$ kept symbolic complex: $\Delta x^2=\Delta p^2=\tfrac12$, so
$\Delta x\,\Delta p=\tfrac12$ exactly on the Heisenberg bound, with $\alpha\to0$ recovering the
vacuum. Close with the classical reading: $\langle\hat n\rangle=\vert\alpha\vert^2$ with relative
spread $\vert\alpha\vert^{-1}$, the packet turning classical at large $\vert\alpha\vert$.

### 4.7 [MSc] How do I build a squeezed state by applying $e^{\frac12(\xi^* a^2-\xi a^{\dagger2})}$ to the ground state and compute its unequal quadrature variances?

Apply $S(\xi)=e^{\frac12(\xi^*a^2-\xi a^{\dagger2})}$ to the vacuum with $\xi=re^{i\pi/3}$ and $r$
symbolic (quantum optics): the phase of $\xi$ kept, not $\theta=0$, so the noise ellipse is
genuinely rotated off the $x,p$ axes. Symbolic side: apply the Bogoliubov transform
$S^\dagger aS=a\cosh r-a^\dagger e^{i\theta}\sinh r$ as replacement rules, reduce the vacuum
moments, and close
$\operatorname{Var}(X_\phi)=\tfrac12\left(\cosh2r-\sinh2r\cos(2\phi-\theta)\right)$ with
FullSimplify: principal axes at $\phi=\theta/2$ carrying $e^{-2r}/2$ and $e^{+2r}/2$, product
exactly $\tfrac14$ there and larger on any other quadrature pair. Numeric side, on independent
machinery: truncated ladder matrices via SparseArray, MatrixExp of the truncated generator applied
to the vacuum vector; truncation and exponentiation do not commute, so sweep the truncation $N$
and compare against the exact even-$n$ amplitudes
$c_{2n}=(-e^{i\theta}\tanh r)^n\sqrt{(2n)!}/(2^n\,n!\,\sqrt{\cosh r})$, a check that refutes
truncation error and a sign-convention error at once; the $(\tanh r)^n$ tail slows convergence at
large $r$. Close with the limits: $r\to0$ recovers the vacuum disk, and large $r$ pushes one
rotated quadrature below shot noise while the conjugate diverges at fixed product.

### 4.8 [MSc] How do I compute the photon-number statistics of a number state and a coherent state (the Poisson distribution)?

Set the number state $n=4$ beside the coherent state $\alpha=2e^{i\pi/3}$, so
$\vert\alpha\vert^2=4$: identical mean energy $\langle\hat H\rangle=\tfrac92$, so only
the statistics discriminate (photon counting). Compute the amplitudes
$c_n=e^{-\vert\alpha\vert^2/2}\alpha^n/\sqrt{n!}$, probabilities with Abs under ComplexExpand,
and close the symbolic Sums $\sum_n P(n)=1$, $\bar n=\vert\alpha\vert^2$, and
$\Delta n^2=\vert\alpha\vert^2$: the check that refutes normalization or weight errors. Compare
against PDF of PoissonDistribution with Mean and Variance as the independent route, and watch the
partial-sum tail numerically. Mandel $Q=(\Delta n^2-\bar n)/\bar n$ lands at $0$ for the coherent
state and at $-1$ for the Fock state, the sub-Poissonian floor, and the large-$\bar n$ Gaussian
envelope of the Poisson is the classical limit. Close on the phase: the phase of $\alpha$ cancels
in every $P(n)$, so photon counting is blind to exactly what the quadratures of 4.6 and 4.7
resolve.
