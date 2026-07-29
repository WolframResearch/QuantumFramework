# What a Continuous-Measurement SDE Tells You Before You Simulate It

**Given only the stochastic equation of a continuously monitored quantum system, a great deal is exactly computable: closed-form moments, the ordinary differential equations a hierarchy closes into, exact rates and steady states, and the one deterministic law that actually tests whether a trajectory simulation got the measurement right. This is a guide to those relations, to the questions they answer, and to how the Wolfram Language noncommutative-algebra system derives and automates them from the equation alone.**

Mads Bahrami (last updated: July 28, 2026)

### Setting the Stage: How This Document Flows

You are handed a stochastic differential equation for a quantum system whose observable is being watched, and you want to simulate it. Before writing a single line of the integrator, you can extract exact statements the simulation must reproduce, and those statements are the sharpest tests it will ever face. This document is about getting them out of the equation and using them.

I strongly believe in a computation-first narrative for learning: in a sense, if I cannot compute it, I cannot claim to understand it. So nothing here is asserted that is not computed on the spot. The engine that does the computing is the noncommutative-algebra system added in Wolfram Language 15.0, and we hand it exactly one physical input, the commutator of the monitored variables, and let it carry the operator algebra.

This is a live notebook-style document. The cells depend on each other and are meant to be run top to bottom, like frames of one continuous movie rather than isolated snapshots. My suggestion is to read each output and its meaning first, and only then worry about the input that produced it. And remember that you are not locked into the code as given: change the Hamiltonian, change the monitored observable, and rerun. Questions and corrections are welcome at mohbahrami@gmail.com.

### Prerequisites and How to Read This

I assume quantum mechanics at the level of states, operators, and the master equation, and enough Wolfram Language to read a cell. No stochastic calculus is needed beyond one fact, introduced where it is used: the noise increment over a step $dt$ has variance $dt$. The pace steepens once, in Part II, where we hand the kernel the canonical commutator and let it do the operator ordering; slow down there. Throughout, $\hbar$, the mass, and the measurement strength $\lambda$ stay symbolic for as long as possible.

Let's start!

## Part I: From a Stochastic Equation to Exact Facts

A system observable $\hat A$ is monitored continuously. Conditioned on the detector record, the state obeys a stochastic differential equation; averaged over records, it obeys the Lindblad master equation. The two live side by side:

$$
\text{averaged:}\quad \dot\rho=-\tfrac{i}{\hbar}[\hat H,\rho]+\lambda\,\mathcal D[\hat A]\rho,
\qquad
\text{conditioned:}\quad d\rho=\Big(-\tfrac{i}{\hbar}[\hat H,\rho]+\lambda\,\mathcal D[\hat A]\rho\Big)dt+\sqrt{\lambda\eta}\,\mathcal H[\hat A]\rho\,dW,
$$

with $\mathcal D[\hat A]\rho=\hat A\rho\hat A^\dagger-\tfrac12\{\hat A^\dagger\hat A,\rho\}$ the dissipator, $\mathcal H[\hat A]\rho=\hat A\rho+\rho\hat A^\dagger-\langle\hat A+\hat A^\dagger\rangle\rho$ the measurement backaction, and $\eta\in(0,1]$ the detector efficiency. In other words, the first equation is what you see if you throw the record away, and the second is what one run of one experiment actually produces.

Here is a habit worth acquiring: before writing a simulation, extract every exact statement you can, because those become the benchmarks that catch the simulation being wrong. A stochastic integrator has many quiet failure modes, a drifting norm, a wrong sign on the noise, an Itô-versus-Stratonovich slip, a grid too coarse to hold the state, and none of them announces itself. An exact relation is a trap laid for each one. The rest of this document is about laying the traps.

## Part II: One Superoperator Runs Everything, the Adjoint Lindbladian

The averaged dynamics of *any* observable $\hat O$ is governed by a single superoperator acting on it:

$$
\frac{d}{dt}\langle\hat O\rangle=\langle\mathcal L^\dagger\hat O\rangle,\qquad
\mathcal L^\dagger\hat O=\tfrac{i}{\hbar}[\hat H,\hat O]+\lambda\big(\hat A^\dagger\hat O\hat A-\tfrac12\{\hat A^\dagger\hat A,\hat O\}\big).
$$

In other words, to know how any average evolves, you apply $\mathcal L^\dagger$ and take the expectation, and $\mathcal L^\dagger$ is nothing but commutators and anticommutators. That is exactly the object a noncommutative-algebra system exists to manipulate.

We will work the textbook case a monitored *position*, $\hat A=\hat x$, of a free particle, $\hat H=\hat p^2/2m$, because everything is checkable there. Version 15.0 lets us declare the two generators and the single relation between them and hand the ordering to the kernel. Define the algebra of $\hat x$ and $\hat p$ from the canonical commutator $[\hat x,\hat p]=i\hbar$, entered leading-monomial-positive, together with two thin wrappers for expanding and for commutators:

```wl
ClearAll[xp, nc, comm];
xp = NonCommutativeAlgebra[<|
    "Generators" -> {{xo, po}},
    "CommutationRelations" -> {po ** xo - xo ** po + I \[HBar]}|>];
nc[e_] := NonCommutativeExpand[e, xp];
comm[a_, b_] := nc[Commutator[a, b, xp]];
ListQ[xp["Generators"]]
```

The bare `True` confirms the object actually built; a relation that failed the algebra's consistency conditions would have returned unevaluated instead. Now the single most important commutator in the whole document, the one the heating rate rests on: the double commutator of $\hat x$ with $\hat p^2$.

```wl
comm[xo, comm[xo, po ** po]]
```

As one can see, $[\hat x,[\hat x,\hat p^2]]=-2\hbar^2$, a pure number with no operator left in it. That it collapses to a c-number is the whole reason the next result is exact. The measurement piece of $\mathcal L^\dagger$ acting on $\hat p^2$ is $-\frac{\lambda}{2}[\hat x,[\hat x,\hat p^2]]$, so compute it directly:

```wl
nc[-(λ/2) comm[xo, comm[xo, po ** po]]]
```

Therefore $\frac{d}{dt}\langle\hat p^2\rangle=\lambda\hbar^2$, exactly, with no approximation and nothing to tune. In words, watching a particle's position heats its momentum at a fixed rate forever, which is the price the uncertainty principle charges for position information. Note what this does *not* depend on: not the state, not the mass, not $\hbar$ except through the overall scale. This will be the single sharpest benchmark in Part V, precisely because it is a straight line with a known slope.

## Part III: When the Hierarchy Closes, You Get Closed Forms

One exact rate is good; a closed *system* is better, because it hands you every moment at once. The question that decides whether you get one is closure: does $\mathcal L^\dagger$ map the span of a finite set of operators back into that span? Equivalently, when you differentiate a moment, do only the moments you already have appear on the right-hand side, or does a new, higher one keep intruding?

Let us just ask the kernel. Encode $\mathcal L^\dagger$ for the free particle, apply it to the five low operators, and read what comes back. Define the generator and the candidate basis:

```wl
ClearAll[ldag, basis];
ldag[o_] := nc[(I/\[HBar]) comm[po ** po/(2 mass), o] - (λ/2) comm[xo, comm[xo, o]]];
basis = {xo, po, xo ** xo, po ** po, (xo ** po + po ** xo)/2};
ldag /@ basis
```

Have you noticed something interesting? Every entry that comes back is a linear combination of operators already in `basis`, together with constants: the $\hat p^2$ row is the pure constant $\lambda\hbar^2$, and the $\hat x^2$ row returns $2\hat x\hat p/m$ shifted by a constant, which is just the symmetrized $(\hat x\hat p+\hat p\hat x)/2$ already in the list (recall $\hat x\hat p=(\hat x\hat p+\hat p\hat x)/2+i\hbar/2$). In other words the list closes on itself: $\mathcal L^\dagger$ never produces a sixth operator. So writing $\vec m=(\langle\hat x\rangle,\langle\hat p\rangle,\langle\hat x^2\rangle,\langle\hat p^2\rangle,\langle(\hat x\hat p+\hat p\hat x)/2\rangle)$, the averages obey a closed linear system

$$
\dot{\vec m}=M\vec m+\vec b,
$$

and once a system is linear and closed, ordinary linear algebra delivers *everything*: the exact time course is $e^{Mt}$, the exact rates are the eigenvalues of $M$, and any steady state is $-M^{-1}\vec b$. Read the reading aloud: the $\langle\hat x\rangle$ row gives $\dot{\langle\hat x\rangle}=\langle\hat p\rangle/m$, the $\langle\hat p\rangle$ row gives $\dot{\langle\hat p\rangle}=0$, and the $\langle\hat p^2\rangle$ row is our heating law. Assemble $M$ and $\vec b$ by reading the coefficients off, and confirm the momentum row is a pure source with no linear part:

```wl
ClearAll[mMat, bVec];
mMat = {{0, 1/mass, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 2/mass},
        {0, 0, 0, 0, 0}, {0, 0, 0, 1/mass, 0}};
bVec = {0, 0, 0, λ \[HBar]^2, 0};
Eigenvalues[mMat]
```

As one can see, every eigenvalue is zero: $M$ is nilpotent. In physical terms the free-particle *ensemble* has no relaxation, it only spreads, and the closed system tells you exactly how. Solve it in closed form for the momentum variance and the position variance, starting from a state at rest with sharp initial moments:

```wl
sol = DSolve[{mx'[t] == mp[t]/mass, mp'[t] == 0, mxx'[t] == 2 s[t]/mass,
    mpp'[t] == λ \[HBar]^2, s'[t] == mpp[t]/mass,
    mx[0] == 0, mp[0] == 0, mxx[0] == x0, mpp[0] == p0, s[0] == 0},
   {mx, mp, mxx, mpp, s}, t];
{mpp[t], mxx[t]} /. First[sol]
```

Therefore $\langle\hat p^2\rangle(t)=p_0+\lambda\hbar^2 t$ grows linearly and $\langle\hat x^2\rangle(t)$ grows like $t^3$, both in exact closed form with no simulation. Recall what bought this: not a clever integrator, just the fact that $\mathcal L^\dagger$ closed on five operators. That closure is the entire mechanism, and it is a property of the equation you can test before simulating anything.

A caution worth stating plainly, because it is where the idea has teeth and limits. Closure here means closure on a finite set of *monomials* in the generators; a system that closed only on some non-monomial set (a projector, a symmetry-adapted combination) would look non-closing to this search. And closure is decided by this computation, not assumed: handing the algebra a quadratic $\hat H$ and a linear $\hat A$ is what makes the low moments close, and a nonlinear $\hat H$ generically breaks it. We return to the broken case at the end.

## Part IV: The Relation That Actually Tests the Measurement

Everything in Part III came from the *averaged* equation, so it holds with the detector unplugged: a purely dissipative particle heats and spreads in exactly the same way. In other words, those benchmarks certify the drift, and a simulator can reproduce every one of them while getting the *measurement* completely wrong. The relation that tests the measurement lives in the conditioned equation, and it is the sharpest of all.

Here is the fact, and it is genuinely surprising. Watch the conditional variances of a single monitored trajectory, the width $V_x$, the momentum spread $V_p$, and the covariance $C$. Each individual run has a randomly wandering center, so you would expect its width to wander too. It does not. The same operator algebra plus the Itô rule (the noise increment squares to $dt$) shows that the random parts of the variance equations cancel identically, leaving *deterministic* ordinary differential equations,

$$
\dot V_x=\frac{2C}{m}-4\lambda\eta\,V_x^2,\qquad
\dot C=\frac{V_p}{m}-4\lambda\eta\,C\,V_x,\qquad
\dot V_p=\lambda\hbar^2-4\lambda\eta\,C^2 .
$$

In other words, every experimental run produces a differently placed wave packet, but all of them have the identical width at the identical time. This is a Riccati system, it is the continuous-measurement instance of the Kalman-Bucy filter, and it is where the efficiency $\eta$ finally appears: the heating term $\lambda\hbar^2$ is the backaction and is $\eta$-independent, while every localization term carries $\eta$, the information you actually collect. Set $\eta$ to zero and the localization vanishes and the width spreads like the ensemble; that is the sanity check that ties this back to Part III.

Find the steady state, keeping everything symbolic, and select the physical (positive) branch:

```wl
ClearAll[fx];
fx = First@Select[
   Solve[{2 Cc/mass - 4 λ η Vx^2 == 0, Vp/mass - 4 λ η Cc Vx == 0,
      λ \[HBar]^2 - 4 λ η Cc^2 == 0}, {Vx, Cc, Vp}],
   TrueQ@Simplify[(Vx /. #) > 0 && (Vp /. #) > 0,
      \[HBar] > 0 && mass > 0 && λ > 0 && η > 0] &];
Simplify[{Sqrt[Vx], Vx Vp - Cc^2} /. fx, \[HBar] > 0 && mass > 0 && λ > 0 && η > 0]
```

As one can see, a watched free particle relaxes to a packet of fixed width $\sigma_x=\frac{1}{\sqrt2}\big(\tfrac{\hbar}{\lambda m}\big)^{1/4}\eta^{-3/8}$, and the second entry says $V_xV_p-C^2=\hbar^2/(4\eta)$. At perfect efficiency $\eta=1$ this is $\hbar^2/4$, a *minimum-uncertainty* packet; below it the product rises, because an inefficient detector leaves the conditional state mixed, and its purity, $\hbar/\big(2\sqrt{V_xV_p-C^2}\big)$, comes out to exactly $\sqrt\eta$. Three things to read from this. The width grows as the efficiency falls, like $\eta^{-3/8}$, so a detector at half efficiency gives a measurably wider steady packet, and that dependence is a fingerprint of the conditioning term. The steady purity $\sqrt\eta$ is a second, cleaner fingerprint, a single number that a drift-only simulation cannot produce. And $\sigma_x$ is precisely the length your spatial grid must resolve, the guideline a simulation needs before it starts.

Now the benchmark this is really for. Build one conditional quantity that depends on the sign, the prefactor, the efficiency, and the Itô-versus-Stratonovich reading of the noise term $\sqrt{\lambda\eta}\,\mathcal H[\hat A]\rho\,dW$, and you have a test the drift benchmarks cannot fake. Two good ones: the steady width above, which pins $\lambda\eta$; and the short-time growth of the conditional purity, which is linear in $\lambda\eta\,dt$ and *changes sign* under the wrong convention. A simulator with the drift right and the conditioning wrong reproduces the heating slope $\lambda\hbar^2$ and fails these. That asymmetry is the whole reason to derive them.

## Part V: The Questions This Answers, and Who Should Care

Let us step back from the free particle and ask what this approach buys anyone with a continuous-measurement code. Each question below is one a simulation cannot answer about itself and an exact relation can.

- **Is my conditioning term correct?** The steady width $\sigma_x\propto\lambda^{-1/4}\eta^{-3/8}$, the steady purity $\sqrt\eta$, and the sign of the short-time purity growth together pin the $dW$ coefficient, its efficiency, and its stochastic-calculus convention. This is the question the ensemble averages cannot touch.
- **Does my code conserve what the physics conserves?** The exact heating rate $\lambda\hbar^2$ is a straight line of known slope; a norm-drifting or over-damped integrator bends it.
- **Can I avoid Monte Carlo entirely?** If $\mathcal L^\dagger$ closes on a finite set, the moments obey a closed ODE system and you never sample at all. Testing closure is itself an exact, pre-simulation question.
- **What grid or cutoff do I need?** The steady width sets the resolution from below; the heating $\sqrt{\lambda\hbar^2 t}$ sets the momentum range from above, and those two requirements pull apart as $\lambda$ grows.
- **Does my ensemble reproduce the master equation?** Averaging many trajectories must return the Lindblad result; the closed moment system gives the exact target to hit.

In short, someone should care because these relations turn "my simulation looks plausible" into "my simulation reproduces a number I derived without it." That is the difference between a picture and a result.

## Part VI: Doing It From the Equation Alone, with Noncommutative Algebra

Notice that nothing in Parts II to IV needed anything but the commutator $[\hat x,\hat p]=i\hbar$ and the two operators $\hat H$ and $\hat A$. In other words, the entire pipeline is: name the algebra the monitored variable lives in, build $\mathcal L^\dagger$ from $\hat H$ and $\hat A$, normal-order with `NonCommutativeExpand`, and read off closure. The algebra is chosen by the system, and the choice is not a matter of taste.

For a bosonic mode it is the canonical commutator we already used. For a spin it is $\mathfrak{su}(2)$, and the same declaration works: hand the kernel the three angular-momentum relations and it builds the enveloping algebra directly. Verify that a spin monitored through $\hat J_z$ lives in an algebra the kernel accepts:

```wl
ClearAll[su2, ncS, commS];
su2 = NonCommutativeAlgebra[<|
    "Generators" -> {{jx, jy, jz}},
    "CommutationRelations" -> {jy ** jx - jx ** jy + I jz,
       jz ** jy - jy ** jz + I jx, jz ** jx - jx ** jz - I jy}|>];
ncS[e_] := NonCommutativeExpand[e, su2];
commS[a_, b_] := ncS[Commutator[a, b, su2]];
{ListQ[su2["Generators"]], commS[jx, jy]}
```

The pair reads `{True, I jz}`: the algebra built, and it reduces $[\hat J_x,\hat J_y]$ to $i\hat J_z$ as it must. So the machinery that gave the position heating law will, unchanged, give the dephasing rate of a monitored spin; only the generators and the one relation differ. The lesson that transfers is the question, not the answer: *which algebra does my monitored observable generate, and does the moment hierarchy close in it?*

One honest note on the reading. Asking the kernel for the Gröbner basis of the three spin relations returns the three relations themselves, which tells you they are already a consistent set (a reduced basis in the free algebra), i.e. the algebra builds with no completion needed. It is a statement about the encoding being well-posed, not yet about which moments close; closure is the separate computation of Part III, applied in this algebra.

## Part VII: Automating It

The pipeline is mechanical, so let us mechanize it. Write one function that takes the Hamiltonian, the monitored observable, an operator basis, and the algebra, and applies $\mathcal L^\dagger$ to every basis element, so that closure can be read straight off the result. Define it and run it on the free particle:

```wl
ClearAll[applyGenerator];
applyGenerator[ham_, aObs_, bas_, alg_] :=
  NonCommutativeExpand[
     (I/\[HBar]) Commutator[ham, #, alg]
     - (λ/2) Commutator[aObs, Commutator[aObs, #, alg], alg], alg] & /@ bas;
applyGenerator[po ** po/(2 mass), xo, basis, xp]
```

This is the same list we read by hand in Part III: every entry sits in the span of `basis`, so the system closes and the moment ODEs are the ones we already solved. Now feed it a monitored anharmonic oscillator, $\hat H=\hat p^2/2m+g\,\hat x^4$, and read the result:

```wl
applyGenerator[po ** po/(2 mass) + g xo ** xo ** xo ** xo, xo, basis, xp]
```

As one can see, monomials of higher degree in $\hat x$ now appear, operators that are not in `basis`. In other words the quartic term pushes each moment onto a higher one, the list does not close, and no finite $M$ exists. That is the honest verdict the tool gives before any simulation: whether a closed benchmark exists at all. Deciding closure by eye from a five-operator list is easy; a fuller tool would compare the operator monomials of the two lists automatically, with one care, that the algebra's normal form writes $\hat p\hat x$ as $\hat x\hat p-i\hbar$, so the comparison must be between normal-ordered monomials rather than raw products.

I want to be equally honest about what this automation is and is not. It is not a new method: symbolic engines already derive these operator equations of motion from a Hamiltonian and its jump operators, QNET over the input-output formalism and QuantumCumulants for mean-field hierarchies among them, and the closure question is the well-studied finite-dimensional-filter problem in disguise. What the noncommutative-algebra route offers is a second, independent symbolic path, run in a different system, that you can set against those and against a concrete matrix truncation. In short, its worth is validation, and a benchmark is worth having even when the method behind it is not new.

## Where This Leaves Us (and What Comes Next)

We now have a computation-first toolkit for reading a continuous-measurement SDE before simulating it: one superoperator $\mathcal L^\dagger$ that governs every average, a heating law $\frac{d}{dt}\langle\hat p^2\rangle=\lambda\hbar^2$ that follows from a single double commutator, a closure test that decides whether the moments obey a closed ODE system and hence a closed-form solution, a deterministic conditional-covariance Riccati whose steady width $\sigma_x\propto(\lambda\eta)^{-1/4}$ is the one relation that actually tests the measurement, and a small function that runs the whole pipeline from the equation alone and reports non-closure honestly.

The natural next step, and the one I would take first, is to build the conditioning benchmark for a concrete system, a monitored qubit or a single quadrature, and set it against a trajectory simulator: draw many records, average, and watch the conditional width sit on the deterministic Riccati value to many digits while the ensemble decoheres around it. That single test exercises the sign, the efficiency, and the stochastic-calculus convention all at once, which is exactly the part of the code that can be wrong. Everything here is one modification away from it: change `ham`, change the monitored observable, and rerun.
