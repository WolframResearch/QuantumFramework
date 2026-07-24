# Part 17 Plan: Continuous-variable quantum optics and information

Seven questions, all MSc. Class census per `Route-Table.md`: C0: 17.1, 17.2, 17.3, 17.4, 17.5,
17.7 (no differential equation; WL machinery named per entry, C0 kernel facts binding); C5: 17.6
(truncated-basis ODE-IVP: `MatrixExp` action form primary, route agreement provably blind to
truncation, $N$-sweep mandatory).

Everything in this part lives on one or two bosonic modes with $[a,a^\dagger]=1$ and the
quadratures $X=(a+a^\dagger)/\sqrt2$, $P=-i(a-a^\dagger)/\sqrt2$, obeying $[X,P]=i$ with vacuum
variances $\tfrac12$ each; the Gaussian unitaries are the displacement
$D(\alpha)=e^{\alpha a^\dagger-\alpha^*a}$ and the squeeze
$S(\xi)=e^{(\xi^*a^2-\xi a^{\dagger2})/2}$, with their two-mode kin, the beam splitter
$e^{\theta(a^\dagger b-ab^\dagger)}$ and the two-mode squeezer $e^{r(ab-a^\dagger b^\dagger)}$
(each phase and sign convention: verify at authoring); phase space carries the Wigner function
$W(x,p)$ and the Husimi $Q(\alpha)=\langle\alpha\vert\rho\vert\alpha\rangle/\pi$. The working
representation throughout is truncated Fock: `SparseArray` ladder matrices, `KroneckerProduct` for
two modes, `MatrixExp` for the unitaries, with the truncation edge exhibited and checked rather
than hidden. The C0 kernel facts bind (in code `E` is Euler's number, so energies and fields need
other names; `Simplify` and `Integrate` need explicit reality and positivity assumptions such as
$r>0$ to close Gaussian moments).

### 17.1 [MSc] How do I define the quadrature operators and the optical phase space of a single mode?

Build the ladder matrix $a$ as a `SparseArray` with superdiagonal $\sqrt1,\dots,\sqrt N$ at $N=12$,
form $X=(a+a^\dagger)/\sqrt2$ and $P=-i(a-a^\dagger)/\sqrt2$ with `ConjugateTranspose`, and compute
the commutator honestly: it equals $i\,\mathbb 1$ everywhere except the corner,
$[X,P]=i(\mathbb 1-(N{+}1)\vert N\rangle\langle N\vert)$, the truncation defect Part 4 establishes
for $\hat x,\hat p$ (4.4), acknowledged here as the price of a finite matrix while every interior
element is exact. Vacuum variances $\langle X^2\rangle=\langle P^2\rangle=\tfrac12$ by `Dot` on the
vacuum unit vector; the refuting cross-check runs on machinery the matrices never touch:
`Integrate` of $x^2\vert\psi_0\vert^2$ for the position-space vacuum
$\psi_0(x)=\pi^{-1/4}e^{-x^2/2}$ must return the same $\tfrac12$, and a wrong $\sqrt2$ convention
in either construction breaks the agreement. Close operationally: the product saturates at
$\Delta X\,\Delta P=\tfrac12$, the vacuum's footprint in the $(x,p)$ plane is the uncertainty disk
of radius $\sqrt{1/2}$, and the variance $\tfrac12$ is the shot-noise floor a balanced detector
reads on an empty port.

### 17.2 [MSc] How do I revisit the displacement and squeeze operators of Part 4 as phase-space transformations and read off their action on the quadratures?

Carry Part 4's parameters $\alpha=2e^{i\pi/3}$ and $\xi=re^{i\pi/3}$ so both phases are live, and
conjugate by BCH: for $D(\alpha)$ the commutator series terminates after one step (the commutator
of $X$ with $\alpha a^\dagger-\alpha^*a$ is a c-number), giving
$D^\dagger(\alpha)XD(\alpha)=X+\sqrt2\operatorname{Re}\alpha$ and $P$ shifted by
$\sqrt2\operatorname{Im}\alpha$; for $S(\xi)$ the nested commutators resum to hyperbolic mixing,
the quadratures along the squeeze axes scaling as $e^{\mp r}$ and mixing $X$ with $P$ whenever
$\phi\neq0$ (squeeze sign convention: verify at authoring). Verify twice, on pictures that must
agree: conjugate the truncated matrices of the 17.1 construction by `MatrixExp` of the truncated
generators and compare interior matrix elements to the closed forms (the corner rows are the known
truncation edge, not evidence); and transform the vacuum Wigner Gaussian, letting `Integrate` under
reality assumptions show the ellipse center lands at
$(\sqrt2\operatorname{Re}\alpha,\sqrt2\operatorname{Im}\alpha)$ and the principal variances scale
by $e^{\mp2r}$ at fixed area. Close with the shared picture: Gaussian unitaries act affinely on
phase space, displacement translates the vacuum disk, squeezing deforms it area-preservingly.

### 17.3 [MSc] How do I model a beam splitter on two modes and exhibit Hong-Ou-Mandel interference?

The beam splitter is $B(\theta)=e^{\theta(a^\dagger b-ab^\dagger)}$ (phase convention: verify at
authoring), realized by `MatrixExp` on the `KroneckerProduct` two-mode Fock space truncated at two
photons per mode; the generator conserves total photon number, so the two-photon sector
$\{\vert2,0\rangle,\vert1,1\rangle,\vert0,2\rangle\}$ closes exactly and the truncation is not an
approximation here. At the balanced point $\theta=\pi/4$ act on $\vert1,1\rangle$: the output is
$(\vert2,0\rangle+\vert0,2\rangle)/\sqrt2$ up to a convention-dependent relative sign, with
coincidence amplitude exactly zero. The refuting cross-check derives the null on different
machinery: transform the creation operators symbolically,
$a^\dagger\to(a^\dagger+b^\dagger)/\sqrt2$, $b^\dagger\to(b^\dagger-a^\dagger)/\sqrt2$, expand
$a^\dagger b^\dagger\vert0,0\rangle$, and watch the $\vert1,1\rangle$ terms cancel bosonically,
with unitarity of $B$ and output probabilities summing to 1 as structural gates. Then the pinned
distinguishability sweep, the point of the example: give the second photon overlap $s\in[0,1]$
with the first by splitting its mode into parallel and orthogonal components, derive the closed
form $P_{\mathrm{coinc}}(s)=\tfrac12(1-s^2)$, and validate the endpoints against a four-mode (two
spatial times two internal) matrix computation; the dip visibility $s^2$ interpolates from the
distinguishable plateau $\tfrac12$ down to the null. Close in the laboratory: the HOM dip is the
standard meter of single-photon indistinguishability, and classical fields cap the dip at half the
plateau, so a deeper dip certifies quantum interference.

### 17.4 [MSc] How do I build the two-mode squeezed vacuum and exhibit its EPR (position-momentum) correlations?

Apply `MatrixExp` in its action form to the two-mode vacuum with the generator
$r(ab-a^\dagger b^\dagger)$ (sign convention: verify at authoring) on a `KroneckerProduct` space
at $N=24$ per mode, keeping $r$ symbolic in the closed forms and $r=1$ for numbers. The state must
be the Schmidt-diagonal $\sqrt{1-\lambda^2}\sum_n\lambda^n\vert n,n\rangle$ with
$\lambda=\tanh r$: `ArrayReshape` the output into its coefficient matrix $c$ and let
`SingularValueList` return a geometric sequence of ratio $\tanh r$, while the reduced state
$cc^\dagger$ is thermal with $\bar n=\sinh^2 r$, the first quantification of the entanglement; the
second is the EPR pair of variances $\Delta(X_1-X_2)^2=e^{-2r}$ and $\Delta(P_1+P_2)^2=e^{-2r}$,
computed with per-mode quadrature matrices on the same state. Truncation is checked by an
$N$-sweep with its tail predicted in advance: the missing weight is
$\lambda^{2(N+1)}\approx10^{-6}$ at $r=1$, $N=24$, and pushing $r$ up at fixed $N$ must visibly
break the variance identity, the honest failure edge. Close on the limit: as $r\to\infty$ the pair
approaches the ideal EPR state, the correlations sharpening without bound while each mode alone
heats toward a featureless thermal state.

### 17.5 [MSc] How do I detect continuous-variable entanglement with the Duan separability criterion?

Evaluate the Duan combination $\Delta(X_1-X_2)^2+\Delta(P_1+P_2)^2$, bounded below by 2 for every
separable state in the $[X,P]=i$ normalization (bound normalization: verify at authoring). On the
two-mode squeezed vacuum, rebuilt inside the entry from the $r(ab-a^\dagger b^\dagger)$ generator
so the answer stays self-contained, the Bogoliubov action sends both combinations to $e^{-r}$
times themselves, so the sum is $2e^{-2r}$, below the bound for every $r>0$; on the coherent
product $D(\alpha)\otimes D(\beta)\vert0,0\rangle$ with $\alpha,\beta$ kept symbolic,
displacement moves means and not variances, so the sum is exactly 2 for all $\alpha,\beta$: one
violating state and one boundary state, both in closed form. The truncated-matrix evaluation must
reproduce both numbers, the $r\to0$ limit of the squeezed sum must land exactly on the
coherent-product value 2 (gluing the two computations together), and a sign slip in the squeeze
generator would squeeze the wrong combinations and push the sum above 2, so the check
discriminates conventions as well as states. Close on tightness: a product state sitting exactly
on the boundary shows the constant 2 cannot be improved, and this same combination is the witness
measured on twin beams in the laboratory.

### 17.6 [MSc] How do I solve the Jaynes-Cummings model (one mode plus a two-level atom, truncated Fock) and exhibit vacuum Rabi oscillations?

C5 per `Route-Table.md`. On resonance the initial state $\vert e,0\rangle$ closes into the
two-dimensional block $\{\vert e,0\rangle,\vert g,1\rangle\}$, and `DSolveValue` on the
two-amplitude system returns the exact pair $\{\cos gt,\,-i\sin gt\}$ (the probed C5 gate), so
$P_e(t)=\cos^2 gt$: vacuum Rabi oscillation at frequency $2g$ out of a field containing zero
photons. Generalize with detuning $\Delta$ through the same `DSolveValue`: oscillation at
$\sqrt{4g^2+\Delta^2}$ with transfer amplitude $4g^2/(4g^2+\Delta^2)$, swept over $\Delta$. The
certified numeric route is the C5 primary: assemble the full truncated JC Hamiltonian with
`SparseArray` and `KroneckerProduct` (mode times atom) and evolve with the `MatrixExp` action
form; check truncation with the mandatory $N$-sweep, because per the C5 verdict route agreement is
provably blind to truncation (two routes agreed to $10^{-7}$ while both were wrong by $10^{-3}$).
Here the sweep also teaches: the rotating-wave Hamiltonian conserves excitation number, so from
$\vert e,0\rangle$ the $N$-sweep must come out exactly flat, and any $N$-dependence refutes the
Hamiltonian assembly itself. Close in the cavity: the $2g$ oscillation is the single-atom vacuum
Rabi frequency of cavity QED, and the dispersive limit $\Delta\gg g$ shrinks the transfer to
$4g^2/\Delta^2$, the regime qubit readout lives in.

### 17.7 [MSc] How do I model balanced homodyne and heterodyne detection, measuring a quadrature and reconstructing the Husimi-$Q$ distribution operationally?

Take the squeezed vacuum at $r=1$ and derive by `Integrate` (assumptions $r>0$, $\theta$ real) the
closed-form marginal of its Wigner Gaussian along the rotated quadrature
$X_\theta=X\cos\theta+P\sin\theta$: a Gaussian of variance
$\sigma_\theta^2=\tfrac12(e^{-2r}\cos^2\theta+e^{2r}\sin^2\theta)$. Homodyne is sampling this
marginal at LO phase $\theta$: draw $10^5$ samples per phase with `RandomVariate` on
`NormalDistribution` of standard deviation $\sigma_\theta$, `Histogram` them against the closed
density, and trace the sample `Variance` around the $\sigma_\theta^2$ curve, sub-vacuum at
$\theta=0$ and anti-squeezed at $\theta=\pi/2$. Heterodyne samples the Husimi
$Q(\alpha)=\vert\langle\alpha\vert\psi\rangle\vert^2/\pi$ instead: realize it operationally by
adding an independent vacuum draw of variance $\tfrac12$ per quadrature to a draw from the state's
(here positive) Wigner Gaussian, exhibiting the convolution $Q=W\ast W_{\mathrm{vac}}$ with
covariances $\sigma_\theta^2+\tfrac12$. The refuting check runs twice: sampled variances must
match the closed forms within the $\sqrt{2/n}$ relative statistical tolerance, and the marginal's
variance must equal the operator variance $\langle X_\theta^2\rangle$ computed on truncated
quadrature matrices built as in 17.1, an independent route through the operators. Close on the
contrast the question asks for: the operational reconstruction (histograms of detector clicks)
converges to the formal objects (the Wigner marginal, the $Q$ function), heterodyne pays exactly
half a vacuum unit per quadrature, the 3 dB penalty for reading both knobs at once, and the
$\theta$-resolved homodyne variances are the raw data of quantum state tomography.
