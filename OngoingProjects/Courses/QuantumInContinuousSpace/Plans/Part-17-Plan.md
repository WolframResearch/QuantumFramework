# Part 17 Plan: Continuous-variable quantum optics and information

Seven questions, all MSc. Class census per `Route-Table.md`: C0: 17.1, 17.2, 17.3, 17.4, 17.5,
17.7 (no differential equation; WL machinery named per entry, C0 kernel facts binding); C5: 17.6
(truncated-basis ODE-IVP: `MatrixExp` action form primary, route agreement provably blind to
truncation, $N$-sweep mandatory).

This part lives on one or two bosonic modes, and in 17.6 on one mode times a two-level atom, with
$[a,a^\dagger]=1$ and the quadratures $X=(a+a^\dagger)/\sqrt2$, $P=-i(a-a^\dagger)/\sqrt2$ obeying
$[X,P]=i$ with vacuum variances $\tfrac12$ each, the rotated one being
$X_\theta=(ae^{-i\theta}+a^\dagger e^{i\theta})/\sqrt2$. The Gaussian unitaries are the
displacement $D(\alpha)=e^{\alpha a^\dagger-\alpha^*a}$; the squeeze
$S(\xi)=e^{(\xi^*a^2-\xi a^{\dagger2})/2}$, the convention pinned by 4.7, under which $S(r)$ with
$r>0$ squeezes $X$ to $\langle X^2\rangle=e^{-2r}/2$ and stretches $P$; the beam splitter
$B(\theta)=e^{\theta(a^\dagger b-ab^\dagger)}$; and the two-mode squeezer
$S_2(r)=e^{r(a^\dagger b^\dagger-ab)}$, written creation-minus-annihilation, the opposite overall
sign from the single-mode pattern, which is exactly what makes $X_1-X_2$ and $P_1+P_2$ the
squeezed pair (each scaling by $e^{-r}$) rather than the anti-squeezed one. Where the
beam-splitter phase sits is a genuine convention (verify at authoring): the coincidence null does
not depend on it, the relative sign between $\vert2,0\rangle$ and $\vert0,2\rangle$ does. The atom
enters through the Jaynes-Cummings Hamiltonian
$H=\omega_c a^\dagger a+\tfrac{\omega_a}{2}\sigma_z+g(a\sigma_++a^\dagger\sigma_-)$ with detuning
$\Delta=\omega_a-\omega_c$; in the frame rotating at $\omega_c$ for both mode and atom this becomes
$\tfrac{\Delta}{2}\sigma_z+g(a\sigma_++a^\dagger\sigma_-)$, where the exact two-amplitude solutions
below are read. Phase space carries the Wigner function $W(x,p)$ and the Husimi
$Q(\alpha)=\vert\langle\alpha\vert\psi\rangle\vert^2/\pi$. The working representation throughout is
truncated Fock: `SparseArray` ladder matrices, `KroneckerProduct` for two modes or mode times atom,
`MatrixExp` for the unitaries, with the truncation edge exhibited and swept rather than hidden. The
C0 kernel facts bind (in code `E` is Euler's number, so energies and fields need other names;
`Simplify` and `Integrate` need explicit reality and positivity assumptions such as $r>0$ to close
Gaussian moments).

### 17.1 [MSc] How do I define the quadrature operators and the optical phase space of a single mode?

Build the ladder matrix $a$ as a `SparseArray` with superdiagonal $\sqrt1,\dots,\sqrt N$ at $N=12$,
form $X=(a+a^\dagger)/\sqrt2$ and $P=-i(a-a^\dagger)/\sqrt2$ with `ConjugateTranspose`, and read the
optical phase space off them: $[X,P]=i$ follows from $[a,a^\dagger]=1$ alone, the vacuum unit vector
gives $\langle X\rangle=\langle P\rangle=0$ and $\langle X^2\rangle=\langle P^2\rangle=\tfrac12$ by
`Dot`, so the vacuum saturates $\Delta X\,\Delta P=\tfrac12$ and occupies a disk of radius
$\sqrt{1/2}$ in the $(x,p)$ plane, and a `Table` of $\langle X_\theta^2\rangle$ over $\theta$ comes
out flat at $\tfrac12$: the vacuum is isotropic, which is what makes $\tfrac12$ one shot-noise
floor rather than a phase-dependent one. On truncated matrices the commutator carries the corner
defect $[X,P]=i(\mathbb 1-(N{+}1)\vert N\rangle\langle N\vert)$, exact in the interior, the same
edge 4.4 exhibits for $\hat x,\hat p$. The refuting cross-check runs on machinery the matrices never
touch: `Integrate` of $x^2\vert\psi_0\vert^2$ for $\psi_0(x)=\pi^{-1/4}e^{-x^2/2}$, and of
$p^2\vert\tilde\psi_0\vert^2$ after `FourierTransform`, must both return $\tfrac12$, where a
misplaced $\sqrt2$ in either quadrature definition would land on $1$ or $\tfrac14$ instead. Close
operationally: $\tfrac12$ is the noise power a balanced detector reads with nothing in the signal
port, and every Gaussian state in this part is this one disk translated, rotated, or deformed at
fixed area.

### 17.2 [MSc] How do I revisit the displacement and squeeze operators of Part 4 as phase-space transformations and read off their action on the quadratures?

Carry Part 4's parameters $\alpha=2e^{i\pi/3}$ and $\xi=re^{i\pi/3}$ so both phases stay live, and
conjugate by BCH: for $D(\alpha)$ the commutator series terminates after one step, since the
commutator of $X$ with $\alpha a^\dagger-\alpha^*a$ is a c-number, giving
$D^\dagger(\alpha)XD(\alpha)=X+\sqrt2\operatorname{Re}\alpha$ and $P$ shifted by
$\sqrt2\operatorname{Im}\alpha$; for $S(\xi)$ the nested commutators resum to hyperbolic mixing,
the quadratures along the squeeze axes scaling as $e^{\mp r}$ and mixing $X$ with $P$ whenever
$\phi\neq0$, with the overall sign fixed by 4.7's convention rather than left open ($S(r)$ shrinks
$X$). Verify twice, on pictures that must agree: conjugate the truncated matrices of the 17.1
construction by `MatrixExp` of the truncated generators at $N=40$ and compare the leading
$20\times20$ block against the closed forms, demanding `Max@Abs` of the difference below $10^{-8}$
and raising $N$ until it is, with the deviation climbing toward the corner as the diagnostic rather
than a failure; and transform the vacuum Wigner Gaussian, letting `Integrate` under reality
assumptions put the ellipse center at
$(\sqrt2\operatorname{Re}\alpha,\sqrt2\operatorname{Im}\alpha)$ and scale the principal variances by
$e^{\mp2r}$ at fixed area. Close with the shared picture: Gaussian unitaries act affinely on phase
space, displacement translating the vacuum disk, squeezing deforming it at constant area, which is
why neither can produce Wigner negativity.

### 17.3 [MSc] How do I model a beam splitter on two modes and exhibit Hong-Ou-Mandel interference?

The beam splitter is $B(\theta)=e^{\theta(a^\dagger b-ab^\dagger)}$, realized by `MatrixExp` on the
`KroneckerProduct` two-mode Fock space truncated at two photons per mode; the generator conserves
total photon number, so the two-photon sector
$\{\vert2,0\rangle,\vert1,1\rangle,\vert0,2\rangle\}$ closes exactly and truncation is not an
approximation here. At the balanced point $\theta=\pi/4$ act on $\vert1,1\rangle$: the substitution
this convention produces, $a^\dagger\to(a^\dagger+b^\dagger)/\sqrt2$ and
$b^\dagger\to(b^\dagger-a^\dagger)/\sqrt2$ in $a^\dagger b^\dagger\vert0,0\rangle$, gives
$(\vert0,2\rangle-\vert2,0\rangle)/\sqrt2$ with the $\vert1,1\rangle$ terms cancelling bosonically,
so the coincidence amplitude is exactly zero; the matrix route returns the same state up to an
overall sign, a global phase carrying no physics, so the comparison is on the coincidence weight and
the relative sign, never the global one. Unitarity of $B$ and output probabilities summing to 1 are
the structural gates. Then the pinned distinguishability sweep, the point of the example: give the
second photon the internal mode $s\,e_1+\sqrt{1-s^2}\,e_2$ on a four-mode (two spatial times two
internal) space, derive the closed form $P_{\mathrm{coinc}}(s)=\tfrac12(1-s^2)$, and evaluate the
four-mode matrix computation at intermediate $s$, notably $s=1/\sqrt2$, where the rival forms
$(1-s)/2$, $(1-s^3)/2$ and $(1-s^2)^2/2$ that agree with it at both endpoints separate from it by
tens of percent: the endpoints alone certify nothing. Close in the laboratory: the dip visibility
$s^2$ is the standard meter of single-photon indistinguishability, and a classical field caps the
dip at half the distinguishable plateau, so a deeper dip certifies quantum interference.

### 17.4 [MSc] How do I build the two-mode squeezed vacuum and exhibit its EPR (position-momentum) correlations?

Apply `MatrixExp` in its action form to the two-mode vacuum with the generator
$r(a^\dagger b^\dagger-ab)$ on a `KroneckerProduct` space at $N=24$ per mode, keeping $r$ symbolic
in the closed forms and $r=1$ for numbers. The state must be the Schmidt-diagonal
$\sqrt{1-\lambda^2}\sum_n\lambda^n\vert n,n\rangle$ with $\lambda=\tanh r$: `ArrayReshape` the
output into its coefficient matrix $c$ and let `SingularValueList` return a geometric sequence of
ratio $\tanh r$, while the reduced state $cc^\dagger$ is thermal with $\bar n=\sinh^2r$, the first
quantification of the entanglement; the second is the EPR pair
$\Delta(X_1-X_2)^2=\Delta(P_1+P_2)^2=e^{-2r}$ against the vacuum value 1, computed with per-mode
quadrature matrices on the same state. That variance identity is also the sign discriminator, and
the only check here that is one: `SingularValueList` cannot see a sign (the flipped generator gives
$(-\lambda)^n$, the same singular values) and neither can the truncation tail, whereas the flipped
generator sends $\Delta(X_1-X_2)^2$ to $e^{+2r}$ and squeezes $X_1+X_2$ instead. Truncation is
checked by an $N$-sweep with its tail predicted in advance, the missing weight
$\lambda^{2(N+1)}\approx10^{-6}$ at $r=1$, $N=24$, and pushing $r$ up at fixed $N$ must visibly
break the variance identity, the honest failure edge. Close on the limit: as $r\to\infty$ the pair
approaches the ideal EPR state, the correlations sharpening without bound while each mode alone
heats toward a featureless thermal state.

### 17.5 [MSc] How do I detect continuous-variable entanglement with the Duan separability criterion?

Evaluate the Duan combination $\Delta(X_1-X_2)^2+\Delta(P_1+P_2)^2$, bounded below by 2 for every
separable state in the $[X,P]=i$ normalization (bound normalization: verify at authoring). On the
two-mode squeezed vacuum, rebuilt inside the entry from the $r(a^\dagger b^\dagger-ab)$ generator
at $N=24$ so the answer stays self-contained, the Bogoliubov action sends both combinations to
$e^{-r}$ times themselves, so the sum is $2e^{-2r}$, below the bound for every $r>0$; on the
coherent product $D(\alpha)\otimes D(\beta)\vert0,0\rangle$ with $\alpha,\beta$ kept symbolic,
displacement moves means and not variances, so the sum is exactly 2 for all $\alpha,\beta$: one
violating state and one boundary state, both in closed form. The truncated-matrix evaluation must
reproduce both numbers, the $r\to0$ limit of the squeezed sum must land exactly on the
coherent-product value 2, gluing the two computations together, and a sign slip in the two-mode
generator would squeeze the complementary combinations and push the sum to $2e^{+2r}>2$, reporting
no violation at any $r$, so the check discriminates conventions as sharply as it discriminates
states. Close on tightness: a product state sitting exactly on the boundary shows the constant 2
cannot be improved, and this same combination is the witness measured on twin beams in the
laboratory.

### 17.6 [MSc] How do I solve the Jaynes-Cummings model (one mode plus a two-level atom, truncated Fock) and exhibit vacuum Rabi oscillations?

C5 per `Route-Table.md`. On resonance the initial state $\vert e,0\rangle$ closes into the
two-dimensional block $\{\vert e,0\rangle,\vert g,1\rangle\}$, and `DSolveValue` on the
two-amplitude system returns the exact pair $\{\cos gt,\,-i\sin gt\}$ (the probed C5 gate), so
$P_e(t)=\cos^2gt$: vacuum Rabi oscillation at frequency $2g$ out of a field containing zero
photons. Generalize with detuning through the same `DSolveValue`,
$P_e(t)=1-\frac{4g^2}{4g^2+\Delta^2}\sin^2(\tfrac12\sqrt{4g^2+\Delta^2}\,t)$, oscillating at
$\sqrt{4g^2+\Delta^2}$ with peak transfer probability $4g^2/(4g^2+\Delta^2)$, swept over $\Delta$.
The certified numeric route is the C5 primary: assemble the full truncated Jaynes-Cummings
Hamiltonian with `SparseArray` and `KroneckerProduct` (mode times atom) and evolve with the
`MatrixExp` action form. From $\vert e,0\rangle$ the rotating-wave coupling conserves excitation
number, so the state never leaves the one-excitation manifold, the truncated ladder is never asked
for $a^\dagger\vert N\rangle$, and the $N$-sweep is exactly flat: keep it, but label it the
assembly check it actually is, since all it can detect is an excitation-non-conserving error. The
truncation trap is exercised by a second run that does leave that manifold, an initial coherent
field $\vert g,\alpha\rangle$ with $\vert\alpha\vert^2\approx9$ showing collapse and revival near
$t\approx2\pi\sqrt{\bar n}/g$, where $N$ must exceed $\vert\alpha\vert^2$ plus several
$\vert\alpha\vert$ and the sweep genuinely moves; per the C5 verdict route agreement never
certifies truncation (two routes agreed to $1.1\times10^{-7}$ while both were wrong by
$5.6\times10^{-3}$, a blindness factor $5.1\times10^4$), so only the sweep or an exact benchmark
can. Close in the cavity: the $2g$ oscillation is the single-atom vacuum Rabi frequency of cavity
QED, the revivals are the discreteness of the photon number made visible in an atomic signal, and
the dispersive limit $\Delta\gg g$ shrinks the transfer to $4g^2/\Delta^2$, the regime qubit
readout lives in.

### 17.7 [MSc] How do I model balanced homodyne and heterodyne detection, measuring a quadrature and reconstructing the Husimi-$Q$ distribution operationally?

Take the squeezed vacuum $S(r)\vert0\rangle$ at $r=1$ on a signal mode truncated at $N_a=40$, where
the neglected weight sits near $10^{-4}$; at $N=12$ it is concentrated at high $n$ and the
anti-squeezed variance falls about 4 percent short of $e^{2r}/2$, ten times the statistical
tolerance quoted below. Derive by `Integrate` (assumptions $r>0$, $\theta$ real) the closed-form
marginal of the state's Wigner Gaussian along $X_\theta$, variance
$\sigma_\theta^2=\tfrac12(e^{-2r}\cos^2\theta+e^{2r}\sin^2\theta)$, sub-vacuum at $\theta=0$ under
4.7's convention, where a flipped squeeze sign would show up loudly as $\langle X^2\rangle=e^{2r}/2$
instead; that marginal is the benchmark, never the sampler. Homodyne is built as apparatus: put a
local oscillator in the coherent state $\beta=\vert\beta\vert e^{i\theta}$ on mode $b$, combine it
with the signal on 17.3's $B(\pi/4)$ by `MatrixExp` in action form, and detect photon number in
both output ports; the difference $n_c-n_d$ pulled back through the splitter is
$n_-=a^\dagger b+ab^\dagger$, whose scaled version $n_-/(\sqrt2\vert\beta\vert)$ tends to $X_\theta$
as the local oscillator grows classical. Sampling is then genuinely operational, `RandomChoice` on
the output Fock weights $\vert\Psi_{\mathrm{out}}\vert^2$ followed by binning
$(n_c-n_d)/(\sqrt2\vert\beta\vert)$, with no appeal to the answer. The refuting check is
quantitative and exact: the scaled variance is $\sigma_\theta^2+\sinh^2r/(2\vert\beta\vert^2)$, so
a sweep over $\vert\beta\vert\in\{2,4,8\}$ (sizing $N_b$ beyond $\vert\beta\vert^2$ plus several
$\vert\beta\vert$, and $N_a$ beyond the $\vert\beta\vert^2/2$ photons each output port then carries)
must show the excess falling as $1/\vert\beta\vert^2$ with that coefficient while the histogram
converges onto the closed marginal, and an independent $N$-sweep at $\theta=\pi/2$, where
truncation bites hardest, must leave the anti-squeezed variance stationary. Heterodyne adds the
vacuum port: splitting the signal against vacuum and reading both quadratures samples
$Q=W\ast W_{\mathrm{vac}}$ with covariance $\sigma_\theta^2+\tfrac12$, cross-checked against
$Q(\alpha)=\vert\langle\alpha\vert\psi\rangle\vert^2/\pi$ evaluated directly from the truncated Fock
amplitudes. Close on the contrast the question asks for: histograms of detector clicks converge to
the formal objects (the Wigner marginal, the $Q$ function), heterodyne adds one full vacuum unit of
variance per quadrature ($\tfrac12$ in these units, half a photon per mode) and so doubles the
vacuum noise, the 3 dB penalty for reading both knobs at once, and the $\theta$-resolved homodyne
variances are the raw data of quantum state tomography.
