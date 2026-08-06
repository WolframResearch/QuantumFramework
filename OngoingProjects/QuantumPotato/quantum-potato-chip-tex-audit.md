# Quantum Potato Chip: full audit of `main.tex` and desk-rejection diagnosis

Audited 2026-08-05. Sources read in full: `quantum_potato_chip.zip/main.tex` (759 lines, submission title "Geometric Factorization of Qubit States into Separable Binary Distributions"), all 15 figure PNGs in the zip, and `QuantumPotato V2.nb` (via nb-reader). Every displayed equation was re-derived independently (sympy/numpy, exact arithmetic; scripts in session scratchpad `audit_checks*.py`). Public record: [arXiv:2411.01082](https://arxiv.org/abs/2411.01082) "Quantum Potato Chips", v1 of 2024-11-01, never revised, 0 citations on INSPIRE as of today.

The single most consequential fact found: **the notebook's computed results are correct, and almost every mathematical error in the TeX is a hand-transcription error from the notebook.** Rebuilding the TeX equations from the notebook output fixes E2 and E4 through E8 below verbatim.

---

## 1. Confirmed mathematical errors

### E1. The fiducial vector does not generate the printed SIC-POVM (main.tex:171)
> "generated from a fiducial vector ... as $\{e^{-i (3\pi)/4}(\sqrt{3} - 1), 1\}$"

Normalizing this vector and forming $\tfrac12|\psi\rangle\langle\psi|$ does not give $\mathcal{Q}_1$; its Bloch direction is not even a vertex of the SIC tetrahedron. The correct fiducial is $\{e^{-3i\pi/4}(\sqrt3-1),\ \sqrt2\}$ (second component $\sqrt2$, not 1). With that vector, the Weyl-Heisenberg orbit $X^pZ^q$, $(p,q)=(0,0),(0,1),(1,0),(1,1)$, reproduces $\mathcal{Q}_1,\mathcal{Q}_2,\mathcal{Q}_3,\mathcal{Q}_4$ exactly in the printed order (verified). The notebook never derives the fiducial (it pulls the POVM from the framework), so this error exists only in the TeX.

### E2. The fourth Wootters basis element is not Hermitian (main.tex:431-436, eq. wootters-matrix)
Printed fourth vector: $(0,\ e^{i\pi/4}/\sqrt2,\ -e^{-i\pi/4}/\sqrt2,\ 1)$. Read as a $2\times2$ matrix this is non-Hermitian, the four elements do not sum to $2\mathbb{I}$, and $\rho=\sum_i w_i A_i$ fails to reproduce the Bloch state. The second component needs a minus sign: $(0,\ -e^{i\pi/4}/\sqrt2,\ -e^{-i\pi/4}/\sqrt2,\ 1)$, i.e. $\tfrac12(\mathbb{I}-\sigma_x+\sigma_y-\sigma_z)$. With that sign, everything checks (sum $=2\mathbb{I}$, reconstruction exact). The notebook's Eq. 24 prints the correct element $(0,\ e^{-3i\pi/4}/\sqrt2,\ e^{3i\pi/4}/\sqrt2,\ 1)$; note $e^{-3i\pi/4}=-e^{i\pi/4}$, so the TeX dropped the sign in transcription.

### E3. Amplitude-damping Kraus operator is wrong (main.tex:533, Table 1)
> $\{\frac{(1-\sqrt{\xi})}{2} \sigma_z+\frac{(1+\sqrt{1-\xi})}{2} \mathbb{I},\ \ldots\}$

With $(1-\sqrt{\xi})/2$ the pair violates $\sum_k K_k^\dagger K_k=\mathbb{I}$ (verified nonzero residual). It must be $(1-\sqrt{1-\xi})/2$, the same pattern the phase-damping row directly below uses correctly, giving $K_0=\mathrm{diag}(1,\sqrt{1-\xi})$.

### E4. Table 2, bit-phase-flip row: wrong sign on the x component (main.tex:553)
Printed $\sqrt3\{-f_q f_\xi,\ f_p f_q,\ f_p f_\xi\}$; the actual channel action on the chip gives $+\sqrt3 f_q f_\xi$ for the first component (verified at four exact values of $\xi$; the notebook's table has the correct $+$ sign). Compare the phase-flip row, whose x component is printed correctly as $+f_qf_\xi$; the $\sigma_y$ channel acts identically on $x$.

### E5. Table 2, amplitude-damping row: wrong damping power on the x component (main.tex:555)
Printed $-f_q(1-\xi)$; the actual transverse contraction is $-f_q\sqrt{1-\xi}$ (the y component right next to it is printed with the correct $\sqrt{1-\xi}$). The notebook's $(1-2q)\sqrt{3-3\xi}$ is right.

### E6. The Bloch-sphere border equation covers half the border and one branch is off the chip (main.tex:286, eq. bloch-chip-border)
Printed: $\vec r=\{-\sqrt{A},\ (\mp1\pm2p)\sqrt{A},\ (2p-1)\sqrt3\}$ with $A = \frac{2}{1+2p(p-1)}-3$.

The x component must carry the $\mp$ correlated with the y sign. As printed, x is always negative, so the $x>0$ half of the border curve is never produced, and the second sign branch satisfies $|\vec r|=1$ but lies on the mirror surface $y=-xz/\sqrt3$, not on the chip (numerically: at $p=0.3$ the printed upper branch is on the unit sphere but violates $y=xz/\sqrt3$). The notebook's Eq. 16 output is correct and complete: $\{\{+\sqrt{A},(1-2p)\sqrt{A},\sqrt3(1-2p)\},\{-\sqrt{A},(2p-1)\sqrt{A},\sqrt3(1-2p)\}\}$. The TeX also flipped the z component's sign relative to the notebook while re-typing. Figure 3(b)'s red curve visibly spans both hemispheres, so the figure contradicts the printed formula.

### E7. The M_x and M_z statements are swapped, and contradict the paper's own worked example
- main.tex:376: "the measurement probabilities for $\mathcal{M}_z$ is $\{q,1-q\}$ and for $\mathcal{M}_x$ is $\{p,1-p\}$". Actual: $P(\mathcal{Q}_1+\mathcal{Q}_2)=pq+p(1-q)=p$, so $\mathcal{M}_z\leftrightarrow\{p,1-p\}$ and $\mathcal{M}_x\leftrightarrow\{q,1-q\}$.
- main.tex:387: "Adding up columns in Eq.(sic-density) return $\mathcal{M}_z$ probabilities, while for rows summation one gets $\mathcal{M}_x$". Actual: column sums give $\tfrac12(1\mp x/\sqrt3)$, the $\mathcal{M}_x$ probabilities; row sums give $\mathcal{M}_z$. (The Wootters analog at main.tex:450 is stated correctly, which makes the SIC section internally inconsistent with it.)
- The worked example at main.tex:403 ($p=1/3$: $\mathcal{M}_z$ gives $\{1/3,2/3\}$) uses the correct assignment, so the draft contradicts itself between lines 376 and 403. The notebook has all of this right; the TeX swapped the labels when transcribing the list `{mx, mz}`.

### E8. The second Lindblad operator is the reciprocal of the correct one (main.tex:673)
Printed: $L_2 = \frac12\sqrt{\frac{p - 3 p^2 + 4 p^3 - 2 p^4}{1-2p}}\, \sigma_x$. The trajectory requires the inverted fraction, $L_2 = \frac12\sqrt{\frac{1-2p}{p - 3 p^2 + 4 p^3 - 2 p^4}}\, \sigma_x$. Verified numerically: with the printed $L_2$ the master equation's right side misses $d\rho/dp$ by a large residual at $p=0.3$ and $p=0.7$; with the inverted fraction the residual vanishes to machine precision at both. The notebook's Eq. 34 stores exactly the inverted (correct) form. Additionally main.tex:668 writes $\gamma_1$ inside the sum over $j$; it must be $\gamma_j$.

### E9. The CNOT-permutation footnote pairs the wrong permutations (main.tex:154)
> "Any two potato chips are related to the third one by two different ways of applying CNOT matrix (corresponding to permutations $\sigma_1$ [swap 1,2] and $\sigma_2$ [swap 3,4])"

Checked over all 24 permutations of the outcome labels: 8 permutations fix the first chip as a set, 8 produce the second printed chip, 8 produce the third. Both printed permutations, (12) and (34), lie in the same class: they produce the same second chip. The third chip requires (13) or (24), i.e. CNOT with control on the other bit. So the correct statement is: standard CNOT gives one of the other chips, the reversed CNOT (control and target exchanged) gives the remaining one.

### E10. Eq. (3new-M) lists two of the three Bloch pairs in the wrong order, and its label is duplicated (main.tex:332-372)
$\mathcal{Q}_1+\mathcal{Q}_3$ has Bloch vector $(-1/\sqrt3,0,0)$ and $\mathcal{Q}_1+\mathcal{Q}_2$ has $(0,0,-1/\sqrt3)$, but the display lists the $+$ vectors first for $\mathcal{M}_x$ and $\mathcal{M}_z$ (while $\mathcal{M}_y$, whose first element genuinely is $+1/\sqrt3$, is listed correctly). The ordering matters because the text immediately identifies first elements with $\{p,\ldots\}$ probabilities. Both displays at main.tex:327 and main.tex:370 carry the same `\label{eq:3new-M}` (LaTeX multiply-defined label; \ref resolution is then arbitrary). Figure 3(c)'s legend also labels the measurements $\mathcal{M}_1,\mathcal{M}_2,\mathcal{M}_3$ while the text uses $\mathcal{M}_x,\mathcal{M}_y,\mathcal{M}_z$.

### E11. The insphere footnote rescales in the wrong direction and conflates two spheres (main.tex:131)
The image of $\{1,0,0,0\}$ has norm $\sqrt3/2$, so normalizing it means dividing by $\sqrt3/2$ (multiplying by $2/\sqrt3$); "rescaled by a factor of $\sqrt{3}/2$" says the opposite. And the radius-1 sphere so obtained passes through the tetrahedron's vertices (circumsphere), while the physical region argued for in the text is the insphere; "analogous to the Bloch sphere" blurs the two. In Bloch terms the simplex vertices sit at radius 3, not 1.

---

## 2. Verified correct (independent re-derivation, exact)

- Eq. (1) is a genuine rotation (orthogonal, det 1) in the plane spanned by $\{1,1,1,1\}$ and $\{1,0,0,0\}$; Eq. (2) is its $\theta=\pi/3$ case; the stated action on $\{p_1,p_2,p_3,1-\Sigma\}$ (main.tex:95) is exact.
- Image tetrahedron centered at the origin; insphere radius $1/(2\sqrt3)$; the insphere is exactly the image of the Bloch ball, so "physical = insphere" is correct for qubits.
- Eq. (kp-vector) maps to Eq. (chip1) exactly; the surface lies in the tetrahedron; Eq. (constraint) solves $|s|^2=1/12$ on both branches; the p-range in Eq. (boundary) is exact.
- Eq. (povm) is a valid SIC: $\sum_i\mathcal{Q}_i=\mathbb{I}$, each rank-1 with trace 1/2, $\mathrm{tr}(\mathcal{Q}_i\mathcal{Q}_j)=1/12$ for $i\ne j$, Bloch vectors as listed at main.tex:196.
- Eq. (qbism): $\mathcal{B}_i=3\Pi_i-\mathbb{I}$ with the duality $\mathrm{tr}(\mathcal{B}_i\mathcal{Q}_j)=\delta_{ij}$; Eq. (density-phase-space) equals $\mathcal{B}\cdot\vec p$; the eigenvalues are $\tfrac12(1\pm\sqrt{3\omega})$ with the printed $\omega$; the replacement rule and its inverse are exact; $3\omega=x^2+y^2+z^2$, so positivity is exactly the Bloch ball.
- The chip's Bloch surface is $\sqrt3\{1-2q,\ (2p-1)(2q-1),\ 1-2p\}$ (main.tex:591), equivalently $y=xz/\sqrt3$ inside the ball.
- Eq. (sic-density) is the probability 4-vector arranged as a matrix; Eq. (stoch-shrink) implements $P\mapsto (P-\tfrac12)/\sqrt3+\tfrac12$; every number in the $p=1/3,q=2/5$ example (main.tex:403) is exact, including the outer-product reconstruction.
- Wootters: Eq. (bloch-in-wooters) with correct Z/X marginals; the Wootters chip is $x=2p-1$, $z=2q-1$, $y=xz$; Eq. (constraint-wootters) is exactly the purity condition; the tetrahedron/ball non-containment statements in the Fig. 5 caption are correct.
- Matthews: Eq. (matthews) follows exactly from the SIC probabilities; $\phi=0$ inside the ball is exactly the chip; extremes $\pm1/\sqrt3$ on the sphere at the $\pm y$ poles; the rescaled form (main.tex:513) is exact under $r\to r/\sqrt3$.
- Channels: bit-flip, phase-flip, depolarizing and phase-damping rows of Table 2 are exact; bit-flip, phase-flip and phase-damping map the chip into itself with exactly the reparametrizations claimed at main.tex:593; bit-phase-flip, depolarizing and amplitude damping leave the chip surface. All Kraus lists except amplitude damping (E3) are complete.
- Liouvillian: $x=1/(1-2p)$ and the printed $y(p)$ are the unique marginal generator rates on both branches; the printed $4\times4$ $\mathcal{L}$ equals the Kronecker sum exactly; the Shannon entropy along the border satisfies $1\le H\le 1.35226$ (grid max 1.352260, endpoints exactly 1).

---

## 3. Conceptual overclaims (what a referee or editor rejects on substance)

### C1. The headline uniqueness claim is false as stated
Abstract (main.tex:54): "States within this region can be fully reconstructed using only two projective measurements, a feature not found elsewhere in the state space." Repeated at main.tex:68, 319, 401.

Any fixed two-parameter family of states with a known chart is determined by two measured numbers. The equatorial disk $y=0$ (real density matrices) is fully determined by the same two measurements, X and Z, given the same kind of promise. So "not found elsewhere" fails in one minute for any editor who knows qubit tomography. What actually is unique to the chip, and what the paper proves without stating it as the theorem: the chip is exactly the locus where the SIC probability 4-vector factorizes into the product of its two binary marginals (equivalently, Matthews correlation zero; equivalently $y=xz/\sqrt3$). Reconstruction-by-outer-product is a corollary that additionally requires the promise that the state lies on a known chip; the abstract omits the promise.

### C2. The two-binary-variables claim is backwards (main.tex:68)
> "any classical problem with two binary variables can be mapped into qubits"

Only independent pairs map onto the chip. A correlated pair, e.g. $\{1/2,0,0,1/2\}$, sits at distance $1/2 > 1/(2\sqrt3)$ from the simplex center: outside the insphere, hence not a physical state in this representation at all. The chip encodes exactly the zero-correlation pairs; correlated distributions are the ones this construction cannot host.

### C3. The $\mathcal{M}_i$ are unsharp POVMs, not projective measurements
Their effects have Bloch length $1/\sqrt3$ (noisy Pauli measurements). The operational protocol is: measure projective X and Z, then rescale by Eq. (stoch-shrink). Calling the $\mathcal{M}_i$ themselves "projective" invites an easy referee objection; the fix is one sentence of wording.

### C4. The Liouvillian section solves a problem whose standard solution is trivial, without saying why
The border consists of pure states, and any differentiable curve of pure states is generated exactly by a time-dependent Hamiltonian, $\vec\omega = \vec r\times\dot{\vec r}$: no dissipators, no negative rates, no $p=1/2$ singularity. Demanding instead that the generator be a Kronecker sum of two classical marginal generators is a structural choice the text never motivates. Negative time-local rates are also standard in the non-Markovian master-equation literature, so "unconventional" both undersells the known context and oversells the novelty. The section additionally switches, without comment, from the SIC chip used everywhere before to the Wootters chip (Eq. constraint-wootters).

### C5. Unsupported computation claims in the conclusion (main.tex:693)
"might provide an advantage for these states as a potential source for any computation" and "how any [universal] quantum computation can be reduced to only probabilistic rules, in a classical way, if possible": speculation with no construction, no bound, no citation; the bracketed "[universal]" survives in the submitted text. Three separate passages defer all applications to future work (main.tex:68, 405 comment, 693).

---

## 4. Prose, LaTeX and packaging defects (what a desk editor sees in five minutes)

- main.tex:66, the thesis sentence of the introduction: "This is manifold looks like a potato chip we shall refer to it as {\it quantum potato chips}." Two grammatical breaks in the sentence that names the paper's object.
- main.tex:116 "Any point sampled from this simplex have the form"; main.tex:130 "a three-dimensional surface" for a 2D surface embedded in 3D; main.tex:314 "free of correlation between its complete projective measurements"; main.tex:452 "Wootters phase-space basis is no longer a proper probability distribution" (a basis is not a distribution); main.tex:649 "The opposite direction of inference is what keeps the overall distribution spread (and purity with von Neumann entropy) constant"; main.tex:591 "as shown in the following equation" with no following equation (the content lives in Table 2).
- Terminology collisions: "separable" for statistical independence of a single qubit's marginals (readers will parse entanglement separability; the title inherits this); "Weyl-Wigner transformation" (main.tex:282) for applying the SIC dual frame; "product state" (main.tex:494) for a product distribution.
- LaTeX mechanics: duplicate `\label{eq:3new-M}`; `$\sqrt(3)$` twice in the Fig. 7 caption (main.tex:506); `revtex4-2` with `prl` option plus a `geometry` a4paper override for PRA/FoP submissions; an IOP-style "Data availability statement" inside a REVTeX/APS manuscript; a MathWorld URL as a footnote (main.tex:390); "Here is the complete Mathematica notebook: wolfr.am/QPC" as the code statement (main.tex:700); 50 lines of commented-out duplicate tables left in the source (main.tex:567-615); `subfigure` (deprecated) alongside `natbib` with a hand-built `thebibliography`; an unused `Pringles.png`, a render of the trademarked product, shipped in the submission bundle.
- Reference list: 10 entries, newest 2017, five of them SIC surveys. Nothing from the three directly adjacent literatures:
  - discrete-Wigner classicality for qubits: [Galvão, PRA 71, 042302 (2005)](https://arxiv.org/abs/quant-ph/0405070) and [Cormick, Galvão, Gottesman, Paz, Pittenger, PRA 73, 012301 (2006)](https://arxiv.org/abs/quant-ph/0506222); the Fig. 5 tetrahedron of nonnegative-Wigner states is exactly their object (the draft cites only Appleby-Ericsson-Fuchs 2011 in a caption);
  - Kirkwood-Dirac positivity, an actively publishing 2024-2026 community whose objects (states with nonnegative joint X-Z quasiprobability, its marginals, classical simulability) overlap the Wootters section directly, e.g. [KD nonpositivity as a resource](https://arxiv.org/abs/2506.08092) and the [KD review](https://en.wikipedia.org/wiki/Kirkwood-Dirac_quasiprobability);
  - probability representations of qubit states (Man'ko school: a qubit as three coin probabilities, e.g. [arXiv:1803.09339](https://arxiv.org/abs/1803.09339)), and minimal-measurement state determination (the Pauli problem literature).

---

## 5. Why it desk-rejects (direct answer to "what is wrong with the draft?")

1. **The abstract's central claim is refutable by a specialist in one minute** (C1). An editor who knows that two projective measurements fix any promised two-parameter family reads "a feature not found elsewhere in the state space" as wrong, and stops there. The genuinely new object (the factorization locus) is never stated as the theorem.
2. **The reference list signals no research conversation.** Ten references, none newer than 2017, zero contact with the three communities that own the neighboring results (discrete-Wigner classicality, KD positivity, probability representations). A desk editor uses the bibliography to find the audience and the referee pool; this one names neither. The arXiv record (Nov 2024, no revisions, no citations) reinforces the impression.
3. **The manuscript reads as a notebook export, not a journal article**: broken thesis sentence in the introduction, eight transcription errors in displayed equations, mixed APS/IOP boilerplate, a marketing name ("quantum potato chip") doing the work of a mathematical one, a URL-shortener link as the code statement. Any spot-check by an editor or referee hits one of these within minutes.
4. **Scope vs venue mismatch.** For PRA: a single-qubit, elementary-geometry observation with no operational task, no bound, no qudit statement, and applications deferred to future work three times does not clear "significant advance". For Foundations of Physics: the paper gestures at foundations (QBism naming, quasi-probabilities) but engages none of the foundations literature (no Spekkens/epistricted theories, no KD, no contextuality), so it does not read as a foundations paper either.

None of this says the content is wrong-headed: the verified core (Sections 1-2 here) is a clean, correct piece of state-space geometry. It is packaged as more than it is, referenced as less than it needs, and typeset with errors the notebook does not have.

## 6. What to change before resubmitting anywhere

1. Regenerate every displayed equation from the notebook; it is the correct source (fixes E2, E4-E8 verbatim; E1, E3, E9-E11 need small derivations).
2. State the real theorem: the chip is exactly the set of qubit states whose SIC (equivalently Wootters) distribution factorizes into independent binary marginals, equivalently the zero set of the Matthews correlation; make two-measurement reconstruction a promise-based corollary and delete "not found elsewhere".
3. Position against the three adjacent literatures (one paragraph each, 25-40 references total), and say explicitly what is new relative to Galvão/Cormick-et-al. and KD positivity: those characterize nonnegativity; this draft characterizes factorizability, a strictly smaller and different locus. That one contrast sentence is the paper's actual novelty claim, and it currently appears nowhere.
4. Cut the Liouvillian section, or open it by conceding the Hamiltonian solution and motivating the classical-marginal-generator constraint as the point.
5. Add the $d\ge3$ question, even as one honest open paragraph (does any factorization surface survive for qutrit SICs?). This is the cheapest way to turn a curiosity into a research direction.
6. One professional language pass; one journal template; "potato chip" demoted to a parenthetical nickname.
7. Venue after the fix: J. Phys. A (geometry of quantum state space is squarely in scope), Physica Scripta, Quantum Studies: Mathematics and Foundations, or Entropy. PRA only if item 5 produces an actual result, or an operational task is added.
