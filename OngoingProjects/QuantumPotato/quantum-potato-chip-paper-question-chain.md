# The quantum potato chip: a formal question chain

*This chain compiles the 2024 source paper (arXiv:2411.01082; `main.tex` in `quantum_potato_chip.zip`) into pinned definitions, lemmas, and problems, in the paper's own conventions. The sibling chain [quantum-potato-chip-question-chain.md](quantum-potato-chip-question-chain.md) formalizes the two-coins concept of [quantum-potato-chip-concept-draft.md](quantum-potato-chip-concept-draft.md).*

## Preamble

This document compiles the quantum-potato-chip idea into named definitions, certified lemmas, and problems whose inputs are entirely contained below. Every formal item uses only the declared primitives and earlier chain items; every convention that changes the object is named. The problems contain no solutions: the central derivation must be executable from this file alone.

**Audit status.** Round 0 reconstructed the algebra independently from the TeX source and checked it with exact symbolic computation. Three fresh-context rounds each verified every certified lemma and solved P1 from this file alone. Round 3 found one normalization ambiguity in P3, repaired below; under the three-round cap, that last repair received a local rule/type check but no fourth independent audit, so the document claims zero known open issues rather than an independently certified clean pass.

---

## Primitives

The ambient scalars are $\mathbb{R}$ and $\mathbb{C}$. Let $M_2(\mathbb{C})$ denote the $2\times2$ complex matrices, with identity $\mathbb{I}$, adjoint $\dagger$, trace $\operatorname{tr}$, determinant $\det$, matrix rank, and positive-semidefinite order $\succeq0$. Let $\lVert\cdot\rVert$ denote the Euclidean norm on $\mathbb{R}^3$.

The Pauli matrices are

$$
\sigma_x=\begin{pmatrix}0&1\\1&0\end{pmatrix},\qquad
\sigma_y=\begin{pmatrix}0&-i\\i&0\end{pmatrix},\qquad
\sigma_z=\begin{pmatrix}1&0\\0&-1\end{pmatrix}.
$$

For $\mathbf r=(x,y,z)\in\mathbb{R}^3$, write

$$
\mathbf r\boldsymbol\cdot\boldsymbol\sigma:=x\sigma_x+y\sigma_y+z\sigma_z.
$$

The binary simplex and four-outcome simplex are

$$
\Delta_1:=\{(u,1-u):0\leq u\leq1\},
$$

$$
\Delta_3:=\left\{(v_1,v_2,v_3,v_4)\in\mathbb{R}_{\geq0}^4:\sum_{i=1}^4v_i=1\right\}.
$$

The rotation group of $\mathbb{R}^3$ is $SO(3)$. A qubit channel means a completely positive trace-preserving linear map on $M_2(\mathbb{C})$. For a state matrix $\rho$ and a finite ordered POVM $\mathsf E=(E_i)_i$, the Born outcome distribution is

$$
\Pr_\rho(\mathsf E):=\bigl(\operatorname{tr}(\rho E_i)\bigr)_i.
$$

*Reading.* The ambient ingredients are one qubit, one tetrahedral four-outcome probability representation, and the classical probability spaces of one bit and two bits. No tensor-product Hilbert space is present.

---

## Question-chain rules

**R1. Declared symbols and exact dependencies.** Every formal symbol must be declared before use. A **Uses** line lists the earlier chain items directly referenced by that item’s definition, lemmas, or problem statement; ambient primitives are omitted, and transitive prerequisites are not repeated.

**R2. Pinned constants and conventions.** Every constant, orientation, outcome order, normalization, and characterization value that changes a result must be explicit.

**R3. Valid quantifier ranges.** Every quantified family must have a stated range, and any restricted existence region must be established before it is used.

**R4. Type correctness.** Scalars, real vectors, probability vectors, matrices, states, measurements, and channels must not be interchanged without a declared map or lemma.

**R5. Collision-free notation.** A symbol has one role throughout the chain; visually confusable or overloaded symbols must be renamed.

**R6. No solutions in the problem set.** Definitions and lemmas may supply inputs, but the document must not contain solutions to P1–P6.

**R7. Jargon stripping.** Every technical label in a question must expand to an explicit mathematical condition or deliverable. Readings are interpretive only and carry no datum needed to solve a problem.

---

## The chain

### D1. Qubit state space

Define

$$
\mathcal S:=\{\rho\in M_2(\mathbb{C}):\rho\succeq0,\ \operatorname{tr}\rho=1\}.
$$

**Lemma 1a (substrate).** The map

$$
\mathbf r\longmapsto\rho(\mathbf r):=\frac12\left(\mathbb I+\mathbf r\boldsymbol\cdot\boldsymbol\sigma\right)
$$

is a bijection from the closed unit ball

$$
\mathbb B:=\{\mathbf r\in\mathbb R^3:\lVert\mathbf r\rVert\leq1\}
$$

onto $\mathcal S$. Its inverse is

$$
\rho\longmapsto
\bigl(\operatorname{tr}(\rho\sigma_x),\operatorname{tr}(\rho\sigma_y),\operatorname{tr}(\rho\sigma_z)\bigr).
$$

**Uses:** none.

*Reading.* A qubit has three real coordinates, but positivity confines them to the Bloch ball. The sphere is the set of pure states; the interior contains mixed states.

### D2. Pinned tetrahedral SIC measurement

Define four unit vectors

$$
\begin{aligned}
\mathbf n_1&=\frac1{\sqrt3}(-1,1,-1),&
\mathbf n_2&=\frac1{\sqrt3}(1,-1,-1),\\
\mathbf n_3&=\frac1{\sqrt3}(-1,-1,1),&
\mathbf n_4&=\frac1{\sqrt3}(1,1,1),
\end{aligned}
$$

and effects

$$
Q_i:=\frac14\left(\mathbb I+\mathbf n_i\boldsymbol\cdot\boldsymbol\sigma\right),
\qquad i\in\{1,2,3,4\}.
$$

**Lemma 2a (certification).** Each $Q_i$ is positive semidefinite with rank one,

$$
\sum_{i=1}^4Q_i=\mathbb I,
$$

and, for $i\neq j$,

$$
\operatorname{tr}(Q_iQ_j)=\frac1{12}.
$$

**Lemma 2b (interpretation).** The directions $\mathbf n_i$ are the vertices of a regular tetrahedron centered at the origin; replacing every $\mathbf n_i$ by $R\mathbf n_i$ for a fixed $R\in SO(3)$ gives another tetrahedral SIC measurement.

**Uses:** none.

*Reading.* This fixes the convention that makes a “potato chip” a definite subset of the Bloch ball. Rotating the tetrahedron rotates the chip; relabeling the four outcomes can change which two classical bits are being declared independent.

### D3. SIC probability representation

For $\rho\in\mathcal S$, define

$$
T(\rho):=\bigl(\operatorname{tr}(\rho Q_1),\operatorname{tr}(\rho Q_2),
\operatorname{tr}(\rho Q_3),\operatorname{tr}(\rho Q_4)\bigr).
$$

**Lemma 3a (certification).** $T(\rho)\in\Delta_3$ for every $\rho\in\mathcal S$.

**Lemma 3b (substrate).** For $\rho=\rho(x,y,z)$,

$$
T(\rho)=\frac14\left(
1+\frac{-x+y-z}{\sqrt3},
1+\frac{x-y-z}{\sqrt3},
1+\frac{-x-y+z}{\sqrt3},
1+\frac{x+y+z}{\sqrt3}
\right).
$$

**Uses:** D1, D2.

*Reading.* The SIC turns a density matrix into an ordinary four-outcome probability vector. Informational completeness means that no information is lost, even though the probabilities occupy only a ball-shaped part of the full tetrahedral simplex.

### D4. Physical probability body

Define

$$
\mathcal K:=T(\mathcal S)\subset\Delta_3.
$$

**Lemma 4a (substrate).** A vector $\mathbf v=(v_1,v_2,v_3,v_4)\in\Delta_3$ belongs to $\mathcal K$ if and only if

$$
\sum_{i=1}^4v_i^2\leq\frac13.
$$

On $\mathcal K$, the inverse of $T$ has Bloch vector

$$
\begin{aligned}
x&=\sqrt3(1-2v_1-2v_3),\\
y&=\sqrt3(1-2v_2-2v_3),\\
z&=\sqrt3(1-2v_1-2v_2).
\end{aligned}
$$

**Lemma 4b (interpretation).** Pure states correspond to $\sum_i v_i^2=1/3$; the maximally mixed state corresponds to $(1/4,1/4,1/4,1/4)$.

**Uses:** D1, D3.

*Reading.* The full simplex contains every classical distribution on four labels, but only the inscribed quadratic body $\mathcal K$ represents qubit states. This physicality restriction is what cuts a finite potato-chip-shaped patch from an otherwise larger product surface.

### D5. The three two-bit pairing determinants

For $\mathbf v=(v_1,v_2,v_3,v_4)\in\mathbb R^4$, define

$$
\begin{aligned}
\delta_1(\mathbf v)&:=v_1v_4-v_2v_3,\\
\delta_2(\mathbf v)&:=v_1v_3-v_2v_4,\\
\delta_3(\mathbf v)&:=v_1v_2-v_3v_4.
\end{aligned}
$$

**Uses:** none.

*Reading.* Four labels can be placed into a $2\times2$ table in three inequivalent ways, after ignoring row swaps, column swaps, and interchange of the two bits. Each $\delta_k$ is the determinant for one such table.

### D6. Independent binary distributions

Define the product map

$$
F:\Delta_1\times\Delta_1\longrightarrow\Delta_3,
$$

$$
F\bigl((p,1-p),(q,1-q)\bigr)
:=\bigl(pq,\ p(1-q),\ (1-p)q,\ (1-p)(1-q)\bigr).
$$

Write $F(p,q)$ for this value when the two binary distributions are clear.

**Lemma 6a (certification).** For $\mathbf v\in\Delta_3$, the following are equivalent:

1. $\mathbf v=F(p,q)$ for some $(p,q)\in[0,1]^2$;
2. $\delta_1(\mathbf v)=0$;
3. the table $\begin{pmatrix}v_1&v_2\\v_3&v_4\end{pmatrix}$ has rank one.

When these conditions hold, the parameters are fixed by the marginals:

$$
p=v_1+v_2,\qquad q=v_1+v_3.
$$

**Uses:** D5.

*Reading.* This is the precise classical content of the construction: the four SIC outcome probabilities behave as the joint distribution of two independent biased bits. The two binary simplices are combined by a Cartesian product; they are not “disjoint,” and they are not two quantum subsystems.

### D7. The quantum potato chips

For $k\in\{1,2,3\}$, define

$$
\mathcal C_k:=\{\rho\in\mathcal S:\delta_k(T(\rho))=0\}.
$$

Call $\mathcal C_1$ the reference quantum potato chip and $\{\mathcal C_1,\mathcal C_2,\mathcal C_3\}$ its pairing orbit.

**Lemma 7a (certification).** The maximally mixed state belongs to all three $\mathcal C_k$, so every $\mathcal C_k$ is nonempty.

**Lemma 7b (interpretation).** Membership in $\mathcal C_k$ asserts classical independence of two labels inside one chosen SIC outcome table. It does not assert separability, entanglement, or absence of quantum correlations between subsystems, because D1 contains only one qubit and no subsystem decomposition.

**Lemma 7c (substrate).** For $\mathbf v\in\Delta_3$, define the three tables

$$
\mathsf J_1(\mathbf v):=\begin{pmatrix}v_1&v_2\\v_3&v_4\end{pmatrix},\qquad
\mathsf J_2(\mathbf v):=\begin{pmatrix}v_1&v_2\\v_4&v_3\end{pmatrix},\qquad
\mathsf J_3(\mathbf v):=\begin{pmatrix}v_1&v_3\\v_4&v_2\end{pmatrix}.
$$

Then $\delta_k(\mathbf v)=\det \mathsf J_k(\mathbf v)$, and $\rho\in\mathcal C_k$ if and only if $\mathsf J_k(T(\rho))$ is the joint distribution of two independent binary variables.

**Uses:** D1, D3, D5, D6.

*Reading.* The chip is the intersection of two constraints: the SIC probabilities must come from a physical qubit, and the chosen $2\times2$ table of those probabilities must factor. The name describes the resulting curved patch inside the Bloch ball.

### P1. Derive the reference chip and its pairing orbit

Starting only from D1–D7, give all of the following in closed form:

1. one polynomial $g_1(x,y,z)$ such that
   $$
   \mathcal C_1=\{\rho(x,y,z):(x,y,z)\in\mathbb B,\ g_1(x,y,z)=0\};
   $$
2. a two-parameter Bloch-vector parametrization obtained from $T(\rho)=F(p,q)$;
3. necessary and sufficient inequalities on $(p,q)\in[0,1]^2$ for the parametrized vector to lie in $\mathbb B$;
4. the two relative-boundary branches $q=q_\pm(p)$ of the admissible product-parameter domain and the exact interval of $p$ on which those branches are real;
5. corresponding polynomials $g_2,g_3$ giving the same typed set description for $\mathcal C_2$ and $\mathcal C_3$;
6. the dimension and topology of each $\mathcal C_k$, a closed-form description of its boundary relative to the zero set $g_k=0$ (equivalently, the boundary inherited from its product-parameter domain), and a determination of which relative-boundary states are pure.

The deliverable must show every substitution from the declared formulas and must identify which statements change under outcome relabeling.

**Solvability note.** This is a finite calculation: insert Lemma 3b into D5, or insert D6 into Lemma 4a and then apply Lemma 4a’s inverse map. No undeclared fact is needed.

**Uses:** D1, D3, D4, D5, D6, D7.

*Reading.* Solving P1 produces the object normally shown in the potato-chip plots: its implicit equation, its coin-bias coordinates, its physical boundary, and the other two chips obtained by changing how the four outcomes are paired.

### D8. The two observable marginals

Define two binary POVMs

$$
\mathsf U_z:=\bigl(Q_1+Q_2,\ Q_3+Q_4\bigr),\qquad
\mathsf U_x:=\bigl(Q_1+Q_3,\ Q_2+Q_4\bigr),
$$

and the sharp Pauli POVMs

$$
\mathsf P_j:=\left(\frac12(\mathbb I-\sigma_j),\frac12(\mathbb I+\sigma_j)\right),
\qquad j\in\{x,z\}.
$$

Let

$$
\eta:=\frac1{\sqrt3},\qquad
S_\eta:=\frac12
\begin{pmatrix}
1+\eta&1-\eta\\
1-\eta&1+\eta
\end{pmatrix}.
$$

**Lemma 8a (substrate).** For every $\rho\in\mathcal S$, the outcome distributions satisfy

$$
\Pr_\rho(\mathsf U_j)=S_\eta\Pr_\rho(\mathsf P_j),
\qquad j\in\{x,z\}.
$$

**Lemma 8b (substrate).** If $T(\rho)=F(p,q)$, then

$$
\Pr_\rho(\mathsf U_z)=(p,1-p),\qquad
\Pr_\rho(\mathsf U_x)=(q,1-q).
$$

**Uses:** D1, D2, D3, D6.

*Reading.* Row and column sums of the four-outcome SIC table are two unsharp binary measurements aligned with $z$ and $x$. Ordinary sharp Pauli measurements contain the same two expectation values, and the fixed stochastic matrix $S_\eta$ converts their outcome probabilities to the SIC marginals.

### P2. Promised two-setting reconstruction

Suppose an unknown state $\rho$ is promised to lie in $\mathcal C_1$ and exact outcome distributions for $\mathsf P_x$ and $\mathsf P_z$ are given. Construct a closed-form reconstruction map for $\rho$ and prove it is unique on $\mathcal C_1$. Then remove the promise $\rho\in\mathcal C_1$ and give the full set of states compatible with the same two projective-measurement distributions.

The deliverable must distinguish:

1. the two independent real numbers obtained from the two settings;
2. the nonlinear chip constraint that supplies the otherwise unmeasured Bloch coordinate;
3. the physicality test that rejects incompatible noisy estimates;
4. reconstruction of the SIC table by multiplying its recovered marginals.

**Solvability note.** Lemma 8a recovers the two marginals, Lemma 6a reconstructs the table under the chip promise, and Lemma 4a reconstructs the density matrix.

**Uses:** D4, D6, D7, D8.

*Reading.* Two settings do not tomograph an arbitrary qubit. They tomograph this promised constrained model because the chip equation fixes the third Bloch coordinate. The saving is therefore model-based tomography, not a violation of the three-parameter dimension of general qubit states.

### D9. Classical correlation diagnostic

For $\mathbf v\in\Delta_3$, define row and column marginals

$$
r_0=v_1+v_2,\quad r_1=v_3+v_4,\quad
c_0=v_1+v_3,\quad c_1=v_2+v_4.
$$

When $r_0r_1c_0c_1>0$, define the binary Matthews coefficient

$$
\phi(\mathbf v):=\frac{\delta_1(\mathbf v)}{\sqrt{r_0r_1c_0c_1}}.
$$

**Lemma 9a (certification).** For every $\rho\in\mathcal S$, all four marginals of $T(\rho)$ are strictly positive, so $\phi(T(\rho))$ is defined.

**Lemma 9b (certification).** For $\rho\in\mathcal S$,

$$
\phi(T(\rho))=0\quad\Longleftrightarrow\quad\rho\in\mathcal C_1.
$$

**Uses:** D1, D3, D5, D7.

*Reading.* The chip is exactly the zero-correlation locus for the two chosen SIC outcome bits. This is ordinary correlation between labels of a single four-outcome measurement, not a measure of entanglement or discord.

### P3. Quantify departure from the chip

Derive $\phi(T(\rho(x,y,z)))$ in closed Bloch-coordinate form. Determine its maximum and minimum over $\mathbb B$ and identify all optimizers. At a Bloch vector $\mathbf r_0$ representing a relative-interior point of $\mathcal C_1$, let $\widehat{\mathbf n}(\mathbf r_0)$ be the unique unit vector normal to P1’s zero set $g_1=0$ whose $y$-component is negative, and define

$$
\mathbf r(\varepsilon):=\mathbf r_0+\varepsilon\widehat{\mathbf n}(\mathbf r_0).
$$

Give the coefficient of $\varepsilon$ in the expansion of $\phi(T(\rho(\mathbf r(\varepsilon))))$ at $\varepsilon=0$.

The deliverable must state how $\phi$ can be estimated from four-outcome SIC counts and which part of that interpretation depends on the fixed pairing D5.

**Solvability note.** Substitute Lemma 3b into D9; compactness of $\mathbb B$ guarantees extrema, and their exact values follow from a constrained finite-dimensional calculation.

**Uses:** D1, D3, D5, D7, D9, P1.

*Reading.* P3 turns the exact chip equation into a diagnostic: not merely “on” or “off” the chip, but a signed and normalized measure of how strongly the chosen outcome bits fail to factor.

### D10. Six pinned noise families

For parameters $\xi,\gamma\in[0,1]$ and $s\in[0,1]$, define six affine Bloch-coordinate maps $f_{\mathsf E}:\mathbb R^3\to\mathbb R^3$ by

$$
\begin{aligned}
f_{\mathsf B_\xi}(x,y,z)&=(x,(1-2\xi)y,(1-2\xi)z),\\
f_{\mathsf Z_\xi}(x,y,z)&=((1-2\xi)x,(1-2\xi)y,z),\\
f_{\mathsf Y_\xi}(x,y,z)&=((1-2\xi)x,y,(1-2\xi)z),\\
f_{\mathsf D_s}(x,y,z)&=(sx,sy,sz),\\
f_{\mathsf A_\gamma}(x,y,z)&=(\sqrt{1-\gamma}\,x,\sqrt{1-\gamma}\,y,(1-\gamma)z+\gamma),\\
f_{\mathsf{PD}_\gamma}(x,y,z)&=(\sqrt{1-\gamma}\,x,\sqrt{1-\gamma}\,y,z).
\end{aligned}
$$

For each displayed label $\mathsf E$, define $\mathsf E:M_2(\mathbb C)\to M_2(\mathbb C)$ to be the unique complex-linear map satisfying

$$
\mathsf E(\rho(\mathbf r))=\rho(f_{\mathsf E}(\mathbf r))
$$

for every $\mathbf r\in\mathbb B$.

The channels $\mathsf B_\xi,\mathsf Z_\xi,\mathsf Y_\xi,\mathsf D_s,\mathsf A_\gamma,\mathsf{PD}_\gamma$ are called, respectively, bit flip, phase flip, bit-phase flip, depolarizing contraction, amplitude damping toward $z=1$, and phase damping in the $z$ basis.

**Lemma 10a (certification).** Each $f_{\mathsf E}$ sends $\mathbb B$ into $\mathbb B$, and each corresponding linear map $\mathsf E$ is a qubit channel for the stated parameter range.

**Uses:** D1.

*Reading.* A curved model is useful only if one knows which laboratory imperfections keep the model closed and which force the state out of it. These formulas pin the noise conventions; channel names alone are insufficient because parameterizations vary.

### P4. Classify chip-preserving noise

For each displayed channel

$$
\mathsf E\in\{\mathsf B_\xi,\mathsf Z_\xi,\mathsf Y_\xi,\mathsf D_s,\mathsf A_\gamma,\mathsf{PD}_\gamma\},
$$

determine all parameter values for which

$$
\mathsf E(\mathcal C_1)\subseteq\mathcal C_1.
$$

For every preserving family, derive its induced update of the product parameters $(p,q)$ from D6. Let $\mathcal V_1\subset\mathbb R^3$ be the full real zero set of P1’s polynomial for $\mathcal C_1$. For every non-preserving family, compute

$$
\mathcal V_1\cap f_{\mathsf E}^{-1}(\mathcal V_1),
$$

decompose it into irreducible real algebraic components, where irreducibility is taken in the real Zariski topology, and then prove the typed equality

$$
\left\{\rho(\mathbf r):
\mathbf r\in\bigl(\mathcal V_1\cap f_{\mathsf E}^{-1}(\mathcal V_1)\bigr)\cap\mathbb B\right\}
=\mathcal C_1\cap\mathsf E^{-1}(\mathcal C_1).
$$

Repeat the global-preservation classification for $\mathcal C_2$ and $\mathcal C_3$.

**Solvability note.** Use P1’s polynomial equation as an identity on its zero set and substitute each map from D10. No master equation or external channel classification is required.

**Uses:** D1, D6, D7, D10, P1.

*Reading.* P4 finds invariant statistical models under noise. If a channel preserves a chip, the evolution of a qubit on that chip reduces to updates of two classical biases while remaining subject to the qubit physicality domain.

### D11. Rotated and relabeled chips

For $R\in SO(3)$ and a permutation $\pi$ of $\{1,2,3,4\}$, define

$$
Q_i^{R,\pi}:=\frac14\left(\mathbb I+(R\mathbf n_{\pi(i)})\boldsymbol\cdot\boldsymbol\sigma\right),
$$

$$
T_{R,\pi}(\rho):=\bigl(\operatorname{tr}(\rho Q_1^{R,\pi}),\ldots,
\operatorname{tr}(\rho Q_4^{R,\pi})\bigr),
$$

and

$$
\mathcal C(R,\pi):=\{\rho\in\mathcal S:\delta_1(T_{R,\pi}(\rho))=0\}.
$$

**Lemma 11a (certification).** Every $\mathcal C(R,\pi)$ is nonempty. Moreover,

$$
T_{I,\mathrm{id}}=T,\qquad \mathcal C(I,\mathrm{id})=\mathcal C_1.
$$

Changing $R$ transports the construction by the corresponding Bloch-ball rotation; changing $\pi$ selects one of the three pairing classes from D5.

**Lemma 11b (certification).** For every $\rho\in\mathcal S$, $T_{R,\pi}(\rho)\in\Delta_3$, and all row and column marginals associated with $\delta_1$ are strictly positive. Hence $\phi(T_{R,\pi}(\rho))$ is defined for every allowed $(R,\pi,\rho)$.

**Uses:** D1, D2, D3, D5, D7, D9.

*Reading.* There is no basis-free distinguished potato chip. The invariant idea is a family: choose a tetrahedral SIC and choose how its four outcomes are read as two bits; the corresponding zero-determinant patch is the chip.

### P5. Audit what survives the convention choice

For this problem, define the frame-dependent correlation value

$$
\phi_{R,\pi}(\rho):=\phi(T_{R,\pi}(\rho)).
$$

Give necessary and sufficient conditions on $(R,\pi)$ and $(R',\pi')$ for

$$
\mathcal C(R,\pi)=\mathcal C(R',\pi').
$$

Then identify, with proofs, which of the following are invariant across the full D11 family:

1. dimension and topology;
2. boundary purity;
3. existence of a two-setting promised-tomography protocol;
4. the physical axes measured in that protocol;
5. the numerical value and sign of $\phi_{R,\pi}(\rho)$ for a fixed state $\rho$;
6. preservation by each fixed laboratory-frame channel in D10.

**Solvability note.** The action of $SO(3)$ is explicit in D11, and the three pairing polynomials are explicit in D5.

**Uses:** D5, D8, D9, D10, D11, P1.

*Reading.* P5 prevents a coordinate artifact from being mistaken for an intrinsic class of qubit states. Geometry survives rotations, but axes, signs, pairings, and channel compatibility generally carry the chosen convention.

### D12. A pinned Wootters-type quasi-probability table

For $\rho(x,y,z)\in\mathcal S$, define the row-major vector

$$
W(\rho):=\frac14\bigl(
1+x+y+z,\ 1+x-y-z,\ 1-x-y+z,\ 1-x+y-z
\bigr).
$$

**Lemma 12a (certification).** The four entries of $W(\rho)$ are real and sum to one, but need not be nonnegative. In the displayed row-major ordering, its row marginals are

$$
\left(\frac{1+x}{2},\frac{1-x}{2}\right),
$$

and its column marginals are

$$
\left(\frac{1+z}{2},\frac{1-z}{2}\right).
$$

**Uses:** D1.

*Reading.* This rival representation makes the two desired marginals sharp without the SIC’s stochastic unsharping. Its price is negativity: it is not an ordinary probability table for every qubit state.

### P6. Compare SIC and quasi-probability factorization

Determine in closed form the set

$$
\mathcal C_W:=\{\rho\in\mathcal S:
\text{every entry of }W(\rho)\text{ is nonnegative and }\delta_1(W(\rho))=0\}.
$$

Give its Bloch equation, parameter domain, boundary relative to the locus $\delta_1(W(\rho))=0$, the purity classification of that relative boundary, and its relation to $\mathcal C_1$. State exactly which two-setting reconstruction steps become simpler and which positivity restrictions become stronger or weaker when $T$ is replaced by $W$.

**Solvability note.** D12 is an explicit row-major $2\times2$ table; its determinant, marginals, and entrywise positivity are finite polynomial inequalities. The maximally mixed state certifies that the quantified set is nonempty.

**Uses:** D1, D5, D7, D8, D12.

*Reading.* P6 separates two ideas that the source sometimes blends: factorization of a four-entry table, and positivity of the representation chosen for that table.

---

## Meta-problems

These duties require judgment, experimental context, or literature search and therefore sit outside the formal chain.

### M1. What is the chip good for?

Assess the practical value of the following uses without assuming an advantage in advance:

1. promised two-setting tomography and model validation;
2. low-dimensional calibration models for devices known to prepare chip states;
3. tracking invariant families under selected noise channels;
4. embedding physically admissible independent two-bit distributions into qubit SIC probabilities;
5. using $\phi$ as a goodness-of-fit statistic for the independence constraint.

For each use, specify the promise, resource count, estimator, failure mode, and comparison class. In particular, separate fewer measurement **settings** from fewer samples, lower variance, or lower total experimental cost.

### M2. Where has this structure appeared before?

Search the literature for the intersection of:

- tetrahedral qubit SIC probability representations;
- rank-one $2\times2$ probability tables and independence varieties;
- discrete Wigner or Wootters phase-space marginals;
- tomography under algebraic or manifold promises;
- invariant statistical models under quantum channels.

The deliverable is a claim-by-claim priority map, not a list of nearby papers: identify which exact definition, equation, reconstruction result, or channel-invariance result is already present and which combination, if any, is new.

### M3. What does the construction not show?

Evaluate and, where necessary, refute the following extrapolations:

1. that a one-qubit chip state is a bipartite quantum product state;
2. that $\phi=0$ proves absence of quantum correlations;
3. that two projective settings determine an arbitrary qubit;
4. that every distribution of two binary variables is represented by a physical chip state;
5. that the construction by itself yields a classical simulation of universal quantum computation;
6. that a time-local equation with a negative rate is automatically a physical Lindblad evolution.

---

## Minimality record

### Deletion test

| Item | Why it survives deletion |
|---|---|
| D1 | Supplies the state type, Bloch coordinates, and physical unit ball. Without it, “physical qubit state” is undefined. |
| D2 | Pins the SIC orientation and outcome order. Without it, the chip is convention-dependent but unnamed. |
| D3 | Connects qubit states to four ordinary probabilities. Without it, the product constraint cannot act on a state. |
| D4 | Distinguishes the physical probability body from the full simplex and supplies the inverse used by P1 and P2. |
| D5 | Makes the three inequivalent two-bit pairings explicit. Without it, multiplicity is hidden in a relabeling convention. |
| D6 | Defines classical independence and its two bias parameters. Without it, the word “factorization” has no operational content. |
| D7 | Names the object under study and quarantines the no-subsystem interpretation. |
| D8 | Connects product marginals to experimentally available projective settings. Without it, the tomography claim is not a protocol. |
| D9 | Turns the determinant constraint into an estimable correlation diagnostic with a pinned domain. |
| D10 | Pins the channel conventions needed for a checkable invariance classification. |
| D11 | Exposes rotation and relabeling dependence instead of presenting one chip as intrinsic. |
| D12 | Pins the quasi-probability rival used in the source and separates sharp marginals from positivity. |

### Collapse test

- D1 and D2 do not collapse: states and measurements have different types.
- D3 and D4 do not collapse: one is a map; the other is its physically allowed image and inverse characterization.
- D5 and D6 do not collapse: D5 supplies all pairing predicates; D6 supplies the specific product parametrization consumed by tomography.
- D7 and D11 do not collapse: D7 pins the reference witness; D11 audits the full convention orbit.
- D8 and D9 do not collapse: marginals enable reconstruction, whereas $\phi$ tests lack of factorization.
- D10 does not collapse into P4 because channel formulas are data and preservation is the requested theorem.
- D12 remains separate from P6 because it supplies a rival representation, certified marginals, and normalization properties before P6 imposes factorization and positivity.

### Primitive consumption

| Primitive | First load-bearing use |
|---|---|
| $M_2(\mathbb C)$, $\succeq0$, trace | D1 state space |
| Pauli matrices | D1 Bloch form and D2 effects |
| Euclidean norm | Lemma 1a physicality |
| $\Delta_1$ | D6 product-map domain |
| $\Delta_3$ | D3 and D6 probability-vector types |
| $SO(3)$ | D2’s rotated witnesses and D11 |
| Qubit channel | Lemma 10a and P4 |

### Per-problem increments

| Problem | Exact increment beyond the preceding problems |
|---|---|
| P1 | Adds only physicality D4, pairing D5, product parameters D6, and the named chips D7. |
| P2 | Adds only the two observable marginals and their fixed stochastic conversion D8. |
| P3 | Adds only the normalized diagnostic D9. |
| P4 | Adds only the six pinned channel families D10. |
| P5 | Adds only the explicit convention family D11. |
| P6 | Adds only the pinned quasi-table D12 and the explicit entrywise-nonnegativity restriction in its own set definition. |

---

## Audit record

### Round 0: source reconstruction and exact checks

**Method.** The TeX source `main.tex` inside `quantum_potato_chip.zip` was read in full. The SIC coordinate map, its inverse, the three pairing determinants, the product parametrization, the physicality inequality, boundary discriminant, worked rational example, and channel-preservation residuals were independently recomputed with exact symbolic algebra.

**Verified source backbone.** The tetrahedral SIC gives an informationally complete four-probability representation; a chosen $2\times2$ arrangement factors exactly when its determinant vanishes; intersecting that product surface with the qubit probability body defines the chip; and a chip-membership promise makes a two-setting reconstruction problem well formed. Exact channel residuals were computed, but their classification is deliberately left to P4.

**Findings and repairs.** Findings are retained here even though the chain above uses corrected statements.

1. **FATAL conceptual overclaim — repaired.** “Product state” and “no quantum correlations” were used for a one-qubit four-outcome table. D6, Lemma 7b, D9, and M3 replace this with classical independence of two outcome labels and explicitly deny a subsystem claim.
2. **FATAL scope overclaim — repaired.** The source suggested that any classical problem with two binary variables can be mapped to qubits. D4 and P1 require the independent distribution to lie in the qubit physicality body; M3 preserves the rejected extrapolation for explicit audit.
3. **MISLEADING tomography claim — repaired.** Two Pauli settings do not determine a general qubit. P2 makes chip membership a promise and asks for the non-unique unrestricted fiber.
4. **IMPRECISE geometric language — repaired.** Two “disjoint simplices” were described where the actual construction is the Cartesian product of two binary distributions. D6 states the product map.
5. **FATAL label mismatch in the source prose — repaired.** For the ordering in D6, the $z$ marginal is $(p,1-p)$ and the $x$ marginal is $(q,1-q)$. One paragraph in the TeX source swaps $p$ and $q$, while its own worked example uses the ordering recorded in Lemma 8b.
6. **MISLEADING measurement language — repaired.** The coarse-grained SIC marginals are unsharp binary POVMs, not projective measurements. D8 distinguishes them from sharp Pauli POVMs and pins the stochastic conversion.
7. **IMPRECISE uniqueness claim — repaired.** There are three chips only after fixing a SIC and quotienting outcome relabelings by row swaps, column swaps, and interchange of the two bits. D5 and D11 expose both choices.
8. **UNCERTIFIED dynamics claim — quarantined.** The TeX source’s time-local boundary evolution uses negative rates and a singular parameterization. It was not promoted to a lemma; M3 requires a separate physicality audit before calling it Lindblad evolution.

**Computation fidelity.** Exact algebra was performed with SymPy because no project-specific verification command was needed for these polynomial identities. No floating-point evidence was used for the certified clauses.

**Open issues after round 0:** 2 independent audits pending; schematic formal-language compilation not attempted.

### Round 1: fresh-context derivability and minimality audit

**Input isolation.** The auditor received this document alone, with no source TeX, other drafts, or web access.

**Checks completed.** Every certified lemma was independently recomputed and survived. P1 was solved privately in all six requested parts from D1–D7 alone; no undeclared input was needed. All dependency references pointed backward and were sufficient.

**Findings and repairs.**

1. **DISCIPLINE — repaired.** R1–R7 were referenced by the audit request but absent from the artifact. The Question-chain rules section now declares them locally and specifies direct-dependency semantics for **Uses** lines.
2. **IMPRECISE — repaired.** “Boundary” was ambiguous because a surface patch has empty interior in $\mathbb R^3$. P1 and P6 now ask for the boundary relative to the factorization surface or parameter domain.
3. **IMPRECISE — repaired.** P4’s “maximal algebraic subset” lacked a base field and order. P4 now names the exact real-algebraic preimage intersection, requests inclusion-maximal components, and only then imposes Bloch-ball physicality.
4. **IMPRECISE — repaired.** P5 used a frame-dependent correlation value without naming it. P5 now defines $\phi_{R,\pi}(\rho)$ explicitly.
5. **DISCIPLINE — repaired.** D13 was a one-use wrapper that failed the collapse test. It was removed, and its entrywise-positivity clause was inlined into P6.

**Machine-readable findings before repair:** {"open_issues":5,"fatal":0,"misleading":0,"imprecise":3,"discipline":2}.

**Open issues after repair:** 1 fresh-context re-audit pending; schematic Lean sketch intentionally remains noncompiling and is labeled as such.

### Round 2: fresh-context re-audit

**Input isolation.** A new auditor received the repaired document alone, without Round 1’s report, the source TeX, other drafts, or web access.

**Checks completed.** Every certified and substrate lemma again survived independent exact algebra or channel-level verification. P1 was again derived privately in all six parts from D1–D7 alone. The deletion/collapse record otherwise survived.

**Findings and repairs.**

1. **FATAL specification type error — repaired.** D10 used the same symbols as Bloch-vector maps and qubit channels, while P4 applied them to both types. D10 now distinguishes each affine coordinate map $f_{\mathsf E}:\mathbb R^3\to\mathbb R^3$ from its complex-linear channel $\mathsf E:M_2(\mathbb C)\to M_2(\mathbb C)$, and P4 uses each at the correct type.
2. **IMPRECISE — repaired.** P4’s algebraic “components” still lacked an irreducibility condition. It now requests irreducible components in the real Zariski topology.
3. **DISCIPLINE — repaired.** P4 directly requested D6’s product-parameter updates but omitted D6 from its **Uses** line.
4. **DISCIPLINE — repaired.** D9 directly quantified over D1’s state space but omitted D1 from its **Uses** line.
5. **DISCIPLINE — repaired.** D7’s interpretation invoked independence without consuming D6. Lemma 7c now pins all three permuted $2\times2$ tables and connects them formally to D6; D7’s **Uses** line includes D6.

**Machine-readable findings before repair:** {"open_issues":5,"fatal":1,"misleading":0,"imprecise":1,"discipline":3}.

**Open issues after repair:** 1 final fresh-context convergence audit pending; schematic Lean sketch intentionally remains noncompiling and is labeled as such.

### Round 3: final fresh-context convergence audit

**Input isolation.** A third new auditor received the twice-repaired document alone, without prior reports, the source TeX, other drafts, or web access.

**Checks completed.** Every certified and substrate lemma survived independent exact symbolic checks, including the SIC identities, probability-body inverse, factorization, observable marginals, channel actions and CPTP realizations, convention covariance, and quasi-probability table. P1 was privately derived in all six parts from D1–D7 alone. Type refinements and deletion/collapse tests passed.

**Finding and repair.**

1. **IMPRECISE — repaired.** P3’s normal displacement did not fix scale or orientation. P3 now defines the unit normal with negative $y$-component, pins the path $\mathbf r(\varepsilon)$, and asks for the coefficient of the declared parameter. Its **Uses** line now explicitly includes P1.

**Machine-readable findings before repair:** {"open_issues":1,"fatal":0,"misleading":0,"imprecise":1,"discipline":0}.

**Open issues after local repair:** 0 known. The repair was checked locally against R1–R7 and for notation/type consistency. In accordance with the three-round cap, it was not sent to a fourth fresh auditor; therefore this record does not call the final state an independently certified clean pass.

---

## Optional Lean-flavored sketch

The following is schematic and is **not claimed to compile**.

```lean
structure BlochState where
  x y z : ℝ
  physical : x^2 + y^2 + z^2 ≤ 1

def sicTable (r : BlochState) : Fin 4 → ℝ :=
  -- the four coordinates in Lemma 3b
  sorry

def pairingDet (k : Fin 3) (v : Fin 4 → ℝ) : ℝ :=
  -- the three polynomials in D5
  sorry

def chip (k : Fin 3) : Set BlochState :=
  {r | pairingDet k (sicTable r) = 0}

theorem reference_chip_closed_form :
  -- P1: derive the implicit equation and the exact parameter domain
  True := by
  sorry
```
