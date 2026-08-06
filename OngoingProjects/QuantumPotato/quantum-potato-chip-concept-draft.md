# Two independent coins inside a qubit

*Concept draft for the rewritten core of the quantum-potato-chip paper. Revision 6: the backbone is six standalone questions, each short, concrete, and executable with no context; every section answers exactly one. All mathematics verified in exact arithmetic (author plus two independent fresh-context reviewers in earlier rounds; the revision-4 wall fact machine-checked); literature grounded in the TeX under `related-papers/`. The appendix frame records what remains open.*

---

## The backbone

> **Q1.** For which qubit states do the four outcome probabilities of the SIC measurement factor into a product of two binary probability distributions, and what geometric object do those states form?
>
> **Q2.** For which qubit states does the $2\times2$ outcome table of a joint unsharp measurement of $\sigma_z$ and $\sigma_x$ factor into the product of its marginals, and which states do so for every such measurement?
>
> **Q3.** Which qubit channels map the set of states with factorizing SIC statistics into itself?
>
> **Q4.** How does the set of states with factorizing SIC statistics relate to the sets of states with nonnegative discrete Wigner or Kirkwood-Dirac distributions?
>
> **Q5.** What is a state with factorizing SIC statistics good for?
>
> **Q6.** Where, in any field's vocabulary, has the factorization of a quantum state's informationally complete statistics been studied?

Each question can be handed to a physicist with nothing else and acted on. Q1 alone generates the object: the SIC is named inside it, so the derivation below involves no choice the question has not already made. Q2 audits that choice against the full joint-measurement family; Q3 adds dynamics; Q4 locates the object against the standard classicality notions; Q5 and Q6 assess. Why these six are minimal, and why the natural-language question behind them ("is any qubit state exactly a pair of independent classical coins?") compiles to Q1, is recorded in the appendix; the one-line reason the SIC appears in Q1 is that a pair of coins has four outcomes, so the completeness witness must be a four-outcome (hence minimal) IC measurement, and the SIC is the unique symmetric sharp one.

## Q1. The chip

**Step 1.** The SIC effects are $\tfrac14(\mathbb{I}+\vec n_{ij}\cdot\vec\sigma)$ with the four unit vectors a regular tetrahedron; its outcome distribution is informationally complete, so the state *is* its table. Taking coordinates so the table's marginal directions are $z$ and $x$ (the tetrahedron fixes them to be orthogonal, and the marginals to be unbiased, automatically), the table of $\rho = \tfrac12(\mathbb{I}+\vec r\cdot\vec\sigma)$, $\vec r=(x,y,z)$, is

$$p_{ij} = \tfrac14\!\left(1 + \tfrac{1}{\sqrt3}\left[(-1)^i z + (-1)^j x + (-1)^{i+j} y\right]\right),$$

affine in $\vec r$, one Pauli component per character slot: for the outcome bits $A=(-1)^i$, $B=(-1)^j$: $\langle A\rangle = z/\sqrt3$, $\langle B\rangle = x/\sqrt3$, $\langle AB\rangle = y/\sqrt3$. The marginals are unsharp $Z$ and $X$ coins, $p=\tfrac12(1+z/\sqrt3)$, $q=\tfrac12(1+x/\sqrt3)$; sharp Pauli data recover them by the fixed rescaling $P\mapsto\tfrac{1}{\sqrt3}(P-\tfrac12)+\tfrac12$.

**Step 2 (elementary $2\times2$ identity).** Product-of-marginals $=$ zero covariance $=$ vanishing determinant, since for $2\times2$ tables the covariance *is* the determinant.

**Step 3 (the equation).** $4\det T = \langle AB\rangle-\langle A\rangle\langle B\rangle$ gives in one line

$$\det T \;=\; \frac{\sqrt3\,\langle\sigma_y\rangle - \langle\sigma_x\rangle\langle\sigma_z\rangle}{12},$$

so the factorizing states of this pairing satisfy $\sqrt3\,y = xz$: a hyperbolic paraboloid, doubly ruled, its rulings the constant-$p$ and constant-$q$ lines. The saddle itself is classical statistics: the independence surface of the $2\times2$ contingency table (the Segre quadric), drawn by Fienberg and Gilbert in 1970. Nothing quantum has happened yet.

**Step 4 (the physical cut).** Positivity of $\rho$ is $|\vec r\,|\le1$, which in table space is exactly the insphere of the probability simplex. Both inclusions are immediate: products have $\det T=0$, and every physical zero of $\det T$ has coins in $(0,1)$ whose product reproduces the table.

**Step 5 (the three pairings).** Four outcomes admit exactly three $2\times2$ arrangements, so "factor into a product of two binary distributions" is the union over the three pairings: the loci $\sqrt3y=xz$, $\sqrt3z=xy$, $\sqrt3x=yz$, cut by the ball. **The answer to Q1**: three congruent saddle-shaped surfaces (the *chips*), each containing two Cartesian axes as rulings and excluding the third, meeting pairwise on the shared axes and jointly at the maximally mixed state.

**Structure (all verified exactly).** The boundary of a chip, its pure states, is the quartic curve through the four eigenstates of its two coin observables, where one coin reaches its extremal allowed bias $\tfrac12\pm\tfrac{1}{2\sqrt3}$ while staying unsharp; coin parametrization $\vec r(p,q)=\sqrt3\,(2q-1,\ (2p-1)(2q-1),\ 2p-1)$. The unitary covariance group of one chip is the Klein group of $\pi$-rotations about the coordinate axes; the 24 outcome relabelings split $8+8+8$ over the three chips (12 unitary), with CNOT-type relabelings permuting them. With $\phi$ the Matthews (Yule) coefficient of the table, a chip is exactly $\{\phi=0\}$; the nearby bounds $|\phi|\le1/\sqrt3$ (equality only at the excluded-axis poles) and $|2p-1|\le1/\sqrt3$ hold for *all* states, so they belong to the state space and the sharpness, not to chip membership, and are claimed for neither. The $\epsilon$-version is the odds-ratio foliation: level sets of $p_{00}p_{11}/(p_{01}p_{10})$ are again quadrics cut by the ball, the chip the unit level. Rearranged, Step 3 reads $\sqrt3\,\langle\sigma_y\rangle = 12\det T + \langle\sigma_x\rangle\langle\sigma_z\rangle$: the Bloch component invisible to the two coins is an affine readout of their covariance.

## Q2. The full family: what is convention, what survives

The joint unsharp measurements of $\sigma_z$ and $\sigma_x$ with fixed marginal sharpnesses $\eta_1, \eta_2$ form, by the joint-measurability literature (Yu-Liu-Li-Oh), a family with a parity scalar and a parity vector beyond the marginals. Organize the answer by symmetry:

- *Pauli-covariant members* (conjugation by $X$ and $Z$ flips the corresponding outcome label; provably equivalent to parity scalar $=0$ and parity vector along $\hat y$): $M_{ij}=\tfrac14(\mathbb{I}+(-1)^i\eta_1\sigma_z+(-1)^j\eta_2\sigma_x+(-1)^{i+j}\gamma\sigma_y)$, $\eta_1^2+\eta_2^2+\gamma^2\le1$. The identity becomes $4\det T=\gamma y-\eta_1\eta_2 xz$; the locus is $y=kxz$ with $k=\eta_1\eta_2/\gamma$ covering all of $\mathbb{R}\setminus\{0\}$ (content: $|k|$). The factorizing-for-every-member states are $\bigcap_k\{y=kxz\}$: **the two coin axes**, the states diagonal in one coin basis.
- *Symmetric members* (no distinguished outcome): depolarized SICs, $k=s$ with effect eigenvalues $(1\pm\sqrt3 s)/4$, so positivity walls the family at $s\le1/\sqrt3$: **the SIC chip is the most curved chip any genuine symmetric measurement allows**, which is the precise sense in which naming the SIC in Q1 was canonical rather than arbitrary.
- *Non-covariant members:* nothing survives. A parity-biased joint measurement of the same coins has $\det T = Z_{\rm parity}/4$ at the maximally mixed state, so even the fair-coin point leaves the locus, and the loci of all members intersect in the empty set. Independence inside one system is measurement-relative all the way down (contrast two qubits, where "local tables factorize for every choice of local IC measurements" picks out the product states declaration-free).
- *Beyond positivity:* the algebra continues into quasi-probability frames ($\eta_1^2+\eta_2^2+\gamma^2>1$; not measurements). The Wootters phase-point frame is $\eta_1=\eta_2=\gamma=1$: its locus $y=xz$ is the factorization locus of the qubit's discrete Wigner function, beyond the positivity wall, which is why it is a frame and not a measurement. The complex Kirkwood-Dirac table of sharp $Z$ and $X$ has $\det Q_{\rm KD}=-(xz+iy)/4$: zero set again the two coin axes.

## Q3. The invariant dynamics

Bit flip, phase flip, and phase damping map each chip into itself, acting as classical noise on a single coin ($p\mapsto p+\xi(1-2p)$; $q\mapsto q+\xi(1-2q)$; $q\mapsto\tfrac12(1-\sqrt{1-\xi}+2q\sqrt{1-\xi})$); bit-phase flip, depolarizing, and amplitude damping generically break factorization. The preservers found so far are exactly classical noises *of the coins themselves*. Open, and posed as the paper's forward question: is "acts as a product of classical channels on the coins" the exact characterization of the chip-preserving channels?

## Q4. Against the positivity notions of classicality

Transverse, not nested, verified in both directions: there is a pure SIC-chip state whose Wootters function is negative (minimum phase-point entry $\approx-0.024$), and a physical Wootters-locus state outside the frame-independent classical octahedron $C_2$ of Galvão and Cormick et al. ($|x|+|y|+|z| = 2\sqrt{\sqrt2-1}+\sqrt2-1\approx1.70$); the one containment that holds (Wootters locus inside Wootters-nonnegative) is a tautology. Against Kirkwood-Dirac positivity the relation is structural: the KD-positive qubit polytope is characterized by an *additive* condition on the same $2\times2$ table, the chip by the *multiplicative* one. Classicality of description (some nonnegative quasi-distribution exists) and classicality of correlation (the jointly measured complementary marginals are independent) are different axes, and the chip is the locus that sets the second to zero.

## Q5. What it is good for

- **A tomography mechanism.** Under the promise of chip membership, two projective measurements do not merely coordinate the state: their rescaled outer product *is* its complete statistics. (Determination under any two-parameter promise is generic; the mechanism is chip-specific.)
- **A foundations statement.** Q4's transversality: correlation-classicality is an axis the positivity programs do not contain.
- **A teaching representation.** The coin table gives the Bloch ball's $y$ direction a classical meaning (the covariance of the two coins), and the chip is the cleanest instance of "state $=$ classical data $+$ a promise".
- **A template.** The same question at $d\ge3$ is genuinely open: a qutrit SIC 9-vector as a $3\times3$ table factorizes iff all $2\times2$ minors vanish (the determinantal Segre variety, dimension $2(d-1)$), met with the qudit body; whether the cut reaches pure states and what replaces the ruled saddle are computable unknowns.
- **Not claimed.** No computational advantage (chip states are two coin flips, the maximally simulable corner of the qubit); no bound as a chip theorem; no measurement-free surface (Q2 is the honest version).

## Q6. Prior art, in every vocabulary found

Grounded in the fetched TeX under `related-papers/` (five sources read in full; three large ones at theorem level with exhaustive term sweeps) plus two reviewers' independent searches:

- **Algebraic statistics** (Fienberg-Gilbert 1970; Slavković-Fienberg; the independence model as Segre variety): owns the saddle and the $\phi$ level sets in the bare simplex; no quantum region to cut with.
- **Joint measurability** (Busch; Busch-Heinonen; Yu-Liu-Li-Oh, arXiv:0805.1538): owns Q2's measurement families including the parity freedom; asks when joint measurements exist, never which states factorize a given table.
- **Discrete-Wigner classicality** (Galvão, quant-ph/0405070; Cormick et al., quant-ph/0506222; read in full): nonnegativity polytopes, stabilizer states, Clifford covariance; the product predicate never posed.
- **Kirkwood-Dirac positivity** (Langrenez et al., 2306.00086; review 2403.18899): the additive-condition polytope; independence absent; the KD determinant connects their object to Q2's core.
- **Probability representations** (Chernega-Man'ko-Man'ko, 1803.09339, 1902.03613; read in full): a qubit as three *separate* coins, explicitly no joint distribution; the factorization question is undefinable there.
- **Spekkens' toy theory** (quant-ph/0401052; consolidation 2105.03277): the combinatorial skeleton (which pure states are coin-like; where correlation is extremal) is answered there *under the epistemic restriction taken as an axiom*; the quantitative content (the quadric with mixed interior, the wall, the Born-rule covariance identity, the transversality) is not.
- **SIC state-space geometry** (Appleby-Ericsson-Fuchs, 0910.2750, read in full; Qplex, 1612.03234): "which probability vectors are quantum"; Step 4 uses their insphere; this work asks the reverse question and composes with theirs.
- **Direct searches** in SIC/Wigner/KD vocabulary for the question or the answer: nothing; arXiv:2411.01082 has no citing follow-ups on the concept.

**What is new, sized honestly.** The question is unasked in quantum vocabularies (confirmed independently twice), with the gloss that the answering identity is elementary once the witness is parametrized; the saddle is classical statistics; the skeleton is the toy bit's under its axiom. Remaining new: the physical cut in closed form; the Q2 audit with its positivity wall at the SIC, its two-axes core, and its emptiness theorem without covariance; the transversality counterexamples; the invariance structure with its open characterization; the $d\ge3$ program inside the qudit body. The novelty move is composition across community boundaries.

---

---

## Appendix: Draft brief and inquiry frame (working scaffolding, drop from any journal version)

```
DRAFT BRIEF
Goal:           Core-concept text whose backbone is six standalone questions, each executable with
                no context, with every section answering exactly one.
Reader:         Expert physicists (foundations / quantum information referees) and co-authors.
Claim:          Q1 alone generates the chip (the SIC is named in the question): sqrt3 y = xz and its
                two relabelings, cut by the Bloch ball. Q2 shows the choice was canonical (positivity
                walls symmetric measurements at the SIC) and finds the measurement-free core (the
                two coin axes; empty once covariance is dropped). Q3: preservers found = classical
                coin noises; exactness open. Q4: independence and positivity are transverse.
Genre/length:   Research-note markdown, seed of the paper; math in LaTeX.
Scope:          in: Q1-Q6 as sectioned. out: Liouvillian section; d>=3 computation (posed open);
                epsilon-version beyond the foliation remark.
Structure:      backbone -> Q1 -> Q2 -> Q3 -> Q4 -> Q5 -> Q6 -> frame. Build order = numbering.
Evidence plan:  All equations verified exactly (author + two fresh-context reviewers); wall fact
                hand- and machine-verified; literature via fetched TeX + reviewers' searches.
Voice:          Direct physicist prose; verbs before nouns; no em dashes; elementary steps labeled
                elementary; every convention declared where used.
Skill invoked:  adversarial-draft (research mode).
Verify first:   [cleared; scopes noted] Q1 Steps 1-5 (SIC identity, insphere both ways, both
                inclusions, three-pairing union); Q2 covariant identity, k onto R\{0}, covariance
                <=> parity-freedom vanishing [hand-proven], parity det = Z/4 at center
                [hand-verified], wall [hand- and machine-verified: eigenvalues (1 +- sqrt3 s)/4;
                s = 1/sqrt3 exactly rank-one; Wootters min eigenvalue -0.183; not yet reviewed],
                KD determinant [hand-verified]; Q3 channel reparametrizations; Q4 counterexamples
                [independently spot-verified]; phi and bias bounds [verified; attributed to state
                space/sharpness]; 8+8+8 split; Klein group.
Constraints:    No em dash; delimited TeX; no bounds-as-chip-theorems; no "refines positivity";
                Wootters labeled frame, not measurement; signs consistent with p = (1+z/sqrt3)/2.
Success test:   (i) Each backbone question is actionable by a physicist given nothing else (the
                author's test: "if I give you these questions with no context, you can do something
                with them"); (ii) Q1 alone yields the chip; (iii) no section is orphaned from the
                backbone; (iv) a hostile referee can neither refute a stated claim nor source what
                Q6 claims as new.
Objections preempted:
                (1) "why the SIC in Q1" -> one line in the backbone (four outcomes force minimal IC;
                    SIC = unique symmetric sharp witness) + Q2's wall makes it canonical;
                (2) "independence is measurement-relative" -> Q2, including the two-qubit contrast;
                (3) "the saddle is classical statistics" -> attributed at Q1 Step 3;
                (4) "the identity is elementary" -> labeled at Q1 Steps 2-3;
                (5) "bounds are not chip theorems" -> stated inside Q1 structure;
                (6) "Spekkens has the skeleton" -> answered-under-rejected-assumptions in Q6;
                (7) "who acts differently" -> Q5's uses + Q3's open question.

INQUIRY FRAME
Backbone:       Q1 [generative; survived deletion trivially; the SIC named inside it, so the
                    derivability test passes from Q1 alone]
                Q2 [survived deletion: without it Q1's choice is convention presented as physics;
                    survived collapse into Q1: different family, different answer type]
                Q3 [survived deletion: standalone dynamics; the forward question]
                Q4 [survived deletion: locates the object against both positivity programs;
                    concrete set relations, not rhetoric]
                Q5, Q6 [assessment; survived deletion as the worth and location duties]
                Self-containment test (the trigger of this revision): each question names its
                objects (SIC, joint unsharp Z-X measurement, channels, Wigner/KD positivity,
                IC statistics) with no reference to any other question or to "the answer".
                Derived, not listed: parametrization, boundary quartic, three-chip orbit, Klein
                group, odds-ratio foliation, reconstruction mechanism, y-as-covariance reading.
Missing-hunt:   Boundary: one qubit. Defended: the two-qubit contrast (inside Q2's answer) shows
                why the single-system question is the one with content to audit. Adjacent probes
                engaged: statistics (saddle), joint measurability (Q2 family incl. parity), toy
                theory (skeleton), CV (closed via Cohen-Zaparovanny + KLM; retained in session
                record). Contingency: outside covariance even the fair-coin point leaves the locus,
                so Q2's finding is load-bearing and stated, not hidden.
Assumptions:    Background QM + qubit SIC existence (constructed). No further assumptions: the
                witness is named in Q1 and audited in Q2, so the former axiom S has become a
                theorem-backed canonicity claim (the wall) instead of an input.
Prior work:     per backbone question: Q6 section.
Translations:   Quantum T1-T12 (round 1; T7 structurally blind, acknowledged); round 2 added:
                independence model / Segre / contingency geometry; joint unsharp qubit observables;
                epistricted toy bit; positive joint quantum distributions; Wigner-representability.
Originality:    Question: unasked in quantum vocabularies (confirmed twice); identity: elementary;
                saddle: classical statistics; skeleton: toy bit under its axiom. Remaining-new:
                physical cut in closed form; Q2 audit (wall, core, emptiness); transversality;
                invariance + open characterization; d>=3 in-body.
Novelty move:   Composition of known objects across community boundaries.
Criterion:      best backbone = fewest standalone questions, each executable cold, jointly
                generating and auditing the object; best derivation = both inclusions proven with
                the measurement dependence exposed.
Approach:       backbone-first with concrete questions. Alternatives against the criterion:
                rev-5's referential chain (failed the self-containment test: "that choice",
                "the answer", "found in Q1"; the trigger of this revision); rev-3 declaration-in-
                prose (failed derivability); parametrization-first (arXiv v1; one inclusion).
Sense check:    Q1 Step 4 reproduces the known d=2 insphere (0910.2750); Q2's k=1 frame reproduces
                the Wootters marginal structure; the KD determinant reproduces the two-axes core;
                the parity member reproduces the Yu-Liu-Li-Oh freedom; the wall is consistent with
                Wootters (k=1) being a non-measurement.
So what:        Q5's four uses + the not-claimed list.
Section map:    backbone -> all; sections Q1..Q6 one-to-one. No orphans.
Revision log:   (1)-(2) [r0] reconstruction demoted; question replaced (procedure -> property).
                (3) [r1->2] representation dependence promoted to a question.
                (4) [r1->2] "refines positivity" -> transversality (counterexamples both ways).
                (5) [r1->2] single chip -> k-family; sqrt3 demoted to convention.
                (6) [r1->2] saddle attributed to classical statistics; novelty re-sized.
                (7) [r1->2] bounds-as-significance added; WITHDRAWN in (10).
                (8) [r2->3] family re-scoped to the covariant slice; covariance <=> slice theorem;
                "dropping covariance: nothing survives" added.
                (9) [r2->3] Wootters reclassified as frame; k range R\{0}; signs fixed.
                (10) [r2->3] bounds withdrawn as chip theorems; operational core -> invariance.
                (11) [r2->3] identity labeled elementary; Spekkens upgraded; CV closed; foliation
                stated; clearance scopes recorded.
                (12) [r3->4] kernel restructuring {Q, S}; trigger: derivability test failed on
                rev-3's questions; rung-1 positivity wall surfaced and verified.
                (13) [r4->5] backbone-first restructuring; trigger: author directive.
                (14) [r5->6] BACKBONE REWRITTEN AS STANDALONE QUESTIONS; trigger: author's
                self-containment test failed on rev-5's set (referential chain: "that choice",
                "the answer", "found in Q1"). Each question now names its objects; the SIC moved
                inside Q1, making Q1 alone generative; the axiom S dissolved into Q2's canonicity
                finding (the wall); sections re-hung one-to-one. Content unchanged except
                connective prose.
OPEN (rounds capped at 2 before revisions 4-6):
                (a) Q2's covariant-core and emptiness statements rest on verified witnesses; short
                proof write-ups due.
                (b) Q3's characterization conjecture unproven.
                (c) d>=3: program only.
                (d) Revisions 4-6 have not been through a fresh adversarial round.
```
