# The chip, as a formal chain

*Sibling chain: [quantum-potato-chip-paper-question-chain.md](quantum-potato-chip-paper-question-chain.md) compiles the 2024 source paper (arXiv:2411.01082) in the paper's conventions; the present chain formalizes the two-coins question and pairs with [quantum-potato-chip-concept-draft.md](quantum-potato-chip-concept-draft.md) and [quantum-potato-chip-derivation-pra.md](quantum-potato-chip-derivation-pra.md).*

*One question drives this document: **is any qubit state exactly a pair of independent classical coins?** The document turns that question into mathematics the way a proof assistant such as Lean would demand it: objects are defined strictly in order, each built only from what is already on the page; facts are stated as lemmas (provable, though not proven here); open questions are stated as problems, posed only about objects already defined. Each item lists the earlier items it depends on ("Uses"), and when a definition quietly claims a fact (for example, that four explicitly written matrices really do form a measurement), that claim belongs to the lemmas attached to the item. Under every item, a* ***Reading*** *paragraph says in plain language what the object is, what each clause of the formula does, and why the item is there; Readings explain, but they never smuggle in assumptions and never reveal the answer to any problem. Lemmas carry one of two tags: **certification** (the lemma is why an object deserves its name) or **interpretation** (the lemma gives meaning; no computation needs it).*

*Checking: two independent reviewers examined earlier versions, each given this file and nothing else. Both solved Problem P1 using only what is written here; both re-derived every lemma by direct computation; the second reviewer found one false claim (an earlier form of Lemma 6a's uniqueness statement, corrected by pinning the constant $\tfrac{1}{12}$) and several places where a symbol was used before being defined, all since repaired. The present version has not been re-checked.*

---

## Primitives (the starting vocabulary, declared once)

$\mathbb{R}$, $\mathbb{C}$; $\mathbb{R}^3$ with Euclidean norm $|\cdot|$; the $2\times2$ complex matrices $M_2$ with adjoint $\dagger$, trace $\mathrm{tr}$, identity $\mathbb{I}$, matrix rank; **positive semidefinite** $\succeq 0$ (Hermitian with nonnegative spectrum); **unitary** ($U \in M_2$ with $U^\dagger U = \mathbb{I}$); the Pauli matrices

$$\sigma_x = \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix}, \qquad
\sigma_y = \begin{pmatrix} 0 & -i \\ i & 0 \end{pmatrix}, \qquad
\sigma_z = \begin{pmatrix} 1 & 0 \\ 0 & -1 \end{pmatrix};$$

finite probability distributions (nonnegative entries summing to $1$); finite index families; the outcome set $\Omega := \{0,1\}\times\{0,1\}$ (the sample space of two bits), elements $(i,j)$.

*Reading: the raw materials: the linear algebra of one qubit plus the sample space of two classical bits. $M_2$ is the qubit's operator space: every state and every measurement operator below is a $2\times2$ matrix; $\{\mathbb{I}, \sigma_x, \sigma_y, \sigma_z\}$ is a basis of its Hermitian part, so any Hermitian matrix is $a\mathbb{I} + \vec b\cdot\vec\sigma$ with four real numbers, which is why three real coordinates will soon carry an entire state. $\Omega$ is the four-point sample space of a pair of bits: every measurement below has its outcomes labeled by $\Omega$, and every probability or quasi-probability object below is a function on $\Omega$; which physical observable each bit tracks is deliberately fixed only at D6. Physics enters this chain exactly twice: through the Born rule (D3) and through completely positive dynamics (D9); everything else is bookkeeping over these primitives.*

## The chain

**D1 (State).**

$$S \;:=\; \{\, \rho \in M_2 \;\mid\; \rho \succeq 0,\ \ \mathrm{tr}\,\rho = 1 \,\}.$$

*Lemma 1a (Bloch form): the map $\rho \mapsto \big(\mathrm{tr}(\rho\sigma_x),\, \mathrm{tr}(\rho\sigma_y),\, \mathrm{tr}(\rho\sigma_z)\big)$ is a bijection from $S$ onto the closed unit ball of $\mathbb{R}^3$, with inverse $(x,y,z) \mapsto \tfrac12\big(\mathbb{I} + x\sigma_x + y\sigma_y + z\sigma_z\big)$.*
*Uses: none.*

*Reading: the states of a qubit: density matrices, pure or mixed. The two clauses are the two probability axioms in operator form: $\rho \succeq 0$ will make every outcome probability nonnegative, and $\mathrm{tr}\,\rho = 1$ will make probabilities total $1$, once the Born rule (D3) converts operators to numbers. Lemma 1a unpacks to the Bloch ball: expanding $\rho$ in the Hermitian basis gives $\rho = \tfrac12(\mathbb{I} + x\sigma_x + y\sigma_y + z\sigma_z)$, and positivity is exactly $x^2+y^2+z^2 \le 1$, with pure states on the sphere and the maximally mixed state at the center. The coordinates are themselves expectation values, $x = \langle\sigma_x\rangle$ and so on, so "a state is three numbers" is literal; every problem below is solved in these coordinates.*

**D2 (Table measurement).**

$$\mathrm{TM} \;:=\; \Big\{\, E : \Omega \to M_2 \;\Big|\; \forall p \in \Omega:\ E(p) \succeq 0, \ \ \text{and}\ \ \sum_{p\in \Omega} E(p) = \mathbb{I} \,\Big\}.$$

*Uses: none.*

*Reading: a four-outcome quantum measurement, one outcome per point of the sample space $\Omega$. $E$ assigns to each sample point $(i,j) \in \Omega$ its **effect** $E(i,j) \in M_2$, the operator whose Born pairing with the state (D3) is that outcome's probability; positivity of the effects makes every probability nonnegative, and the sum condition makes the four probabilities total $1$, so $\mathrm{TM}$ is exactly the four-outcome POVMs arranged as a $2\times2$ array. A projective measurement of a qubit has at most two outcomes, so any member with four nonzero effects is necessarily unsharp (its effects are not projectors); that unsharpness is what makes four informative outcomes possible on a two-dimensional system. Nothing else is assumed; the two bits are what will play the role of the two coins.*

**D3 (Table of a state).**
For $\rho \in S$ and $E \in \mathrm{TM}$:

$$T(\rho,E) : \Omega \to \mathbb{C}, \qquad T(\rho,E)(p) \;:=\; \mathrm{tr}\big(\rho\, E(p)\big).$$

*Lemma 3a: the values of $T(\rho,E)$ are real, lie in $[0,1]$, and sum to $1$: $T(\rho,E)$ is a probability distribution on $\Omega$.*
*Uses: D1, D2.*

*Reading: the Born rule, the first of the chain's two physics inputs: $T(\rho,E)(p)$ is the probability of outcome $p$ when the measurement $E$ is performed on the state $\rho$, so $T(\rho,E)$ is the joint distribution of the two outcome bits. Lemma 3a is where D1's and D2's clauses pay off in pairs: positivity of $\rho$ against positivity of each $E(p)$ yields nonnegative values, and $\mathrm{tr}\,\rho = 1$ against $\sum_p E(p) = \mathbb{I}$ yields total $1$. The structural fact hiding here, and doing all the work later: $T$ is linear in $\rho$, so in Bloch coordinates every table entry is an affine function of $(x,y,z)$, which is why every locus in this chain comes out as a low-degree algebraic surface.*

**D4 (Determinant and product).**
For any $v : \Omega \to \mathbb{C}$:

$$\det v \;:=\; v(0,0)\,v(1,1) \;-\; v(0,1)\,v(1,0).$$

A probability distribution $v$ on $\Omega$ is **product** iff $\det v = 0$.
*Lemma 4a: a distribution $v$ is product iff $v(i,j) = a(i)\,b(j)$, where $a(i) := v(i,0)+v(i,1)$ and $b(j) := v(0,j)+v(1,j)$ are its marginals. (Gloss, in fresh symbols: with $\langle f\rangle_v := \sum_{p\in \Omega} v(p)f(p)$ and $\alpha(i,j) := (-1)^i$, $\beta(i,j) := (-1)^j$: $\ 4\det v = \langle\alpha\beta\rangle_v - \langle\alpha\rangle_v\langle\beta\rangle_v$, the covariance of the two outcome bits.)*
*Uses: none.*

*Reading: statistical independence of the two bits, in three equivalent guises. $\det v$ treats the table as a $2\times2$ matrix; Lemma 4a says its vanishing is the same as the table being the outer product of its own two marginals, and the gloss says the same again probabilistically: $4\det v$ is the covariance of the $\pm1$-valued bits, so "product" $=$ "zero covariance" $=$ "zero determinant". The special luck of the $2\times2$ case, which the whole chain exploits: for two binary variables zero covariance already forces full independence, so one scalar equation carries the entire condition (for larger alphabets it would not). The determinant is deliberately defined on complex-valued tables: the same quantity will later be applied unchanged to a quasi-probability table (D11) where the entries are not real.*

**D5 (Re-indexing).**
For a bijection $\pi : \Omega \to \Omega$ and $E \in \mathrm{TM}$: $\ E\circ\pi \in \mathrm{TM}$.
*Uses: D2.*

*Reading: relabeling the outcomes: $E\circ\pi$ is the same apparatus with the same four effects, just with the sample points renamed. Relabeling is physically free, and it is the only source of multiplicity in this chain: the $24$ bijections of $\Omega$ regroup the four outcomes into exactly three inequivalent ways of splitting them into a row bit and a column bit. This is why whatever locus a problem produces will come in three copies, and it is the operation that Lemma 6a's uniqueness clause and P1's list of relabelings range over.*

**D6 (The reference measurement).**
$E^\star \in \mathrm{TM}$,

$$E^\star(i,j) \;:=\; \tfrac14\left( \mathbb{I} \;+\; \tfrac{1}{\sqrt3}\Big[\, (-1)^j\, \sigma_x \;+\; (-1)^{i+j}\, \sigma_y \;+\; (-1)^i\, \sigma_z \,\Big] \right).$$

*Lemma 6a (certification): $E^\star \in \mathrm{TM}$; each $E^\star(p)$ has rank one; $\mathrm{tr}\big(E^\star(p)\,E^\star(q)\big) = \tfrac{1}{12}$ for all $p \ne q$; and any $E \in \mathrm{TM}$ with rank-one effects and $\mathrm{tr}\big(E(p)E(q)\big) = \tfrac{1}{12}$ for all $p \ne q$ equals $U\, E^\star(\pi(\cdot))\, U^\dagger$ for some unitary $U$ and re-indexing $\pi$. (The pinned value matters: rank-one effects with some common pairwise overlap do not suffice; families with weight multiset $\{s,s,1-s,1-s\}$, $s \ne \tfrac12$, and common overlap $s(1-s)/3$ exist.)*
*Lemma 6b (interpretation; consumed by no problem's computation): $\rho \mapsto T(\rho,E^\star)$ is injective on $S$: the table is a complete description of the state.*
*Lemma 6c (marginals): the marginals of $T(\rho,E^\star)$ are*

$$a(i) = \tfrac12\big(1 + (-1)^i\, \eta\, \mathrm{tr}(\rho\sigma_z)\big), \qquad b(j) = \tfrac12\big(1 + (-1)^j\, \eta\, \mathrm{tr}(\rho\sigma_x)\big), \qquad \eta = \tfrac{1}{\sqrt3}.$$

*(Informal gloss: two unbiased coins tied to $\sigma_z$ and $\sigma_x$, of sharpness the coefficient $\eta$.)*
*Uses: D1, D2, D3, D4, D5.*

*Reading: the reference measurement of the whole chain: the qubit SIC-POVM, written out. Each effect is $\tfrac14(\mathbb{I} + \vec n\cdot\vec\sigma)$ with $\vec n$ one of the four unit vectors of a regular tetrahedron in the Bloch ball; its eigenvalues are $\{0, \tfrac12\}$, so each effect is half a pure-state projector, and the sign pattern $((-1)^j, (-1)^{i+j}, (-1)^i)/\sqrt3$ fixes both the tetrahedron and the grouping into two bits at once, so no silent choice remains. Lemma 6a certifies the name (a measurement, rank-one effects, all pairwise overlaps $\tfrac{1}{12}$) and records a warning earned during checking: the pinned value is essential, because equal overlaps at other values admit unequal-weight families that are not the SIC. Lemma 6b is informational completeness: a qubit has $4 = d^2$ real parameters counting normalization, so four is the minimum outcome count for a complete measurement, and for this one the four probabilities determine $\rho$: the table* is *the state. Lemma 6c identifies the two bits physically at last: the row bit is an unsharp $\sigma_z$ measurement and the column bit an unsharp $\sigma_x$ measurement, both unbiased with sharpness $1/\sqrt3$, and this orthogonal, unbiased coin pair is forced by the tetrahedron, not assumed: these are the two coins of the title.*

**D7 (Chip).**

$$\mathrm{Chip} \;:=\; \{\, \rho \in S \;\mid\; T(\rho,E^\star)\ \text{is product} \,\}.$$

*Uses: D1, D3, D4, D6.*

*Reading: the object of the paper. By Lemma 6b the table $T(\rho,E^\star)$ is the complete information content of $\rho$, and by D4 "product" says its two coins are statistically independent; so $\mathrm{Chip}$ is the set of qubit states that are, informationally, nothing but a pair of independent biased coins. A dimension count says what to expect before computing: states carry three parameters (Lemma 1a), coin pairs carry two, and D4 is one scalar equation, so if the set is nonempty it should generically be a surface in the Bloch ball; existence is therefore cheap, and the exact surface is the content. This named set is what P1 computes, P3 evolves, and P4 locates.*

> **P1.** Give $\mathrm{Chip}$ in closed form, and give the sets $\{\, \rho \in S \mid T(\rho,\, E^\star\!\circ\pi)\ \text{is product} \,\}$ for all bijections $\pi : \Omega \to \Omega$.
> *Uses: D5, D7 (hence D1-D6).*
> *Solvability note: P1 is a finite computation from this file alone: parametrize $S$ by Lemma 1a, evaluate D3 on D6, impose D4.*

*Reading: derive the chip. The deliverable is an explicit equation in Bloch coordinates together with the ball constraint, and then the second task: by D5's reading the relabelings produce the sibling surfaces of the other groupings, so the task is to exhibit all of them and count the distinct ones. The solvability note is not a hope but a record: both independent reviewers carried out exactly this three-step computation from the file alone and agreed; the answer is deliberately not printed anywhere in this document, so it remains a problem set.*

**D8 (Coin-pair family).**
For $\eta_1, \eta_2 \in (0,1]$:

$$J(\eta_1,\eta_2) \;:=\; \Big\{\, E \in \mathrm{TM} \;\Big|\; \forall i:\ E(i,0)+E(i,1) = \tfrac12\big(\mathbb{I} + (-1)^i \eta_1 \sigma_z\big), \ \ \forall j:\ E(0,j)+E(1,j) = \tfrac12\big(\mathbb{I} + (-1)^j \eta_2 \sigma_x\big) \,\Big\}.$$

*Lemma 8a (certification): $E^\star \in J\big(\tfrac{1}{\sqrt3}, \tfrac{1}{\sqrt3}\big)$.*
*Lemma 8b (nonemptiness): $J(\eta_1,\eta_2) \ne \emptyset$ iff $\eta_1^2 + \eta_2^2 \le 1$. (So $E^\star$ sits strictly inside the admissible disk: $\tfrac13 + \tfrac13 < 1$.)*
*Uses: D2, D6.*

*Reading: all joint measurements of the same two coins. The first constraint line fixes the row marginal to be the unsharp $\sigma_z$ coin of sharpness $\eta_1$, the second fixes the column marginal to the unsharp $\sigma_x$ coin of sharpness $\eta_2$; everything the marginals do not fix stays free, so $J(\eta_1,\eta_2)$ is the full set of four-outcome measurements that jointly realize this coin pair. Lemma 8b is the known trade-off for measuring two incompatible qubit observables together, in its disk form $\eta_1^2+\eta_2^2 \le 1$: sharpening one coin forces blurring the other, and beyond the disk no joint apparatus exists at all, which is why P2's range stops there. Lemma 8a places the SIC inside this family, strictly interior to the disk; the family is P2's arena, the space of conventions against which the chip's chosen reference measurement is examined.*

> **P2.** For $\eta_1^2 + \eta_2^2 \le 1$ and each $E \in J(\eta_1,\eta_2)$, give $\{\, \rho \in S \mid T(\rho,E)\ \text{is product} \,\}$ in closed form as a function of the member, through any in-file parametrization of $J(\eta_1,\eta_2)$; then give
> $$\bigcap_{E \,\in\, J(\eta_1,\eta_2)} \{\, \rho \in S \mid T(\rho,E)\ \text{is product} \,\}.$$
> *Uses: D1, D3, D4, D8.*

*Reading: the examination of convention. The first task tracks how the independence locus moves as the joint measurement ranges over every apparatus realizing the same two coins: each member has its own version of the chip, and the closed form must show it as a function of the member's free parameters. The second task, the intersection, is the sharper question: which states, if any, have independent coins no matter who does the joint measuring. Whatever survives the intersection is measurement-independent physics; whatever does not is a property of the chosen reference measurement, and this problem is what entitles the chain to tell those apart.*

**D9 (Channel).**

$$\mathrm{Ch} \;:=\; \Big\{\, \Phi : M_2 \to M_2 \;\Big|\; \Phi(\rho) = \sum_{a\in F} K_a\, \rho\, K_a^\dagger \ \text{for some finite family } (K_a)_{a\in F} \text{ in } M_2 \text{ with } \sum_{a\in F} K_a^\dagger K_a = \mathbb{I} \,\Big\}.$$

*Lemma 9a: every $\Phi \in \mathrm{Ch}$ maps $S$ to $S$ and, through Lemma 1a, induces a map $r \mapsto Ar + c$ of the closed unit ball into itself, for some $A \in \mathbb{R}^{3\times3}$ and $c \in \mathbb{R}^3$; the pair $(A,c)$ determines $\Phi$.*
*Uses: D1.*

*Reading: qubit dynamics, the chain's second and last physics input: quantum channels in Kraus form. The sum $\sum_a K_a \rho K_a^\dagger$ is the general completely positive evolution (a single Kraus term is a unitary; several terms are noise), and the normalization $\sum_a K_a^\dagger K_a = \mathbb{I}$ is conservation of probability, dual to D2's normalization of effects. Lemma 9a converts this to geometry: through the Bloch coordinates of Lemma 1a a channel acts as an affine contraction $r \mapsto Ar + c$ of the ball, and the pair $(A,c)$ carries all of $\Phi$; P3 is phrased directly on $(A,c)$ for exactly this reason.*

> **P3.** Give necessary and sufficient conditions on $(A,c)$ of Lemma 9a for $\Phi(\mathrm{Chip}) \subseteq \mathrm{Chip}$.
> *Uses: D7, D9.*

*Reading: which noises keep two independent coins independent. The deliverable is a complete test on the affine data $(A,c)$: the image of every chip state must again be a chip state. The known examples of preservers act as purely classical noise on a single coin (flip one bit with some probability, or damp one bias toward fair); the problem asks whether classical coin noise is the whole class, which would make the chip a set distinguished by its dynamics and not only by its statics.*

**D10 (Wigner-nonnegative set).**

$$G(i,j) \;:=\; \tfrac12\Big( \mathbb{I} + (-1)^j \sigma_x + (-1)^{i+j} \sigma_y + (-1)^i \sigma_z \Big), \qquad
W \;:=\; \{\, \rho \in S \;\mid\; \forall p \in \Omega:\ \mathrm{tr}\big(\rho\, G(p)\big) \ge 0 \,\}$$

(the value is real, since $\rho$ and $G(p)$ are Hermitian; the letter $G$ avoids collision with Lemma 9a's $A$).
*Lemma 10a (bridge): $G(p) = 2\sqrt3\; E^\star(p) + \tfrac{1-\sqrt3}{2}\,\mathbb{I}$, so membership in $W$ is a floor on the $E^\star$ table: $\rho \in W$ iff $T(\rho,E^\star)(p) \ge \tfrac{\sqrt3-1}{4\sqrt3}$ for all $p \in \Omega$.*
*Uses: D1, D3, D6.*

*Reading: the first rival notion of classicality: nonnegativity of a discrete Wigner function. The $G(p)$ are the phase-point operators of a qubit discrete Wigner function (several conventions exist for building one; this file fixes one), carrying the same sign pattern as $E^\star$ but with Bloch vectors of length $\sqrt3$ instead of $1$, which makes them Hermitian yet not positive: the numbers $\mathrm{tr}(\rho\,G(p))$ are quasi-probabilities that can go negative, and $W$ is the region where none does. Lemma 10a is the bridge that keeps D10 inside the chain: $G$ is an affine rescale of $E^\star$, so the Wigner values are a shifted rescale of the very SIC table everything else lives on, and Wigner classicality becomes a* floor *on that table's entries. The contrast it sets up for P4: the chip is the vanishing of the table's determinant, $W$ is a lower bound on its entries: two different functionals of one and the same object.*

**D11 (Kirkwood-Dirac-positive set).**

$$\Pi_z(i) := \tfrac12\big(\mathbb{I} + (-1)^i \sigma_z\big), \qquad \Pi_x(j) := \tfrac12\big(\mathbb{I} + (-1)^j \sigma_x\big),$$

$$Q(\rho)(i,j) \;:=\; \mathrm{tr}\big( \Pi_x(j)\, \Pi_z(i)\, \rho \big), \qquad
K \;:=\; \{\, \rho \in S \;\mid\; \forall (i,j) \in \Omega:\ Q(\rho)(i,j) \in \mathbb{R} \ \text{and}\ Q(\rho)(i,j) \ge 0 \,\}.$$

*Uses: D1.*

*Reading: the second rival notion: Kirkwood-Dirac positivity for the sharp $\sigma_z$/$\sigma_x$ pair. $\Pi_z(i)$ and $\Pi_x(j)$ are the sharp eigenprojectors of the two observables the coins unsharply measure, and $Q(\rho)(i,j)$ pairs the state with the ordered product $\Pi_x\Pi_z$; because the two projectors do not commute, that product is not Hermitian and $Q$ is in general complex, which is exactly why D4 defined $\det$ on complex tables. Still, $Q$'s row and column marginals are the genuine sharp $\sigma_z$ and $\sigma_x$ statistics, so $Q$ is the sharp-limit sibling of the coin table, with quasi-ness as the price of sharpness. $K$ is the set where all four entries are real and nonnegative, the KD-positive states; P4 will place the chip against both this and $W$, and apply the same independence quantity $\det$ to $Q$ itself.*

> **P4.** Determine the inclusion and intersection relations among $\mathrm{Chip}$, $W$, and $K$, and the locus $\{\, \rho \in S \mid \det Q(\rho) = 0 \,\}$.
> *Uses: D4, D7, D10, D11.*

*Reading: locate the chip on the map of classicality. The first task is set relations: whether zero correlation of the coins is contained in, contains, or cuts transversally across each of the two standard positivity notions, decided by inclusions, intersections, and witness states. The second task applies the chip's own defining quantity to the rival object: $\det Q = 0$ is the independence condition imposed on the sharp-pair quasi-table, the exact analogue of D7 one level of sharpness up, and its locus completes the picture of what the independence condition does across the file's three tables (probability, floor-bounded, quasi).*

## Two questions that cannot be formalized (outside the chain, listed for completeness)

> **M1 (worth).** What tasks do members of $\mathrm{Chip}$ enable that generic states do not?
> **M2 (priority).** In which literatures, under which vocabularies, do D4 applied to quantum tables, or P1 itself, already appear?

*Reading: the assessment duties of a research paper, deliberately kept outside the formal chain: M1 is the operational-payoff question and M2 the prior-art question. Neither is a mathematical statement about the objects above: no computation on D1-D11 can answer them, which is precisely the test that expelled them from the formal list; in a paper they become the applications and related-work sections rather than theorems.*

## Why nothing can be removed

- **Foundation:** D1-D4 (state, measurement, Born table, determinant/product). Remove any one and D7/P1 can no longer be stated. Lemmas 1a, 3a, 4a establish that the objects have the properties their use assumes and supply the coordinates the solvability note relies on; D4's $\det$ is needed twice (through "product" in D7, directly in P4).
- **Reference-measurement block:** D5 (relabeling; needed by Lemma 6a's uniqueness clause and by P1's second task), D6 (the reference measurement; its formula fixes both the tetrahedron and the grouping into two bits, so no silent choice remains), D7 (gives the object a name, so P3 can refer to it directly instead of by pronoun). The unitary notion in the Primitives is needed only by Lemma 6a's uniqueness clause; delete that clause and the notion goes with it.
- **Lemma roles:** 6a, 8a, 8b, 10a are certification: they are why $E^\star$ may be called "the SIC", why it belongs to D8's family, and how D10 connects back to the chain. 6b and 6c are interpretation: they make "the state is two coins" the meaning of P1. No problem's computation uses 6a, 6b, or 6c, given the preamble's rule that a definition's membership claims are covered by its own lemmas.
- **What each problem adds:** P2 needs only D8 beyond the foundation; P3 only D9; P4 only D10 and D11. Each added block is used by its problem and by nothing else, and each problem fails without its block.
- **M1, M2** are not mathematical statements about the objects above; keeping them inside the formal list was a mistake of earlier drafts, recorded here so it is not repeated.

## Record of the checks

- Round 1 (nine findings, repaired): Dirac bra-ket notation used without being declared; "Bloch angle" and "unitary" used without definitions; incomplete dependency lists; the coin-pair family ranged over sharpness values for which no joint measurement exists (the repair is Lemma 8b's disk $\eta_1^2+\eta_2^2 \le 1$); real-versus-complex mismatches; an answer leaked in the sketch.
- Round 2 (seven findings, repaired): the uniqueness claim of Lemma 6a was **false** as then written, with an exact counterexample (rank-one measurements with all pairwise overlaps equal to $s(1-s)/3$ and weights $\{s,s,1-s,1-s\}$, $s \ne \tfrac12$); the determinant needed to accept complex tables for P4 to make sense; an expectation symbol collided with $E$; dependency lists were inconsistent; the pair $(A,c)$ was never quantified.
- Both reviewers independently solved P1 and computed the locus in P4 using this file alone; their answers agree and are deliberately not recorded here, so the file remains a clean problem set.

## Sketch in Lean notation (shape only; this code does not compile)

```lean
-- schematic: names and shapes, not working mathlib code
def State := { ρ : Matrix (Fin 2) (Fin 2) ℂ // ρ.PosSemidef ∧ ρ.trace = 1 }
def TableMeas := { E : Fin 2 × Fin 2 → Matrix (Fin 2) (Fin 2) ℂ // (∀ p, (E p).PosSemidef) ∧ ∑ p, E p = 1 }
def table (ρ : State) (E : TableMeas) : Fin 2 × Fin 2 → ℝ := fun p => ((ρ.val * E.val p).trace).re
def IsProduct (v : Fin 2 × Fin 2 → ℝ) : Prop := v (0,0) * v (1,1) = v (0,1) * v (1,0)
def Estar : TableMeas := ⟨fun (i,j) => (1/4 : ℝ) • (1 + (Real.sqrt 3)⁻¹ • ((-1:ℝ)^j.val • σx + (-1:ℝ)^(i.val+j.val) • σy + (-1:ℝ)^i.val • σz)), by admit⟩
def Chip : Set State := { ρ | IsProduct (table ρ Estar) }
-- P1 : give Chip in closed form (deliberately not stated here)
```
