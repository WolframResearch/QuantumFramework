# Qubit states whose SIC-POVM statistics factorize

**Abstract.** A single run of the qubit symmetric informationally complete measurement returns two bits, and the joint distribution of these bits carries the complete quantum state. We determine the states for which the two bits are statistically independent, so that the entire information content of the state is a pair of independent biased coins. In Bloch coordinates the set is the surface $\sqrt{3}\,\langle\sigma_y\rangle = \langle\sigma_x\rangle\langle\sigma_z\rangle$; relabeling the four outcomes produces two further copies, and the three surfaces are permuted by the tetrahedral symmetry of the measurement. Each surface is a doubly ruled saddle clipped by the Bloch ball: its straight rulings hold one coin fixed while the other varies, and its rim of pure states is a quartic curve through the four eigenstates of $\sigma_x$ and $\sigma_z$.

## I. Setup

Measured once, an informationally complete measurement turns a quantum state into a probability distribution that determines the state. For a qubit the minimal such measurement has four outcomes, matching the four real parameters of the state counting normalization, and four outcomes are naturally labeled by two bits. A sharp question follows: for which states are the two bits statistically independent, so that the state is, informationally, nothing more than a pair of independent biased coins? We answer this in closed form and describe the geometry of the answer.

A qubit state is a $2\times2$ complex matrix $\rho$, positive semidefinite with unit trace. Writing $\mathbb{I}$ for the identity, $\vec\sigma = (\sigma_x, \sigma_y, \sigma_z)$ for the Pauli matrices, $\mathrm{tr}$ for the trace, and $\langle M\rangle := \mathrm{tr}(\rho M)$ for the expectation of an observable $M$, every state takes the Bloch form

$$\rho = \tfrac{1}{2}\big(\mathbb{I} + \vec r\cdot\vec\sigma\big), \qquad \vec r = (x, y, z) := \big(\langle\sigma_x\rangle,\ \langle\sigma_y\rangle,\ \langle\sigma_z\rangle\big),$$

and the eigenvalues of $\rho$ are $\tfrac12(1 \pm |\vec r|)$: positivity of the state is exactly the closed unit ball $|\vec r| \le 1$, with the pure states on the sphere.

A measurement with finitely many outcomes assigns to each outcome an effect, a positive semidefinite matrix, the effects summing to $\mathbb{I}$; such a family is a positive operator-valued measure (POVM), and the Born rule gives the probability of the outcome with effect $E$ as $\mathrm{tr}(\rho E)$. We take the four-outcome measurement whose outcomes carry a pair of bits $(i,j) \in \{0,1\}^2$, with effects

$$E_{ij} = \tfrac14\big(\mathbb{I} + \vec n_{ij}\cdot\vec\sigma\big), \qquad \vec n_{ij} = \tfrac{1}{\sqrt3}\big((-1)^j,\ (-1)^{i+j},\ (-1)^i\big).$$

The four unit vectors $\vec n_{ij}$ are the vertices of a regular tetrahedron: writing $p$ and $q$ for outcome labels $(i,j)$, one has $\vec n_p\cdot\vec n_q = -\tfrac13$ for $p \ne q$, and the four vectors sum to zero, so the effects sum to $\mathbb{I}$. Each effect is half of a rank-one projector, and all pairwise overlaps are equal, $\mathrm{tr}(E_p E_q) = \tfrac1{12}$: these are the defining properties of the qubit symmetric informationally complete POVM (SIC-POVM). The bit labels record the signs of the vertex components: $i$ is the sign bit of the $z$ component, $j$ that of the $x$ component, and the mod-2 sum $i \oplus j$ that of the $y$ component.

By the Born rule a state $\rho$ produces the $2\times2$ probability table $T(\rho)$ with entries

$$p_{ij} = \mathrm{tr}\big(\rho\, E_{ij}\big) = \frac14\left[\,1 + \frac{(-1)^j\,x + (-1)^{i+j}\,y + (-1)^i\,z}{\sqrt3}\,\right],$$

the bit $i$ labeling the row and $j$ the column. One run of the measurement is thus a pair of random bits, in general correlated, and the opening question becomes: for which $\vec r$ does this table factorize?

## II. Two coins in one measurement, and a third

Summing the table over one bit leaves the distribution of the other. Because the sign patterns of the $x$ and $y$ terms average away across a row and those of the $y$ and $z$ terms down a column, the marginals are

$$a_i := p_{i0} + p_{i1} = \tfrac12\big(1 + (-1)^i\,\eta\, z\big), \qquad b_j := p_{0j} + p_{1j} = \tfrac12\big(1 + (-1)^j\,\eta\, x\big), \qquad \eta := \tfrac1{\sqrt3}.$$

The same statement holds at the level of effects,

$$E_{i0} + E_{i1} = \tfrac12\big(\mathbb{I} + (-1)^i\,\eta\,\sigma_z\big) = \eta\,\Pi^z_i + (1-\eta)\,\tfrac{\mathbb{I}}{2},$$

where $\Pi^z_i := \tfrac12(\mathbb{I} + (-1)^i\sigma_z)$ projects onto the $\sigma_z$ eigenstate of eigenvalue $(-1)^i$. The row bit is therefore an unsharp $\sigma_z$ measurement: with probability $\eta$ the projective measurement is performed, and otherwise a fair coin is returned. It is unbiased, both outcomes being equally likely at the maximally mixed state, and its sharpness is $\eta = 1/\sqrt3$: the outcome statistics respond to the state only through $\eta\langle\sigma_z\rangle$. The column bit is in the same way an unbiased unsharp $\sigma_x$ measurement of the same sharpness. We call a binary measurement of this kind a coin, and the mean of its $\pm1$-valued outcome its bias. The sharpness pair $(\eta,\eta)$ obeys $\eta^2 + \eta^2 = \tfrac23 < 1$, placing it strictly inside the region $\eta_1^2 + \eta_2^2 \le 1$ known to admit joint measurement of unsharp $\sigma_z$ and $\sigma_x$ coins; the SIC-POVM is one such joint apparatus.

Attach to the outcome the two $\pm1$-valued variables $A := (-1)^i$ and $B := (-1)^j$, and for a function $f$ of the outcome write $\langle f\rangle := \sum_{ij} p_{ij}\, f(i,j)$; the two uses of $\langle\cdot\rangle$ are told apart by their argument. The product $AB = (-1)^{i+j}$ is a third binary variable, and it has its own effect pair, $E_{00} + E_{11} = \tfrac12(\mathbb{I} + \eta\,\sigma_y)$ and $E_{01} + E_{10} = \tfrac12(\mathbb{I} - \eta\,\sigma_y)$: an unbiased unsharp $\sigma_y$ coin, again of sharpness $\eta$. Altogether

$$\langle A\rangle = \eta\,\langle\sigma_z\rangle, \qquad \langle B\rangle = \eta\,\langle\sigma_x\rangle, \qquad \langle AB\rangle = \eta\,\langle\sigma_y\rangle .$$

A single run of the SIC-POVM thus serves three coins at once, one per Pauli direction, all of the same sharpness, subject to one multiplicative constraint: the third variable is, outcome by outcome, the product of the first two. Inverting,

$$(x,\ y,\ z) = \sqrt3\,\big(\langle B\rangle,\ \langle AB\rangle,\ \langle A\rangle\big),$$

so the table determines the state completely; this is the informational completeness of the measurement. A qubit state is exactly two coin biases together with the correlation of the coins.

## III. The independence condition

Let $\det T := p_{00}\,p_{11} - p_{01}\,p_{10}$. This determinant is the entire obstruction to factorization for any $2\times2$ array of reals with unit sum: writing $a_i$ and $b_j$ for the marginals of the array as above, elimination of $p_{11} = 1 - p_{00} - p_{01} - p_{10}$ verifies the identity

$$p_{ij} = a_i\, b_j + (-1)^{i+j}\,\det T ,$$

so the table equals the product of its own marginals if and only if $\det T = 0$. In the same notation,

$$4\,\det T = \langle AB\rangle - \langle A\rangle\langle B\rangle ,$$

the covariance of the two bits. The two identities express the special luck of binary alphabets: for a pair of two-valued variables, vanishing covariance is not merely necessary for independence but equivalent to it, and the entire independence condition is a single scalar equation. Evaluating on the SIC table through the coin expectations of the previous section,

$$\det T(\rho) = \frac{\langle AB\rangle - \langle A\rangle\langle B\rangle}{4} = \frac{\sqrt3\,\langle\sigma_y\rangle - \langle\sigma_x\rangle\langle\sigma_z\rangle}{12} .$$

The two bits of the SIC-POVM are therefore independent precisely on the locus

$$\sqrt3\, y = x\, z, \qquad x^2 + y^2 + z^2 \le 1 .$$

Read through the coins the condition is transparent: independence forces $\langle AB\rangle = \langle A\rangle\langle B\rangle$, so the bias of the $\sigma_y$ coin is not free but equal to the product of the other two biases. The locus is far from empty. It contains the maximally mixed state at the origin, the full $x$ and $z$ diameters of the ball, and, as shown below, a closed curve of pure states; being one equation inside a three-parameter ball, it is a surface.

## IV. Relabeling the outcomes: three surfaces

Grouping the four outcomes into a row bit and a column bit was a choice, and relabeling outcomes is physically free, so the locus must be examined under the $24$ bijections of the outcome set. After a relabeling, each of the two new bits is a binary function of the old outcome pair that takes each value twice, and the balanced binary functions on four outcomes are exactly the three parities $i$, $j$, $i\oplus j$ and their negations. Independence of two bits is unchanged by negating either bit or by exchanging the two, so the relabeled independence condition depends only on the unordered pair of distinct parities that the new bits realize: three possibilities, and the count $24 = 3 \times 2 \times 4$ over pairs, orders, and negations shows every bijection is accounted for. Each possibility declares two of the three coins $A$, $B$, $AB$ to be the bit pair; since the pointwise product of any two of these variables is the third, each condition again states that the omitted coin's bias equals the product of the two declared biases:

$$\sqrt3\, y = xz \quad (i \text{ with } j), \qquad \sqrt3\, x = yz \quad (i \text{ with } i\oplus j), \qquad \sqrt3\, z = xy \quad (j \text{ with } i\oplus j),$$

always within the unit ball. The three surfaces are congruent: the rotation by $2\pi/3$ about the tetrahedron vertex $\vec n_{00} = (1,1,1)/\sqrt3$, implemented on states by a unitary conjugation, permutes the three Pauli directions cyclically, maps the SIC-POVM onto a relabeling of itself, and carries the three surfaces onto one another. Any two of them meet in a coordinate diameter of the ball: for the first two surfaces, combining the equations gives $x\,(3 - z^2) = 0$ and $y\,(3 - z^2) = 0$, and $z^2 < 3$ throughout the ball, so the intersection is the $z$ axis. All three surfaces meet only at the maximally mixed state, the one state whose three coins are pairwise independent under every grouping.

## V. Geometry: rulings, rim, and the saddle

On the locus, take as coordinates the two coin biases themselves, $\mu := \langle A\rangle = \eta z$ and $\nu := \langle B\rangle = \eta x$. Then

$$(x,\ y,\ z) = \sqrt3\,\big(\nu,\ \mu\nu,\ \mu\big),$$

and the ball constraint factorizes:

$$|\vec r|^2 = 3\big(\mu^2 + \nu^2 + \mu^2\nu^2\big) \le 1 \iff \big(1+\mu^2\big)\big(1+\nu^2\big) \le \tfrac43 .$$

The correspondence between $(\mu,\nu)$ and the state is one to one: the surface is globally a graph over this region of the coin-bias plane, and every admissible bias pair is realized by exactly one state. The rim of the surface, where the state is pure, is the curve $(1+\mu^2)(1+\nu^2) = \tfrac43$; in Bloch coordinates, substituting $y = xz/\sqrt3$ into $|\vec r| = 1$ gives the quartic

$$\big(3 + x^2\big)\big(3 + z^2\big) = 12, \qquad y = \frac{x z}{\sqrt3},$$

a closed curve on the Bloch sphere passing through the four eigenstates of $\sigma_z$ and $\sigma_x$, the points $(x,z) = (0,\pm1)$ and $(\pm1,0)$. These four states are the extreme points of the rim in the biases: on the rim, $|\mu|$ attains its maximum $\eta$ only at $\nu = 0$, and symmetrically, so a coin reaches its sharpness-limited extreme bias only when its partner is fair.

Slicing instead of substituting exposes the internal structure. Fix the $\sigma_x$ coin, $\nu = \text{const}$: the slice of the surface by the plane $x = \sqrt3\,\nu$ is the segment $\mu \mapsto \sqrt3\,(\nu,\ \mu\nu,\ \mu)$, affine in $\mu$, hence straight; fixing $\mu$ works the same way. Through every point of the surface pass two straight lines lying entirely on it. The lines have a statistical meaning. A straight segment in the ball is a family of mixtures of its endpoint states, and if two tables $a^{(1)}_i b_j$ and $a^{(2)}_i b_j$ factorize with a common column marginal $b_j$, their mixture with weight $\lambda \in [0,1]$ equals $\big[\lambda\, a^{(1)}_i + (1-\lambda)\, a^{(2)}_i\big]\, b_j$ and factorizes as well: independence survives mixing exactly when one coin's bias is held fixed. It does not survive general mixing: the equal mixture of the pure states at $z = 1$ and at $x = 1$ leaves the surface, which is not a convex set.

Finally, rotate the frame by $\pi/4$ about the $y$ axis, $u := (x+z)/\sqrt2$ and $w := (z-x)/\sqrt2$. The equation $\sqrt3\,y = xz$ becomes

$$y = \frac{u^2 - w^2}{2\sqrt3},$$

the normal form of a hyperbolic paraboloid: a saddle, one of the two doubly ruled quadric surfaces, with the two families of straight lines found above as its rulings. A remark on nomenclature is now in order. A hyperbolic paraboloid clipped to a bounded patch is the familiar saddle shape of a fried potato chip, and we call the locus the quantum potato chip; the relabelings of the previous section then place three congruent chips inside one Bloch ball.

## VI. Discussion

On the chip the third coin is redundant, its bias being the product of the other two, so the bias pair $(\mu,\nu)$ is a complete description: a chip state is, informationally, a pair of independent biased coins, and every pair in the region $(1+\mu^2)(1+\nu^2) \le \tfrac43$ occurs exactly once. The shape of that region is itself physics. Allowing correlation, a qubit can display any bias pair in the disk $3(\mu^2 + \nu^2) \le 1$, which is strictly larger: there are bias pairs, for instance $\mu = \nu = \tfrac25$, that a qubit can carry only by correlating its coins. At the opposite extreme sit the two eigenstates of $\sigma_y$: both coins are individually fair there, yet the covariance $\langle AB\rangle - \langle A\rangle\langle B\rangle = \tfrac13(\sqrt3\, y - xz)$ is extremal there over the whole ball, reaching $\pm\eta$: these are the states as far from coin independence as the sharpness permits.

Two companion questions are left open. The first is dynamical: which quantum channels map the chip into itself. Classical noise on a single coin does; flipping the $\sigma_z$ coin with probability $p$, enacted by conjugation with $\sigma_x$ mixed into the identity, damps that coin's bias by $1-2p$, leaves the $\sigma_x$ coin untouched, and preserves the factorization of every chip state. Whether the preservers go beyond classical coin noise is the question. The second is comparative: the chip is the vanishing of a determinant of the SIC table, while the standard positivity notions of qubit classicality constrain close relatives of the same object, nonnegativity of a discrete Wigner function built on the same tetrahedron being a floor on the entries of that very table, and Kirkwood-Dirac positivity for the sharp $\sigma_z$, $\sigma_x$ pair being positivity of the table's sharp-marginal, generally complex sibling. How the chip sits among those sets deserves a treatment of its own.
