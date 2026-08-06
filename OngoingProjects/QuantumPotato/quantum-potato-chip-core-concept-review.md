# The core concept: is the chip derived and explained to world-class standard?

Scope: only the core concept (the chip itself), judged on two axes: quality of derivation, quality of exposition. The error-by-error audit lives in [quantum-potato-chip-tex-audit.md](quantum-potato-chip-tex-audit.md). Every mathematical statement below was verified in a kernel this session (exact arithmetic).

**Verdict.** The computations behind the chip are correct, but the chip is never actually *derived*: it is parametrized. The paper builds the surface as the image of a construction and never states, let alone proves, the invariant condition that defines it. And it is not explained: the object is named before it is defined, motivated by the wrong property, and narrated cell-by-cell in notebook voice. Both faults are repairable, because the concept underneath is genuinely good and has a ten-line world-class derivation the paper walks past.

---

## 1. The derivation as it stands

The paper's route (Sec. 2): represent the state as a SIC probability 4-vector; rotate R⁴ by the θ=π/3 matrix of Eq. (1)-(2) and drop a coordinate to get a tetrahedron in R³; push the outer product {p,1-p}⊗{q,1-q} through the rotation to get a parametrized surface; assert (footnote) that physical states form the insphere; intersect. Then Sec. 2.2 repeats the construction in Bloch coordinates and calls it a second method.

Seven structural judgments:

**D1. The rotation is scaffolding doing no conceptual work, and its one relevant property is never named.** The factorization condition and the physicality condition are both statements about the 4-vector itself; neither needs R³. What the rotation actually provides is an *isometry* of the simplex's affine hull onto coordinate R³ (isometry is what keeps the physical ball round), any such map works, and θ=π/3 is in no way canonical. The text introduces it as "Consider the following rotation matrix", which is the direct transcript of the notebook's `RotationMatrix[θ,{{1,1,1,1},{1,0,0,0}}]` cell. A derivation that opens with unexplained machinery whose purpose is visualization has inverted its priorities.

**D2. The object is only ever parametrized, never defined.** Section 2 proves "every product vector lands on this surface" (an inclusion) and never asks the converse question "which states have factorizing statistics?". The defining equation is one line and appears nowhere:

> v equals the product of its own two binary marginals ⇔ det [[p₁,p₂],[p₃,p₄]] = 0,

because for a 2×2 probability table the determinant *is* the covariance of the two marginal bits. In Bloch coordinates (verified exact):

> det = ( √3⟨σ_y⟩ − ⟨σ_x⟩⟨σ_z⟩ ) / 12,

so the chip is exactly {√3 y = xz} ∩ {|r| ≤ 1}, and the parametrized image coincides with this zero set inside the ball (both inclusions checked). The paper's own Matthews section secretly contains this fact as φ = 0, but files it as a downstream "property", and establishes exactness by pointing at a plot: "As shown in Fig. 7, the only region with φ=0 is the quantum potato chip". A figure stands where a one-line proof belongs.

**D3. Physicality is used before it is derived.** The intersection step in Sec. 2.1 rests on "the only region ... is the sphere that is inscribed within it" plus a garbled rescaling footnote (audit E11). The genuine proof (eigenvalues, 3ω = x²+y²+z²) arrives one section later and is never connected back. The logic of the central construction is therefore circular-by-forward-reference at the point of first use.

**D4. The "two distinct methods" are one method.** Sec. 2.2 conjugates the same construction by the linear map between probability space and Bloch space. Genuinely distinct derivations were available: algebraic (the determinant condition), operational (independence of the two unsharp Pauli marginals), group-theoretic (the three pairings and their relabeling orbit). The paper has fragments of all three and derives none as such.

**D5. Only the constructive direction is run.** World-class treatments of a locus run both directions: build it, then characterize it. Running only "build it" is precisely why the paper ends up advertising the wrong uniqueness property (two-measurement reconstruction, which any promised 2-parameter family enjoys) instead of the right one (zero coin covariance, which is exactly the chip).

**D6. The object's mathematical identity is never extracted.** The surface y = xz/√3 is a hyperbolic paraboloid: a doubly ruled quadric, and literally the shape of a Pringle, so the paper's name is mathematically earned and the paper never says it. The two ruling families are exactly the constant-p and constant-q lines: the chip is woven from straight lines of constant coin bias. The boundary is the quartic curve (sphere ∩ quadric) through |0⟩, |1⟩, |+⟩, |−⟩. None of this appears.

**D7. Structure that falls out in lines is missing** (all verified this session):
- The covariance identity above, which is the physical headline: in the SIC frame, the y coherence measures, up to the ⟨σ_x⟩⟨σ_z⟩ product term, the *classical correlation between the qubit's Z coin and X coin*. The chip is the zero-correlation locus. Stated nowhere.
- Three chips = the three ways to arrange four outcomes as a 2×2 product table; the 24 relabelings split 8 (stabilizer) + 8 + 8; a CNOT relabeling gives the second chip and a reversed CNOT the third (the printed footnote pairs two relabelings from the same coset, audit E9).
- Each Cartesian axis lies on exactly two chips: each chip contains its own two coin axes (they are its rulings through the center) and excludes the third axis, its correlation direction, which is where |φ| is maximal on the sphere (±1/√3). All three chips intersect only at the maximally mixed state: the fair coin pair.

On top of these structural faults sit the literal transcription errors in the derivation chain (wrong fiducial E1, wrong border equation E6, swapped marginal labels E7); those are typos, D1-D7 are not.

---

## 2. The derivation as it should be (ten lines, no rotation)

1. **Frame.** Qubit SIC: Q_i = ¼(I + n_i·σ) with tetrahedral unit vectors, n_i·n_j = −1/3 (i≠j). State ρ = ½(I + r·σ) ↔ p_i = ¼(1 + r·n_i). Every qubit SIC is a rotation of this one, so nothing is lost.
2. **Physicality, both directions.** |p − c|² = |r|²/12 for c the simplex center, and the insphere radius of the 3-simplex is 1/(2√3); hence physical states ⇔ the insphere, exactly, with proof, before first use.
3. **The table.** Arrange outcomes as T = [[p₁,p₂],[p₃,p₄]]. The row coarse-graining {Q₁+Q₂, Q₃+Q₄} is the unsharp Pauli-Z measurement with outcome probability p = ½(1 − z/√3); the columns give unsharp Pauli-X with q = ½(1 − x/√3). The two coins are the qubit's Z and X statistics from the first minute, and the sharp Pauli data recover them through the fixed doubly stochastic map S.
4. **Definition.** The chip is the set of states for which T is the product of its own marginals.
5. **Theorem.** For 2×2 tables, product ⇔ det T = 0 ⇔ zero covariance of the two bits; and det T = (√3 y − xz)/12. Hence the chip is exactly the hyperbolic-paraboloid section {√3 y = xz} of the Bloch ball. Both inclusions are immediate; the parametrization r(p,q) = √3(1−2q, (2p−1)(2q−1), 1−2p) and its (p,q)-rulings are corollaries, as is the pure boundary quartic through |0⟩,|1⟩,|±⟩.
6. **Structure.** Three pairings → three chips, permuted by CNOT-type relabelings (S₄ = 8+8+8); pairwise intersections are the shared coin axes; triple intersection is the fair-coin state. Wootters twin: same construction on the phase-point frame gives y = xz with no √3 and automatic nonnegativity.
7. **Operational corollary, honestly stated.** Under the promise that the state lies on a known chip, sharp X and Z measurements determine it (outer product after S). The promise clause is part of the statement; without a promise two measurements never determine a qubit, and with some promise any 2-parameter family is determined.

Everything the paper computes survives; the difference is that the object now has a definition, a theorem, a shape, and a meaning, in that order.

---

## 3. The exposition as it stands

**E1. The motivating question is never asked.** The reader meets a rotation matrix before any physics. The question the whole paper answers, "which qubit states carry the statistics of two independent classical coins?", appears nowhere in the introduction; its two "key questions" describe simplex operations, not meaning.

**E2. The name arrives before the object.** The christening sentence is also the broken one ("This is manifold looks like a potato chip we shall refer to it as quantum potato chips"), and the shape claim is never cashed out (the surface really is the Pringle quadric, which would earn the name mathematically; unsaid).

**E3. Four representation layers with no dictionary.** Probability simplex, rotated tetrahedron, Bloch ball, and Wootters quasi-probabilities, with three different radius conventions (1, 1/(2√3), 1/√3), are interleaved with no table telling the reader where each equation lives. The insphere footnote confusion and the "radius 1/√3" remark before Eq. (19) are symptoms. The Matreshka figure is the one place the layers meet, in a caption.

**E4. The coins' identity is withheld, then swapped.** That p and q are the (rescaled) Z and X statistics is revealed only in Sec. 3, and there with the labels swapped relative to the paper's own worked example (audit E7). The single best explanatory sentence in the paper, "much like flipping a pair of independent biased coins" (main.tex:494), is buried mid-Sec. 3.2; it belongs in the abstract.
 
**E5. The wrong "why it matters".** The exposition hangs the chip's significance on two-measurement reconstruction (refutable, see D5) rather than on the property that is actually unique and physical: zero correlation between the two coins, φ = 0. The correct headline is demoted to a "property" section.

**E6. Surface and boundary blur.** The chip is a 2D surface of mostly mixed states; its boundary is the pure curve. The Liouvillian section's title says "boundary of quantum potato chips" while the introduction speaks of states "within" the chip; the reader is never told which object is in play, or that the boundary is where the eigenvalue vanishes.

**E7. Notebook voice throughout.** "Consider the following rotation matrix", "one gets", "one can drop the redundant first element": each cell's operation is narrated, the arc is not. Sentences at definitional moments are the broken ones; terminology collides with standard usage ("separable" for independence, "product state" for a single-qubit distribution, "Weyl-Wigner transformation" for applying the dual frame).

**Assets already in hand.** The figures are genuinely good (the tetrahedron-with-saddle, the chip in the Bloch ball, the Matreshka, the φ contours); the p=1/3, q=2/5 worked example is exactly the right pedagogy; the coin sentence exists. The world-class exposition is a reordering plus a definition, not a rewrite from nothing: question → table and coins → covariance = determinant → chip as its zero set (shape named) → structure (three chips, axes, rulings) → promise-based tomography corollary → contrast with the Wigner-nonnegativity and Kirkwood-Dirac positivity literatures.
