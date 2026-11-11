# Bosonic minimal design — status

Tracks Mads's `Minimal-Bosonic-Design.md` in this folder. Every ✅ below was checked with the declared-vars API by `tracker-probe.wls` (run with `wolframscript -file tracker-probe.wls`); each check compares against the formula in Mads's table, not just "it evaluated".

Legend: ✅ works · 🟡 partly works · ⬜ needs development (idea number) · ➖ out of scope for an exact operator layer (his class d)

## Summary

| | Rows |
|---|---|
| ✅ works now | ~~44~~ ~~47~~ ~~49~~ 51 (his 29 class a, plus B-053 via Idea 7, B-020 via Ideas 2 and 5, G-012, G-038, G-041 via Idea 5, B-026, B-027, P-020, C-031, G-024 via Idea 2, G-040, B-106, X-003 via Idea 4, P-017 and G-027 from class d via Idea 2, C-029, C-032, C-050, X-006, P2-024 via Idea 1, C-022, C-054 via Ideas 1 and 2) |
| 🟡 partly | 3 (P-026, P-030, X-007: one part works) |
| ⬜ needs an idea | ~~16~~ ~~13~~ ~~11~~ 9 |
| ➖ out of scope | 66 |

## Ideas

| Idea | What | Status | Rows it unlocks |
|---|---|---|---|
| 1 | one rule for `Exp` of any quadratic (2N×2N matrix, Gauss factors) | ~~🟡~~ ✅ single mode done (~~`Exp[x b**b + y bd**bd + w bd**b]` returned unchanged~~), quadratic plus linear: `quadraticExp` reads the normal-ordered exponent, disentangles through the su(1,1) matrix and shifts the linear terms away (closed-form u, v and scalar, regular in λ = √(4pq − w²), so the parabolic λ = 0 case such as e^{t x²} is its `Series` limit); antinormal is the normal formula under a → a†, a† → −a. Replaces the pure-squeeze rule with the same output. `adFamilyRules` braids any single-mode generator of degree ≤ 2 (one affine rule replacing the linear and squeeze cases). Checked by exact series of symbols and of Fock elements. ~~Multimode ⬜~~ ~~Multimode 🟡~~ Multimode ✅: coupled pure-quadratic generators (`multimodeQuadraticExp`): Gauss factors from the blocks {{A, B}, {C, D}} of exp[{{W, P}, {−Q, −Wᵀ}}], e^{−Tr W/2} det D^{−1/2} e^{a†BD⁻¹a†/2} :e^{a†(D⁻ᵀ−1)a}: e^{−aD⁻¹Ca/2}, and :e^{a†(e^W−1)a}: for P = Q = 0; antinormal the same under a → a†, a† → −a; uncoupled modes split into single-mode factors. The root det D^{−1/2} takes its winding phase from `rootPhase` (dominant Cauchy–Binet term of det D over subsets of the eigenvalues), continuous when the generator is t × a numeric matrix, principal otherwise. Checked by exact series to t³ in both orders, and the root against the closed form and 30-level evolutions to t = 12. ~~Linear terms over coupled modes ⬜~~ Linear terms over coupled modes ✅ (`affineShift`): e^{G₂+L} = e^{κ} e^{L̃} e^{G₂} with L̃ = (a†, a)·φ(M)ℓ, κ = ½ ℓᵀJψ(M)ℓ, M = S Xᵀ S the adjoint matrix on (a†, a), both from one augmented MatrixExp (regular for singular X, so a pure linear G is included), then e^{L̃} braided into the wings. Checked by exact series to t³ in both orders (hopping, squeezing with and without number terms, a constant), and a displaced two-mode squeezer against `MatrixExp` on 22 and 25 levels (≤ 2·10⁻¹⁴) | 10 |
| 2 | coherent elements of ordered products, ⟨α\|e^{P(a†)} f(N) e^{Q(a)}\|β⟩ | ✅ done 2026-11-04, committed `62a54595` (`normalProductElement` with a coherent flag; checked symbolically against D(g)\|β⟩ = e^{(gβ̄−ḡβ)/2}\|β+g⟩ and numerically against 60-level matrices for D, normal and antinormal products, squeeze, wing–f(N)–wing) | 9 |
| 3 | s-ordered (Weyl) symbol, `NormalOrdered[expr, vars]["Symbol", s]` | 🟡 prototype outside QF (Mads's choice), `s-ordered-symbol.wls`: e^{−(1−s)/2 Σ∂_{A_i}∂_{B_i}} on the normal symbol, exact on polynomial × Gaussian terms (a linear term j added to the exponent, differentiated and set to 0, plus the Gaussian heat formula), any number of modes. Checked: D(γ) for all s, s = ±1 against normal and antinormal order, C-053, Weyl symbols of \|n⟩⟨n\| for n ≤ 4, C-016 (quadratic H agrees with the Poisson bracket, quartic does not). Returns `$Failed` when `BosonicExpOrder` leaves a mixed exponential unordered. Not in QF: a property cannot attach to polynomials, since `NormalOrdered` of a polynomial evaluates away | 5 |
| 4 | several modes in the exponential branch of `BosonicMatrixElement` | ✅ functions of the N_k alone (Fock basis: ⟨m\|f(N₁, …)\|n⟩ = f(n) δ_{mn}, substituting each cr_k ** an_k by n_k). Wings ✅ (~~`⟨1,1\|Exp[z bd1**bd2]\|0,0⟩` unevaluated~~): `multimodeProductElement` turns creators, annihilators and N_k into commuting x, y, k, reads F(x) f(k) G(y), and sums the multivariate Taylor coefficient of F G over 0 ≤ k ≤ min(m, n); integer indices, Fock basis. Coherent basis with several modes and symbolic multimode indices not covered | 4 |
| 5 | f(N) beside bare ladder operators | ✅ Fock/coherent part committed `f76ab7da` (`"Creators"`/`"Annihilators"` tags in `normalFactor`; ⟨m\|a†^p f(N) a^q\|n⟩ in closed form, zero for n < q; checked against 50-level matrices, also between wings and in the coherent basis). Shift rule committed `1d307d37` (`shiftRules` in `BosonicExpSimplify`: a^q f(N) = f(N+q) a^q, f(N) a†^p = a†^p f(N+p), functions of N merge with `Together`, a balanced word beside f(N) becomes a polynomial in N; a† f(N) a left alone for the vacuum). Checked: a f(N) = f(N+1) a against 40-level matrices, ⟨1\|b ** f(N)\|2⟩ = √(2/3), E E† = E² E†² = 1, f-deformed A A† = (N+1) g(N+1)², Holstein–Primakoff S₊S₋ = (2j−N)(N+1). Not covered: another mode between the ladder word and f(N) | 4 |
| 6 | su(1,1)/su(2) with free label k | ⬜ recipe only — `NonCommutativeAlgebra` with the commutation relations already works in plain WL; his note builds on `BosonicRelations`, which was removed | 3 |
| 7 | conjugation linear over sums and scalars | ✅ done 2026-11-03 (BosonicExp.m `conjugateParts`, committed `1a4ca028`) | 2 |

## ✅ Works now (his class a, re-checked)

| Row | Formula | Call |
|---|---|---|
| ✅ P-014 / P2-050 | D(u)D(w) vs D(w)D(u), ratio e^{2i Im(u w̄)} | `BosonicExpSimplify[D ** D, v]` both orders |
| ✅ P-006 | displacement as e^{iPq} e^{-iQp} | `BosonicExpSimplify` with q, p written in b, b† |
| ✅ C-004 | T q T⁻¹ = q − Q | `BosonicExpSimplify[Exp[-I Q p] ** q ** Exp[I Q p], v]` |
| ✅ C-061 | d/dt T(f_t) | ε-coefficient of the merged exponent |
| ✅ G-001 | [q, p] = iħ | `BosonicNormalOrder` (h, w are scalars automatically; no `"Scalars"`) |
| ✅ C-046 | two-mode su(1,1) Casimir = ((N₁−N₂)² − 1)/4 | `BosonicNormalOrder[…, v2]` |
| ✅ G-068 | S₂† a S₂ = a cosh r − b† sinh r | `BosonicExpSimplify[…, v2]` |
| ✅ G-071 | S₂,ab = S_c(ξ) S_d(−ξ) | merged exponent normal-ordered |
| ✅ G-075 | amplifier gives S₂(ηt) | same call as G-068 |
| ✅ G-095 | low-gain SPDC elements | `BosonicMatrixElement[{{1,1},{0,0}}, …, v2]` = −iηt |
| ✅ G-037 | σ algebra, [H_JC, N] = 0 | `BosonicNormalOrder[…, v3]` |
| ✅ G-046 | Schwinger J's and the J₃ phase action | `BosonicNormalOrder`, `BosonicExpSimplify[…, v2]` |
| ✅ C-024 | u, v from quadrature blocks | scalar expansion |
| ✅ C-018 | e^G a e^{−G} for G = x a² + y a†² + w a†a | `BosonicExpSimplify[…, v]` |
| ✅ P2-010 / G-022 | forced oscillator phase | `BosonicExpSimplify[D ** D, v]` |
| ✅ P2-012 | K₀, K₂ in quadratures | `BosonicNormalOrder[…, v]` |
| ✅ P2-042 | [a, H_d] = ε | `BosonicNormalOrder[…, v]` |
| ✅ B-075 | resolvent / coherent kernel | `BosonicMatrixElement` Fock and coherent |
| ✅ G-025 | ⟨β\|a†\|α⟩ = β̄⟨β\|α⟩ | `BosonicMatrixElement[…, v, "Basis" -> "Coherent"]` |
| ✅ P2-034 | thermal adjoint Lindbladian on a, a†a | `BosonicNormalOrder[…, v]` |
| ✅ P2-048 | Tr e^{−s(2N+1)} = 1/(2 sinh s) | `Sum` of `BosonicMatrixElement` |
| ✅ B-097 | retarded contraction | `BosonicNormalOrder[…, v]` |
| ✅ C-001 | [P, Q] = −i | `BosonicNormalOrder[p ** q - q ** p, v]` |
| ✅ C-057 | Segal [a(f), a†(g)] = Σ f̄g, 3 modes | `BosonicNormalOrder[…, v3]` |
| ✅ G-094 | broadband commutator, finite band | same call as C-057 |
| ✅ X-004 / P-025 / C-041 | SU(2) disentangling | `BosonicExpOrder[Exp[x bd1**b2 - y bd2**b1], v2]` |
| ✅ P-024 / C-038 | su(2) ladder element, symbolic j, m | `BosonicMatrixElement[{{j+m+1, j-m-1}, {j+m, j-m}}, bd1 ** b2, v2]` |
| ✅ P-026 / C-042 | SU(2) annihilator conjugation | `BosonicExpSimplify[…, v2]` |
| ✅ C-012 (QF part) | coherent Mehler kernel e^{wᾱβ}⟨α\|β⟩ | `BosonicMatrixElement[…, "Basis" -> "Coherent"]` |
| ✅ B-053 | passive transformation of normal symbols | **new with Idea 7** |

## 🟡 / ⬜ Needs development (his class b and c)

| Row | Formula | Status | Needs |
|---|---|---|---|
| ✅ X-003 / G-072 | two-mode squeezed vacuum | nullifier ✅; Fock series ✅: `BosonicExpOrder` disentangles S₂(re^{iθ}) into e^{e^{iθ}tanh r a₁†a₂†} sech r^{N₁+N₂+1} e^{−e^{−iθ}tanh r a₁a₂}, and the multimode wings give ⟨k,k\|S₂\|0,0⟩ = e^{ikθ} tanh^k r sech r, k ≤ 3, off-diagonal 0 (`Assumptions -> r > 0 && θ ∈ Reals`) | — |
| ⬜ C-015 | symplectic covariance of Weyl symbols | Idea 7 done | 3 |
| ✅ X-006 / P-036 / B-048 | multimode Gaussian vacuum, det(I − ZZ†)^{1/4} | via Idea 1 (multimode): for U = e^{−itH}, H = N₁ + 2N₂ + (½ + i/5)a₁†a₂† + h.c., the ordered U\|0⟩ = s e^{a†Za†/2}\|0⟩ gives \|s\|⁴ = det(I − ZZ†) exactly for real t, Z ≠ 0. Z here has zero diagonal; a generator with a beam splitter too (all of Z nonzero) times out at symbolic t | — |
| ✅ P2-024 | multimode pair Riccati | via Idea 1 (multimode): for H = N₁ + N₂/2 + αa₁†a₂† + ᾱa₁a₂, α = 1/3 + i/4, the a₁†a₂† coefficient z of the ordered e^{−itH} solves iż = α + ᾱz² + (1 + ½)z (the atlas iŻ = A + ZA†Z + CZ + ZDᵀ for one pair) and the vacuum scalar s solves ṡ/s = −iᾱz, symbolic t | — |
| ⬜ C-019 | Gaussian propagation by linear-fractional map | Idea 1 done for one mode, but the row is a position-space wavefunction and QF has no exact position representation | 1 |
| ✅ C-029 | squeeze-rotation factorization | via Idea 1: U_t = e^{−itH}, H = (αq² + βp² + γ(qp + pq))/2, ordered by `BosonicExpOrder`; ⟨0\|U_t\|0⟩ = ϑ_t⟨0\|S(B_t)\|0⟩ = ū_t^{−1/2} with u_t = (A + D + i(C − B))/2 from the classical flow e^{t[[γ, β], [−α, −γ]]}, checked as ⟨0\|U_t\|0⟩⁻² = ū_t for symbolic α, β, γ, t, stable and unstable | — |
| ✅ C-032 | driven quadratic propagator | via Idea 1 (affine; ~~the audit's `BosonicExpOrder[Exp[-I t (ad**a + 1/2 + f (a + ad)/Sqrt[2])]]` returned unfactored~~): Fock elements m, n ≤ 2 equal those of e^{−it/2 + itf²/2} D(−f/√2) e^{−itN} D(f/√2), symbolic f and all real t (without squeeze terms the ordered form has no root) | — |
| ✅ C-050 | time-dependent SU(1,1) Hamiltonian | via Idea 1, constant coefficients as in the audit: H = λK₊ + λ̄K₋ + μK₀ (λ̄ an independent symbol); from the ordered e^{−itH}, ω = 2·[a†² coefficient] solves iω̇ = λ̄ω² + μω + λ, ω(0) = 0, and γ, ω̄ read from the a² and N factors solve γ̇ = −λω̄ − λ̄ω − μ, symbolically. ~~Known limit: with squeeze terms the scalar D^{−1/2} is the principal root, right until D first crosses the negative real axis (oscillating evolutions past half a period)~~ Fixed 2026-11-09 after Mads's long-time check: the root is written e^{−il/2}(e^{−il}D)^{−1/2} ~~with the branch l = iw√(1 − 4pq/w²), continuous for every Hermitian generator at real time and every real squeeze~~ (that l is 0·∞ at w = 0, so ω → 0 of a symbolic ω gave `Indeterminate`, Mads 2026-11-10); since 2026-11-10 l comes from `rootPhase`: continuous when the coefficients are numeric times t, ~~principal root when they are symbolic (regular at ω = 0, right up to the first crossing)~~ and since 2026-11-11 (Mads: symbolic coefficients substituted after ordering flipped sign past the first crossing) for symbolic coefficients l = i sign(Re iw)√(4pq − w²) in the oscillating regime (4pq − w² real positive) and 0 elsewhere, so continuous in t for both rotation senses, principal in the hyperbolic regime and at ω = 0 (regular by substitution and by `Limit`); checks in `Downloads\RootBranchTests\RootBranches.wlt`. Probe `root continuous` checks one limit from both sides of the first crossing; the Riccati check here does not depend on the branch | — |
| ✅ C-022 | coherent elements of a metaplectic operator | via Ideas 1 and 2: for H = (αq² + βp²)/2, the coherent element of the ordered e^{−itH} (labels α = (x + iξ)/√2 real, ᾱ written out) equals 2^{1/2}Y^{−1/2}e^{−iσ(Fz′,z)/2}e^{−(x²+ξ²)/4}e^{−(ξ+ix)²W/4}: label derivatives of the logs agree, prefactor² = 1, ratio → 1 as t → 0, oscillating and squeezing flows. σ(z, z′) = x·ξ′ − ξ·x′ is the atlas's own (from its ⟨φ_z, φ_z′⟩ = e^{iσ(z,z′)/2}e^{−\|z−z′\|²/4}); the opposite sign fails | — |
| ✅ C-054 | Herman–Kluk for quadratics | via Ideas 1 and 2: (2π)⁻¹∫ u^{1/2}⟨α_w\|φ_{Fz}⟩⟨φ_z\|α_{w′}⟩dz, a 2-D Gaussian in closed form (2π e^{½bᵀM⁻¹b}/√det M), equals the same coherent element, same checks | — |
| ⬜ C-030 / C-053 | Weyl symbol of a metaplectic operator | | 1, 3 |
| ✅ B-026 | kernel acting on a vector | Bargmann kernel of D(h) from the coherent element, integrated against Φ = \|1⟩ with weight e^{−ᾱα} (Gaussian formula with a source sᾱ, then ∂_s; A = 1), gives (z̄ − h̄) e^{−\|h\|²/2} e^{z̄h} for symbolic z, h | — |
| ✅ B-027 | product of kernels | the weighted integral of K_{D(h)} K_{D(k)} (Gaussian formula, A = 1) equals the kernel of D(h) ** D(k) for symbolic z, w, h, k | — |
| ✅ P-020 / G-029 | normally ordered moments from χ_N | χ_N = ⟨ψ\|e^{γa†} e^{−γ̄a}\|ψ⟩; (−1)^n ∂_γ^m ∂_γ̄^n χ_N at 0 gives ⟨a†^m a^n⟩ for \|β⟩ and \|2⟩, m, n ≤ 2. Class-d line ρ = ∫d²γ/π χ_W(γ) D(−γ) ✅ too: on coherent elements ⟨z\|·\|w⟩ it rebuilds \|β⟩⟨β\| (`Integrate`, symbolic β) and \|m⟩⟨m\|, m ≤ 3 (χ_W = e^{−\|γ\|²/2} L_m from `BosonicMatrixElement`; polynomial × Gaussian by source derivatives of the Gaussian formula) | — |
| ✅ G-027 | P(α) = δ²(α−β) for \|β⟩, P_m = e^{\|α\|²}/m! ∂_α^m ∂_ᾱ^m δ²(α) for \|m⟩ | moved from class d. \|β⟩: `FourierTransform` of the coherent χ_N gives `DiracDelta[ax - bx] DiracDelta[ay - by]`. \|m⟩, weakly: ∫d²α P_m e^{γᾱ−γ̄α} (`Integrate` with `DiracDelta` derivatives, ∂_α∂_ᾱ = Δ/4) equals the Fock χ_N, m ≤ 3 | — |
| ✅ C-031 | translation coefficient of a metaplectic operator | oscillator case: `BosonicMatrixElement` reads ⟨α\|D(−α₀) e^{−itN}\|α⟩ directly (since `BosonicExpOrder` moves e^{ᾱa} past the N factor); the coherent trace (Gaussian ∫d²α/π e^{−A\|α\|²+Bα+Cᾱ} = e^{BC/A}/A, Re A = 1 − cos t > 0) gives e^{i\|α\|² cot(t/2)/2}/(2i sin(t/2)) for symbolic α, t; general Cayley-chart formula not checked | — |
| ✅ G-024 | Fock elements from the coherent diagonal symbol | ∂_ᾱ^m ∂_α^n (e^{\|α\|²}⟨α\|F\|α⟩)/√(m!n!) at 0 equals ⟨m\|F\|n⟩ for F = D(g) and e^{cN}, m, n ≤ 3 | — |
| ✅ P-017 | Tr D(γ) = πδ²(γ), Q and P symbols of D(γ), Tr[D(α)D(β)†] = πδ²(α−β), Â = ∫d²γ/π Tr[D(γ)Â] D(−γ) | moved from class d. Trace as ∫d²β/π ⟨β\|·\|β⟩: the coherent element is a pure phase in β, so `FourierTransform` at frequency 0 gives `Pi DiracDelta[p] DiracDelta[q]`; Tr[D(α)D(β)†] gives πδ(p₁−p₂)δ(q₁−q₂) with its phase 1 on the support. Q = ⟨α\|D\|α⟩; P read from the antinormal D and checked through ∫d²α/π P(α)⟨z\|α⟩⟨α\|w⟩ = ⟨z\|D\|w⟩. Reconstruction checked on coherent elements for Â = a, a†, a a†, with Tr[D Â] read directly (~~through Tr[a† D], Tr[a† D a]~~ until the D ** a† gap was fixed) | — |
| ✅ B-020 | rank-one projector symbols | \|1⟩⟨1\| = a† :e^{−a†a}: a: Fock elements δ_{m1}δ_{n1} and coherent ᾱβ⟨α\|0⟩⟨0\|β⟩ checked (`f76ab7da`) | — |
| ⬜ C-016 | Moyal product, exact for deg ≤ 2 | | 3 |
| ⬜ C-053 | Weyl symbol of the oscillator exponential | | 3 |
| ⬜ C-049 | Wigner function of a squeezed Gaussian | Idea 2 done | 3 |
| 🟡 X-007 / B-107 | Tr_B Γ(M) = det(I − M)⁻¹ | diagonal M ✅: Σ ⟨k₁,k₂\|m₁^{N₁} m₂^{N₂}\|k₁,k₂⟩ = det(I − diag(m₁, m₂))⁻¹ for symbolic 0 < mᵢ < 1; general M ⬜: Idea 1 (multimode) now orders Γ(M) = e^{a† log M a} as :e^{a†(M − 1)a}:, but `BosonicMatrixElement` does not read a mode-mixing normal-ordered middle, so the trace waits for the normal symbol | ~~1 + 4 (general)~~ normal symbol |
| ✅ B-106 | grand partition function, free limit | Σ over two modes of e^{−β Σ(εᵢ−μ)Nᵢ} = Π 1/(1 − e^{−β(εᵢ−μ)}), symbolic β, εᵢ, μ | — |
| ✅ G-040 | dispersive conditional rotation | e^{−iχt(N_e + N_a(N_e − N_g))}: e^{−iχt(n+1)} on the excited branch, e^{iχtn} on the ground branch, 0 between, symbolic n | — |
| ✅ G-012 | Susskind–Glogower E E† = 1, E†E = 1 − \|0⟩⟨0\| | E†E ✅ (diagonal {0, 1, 1, …}, symbolic `Piecewise[{{1, n >= 1}}, 0]`); E E† = 1 by `BosonicExpSimplify` (`1d307d37`) | — |
| ✅ G-038 | resonant Jaynes–Cummings propagator | his 2×2 block U (cos λt√(N+1), −i a sin λt√N/√N, −i a† sin λt√(N+1)/√(N+1), cos λt√N): `BosonicExpSimplify` gives H U entries −i√(N+1) sin λt√(N+1), cos λt√(N+1) a, a† cos λt√(N+1), −i√N sin λt√N; H U − i dU/dt has zero Fock elements m, n ≤ 3; U(0) = 1 | — |
| ✅ G-041 | Raman number-controlled rotation | U = cos λtN − i sin λtN σ₁ with H = λN σ₁: H U − i dU/dt zero block, U(0) = 1 | — |
| 🟡 P-026 | SU(2) coherent state expansion | fixed j ✅; symbolic j ⬜ | 4 |
| 🟡 P-030 / C-047 | SU(1,1) coherent expansion and overlap | k = 1/4 annihilator `(b**b)/2` ✅ (Idea 7); class-d measure (2k−1)/π ∫_{\|ζ\|<1} d²ζ/(1−\|ζ\|²)² \|ζ;k⟩⟨ζ;k\| = 1 ✅ at k = 3/4, the odd Fock sector of K₊ = a†²/2 (the only oscillator sector with k > 1/2): exact polar `Integrate` of ⟨2m+1\|e^{ζa†²/2}\|1⟩ products gives `IdentityMatrix[3]`, m, n ≤ 2; general k ⬜ | 6 |
| ⬜ P-028 / C-045 | su(1,1) with label k | | 6 |
| ⬜ P2-013 | SU(1,1) Riccati equation and phase | | 6 |

## Extra gaps from the audit

From the "QF problems found (exact layer)" section of `Operator-Exponentials-QF-Audit.md`, which points to the design note for the fixes. Rechecked with `audit-gaps.wls` in this folder.

| Gap | Status | How |
|---|---|---|
| ✅ ⟨0\|0^(−n)\|0⟩ and ⟨1\|0^(n−1)\|1⟩ were `Indeterminate` | fixed `1a4ca028` | `numberDiagonal`: zero base with an n-dependent exponent → `Piecewise[{{1, p == 0}}, 0^p]`; ⟨1\|0^(−n)\|1⟩ still `ComplexInfinity` |
| ✅ list as first argument, default (Weyl) method | worked already | his bug was the Gröbner path, removed here |
| ✅ list with `Method -> "Blasiak"` | fixed `1a4ca028` | list rule in `MultiModeBlasiakOrder` |
| ✅ list in `BosonicVEV` | fixed `1a4ca028` | list rule in `BosonicVEV` |
| ✅ `BosonicExpOrder` of e^{−(w+μ) n} (sum as coefficient) | fixed `1a4ca028` | `canonicalizeGenerators` collects the exponent by monomial (`exponentTerm`) instead of fully expanding it |
| ✅ `BosonicExpSimplify` of S(ξ) ** S(−ξ) left as `E^(X/2 − X/2)` | fixed `1a4ca028` | merged exponent is expanded |
| ✅ `BosonicExpSimplify` of S (a + f) S⁻¹ | fixed by Idea 7 | |
| ✅ `BosonicExpSimplify` with a mixed quadratic + linear generator | works | gives b − 2x b† − y |
| ✅ `BosonicExpOrder` of e^{c n + linear}, general quadratic + linear | ~~open~~ fixed by Idea 1 (single mode) | e^{wN + μa† + νa} matches e^{μν(eʷ−1−w)/w²} e^{μ(eʷ−1)/w a†} e^{wN} e^{ν(eʷ−1)/w a} on Fock elements (probe) |
| ✅ two-mode general complex-symmetric squeeze in `BosonicExpSimplify` | ~~open~~ fixed 2026-11-10 | ~~Idea 1 (multimode) orders it, but `BosonicExpSimplify` of S b₁ S⁻¹ still returns unchanged (rechecked 2026-11-10, symbolic and numeric z): the conjugation rules do not use the multimode matrix~~ `adRotate` only sums the series when the commutators of op close after two steps; when they do not, `linearConjugate` conjugates a linear op by any B of degree ≤ 2 over any modes as (a†, a)·e^{M}ℓ plus the constant from B's linear part, M = {{W, −P}, {Q, −Wᵀ}}; for W = 0 and no linear part, e^{M} is written as cosh√(−PQ) and sinh√(−PQ)/√(−PQ) blocks (`bogoliubovExp`). Sums, products and exponentials of linear ops follow through `conjugateParts`. Checked: S a S⁻¹ and S b† S⁻¹ for symbolic complex-symmetric Z, single-mode plus pair, a B with hopping and linear terms on a + 2b† + 3, a singular −PQ (squeeze on a, op a + b), all equal the commutator series to t³; numeric Z against 25-level matrices; G-068 output unchanged |
| ✅ linear exponential over several modes, e^{γa₁† − γ̄a₁ + ηa₂† − η̄a₂} | ~~open (found 2026-11-08)~~ fixed by Idea 1 (multimode) | ~~`BosonicExpOrder` returns it unchanged and `BosonicVEV` stays unevaluated~~ uncoupled modes split into single-mode factors, e^{−\|γ\|²/2−\|η\|²/2} e^{γa₁†} e^{ηa₂†} e^{−γ̄a₁} e^{−η̄a₂}; `BosonicVEV` and ⟨1,0\|·\|0,0⟩ equal their closed forms |
| ✅ power of a sum above 30, (a + a†)³¹ | fixed 2026-11-08, uncommitted | `NonCommutativeExpand` refuses exponents above its `"ExpandMaxExponent"` (30), so `BosonicNormalOrder` came back unevaluated and `BosonicVEV` errored with private symbols in the output. `ladderOrder` (BosonicAlgebra.m) now builds such powers one factor at a time, normal-ordering after each product (exponent 60 in about 0.5 s, (a₁ + a₁† + a₂ + a₂†)³² in about 14 s); ⟨0\|(a + a†)ⁿ\|0⟩ = (n−1)!! for n = 28 to 34, ⟨1\|(xa + ya†)³³\|0⟩ = 33!! x¹⁶y¹⁷, nested and two-mode powers. Two tests in BosonicMatrixElement.wlt |
| ✅ `BosonicVEV` of a product of exponentials | fixed | ~~`BosonicVEV` itself still unevaluated~~; ~~a second `BosonicVEV` definition routes non-polynomials~~ `BosonicVEV` is now one definition: `BosonicMatrixElement` between vacuum Fock states, any number of modes, passing `Assumptions` (`Method` dropped) |
| ✅ ⟨m\|D(α)\|m⟩, ⟨2m\|S(ξ)\|0⟩ with symbolic m | ~~open~~ fixed | ~~not small: closed forms need a `Piecewise` on the sign of m − n~~ the D and S closed forms take symbolic indices: m − n picks the branch when it is a number (⟨n\|D\|n⟩ = e^{−\|α\|²/2}L_n(\|α\|²), ⟨n+2\|S\|n⟩, odd differences 0), `Piecewise` otherwise; ⟨2m\|S\|0⟩ = √((2m)! sech r)(−e^{iθ}tanh r)^m/(2^m m!) under m ∈ ℤ, m ≥ 0 |
| ➖ Kerr e^{λ n²} in `BosonicExpOrder` | out of reach | not Gaussian; no finite disentangling |
| ✅ `BosonicExpOrder` of D(−α) ** e^{−itN} | fixed (found checking C-031) | ~~returned e^{−αa†} ** e^{ᾱa} ** :e^{(e^{−it}−1)N}:, not normal-ordered~~. Two rules per ordering in `modeRulesNO`/`modeRulesAN` move a ladder exponential past the number factor, from e^{−sN} a e^{sN} = eˢ a: normal e^{βa} :e^{cN}: = :e^{cN}: e^{β(1+c)a}, antinormal e^{βa†} ⋮e^{caa†}⋮ = ⋮e^{caa†}⋮ e^{β(1−c)a†}, and the mirrored pairs; checked against hand-written Fock elements |
| ⬜ symbolic Fock index with both wings, ⟨n\|e^{γa†} e^{−γ̄a}\|n⟩ | open (found checking P-020) | unevaluated: `normalProductElement` sums the middle only for integer m, n (BosonicAlgebra.m:324); the closed form is L_n(γγ̄). Polynomials already take symbolic n (⟨n\|a†² a²\|n⟩ = n(n−1)) |
| ✅ coherent element of D(γ) ** a† | fixed (found checking P-017) | ~~`BosonicMatrixElement[{β, β}, disp ** bd, v, "Basis" -> "Coherent"]` unevaluated~~. `BosonicExpOrder` already gave e^{γa†} ** (a† − γ̄) ** e^{−γ̄a}; the element could not read the sum factor. `sumProductElement` distributes `**` over sums, pulls scalars out, and adds the elements of the terms, so D ** a†, D ** a ** a†, e^{ca} ** a† and sums like e^{ca} + a† now evaluate, Fock and coherent |

## ➖ Out of scope (his class d, ~~68~~ ~~67~~ 66 rows)

P-017 and G-027 moved to ✅ (2026-11-05): their distributions come out of `FourierTransform` of coherent elements, or weakly from `Integrate` with `DiracDelta`. The class-d parts of P-020 (ρ from χ_W) and P-030 (measure, k = 3/4) are ✅ in their rows. Tried and left out: C-014 (needs Weyl symbols, Idea 3), C-097 (disk symbol calculus with free label, Idea 6). The rest unchanged from his note: scalar or matrix identities plain WL already does, infinite-mode and continuum statements, distributions, inequalities, impossibility and existence statements, spectral classifications and cited theorems, half-line differential operators, position-space kernels, untranscribed definitions, phase-space integrals, approximations. See the last section of `Minimal-Bosonic-Design.md` for the row lists.

## Coverage by component (2026-11-08)

Each in-scope row (his class a, b and c, plus P-017 and G-027 moved from class d) checked on paper against five components:

- **Ordering**: polynomials in the ladder operators brought to normal order by `NonCommutativeExpand` over the declared algebra; this includes su(1,1) with a free label k (commutation relations, the Casimir set to k(k−1), the lowest-weight module as a left ideal) and N declared as a generator with N = a†a.
- **Exponential rule**: Idea 1 written on K₊, K₀, K₋ (2N×2N matrix, affine with linear terms). It both disentangles e^{G} and conjugates linear operators by e^{G}, since conjugation uses the same matrix; conjugation of a sum or a scalar multiple is Idea 7.
- **Shift rule**: Idea 5's a^q f(N) = f(N+q) a^q, also for f an exponential (moving e^{βa} past e^{cN}).
- **Elements**: reading Fock elements (Taylor coefficients, n in place of N) and coherent elements (β for a, ᾱ for a†, times the overlap) from the ordered form; Ideas 2 and 4 and the Fock part of Idea 5.
- **s-ordered symbol**: Idea 3, the only component with no code yet.

Only C-019 fits none: it is a position-space wavefunction. This is a paper check; it has not been run.

| Row | Ordering | Exponential rule | Shift rule | Elements | s-ordered symbol |
|---|---|---|---|---|---|
| P-014 / P2-050 | | ✓ | | | |
| P-006 | | ✓ | | | |
| C-004 | | ✓ | | | |
| C-061 | | ✓ | | | |
| G-001 | ✓ | | | | |
| C-046 | ✓ | | | | |
| G-068 | | ✓ | | | |
| G-071 | ✓ | ✓ | | | |
| G-095 | ✓ | | | ✓ | |
| G-037 | ✓ | | | | |
| G-046 | ✓ | ✓ | | | |
| C-024 | ✓ | | | | |
| C-018 | | ✓ | | | |
| P2-010 / G-022 | | ✓ | | | |
| P2-012 | ✓ | | | | |
| P2-042 | ✓ | | | | |
| B-075 | | | | ✓ | |
| G-025 | | | | ✓ | |
| P2-034 | ✓ | | | | |
| P2-048 | | | | ✓ | |
| B-097 | ✓ | | | | |
| C-001 | ✓ | | | | |
| C-057 | ✓ | | | | |
| G-094 | ✓ | | | | |
| X-004 / P-025 / C-041 | | ✓ | | | |
| P-024 / C-038 | ✓ | | | ✓ | |
| P-026 / C-042 | | ✓ | | ✓ (symbolic j open) | |
| C-012 | | | | ✓ | |
| B-053 | | ✓ | | | |
| X-003 / G-072 | | ✓ | | ✓ | |
| C-015 | | ✓ | | | ✓ |
| X-006 / P-036 / B-048 | | ✓ | | ✓ | |
| P2-024 | | ✓ | | | |
| C-019 | | | | | |
| C-029 | | ✓ | | ✓ | |
| C-032 | | ✓ | | ✓ | |
| C-050 | | ✓ | | | |
| C-022 | | ✓ | | ✓ | |
| C-054 | | ✓ | | ✓ | |
| C-030 | | ✓ | | | ✓ |
| C-053 | | ✓ | | | ✓ |
| B-026 | | | | ✓ | |
| B-027 | | | | ✓ | |
| P-020 / G-029 | | | | ✓ | |
| G-027 | | | | ✓ | |
| C-031 | | | ✓ | ✓ | |
| G-024 | | | | ✓ | |
| P-017 | | | | ✓ | |
| B-020 | | | | ✓ | |
| C-016 | ✓ | | | | ✓ |
| C-049 | | | | ✓ | ✓ |
| X-007 / B-107 | | ✓ (general M) | | ✓ | |
| B-106 | | | | ✓ | |
| G-040 | | | | ✓ | |
| G-012 | | | ✓ | ✓ | |
| G-038 | ✓ | | ✓ | | |
| G-041 | ✓ | | ✓ | | |
| P-030 / C-047 | ✓ | | | | |
| P-028 / C-045 | ✓ | | | | |
| P2-013 | | ✓ | | | |

For P-028, P-030 and P2-013 no bosonic rule is used: conjugation by e^{ζK₊} ends after three terms, the coherent-state norms come from the left ideal, and the exponential rule needs only the 2×2 group matrix, not k. The claim that the rule written on K₊, K₀, K₋ reproduces the bosonic results is the weakest premise here; it is checked by the exact differential operators a = ∂_z, a† = z, N = z∂_z, K₊ = z, K₀ = z∂_z + k, K₋ = z∂_z² + 2k∂_z acting on z^m.

## Log

- 2026-11-03: tracker created; 31 probe checks pass at `ffeb814d` plus the uncommitted Idea 7 edit; Ideas 1–5 confirmed missing.
- 2026-11-03: `1a4ca028` commits Idea 7 and five audit gaps (vacuum powers, Blasiak and VEV lists, summed exponent coefficient, S(ξ)S(−ξ)); "Extra gaps from the audit" section added.
- ~~Test files still using old call forms (not yet converted): 9 lines across BosonicExpOrder.wlt, BosonicExpOrderAntinormal.wlt, BosonicExpBraid.wlt; SecondQuantization.wlt (12 tests).~~ Superseded by the entry below.
- 2026-11-04: `62a54595` commits Idea 2 (coherent elements of ordered products) and the annihilator shorthand: vars may be `{a, b}`, creators `SuperDagger` of each name; `{{c…}, {a…}}` and `WeylAlgebra` unchanged. The shorthand makes the old `ToBosonicOperator[…, {𝕒, 𝕓}, s]` test lines valid again.
- 2026-11-04: `f76ab7da` commits the Fock/coherent part of Idea 5 (bare ladder powers beside f(N)); B-020 ✅, G-012 🟡 (E†E done, E E† needs the shift rule). Probe now 35 checks.
- ~~Test lines still in the old form: 2 in BosonicExpOrder.wlt (`First @ FieldVariables[{1}]`, `BosonicMatrixElement` without vars); SecondQuantization.wlt (12 tests).~~ Converted in `9ea7d75f`; all 274 SecondQuantization tests pass. Usage strings and docs do not yet mention the shorthand.
- 2026-11-04: `1d307d37` commits the shift rule (second half of Idea 5); Idea 5 ✅, G-012 ✅, G-038 and G-041 🟡 (idea done, rows not checked). 274 SecondQuantization tests pass. Probe now 43 checks, all pass: shift identities, 40-level matrix check, ⟨1\|a f(N)\|2⟩ = √(2/3), E E† and E² E†² = 1, E†E kept, f-deformed, Holstein–Primakoff, `b ** b†` untouched. Open section adds the direct `BosonicMatrixElement` of `b ** f(N)` (needs `BosonicExpSimplify` first) and a second mode between the word and f(N).
- 2026-11-04: G-038 (resonant Jaynes–Cummings) and G-041 (Raman) checked on `1d307d37`, both ✅; all four Idea 5 rows done. Probe now 45 checks, all pass. Seen while checking: a scalar inside a ladder factor, `(λ b) ** f(N)`, blocks the rules; write `λ (b ** f(N))`.
- 2026-11-05: B-026, B-027, P-020/G-029, C-031, G-024 checked, all ✅ with no kernel change; Idea 2's rows done. Probe now 51 checks, all pass, all exact. New open gap: `BosonicExpOrder` leaves e^{ᾱa} left of a function of N.
- 2026-11-05: first part of Idea 4: `BosonicMatrixElement` reads a multimode function of the number operators in the Fock basis, four lines in the dispatcher. G-040 and B-106 ✅, X-007 🟡 (diagonal M). 276 SecondQuantization tests pass (2 new in BosonicMatrixElement.wlt). Probe now 54 checks, all pass.
- 2026-11-05: `BosonicExpOrder` moves ladder exponentials past a number exponential, normal and antinormal; the D(−α) e^{−itN} gap is fixed and C-031 is read without the manual braid. 279 SecondQuantization tests pass (3 new). Probe now 55 checks, all pass.
- 2026-11-05: P-017 moved from ➖ to ✅ with no kernel change: traces of D(γ) as `FourierTransform` of coherent elements give the deltas exactly; Q, P symbols and the reconstruction for a, a†, a a† checked. Probe now 59 checks, all pass. New open gap: coherent element of D ** a†.
- 2026-11-05: G-027 moved from ➖ to ✅; class-d parts of P-020 (ρ from χ_W) and P-030 (measure at k = 3/4) ✅. No kernel change. Probe now 64 checks, all pass, all exact.
- 2026-11-05: D ** a† gap fixed: `BosonicMatrixElement` is linear over sums in the exponential branch (`sumProductElement`), so shifted ladder factors left by `BosonicExpOrder` are read. P-017's reconstruction now uses D ** a† directly. 283 SecondQuantization tests pass (2 new in BosonicMatrixElement.wlt). Probe now 65 checks, all pass.
- 2026-11-06: Idea 4 done: multimode wings in `BosonicMatrixElement` (`multimodeProductElement`, sharing the `BosonicExpOrder` retry with the single-mode path), X-003 ✅. `BosonicVEV` routes exponentials to the vacuum element. 286 SecondQuantization tests pass (3 new; `BME-Modes-ExponentialUnevaluated` became `BME-Modes-ExponentialOneMode` since the element now evaluates). Probe now 67 checks, all pass.
- 2026-11-06: cleanup: `BosonicVEV` is one definition through `BosonicMatrixElement`; one `dequantize` helper (Utils.m) replaces `dequantizeLadder`, `vevCollapse` and the hand-written collapses. Idea 1 🟡: single-mode quadratic disentangling and braid (see the Ideas table). 288 SecondQuantization tests pass (2 new). Probe: `Idea1 quad3` moved out of the open section, 68 checks, all pass.
- 2026-11-06: Idea 1 affine step: quadratic plus linear in one mode, ordering and braid; the e^{c n + linear} gap is fixed. `BosonicExpBraid-MixedFamilyIsInert` removed (number + linear now braids; covered by `BosonicExpBraid-QuadraticPastLadder`). 287 SecondQuantization tests pass. Probe now 69 checks, all pass.
- 2026-11-06: C-029, C-032, C-050 ✅ via Idea 1, no further kernel change except that a number-plus-linear exponent (no squeeze terms) is now written with e^{wN} exactly, so C-032 holds for all real t. C-019 left ⬜: a position-space statement, and QF has no exact position representation. Probe now 72 checks, all pass.
- 2026-11-06: C-022 and C-054 ✅ via Ideas 1 and 2, no kernel change: coherent elements of the ordered quadratic propagator against the closed Gaussian and the Herman–Kluk integral (closed-form 2-D Gaussian), oscillating and squeezing. Probe now 74 checks, all pass.
- 2026-11-08: "Coverage by component" section added (Mads): every in-scope row placed against ordering, exponential rule, shift rule, elements and s-ordered symbol; only C-019 fits none. Paper check, not run.
- 2026-11-08: exponential rule written on K₊, K₀, K₋ checked (`exponential-rule-check.wls`, `exponential-rule-vs-qf.wls`): with K₋ = d/dz, K₀ = z d/dz + k, K₊ = z² d/dz + 2kz, the product e^{pK₊}e^{qK₀}e^{rK₋} read from the 2×2 matrix solves ∂ₜu = Gu exactly for symbolic k, coefficients and f; at k = 1/4 it equals QF's ordered quadratic exponential (Fock m, n ≤ 3, three points). Its one flaw is the scalar g₂₂^{−2k} on the principal branch: past the time g₂₂ crosses the negative real axis it is off by e^{4πik} for every k (−1 at k = 1/4, the known sign bug), so the fix belongs in the rule, not in its bosonic use.
- 2026-11-08: power-above-30 fix in `ladderOrder` (uncommitted), passed /wl-verify in two rounds (round 1 replaced repeated squaring, too slow, by one factor at a time); all SecondQuantization tests pass with 2 new; probe unchanged (the 3 open checks still open). Idea 3 prototype kept outside QF in `s-ordered-symbol.wls`.
- 2026-11-10: Idea 1 multimode 🟡 (pure quadratic, committed `22f6ceb2`, `4fa82ef7`); X-006 / P-036 / B-048 and P2-024 ✅, X-007 general M now waits only on the normal symbol. Probe adds `P2-024 pair Riccati` and `X-006 Gaussian vacuum`; now 77 checks, all pass (the 3 open checks still open).
- 2026-11-10: Idea 1 ✅: linear terms over coupled modes (`affineShift`); the several-mode linear exponential gap is closed. The two-mode complex-symmetric squeeze in `BosonicExpSimplify` stays open. Probe adds `Idea1 multimode linear`; now 78 checks, all pass.
- 2026-11-10: the two-mode complex-symmetric squeeze gap in `BosonicExpSimplify` is closed (`linearConjugate`, `bogoliubovExp`).
