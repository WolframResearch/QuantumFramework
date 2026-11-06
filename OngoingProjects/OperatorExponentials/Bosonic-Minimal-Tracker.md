# Bosonic minimal design — status

Tracks Mads's `Minimal-Bosonic-Design.md` in this folder. Every ✅ below was checked with the declared-vars API by `tracker-probe.wls` (run with `wolframscript -file tracker-probe.wls`); each check compares against the formula in Mads's table, not just "it evaluated".

Legend: ✅ works · 🟡 partly works · ⬜ needs development (idea number) · ➖ out of scope for an exact operator layer (his class d)

## Summary

| | Rows |
|---|---|
| ✅ works now | 44 (his 29 class a, plus B-053 via Idea 7, B-020 via Ideas 2 and 5, G-012, G-038, G-041 via Idea 5, B-026, B-027, P-020, C-031, G-024 via Idea 2, G-040, B-106, X-003 via Idea 4, P-017 and G-027 from class d via Idea 2) |
| 🟡 partly | 3 (P-026, P-030, X-007: one part works) |
| ⬜ needs an idea | 16 |
| ➖ out of scope | 66 |

## Ideas

| Idea | What | Status | Rows it unlocks |
|---|---|---|---|
| 1 | one rule for `Exp` of any quadratic (2N×2N matrix, Gauss factors) | ⬜ — `Exp[x b**b + y bd**bd + w bd**b]` returned unchanged | 10 |
| 2 | coherent elements of ordered products, ⟨α\|e^{P(a†)} f(N) e^{Q(a)}\|β⟩ | ✅ done 2026-11-04, committed `62a54595` (`normalProductElement` with a coherent flag; checked symbolically against D(g)\|β⟩ = e^{(gβ̄−ḡβ)/2}\|β+g⟩ and numerically against 60-level matrices for D, normal and antinormal products, squeeze, wing–f(N)–wing) | 9 |
| 3 | s-ordered (Weyl) symbol, `NormalOrdered[expr, vars]["Symbol", s]` | ⬜ — no such property | 5 |
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
| ⬜ X-006 / P-036 / B-048 | multimode Gaussian vacuum, det(I − ZZ†)^{1/4} | | 1 |
| ⬜ P2-024 | multimode pair Riccati | | 1 |
| ⬜ C-019 | Gaussian propagation by linear-fractional map | | 1 |
| ⬜ C-029 | squeeze-rotation factorization | | 1 |
| ⬜ C-032 | driven quadratic propagator | | 1 (affine) |
| ⬜ C-050 | time-dependent SU(1,1) Hamiltonian | | 1 |
| ⬜ C-022 | coherent elements of a metaplectic operator | Idea 2 done | 1 |
| ⬜ C-054 | Herman–Kluk for quadratics | Idea 2 done | 1 |
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
| 🟡 X-007 / B-107 | Tr_B Γ(M) = det(I − M)⁻¹ | diagonal M ✅: Σ ⟨k₁,k₂\|m₁^{N₁} m₂^{N₂}\|k₁,k₂⟩ = det(I − diag(m₁, m₂))⁻¹ for symbolic 0 < mᵢ < 1; general M ⬜ | 1 + 4 (general) |
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
| ⬜ `BosonicExpOrder` of e^{c n + linear}, general quadratic + linear | open | Idea 1 (affine form) |
| ⬜ two-mode general complex-symmetric squeeze in `BosonicExpSimplify` | open | Idea 1 (multimode) |
| ✅ `BosonicVEV` of a product of exponentials | fixed | ~~`BosonicVEV` itself still unevaluated~~; a second `BosonicVEV` definition routes non-polynomials to `BosonicMatrixElement` between vacuum Fock states, any number of modes, passing `Assumptions` |
| ⬜ ⟨m\|D(α)\|m⟩, ⟨2m\|S(ξ)\|0⟩ with symbolic m | open | not small: closed forms need a `Piecewise` on the sign of m − n |
| ➖ Kerr e^{λ n²} in `BosonicExpOrder` | out of reach | not Gaussian; no finite disentangling |
| ✅ `BosonicExpOrder` of D(−α) ** e^{−itN} | fixed (found checking C-031) | ~~returned e^{−αa†} ** e^{ᾱa} ** :e^{(e^{−it}−1)N}:, not normal-ordered~~. Two rules per ordering in `modeRulesNO`/`modeRulesAN` move a ladder exponential past the number factor, from e^{−sN} a e^{sN} = eˢ a: normal e^{βa} :e^{cN}: = :e^{cN}: e^{β(1+c)a}, antinormal e^{βa†} ⋮e^{caa†}⋮ = ⋮e^{caa†}⋮ e^{β(1−c)a†}, and the mirrored pairs; checked against hand-written Fock elements |
| ⬜ symbolic Fock index with both wings, ⟨n\|e^{γa†} e^{−γ̄a}\|n⟩ | open (found checking P-020) | unevaluated: `normalProductElement` sums the middle only for integer m, n (BosonicAlgebra.m:324); the closed form is L_n(γγ̄). Polynomials already take symbolic n (⟨n\|a†² a²\|n⟩ = n(n−1)) |
| ✅ coherent element of D(γ) ** a† | fixed (found checking P-017) | ~~`BosonicMatrixElement[{β, β}, disp ** bd, v, "Basis" -> "Coherent"]` unevaluated~~. `BosonicExpOrder` already gave e^{γa†} ** (a† − γ̄) ** e^{−γ̄a}; the element could not read the sum factor. `sumProductElement` distributes `**` over sums, pulls scalars out, and adds the elements of the terms, so D ** a†, D ** a ** a†, e^{ca} ** a† and sums like e^{ca} + a† now evaluate, Fock and coherent |

## ➖ Out of scope (his class d, ~~68~~ ~~67~~ 66 rows)

P-017 and G-027 moved to ✅ (2026-11-05): their distributions come out of `FourierTransform` of coherent elements, or weakly from `Integrate` with `DiracDelta`. The class-d parts of P-020 (ρ from χ_W) and P-030 (measure, k = 3/4) are ✅ in their rows. Tried and left out: C-014 (needs Weyl symbols, Idea 3), C-097 (disk symbol calculus with free label, Idea 6). The rest unchanged from his note: scalar or matrix identities plain WL already does, infinite-mode and continuum statements, distributions, inequalities, impossibility and existence statements, spectral classifications and cited theorems, half-line differential operators, position-space kernels, untranscribed definitions, phase-space integrals, approximations. See the last section of `Minimal-Bosonic-Design.md` for the row lists.

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
