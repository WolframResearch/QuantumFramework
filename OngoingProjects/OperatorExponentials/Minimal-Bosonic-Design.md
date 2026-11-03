# Minimal additions to QF's exact bosonic layer for the missing atlas formulas

Source: [Operator-Exponentials-QF-Audit.md](Operator-Exponentials-QF-Audit.md). Kernel: Wolfram 15.0.1, QF loaded from this checkout at `fbebae87`; no QF source was edited. File references are to `QuantumFramework/Kernel/SecondQuantization/`.

## Which entries are in scope

The audit has 326 formula entries. 182 of them say "No QF code" (106 as a bare line, 76 as a `QF code:` line that opens with "No" or "none"); the audit's own header counts 175, the difference being entries whose code line names a closest tool but no call. Dropping 55 fermionic, Grassmann, Pfaffian, Clifford and Thirring entries leaves **127 bosonic or general entries**, all read in full. A further 94 entries were checked only with truncated matrices; they are outside this brief and not classified here.

Of the 127:

- **class a, 29**: computable today with the right call; the audit never tried the symbolic layer on them.
- **class b, 27**: one missing rule or branch in an existing function (a few need two such ideas).
- **class c, 3**: need su(1,1) with an arbitrary representation label, which no bosonic realization gives; `NonCommutativeAlgebra` handles them.
- **class d, 68**: nothing for a finite exact operator layer to compute: scalar or matrix identities that plain Wolfram Language already does, infinite-mode or distributional statements, inequalities, cited theorems, position-space kernels.

The 27 class b rows and 3 class c rows fall under seven ideas: one rule for any quadratic exponential (10 rows), coherent-basis elements of ordered products (9), the s-ordered symbol of an operator (5), several modes in the exponential element branch (4), functions of N beside bare ladder operators (4), su(1,1) with a free label through `NonCommutativeAlgebra` (3), and conjugation that is linear over sums (2). Rows needing two ideas are counted under both.

## Table

"Today" means the call returns the atlas result in the current kernel; the probe that shows it is named in brackets. Ideas 1 to 7 are defined in the next section.

| Formula | Atlas ID | Class | Where QF stops | Smallest fix |
|---|---|---|---|---|
| Duhamel derivative of e^{L(t)} | C-051 | d | n/a | none: an integral over s of products for an arbitrary L(t) |
| [D(ω1), D(ω2)] = 0 iff Im(ω1 ω̄2) ∈ πZ; lattice completeness | P-014, P2-050 | a | completeness is d | today: `BosonicExpSimplify` of D1 ** D2 and D2 ** D1 gives prefactors whose ratio is e^{2i Im(u v̄)} [probe1 "11"] |
| Displacement as e^{-iPQ/2} e^{iPq} e^{-iQp} | P-006 | a | position action is d | today: `BosonicExpSimplify[Exp[I P q] ** Exp[-I Q p]]` with q, p written in a, a† returns e^{iPQ/2} e^{i(Pq - Qp)} [probe1 "12"] |
| T(z) Q T(z)⁻¹ = Q - q, coordinate factorization | C-004 | a | position action is d | today: `BosonicExpSimplify[Exp[B] ** q ** Exp[-B]]` returns q - Q [probe1 "13"] |
| d/dt T(f_t) = T(f_t)[a†(ḟ) - a(ḟ) + i Im⟨f, ḟ⟩] | C-061 | a | | today: `BosonicExpSimplify[Exp[-(f ad - fc a)] ** Exp[(f+ε g) ad - (fc+ε gc) a]]`; the ε coefficient of the exponent is g a† - gc a + (fc g - f gc)/2 [probe1 "15"] |
| [q, p] = iħ, a = (ωq + ip)/√(2ħω), H = ħω(n + 1/2) | G-001 | a | | today: `BosonicNormalOrder` with `"Scalars" -> {h, w}` gives 0 for H and i(√(h/w)√(hw) - h) for the commutator, zero for positive h, w [probe1 "24"] |
| a† has no normalizable eigenvector | G-014 | d | n/a | none: a statement about every vector of Fock space |
| Two-mode squeezed vacuum, Fock series and nullifiers | X-003, G-072 | b | BosonicAlgebra.m:404 rejects multimode indices for exponentials | Idea 4; the nullifier part is today: `BosonicExpSimplify[(a1 + z ad2) ** Exp[-z ad1 ** ad2]]` gives `Exp[...] ** a1` [probe1 "47null"] |
| su(1,1) realizations and Casimirs | C-046 | a | | today: `BosonicNormalOrder` of the two-mode Casimir minus ((N1 - N2)² - 1)/4 is 0 [probe1 "52"] |
| S2† a S2 = μ a - ν b† | G-068 | a | | today: `BosonicExpSimplify[Exp[-B] ** a1 ** Exp[B]]` gives a1 cosh r - e^{iθ} sinh r a2† [probe1 "53"] |
| S2,ab(ξ) = S_c(ξ) S_d(-ξ) | G-071 | a | | today: `BosonicExpSimplify` merges the two commuting exponentials; `BosonicNormalOrder` of the merged exponent is x a1 a2 - y a1† a2† [probe2 "54exp"] |
| Pair coherent states | G-090 | d | n/a | none: a state defined by its Fock series; the normalization is a scalar Bessel sum |
| Low-gain down-conversion to second order | G-095 | a | | today: `BosonicMatrixElement[{{n,n},{0,0}}, 1 + X + X ** X/2]` gives 1 - |η|²t²/2, -iηt, -η²t² [probe1 "57"] |
| Passive transformation of normal symbols | B-053 | b | BosonicExp.m:410-416, conjugation of a sum fails | Idea 7; monomials work today [probe1 "58", "58sum"] |
| e^{tA} for a 2×2 su(2) matrix | C-039 | d | n/a | none: a scalar 2×2 `MatrixExp` |
| σ algebra and [H_JC, n + |e⟩⟨e|] = 0 | G-037 | a | | today: write σ+ = c_e† c_g, σ- = c_g† c_e, σ3 = N_e - N_g; all three commutators reduce to 0 with `BosonicNormalOrder` [probe2 "61"] |
| Detuned Jaynes-Cummings doublet propagator | G-039 | d | n/a | none: a scalar 2×2 matrix per sector; the operator form is the G-038 row |
| Schwinger J's, Casimir, J3 phase action | G-046 | a | | today: `BosonicNormalOrder` gives 0 for both relations; `BosonicExpSimplify` gives e^{iφ/2} a1 (as cosh of √(-φ²)/2) [probe1 "66"] |
| M^T J M = J and its inverse | B-044 | d | n/a | none: matrix identity |
| General transformation kernel for normal symbols | B-054 | d | n/a | none: the kernel is not transcribed in the atlas |
| sp(2N,R) block form of quadratic generators | P-035 | d | n/a | none: matrix condition |
| Fractional-linear action on Z | P-038 | d | n/a | none: matrix map |
| Siegel half-space action | P2-002 | d | n/a | none: matrix map |
| Particle-antiparticle blocks | P2-020 | d | n/a | none: matrix identity |
| Symplectic covariance of Weyl symbols | C-015 | b | no Weyl symbol of an operator anywhere in the layer | Idea 3 (with Idea 7 for sums) |
| S_t = -J Ḟ F⁻¹, metaplectic sign | C-023 | d | n/a | none: matrix ODE and a covering-group sign |
| u, v from the quadrature blocks A, B, C, D | C-024 | a | | today: expand ((AQ + BP) + i(CQ + DP))/√2 with Q, P written in a, a†; the a coefficient is (A + D + i(C - B))/2 (the BosonicExpSimplify conjugation of a in probe2 "132" gives the same row for a one-mode quadratic) |
| Multimode Gaussian vacuum, normalization det(I - ZZ†)^{1/4} | X-006, P-036, B-048 | b | BosonicExp.m:85-89 is one mode only | Idea 1 (multimode form) [probe3 "109"] |
| Coordinate form of the squeezed vacuum | C-028, C-049 | d | n/a | none for the layer: a position-space wavefunction; ⟨α|S|0⟩ would come from Idea 2, the rest is an `Integrate` against the P-013 kernel |
| Reduced Gaussian channel from joint evolution | C-056 | d | n/a | none: covariance-matrix identity |
| d/dt log det M = Tr(M⁻¹ Ṁ) | B-073 | d | n/a | none: matrix identity |
| Riccati fixed point Z for the ground state | B-085 | d | n/a | none: matrix square-root identity |
| SU(1,1) Riccati equation and phase | P2-013 | c | no abstract su(1,1) | Idea 6: the holomorphic carrier of the NonCommutativeAlgebra encoding |
| Floquet classification | P2-014 | d | n/a | none: spectral classification |
| Multimode pair Riccati and vacuum phase | P2-024 | b | BosonicExp.m:123-128 is the two-mode one-pair case | Idea 1 (multimode form) |
| Mehler kernel | C-012 | d | n/a | none for the layer: a position-space kernel. QF already gives its coherent kernel, `BosonicMatrixElement[{al, be}, Exp[Log[w] ad ** a], "Basis" -> "Coherent"]` = e^{wᾱβ}⟨α|β⟩ [probe3 "131"]; the Hermite-sum route to the position kernel does not close in `Sum` |
| Linear flow and Heisenberg action F_t = e^{tJS} | C-018 | a | multimode case not probed | today, one mode: `BosonicExpSimplify[Exp[G] ** a ** Exp[-G]]` for G = x a² + y a†² + w a†a gives a cosh Ω - (w a + 2y a†) sinh Ω/Ω, Ω² = w² - 4xy [probe2 "132"] |
| Gaussian propagation by a linear-fractional map | C-019 | b | BosonicExpOrder has no three-term quadratic rule (BosonicExp.m:85-97) | Idea 1: the map Γ → Γ_t is the Gauss factor of the 2×2 matrix product |
| Siegel invariance and determinant cocycle | C-020 | d | n/a | none: matrix identities |
| Matrix Riccati and amplitude equations | C-021 | d | n/a | none: matrix ODE |
| Coherent elements of a metaplectic operator | C-022 | b | BosonicAlgebra.m:406 and BosonicExp.m:85-97 | Ideas 1 and 2 |
| Squeeze-rotation factorization and vacuum phase | C-029 | b | BosonicExp.m:85-97 | Idea 1; the phase is the {0, 0} Fock element of the ordered form |
| Herman-Kluk, exact for quadratics | C-054 | b | BosonicAlgebra.m:406 | Ideas 1 and 2, then a Gaussian `Integrate` over the label |
| Regularized forced Hamiltonian | B-083 | d | n/a | none: a continuum of modes; the one-mode algebra is one `BosonicNormalOrder` |
| Forced oscillator: displacement plus accumulated phase | P2-010, G-022 | a | | today: `BosonicExpSimplify[Exp[b dt ad - bc dt a] ** Exp[al ad - alc a]]` carries the phase (alc b - al bc) dt/2 = i Im(ᾱβ) dt [probe3 "141"] |
| Parametric oscillator K0, K2 in quadratures | P2-012 | a | | today: `BosonicNormalOrder` of K0 - (N + 1/2)/2 and of K2 - i(a†² - a²)/4, with q, p written in a, a†, gives 0 and 0 [probe3 "142"] |
| Singular oscillator ground state, label k | P2-015 | d | n/a | none: a differential operator with 1/x² on the half line |
| Forced parametric generating kernel | P2-018 | d | n/a | none: the parameter map is a source lookup (Unchecked in the audit) |
| Drive equation and displacement channel | P2-042 | a | | today: `BosonicNormalOrder[a ** Hd - Hd ** a]` is ε; the rest is a scalar ODE |
| Propagator of a driven quadratic Hamiltonian | C-032 | b | BosonicExp.m:72-133, no rule for linear plus quadratic | Idea 1 with the linear terms (the 3×3 affine form of the same matrix) |
| Periodic quadratic evolution, quasi-energies | C-035 | d | n/a | none: needs Mathieu stability data |
| Time-dependent SU(1,1) Hamiltonian | C-050 | b | BosonicExp.m:85-97 | Idea 1: products of one-mode quadratic exponentials are 2×2 matrix products |
| Nondegenerate amplifier gives S2(ηt) | G-075 | a | | today: -iH_I t/ħ is the S2 exponent term by term; the action is the G-068 row |
| Rank-one projector symbols | B-020 | b | BosonicAlgebra.m:406 | Idea 2, with the bare ladder factors of Idea 5 for |1⟩⟨1| = a† :e^{-a†a}: a |
| Kernel acting on a vector | B-026 | b | BosonicAlgebra.m:406 | Idea 2, then a Gaussian `Integrate` with weight e^{-ᾱα} |
| Product of matrix kernels | B-027 | b | BosonicAlgebra.m:406 | Idea 2, then the same `Integrate` |
| Resolvent kernel of a number-type H | B-075 | a | | today: `BosonicMatrixElement[{al, be}, Exp[c Log[xi] ad ** a], "Basis" -> "Coherent"]` gives the kernel exp(ᾱβ ξ^c) (times the overlap normalization), and `BosonicMatrixElement[{n, n}, (c ad ** a - z)^-1]` gives 1/(cn - z) [probe1 "161"] |
| Disk star product to first order in h | P2-006 | d | n/a | none: asymptotic, Unchecked |
| Tr D(γ) = πδ²(γ), operator reconstruction | P-017, P-020 | d | n/a | none: distributional |
| Normally ordered moments from χ_N | P-020, G-029 | b | BosonicAlgebra.m:406 | Idea 2: χ_N of a given state is a coherent-basis trace; the derivative rule itself is a Taylor statement |
| Exact attenuation solution for χ_k | P2-036, P2-039 | d | n/a | none: a scalar PDE in phase space |
| Gaussian convolution kernel for attenuation | P2-038 | d | n/a | none: scalar Gaussian integral |
| Weyl quantization kernel convention | C-013 | d | n/a | none: a definition in position space |
| Symbol traces and Wigner overlaps | C-014 | d | n/a | none: phase-space integrals of Hilbert-Schmidt operators |
| Moyal product, exact commutator for deg H ≤ 2 | C-016 | b | no Weyl symbol | Idea 3 [protoB "moyal183"] |
| Weyl symbol of a metaplectic operator (Cayley) | C-030, C-053 | b | no Weyl symbol; BosonicExp.m:85-97 | Ideas 1 and 3 |
| Translation coefficient of a metaplectic operator | C-031 | b | `Sum` over QF's diagonal Laguerre elements of D(-α) stays unevaluated [probe3 "185"] | Idea 2: Tr = ∫ d²β/π ⟨β| ... |β⟩ after `BosonicExpBraid[..., "Direction" -> "Left"]` moves e^{ᾱa} past e^{-itN} (not probed) |
| Wigner function of a one-mode squeezed Gaussian | C-049 | b | no Weyl symbol | Ideas 2 and 3: the normal symbol of |ψ⟩⟨ψ| is ⟨α|ψ⟩⟨ψ|α⟩ with α, ᾱ independent, then the Gaussian heat formula |
| Weyl symbol of the oscillator exponential | C-053 | b | no Weyl symbol | Idea 3 [protoC "weyl187"] |
| Coherent-projector differential correspondences | G-025 | a | | today: `BosonicMatrixElement[{be, al}, ad, "Basis" -> "Coherent"]` gives β̄⟨β|α⟩ [probe3 "195"]; the identity is the Wirtinger α-derivative of the overlap |
| Glauber P of a number state | G-027 | d | n/a | none: a distribution (QF's `SOrderedFunction` correctly refuses s = 1) |
| Fourier pair between P and χ_N | G-031 | d | n/a | none: scalar Fourier pair |
| Real Gaussian integral | B-022 | d | n/a | none: scalar |
| Doubled complex Gaussian with sources | B-023 | d | n/a | none: scalar |
| det = e^{Tr log}, log-det series | P-048 | d | n/a | none: matrix |
| Commuting Gaussian integral, Fresnel factor | C-091 | d | n/a | none: scalar |
| Tr_B Γ(M) = det(I - M)⁻¹ | X-007, B-107 | b | BosonicAlgebra.m:404 (multimode); BosonicExp.m:79-80 is one mode | Idea 4 for diagonal M; Idea 1 (multimode e^{a†ha} = :e^{a†(e^h - 1)a}:) for general M, then a Gaussian integral |
| Grand partition functions | B-106 | b | BosonicAlgebra.m:404 | Idea 4 for the free limit; the interacting case is d |
| Berezin-Lieb and Golden-Thompson bounds | P-016 | d | n/a | none: inequalities over all H |
| Quasiclassical partition bounds | P2-032 | d | n/a | none: cited theorem |
| Thermal oscillator master equation | P2-034 | a | | today: `BosonicNormalOrder` of the adjoint Lindblad generator applied to a and to a†a gives -(γ/2) a and γν - γ a†a [probe3 "236"] |
| Tr e^{-s(2N+1)} = 1/(2 sinh s) | P2-048 | a | | today: `Sum[BosonicMatrixElement[{n, n}, Exp[-s (2 ad ** a + 1)]], {n, 0, Infinity}]` gives e^s/(e^{2s} - 1) [probe1 "240"] |
| Berezin-Lieb bounds for spin | C-043 | d | n/a | none: inequality |
| ‖e^{itA}‖ ≤ e^{|t|‖T‖} | B-068 | d | n/a | none: norm bound |
| Pair potential on the n-particle sector | B-093 | d | n/a | none: field theory in first quantization |
| Dyson series, e^{+itH} convention | B-094 | d | n/a | none: formal time-ordered series |
| Retarded contraction, T[a*(1) a(2)] | B-097 | a | many-mode labels are d | today, one mode: `BosonicNormalOrder[(e^{-iωt2} a) ** (e^{iωt1} a†)]` gives e^{-iω(t2 - t1)}(a†a + 1), the contraction D(2|1) beside the normal product [probe3 "256"] |
| Wick's theorem as e^{D} | B-099 | d | n/a | none: functional derivatives in field labels |
| Normal symbol of interaction-picture evolution | B-100 | d | n/a | none: formal |
| Contractions on an exponential source | B-101 | d | n/a | none: scalar calculus on sources |
| Determinant identity for a quadratic symbol | B-102 | d | n/a | none: checked in the fermionic reading |
| Auxiliary-field decoupling | B-103 | d | n/a | none: coefficients not transcribed |
| Retarded loops vanish | B-104 | d | n/a | none: formal |
| Second-order Dyson elimination | G-102 | d | n/a | none: an approximation |
| No Hermitian phase operator | G-011 | d | n/a | none: impossibility |
| Susskind-Glogower E E† = 1, E† E = 1 - |0⟩⟨0| | G-012 | b | BosonicAlgebra.m:357-369 sends a bare a to `numberDiagonal`, which fails at line 342 | Idea 5 [protoC "sgEEd", "sgEdE"] |
| Exact resonant Jaynes-Cummings propagator | G-038 | b | BosonicExp.m:418-443 has no rule moving a past f(N) | Idea 5 [protoD "jc"] |
| Dispersive U as conditional number rotation | G-040 | b | BosonicAlgebra.m:404 | Idea 4: per-mode substitution of each number operator [protoD "270"] |
| Raman number-controlled rotation | G-041 | b | BosonicExp.m:418-443 | Idea 5 [protoD "raman"] |
| Coordinate to Bargmann kernel | P-013 | d | n/a | none: Hermite generating function, scalar |
| Laguerre generating function, lower index | P2-053 | d | n/a | none: scalar `Sum` |
| Laguerre generating function along m - n | P2-053 | d | n/a | none: scalar `Sum` |
| Fock elements from the coherent diagonal symbol | G-024 | b | BosonicAlgebra.m:406 | Idea 2, then Wirtinger derivatives of the symbol |
| Weighted conditions for a quadratic normal form | B-078 | d | n/a | none: infinite-mode theorem |
| Square-root perturbation stays in the ideal | B-087 | d | n/a | none: operator square root and an ideal theorem |
| Calogero degeneracy generating function | P2-016 | d | n/a | none: counting |
| [P_j, Q_k] = -iħ δ_jk | C-001 | a | | today: `BosonicNormalOrder[p ** q - q ** p]` with q, p written in a, a† is -i, exact (the truncated `QuadratureOperators` fail in the last level) [probe3 "303"] |
| Periodic-orbit determinant | C-036 | d | n/a | none: semiclassical asymptotics |
| Segal fields [a(f), a†(g)] = ⟨f, g⟩ | C-057 | a | infinite-mode form is d | today for N modes: `BosonicNormalOrder` gives Σ f̄_k g_k (three modes) [probe3 "305"] |
| Collective broadband commutator | G-094 | a | infinite band is d | today for a finite band: the C-057 call with f = g = l gives C = Σ|l_k|² [probe3 "305"] |
| Continuum CCR, pair-creation kernel | G-105 | d | n/a | none: continuum modes |
| SU(2) disentangling | X-004, P-025, C-041 | a | | today: `BosonicExpOrder[Exp[x ad1 ** a2 - y ad2 ** a1]]` in the Schwinger realization gives e^{ζJ+} (cos r)^{-2J0} e^{-ζ̄J-} with ζ = e^{iφ} tan r (written with Tanh[√(-xy)]/√(-xy)) [probe1 "309"] |
| su(2) commutators, ladder action, holomorphic form | P-024, C-038 | a | holomorphic form is Idea 6's carrier | today: commutators by `BosonicNormalOrder` (as in G-046); `BosonicMatrixElement[{{j+m+1, j-m-1}, {j+m, j-m}}, ad1 ** a2]` gives √(j-m) √(j+m+1) with symbolic j, m [probe3 "311"] |
| SU(2) coherent state expansion | P-026 | a | symbolic j needs Idea 4 | today for fixed j: `BosonicExpSimplify[Exp[z ad1 ** a2] ** ad2 ** Exp[-z ad1 ** a2]]` is ad2 + z ad1, so e^{ζJ+}|0, 2j⟩ is (a2† + ζ a1†)^{2j}|0⟩/√(2j)!, whose elements the multimode polynomial branch gives |
| SU(2) overlap, annihilator, resolution | P-026, C-042 | a | overlap and resolution are scalar sums and integrals | today: `BosonicExpSimplify[Exp[z ad1 ** a2] ** (ad2 ** a1) ** Exp[-z ad1 ** a2]]` gives (z a1† + a2†) ** (a1 - z a2) = J- + 2zJ0 - z²J+ [probe1 "313"] |
| su(1,1) algebra, Casimir, ladder, holomorphic form | P-028, C-045 | c | QF realizes only k = 1/4, 3/4 (one mode) and k = (1 + |n1 - n2|)/2 (two modes) | Idea 6 [protoE] |
| SU(1,1) coherent expansion and overlap | P-030, C-047 | c | as above | Idea 6: the norms n(2k + n - 1) give the coefficients √(Γ(n+2k)/(n! Γ(2k))); the oscillator case k = 1/4 of the annihilator is today, `BosonicExpSimplify[Exp[z ad ** ad/2] ** a ** a ** Exp[-z ad ** ad/2]]` = (a - z a†) ** (a - z a†) [probe3 "315"] |
| SU(1,1) resolution of identity | P-030 | d | n/a | none: scalar radial integral once the coefficients are known |
| Root-space coherent states | P-047 | d | n/a | none: general semisimple group, Unchecked |
| Compact Kähler kernel | P2-007 | d | n/a | none: scalar |
| Jacobians of classical domains | P2-056 | d | n/a | none: Unchecked convention |
| Completeness for discrete SU(1,1) subsystems | P2-004 | d | n/a | none: cited theorem |
| Clock and shift matrices | C-037 | d | n/a | none: finite N × N matrices |
| Gauss decomposition of SL(2,C) | C-040 | d | n/a | none: scalar 2×2 identity (Idea 1 uses it) |
| Spin-to-oscillator contraction | C-044 | d | n/a | none: scalar limit |
| Singular oscillator lowest weight | C-052 | d | n/a | none: half-line differential operator |
| Berezin symbol product on the disk | C-097 | d | n/a | none: integral over the disk |

## Ranked design ideas

Counts are table rows the idea is needed for. "Alone" means the idea by itself makes the row computable; "with" means the row needs a second idea too. Ideas 1 to 5 and 7 each extend one existing function; none adds a symbol. Idea 6 needs no QF change.

### Idea 1. One rule for the exponential of any quadratic (10 rows: 6 alone, 4 with Ideas 2, 3 or 4)

What it does: `BosonicExpOrder` normal-orders e^{G} for every quadratic G, including the three-term one-mode case x a² + y a†² + w a†a, which it now returns untouched [probe2 "quad3"], and products of such exponentials.

Where: `BosonicExpOrder`, rulesNO at BosonicExp.m:85-128. The one-mode squeeze, two-mode squeeze and beam-splitter rules there (lines 85-97, 115-128) are special cases and would become instances of this rule.

Rule, one mode: with K+ = a†²/2, K- = a²/2, K0 = (a†a + 1/2)/2, send G to the 2×2 matrix 2y K+ + 2x K- + 2w K0 in the representation K+ = ((0,1),(0,0)), K- = ((0,0),(-1,0)), K0 = diag(1/2, -1/2), take g = `MatrixExp` of it, and read the Gauss factors of C-040: p = g12/g22, q = -2 log g22, r = -g21/g22. Then

e^{G} = e^{-w/2} e^{p a†²/2} e^{q (a†a + 1/2)/2} e^{r a²/2}.

A product of quadratic exponentials maps to the product of their matrices, which is how C-019, C-029 and C-050 come out. For N modes the same reading of the 2N × 2N matrix of the adjoint action gives e^{½ a†Z a†} :e^{a†(Y - 1)a}: e^{½ aW a} with the scalar det(Y)^{-1/2}, and with linear terms added it is the 3×3 (affine) version of the same matrix. The branch of the scalar is fixed by continuity from G = 0.

Unlocks: C-019, C-029, C-050, C-032, X-006/P-036/B-048, P2-024 alone; C-022 and C-054 with Idea 2; C-030 with Idea 3; X-007/B-107 for a general M with Idea 4.

Evidence: (protoB) the Fock elements of the ordered form, `BosonicMatrixElement[{m, n}, ordered[t x, t y, t w]]`, agree with QF's exact elements of the polynomial series of e^{t(x a² + y a†² + w a†a)} for every m, n ≤ 3 through order t⁴ (16 differences, all 0). (probe3 "109") for two modes with symbolic symmetric A and B, the vacuum element ⟨0| e^{½ aAa} e^{½ a†Ba†} |0⟩ computed from QF's multimode polynomial branch agrees with det(I - tAB)^{-1/2} through order t³ (difference 0); that determinant is the scalar the multimode rule must produce, and with A = -Z̄, B = -Z it is the X-006 normalization.

### Idea 2. Coherent-basis elements of ordered exponential products (9 rows: 5 alone, 4 with Ideas 1, 3 or 5)

What it does: lets `BosonicMatrixElement[{α, β}, expr, "Basis" -> "Coherent"]` evaluate any product that `BosonicExpOrder` brings to the form e^{P(a†)} f(N) e^{Q(a)}, not only functions of N.

Where: `BosonicMatrixElement`, coherent branch at BosonicAlgebra.m:406-409. After `numberDiagonal` fails, retry with `BosonicExpOrder` as the Fock branch already does (lines 412-414), and tag the factors with the existing `normalFactor` (line 357).

Rule, with A = ᾱ and B = β:

⟨α| e^{P(a†)} f(N) e^{Q(a)} |β⟩ = e^{P(A) + Q(B)} ⟨α| f(N) |β⟩,

because a|β⟩ = β|β⟩, ⟨α|a† = ᾱ⟨α|, and both exponential series are entire. The last factor is what line 409 already returns. No Taylor coefficient is taken, so the integer-index restriction of the Fock branch (line 394) does not arise; the Up, Diagonal, Down order check of line 384 must stay. Bare ladder factors beside f(N) (B-020) multiply by A^p B^q.

Unlocks: B-026, B-027, P-020/G-029, C-031, G-024 alone; C-022 and C-054 with Idea 1; C-049 with Idea 3; B-020 with Idea 5.

Evidence (protoA, a scratch function outside QF that calls QF's own `BosonicExpOrder` and `numberDiagonal`): for the displacement it returns exp(-g gc/2 + gA - gc B), exactly the closed form in the earlier coherent-basis finding; for the squeeze e^{(x a² - y a†²)/2} it returns √(sech √(xy)) exp(AB(sech √(xy) - 1) + (xB² - yA²) tanh √(xy)/(2√(xy))); it also handles e^{g a†} e^{cN} e^{-gc a}. In all three cases e^{AB} times the result agrees, coefficient by coefficient through total order 4 or 5 in A and B, with the generating function Σ A^m B^n ⟨m|expr|n⟩/√(m! n!) built from QF's exact Fock elements (differences simplify to 0).

### Idea 3. The s-ordered symbol of a normal-ordered expression (5 rows: 2 alone, 3 with Ideas 1, 2 or 7)

What it does: returns the Weyl symbol, or any s-ordered symbol, of an operator. Nothing in the layer produces one today: the functions in PhaseSpaceRepresentations.m act on truncated states.

Where: a property beside `NormalOrdered[expr]["Series"]` (BosonicExp.m:242-250), for example `NormalOrdered[expr]["Symbol", s]`, fed by `BosonicNormalOrder` for polynomials and by `BosonicExpOrder` for exponentials.

Rule: with A = ᾱ and B = α, the s-ordered symbol is e^{-(1-s)/2 ∂_A ∂_B} applied to the normal symbol. On a polynomial the series stops. On a Gaussian it has a closed form, e^{λ∂_A∂_B} e^{cAB} = (1 - λc)^{-1} e^{cAB/(1 - λc)}, and on a general Gaussian it is a Gaussian integral.

Unlocks: C-016 and C-053 alone; C-015 with Idea 7; C-030 with Idea 1; C-049 with Idea 2.

Evidence (protoB, protoC): for C-053, `BosonicExpOrder[Exp[-tau ad ** a]]` gives `NormalOrdered[Exp[(E^-tau - 1) ad ** a]]`, and the Gaussian formula at s = 0 turns e^{-τ/2} times it into sech(τ/2) e^{-2 tanh(τ/2)|α|²}, the atlas symbol with |X|² = 2|α|² (difference 0). For C-016, with H = a†a + κ(a² + a†²) + μ a†² and B = a†²a + ν a³, the Weyl symbol of `BosonicNormalOrder[H ** B - B ** H]` equals the Poisson bracket of the Weyl symbols exactly; with the quartic H = a†²a² the same comparison fails, as it must.

### Idea 4. Several modes in the exponential branch of BosonicMatrixElement (4 rows: 3 alone, 1 with Idea 1; also symbolic j in P-026)

What it does: lets the Fock branch for exponentials take index vectors, as the polynomial branch already does.

Where: `BosonicMatrixElement`, BosonicAlgebra.m:401-404, which returns `$Failed` for list indices and for more than one field variable [probe1 "47fock", "228"].

Rule: a function of the number operators is diagonal, so substitute each mode's number (`numberDiagonal` once per mode). An exponential of creation operators alone, e^{P(a1†, ..., ak†)}, has the element Π_i √(m_i!/k_i!) times the coefficient of x^{m-k} in e^{P(x)} (the several-variable form of the helper `wing`, line 349), and the element of a product sums over the intermediate occupation vector k ≤ min(m, n).

Unlocks: X-003/G-072, B-106 (free limit), G-040 alone; X-007/B-107 with Idea 1 (diagonal M needs only this idea); the symbolic-j case of P-026.

Evidence (protoB, protoC, protoD): `BosonicExpOrder` already writes the two-mode squeezer as sec Ω e^{c a1†a2†} ** (sec Ω)^{N1+N2} ** e^{c' a1 a2}; reading its vacuum column with the coefficient of (x1 x2)^n gives sech r (-e^{iθ} tanh r)^n for symbolic n (difference 0) and for n = 0 to 4. For G-040, substituting each mode's number in e^{-iχt(N_e + N_a(N_e - N_g))} gives e^{-iχt(n+1)} on the excited branch and e^{iχtn} on the ground branch, and QF's `numberDiagonal` already reads each one-mode factor.

### Idea 5. Functions of N beside bare ladder operators (4 rows: 3 alone, 1 with Idea 2)

What it does: handles a f(N), f(N) a† and a† f(N) a, the products in the Jaynes-Cummings and phase-operator rows.

Where, two places:

- Fock elements: `normalFactor` (BosonicAlgebra.m:357-369) passes a bare a to `numberDiagonal`, which fails at line 342 [probe1 "268"]. Tag a^q as a Down monomial and a†^p as an Up monomial, with ⟨m| a†^p f(N) a^q |n⟩ = δ(m - p, n - q) √(m!/(m-p)!) √(n!/(n-q)!) f(n - q) for n ≥ q and 0 for n < q. The same branch can drop the integer requirement of line 394 when one index is 0, since an exponential of annihilators acting on the vacuum (or of creators met by the vacuum bra) contributes only its constant term: today `{4, 0}` evaluates and `{m, 0}` does not [probe1 "121edge"].
- Operator identities: a rule in `BosonicExpSimplify` (expReduceRules, BosonicExp.m:418-443) that moves a past f(N): a f(N) = f(N + 1) a and f(N) a† = a† f(N + 1). Today `BosonicExpSimplify[a ** (ad ** a + 1)^(-1/2)]` comes back unchanged [probe1 "shift"].

A trap the design has to respect: a† f(N) a = N f(N - 1), and Wolfram Language cancels N · N⁻¹ to 1 on sight, so the operator rule alone turns E†E into 1 and loses the vacuum. That product must go through the Fock rule, whose n ≥ q condition keeps it.

Unlocks: G-012, G-038, G-041 alone; B-020 with Idea 2.

Evidence (protoC, protoD): with the shift rule, E E† = (N+1)^{-1/2} a a† (N+1)^{-1/2} reduces to 1; with the Fock rule, ⟨n|E†E|n⟩ is 1 for n ≥ 1 and 0 at n = 0, so E†E = 1 - |0⟩⟨0|. For the resonant Jaynes-Cummings propagator in the atlas's 2×2 block form (entries cos(λt√(N+1)), -i a sin(λt√N)/√N, -i a† sin(λt√(N+1))/√(N+1), cos(λt√N)), H U - i dU/dt reduces to the zero block and U(0) is the identity; the Raman propagator cos(λtN) - i sin(λtN) σ1 also gives the zero block.

### Idea 6 (class c). su(1,1) and su(2) with an arbitrary label, by NonCommutativeAlgebra (3 rows; also the holomorphic part of P-024/C-038)

What it does: covers the rows that need a representation label k (or j) that no bosonic realization reaches. QF's realizations give k = 1/4, 3/4 (one mode) and k = (1 + |n1 - n2|)/2 (two modes) only, which is why the audit's P-028 check used k = 0.7.

Where: no change to a QF function. A short recipe beside `BosonicRelations` (BosonicAlgebra.m:25), building the algebra the way the layer already calls `NonCommutativePolynomialReduce` (line 47).

Encoding, by the `nca-route` procedure: generators {K+, K0, K-} with [K-, K+] = 2K0 and [K0, K±] = ±K± typed as `"CommutationRelations"` (finite Gröbner bases, a canonical normal form, left ideals). The lowest-weight module is the quotient by the left ideal ⟨K-, K0 - k⟩, with k a scalar. Rejected: the free algebra knows no relations; the same relations typed as `"Relations"` can return an infinite basis truncated without warning; no named structure fits su(1,1); a matrix carrier needs a finite dimension, which the discrete series does not have.

Rule the rows use: K- K+^n ≡ n(2k + n - 1) K+^{n-1} modulo the left ideal gives the norms, hence the coherent-state coefficients √(Γ(n + 2k)/(n! Γ(2k))); ad_{K+}³(K-) = 0 makes the annihilator series e^{ζK+} K- e^{-ζK+} finite.

Unlocks: P-028/C-045, P-030/C-047, P2-013; the holomorphic form in P-024/C-038 (the same recipe with su(2) relations).

Evidence (protoC, protoE): the commutation set comes back as its own Gröbner basis, so it is admissible; the Casimir K0² - (K+K- + K-K+)/2 commutes with all three generators; ad_{K+} applied to K- gives -2K0, then 2K+, then 0; the left Gröbner basis of ⟨K-, K0 - k⟩ has length 2 at caps 10, 20 and 40, so it is finite and does not track the cap; K- K+^n - n(2k + n - 1) K+^{n-1} reduces to 0 for n = 1 to 6; the Casimir reduces to k(k - 1), and the wrong value k(k + 1) does not reduce to 0. Independent carrier: the holomorphic operators K- = d/dz, K+ = z² d/dz + 2kz, K0 = z d/dz + k satisfy both relations, the same ladder identity on the constant function for n = 1 to 4, and the Casimir value k(k - 1), all exactly.

### Idea 7. Conjugation that is linear over sums and scalars (2 rows: 1 alone, 1 with Idea 3)

What it does: makes `BosonicExpSimplify[Exp[B] ** X ** Exp[-B]]` work when X is a sum of monomials or a scalar times a monomial.

Where: `expConjugate`, BosonicExp.m:410-416. A sum goes whole into `adAction` (line 396), whose test (the second commutator a scalar multiple of X) fails for most sums, so the expression comes back unchanged.

Rule: map over the terms of a `Plus` and pull scalars out first; each monomial is already conjugated factor by factor by the rule at lines 426-430.

Unlocks: B-053 alone; C-015 with Idea 3. It also lets the k = 1/4 annihilator of P-030 be written as (a ** a)/2.

Evidence: conjugating `ad1 ** a1 ** a2` by a beam splitter returns the product of the three rotated factors, while `ad1 ** a1 + a2` returns unchanged [probe1 "58", "58sum"]; mapping the same call over the two terms returns both rotated terms [probe3 "58map"]. Likewise `a ** a` conjugated by e^{z a†²/2} gives (a - z a†) ** (a - z a†) [probe3 "315"], while (a ** a)/2 comes back unchanged [probe1 "315"].

## Entries out of reach for a finite exact layer (class d)

Grouped by reason; each reason applies to every entry listed after it.

- Scalar or matrix identities with no operator in them; plain Wolfram Language (`MatrixExp`, `Sum`, `Integrate`, `Det`) already does them: C-039, G-039, B-044, P-035, P-038, P2-002, P2-020, C-023, C-056, B-073, B-085, C-020, C-021, P2-038, B-022, B-023, P-048, C-091, P-013, P2-053 (both Laguerre sums), P2-016, P2-007, C-037, C-040, C-044, G-031, P2-036/P2-039.
- Infinite-mode, continuum or field-theory statements: B-083, B-078, B-087, G-105, B-093, B-099, B-100, B-101, B-103, B-104, B-094.
- Distributions: P-017/P-020 (Tr D(γ) = πδ²), G-027 (P function of a number state).
- Inequalities and norm bounds: P-016, P2-032, C-043, B-068.
- Existence or impossibility statements: G-014, G-011.
- Spectral or stability classifications and cited theorems: P2-014, C-035, P2-004, P2-056, P-047, C-036, P2-006.
- Position-space differential operators on the half line: P2-015, C-052.
- Position-space wavefunctions and kernels, which the Fock-and-ladder layer does not represent: C-028/C-049, C-012.
- Definitions or kernels not transcribed in the atlas, or source lookups: B-054, C-013, P2-018, B-102.
- Phase-space integrals of Hilbert-Schmidt operators and integral resolutions: C-014, P-030 (resolution of identity), C-097.
- Approximations: G-102.
- Fock-series states with no exponential to evaluate: G-090.
- Duhamel's formula for an arbitrary L(t): C-051.
