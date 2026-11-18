---
Template: TechNote
Name: BosonicAlgebra
Title: Bosonic Operator Algebra
Context: Wolfram`QuantumFramework`
CellContext: Global`
ContextPath: [Wolfram`QuantumFramework`SecondQuantization`]
Paclet: Wolfram/QuantumFramework
URI: Wolfram/QuantumFramework/tutorial/BosonicAlgebra
Keywords: [bosonic, ladder operators, normal order, anti-normal order, disentangling, Gauss decomposition, Baker-Campbell-Hausdorff, Zassenhaus, displacement operator, squeeze operator, beam splitter, matrix element, coherent state, normal symbol]
RelatedGuides: [WolframQuantumComputationFramework]
RelatedTutorials: [SecondQuantization, TimeEvolution]
Typeset: _SuperDagger -> StandardForm
---

This tech note introduces the symbolic bosonic algebra of the QuantumFramework. The ladder operators $a$ and $a^\dagger$ are treated as symbols obeying $[a, a^\dagger] = 1$, so ordering identities, exponentials of quadratic generators and their matrix elements come out exact, with no Fock-space truncation. A truncation enters only at the end, when a symbolic result is realized as a matrix.

Install and load the QuantumFramework paclet and the second-quantization subpackage (skip the install step if it is already available):

```wl
#| eval: false
PacletInstall["https://www.wolfr.am/DevWQCF", ForceVersionInstall -> True]
```

```wl
<< Wolfram`QuantumFramework`
<< Wolfram`QuantumFramework`SecondQuantization`
```

<!-- #| tags: core -->
## Bosonic algebra functions

### Ladder variables and polynomials

<!-- #| style: DefinitionBox -->
|   |   |
|---|---|
| [FieldVariables]()[] | Ladder variables $\{\{a^\dagger\}, \{a\}\}$ of one mode; [FieldVariables]()[{*label*₁, ...}] gives one mode per label |
| [BosonicNormalOrder]()[*expr*, *vars*] | *expr* rewritten with every creator to the left of every annihilator |
| [BosonicAntinormalOrder]()[*expr*, *vars*] | *expr* rewritten with every annihilator to the left of every creator |
| [BosonicVEV]()[*expr*, *vars*] | Vacuum expectation value ⟨0\|*expr*\|0⟩ |
| [BosonicBCHTerms]()[{*A*, *B*}, *n*, *vars*] | *n*-th order term of $\log(e^{A} e^{B})$ |
| [BosonicZassenhausTerms]()[{*A*, *B*}, *n*, *vars*] | *n*-th term of the Zassenhaus factorization of $e^{A+B}$ |

### Exponentials

<!-- #| style: DefinitionBox -->
|   |   |
|---|---|
| [BosonicExpOrder]()[*expr*, *vars*] | Product of exponentials of ladder operators brought into normal order |
| [NormalOrdered]()[*expr*, *vars*] | Normal-ordered product :*expr*: |
| [AntinormalOrdered]()[*expr*, *vars*] | Anti-normal-ordered product of *expr* |
| [BosonicExpSimplify]()[*expr*, *vars*] | Merges products of exponentials and functions of number operators |
| [BosonicExpBraid]()[*expr*, *vars*] | Moves the left exponential of a pair past the right one, $e^{A} e^{B} = e^{e^{A} B e^{-A}} e^{A}$ |

### Matrix elements and matrices

<!-- #| style: DefinitionBox -->
|   |   |
|---|---|
| [BosonicMatrixElement]()[{*m*, *n*}, *expr*, *vars*] | Fock element ⟨*m*\|*expr*\|*n*⟩, several modes with {{*m*₁, ...}, {*n*₁, ...}} |
| [BosonicMatrixElement]()[{α, β}, *expr*, *vars*, "Basis" -> "Coherent"] | Coherent element ⟨α\|*expr*\|β⟩ |
| [ToBosonicOperator]()[*expr*, *vars*, *size*] | *expr* realized as a QuantumOperator with each mode truncated to dimension *size* |

<!-- #| tags: details -->
## Details and options

### Declaring the ladder operators

<!-- #| style: DefinitionBox -->
|   |   |
|---|---|
| {{$c_1$, ...}, {$a_1$, ...}} | Creators first, then annihilators, one entry per mode |
| {$a_1$, ...} | Annihilators alone, the creators are their daggers |
| WeylAlgebra[...] | A Weyl algebra built elsewhere |

### Options

<!-- #| style: DefinitionBox -->
|   |   |
|---|---|
| "Ordering" -> *st* | For [BosonicExpOrder](), *st* can be "Normal" (default), "Antinormal" or "Disentangled" |
| Assumptions -> *assum* | Assumptions on the coefficients, used to pick branches and simplify |
| TimeConstraint -> *t* | Seconds allowed for ordering (default 30); past it the input is returned unordered with a message |
| "Basis" -> "Coherent" | Matrix elements between coherent states instead of Fock states |
| Method -> "Blasiak" | [BosonicNormalOrder]() through the Blasiak closed-form formula instead of the default reduction |

### Ordering identities

<!-- #| style: DefinitionBox -->
|   |   |
|---|---|
| Displacement | $e^{\alpha a^\dagger - \alpha^{*} a} = e^{-\|\alpha\|^2/2} e^{\alpha a^\dagger} e^{-\alpha^{*} a}$ |
| Number exponential | $e^{\lambda a^\dagger a} = \colon e^{(e^{\lambda}-1) a^\dagger a} \colon$ |
| Passive mode mixing | $e^{a^\dagger W a} = \colon e^{a^\dagger (e^{W}-1) a} \colon$ |
| Gauss form ("Disentangled") | $\colon e^{a^\dagger (M-1) a} \colon = e^{a^\dagger (\log U) a} \Lambda^{N} e^{a^\dagger (\log L) a}$ for $M = U \Lambda L$, $U$ upper and $L$ lower unitriangular |

<!-- #| tags: basic1 -->
## Basic usage

### Ladder variables

Ladder variables of one mode:

```wl
FieldVariables[]
```

Two modes labeled 1 and 2:

```wl
FieldVariables[{1, 2}]
```

Name the one-mode variables used below:

```wl
{{ad}, {a}} = v = FieldVariables[];
```

### Normal and anti-normal order

The commutator leaves a c-number when $a a^\dagger$ is brought into normal order:

```wl
BosonicNormalOrder[a ** ad, v]
```

Higher powers collect the combinatorial factors:

```wl
BosonicNormalOrder[GeneralizedPower[NonCommutativeMultiply, a, 2] ** GeneralizedPower[NonCommutativeMultiply, ad, 2], v]
```

Anti-normal order picks up the opposite sign:

```wl
BosonicAntinormalOrder[ad ** a, v]
```

The Blasiak formula and the default reduction agree:

```wl
With[{w = GeneralizedPower[NonCommutativeMultiply, a, 3] ** GeneralizedPower[NonCommutativeMultiply, ad, 3]},
 BosonicNormalOrder[w, v, Method -> "Blasiak"] === BosonicNormalOrder[w, v]]
```

### Vacuum expectation values

Vacuum expectation value of $a a^\dagger$:

```wl
BosonicVEV[a ** ad, v]
```

The moments $\langle 0 | a^m a^{\dagger m} | 0 \rangle = m!$ check the algebra at every degree:

```wl
Table[BosonicVEV[GeneralizedPower[NonCommutativeMultiply, a, m] ** GeneralizedPower[NonCommutativeMultiply, ad, m], v], {m, 5}]
```

### Baker-Campbell-Hausdorff and Zassenhaus

Second-order term of $\log(e^{\alpha a^\dagger} e^{\beta a})$:

```wl
BosonicBCHTerms[{\[Alpha] ad, \[Beta] a}, 2, v]
```

The commutator $[a, a^\dagger]$ is a number, so every higher term vanishes:

```wl
BosonicBCHTerms[{\[Alpha] ad, \[Beta] a}, 3, v]
```

Second Zassenhaus term, with the opposite sign:

```wl
BosonicZassenhausTerms[{\[Alpha] ad, \[Beta] a}, 2, v]
```

## Ordering exponentials

An exponential of a sum of non-commuting operators does not factor for free. For generators of degree at most two, the factors are fixed by a $2N \times 2N$ matrix built from the coefficients (Wilcox 1967), and [BosonicExpOrder]() computes them exactly.

### One mode

Displacement operator in normal order:

```wl
BosonicExpOrder[Exp[\[Alpha] ad - Conjugate[\[Alpha]] a], v]
```

Squeeze operator, with $r > 0$ declared so the hyperbolic functions come out directly:

```wl
BosonicExpOrder[Exp[r/2 (GeneralizedPower[NonCommutativeMultiply, a, 2] - GeneralizedPower[NonCommutativeMultiply, ad, 2])], v, Assumptions -> r > 0]
```

An exponential of the number operator has no finite factorization, and comes back as a normal-ordered product:

```wl
BosonicExpOrder[Exp[\[Lambda] ad ** a], v]
```

A driven and squeezed oscillator, $H = a^\dagger a + \tfrac{1}{3}(a^2 + a^{\dagger 2}) + \tfrac{1}{5}(a + a^\dagger)$, evolved for a symbolic time:

```wl
h1 = ad ** a + 1/3 (GeneralizedPower[NonCommutativeMultiply, a, 2] + GeneralizedPower[NonCommutativeMultiply, ad, 2]) + 1/5 (ad + a);
u1 = BosonicExpOrder[Exp[-I t h1], v];
```

The ordered form carries a square root of a time-dependent quantity, and the principal branch of that root flips sign every half period. The engine follows the branch continuously in $t$, so the vacuum amplitude stays correct at long times. Against a truncated evolution at $t = 7.3$ with 40, 60 and 80 levels, the difference is the truncation error and goes to zero:

```wl
Table[With[{low = N@DiagonalMatrix[Sqrt[Range[s - 1]], 1]},
  Abs[N[BosonicMatrixElement[{0, 0}, u1, v] /. t -> 73/10] -
    MatrixExp[-7.3 I (h1 /. NonCommutativeMultiply -> Dot /. GeneralizedPower[Dot, m_, k_] :> MatrixPower[m, k] /. {ad -> Transpose[low], a -> low})][[1, 1]]]], {s, {40, 60, 80}}]
```

### Several modes

Two modes for the examples below:

```wl
{{bd1, bd2}, {b1, b2}} = v2 = FieldVariables[{1, 2}];
```

Two-mode squeezer, with wings creating and annihilating pairs and a number middle:

```wl
BosonicExpOrder[Exp[r (bd1 ** bd2 - b1 ** b2)], v2, Assumptions -> r > 0]
```

A beam splitter conserves the photon number, and its normal-ordered form is a single middle that mixes the modes, $:e^{a^\dagger (e^{W}-1) a}:$:

```wl
bs = BosonicExpOrder[Exp[-I t (bd1 ** b2 + bd2 ** b1)], v2, Assumptions -> Element[t, Reals]]
```

Writing out [Sinc]() gives the familiar $\cos t$ and $\sin t$ of the beam splitter:

```wl
FullSimplify[FunctionExpand[bs], Element[t, Reals]]
```

Coupled evolution with unequal frequencies, a beam splitter and two-mode squeezing:

```wl
h2 = bd1 ** b1 + 3/2 bd2 ** b2 + 1/4 (bd1 ** b2 + bd2 ** b1) + 1/5 (bd1 ** bd2 + b1 ** b2);
BosonicExpOrder[Exp[-0.7 I h2], v2]
```

### Anti-normal order

Displacement operator in anti-normal order, the Gaussian factor changes sign:

```wl
BosonicExpOrder[Exp[\[Alpha] ad - Conjugate[\[Alpha]] a], v, "Ordering" -> "Antinormal"]
```

### Disentangled form

The normal-ordered middle of a passive evolution is the second-quantized image $\Gamma(M)$ of a matrix $M$, and $\Gamma(M_1 M_2) = \Gamma(M_1)\Gamma(M_2)$. Factoring $M = U \Lambda L$ therefore factors the operator into raising, diagonal and lowering exponentials. For a beam splitter this is the SU(2) disentangling (Truax 1985):

```wl
BosonicExpOrder[Exp[x bd1 ** b2 - y bd2 ** b1], v2, "Ordering" -> "Disentangled"]
```

The factorization needs a non-vanishing pivot, and the Gauss form itself is singular at the 50:50 point, where $\tan$ diverges. There the middle stays normal ordered:

```wl
BosonicExpOrder[Exp[-I Pi/2 (bd1 ** b2 + bd2 ** b1)], v2, "Ordering" -> "Disentangled"]
```

## Ordered products

On a polynomial, [NormalOrdered]() is the normal-ordered product $:\cdot:$ of field theory, which moves creators left without applying the commutator. It is a different operator from [BosonicNormalOrder](), and the two differ by the vacuum term:

```wl
NormalOrdered[a ** ad, v]
```

```wl
BosonicNormalOrder[a ** ad, v]
```

## Merging and braiding

[BosonicExpSimplify]() goes the other way, merging exponentials whenever the algebra closes. For linear generators Baker-Campbell-Hausdorff terminates and the merge is exact:

```wl
BosonicExpSimplify[Exp[\[Alpha] ad] ** Exp[\[Beta] a], v]
```

A conjugation by a displacement shifts the mode:

```wl
BosonicExpSimplify[Exp[\[Alpha] ad - \[Alpha] a] ** a ** Exp[\[Alpha] a - \[Alpha] ad], v]
```

[BosonicExpBraid]() only moves one exponential past the other, transporting its generator:

```wl
BosonicExpBraid[Exp[\[Alpha] ad] ** Exp[\[Beta] a], v]
```

## Matrix elements

### Fock basis

The displacement elements are Laguerre polynomials, also for a symbolic index:

```wl
BosonicMatrixElement[{n, n}, DisplacementOperator[\[Alpha]]]
```

Squeezed vacuum amplitude $\langle 2 | S | 0 \rangle$ from the ordered exponential:

```wl
BosonicMatrixElement[{2, 0}, Exp[r/2 (GeneralizedPower[NonCommutativeMultiply, a, 2] - GeneralizedPower[NonCommutativeMultiply, ad, 2])], v, Assumptions -> r > 0]
```

The two-mode squeezer creates photons in pairs, $\langle k, k | S_2 | 0, 0 \rangle = \tanh^k r / \cosh r$:

```wl
Table[Simplify[BosonicMatrixElement[{{k, k}, {0, 0}}, Exp[r (bd1 ** bd2 - b1 ** b2)], v2, Assumptions -> r > 0], r > 0], {k, 0, 3}]
```

Photon hopping $\langle 1, 0 | U | 0, 1 \rangle$ for the coupled evolution at $t = 1, 2, 3$:

```wl
uh = BosonicExpOrder[Exp[-I t h2], v2];
hop = Table[BosonicMatrixElement[{{1, 0}, {0, 1}}, uh, v2] /. t -> tt, {tt, 1, 3}] // N
```

Check against a truncated evolution with 22 levels per mode:

```wl
With[{low = N@DiagonalMatrix[Sqrt[Range[21]], 1]},
 With[{m1 = KroneckerProduct[low, IdentityMatrix[22]], m2 = KroneckerProduct[IdentityMatrix[22], low]},
  hop - Table[MatrixExp[-I tt (h2 /. NonCommutativeMultiply -> Dot /. {bd1 -> Transpose[m1], bd2 -> Transpose[m2], b1 -> m1, b2 -> m2})][[23, 2]], {tt, 1., 3.}]]] // Chop
```

A function of the number operator between ladder operators, $a_1^\dagger f(N_1) a_1 | n \rangle = n f(n-1) | n \rangle$, with $f(N) = e^{N} - e$ vanishing at the index but not one below it:

```wl
BosonicMatrixElement[{{1, 0}, {1, 0}}, bd1 ** (Exp[bd1 ** b1] - E) ** b1, v2]
```

### Coherent basis

Coherent amplitude of the two-mode squeezed vacuum:

```wl
amp = Simplify[BosonicMatrixElement[{{\[Alpha]1, \[Alpha]2}, {0, 0}}, Exp[r (bd1 ** bd2 - b1 ** b2)], v2, "Basis" -> "Coherent", Assumptions -> r > 0], r > 0]
```

Same result as the textbook form $e^{(\tanh r) \alpha_1^{*} \alpha_2^{*}} \langle \alpha | 0 \rangle / \cosh r$:

```wl
Simplify[amp - Sech[r] Exp[Tanh[r] Conjugate[\[Alpha]1] Conjugate[\[Alpha]2]] Exp[-(\[Alpha]1 Conjugate[\[Alpha]1] + \[Alpha]2 Conjugate[\[Alpha]2])/2], r > 0]
```

### Normal symbol

The normal symbol $F(\alpha^{*}, \beta)$ of an operator is defined by $\langle \alpha | A | \beta \rangle = \langle \alpha | \beta \rangle F(\alpha^{*}, \beta)$ (Cahill and Glauber 1969), so it is the coherent element divided by the overlap. For a number exponential it is $e^{(e^{\lambda}-1)\alpha^{*} \beta}$:

```wl
Simplify[BosonicMatrixElement[{\[Alpha], \[Beta]}, Exp[\[Lambda] ad ** a], v, "Basis" -> "Coherent"] /
  Exp[Conjugate[\[Alpha]] \[Beta] - (\[Alpha] Conjugate[\[Alpha]] + \[Beta] Conjugate[\[Beta]])/2]]
```

### Traces

The trace of $x_1^{N_1} x_2^{N_2}$ is a product of geometric series, $1/((1 - x_1)(1 - x_2))$:

```wl
Simplify[Sum[BosonicMatrixElement[{{k1, k2}, {k1, k2}}, Exp[Log[x1] bd1 ** b1 + Log[x2] bd2 ** b2], v2], {k1, 0, Infinity}, {k2, 0, Infinity}], 0 < x1 < 1 && 0 < x2 < 1]
```

## From symbols to matrices

[ToBosonicOperator]() realizes a symbolic result as a QuantumOperator at a chosen truncation. The normal-ordered displacement, realized on 10 levels:

```wl
d = ToBosonicOperator[BosonicExpOrder[Exp[0.4 (ad - a)], v], v, 10]
```

On the lowest levels the truncated matrix agrees with $e^{0.4(a^\dagger - a)}$ computed on 40 levels:

```wl
With[{low = N@DiagonalMatrix[Sqrt[Range[39]], 1]},
 Max[Abs[Normal[d["Matrix"]][[;; 5, ;; 5]] - MatrixExp[0.4 (Transpose[low] - low)][[;; 5, ;; 5]]]]]
```

## Notes and limits

Results are written with [Sinc]() so they stay finite where a frequency vanishes; [FunctionExpand]() turns them into sines and cosines. With machine numbers any number of coupled modes orders quickly. With exact numbers, or a symbolic time, two coupled modes order in a few seconds, while three coupled modes bring the roots of a cubic into every coefficient and run past the default TimeConstraint, so give the coefficients as machine numbers there. Symbolic coefficients work for one mode, the beam splitter and the two-mode squeezer with frequencies. For general coupled modes the continuous branch of the root depends on which stability region the coefficients are in, and symbolic coefficients alone do not decide that.

## References

- R. M. Wilcox, "Exponential operators and parameter differentiation in quantum physics", J. Math. Phys. 8, 962 (1967).
- K. E. Cahill and R. J. Glauber, "Ordered expansions in boson amplitude operators", Phys. Rev. 177, 1857 (1969).
- D. R. Truax, "Baker-Campbell-Hausdorff relations and unitarity of SU(2) and SU(1,1) squeeze operators", Phys. Rev. D 31, 1988 (1985).
- X. Ma and W. Rhodes, "Multimode squeeze operators and squeezed states", Phys. Rev. A 41, 4625 (1990).
- P. Blasiak, A. Horzela, K. A. Penson, A. I. Solomon and G. H. E. Duchamp, "Combinatorics and boson normal ordering: a gentle introduction", Am. J. Phys. 75, 639 (2007).
