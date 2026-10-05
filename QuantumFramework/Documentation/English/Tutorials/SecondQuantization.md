---
Template: TechNote
Name: SecondQuantization
Title: Second Quantization Functions
Context: Wolfram`QuantumFramework`
CellContext: Global`
Paclet: Wolfram/QuantumFramework
URI: Wolfram/QuantumFramework/tutorial/SecondQuantization
Keywords: [Second Quantization, Fock Space, Quantum Optics]
RelatedGuides: [WolframQuantumComputationFramework]
RelatedTutorials: [TransmonGates]
Typeset: _SuperDagger -> StandardForm
---

This tech note introduces the QuantumFramework implementation of bosonic second quantization on a truncated Fock space. The truncation provides a finite-dimensional representation that is convenient for computation while retaining the structure of common states and operators. The examples below focus on practical workflows for building states, applying operators, and visualizing results in quantum mechanics and quantum optics.

Install and load the QuantumFramework paclet and the second-quantization subpackage (skip the install step if it is already available):

```wl
#| eval: false
PacletInstall["https://www.wolfr.am/DevWQCF", 
 ForceVersionInstall -> True]
```

```wl
<< Wolfram`QuantumFramework`
<< Wolfram`QuantumFramework`SecondQuantization`
```

<!-- #| tags: core -->
## Bosonic states and operators

### Fock space size and utilities

<!-- #| style: DefinitionBox -->
|   |   |
|---|---|
| [SetFockSpaceSize]()[*size*] | Sets *size *as the default truncation dimension , all the other definitions can use this default *size* |
| [OperatorVariance]()[ψ,*op*] | Variance of the operator *op *for a state ψ |

### States

<!-- #| style: DefinitionBox -->
|   |   |
|---|---|
| [FockState]()[*n*,*size*] | *n*-th Fock state in a basis of dimension *size* |
| [FockState]()[*{$n_1$,$n_2$,..}*,*size*] | Multimode Fock state with occupation numbers $n_{i}$ , with each mode in a basis of dimension *size* |
| [CoherentState]()[*size*] [CoherentState]()[*size*][α] | Parametric state representing a normalized coherent state in a basis of dimension *size *, the complex amplitude is a parameter |
| [ThermalState]()[*nbar*,*size*] | Thermal mixed state with the average number of photons *nbar* in a basis of dimension *size* |
| [CatState]()[*size*] [CatState]()[*size*][α,ϕ] | Parametric general cat state representing a superposition of coherent states in a basis of dimension *size*, the complex amplitude and the phase are parameters |

### Operators

<!-- #| style: DefinitionBox -->
|   |   |
|---|---|
| [AnnihilationOperator]()[*size, order*] | Bosonic annihilation operator in a truncated Fock space with basis dimension *size* and qudit order *order* |
| [DisplacementOperator]()[α,*size*, *order*] | Phase space displacement operator of a single mode with complex amplitude *alpha* , basis dimension *size* and qudit order *order* |
| [SqueezeOperator]()[ξ,*size*,*order*] | Squeeze operator of a single mode with complex parameter *xi* , basis dimension *size* and qudit order *order* |
| [QuadratureOperators]()[*size,*order**] | $X_{1}$ and $X_{2}$ position and momentum quadrature operators with basis dimension *size* and qudit order *order* |
| [PhaseShiftOperator]()[*θ,size,*order**] | Phase space rotation operator with rotation angle *θ*, basis dimension *size* and qudit order *order* |
| [BeamSplitterOperator]()[*{θ,ϕ},size,*order**] | Two mode beam-splitter operator, *θ* indicates the reflectivity, *ϕ* the relative phase in a basis of dimension *size* and qudit order *order* |

In the sections below, we demonstrate how to use these functions in the QuantumFramework.

<!-- #| tags: details -->
## Details and options

### Operators Mathematical Definitions

<!-- #| style: DefinitionBox -->
|   |   |
|---|---|
| `DisplacementOperator` | $\mathrm{Exp}[\alpha \, a^{\dagger }\, -\alpha^{*}a]$ |
| `SqueezeOperator` | $\exp[\tfrac{1}{2}(\xi^{*} a^{2} - \xi\, a^{\dagger 2})]$ |
| `PhaseShiftOperator` | [Exp]()[i θ $a^{\dagger }$a] |
| `BeamSplitterOperator` | [Exp]()[θ($e^{i \phi }$$a_{1}$$a_{2}^{\dagger }$-$e^{-i \phi }$$a_{1}^{\dagger }$$a_{2}$)] |
| `QuadratureOperators` | $X_1$:$1/2$(a+$a^\dagger $) $X_2$:$1/2\, i$(a-$a^\dagger $) |

#### Displacement operator ordering

By default, when "Ordering" is not specified, normal ordering is used. The formulas for the three orderings are shown below:

<!-- #| style: DefinitionBox -->
|   |   |
|---|---|
| [DisplacementOperator]()[..,"Ordering"->*st*] | *st* can be "Normal", "Weak" or "Antinormal". For "Normal", one gets $e^{-|\alpha|^2/2}\, e^{\alpha a^\dagger}\, e^{-\alpha^* a}$; for "Weak", $e^{\alpha a^\dagger - \alpha^* a}$; for "Antinormal", $e^{|\alpha|^2/2}\, e^{-\alpha^* a}\, e^{\alpha a^\dagger}$. |

#### Squeeze operator ordering

By default, when "Ordering" is not specified, normal ordering is used. The formulas for the three orderings are shown below:

<!-- #| style: DefinitionBox -->
|   |   |
|---|---|
| [SqueezeOperator]()[…,"Ordering"->*st*] | *st* can be "Normal", "Antinormal" or "Weak". Given $\xi = r e^{i\theta}$, for "Normal", one gets $e^{-\frac{1}{2} e^{i\theta} \tanh r\, a^{\dagger 2}}\, e^{-\ln(\cosh r)\,(a^\dagger a + \frac{1}{2})}\, e^{\frac{1}{2} e^{-i\theta} \tanh r\, a^{2}}$; for "Weak", $e^{\frac{1}{2}(\xi^* a^{2} - \xi a^{\dagger 2})}$; for "Antinormal", $e^{\frac{1}{2} e^{-i\theta} \tanh r\, a^{2}}\, e^{\ln(\cosh r)\,(a^\dagger a + \frac{1}{2})}\, e^{-\frac{1}{2} e^{i\theta} \tanh r\, a^{\dagger 2}}$. |

#### Coherence and correlation functions

<!-- #| style: DefinitionBox -->
|   |   |
|---|---|
| [G2Coherence]()[ψ] | Second-order coherence $g^{(2)}$(0) for state ψ |
| [G1Correlation]()[ψ,{{$r_{1}$,$t_{1}$},{$r_{2}$,$t_{2}$}}] | First-order correlation function $G^{(1)}(x_{1},x_{2})$ for state ψ at the space and time coordinates $r_{1},t_{1},r_{2},t_{2}$ |

<!-- #| tags: basic1 -->
## Basic usage

### Setting the space size

The global variable `$FockSize` stores the default truncation size (dimension) of the Fock space. It is used by all second-quantization functions when the *size* parameter is not explicitly specified. • Default value: 16 • Use $\mathrm{SetFockSpaceSize}[n]$ to change the default. • Use $\mathrm{SetFockSpaceSize}[]$ to reset to 16.

**Note:** Larger truncation sizes improve accuracy for high photon-number states, but increase computational cost.

Show the default truncation size:

```wl
$FockSize
```

Set the truncation size to 20:

```wl
SetFockSpaceSize[20];
```

Show new truncation size:

```wl
$FockSize
```

Set the default size again:

```wl
SetFockSpaceSize[];
```

Show the current truncation size:

```wl
$FockSize
```

### Fock States

Create a single-mode Fock state with occupation number 5 using the default truncation size:

```wl
FockState[5]
```

Create a single-mode Fock state with occupation number 5 using truncation size 10:

```wl
FockState[5, 10]
```

Create a two-mode Fock state with occupation numbers 2 and 12:

```wl
FockState[{2, 12}]
```

Verify the shorthand string representation matches the two-mode state:

```wl
FockState[{6, 6}] == QuantumState["66", $FockSize]
```

This shorthand does not work for numbers ≥ 10 because the digits are interpreted as separate qudits:

```wl
FockState[{10, 10}]
QuantumState["1010", $FockSize]
% == %%
```

### Annihilation and creation operators

Create a single-mode annihilation operator with the default truncation size:

```wl
AnnihilationOperator[]
```

Create an annihilation operator with truncation size 8 (an 8×8 matrix) using the default mode order {1}:

```wl
AnnihilationOperator[8]
```

Create an annihilation operator acting on mode 2 using the default truncation size:

```wl
AnnihilationOperator[{2}]
```

Create an annihilation operator with size 5 acting on mode 3:

```wl
AnnihilationOperator[5, {3}]
```

Apply the annihilation operator to a Fock state:

```wl
AnnihilationOperator[][FockState[2]]["Formula"]
```

Define the creation operator by taking the adjoint (SuperDagger) of the annihilation operator, then apply it to a Fock state:

```wl
(SuperDagger[AnnihilationOperator[]]@ FockState[1])["Formula"]
```

Or use the "Dagger" property directly:

```wl
(SuperDagger[AnnihilationOperator[]]@FockState[1])["Formula"]
```

### Coherent state

Create a symbolic coherent state with the default truncation size:

```wl
CoherentState[][\[Beta]]
```

Create a symbolic coherent state with truncation size `20`, it is possible to specify that the state is not normalized to avoid the truncated space normalization factor (can be impractical for symbolic manipulations):

```wl
CoherentState[20, "Normalized" -> False][\[Beta]]["Formula"]
```

Photon-number distribution for a coherent state with $\alpha \, =\, 1.3+2i$ using truncation size `25`:

```wl
CoherentState[25][1.3 + 2 I]["ProbabilityPlot", LabelStyle -> 8, 
 AspectRatio -> 1/4]
```

### Cat State

Create a cat state $\mathcal{N}_{\alpha }(|\alpha \rangle +e^{i\phi }|-\alpha \rangle )$ with amplitude $\alpha =2+i$ and phase $\phi =3\pi /4$:

```wl
CatState[25][2 + I, 3 \[Pi]/4]
```

Photon distribution of the same state:

```wl
CatState[25][2. + I, 3 \[Pi]/4]["ProbabilityPlot", LabelStyle -> 6, 
 AspectRatio -> 1/4]
```

### Thermal state

Create a thermal state with mean photon number $\bar{n}=2$ using truncation size `40`:

```wl
ThermalState[2., 40]
```

Plot the photon-number distribution (diagonal of the density matrix), for $\bar{n}=4$ and truncation dimension `40`:

```wl
ListPlot[Diagonal[ThermalState[4, 40]["DensityMatrix"]], 
 PlotRange -> All, Filling -> Axis, 
 AxesLabel -> {"Fock States", "Probability"}]
```

Define a parametric thermal state in terms of the mean photon number:

```wl
state = QuantumState[ThermalState[x, 40], "Parameters" -> x]
```

Plot the entropy as a function of the mean photon number:

```wl
Plot[state[x]["Entropy"], {x, 0, 5}, 
 AxesLabel -> {"\!\(\*OverscriptBox[\(n\), \(\[LongDash]\)]\)", 
   "Entropy"}]
```

### Displacement operator

Create a displacement operator with truncation size 40 and amplitude 3+i:

```wl
DisplacementOperator[3 + I, 40]
```

Verify that displacing the vacuum produces a coherent state:

```wl
DisplacementOperator[\[Alpha], 40][] == CoherentState[40][\[Alpha]]
```

Specify the mode (order) of the operator with the default size, mode 2 in this case:

```wl
DisplacementOperator[0.5, {2}]
```

### Phase shift operator

Matrix form of the phase-shift operator for truncation size 10:

```wl
PhaseShiftOperator[\[Theta], 10]["Matrix"] // MatrixForm
```

Effect of phase shift on a single-mode symbolic state:

```wl
PhaseShiftOperator[\[Theta]][
  QuantumState[
   Array[Subscript[\[FormalC], #] &, $FockSize, 
    0], $FockSize]] // TraditionalForm
```

Effect of phase shift on a two-mode superposition state:

```wl
PhaseShiftOperator[\[Theta], {2}][
  FockState[{5, 3}] + FockState[{5, 6}]]["Formula"]
```

Phase shift operator coincides with the "Z" gate for qudits when the angle is 2π/n (up to conjugation):

```wl
With[{n = RandomInteger[{1, $FockSize - 1}]}, 
 PhaseShiftOperator[2 \[Pi]/n, n] == 
  QuantumOperator["Z"[n]]["Conjugate"]]
```

### Squeeze operator

Create a squeeze operator with complex parameter 1+i using the default size:

```wl
SqueezeOperator[1. + I]
```

Create a two-mode state with mode 2 in a squeezed vacuum |0〉 ⊗ |ξ〉 ($\xi \, =\, 0.9$) using truncation size 15:

```wl
state = SqueezeOperator[0.9, 15, {2}]@FockState[{0, 0}, 15]
```

Plot the non-zero amplitudes of the state:

```wl
state["AmplitudePlot", ChartLegends -> Automatic, ImageSize -> 200]
```

### Beam-Splitter Operator

Set a symmetric beam splitter with $\theta \, =\, \pi /4$ and $\phi \, =\, \pi /2$:

```wl
symmetricBS = BeamSplitterOperator[{\[Pi]/4., \[Pi]/2.}]
```

Transform the state |1,1〉 with the symmetric beam splitter:

```wl
Chop[symmetricBS[FockState[{1, 1}]]]["Formula"]
```

Specify the mode order (where the operator acts) and the truncation size:

```wl
BeamSplitterOperator[{0.2, 0.9}, 6, {1, 3}]
```

#### BeamSplitter operator: Method option

BeamSplitterOperator supports two computation methods via the Method option. The default is *MatrixExp *which uses direct matrix exponentiation; *Recurrence *uses an efficient recurrence relation for the beam splitter matrix elements. Performance can vary depending on whether the angle arguments are numeric or exact/symbolic.

Use the "Recurrence" method with exact parameters:

```wl
BeamSplitterOperator[{\[Pi]/4, \[Pi]/2}, Method -> "Recurrence"][
  FockState[{1, 1}]]["Formula"]
```

### OperatorVariance

Calculate the variance of the quadrature operators for the ground state:

```wl
OperatorVariance[FockState[0], #] & /@ QuadratureOperators[]
```

A Fock state has a definite number of photons, so the variance of the number operator is zero:

```wl
With[{a = AnnihilationOperator[]}, 
 OperatorVariance[FockState[2], SuperDagger[a]@a]]
```

### G2Coherence (Second-Order Coherence)

The function G2Coherence computes the equal-time second-order coherence $g^{(2)}(0)$ which characterizes the photon statistics of a quantum state. It is defined as:

$g^{(2)}(0)=\frac{\langle a^{\dagger }a^{\dagger }a\, a\rangle }{\langle a^{\dagger }a\rangle^{2}}$ Physical interpretation: • $g^{(2)}=1$: Poissonian statistics (coherent light) • $g^{(2)}<1$: Sub-Poissonian (antibunching, nonclassical light) • $g^{(2)}>1$: Super-Poissonian (bunching, thermal light) • $g^{(2)}=0$: Single-photon state

For a coherent state, $g^{(2)}=1$:

```wl
G2Coherence[CoherentState[][2.]]
```

For a Fock state |n〉 with n > 0, $g^{(2)}=1-1/n$:

```wl
G2Coherence[FockState[3]]
```

For a thermal state with mean photon number $\overline{n}$, $g^{(2)}=2$:

```wl
G2Coherence[ThermalState[2.]]
```

Improve the result increasing the truncation size:

```wl
G2Coherence[ThermalState[2., 40]]
```

It is not exactly 2 because of truncation: the thermal distribution has some weight at photon numbers beyond the cutoff.

The single photon state |1〉 has $g^{(2)}=0$ (perfect antibunching):

```wl
G2Coherence[FockState[1]]
```

Consider the mixed state (1-ϵ) |0〉 〈0| + ϵ |2〉 〈2|:

```wl
\[Rho] = (1 - \[Epsilon]) FockState[0][
    "MatrixState"] + \[Epsilon] FockState[2]["MatrixState"]
```

Extreme bunching when ϵ -> 0:

```wl
G2Coherence[\[Rho]]
```

If your detector clicks, it means you caught a photon. Because photons only exist in pairs in this specific state, the probability that a second photon is right there with it is virtually 100%.Therefore, $g^{2}$ diverges because the conditional probability of detecting a second photon (given that you just detected a first) is massively disproportionate to the absolute probability of detecting a photon randomly in the dark field (ϵ ->0)

### G1Correlation (first - order correlation)

The first-order correlation is defined for some state |ψ〉 as

$G^{(1)}(x_{1},x_{2})=\langle \, \psi \, |\, E^{(-)}(r_{1},t_{1})\, E^{(+)}(r_{2},t_{2})\, |\psi \rangle $

where the $E^{(+)},E^{(-)}$ are the positive and negative frequency components of a single-mode field

Calculate $G^{(1)}(x,x)$ for |4〉:

```wl
G1Correlation[FockState[4], {{r, t}, {r, t}}]
```

Calculate $G^{(1)}(x_{1},x_{2})$ for |4〉:

```wl
G1Correlation[FockState[4], {{r1, t1}, {r2, t2}}]
```

### Quadrature Operators

Compute the commutator of the quadrature operators. It should be $\frac{i}{2}1$; any small deviation in the last term is expected due to truncation:

```wl
(Commutator @@ QuadratureOperators[]) // TraditionalForm
```

Expectation value of $\langle X_{1}^{2}\rangle $ for Fock state $|3\rangle $:

```wl
With[{s = FockState[3]}, 
  SuperDagger[s]@(QuadratureOperators[][[1]]^2)@s]["Scalar"]
```

Or alternatively because the operator is Hermitian:

```wl
(QuadratureOperators[][[1]]@FockState[3])["Norm"]^2
```

### Mode Order (Multimode Systems)

Many operators accept an optional *order* parameter that specifies which mode(s) the operator acts on in a multimode system. The order is a list of positive integers. • Default: `{1}` (single mode, acts on mode 1) • Two modes: `{1,2}` or `{2,3}` etc. The order parameter is essential for: • Creating operators that act on specific modes in a tensor-product space. • Building multimode quantum circuits (e.g., beam splitters coupling modes 1 and 2).

Single-mode operator with explicit order:

```wl
AnnihilationOperator[{2}]
```

$\langle n_{1},n_{2}|a_{2}|n_{1},n_{2}\rangle =\sqrt{n_{2}}|n_{1},n_{2}-1\rangle $

Verify the property:

```wl
AnnihilationOperator[{2}][FockState[{3, 5}]] == 
 Sqrt[5] FockState[{3, 5 - 1}]
```

Two-mode beam splitter acting on modes 1 and 3:

```wl
BeamSplitterOperator[{\[Pi]/4, \[Pi]/2}, {1, 3}, 
   Method -> "Recurrence"][FockState[{1, 3, 1}]]["Formula"]
```

Displacement operator on mode 3:

```wl
DisplacementOperator[1.5, 8, {3}]
```

Note: When using mode order, all operators in a calculation should use consistent mode assignments. The *order* parameter does NOT change the Fock space dimension—it only labels which mode the operator acts on.

### Error Messages

The following error messages may be generated by second-quantization functions:

<!-- #| style: DefinitionBox -->
|   |   |
|---|---|
| `FockState::clip` | Issued when the Fock state index is outside the valid range [0, size-1]. The index is clipped to the nearest valid value. |
| `FockVals::len` | Issued when multimode Fock state values contain invalid (negative or too large) integers. |
| `DisplacementOperator::invalidorder` | Issued when the "Ordering" option is not "Normal", "Weak", or "Antinormal". |
| `SqueezeOperator::invalidorder` | Issued when the "Ordering" option is not "Normal", "Weak", or "Antinormal". |
| `BeamSplitterOperator::badmethod` | Issued when the Method option is not "MatrixExp" or "Recurrence". |

Example: FockState index clipping

```wl
FockState[100, 16]
```

### Ordering option

Displacement and squeeze operators can be defined with different orderings. Use the "Ordering" option to choose "Weak", "Normal", or "Antinormal". In the infinite-dimensional limit these are equivalent, but in a truncated space the numerical error depends on the ordering, so one choice may be more accurate for a given calculation.

Check unitarity of the displacement operator using "UnitaryQ" for $\alpha \, =\, 0.5+i$. For this case, weak ordering gives the expected result:

```wl
DisplacementOperator[0.5 + I]["UnitaryQ"]
```

Repeat the check using weak ordering:

```wl
DisplacementOperator[0.5 + I, "Ordering" -> "Weak"]["UnitaryQ"]
```

For a squeezed vacuum state with $\xi \, =\, 1.5\, +0.5\, i$ and `size = 20` compare the numerical error of "Weak" and "Normal" orderings against the analytic expression. For this value of ξ, the truncation size is not sufficient for the weak-ordered operator, and the normal-ordered one is more accurate.

Amplitudes obtained from a known analytic formula:

```wl
analytic = 
  1/Sqrt[Cosh[Abs[\[Xi]]]] Table[
     Sqrt[(2 m)!]/(2^m m!) (-1)^m E^(I m Arg[\[Xi]]) Tanh[Abs[\[Xi]]]^
      m, {m, 0, 9}] /. \[Xi] -> (1.5 + 0.5 I);
```

Compute absolute errors:

```wl
absoluteErrors = (SqueezeOperator[1.5 + 0.5 I, 20, "Ordering" -> #][][
       "AmplitudeList"] - analytic) & /@ {"Normal", "Weak"};
```

Show the results:

```wl
TableForm[absoluteErrors // Transpose, 
 TableHeadings -> {Table[Ket[{i}], {i, 0, 20, 2}], {"Normal", "Weak"}}]
```

### Applications: States and operators

### Mean number of photons

Define the annihilation operator and create a random pure state:

```wl
a := AnnihilationOperator[]
\[Psi] = QuantumState["RandomPure", $FockSize];
```

Compute the mean number of photons using QuantumMeasurementOperator:

```wl
(QuantumMeasurementOperator[SuperDagger[a]@a][\[Psi]])["Mean"]
```

Equivalent bracket notation $\langle \, \psi \, |a^{\dagger }a|\psi \rangle \, $:

```wl
Chop[( SuperDagger[\[Psi]]@(SuperDagger[a]@a)@\[Psi] )["Scalar"]]
```

Number operator is Hermitian:

```wl
a[\[Psi]]["Norm"]^2
```

Using the probabilities of the diagonal of the density matrix:

```wl
Chop@Dot[Range[0, $FockSize - 1], Diagonal@\[Psi]["DensityMatrix"]]
```

### Heisenberg evolution of annihilation operator 

In the Heisenberg picture, the annihilation operator evolves as $\hat{a}(t)=\hat{a}\, e^{-i\, \omega \, t}$. We can reproduce this using super-operators in the framework.

Define the Hamiltonian:

```wl
H = \[HBar] \[Omega] (SuperDagger[a]@a + 1/2);
```

Construct the evolution super-operator:

```wl
exp\[ScriptCapitalH] = 
  MatrixExp[-I QuantumOperator["Hamiltonian"[H/\[HBar]]] t];
```

Evolve the operator:

```wl
heisenbergAnnihilation = 
  Simplify[
   SuperDagger[exp\[ScriptCapitalH]][a["MatrixQuantumState"]][
    "Operator"], {t, \[Omega]} \[Element] Reals];
```

Compare the result to the known result of quantum mechanics:

```wl
a E^(-I \[Omega] t) == heisenbergAnnihilation
```

<!-- #| tags: quasi -->
## Quasi-probability Representations 

Two phase-space quasi-probability representations are included: Wigner and HusimiQ. They are used extensively in quantum mechanics and quantum optics. Each is computed numerically on a grid of phase-space points and returns an interpolating function. Glauber representation is not included because, for a general state, a stable numerical algorithm is not feasible due to highly singular behavior.

<!-- #| style: DefinitionBox -->
|   |   |
|---|---|
| [WignerRepresentation]()[ψ,{*xmin*,*xmax*},{*pmin*,*pmax*},*opts*] | Computes the Wigner quasi-probabilty representation of a quantum state ψ , in a phase space rectangular region defined by the given limits of *x* and *p* |
| [HusimiQRepresentation]()[ψ,{*xmin*,*xmax*},{*pmin*,*pmax*},*opts*] | Computes the Husimi Q function representation of a quantum state *ψ *, in a phase space rectangular region defined by the given limits of *x* and *p* |

In the sections below, we demonstrate how to use these functions in the QuantumFramework.

### Examples

#### Ground state Wigner Representation:

Compute the Wigner representation on the square region [-4, 4] × [-4, 4]:

```wl
wig = WignerRepresentation[FockState[0], {-4, 4}, {-4, 4}];
```

Visualize the function in 3D:

```wl
Plot3D[
 wig[x, p], {x, -4, 4}, {p, -4, 4},
 ColorFunction -> "Rainbow",
 PlotRange -> All,
 PlotLegends -> Automatic,
 PlotPoints -> 80
 ]
```

#### HusimiQ Representation of a Fock state

Compute the HusimiQ representation for a Fock state with `n = 4`:

```wl
husi = HusimiQRepresentation[FockState[4], {-4, 4}, {-4, 4}];
```

Visualize the density plot of the representation:

```wl
DensityPlot[
 husi[x, p], {x, -4, 4}, {p, -4, 4},
 ColorFunction -> "Rainbow",
 PlotRange -> All,
 PlotLegends -> Automatic,
 PlotPoints -> 80,
 FrameLabel -> {"x", "p"}
 ]
```

## Quasi-probability representations: Applications and options

### Applications

### Expectation values

Compute the Wigner function for the Fock state |2〉:

```wl
wig = WignerRepresentation[FockState[2], {-5, 5}, {-5, 5}];
```

Phase space integral (normalization):

```wl
NIntegrate[wig[x, p], {x, -4, 4}, {p, -4, 4}, AccuracyGoal -> 4]
```

Compute the expectation value $\langle x^{4}\rangle $, for a purely position operator $x^{4}$ its Weyl symbol is $x^{4}$, hence the expectation value is $\int \, \int x^{4}W_{2}(x,p)\, dx\, dp$:

```wl
NIntegrate[x^4 wig[x, p], {x, -4, 4}, {p, -4, 4}, 
 Method -> "MultiPeriodic"]
```

Alternatively using $|x^{2}|2\rangle |^{2}$ :

```wl
(((Sqrt[2] QuadratureOperators[][[1]])^2)@FockState[2])["Norm"]^2 // N
```

#### Marginal distributions

Compute the Wigner function for a coherent state with $\alpha =1$:

```wl
wig = WignerRepresentation[CoherentState[][1.], {-4, 4}, {-4, 4}]
```

Visualize momentum distribution:

```wl
Plot[Evaluate[Integrate[wig[x, p], {x, -4, 4}]], {p, -4, 4}]
```

Visualize position distribution:

```wl
Plot[Evaluate[Integrate[wig[x, p], {p, -4, 4}]], {x, -4, 4}]
```

### Squeezed vacuum Wigner Representation

Wigner representation of a squeezed vacuum state with parameter $\xi \, =\, 0.5$:

```wl
squeezedVacuumFunc = 
  WignerRepresentation[SqueezeOperator[0.5][], {-4, 4}, {-4, 4}];
```

Plot the function, the nonclassical quadrature behavior of squeezed states is visible:

```wl
Plot3D[
 squeezedVacuumFunc[x, p], {x, -4, 4}, {p, -4, 4},
 ColorFunction -> "SunsetColors",
 PlotRange -> All,
 PlotLegends -> Automatic,
 PlotPoints -> 80
 ]
```

### Options

<!-- #| style: DefinitionBox3Col -->
|   |   |   |
|---|---|---|
| `"GaussianScaling"` | $\sqrt{2}$ | Scaling factor of the representation |
| `"GridSize"` | `100` | Granularity of the sampling |

#### Scaling of the representation for ground state

Set the Gaussian scaling parameter to 0.5:

```wl
wig = WignerRepresentation[FockState[0], {-4, 4}, {-4, 4}, 
   "GaussianScaling" -> 0.5];
```

Visualize the result:

```wl
Plot3D[wig[x, p], {x, -4, 4}, {p, -4, 4},
 ColorFunction -> "Rainbow",
 PlotRange -> All,
 PlotPoints -> 80,
 PlotLegends -> Automatic,
 ImageSize -> 200]
```

#### Using GridSize

Compute the Wigner function value at the origin for the Fock state |2〉 using a grid size of 140:

```wl
wig = WignerRepresentation[FockState[2], {-6, 6}, {-6, 6}, 
   "GridSize" -> 140];
wig[0, 0]
```

A grid size of 20 is insufficient for this region and yields noticeable numerical error:

```wl
wig2 = WignerRepresentation[FockState[2], {-6, 6}, {-6, 6}, 
   "GridSize" -> 20];
wig2[0, 0]
```

The exact value from the analytic expression is:

```wl
exact = 1/\[Pi] E^-(x^2 + p^2) LaguerreL[2, 2 (x^2 + p^2)] /. {x -> 
     0 \.08, p -> 0};
```

Compare absolute errors:

```wl
{exact - wig[0, 0], exact - wig2[0, 0]}
```

<!-- #| tags: advanced -->
## Advanced examples

### Decay of a coherent state

Consider $\mathcal{H}=\, \omega \, \hat{a}^{\dagger }\hat{a}$ with Lindblad jump operator $L_{1}\, =\, \sqrt{\gamma }\hat{a}$ (single-photon loss channel), and set $\omega \, =\, 1$

Define the annihilation operator and use a coherent state with $|\alpha |\, =\, 2$ as the initial state:

```wl
a := AnnihilationOperator[];
SetFockSpaceSize[15];
\[Alpha] = 2 Exp[I RandomReal[{0, 2 Pi}]];
\[Rho]0 = CoherentState[][\[Alpha]];
```

Set the decay rate γ:

```wl
\[Gamma] = RandomReal[{1, 2}];
```

Define the Hamiltonian including photon loss:

```wl
\[ScriptCapitalH] = 
  QuantumOperator["Hamiltonian"[SuperDagger[a]@a, {a}, {\[Gamma]}]];
```

State evolution, we can do it symbolically:

```wl
\[Rho]t = QuantumEvolve[\[ScriptCapitalH], \[Rho]0, t];
```

Plot the mean photon number as a function of time:

```wl
meanPhotonPlot = 
  Plot[Chop[
    Range[0, $FockSize - 1] . 
     Diagonal[\[Rho]t[t]["DensityMatrix"]]], {t, 0, 1}, Frame -> True, 
   FrameLabel -> {"t", "Mean photon number"}, GridLines -> Automatic, 
   LabelStyle -> 10, ImageSize -> 220, PlotRange -> {0, 4}];
```

Plot the fidelity with respect to the analytic expression:

```wl
fidelityPlot = 
  Plot[QuantumSimilarity[CoherentState[][\[Alpha] \!\(TraditionalForm\`
\*SuperscriptBox[\(E\), \(\(-I\)\ t\)]\ 
\*SuperscriptBox[\(E\), \(\(-
\*FractionBox[\(1\), \(2\)]\)\ t\ \[Gamma]\)]\)], \[Rho]t[t]], {t, 0, 
    1}, Frame -> True, FrameLabel -> {"t", "State similarity"}, 
   GridLines -> Automatic, LabelStyle -> 10, PlotRange -> {0, 1.1}, 
   ImageSize -> 220];
```

Show the results:

```wl
Labeled[GraphicsRow[{meanPhotonPlot, fidelityPlot}], 
 Style["Evolution statistics of a coherent state", 
  FontFamily -> "Arial"], Top]
```

### Optical balance

Consider the following master equation, which balances a driving field Hamiltonian and single-photon loss:

$\frac{d\, \rho }{dt}=\, -ig[\hat{a}+\hat{a}^{\dagger },\hat{\rho }]+\frac{\gamma }{2}(2\hat{a}\hat{\rho }\hat{a}^{\dagger }-\hat{a}^{\dagger }\hat{a}\hat{\rho }-\hat{\rho }\hat{a}^{\dagger }\hat{a})$

Set the drive strength g:

```wl
g = RandomReal[{1, 2}];
```

Random initial state:

```wl
\[Rho]0 = QuantumState["RandomPure", $FockSize];
```

Set up the master equation:

```wl
\[ScriptCapitalH] = 
  QuantumOperator["Hamiltonian"[g (a + SuperDagger[a]), {a}, {\[Gamma]}]];
```

Numerical evolution from $t_{i}=0\, $ to $t_{f}=3$ :

```wl
\[Rho]t = QuantumEvolve[ \[ScriptCapitalH], \[Rho]0, {t, 0, 3}];
```

Steady state solution is a coherent state with amplitude $\alpha \, =\, -2i\, g/\gamma $ via phase-space integration:

```wl
Plot[QuantumSimilarity[\[Rho]t[t], 
  CoherentState[][-I 2 g/\[Gamma]]], {t, 0, 3}, PlotRange -> {0, 1}, 
 Frame -> True, 
 PlotLabel -> 
  "Convergence to a coherent state |-2\[ImaginaryI]G/\[Gamma]\[RightAngleBracket]", 
 FrameLabel -> {"Time(t)", "Fidelity"}, ImageSize -> 300]
```

For the initial state $\hat{\rho }$(0) = |0〉 〈0|, the system evolves to the coherent state $|{-i g t}\rangle\langle {-i g t}|$ for short times; the complex amplitude is independent of γ.

Initial state:

```wl
\[Rho]0 = QuantumState["0", $FockSize];
```

Evolved state:

```wl
\[Rho]t = QuantumEvolve[\[ScriptCapitalH], \[Rho]0, {t, 0, 1}];
```

The integrator leaves the evolved density matrix with eigenvalues a little below zero, so the fidelity is taken with its repaired form, the "Physical" property. Photon number expectations and fidelity:

```wl
time = Range[0, 0.8, 0.02];
qs = QuantumSimilarity[CoherentState[][-I g #], \[Rho]t[#]["Physical"]] & /@ time;
meanN = Chop[
     Range[0, $FockSize - 1] . 
      Diagonal[\[Rho]t[#]["DensityMatrix"]]] & /@ time ;
```

Plot both quantities versus time:

```wl
ListLinePlot[{Thread[{time, qs}], Thread[{time, meanN}]}, 
 PlotRange -> {0, 1.1}, Frame -> True, FrameLabel -> {"t (Time)", ""}, 
 LabelStyle -> 10, 
 PlotLegends -> {"Fidelity", "Mean number of photons"}, 
 PlotLabel -> "Optical balance for short times", ImageSize -> 300]
```

### Jaynes-Cummings Model: Calculating the unitary

```wl
ClearAll[g]
```

Define the basis and the relevant operators:

```wl
SetFockSpaceSize[32];
atomBasis = QuantumBasis[{"e", "g"}];
```

Define the atomic operators and the annihilation operator in second mode:

```wl
\[Sigma]plus =  QuantumOperator["+", atomBasis];
\[Sigma]minus = QuantumOperator["-", atomBasis];
a := AnnihilationOperator[{2}];
```

Free Hamiltonian $\mathcal{H}_{0}=\, \hbar \nu \hat{a}^{\dagger }\hat{a}+\frac{1}{2}\hbar \omega \sigma_{z}$:

```wl
H0 = QuantumOperator[
   h \[Nu] SuperDagger[a]@a + 
    1/2 h \[Omega] QuantumOperator["Z", atomBasis], 
   "Parameters" -> {\[Nu], \[Omega]}];
```

Atom-field interaction term $\mathcal{H}_{I}=\, \hbar g(\sigma_{+}\hat{a}+\sigma_{-}\hat{a}^{\dagger })$:

```wl
H1 = h g (\[Sigma]plus@a + \[Sigma]minus@SuperDagger[a]);
```

Going to the interaction picture at resonance, we get a time-independent Hamiltonian and from it we get the unitary operator to evolve states:

```wl
HInteract = MatrixExp[I H0/h t]@H1@MatrixExp[-I H0/h t];
U = MatrixExp[-I HInteract[\[Nu], \[Nu]]/h t];
```

### Jaynes-Cummings Model: Evolving states

#### Field in a Fock State

Rabi oscillations when the field is initially in a Fock state

```wl
\[Psi]0 = 
  QuantumTensorProduct[QuantumState["1", atomBasis], FockState[1]];
\[Psi]t = U @ \[Psi]0
```

Find the probability of the atom being in the excited state, using the projector as a measurement operator:

```wl
proj = QuantumState["0", atomBasis]["Operator"];
probE = (QuantumMeasurementOperator[proj]@ \[Psi]t)["Mean"];
```

Plot the result:

```wl
Plot[probE /. g -> 1 , {t, 0, 3 \[Pi]}, Frame -> True, 
 FrameLabel -> {"T", "Probability"}, LabelStyle -> 10, 
 PlotLabel -> "Resonant Rabi oscillations", ImageSize -> 300]
```

Calculate the evolution for the states |g〉 ⊗ |n〉, varying the number state n:

```wl
TableForm[
  Table[s = 
    QuantumTensorProduct[QuantumState["1", atomBasis], 
     FockState[x]]; {s["Formula"], U@s}, {x, 0, 10}], 
  TableHeadings -> {None, {"Initial State", 
     "Evolved State"}}] // TraditionalForm
```

#### Coherent state: Collapse and revivals

Define the initial state $|\psi_{0}\rangle \, =\, |g\rangle |\alpha \rangle \, with\, |\alpha |=4$ and evolve it (we use a random α):

```wl
\[Alpha] = 4 Exp[I RandomReal[{0, 2 \[Pi]}]];
\[Psi]0 = 
  QuantumTensorProduct[QuantumState["1", atomBasis], 
   CoherentState[][\[Alpha]]];
\[Psi]t = QuantumState[U @ \[Psi]0, "Parameters" -> {g, t}];
```

Extracting the coefficients $c_{e}$ and $c_{g}$ for all n as a function of time:

```wl
{cet[t_], cgt[t_]} = (SuperDagger[#]@\[Psi]t[1, t])["AmplitudesList"] & /@ 
   atomBasis["BasisStates"];
```

Plotting the population inversion $W(t)=\sum_{\, n\, }|c_{e,n}(t)|^{2}-|c_{g,n}(t)|^{2}$

```wl
Plot[Total @(Abs[cet[t]]^2 - Abs[cgt[t]]^2), {t, 0, 40}, 
 PlotRange -> All, Frame -> True, FrameLabel -> {"gt", "W(t)"}, 
 PlotTheme -> "Detailed", PlotLegends -> None, ImageSize -> 400, 
 LabelStyle -> 10, 
 PlotLabel -> "Population inversion W(t) for the atom"]
```

#### Atom "thermalization" (atom interacting with a thermal state field)

Atom initially in $\frac{1}{\sqrt{2}}(|g\rangle +|e\rangle )$ and field in a thermal state with $\overline{n}=5$:

```wl
\[Rho]0 = 
  QuantumTensorProduct[QuantumState["+", atomBasis], ThermalState[5.]];
```

Evolve the density matrix $\rho (t)\, =\hat{U}(t)\, \rho (0)\, \hat{U}^{\dagger }(t)\, $. We use shortcut notation applied to matrix states:

```wl
\[Rho]t  = 
  QuantumState[SuperDagger[U][\[Rho]0], "Parameters" -> {g, t}];
```

Reduced density matrix of the atom:

```wl
atomMatrix = 
  Normal[QuantumPartialTrace[\[Rho]t[1, t], {2}]["DensityMatrix"]] // 
    ComplexExpand // Simplify;
```

Plot two components of the matrix:

```wl
Plot[Evaluate[atomMatrix[[1, {1, 2}]]], {t, 0, 200}, ImageSize -> 350, 
 PlotRange -> {-0.3, 1}, Frame -> True, LabelStyle -> 10, 
 FrameLabel -> {"t", ""}, 
 PlotLegends -> {"\!\(\*SubscriptBox[\(\[Rho]\), \(EE\)]\)", 
   "\!\(\*SubscriptBox[\(\[Rho]\), \(EG\)]\)"}, AspectRatio -> 0.5, 
 PlotLabel -> "Evolution of matrix elements of the atomic part"]
```

Atom gets "thermalized"; the time-average of the coherences → 0:

```wl
timeAverages = Integrate[(atomMatrix[[1, 2]]), {t, 0, tf}]/tf;
Plot[timeAverages, {tf, 5, 200}, PlotRange -> All, 
 AxesLabel -> {"g t", "Time Average"}, LabelStyle -> 10]
```

### Jaynes-Cummings with dissipation

$$\hat{H} = \hbar\nu\, \hat{a}^{\dagger}\hat{a} + \tfrac{1}{2}\hbar\omega\, \sigma_{z} + \hbar g\,(\sigma_{+}\hat{a} + \sigma_{-}\hat{a}^{\dagger}), \qquad \frac{d\rho}{dt} = -\frac{i}{\hbar}[\hat{H}, \rho] + \mathcal{D}_{\mathrm{decay}}(\rho) + \mathcal{D}_{\mathrm{loss}}(\rho)$$

where $\mathcal{D}_{decay},\mathcal{D}_{loss}$ are the master equation dissipator terms related to Lindblad jump operators $L_{1}=\sqrt{\gamma_{1}}(\sigma_{-}\otimes \, \mathcal{I})$ and $L_{2}=\sqrt{\gamma_{2}}(\mathcal{I}\, \otimes \, \hat{a})$, consider the initial state |e〉|n〉 with `n=2`.

Definition of the basis and the necessary operators:

```wl
atomBasis = QuantumBasis[{"e", "g"}];
\[Sigma]plus =  QuantumOperator["+", atomBasis];
\[Sigma]minus = QuantumOperator["-", atomBasis];
a := AnnihilationOperator[{2}];
```

Parameter setting:

```wl
SetFockSpaceSize[8];
\[Nu] = 10.;   (* field frequency *)
\[Omega] = 10.;   (* atomic frequency *)
gc = 15.; (* interaction coupling strength *)
tf = 1;        (* Final time *)
\[Psi]0 = 
  QuantumTensorProduct[QuantumState["0", atomBasis], FockState[2]];
```

Hamiltonian definition:

```wl
H =  \[Nu] SuperDagger[a]@a + 
   1/2 \[Omega] QuantumOperator["Z", atomBasis] + 
   gc (\[Sigma]plus@a + \[Sigma]minus@SuperDagger[a]);
```

Jump operators:

```wl
decayAtom = \[Sigma]minus@QuantumOperator["I"[$FockSize], {2}];
lossField = QuantumOperator["I", atomBasis]@a;
```

Helper definition to get the probability of finding the atom in the excited state in a time window:

```wl
atomProb[evolved_, tspecs_] := 
 Table[{t, (SuperDagger[QuantumState["0", atomBasis]]@evolved[t])[
    "Norm"]}, Evaluate[{t, Sequence @@ tspecs}]]
```

Compute the evolution for different damping parameters $\gamma =\gamma_{1}=\gamma_{2}$, and record the probability of finding the atom in the excited state:

```wl
probsE = Table[\[ScriptCapitalH] = 
    QuantumOperator["Hamiltonian"[H, {decayAtom, lossField}, {x, x}]];
      evolved = QuantumEvolve[\[ScriptCapitalH], \[Psi]0, {t, 0, 1}];
      atomProb[evolved, {0, 1, 0.01}],
      {x, {0.8, 0.4, 0.2}}];
```

Damped Rabi oscillations:

```wl
ListLinePlot[probsE, PlotRange -> All, Mesh -> All, Frame -> True, 
 PlotStyle -> Dashed, 
 FrameLabel -> {"t", 
   "\!\(\*SubscriptBox[\(P\), \
\(\(\(|\)\(e\)\)\(\[RightAngleBracket]\)\)]\)(t)"}, LabelStyle -> 9, 
 ImageSize -> 350, 
 PlotLegends -> {"\[Gamma],\[Kappa] = 0.8", "\[Gamma],\[Kappa] = 0.4", 
   "\[Gamma],\[Kappa] = 0.2"}, 
 PlotLabel -> 
  "Damped Rabi oscillations for different dissipation parameters"]
```

Now let's evolve different states of the form |e〉 ⊗ |α〉, for decay parameters $\gamma ,\kappa =\, 0.5$

Define the states:

```wl
alphas = {0.3, 0.6, 1.2, 1.8};
states0 = 
  Table[QuantumState["0", atomBasis] // 
    CoherentState[][\[Alpha]], {\[Alpha], alphas}];
```

Define the master equation super-operator:

```wl
\[ScriptCapitalH] = 
  QuantumOperator["Hamiltonian"[H, {decayAtom, lossField}, {0.5, 0.5}]];
```

Calculate the probabilities of the evolved states:

```wl
probsE = 
  atomProb[
     QuantumEvolve[\[ScriptCapitalH], #, {t, 0, 2}], {0, 2, 0.01}] & /@
    states0;
```

Plot the results:

```wl
ListLinePlot[probsE, InterpolationOrder -> 2, ImageSize -> 350, 
 AspectRatio -> 1/2, 
 PlotLegends -> (StringForm["\[Alpha] = ``", #] & /@ alphas), 
 PlotRange -> All, 
 AxesLabel -> {"t", 
   "\!\(\*SubscriptBox[\(P\), \(\( | \
\(e\)\)\(\[RightAngleBracket]\)\)]\)"}, LabelStyle -> 13, 
 PlotStyle -> Thickness[0.005]]
```

### Driven-Dissipative Kerr Nonlinear Oscillator

The Hamiltonian is:

$\hat{H}=\, -\Delta \hat{a}^{\dagger }\hat{a}+\frac{K}{2}(\hat{a}^{\dagger 2}\hat{a}^{2})+\varepsilon (\hat{a}^{\dagger }+\hat{a})$

where Δ is the detuning, K, is the Kerr nonlinearity and ε is the drive strength. Dissipation can be added with the photon loss channel with Lindblad operator $L_{loss}=\sqrt{\kappa }\hat{a}$

```wl
a := AnnihilationOperator[];
SetFockSpaceSize[25];
```

Define the Hamiltonian:

```wl
H = QuantumOperator[-\[CapitalDelta] SuperDagger[a]@a + 
    K/2 ((SuperDagger[a]^2)@(a^2)) + \[CurlyEpsilon] (SuperDagger[a] +
        a), "Parameters" -> {\[CapitalDelta], K, \[CurlyEpsilon]}];
```

We take the vacuum as the initial state:

```wl
\[Rho]0 = FockState[0];
```

Set up the Lindbladian with $\Delta =3,K=1,\varepsilon =1.5,\kappa =0.5\, $:

```wl
\[ScriptCapitalL] = 
  QuantumOperator["Liouvillian"[H[3, 1, 1.5], {a}, {0.5}]];
```

Calculate the evolution numerically , with final time $\, t_{f}=15$:

```wl
evolved = QuantumEvolve[I \[ScriptCapitalL], \[Rho]0, {t, 0, 15}];
```

Wigner function snapshots:

```wl
Partition[
  Table[wig = WignerRepresentation[evolved[t], {-5, 5}, {-5, 5}]; 
   DensityPlot[wig[x, p], {x, -5, 5}, {p, -5, 5}, PlotPoints -> 60, 
    ColorFunction -> "DeepSeaColors", FrameLabel -> {"x", "p"}, 
    PlotRange -> All, 
    PlotLabel -> StringForm["t = ``", t]], {t, {0, 1.25, 2.5, 5, 7.5, 
     10}}], {2}] // GraphicsGrid
```

We see the steady state resembles a coherent state up to distortion that depends of the parameters, it is known that evolution in Kerr nonlinear oscillators can produce coherent state manifold convergence

Set up a numerical function for comparing the fidelity of the evolved state at $t_{f}$ and a coherent state:

```wl
f[x_?NumericQ] := QuantumSimilarity[evolved[15], CoherentState[][x]]
```

Look for a coherent state that maximizes the fidelity with $\alpha \, =x+i\, y$, we limit the region from $1/4\le |\alpha |\le 1$:

```wl
sol = NMaximize[
  f[x + I  y], {x, y} \[Element] Annulus[{0, 0}, {1/4, 1}], 
  Method -> "SimulatedAnnealing"]
```

Following a semi-classical approach such the parameters are in the stability regime the coherent state amplitude of convergence can be derived as $\alpha \, =\, \frac{\varepsilon }{(\Delta \, -\, K\, |\alpha |^{2})+i\, \frac{\kappa }{2}}$

Taking the root with less energy:

```wl
Clear[\[Alpha]];
With[{\[CapitalDelta] = 3, \[CurlyEpsilon] = 3/2, \[Kappa] = 1/2, 
  K = 1}, First@
  SortBy[NSolve[{\[Alpha] == \[CurlyEpsilon]/((\[CapitalDelta] - 
         K Abs[\[Alpha]]^2) + I (\[Kappa]/2) )}, \[Alpha]], 
   Abs[#[[1, 2]]] &]]
```
