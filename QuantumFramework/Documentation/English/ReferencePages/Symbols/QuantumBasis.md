---
Template: Symbol
Name: QuantumBasis
Context: Wolfram`QuantumFramework`
Paclet: Wolfram/QuantumFramework
URI: Wolfram/QuantumFramework/ref/QuantumBasis
Keywords: [quantum basis, basis transformation, change of basis, output basis, input basis, picture, Schrodinger picture, Heisenberg picture, interaction picture, phase space]
SeeAlso: [QuantumState, QuantumTensorProduct, QuantumOperator, QuantumPhaseSpaceTransform]
RelatedGuides: [WolframQuantumComputationFramework]
RelatedTutorials: ["[Wolfram Quantum Framework Tutorial](Tutorial)"]
---

<!--
Faithful md source for Documentation/English/ReferencePages/Symbols/QuantumBasis.nb,
recovered with NotebookToMarkdown and re-verified against the repo kernel at ebc91a15
(2026-07-24). Deliberate deviations from the .nb, each verified against the live kernel:
  1. New picture coverage: the Details bullets and picture-name table, the options
     table, the Basic Examples unit on the picture, and the whole "Options" example
     section. The .nb documented no picture at all, though every QuantumBasis carries
     one; the forms shown here were fixed across 3a830e06, 669173d2, cb40c1e4 and ebc91a15, and
     are pinned by the "QuantumBasis - pictures" section of Tests/QuantumBasis.wlt.
     First-picture-wins precedence is documented from that section's measurement, not
     from the argument order it might suggest.
  2. Two Usage lines are repaired, not rewritten. NotebookToMarkdown recovered the
     first one's description as a literal StyleBox box-string, and the second lost its
     association braces and arrows to the TraditionalForm FormBox math route (literal
     {} and -> are invisible TeX grouping there). Both now say what the .nb rendered.
  3. Two example captions ended with a stray double quote instead of a colon
     ("Transform the state into the Wootters basis"") and are corrected; a caption must
     end in ":" or CheckDefinitionNotebook flags ExampleTextLastCharacter.
  4. Seven cells referred back with % and %%, which only resolve in an interactive
     session; on a rebuild they came back as Out[0][...]. Each now binds its
     intermediate to a name, which also makes those units self-contained.
  5. RandomSeed[0] was replaced by SeedRandom[0] in three places. RandomSeed is a bare
     System symbol with no down-values: it returns unevaluated and never seeded, so
     those examples were not reproducible. The affected output hints were re-taken
     from the reseeded run.
  6. Three pre-existing example defects, each refuted by the kernel and by the cell
     following it: QuantumBasis[{3,2}] was captioned as having 2-dimensional input
     and 3-dimensional output, but it is two OUTPUT qudits and no input
     ({"OutputQudits","InputQudits"} is {2,0}); basis["Split",2] was a no-op on a
     basis that already had 2 output qudits, and is now Split,1, which does move one
     to the input; and the purely-symbolic basis reused the symbol w, still bound
     from the Wootters cell, so a whole QuantumState was appearing inside it, and its
     free symbols are now Greek.
Note "HasseSIC" in the informationally-complete table is carried over unchanged and is
a pre-existing defect: the kernel's $QuditPhaseSpaceBasisNames spells it "HesseSIC".
So are the two N::meprec messages under Scope, which come from
PositiveSemidefiniteMatrixQ on the exact Tetrahedron POVM elements and would need
Quiet to hide. Both are tracked separately and deliberately not touched here.
-->

## Usage

<code>[QuantumBasis]()["*name*"]</code> represents a named quantum basis *name*.

<code>[QuantumBasis]()[<|$name_{1}$ -> $b_{1}$, $name_{2}$ -> $b_{2}$, …|>]</code> represents a quantum basis with basis elements $b_{i}$, having names $\mathit{name}_{i}$.

<code>[QuantumBasis]()[{$n_{1}$,$n_{2}$,…}]</code> represents a $n_{1}\times n_{2}\times \ldots $ dimensional computational basis of a composite system (many qudits).

<code>[QuantumBasis]()[*n*,*m*]</code> represents a $n^{m}$ dimensional computational basis of a composite system (*m* qudits, each one, *n*-dimensional).

<code>[QuantumBasis]()[{{$n_{1}$,$n_{2}$,…},{$m_{1}$,$m_{2}$,…}}]</code> represents a $n_{1}\times n_{2}\times \ldots $ dimensional computational basis output qudits, and $m_{1}\times m_{2}\times \ldots $ dimensional of the input qudits. Instead of dimension, one can add named basis, too.

## Details & Options

- A quantum basis may consist of distinct input and output qudits. For instance, a [QuantumOperator]() operates with both input and output qudits, while a typical [QuantumState]() only pertains to output. Therefore, it's important to recognize that the concept of a [QuantumBasis]() differs fundamentally from that of a basis in a Hilbert space. While users typically won't need to interact with this distinction directly, it is automatically handled behind the scenes for most practical purposes.
- The keys names can be any expression.
- Built-in names in the Hilbert space include:

|   |   |
|---|---|
| "Computational" or "I" | 2-dimensional computational basis |
| "Computational"[*d*] or "I"[*d*] | *d*-dimensional computational basis |
| `"Bell"` | Bell basis of 4-dimensional qudit |
| `"Pauli"` | basis of the Pauli matrices |
| "PauliX" or "X" | Pauli-X (Hadamard) basis |
| "PauliY" or "Y" | Pauli-Y basis |
| "PauliZ" or "Z" | Pauli-Z (computational) basis |
| "X"[*d*] or "Y"[*d*] or "Z"[*d*] | Generalized Pauli basis in *d* dimension (for qudits) |
| "JX" or "JY" or "JZ" | Generalized Spin basis with spin $\frac{1}{2}$ (or 2 dimension) |
| "JX"[*j*] or "JY"[*j*] or "JZ"[*j*] | Generalized Spin basis with spin *j* (or $2j+1$ dimension) |
| `"Fourier"` | 2-dimensional basis of the quantum Fourier transform |
| `"Fourier"[d]` | *d*-dimensional basis of the quantum Fourier transform |
| `"Schwinger"` | 2-dimensional Schwinger basis |
| `"Schwinger"[d]` | *d*-dimensional Schwinger basis |
| `"Dirac"` | basis of the Dirac (gamma) matrices |

- Built-in names in the phase space that are informationally complete, include:

|   |   |
|---|---|
| `"WignerMIC"` | Phase space basis corresponding to the WignerMICPOVM measurement (Minimal Information-ally Complete) |
| `"WignerMIC"[d]` | WignerMIC basis in the *d*-dimension |
| `"GellMannMIC"` | Phase space basis corresponding to the GellMannMICPOVM (generalized Pauli) measurements |
| `"GellMannMIC"[d]` | GellMann basis in *d*-dimension |
| `"Tetrahedron"` | Phase space basis corresponding to the TetrahedronSICPOVM measurements (Symmetric Information-ally Complete), which eigenstates form vertices of a [Tetrahedron]() |
| `"Tetrahedron"[a,b,c]` | A tetrahedron rotated by <code>"U"[*a*,*b*,*c*]</code> gate |
| `"HasseSIC"` | 3-dimensional [Hasse](https://arxiv.org/abs/1609.03075) basis corresponding to the HasseSICPOVM |
| `"HoggarSIC"` | 8-dimensional [Hoggar](https://arxiv.org/abs/1609.03075) basis corresponding to the HoggarSICPOVM |
| `"QBismSIC"[d]` | numeric SIC from [QBism](https://github.com/heyredhat/qbism/tree/master/qbism/sic_povms) research program up-to dimension $d=151$ |
| `"RandomHaarMIC"` | A random MIC basis based on Haar method |
| `"RandomBlochMIC"` | A MIC basis based on random Bloch sphere affine transformation |

- Built-in names in the phase space that are not informationally complete, include:

|   |   |
|---|---|
| `"GellMann"` | basis of the Gell-Mann matrices |
| `"GellMann"[d]` | basis of generalized Gell-Mann matrices for dimension *d* |
| `"Wigner"` | 2-dimensional basis of the Wigner phase space operators with respect to the computational basis |
| `"Wigner"[d]` | *d*-dimensional basis of the Wigner phase space operators with respect to the computational basis |
| `"Wootters"` | 2-dimensional phase space basis |
| `"Wootters"[p]` | *p*-dimensional phase space basis where *p* is a prime number |
| `"Bloch"` | Phase space basis of generalized Bloch coordinates |

- A few properties of [QuantumBasis]() are:

|   |   |
|---|---|
| `"Association"` | An association with keys the formatted symbols of basis elements, and values the corresponding tensor |
| `"Elements"` | Corresponding elements of a basis |
| `"OrthogonalElements"` | Corresponding orthogonal elements of a basis |

- Besides its input and output qudits, a [QuantumBasis]() carries a **picture**, recovered with `qb["Picture"]`. The picture records where the time dependence of a system written in this basis is kept; it is part of the basis identity but says nothing about the basis elements themselves. Built-in pictures are:

|   |   |
|---|---|
| `"Schrodinger"` | states carry the time dependence and operators are fixed; the default |
| `"Heisenberg"` | operators carry the time dependence and states are fixed |
| `"Interaction"` | the time dependence is divided between states and operators |
| `"PhaseSpace"` | the tag a phase-space basis carries; a basis whose *elements* form a phase-space frame sets it for itself, and a state expanded in such a basis is a quasi-probability distribution rather than an amplitude vector |

- A picture may be named positionally, on either side of the basis specification, or given as the `"Picture"` option, and the three forms agree. So <code>[QuantumBasis]()["Heisenberg"]</code> is the default basis in the Heisenberg picture, and <code>[QuantumBasis]()["Schrodinger"]</code> is <code>[QuantumBasis]()[]</code>. A picture may also trail a two-sided specification, as in <code>[QuantumBasis]()["X", "Y", "Heisenberg"]</code>.
- When more than one picture is written, the one written **first** wins, positional and option alike: the arguments are folded in reverse, so the earliest is applied last. <code>[QuantumBasis]()["Heisenberg", "Picture" -> "Interaction"]</code> is in the Heisenberg picture.
- A trailing integer is read as a multiplicity rather than a dimension, so <code>[QuantumBasis]()["Heisenberg", 3]</code> is three qubits in the Heisenberg picture.
- A phase-space basis name sets `"PhaseSpace"` on its own, so a picture normally only needs naming for `"Heisenberg"` and `"Interaction"`.
- Such a basis will not be renamed to another picture positionally: its elements, not a label, are what put it in phase space, so <code>[QuantumBasis]()["Wigner", "Heisenberg"]</code> is rejected with `QuantumBasis::phaseSpacePicture`. The `"Picture"` option is deliberately left open for code that has actually undone the transform.
- The picture travels with the basis into whatever is built on it, so <code>[QuantumState]()[*spec*, "Heisenberg"]</code> and <code>[QuantumOperator]()[*spec*, "Heisenberg"]</code> both land in the Heisenberg picture.
- Any picture outside the four above is rejected with `QuantumBasis::picture`. A picture followed by a tail that is not a basis specification at all, such as a bare symbol or a non-integer number, is rejected with `QuantumBasis::invalidSpec` naming that tail; a tail that fails on its own terms reports its own reason instead, so an unrecognized basis *name* gives `QuditBasis::invalidName` and a negative integer gives `QuditBasis::invalidArgs`.

- Options for [QuantumBasis]() are:

| option | default | effect |
|---|---|---|
| <code>"Picture"</code> | <code>"Schrodinger"</code> | the picture the basis is written in, one of the four names above |
| <code>"Label"</code> | <code>None</code> | the label the basis displays under; a named basis, or one given a single dimension, fills this in for itself, so <code>[QuantumBasis]()["PauliX"]</code> is labelled `"PauliX"` and <code>[QuantumBasis]()[3]</code> is labelled `"I"[3]`, while a list of dimensions or an explicit [QuditBasis]() leaves it `None` |
| <code>"ParameterSpec"</code> | <code>{}</code> | symbolic parameters the basis depends on, given as `{p, min, max}` triples; `"Parameter"` and `"Parameters"` are accepted as spellings of the same option, and a bare symbol is completed to a unit range |

## Basic Examples

Create a 2-dimensional basis:

```wl
QuantumBasis[2]
```

Note with no input, the basis is automatically set to 2D by default:

```wl
QuantumBasis[2] == QuantumBasis[]
```
<!-- => True -->

---

Create a 3-dimensional basis:

```wl
QuantumBasis[3]
```

---

Create a 2×2×2 dimensional basis (three qubits):

```wl
QuantumBasis[2, 3]
```

---

Create a composite basis of 2- and 3-dimensional qudits (a qubit with a qutrit):

```wl
QuantumBasis[{2, 3}]
```

---

Create a 2-dimensional basis using an explicit element representation with arbitrary names:

```wl
QuantumBasis[<|"\[DownArrow]" -> {1, 2 I}, "\[UpArrow]" -> {0, 1}|>]
```

---

Construct a Pauli-Y basis:

```wl
basis = QuantumBasis["Y"]
```

---

Construct a Bell basis for a single 4-dimensional qudit (quqrit):

```wl
basis = QuantumBasis["Bell"]
```

Return its association (basis names with corresponding representations):

<!-- #| screenshot: true -->
```wl
basis["Association"]
```

---

A quantum basis with Pauli-X as the output and computational 2-dimensional basis as the input:

```wl
QuantumBasis[{"X"}, {2}]
```

Show its elements:

```wl
QuantumBasis[{"X"}, {2}]["Association"]
```

---

The elements of a 3-dimensional [QuantumBasis]() can be arbitrary (potentially non-orthonormal) vectors:

```wl
skewed = QuantumBasis[<|"a" -> {1, -3, 0}, "b" -> {3, 1, 2}, "c" -> {2, -1, 4}|>]
```

Return the list of its orthogonalized basis elements:

```wl
skewed["OrthogonalElements"]
```
<!-- => {{1/Sqrt[10], -(3/Sqrt[10]), 0}, {3/Sqrt[14], 1/Sqrt[14], Sqrt[2/  7]}, {-(3/Sqrt[35]), -(1/Sqrt[35]), Sqrt[5/7]}} -->

---

Represent the 2-dimensional Schwinger basis of rank-2 (matrix) elements:

```wl
schwinger = QuantumBasis["Schwinger"]
```

<!-- #| screenshot: true -->
```wl
schwinger["Association"]
```

---

Every basis is also written in a picture, which defaults to the Schrodinger one:

```wl
QuantumBasis[]["Picture"]
```
<!-- => "Schrodinger" -->

Name a picture to get the default basis written in it:

```wl
QuantumBasis["Heisenberg"]["Picture"]
```
<!-- => "Heisenberg" -->

The picture may sit on either side of the basis specification, or be given as an option, and the three agree:

```wl
QuantumBasis["Heisenberg", "PauliX"] === QuantumBasis["PauliX", "Heisenberg"] === 
 QuantumBasis["PauliX", "Picture" -> "Heisenberg"]
```
<!-- => True -->

## Scope

### Tetrahedron basis:

Tetrahedron basis corresponds to the coordinates of a tetrahedron (after some normalization);

```wl
pnts1 = PolyhedronCoordinates@Tetrahedron[2 Sqrt[2/3]] // Simplify;
```

Show the tetrahedron and Bloch sphere:

```wl
Show[QuantumState["UniformMixture"]["BlochPlot"], 
 Graphics3D[{Opacity[.3], Red, Tetrahedron[pnts1]}]]
```

Given the elements of tetrahedron basis, find the corresponding quantum states and their Bloch vectors:

```wl
pnts2 = Normalize@QuantumState[#]["BlochCartesianCoordinates"] & /@ 
   QuditBasis["Tetrahedron"]["Elements"];
```

Show the Bloch vectors:

```wl
Show[QuantumState["UniformMixture"]["BlochPlot"], 
 Graphics3D[{Opacity[.3], Red, Tetrahedron[pnts2]}]]
```

From TetrahedronSICPOVM which is a SIC-POVM, find the corresponding quantum states and their Bloch vectors and classify into two distinct categories based on their eigenvalues: those with an eigenvalue of 0 and those with an eigenvalue of 1/2:

```wl
pnts3 = Map[QuantumState[#]["BlochCartesianCoordinates"] &] /@ 
    Eigenvectors /@ 
     QuantumMeasurementOperator["TetrahedronSICPOVM"][
      "POVMElements"] // Transpose;
```

To visualize these states, plot them on the Bloch sphere, with a red tetrahedron representing states with eigenvalue 1/2, and a blue tetrahedron representing states with eigenvalue 0.

```wl
Show[QuantumState["UniformMixture"]["BlochPlot"], 
 Graphics3D[{Opacity[.5], Red, Tetrahedron[pnts3[[1]]], Blue, 
   Tetrahedron[pnts3[[2]]]}]]
```

Verify that the coordinates of the previously obtained points are identical.

```wl
pnts1 == pnts2 == pnts3[[1]]
```
<!-- => True -->

### Gell-Mann basis:

Gell-Mann matrices are traceless:

```wl
Tr /@ Rest@QuantumBasis["GellMann"[5]]["Elements"]
```
<!-- => {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0} -->

Gell-Mann matrices are Hermitian:

```wl
HermitianMatrixQ /@ Rest@QuantumBasis["GellMann"[5]]["Elements"]
```
<!-- => {True, True, True, True, True, True, True, True, True, True, True, \ True, True, True, True, True, True, True, True, True, True, True, \ True, True} -->

There are $d-1$ diagonal elements:

```wl
DiagonalMatrixQ /@ 
  Rest@QuantumBasis["GellMann"[7]]["Elements"] // Counts
```
<!-- => <|False -> 42, True -> 6|> -->

The rest are symmetric or anti-symmetric matrices:

```wl
SymmetricMatrixQ /@ 
  Rest@QuantumBasis["GellMann"[3]]["Elements"] // Counts
```
<!-- => <|True -> 5, False -> 3|> -->

```wl
AntisymmetricMatrixQ /@ 
  Rest@QuantumBasis["GellMann"[3]]["Elements"] // Counts
```
<!-- => <|False -> 5, True -> 3|> -->

For qudits, particularly qutrits, the POVM elements can be constructed using the Gell-Mann matrices as a basis.

Generate a random mixed state in 3D Hilbert space:

```wl
SeedRandom[0];
r = QuantumState["RandomMixed", {d = 3}]
```

Find the amplitudes (i.e. the vectorized density matrix) in the phase space of GellMann basis:

```wl
amplitudes = QuantumPhaseSpaceTransform[r, "GellMann"[d]]["AmplitudesList"]
```

Given the GellMann basis elements $E_{i}$, find <code>[Tr]()[ρ.$E_{i}$]</code>:

```wl
traces = Tr[r["DensityMatrix"] . #] & /@ 
   QuantumBasis["GellMann"[3]]["Elements"] // Chop
```

Check the density matrix for a qutrit is expressed as $\rho =\frac{1}{3}I+\frac{1}{2}\sum_{i}r_{i}E_{i}$

```wl
traces/amplitudes
```
<!-- => {3., 2., 2., 2., 2., 2., 2., 2., 2.} up to floating-point noise -->

### Wootters basis:

Create a random pure state:

```wl
SeedRandom[0];
\[Rho] = QuantumState["RandomPure"]
```

Represent the Wootters basis elements:

<!-- #| screenshot: true -->
```wl
QuantumBasis["Wootters"]["Association"]
```

Transform the state into the Wootters basis:

```wl
w = QuantumPhaseSpaceTransform[\[Rho], "Wootters"]
```

Show the corresponding amplitudes in the Wootters basis:

<!-- #| screenshot: true -->
```wl
w["Amplitudes"]
```

Reshape the amplitudes into 2×2 matrix and check if it is the same as "PrimeQuasiProbability" of the state:

```wl
Partition[w["AmplitudesList"], {2}] == \[Rho]["PrimeQuasiProbability"]
```
<!-- => True -->

Verify that the sum of each column matches the probabilities in the computational basis (i.e., position basis):

```wl
Total@\[Rho]["PrimeQuasiProbability"] == \[Rho]["ProbabilitiesList"]
```
<!-- => True -->

Verify that the sum of each row matches the probabilities in the Pauli-X basis (i.e., momentum basis):

```wl
Total /@ \[Rho]["PrimeQuasiProbability"] == 
 QuantumState[\[Rho], "X"]["ProbabilitiesList"]
```
<!-- => True -->

### From SIC-POVM to corresponding phase space basis:

Every POVM corresponds to a basis of the doubled state representation. Consider Tetrahedron:

```wl
\[ScriptCapitalM] = QuantumMeasurementOperator["TetrahedronSICPOVM"]
```

Check POVM elements are explicitly positive semidefinite, and sum up to identity:

```wl
And @@ (PositiveSemidefiniteMatrixQ /@ \[ScriptCapitalM][
    "POVMElements"])
Total[\[ScriptCapitalM]["POVMElements"]] // Normal // N // Chop
```
<!-- => True -->

<!-- => {{1.`, 0}, {0, 1.`}} -->

Gram matrix is used to transform between basis and its corresponding POVM:

```wl
\[ScriptCapitalG] = 
  Inverse[Outer[Tr@*Dot, #, #, 1]] &@\[ScriptCapitalM][
     "POVMElements"] // Simplify;
```

Check if Gram matrix is explicitly positive semidefinite:

```wl
PositiveSemidefiniteMatrixQ[\[ScriptCapitalG]]
```
<!-- => True -->

Create a basis based on above POVM elements:

```wl
basis = QuantumBasis[\[ScriptCapitalG] . \[ScriptCapitalM][
     "POVMElements"] // Simplify]
```

Create a random mixed state:

```wl
SeedRandom[0];
\[Rho] = QuantumState["RandomMixed"];
```

Apply the measurement operator to obtain the corresponding probabilities.

<!-- #| screenshot: true -->
```wl
mea = \[ScriptCapitalM][\[Rho]]["Probabilities"]
```

Find the phase space representation of the state in the given basis and return its amplitudes:

```wl
ps = QuantumPhaseSpaceTransform[\[Rho], basis]["Amplitudes"]
```
<!-- => values {0.392528, 0.300255, 0.114132, 0.193084} keyed by QuditName[0..3] -->

Double the state (analogous to vectorizing a density matrix) and transform it into the new basis:

```wl
db = QuantumState[\[Rho]["Double"], basis]["Amplitudes"]
```
<!-- => values {0.392528, 0.300255, 0.114132, 0.193084} keyed by QuditName[0..3] -->

Calculate probabilities using the geometric relationship between the Tetrahedron corners and the Bloch vector:

```wl
mn = (\[Rho]["BlochVector"] . Normalize[#] + 1)/4 & /@ 
  PolyhedronCoordinates@Tetrahedron[]
```
<!-- => {0.392528, 0.300255, 0.114132, 0.193084} -->

Check the amplitudes are the same:

```wl
mn == Values[db] == Values[ps] == Values[mea]
```
<!-- => True -->

[QuantumWeylTransform]() can be used to undouble the state:

```wl
QuantumWeylTransform[QuantumState[\[Rho]["Double"], basis]] == \[Rho]
```
<!-- => True -->

## Generalizations and Extensions

Create a composite basis of a qutrit and a qubit, both of them output qudits:

```wl
basis = QuantumBasis[{3, 2}]
```

<!-- #| screenshot: true -->
```wl
Normal /@ basis["ElementAssociation"]
```

```wl
{basis["OutputQudits"], basis["InputQudits"]}
```
<!-- => {2, 0} -->

Set the number of output qudits, which sends the remaining ones to the input:

```wl
basis2 = basis["Split", 1]
```

```wl
{basis2["OutputQudits"], basis2["InputQudits"]}
```
<!-- => {1, 1} -->

Input or output can be empty (consisting of only 1-dimensional qudits), as the input of the unsplit basis is:

```wl
basis["Input"]
```

## Options

### "Picture"

All four pictures are accepted as a bare argument:

```wl
AssociationMap[QuantumBasis[#]["Picture"] &, 
 {"Schrodinger", "Heisenberg", "Interaction", "PhaseSpace"}]
```
<!-- => <|"Schrodinger" -> "Schrodinger", "Heisenberg" -> "Heisenberg", "Interaction" -> "Interaction", "PhaseSpace" -> "PhaseSpace"|> -->

---

`"Schrodinger"` is the default, so naming it changes nothing:

```wl
QuantumBasis["Schrodinger"] === QuantumBasis[]
```
<!-- => True -->

---

A picture chooses a convention, not a basis, so the qudits are left at the default:

```wl
QuantumBasis["Heisenberg"]["Dimensions"]
```
<!-- => {2} -->

A trailing integer is still a multiplicity, so this is three qubits rather than one qutrit:

```wl
QuantumBasis["Heisenberg", 3]["Dimensions"]
```
<!-- => {2, 2, 2} -->

---

A picture may also trail a two-sided specification, leaving both qudits in place:

```wl
QuantumBasis["X", "Y", "Heisenberg"]["Picture"]
```
<!-- => "Heisenberg" -->

```wl
QuantumBasis["X", "Y", "Heisenberg"]["Dimensions"]
```
<!-- => {2, 2} -->

---

When two pictures are written, the first one wins, whichever spelling each uses:

```wl
{QuantumBasis["Heisenberg", "Interaction"]["Picture"], 
 QuantumBasis["Heisenberg", "Picture" -> "Interaction"]["Picture"], 
 QuantumBasis["Picture" -> "Interaction", "Heisenberg"]["Picture"]}
```
<!-- => {"Heisenberg", "Heisenberg", "Interaction"} -->

---

A phase-space basis name sets its own picture, without being asked:

```wl
AssociationMap[QuantumBasis[#]["Picture"] &, 
 {"Computational", "PauliX", "Fourier", "Wigner", "Pauli", "Wootters"}]
```
<!-- => <|"Computational" -> "Schrodinger", "PauliX" -> "Schrodinger", "Fourier" -> "Schrodinger", "Wigner" -> "PhaseSpace", "Pauli" -> "PhaseSpace", "Wootters" -> "PhaseSpace"|> -->

---

The picture carries into whatever is built on the basis:

```wl
QuantumState["0", "Heisenberg"]["Picture"]
```
<!-- => "Heisenberg" -->

```wl
QuantumOperator["X", "Heisenberg"]["Picture"]
```
<!-- => "Heisenberg" -->

### "Label"

A bare basis has no label:

```wl
QuantumBasis[]["Label"]
```
<!-- => None -->

---

A named basis labels itself:

```wl
QuantumBasis["PauliX"]["Label"]
```
<!-- => "PauliX" -->

---

Give `"Label"` to override it, alongside a picture if wanted:

```wl
QuantumBasis["Heisenberg", "Label" -> "\[ScriptCapitalH]"]["Label"]
```
<!-- => "\[ScriptCapitalH]" -->

### "ParameterSpec"

A basis with no symbolic parameters has an empty specification:

```wl
QuantumBasis[]["ParameterSpec"]
```
<!-- => {} -->

---

A bare symbol is completed to a unit range, and `"Parameter"` is accepted for the same option:

```wl
QuantumBasis["Parameter" -> t]["ParameterSpec"]
```
<!-- => {{t, 0, 1}} -->

---

Give the range explicitly to control it; `"Parameters"` returns just the symbols:

```wl
QuantumBasis["ParameterSpec" -> {{t, 0, 2 Pi}}]["Parameters"]
```
<!-- => {t} -->

## Applications

Visualize a quantum state in terms of association between basis states and the corresponding complex amplitudes:

```wl
state = QuantumState[{1/Sqrt[2], -1/Sqrt[2]}, "Z"]
```

```wl
state["Formula"]
```
<!-- => the typeset formula (1/Sqrt[2])|z+> - (1/Sqrt[2])|z->; "Formula" returns RawBoxes, so there is no literal value to compare against -->

Change the basis into the Pauli-X basis:

```wl
newState = QuantumState[state, "X"]
```

```wl
newState["Formula"]
```
<!-- => the typeset formula |x->; "Formula" returns RawBoxes, so there is no literal value to compare against -->

---

[QuantumBasis]() objects can be constructed purely symbolically (without explicit vector or matrix elements):

```wl
basis = QuantumBasis[<|"\[UpArrow]" -> {\[Alpha], \[Beta]}, 
   "\[DownArrow]" -> {\[Gamma], \[Delta]}|>]
```

Represent a collection of 2 qudits:

```wl
pair = QuantumBasis[basis, 2]
```

```wl
pair["Names"]
```
<!-- => {Wolfram`QuantumFramework`QuditName[{"\[UpArrow]", "\[UpArrow]"}, "Dual" -> False], Wolfram`QuantumFramework`QuditName[{"\[UpArrow]", "\[DownArrow]"}, "Dual" -> False], \ Wolfram`QuantumFramework`QuditName[{"\[DownArrow]", "\[UpArrow]"}, "Dual" -> False], \ Wolfram`QuantumFramework`QuditName[{"\[DownArrow]", "\[DownArrow]"}, "Dual" -> False]} -->

---

Create a 3-qubit version of the same basis by taking a tensor product of the basis with itself 3 times:

```wl
threeX = QuantumBasis["X", 3]
```

View the element names of the 3-qubit system:

<!-- #| screenshot: true -->
```wl
threeX["Names"]
```

View the dual element names:

<!-- #| screenshot: true -->
```wl
threeX["Dual"]["Names"]
```

---

Winger basis of a qudit:

```wl
wigner = QuantumBasis["Wigner"[3]]
```

Note that the basis names (i.e. [QuditName]()) has the format as $\mathcal{W}_{p,q}$ with *p* and *q* runs from 0 to $d-1$. See [QuantumWignerTransform]() for more details.

<!-- #| screenshot: true -->
```wl
MatrixForm /@ wigner["Association"]
```

## Properties and Relations

[QuantumBasis]() is often used as an argument of [QuantumState](), [QuantumOperator](), and [QuantumMeasurementOperator]():

<!-- #| screenshot: true -->
```wl
QuantumState[{1, 0}, "Y"]
```

<!-- #| screenshot: true -->
```wl
QuantumOperator[{{1, 2}, {3, 4}}, "X"]
```

<!-- #| screenshot: true -->
```wl
QuantumMeasurementOperator["X", {1}]
```

---

[QuantumTensorProduct]() can be used to tensor product different bases:

```wl
basisTensor = 
 QuantumTensorProduct[{QuantumBasis["X"], QuantumBasis["Fourier"]}]
```

<!-- #| screenshot: true -->
```wl
basisTensor["Names"]
```

```wl
basisTensor == QuantumBasis[{"PauliX", "Fourier"}]
```
<!-- => True -->

The first basis can then be traced out, leaving only the computational basis:

```wl
basisTrace = QuantumPartialTrace[basisTensor, {1}]
```

```wl
basisTrace["Names"]
```
<!-- => {Wolfram`QuantumFramework`QuditName[ Subscript["F", 1], "Dual" -> False], Wolfram`QuantumFramework`QuditName[ Subscript["F", 2], "Dual" -> False]} -->

---

<code>[QuantumBasis]()[…,*n*]</code> is equivalent to the [QuantumTensorProduct]() of multiple copies:

```wl
QuantumBasis["X", 3] == 
 QuantumTensorProduct[QuantumBasis["X"], QuantumBasis["X"], 
  QuantumBasis["X"]]
```
<!-- => True -->

## Interactive Examples

Geometrical representations of qubit pure states in different bases:

<!-- #| screenshot: true -->
```wl
DynamicModule[{
  refpoints = 
   1/Sqrt[2] Prepend[#, Norm[#]] &@
      QuantumState[#]["BlochVector"] & /@ {"1", "0", "+", "R"},
  pointNames = {"1", "0", "+", "R"},
  spheres = {"Bloch", "Wigner", "WignerMIC", "Wootters", "Tetrahedron",
     "Pauli", "GellMann", "RandomMIC", "RandomMIC"[Method -> "Bloch"]},
   colors = <|"Bloch" -> LightGray, "Wigner" -> Red, 
    "WignerMIC" -> Blue, "Wootters" -> Orange, "Tetrahedron" -> Green,
     "Pauli" -> Black, "GellMann" -> Magenta, "RandomMIC" -> Cyan, 
    "RandomMIC"[Method -> "Bloch"] -> LightCyan|>, showLabels = False,
   vectors},
 Panel@Column[{
    CheckboxBar[Dynamic[spheres], spheres, 
     Appearance -> "Vertical" -> {3, 4}],
    spheres = 
     DeleteElements[
      spheres, {"RandomMIC", "RandomMIC"[Method -> "Bloch"]}];
    Row[{"Labels: ", Checkbox[Dynamic[showLabels]]}],
    Dynamic[
     Panel[Legended[Graphics3D[{
         Opacity[.3],
         
         If[MemberQ[spheres, "Bloch"], {LightGray, Sphere[], 
           If[showLabels, 
            MapThread[
             Text[Style[Ket[{#1}], Black], Rest[#2]] &, {pointNames, 
              refpoints}], Nothing]}, Nothing],
         
         Map[basisName |-> 
           Block[{basis = QuditBasis[basisName], center,
             points
             },
            
            center = 
             QuantumState[QuantumState["UniformMixture"]["Double"], 
               basis]["StateVector"];
            
            points = 
             Chop[QuantumState[QuantumState[#]["Double"], basis][
                   "StateVector"] - center] & /@ pointNames // Normal;
            {colors[basisName], PointSize[0.03], 
             If[showLabels, 
              MapThread[
               Text[Style[Ket[{#1}], Black], Rest[#2]] &, {pointNames,
                 Normal[points + Threaded[center]]}], Nothing],
             
             ConvexHullMesh[
              Chop@*Rest /@ (Threaded[center] + 
                 AffineTransform[
                   Transpose[points] . Inverse[Transpose[refpoints]]]@
                  RandomPoint[Sphere[{0, 0, 0, 0}], 6000]), 
              MeshCellStyle -> {1 -> None}]
             }],
          DeleteCases["Bloch"]@spheres
          ]}, Axes -> True, PlotRange -> {{-1, 1}, {-1, 1}, {-1, 1}}],
        SwatchLegend[Values@colors, Keys@colors]], 
      Background -> White
      ]]
    }]
 ]
```
