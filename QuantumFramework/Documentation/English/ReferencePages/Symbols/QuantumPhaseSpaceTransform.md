---
Template: Symbol
Name: QuantumPhaseSpaceTransform
Context: Wolfram`QuantumFramework`
Paclet: Wolfram/QuantumFramework
URI: Wolfram/QuantumFramework/ref/QuantumPhaseSpaceTransform
Keywords: [phase space, Wigner function, Husimi Q, Glauber-Sudarshan P, quasi-probability]
SeeAlso: [QuantumWignerTransform, QuantumBasis]
RelatedGuides: [WolframQuantumComputationFramework]
RelatedTutorials: ["[Wolfram Quantum Framework Tutorial](Tutorial)"]
---

<!--
Faithful md source for Documentation/English/ReferencePages/Symbols/QuantumPhaseSpaceTransform.nb,
recovered with NotebookToMarkdown and re-verified against the repo kernel at d0c3b886
(2026-07-24). Deliberate deviations from the .nb, each verified against the live kernel:
  1. The informationally-complete table named a "HasseSIC" basis, a "Hasse" basis link
     and a HasseSICPOVM measurement. None of the three exists: both
     QuantumBasis["HasseSIC"] and QuantumMeasurementOperator["HasseSICPOVM"] fail with
     QuditBasis::invalidName. The kernel spells the basis "HesseSIC"
     ($QuditPhaseSpaceBasisNames, NamedBases.m:22) and the measurement "HesseSICPOVM"
     ($QuantumMeasurementOperatorNames, NamedMeasurementOperators.m:7), after the Hesse
     configuration behind the exceptional d=3 SIC; the fiducial is {0,1,-1}/Sqrt[2]
     (MIC.m:168). The linked paper spells it Hesse throughout and Hasse nowhere. The row
     keeps its "3-dimensional" wording, which names the Hilbert space the SIC lives in,
     as the neighbouring 8-dimensional Hoggar row does; the basis itself has 9 elements.
     Every other name in the table was re-run and constructs.
  2. The cell reading "% == %%" compared the two results above it, which only resolve in
     an interactive session; on a rebuild it comes back as Out[0][...]. The two cells now
     bind "amplitudes" and "probabilities", and the check is "probabilities == amplitudes",
     preserving the original operand order (% was the probabilities, %% the amplitudes).
     That also makes the unit self-contained.
  3. Two example captions ended in a word rather than a colon ("...and return amplitudes",
     "Check they are the same") and now end in ":", the convention for a caption that
     introduces an input. The two cells that close a thought rather than introduce code
     keep their full stop. Note that CheckDefinitionNotebook cannot be used to police this
     here: it rejects a paclet reference page with ::nbortype, "not a valid resource
     definition notebook", and does so for the original page as much as for this rebuild,
     so the ExampleTextLastCharacter hint never fires either way.
  4. Two output hints under Scope were re-taken. SeedRandom[0]; QuantumState["RandomPure"]
     is reproducible, but it no longer returns the state the stored outputs came from: the
     magnitudes still agree, so ProbabilitiesList matched to its last two digits, while the
     phase differs, so the QBismSIC amplitudes are a different probability vector entirely.
     The example's own claim is unaffected and was re-run: trans["Matrix"] . ps - prob
     chops to {0, 0} for the state the kernel actually produces. The RandomMixed draws in
     Basic Examples still reproduce their stored outputs to all 17 digits and are untouched.
  5. The cell after "Inverse basis means that the direction of input is flipped:" is an
     Output cell with no Input: a small ProcessTheory Diagram of a box with its ports, whose
     generating input had been deleted from the notebook. NotebookToMarkdown drops it, and
     "PreserveOutputs" -> True does not rescue it either (it skips graphics), so the figure
     would have been lost on the rebuild. It is pinned here as a "#| boxes: true" Output
     cell carrying the notebook's own box expression, extracted programmatically and checked
     to round-trip through the parse the forward path applies to a boxes fence.
  6. The .nb's empty scaffolding sections (Generalizations & Extensions, Options with two
     XXXX subsections, Applications, Properties & Relations, Possible Issues, Interactive
     Examples, Neat Examples) carry no content and are not carried here. MarkdownToNotebook
     emits that scaffolding itself, which is why the QuantumBasis rebuild kept all of them
     while its .md names none; the same was confirmed on the notebook built from this file.
  7. Output hints are written on one line. The recovery hard-wrapped the longer ones with
     backslash continuations inside the comment, which is only an annotation artefact.
  8. The "GellMannMIC"[d] row described a "GellMann basis", naming a different basis that
     this same table documents separately. The two are not one object: "GellMann"[d] is
     the identity prepended to the Gell-Mann matrices (NamedBases.m:189), while
     "GellMannMIC"[d] is the Gram dual of the GellMann MIC POVM (NamedBases.m:194,
     MIC.m:113). Their traces separate them exactly, {d, 0, ...} against {1, 1, ...},
     and at d = 3 the element lists were run and are unequal. The row now names
     GellMannMIC, the basis it is about, as the parallel WignerMIC row above it does.
  9. The table documented "QBismSIC"[d] but not the bare "QBismSIC", which three of the
     page's own examples use. The bare form is the d = 2 default of that same definition
     (NamedBases.m:274, over QBismSICPOVM[d : _Integer : 2] at MIC.m:156), and its
     elements were confirmed identical to "QBismSIC"[2]: four 2 x 2 elements, Dimension
     4. A bare row is added ahead of the parameterized one, worded as the other
     bare/parameterized pairs here are, so the examples run on a documented form.
 10. Two typographic repairs, neither of which changes what the page claims. A caption
     under Basic Examples read "returns the same stat as the one in the Hilbert space"
     and now reads "state". Two rows read "Information-ally Complete", the residue of a
     hyphenated line break in the notebook, and now read "Informationally Complete".
Three further differences are the forward path's doing, not this file's, and were measured
by converting the old and the new .nb with nb-reader and diffing:
  a. The Keywords cell is not emitted. The keywords survive in the frontmatter above and
     the KeywordsSection header is still built, but the cell listing them is dropped. The
     QuantumBasis rebuild did the same, so this is MarkdownToNotebook's behaviour rather
     than anything lost in recovery.
  b. The Related Guides link is labelled with the URI slug, WolframQuantumComputationFramework,
     where the .nb had "Wolfram Quantum Computation Framework". The cell is still the plain
     Link ButtonBox that DocumentationBuild's harvester requires and the paclet: URI is
     unchanged, so only the visible label differs; guideLinkContent derives both the label
     and the URI from one bare frontmatter name, and RelatedGuides has no labelled form the
     way RelatedTutorials does. The QuantumBasis rebuild shows the same label.
  c. "trans . ps" in a caption was an InlineFormula in the .nb and is an InlineCode span
     here, because its box tree is neither a call form nor 2D math and the walker renders
     that as a code span.
-->

## Usage

<code>[QuantumPhaseSpaceTransform]()[*object*,*basis*]</code> represents the transformation of quantum object into the phase space basis

## Details & Options

- Built-in names for the phase space include:

|   |   |
|---|---|
| `"GellMann"` | basis of the Gell-Mann matrices |
| `"GellMann"[d]` | basis of generalized Gell-Mann matrices for dimension *d* |
| `"Wigner"` | 2-dimensional basis of the Wigner phase space operators with respect to the computational basis |
| `"Wigner"[d]` | *d*-dimensional basis of the Wigner phase space operators with respect to the computational basis |
| `"Wootters"` | 2-dimensional phase space basis |
| `"Wootters"[p]` | *p*-dimensional phase space basis where *p* is a prime number |
| `"WignerMIC"` | Phase space basis corresponding to the WignerMICPOVM measurement (Minimal Informationally Complete) |
| `"WignerMIC"[d]` | WignerMIC basis in the *d*-dimension |
| `"Bloch"` | Phase space basis of generalized Bloch coordinates |
| `"GellMannMIC"` | Phase space basis corresponding to the GellMannMICPOVM (generalized Pauli) measurements |
| `"GellMannMIC"[d]` | GellMannMIC basis in *d*-dimension |
| `"Tetrahedron"` | Phase space basis corresponding to the TetrahedronSICPOVM measurements (Symmetric Informationally Complete), which eigenstates form vertices of a [Tetrahedron]() |
| `"Tetrahedron"[a,b,c]` | A tetrahedron rotated by <code>"U"[*a*,*b*,*c*]</code> gate |
| `"HesseSIC"` | 3-dimensional [Hesse](https://arxiv.org/abs/1609.03075) basis corresponding to the HesseSICPOVM |
| `"HoggarSIC"` | 8-dimensional [Hoggar](https://arxiv.org/abs/1609.03075) basis corresponding to the HoggarSICPOVM |
| `"QBismSIC"` | 2-dimensional numeric SIC from [QBism](https://github.com/heyredhat/qbism/tree/master/qbism/sic_povms) research program |
| `"QBismSIC"[d]` | numeric SIC from [QBism](https://github.com/heyredhat/qbism/tree/master/qbism/sic_povms) research program up-to dimension $d=151$ |
| `"RandomMIC"` | A random MIC basis based on Haar method |
| `"RandomMIC"[Method->"Bloch"]` | A MIC basis based on random Bloch sphere affine transformation |

## Basic Examples

Create a random mixed state:

```wl
SeedRandom[0];
\[Rho] = QuantumState["RandomMixed"]
```

Transform into Tetrahedron basis and return amplitudes:

<!-- #| screenshot: true -->
```wl
amplitudes = QuantumPhaseSpaceTransform[\[Rho], "Tetrahedron"]["Amplitudes"]
```

Show probabilities from TetrahedronSICPOVM measurement:

<!-- #| screenshot: true -->
```wl
probabilities = QuantumMeasurementOperator[
   "TetrahedronSICPOVM"][\[Rho]]["Probabilities"]
```

Check they are the same:

```wl
probabilities == amplitudes
```
<!-- => True -->

---

Transform a quantum circuit:

```wl
QuantumPhaseSpaceTransform[QuantumCircuitOperator["CHSH"]]["Diagram", 
 "ShowWireDimensions" -> True]
```

Pay attention to the dimension of wire (which are doubled, compared to the conventional Hilbert space representation).

---

Generate a random mixed state:

```wl
SeedRandom[0];
\[Rho] = QuantumState["RandomMixed"]
```

Show Hadamard operator in the QBismSIC basis:

<!-- #| screenshot: true -->
```wl
QuantumPhaseSpaceTransform[QuantumOperator["H"], "QBismSIC"]["Table"]
```

Show the state (vectorized) in the QBismSIC basis:

```wl
QuantumPhaseSpaceTransform[\[Rho], "QBismSIC"]["AmplitudesList"]
```
<!-- => {0.39233282167545186, 0.27224439271969564, 0.09614688168225288, 0.23927590392259973} -->

Check that the transformation in the phase space returns the same state as the one in the Hilbert space:

```wl
QuantumWeylTransform@
  QuantumPhaseSpaceTransform[QuantumOperator["H"], "QBismSIC"]@
   QuantumPhaseSpaceTransform[\[Rho], "QBismSIC"] == 
 QuantumOperator["H"][\[Rho]]
```
<!-- => True -->

## Scope

### Calculate quantum probabilities from the phase space representation

Generate a random pure state in 2D:

```wl
SeedRandom[0];
\[Rho] = QuantumState["RandomPure"]
```

Calculate corresponding quantum probabilities in the computational basis:

```wl
prob = \[Rho]["ProbabilitiesList"]
```
<!-- => {0.8643837673262976, 0.13561623267370243} -->

Find the phase space representation in QBismSIC:

```wl
ps = QuantumPhaseSpaceTransform[\[Rho], "QBismSIC"]["AmplitudesList"]
```
<!-- => {0.3504166746438141, 0.3599603915103534, 0.28450628089946295, 0.0051166529463696694} -->

Find the transformation from QBismSIC to computational basis:

<!-- #| screenshot: true -->
```wl
trans = QuantumOperator["Measure", QuantumBasis["I", "QBismSIC"]];
trans["Table"] // Chop
```

Check `trans.ps` returns the quantum probabilities:

```wl
trans["Matrix"] . ps - prob // Chop
```
<!-- => {0, 0} -->

### Inverse basis in the phase space

Operators have input and output. When transforming an operator such as Hadamard into a phase space basis, both input and output will be the same as the given basis.

Transform Hadamard gate into QBismSIC basis:

<!-- #| screenshot: true -->
```wl
QuantumPhaseSpaceTransform[QuantumOperator["H"], "QBismSIC"]["Table"]
```

Some elements are negative, implying negative joint probability if one takes the transformation more seriously. To turn those negative (quasi) probabilities into positive values, one can use inverse basis:

<!-- #| screenshot: true -->
```wl
inverseQB = QuantumBasis["QBismSIC", QuditBasis["QBismSIC"]["Inverse"]];
QuantumPhaseSpaceTransform[QuantumOperator["H"], inverseQB]["Table"]
```

Inverse basis means that the direction of input is flipped:

```wl
#| style: Output
#| boxes: true
InterpretationBox[GraphicsBox[{EdgeForm[GrayLevel[0]], FaceForm[{GrayLevel[0], Opacity[0]}], GeometricTransformationBox[RectangleBox[NCache[{Rational[-1, 2], Rational[-1, 2]}, {-0.5, -0.5}], RoundingRadius -> {{Left, Top} -> 0.2}], {{{1, 0}, {0, 1}}, {0, 0}}], InsetBox[TemplateBox[{TagBox["\"H\"", HoldForm], RowBox[{"Diagram", "[", RowBox[{"\"H\"", ",", "a", ",", SuperscriptBox["b", "*"]}], "]"}]}, "ClickToCopy2"], {0, 0}], {Arrowheads[Small], {ArrowBox[{{0., 0.5}, {0., 0.75}}], InsetBox[TemplateBox[{InterpretationBox[TagBox["a", HoldForm], ProcessTheory`Port`Port[Association["Expression" :> $CellContext`a, "Type" -> ]], Tooltip -> ""], RowBox[{"ProcessTheory`Port`Port", "[", "a", "]"}]}, "ClickToCopy2"], {0., 1.}]}, {ArrowBox[{{0., -0.5}, {0., -0.75}}], InsetBox[TemplateBox[{InterpretationBox[TagBox["b", HoldForm], ProcessTheory`Port`Port[Association["Expression" :> $CellContext`b, "Type" -> ]], Tooltip -> ""], RowBox[{"ProcessTheory`Port`Port", "[", "b", "]"}]}, "ClickToCopy2"], {0., -1.}]}}}, BaseStyle -> {GraphicsHighlightColor -> RGBColor[1, 0, 1]}, FormatType -> StandardForm, ImageSize -> Tiny], Diagram[Association["DiagramOptions" -> {}, "Expression" :> "H", "InputPorts" -> {ProcessTheory`Port`Port[Association["Expression" :> $CellContext`b, "Type" -> ]]}, "OutputPorts" -> {ProcessTheory`Port`Port[Association["Expression" :> $CellContext`a, "Type" -> ]]}]]]
```

Show how the inverse basis is related the original basis:

<!-- #| screenshot: true -->
```wl
QuantumPhaseSpaceTransform[QuantumOperator["I"], inverseQB]["Table"]
```