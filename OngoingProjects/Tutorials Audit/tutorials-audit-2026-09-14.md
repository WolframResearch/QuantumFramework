# QuantumFramework tutorials audit

Date: 2026-09-14. Anchor: repo `main` at `d874f7db`, paclet 2.1.1. Every kernel check below ran against the repo paclet (`PacletDirectoryLoad` on the working tree), not the installed copy.

Scope: the 15 pages under `QuantumFramework/Documentation/English/Tutorials/`. All 15 notebooks, the 5 existing markdown sources and the one image were read in full. The 10 notebooks that had no markdown source now have one beside them (NotebookToMarkdown twins with tech-note frontmatter; details in section 5).

## 1. Bottom line

- Two tutorial generations coexist. Five pages (TimeEvolution, QuantumMachineLearning, SendingQueriesToIBMQPUs, IBMQuantumErrorMap, QPUServiceConnect) are current md2nb sources from 2026. The other ten are hand-authored notebooks last touched between 2025-09 and 2026-06, and most of them predate the 2.0 syntax change.
- The old list-form gate and state specs are not rejected on 2.1.1, they are parsed as literal arrays. `QuantumOperator[{"RZ", Pi/8}]` now has the matrix `{{"RZ"}, {Pi/8}}` and `QuantumState[{"Register", 2}]` has the amplitudes `{"Register", 2}`. Every cell that uses these forms produces a wrong object with no message. Bellstheorem and ExploringFundamentalsOfQuantumTheory are built on them; the flagship circuit of Tutorial and the QAOA cost layers of QuantumOptimization use them too.
- The TensorNetwork tutorial documents an API that no longer exists: the `"TensorNetwork"` property now returns a `Wolfram\`TensorNetworks\`` object and rejects graph options, and the four functions the page is built on are undefined symbols.
- QuantumObjectComposition cannot be run by any reader: it loads an unpublished helper paclet from a developer path that is not in the repository.
- The guide page lists only two tech notes (Tutorial, GettingStarted). The other 13 are reachable only through search or cross-links.
- Merges that make sense: GettingStarted into Tutorial; IBMQuantumErrorMap into SendingQueriesToIBMQPUs (the parameter table belongs in a service-connection reference page); QPUServiceConnect keeps its umbrella role but loses its verbatim copy of the IBM walkthrough. Splits that make sense: QuantumOptimization (491 code cells, eight topics, two embedded reference pages) and SecondQuantization (the only documentation of a subpackage with no reference pages).

## 2. Per-tutorial verdicts

| Tutorial | Purpose | Serves it | Runs on 2.1.1 | Linked from guide | Verdict |
|---|---|---|---|---|---|
| Tutorial ("Quantum Computation") | framework overview, hub linked from 17 reference pages | mostly | one silently wrong cell | yes | keep as hub; absorb GettingStarted; fix syntax |
| GettingStarted | install, load, one example | duplicate of Tutorial's last section | one failing cell | yes | merge into Tutorial |
| Quantumobjectabstraction | internals: the object hierarchy | yes (reference-like) | one recursion failure, one undefined variable | no | keep as "Framework internals" note |
| QuantumObjectComposition | function-composition graph via a helper paclet | not for readers | no (missing paclet, `NotebookDirectory[]`) | no | remove from the shipped docs |
| CircuitDiagram | the 24 "Diagram" options | yes, and it is the only place they are documented | yes | no | keep; copy the option tables into the QuantumCircuitOperator page |
| TensorNetwork | circuit as tensor network, contraction | no | no (API removed) | no | retire, or rewrite on `Wolfram\`TensorNetworks\`` |
| Bellstheorem | CHSH derivation and circuits | good pedagogy | silently wrong (8 list-form gates) | no | keep as application note; fix syntax |
| ExploringFundamentalsOfQuantumTheory | eraser, bomb, Hardy, switch | good pedagogy | silently wrong (20 list-form specs) | no | keep as application note; fix syntax |
| SecondQuantization | bosonic subpackage, all functions | yes, but it is carrying reference duty | yes; installs a dev paclet in cell 1 | no | split: reference pages plus a shorter note |
| QuantumOptimization | eight optimization topics | too much for one page | mixed; QAOA layers silently wrong; cells of 7 and 38 minutes | no | split into four notes plus reference pages |
| TimeEvolution | QuantumEvolve end to end | yes | one plotting bug (`Mesh`) | no | keep; fix `Mesh` |
| QuantumMachineLearning | circuit compiled to a trainable net | yes | yes | no | keep |
| QPUServiceConnect | OpenQASM, IBM, Braket umbrella | partly; IBM section duplicates SendingQueriesToIBMQPUs; notebook rebuilt from the .md on 2026-09-14 | Braket path unverified and known broken | no | keep as umbrella; de-duplicate |
| SendingQueriesToIBMQPUs | the IBM hardware walkthrough | yes | hardware cells need an account | no | keep; absorb the error-map prose |
| IBMQuantumErrorMap | one service request and its 9 parameters | it is a reference page in tutorial clothing | needs an account | no | merge into SendingQueriesToIBMQPUs; table into a service-connection reference page |

Notes per page, in the same order.

**Tutorial.** The circuit `{"X", 1, "CNOT" -> {3, 2}, {"R", θ, "YY" -> {2, 3}}, "SWAP", "SX" -> 3, ...}` evaluates without a message, but the `{"R", θ, "YY"}` element is now a 3-by-1 literal array, so the displayed circuit and the measurement that follows are wrong. It links `QuantumPartialTranspose`, which has no reference page. It has one `%` chain, `XXXX` placeholders in every metadata slot, and stale outputs. The magic-basis example at its end is a byte-for-byte duplicate of GettingStarted.

**GettingStarted.** Eleven cells. The `qc["TensorNetwork", EdgeLabels -> Automatic]` cell now fails (`OptionValue::nodef`), and the stored `Names["Quantum*"]` output lists symbols that no longer exist (QuantumDiagramProcess, QuantumLabelName, QuantumStateSampler) while missing the ones that do. Its install instruction is the only content not already in Tutorial.

**Quantumobjectabstraction.** `QuantumChannel[{"BitFlip", p}, {2, 3}]` hits the recursion limit; the call form `"BitFlip"[p]` works. The cell `QuantumCircuitOperator[<|"Elements" -> ops, ...|>]` refers to an undefined `ops`. The categorization cell names the context `Wolfram\`QuantumFrameworkLoader\``. Otherwise the InputForm walkthrough still runs and is genuinely useful as an internals note.

**QuantumObjectComposition.** `PacletDirectoryLoad["~/QuantumFramework/OngoingProjects/Improving doc pages/PacletFunctionGraph"]` points at a path that exists on no reader's machine and is not tracked in this repository; `Names["Wolfram\`PacletFunctionGraph\`*"]` is empty on the repo kernel. The `$Meta` output says paclet 2.0.0. This is a design document for a tool, not a tutorial for the paclet.

**CircuitDiagram.** The QuantumCircuitOperator reference page mentions one of the 24 diagram options ("WireLabels"), so this note is the de facto option reference. All 20 cells run. The notebook's own metadata has "New in: ??".

**TensorNetwork.** On 2.1.1, `circuit["TensorNetwork"]` returns a `TensorNetwork` object from the TensorNetworks paclet; `GraphLayout` and `EdgeLabels` are rejected; `TensorNetworkQ`, `ContractTensorNetwork`, `TensorNetworkIndexGraph`, `TensorNetworkFreeIndices`, `TensorNetworkData` and `GraphTensorNetwork` return unevaluated. Its `QuantumState[{"RandomPure", 3}]` is a two-amplitude garbage state. The current API is the one QuantumMachineLearning already shows (`"TensorNetwork"` property, `GreedyContractionPath`, `TensorNetworkContraction`).

**Bellstheorem.** Eight `{"RZ", ...}` specs, two `%` chains, a first cell that fetches `EntityValue` over the network, links to `resources.wolframcloud.com` instead of paclet links, and the typo "voilated". The physics and the circuit-versus-analytic comparison are good and worth keeping.

**ExploringFundamentalsOfQuantumTheory.** Six `{"RY", θ}`, six `{"YRotation", θ}` or `{"Phase", ...}`, and eight `{"Register", n}` or similar specs, all silently wrong now; `{"Switch", A, B}` still works. The prose formulas were typed with the TeX assistant, one of them with a stray brace (`R_{y}(\theta})`), now fixed in the twin.

**SecondQuantization.** Its first cell is `PacletInstall["https://www.wolfr.am/DevWQCF", ForceVersionInstall -> True]`, which replaces the reader's installed paclet with a development build; that line must not ship. The function tables link `AnnihilationOperator`, `FockState`, `SetFockSpaceSize` and the rest to reference pages that do not exist (none of the subpackage's exports has one). One input cell carries a stray backspace character (`0 \.08`). Bullet lists were typed inside single text cells and come out as one paragraph. The applications (coherent-state decay, optical balance, Jaynes-Cummings with and without dissipation, driven Kerr oscillator) are the strongest part.

**QuantumOptimization.** Eight topics in one page: variational circuits and layer helpers, VQE with a PennyLane-generated H3+ Hamiltonian, adiabatic evolution (with a full Usage / Details / Scope / Applications / Possible Issues block for `QuantumAdiabaticEvolve`), QAOA plus QAOA-in-QAOA, parameter-shift rules and natural gradient, HHL, three VQLS examples, `QuantumLinearSolve` (again with Details / Options / Properties), and a Classiq setup section. An editorial note survives in the text: "(remove this, focus only on our convention, we do not want to say we reproduce exactly a specific result etc)". Three cells use `Quiet`. Seven `{"R", γ, "ZZ"}` specs feed the QAOA cost layers, so those sections compute with a wrong operator on 2.1.1. Stored timings show one gradient-descent cell at 446 s and another at 2286 s. Eighteen `%` chains. `QuantumLinearSolve` and every other export of the QuantumOptimization context lacks a reference page.

**TimeEvolution.** Current and thorough. One defect: `Plot[..., Mesh -> First[points], ...]` passes a flat list, which `Plot` rejects (`Mesh::ilevels`, visible in the built notebook's output); `Mesh -> {First[points]}` works.

**QuantumMachineLearning.** Current, focused, no findings.

**QPUServiceConnect.** The IBM section is a near-verbatim copy of SendingQueriesToIBMQPUs. The built notebook predates the last source edit: it lacks the "Inspecting a backend" section and still holds leftover debugging cells (`braket["Requests"]["SearchDevices"]["GetDevice"]`) with failed outputs. The Braket section is marked non-evaluatable and the end-to-end Braket path was found broken in an earlier audit.

**SendingQueriesToIBMQPUs.** Current. It is the page a reader with an IBM account actually needs.

**IBMQuantumErrorMap.** Five non-evaluatable cells and a static image. Two of its four sections are a parameter table and an accessor list for one `ServiceConnect` request; the other two ("Reading the map", "Why the bonds are bidirectional") are two paragraphs of prose that belong next to the submission walkthrough.

## 3. Cross-cutting findings

Syntax drift in the legacy pages, counted in the code cells of the twins.

| Tutorial | `{"RX/RY/RZ", ..}` | `{"R", ..}` | `{"YRotation"/"Phase", ..}` | `{"Register"/"RandomPure", ..}` | channel list form | `%` chains | `Quiet` |
|---|---|---|---|---|---|---|---|
| Tutorial | 0 | 1 | 0 | 0 | 0 | 1 | 0 |
| Quantumobjectabstraction | 0 | 0 | 0 | 0 | 1 | 0 | 0 |
| TensorNetwork | 0 | 0 | 0 | 1 | 0 | 1 | 0 |
| Bellstheorem | 8 | 0 | 0 | 0 | 0 | 2 | 0 |
| ExploringFundamentalsOfQuantumTheory | 6 | 0 | 6 | 8 | 0 | 3 | 0 |
| SecondQuantization | 0 | 0 | 0 | 0 | 0 | 1 | 0 |
| QuantumOptimization | 0 | 7 | 0 | 0 | 0 | 18 | 3 |

The other eight tutorials have none of these.

Kernel probe of the forms in question (repo paclet, fresh kernel).

| Form | Result on 2.1.1 |
|---|---|
| `QuantumOperator[{"RZ", Pi/8}]["Matrix"]` | `{{"RZ"}, {Pi/8}}`, no message |
| `QuantumOperator["RZ"[Pi/8]]["Matrix"]` | correct diagonal phases |
| `{"RZ", Pi/8} -> 2` versus `"RZ"[Pi/8] -> 2` inside a circuit | not equal |
| `QuantumCircuitOperator[{{"R", g, "ZZ"} -> {1, 2}}]["Matrix"]` | 3-by-1 array containing a QuantumOperator |
| `QuantumCircuitOperator[{"R"[g, "ZZ"] -> {1, 2}}]["Matrix"]` | 4-by-4 |
| `QuantumOperator[{"YRotation", t}]`, `QuantumOperator[{"Phase", 0.3}]` | 2-by-1 literal arrays |
| `QuantumState[{"Register", 2}]["AmplitudesList"]` | `{"Register", 2}` |
| `QuantumState["Register"[2]]["AmplitudesList"]` | `{1, 0, 0, 0}` |
| `QuantumState[{"RandomPure", 3}]` | two-amplitude state `{"RandomPure", 3}` |
| `QuantumOperator[{"Switch", A, B}]` | works (4-by-4, labelled) |
| `QuantumChannel[{"BitFlip", p}, {2, 3}]` | recursion-limit failure; `"BitFlip"[p]` works |
| `qc["TensorNetwork", EdgeLabels -> Automatic]`, `..., GraphLayout -> ...` | `OptionValue::nodef`, failed |
| `ContractTensorNetwork`, `TensorNetworkIndexGraph`, `TensorNetworkFreeIndices`, `TensorNetworkQ` | undefined, return unevaluated |
| `Names["Wolfram\`PacletFunctionGraph\`*"]` | empty |
| `Plot[..., Mesh -> {1., 2.}]` versus `Mesh -> {{1., 2.}}` | first fails (`Mesh::ilevels`), second works |

Other cross-cutting points.

- Reachability: the guide's Tech Notes section lists Tutorial and GettingStarted only. Symbol pages link Tutorial (17 pages) and IBMQuantumErrorMap (IBMJobSubmit). Nothing links TimeEvolution, QuantumMachineLearning, SecondQuantization or QuantumOptimization except each other.
- Metadata: seven legacy pages have `XXXX` keywords, five carry the stale context `Wolfram\`QuantumFrameworkLoader\``, all ten have `XXXX` in their Related Guides and Related Tech Notes slots, and the categorization alternates between "Tech Note" and "Tutorial".
- Dead links: the reference-page links from SecondQuantization (whole subpackage), QuantumOptimization (`QuantumLinearSolve`, the `ExampleRepositoryFunctions` tutorial, a Wolfram Cloud paclet-repository URL for itself), TensorNetwork (three removed functions) and Tutorial (`QuantumPartialTranspose`). The classification of every symbol link is in section 6.
- Reference-page material sitting inside tutorials: CircuitDiagram (Diagram options), SecondQuantization (function tables, orderings, error messages), QuantumOptimization (`QuantumAdiabaticEvolve` and `QuantumLinearSolve` blocks), IBMQuantumErrorMap (request parameters).
- Two pages with authoring residue: the editorial note in QuantumOptimization, and the debugging cells in the built QPUServiceConnect notebook.

## 4. Recommended target set

| Today (15) | Proposed | Action |
|---|---|---|
| Tutorial + GettingStarted | Getting Started (one page) | merge; install section first, then the overview; delete the duplicated magic-basis block; replace the list-form circuit |
| Quantumobjectabstraction | Framework Internals | keep; fix the channel form and the `ops` cell; fill metadata |
| QuantumObjectComposition | (none) | move out of the shipped docs until the helper paclet ships |
| CircuitDiagram | Circuit Diagrams | keep; copy the option tables into the QuantumCircuitOperator page's Details |
| TensorNetwork | Tensor Networks | rewrite on the TensorNetworks paclet API, else retire |
| Bellstheorem, ExploringFundamentalsOfQuantumTheory | two application notes under one "Applications" group | keep both; convert every list-form spec to the call form |
| SecondQuantization | Second Quantization note + reference pages for its exports | split; drop the `PacletInstall` cell; move tables, orderings and messages to the reference pages |
| QuantumOptimization | Variational Circuits and Optimizers; VQE and QAOA; Adiabatic Quantum Computing; Quantum Linear Solvers | split; reference pages for `QuantumLinearSolve`, `QuantumAdiabaticEvolve`, `GradientDescent`, `QuantumNaturalGradientDescent`, `FubiniStudyMetricTensor`, the layer helpers and the shift-rule functions; remove the editorial note, the `Quiet` calls and the 38-minute cell; Classiq setup becomes a reference page for `ClassiqSetup` |
| TimeEvolution | Time Evolution | keep; fix `Mesh` |
| QuantumMachineLearning | Quantum Machine Learning in Phase Space | keep |
| QPUServiceConnect | Running on Quantum Hardware (umbrella) | keep the OpenQASM section; replace the IBM section by a short pointer to the IBM page; label Braket as unverified or remove it until the pipeline works |
| SendingQueriesToIBMQPUs + IBMQuantumErrorMap | Sending Queries to IBM QPUs | merge the map prose in as "Inspecting a backend before you submit"; the parameter table and accessor list go to a `Template: ServiceConnection` reference page for `IBMQuantumPlatform` |

Result: 10 tutorial pages plus new reference pages, every page reachable from the guide's Tech Notes list.

## 5. What is now on disk

Ten new markdown sources beside their notebooks, produced with `UsingFrontEnd @ NotebookToMarkdown[nb, md]` (faithful front-end input text for every code cell; code-cell counts match the notebooks cell for cell, the six image-only cells of QuantumOptimization becoming figure links to six `QuantumOptimization-fig-N.png` sidecars).

Added to each twin, since the walker recovers no frontmatter for tech-note notebooks: `Template: TechNote`, `Name`, `Title` (the notebook's title cell), `Context: Wolfram\`QuantumFramework\``, `CellContext: Global\``, `Paclet`, `URI` (the notebook's own), `Keywords` where the notebook had any, `RelatedGuides: [WolframQuantumComputationFramework]`, and `RelatedTutorials` where the notebook had links (QuantumObjectComposition). The duplicate H1 and the unfilled "Related Guides / Related Tech Notes / XXXX" scaffolding were dropped. Everything else is the walker's output.

Hand edits, each asserted by the finalizer so a regenerated twin that no longer matches is noticed:

- Quantumobjectabstraction: one inline signature that the walker dumped as a `ButtonBox`.
- QuantumObjectComposition: the Kind table (a grid typed inside a text cell) rewritten as a pipe table.
- ExploringFundamentalsOfQuantumTheory: the author's `R_{y}(\theta})` corrected to `R_{y}(\theta)` (5 places).
- Bellstheorem: raw `\[Implies]` glyphs in prose written as `$\Rightarrow$`.
- SecondQuantization: the two multi-paragraph "Ordering" table cells and the Jaynes-Cummings display formula rewritten as TeX.
- QuantumOptimization: the two VQLS cost-function display formulas rewritten as TeX.
- After a KaTeX pass over every math span of all fifteen files (683 spans): in QuantumOptimization, the author's typed TeX carried a stray `&` in a two-part display formula (now a comma and `\quad`), `\ij` for a subscript `ij`, `\ψ` for `\psi`, two raw backspace characters inside the QAOA oracle formula, unmapped `\[VerticalSeparator]` bars in the HHL amplitude ratio (now `|x_1|^2:|x_2|^2`), `\overset{^\b}{x}` for `\hat{x}`, and a Unicode `√` for `\sqrt` (twice); in SecondQuantization, a double superscript `a^{\dagger}^{2}` (now `a^{\dagger 2}`) and an em-dash overscript standing in for `\bar{n}` (twice). The formal-symbol default `s` in the `QuantumAdiabaticEvolve` option table is written `\[FormalS]`; `$Failed` and `$Meta` in prose are backticked so no markdown parser reads them as math delimiters; a page-range en dash in a Bellstheorem reference is a hyphen.
- After a line-by-line read of the QuantumOptimization twin: two images that sat inside text cells (the natural-gradient update rule and the multiplexer-solver flow chart) are dropped by the walker, which exports only images in input cells; they were extracted from the notebook and are now `QuantumOptimization-fig-7.png` and `QuantumOptimization-fig-8.png`. Three Python cells (the PennyLane Hamiltonian, the Classiq authentication, the session example) came out as one-line prose with a style directive; they are `python` fences now (md2nb builds a foreign-language fence as a plain program cell, which is right for cells that need Python anyway). One sentence carried four `\[VerticalSeparator]` bars as raw named characters (now math), two formulas carried WL precision marks copied from output, an empty input cell was dropped, and a code-styled comma between two formulas is a plain comma. In SecondQuantization the squeeze-operator definition row and a coherent-state bra-ket were typed with raw named characters (now math). In TensorNetwork the package-scope context was typed inside math delimiters (now a code span).
- Left as the notebook has them: twenty-seven `$CellContext\`` prefixes inside QuantumOptimization input cells (the front end's own input text for symbols pasted from evaluated graphics; they resolve inside a notebook), and one of those cells references an undefined `qparameters` where `quantump` was meant, an author error the notebook already carries.
- Processing directives, not content: `#| eval: false` on the cells a rebuild must not run, namely the two `PacletInstall` cells (GettingStarted, SecondQuantization), the missing-paclet `PacletDirectoryLoad` cell (QuantumObjectComposition) and the six `ClassiqSetup` cells (QuantumOptimization), so `build_docs.wls` cannot install anything on a reader's machine.

Checks run on all fifteen `.md` files after these edits: code-cell count against each notebook (equal, the six image-only cells aside), no box dumps, no private-use or control characters, no `XXXX` placeholders, every inline and display math span renders under KaTeX 0.16 with `throwOnError`, and the QuantumOptimization twin was also read end to end after the scans. Still present: em dashes in the authors' prose (QuantumObjectComposition 10, QuantumOptimization 25, SecondQuantization 1), left as written pending a decision.

Known residues in the twins: `<!-- #| style: MathCaption -->` directives (they round-trip the caption style, harmless); bullet lists typed inside single text cells flattened into one paragraph (SecondQuantization); four formulas the author typed inside parentheses-grids kept as `\begin{matrix}` (QuantumOptimization); typeset kets inside inputs appear as `Wolfram\`QuantumFramework\`QuditName[...]` (that is the front end's own input text for them); the stray backspace in a SecondQuantization input; `%` chains that md2nb cannot rebuild (counts above).

Rebuild (2026-09-14, after the first version of this report): all fifteen notebooks were rebuilt from their `.md` with `MarkdownToNotebook[md, nb, "Evaluate" -> False]`, one call per page in a single kernel (the local `MarkdownToNotebook.wl`, front-end appearance pinned to Light as `build_docs.wls` does). The notebooks are input-only: no stored outputs, nothing evaluated, and no verification loop was run, on instruction. Without outputs they are a fraction of their former size (QuantumOptimization from 12 MB to 1.9 MB, SecondQuantization from 7.7 MB to 239 KB). Fourteen notebooks changed; IBMQuantumErrorMap.nb came out byte-identical to the build already committed, so the build is reproducible. Three `FrontEndObject::notavail` messages appeared during the QPUServiceConnect build; the same page rebuilt with the whole conversion under `UsingFrontEnd` is byte-identical and message-free, so the messages changed nothing. The `%` chains (counts above) are now literal `%` in input cells and evaluate correctly only in document order. The pre-rebuild notebooks, with their stored outputs, are the copies in the backup folder (section 7).

Converter changes in `MarkdownToNotebook/NotebookToMarkdown.wl` (uncommitted, in a working tree that already carried other sessions' uncommitted edits):

- `DefinitionBox3Col` cells walk to three-column pipe tables (they fell through to matrix math).
- TeX-assistant formulas splice the author's typed TeX (the box tree was dumped; rendering the boxes instead turned `\mathrm{}` into bogus commands); the placeholder tokens are zero-padded after a verification round found that ten or more boxes in one formula collided.
- Nested style-less `Cell` leaves and `HyperlinkDefault`-wrapped links inside signatures render clean.
- `HyperlinkTemplate` web links walk to markdown links.

Five regression tests were added to `NotebookToMarkdown.md`. Test-suite result: 361 of 362 pass. The one failure, "nested Details heading groups notes (issue #77)", exercises the forward converter only (`MarkdownToNotebook`, whose working copy carried earlier uncommitted edits) and predates this session..

## 6. Symbol links without a reference page

Every `[Symbol]()` link in the fifteen pages was classified on the repo kernel. Links to existing paclet reference pages and to System symbols are omitted here.

| Tutorial | Linked symbol | Status |
|---|---|---|
| SecondQuantization | AnnihilationOperator, BeamSplitterOperator, CatState, CoherentState, DisplacementOperator, FockState, G1Correlation, G2Coherence, HusimiQRepresentation, OperatorVariance, PhaseShiftOperator, QuadratureOperators, SetFockSpaceSize, SqueezeOperator, ThermalState, WignerRepresentation | paclet symbols with no reference page |
| QuantumOptimization | QuantumLinearSolve | paclet symbol with no reference page |
| Tutorial | QuantumPartialTranspose | paclet symbol with no reference page |
| TensorNetwork | ContractTensorNetwork | undefined symbol |
| TensorNetwork | TensorNetworkQ, TensorNetworkFreeIndices | names survive in the paclet context but carry no definition |
| TensorNetwork | TensorNetworkIndexGraph | now lives in the TensorNetworks paclet, so the paclet-local link is dead |
| QuantumMachineLearning | TensorNetworkContraction | TensorNetworks paclet symbol; the link resolves only if md2nb targets that paclet |
| TimeEvolution | PauliX, PauliZ | undefined symbols; the source writes `[PauliX]()` for what is the operator name string `"PauliX"` |
| QuantumObjectComposition | QuantumFramework | a guide link written as a symbol link |

## 7. Backup

All fifteen notebooks were copied, unchanged and checksum-verified, to `/Users/mohammadb/Documents/GitHub/Old QF notebooks/` before the rebuild; they are the pre-rebuild notebooks with their stored outputs.

## 8. Open items

- Rebuild done, unverified: all fifteen notebooks are `"Evaluate" -> False` builds of their `.md`; none has been opened, evaluated, or round-tripped, and the fourteen changed notebooks plus the ten twins and eight figures are uncommitted.
- QPUServiceConnect: its notebook was behind its `.md` (missing the "Inspecting a backend" section, still holding debugging cells); the rebuild closed that gap.
- The Braket path in QPUServiceConnect and the hardware cells of the IBM pages cannot be checked without credentials.
- The converter edits went through `/wl-verify`: round 1 found three issues (the suite count stated here was off by one, placeholder tokens collided once a formula held ten or more TeX-assistant boxes, and a non-string `input` produced a `StringExpression`); all three were fixed, a fifth regression test was added, and round 2 reproduced every claim with no open issues: the five tests pass from the shipped text, the suite is 361 of 362 with only the pre-existing forward-path failure, non-TeX formulas serialize byte-identically to the committed walker, and both regenerated pages (CircuitDiagram, ExploringFundamentalsOfQuantumTheory) carry no box dumps and raise no messages.
- Non-blocking observations left as they are: the token width wraps only past 10000 boxes in one formula; a `$` inside an author's typed TeX would close the math span early; the definition-table rule has no divider guard, so a bare two-by-two grid in a `DefinitionBox` becomes a table rather than a matrix; the two new inline-prose `Cell` rules rely on dispatch order (not an explicit guard) to leave decoration cells alone.
- The diagnostic probe scripts in the session scratchpad were not verified and are throwaway.
