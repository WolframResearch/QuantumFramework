# AI-Actionable Audit of the QuantumInDiscreteSpace Answers

Date: July 24, 2026

Status: Audit only. The ten `Part-*.md` source files were read in full. No source file was revised as part of this audit.

## AI execution brief

Goal: Revise the answer collection so each question is clear, each answer is concise and precise, jargon is explained where it first appears, and the pedagogy follows a computation-first rhythm.

Reader: A future AI agent editing the files in this directory.

Central instruction: Preserve the paired bare Wolfram Language and QuantumFramework computations, but make every answer solve exactly the task named in its question before offering any extension.

Scope:

- In scope: question wording, answer scope, conceptual precision, inline terminology, prose-to-code bridges, interpretation of results, code readability, and consistency across the ten Markdown files.
- Out of scope: redesigning the course, adding unrelated topics, replacing the WL/QF comparison, or changing notebook-generation infrastructure.
- This audit did not execute the code. Any agent that changes Wolfram Language or QuantumFramework code must run the relevant cells and verify their outputs.

Required skills for a revision pass:

1. Use `learning-by-computing` for pedagogy and prose structure.
2. Use `qf` when changing or assessing QuantumFramework code.
3. Use `wl-verify` after modifying Wolfram Language or QuantumFramework code.

Voice:

- Direct, mentor-like, and computation-first.
- Short paragraphs.
- Plain language before specialized terminology.
- First person only for genuine authorial choices; use "we" for the shared computation.
- No em dash or en dash.

Success test:

- A reader can identify the direct answer within the first two sentences.
- Every question is fully answered by an explicit computation.
- No answer silently expands into a second, harder topic.
- Every technical term needed by a BSc reader is glossed inline.
- Each code cell has a clear prose bridge.
- Nontrivial output is interpreted by meaning, not repeated verbatim.
- Changed cells evaluate successfully and support the surrounding claims.

## What must be preserved

Do not remove the collection's strongest features:

- The side-by-side WL and QF routes.
- Exact symbolic computations where they improve understanding.
- Explicit checks that the two routes agree.
- Concrete examples before generalization.
- The distinction between a bare matrix and a structured quantum object, explained once and then used without repeated promotional boilerplate.
- The strongest compact answers, especially:
  - `Part-04.md`, section 4.1.
  - `Part-05.md`, section 5.3.
  - `Part-06.md`, sections 6.3 and 6.4b.
  - `Part-07.md`, sections 7.2 through 7.5a.
  - `Part-10.md`, sections 10.2a, 10.2b, and 10.3.

## Revision rules for every answer

Apply these rules before addressing the file-specific backlog.

### One question, one primary task

The first example must completely answer the question. Move optional qudit, heterogeneous-register, degenerate-case, or group-theory extensions into a new question with their own level label.

Do not let a BSc answer become an MSc answer halfway through.

### Use the five-beat computation-first rhythm

For each nontrivial idea:

1. State the direct answer in two to four sentences.
2. Restate the difficult point plainly when it materially helps.
3. Introduce the computation with an imperative bridge ending in a colon.
4. Compute one clear fact.
5. Interpret what the result means and connect it to the question.

A bare Boolean verification may stand without further explanation when the bridge already states the claim.

### Define jargon inline

On first use, gloss terms by what they do. Terms that currently need attention include:

- sesquilinear
- gauge freedom
- affine
- extreme point
- commutant
- Bargmann invariant
- Casimir
- adjoint representation
- Haar-random
- arity
- Fubini-Study geodesic

Expand `WL` as Wolfram Language and `QF` as QuantumFramework once in every file intended to stand alone, or state clearly that a shared introduction defines them.

### Keep code cells self-contained and readable

- Add `ClearAll[name]` before defining or redefining reusable symbols.
- Prefer `With` or `Module` for local teaching examples.
- Do not use `%` to refer to a previous output.
- Break very long one-line inputs into readable multi-line expressions.
- Give each code cell one prose bridge.
- Define a reusable object once, then test it immediately.
- State symbolic assumptions explicitly.
- Explain `Chop` the first time it appears in the collection.

### Interpret meaning, not raw output

Avoid prose such as:

- "Both return `True`."
- "The cell returns `False`."
- "Both return $0$."

Prefer the physical or mathematical conclusion:

- "This confirms that the probabilities are unchanged by a global phase."
- "No physical Bloch vector violates the inequality."
- "The residual vanishes, so both sides of Ehrenfest's equation agree identically."

### Remove repeated framework boilerplate

Explain the following ideas once, early:

- A summary box is the live object, not a printed matrix.
- A `QuantumState` or `QuantumOperator` stores dimensions, basis, and subsystem order.
- A matrix is one representation read from the object.

After that introduction, do not repeat phrases such as "not a bare matrix," "one representation read from the object," or "the object knows its structure" unless that distinction directly answers the current question.

## Priority backlog

Complete P0 before P1, and P1 before P2.

## P0: Promise, precision, and evidence problems

### P0-01: Compute the collapsed states requested by Part 5, section 5.2

Location: `Part-05.md:114-200`

Problem: The question asks for post-measurement collapsed states. The answer computes amplitudes and probabilities, then asserts that the state becomes $|+\rangle$ or $|-\rangle$. It never retrieves or verifies the conditional post-measurement states.

Required action:

- Add an explicit WL computation of the normalized projected states.
- Add the corresponding QF property or operation for outcome-conditioned states.
- Verify that the states are the $X$ eigenstates, up to global phase.
- Keep the qutrit Fourier-basis extension only if it becomes a separate question.

Done when: The answer visibly computes both the outcome probabilities and the state associated with each outcome.

### P0-02: Match the scope of the driven-system question in Part 7, section 7.8

Location: `Part-07.md:291-329`

Problem: The question says "solve" the classically driven two-level system. The answer derives the rotating-wave Hamiltonian in WL, but the QF route only computes the resonant population. The general detuned solution is not shown.

Required action: Choose one of these scopes and make the question and answer agree:

1. Derive the rotating-frame Hamiltonian under the rotating-wave approximation.
2. Solve the effective detuned Hamiltonian and compute the transition probability.

If both are retained, split them into two questions.

Done when: Every operation promised in the heading is explicitly computed in both the main derivation and the final conclusion.

### P0-03: Recover the global phase in Part 9, section 9.4

Location: `Part-09.md:119-153`

Problem: The answer states the decomposition

$$
U=e^{i\phi}R_z(\alpha)R_y(\beta)R_z(\gamma)
$$

for an arbitrary $U(2)$ matrix, but samples an $SU(2)$ matrix and extracts only the three Euler angles. The global phase $\phi$ is never computed.

Required action:

- Either restrict the question to an $SU(2)$ matrix, or compute $\phi$ from a general $U(2)$ matrix before extracting the $SU(2)$ part.
- Reconstruct the original matrix from the extracted parameters and verify equality.
- State the Euler-angle branch convention.

Done when: The computed parameters reconstruct the same input unitary, including its global phase.

### P0-04: Correct the Hermiticity claim in Part 2, section 2.3b

Location: `Part-02.md:129-151`

Problem: "Hermiticity is exactly the condition that forces those eigenvalues to be real" is too strong. Hermiticity guarantees a real spectrum, but a non-Hermitian matrix can also have real eigenvalues.

Required action:

- State Hermiticity as a sufficient structural condition that guarantees real eigenvalues and an orthonormal eigenbasis.
- Do not describe real eigenvalues alone as equivalent to Hermiticity.
- Keep the existing symbolic eigenvalue verification if it remains clear.

Done when: The answer no longer implies that every matrix with real eigenvalues is Hermitian.

### P0-05: State the sign assumption or use the absolute value in Part 2, section 2.13

Location: `Part-02.md:640-679`

Problem: The minimum gap is described as $2\Delta$ while the code assumes only that $\Delta$ is real.

Required action:

- Use $2|\Delta|$, or explicitly assume $\Delta\ge0$ in both prose and computation.

Done when: The stated minimum gap is correct over the declared parameter domain.

### P0-06: Normalize the general state in Part 1, section 1.8

Location: `Part-01.md:252-297`

Problem: The general vector $\{\alpha,\beta\}$ is used directly in Pauli expectation values, although the result is called a Bloch vector of a pure qubit. This requires normalization.

Required action:

- State $|\alpha|^2+|\beta|^2=1$, or normalize the vector in the computation.
- Verify that the Bloch-vector norm is one for the pure normalized state.

Done when: The calculation cannot be misread as assigning a unit Bloch vector to an arbitrary unnormalized pair.

### P0-07: Qualify the mixture claim in Part 3, section 3.7

Location: `Part-03.md:232-280`

Problem: The conclusion says every mixture has radius and purity below one. Trivial mixtures and mixtures of identical pure states remain pure.

Required action:

- Restrict the statement to nontrivial mixtures of distinct pure states.
- Make the conditions $0<p<1$ and $\rho_1\ne\rho_2$ explicit.
- Connect radius-one boundary points to extremality rather than treating radius one alone as a proof.

Done when: The conclusion handles trivial and identical-state mixtures correctly.

### P0-08: Make the Berry-phase branch convention explicit in Part 7, section 7.13

Location: `Part-07.md:504-530`

Problem: The heading defines $-\arg\prod_k z_k$, while the code computes $-\sum_k\arg z_k$. These agree only modulo $2\pi$ after branch handling.

Required action:

- Compute the argument of the complete overlap product, or explicitly reduce the summed phase modulo $2\pi$.
- State the chosen branch or phase interval.
- Explain in plain language that the result is a phase defined modulo $2\pi$.

Done when: The formula, implementation, and reported phase use the same convention.

### P0-09: Use or remove the magnetic-field model in Part 8, section 8.1

Location: `Part-08.md:7-50`

Problem: The question asks to model a spin in a magnetic field. The answer introduces $\hat H=-\gamma B S_z$ but never uses the Hamiltonian.

Required action: Either compute its energies or time evolution, or narrow the question to Stern-Gerlach measurement of $S_z$.

Done when: Every major object introduced in the opening is used to answer the heading.

### P0-10: Compute the positivity condition in Part 3, section 3.6

Location: `Part-03.md:198-231`

Problem: The answer asserts that positivity is exactly $|\vec r|\le1$, but the computation shown recovers the Bloch vector and purity, not the eigenvalues that establish positivity.

Required action:

- Compute the eigenvalues $\frac12(1\pm|\vec r|)$.
- Use them to derive the positivity condition.
- Then connect $|\vec r|=1$ to purity one.

Done when: Positivity and the pure-state boundary both follow from computations shown in the answer.

## P1: Overscoped or level-mismatched answers

### P1-01: Split Part 8, section 8.2

Location: `Part-08.md:51-148`

Problem: A BSc construction of spin-$j$ matrices expands into spherical tensors, Wigner-Eckart, Clebsch-Gordan coefficients, product representations, and direct-sum decompositions.

Required action:

- Keep section 8.2 focused on building $J_x,J_y,J_z$, verifying the commutators, and verifying $J^2=j(j+1)I$.
- Move lines 78-127 into a separate MSc question on recovering angular-momentum matrix elements from Wigner-Eckart or Clebsch-Gordan coefficients.
- Introduce "Casimir" plainly as the operator $J^2$ that has the same value $j(j+1)$ throughout one spin-$j$ multiplet.

Done when: A BSc reader can complete section 8.2 without first learning spherical-tensor representation theory.

### P1-02: Split the heterogeneous-register extension from Part 4, section 4.2

Location: `Part-04.md:63-224`

Problem: The Bell-state partial trace is clear and sufficient, but the answer then adds a long $(3,2,5)$ register derivation.

Required action:

- End the BSc answer after comparing Bell and product reductions.
- Move the heterogeneous register to a separate advanced question about partial traces of qudits with unequal dimensions.

Done when: Section 4.2 answers the two-party reduced-state question in roughly one conceptual arc.

### P1-03: Simplify Part 1, section 1.5

Location: `Part-01.md:108-177`

Problem: The direct count $2d-2$ is obscured by Jacobian ranks, real embeddings, density-matrix maps, and gauge terminology.

Required action:

- Lead with the degree-of-freedom count: $2d$ real amplitude parameters, minus one for normalization, minus one for global phase.
- Verify the qubit count with the Bloch-angle parametrization.
- Move the Jacobian-rank derivation to a separate MSc question if it is pedagogically important.

Done when: A reader encounters $2d-2$ and understands the two removed parameters before any differential-geometric machinery appears.

### P1-04: Split ordinary and degenerate compatibility in Part 2, section 2.8

Location: `Part-02.md:316-407`

Problem: One BSc answer handles commuting observables, non-degenerate spectra, degenerate eigenspaces, generic combinations, and the Bell basis.

Required action:

- Keep the main answer focused on checking $[A,B]=0$ and verifying a shared eigenbasis in the non-degenerate case.
- Move the degenerate case and the $A+cB$ construction into a separate MSc question.
- Define "degenerate" immediately as "more than one independent eigenvector with the same eigenvalue."

Done when: The basic compatibility test can be learned without the degenerate-case construction.

### P1-05: Split the qutrit extensions from Part 5, sections 5.1 and 5.2

Locations:

- `Part-05.md:3-113`
- `Part-05.md:114-200`

Problem: Each qubit answer becomes a second answer about qutrits.

Required action:

- Keep 5.1 focused on outcomes, probabilities, and the mean for one qubit measurement.
- Keep 5.2 focused on basis change, probabilities, and collapsed qubit states.
- Move dimension-$d$ or qutrit generalizations into separate questions.

Done when: Each section has one example family and one conclusion.

### P1-06: Split gate definitions from wire placement in Part 10, section 10.1

Location: `Part-10.md:3-96`

Problem: One answer covers controlled-$U$, four named gates, block decompositions, non-adjacent placement, spectator identities, and register scaling.

Required action:

- Keep 10.1 focused on constructing and comparing the named matrices.
- Create a separate question about placing a gate on non-adjacent wires.
- State the wire-ordering convention before displaying matrices.

Done when: The gate-construction answer no longer contains a second tutorial on operator placement.

### P1-07: Add plain-language orientation to the densest Part 6 sections

Locations:

- `Part-06.md:3-70`
- `Part-06.md:71-134`
- `Part-06.md:247-322`
- `Part-06.md:323-375`

Problem: The answers begin at a high symbolic level and ask the reader to hold many definitions before seeing a concrete case.

Required action:

- Start each section with one concrete pair such as $X$ and $Z$ before moving to general Bloch directions.
- Define "slack" as the left side minus the right side.
- Define "saturate" as "reach equality."
- Separate the general proof from the concrete example when both are retained.

Done when: The first computation can be understood without parsing the fully general six-angle expression.

### P1-08: Reformat the long Part 7 computations

Location: `Part-07.md`

Problem: Part 7 contains 36 lines longer than 160 characters, including one 694-character line. These cells are difficult to inspect, teach from, or debug.

Required action:

- Break nested `With` expressions across lines.
- Name Hamiltonians, states, observables, and intermediate quantities.
- Keep each line to one logical operation where practical.
- Do not change the mathematics merely to shorten a line.

Done when: A reader can identify the state, operator, assumptions, and returned quantity without horizontally scanning a single large expression.

## P2: Collection-wide consistency

### P2-01: Remove output-history dependencies

Locations:

- `Part-05.md:253`
- `Part-09.md:44`
- `Part-09.md:56`
- `Part-10.md:28`

Required action: Bind each prior result to a descriptive symbol and refer to that symbol instead of `%`.

Done when: Every answer still works if its cells are copied into a fresh kernel and evaluated in the displayed order.

### P2-02: Add symbol hygiene

Location: All ten files.

Problem: Reusable global symbols are repeatedly assigned without `ClearAll`.

Required action:

- Add `ClearAll` before reusable function definitions.
- Prefer local scoping for one-off examples.
- Avoid allowing definitions from an earlier question to affect a later one.

Done when: Evaluating one answer cannot silently change the result of another answer.

### P2-03: Stop narrating literal Boolean and zero outputs

Representative locations:

- `Part-01.md:104`
- `Part-02.md:106`
- `Part-02.md:127`
- `Part-02.md:149`
- `Part-02.md:314`
- `Part-03.md:95`
- `Part-06.md:60`
- `Part-07.md:216`
- `Part-07.md:283`
- `Part-08.md:144`

Required action: Replace output-token narration with the mathematical or physical meaning, unless the Boolean stands silently after an explicit verification bridge.

Done when: The prose does not say "returns `True`," "returns `False`," or "returns $0$" merely to repeat visible output.

### P2-04: Add missing bridges before code cells

Representative locations:

- `Part-08.md:14`
- `Part-08.md:38`
- `Part-08.md:42`
- `Part-10.md:46`

Required action: Add a terse imperative sentence naming the exact computation and ending with a colon.

Done when: No code cell appears without a reader knowing what fact it is meant to compute.

### P2-05: Fix the Part 6 numbering

Location: `Part-06.md:156-246`

Problem: Section 6.4 is followed by section 6.4b, with no 6.4a label.

Required action: Rename 6.4 to 6.4a or renumber the pair consistently.

Done when: Cross-references can identify both sections unambiguously.

### P2-06: Define the WL/QF comparison once

Location: All ten files.

Required action:

- If the files are standalone, add one compact local definition of `WL` and `QF`.
- If they are chapters of one assembled document, add the definition to the shared introduction and avoid repeating it in every answer.

Done when: A reader never has to infer what `WL` or `QF` means.

## File-by-file action map

Use this as a completion checklist after working through the priority backlog.

### Part 1

- Preserve sections 1.1 through 1.4 and 1.9a through 1.9b as compact examples.
- Simplify 1.5 and move its Jacobian treatment to an advanced question.
- Gloss "sesquilinear" in 1.6 or remove it if it adds no practical understanding.
- Add normalization to 1.8.
- Reduce repeated summary-box explanations after their first introduction.

### Part 2

- Shorten the operator-object tour in 2.1 after the commutator and anticommutator have been computed.
- Correct the Hermiticity claim in 2.3b.
- Preserve the two-method agreement in 2.5.
- Split 2.8 into ordinary and degenerate cases.
- Keep 2.9a and 2.9b linked, but reduce literal-output narration.
- State the domain of $\Delta$ correctly in 2.13.

### Part 3

- Preserve the focused construction in 3.1a and the direct purity computation in 3.1b.
- Make the pure-state criteria in 3.2 visibly equivalent without overexplaining each Boolean.
- Compute the eigenvalues needed for positivity in 3.6.
- Qualify the mixture statement and strengthen the extremality argument in 3.7.

### Part 4

- Preserve the opening and section 4.1 as models for the desired voice.
- End the basic partial-trace answer after the Bell and product-state comparison.
- Move the $(3,2,5)$ example into its own advanced section.

### Part 5

- Split qutrit extensions from 5.1 and 5.2.
- Add explicit conditional collapsed states to 5.2.
- Preserve 5.3 as a concise model.
- Replace `%` in 5.4 with a named result.

### Part 6

- Add concrete $X,Z$ examples before the fully general formulas.
- Define "slack," "saturate," and covariance plainly.
- Fix the 6.4 and 6.4b numbering.
- Preserve 6.3 and 6.4b as comparatively focused answers.
- Shorten the setup of 6.5 and 6.6 without removing the essential computation.

### Part 7

- Preserve the concise physical interpretations in 7.2 through 7.5a.
- Reformat long cells throughout the file.
- Align the scope of 7.8 with what is actually solved.
- Make the adiabatic rate condition in 7.12 explicit rather than relying only on three sample rates.
- Fix phase-branch handling in 7.13.
- Define specialized terms such as "commutant," "Bargmann invariant," and "Fubini-Study geodesic" on first use.

### Part 8

- Use or remove the Hamiltonian in 8.1.
- Split 8.2 into BSc matrix construction and MSc representation-theory material.
- Add prose bridges before the consecutive cells in 8.1.
- Gloss "adjoint representation" in 8.3.
- Keep 8.4 through 8.7 focused on the operation named in each heading.

### Part 9

- Replace all `%` dependencies.
- Define "Haar-random" and "arity" or use plainer language.
- Recover the global phase in 9.4 or narrow the question to $SU(2)$.
- In 9.5, map the actual $2\pi$ unitary into $SO(3)$ so the computation mirrors the prose directly.
- Preserve the explicit statement in 9.6 that the search demonstrates density, not the efficient Solovay-Kitaev recursion.

### Part 10

- Split gate construction from non-adjacent wire placement in 10.1.
- Add a bridge before the named block-matrix display.
- Replace `%` with a named controlled-$U$ matrix.
- Preserve 10.2a, 10.2b, and 10.3 as models of concise task-answer alignment.

## Recommended execution order for the next AI agent

1. Record the existing Git diff and do not overwrite unrelated user changes.
2. Read `learning-by-computing`, `qf`, and `wl-verify` completely.
3. Fix P0 issues one file at a time.
4. Execute and verify every changed code cell before changing prose that depends on its output.
5. Split the P1 sections without deleting advanced material; move it into clearly labeled new questions.
6. Apply the P2 consistency rules mechanically across all ten files.
7. Re-read every question followed only by its first and last paragraph. Confirm that they answer the same task.
8. Re-read every BSc section as a reader who knows linear algebra but not the advanced terminology.
9. Run a final search for the recurring anti-patterns:

   ```bash
   rg -n '%|returns `True`|returns `False`|return(s)? \$0\$|summary box|bare matrix|\x{2014}|\x{2013}' Part-*.md
   ```

10. Run Wolfram verification on all modified cells.
11. Render or rebuild the notebook outputs if these Markdown files are notebook sources.
12. Report changed sections, verification results, and any issue deliberately deferred.

## Final acceptance checklist

- [ ] Every heading names one primary task.
- [ ] Every answer gives the direct answer within its first two sentences.
- [ ] Every requested quantity or state is explicitly computed.
- [ ] Every nontrivial claim is either computed or carefully qualified.
- [ ] Every technical term needed for comprehension is glossed inline.
- [ ] Every code cell has a bridge.
- [ ] No cell depends on `%`.
- [ ] Reused symbols are cleared or locally scoped.
- [ ] Long cells are formatted for inspection.
- [ ] Result prose interprets meaning instead of echoing tokens.
- [ ] WL and QF are defined for the intended reading context.
- [ ] The paired WL/QF structure is preserved.
- [ ] BSc answers do not silently become MSc treatments.
- [ ] All changed code has been executed and verified.
- [ ] No em dash or en dash appears in prose.

## Audit limitation

This report evaluates question clarity, answer clarity, jargon, concision, conceptual precision, and computation-first pedagogy. It identifies several mathematical claims that require correction or stronger evidence, but it is not a complete runtime or numerical verification of the Wolfram Language and QuantumFramework code.
