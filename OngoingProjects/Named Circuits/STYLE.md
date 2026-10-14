# Style contract for the named-circuit pages

Every family page follows this. It is the spec behind the Grover page.

## Format
- `Template: Default` md2nb, with `Context: Wolfram`QuantumFramework`` so cells
  need no `Needs`.
- Every example cell is `QuantumCircuitOperator[{"Name"[...]}]["Diagram"]`, one
  diagram per cell, siblings separated by `---`. Operators that are not circuits
  are embedded in the one-element list so the whole page is one object type.
- One-line intro; short section headers naming the layers/roles of the family,
  not the call-forms. Group signatures by role.

## Captions
- Exactly one line, stating what *this* signature does as distinct from its
  siblings: the differentiator (sets the wire count, the target wire, the gate,
  names the variables so no ancilla, picks the about-zero diffusion), never a
  generic restatement.
- Concrete, never vague. Name the actual wires, the actual gate, the actual
  width change. Banned: "reserving ancillas", "wrapped as given", "shifts every
  wire" and any phrase that hides the mechanic.
- No claimed identity that is not literally true. State only what the kernel
  shows (matrix, wire layout), not what the name "should" mean.

## What goes in / what is cut
- Exhaust the real signatures and downvalue branches, but every cell must build.
  Anything that returns a `Failure`, or a half-working alias, is cut from the
  page and logged in `ISSUES.md`.
- Prefer the clean illustrative input over the valid-but-murky one.
- Terse idiomatic calls: drop defensive placeholders (e.g. `Automatic`) when the
  shorter form is identical.

## Dual-registry names
- `Fourier`, `InverseFourier`, `Toffoli`, `Fredkin`, `Switch`, `Multiplexer`,
  `Deutsch` are in both registries. Under the `QuantumCircuitOperator[{...}]`
  convention they resolve to the circuit decomposition, so they are documented
  once (as circuits) and the caption says the diagram is the decomposition, not
  an atomic box.

## Discipline
- Every claim verified by running the kernel, not recalled. Pull each name's
  signatures from the kernel source (`NamedCircuits.m` / `NamedOperators.m`),
  then one verification script per family builds every example and prints width,
  wire layout (`["Flatten"]["Operators"]` orders), and matrix identity where it
  matters. One kernel boot per family, not per cell.
- Self-test each caption: if the differentiator cannot be stated in one true
  concrete sentence, the example is the wrong one to include.
