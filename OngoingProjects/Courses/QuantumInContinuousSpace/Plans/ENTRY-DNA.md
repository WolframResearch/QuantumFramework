# Entry DNA: how a plan entry shows the computation

Governs the shape of every per-question entry in `Part-NN-Plan.md`. Distilled from a single worked
case, 4.1, and grounded in a live kernel (WL 15.0.0) rather than recalled. `PIPELINE.md` remains the
standard for the finished answers; this file is the standard for the plans that generate them, and
the two now have the same shape.

## The principle

**Every physical conclusion is a return value, not a sentence.**

If an entry states that a constant vanishes, a branch is inadmissible, a spectrum is
$E_n=n+\tfrac12$, or a state is normalized, then some cell in the entry returns that fact. Prose
introduces and interprets; it never supplies. The test is mechanical: take any physics claim in the
entry and point at the binding that produced it. If the answer is "the author works it out", the
entry is not finished.

This is `PIPELINE.md` section 0 ("earn every claim with a visible computation") applied one level
up, to the plan.

## The failure this replaces

The first draft of 4.1 read, in part: "select the decaying branch by hand and read the quantization
off as the manual termination condition, Solve giving $E_n=n+\tfrac12$."

Three defects, all of one kind. "By hand" hands the load-bearing step back to the reader. "Solve
giving" never says what equation `Solve` receives. And $E_n=n+\tfrac12$ arrives as an assertion in
prose. The entry named seven Wolfram functions and showed no computation, so a cold authoring
session would have had to invent the method at exactly the point where the physics lives.

## The pipeline, six moves

1. **Make the defining object applicable.** Bind the operator as a function of an expression, so
   every later cell can apply it: `h[f_] := -D[f, {x, 2}]/2 + x^2 f/2;`. One line, reused
   everywhere, never re-typed.
2. **Generic first.** Solve with the eigenvalue symbolic. The output is a family carrying free
   constants, which is the honest state of knowledge before boundary conditions.
3. **Admissible second, by independent machinery.** A route whose domain restricts the solution
   set, returning concrete pairs. This is where quantization becomes a return value. Check what the
   domain actually imposes: an unspecified boundary condition is not the same as decay.
4. **Confront the two, and let the confrontation return the reduction.** A computation whose output
   says which constants die. This is the move that converts "select the branch by hand" into a
   result.
5. **Infer the law from the computed data.** Turn the returned list into a closed form in $n$
   instead of asserting one, and state the index convention explicitly.
6. **Refute.** A residual on a nontrivial carrier, plus the asymptotic that shows why the rejected
   branch had to be rejected.

Moves 2 and 3 are co-generators, not primary-and-witness. This is the distinction
`Route-Table.md` does not carry: its verdicts pick a primary route plus an independent cross-check,
which is the right frame for *trusting* an answer and the wrong frame for *deriving* one. When a
second route already encodes a constraint, feed its output into a third computation rather than
saving it for confirmation at the end.

## The worked exemplar: 4.1

Every output below was reproduced in WL 15.0.0 by independent fresh-context verification rounds that
wrote their own scripts rather than rerunning the author's. The constant names depend on the
`GeneratedParameters -> C` default; under `GeneratedParameters -> K` the `Cases` pattern in move 3
returns nothing.

The oscillator Hamiltonian as an applicable operator, then the general solution at symbolic energy:

```wl
h[f_] := -D[f, {x, 2}]/2 + x^2 f/2;
generic = DSolveValue[h[u[x]] == en u[x], u, x]
```

Returns a `Function` whose body is
`C[2] ParabolicCylinderD[(-1-2 en)/2, I Sqrt[2] x] + C[1] ParabolicCylinderD[(-1+2 en)/2, Sqrt[2] x]`:
two constants, no quantization, every energy still allowed. Do not attach decay conditions at
$\pm\infty$ here. Ask what happens rather than being told:

```wl
AbsoluteTiming[DSolve[{h[u[x]] == en u[x], u[-Infinity] == 0, u[Infinity] == 0}, u, x]]
```

Returns `{27.8, {}}`: an empty solution set, so this fails rather than hangs, and the cost lands
between 27 and 30 seconds across isolated runs (the `u[x]` second-argument form costs about 43).
`DSolveValue` on the same system returns unevaluated. The cost is paid once per kernel by whichever
of the two runs first, after which the other returns in under a tenth of a second and
`ClearSystemCache[]` restores the full price, so neither form is the cheap escape. Decay has to
enter through the eigensolver's domain instead (C1 verdict, probe 3).

The recorded *hang* belongs to a different equation, the Coulomb radial problem with conditions at
$0$ and $\infty$, which does time out. Transplanting that verdict to this equation was an error this
document made and its own verification caught, which is the failure mode the checklist below exists
to prevent.

The admissible pairs, from the domain form. Note what the domain does and does not do: with no
boundary condition supplied, the eigensolver applies Neumann-zero on $\partial\Omega$, which is
harmless on the whole line but is not decay, so this phrasing must not be transplanted to a finite
domain without an explicit condition:

```wl
{energies, eigenfunctions} =
  DEigensystem[h[u[x]], u[x], {x, -Infinity, Infinity}, 8]
```

Returns `{1/2, 3/2, 5/2, 7/2, 9/2, 11/2, 13/2, 15/2}` with $e^{-x^2/2}$, $2x\,e^{-x^2/2}$, and so on.

The free parameters, found rather than counted:

```wl
constants = DeleteDuplicates@Cases[generic[x], C[_], Infinity]
```

Returns `{C[2], C[1]}`, in traversal order, not the order they were named.

Now the move that earns the physics. Match the generic family to each admissible pair at a point, in
value and slope, and let `Solve` report what survives:

```wl
matchingRules = Table[
   Solve[{(generic[0] /. en -> energies[[j]]) == (eigenfunctions[[j]] /. x -> 0),
          (D[generic[x], x] /. {x -> 0, en -> energies[[j]]}) ==
            (D[eigenfunctions[[j]], x] /. x -> 0)},
     constants],
   {j, Length[energies]}]
```

Returns `{{C[2] -> 0, C[1] -> 1}}`, `{{C[2] -> 0, C[1] -> Sqrt[2]}}`, `{{C[2] -> 0, C[1] -> 2}}`,
and on through `{{C[2] -> 0, C[1] -> 8 Sqrt[2]}}` at the eighth level, one solution per level, with
$C_1=2^{(j-1)/2}$ throughout. `C[2] -> 0` at every level, computed, not assumed: the branch that
grows like $e^{+x^2/2}$ is excluded by admissibility.

Two things this move does not do, both worth stating because the document's thesis invites the
overclaim. First, the surviving constant is not fixed by normalization. The eigensolver does not
normalize its output, so $C_1$ is fixed by matching an arbitrary scale. Do not take that on my word,
which is the whole point of this document; ask:

```wl
Table[Integrate[eigenfunctions[[j]]^2, {x, -Infinity, Infinity}], {j, 8}]
```

Returns `{Sqrt[Pi], 2 Sqrt[Pi], 8 Sqrt[Pi], 48 Sqrt[Pi], ...}`, that is $\|f_n\|^2=2^n n!\sqrt{\pi}$,
squared norms and not norms. `Method -> "Normalize"` is the documented option that would make a
normalization reading true, and it must be confirmed by the integral rather than by its own
presence: the eigensolver does not validate `Method` strings, so `Method -> "Normalise"` or any
other misspelling returns unnormalized output silently, with no message. Second, the matching
substitutes `energies[[j]]` before solving, so it
cannot return discreteness; it returns branch death at energies that were already quantized one move
earlier. Quantization is a return value of move 3, not of move 4.

The spectrum as a law rather than a claim:

```wl
FindSequenceFunction[energies, n]
```

Returns `(-1 + 2 n)/2`, that is $n-\tfrac12$ on the 1-indexed list. State the offset in the same
breath: the Wolfram index is the quantum number plus one, so $E_n=n+\tfrac12$ for
$n=0,1,2,\dots$. An entry that quotes $n+\tfrac12$ without the offset has skipped a real step.

Refutation, on a carrier that is not the ground state:

```wl
FullSimplify[h[(2 x^3 - 3 x) Exp[-x^2/2]] - (7/2) (2 x^3 - 3 x) Exp[-x^2/2]]
```

Returns the exact integer `0`, not a residue that merely vanishes numerically. `FullSimplify` is
doing real work here: the unsimplified difference is a nonzero expression, so the wrapper is part of
the claim rather than decoration. Pair it with the large-$x$ behavior of the rejected branch, which
is why admissibility killed it:

```wl
Asymptotic[ParabolicCylinderD[-1, I Sqrt[2] x], x -> Infinity]
```

Returns `-I E^(x^2/2)/(Sqrt[2] x)`, against $e^{-x^2/2}$ for the branch that survived. Use
`Asymptotic` rather than `Series` here: raw `Series` returns a `Plus` of `E^SeriesData[...]`
products in which the growth factor is not literally present until `Normal` is applied, so the
one-step form is both shorter and honest about what it shows.

## Checklist for an entry

- [ ] Every physics claim in the prose points at a cell that returned it. No "by hand", no "one
      then finds", no formula stated but not computed.
- [ ] Each cell consumes a binding from an earlier cell. No re-typed literals, no wavefunction
      written twice.
- [ ] Named functions appear as calls with their actual arguments, not as bare names in a sentence.
- [ ] Where two routes are available, it is stated whether the second is a co-generator (its output
      feeds a derivation) or a cross-check (it independently confirms). Do not default to the second.
- [ ] Any general law is inferred from computed data, or derived, never asserted.
- [ ] Index and sign conventions that differ between Wolfram and physics are stated where they occur.
- [ ] The refuting check is a computation whose failure mode is loud.

## Skills to invoke when drafting entries

`learning-by-computing` is the genre skill for this shape (short prose, then working code) and must
be loaded before drafting. `wl-native-style` and `wl-style-guide` govern the idiom. `pde-route`
decides routes and is not a drafting skill: its primary-plus-cross-check frame answers "can I trust
this result", not "how is this computed", and using it alone is what produced the prose-only first
draft. Before settling on an idiom, enumerate at least two and read the primitive's doc page; the
introspection and inference primitives (`Cases` for free parameters, `FindSequenceFunction` for a
law) are exactly what that search surfaces and what assertion-shaped prose hides.
