# Independent quality review, round 3/3

Reviewed the complete requested source, documentation and saved JSON artifacts under
`OngoingProjects/VerifiableCircuitCompilation/dfthc_gram`, plus the explicitly permitted
resource definitions, pinned inputs, official mathematical routines and saved candidates.
The reviewed `GramCores.wl` SHA-256 is
`fb1f66cf2373916551906e5b322930ffcda04980ac3a01171e69721fae5143b7`.
No prior reviews were consulted. Evidence commands ran on isolated copies; source artifacts,
incumbents and candidate files were not modified. No optimization pilots or historical
circuit runs were repeated.

Kernel: Wolfram 15.0.0 for Mac OS X ARM (64-bit), February 17, 2026.
`WolframKernel -noinit -noprompt -script run.wls` returned PASS, all four all-data
checks true, and SDP minimum `1.399999996255239` against exact `7/5`.
`WolframKernel -noinit -noprompt -script verify-data.wls` returned PASS on all ten
checks, 21 reference resource rows and one unchanged-loader row.
`python -m unittest test_study -v` passed all five tests.

**A. Documented symbols — PASS.** The only non-System WL APIs are the six exports defined
with usage strings and explicit definitions in `GramCores.wl:1–21,23–117`, and
`DFTHCFeasibility`PreparationCost` / `LookupCost` from
`../dfthc_feasibility/RankCost.wl`. The current source establishes arguments, guards,
return associations and failure shapes; executable probes returned core dimensions `{3,3}`,
rectangular projector-Gram dimensions `{2,3}`, gradient dimensions `{2,2}`, and a malformed
core returned `Failure`. All eight exported functions have empty `Options`, consistent with
the calls. The source-pinned resource definitions reproduced the complete official-model
rows. The primary official `dfthc_math.py` definitions for normalized projectors, native
offsets, one-body correction and LCU normalization agree with `README.md:10–31,40–59`.
No QuantumFramework symbols or undocumented package APIs occur.

**B. Idiomatic functional WL — PASS.** A search of all three WL targets found no
`Do`, `For`, `While`, `AppendTo` or `Quiet`. `MapThread` constructs block records
(`run.wls:17–29`); `Map` / `AllTrue` express checks and reductions. The `Table` uses
(`GramCores.wl:31`, `run.wls:37–38`) enumerate principal-minor sizes, not indices over
an existing data list. The algebraic `Module` scopes isolate free symbols as well as
intermediate expressions; they are not wrappers around lexical constants alone.
The dense matrices here are small symbolic fixtures or genuinely dense measured arrays.
`Attributes[Expand]` returned `{Protected}`, so the explicit expansion maps are not maps
of a Listable function.

**C. Terse physics reasoning — PASS.** Core blocks, projected residuals, copy-space
projectors, offset contractions and rescaling Hessians form a readable mathematical chain
(`GramCores.wl:27–69,76–117`). The Python boundary labels the factor-residual pilot,
global coupled residual, covariance and final official metrics separately
(`study.py:116–138,205–267,306–375`). The reports state the failed feasibility result,
iteration cap and conditional normalization directly; they do not promote rank reduction
to an objective improvement. No material dead scaffold or repeated physics derivation
requires a fix.

**D. Physics invariants — PASS.** The rerun passed all 37 exact identities, including PSD
quadratic form, two-body reconstruction, one-body contraction, offset norm and copy-space
invariants (`GramCores.wl:47–69,90–98,111–117`). Independent reads of all eight SHA-matched
saved candidates confirmed real finite factors, declared shapes and exact equality to the
vectors decoded from their saved angle words. Maximum unit-row error was `4.45e-16`,
maximum factor antisymmetry `4.17e-17`, and augmented-core antisymmetry zero. The most negative
computed core eigenvalue was `-1.25e-14`, consistent with roundoff in the exactly PSD
construction `A.Transpose[A]`. Independent C29 offset and norm drifts were
`8.881784197001252e-16` and `2.842170943040401e-14`; its complete word table equals the
incumbent's. `study.py:256–258,336–337` also enforces offset preservation and agreement
between the coupled physical residual and official error. Evolution, unitarity and cleanup
are correctly excluded by the brief.

**E. Nontrivial examples — PASS.** The trivial endpoint would be a single scalar or
commuting projector and a redundant factor rotation. The exact two-projector fixture has
commutator `{{0,12/25},{-12/25,0}}` (`GramCores.wl:64`), the same-core example changes
normalization from `1/2` to `1` (`:43–44,63`), and the singular Schur fixture rejects an
out-of-range vector (`:30–33,52–53`). Coupled multi-block removal and rectangular
variable-projection finite differences are tested (`test_study.py:14–42`). All 65 actual
blocks, eight failed candidate metrics and a rejected fixed-dictionary global bound provide
structure beyond a size-one demonstration.

**F. Symbolic-first derivations — FINDING (one).** The claimed general projector-gradient
identity is checked only after fixing both weight matrices to numbers
(`GramCores.wl:34–38,59`; `README.md:46–53,139`). In particular, the selected matrices give

```wl
{{1,2},{3,-1}} . Transpose[{{2,-1},{1,3}}]
(* {{0,7},{7,0}} *)
```

Thus the saved exact check never exercises the diagonal cross-weight terms; the similarly
weighted WL fixture in `run.wls:42–44` does not fill this symbolic gap. This is a coverage
finding, not evidence that the formula is incorrect: an independent derivation using fully
symbolic `u` of shape `{2,3}`, `u0` of shape `{3,3}`, `w` of shape `{2,2}` and `w0` of shape
`{3,2}`, differentiating the directly flattened projector residual, returned zero for every
component of `derivative - formula`.

Smallest meaningful fix: retain symbolic weight entries in `GramSymbolicChecks[]` instead
of the numeric assignments at line 34, localize their coefficient heads, and rerun the same
exact derivative assertion. Keep the rational noncommuting numeric fixture as a separate
reference. This adds no optimization or circuit work. The remaining core, compression and
rescaling derivations retain their physical parameters symbolically; full-size numerical
linear algebra is an appropriate final data check.

**G. Regimes — PASS, with the general-gradient coverage defect counted once under F.**
The explicit table is `README.md:137–143`. General symbolic core/Schur/compression/rescale
identities execute; exact noncommutation, normalization and SDP cases execute; the identity
rotation and vanishing-residual limits execute (`GramCores.wl:65–66`). As an independent
zero-scale check, set `o={0,0,1}`, `v={2,3,0}`, `Q=I-o o^T`,
`Pv=v v^T/(v.v)`, and retain a symbolic two-row `w`: at scale zero,
`w.Pv.Q.Transpose[Pv].Transpose[w]-w.Pv.Transpose[w]`,
`w.Pv.Q.v-w.v`, and `v.Q.v-v.v` all expand to zero. All-data numeric references and
resource checks reran successfully. Invalid WL shapes, out-of-range Schur vectors, zero
removal direction and nonnegative-fit boundary execute; an independent call to
`weights_for` with duplicate projector rows raised `ValueError: Singular or unresolved
dictionary Gram.` The saved metric and rank screens include actual rejection boundaries.
No additional missing physical regime was found.

OPEN ISSUES: 1
