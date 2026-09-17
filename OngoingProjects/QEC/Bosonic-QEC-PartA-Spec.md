# Part (a) spec: the bosonic code object and symbolic Knill-Laflamme

*One-page implementation spec for part (a) of `Bosonic-QEC-Plan.md`. Concrete enough to
start coding from. The design keeps a single public symbol with a property surface, so the
API stays small by construction.*

## The one public symbol

`QECBosonicCode` is the object. Everything is a property on it; there are no free
functions to sprawl. It is the bosonic sibling of the qubit `QECCode`.

### Constructors

```
QECBosonicCode[{w0, w1}]                    (* from two codewords *)
QECBosonicCode["Binomial", N, S]            (* named family, parameters (N, S) *)
QECBosonicCode["Cat", legs, \[Alpha]]       (* 2- or 4-component cat, amplitude \[Alpha] *)
QECBosonicCode[qs_QuantumState, ...]        (* codewords as QuantumFramework states *)
```

A codeword is given in the Fock basis as an association `<|n -> amplitude, ...|>` of
occupation number to (possibly symbolic) amplitude, as a `CoherentState` combination, or
as a `QuantumState`. The named constructors build the association. The object stores the
codewords, the mode count, and the assumptions on any symbolic parameters.

### Direct properties

```
code["Codewords"]           (* the two logical codewords *)
code["Modes"]               (* number of bosonic modes (1 for single-mode) *)
code["MeanPhotonNumber"]    (* {<W0|n|W0>, <W1|n|W1>}, symbolic *)
code["Parameters"]          (* e.g. <|"N" -> 1, "S" -> 1|> for a named code *)
code["Properties"]          (* the list *)
```

### The Knill-Laflamme surface (the payoff)

An error channel is named as `"Loss"[L]` (error set `{I, a, ..., a^L}`), `"Dephasing"[D]`
(`{I, n, ..., n^D}`), or an explicit list of bosonic operators built from QuantumFramework
primitives.

```
code["KnillLaflammeMatrix", channel]   (* the Hermitian coefficient matrix h_ab, symbolic *)
code["CorrectableQ", channel]          (* exact True/False, by symbolic identity *)
code["CorrectionOrder", "Loss"]        (* largest L with exact correction *)
code["Signature", channel]             (* the local-unitary invariant read off h (Du et al.) *)
```

`"CorrectableQ"` returns `True` only when the overlap matrix `<W_i|E_a\[Dagger] E_b|W_j>`
equals `h_ab \[Delta]_ij` as a symbolic identity under the code's assumptions, decided with
`PossibleZeroQ` on each residual, never a numerical tolerance.

## Worked example (the first milestone)

The smallest binomial code, parameters `N = 1`, `S = 1`, codewords
`(|0> + |4>)/Sqrt[2]` and `|2>`:

```
code = QECBosonicCode["Binomial", 1, 1];

code["Codewords"]              (* {<|0 -> 1/Sqrt[2], 4 -> 1/Sqrt[2]|>, <|2 -> 1|>} *)
code["MeanPhotonNumber"]       (* {2, 2}  -- equal, which is why loss is correctable *)

code["KnillLaflammeMatrix", "Loss"[1]]
(* {{1, 0}, {0, 2}}  -- exact:
   h00 = <W|I|W> = 1;  h11 = <W|a\[Dagger]a|W> = <W|n|W> = 2, same for both codewords;
   h01 = <W|a|W> = 0 because a flips photon-number parity (even codewords -> odd image). *)

code["CorrectableQ", "Loss"[1]]   (* True  *)
code["CorrectableQ", "Loss"[2]]   (* False -- one binomial mode does not correct 2 losses *)
code["CorrectionOrder", "Loss"]   (* 1 *)
```

Then the same three checks on `QECBosonicCode["Cat", 4, \[Alpha]]` (the four-component cat,
`\[Alpha]` symbolic), which corrects single loss exactly in the ideal limit, and on
`QECBosonicCode["Cat", 2, \[Alpha]]`, which corrects it only approximately, so
`"CorrectionOrder"` returns the leading order in `\[Alpha]` rather than an exact integer.

## How the matrix element is computed (reuse, do not reinvent)

`<W_i| E_a\[Dagger] E_b |W_j>` is assembled by the `SecondQuantization` machinery Bruno
already owns:

1. Form `E_a\[Dagger] E_b` (a product of `a`, `a\[Dagger]` powers) and normal-order it with
   `BosonicNormalOrder` into a sum of `a\[Dagger]^m a^n`.
2. Evaluate each `<W_i| a\[Dagger]^m a^n |W_j>` in closed form: for Fock-sum codewords this
   is a finite symbolic sum over the shared occupations (each term a ratio of factorials);
   for coherent-state (cat) codewords it is
   `\[Alpha]*^m \[Beta]^n <\[Alpha]|\[Beta]>`, a Gaussian exponential.
3. Reduce the resulting matrix with `FullSimplify` under the parameter assumptions;
   subtract `h_ab \[Delta]_ij` and decide each entry with `PossibleZeroQ`.

No Fock truncation and no numerics enter: the sums are finite for Fock codewords and
closed-form for coherent-state codewords. `HypergeometricPFQ` and the theta functions close
the cat and binomial overlaps.

## How it interacts with QuantumFramework

- Codewords are `QuantumState` objects (or convert to them), so `code["Codewords"]` composes
  with the rest of QuantumFramework and with the phase-space tools (`WignerFunction`,
  `HusimiQFunction`).
- Error operators come from the QuantumFramework bosonic layer (`AnnihilationOperator` and
  its powers), so an explicit error set is native QuantumFramework, not strings.
- The object is forward-compatible with the qubit `QECCode`: both answer "is this a code,
  and does it correct this error set", one over the Pauli group and GF(2), the other over
  the oscillator and its Fock or coherent-state overlaps.

## Scope of part (a)

One public symbol, roughly four direct properties and four Knill-Laflamme properties, the
binomial and cat milestones, and the reuse of `BosonicNormalOrder` and the coherent-state
overlaps. Code discovery (part b) and channel-adapted recovery (part c) are separate and
follow.
