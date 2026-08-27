# PauliStabilizer named-codes audit

A Markdown reference on named stabilizer codes for the Wolfram
`QuantumFramework` `PauliStabilizer` kernel (paclet 2.1.1).

## What this is

`main.md` is a self-contained reference document with four parts:

1. **How stabilizer codes are named in the literature**: the `[[n,k,d]]`
   convention, a table of folklore names and their aliases (Shor, Steane,
   perfect/Laflamme, bit-flip/phase-flip, C4/iceberg, Bacon-Shor, surface/toric,
   color, quantum Reed-Muller, quantum Hamming, Golay, LDPC families), the
   distinction between a code (a subspace) and the codeword *state* that
   `PauliStabilizer` actually stores, and the exact generator convention used in
   `Kernel/Stabilizer/NamedCodes.m` (signed Pauli strings, qubit 1 leftmost,
   parity checks then logical Z-bar, `|0_L>` base with a sign-flipped `|1_L>`
   variant).

2. **A proposed catalog** of codes worth adding, in priority order, each with
   name(s) and aliases, `[[n,k,d]]`, explicit stabilizer generators as Pauli
   strings, a one-line description, and a primary reference. Every generator set
   was checked against the live kernel as a valid stabilizer state; distances
   are cited, not recomputed. Codes whose standard presentation is not a single
   unique generator set (rotated surface code layout, Golay basis, k > 1 logical
   choice) are flagged.

3. **A contributor checklist** for adding a named code: register the name(s) in
   `$PauliStabilizerNames`, add a `PauliStabilizer["name"] := PauliStabilizer[{
   generators }]` rule, comment with the physics, add tests under
   `Tests/Stabilizer/`, and update the doc page.

4. **References**, with full bibliographic data and arXiv links, in the
   document's final section.

The document proposes documentation and physics; it changes no kernel code. The
Pauli-string generators live in the Markdown as documentation only.

## Files

- `main.md`: the reference document, bibliography included.
- `README.md`: this file.
