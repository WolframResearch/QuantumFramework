# References

Source material for the QEC project. Each PDF has a plain-text twin extracted with
`pdftotext -layout`: **grep the `.txt`, read the `.pdf`** — the text version is far
cheaper to search and keeps section numbering intact.

## The files are deliberately not in the repository

Only this README is tracked; `.gitignore` excludes everything else here. Gottesman's 2026
book is an unpublished draft from the author's own page and is not ours to redistribute,
and none of the material needs to live in version control for the layer to build or its
tests to pass — nothing in `QECCore/`, `Tests/` or `docs/` reads these files. What is
worth keeping under version control is the index below, which records exactly which
sections each piece of the implementation rests on.

To populate the folder, download from the links in the table and extract the text twins:

```
pdftotext -layout Gottesman-1997-Thesis.pdf Gottesman-1997-Thesis.txt
pdftotext -layout Gottesman-2026-Book.pdf   Gottesman-2026-Book.txt
```

## In this folder

| File | What it is |
|---|---|
| `Gottesman-1997-Thesis.pdf` / `.txt` | Daniel Gottesman, *Stabilizer Codes and Quantum Error Correction*, Caltech PhD thesis, [arXiv:quant-ph/9705052](https://arxiv.org/abs/quant-ph/9705052) |
| `Gottesman-2026-Book.pdf` / `.txt` | Daniel Gottesman, *Surviving as a Quantum Computer in a Classical World*, 2026 draft, [source](https://www.cs.umd.edu/~dgottesm/QECCbook-2026.pdf). Cite as the 2026 draft |
| `Bahrami-StabilizerFormalism.nb` | Mohammad Bahrami, *Stabilizer formalism via computation*, [Wolfram Community 3738980](https://community.wolfram.com/groups/-/m/t/3738980) — the framework's own stabilizer tour |

## Sections actually used

**Thesis** — §3.2 general stabilizer codes (normalizer, degeneracy, syndromes,
logical operators as N(S)/S); §3.4 the binary symplectic and GF(4) languages;
§3.5 building codes from old codes (qubit removal, pasting, concatenation);
§4.1 standard form and the closed formulas for logical operators; §4.2 encoding
networks; §8.1 distance-two codes; §8.4 perfect codes; §8.6 CSS codes.

**Book** — ch. 4 classical linear codes (prerequisite for CSS); §5.1 CSS codes;
§5.2 GF(4); §6.3 Clifford generators and **Table 6.1**, the gate/row-operation
dictionary; §6.4 encoding circuits (**Procedure 6.6**, the approach the encoder
here follows — it handles non-CSS codes where the thesis shortcut fails);
ch. 7 bounds; **ch. 8 qudits** (§8.1.1 prime dimension, §8.1.2 prime powers,
§8.2 qudit stabilizer codes, §8.3 qudit CSS, §8.4 qudit Clifford);
§9.1.2 concatenated stabilizer codes; **ch. 17 toric and surface codes**.

**Book, for the noisy-extraction half (roadmap item 2)** — §10.1.1 the basic error model
(Definitions 10.1-10.3: locations, faults, one rate per location type), which the
circuit-level noise model is a Pauli specialisation of; §10.1.2 error propagation;
**§10.2** the formal definition of fault tolerance, whose propagation properties are what
this extraction circuit fails; §10.4 the phenomenological model, and the warning that
thresholds are not comparable across noise levels; **§12.1.1 and figures 12.1a/12.1b**
the extraction circuit itself and its non-fault-tolerance; **§12.2.2 with figures 12.6
and 12.7** measuring the whole syndrome and then repeating, and how many repetitions a
guarantee needs; §12.1.2-12.1.3, §12.2, §12.3 (incl. §12.3.3), §12.4, §12.5 the verified-
ancilla protocols that fix it, and §12.5.1 for why surface-code practice uses the
non-fault-tolerant circuit anyway; §15.4 and §15.4.3 resetting and reusing ancillas;
**§15.5.1 and §15.5.2** parallelism, waiting, and why serial extraction makes idle noise
matter more rather than less.

**Not in Gottesman**, and cited in the source where used: detectors as differences of
repeated noisy syndromes are Dennis, Kitaev, Landahl and Preskill, *Topological quantum
memory*, [arXiv:quant-ph/0110143](https://arxiv.org/abs/quant-ph/0110143); the detector
error model as an explicit object, and the circuit syntax this layer exports, are Gidney,
*Stim: a fast stabilizer circuit simulator*,
[arXiv:2103.02202](https://arxiv.org/abs/2103.02202). §12.2.2 raises the differencing
strategy and sets it aside as hard to analyse in general, so the book is not its source.

Where both cover the same ground, the book is the better guide — its encoder
construction is uniform for CSS and non-CSS codes, and its notation is current.
Note that the book's appendices B (group theory), C (finite fields) and D (linear
algebra) are **empty placeholders** in the 2026 draft, and that its concatenation
section uses "inner"/"outer" in the opposite sense to the classical convention
used in this package.

## Not stored here (fetch if needed)

- Aaronson & Gottesman, *Improved Simulation of Stabilizer Circuits*,
  [arXiv:quant-ph/0406196](https://arxiv.org/abs/quant-ph/0406196) — the CHP tableau
  with destabilizers, the basis of the framework's engine.
- Gidney, *Stim: a fast stabilizer circuit simulator*,
  [arXiv:2103.02202](https://arxiv.org/abs/2103.02202) — detector error models,
  Pauli-frame sampling; the cross-check target.
- Higgott & Gidney, PyMatching / sparse blossom — [GitHub](https://github.com/oscarhiggott/PyMatching).
- Delfosse & Nickerson, union-find decoder,
  [arXiv:1709.06218](https://arxiv.org/abs/1709.06218) — the simplest decoder worth
  implementing after lookup.
- Dennis, Kitaev, Landahl & Preskill, *Topological quantum memory*,
  [arXiv:quant-ph/0110143](https://arxiv.org/abs/quant-ph/0110143); Fowler et al.,
  [arXiv:1208.0928](https://arxiv.org/abs/1208.0928) — surface codes.
- Bravyi et al., bivariate-bicycle qLDPC codes,
  [arXiv:2308.07915](https://arxiv.org/abs/2308.07915) — the [[144,12,12]] Gross code.
- Forney, *Introduction to Finite Fields*, MIT 6.451 ch. 7 —
  [PDF](https://ocw.mit.edu/courses/6-451-principles-of-digital-communication-ii-spring-2005/8921c013c1d7802bb197572c52442b97_chap7.pdf).
  Needed only for prime-power qudit dimensions (book §8.1.2); prime dimensions
  need nothing beyond arithmetic mod p.
- Gottesman's UMD course CMSC 858G (Spring 2026) —
  [problem sets](https://www.cs.umd.edu/class/spring2026/cmsc858G/), including qudit
  exercises usable as test cases.
