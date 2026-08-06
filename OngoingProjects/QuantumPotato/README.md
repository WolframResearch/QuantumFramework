# Quantum Potato Chip: folder guide

The object of the project: for a fixed tetrahedral SIC measurement, the qubit states whose four outcome probabilities factor into a product of two binary distributions form the surface $\sqrt{3}\,\langle\sigma_y\rangle=\langle\sigma_x\rangle\langle\sigma_z\rangle$ cut by the Bloch ball, a doubly ruled saddle (the "chip"); relabeling the outcomes gives two congruent siblings. The folder holds two generations of writing about it: the 2026-08-05 postmortem of the 2024 submission ([arXiv:2411.01082](https://arxiv.org/abs/2411.01082), "Quantum Potato Chips"), and the rebuilt concept ("a qubit state that is exactly a pair of independent coins") developed 2026-08-05/06.

**Two sign conventions coexist.** Files that pin the paper's tetrahedron ($\mathbf n_1=\tfrac{1}{\sqrt3}(-1,1,-1)$, ...) get row coin $p=\tfrac12(1-z/\sqrt3)$; files that index effects by sign patterns ($\vec n_{ij}=\tfrac{1}{\sqrt3}((-1)^j,(-1)^{i+j},(-1)^i)$) get $p=\tfrac12(1+z/\sqrt3)$. Formulas moved between the two lines differ by signs, not by errors.

## The documents

| File | Role | Convention |
|---|---|---|
| [quantum-potato-chip-tex-audit.md](quantum-potato-chip-tex-audit.md) | Audit of record for the submission: transcription errors E1-E11, overclaims C1-C5, venue diagnosis | paper |
| [quantum-potato-chip-core-concept-review.md](quantum-potato-chip-core-concept-review.md) | Judgment of the paper's derivation (D1-D7) and exposition (E1-E7), with the ten-line target derivation | paper |
| [quantum-potato-chip-paper-question-chain.md](quantum-potato-chip-paper-question-chain.md) | The paper compiled into a typed question chain (D1-D12, P1-P6, M1-M3; 4 audit rounds) | paper |
| [quantum-potato-chip-paper-chain-answers.md](quantum-potato-chip-paper-chain-answers.md) | Closed-form answer sheet to that chain: P1, the P2 reconstruction, $\phi$, chip-preserving noise | paper |
| [quantum-potato-chip-prl-derivation.md](quantum-potato-chip-prl-derivation.md) | Compact 8-step derivation, independent sibling from a parallel session; the verified current piece is derivation-pra | paper |
| [quantum-potato-chip-concept-draft.md](quantum-potato-chip-concept-draft.md) | Rev 6 of the rebuilt core: six standalone questions Q1-Q6, with draft brief and inquiry frame | coins |
| [quantum-potato-chip-question-chain.md](quantum-potato-chip-question-chain.md) | The coins concept as a formal chain with Readings (D1-D11, P1-P4, M1-M2; two cold reviews) | coins |
| [quantum-potato-chip-question-chain.pdf](quantum-potato-chip-question-chain.pdf) | 2026-08-05 render of the chain above; regenerate after editing the md | coins |
| [quantum-potato-chip-derivation-pra.md](quantum-potato-chip-derivation-pra.md) | Verified article-form derivation: the two coins, the independence locus, three chips, rulings, rim | coins |
| [question-chain-method-prompt.md](question-chain-method-prompt.md) | Portable question-chain method (no chip content) | n/a |

## Sources

- **`QuantumPotato.nb` / `QuantumPotato V2.nb`**: the computations behind the paper; V2 is the notebook published at [wolfr.am/QPC](https://wolfr.am/QPC), while V1 additionally holds the introduction's simplex figures and a closing uncertainty-analysis section absent from the submission. The notebook results are correct; most printed equation errors entered during hand-transcription into TeX (tex-audit, headline finding).
- **`quantum_potato_chip.zip`**: the submission source (`main.tex` + 15 figure PNGs) that the audits reference.
- **`related-papers/`**: TeX sources of the eight papers mined for prior art (0910.2750 QBist state spaces; 1612.03234 qplex; 1803.09339 and 1902.03613 probability representation; 2306.00086 KD geometry; 2403.18899 KD review; quant-ph/0405070 and quant-ph/0506222 discrete-Wigner classicality). `q_T1..T12.xml` and `q_V1..V2.xml` are 2026-08-05 arXiv API snapshots of the prior-art queries; `search_arxiv.py` ran the T-battery; `parse_atom.py` pretty-prints a saved feed.

## Status

- Concept-draft revisions 4-6 have not had a fresh adversarial round; short proof write-ups for Q2's covariant-core and emptiness statements are due.
- The Q3/P3 characterization ("are the chip-preserving channels exactly classical noises of the coins?") is open.
- Each chain's final repair was rule-checked locally but not re-audited (round caps); neither claims an independently certified clean pass.
