# Stabilizer Papers — arXiv Fetch Report

**Source list:** `paper-bibliography.md` (33 papers)  
**Output folder:** `/Users/mohammadb/Documents/GitHub/Planning for research/Stabilizer`  
**Generated:** 2026-04-30 13:20:47

## Summary

- **Total papers in source list:** 33
- **arXiv ID resolved:** 29/33 (unique IDs: 28; #14 is a duplicate of #13)
- **TeX archives downloaded:** 27 unique archives covering 29 papers
- **Papers not on arXiv:** 4 (industry conferences, Springer chapters, lecture notes — see bottom of report)

### Resolution method breakdown

| Method | Count | Description |
|---|---:|---|
| `S2-doi` | 21 | Semantic Scholar API lookup by DOI |
| `S2-ssId` | 2 | Semantic Scholar API lookup by paper ID |
| `S2-retry` | 1 | S2 lookup that initially returned HTTP 429, succeeded on retry |
| `arxiv-search` | 2 | arXiv API search by author + title (after S2 returned no_arxiv) |
| `known` | 2 | arXiv ID supplied directly (Gottesman thesis is widely cited; PauliEngine has arXiv DOI) |
| `duplicate` | 1 | Same arXiv preprint as another entry |
| `not-on-arxiv` | 4 | Confirmed no arXiv version (S2 + arXiv title search both empty) |

## Per-paper status

| # | Tag | arXiv ID | Method | Archive | Title |
|---:|---|---|---|---|---|
| 1 | `Gid21_Stim` | `2103.02202` | S2-doi | ✅ 659,075 bytes | Stim: a fast stabilizer circuit simulator |
| 2 | `AarGot04` | `quant-ph/0406196` | S2-doi | ✅ 51,341 bytes | Improved Simulation of Stabilizer Circuits |
| 3 | `Got98_Heisenberg` | `quant-ph/9807006` | S2-retry | ✅ 15,849 bytes | The Heisenberg Representation of Quantum Computers |
| 4 | `DehMoo03` | `quant-ph/0304125` | S2-doi | ✅ 15,628 bytes | Clifford group, stabilizer states, and linear and quadratic operations over GF(2) |
| 5 | `AndBri05` | `quant-ph/0504117` | S2-doi | ✅ 25,827 bytes | Fast simulation of stabilizer circuits using a graph-state representation |
| 6 | `EllEasCav08` | `0806.2651` | S2-doi | ✅ 98,752 bytes | Graphical description of Pauli measurements on stabilizer states |
| 7 | `KocHuaLov17` | `1703.04630` | S2-doi | ✅ 43,371 bytes | Discrete Wigner Function Derivation of the Aaronson-Gottesman Tableau Algorithm |
| 8 | `Reid24` | `2404.19408` | S2-doi | ✅ 66,698 bytes | A Simple Method for Compiling Quantum Stabilizer Circuits |
| 9 | `deSilSalYin23` | `2311.10357` | S2-doi | ✅ 72,814 bytes | Fast algorithms for classical specifications of stabiliser states and Clifford gates |
| 10 | `RuhDev25` | `2405.03970` | arxiv-search | ✅ 1,055,755 bytes | Quantum Circuit Optimization and MBQC Scheduling With a Pauli Tracking Library |
| 11 | `Paler14` | `1401.5872` | arxiv-search | ✅ 71,129 bytes | Software-based Pauli tracking in fault-tolerant quantum circuits (arXiv: 'Software Pauli Tracking for Quantum Computation') |
| 12 | `Riesebos17` | `—` | not-on-arxiv | — | Pauli frames for quantum computer architectures |
| 13 | `PatGuh26` | `2312.02377` | S2-doi | ✅ 3,897,032 bytes | A Graphical Rule Book for Clifford Manipulations of Stabilizer States |
| 14 | `PatGuh23` | `2312.02377` | duplicate | — | Clifford Manipulations of Stabilizer States: A graphical rule book... (preprint of #13) |
| 15 | `Yashin25` | `2504.14101` | S2-doi | ✅ 40,653 bytes | A Streamlined Demonstration That Stabilizer Circuits Simulation Reduces to Boolean Linear Algebra |
| 16 | `Winderl23` | `2309.08972` | S2-ssId | ✅ 11,882,419 bytes | Architecture-Aware Synthesis of Stabilizer Circuits from Clifford Tableaus |
| 17 | `KoeSmo14` | `1406.2170` | S2-doi | ✅ 14,415 bytes | How to efficiently select an arbitrary Clifford group element |
| 18 | `HosDehMoo04` | `quant-ph/0408190` | S2-doi | ✅ 17,322 bytes | Stabilizer states and Clifford operations for systems of arbitrary dimensions and modular arithmetic |
| 19 | `Beaudrap11` | `1102.3354` | S2-doi | ✅ 58,420 bytes | A linearized stabilizer formalism for systems of finite dimension |
| 20 | `Got97_Thesis` | `quant-ph/9705052` | known | ✅ 97,560 bytes | Stabilizer Codes and Quantum Error Correction (Gottesman PhD thesis) |
| 21 | `Got00_Intro` | `quant-ph/0004072` | S2-doi | ✅ 15,836 bytes | An Introduction to Quantum Error Correction |
| 22 | `Got09_Intro` | `0904.2557` | S2-doi | ✅ 51,165 bytes | An Introduction to Quantum Error Correction and Fault-Tolerant Quantum Computation |
| 23 | `MonPar23` | `2309.11793` | S2-doi | ✅ 739,585 bytes | Quantum Circuits for Stabilizer Error Correcting Codes: A Tutorial |
| 24 | `Biswas24` | `2405.13590` | S2-doi | ✅ 81,857 bytes | On Classical Simulation of Quantum Circuits Composed of Clifford Gates |
| 25 | `Fujii15` | `—` | not-on-arxiv | — | Stabilizer Formalism and Its Applications |
| 26 | `Browne14` | `—` | not-on-arxiv | — | Quantum Error Correction and The Stabilizer Formalism (Revision) |
| 27 | `Arab24` | `—` | not-on-arxiv | — | Lecture notes on quantum entanglement: From stabilizer states to stabilizer channels |
| 28 | `PayWin24a` | `2407.16801` | S2-doi | ✅ 132,328 bytes | Qudit Quantum Programming with Projective Cliffords |
| 29 | `WinPay24b` | `2407.16861` | S2-doi | ✅ 43,457 bytes | Condensed encodings of projective Clifford operations in arbitrary dimension |
| 30 | `Mueller26_PauliEngine` | `2601.02233` | known | ✅ 81,583 bytes | PauliEngine: High-Performant Symbolic Arithmetic for Quantum Operations |
| 31 | `FangYing23_SymPhase` | `2311.03906` | S2-doi | ✅ 1,090,787 bytes | SymPhase: Phase Symbolization for Fast Simulation of Stabilizer Circuits |
| 32 | `GarMar15_Frames` | `1712.03554` | S2-doi | ✅ 838,465 bytes | Simulation of Quantum Circuits via Stabilizer Frames |
| 33 | `GarMarCro12_InnerProduct` | `1210.6646` | S2-ssId | ✅ 174,349 bytes | Efficient Inner-product Algorithm for Stabilizer States |

## Folder layout

```
Stabilizer/
├── eprint/                            # raw .tar.gz e-print archives
│   └── <Tag>_<safe-arxiv-id>.tar.gz
├── tex/                               # extracted LaTeX projects
│   └── <Tag>_<safe-arxiv-id>/
│       ├── <main>.tex
│       └── ...
└── report.md                          # this file
```

Old-style arXiv IDs use `_` instead of `/` in filenames (e.g. `quant-ph_9705052`).

## Papers with no arXiv version

### #12  `Riesebos17` — Pauli frames for quantum computer architectures
- **Year:** 2017  
- **Venue:** DAC  
- **Status:** S2 returned `no_arxiv`; arXiv title-search returned no matching paper.  
  These are typically industry-conference papers, Springer book chapters, or local lecture notes that bypass arXiv.

### #25  `Fujii15` — Stabilizer Formalism and Its Applications
- **Year:** 2015  
- **Venue:** Springer book chapter  
- **Status:** S2 returned `no_arxiv`; arXiv title-search returned no matching paper.  
  These are typically industry-conference papers, Springer book chapters, or local lecture notes that bypass arXiv.

### #26  `Browne14` — Quantum Error Correction and The Stabilizer Formalism (Revision)
- **Year:** 2014  
- **Venue:** lecture notes  
- **Status:** S2 returned `no_arxiv`; arXiv title-search returned no matching paper.  
  These are typically industry-conference papers, Springer book chapters, or local lecture notes that bypass arXiv.

### #27  `Arab24` — Lecture notes on quantum entanglement: From stabilizer states to stabilizer channels
- **Year:** 2024  
- **Venue:** Frontiers of Physics  
- **Status:** S2 returned `no_arxiv`; arXiv title-search returned no matching paper.  
  These are typically industry-conference papers, Springer book chapters, or local lecture notes that bypass arXiv.

## Notes on duplicates and overlaps

**#13 PatGuh26** and **#14 PatGuh23** map to the same arXiv preprint, `2312.02377`. The Undermind catalog separated them because Semantic Scholar tracked the December-2023 arXiv preprint and the 2026 IEEE TQE journal publication as distinct records, but only one TeX source exists. The archive is saved under the `PatGuh26` tag.

**#11 Paler14** notice: the conference paper title is *“Software-based Pauli tracking in fault-tolerant quantum circuits”* (DATE 2014). The arXiv preprint title is *“Software Pauli Tracking for Quantum Computation”* (arXiv:1401.5872, Jan 2014). Same authors (Paler, Devitt, Nemoto, Polian) and same content; the journal title was reworked for the conference.

**#3 Got98_Heisenberg** notice: there is no formal venue for this paper — it is a 1998 talk-to-arXiv submission `quant-ph/9807006`. The Semantic Scholar record had no DOI; resolved via the Semantic Scholar paper ID, which initially hit HTTP 429 and succeeded on retry.

## How this list was produced

1. Read source markdown, extracted 33 entries (title + DOI or Semantic Scholar paper ID).
2. **Phase 2 — Resolve arXiv IDs:** queried `api.semanticscholar.org/graph/v1/paper/{DOI:...|paperId}?fields=externalIds`; pulled `externalIds.ArXiv` when present. Rate-limited to 1.6s per call.
3. **Phase 2b — Recovery:**
   - Retried HTTP 429s with longer backoff (one resolved: Got98_Heisenberg).
   - For the 8 still-unresolved papers, ran an arXiv API search (`export.arxiv.org/api/query`) using `au:<author> AND ti:<keywords>` queries.
   - Verified each candidate by fetching the abstract page and reading the `<meta name="citation_title">` tag — only accepted matches where the metadata title clearly matched the source title (this caught 3 false positives from a naive set-intersection similarity score above 0.5 that were actually unrelated papers).
   - For Gottesman's PhD thesis, supplied the well-known arXiv ID `quant-ph/9705052` directly (the source DOI `10.7907/RZR7-DT72` is a Caltech thesis ID, not in S2).
4. **Phase 3 — Download:** `curl -sSL https://arxiv.org/e-print/<id>` per paper, saved to `eprint/`. Extracted with `tar -xzf` (most papers) or `gunzip` for single-file submissions (Dehaene-De Moor 2003, Hostens 2004, Gottesman 2000, Gottesman 2009, Gottesman 1998).

All HTTP traffic ran through `/usr/bin/curl`; no Python TLS dependencies were required.