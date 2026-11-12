# arXiv 2512.13908 (Google Willow, magic state cultivation) and arXiv 2404.18828 (HQS on IBM-Q, system-bath physics)

Sources:
- `../tex/2512.13908/main.tex` (501 lines) and `supplemental_material.tex` (483 lines)
- `../tex/2404.18828/demonstrationBathMapping_resubmission.tex` (1583 lines)

All three files were read in full. Neither paper ships .csv or .txt data.

2512.13908 has no table in its TeX. Its per-qubit values are in vector figure PDFs (`figures_for_SI/qubit_specs_105.pdf`, `2q_specs_105.pdf`, `ro_specs_105.pdf`), extracted with `pdftotext -bbox`. Statistics computed from those figures are marked **[fig-derived]**; they are not prose statements of the paper.

## arXiv 2512.13908: Google Willow

**Device.** "a Willow processor ..., which contains 105 tunable transmon qubits in a square lattice, coupled by tunable transmon couplers" (supp:38).
- Specs cover the 54 qubits of Fig. 4: "performance numbers are specific to the indicated operating frequencies for Fig. 4" (supp:15).
- No date is given.

**1Q gate duration (verbatim).** "The single-qubit gates all had duration 35 ns." (Fig. S1 caption, supp:15).
- CZ gates are the same length: "The CZ gates all had duration 35 ns." (supp:25).

**Pulse shape, DRAG, drive strength.** Not stated.

**Qubit frequency** [fig-derived, Fig. S1a/c, n=54]:
- Operating frequency 6.05 to 6.74 GHz, median 6.46 GHz.
- Maximum frequency 6.60 to 6.89 GHz, median 6.725 GHz.

**Anharmonicity** [fig-derived, Fig. S1b, "Qubit f10 - f21 [MHz]"]: 189 to 202 MHz, median 196 MHz. That is α/2π ≈ −196 MHz.

**T1** [fig-derived, Fig. S1d]: 27 to 112 µs, median 89 µs.

**T2.** "The coherence T_2 is from a CPMG experiment" (supp:15). [fig-derived, Fig. S1e]: 19 to 113 µs, median 64.5 µs.

**T_phi.** Not stated. The model uses a quasi-static detuning drawn uniformly from −200 to 200 kHz (supp:136, 154).

**1Q gate error.** "Single-qubit gate error is from Clifford randomized benchmarking (including both π and π/2 rotations), and we report 'average' error." (supp:15).
- [fig-derived, Fig. S1f]: 2.5 to 49 ×10^-4, median 6.25×10^-4, mean 9.7×10^-4.

**Leakage.** Qualitative only:
- Leakage is likely why the logical error was "about seven times higher than expected" (supp:459).
- "Our circuits did not contain the leakage removal units" (supp:463).
- Possible cause: "frequency collisions involving qubit |2⟩ states and subsequent leakage" (supp:463).

**What limits gate speed.** Not stated.

## arXiv 2404.18828: HQS on IBM-Q

**Devices.** IBM-Q Nairobi (Fig. 4, Fig. 6) and IBM-Q Lagos (Fig. 5).
- No processor family or date is given.
- The calibration data was "obtained at the time the simulations were run" (line 1146).

**1Q gate duration (verbatim).** "All single-qubit gate times were 35.5 ns, except the (virtual) rotate-Z gate whose gate time is zero." (Table I caption, line 1148).
- Identity: "35 ns qubit-sleep" (lines 1199, 1236).
- Gate set: "√X, X, CNOT, virtual R_z, and identity" (line 1442).

**Pulse shape, DRAG, drive strength.** Not stated.

**Frequency and anharmonicity.** "had commonly frequencies close to 5 GHz with anharmonicities of the order -0.3 GHz" (line 1120).

**T1 and T2.** Table I (lines 1124-1140), 13 entries:
- T1: 54.5 to 126 µs, median 101 µs (median computed by the reading agent).
- T2: 21.5 to 220 µs, median 84.1 µs. The type of T2 measurement is not stated.
- The text gives "T_1 ~ 100 µs" (1194) and "T_2 ~ 80 µs" (1197).
- Drift after calibration is noted (673).

**1Q gate error.** Table I: X or √X error from 1.43 to 4.01 ×10^-4, median 2.79×10^-4 (median computed by the reading agent).

**Coherent errors (measured).**
- √X overrotation: δφ = −0.01 and 0.003 rad (line 1201).
- Idle Z error: −0.03 and 0.01 rad per 35 ns identity (lines 1199, 1236).
- Overrotations "could be due to the influence of higher levels of the qubits" (1258-1259).

**CNOT.** Gate time 249 to 960 ns; error 0.48 to 2.26 ×10^-2 (Table II).

## Summary table

| Parameter | Willow (2512.13908) | IBM-Q (2404.18828) |
|---|---|---|
| 1Q gate duration | 35 ns (verbatim) | 35.5 ns (verbatim); virtual R_z 0 ns |
| α/2π | −196 MHz median, \|α\| from 189 to 202 MHz [fig-derived] | "of the order −0.3 GHz" |
| ω_q/2π | 6.46 GHz median, from 6.05 to 6.74 [fig-derived] | "close to 5 GHz" |
| T1 | 89 µs median, from 27 to 112 [fig-derived] | 101 µs median, from 54.5 to 126; "~100 µs" |
| T2 | 64.5 µs median (CPMG), from 19 to 113 [fig-derived] | 84.1 µs median, from 21.5 to 220; "~80 µs" |
| 1Q error | RB average 6.25e-4 median, from 2.5e-4 to 4.9e-3 [fig-derived] | 2.79e-4 median, from 1.43e-4 to 4.01e-4 |
| Pulse shape / DRAG / drive | not stated | not stated |
| Leakage | qualitative only | qualitative (higher levels) |
