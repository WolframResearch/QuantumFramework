# Transmon single-qubit gate parameters: what the cited papers state

Compiled 2026-10-12 for `TransmonGates.md` and `TransmonGatesV2.md`. Every TeX source below was read in full; nothing was sampled.

Layout:
- `tex/<arXiv id>/`: the TeX sources and the original `source.tar.gz` of each paper.
- `readings/`: one file per reading, with verbatim quotes and file:line anchors.

| arXiv | Paper | Reading |
|---|---|---|
| cond-mat/0703002 | Koch et al., PRA 76, 042319 (2007) | readings/Koch2007-Krantz2019.md |
| 1904.06560 | Krantz et al., Appl. Phys. Rev. 6, 021318 (2019) | readings/Koch2007-Krantz2019.md |
| 0901.0534 | Motzoi et al., PRL 103, 110501 (2009) | readings/DRAG-…md |
| 1011.1949 | Gambetta et al., PRA 83, 012308 (2011) | readings/DRAG-…md |
| 2511.12799 | Malarchick (2025), simulation-only calibration guidance | readings/DRAG-…md |
| 2512.19919 | Da Costa Jesus et al. (2025), DRAG under strong driving | readings/DRAG-…md |
| 1910.11333 | Google Sycamore, Supplementary Information (2019) | readings/1910.11333-Sycamore-SI.md |
| 2512.13908 | Google, magic state cultivation on Willow (2025) | readings/Willow2512.13908-IBM2404.18828.md |
| 2404.18828 | HQS, system-bath physics on IBM-Q (2024/25) | readings/Willow2512.13908-IBM2404.18828.md |

Provenance labels: **meas** = measured on a device; **sim** = simulation; **design** = design target or "typical" statement; **est** = theoretical estimate.

## Ranges per parameter

| Parameter | What the papers state | TransmonGates value |
|---|---|---|
| **1Q gate length** | **meas:** 35 ns (Willow, all 1Q gates), 35.5 ns (IBM-Q Nairobi and Lagos, all 1Q gates). **design/sim:** DRAG suffices for ≳ 20 ns when T2/T ≳ 500 (Malarchick 2025, sim); 1-10 ns "usually" and 10 ns "typical" (Koch 2007); "several tens of ns" reasonable (Koch). **research, fast gates:** 6.8 ns at 99.98% (Gao et al., cited in J25); simulated minimum ≈ 4.5 ns. **Not stated:** Sycamore's 1Q gate length. Its SI gives only the design example "25 μs coherence implies a 25 ns gate", and the Nature main text was not retrieved. | 8 ns (Leakage section and 8 ns calibrations); 10.2 ns (best DRAG pulse); 22.4 ns (best without DRAG) |
| **Anharmonicity α/2π** | **meas:** −208 MHz median (Sycamore, 53 qubits, stdev 4.7); −196 MHz median, \|α\| from 189 to 202 MHz (Willow, read from figures); "of the order −0.3 GHz" (IBM-Q). **design:** −100 to −300 MHz designed, −200 to −300 MHz typical (Krantz). **sim:** −200 (Krantz Fig. 14, Malarchick), −225 (J25), −400 MHz (M09) | −200 MHz |
| **Qubit frequency ω_q/2π** | **meas:** 6.66 GHz idle median (Sycamore); 6.46 GHz median (Willow); "close to 5 GHz" (IBM-Q). **design:** 3-6 GHz, ≈ 5 GHz typical (Krantz); 1-10 GHz (Koch 2007) | not used (rotating frame) |
| **T1** | **meas:** 15.5 µs median (Sycamore 2019); 85 µs (Krantz, citing Kjaergaard 2017); 89 µs median, 27 to 112 (Willow); ~100 µs, median 101 (IBM-Q). **trend:** "above 100 µs" (Krantz). **sim:** 37 µs (Garnet-like), 40/350/500 µs "achievable today" (J25) | 100 µs |
| **T2** | **meas:** 64.5 µs median CPMG, 19 to 113 (Willow); ~80 µs, median 84 (IBM-Q); T2* 95 or 98 µs and T2E 120 µs (Krantz); not stated for Sycamore. **bound:** T2 ≤ 2T1 (Krantz) | 100 µs, at the upper end of measured medians |
| **Pure dephasing** | Γ2 = Γ1/2 + Γφ (Krantz). Malarchick and J25 put the rate 1/T2 − 1/(2T1) on D[n̂], which damps the 0-1 coherence at half that rate (reading-agent note). | rate 2/T_φ on b†b, which gives coherence decay at 1/T2 exactly (kernel-checked) |
| **DRAG coefficient** (our λ, where the Q quadrature is λΩ̇/α) | **theory:** λ = 1 cancels leakage to first order (M09, Krantz); λ = 0.5 for phase error (Krantz, citing G11); first-order family 0, 0.5, 0.707, 1 (G11), with M09's λ = 1 best at second order. **sim:** optimized 1.0 to 1.3, analytic ≈ 1.2 (J25: Hann pulse, −225 MHz, 7-15 ns). **meas, cited:** δ_c = 0 at λ ≈ 0.5 (PRL 116, 020501, via J25); "in practice ... deviation from these two optimal values" (Krantz) | 1.41 at 8 ns; 1.49 at 10.2 ns |
| **Detuning** | Time-dependent δ1 = (λ1² − 4)Ω²/(4Δ) (M09). A constant δ0 is as good (G11 conjecture; J25: "constant detuning can outperform"). Constant-δ formula for a Hann π pulse: δ_c ≈ 0.712(λ2² − 4α)π²/(Δ2T²) (J25), with δ_c ∝ 1/T². Calibrated in experiments (Sycamore "drive detuning"; Krantz δf); no numbers given. | +0.24 rad/ns at 8 ns; +0.175 at 10.2 ns. Same sign as J25's formula. |
| **Pulse shape** | Gaussian (M09, G11, Krantz, Malarchick); Hann sin² = raised cosine (J25); tanh (M09). "DRAG is largely impervious to envelope shape, provided it has the right area and starts and ends at 0" (M09) | raised cosine |
| **1Q gate error** | **meas:** 0.13 to 0.19% Pauli (Sycamore); RB average 6.3e-4 median (Willow); 2.8e-4 median, from 1.4e-4 to 4.0e-4 (IBM-Q); F ≳ 0.99 routine (Krantz). **sim:** 5.3e-5 with T1 = 40 µs (M09) | 5.5e-5 (best DRAG) and 1.3e-4 (best without DRAG), gate error |
| **Leakage per gate** | Not quantified in any of the nine papers. Willow cites leakage qualitatively as the likely cause of excess logical error. | 2.6e-6 (best DRAG), 1.6e-5 (best without DRAG), mean over six inputs |
| **Levels kept** | 3 (M09, Malarchick); 3-4 (Koch); d = 5 with RWA limit d < 7 (G11); 4 analytic and 6 simulated (J25) | 5, with a 6-level check |
| **Speed limit** | Anharmonicity, through pulse bandwidth: leakage plus phase error (Koch, Krantz, M09). 2Q gate: "too fast" means leakage (Sycamore). | same physics |

## Consequences for the tutorial

- **Gate length.** Measured single-qubit gates on current devices are 35 to 35.5 ns. The tutorial's 8 to 22 ns lengths are fast-gate territory. They lie within the research range (6.8 ns at 99.98%) and below the ≈ 20 ns where calibrated DRAG suffices in Malarchick's simulations. A calibration length of about 20 to 35 ns would match devices.
- **Anharmonicity.** −200 MHz is well supported: Sycamore −208, Willow −196, Krantz "200 to 300".
- **Coherence.** T1 = 100 µs is supported (IBM median 101, Willow 89, Krantz "above 100 µs"). T2 = 100 µs is at the upper end of the measured medians (Willow 64.5, IBM 84). T2 slightly below T1, around 60 to 85 µs, would be more representative.
- **DRAG coefficient.** λ = 1 is the textbook leakage value. Our 1.41 to 1.49 at 8 to 10 ns lies above J25's optimized 1.0 to 1.3, which J25 found for 7 to 15 ns Hann pulses at −225 MHz with a different objective. A longer calibration length should bring λ toward 1.
- **Symbol clash in the literature.** In M09, G11 and J25, λ is the √2 coupling ratio, not the DRAG coefficient.
