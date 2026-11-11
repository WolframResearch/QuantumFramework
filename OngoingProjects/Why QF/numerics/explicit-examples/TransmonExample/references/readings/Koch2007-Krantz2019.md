# Koch et al. 2007 (cond-mat/0703002) and Krantz et al. 2019 (arXiv 1904.06560)

Sources: TeX in `../tex/cond-mat_0703002/transmon.tex` (1342 lines) and `../tex/1904.06560/Manuscript_Krantz_ArXiv_v4b_blackline.tex` (2534 lines). Both were read in full.
- The Krantz file is a "blackline" version. Its only change markup is `%`-commented `\textcolor{red}` blocks.
- No on-topic value appears only in deleted text. The commented numbers are Krantz 830 (T_phi,G = 98 µs, the same as the final text) and Krantz 2343 (an off-topic Purcell example).

Labels used below:
- **M**: measured on a device (from a cited experiment)
- **S**: simulation
- **E**: theoretical estimate
- **D**: design target or typical value
- **I**: illustrative example

## 1. Single-qubit gate duration and pulse shape

**Koch (2007):**
- 243, **E**: "minimum pulse duration to be τ_p ~ |ω01 α_r|^-1", from the transform-limited bandwidth.
- 243, **E**: with dephasing times of a few hundred ns, "reasonable pulse durations would be in the range of several tens of nanoseconds".
- 243: "significantly shorter microwave pulses are difficult to achieve" (2007 hardware).
- 243-245, **D/I**: "Using a typical pulse length of 10 ns" gives |α_r^min| ~ 1/200π at ω01/2π = 10 GHz.
- 247, **E**: with dephasing "of the order of 20 µs ... pulse durations much longer then 10 ns could be used".
- 661 (Table I caption), **D** for 2007 experiments: "pulse durations usually range between 1 and 10 ns".
- No pulse shape is named.

**Krantz (2019):**
- 1530, **I**: "Gaussian pulses with standard deviations σ = {1,2,5} ns" overlap spectrally with the 1-2 transition.
- 1534 (Fig. 14), **S/I**: an Xπ pulse with and without DRAG, simulated with QuTiP mesolve. Only σ is given, not the total length.
- 1539-1547: the DRAG envelope is s(t) on I and λ ṡ(t)/α on Q. This is the same form as our H_lab.
- A single-qubit gate time in ns is not stated.
- Cosine envelopes are not stated. Slepian envelopes appear only for noise spectroscopy (1123) and the CPHASE trajectory (1767).

## 2. Anharmonicity, E_C, qubit frequency, E_J/E_C

**Koch (2007):**
- Definitions: α ≡ E12 − E01; α_r ≡ α/E01 (217-218).
- Model, 238-239: α ≃ −E_C; α_r ≃ −(8E_J/E_C)^{-1/2}.
- α_r changes sign at E_J/E_C ≈ 9 and has a local minimum near 17.5 (221).
- Qubit frequency:
  - "ω01/2π ≈ 10 GHz typical for experiments" (243, 2007)
  - 7 GHz measured (404, Schuster 2007)
  - 8 GHz illustrative (428)
  - "Typical qubit frequencies are of the order of 1-10 GHz" (661)
- E_C/h:
  - 0.35 GHz (Fig. 4 at 196; 521, 591, 620, 641)
  - 300 MHz, a reachable design value (404)
- E_J/E_C:
  - "20 ≲ E_J/E_C ≪ 5×10^4" (247)
  - "several tens up to several hundreds" (101)
  - ~10^2 (252, 667)
  - 60 (Fig. 10) and 85 (Table I)
- Junction asymmetry d ≃ ±10% (273).

**Krantz (2019):**
- 341, **D**: "α = −E_C is usually designed to be 100-300 MHz", with ω_q = (√(8E_J E_C) − E_C)/ħ = 3-6 GHz and "E_J/E_C ≥ 50".
- 1530, **D**: α "is negative and typically around 200 to 300 MHz". This differs from line 341.
- 346: Duffing model H = ω a†a + (α/2) a†a†aa.
- 848, **D**: "Typical qubits are designed at frequency ω/2π ≈ 5 GHz", at T ≈ 20 mK.
- 936: charge dispersion matters if E_J/E_C is "smaller than around 60".
- 1534 (Fig. 14, **S/I**): α = −200 MHz, ω_q = 4 GHz.
- 2026 (Fig. 20, **I**): α/2π = −300 MHz.
- 1146: tunability "typically around 1 GHz".

## 3. T1, T2, T2*, T2 echo, T_phi

**Koch (2007).** Estimates (**E**), by channel:
- T1:
  - radiative 0.3 ms at 8 GHz (428)
  - Purcell 16 µs (460)
  - quasiparticles ~1 s (480)
  - flux line 20 ms to 1 s (521)
  - flux coupling to the whole circuit 70 ms (527)
- T2 (pure dephasing):
  - charge ~8 s (591); slow charge noise 0.4 ms (606)
  - flux 1 µs at Φ0/4 (620), 3.6 ms at the sweet spot (626)
  - critical current 35 µs, "likely to be the limiting dephasing mechanism" (641)
- Projected dephasing "of the order of 20 µs" (247).
- 661: transmon coherence "should be limited by relaxation (T1) processes only".
- CPB context: T2 of "a few hundreds of nanoseconds" (243, **M**).

**Krantz (2019).** Device data from Kjaergaard 2017, Fig. 5:
- T1 = 85 µs (815, 852, **M**).
- T2* = 95 µs (Fig. 5 caption, 819) but 98 µs in the text (868). The paper is inconsistent here.
- T2E = 120 µs (823, 870, **M**).
- Gaussian T_phi,G = 98 µs (902, **M**).
- Relations: Γ2 = Γ1/2 + Γφ (776); "T2 can never be larger than 2T1" (891).
- Trend, 2341: "intrinsic coherence times for superconducting qubits have reached above 100 µs".

## 4. Drive strength, DRAG coefficient, detuning

- Rabi frequency compared with α: not stated in either paper.
- Krantz 1547, **E**: "The theoretically optimal choice for reducing dephasing error is λ = 0.5 and an optimal choice for reducing leakage error is λ = 1" (citing Gambetta 2011 and Motzoi 2013). Fig. 14 uses λ = 0.5.
- Krantz 1549:
  - "In practice there can be a deviation from these two optimal values, often due to pulse distortions".
  - λ is calibrated by RB with single-shot measurement of |2⟩.
  - A detuning δf is added, s′(t) e^{i2πδf t} (1551). No value is given.
- Koch: DRAG is not discussed.

## 5. Gate fidelity and leakage

- Krantz 1553: "Modern single-qubit gates using DRAG pulsing now routinely reach fidelities F_1qb ≳ 0.99".
- Leakage per gate: not stated in either paper.

## 6. Levels kept

- Koch: 3 levels (Fig. 2, 138), 4 levels (Fig. 4, 196), and the lowest three states in the transformation (762).
  - 362: "virtual transitions through excited transmon states have to be taken into account".
- Krantz: the 2014 footnote says higher levels "need to be taken into account for most qubit simulations". The Fig. 14 level count is not stated.

## 7. What limits gate speed

- Koch 97: anharmonicity "limits the speed of qubit operations".
- Koch 214: it "ultimately sets a lower bound on the duration of control pulses".
- Koch 243: τ_p ~ |ω01 α_r|^-1.
- Krantz 314: "the amount of anharmonicity sets a limit on how short the pulses ... can be".
- Krantz 1530-1538: the spectral overlap of short pulses with the 1-2 transition causes leakage plus phase error from level repulsion.

## Summary table

| Parameter | Values stated | Where |
|---|---|---|
| π-pulse length | 1-10 ns typical (2007); 10 ns "typical"; several tens of ns reasonable; τ_p ~ \|ω01 α_r\|^-1 | Koch 661, 243-245 |
| Pulse shape | Gaussian σ = 1, 2, 5 ns (illustrative); DRAG I = s, Q = λṡ/α | Krantz 1530, 1539-1547 |
| α/2π | −(100 to 300) MHz design; −(200 to 300) MHz typical; −200 MHz (sim); α ≃ −E_C | Krantz 341, 1530, 1534; Koch 238 |
| ω_q/2π | 3-6 GHz design; ≈5 GHz typical; 1-10 GHz (2007) | Krantz 341, 848; Koch 661 |
| E_J/E_C | ≥ 50; charge-sensitive below ~60; tens to hundreds | Krantz 341, 936; Koch 101 |
| T1 | 85 µs measured; > 100 µs trend | Krantz 852, 2341 |
| T2* / T2E / T_phi | 95 (or 98) µs / 120 µs / 98 µs measured; T2 ≤ 2T1 | Krantz 819, 868, 823, 902, 891 |
| DRAG λ | 1 for leakage, 0.5 for phase error; deviates in practice | Krantz 1547-1549 |
| DRAG detuning | introduced, no value | Krantz 1551 |
| 1Q fidelity | ≳ 0.99 routine | Krantz 1553 |
| Rabi Ω vs α; leakage per gate | not stated | none |
| Speed limit | anharmonicity, through bandwidth: leakage and phase error | Koch 97, 214, 243; Krantz 314, 1538 |

Our doc's convention matches Krantz's DRAG form, Q = λṡ/α, where λ = 1 cancels leakage to first order.
