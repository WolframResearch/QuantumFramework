# DRAG papers

| Code | arXiv | Paper |
|---|---|---|
| M09 | 0901.0534 | Motzoi, Gambetta, Rebentrost, Wilhelm, PRL 103, 110501 (2009) |
| G11 | 1011.1949 | Gambetta, Motzoi, Merkel, Wilhelm, PRA 83, 012308 (2011) |
| R25 | 2511.12799 | R. Malarchick, "When does numerical pulse optimization actually help? Error budgets, robustness tradeoffs, and calibration guidance for transmon single-qubit gates" (2025). Simulation only: "no pulses were executed on physical hardware" (line 251). |
| J25 | 2512.19919 | Da Costa Jesus, Li, Gao, Barends, Cárdenas-López, Motzoi, "Analytical blueprint for 99.999% fidelity X-gates on present superconducting hardware under strong driving" (2025, FZ Jülich) |

Sources: `../tex/0901.0534/drag.tex` (366 lines), `../tex/1011.1949/Drag_v2.tex` (945), `../tex/2511.12799/main.tex` (276), `../tex/2512.19919/main.tex` (849). All read in full.

Labels:
- **SIM**: simulation result
- **DES**: design or theory choice
- **MEAS**: experimental value the paper cites
- **ILL**: illustrative

## Points that matter for TransmonGates

1. **The "DRAG suffices at ≳20 ns" claim** is R25's abstract, verbatim (line 41): "(1) properly calibrated DRAG is sufficient for gate times ≳ 20 ns on hardware with T_2/T ≳ 500". It is a recommendation drawn from simulation.
   - It reappears at lines 235 and 258.
   - Their own parameters give T2/T = 9.6 µs / 20 ns = 480 (reading agent's arithmetic).
2. **R25's "DRAG" is λ = 0.5 in our convention, with no detuning.**
   - R25 uses β = −1/(2α) and attributes it to M09. M09's correction is actually λ = 1. The ½ matches G11's half-derivative "Y-only" solution.
   - So R25's numbers describe an uncalibrated half-DRAG: crossover with GRAPE near 15 ns, DRAG F = 0.967 at 10 ns.
3. **λ name clash.** In M09, G11 and J25, λ (λ1, λ2) is the 1-2/0-1 coupling ratio, ≈√2. Our λ (the DRAG coefficient) corresponds to G11's ansatz β (line 405) and to J25's α (line 421).
4. **J25's constant-detuning formula fits our pulse.** For a Hann (raised-cosine) π pulse, line 482: δ_c ≈ 0.712 (λ2² − 4α) π² / (Δ2 T²).
   - At λ2² = 2, Δ2 = −2π·0.2 rad/ns, T = 10 ns and α = 1.49, this gives δ_c ≈ +0.22 rad/ns.
   - Our calibrated δ = +0.175 rad/ns corresponds to α ≈ 1.28 by this formula. The sign agrees with ours.
   - This is the reading agent's arithmetic. The paper says agreement with numerics improves with gate time (486).
5. **Dephasing-rate convention.** R25 (108) and J25 (605) put γ_φ = 1/T2 − 1/(2T1) on D[n̂].
   - A rate γ on D[n̂] damps the 0-1 coherence at γ/2, so their realized T2 is longer than the T2 they quote.
   - TransmonGates uses L_φ = √(2/T_φ) b†b, rate 2/T_φ. That gives 0-1 coherence decay at 1/T2 exactly. This was checked in a kernel earlier in the TransmonGates work.

## M09: Motzoi et al. 2009

- **Levels:** "We consider the lowest three levels" (94).
- **Parameters:** "λ = √2, Δ = 2π(−400 MHz) (typical for the Phase ... and Transmon ... qubits)" (157).
- **Pulse shapes:**
  - Shifted Gaussian E_G and tanh "tangential" E_T (136-144).
  - Fig. 1: Gaussian σ = 0.5 t_g; tangential σ = 0.6 t_g − 1.4 (157).
  - Fig. 2: "6ns analog Gaussian (σ = t_g/2) with DRAG" (203).
- **Without DRAG:** "the error is 1% at a gate time of 6 ns", and it scales as (E_max^x)²/Δ (144). SIM, at −400 MHz, 3 levels.
- **With DRAG:**
  - "DRAG reduces the gate error by orders of magnitude (e.g. 4 orders at 6ns)" (193).
  - "DRAG is largely impervious to envelope shape, provided it has the right area and starts and ends at 0" (193).
- **First-order DRAG** (175): E^y = −Ė^x/Δ, with δ1 = (λ² − 4)E_x²/(4Δ), removing leakage to order E_x⁴/Δ³ (177). In our convention this is the full correction, λ = 1, with a time-dependent detuning.
- **Higher order:** fifth-order formulas at 182-185.
- **Phase ramping:** use the mean detuning δ̄1 and rotate the quadratures (195-197).
- **GRAPE:**
  - 1 ns pixels (218); converges to 1 − Φ2 < 10^-5 (213).
  - With one control the minimum time saturates "at a time of 2π/Δ"; two controls allow arbitrarily short pulses (220).
- **Earlier pulses had a minimum time** (88): "these pulses are limited to a minimal time, given by the inverse of the anharmonicity".
- **Decoherence:**
  - T1 = 40 µs, with pure dephasing set to zero "consistent with the Transmon qubit" (238).
  - Optimal errors 5.3e-5 (DRAG) and 3.4e-5 (GRAPE), SIM.
  - "current experimental values T_1 ≈ (1-4) µs" (239, 2009, MEAS cited).
- **Not stated:** qubit frequency, T2.

## G11: Gambetta et al. 2011

- **Model:** Δ_j = Δ2 (j−1)j/2, a Duffing oscillator in the RWA (116), with λ_j = √j (141).
  - The appendix says cQED drives can differ, so λ_j should be "used as a fitting parameter in any experiment" (625).
- **Levels:**
  - RWA limit d < √(2ω/Δ), specifically "d = 7" (165).
  - Simulations use d = 5 (194, 208, 297, 408); extension examples use d = 6.
- **Gaussian NOT pulse, no DRAG:** offset Gaussian, A = π, t_g = 4σ (186-189).
  - σ = {1/3, 2/3, 3/2}·2π/Δ2 give errors 0.198, 0.0160 and 0.0030 (208).
  - At |Δ|/2π = 200 MHz that is t_g ≈ 6.7, 13.3 and 30 ns (reading agent's conversion).
- **Small parameter:** ε = 1/(t_g Δ) (248).
- **First-order family** (our λ in the last column, assuming λ1 = √2):

| Solution | Ω_y | δ(t) | our λ | Lines |
|---|---|---|---|---|
| Z-only | 0 | λ1²Ω_G²/(4Δ2) | 0 | 331 |
| Y-only (Chow 2010 "simple DRAG", Lucero 2010 "half derivative") | −λ1²Ω̇_G/(4Δ2) | 0 | 0.5 | 344-347 |
| Optimal first order | −Ω̇_G λ1/(2Δ2) | Ω_G²(λ1² − 2λ1)/(4Δ2) | ≈ 0.707 | 361-362 |
| M09 | −Ω̇_G/Δ2 | Ω_G²(λ1² − 4)/(4Δ2) | 1 | 371-372 |

- **M09's first-order solution** is "not optimal" (374), but becomes "the best solution" when corrected to second order (388). The second-order amplitude is Ω_x = Ω_G + (λ1² − 4)Ω_G³/(8Δ2²) (387).
- **Numerical fit** (404-409): Ω_x = αΩ_G, Ω_y = −βΩ̇_G/Δ2, δ = γΩ_G²/Δ2 + δ0, in d = 5.
  - A time-varying δ "does not lead to any improvements" once a constant δ0 is allowed.
  - Conjecture: the optimum is a pulse "with a frequency that is not equal to the transition frequency of qubit". This supports TransmonGates' constant δ.
- **Phase ramping:** Appendix B (636-663).
- **Commented-out table, not published:** transmon ω/2π 6-10 GHz, Δ/2π from −200 to −400 MHz.
- **Not stated:** T1, T2, T_φ, gate times in ns.

## R25: Malarchick 2025 (simulation only)

- **Parameters** (112): "representative of IQM's Garnet 20-qubit processor: ω_q/2π = 5.0 GHz, α/2π = −200 MHz, T_1 = 37 µs, T_2 = 9.6 µs". Baseline T = 20 ns.
- **Model:** three levels (68). H_d = (α/2) n̂(n̂ − I), H_x = ½(a + a†), H_y = ½ i(a† − a). This is our operator form.
- **Pulses:**
  - Gaussian σ = T/8 (84).
  - DRAG β = −1/(2α) = 0.398 ns (86), which is λ = 0.5 in our convention, with no detuning.
  - GRAPE with 50 slices.
- **Closed-system three-level X gate, Table I (137-143).** F is fidelity, P2 the population left in level 2.

| T (ns) | Gaussian F | Gaussian P2 | DRAG F | DRAG P2 |
|---|---|---|---|---|
| 10 | 0.786 | 1e-1 | 0.967 | 3e-2 |
| 15 | 0.939 | 1e-2 | 0.996 | 3e-3 |
| 20 | 0.972 | 6e-4 | 0.9995 | 1e-4 |
| 25 | 0.982 | 1e-5 | 0.99986 | 1e-6 |
| 30 | 0.988 | 3e-8 | 0.99993 | 6e-8 |
| 50 | 0.996 | 1e-9 | 0.999991 | 1e-9 |
| 100 | 0.999 | 3e-10 | 0.999999 | 3e-10 |

- **Crossover:** around 15 ns. GRAPE is "necessary at short gate times (≲ 15 ns)" (154, 237).
- **Error budget at 20 ns** (169-174): DRAG coherent 4.9e-4; T1 + T2 8.4e-4; DRAG only 1.2 times worse than GRAPE.
- **Analytic floor:** ε_T1 ≈ T/(2T1) and ε_φ ≈ T/T2 − T/(2T1) (193).
- **Robustness:** sweeps of ±5 MHz and ±5% (200-213). "frequency drift of O(1 MHz) is common between calibration cycles" (239).
- **Internal inconsistency:** the two-level DRAG result is F = 0.993 (122), but Table II gives coherent error 2.8e-2 (169).

## J25: Da Costa Jesus et al. 2025

- **Parameters:**
  - Δ2/2π = −0.225 GHz "used for all simulations" (177); Δ3 = 3Δ2; λ_j ≈ √j.
  - Levels: 4 for the analytic derivations, 6 in simulation (177).
  - ω_q = √(8E_J E_C) − E_C and Δ = −E_C (649).
- **Pulses:**
  - Hann Ω0 sin²(πt/T), the raised cosine, the same as TransmonGates (177).
  - Ω_DRAG = Ω_Hann − iΩ̇_Hann/Δ2, which is λ = 1.
  - R1D uses sin³, R2D sin⁴, plus a cosine series (233, 301, 309).
  - The sign convention for σ^y is opposite to M09 and G11 (124).
- **Detuning:**
  - δ(t) = −Ω_x²(4 − λ2²)/(4Δ2), "often approximated by a constant detuning shift" (153).
  - "counterintuitively, a constant detuning can outperform time-dependent detuning in relevant regimes" (95).
- **Uncalibrated linear DRAG** (179, SIM):
  - Errors below 10^-4 for t > ≈12.7 ns.
  - Fidelity below 99.9% for t < ≈10 ns.
- **Prefactor triad** (421): Ω → β(Ω_x + αΩ_y), δ → δ_c. Here α is our λ.
  - α = 1 "should completely compensate" 1-2 leakage.
  - The analytic estimate is ≈ 1.2; numerically optimized values run from 1.0 to 1.3 (536-537).
- **Experiment cited** (PRL 116, 020501): "δ_c = 0 roughly at α = 0.5", with δ_c ∝ 1/T² (486).
- **Constant-detuning formula** (482): δ_c ≈ 0.712 (λ2² − 4α) π²/(Δ2 T²).
- **Optimized DRAG:**
  - Error below 1e-5 between 8.93 and 9.33 ns, through destructive interference of leakage channels (426).
  - The paper warns that relying on such a specific duration is problematic.
- **Minimum gate time:** T_min = √(2n) π/Δ2 (394); about 4.45 ns for the cosine pulse (403).
- **Decoherence** (605): (T1, T2) = (40, 56), (350, 350) and (500, 500) µs, described as "achievable today"; (1000, 1000) µs for the future.
- **Experiment by Gao et al.** (620): 6.8 ns, 99.98%.

## Summary table

| Parameter | Values stated | Where |
|---|---|---|
| Δ/2π | −400 MHz ("typical", 2009); −200 MHz (Garnet-like); −225 MHz; G11 table −200 to −400 MHz | M09:157; R25:112; J25:177; G11:124-136 (commented) |
| ω_q/2π | 5.0 GHz; G11 commented table 6-10 GHz | R25:112 |
| T1 | 40 µs (sim); 1-4 µs (2009); 37 µs (Garnet-like); 40/350/500 µs "achievable today" | M09:239; R25:112; J25:605 |
| T2 | 9.6 µs; 56/350/500 µs | R25:112; J25:605 |
| π-pulse time | 6 ns (M09); 10-100 ns sweep, 20 ns baseline; 6.75-12.7 ns; 6.8 ns experiment; min ≈ 4.45 ns | M09:144,203; R25:112,137; J25:179,569,620,403 |
| Pulse shapes | Gaussian, tanh, offset Gaussian, Hann sin² (= raised cosine), sin³, sin⁴ | M09:139; G11:187; R25:84; J25:177,233,301 |
| DRAG coefficient (our λ) | 1 (M09); 0, 0.5, 0.707, 1 (G11 family); 0.5 (R25); 1.0-1.3 optimized, ≈ 1.2 analytic (J25) | M09:175; G11:331-372; R25:86; J25:536 |
| Detuning | time-dependent δ1 (M09); constant δ0 suffices (G11); δ_c ∝ 1/T² formula (J25) | M09:175; G11:409; J25:482 |
| Levels kept | 3; d = 5 (RWA limit 7); 3; 4 analytic and 6 simulated | M09:94; G11:165,208; R25:68; J25:177 |
| DRAG ≥ 20 ns sufficient | R25 recommendation (simulation, T2/T ≳ 500) | R25:41,235,258 |
