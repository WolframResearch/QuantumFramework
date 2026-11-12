# Four Hardware Cases for Qolab: Formulation and Feasibility

Target hardware assumed throughout: **Arbel**, the 21-qubit IQCC backend that QUALink already runs on (Quantum Machines OPX1000; drive and readout on MW-FEM outputs, flux on LF-FEM outputs). Device values below come from the QUALink fixture snapshot of 2026-07-07 and must be re-read from the live calibration before they are used. If Qolab means a different chip, the formulations and the simulation plans stand; the hardware notes do not.

**Ranking:** DRAG, then charge dispersion with three-state readout, then cryoscope, then sweet spot / CPMG. All four can be simulated in QF+WL, and none needs TransmonLab. The hardware route for all four is QUALink, which needs a different small addition for each.

| Case | Simulation in QF+WL | TransmonLab | QUALink today | Main risk |
|---|---|---|---|---|
| 1. Charge dispersion, parity, three-state readout | yes; the spectrum is exact in plain WL | not needed: it has the spectrum, not three-state readout | flux bias and per-shot I/Q work; no EF pulse, no GEF readout | readout and drive at the low-frequency flux point |
| 2. Cryoscope | yes; mostly WL system identification, QF minor | cannot: no distortion stage in its pipeline | plays flux pulses; emits no output filter | time cost; 1 ns steps on a 4 ns clock |
| 3. Sweet spot / CPMG | yes, with filter functions; QF's Lindblad equation cannot produce 1/f decay | has the noise formulas; optional | Ramsey, echo, DC bias work | wiring noise during APS |
| 4. DRAG | yes; TransmonGatesV3 is most of it | wrong engine | explicit {I, Q} pulses work; no GEF readout | at reachable lengths the leakage is too small to see in fidelity |

---

## 1. Circuit to Spectrum: Charge Dispersion, Parity and Three-State Readout

### Formulation

The transmon is a Cooper-pair box with $E_J \gg E_C$:

$$H = 4E_C(\hat n - n_g)^2 - E_J\cos\hat\varphi .$$

Its levels are Mathieu characteristic values: the $m$-th smallest of

$$E_C\, a_{2(n_g+k)}\!\left(-\frac{E_J}{2E_C}\right),\qquad k\in\mathbb Z,$$

is $E_m(n_g)$. Each level oscillates with the offset charge,

$$E_m(n_g) \simeq E_m(\tfrac14) - \frac{\epsilon_m}{2}\cos 2\pi n_g,\qquad
\epsilon_m \simeq (-1)^m E_C\,\frac{2^{4m+5}}{m!}\sqrt{\frac{2}{\pi}}\left(\frac{E_J}{2E_C}\right)^{\frac m2+\frac34} e^{-\sqrt{8E_J/E_C}},$$

so the band width grows with $m$ and falls exponentially in $\sqrt{8E_J/E_C}$, and the transition from $|1\rangle$ to $|2\rangle$ is far more charge sensitive than the qubit transition.

A SQUID makes $E_J$ flux tunable,

$$E_J(\Phi) = E_{J\Sigma}\left|\cos\frac{\pi\Phi}{\Phi_0}\right|\sqrt{1+d^2\tan^2\frac{\pi\Phi}{\Phi_0}},\qquad hf_{01}\simeq\sqrt{8E_J(\Phi)E_C}-E_C,$$

so lowering the qubit frequency lowers $E_J/E_C$ and widens the band.

One quasiparticle tunneling across the junction changes $n_g$ by $1/2$. When parity switches faster than one spectroscopy trace averages, the 1→2 line appears at both $f_{12}(n_g)$ and $f_{12}(n_g+\tfrac12)$, split by

$$\Delta f_{12} = \frac{|\epsilon_2-\epsilon_1|}{h}\,|\cos 2\pi n_g| .$$

$n_g$ drifts between traces, so repeated traces sample the splitting from zero up to $|\epsilon_2-\epsilon_1|/h$, the full charge dispersion of the transition.

To read out three states, note that with the transmon in $|j\rangle$ the resonator sits at

$$\omega_r^{(j)} = \omega_r + \chi_{j-1,j} - \chi_{j,j+1},\qquad \chi_{ij} = \frac{g_{ij}^2}{\omega_{ij}-\omega_r},\qquad \chi_{-1,0}=0,$$

and a readout drive $\varepsilon$ at $\omega_d$ leaves it in the coherent state

$$\alpha_j = \frac{-i\varepsilon}{\kappa/2 + i\big(\omega_r^{(j)}-\omega_d\big)} .$$

The three pointer states give three IQ blobs. Two blobs of width $\sigma$ whose centers lie $d$ apart are confused with probability $\tfrac12\,\mathrm{erfc}\big(d/2\sqrt2\,\sigma\big)$. $\chi_{ij}$ falls as the qubit is tuned away from the resonator, so readout degrades exactly where the dispersion is largest.

### Scale on Arbel

qA1 ($f_{01}$ = 5.08 GHz, $|\alpha|/2\pi$ = 213 MHz) sits near $E_J/E_C \approx 90$, $E_C/h \approx 195$ MHz, where the 1→2 band is far too narrow to resolve. For $E_C/h$ = 250 MHz the exact spectrum gives a 1→2 dispersion of about 0.3 MHz at $E_J/E_C$ = 50, 1.7 MHz at 40, 10 MHz at 30 and 26 MHz at 25. On qA1 that means tuning down to roughly 2.5-3.3 GHz. The stored flux curvature and flux period match a nearly symmetric SQUID, so that range should exist.

### Feasibility

- **QF+WL.** The spectrum comes from the charge-basis matrix or `MathieuCharacteristicA`. The drift of $n_g$ and the parity switching are WL random processes (`TelegraphProcess` for parity, a random walk for $n_g$). The line shape is a three-level driven Lindblad equation in QF. The readout is the three pointer states plus amplifier noise and decay during readout, classified with `ClusterClassify` / `Classify` and compared with Arbel's measured 3×3 assignment matrix.
- **TransmonLab.** Has `CPBHamiltonian`, `TransmonSpectrum`, `ChargeDispersion`, `CPBSpectrumPlot` and a flux model. All its IQ tools are two-state. It saves a few lines on the spectrum and gives nothing for the readout.
- **Hardware.** Arbel calibrates every qubit's `EF_x180` pulse and `readout_GEF` operation, with a 3×3 assignment matrix (diagonal about 0.99, 0.93-0.95, 0.90). QUALink sets the DC flux bias and saves per-shot I/Q, but it does not play `EF_x180` or `readout_GEF`. Adding those two from the calibration state is the hardware work.
- **Prior work.** TransmonLab has only the static picture (`TransmonSpectrum` with `"ChargeOffset"`, `CPBSpectrumPlot`, Demo 1 of `Examples/TransmonLab_Examples.md`, a CPB-versus-$n_g$ toy notebook); nothing on the 1→2 line, parity or three-state readout. QUALink has run a DC flux sweep on qC2 (`LIVE_08_FluxSweep_Arbel.wl`), resonator response only. Flagged: TransmonLab's `Docs/WhatIsTransmonCQED-v2.md` puts the ground-level charge dispersion at $E_J/E_C$ = 50 near $10^{-9}E_C$ (its own `ChargeDispersion` page gives $5.6\times10^{-7}E_C$) and concludes that charge noise is irrelevant; the 1→2 transition disperses by $|\epsilon_2-\epsilon_1|\approx1.3\times10^{-3}E_C$, the splitting this case measures.
- **Opinion.** Most physics per unit of effort, and the parity doubling is real-device physics a paper does not give students. The risk sits at the low-frequency flux point: weaker dispersive shift, an xy chain that must pass about 2.5-3.5 GHz, and an unknown $T_1$ there. Ask Qolab what frequency range their xy chain covers before committing.

---

## 2. Cryoscope

### Formulation

Near its parking point the qubit frequency is quadratic in flux,

$$f(\Phi) \simeq f_{\max} + c\,\Phi^2,\qquad c<0,$$

and Arbel stores $c$ per qubit (`freq_vs_flux_01_quad_term`, in Hz/V²). A Ramsey sequence brackets a flux pulse cut off at time $\tau$: $\pi/2$, flux pulse, then $\pi/2$ about $x$ or $y$. It measures

$$\langle\sigma_x\rangle + i\langle\sigma_y\rangle = e^{-\Gamma\tau}\,e^{i\varphi(\tau)},\qquad \varphi(\tau) = 2\pi\int_0^\tau \Delta f(t)\,dt,$$

so the detuning the qubit saw at time $\tau$, and the flux that reached the chip, are

$$\Delta f(\tau) = \frac{1}{2\pi}\frac{d\varphi}{d\tau},\qquad \Phi(\tau) = \sqrt{\Delta f(\tau)/c}\, .$$

The lines between the controller and the chip act as a linear time-invariant filter, $\Phi_{\rm chip} = h * \Phi_{\rm AWG}$. Bias tees, cable skin effect and on-chip lines give a step response well described by a sum of exponentials,

$$s(t) = \Big(1+\sum_k A_k e^{-t/\tau_k}\Big)\Theta(t),\qquad H(s) = 1+\sum_k \frac{A_k\, s}{s+1/\tau_k},$$

with $A_k>0$ for overshoot and $A_k<0$ for slow undershoot. Pre-distortion applies $G = 1/H$ before the line, so $GH = 1$ and the step arrives flat. For one exponential the inverse has a single pole,

$$G(s) = \frac{s+1/\tau}{(1+A)\,s + 1/\tau},$$

which is stable when $A>-1$. The controller implements $G$ as a digital filter,

$$y[n] = \sum_k b_k\, x[n-k] + \sum_m a_m\, y[n-m],$$

with FIR taps $b_k$ and IIR taps $a_m$.

### Feasibility

- **QF+WL.** The quantum part is one qubit with a time-dependent detuning, and its phase is the integral above, so QF adds little. The substance is WL system identification: `TransferFunctionModel`, `OutputResponse`, `ToDiscreteTimeModel`, `RecurrenceFilter` for IIR, `ListConvolve` for FIR, `NonlinearModelFit` for $(A_k,\tau_k)$, and a stability check of $G$. The notebook should implement Quantum Machines' exact filter equations so the fitted taps paste into the configuration unchanged.
- **TransmonLab.** Cannot. Its own QUA audit (Case 16) marks output filtering as not implementable: nothing sits between the designed envelope and the solver.
- **Hardware.**
  - Arbel's z lines are LF-FEM outputs, which support filters (one FIR with 48 taps at 0.5 ns, six single-pole IIR filters, configuration fields `feedforward` and `exponential`). MW-FEM outputs have none.
  - Every z line in the fixture has `filter_fir_taps = None` and `filter_iir_taps = None`, and Arbel plays its CZ gates as flux pulses on these lines. The uncorrected pulse Qolab wants students to see is therefore what the chip plays now, and a working pre-distortion would help Arbel's CZ, not only the course.
  - QUALink emits no `filter` field. The $\tau$ sweep needs 1 ns steps on a 4 ns clock, meaning one precomputed waveform per $\tau$; QUALink's custom sweep makes per-value waveforms but is untested for this.
- **Prior work.** None on the physics. TransmonLab's `Case16_OutputFilter_Improvements.md` is a design sketch for a filter stage in `BuildTimeline`, and its `Case15` notes that the solver has no time-dependent qubit frequency. QUALink already plays flux pulses at their written amplitude (`ClaudeNotes/why symbolic/symbolic-wl-for-hardware.md`), so the step can be authored; the taps cannot be emitted.
- **Opinion.** Correct but slow, as Qolab says:
  - every $\tau$ point is a two-quadrature Ramsey with many shots;
  - long time constants need long windows;
  - the filter must be iterated.

  Quantum Machines already ships filter-tap tools, so the Wolfram value must be the modeling: identify $H$, invert it, prove $G$ stable. Filter state persists across programs, so loading filters on a shared z line needs IQCC's agreement. Good second-phase project, weak live demo.

---

## 3. Sweet-Spot Mapping and CPMG

### Formulation

From the SQUID formula of case 1, $f_{01}(\Phi)$ has a maximum at $\Phi=0$, where the flux sensitivity

$$D_\Phi \equiv 2\pi\,\frac{\partial f_{01}}{\partial\Phi}$$

vanishes. Flux noise $\delta\Phi(t)$ with spectrum $S_\Phi(\omega)$ shifts the qubit frequency by $D_\Phi\,\delta\Phi$. A sequence of $\pi$ pulses flips the sign of the accumulated phase with a switching function $y_N(t)=\pm1$, and Gaussian noise gives coherence $W(\tau) = e^{-\chi_N(\tau)}$ with

$$\chi_N(\tau) = \frac{D_\Phi^2}{2}\int_{-\infty}^{\infty}\frac{d\omega}{2\pi}\,S_\Phi(\omega)\,\big|Y_N(\omega,\tau)\big|^2,\qquad Y_N(\omega,\tau)=\int_0^\tau y_N(t)\,e^{i\omega t}\,dt .$$

For Ramsey ($N=0$) and Hahn echo ($N=1$),

$$|Y_0|^2 = \tau^2\,\mathrm{sinc}^2\frac{\omega\tau}{2},\qquad |Y_1|^2 = \frac{16}{\omega^2}\sin^4\frac{\omega\tau}{4},$$

and for CPMG with $\pi$ pulses at $t_j = (j-\tfrac12)\tau/N$,

$$|Y_N|^2 = \frac{1}{\omega^2}\Big|1 + (-1)^{N+1}e^{i\omega\tau} + 2\sum_{j=1}^{N}(-1)^j e^{i\omega t_j}\Big|^2 .$$

For $1/f$ flux noise, $S_\Phi = A_\Phi/|\omega|$, the Ramsey decay is Gaussian with an infrared logarithm, and the echo removes the logarithm:

$$\chi_0 \simeq \frac{D_\Phi^2 A_\Phi\,\tau^2}{2\pi}\,\ln\frac{1}{\omega_{\rm ir}\tau},\qquad \chi_1 = \frac{D_\Phi^2 A_\Phi\,\tau^2}{2\pi}\,\ln 2 .$$

For large $N$ the switching function is a square wave of angular frequency $\pi N/\tau$, which carries $8/\pi^2$ of its power in the fundamental, so the decay reads off the noise spectrum at that one frequency:

$$\chi_N(\tau) \simeq \frac{4}{\pi^2}\,D_\Phi^2\,\tau\,S_\Phi\!\left(\frac{\pi N}{\tau}\right).$$

At the sweet spot $D_\Phi=0$, and dephasing comes from the second-order term $\tfrac12\,\partial^2\omega/\partial\Phi^2\,\delta\Phi^2$, which is not Gaussian.

### Feasibility

- **QF+WL.** Possible, with one limit: QF's Lindblad equation is Markovian, so it gives exponential decay and cannot reproduce the $1/f$ Gaussian decay or echo refocusing. The formulation above is classical and closes in WL: the filter functions, $\chi_N$ by `Integrate` or `NIntegrate`, and a Monte Carlo over noise realizations as a cross-check. QF would mostly sit idle.
- **TransmonLab.** The one case where it has ready pieces: `NoisePSD`, `FilterFunction` (Ramsey, echo, CPMG: the $|Y_N|^2$ above divided by $\tau^2$), `TphiFromNoise`, `SweetSpotAnalysis`. Useful but optional, since the same formulas are a few lines of WL.
- **Hardware.** QUALink has Ramsey and echo presets and DC flux bias; CPMG is a pulse train in an explicit sequence. Arbel already stores its sweet-spot curvature as a routine calibration.
- **Prior work.** The furthest along of the four. TransmonLab's examples (Demo 2) already map $T_\varphi$ against flux, and Demo 12 plots the CPMG filter functions. On Arbel qA1 a QUALink Hahn echo matched the native QUA SDK (41.85 vs 39.05 μs, `HahnEcho_WL_QUALink_vs_Native.md`; it needs the `"Thermalization"` reset or the trace comes out flat), and Ramsey gave $T_2^*$ = 34 μs (`Ramsey.md`). The DOE Genesis dynamical-decoupling proposal (MSU lead, Qolab a partner) already assigns TransmonLab filter-function scoring of CPMG, XY8, UDD and KDD, so this case should align with it. Flagged: QUALink's `11_CPMG.wl` only compiles and is wrong (its π/2 and π are one waveform, and it plays five π pulses for every N); TransmonLab's `SweetSpotAnalysis` page gives $T_\varphi\approx$ 10 ns at 0.3 Φ₀ where its `TphiFromNoise` page gives 1.8 μs for the same system and noise.
- **Opinion.** Lowest priority. Qolab's worry about wiring noise during a live APS demo is right, and the sweet-spot map reproduces a calibration Arbel already runs. Use recorded data, or fold the flux arc into case 1, which needs the same arc.

---

## 4. DRAG

### Formulation

The five-level Duffing transmon of TransmonGatesV3, in the frame of the drive:

$$H(t) = \frac{\alpha}{2}b^\dagger b^\dagger bb + \delta\,b^\dagger b + \frac{\Omega(t)}{2}\big(b+b^\dagger\big) - \frac{\lambda\,\dot\Omega(t)}{2\alpha}\,i\big(b^\dagger-b\big),\qquad
\Omega(t) = A\,\frac{\pi}{\tau}\Big(1-\cos\frac{2\pi t}{\tau}\Big).$$

Leakage is driven by the pulse's Fourier component at the 1→2 transition. Integrating the DRAG term by parts, with $\Omega(0)=\Omega(\tau)=0$,

$$\int_0^\tau\Big(\Omega - \frac{i\lambda\dot\Omega}{\alpha}\Big)e^{i\omega t}\,dt = \Big(1-\frac{\lambda\omega}{\alpha}\Big)\tilde\Omega(\omega),$$

which vanishes at $\omega=\alpha$ for $\lambda=1$. The raised cosine has

$$\big|\tilde\Omega(\omega)\big| = \frac{2\pi A}{\tau}\;\frac{(2\pi/\tau)^2\,\big|\sin(\omega\tau/2)\big|}{|\omega|\,\big|(2\pi/\tau)^2-\omega^2\big|}\;\propto\;\omega^{-3},$$

so with amplitude and detuning calibrated, the error left without DRAG falls as

$$\epsilon_{\rm leak}\simeq\frac{8\pi^6}{(|\alpha|\tau)^6}.$$

The drive also shifts the qubit transition, by $\Omega(t)^2/2|\alpha|$ to second order, which leaves an uncalibrated pulse off target; calibrating $\delta$ and $A$ removes this. Decay and dephasing add, to first order along the rotation,

$$\epsilon_{\rm dec}\simeq K\tau,\qquad K = \frac{1+J_1(1)}{4T_1}+\frac{1-J_1(1)}{4T_2},$$

so without DRAG the best length solves $d(\epsilon_{\rm leak}+\epsilon_{\rm dec})/d\tau=0$:

$$\tau_\star^7 = \frac{48\pi^6}{K\,|\alpha|^6},$$

which reduces to $(96\pi^6T_1/|\alpha|^6)^{1/7}$ for $T_1=T_2$: 21.8 ns for V3's $\alpha/2\pi=-200$ MHz and $T_1$ = 100 μs. At the optimum $\epsilon_{\rm dec} = 6\,\epsilon_{\rm leak}$.

On the hardware, Quantum Machines' DragCosine pulse of $L$ samples at 1 ns is

$$I(t)=\frac{a}{2}\Big(1-\cos\frac{2\pi t}{L-1}\Big),\qquad Q(t)=\frac{\alpha_{\rm QM}}{2\pi\,(f_{\rm anh}-f_{\rm det})}\,\dot I(t),$$

then rotated at the detuning frequency $f_{\rm det}$. This is V3's pulse family: $\tau = L-1$ ns, $f_{\rm det}\leftrightarrow\delta$, and

$$\alpha_{\rm QM} = \pm\lambda\Big(1-\frac{f_{\rm det}}{f_{\rm anh}}\Big),$$

with the sign set by the IQ convention, which one hardware sweep fixes. The pulse area is fixed, so a qubit calibrated at length $\tau_{\rm cal}$ and amplitude $a_{\rm cal}$ (as a fraction of full scale) reaches at best $\tau_{\min}\approx \tau_{\rm cal}\,a_{\rm cal}$, if the output stays linear to full scale.

### Scale on Arbel

- **Pulses.** Arbel's x180 pulses are 48 ns DragCosine pulses (one is 52 ns, one 144 ns) at 17-96% of full scale, which gives $\tau_{\min}$ of about 15 ns on qA1, about 8 ns on qC5, and no headroom on qD5. Quantum Machines schedules on a 4 ns clock and QUALink flags durations under 16 ns, so shorter pulses go zero-padded into a longer slot.
- **Anharmonicity.** $|\alpha|/2\pi$ = 197-215 MHz, close to V3's choice.
- **Best length.** With qA1's $\alpha$, $T_1$ and $T_2$, $\tau_\star \approx 18$ ns at a total error near $3\times10^{-4}$, against $\epsilon_{\rm dec}\approx 7\times10^{-4}$ at today's 48 ns.

A closed five-level run at $|\alpha|/2\pi$ = 213 MHz (amplitude and detuning calibrated, with and without $\lambda$, local optimum) gives:

| $\tau$ (ns) | uncalibrated | calibrated, no DRAG | calibrated, DRAG | $8\pi^6/(\lvert\alpha\rvert\tau)^6$ |
|---|---|---|---|---|
| 12 | 3.0e-2 | 2.5e-4 | 2.3e-5 | 4.5e-4 |
| 16 | 1.7e-2 | 7.4e-5 | 1.3e-5 | 8.0e-5 |
| 20 | 1.1e-2 | 2.2e-5 | 2e-8 | 2.1e-5 |
| 24 | 7.3e-3 | 7.4e-6 | 7.5e-7 | 7.0e-6 |
| 48 | 1.8e-3 | 1.1e-7 | 1.0e-8 | 1.1e-7 |

The uncalibrated error falls as $1/\tau^2$ and is almost all drive-induced frequency shift; the leakage it contains is $10^{-4}$ or less. The closed form tracks the calibrated no-DRAG error from 16 ns up.

### Feasibility

- **QF+WL.** The physics is TransmonGatesV3. The translation is WL work:
  - sample at 1 ns;
  - put the length on the 4 ns grid with the 16 ns floor;
  - cap the amplitude;
  - map $(A,\lambda,\delta,\tau)$ to $(a,\alpha_{\rm QM},f_{\rm det},L)$ with the formulas above;
  - predict what the hardware measures: repeated-$\pi$ error amplification and $\lambda$ calibration sequences (compose the pulse's channel $N$ times), and the $|2\rangle$ population read by GEF readout.
- **TransmonLab.** Wrong engine:
  - its `FastDuffing` model drops the DRAG quadrature;
  - its Duffing dephasing is 4× too weak;
  - its multilevel decay skips coherence transfer;
  - its kernel audit reports a leakage-optimal $\beta\approx0$, contrary to V3.

  QUALink's `QUAPhysics` calls these solvers, so predictions must come from QF.
- **Hardware.** QUALink sends an explicit {I, Q} pair (commit `ebfe81a`), so a WL-computed waveform reaches the chip as computed. It needs GEF readout added. The xy lines are MW-FEM, which has no output filters, so any line correction has to be built into the waveform; at these bandwidths it should be small (unchecked).
- **Prior work.** QUALink already builds Arbel's calibrated DragCosine, Q quadrature included, and has played it live in T₁, echo and Ramsey runs; a symbolic sweep of the DRAG coefficient compiles, but no DRAG sweep or leakage measurement has run (`14_DRAGCalibration.wl` only compiles, and the QUA it emits is a plain amplitude sweep). TransmonLab's public APS 2026 notebook shows leakage falling from 1.4×10⁻³ to 3×10⁻⁴ (15 ns, $\alpha/2\pi$ = −150 MHz), and QUALink's APS 2027 launch case names closed-form DRAG on Arbel as its flagship. Flagged: the earlier APS state-preparation notebooks compare pulses of different rotation angle, or report $P(|1\rangle)$ instead of leakage; TransmonLab's Case 11 and Demo 8 use the two-level model, which has no $|2\rangle$; the launch case writes $Q=-(\alpha/\Delta)\,\dot I$ with the qubit-cavity detuning $\Delta$, which does not belong, and the analytic-versus-finite-difference gain it plans to show is far below the $10^{-7}$ leakage at 48 ns in the table above. The `FastDuffing` point under TransmonLab follows the June kernel audit; the current `TransmonPhysics` page says both quadratures are now propagated, so re-check it before relying on it.
- **Opinion.** The best case, with the story set by the numbers above. At the lengths Arbel reaches, DRAG's leakage benefit is invisible in gate fidelity, because $\epsilon_{\rm leak}$ sits below $\epsilon_{\rm dec}$ and Arbel's measured gate errors are about $2\times10^{-3}$. What the hardware can show:
  1. **Calibration prediction.** The model predicts the calibrated $(A,\lambda,\delta)$ at Arbel's $\alpha$, and the hardware's own calibration lands there; this is where DRAG acts as a phase correction.
  2. **Leakage amplification.** $|2\rangle$ builds up over $N$ repeated $\pi$ pulses, and the GEF readout reads it above its floor of about 1%.
  3. **Shortened pulse.** Without DRAG, the pulse shortens from 48 ns toward $\tau_\star\approx18$ ns, on a qubit with the headroom.
  4. **Strong drive.** Below about 12 ns DRAG matters for fidelity too, but only qC5-like qubits get there.

---

## How the Cases Connect

Case 1's three-state readout is the instrument case 4 uses to measure leakage, and case 1's flux arc is case 3's sweet-spot map. One QUALink extension (EF pulse and GEF readout from the calibration state) serves cases 1 and 4; filter emission serves case 2 only.

## References

- J. Koch et al., "Charge-insensitive qubit design derived from the Cooper pair box," [Phys. Rev. A 76, 042319 (2007)](https://doi.org/10.1103/PhysRevA.76.042319).
- F. Motzoi, J. M. Gambetta, P. Rebentrost and F. K. Wilhelm, "Simple pulses for elimination of leakage in weakly nonlinear qubits," [Phys. Rev. Lett. 103, 110501 (2009)](https://doi.org/10.1103/PhysRevLett.103.110501).
- P. Krantz et al., "A quantum engineer's guide to superconducting qubits," [Appl. Phys. Rev. 6, 021318 (2019)](https://doi.org/10.1063/1.5089550).
- M. A. Rol et al., "Time-domain characterization and correction of on-chip distortion of control pulses in a quantum processor," [arXiv:1907.04818](https://arxiv.org/abs/1907.04818).
- Quantum Machines, [output filter guide (QOP 1.3.1)](https://docs.quantum-machines.co/1.3.1/docs/Guides/output_filter/).
- Still to be checked against full text before citing: parity switching in transmons (Ristè et al. 2013; Serniak et al. 2018), CPMG noise spectroscopy (Bylander et al. 2011), echo decay under $1/f$ noise (Ithier et al. 2005).
