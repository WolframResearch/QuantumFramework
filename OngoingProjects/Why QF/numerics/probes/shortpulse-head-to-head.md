# Short pulses at default settings: QuTiP 5.3.1 against QF

Does QuTiP 5.3.1 still step over short pulses with every solver option at its default, and does QF's `QuantumEvolve`, which since 5 October 2026 restarts its integration at times it reads off the symbolic Hamiltonian, get the same problems right at its defaults? Each case below runs in both packages on the same Hamiltonian and is scored against an exact value.

## Verdict

QuTiP 5.3.1 at its defaults still fails every reported case. It returns the state the pulse never touched for its documentation's example, for #1265 once the wait between the two pulses is 20 or longer (10 still passes; QuTiP 4.3 failed above 23), for #1945 at 268 output times, for #2552 index93, for the #2814-style train, for every shape and width of the isolated pi pulse (the widest included, with 2 or 201 output times), and for the transmon pulse after an idle gap; the untruncated Gaussian fails even without a gap, because its drive is negligible for the first 48 ns. While the state sits still under a zero drive the local error estimate vanishes, the step grows without bound, and the integrator never samples the pulse. The `max_step` that QuTiP's docstring advises fixes every case, but at widths of 0.04 and below QuTiP then raises `IntegratorException` until `nsteps` is raised too. QF at its defaults passes every reported case within its default tolerance. The restarts are what make it pass for the Gaussian, sech and narrow `Piecewise` square pulses, the train and the transmon after a gap: with them removed, NDSolve at its defaults skips those as QuTiP does. NDSolve alone already captures `UnitStep` and `Boole` edges and the wider `Piecewise` squares, which its own discontinuity handling reads off the symbolic Hamiltonian, and the sampled #2552 pulses, whose drive QF does not read. QF's pass is not unconditional. When one step-function factor has no edge strictly inside the time range, as `Boole[t >= 0]` on [0, 100] has none, the landmark search fails and returns no restart times at all. So a square pulse from t = 0 written the way #1265 writes its windows makes QF skip a later Gaussian exactly as QuTiP does. Where QuTiP needs no extra steps it is the faster package, by more than an order of magnitude on these two- and three-level problems; at width 0.001 the advised `max_step` costs QuTiP hundreds of thousands of steps and QF becomes the faster one.

## How to run

From this folder, with a Python venv holding `qutip==5.3.1`:

```bash
python3.13 -m venv /tmp/qutip531 && /tmp/qutip531/bin/pip install qutip==5.3.1
```

```bash
/tmp/qutip531/bin/python -I shortpulse_qutip.py
```

```bash
wolframscript -file shortpulse_reference.wls
```

```bash
wolframscript -file shortpulse_qf.wls
```

```bash
python3 -I shortpulse_table.py
```

`shortpulse_qutip.py`, `shortpulse_qf.wls` and `shortpulse_reference.wls` write `shortpulse_qutip.tsv`, `shortpulse_qf.tsv` and `shortpulse_reference.tsv`; `shortpulse_table.py` prints the tables below from them. `shortpulse_data2552/` holds the six data files attached to #2552.

## Cases

1. **QuTiP documentation**, "Working with pulses" (`doc/guide/dynamics/dynamics-time.rst`, the same text at tag v5.3.1): `sigmaz() + 10 pi (0.7 < t < 0.75) sigmax()` from |1>, `<n>` at t = 1. The pulse is a rotation about Z + aX with a = 10 pi for 0.05, so the exact value is cos^2(0.05 W) + sin^2(0.05 W)/W^2 with W = sqrt(1 + a^2).
2. **Reported issues.**
   - **#1265**: the issue's code, with the three edits QuTiP 5 requires (`progress_bar` dropped, `e_ops` passed by keyword, `lambda t, arg` renamed `lambda t, args`). Two raised-cosine pi/2 pulses of width 20 about X, `mesolve` with T1 = T2 = inf (zero collapse operators), Pe at t = 99.5; exact 1. The wait between the pulses runs from 0 (the issue's `c_x`) through 23.0018 (the threshold found on QuTiP 4.3) to 40 (`c_x_bug`) and 50. QF evolves the density matrix.
   - **#1945**: the issue's code (`mesolve`, a square pulse of 2 pi 0.002 sigmax on [100, 600]), scored by the largest error of `<sz>` over each output grid against cos(2 Omega (t - 100)) inside the pulse.
   - **#2552**: the issue's code and its attached data, two sampled pi/2 pulses about X. The drive commutes with itself (its imaginary part is below 1e-14 of its real part), so the exact value of |<0|psi(T)>|^2 is cos^2 Theta with Theta the integral of the interpolated drive: QuTiP's cubic spline (Simpson's rule per knot interval, exact for a cubic) for the QuTiP rows, and the trapezoid rule for the QF rows, because `QuantumEvolve`'s default `"MergeInterpolatingFunctions" -> True` joins the two interpolations into a linear one.
   - **#2814** has no code, only a sketch of a pulse train with long waits. Five Gaussian pi/5 pulses (sigma = 0.01) centred at 10, 30, 50, 70, 90 in [0, 100] give P1 = sin^2(k pi/10) with k pulses captured, so P1 = 1 needs all five and P1 identifies k.
3. **A pi pulse at t = 50 in [0, 100]**, H = Omega(t) X/2 from |0>, as in P1 item 4 of `QF-Numeric-Inventory.md`: square (`Piecewise` and a product of two `UnitStep`), Gaussian with sigma = w/4 and sech with tau = w/4, all of area pi, at widths 1, 0.1, 0.01, 0.001. Exact P1(100) = 1.
4. **Three-level transmon**, alpha/2pi = -200 MHz, time in ns, rotating frame, H = alpha/2 a+a+aa + Omega(t)/2 (a + a+), no DRAG, from |0>. "Windowed": a Gaussian with sigma = T/4, T = 16, truncated to [gap, gap + T] with area pi over the window, ending at gap + T (the inventory's transmon probe, at three levels). "Untruncated": sigma = 4 centred at gap + 48, ending at gap + 96. Gaps of 0, 100, 1000 and 10,000 ns. |0> is stationary while the drive is off, so the exact P1 and P2 do not depend on the gap; `shortpulse_reference.wls` computes them with NDSolveValue at `WorkingPrecision` 40 and a step bound of 1/20 ns over the pulse alone, converged to 1e-19 against a run at `WorkingPrecision` 50, and SciPy's DOP853 at rtol 1e-13 agrees to 5e-15.
5. **A square pulse from t = 0, then a Gaussian**: pi/2 on [0, 1], then a Gaussian pi/2 pulse (sigma = 0.025) at t = 50; exact P1(100) = 1, and a skipped Gaussian leaves 1/2. Built to test QF's landmark search; the square is written three ways that are the same function.

**Settings.** QuTiP default: `sesolve`/`mesolve` with no options (method `adams`, atol 1e-8, rtol 1e-6, nsteps 2500, max_step 0). QuTiP with the advised step: `max_step` = 0.49 w, under the docstring's "less than half the width of the thinnest pulse"; where that exceeds the default `nsteps` and QuTiP raises `IntegratorException`, the run is repeated with `nsteps` = 2 x span / max_step and marked + in the table. The width w is the full width of a square pulse, 4 sigma for a Gaussian, 4 tau for sech, 20 for #1265, 500 for #1945, 9.95e-8 s for #2552 and 16 ns for the transmon. QF default: `QuantumEvolve` with no options. QF control: NDSolveValue at its defaults on the equations `QuantumEvolve` builds (`"ReturnEquations" -> True`) with its restart events removed; `QuantumEvolve` passes NDSolveValue no options of its own, so this is QF without the landmarks.

**Timing.** Median wall time of three runs: QuTiP's `time.perf_counter` around one solver call, QF's `EvaluationData` around `QuantumEvolve` and the read at the final time. Apple M2 Pro, 32 GB, load average about 4 from other processes. The QF time for #1945 is mostly the read of 268 states through `psi[t]`; for #2552 it is mostly the merged 3406-piece linear interpolation of the drive, evaluated at every step.

## Results

Versions:

- QuTiP probe: QuTiP 5.3.1, SciPy 1.18.1, NumPy 2.5.4, Python 3.13.1, macOS-14.2-arm64-arm-64bit-Mach-O, 2026-10-19 13:17
- QF probe: Wolfram Language 15.0.1 for Mac OS X ARM (64-bit) (July 2, 2026), QuantumFramework 2.1.3 at b92953cd, uncommitted kernel files: {M QuantumFramework/Kernel/QuantumState/Properties.m}, Arrays 1.4.1, TensorNetworks 1.2.0, 2026-10-19 13:19
- Transmon reference: Wolfram Language 15.0.1 for Mac OS X ARM (64-bit) (July 2, 2026), 2026-10-19 13:09

| Case | Exact | QF default | error | time | QuTiP default | error | time | QuTiP, max_step 0.49 w | error | time |
|---|---|---|---|---|---|---|---|---|---|---|
| QuTiP docs example, <n>(1) | 0.00101282 | 0.00101282 | 4.8e-11 | 17 ms | 1 | 1.0e+00 | 1.2 ms | 0.00101281 | 1.3e-08 | 1.4 ms |
| #1265, pi/2 + pi/2, no wait, Pe | 1 | 1 | 1.9e-08 | 26 ms | 0.999999 | 5.4e-07 | 3.1 ms | 0.999999 | 5.4e-07 | 3.2 ms |
| #1265, pi/2, wait 40, pi/2, Pe | 1 | 1 | 6.9e-08 | 27 ms | 0.499994 | 5.0e-01 | 2.6 ms | 1 | 5.9e-07 | 3.3 ms |
| #1945, 268 output times, max error of <sz> | 0 | 1.60147e-07 | 1.6e-07 | 6.1e+02 ms | 1.99993 | 2.0e+00 | 1.5 ms | 2.59226e-06 | 2.6e-06 | 1.7 ms |
| #2552 index93, sampled pulses, P0 | 5.72404e-16 | 6.6926e-14 | 6.6e-14 | 8.9e+02 ms | 0.499994 | 5.0e-01 | 18 ms | 4.75647e-11 | 4.8e-11 | 19 ms |
| #2814-style train of 5 pulses, P1 | 1 | 1 | 3.3e-12 | 21 ms | 0 | 1.0e+00 | 0.2 ms | 1 + | 4.8e-10 | 21 ms |
| pi pulse, square Piecewise, w = 1.0 | 1 | 1 | 8.9e-16 | 13 ms | 0 | 1.0e+00 | 0.2 ms | 1 | 1.2e-13 | 0.5 ms |
| pi pulse, square UnitStep, w = 1.0 | 1 | 1 | 1.4e-15 | 18 ms | 0 | 1.0e+00 | 0.1 ms | 1 | 1.2e-13 | 0.9 ms |
| pi pulse, Gaussian sigma=w/4, w = 1.0 | 1 | 1 | 5.2e-14 | 11 ms | 0 | 1.0e+00 | 0.1 ms | 1 | 2.2e-16 | 0.7 ms |
| pi pulse, sech tau=w/4, w = 1.0 | 1 | 1 | 3.0e-13 | 12 ms | 6.58121e-51 | 1.0e+00 | 0.1 ms | 1 | 7.4e-13 | 0.9 ms |
| pi pulse, square Piecewise, w = 0.1 | 1 | 1 | 1.4e-15 | 13 ms | 0 | 1.0e+00 | 0.1 ms | 1 | 1.7e-12 | 1.9 ms |
| pi pulse, square UnitStep, w = 0.1 | 1 | 1 | 1.4e-15 | 18 ms | 0 | 1.0e+00 | 0.1 ms | 1 | 1.7e-12 | 4.5 ms |
| pi pulse, Gaussian sigma=w/4, w = 0.1 | 1 | 1 | 2.4e-14 | 12 ms | 0 | 1.0e+00 | 0.1 ms | 1 | 3.0e-14 | 2.9 ms |
| pi pulse, sech tau=w/4, w = 0.1 | 1 | 1 | 3.4e-13 | 11 ms | 0 | 1.0e+00 | 0.1 ms | 1 | 2.1e-12 | 3.1 ms |
| pi pulse, square Piecewise, w = 0.01 | 1 | 1 | 1.6e-15 | 13 ms | 0 | 1.0e+00 | 0.1 ms | 1 + | 6.8e-12 | 18 ms |
| pi pulse, square UnitStep, w = 0.01 | 1 | 1 | 1.6e-15 | 18 ms | 0 | 1.0e+00 | 0.1 ms | 1 + | 6.8e-12 | 41 ms |
| pi pulse, Gaussian sigma=w/4, w = 0.01 | 1 | 1 | 3.7e-14 | 13 ms | 0 | 1.0e+00 | 0.2 ms | 1 + | 6.9e-13 | 26 ms |
| pi pulse, sech tau=w/4, w = 0.01 | 1 | 1 | 3.7e-13 | 13 ms | 0 | 1.0e+00 | 0.1 ms | 1 + | 2.1e-12 | 27 ms |
| pi pulse, square Piecewise, w = 0.001 | 1 | 1 | 1.8e-08 | 13 ms | 0 | 1.0e+00 | 0.1 ms | 1 + | 4.4e-13 | 1.6e+02 ms |
| pi pulse, square UnitStep, w = 0.001 | 1 | 1 | 1.9e-15 | 18 ms | 0 | 1.0e+00 | 0.1 ms | 1 + | 4.4e-13 | 3.9e+02 ms |
| pi pulse, Gaussian sigma=w/4, w = 0.001 | 1 | 1 | 5.2e-14 | 13 ms | 0 | 1.0e+00 | 0.1 ms | 1 + | 6.9e-13 | 2.5e+02 ms |
| pi pulse, sech tau=w/4, w = 0.001 | 1 | 1 | 3.2e-13 | 12 ms | 0 | 1.0e+00 | 0.1 ms | 1 + | 4.3e-13 | 2.6e+02 ms |
| transmon, windowed Gaussian, gap 0 ns | P1 0.987074, P2 0.000549 | P1 0.987074, P2 0.000549 | 1.4e-08 | 19 ms | P1 0.987074, P2 0.000549 | 4.8e-09 | 0.4 ms | P1 0.987074, P2 0.000549 | 4.8e-09 | 0.3 ms |
| transmon, windowed Gaussian, gap 1000 ns | P1 0.987074, P2 0.000549 | P1 0.987074, P2 0.000549 | 1.0e-08 | 23 ms | P1 0, P2 0 | 9.9e-01 | 0.1 ms | P1 0.987074, P2 0.000549 | 4.8e-08 | 0.6 ms |
| transmon, untruncated Gaussian, gap 0 ns | P1 0.989075, P2 1.83e-08 | P1 0.989075, P2 1.8e-08 | 6.9e-09 | 16 ms | P1 1.78271e-27, P2 3.54e-63 | 9.9e-01 | 0.1 ms | P1 0.989075, P2 1.81e-08 | 8.0e-08 | 1 ms |
| transmon, untruncated Gaussian, gap 1000 ns | P1 0.989075, P2 1.83e-08 | P1 0.989075, P2 1.8e-08 | 4.8e-09 | 13 ms | P1 0, P2 0 | 9.9e-01 | 0.1 ms | P1 0.989075, P2 1.77e-08 | 3.2e-08 | 1.2 ms |
| square pi/2 from t = 0 as Boole[t>=0] Boole[t<=1], Gaussian pi/2 at 50 | 1 | 0.5 | 5.0e-01 | 14 ms | 0.500043 | 5.0e-01 | 0.2 ms | 1 | 4.3e-12 | 3.3 ms |
| square pi/2 from t = 0 as UnitStep[t] UnitStep[1-t], Gaussian pi/2 at 50 | 1 | 0.5 | 5.0e-01 | 13 ms | 0.500043 | 5.0e-01 | 0.3 ms | 1 | 4.3e-12 | 5.7 ms |
| square pi/2 from t = 0 as Boole[0<=t<=1], Gaussian pi/2 at 50 | 1 | 1 | 1.1e-13 | 14 ms | 0.500043 | 5.0e-01 | 0.2 ms | 1 | 4.3e-12 | 3.3 ms |

QF with and without its restart events (NDSolveValue at its defaults on the same equations):

| Case | Restart times | QF default | error | QF, restart events removed | error |
|---|---|---|---|---|---|
| QuTiP docs example, <n>(1) | 2 | 0.00101282 | 4.8e-11 | 0.00101282 | 8.0e-11 |
| #1265, pi/2 + pi/2, no wait, Pe | 0 | 1 | 1.9e-08 | 1 | 1.9e-08 |
| #1265, pi/2, wait 40, pi/2, Pe | 0 | 1 | 6.9e-08 | 1 | 6.9e-08 |
| #1945, 268 output times, max error of <sz> | 2 | 1.60147e-07 | 1.6e-07 | 3.56772e-06 | 3.6e-06 |
| #2552 index93, sampled pulses, P0 | 0 | 6.6926e-14 | 6.6e-14 | 6.6926e-14 | 6.6e-14 |
| #2814-style train of 5 pulses, P1 | 25 | 1 | 3.3e-12 | 0 | 1.0e+00 |
| pi pulse, square Piecewise, w = 1.0 | 2 | 1 | 8.9e-16 | 1 | 8.3e-09 |
| pi pulse, square UnitStep, w = 1.0 | 2 | 1 | 1.4e-15 | 1 | 1.6e-08 |
| pi pulse, Gaussian sigma=w/4, w = 1.0 | 5 | 1 | 5.2e-14 | 5.89989e-42 | 1.0e+00 |
| pi pulse, sech tau=w/4, w = 1.0 | 5 | 1 | 3.0e-13 | 1 | 7.9e-08 |
| pi pulse, square Piecewise, w = 0.1 | 2 | 1 | 1.4e-15 | 1 | 8.3e-09 |
| pi pulse, square UnitStep, w = 0.1 | 2 | 1 | 1.4e-15 | 1 | 1.4e-08 |
| pi pulse, Gaussian sigma=w/4, w = 0.1 | 5 | 1 | 2.4e-14 | 0 | 1.0e+00 |
| pi pulse, sech tau=w/4, w = 0.1 | 5 | 1 | 3.4e-13 | 2.54648e-82 | 1.0e+00 |
| pi pulse, square Piecewise, w = 0.01 | 2 | 1 | 1.6e-15 | 0 | 1.0e+00 |
| pi pulse, square UnitStep, w = 0.01 | 2 | 1 | 1.6e-15 | 1 | 1.4e-08 |
| pi pulse, Gaussian sigma=w/4, w = 0.01 | 5 | 1 | 3.7e-14 | 3.24751e-23 | 1.0e+00 |
| pi pulse, sech tau=w/4, w = 0.01 | 5 | 1 | 3.7e-13 | 0 | 1.0e+00 |
| pi pulse, square Piecewise, w = 0.001 | 2 | 1 | 1.8e-08 | 0 | 1.0e+00 |
| pi pulse, square UnitStep, w = 0.001 | 2 | 1 | 1.9e-15 | 1 | 1.4e-08 |
| pi pulse, Gaussian sigma=w/4, w = 0.001 | 5 | 1 | 5.2e-14 | 0 | 1.0e+00 |
| pi pulse, sech tau=w/4, w = 0.001 | 5 | 1 | 3.2e-13 | 7.5274e-63 | 1.0e+00 |
| transmon, windowed Gaussian, gap 0 ns | 0 | P1 0.987074, P2 0.000549 | 1.4e-08 | P1 0.987074, P2 0.000549 | 3.0e-08 |
| transmon, windowed Gaussian, gap 1000 ns | 4 | P1 0.987074, P2 0.000549 | 1.0e-08 | P1 0, P2 0 | 9.9e-01 |
| transmon, untruncated Gaussian, gap 0 ns | 5 | P1 0.989075, P2 1.8e-08 | 6.9e-09 | P1 0.989075, P2 1.81e-08 | 3.5e-08 |
| transmon, untruncated Gaussian, gap 1000 ns | 5 | P1 0.989075, P2 1.8e-08 | 4.8e-09 | P1 2.12914e-61, P2 1.2e-160 | 9.9e-01 |
| square pi/2 from t = 0 as Boole[t>=0] Boole[t<=1], Gaussian pi/2 at 50 | 0 | 0.5 | 5.0e-01 | 0.5 | 5.0e-01 |
| square pi/2 from t = 0 as UnitStep[t] UnitStep[1-t], Gaussian pi/2 at 50 | 0 | 0.5 | 5.0e-01 | 0.5 | 5.0e-01 |
| square pi/2 from t = 0 as Boole[0<=t<=1], Gaussian pi/2 at 50 | 6 | 1 | 1.1e-13 | 0.5 | 5.0e-01 |
| #1265, wait 10 | 0 | 1 | 6.5e-08 | 1 | 6.5e-08 |
| #1265, wait 20 | 0 | 1 | 7.3e-09 | 1 | 7.3e-09 |
| #1265, wait 23 | 0 | 1 | 7.1e-08 | 1 | 7.1e-08 |
| #1265, wait 23.0017 | 0 | 1 | 6.7e-08 | 1 | 6.7e-08 |
| #1265, wait 23.0018 (BoxiLi threshold) | 0 | 1 | 6.7e-08 | 1 | 6.7e-08 |
| #1265, wait 24 | 0 | 1 | 7.0e-08 | 1 | 7.0e-08 |
| #1265, wait 30 | 0 | 1 | 6.9e-08 | 1 | 6.9e-08 |
| #1265, wait 50 | 0 | 1 | 6.9e-08 | 1 | 6.9e-08 |
| transmon, windowed Gaussian, gap 100 ns | 4 | P1 0.987074, P2 0.000549 | 1.0e-08 | P1 0.987074, P2 0.000549 | 3.0e-08 |
| transmon, windowed Gaussian, gap 10000 ns | 4 | P1 0.987074, P2 0.000549 | 2.8e-08 | P1 0, P2 0 | 9.9e-01 |
| transmon, untruncated Gaussian, gap 100 ns | 5 | P1 0.989075, P2 1.8e-08 | 4.7e-09 | P1 0.989075, P2 1.8e-08 | 5.2e-08 |
| transmon, untruncated Gaussian, gap 10000 ns | 5 | P1 0.989075, P2 1.8e-08 | 3.7e-09 | P1 1.80668e-59, P2 0 | 9.9e-01 |

#1265, Pe at t = 99.5 against the wait between the two pi/2 pulses (exact 1):

| Wait | QF default | error | QuTiP default | error | QuTiP, max_step 9.8 | error |
|---|---|---|---|---|---|---|
| 0 (c_x) | 1 | 1.9e-08 | 0.999999 | 5.4e-07 | 0.999999 | 5.4e-07 |
| 10 | 1 | 6.5e-08 | 1 | 2.7e-07 | 0.999999 | 1.4e-06 |
| 20 | 1 | 7.3e-09 | 0.499994 | 5.0e-01 | 0.999999 | 1.4e-06 |
| 23 | 1 | 7.1e-08 | 0.499994 | 5.0e-01 | 1 | 2.7e-07 |
| 23.0017 | 1 | 6.7e-08 | 0.499994 | 5.0e-01 | 1 | 3.0e-07 |
| 23.0018 (BoxiLi threshold) | 1 | 6.7e-08 | 0.499994 | 5.0e-01 | 1 | 3.0e-07 |
| 24 | 1 | 7.0e-08 | 0.499994 | 5.0e-01 | 0.999998 | 2.0e-06 |
| 30 | 1 | 6.9e-08 | 0.499994 | 5.0e-01 | 1 | 1.5e-06 |
| 40 (c_x_bug) | 1 | 6.9e-08 | 0.499994 | 5.0e-01 | 1 | 5.9e-07 |
| 50 | 1 | 6.9e-08 | 0.499994 | 5.0e-01 | 1 | 2.1e-07 |

Other rows:

| Case | Package, setting | Result | Error | Time | Note |
|---|---|---|---|---|---|
| 1 docs, pulse 0.7<t<0.75 | QuTiP 5.3.1, max_step 0.01 (docs example) | 0.00101282 | 6.5e-10 | 1.4 ms |  |
| 2b #1945, tlist linspace(0,700,267) | QuTiP 5.3.1, default | 1.01847e-05 | 1.0e-05 | 1.8 ms | FutureWarning: The signature f(t, args) is deprecated and will be removed in QuTiP 5. |
| 2b #1945, tlist linspace(0,700,267) | QuTiP 5.3.1, max_step 245 | 1.01847e-05 | 1.0e-05 | 1.7 ms | FutureWarning: The signature f(t, args) is deprecated and will be removed in QuTiP 5. |
| 2b #1945, tlist linspace(0,700,800) | QuTiP 5.3.1, default | 2.17971e-06 | 2.2e-06 | 4.5 ms | FutureWarning: The signature f(t, args) is deprecated and will be removed in QuTiP 5. |
| 2b #1945, tlist linspace(0,700,800) | QuTiP 5.3.1, max_step 245 | 9.76273e-06 | 9.8e-06 | 4.4 ms | FutureWarning: The signature f(t, args) is deprecated and will be removed in QuTiP 5. |
| 2c #2552, index92 | QuTiP 5.3.1, default | 2.26754e-15 | 3.4e-14 | 19 ms | Theta - pi/2 = -1.914e-07 |
| 2c #2552, index92 | QuTiP 5.3.1, max_step 4.875e-08 | 4.26145e-11 | 4.3e-11 | 18 ms | Theta - pi/2 = -1.914e-07 |
| 2c #2552, index94 | QuTiP 5.3.1, default | 7.58926e-15 | 2.9e-14 | 19 ms | Theta - pi/2 = -1.914e-07 |
| 2c #2552, index94 | QuTiP 5.3.1, max_step 4.875e-08 | 3.41353e-11 | 3.4e-11 | 19 ms | Theta - pi/2 = -1.914e-07 |
| 3 pi pulse, square Piecewise, w 1.0 | QuTiP 5.3.1, default, tlist linspace(0,100,201) | 0 | 1.0e+00 | 1.2 ms |  |
| 3 pi pulse, square UnitStep, w 1.0 | QuTiP 5.3.1, default, tlist linspace(0,100,201) | 0 | 1.0e+00 | 1.3 ms |  |
| 3 pi pulse, Gaussian sigma=w/4, w 1.0 | QuTiP 5.3.1, default, tlist linspace(0,100,201) | 0 | 1.0e+00 | 1.1 ms |  |
| 3 pi pulse, sech tau=w/4, w 1.0 | QuTiP 5.3.1, default, tlist linspace(0,100,201) | 1.37601e-104 | 1.0e+00 | 1.2 ms |  |
| 3 pi pulse, square Piecewise, w 0.1 | QuTiP 5.3.1, default, tlist linspace(0,100,201) | 0 | 1.0e+00 | 1.1 ms |  |
| 3 pi pulse, square UnitStep, w 0.1 | QuTiP 5.3.1, default, tlist linspace(0,100,201) | 0 | 1.0e+00 | 1.2 ms |  |
| 3 pi pulse, Gaussian sigma=w/4, w 0.1 | QuTiP 5.3.1, default, tlist linspace(0,100,201) | 0 | 1.0e+00 | 1.1 ms |  |
| 3 pi pulse, sech tau=w/4, w 0.1 | QuTiP 5.3.1, default, tlist linspace(0,100,201) | 0 | 1.0e+00 | 1.1 ms |  |
| 3 pi pulse, square Piecewise, w 0.01 | QuTiP 5.3.1, default, tlist linspace(0,100,201) | 0 | 1.0e+00 | 1.1 ms |  |
| 3 pi pulse, square UnitStep, w 0.01 | QuTiP 5.3.1, default, tlist linspace(0,100,201) | 0 | 1.0e+00 | 1.1 ms |  |
| 3 pi pulse, Gaussian sigma=w/4, w 0.01 | QuTiP 5.3.1, default, tlist linspace(0,100,201) | 0 | 1.0e+00 | 1.1 ms |  |
| 3 pi pulse, sech tau=w/4, w 0.01 | QuTiP 5.3.1, default, tlist linspace(0,100,201) | 0 | 1.0e+00 | 1.1 ms |  |
| 3 pi pulse, square Piecewise, w 0.001 | QuTiP 5.3.1, default, tlist linspace(0,100,201) | 0 | 1.0e+00 | 1.1 ms |  |
| 3 pi pulse, square UnitStep, w 0.001 | QuTiP 5.3.1, default, tlist linspace(0,100,201) | 0 | 1.0e+00 | 1.2 ms |  |
| 3 pi pulse, Gaussian sigma=w/4, w 0.001 | QuTiP 5.3.1, default, tlist linspace(0,100,201) | 0 | 1.0e+00 | 1.2 ms |  |
| 3 pi pulse, sech tau=w/4, w 0.001 | QuTiP 5.3.1, default, tlist linspace(0,100,201) | 0 | 1.0e+00 | 1.2 ms |  |
| 4 transmon, windowed | SciPy DOP853 reference, rtol 1e-13, max_step 0.01, pulse only | P1 0.987073888304, P2 0.000549442 | 6.7e-16 |  |  |
| 4 transmon, untruncated | SciPy DOP853 reference, rtol 1e-13, max_step 0.01, pulse only | P1 0.989075103923, P2 1.83134e-08 | 4.6e-15 |  |  |
| 2a #1265, wait 40 (c_x_bug) | QF, diagnostic: MaxStepFraction -> 1 | 1 | 6.6e-08 | 25 ms | Catenate::invrp; Join::heads |
| 2b #1945, tlist linspace(0,700,267) | QF, default | 1.74306e-07 | 1.7e-07 | 6e+02 ms | restart times 2 |
| 2b #1945, tlist linspace(0,700,800) | QF, default | 1.73237e-07 | 1.7e-07 | 1.89 s | restart times 2 |
| 2b #1945, tlist linspace(0,700,268) | QF, default, read via psi["StateVector"] /. t -> time | 4.0435e+08 | 4.0e+08 | 16 ms |  |
| 4 transmon, windowed, gap 0 | QF, diagnostic: MaxStepFraction -> 1 | P1 0.987074, P2 0.000549 | 1.4e-08 | 18 ms | Catenate::invrp; Join::heads |
| 2c #2552, index92 | QF, default | 2.81681e-14 | 1.8e-14 | 9.5e+02 ms | restart times 0 |
| 2c #2552, index93 | QF, diagnostic: MaxStepFraction -> 1 | 2.36902e-13 | 2.4e-13 | 9.2e+02 ms |  |
| 2c #2552, index94 | QF, default | 7.20841e-15 | 3.9e-14 | 9.2e+02 ms | restart times 0 |

## QF defects this run found

1. **The landmark search returns nothing when one step-function factor has no edge inside the time range.** `pulseLandmarks` (`QuantumEvolution.m`) takes `t /. Solve[FunctionDiscontinuities[f, t] && t0 < t < t1, t, Reals]` for each step function f. For `Boole[t >= 0]` on [0, 100] the solution list is empty, `t /. {}` leaves the bare `t`, `Catenate` fails (`Catenate::invrp`, then `Join::heads`), and every landmark is lost, Gaussian ones included. Case 5 shows the consequence: written as `Boole[t >= 0] Boole[t <= 1]` or `UnitStep[t] UnitStep[1 - t]`, the square pulse from t = 0 hides the later Gaussian and QF returns P1 = 0.5; written as `Boole[0 <= t <= 1]` it returns 1. #1265 (whose windows are such products) and the windowed transmon at gap 0 lose their landmarks the same way but pass, the first through NDSolve's own handling of `Boole` and the second because the drive starts at t = 0. A fix: read the solutions with `Cases[Solve[...], {t -> x_} :> x]`, which gives `{}` for an empty solution.
2. **Reading an evolved state before its time is bound is still wrong near a restart.** For #1945's pulse on [100, 600] in [0, 700], `Normal[psi["StateVector"]] /. t -> time` departs from normalization by up to about 5e8 between roughly t = 588 and t = 600, while `psi[t]` stays normalized to 2e-16. The solution grid holds the restart a hair before t = 600 next to NDSolve's own grid point at 600. Inventory item 15 records this read as resolved by Arrays 1.4.1; it is not resolved for this pulse.

## Proposed edits to `QF-Numeric-Inventory.md` (not applied)

Replace P1 item 4 with:

> 4. **Short pulses**, the most reported failure across QuTiP users (issue #2814), resolved on 5 October 2026 for pulses written as formulas. A numeric `QuantumEvolve` reads landmarks off the symbolic Hamiltonian and restarts its integration at each one: every discontinuity, and two and four widths either side of the centre of every Gaussian or exponential profile, which covers square, Gaussian, sech and tanh pulses. A pi pulse at t = 50 in [0, 100], with no options set, is captured at every width and shape below; "skipped without restarts" marks the pulses that NDSolve at its defaults skips on the same equations once the restarts are removed. QuTiP 5.3.1 at its defaults skips all sixteen (probes `shortpulse_qf.wls` and `shortpulse_qutip.py`, 19 October 2026).
>
>    | Width | Square `Piecewise` | Square `UnitStep` | Gaussian | sech |
>    |---|---|---|---|---|
>    | 1.0 | captured | captured | captured (skipped without restarts) | captured |
>    | 0.1 | captured | captured | captured (skipped without restarts) | captured (skipped without restarts) |
>    | 0.01 | captured (skipped without restarts) | captured | captured (skipped without restarts) | captured (skipped without restarts) |
>    | 0.001 | captured (skipped without restarts) | captured | captured (skipped without restarts) | captured (skipped without restarts) |
>
>    Each solve takes about 0.01 to 0.02 s; a train of 50 Gaussian pi/2 pulses of width 0.04 lands on the exact population to 2e-10 in 0.07 s. On the cases QuTiP users reported (the example in QuTiP's documentation, #1265, #1945, #2552 and a #2814-style train) and on a three-level transmon pi pulse after an idle gap, QF at its defaults lands within its default tolerance of the exact value wherever QuTiP 5.3.1 at its defaults returns the state the pulse never touched ([head-to-head](shortpulse-head-to-head.md)). A smooth drive has no landmarks and runs as before. A pulse given as data (an `InterpolatingFunction`) is not read: the sampled pulses of #2552 are captured only because NDSolve's step growth does not outrun them, so `MaxStepSize` remains the control there.
>
>    Open: the landmark search returns no landmarks at all when one step-function factor has no edge strictly inside the time range. An empty solution leaves the bare time symbol where `Catenate` expects a list (`Catenate::invrp`, `Join::heads`). A pi/2 square pulse from t = 0 written as `Boole[t >= 0] Boole[t <= 1]`, the way #1265 writes its windows, followed by a Gaussian pi/2 pulse at t = 50, gives P1 = 0.5 instead of 1; written as `Boole[0 <= t <= 1]` it gives 1.

Replace P2 item 15 with:

> 15. An evolved state read before its time is bound (`psi["StateVector"]` rather than `psi[t]`) failed with `Interpolation::inddp` when NDSolve repeats a grid point, which it does at a discontinuity of the Hamiltonian, because the Wolfram/Arrays lazy container rebuilt its interpolations over the whole grid. Arrays 1.4.1, which QF now requires, removes the failure, but the read can still be wrong beside a restart. For a square pulse on [100, 600] in [0, 700], the grid holds the restart a hair before t = 600 next to NDSolve's own point at 600, and the read departs from normalization by up to about 5e8 over the last dozen time units of the pulse, while `psi[t]` stays normalized to machine precision (probe `shortpulse_qf.wls`, 19 October 2026). Read evolved states through `psi[t]`.
