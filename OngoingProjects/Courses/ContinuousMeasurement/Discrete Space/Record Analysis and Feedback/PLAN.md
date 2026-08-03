# Plan: Record Analysis and Feedback on the Continuous-Measurement SME Integrator

Scope of this document: a plan (not the work) for "item 4" of the SDE audit. We have a
verified diffusive-SME integrator and its readout; everything downstream of "simulate one
trajectory" is unbuilt. This plan says what to build on top of the two functions we have,
the questions worth settling, and exactly where the three reference texts treat each piece.

The work itself is to be done in a separate session, driven by this plan.

## 0. The gap (item 4, verbatim from the audit)

> Record analysis and feedback: the readout is `Sow`-ed but nothing is built on it, no
> spectrum / estimation / smoothing, and no real-time feedback loop, which is the very
> setting the scheme was designed for.

The integrator turns a detector record into a conditioned state (it *is* the filter), and
it Sow-s the record. What it does not do: analyze the record (spectrum, parameter
estimation, smoothing) or close a loop (drive the Hamiltonian from the record or the
filtered state). Rouchon and Ralph built the scheme for exactly the last item.

## 1. What we have: the two functions everything is built on

Both live in the essay `../Manual StepLike Ito/continuous-measurement-trajectories.md`
(and QF-free in `../Manual StepLike Ito/manual-implementation-ito-qf-independent.wl` as
`\[ScriptCapitalT]` / `\[ScriptCapitalR]`). They are plain dense-matrix functions, no
QuantumFramework dependency.

**`\[ScriptCapitalR]` (readout increments of one step)** returns the record vector `dy`,
one entry per monitored operator:

```wl
\[ScriptCapitalR][ρ_, Ls_, dt_, dw_, η_] :=
   MapThread[Sqrt[#3] Tr[(#1 + ConjugateTranspose[#1]) . ρ] dt + #2 &, {Ls, dw, η}];
```

i.e. `dy_r = Sqrt[η_r] Tr[(L_r + L_r^†) ρ] dt + dW_r`, signal (the mean) plus shot noise.

**`\[ScriptCapitalT]` (the trajectory map)** folds the positivity-preserving update over the
noise and `Sow`s the record each step:

```wl
\[ScriptCapitalT][ρ0_, H_, Ls_, η:(None|_List):None, V:(None|_List):None, dt_, tf_]
```

returns the list of conditioned density matrices `{ρ0, ρ1, ...}`, and `Reap` around the
call recovers `{trajectory, {records}}` where `records[[n]]` is the step-`n` `dy`.

Everything item 4 needs consumes exactly these two outputs: the **record** (`dy` stream)
and the ability to **re-run `\[ScriptCapitalT]` with a modified `H`** (for feedback). No
new physics is required, only tooling around them, plus one small extension (a filter mode
of `\[ScriptCapitalT]` that ingests an externally supplied record instead of drawing its
own noise).

## 2. How to address item 4: the five pieces to build

Each piece names the new WL surface, what it consumes from `\[ScriptCapitalT]`/`\[ScriptCapitalR]`,
and its reference anchor. Pieces 2.1 and 2.5(a) are the highest value.

### 2.1 Power spectrum of the record, `S(ω)`
The Fourier transform of the record's two-point correlation (Wiener-Khinchin), or the
Welch-averaged periodogram of the current. Already prototyped in the essay's Part VI
(`powerSpectrum`, verified). What it reveals: a flat shot-noise floor with peaks at the
system's transition frequencies (Rabi peak, Mollow triplet), peak widths giving decay
rates, and sub-shot-noise dips signaling squeezing.
- **Consumes:** the `records` stream from one long `\[ScriptCapitalT]` run.
- **Refs:** Jacobs §3.1.4, §3.6, §3.7.2, App F.8 (the record spectrum, incl. linear/Heisenberg
  route); Barchielli-Gregoratti §4.5 (rigorous spectrum of the output process) and §9.2
  (two-level-atom homodyne spectrum, §9.2.1 the homodyne-current spectrum, §9.2.2 squeezing);
  Wiseman-Milburn Ch 4 (homodyne detection).

### 2.2 State estimation (filtering)
Running `\[ScriptCapitalT]` on a record *is* the optimal (quantum Kalman-type) filter: it
returns the observer's conditional state. The one missing capability is a **filter mode**
that ingests an *external* record (a real or independently simulated `dR`) rather than
drawing its own `dW`. Small change: invert `dR = Sqrt[η] <L+L^†> dt + dW` for `dW` at each
step, or accept `dW` directly, and fold.
- **Consumes:** an external `dR`/`dW` stream; returns the conditional-state trajectory.
- **Refs:** Jacobs §3.1.3 (in the linear-Gaussian case the SME *is* the classical
  Kalman-Bucy filter), §3.2; Barchielli-Gregoratti §2.4, §3.4, §5 (a posteriori states);
  Wiseman-Milburn Ch 4.

### 2.3 Parameter estimation
Infer model parameters (Rabi `Ω`, rate `γ`, efficiency `η`) from the record. Two routes:
(a) read them off the spectrum (peak position = `Ω`, width = `γ`); (b) likelihood/Bayesian
estimation over the record, i.e. the hybrid master equation carrying an unknown parameter.
- **Consumes:** the record; and, for (b), a `\[ScriptCapitalT]` variant carrying a
  parameter distribution.
- **Refs:** Jacobs §3.8 (parameter estimation via the hybrid master equation), Ch 6
  (metrology, Cramér-Rao); Wiseman-Milburn Ch 2 (quantum parameter estimation).

### 2.4 Smoothing (retrodiction)
Filtering uses only the past record; the *future* record further constrains the state at an
earlier time. A smoother combines past and future for a strictly sharper estimate. Build a
forward-backward pass on top of `\[ScriptCapitalT]`.
- **Consumes:** the full record, forward and backward.
- **Refs:** Barchielli-Gregoratti §5.2 (a posteriori states, purification) for the
  a-posteriori structure; primary sources for quantum smoothing are Tsang (2009) and
  Guevara-Wiseman (2015), to be cited by the working session; Jacobs / Wiseman-Milburn for
  the filtering baseline.

### 2.5 Feedback control (close the loop)
The payoff. Drive `H(t)` from the measurement in real time. Two modes:

(a) **Markovian (current-based) feedback**, Wiseman-Milburn: add a feedback Hamiltonian
proportional to the instantaneous current, `H_fb ∝ dR/dt`, giving a modified (feedback)
SME. Simple, no state estimate needed; the canonical route for spin-squeezing and
stabilization. Requires checking that the positivity-preserving step survives the extra
term.
- **Refs:** Wiseman-Milburn §5.5 (homodyne-mediated feedback), §5.6 (Markovian feedback in
  a linear system), §5.7 (deterministic spin-squeezing), App H.3 (derivation of the
  Markovian feedback SME); Jacobs §5.3, §5.4; Barchielli-Gregoratti §10.2 (the WM scheme),
  §10.5 (squeezing of fluorescence), §10.6 (line narrowing).

(b) **State-based (Bayesian) feedback**: compute the control from the *filtered* conditional
state (LQG-optimal for linear systems). More powerful, needs the filter (2.2).
- **Refs:** Wiseman-Milburn Ch 6 (state-based feedback: §6.2 freezing a conditional state,
  §6.6 linear quantum systems = LQG); Jacobs §5.4.1 (rapid purification), §5.4.3
  (near-optimal single-qubit feedback), §5.5 (HJB / optimal control for linear systems).

Canonical demonstrations to reproduce: rapid purification (Jacobs §5.4.1), deterministic
spin-squeezing (WM §5.7), freezing / stabilizing a target state (WM §6.2), line narrowing
and squeezing of fluorescence (Barchielli-Gregoratti §10.5-10.6).

## 3. The most important questions to discuss

1. **Observable vs inferred.** The record is the only directly measured object; the state
   is the filter output. Every downstream tool is either record statistics (spectrum,
   parameter estimation) or an estimate (filtering, smoothing). Keep this line sharp.
2. **Filtering vs smoothing.** What does the future record buy over the past-only estimate,
   and when is smoothing worth its non-causal cost (post-processing, not real-time control)?
3. **Markovian vs state-based feedback.** Current-driven (simple, no estimator, WM) versus
   state-estimate-driven (LQG-optimal, needs the filter). When is each the right choice,
   and what is the performance gap?
4. **Does feedback preserve structure?** The bare step is positivity-preserving at any `dt`.
   Does the Markovian feedback term keep the numerator a sum of `A ρ A^†`? (See the WM
   feedback SME, App H.3.) If not, what is the safe step size?
5. **Inefficiency `η < 1`.** A lossy detector makes the conditional state mixed and degrades
   estimation and feedback. Quantify how spectrum contrast, filter sharpness, and control
   performance scale with `η`.
6. **The linear-Gaussian special case.** For quadratic `H` + linear coupling the SME reduces
   to the Kalman-Bucy filter (Jacobs §3.1.3) and optimal control is LQG (WM §6.6): the
   filter closes in a few moments, `O(1)` in Hilbert dimension. Where does closure hold and
   where does it fail (nonlinear `H`, non-Gaussian state)?
7. **Delay and noise in a real loop.** A physical feedback loop has finite delay and a noisy
   current. What are the stability limits, and how do they enter the discrete step?
8. **What the spectrum can and cannot show.** Transition frequencies, linewidths, squeezing
   (sub-shot-noise) yes; but weak measurement gives modest peak-to-floor contrast. When is
   the time-domain estimate or a direct fit better?

## 4. What each reference text says (precise pointers)

**Wiseman & Milburn, *Quantum Measurement and Control* (CUP, 2010)**: the feedback text.
- Ch 4 Quantum trajectories: §4.3 photodetection, §4.4 homodyne detection, §4.5 heterodyne,
  §4.8 imperfect detection (`η<1`).
- Ch 5 Quantum feedback control: §5.4 feedback control of a monitored system, §5.5
  homodyne-mediated feedback control, §5.6 Markovian feedback in a linear system, §5.7
  deterministic spin-squeezing.
- Ch 6 State-based quantum feedback control: §6.2 freezing a conditional state, §6.4 linear
  classical systems, §6.6 linear quantum systems (LQG).
- App H.3 Derivation of the Wiseman-Milburn Markovian feedback SME.

**Jacobs, *Quantum Measurement Theory and its Applications* (CUP, 2014)**: the most
practical/computational for item 4.
- Ch 3 Continuous measurement: §3.1.3 when the SME is the classical Kalman-Bucy filter,
  §3.1.4 the power spectrum of the measurement record, §3.6 inputs/outputs/spectra, §3.7.2
  calculating the power spectrum for linear systems, §3.8 parameter estimation: the hybrid
  master equation.
- Ch 5 Quantum feedback control: §5.3 explicit continuous-time feedback, §5.4 feedback via
  continuous measurement (§5.4.1 rapid purification, §5.4.2 control via back-action, §5.4.3
  near-optimal single-qubit feedback), §5.5 optimization (HJB; §5.5.2 optimal control for
  linear quantum systems = LQG).
- Ch 6 Metrology (§6.1 Cramér-Rao). App F.8 spectrum of the measurement signal.

**Barchielli & Gregoratti, *Quantum Trajectories and Measurements in Continuous Time: The
Diffusive Case* (Springer LNP 782, 2009)**: the rigorous mathematical foundation.
- Ch 4 Continuous measurements and instruments: §4.4 classical post-measurement processing,
  §4.5 autocorrelation and spectrum of the output process (§4.5.3 the spectrum of the output
  of the continuous measurement).
- Ch 5 SME Part II: §5.2 purification of a posteriori states (basis for smoothing/estimation).
- Ch 6 Mutual entropies and information gain in quantum continuous measurements (how much the
  record tells you about the state, an information-theoretic bound).
- Ch 9 A two-level atom, heterodyne and homodyne spectra: §9.1 heterodyne detection and the
  Mollow spectrum, §9.2 homodyne detection (§9.2.1 the spectrum of the homodyne current,
  §9.2.2 examples of spectra and squeezing).
- Ch 10 Feedback: §10.2 the feedback scheme of Wiseman and Milburn, §10.4 control of the
  atomic state, §10.5 control of the squeezing of fluorescence light, §10.6 control and line
  narrowing.

Division of labor across the three: Wiseman-Milburn for the feedback constructions,
Jacobs for the computational spectrum/estimation/optimization recipes and the Kalman-Bucy
connection, Barchielli-Gregoratti for the rigorous output-spectrum theory and worked
two-level-atom spectra and feedback.

## 5. Suggested build order

1. **Spectrum module** (extend the essay's `powerSpectrum`): `S(ω)`, squeezing check,
   Mollow triplet on a strongly driven atom. [prototype exists in Part VI]
2. **Filter mode of `\[ScriptCapitalT]`** ingesting an external record; verify it returns the
   same conditional state the generating run did.
3. **Parameter estimation** demo: recover `Ω`, `γ`, and `η` from a record (spectrum peak +
   width; likelihood for `η`).
4. **Markovian feedback** (WM): add `H_fb ∝` current to the step; reproduce spin-squeezing
   or state stabilization; check positivity of the feedback step.
5. **State-based / LQG feedback** in the linear-Gaussian case: Gaussian moment filter +
   optimal gain; freeze a conditional state.
6. **(Stretch) Smoothing**: forward-backward pass; quantify the gain over filtering.

## 6. Deliverable shape

A learning-by-computing essay + notebook (same voice and md2nb pipeline as
`continuous-measurement-trajectories.md`), QF-free, built entirely on `\[ScriptCapitalT]` /
`\[ScriptCapitalR]`. Each concept: short prose, then a cell that computes it, then a check.
Verify every cell with `wolframscript`, then `/wl-verify` and `/wl-quality`. The three books
above are the sources; cite specific sections, and add the primary papers the working
session finds (Wiseman-Milburn 1993/1994 feedback, Tsang/Guevara-Wiseman smoothing,
Doherty-Jacobs state-estimation feedback).
