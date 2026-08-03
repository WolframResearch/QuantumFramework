# Physics & Symbolic Brief: record-analysis-and-feedback.md

The design brief the essay was written against, and the map from each required
regime to the cells that discharge it. (Meta-artifact, not part of the
reader-facing notebook.)

## 1. General / symbolic object

Dimension-general plain dense complex matrices, no QuantumFramework: a `d x d`
conditional state `ρ`, monitored operators `L_r`, unmonitored `V_j`, detector
efficiencies `η_r ∈ (0,1]`, and a Hermitian feedback operator `F`. What stays
symbolic: the Hilbert dimension `d` (qubit, spin-J, oscillator), the number of
channels, the efficiency `η`, and the drive/rate/frequency `(Ω, γ, ω)` in the
closed forms (steady state, Mollow eigenvalues, squeeze minimum, DARE variances).
The single positivity-preserving step `rouchonStep` is representation-independent;
the feedback identity `c -> c - iF` is operator-level.

## 2. Doc grounding

Foundation `\[ScriptCapitalT]` / `\[ScriptCapitalR]`: the audited QF-free
`../Manual StepLike Ito/manual-implementation-ito-qf-independent.wl`. Physics
anchors: Wiseman-Milburn *Quantum Measurement and Control* (feedback ME
Eqs. 5.142-5.143, spin-squeezing 5.204-5.211, homodyne spectrum Ch 4); Jacobs
*Quantum Measurement Theory and its Applications* (Kalman-Bucy §3.1.3, record
spectrum §3.1.4/§3.7.2, hybrid ME §3.8, qubit feedback §5.4, Wineland witness);
Barchielli-Gregoratti *Quantum Trajectories ... The Diffusive Case* (rigorous
output spectrum §4.5, two-level-atom homodyne spectrum and squeezing §9.2).
Non-core WL symbols checked in a kernel before use: `Fourier`/`FourierParameters`,
`LinearSolve`/`NullSpace` (resolvent, stationary state), `MatrixExp`, `Minimize`,
`Eigenvalues`, `ResourceFunction["LindbladSolve"]`, `PositiveSemidefiniteMatrixQ`.

## 3. Physics invariants preserved and checked

Unit trace, Hermiticity, and positive-semidefiniteness of every conditioned state
at every step (the Kraus-sum guarantee, preserved by unitary feedback); ensemble
mean of trajectories = the (feedback-modified) Lindblad state; shot-noise floor
normalized to 1, with squeezing = spectrum below 1; `ℒ ρ_ss = 0`, `Tr ρ_ss = 1`,
and a unique steady state (guarded); `ξ² = 1` at the standard quantum limit for a
coherent spin state; the filter round-trip is exact on its own record.

## 4. Idiomatic plan

Functional WL: one shared captured-`Function` engine (`rouchonPieces` +
`rouchonStep`) driven by `FoldList`; `Map`/`MapThread`/`Reap`/`Sow`; `With`/`Block`
over `Module`; pattern dispatch. No `Do`/`For`/`While`/`AppendTo`/`Module`.
Symbolic / closed-form first (stationary state, Mollow eigenvalues, exact squeeze
rational, DARE variances) before Monte-Carlo estimates, with every numeric estimate
checked against a closed form.

## 5. Required regimes and limits (mapped to cells)

The trivial base case refused throughout: a single undriven qubit, unit-efficiency
`Z`-measurement, no feedback (the companion essay's starting point).

| Regime | Where it is discharged |
|---|---|
| **General symbolic case** | Dimension-general `rouchonPieces`/`rouchonStep`/`\[ScriptCapitalT]`; homodyne spectrum as the resolvent `(iω - ℒ)^{-1}` for any Lindblad generator; symbolic steady state `Ω²/(γ²+2Ω²)` (Part II); symbolic Mollow eigenvalues `{0, -γ/2, (-3γ ± √(γ²-16Ω²))/4}` (Part II); DARE filter/smoother variances symbolic in `(A,Q,R)` (Part VIII); the feedback identity `c -> c - iF` operator-level (Part VI). |
| **Exactly solvable special case** | Closed-form homodyne spectrum vs Welch periodogram (Part II); exact squeeze minimum `147779/205379` vs Barchielli-Gregoratti 0.7195 (Part III); feedback-ensemble mean vs `LindbladSolve` of the WM feedback ME (Part VI); filter round-trip exact to machine precision (Part IV); exact DARE reduction 46.25% vs the single-run smoother (Part VIII). |
| **Asymptotic / limiting case** | `η -> 1` pure conditional state and `η -> 0` squeezing vanishing (Part III/VII); no-feedback limit `F -> 0` recovers the bare SME (Part VI); strong-drive Mollow triplet, sidebands `-3γ/4 ± iΩ'` (Part II); large-record Fisher-information sharpening toward Cramér-Rao (Part V); linear-Gaussian Kalman-Bucy collapse to `O(1)` moments (Part VII). |
| **Numerical reference case** | Periodogram vs resolvent (Part II); pure-Wiener-record control for the shot-noise floor (Part II); Monte-Carlo ensemble vs `LindbladSolve` (Parts VI, VII); single-run smoother vs the exact DARE (Part VIII). |
| **Failure / edge case** | `η < 1` mixed conditional state degrades filter/spectrum/squeezing/feedback (Parts III, IV, VI); the degenerate (pure-dephasing) Liouvillian returns `Missing["NotUnique", 2]` from the uniqueness guard (Part II); a Kerr nonlinearity breaks Gaussian closure, the conditional variance spread jumping two orders of magnitude (Part VII); the feedback rotation contracts the mean spin `⟨Jx⟩`, so the honest Wineland witness `ξ²_R` exceeds the variance parameter (Part VI). |
