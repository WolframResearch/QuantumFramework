# Physics & Symbolic Brief: Smooth-Noise-NDSolve-Monitored-Particle

Committed before the essay was written; the essay's cells are graded against this.

## 1. General / symbolic object

The Stratonovich conversion and the Lindblad identity are certified on a symbolic grid
(three points x1, x2, x3, symbolic strength lam, symbolic hbar), the diagonality of the
measurement operator carrying the result to any n. The covariance references keep lam, m,
hbar, the trap frequency, the starting width, and time symbolic: the steady triple from
the fixed-point equations, the full flow in closed form via the linearized matrix Riccati
equation (MatrixExp of the 4x4 block matrix, fractional-linear map), the unconditional
moment chain from DSolve. Numbers enter only as chosen fixture inputs.

## 2. Doc grounding

NDSolve MethodOfLines / TensorProductGrid ("Coordinates", "DifferenceOrder"), the
matched-endpoint periodic idiom, and PeriodicBoundaryCondition's finite-element routing
were probed in a live kernel this session (step0/step0b/step0c probe scripts), not
recalled; StratonovichProcess/ItoProcess usage mirrors the verified conversion cell of
the qubit sibling essay; Interpolation, ListInterpolation periodic endpoint convention,
RandomChoice weighted form, DSolveValue, LaunchKernels/DistributeDefinitions probed
likewise. The Wolfram MCP was not used.

## 3. Physics invariants

Noise off: unit norm and the exact free/coherent evolutions. Noise on: the conditional
width rides the deterministic Riccati flow pathwise; the covariance triple stays on the
purity shell SigmaX SigmaP - C^2 = hbar^2/4 (invariant and attracting at rate
4 lam SigmaX, verified symbolically and on the grid); the ensemble mean squared norm is
pinned at one (the likelihood is a martingale); unweighted quadratic-form averages obey
the Lindblad moment chain in both position and momentum; the law of total variance links
conditional and unconditional variances.

## 4. Idiomatic plan

Functional WL only (FoldList, Map, MapThread, With/Block; Module solely to mint fresh
NDSolve symbols); packed-array contractions for ensemble readouts; no Do/For/While/
AppendTo, no Quiet, no Print. Symbolic first: every reference is derived (Solve,
DSolve, MatrixExp closed form) before any NDSolve integration, and NDSolve appears only
where a rough driver makes closed forms impossible (the driven means, the field itself).

## 5. Required regimes and limits

| Regime | Where the essay covers it |
|---|---|
| General symbolic case | three-point symbolic Lindblad + conversion cells; symbolic steady triple; symbolic closed-form covariance flow; symbolic unconditional chain |
| Exactly solvable special case | free Gaussian dispersion; coherent state in the trap; exact ring propagator; steady state; law of total variance |
| Asymptotic / limiting case | t -> Infinity limit of the closed flow equals the steady triple; hbar -> 0 collapse to the classical Kalman limit; lambda-scaling of the steady width (fourth root); noise-step and record-step refinement limits |
| Numerical reference case | grid runs vs the exact covariance flow (pathwise); cross-engine filtering vs the Bayes/split-step engine, in the trap and in the double well; paired same-path noise-step sweep |
| Failure / edge case | the nonlocal-coefficient refusal; the Ito-damping negative control; the seam-crossing packet firing the edge guard with the boundary certified under strain; the effective-sample-size erosion at long horizon; the double well outside every closure |

Trivial base case refused: a single free Gaussian trajectory with the measurement off,
compared to nothing.
