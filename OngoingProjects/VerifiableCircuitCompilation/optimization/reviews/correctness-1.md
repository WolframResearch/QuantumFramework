# Independent width-certification verification — round 1

Verdict: PASS. `OPEN ISSUES: 0`.

The fresh verifier read the complete WL artifacts, coefficient certificate,
generator, analysis and dependencies without the author's rationale. It used
a real Wolfram 15.0.0 ARM kernel and exact Python arithmetic. No tracked
artifacts were edited by the reviewer.

Reproduced evidence:

- `verify-widths.wls`: exactly 21 of 21 checks pass.
- `run-widths.wls`: complete stored output reproduced byte for byte; four
  candidate widths certify and four are rejected.
- Selected fractional widths `(30,23,16)` and physical widths `(30,22,15)`
  certify every one of 192 intervals, including positivity and capacity,
  at `eta = 15/268435456`.
- All saved coefficient tables match exact Python `Fraction` ties-to-even
  rounding. The source table and unchanged incumbent generator retain their
  published hashes.
- General signed error envelope, floor residual, exact special case,
  infinite-precision limit, normalization identity and phase modulus
  reproduce. Independent polynomial checks confirm the squared-envelope
  equivalence under the stated positivity assumptions.
- Malformed rows, widths, tolerances and arities return failures. Width
  endpoints 1 and 36 preserve table shape. Injected `False` and `"TIMEOUT"`
  both reject certification; the actual timeout fallback was also exercised.
- Independent integer arithmetic checks every one of 3,133,452 inputs:
  the identity `e0-t*e1`, row-error bounds and positivity/capacity of every
  partial accumulator hold.
- Analysis reproduces the complete screening report, resource estimates and
  16,380 compiled boundary shots across 4,095 inputs, including restoration,
  cleanup, exact phase agreement and priced counts.

No findings or fixes were required. The verdict covers the WL certification
and screening scope. It does not claim full official circuit verification,
speed benchmarking or global optimality.

OPEN ISSUES: 0
