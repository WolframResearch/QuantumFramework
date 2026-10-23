# Reduced-width quadratic phase synthesis

The source remains the hash-pinned, certified `B-quadratic-192.json` table.
No minimax coefficients are refitted. Its coefficient magnitudes are rounded
from 36 fractional bits to `(30, 23, 16)`, using exact ties-to-even rounding.
The three physical registers have `(30, 22, 15)` bits because the last two
coefficients are below one half. These rounded representations receive new
certificates; the old certificate alone does not establish their accuracy.

For interval origin `a = 1 + i/64`, let `t = M-a` lie in `[0,1/64]` and
let `p(t) = c0 - t(c1-t c2)` use the rounded dyadic coefficients. Fourteen
binary input fractions, numbered 7 through 20, contribute partial products
to each Horner stage. Each row drops less than one output unit. Thus
`0 <= e1 <= 14*2^-23` and `0 <= e0 <= 14*2^-30`, and the implemented value is

```
p(t) + e0 - t*e1.
```

`WidthOptimization.wl` keeps this signed dependence. `Resolve` proves, over
each complete real interval, that both endpoints of the error envelope obey
relative error `eta = 15/(16*2^24)`, with a positive lower endpoint.
Squaring is valid because positivity is proved before converting
`|sqrt(M)*y - 1| <= eta` into polynomial inequalities. Every solver result
must literally be `True`; `False` and timeouts both reject a design.

The first intermediate is positive; subtraction of nonnegative partial
products makes each accumulator decrease monotonically from its loaded
value to its positive final value. Bounds on the loaded integers therefore
also cover every partial accumulator. Inverse stages undo the same integer
translations. Normalization uses `M = x*4^s/2^20`; denormalization multiplies
the normalized result by `2^(s-10)`, preserving the relative error exactly.
Ten zero extension bits suffice for `0 <= s <= 10`. Phase extraction uses
unit-modulus diagonal gates. Computational-basis preservation, measurement
corrections, input restoration and scratch cleanup require compiled circuit
checks in addition to these scalar certificates.

The register schedule releases each product scratch register after its
stage. Phase extension bits are allocated only around phase extraction.
The selected schedule additionally measurement-uncomputes `c2` after the
first stage, reuses its clean storage and reloads it before the inverse
first stage. The lookup overhead is included in the compiled resource count.
The analysis compares this schedule against retaining `c2` and against
direct controlled addition, plus bounded adder-workspace choices. It makes
no claim of global optimality. The compiler's adder budget is a resource
probe, not an enforced bound on every internal primitive; reported qubits
always come from the finished circuit.

| Regime | Executable evidence |
|---|---|
| General symbolic | `HornerError[t,e0,e1]`, arbitrary row-count bound and quantified signed envelope |
| Exact special case | Zero row errors; interval origin; ties-to-even half-integers |
| Limit | Row bound vanishes as fractional precision tends to infinity |
| Numeric reference | All 3,133,452 inputs in an integer model; compiled boundary simulations |
| Failure and edge | Rejected narrow designs; all interval endpoints; empty table, invalid widths and inexact tolerance |

## Documentation consulted

- [Resolve](https://reference.wolfram.com/language/ref/Resolve.html): installed WL 15 documentation; polynomial quantifier elimination over `Reals`.
- [ForAll](https://reference.wolfram.com/language/ref/ForAll.html): the three-argument conditional form is also documented on the installed Resolve page.
- [Round](https://reference.wolfram.com/language/ref/Round.html): exact rounding is checked against half-integer cases in a real kernel and Python.
- Kickmix at commit `638777fac4379bc62e2a1e4a52cf1a1753895d48`: `iadd.pybind.cc`, `gen_iadd.cc` and `circuit_builder.pybind.cc` establish quantum offset alignment, scratch allocation, exact modulo arithmetic, lookup cleanup and `set_max_qubits` behavior.

## Reproduction

Run from this project directory with a configured Wolfram kernel and the
kickmix Python environment described in the parent README. Set
`QF_INVSQRT_DATA` to the external downloaded competition folder.

```sh
WolframKernel -noinit -noprompt -script optimization/run-widths.wls
WolframKernel -noinit -noprompt -script optimization/verify-widths.wls
python optimization/analyze.py --data "$QF_INVSQRT_DATA"
python -m unittest test_optimized_quadratic
python verify_optimized_quadratic.py --data "$QF_INVSQRT_DATA" --workers 4
```

`width-certificates.json` stores all accepted and rejected interval results.
`screening.json` contains model screens and resource estimates, which are
not official verification results. The full verifier's report, when PASS,
is `../results/optimized-quadratic-verified.json`. The incumbent generator
and its full report remain intact for comparison.
