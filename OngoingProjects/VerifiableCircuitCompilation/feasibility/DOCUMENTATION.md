# Documentation and runtime evidence

The installed WL 15.0 reference notebooks and package sources were read locally.
The authoritative online links below identify the same documentation. Temporary
reading copies and local installation paths are not included in this repository.

| Symbol/source | Material read | Consequence exercised in this study |
|---|---|---|
| [MiniMaxApproximation](https://reference.wolfram.com/language/FunctionApproximations/ref/MiniMaxApproximation.html) | Complete symbol page; installed Approximations.m header, limitations, examples, usages and options | Load FunctionApproximations`; use degree `{degree,0}`, WorkingPrecision 50, extract `{abscissae,{polynomial,error}}`; reject unrecognized output |
| [Function Approximations tutorial](https://reference.wolfram.com/language/FunctionApproximations/tutorial/FunctionApproximations.html) | Rational interpolation/minimax section, convergence, zero/denominator failures, Brake, derivatives, and restart examples | Relative minimax is discovery only; quantize coefficients and independently prove bounds. No claim of globally optimal quantized coefficients |
| [Resolve](https://reference.wolfram.com/language/ref/Resolve.html) | Complete installed symbol page, including exact domains, quantifiers, WorkingPrecision and possible issues | Rational polynomial inequalities, explicit Reals, no finite WorkingPrecision. Timeout is not True |
| [BooleanConvert](https://reference.wolfram.com/language/ref/BooleanConvert.html) | Complete installed symbol page | ANF and ESOP exist; conditional simplifications require proof of their care set |
| [BooleanMinimize](https://reference.wolfram.com/language/ref/BooleanMinimize.html) | Complete installed symbol page | Minimize a full-adder Boolean expression in ANF, then verify with TautologyQ. Boolean expression length is not a quantum resource objective |
| [PolynomialReduce](https://reference.wolfram.com/language/ref/PolynomialReduce.html) | Complete installed symbol page | Verify bit idempotence over the rationals; Modulus changes coefficient arithmetic and must not be imposed on arbitrary-angle phases |
| [LinearOptimization](https://reference.wolfram.com/language/ref/LinearOptimization.html) | Details/options, basic examples, Integer Variables, Possible Issues | Integer-domain syntax is exercised. Mixed-integer solves use machine precision and are not exact proof certificates. No MILP optimum is claimed or used for architecture ranking |
| Official kickmix installed docstrings | init_lookup, del_lookup, unary_iteration, iadd, z_pow, mark, operation_statistics | Joint coefficient lookup and deletion are charged; Z_POW is in half-turns; static counts include classical conditions |
| Official invsqrt files | Entire nr_reference_circuit.py and verify_submission_invsqrt.py | Restore x and scratch; all x; modulo-one phase; per-partial-product truncation; margin-dependent synthesis price |

The kernel runner and independent reviewer reproduce relevant signatures and
return shapes. No new QuantumFramework API is introduced here: these are core-WL
symbolic certificates and official kickmix resource probes. The earlier QF bridge
is left unchanged.
