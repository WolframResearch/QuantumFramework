"""Phase-resolved Stim benchmark. Companion to bench_phases_qf.wls / bench_phases_qc.jl;
shared output schema  RESULT|stim|<phase>|n=<n>|k=<k>|<ms>.

Phases (each timed separately):
  ingest_py    Python append loop: op tuples -> stim.Circuit. One pybind11 call per
               gate, so this is the slow ingest route; reported because it is what a
               Python-driven pipeline actually pays.
  ingest_text  stim.Circuit(text) parse of the circuit's own serialization: the C++
               bulk ingest route (Stim's file-based workflow).
  simulate     TableauSimulator().do(circuit) on a pre-built circuit: the engine.
  materialize  current_inverse_tableau().to_numpy(): dense bit arrays of the final
               Clifford data in Stim's native (inverse-map) convention.
  mat_forward  canonical_stabilizers(): forward stabilizers via Gaussian elimination.
               MORE work than QF's unpack (which does no elimination); reported so
               the reader can pick either boundary, never silently mixed in.
  e2e          ingest_text + simulate + materialize summed from the phase rows.
  measure1     one qubit measurement on the evolved state (sampling collapse; QF's
               analogue returns both outcome branches, which is more work).

The neutral JSON op list is parsed outside all timed regions, as on every engine.

Estimator: every phase row is min over 7 evaluations, a budget pinned across all
three engines (bench_phases_qf.wls, bench_phases_qc.jl): the min estimator is
downward-biased by sample count, so asymmetric rep budgets would bias sub-ms
cross-engine ratios for reasons unrelated to the engines.
"""
import json, os, time, stim
from stim._detect_machine_architecture import _UNSTABLE_detect_march

HERE = os.path.dirname(os.path.abspath(__file__))


def timeit(f, reps=7):
    best = 1e9
    for _ in range(reps):
        t = time.perf_counter(); f(); best = min(best, time.perf_counter() - t)
    return best * 1000  # ms


def pr(phase, n, k, ms):
    print(f"RESULT|stim|{phase}|n={n}|k={k}|{ms:.4f}")


def build_py(ops):
    c = stim.Circuit()
    for op in ops:
        if op[0] == 'CNOT':
            c.append('CNOT', [op[1], op[2]])
        else:
            c.append(op[0], [op[1]])
    return c


print(f"stim {stim.__version__} march={_UNSTABLE_detect_march()}")
for n in (100, 500, 1000):
    with open(os.path.join(HERE, f"stab_ops_{n}.json")) as fh:
        ops = json.load(fh)
    k = len(ops)
    circuit = build_py(ops)
    assert len(circuit) == k, "instruction count mismatch: fusion would skew the comparison"
    text = str(circuit)

    pr("ingest_py", n, k, timeit(lambda: build_py(ops)))
    pr("ingest_text", n, k, timeit(lambda: stim.Circuit(text)))
    pr("simulate", n, k, timeit(lambda: stim.TableauSimulator().do(circuit)))

    s = stim.TableauSimulator(); s.do(circuit)
    pr("materialize", n, k, timeit(lambda: s.current_inverse_tableau().to_numpy()))
    pr("mat_forward", n, k, timeit(lambda: s.canonical_stabilizers()))

    # single measurement on the evolved (entangled) state: rebuild each rep so the
    # collapse is never measured on an already-collapsed tableau
    def measure_once():
        sm = stim.TableauSimulator(); sm.do(circuit)
        t0 = time.perf_counter(); sm.measure(0)
        return (time.perf_counter() - t0) * 1000
    pr("measure1", n, k, min(measure_once() for _ in range(7)))
print("DONE|stim")
