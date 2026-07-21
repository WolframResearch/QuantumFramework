"""Superseded by bench_stim_all.py, which reports these sizes plus n=1000 and the per-gate control.

Kept because published tables cite it by name. Like bench_stim_all.py it batches the
stream into one stim.Circuit and issues a single sim.do(circuit): driving Stim one gate
per Python call measures the pybind11 boundary (~0.5 us/gate, flat in n) instead of the
tableau engine, which is O(n) per gate and cannot be flat.
"""
import json, time, stim

def timeit(f, reps=5):
    best = 1e9
    for _ in range(reps):
        t = time.perf_counter(); f(); best = min(best, time.perf_counter() - t)
    return best

for n in (100, 500):
    ops = json.load(open(f"/tmp/stab_ops_{n}.json"))
    circuit = stim.Circuit()
    for op in ops:
        circuit.append('CNOT', [op[1], op[2]]) if op[0] == 'CNOT' else circuit.append(op[0], [op[1]])
    ms = timeit(lambda: stim.TableauSimulator().do(circuit)) * 1000
    print(f"stim n={n} m={len(ops)}: {ms:.3f} ms ({ms * 1e6 / len(ops):.1f} ns/gate)")
