"""Stim tableau throughput on the shared Clifford streams in stab_ops_{100,500,1000}.json.

Reports two numbers per size:

  batched  -- the whole op stream is compiled into one stim.Circuit and handed to the
              simulator in a single sim.do(circuit) call, so the gate loop runs inside
              C++. This is the engine measurement, and the apples-to-apples match for
              QF's ps["ApplyCircuit", specs] (qf_final_all.wl), which likewise builds
              its spec list outside the timed region and applies the stream in one call.

  per_gate -- one Python->pybind11 call per gate. This is bounded by call overhead, not
              by the tableau engine: cost per gate sits flat near 0.5 us for every n up
              to ~2000, whereas a genuine tableau update is O(n) per gate. Kept as a
              labelled control so the two can be told apart, never as the Stim figure.

Circuit construction is outside the timed region, matching the QF harness. The circuit
holds one instruction per gate (the H/S/CNOT stream never repeats a gate type back to
back, so stim performs no adjacent-instruction fusion): batching changes the call
boundary, not the work.
"""
import json, os, time, stim
from stim._detect_machine_architecture import _UNSTABLE_detect_march

HERE = os.path.dirname(os.path.abspath(__file__))


def timeit(f, reps=5):
    best = 1e9
    for _ in range(reps):
        t = time.perf_counter(); f(); best = min(best, time.perf_counter() - t)
    return best


def build_circuit(ops):
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
    m = len(ops)
    circuit = build_circuit(ops)
    assert len(circuit) == m, f"expected {m} instructions, got {len(circuit)}"

    def batched():
        stim.TableauSimulator().do(circuit)

    def per_gate():
        s = stim.TableauSimulator()
        for op in ops:
            if op[0] == 'H': s.h(op[1])
            elif op[0] == 'S': s.s(op[1])
            else: s.cnot(op[1], op[2])

    tb, tp = timeit(batched), timeit(per_gate)
    print(f"RESULT|stim_n{n}_m{m}_ms|{tb * 1000:.4f}|ns_per_gate|{tb / m * 1e9:.1f}")
    print(f"RESULT|stim_pergate_n{n}_m{m}_ms|{tp * 1000:.4f}|ns_per_gate|{tp / m * 1e9:.1f}")
