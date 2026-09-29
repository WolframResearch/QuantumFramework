# Gate-class fusion prototype: benchmark results (2026-09-26)

WL 15.0.1, QF working tree, one kernel. Every time is one run, measured after `ClearSystemCache[]`, applied to a random normalized state. Every arm reproduces QF's output state to within 1.6e-15, and the run raised no messages. Raw log: `bench-bench.log`; code: `fusion.wl`, `bench.wls`. During the run one unrelated kernel (PID 47529) and the separate `Cos[qo]` fix session were also using the CPU.

Arms:
- **QF**: QF's default tensor-network contraction, and its `"Greedy"` contraction path.
- **B**: the circuit applied gate by gate to the state vector, each gate along its own qubit axes, with no gate classes and no fusion.
- **C**: gates classified by their nonzero pattern; consecutive diagonal and permutation gates fused into structured `PermutationMatrix` / `DiagonalMatrix` pairs with `Dot`; blocks applied with `Dot`.
- **D**: the same fusion stored as plain vectors (a gather list and a factor list); blocks applied with `Part` and `Times`.
- **E**: D, plus QFT and inverse-QFT blocks applied as an FFT.

"Build" is the time to classify and fuse the gates. "Apply" is the time to run the fused steps on the state.

| Circuit | gates | share fusable | QF default | QF greedy | B | C build + apply | D build + apply | E total |
|---|---|---|---|---|---|---|---|---|
| QFT, 14 qubits | 112 | 0.88 | 1.84 | 0.80 | 0.100 | 1.28 + 0.014 | 0.29 + 0.081 | **0.003** |
| QFT, 16 qubits | 144 | 0.89 | 3.55 | 3.71 | 0.162 | 3.70 + 0.045 | 1.41 + 0.397 | **0.004** |
| QPE, 12 + 2 qubits | 110 | 0.67 | 1.71 | 0.85 | 0.090 | 0.67 + 0.018 | 0.42 + 0.074 | **0.020** |
| QPE, 14 + 2 qubits | 142 | 0.70 | 3.72 | 4.07 | 0.159 | 2.95 + 0.073 | 1.08 + 0.419 | **0.056** |
| Bernstein-Vazirani, 17 qubits | 50 | 0.34 | 3.06 | 3.36 | **0.083** | 0.78 + 0.063 | 0.20 + 0.054 | 0.233 |
| QAOA p=2, 12 qubits | 168 | 0.79 | 1.69 | 0.48 | **0.056** | 0.35 + 0.002 | 0.09 + 0.001 | 0.101 |
| QAOA p=2, 14 qubits | 224 | 0.81 | 2.90 | 1.63 | **0.119** | 1.77 + 0.007 | 0.27 + 0.007 | 0.285 |
| Clifford+T, 12 qubits | 500 | 0.80 | 6.67 | 2.16 | **0.159** | 1.10 + 0.043 | 0.40 + 0.089 | 0.490 |
| Clifford+T, 14 qubits | 500 | 0.79 | 7.87 | 3.96 | **0.245** | 3.48 + 0.181 | 1.17 + 0.456 | 1.576 |
| dense 2-qubit gates, 12 qubits | 150 | 0 | 1.93 | 0.43 | **0.053** | 0.05 + 0.008 | 0.05 + 0.008 | 0.062 |
| dense 2-qubit gates, 14 qubits | 150 | 0 | 2.32 | 1.24 | **0.087** | 0.05 + 0.043 | 0.05 + 0.048 | 0.095 |

(all times in seconds; bold = fastest for a single run of the circuit)

## What the numbers say

1. **Most of QF's time is its own per-gate cost, not missing structure.** Arm B applies every gate the ordinary way along its qubit axes, with no structure at all, and is 8-37 times faster than QF's best method on every circuit, the dense control included.
2. **Fusion does not pay for a circuit run once.** Embedding a gate in the full register costs about as much as applying it, so building fused blocks costs more than it saves. Stored as structured arrays (C), the build is 1.6-6.6 times slower than as plain vectors (D).
3. **Once built, the fused blocks apply faster than gate-by-gate application.** C's apply step is the fastest or tied with D's on every structured circuit. For QAOA at 14 qubits it is 17 times faster than B, and C breaks even after about 16 runs. D breaks even after about 2-6 runs, depending on the circuit. Fusion is worth building when the same circuit is applied many times: many input states, repeated Trotter steps, or repeated shots of the same layer.
4. **Recognizing a named block is the one large single-run win.** The QFT as an FFT (E) beats gate-by-gate application by 40 times at 16 qubits, and QF by about 900 times. In phase estimation it removes the inverse-QFT share of the circuit.
5. **Fusion costs nothing when there is no structure.** On the dense circuits C, D and E do the same work as B.

The build cost measured here belongs to this prototype (it embeds each gate through a table of the bits of every basis state). A build written with integer bit operations could be much cheaper, which would move the break-even point for C and D down.
