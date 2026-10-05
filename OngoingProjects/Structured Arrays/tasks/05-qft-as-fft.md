# Brief 5: the QFT applied to a state as a fast Fourier transform (investigate first)

Read `README.md` in this folder first: it has the workflow, the house rules and the traps every brief shares.

## Goal

The quantum Fourier transform on n qubits, applied to a state vector, is the discrete Fourier transform of its 2ⁿ amplitudes: one call to `Fourier`, which runs a fast Fourier transform. QF instead contracts the gates of its `"Fourier"` circuit, n(n+1)/2 Hadamard and controlled-phase gates and ⌊n/2⌋ swaps, as a tensor network. The task: decide whether QF should recognize a QFT block when it applies a circuit to a numeric state, and if so, apply that block as an FFT.

## What QF does today

Measured on main at `f6299b43` by `baseline.wls qft` and `supplement2.wls` (raw output in `baseline-out/`), on a random normalized state:

| n | `qc[qs]`, default contraction | `qc[qs]`, greedy contraction path | `FourierMatrix[2^n] . v` | `Fourier[v]` |
|---|---|---|---|---|
| 14 | 3.0 s | 0.95 s | 0.003 s | 0.003 s |
| 16 | 3.9 s | 3.8 s | 0.004 s | 0.004 s |

QF's result equals `FourierMatrix[2^n] . v` to 1.6 × 10⁻¹⁵, and `Fourier[v]` with its default `FourierParameters` equals `FourierMatrix[2^n] . v` exactly (checked at 2¹⁰), so `Fourier` follows QF's convention with no change of parameters.

The circuit is `QuantumCircuitOperator["Fourier"[n, m]]` in `QuantumFramework/Kernel/QuantumCircuitOperator/NamedCircuits.m:517-529`: a Hadamard and controlled phase shifts on each qubit, then swaps, on qubits m + 1 to m + n; `"InverseFourier"` is its `"Dagger"`. Applying a circuit to a state (`QuantumCircuitOperator.m:106`) flattens it first, `TensorNetworkApply[qco["Flatten"], qs, ...]`, after which the gates no longer form a named block.

The gate-fusion prototype (`../prototype/RESULTS.md`, arm E) applied QFT and inverse-QFT blocks as FFTs: 40 times faster than applying the gates one by one at 16 qubits and about 900 times faster than QF, and in phase estimation it removed the inverse QFT's share of the time.

## Investigation (do this first, report before changing code)

1. **Can a QFT block be found before the circuit is flattened?** A `QuantumCircuitOperator` can hold another one as an element (phase estimation holds its inverse QFT this way; the survey read it as the last element). Find what identifies a `"Fourier"` or `"InverseFourier"` block reliably: its construction, its label, its qubits. A label alone is not enough, since users relabel circuits.
2. **Where would the FFT step go?** In the apply path, before `"Flatten"`: split the circuit into the parts before the block, the block and the parts after, apply the block to the state vector as a Fourier transform along the axis of its qubits (reshape the amplitudes so those qubits form one axis, `Fourier` along it, reshape back; `invQFT` in `../structured-arrays-in-qf.md`, section 10, does this for the first qubits), and the other parts as today.
3. **Scope.** Only a pure state with numeric amplitudes, and only a block on contiguous qubits (the `"Fourier"[n, m]` circuit's are). A symbolic state, a mixed state, or a circuit returned as an operator keeps the contraction.
4. **Measure** the gain on the QFT alone and on `"PhaseEstimation"` at 12 to 16 counting qubits, against both the default and the greedy contraction.

## Bigger alternative to weigh first

The same prototype found that most of QF's time on every circuit, dense ones included, is its cost per gate: applying each gate along its own qubit axes, with no structure at all, was 8 to 37 times faster than QF's best method (`../prototype/RESULTS.md`, point 1). Speeding up how QF applies any circuit would help far more circuits than recognizing the QFT. That is a separate, larger task; the report from this investigation should say which of the two Mads should fund first.

## What must not change

- Results of every circuit to roundoff (compare against today's result to a stated bound), including circuits with a QFT on part of the register.
- Symbolic and mixed states, and circuits used as operators (`qc["QuantumOperator"]`): today's route.

## Tests to write

- `QuantumCircuitOperator["Fourier"[16]]` on a random numeric state inside a time limit, equal to `Fourier[v]` to a stated bound.
- The QFT on qubits 3 to 8 of a 10-qubit state, and its inverse, against the explicit matrix product.
- `"PhaseEstimation"` with 12 counting qubits: the same state as today, to a stated bound, inside a time limit.
- Regression: a symbolic state through `"Fourier"[3]` keeps its exact result.
