# Named Quantum Circuits and Operators: index

A catalog of every named construction `QuantumCircuitOperator` and
`QuantumOperator` accept, verified against the kernel (QF 2.1.3). 43 named
circuits + 108 named operators, 7 shared. Pages are ordered by impact for a
quantum audience, not by the registry's internal order. See `STYLE.md` for the
authoring contract and `ISSUES.md` for the running defect log.

## Tier 1 - Complete algorithms (the hook)
- [Grover](named-quantum-circuits.md)
- [Deutsch and Deutsch-Jozsa](deutsch-jozsa.md)
- [Simon's algorithm](simon.md)
- [Bernstein-Vazirani](bernstein-vazirani.md)
- [Fourier transform and phase estimation](phase-estimation.md)

## Tier 2 - Oracles (the engine the algorithms compile to)
- [Oracles](oracles.md) - `BooleanOracle`, `BooleanOracleR`, `PhaseOracle`, `GrayOracle`

## Tier 3 - Composition and simulation
- [Composition and simulation](composition.md) - `Controlled`, `Switch`, `Multiplexer`, `ControlledMultiplexer`, `Trotterization`

## Tier 4 - Beyond unitary gates (the differentiation)
- [Nonclassicality tests](foundations.md) - `CHSH`, `LeggettGarg`, `WignerCHSH`
- [Non-unitary elements](open-systems.md) - `Measurement`, `Reset`, `Trace`, `Discard` (and a note on the dynamics generators)
- [Tensor-network atoms](zx-tensor.md) - `Cup`, `Cap`, `Braid`, `Copy`, `ZSpider`
- [Qudits and angular momentum](qudits.md) - `Shift`, `ShiftPhase`, `SUM`, `FlipSign`, `JZ`, `J+`, `WignerD`

## Tier 5 - Gate gadgets and encoders
- [Gate gadgets and encoders](gadgets.md) - `Bell`, `GHZ`, `Toffoli`, `Fredkin`, `Magic`, `Number`, `PhaseNumber`, `Graph`, `QuantumState`

## Tier 6 - Primitive gate set (reference)
- [Primitive gates](primitives.md) - Paulis, `H`, `S`, `T`, `V`, rotations, `U`, `P`, `CNOT`, `CZ`, `SWAP`

## Not documented as circuits (verified not to fit the frame)
Logged in `ISSUES.md`: `Channel`, `Decohere`, `Hamiltonian`, `Liouvillian` are
open-system dynamics objects for `QuantumChannel` / `QuantumEvolve`, not circuit
gates; `JX` / `JY` return the `JZ` matrix (bug); several forms crash the kernel
(`Decohere[p]`, `Spider[phase]`) or reject the association variable-spec.
