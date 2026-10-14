---
Template: Default
Title: Nonclassicality Tests
Context: Wolfram`QuantumFramework`
---

# Nonclassicality Tests

These circuits stage the standard tests that separate quantum correlations from
classical ones: a spatial Bell test, its temporal analogue, and a phase-space
version. Each is parametrized by a single measurement angle, with the Alice,
Bob, and Charlie roles laid out as labelled sub-circuits, barrier-separated in
the spatial and temporal tests. The page builds the measurement circuits; the
correlation value and its classical bound (for CHSH, the $2\sqrt2$ Tsirelson
violation of the bound 2) are a classical post-processing step on the outcome
statistics, not computed here.

## Spatial Bell test: `CHSH`

The CHSH Bell test: a shared entangled pair, Charlie's preparation, and Alice and Bob's measurements, with Bob's rotation carrying the angle, ending in four readouts:

```wl
QuantumCircuitOperator[{"CHSH"[]}]["Diagram"]
```

---

The measurement angle is the single knob; setting it changes Bob's rotation:

```wl
QuantumCircuitOperator[{"CHSH"[Pi/3]}]["Diagram"]
```

## Temporal Bell test: `LeggettGarg`

The temporal analogue of CHSH, correlating measurements of one system at different times rather than two systems in space:

```wl
QuantumCircuitOperator[{"LeggettGarg"[]}]["Diagram"]
```

## Phase-space Bell test: `WignerCHSH`

The CHSH test built on the Wigner minimal informationally complete basis, a discrete phase-space *representation* of the qubit, rather than a physical measurement of phase space:

```wl
QuantumCircuitOperator[{"WignerCHSH"[]}]["Diagram"]
```
