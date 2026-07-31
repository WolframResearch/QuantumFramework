# Confirmed tier-mixing defects (sweep of 2026-08-01)

310 probes across 6 areas; 24 claims adversarially confirmed; the three root
causes below were fixed the same day (operator constructor Diagonal catch-all,
raw ArrayReshape in QuantumPartialTrace, ArrayConjugate deferred regression).
Still open, each verified with a repro in the sweep transcript
(subagents/workflows/wf_2461c4e2-4af/journal.jsonl):

1. Measurement on a symbolic state FABRICATES probabilities: Z/X/POVM on
   QuantumState[VectorSymbol["v",2]] messages TensorContract::contr then
   returns {1,0} (or {1/2,1/2}), independent of the state; Mean/Entropy
   inherit. Site: measurement contraction, same scalar blind spot.
2. Symbolic operator APPLIED to a state (any tier): amplitude is
   1.*HoldForm[...]; the operator now constructs correctly but application
   goes through a basis-split path that re-wraps.
3. Multi-qudit gates on non-explicit states (CNOT on VectorSymbol["w",4],
   H on one qudit of a 2-qubit evolved state): AssociationThread::idim message
   storm; dims wrong (16 for 2 qubits) or ParameterSpec destroyed.
4. QuantumPartialTrace drops ParameterSpec of a lazy evolved state even when
   values are right; rt[1.] then stays inert.
5. QuantumTensorProduct with a symbolic factor: KroneckerProduct unevaluated,
   ProbabilitiesList fabricated {1,0,0,0}.
6. VectorState/StateVector on a symbolic MATRIX state returns the dim-4 Bend
   (PureStateQ is False for Unknown type) - silent dimension change.
7. Concurrence of a symbolic 2-qubit state: elementwise square, no trace.
8. QuantumOperator[Exp[I phi "PauliX"]] hits $RecursionLimit (pre-existing,
   breaks SPSRGradientValues' documented generator form).
