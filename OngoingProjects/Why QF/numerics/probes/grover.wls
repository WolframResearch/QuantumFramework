(* The QF-Showcase Grover example: 8 qubits, 12 GroverPhase iterations on the uniform state, applied
   with Method -> {"TensorNetwork", "Path" -> "Greedy"} and with the default method, in a fresh kernel.
   Usage: grover.wls <QF paclet dir> greedy|greedy-binary|default|default-binary *)
{qfDir, which} = $ScriptCommandLine[[2 ;; 3]];
(* a TensorNetworks paclet directory in the environment variable TN_DIR is loaded ahead of the installed one *)
With[{tnDir = Environment["TN_DIR"]}, If[StringQ[tnDir], PacletDirectoryLoad[tnDir]]];
PacletDirectoryLoad[qfDir];
Needs["Wolfram`QuantumFramework`"];
P[l_, x_] := Print[StringPadRight[l, 40], " ", ToString[x, InputForm]];
SetAttributes[sec, HoldFirst];
sec[e_, limit_ : 300] := MemoryConstrained[
    TimeConstrained[With[{t = AbsoluteTiming[e]}, {N[Round[First[t], 1 / 1000]], Last[t]}], limit, {"timeout", None}],
    6 10^9, {"memory", None}];
P["TN / Arrays", {First[PacletFind["Wolfram/TensorNetworks"]]["Version"], First[PacletFind["Wolfram/Arrays"]]["Version"]}];
n = 8; marked = 200;
oracle = BooleanFunction[2^marked, n];
{tStep, step} = sec[QuantumCircuitOperator["GroverPhase"[oracle]]];
kopt = Round[Pi / (4 ArcSin[Sqrt[1 / 2^n]]) - 1 / 2];
{tGrover, grover} = sec[QuantumCircuitOperator[Join[Thread["H" -> Range[n]], ConstantArray[step, kopt]]]];
P["build step, grover s", {tStep, tGrover}];
{t, found} = sec[Switch[which,
    "greedy", grover[Method -> {"TensorNetwork", "Path" -> "Greedy"}],
    "greedy-binary", grover[Method -> {"TensorNetwork", "Path" -> "Greedy", "DiagonalHyperedges" -> False}],
    "default-binary", grover[Method -> {"TensorNetwork", "DiagonalHyperedges" -> False}],
    _, grover[]]];
P[which <> " apply s", t];
If[MatchQ[found, _QuantumState],
    P["max probability", N[Max[found["Probabilities"]]]];
    P["argmax", First[Ordering[found["Probabilities"], -1]] - 1]
];
