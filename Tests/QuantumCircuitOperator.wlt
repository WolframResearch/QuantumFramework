(* QuantumCircuitOperator tests. *)

BeginTestSection["QuantumCircuitOperator - constructors"]

VerificationTest[QuantumCircuitOperator["GHZ"[3]]["Arity"], 3, TestID -> "GHZ-3"];

VerificationTest[QuantumCircuitOperator["Bell"]["Arity"], 2, TestID -> "Bell-bare"];

VerificationTest[QuantumCircuitOperator["Toffoli"]["Arity"], 3, TestID -> "Toffoli-bare"];

VerificationTest[QuantumCircuitOperator["Fredkin"]["Arity"], 3, TestID -> "Fredkin-bare"];

VerificationTest[QuantumCircuitOperator["Fourier"[3]]["Arity"], 3, TestID -> "Fourier-3"];

VerificationTest[QuantumCircuitOperator["InverseFourier"[3]]["Arity"], 3, TestID -> "InverseFourier-3"];

(* sequence-of-gates list (semantic composition - kept as List) *)
VerificationTest[
  QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}}]["Arity"],
  2,
  TestID -> "Sequence-of-gates-list"
];

(* Pauli-string circuit shortcut *)
VerificationTest[
  QuantumCircuitOperator["XYZ"]["Arity"],
  3,
  TestID -> "PauliString-XYZ"
];

EndTestSection[]


BeginTestSection["QuantumCircuitOperator - basic action"]

(* Bell circuit takes |00> -> Bell pair *)
VerificationTest[
  Chop[Norm[
    QuantumCircuitOperator["Bell"][QuantumState["00"]]["StateVector"] -
      QuantumState["Bell"]["StateVector"]
  ]],
  0,
  TestID -> "Bell-circuit-on-00"
];

VerificationTest[QuantumCircuitOperator["GHZ"[3]]["Elements"] // Length, 3, TestID -> "GHZ-gate-count"];

EndTestSection[]


BeginTestSection["QuantumCircuitOperator - shortcut roundtrip"]

(* QuantumShortcut compresses a circuit to its named-shorthand form;
   re-feeding it through the QuantumCircuitOperator constructor reconstructs
   an equivalent circuit. *)

VerificationTest[
  QuantumCircuitOperator @ QuantumShortcut[QuantumCircuitOperator["Bell"]] == QuantumCircuitOperator["Bell"],
  True,
  TestID -> "Roundtrip-Bell"
];

VerificationTest[
  QuantumCircuitOperator @ QuantumShortcut[QuantumCircuitOperator["GHZ"[3]]] == QuantumCircuitOperator["GHZ"[3]],
  True,
  TestID -> "Roundtrip-GHZ-3"
];

VerificationTest[
  QuantumCircuitOperator @ QuantumShortcut[QuantumCircuitOperator["Toffoli"]] == QuantumCircuitOperator["Toffoli"],
  True,
  TestID -> "Roundtrip-Toffoli"
];

VerificationTest[
  QuantumCircuitOperator @ QuantumShortcut[QuantumCircuitOperator["Fredkin"]] == QuantumCircuitOperator["Fredkin"],
  True,
  TestID -> "Roundtrip-Fredkin"
];

(* hand-built circuit *)
VerificationTest[
  With[{qc = QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}, "Phase"[Pi/4] -> 2}]},
    QuantumCircuitOperator @ QuantumShortcut[qc] == qc
  ],
  True,
  TestID -> "Roundtrip-handbuilt"
];

EndTestSection[]


BeginTestSection["QuantumCircuitOperator - Qiskit export"]

(* Qiskit export must hand Python a tuple (name, args), not a WL function call.
   The bug surfaced as `'Measure'[(1,)]` on the Python side; lock it down. *)

VerificationTest[
  Quiet @ Head @ QuantumCircuitOperator[{{1}}]["Qiskit"],
  QiskitCircuit,
  TestID -> "Qiskit-measurement-from-list"
];

VerificationTest[
  Quiet @ Head @ QuantumCircuitOperator["Bell"]["Qiskit"],
  QiskitCircuit,
  TestID -> "Qiskit-Bell"
];

VerificationTest[
  Quiet @ Head @ QuantumCircuitOperator[{"U"[Pi/3, Pi/4, Pi/5] -> 1}]["Qiskit"],
  QiskitCircuit,
  TestID -> "Qiskit-U3"
];

EndTestSection[]


BeginTestSection["QuantumCircuitOperator - failure"]

VerificationTest[
  Quiet @ QuantumCircuitOperator["NotACircuit"[2]],
  Failure["InvalidName", _],
  SameTest -> MatchQ,
  TestID -> "InvalidName-call-form"
];

VerificationTest[
  Quiet @ QuantumCircuitOperator["NotACircuit"],
  Failure["InvalidName", _],
  SameTest -> MatchQ,
  TestID -> "InvalidName-bare"
];

EndTestSection[]
