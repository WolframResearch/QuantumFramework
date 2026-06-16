#!/usr/bin/env wolframscript
PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"];
sep[s_] := Print["==== ", s, " ===="];

(* A. Symbolic signs / SymbolicMeasure *)
sep["A symbolic measure"];
Module[{bell, sm},
  bell = PauliStabilizer[{"XX", "ZZ"}];
  sm = bell["SymbolicMeasure", 1];
  Print["bell stabilizers: ", bell["Stabilizers"]];
  Print["after SymbolicMeasure[1] signs: ", sm["Signs"]];
  Print["SampleOutcomes 4: ", (#["Signs"] & /@ sm["SampleOutcomes", 4])];
];

(* B. Exact arithmetic + bit-reproducible seeded random Clifford *)
sep["B exact + reproducible"];
Module[{a, b},
  SeedRandom[42]; a = PauliStabilizer["Random"[6]];
  SeedRandom[42]; b = PauliStabilizer["Random"[6]];
  Print["seeded identical tableau: ", a["Tableau"] === b["Tableau"]];
  Print["signs exact integers: ", MatchQ[a["Signs"], {(-1|1)..}]];
  Print["GHZ expectation ZZZ exact (head): ", Head[PauliStabilizer[{"XXX","ZZI","IZZ"}]["Expectation","ZZZ"]]];
];

(* C. Closed-form entropy at large n, no OOM *)
sep["C entropy large n"];
Module[{ghz, n = 200, specs, t},
  specs = Join[{"H" -> 1}, Table["CNOT" -> {i, i+1}, {i, 1, n-1}]];
  ghz = PauliStabilizer[n]["ApplyCircuit", specs];
  t = First @ AbsoluteTiming[ghz["Entropy", Range[100]]];
  Print["GHZ n=200 half-chain entropy: ", ghz["Entropy", Range[100]], "  (", 1000. t, " ms)"];
  Print["single-qubit cut entropy: ", ghz["Entropy", 1]];
];

(* D. StabilizerFrame stabilizer rank / magic *)
sep["D stabilizer frame / rank"];
Module[{plus, f, f2},
  plus = PauliStabilizer[1]["H", 1];
  f = plus["T", 1];
  Print["T on |+>: frame length (stabilizer rank): ", f["Length"]];
  f2 = f["T", 1];   (* note: T on a frame doubles again only if it hits a stab *)
  Print["coefficients: ", f["Coefficients"]];
  Print["materialized state vector: ", Chop @ f["StateVector"]];
];

(* E. CliffordChannel as first-class algebraic object *)
sep["E clifford channel"];
Module[{cc, comp},
  cc = CliffordChannel["Identity", 2];
  Print["identity channel rank: ", cc["Rank"], " in/out: ", {cc["InputQubits"], cc["OutputQubits"]}];
  comp = cc[cc];
  Print["self-composition rank: ", comp["Rank"], " source: ", comp["Source"]];
  Print["from stabilizer (Bell) Choi rows: ", CliffordChannel[PauliStabilizer[{"XX","ZZ"}]]["Rank"]];
];

(* F. GraphState / LocalComplement *)
sep["F graph state"];
Module[{g, gs, lc},
  g = CycleGraph[5];
  gs = GraphState[g];
  Print["graph-state stabilizers: ", gs["Stabilizers"]];
  lc = LocalComplement[gs, 1];
  Print["LocalComplement[.,1] edge count: ", lc["EdgeCount"], " (was ", gs["EdgeCount"], ")"];
];

(* G. Conversion into the QF object model *)
sep["G object-model conversion"];
Module[{ghz},
  ghz = PauliStabilizer[{"XXX","ZZI","IZZ"}];
  Print["State head: ", Head[QuantumState[ghz]]];
  Print["Operator head: ", Head[ghz["QuantumOperator"]]];
  Print["Circuit head: ", Head[ghz["QuantumCircuit"]]];
  Print["State vector: ", Chop @ Normal @ QuantumState[ghz]["StateVector"]];
];

(* H. Pauli-string measurement + expectation *)
sep["H pauli-string measure / expectation"];
Module[{ghz},
  ghz = PauliStabilizer[{"XXX","ZZI","IZZ"}];
  Print["<ZZZ> = ", ghz["Expectation","ZZZ"], "  <XXX> = ", ghz["Expectation","XXX"], "  <ZII> = ", ghz["Expectation","ZII"]];
  Print["measure XXX (deterministic): ", Keys @ ghz["M","XXX"]];
];

(* TRAP 1: PauliStabilizer[circuit] O(4^n) vs ApplyCircuit *)
sep["TRAP PauliStabilizer[circuit] vs ApplyCircuit"];
Module[{n = 12, specs, qc, t1, t2, r1, r2},
  specs = Join[{"H" -> 1}, Table["CNOT" -> {i, i+1}, {i, 1, n-1}]];
  qc = QuantumCircuitOperator[specs];
  t1 = First @ AbsoluteTiming[r1 = PauliStabilizer[qc]];
  t2 = First @ AbsoluteTiming[r2 = PauliStabilizer[n]["ApplyCircuit", specs]];
  Print["PauliStabilizer[qc]  n=12: ", 1000. t1, " ms"];
  Print["ps[ApplyCircuit]     n=12: ", 1000. t2, " ms"];
  Print["same state (Stabilizers match): ", r1["Stabilizers"] === r2["Stabilizers"]];
];

(* TRAP 2: global-phase contract on gate-updated tableau *)
sep["TRAP global-phase contract"];
Module[{qs, ps, viaGate, viaCtor},
  qs = QuantumState["0"];
  ps = PauliStabilizer[qs];
  viaGate = Chop @ Normal @ ps["Y", 1]["State"]["StateVector"];
  viaCtor = Chop @ Normal @ PauliStabilizer[QuantumOperator["Y"][qs]]["State"]["StateVector"];
  Print["gate-updated  Y|0> state: ", viaGate];
  Print["constructor   Y|0> state: ", viaCtor];
  Print["equal up to global phase (mag): ", Chop[Abs[viaGate] - Abs[viaCtor]]];
];

Print["DONE-VERIFY"];
