(* Experiment 3: the precise analogy. SymPhase outcomes are degree-1 over    *)
(* F_2 (Z_2 signs); diagonal Clifford+T phases are degree-<=3 over Z_8.       *)
(* And the crux: measurement = classical DISJUNCTION (substitute -> one       *)
(* physical state); T = coherent SUPERPOSITION (substitute reaches no state). *)
(* Run: wolframscript -file exp3-analogy-table.wls                            *)

PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"];

(* ---- SymPhase side: outcomes are AFFINE-LINEAR (degree 1) over F_2. ----- *)
(* (Fresh kernel: counter starts at 0, so the two free bits are s_1, s_2.)   *)
parityCheck = PauliStabilizer[3][{"H" -> 1, "H" -> 2, "CNOT" -> {1, 3}, "CNOT" -> {2, 3}}];
smPar = parityCheck["SymbolicMeasure", {1, 2, 3}];
free = DeleteDuplicates @ Cases[smPar["Phase"], _\[FormalS], Infinity];
(* The signs are 1 - 2 s, i.e. each sign entry is affine-linear in the bits. *)
Echo[smPar["Phase"], "G1 symbolic phase column (linear F_2 forms in s_k):"];
readBits[ps_] := IntegerDigits[FirstPosition[Normal[ps["State"]["StateVector"]], x_ /; x != 0][[1]] - 1, 2, ps["Qubits"]];
table = Flatten[Table[readBits[smPar["SubstituteOutcomes", Thread[free -> {a, b}]]], {a, 0, 1}, {b, 0, 1}], 1];
Echo[table, "G1 (m1,m2,m3) over the 4 free-bit assignments:"];
Echo[AllTrue[table, Mod[Total[#], 2] == 0 &], "G1 m1 XOR m2 XOR m3 = 0 always (a DEGREE-1 relation):"];

(* Substituting a measurement object yields ONE basis state (disjunction).   *)
oneShot = smPar["SubstituteOutcomes", Thread[free -> {1, 0}]];
sv = Chop[oneShot["State"]["StateVector"] // Normal];
Echo[Count[sv, x_ /; Abs[x] > 10^-9], "G2 nonzero amplitudes after substitution (=1 => a single basis state):"];

(* ---- T side: phases are degree-<=3 over Z_8, and the frame is a SUM. ---- *)
w = Exp[I Pi/4];
(* Degree ladder over Z_8: T (deg 1), CZ (deg 2 via 4 x1 x2), CCZ (deg 3).   *)
degrees = <|"T_j (single)" -> "x_j", "CZ (pair)" -> "4 x_i x_j", "CCZ (triple)" -> "4 x_i x_j x_k"|>;
Echo[degrees, "H1 Z_8 phase-polynomial terms by gate (degree 1,2,3):"];

(* The T frame is a coherent sum: substituting its 'sign' gives a SUMMAND,    *)
(* never the physical state. Show the frame state has TWO nonzero amplitudes  *)
(* (a superposition), while each component has its own support.               *)
fr = PauliStabilizer[1][{"H", "T"}];
svFrame = Chop[fr["StateVector"] // Normal];
Echo[Count[svFrame, x_ /; Abs[x] > 10^-9], "H2 nonzero amplitudes of T|+> (=2 => genuine superposition):"];
Echo[fr["Coefficients"] // Simplify, "H2 the two complex weights that are SUMMED, not selected:"];

(* ---- Side-by-side ring/degree summary. ---- *)
summary = Dataset[<|
  "carrier"      -> <|"measurement/Pauli" -> "sign column R", "diagonal Clifford+T" -> "phase column"|>,
  "ring"         -> <|"measurement/Pauli" -> "Z_2 (F_2 bits)", "diagonal Clifford+T" -> "Z_8 (8th roots)"|>,
  "poly degree"  -> <|"measurement/Pauli" -> "1 (affine XOR)", "diagonal Clifford+T" -> "<= 3 (CNOT-dihedral)"|>,
  "combination"  -> <|"measurement/Pauli" -> "disjunction (OR, probabilities)", "diagonal Clifford+T" -> "still one branch (phase absorbed)"|>,
  "stays rank 1" -> <|"measurement/Pauli" -> "yes (one tableau)", "diagonal Clifford+T" -> "yes (one phase-poly state)"|>
|>];
Echo[summary, "I1 analogy table:"];

Echo["exp3 done", "STATUS:"];
