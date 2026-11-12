(* ::Package:: *)

Package["Wolfram`QuantumFramework`QEC`"]

PackageScope[magicStateInstructions]
PackageScope[magicTestInstructions]
PackageScope[magicTestHeraldChecks]
PackageScope[instructionStateVector]
PackageScope[$nonCliffordOps]


(* ============================================================================ *)
(* The pi/8 gate: a magic state, tested, and injected (Got26 sec. 13.3-13.4,    *)
(* figs. 13.7, 14.8b, 14.9).                                                    *)
(*                                                                              *)
(* Transversal gates are Clifford, and the Clifford group is not universal.     *)
(* The missing gate is R_pi/8 = T = diag(1, e^(i pi/4)), and it is done without  *)
(* a transversal T at all: prepare an encoded MAGIC STATE |A_L> = T_L |+_L>,     *)
(* check it, and spend it.  Spending it is fig. 13.7 -- a transversal CNOT from  *)
(* the magic block into the data, a transversal Z measurement of the data, and  *)
(* when the outcome is 1 the Clifford correction S X on the output -- so every  *)
(* gate the injection uses is one the protocol already has.                     *)
(*                                                                              *)
(* PREPARING |A_L> (fig. 14.8b) is the CSS encoder of |0_L> with one more wire:  *)
(* a qubit carrying a representative of logical X, chosen to avoid every pivot  *)
(* of the X checks, is put in T|+> = T H|0> and fanned out over that            *)
(* representative's support, which makes a|0_L> + b|1_L> from a|0> + b|1>.  For   *)
(* the 7-qubit code: 4 Hadamards, one T, 11 CNOTs, the book's counts.  It is not *)
(* fault tolerant, so the state is CHECKED rather than trusted.                 *)
(*                                                                              *)
(* TESTING it (fig. 14.9, sec. 13.4).  |A> is the +1 eigenstate of             *)
(* U = T X T^dagger, a Hermitian Clifford up to phase, so it can be measured    *)
(* with a cat state: a controlled-U from each cat qubit to its block qubit,      *)
(* broken into T^dagger = T S^dagger, CNOT, T as the book does, then the cat is   *)
(* read in X and its parity is the eigenvalue.  The run is kept when the cat's  *)
(* own checks pass and the eigenvalue is +1.  An EC step on the magic block      *)
(* comes first, because a single error on the block can fool the measurement.   *)
(*                                                                              *)
(* WHAT THE PAULI FRAME CANNOT DO.  A Pauli pushed through T is not a Pauli, so  *)
(* framePropagate refuses T, and the pi/8 exRec is counted -- bad sets, the     *)
(* book's own treatment in sec. 14.5.3 -- but not searched for malignant sets.  *)
(* Its correctness is checked instead on state vectors: the encoder makes        *)
(* exactly T_L|+_L>, and that state passes the test with certainty.             *)
(* ============================================================================ *)

$nonCliffordOps = {"T", "Tdg"};


(* ---- the magic state ---- *)

(* A logical X representative with no support on the pivots of the row-reduced X
   checks, so the fan-out never writes into a qubit the checks still need in |0>. *)
magicWire[a_Association] := magicWire[a] = With[{n = a["Qubits"]},
    With[{basis = gf2Basis[Select[a["CheckMatrix"], Total[#[[n + 1 ;; 2 n]]] === 0 &][[All, 1 ;; n]]],
          xbar = First[codeLogicalVectors[a]["X"]][[1 ;; n]]},
        Flatten @ Position[gf2Reduce[basis, xbar], 1]
    ]
]

(* T|+> (or T^dagger|+> for "Tdg") on the first qubit of that support, fanned out,
   then the X checks added as in cssZeroInstructions.  The two commute: the fan-out
   touches no pivot, and each check's CNOTs copy a pivot that nothing else writes. *)
magicStateInstructions[a_Association, offset_Integer, gate_String : "T"] := With[
    {support = magicWire[a]},
    Join[
        {{"H", offset + First[support]}, {gate, offset + First[support]}},
        Table[{"CNOT", offset + First[support], offset + q}, {q, Rest[support]}],
        cssZeroInstructions[a, offset]
    ]
]


(* ---- the test ---- *)

(* The cat, its checks, the controlled-U's and the readout, on block qubits
   block + 1 .. block + n, cat qubits cat + 1 .. cat + n and one check qubit.  The
   readout is n heralds whose parity is the eigenvalue. *)
magicTestInstructions[a_Association, block_Integer, cat_Integer, check_Integer] := With[
    {n = a["Qubits"], cq = Range[cat + 1, cat + a["Qubits"]]},
    Join[
        catInstructions[cq, check, catCheckPairs[n], 1],
        Catenate @ Table[
            {{"Sdg", block + i}, {"T", block + i}, {"CNOT", cat + i, block + i}, {"T", block + i}},
            {i, n}
        ],
        Table[{"H", q}, {q, cq}],
        Table[{"MH", q}, {q, cq}]
    ]
]

(* The cat's checks read zero one by one; the readout's parity must be even, i.e.
   the eigenvalue +1. *)
magicTestHeraldChecks[a_Association] := With[{n = a["Qubits"], k = Length[catCheckPairs[a["Qubits"]]]},
    Join[
        Table[{{j}, {{1}}}, {j, k}],
        {{Range[k + 1, k + n], {ConstantArray[1, n]}}}
    ]
]


(* ---- a state-vector check, for small circuits only ---- *)

(* The state an instruction list leaves, from |0...0>, as a vector over the 2^nq
   basis states, qubit 1 most significant.  Resets at the start are no-ops;
   measurements are not supported.  It exists to check what the frame cannot. *)
instructionStateVector[instr_List, nq_Integer] := Module[
    {v = UnitVector[2^nq, 1] + 0. I, one, bit, idx0, idx1, apply1},
    bit[q_] := 2^(nq - q);
    apply1[m_, q_] := (
        idx0 = Select[Range[0, 2^nq - 1], BitAnd[#, bit[q]] === 0 &] + 1;
        idx1 = idx0 + bit[q];
        {v[[idx0]], v[[idx1]]} = {m[[1, 1]] v[[idx0]] + m[[1, 2]] v[[idx1]], m[[2, 1]] v[[idx0]] + m[[2, 2]] v[[idx1]]}
    );
    one = <|
        "H" -> {{1, 1}, {1, -1}} / Sqrt[2.], "X" -> {{0, 1}, {1, 0}}, "Z" -> {{1, 0}, {0, -1}},
        "Y" -> {{0, -I}, {I, 0}}, "S" -> {{1, 0}, {0, I}}, "Sdg" -> {{1, 0}, {0, -I}},
        "T" -> {{1, 0}, {0, Exp[I Pi / 4.]}}, "Tdg" -> {{1, 0}, {0, Exp[-I Pi / 4.]}}
    |>;
    Do[
        Switch[First[g],
            "R", Null,
            "CNOT", With[{c = bit[g[[2]]], t = bit[g[[3]]]},
                v = v[[Table[If[BitAnd[i, c] =!= 0, BitXor[i, t], i], {i, 0, 2^nq - 1}] + 1]]],
            "CZ", With[{c = bit[g[[2]]], t = bit[g[[3]]]},
                v = v Table[If[BitAnd[i, c] =!= 0 && BitAnd[i, t] =!= 0, -1, 1], {i, 0, 2^nq - 1}]],
            _, apply1[one[First[g]], g[[2]]]
        ],
        {g, instr}
    ];
    v
]
