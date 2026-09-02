(* ::Package:: *)

Package["Wolfram`QuantumFramework`QEC`"]

PackageScope[codeReductionGates]
PackageScope[invertGate]


(* ============================================================================ *)
(* The encoder.                                                                 *)
(*                                                                              *)
(* Gottesman, QECC book sec. 6.4.1 (Procedure 6.6) with Table 6.1.  The m        *)
(* generators are the columns of a 2n x m symplectic matrix G, rows 1..n the X   *)
(* part and rows n+1..2n the Z part.  Left multiplication by a Clifford gate is  *)
(* a row operation on G (H swaps the two rows of a qubit, S adds the X row into  *)
(* the Z row, CNOT and CZ combine rows across qubits, SWAP exchanges qubits);    *)
(* relabelling generators is a free column operation.  Reduce G to the initial   *)
(* all-Z stabilizer form while recording the gates, and the encoder is the       *)
(* inverse of that sequence.                                                     *)
(* ============================================================================ *)


(* ---- applying Paulis and gate lists to an engine state ---- *)

(* The phase of v is dropped: a stabilizer state carries no global phase, and a
   Pauli's overall factor is exactly that. *)
applyPauliVector[ps_Wolfram`QuantumFramework`PauliStabilizer, v_List] := Module[{n = pauliQubits[v]},
    Fold[
        Replace[#2, {{"I", _} -> #1, {letter_, q_} :> #1[letter, q]}] &,
        ps,
        Table[{Lookup[<|{0, 0} -> "I", {1, 0} -> "X", {1, 1} -> "Y", {0, 1} -> "Z"|>, Key[{v[[q]], v[[n + q]]}]], q}, {q, n}]
    ]
]

applyPauliVector[ps_Wolfram`QuantumFramework`PauliStabilizer, s_String] := applyPauliVector[ps, QECPauliVector[s]]

(* The engine's compiled bulk fold, not a gate-by-gate object fold.  Every gate the
   reduction emits is an arrow `name -> order` over {H, S, X, Y, Z, CNOT, CZ, SWAP},
   which is exactly the spec shape "ApplyCircuit" takes, so this is a straight
   delegation.  Folding one object-gate at a time pays object wrap/unwrap and measures
   ~25x slower than the compiled route; the difference is invisible on a 7-qubit
   encoder and decisive once a noise model runs the same circuit thousands of times. *)
applyGates[ps_Wolfram`QuantumFramework`PauliStabilizer, {}] := ps

applyGates[ps_Wolfram`QuantumFramework`PauliStabilizer, gates_List] := ps["ApplyCircuit", gates]


(* ---- reduction ---- *)

codeReductionGates[a_Association] := Module[
    {n = a["Qubits"], m = codeStabilizerCount[a], g, gates, hh, ss, cnot, cz, swap, coladd, pivot},

    g = Transpose[a["CheckMatrix"]];
    (* Internal`Bag rather than AppendTo: the gate list is built by a greedy sweep and
       AppendTo copies it on every push, which is the pattern the Stabilizer subsystem
       moved away from in its 2026 idiom audit. *)
    gates = Internal`Bag[];

    hh[q_] := (g[[{q, n + q}]] = g[[{n + q, q}]]; Internal`StuffBag[gates, "H" -> q]);
    ss[q_] := (g[[n + q]] = Mod[g[[n + q]] + g[[q]], 2]; Internal`StuffBag[gates, "S" -> q]);
    cnot[c_, t_] := (
        g[[t]] = Mod[g[[t]] + g[[c]], 2];
        g[[n + c]] = Mod[g[[n + c]] + g[[n + t]], 2];
        Internal`StuffBag[gates, "CNOT" -> {c, t}]
    );
    cz[c_, t_] := (
        g[[n + t]] = Mod[g[[n + t]] + g[[c]], 2];
        g[[n + c]] = Mod[g[[n + c]] + g[[t]], 2];
        Internal`StuffBag[gates, "CZ" -> {c, t}]
    );
    swap[u_, v_] := (
        g[[{u, v}]] = g[[{v, u}]];
        g[[{n + u, n + v}]] = g[[{n + v, n + u}]];
        Internal`StuffBag[gates, "SWAP" -> {u, v}]
    );
    coladd[src_, dst_] := (g[[All, dst]] = Mod[g[[All, dst]] + g[[All, src]], 2]);

    Do[
        pivot = SelectFirst[Range[p, n], g[[#, p]] === 1 &];
        If[ MissingQ[pivot],
            pivot = SelectFirst[Range[p, n], g[[n + #, p]] === 1 &];
            hh[pivot]
        ];
        If[pivot =!= p, swap[p, pivot]];
        Do[If[i =!= p && g[[i, p]] === 1, cnot[p, i]], {i, n}];
        If[g[[n + p, p]] === 1, ss[p]];
        Do[If[i =!= p && g[[n + i, p]] === 1, cz[p, i]], {i, n}];
        Do[If[q =!= p && g[[p, q]] === 1, coladd[p, q]], {q, m}],
        {p, m}
    ];

    Do[hh[q], {q, m}];

    Internal`BagPart[gates, All]
]

(* H, CNOT, CZ and SWAP are their own inverses; S^-1 = S^3. *)
invertGate["H" -> q_] := {"H" -> q}
invertGate["CNOT" -> qs_] := {"CNOT" -> qs}
invertGate["CZ" -> qs_] := {"CZ" -> qs}
invertGate["SWAP" -> qs_] := {"SWAP" -> qs}
invertGate["S" -> q_] := {"S" -> q, "S" -> q, "S" -> q}


(* ---- the encoder itself ---- *)

(* The reduction run backwards, plus a Pauli fixup.  Inverting the reduction
   prepares a state stabilized by +-each generator; the fixup solves over GF(2)
   for the Pauli whose symplectic products with the completed generators match the
   pattern of wrong signs, and applying it flips exactly those. *)
codeEncodingGates[a_Association] := codeEncodingGates[a] = Module[
    {n = a["Qubits"], encoder, prepared, completed, system, wrong, fix},

    encoder = Flatten[invertGate /@ Reverse[codeReductionGates[a]], 1];
    prepared = applyGates[Wolfram`QuantumFramework`PauliStabilizer[n], encoder];

    completed = codeCompletedGenerators[a];
    If[completed === $Failed, Return[$Failed]];

    wrong = Boole[prepared["Expectation", #] === -1] & /@ (QECPauliString /@ completed);
    If[Total[wrong] === 0, Return[encoder]];

    (* Row i of the system is generator i with its halves swapped, so that
       system . p = the vector of symplectic products of p with the generators. *)
    system = With[{v = symplecticPart[#]}, Join[v[[n + 1 ;; 2 n]], v[[1 ;; n]]]] & /@ completed;
    fix = gf2Solve[system, wrong];
    If[fix === $Failed, Return[encoder]];

    Join[
        encoder,
        Flatten[Table[{If[fix[[q]] === 1, "X" -> q, Nothing], If[fix[[n + q]] === 1, "Z" -> q, Nothing]}, {q, n}]]
    ]
]

codeEncodingValidQ[a_Association] := With[
    {prepared = applyGates[Wolfram`QuantumFramework`PauliStabilizer[a["Qubits"]], codeEncodingGates[a]]},
    AllTrue[QECPauliString /@ codeVectors[a], prepared["Expectation", #] === 1 &]
]
