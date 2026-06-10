Package["Wolfram`QuantumFramework`"]



(* ============================================================================ *)
(* Clifford-generator gate updates on the AG tableau.                           *)
(* References: AarGot04 \[Section]3 (canonical update rules), Biswas24          *)
(* \[Section]3 (pedagogical derivation), PatGuh26 \[Section]3.3 (Karnaugh-map  *)
(* derivation).                                                                 *)
(*                                                                              *)
(* GlobalPhase contract: every rule below builds a fresh association without   *)
(* forwarding "GlobalPhase". This is intentional. The phase factor a Clifford  *)
(* gate U picks up depends not just on U but on the input state's tableau     *)
(* (e.g., Z\[VerticalSeparator]0\[RightAngleBracket] = \[VerticalSeparator]0\[RightAngleBracket] vs Z\[VerticalSeparator]1\[RightAngleBracket] = -\[VerticalSeparator]1\[RightAngleBracket]; same Z gate, different    *)
(* phase). Recovering the new GlobalPhase exactly requires materializing the   *)
(* state vector at every gate update, which is O(2^n) and defeats the AG      *)
(* complexity advantage.                                                        *)
(*                                                                              *)
(* Contract: PauliStabilizer[qs][...gates...]["State"] is correct *up to a     *)
(* global phase*, not exactly. For exact equality, re-run the constructor on   *)
(* the post-gate state: PauliStabilizer[qs[gate]]["State"] === qs[gate].       *)
(* ============================================================================ *)

(* H_j: swap X[j]<->Z[j]; flip sign on rows where x[j] == z[j] == 1 *)
ps_PauliStabilizer["H", j_Integer] := With[{t = ps["Tableau"]},
    PauliStabilizer[<|
        "Signs" -> MapIndexed[If[t[[1, j, #2[[1]]]] == t[[2, j, #2[[1]]]] == 1, - #, #] &, ps["Signs"]],
        "Tableau" -> ReplacePart[t, Thread[{{1, j}, {2, j}} -> Extract[t, {{2, j}, {1, j}}]]]
    |>]
]

(* S_j: Z[j] := Z[j] XOR X[j]; flip sign on rows where x[j] == z[j] == 1 *)
ps_PauliStabilizer["S", j_Integer] := With[{t = ps["Tableau"]},
    PauliStabilizer[<|
        "Signs" -> MapIndexed[If[t[[1, j, #2[[1]]]] == t[[2, j, #2[[1]]]] == 1, - #, #] &, ps["Signs"]],
        "Tableau" -> ReplacePart[t, {2, j} -> MapThread[BitXor, Extract[ps["Tableau"], {{1, j}, {2, j}}]]]
    |>]
]

(* S^\[Dagger]_j: same Z-update; sign flip when x[j] == 1 and z[j] == 0 *)
ps_PauliStabilizer[SuperDagger["S"], j_Integer] := With[{t = ps["Tableau"]},
    PauliStabilizer[<|
        "Signs" -> MapIndexed[If[t[[1, j, #2[[1]]]] == 1 - t[[2, j, #2[[1]]]] == 1, - #, #] &, ps["Signs"]],
        "Tableau" -> ReplacePart[t, {2, j} -> MapThread[BitXor, Extract[ps["Tableau"], {{1, j}, {2, j}}]]]
    |>]
]

(* CNOT j -> k: X[k] := X[k] XOR X[j]; Z[j] := Z[j] XOR Z[k]; sign flip when AG conditions met *)
ps_PauliStabilizer["CNOT" | "CX", j_Integer, k_Integer] := With[{t = ps["Tableau"]},
    PauliStabilizer[<|
        "Signs" -> MapIndexed[If[t[[1, j, #2[[1]]]] == t[[2, k, #2[[1]]]] == 1 && t[[1, k, #2[[1]]]] == t[[2, j, #2[[1]]]], - #, #] &, ps["Signs"]],
        "Tableau" -> ReplacePart[t, {
            {1, k} -> MapThread[BitXor, Extract[t, {{1, j}, {1, k}}]],
            {2, j} -> MapThread[BitXor, Extract[t, {{2, j}, {2, k}}]]
        }]
    |>]
]

(* X_j: phase XOR with Z[j] row -- tableau unchanged *)
ps_PauliStabilizer["X", j_Integer] := PauliStabilizer[<|"Phase" -> BitXor[ps["Phase"], ps["Z"][[j]]], "Tableau" -> ps["Tableau"]|>]

(* Y_j: phase XOR with X[j] XOR Z[j] *)
ps_PauliStabilizer["Y", j_Integer] := PauliStabilizer[<|"Phase" -> BitXor[ps["Phase"], ps["X"][[j]], ps["Z"][[j]]], "Tableau" -> ps["Tableau"]|>]

(* Z_j: phase XOR with X[j] row *)
ps_PauliStabilizer["Z", j_Integer] := PauliStabilizer[<|"Phase" -> BitXor[ps["Phase"], ps["X"][[j]]], "Tableau" -> ps["Tableau"]|>]

(* CZ = H_k . CNOT(j,k) . H_k *)
ps_PauliStabilizer["CZ", j_Integer, k_Integer] := ps["H", k]["CNOT", j, k]["H", k]

(* SWAP via column permutation *)
ps_PauliStabilizer["SWAP", j_Integer, k_Integer] := ps["PermuteQudits", Cycles[{{j, k}}]]


(* ============================================================================ *)
(* V = sqrt(X) and V^\[Dagger] -- defined via PauliStabilizer string composition *)
(* ============================================================================ *)

ps_PauliStabilizer["V", j_Integer] := PauliStabilizer[{"-Y"}, {"X"}]["PadLeft", j] @ ps
ps_PauliStabilizer[SuperDagger["V"], j_Integer] := PauliStabilizer[{"Y"}, {"X"}]["PadLeft", j] @ ps


(* ============================================================================ *)
(* Permute (rows) and PermuteQudits (columns)                                   *)
(* ============================================================================ *)

ps_PauliStabilizer["Permute", perm_] := With[{n = ps["Qudits"]},
    PauliStabilizer[<|
        "Signs" -> Catenate[Permute[#, perm] & /@ Partition[ps["Signs"], n]],
        "Tableau" -> Map[Map[Catenate[Permute[#, perm] & /@ Partition[#, n]] &, Permute[#, perm]] &, ps["Tableau"]]
    |>]
]

ps_PauliStabilizer["PermuteQudits", perm_] := With[{n = ps["Qudits"]},
    PauliStabilizer[<|
        "Signs" -> ps["Signs"],
        "Tableau" -> Map[Permute[#, perm] &, ps["Tableau"]]
    |>]
]


(* ============================================================================ *)
(* Dagger / Inverse                                                             *)
(* ============================================================================ *)

PauliStabilizer::singular = "Cannot compute Dagger: symplectic matrix is singular mod 2 (this happens when the input PauliStabilizer was built from a stabilizer-only list whose Reverse-padded destabilizers coincide with the stabilizer rows). Use a circuit-built fixture or supply explicit destabilizers."

(* Notes on Dagger:                                                            *)
(*   (1) the inner inverse-matrix PS uses all-+1 signs so the composition's   *)
(*       phase term is well-defined and involutive on the matrix component   *)
(*       for AG-canonical inputs.                                              *)
(*   (2) a singular symplectic matrix (Det == 0 mod 2) emits                  *)
(*       PauliStabilizer::singular with $Failed instead of cascading errors.  *)
ps_PauliStabilizer["Dagger" | "Inverse"] := Block[{n, m, det, mat},
    n = ps["GeneratorCount"];
    m = ps["Matrix"];
    det = Mod[Det[m], 2];
    If[det === 0,
        Message[PauliStabilizer::singular];
        Return[$Failed]
    ];
    mat = Inverse[m, Modulus -> 2];
    PauliStabilizer @ <|
        "Phase" -> BitXor[
            ps["Phase"],
            ps[PauliStabilizer[<|
                "Matrix" -> mat,
                "Signs" -> ConstantArray[1, 2 n]
            |>]]["Phase"]
        ],
        "Matrix" -> mat
    |>
]


(* ============================================================================ *)
(* Pad to a target qubit count                                                  *)
(* ============================================================================ *)

ps_PauliStabilizer["PadRight", n_] := QuantumTensorProduct[ps, PauliStabilizer[Max[n - ps["Qubits"], 0]]]
ps_PauliStabilizer["PadLeft", n_]  := QuantumTensorProduct[PauliStabilizer[Max[n - ps["Qubits"], 0]], ps]


(* ============================================================================ *)
(* Non-Clifford boundary: P[\[Theta]] / T / T^\[Dagger]                         *)
(* Returns a `StabilizerFrame` (GarMar15 \[Section]3) that closes under further *)
(* Clifford updates.                                                            *)
(* ============================================================================ *)

ps_PauliStabilizer["P"[phase_], j_Integer] := With[{c = Exp[I phase / 2]},
    StabilizerFrame[{
        {(1 + c) / 2, ps},
        {(1 - c) / 2, ps["Z", j]}
    }]
]

ps_PauliStabilizer["T", j_Integer] := ps["P"[Pi / 2], j]
ps_PauliStabilizer[SuperDagger["T"], j_Integer] := ps["P"[- Pi / 2], j]
