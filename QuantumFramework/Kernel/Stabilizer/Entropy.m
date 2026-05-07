Package["Wolfram`QuantumFramework`"]

PackageScope[stabilizerEntropy]



(* ============================================================================ *)
(* A.7 (2026-05-07): closed-form polynomial-time stabilizer-state entropy.      *)
(*                                                                              *)
(* Reference: Fattal, Cubitt, Yamamoto, Bravyi, Chuang 2004                    *)
(* (arxiv:quant-ph/0406168) Theorem 1: for an n-qubit stabilizer state psi    *)
(* and a partition A subset {1, ..., n}, the bipartite entanglement entropy   *)
(* across the cut is                                                            *)
(*                                                                              *)
(*   S(A) = rank_{F_2}(stabilizer generators restricted to A) - |A|            *)
(*                                                                              *)
(* where the "restriction to A" projects each generator's symplectic vector   *)
(* onto the qubits in A (keep columns A and (A + n) of the symplectic         *)
(* matrix). The rank is over F_2.                                              *)
(*                                                                              *)
(* Sanity checks (in units of log_2):                                          *)
(*   |00>          A={1}    stabs ZI, IZ -> {Z, I}  rank=1, |A|=1, S=0  (separable) *)
(*   |Bell>        A={1}    stabs XX, ZZ -> {X, Z}  rank=2, |A|=1, S=1  (max entangled) *)
(*   |GHZ_3>       A={1}    stabs XXX, ZZI, IZZ -> {X, Z, I}  rank=2, |A|=1, S=1 *)
(*   |Cluster_4>   A={1,2}  rank=3, |A|=2, S=1                              *)
(*                                                                              *)
(* This is `O(n^3)` (rank computation), independent of `|A|` -- closed form   *)
(* succeeds at any qubit count, unlike the Schmidt-rank path which OOMs at    *)
(* n ~ 12.                                                                     *)
(* ============================================================================ *)

stabilizerEntropy[ps_PauliStabilizer ? PauliStabilizerQ, partition_List] := Module[{
    n, gens, partitionA, indicesA, restricted, rankRestricted
},
    n = ps["Qubits"];
    gens = ps["Matrix"][[ps["GeneratorCount"] + 1 ;; 2 ps["GeneratorCount"]]];
    partitionA = Sort @ DeleteDuplicates @ partition;

    If[!AllTrue[partitionA, IntegerQ[#] && 1 <= # <= n &],
        Message[PauliStabilizer::partition, partition, n];
        Return[$Failed]
    ];

    (* Restrict each generator's symplectic vector to qubits in A. The bits  *)
    (* layout for a row is [x_1 ... x_n | z_1 ... z_n]; selecting qubit set  *)
    (* A means keeping columns A and (n + A).                                *)
    indicesA = Join[partitionA, partitionA + n];
    restricted = gens[[All, indicesA]];

    rankRestricted = MatrixRank[restricted, Modulus -> 2];

    rankRestricted - Length[partitionA]
]


(* Whole-list partition convenience: ps["Entropy", {1, 2}] *)
ps_PauliStabilizer["Entropy", partition_List] := stabilizerEntropy[ps, partition]

(* Single-qubit shortcut: ps["Entropy", q] *)
ps_PauliStabilizer["Entropy", q_Integer] := stabilizerEntropy[ps, {q}]


PauliStabilizer::partition = "Partition `1` is not a subset of {1, ..., `2`}."
