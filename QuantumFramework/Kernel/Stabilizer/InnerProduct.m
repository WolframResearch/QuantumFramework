Package["Wolfram`QuantumFramework`"]

PackageExport[StabilizerInnerProduct]
PackageExport[StabilizerExpectation]



(* ============================================================================ *)
(* Phase 4 \[Dash] inner products and expectation values via direct vector       *)
(* materialization. Phase 5+ may add the closed-form O(n^3) algorithm of        *)
(* Garc\[IAcute]a-Markov-Cross 2012 (arxiv:1210.6646) -- TODO: add the GMC      *)
(* algorithm for n > 8 where state materialization OOMs.                        *)
(* ============================================================================ *)


(* ============================================================================ *)
(* StabilizerInnerProduct[\[Psi], \[Phi]]: \[LeftAngleBracket]\[Psi]|\[Phi]\[RightAngleBracket]  *)
(*                                                                              *)
(* For two stabilizer states |\[Psi]>, |\[Phi]>:                                *)
(*   - If S_\[Psi] and S_\[Phi] disagree on any sign, returns 0.                *)
(*   - Otherwise returns 2^(-s/2) where s = n - dim(S_\[Psi] \[Intersection] S_\[Phi]).         *)
(*                                                                              *)
(* Phase 4 implementation: materialize state vectors and compute directly.      *)
(* Phase 5+ will use the GarMarCro12 \[Section]3 closed-form O(n^3) algorithm.  *)
(* ============================================================================ *)

StabilizerInnerProduct[psi_PauliStabilizer ? PauliStabilizerQ, phi_PauliStabilizer ? PauliStabilizerQ] /;
    psi["Qubits"] == phi["Qubits"] := Module[{vPsi, vPhi},
    vPsi = psi["State"]["StateVector"];
    vPhi = phi["State"]["StateVector"];
    Conjugate[vPsi] . vPhi
]


(* StabilizerFrame inner product: Sum_ij conj(c_i) c_j <psi_i | psi_j> *)
StabilizerInnerProduct[fA_StabilizerFrame, fB_StabilizerFrame] /;
    fA["Qubits"] == fB["Qubits"] := Module[{vA, vB},
    vA = fA["StateVector"];
    vB = fB["StateVector"];
    Conjugate[vA] . vB
]


(* Mixed cases *)
StabilizerInnerProduct[psi_PauliStabilizer ? PauliStabilizerQ, fB_StabilizerFrame] :=
    StabilizerInnerProduct[StabilizerFrame[psi], fB]

StabilizerInnerProduct[fA_StabilizerFrame, phi_PauliStabilizer ? PauliStabilizerQ] :=
    StabilizerInnerProduct[fA, StabilizerFrame[phi]]


(* ============================================================================ *)
(* StabilizerExpectation[ps, P]: \[LeftAngleBracket]\[Psi]| P |\[Psi]\[RightAngleBracket] *)
(*                                                                              *)
(* For stabilizer state |\[Psi]> with stabilizer group S, and a Pauli string P:  *)
(*   - If P or -P is in S \[Rule] expectation = +1 or -1 respectively           *)
(*   - Otherwise (P anticommutes with at least one stabilizer) \[Rule] 0        *)
(*                                                                              *)
(* Reference: AarGot04 \[Section]3, OngoingProjects/Stabilizer/paulistabilizer-source-audit.md \[Section]6.5 *)
(* ============================================================================ *)

(* Pauli string -> binary symplectic vector + sign *)
pauliStringToVec[s_String] := Module[{sign, body, n, xz},
    {sign, body} = If[StringStartsQ[s, "-"], {-1, StringDrop[s, 1]}, {1, s}];
    n = StringLength[body];
    xz = Replace[Characters[body],
        {"I" -> {0, 0}, "X" -> {1, 0}, "Y" -> {1, 1}, "Z" -> {0, 1}}, {1}];
    {sign, n, Join[xz[[All, 1]], xz[[All, 2]]]}
]


StabilizerExpectation[ps_PauliStabilizer ? PauliStabilizerQ, P_String] /;
    StringMatchQ[P, RegularExpression["^-?[IXYZ]+$"]] := Module[{sign, n, pVec, gens, omega, anticomm},
    {sign, n, pVec} = pauliStringToVec[P];
    If[n != ps["Qubits"],
        Message[StabilizerExpectation::dim, n, ps["Qubits"]];
        Return[$Failed]
    ];
    (* Check anticommute with any stabilizer *)
    omega = ArrayFlatten[{
        {ConstantArray[0, {n, n}], IdentityMatrix[n]},
        {IdentityMatrix[n], ConstantArray[0, {n, n}]}
    }];
    gens = ps["Matrix"][[ps["GeneratorCount"] + 1 ;; 2 ps["GeneratorCount"]]];
    anticomm = Mod[gens . omega . pVec, 2];
    If[AnyTrue[anticomm, # =!= 0 &],
        Return[0]   (* P anticommutes with some stabilizer => <P> = 0 *)
    ];
    (* P commutes with all stabilizers. P is in the closure (N(S)) but the SIGN of <P>
       depends on whether P or -P is in the stabilizer group, which requires tracking
       the i-factor from Y = iXZ when decomposing P over the generators.
       PHASE 4 v1: direct vector materialization for n <= 8 (always correct).
       TODO Phase 5+: replace with AG closed-form algorithm tracking the cumulative
       phase via agPhase[x1,z1,x2,z2] over the linear-combo recovered by LinearSolve. *)
    With[{vec = ps["State"]["StateVector"], pauliMat = pauliStringMatrix[P]},
        Re[Conjugate[vec] . (pauliMat . vec)]
    ]
]

StabilizerExpectation::dim = "Pauli string `1` qubits does not match PauliStabilizer `2` qubits."


(* Helper: build the matrix form of a Pauli string for fallback. *)
pauliStringMatrix[s_String] := Module[{sign, body, mats},
    {sign, body} = If[StringStartsQ[s, "-"], {-1, StringDrop[s, 1]}, {1, s}];
    mats = Replace[Characters[body], {
        "I" -> PauliMatrix[0],
        "X" -> PauliMatrix[1],
        "Y" -> PauliMatrix[2],
        "Z" -> PauliMatrix[3]
    }, {1}];
    sign * KroneckerProduct @@ mats
]
