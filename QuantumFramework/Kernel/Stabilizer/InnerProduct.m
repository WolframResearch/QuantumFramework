Package["Wolfram`QuantumFramework`"]

PackageScope[stabilizerInnerProduct]
PackageScope[stabilizerExpectation]



(* ============================================================================ *)
(* Inner products and expectation values.                                       *)
(*                                                                              *)
(* Two paths are exposed via Method:                                            *)
(*   "Direct"     -- materialize state vectors (O(2^n); recovers complex phase) *)
(*   "ClosedForm" -- García-Markov-Cross 2012 §3 (O(n^3); magnitude only)       *)
(*                                                                              *)
(* Method-grade entry points:                                                   *)
(*   ps1["InnerProduct", ps2]                                                   *)
(*   ps["Expectation", pauli]                                                   *)
(* Implementations live here as PackageScope helpers; method dispatch is in     *)
(* Stabilizer/Properties.m and StabilizerFrame.m.                               *)
(* ============================================================================ *)


(* ============================================================================ *)
(* stabilizerInnerProduct[\[Psi], \[Phi]]: \[LeftAngleBracket]\[Psi]|\[Phi]\[RightAngleBracket]  *)
(*                                                                              *)
(* For two stabilizer states |\[Psi]>, |\[Phi]>:                                *)
(*   - If S_\[Psi] and S_\[Phi] disagree on any sign, returns 0.                *)
(*   - Otherwise returns 2^(-s/2) where s = n - dim(S_\[Psi] \[Intersection] S_\[Phi]).         *)
(*                                                                              *)
(* Default implementation: materialize state vectors and compute directly.     *)
(* ============================================================================ *)

(* Closed-form O(n^3) inner-product MAGNITUDE per GarMarCro12 §3.              *)
(*                                                                              *)
(* For two stabilizer states |psi>, |phi> on n qubits with stabilizer groups  *)
(* S_psi, S_phi:                                                                *)
(*   |<psi | phi>| = 0           if S_psi and S_phi disagree on any sign in    *)
(*                                S_psi cap S_phi.                              *)
(*                = 2^(-s/2)     otherwise, where s = n - dim(S_psi cap S_phi).*)
(*                                                                              *)
(* The COMPLEX phase of <psi | phi> requires Gaussian-sum bookkeeping over the *)
(* intersection generators (TODO). For the magnitude (and ±1 sign for         *)
(* perfectly-aligned states), the closed form is `O(n^3)` and works at any n.  *)
(*                                                                              *)
(* Available as `ps["InnerProduct", other, Method -> "ClosedForm"]`. The       *)
(* default `ps["InnerProduct", other]` uses direct-vector materialization      *)
(* (`O(2^n)` time + memory; recovers the full complex phase). For n > 8 where *)
(* materialization OOMs, use the closed-form Method.                            *)

stabilizerInnerProductClosedForm[psi_PauliStabilizer ? ConcretePauliStabilizerQ, phi_PauliStabilizer ? ConcretePauliStabilizerQ] /;
    psi["Qubits"] == phi["Qubits"] := Module[{
    n, gPsi, gPhi, phasePsi, phasePhi,
    joint, rankJoint, dimInt, s, nullVecs,
    intersectionSign
},
    n = psi["Qubits"];
    gPsi = psi["Matrix"][[psi["GeneratorCount"] + 1 ;; 2 psi["GeneratorCount"]]];
    gPhi = phi["Matrix"][[phi["GeneratorCount"] + 1 ;; 2 phi["GeneratorCount"]]];
    phasePsi = psi["Phase"][[psi["GeneratorCount"] + 1 ;; 2 psi["GeneratorCount"]]];
    phasePhi = phi["Phase"][[phi["GeneratorCount"] + 1 ;; 2 phi["GeneratorCount"]]];

    joint = Mod[Join[gPsi, gPhi], 2];
    rankJoint = MatrixRank[joint, Modulus -> 2];
    dimInt = 2 n - rankJoint;
    s = n - dimInt;

    If[dimInt == 0, Return[2^(-s/2)]];

    nullVecs = NullSpace[Transpose[joint], Modulus -> 2];

    intersectionSign[a_, b_] := Module[{phasePsiBit, phasePhiBit},
        phasePsiBit = Mod[
            a . phasePsi + Quotient[stabilizerRowSumAGPhase[gPsi, a, n], 2],
            2
        ];
        phasePhiBit = Mod[
            b . phasePhi + Quotient[stabilizerRowSumAGPhase[gPhi, b, n], 2],
            2
        ];
        phasePsiBit == phasePhiBit
    ];

    If[
        AllTrue[nullVecs, intersectionSign[Take[#, n], Drop[#, n]] &],
        2^(-s/2),
        0
    ]
]


(* Default: direct vector materialization. Practical for n <= ~8. *)
stabilizerInnerProduct[psi_PauliStabilizer ? PauliStabilizerQ, phi_PauliStabilizer ? PauliStabilizerQ, OptionsPattern[{Method -> "Direct"}]] /;
    psi["Qubits"] == phi["Qubits"] := Switch[OptionValue[Method],
    "ClosedForm",
        If[ConcretePauliStabilizerQ[psi] && ConcretePauliStabilizerQ[phi],
            stabilizerInnerProductClosedForm[psi, phi],
            Message[PauliStabilizer::expectationdim,
                "InnerProduct[..., Method -> \"ClosedForm\"] requires concrete (numeric) signs."];
            $Failed
        ],
    _,
        Module[{vPsi, vPhi},
            vPsi = psi["State"]["StateVector"];
            vPhi = phi["State"]["StateVector"];
            Conjugate[vPsi] . vPhi
        ]
]


(* StabilizerFrame inner product: Sum_ij conj(c_i) c_j <psi_i | psi_j> *)
stabilizerInnerProduct[fA_StabilizerFrame, fB_StabilizerFrame] /;
    fA["Qubits"] == fB["Qubits"] := Module[{vA, vB},
    vA = fA["StateVector"];
    vB = fB["StateVector"];
    Conjugate[vA] . vB
]


(* Mixed cases *)
stabilizerInnerProduct[psi_PauliStabilizer ? PauliStabilizerQ, fB_StabilizerFrame] :=
    stabilizerInnerProduct[StabilizerFrame[psi], fB]

stabilizerInnerProduct[fA_StabilizerFrame, phi_PauliStabilizer ? PauliStabilizerQ] :=
    stabilizerInnerProduct[fA, StabilizerFrame[phi]]


(* ============================================================================ *)
(* stabilizerExpectation[ps, P]: \[LeftAngleBracket]\[Psi]| P |\[Psi]\[RightAngleBracket] *)
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


(* Closed-form expectation via AG-phase i-factor tracking. Algorithm:         *)
(*   1. P anticommutes with any stabilizer => <P> = 0.                        *)
(*   2. P commutes with all stabilizers => decompose pVec over generators    *)
(*      via LinearSolve over F_2.                                              *)
(*   3. If LinearSolve fails (P in N(S) \ S, the larger normalizer minus the *)
(*      stabilizer subgroup), <P> = 0.                                        *)
(*   4. Otherwise compute the cumulative AG i-power of the F_2 row sum and    *)
(*      combine with the input sign and the per-row signs to get +-1.        *)

stabilizerExpectation[ps_PauliStabilizer ? PauliStabilizerQ, P_String] /;
    StringMatchQ[P, RegularExpression["^-?[IXYZ]+$"]] := Module[{
    sign, n, pVec, gens, signs, phases, omega, anticomm, coeffs,
    rowSumPhase, signProduct, agSignBit
},
    {sign, n, pVec} = pauliStringToVec[P];
    If[n != ps["Qubits"],
        Message[PauliStabilizer::expectationdim, n, ps["Qubits"]];
        Return[$Failed]
    ];
    omega = ArrayFlatten[{
        {ConstantArray[0, {n, n}], IdentityMatrix[n]},
        {IdentityMatrix[n], ConstantArray[0, {n, n}]}
    }];
    gens = ps["Matrix"][[ps["GeneratorCount"] + 1 ;; 2 ps["GeneratorCount"]]];
    signs = ps["StabilizerSigns"];
    phases = ps["Phase"][[ps["GeneratorCount"] + 1 ;; 2 ps["GeneratorCount"]]];
    anticomm = Mod[gens . omega . pVec, 2];
    If[AnyTrue[anticomm, # =!= 0 &],
        Return[0]   (* P anticommutes with some stabilizer => <P> = 0 *)
    ];
    (* P commutes with all stabilizers. Try to decompose pVec = sum c_i g_i mod 2. *)
    coeffs = Quiet @ Check[LinearSolve[Transpose[gens], pVec, Modulus -> 2], $Failed];
    If[coeffs === $Failed || !VectorQ[coeffs, IntegerQ],
        Return[0]   (* P is in N(S) \ S; orthogonal to ψ. *)
    ];
    (* AG cumulative i-power for the F_2 row sum, and the per-row sign        *)
    (* product. The total sign of <P> is targetSign * (-1)^(rowSumPhase/2) *  *)
    (*   prod_i signs[i]^c_i.                                                  *)
    rowSumPhase = stabilizerRowSumAGPhase[gens, coeffs, n];
    (* For valid stabilizer-group members, the cumulative phase mod 4 is in   *)
    (* {0, 2}. If odd, P is not in S (caught by LinearSolve normally; this is *)
    (* a defensive check).                                                     *)
    If[OddQ[rowSumPhase], Return[0]];
    agSignBit = Quotient[rowSumPhase, 2];
    signProduct = Times @@ (signs^coeffs);
    sign * (1 - 2 agSignBit) * signProduct
]

PauliStabilizer::expectationdim = "Pauli string `1` qubits does not match PauliStabilizer `2` qubits."


(* Helper: build the matrix form of a Pauli string for fallback.              *)
(* The 1-qubit case is handled explicitly because KroneckerProduct with a    *)
(* single argument throws KroneckerProduct::argmu and produces an            *)
(* unevaluated expression that breaks downstream `Re[...]` of the expectation.*)
pauliStringMatrix[s_String] := Module[{sign, body, mats},
    {sign, body} = If[StringStartsQ[s, "-"], {-1, StringDrop[s, 1]}, {1, s}];
    mats = Replace[Characters[body], {
        "I" -> PauliMatrix[0],
        "X" -> PauliMatrix[1],
        "Y" -> PauliMatrix[2],
        "Z" -> PauliMatrix[3]
    }, {1}];
    sign * If[Length[mats] == 1, First[mats], KroneckerProduct @@ mats]
]
