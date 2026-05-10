Package["Wolfram`QuantumFramework`"]



(* ============================================================================ *)
(* Hybrid interop UpValues.                                                     *)
(*                                                                              *)
(* Cross-head dispatch so QF measurement / channel operators consume a          *)
(* PauliStabilizer or StabilizerFrame natively without forcing the caller to    *)
(* materialize via ps["State"] (which costs O(2^n)). Mirrors the existing       *)
(* PauliStabilizer integration UpValues in Stabilizer/Conversions.m.            *)
(*                                                                              *)
(* Design rationale: UpValues attached to PauliStabilizer / StabilizerFrame --  *)
(* not a Picture flag, not a QuantumBasis wrapper. Both alternatives would      *)
(* route through QuantumBasis machinery (KroneckerProduct[Output, Input] +      *)
(* MatrixInverse[ReducedMatrix]) and pay O(2^n)..O(8^n), defeating the          *)
(* formalism's O(n^2) advantage. UpValues stay in the tableau when they can.    *)
(*                                                                              *)
(* Dispatch ladder for qmo[ps]:                                                 *)
(*   1. QMO basis is a Pauli string  -> ps["M", pauli], stays in tableau.       *)
(*   2. QMO basis is non-Pauli       -> emits PauliStabilizer::nonpaulibasis    *)
(*                                       and falls back to                      *)
(*                                       PauliStabilizerApply[                  *)
(*                                         QuantumCircuitOperator[qmo], ps].   *)
(* ============================================================================ *)


(* ============================================================================ *)
(* Helper: extract a Pauli string label from a QMO, when possible.              *)
(*                                                                              *)
(* Returns the string for a Pauli-eigenbasis measurement (e.g. "XZZXI"),        *)
(* or Missing[] otherwise. Used to gate the fast path.                          *)
(* ============================================================================ *)

PackageScope[stabilizerPauliLabelFromQMO]
PackageScope[stabilizerPauliFromMatrix]
PackageScope[$stabilizerPauliMatrixSearchMaxQubits]

(* Recognized label forms:                                                    *)
(*   1. String matching ^-?[IXYZ]+$            -> direct Pauli string         *)
(*   2. -Superscript[X|Y|Z|I, CircleTimes[m]]  -> "-XXXX..." (m copies)       *)
(*   3. Superscript[X|Y|Z|I, CircleTimes[m]]   -> "XXXX..."  (m copies)       *)
(*   4. Times[-1, str] with str a Pauli string -> "-" <> str                  *)
(* Matrix-iteration fallback: for n <= cap qubits, also iterates 4^n*{+1,-1}  *)
(* Pauli candidates against the QMO's MatrixRepresentation. This catches QMOs *)
(* built directly from explicit matrices (QuantumOperator[matrix, ...]) where *)
(* the symbolic Label is None. For larger n the cost is exponential in n;     *)
(* the cap keeps the search bounded (4^4 * 2 = 512 candidates).               *)
(* ============================================================================ *)

(* Cap on the matrix-iteration detector. PackageScope so tests/users may     *)
(* override; default 4 keeps the search cost ~O(4^n) bounded.                *)
$stabilizerPauliMatrixSearchMaxQubits = 4;


(* ============================================================================ *)
(* Helper: build a 2^n x 2^n matrix for a Pauli string.                        *)
(*                                                                              *)
(* Internal helper for the matrix-iteration detector (handles the n=1          *)
(* KroneckerProduct edge case).                                                 *)
(* ============================================================================ *)

stabilizerPauliMatrixFromString[s_String] := Module[{sign, body, mats},
    {sign, body} = If[StringStartsQ[s, "-"], {-1, StringDrop[s, 1]}, {1, s}];
    mats = Replace[Characters[body], {
        "I" -> PauliMatrix[0],
        "X" -> PauliMatrix[1],
        "Y" -> PauliMatrix[2],
        "Z" -> PauliMatrix[3]
    }, {1}];
    sign * If[Length[mats] == 1, First[mats], KroneckerProduct @@ mats]
]


(* ============================================================================ *)
(* Matrix-iteration Pauli detector.                                             *)
(*                                                                              *)
(* Given an explicit 2^n by 2^n matrix and a qubit count n, iterate over all   *)
(* signed Pauli strings of length n and return the first one whose matrix      *)
(* equals the input. Cost: O(4^n * 2 * 4^n) = O(16^n) for the comparisons.    *)
(* Capped at $stabilizerPauliMatrixSearchMaxQubits (default 4).                *)
(*                                                                              *)
(* Returns:                                                                     *)
(*   - String form like "X", "-XYZ", etc., if a match is found.                *)
(*   - Missing["TooManyQubits"] if n > cap.                                    *)
(*   - Missing["DimMismatch"] if the matrix doesn't have shape {2^n, 2^n}.    *)
(*   - Missing["NoPauliMatch"] if no candidate matches.                        *)
(* ============================================================================ *)

stabilizerPauliFromMatrix[mat_, n_Integer ? Positive] /; n > $stabilizerPauliMatrixSearchMaxQubits :=
    Missing["TooManyQubits"]

stabilizerPauliFromMatrix[mat_, n_Integer ? Positive] := Module[{
    dim = 2^n, normalized, candidates, result
},
    If[!MatrixQ[Normal[mat]] || Dimensions[Normal[mat]] =!= {dim, dim},
        Return[Missing["DimMismatch"]]
    ];

    normalized = Normal[mat];
    candidates = Tuples[{"I", "X", "Y", "Z"}, n];

    result = Catch[
        Do[
            With[{str = StringJoin[c]},
                With[{m = stabilizerPauliMatrixFromString[str]},
                    If[normalized === Normal[m], Throw[str]];
                    If[normalized === Normal[-m], Throw["-" <> str]]
                ]
            ],
            {c, candidates}
        ];
        $unmatched
    ];

    If[StringQ[result], result, Missing["NoPauliMatch"]]
]


stabilizerPauliLabelFromQMO[qmo_QuantumMeasurementOperator] := Module[{
    op, label, matched, matrix, target, n, matSearch
},
    op = qmo["Operator"];
    label = op["Label"];

    (* Form 1: plain Pauli string. *)
    If[StringQ[label] && StringMatchQ[label, RegularExpression["^-?[IXYZ]+$"]],
        Return[label]
    ];

    (* Forms 2/3/4: pattern-replace common label expressions. *)
    matched = Replace[label, {
        Times[-1, Superscript[letter : "X" | "Y" | "Z" | "I", CircleTimes[m_Integer ? Positive]]] :>
            "-" <> StringJoin[ConstantArray[letter, m]],
        Superscript[letter : "X" | "Y" | "Z" | "I", CircleTimes[m_Integer ? Positive]] :>
            StringJoin[ConstantArray[letter, m]],
        Times[-1, str_String /; StringMatchQ[str, RegularExpression["^[IXYZ]+$"]]] :>
            "-" <> str,
        _ :> $unmatchedPauliLabel
    }];

    If[matched =!= $unmatchedPauliLabel && StringQ[matched],
        Return[matched]
    ];

    (* Matrix-iteration fallback for n <= cap. We require a square              *)
    (* 2^n x 2^n matrix; non-square shapes (e.g. computational-basis projector  *)
    (* stacks of shape {2^(n+1), 2^n}) signal a multi-Kraus QMO and fall        *)
    (* through to the legacy path.                                              *)
    target = qmo["InputOrder"];
    n = Length[target];
    If[1 <= n <= $stabilizerPauliMatrixSearchMaxQubits,
        matrix = op["MatrixRepresentation"];
        If[MatrixQ[Normal[matrix]] && Dimensions[Normal[matrix]] === {2^n, 2^n},
            matSearch = stabilizerPauliFromMatrix[matrix, n];
            If[StringQ[matSearch], Return[matSearch]]
        ]
    ];

    Missing["NonPauliBasis"]
]


(* ============================================================================ *)
(* PauliStabilizer::nonpaulibasis -- fired when a non-Pauli QMO/channel acts    *)
(* on a stabilizer-form input. Currently falls back to the generic              *)
(* PauliStabilizerApply circuit-conversion path. A future StabilizerFrame       *)
(* decomposition would keep the cost at O(rank * n^2) rather than O(2^n).      *)
(* ============================================================================ *)

PauliStabilizer::nonpaulibasis = "Hybrid interop: the measurement/channel basis is non-Pauli; falling back to the generic circuit path."


(* ============================================================================ *)
(* QuantumMeasurementOperator[qmo][ps_PauliStabilizer] (UpValue on              *)
(*   PauliStabilizer)                                                           *)
(*                                                                              *)
(* Fast path: when the QMO's Operator label parses as a Pauli string, route    *)
(* directly to ps["M", pauliString] which is the existing AG measurement       *)
(* primitive (Stabilizer/PauliMeasure.m). Stays in the tableau, O(n^2).        *)
(*                                                                              *)
(* Fallback: emit ::nonpaulibasis info message and use the generic path        *)
(* PauliStabilizerApply[QuantumCircuitOperator[qmo], ps] which converts        *)
(* the QMO to a circuit and folds gates over the tableau.                      *)
(* ============================================================================ *)

qmo_QuantumMeasurementOperator[ps_PauliStabilizer ? ConcretePauliStabilizerQ] ^:= Module[{label},
    label = stabilizerPauliLabelFromQMO[qmo];
    If[StringQ[label],
        ps["M", label],
        Message[PauliStabilizer::nonpaulibasis];
        PauliStabilizerApply[QuantumCircuitOperator[qmo], ps]
    ]
]


(* StabilizerFrame: same dispatch shape but no native "M" method on a frame    *)
(* yet. Materialize the frame and apply the qmo on the materialized state.     *)

qmo_QuantumMeasurementOperator[sf_StabilizerFrame] ^:= (
    Message[PauliStabilizer::nonpaulibasis];
    qmo[sf["State"]]
) /; StabilizerFrameQ[sf]


(* ============================================================================ *)
(* QuantumChannel[qc][ps_PauliStabilizer] / [sf_StabilizerFrame]                *)
(*                                                                              *)
(* Detect named Pauli channels (BitFlip, PhaseFlip, BitPhaseFlip,              *)
(* Depolarizing) and return a probabilistic-mixture list                        *)
(* {{probability, ps_after_pauli}, ...} where each ps_after_pauli is the        *)
(* original ps with the corresponding Pauli gate applied via tableau update     *)
(* (cost O(n) per branch). This is the natural form for tableau-level Clifford- *)
(* channel application; the user can post-process into a mixed state via QF's   *)
(* QuantumState mixture semantics, or sample stochastically.                    *)
(*                                                                              *)
(* Non-Clifford channels (AmplitudeDamping, PhaseDamping, GeneralizedAmplitude- *)
(* Damping, ResetError) still fall back to dense materialization with the       *)
(* ::nonpaulibasis info message; their Kraus operators are not all Paulis.      *)
(* ============================================================================ *)


(* Helper: extract the Clifford-channel Pauli mixture from qc, when            *)
(* applicable. Returns a list of {probability, pauli_string, qubit_index}       *)
(* triples (one per Kraus branch), or Missing[] for non-Clifford channels.     *)
(* The qubit_index is the channel's input order (single-qubit channels only).  *)

PackageScope[stabilizerCliffordChannelMixture]

stabilizerCliffordChannelMixture[qc_QuantumChannel] := Module[{label, target},
    label = qc["Label"];
    target = qc["InputOrder"];
    (* Single-qubit named Pauli channels. Multi-qubit / non-named channels are  *)
    (* handled by the dense fallback below.                                     *)
    If[Length[target] =!= 1, Return[Missing["MultiQubitChannel"]]];
    target = First[target];
    Switch[label,
        "BitFlip"[_],         {{1 - label[[1]], "I", target}, {label[[1]], "X", target}},
        "PhaseFlip"[_],       {{1 - label[[1]], "I", target}, {label[[1]], "Z", target}},
        "BitPhaseFlip"[_],    {{1 - label[[1]], "I", target}, {label[[1]], "Y", target}},
        "\[CapitalDelta]"[_], (* Depolarizing[p] -- internal label is CapitalDelta[p] *)
            {
                {1 - 3 label[[1]] / 4, "I", target},
                {label[[1]] / 4, "X", target},
                {label[[1]] / 4, "Y", target},
                {label[[1]] / 4, "Z", target}
            },
        _, Missing["NotClifford"]
    ]
]


qc_QuantumChannel[ps_PauliStabilizer ? ConcretePauliStabilizerQ] ^:= Module[{mixture},
    mixture = stabilizerCliffordChannelMixture[qc];
    If[ListQ[mixture],
        (* Fast path: tableau-level Pauli-mixture application. Each branch is a  *)
        (* (probability, post-state) pair; "I" branch leaves ps unchanged.       *)
        {#1, If[#2 === "I", ps, ps[#2, #3]]} & @@@ mixture,
        (* Fallback: non-Clifford channel; materialize. *)
        Message[PauliStabilizer::nonpaulibasis];
        qc[ps["State"]]
    ]
]

qc_QuantumChannel[sf_StabilizerFrame] ^:= (
    Message[PauliStabilizer::nonpaulibasis];
    qc[sf["State"]]
) /; StabilizerFrameQ[sf]
