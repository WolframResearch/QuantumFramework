(* ::Package:: *)

Package["Wolfram`QuantumFramework`QEC`"]

PackageScope[codeSwapped]
PackageScope[codeCorrectableWeight]
PackageScope[codeDecoderToWeight]
PackageScope[codeDecoderReach]


(* ============================================================================ *)
(* Syndromes and decoding.                                                      *)
(*                                                                              *)
(* The syndrome of an error is the pattern of anticommutations with the checks,  *)
(* which over GF(2) is one matrix-vector product: with the halves of the check    *)
(* matrix swapped, syndrome = E . swap(M)^T mod 2.  Whole weight classes go       *)
(* through as a single matrix product.                                           *)
(* ============================================================================ *)

(* Check matrix with its X and Z halves exchanged: multiplying by its transpose
   is the symplectic form. Memoised per code. *)
codeSwapped[a_Association] := codeSwapped[a] = With[{mat = a["CheckMatrix"], n = a["Qubits"]},
    Join[mat[[All, n + 1 ;; 2 n]], mat[[All, 1 ;; n]], 2]
]

codeSyndromeVector[a_Association, v_List] := Mod[symplecticPart[v] . Transpose[codeSwapped[a]], 2]

QECCode::errsize = "The error must act on `1` qubits.";

codeSyndrome[a_Association, err_] := With[{v = QECPauliVector[err]},
    Which[
        v === $Failed, $Failed,
        Length[v] =!= 2 a["Qubits"] + 1, Message[QECCode::errsize, a["Qubits"]]; $Failed,
        True, codeSyndromeVector[a, v]
    ]
]

(* Syndromes of every weight-one error, keyed by the error. *)
codeSyndromeTable[a_Association] := codeSyndromeTable[a] = AssociationMap[
    codeSyndromeVector[a, QECPauliVector[#]] &,
    QECPauliString /@ weightOneVectors[a["Qubits"]]
]


(* ---- the lookup decoder ---- *)

(* The weight the code is guaranteed to correct, t = floor((d-1)/2). *)
codeCorrectableWeight[a_Association] := With[{d = codeDistance[a]},
    If[d === Infinity, 0, Floor[(d - 1) / 2]]
]

(* Minimum-weight lookup: walk the error weights upwards and keep the first error
   found for each syndrome, so every stored correction has least weight for it.

   How far up to walk is a real choice.  Correctness is only guaranteed to
   t = floor((d-1)/2), but a table built past t is still the maximum-likelihood
   guess under independent low-probability noise, and stopping at t would gut the
   distance-one teaching codes: the bit-flip code has d = 1, hence t = 0, yet
   correcting its single X errors is the whole point of it.  So the default reach
   is Max[t, 1] -- never less than the weight-one table the prototype always
   built, and more whenever the code earns it.  code["Decoder", w] asks for any
   other reach.

   This is a table, i.e. exact but exponential in the reach.  The matching and
   belief-propagation decoders of roadmap item 2 replace it for the codes where a
   table is not an option. *)
codeDecoderToWeight[a_Association, maxWeight_Integer] := codeDecoderToWeight[a, maxWeight] = Module[{n = a["Qubits"], table},
    table = <|codeSyndromeVector[a, pauliIdentity[n]] -> pauliIdentity[n]|>;
    Do[
        Do[
            With[{syn = codeSyndromeVector[a, v]},
                If[! KeyExistsQ[table, syn], table[syn] = v]
            ],
            {v, weightKVectors[n, w]}
        ],
        {w, 1, maxWeight}
    ];
    table
]

codeDecoderReach[a_Association] := Max[codeCorrectableWeight[a], 1]

codeDecoder[a_Association] := codeDecoderToWeight[a, codeDecoderReach[a]]

QECCode::badsyn = "The syndrome must be a binary vector of length `1`.";

codeDecode[a_Association, syn_List] := With[{m = codeStabilizerCount[a]},
    If[ Length[syn] =!= m || ! SubsetQ[{0, 1}, DeleteDuplicates[syn]],
        Message[QECCode::badsyn, m]; $Failed,
        Replace[Lookup[codeDecoder[a], Key[syn], Missing["UndecodableSyndrome", syn]], v : {__Integer} :> QECPauliString[v]]
    ]
]


(* ---- the cycle ---- *)

(* Error in, syndrome out, correction inferred, residual reported as an actual
   Pauli including its phase.
   The residual always commutes with every check, so it is either in the
   stabilizer group (up to the unobservable global phase: the correction worked)
   or in N(S) \ S (a silent logical error).  Those two are reported apart, in
   "Outcome": a decoder that returns the wrong coset is not the same failure as a
   syndrome it cannot decode at all, and the prototype's single Success boolean
   could not tell them apart. *)
codeCorrectionCycle[a_Association, err_] := Module[{v, syn, corr, residual, inStabilizer},
    v = QECPauliVector[err];
    If[v === $Failed || Length[v] =!= 2 a["Qubits"] + 1, Message[QECCode::errsize, a["Qubits"]]; Return[$Failed]];

    syn = codeSyndromeVector[a, v];
    corr = Lookup[codeDecoder[a], Key[syn], Missing["UndecodableSyndrome", syn]];

    If[ MissingQ[corr],
        Return[<|
            "Error" -> QECPauliString[v], "Syndrome" -> syn, "Correction" -> corr,
            "Residual" -> Missing["Undecodable"], "Outcome" -> "Undecodable", "Success" -> False
        |>]
    ];

    residual = QECPauliProduct[v, corr];
    inStabilizer = gf2MemberQ[codeBasis[a], symplecticPart[residual]];

    <|
        "Error" -> QECPauliString[v],
        "Syndrome" -> syn,
        "Correction" -> QECPauliString[corr],
        "Residual" -> QECPauliString[residual],
        "Outcome" -> If[inStabilizer, "Corrected", "LogicalError"],
        "Success" -> inStabilizer
    |>
]

(* The same cycle run on an actual encoded state through the engine: the error and
   the correction are applied as gates and the stabilizers are compared with their
   signs.  This checks restoration of that one fiducial state, so a residual that
   fixes it (a logical Z on the |0_L> completion) counts as success here and as a
   logical error in the symplectic cycle above.  The two are meant to disagree
   there; the difference is the point of running both. *)
codePhysicalCorrectionCycle[a_Association, err_] := Module[{gens, reference, damaged, syn, corr, final},
    gens = codeCompletedGenerators[a];
    If[gens === $Failed, Return[$Failed]];

    reference = Wolfram`QuantumFramework`PauliStabilizer[QECPauliString /@ gens];
    damaged = applyPauliVector[reference, QECPauliVector[err]];
    syn = (1 - (damaged["Expectation", #] & /@ (QECPauliString /@ codeVectors[a]))) / 2;
    corr = Lookup[codeDecoder[a], Key[syn], Missing["UndecodableSyndrome", syn]];

    If[ MissingQ[corr],
        Return[<|"Error" -> QECPauliString[QECPauliVector[err]], "Syndrome" -> syn, "Correction" -> corr, "Success" -> False|>]
    ];

    final = applyPauliVector[damaged, corr];

    <|
        "Error" -> QECPauliString[QECPauliVector[err]],
        "Syndrome" -> syn,
        "Correction" -> QECPauliString[corr],
        "Success" -> final["Stabilizers"] === reference["Stabilizers"]
    |>
]


(* ---- the quantum Hamming bound ---- *)

(* A code is perfect when the errors it corrects, translated by the codespace,
   exactly fill the Hilbert space: sum_(w<=t) 3^w Binomial[n,w] * 2^k == 2^n.
   The prototype hardcoded t = 1; with the distance available this is the general
   statement. *)
codePerfectQ[a_Association] := With[{n = a["Qubits"], k = codeLogicalQubits[a], t = codeCorrectableWeight[a]},
    Sum[3^w Binomial[n, w], {w, 0, t}] 2^k === 2^n
]
