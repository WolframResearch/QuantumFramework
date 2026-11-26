(* ::Package:: *)

Package["Wolfram`QuantumFramework`QEC`"]

PackageExport[QECCode]

PackageScope[qecFailure]
PackageScope[qecMessageFailure]

PackageScope[codeData]
PackageScope[codeCheckMatrix]
PackageScope[codeBasis]
PackageScope[codeQubits]
PackageScope[codeStabilizerCount]
PackageScope[codeLogicalQubits]
PackageScope[codeCheckOrder]
PackageScope[codeVectors]
PackageScope[fromVectors]

(* Declared here, defined in the files that follow, so that every file resolves
   them to the same symbol regardless of load order. *)
PackageScope[codeStandardForm]        (* Structure.wl *)
PackageScope[codeLogicalOperators]
PackageScope[codeDistance]
PackageScope[codeMinimumWeightLogical]
PackageScope[codeLogicalPauliQ]
PackageScope[codeStabilizerMemberQ]
PackageScope[codeCompletedGenerators]
PackageScope[codeCSSQ]

PackageScope[codeSyndrome]            (* Syndrome.wl *)
PackageScope[codeSyndromeVector]
PackageScope[codeSyndromeTable]
PackageScope[codeDecoder]
PackageScope[codeDecode]
PackageScope[codeCorrectionCycle]
PackageScope[codePhysicalCorrectionCycle]
PackageScope[codePerfectQ]

PackageScope[codeEncodingGates]       (* Encoder.wl *)
PackageScope[codeEncodingValidQ]
PackageScope[applyPauliVector]
PackageScope[applyGates]

PackageScope[codeDetectorModel]       (* DetectorModel.wl *)
PackageScope[defaultRounds]

PackageScope[demRate]                 (* Memory.wl *)
PackageScope[demGroups]
PackageScope[demRowKeys]

PackageScope[codeStimCircuit]         (* Stim.wl *)
PackageScope[recordDetectors]
PackageScope[faultEffect]

PackageScope[noiseLevel]              (* Noise.wl *)
PackageScope[noiseMeasurementError]
PackageScope[noiseCircuitRates]
PackageScope[$QECCircuitLocations]

PackageScope[codeCircuitInstructions] (* Circuit.wl *)
PackageScope[generatorInstructions]
PackageScope[circuitInstructions]
PackageScope[circuitQubitCount]
PackageScope[framePropagate]
PackageScope[$circuitOneQubitOps]
PackageScope[frameZero]
PackageScope[instructionEngineGates]

PackageScope[$QECNamedCodes]          (* Families.wl *)

PackageScope[pauliProfile]            (* Noise.wl *)

PackageScope[codeMaximumLikelihoodDecoder]   (* ErrorRate.wl *)
PackageScope[codeCosetRepresentatives]
PackageScope[exactLogicalErrorRate]
PackageScope[sampledLogicalErrorRate]
PackageScope[logicalErrorRate]
PackageScope[exactEnumerableQ]


(* ============================================================================ *)
(* The code object.                                                             *)
(*                                                                              *)
(* A stabilizer code is stored the way the engine stores a tableau: a check      *)
(* matrix of symplectic rows plus a separate phase list, never as Pauli strings. *)
(* Strings are produced on demand by the "Generators" property.                  *)
(*                                                                              *)
(*   "CheckMatrix"  m x 2n integer matrix, row i = {x1..xn, z1..zn}             *)
(*   "Phases"       length-m list in Z4; Hermitian generators have 0 or 2       *)
(*   "Qubits"       n                                                           *)
(*                                                                              *)
(* Everything else -- distance, logical operators, standard form, decoder,       *)
(* encoder -- is derived and memoised on the data association, so a property is  *)
(* computed at most once per code no matter how often it is asked for.  The      *)
(* prototype recomputed the standard form inside LogicalOperators, inside        *)
(* CodeDistance, and again inside every summary-box redraw.                      *)
(* ============================================================================ *)

QECCode::usage = "QECCode[{gen1, gen2, ...}] represents a stabilizer code with the given Pauli string generators.\nQECCode[name] builds a named code, such as \"SteaneCode\" or \"5QubitCode\", or a code family at its default size.\nQECCode[name[params]] builds a member of a code family, such as QECCode[\"Surface\"[3]], QECCode[\"Repetition\"[5]] or QECCode[\"CSS\"[hx, hz]].\nQECCode[ps] converts a Wolfram`QuantumFramework`PauliStabilizer into a code.\ncode[prop] gives a property; code[\"Properties\"] lists them; $QECCodeNames lists the names.";


(* ---- equality ---- *)

(* Two codes are equal when they are the same code: the same number of qubits and the
   same stabilizer group, signs included, whatever generators were chosen for it --
   QuantumState and QuantumOperator compare what they represent the same way. *)
codeEqualQ[a1_Association, a2_Association] := a1 === a2 || (
    a1["Qubits"] === a2["Qubits"] &&
    Length[a1["CheckMatrix"]] === Length[a2["CheckMatrix"]] &&
    AllTrue[codeVectors[a2], codeStabilizerMemberQ[a1, #] &]
)

QECCode /: Equal[cs__QECCode] := AllTrue[Partition[First /@ {cs}, 2, 1], Apply[codeEqualQ]]
QECCode /: Unequal[cs__QECCode] := ! Equal[cs]


(* ---- failures, the framework's way ---- *)

(* A constructor that cannot build what it was asked for says why and returns a
   Failure carrying the same message, as QuantumState and QuantumChannel do:
   "InvalidName" for a name it does not know, "InvalidArguments" for a name it knows
   with parameters it cannot use. *)
SetAttributes[qecFailure, HoldRest]
qecFailure[tag_String, msg_MessageName, params___] := With[{ps = {params}},
    Message[msg, Sequence @@ ps];
    Failure[tag, <|"MessageTemplate" :> msg, "MessageParameters" -> ps|>]
]

(* The same for every other refusal: the message's own tag names the failure, so
   QECCode[{"XX", "XZ"}] gives Failure["noncomm", ...] and says why. *)
SetAttributes[qecMessageFailure, HoldAll]
qecMessageFailure[msg : MessageName[_, tag_String], params___] := qecFailure[tag, msg, params]


(* ---- accessors on the data association ---- *)

codeData[QECCode[a_Association]] := a
codeCheckMatrix[a_Association] := a["CheckMatrix"]
codeQubits[a_Association] := a["Qubits"]
codeStabilizerCount[a_Association] := Length[a["CheckMatrix"]]
codeLogicalQubits[a_Association] := a["Qubits"] - Length[a["CheckMatrix"]]

(* The order in which the bare-ancilla circuit touches each generator's qubits.  It
   is not part of the code -- any order gives the same stabilizer group -- but it
   decides which pairs of data errors one ancilla fault leaves behind (the hook error
   of Got26 sec. 12.1.1), so a family whose geometry fixes a good order carries it as
   "CheckOrder"; any other code gets its supports in qubit order. *)
codeCheckOrder[a_Association] := Lookup[a, "CheckOrder",
    Flatten[Position[BitOr @@ Partition[#, a["Qubits"]], 1]] & /@ a["CheckMatrix"]
]

(* Generator rows in full Pauli-row form {x|z|e}. *)
codeVectors[a_Association] := MapThread[Append, {a["CheckMatrix"], a["Phases"]}]

(* Echelon basis of the check matrix, memoised: every stabilizer-membership
   question reduces against this instead of re-running elimination. *)
codeBasis[a_Association] := qecMemoKey[codeBasis[a]] = gf2Basis[a["CheckMatrix"]]


(* ---- construction ---- *)

QECCode::gens = "Generators must be valid Pauli strings or Pauli rows, all acting on the same number of qubits.";
QECCode::noncomm = "Generators `1` and `2` do not commute.";
QECCode::dep = "Generators are not independent: the check matrix has rank `1` but `2` generators were given.";
QECCode::overcomplete = "`1` generators were given on `2` qubits; a stabilizer code needs at most n.";

(* The first pair of rows that anticommute: the symplectic Gram matrix (entry (i,j) is 0
   exactly when generators i and j commute), by the same rotation of the halves as
   symplecticProduct. *)
anticommutingPair[mat_] := FirstPosition[
    Mod[mat . RotateLeft[Transpose[mat], Length[First[mat]] / 2], 2], 1, Missing["None"], {2}
]

fromVectors[vecs : {__List}] /; ! (AllTrue[vecs, QECPauliQ] && SameQ @@ (Length /@ vecs)) :=
    qecMessageFailure[QECCode::gens]

fromVectors[vecs : {__List}] /; Length[vecs] > pauliQubits[First[vecs]] :=
    qecMessageFailure[QECCode::overcomplete, Length[vecs], pauliQubits[First[vecs]]]

fromVectors[vecs : {__List}] := With[{mat = symplecticPart /@ vecs},
    Replace[anticommutingPair[mat], {
        {i_, j_} :> qecMessageFailure[QECCode::noncomm, pauliString[vecs[[i]]], pauliString[vecs[[j]]]],
        _Missing :> With[{rank = gf2Rank[mat]},
            If[ rank === Length[vecs],
                QECCode[<|"CheckMatrix" -> mat, "Phases" -> phasePart /@ vecs, "Qubits" -> pauliQubits[First[vecs]]|>],
                qecMessageFailure[QECCode::dep, rank, Length[vecs]]
            ]
        ]
    }]
]

QECCode[gens : {__String ? QECPauliQ}] := fromVectors[pauliVector /@ gens]

QECCode[gens : {__String}] := qecMessageFailure[QECCode::gens]

QECCode[vecs : {{__Integer} ..}] := fromVectors[vecs]

QECCode[ps_Wolfram`QuantumFramework`PauliStabilizer] := QECCode[ps["Stabilizers"]]


(* ---- properties ---- *)

$codeDirectProperties = {
    "CheckMatrix", "Phases", "Qubits", "Generators", "GeneratorVectors", "Signs",
    "StabilizerCount", "LogicalQubits"
};

$codeDerivedProperties = {
    "Parameters", "Distance", "MinimumWeightLogical", "LogicalOperators", "LogicalX", "LogicalZ",
    "StandardForm", "CompletedGenerators", "SyndromeTable", "Decoder", "PerfectQ", "CSSQ", "SyndromeCircuit", "CheckOrder",
    "EncodingGates", "EncodingCircuit", "EncodingCircuitValidQ", "PauliStabilizer", "State",
    "Encoder", "Codewords", "Codespace", "SyndromeMeasurement",
    "WeightEnumerator", "SignedWeightEnumerator", "DualWeightEnumerator"
};

$codeParametrizedProperties = {
    "Syndrome", "Decode", "Decoder", "LogicalErrorRate", "CorrectionCycle", "PhysicalCorrectionCycle",
    "LogicalPauliQ", "StabilizerMemberQ",
    "Recovery", "LogicalChannel", "LogicalPauliProbabilities", "KnillLaflammeMatrix", "CorrectableQ",
    "Concatenate", "RemoveQubit", "Paste", "Stim", "CodespaceProbability"
};


(* "Decoder" belongs to two of the three lists, because code["Decoder"] (the lookup
   table) and code["Decoder", noise] (the coset map) are two call shapes of one name.
   Both dispatches stay; the listing names it once. *)
QECCode[_Association]["Properties"] :=
    DeleteDuplicates @ Join[$codeDirectProperties, $codeDerivedProperties, $codeParametrizedProperties]

QECCode[a_Association]["CheckMatrix"] := a["CheckMatrix"]
QECCode[a_Association]["Phases"] := a["Phases"]
QECCode[a_Association]["Qubits"] := a["Qubits"]
QECCode[a_Association]["StabilizerCount"] := codeStabilizerCount[a]
QECCode[a_Association]["LogicalQubits"] := codeLogicalQubits[a]
QECCode[a_Association]["GeneratorVectors"] := codeVectors[a]
QECCode[a_Association]["Generators"] := pauliString /@ codeVectors[a]
QECCode[a_Association]["Signs"] := I ^ a["Phases"]

QECCode[a_Association]["Parameters"] := {a["Qubits"], codeLogicalQubits[a], codeDistance[a]}
QECCode[a_Association]["Distance"] := codeDistance[a]
QECCode[a_Association]["MinimumWeightLogical"] := codeMinimumWeightLogical[a]
QECCode[a_Association]["LogicalOperators"] := codeLogicalOperators[a]
QECCode[a_Association]["LogicalX"] := codeLogicalOperators[a]["X"]
QECCode[a_Association]["LogicalZ"] := codeLogicalOperators[a]["Z"]
QECCode[a_Association]["StandardForm"] := codeStandardForm[a]
QECCode[a_Association]["CompletedGenerators"] := codeCompletedGenerators[a]

QECCode[a_Association]["SyndromeTable"] := codeSyndromeTable[a]
QECCode[a_Association]["Decoder"] := codeDecoder[a]

(* Given a noise model, the decoder stops being a minimum-weight table and becomes
   an inference: for each syndrome, the coset carrying the most probability.

   Weighing cosets means enumerating them, so this refuses rather than grinding on a
   code that cannot be enumerated; code["Decoder", w] is the cheap table for those. *)
QECCode[a_Association]["Decoder", QECNoiseModel[noise_Association]] /; exactEnumerableQ[a] :=
    pauliString /@ codeCosetRepresentatives[a, noise]

QECCode[a_Association]["Decoder", _QECNoiseModel] :=
    qecMessageFailure[QECLogicalErrorRate::mlneedsall, a["Qubits"], 4^a["Qubits"], $QECExactEnumerationLimit]

QECCode[a_Association]["LogicalErrorRate", noise : _QECNoiseModel | _Wolfram`QuantumFramework`QuantumChannel, rest___] :=
    QECLogicalErrorRate[QECCode[a], noise, rest]
QECCode[a_Association]["Decoder", maxWeight_Integer] := codeDecoderToWeight[a, maxWeight]
QECCode[a_Association]["PerfectQ"] := codePerfectQ[a]
QECCode[a_Association]["CSSQ"] := codeCSSQ[a]

QECCode[a_Association]["CheckOrder"] := codeCheckOrder[a]
QECCode[a_Association]["SyndromeCircuit"] := QECSyndromeCircuit[QECCode[a]]
QECCode[a_Association]["SyndromeCircuit", rounds_] := QECSyndromeCircuit[QECCode[a], rounds]

QECCode[a_Association]["EncodingGates"] := codeEncodingGates[a]
QECCode[a_Association]["EncodingCircuit"] := Wolfram`QuantumFramework`QuantumCircuitOperator[codeEncodingGates[a]]
QECCode[a_Association]["EncodingCircuitValidQ"] := codeEncodingValidQ[a]

(* The engine holds states, and m < n generators do not fix one: they fix a
   2^(n-m)-dimensional space.  So the code is handed over completed to n generators
   (codeCompletedGenerators), which pins one fiducial codeword, as the tech note
   describes.  Handing over the bare generators used to be accepted by the engine and
   built a tableau whose measurements were wrong; since November 2026 the engine
   refuses strings of rank below n, which is what exposed this. *)
QECCode[a_Association]["PauliStabilizer"] :=
    Replace[codeCompletedGenerators[a], gens_List :> Wolfram`QuantumFramework`PauliStabilizer[pauliString /@ gens]]
QECCode[a_Association]["State"] := QECCode[a]["PauliStabilizer"]["State"]

QECCode[a_Association]["Syndrome", err_] := codeSyndrome[a, err]

(* A code applied to a Pauli gives its syndrome: code[QECPauli["XIIIIII"]].  A bare
   string cannot be the argument, since code["..."] already asks for a property. *)
QECCode[a_Association][p_QECPauli] := codeSyndrome[a, pauliVector[p]]
QECCode[a_Association]["Decode", syn_] := codeDecode[a, syn]
QECCode[a_Association]["CorrectionCycle", err_] := codeCorrectionCycle[a, err]
QECCode[a_Association]["PhysicalCorrectionCycle", err_] := codePhysicalCorrectionCycle[a, err]
QECCode[a_Association]["LogicalPauliQ", p_] := codeLogicalPauliQ[a, p]
QECCode[a_Association]["StabilizerMemberQ", p_] := codeStabilizerMemberQ[a, p]



(* ---- formatting ---- *)

(* Collapsed rows are cheap data only; the expensive properties (distance, logical
   operators) sit behind Dynamic so opening the box is what pays for them. *)
QECCode /: MakeBoxes[obj : QECCode[a_Association] /; KeyExistsQ[a, "CheckMatrix"], form : (StandardForm | TraditionalForm)] :=
    BoxForm`ArrangeSummaryBox[
        QECCode,
        obj,
        ArrayPlot[a["CheckMatrix"], Mesh -> All, MeshStyle -> GrayLevel[0.8],
            ColorRules -> {0 -> White, 1 -> RGBColor[0.15, 0.5, 0.65]},
            ImageSize -> {Automatic, 34}, Frame -> False],
        {
            BoxForm`SummaryItem[{"Qubits: ", a["Qubits"]}],
            BoxForm`SummaryItem[{"Logical qubits: ", codeLogicalQubits[a]}],
            BoxForm`SummaryItem[{"Generators: ", codeStabilizerCount[a]}]
        },
        {
            BoxForm`SummaryItem[{"Stabilizers: ", Row[pauliString /@ codeVectors[a], ", "]}],
            BoxForm`SummaryItem[{"Parameters: ", Dynamic[Row[{"[[", Row[obj["Parameters"], ","], "]]"}]]}],
            BoxForm`SummaryItem[{"Logical X: ", Dynamic[obj["LogicalX"]]}],
            BoxForm`SummaryItem[{"Logical Z: ", Dynamic[obj["LogicalZ"]]}],
            BoxForm`SummaryItem[{"CSS: ", Dynamic[obj["CSSQ"]]}]
        },
        form,
        "Interpretable" -> False
    ]
