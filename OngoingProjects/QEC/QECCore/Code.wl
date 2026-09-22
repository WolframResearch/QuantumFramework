(* ::Package:: *)

Package["Wolfram`QuantumFramework`QEC`"]

PackageExport[QECCode]

PackageScope[codeData]
PackageScope[codeCheckMatrix]
PackageScope[codeBasis]
PackageScope[codeQubits]
PackageScope[codeStabilizerCount]
PackageScope[codeLogicalQubits]
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

QECCode::usage = "QECCode[{gen1, gen2, ...}] represents a stabilizer code with the given Pauli string generators.\nQECCode[name] builds a named code, such as \"SteaneCode\" or \"5QubitCode\".\nQECCode[name, args] builds a member of a code family, such as QECCode[\"Repetition\", 5] or QECCode[\"CSS\", hx, hz].\nQECCode[ps] converts a Wolfram`QuantumFramework`PauliStabilizer into a code.\ncode[prop] gives a property; code[\"Properties\"] lists them.";


(* ---- accessors on the data association ---- *)

codeData[QECCode[a_Association]] := a
codeCheckMatrix[a_Association] := a["CheckMatrix"]
codeQubits[a_Association] := a["Qubits"]
codeStabilizerCount[a_Association] := Length[a["CheckMatrix"]]
codeLogicalQubits[a_Association] := a["Qubits"] - Length[a["CheckMatrix"]]

(* Generator rows in full Pauli-row form {x|z|e}. *)
codeVectors[a_Association] := MapThread[Append, {a["CheckMatrix"], a["Phases"]}]

(* Echelon basis of the check matrix, memoised: every stabilizer-membership
   question reduces against this instead of re-running elimination. *)
codeBasis[a_Association] := codeBasis[a] = gf2Basis[a["CheckMatrix"]]


(* ---- construction ---- *)

QECCode::gens = "Generators must be valid Pauli strings or Pauli rows, all acting on the same number of qubits.";
QECCode::noncomm = "Generators `1` and `2` do not commute.";
QECCode::dep = "Generators are not independent: the check matrix has rank `1` but `2` generators were given.";
QECCode::overcomplete = "`1` generators were given on `2` qubits; a stabilizer code needs at most n.";

fromVectors[vecs : {__List}] := Module[{n, m, mat, x, z, gram, bad, rank},
    If[ ! AllTrue[vecs, QECPauliQ] || Length[DeleteDuplicates[Length /@ vecs]] =!= 1,
        Message[QECCode::gens]; Return[$Failed]
    ];
    m = Length[vecs];
    n = pauliQubits[First[vecs]];
    If[m > n, Message[QECCode::overcomplete, m, n]; Return[$Failed]];
    mat = symplecticPart /@ vecs;
    x = mat[[All, 1 ;; n]];
    z = mat[[All, n + 1 ;; 2 n]];
    (* Symplectic Gram matrix: entry (i,j) is 0 exactly when generators i and j commute. *)
    gram = Mod[x . Transpose[z] + z . Transpose[x], 2];
    bad = FirstPosition[gram, 1, Missing["None"], {2}];
    If[ ! MissingQ[bad],
        Message[QECCode::noncomm, QECPauliString[vecs[[bad[[1]]]]], QECPauliString[vecs[[bad[[2]]]]]];
        Return[$Failed]
    ];
    rank = gf2Rank[mat];
    If[rank =!= m, Message[QECCode::dep, rank, m]; Return[$Failed]];
    QECCode[<|
        "CheckMatrix" -> mat,
        "Phases" -> (phasePart /@ vecs),
        "Qubits" -> n
    |>]
]

QECCode[gens : {__String}] := If[
    AllTrue[gens, QECPauliQ],
    fromVectors[QECPauliVector /@ gens],
    Message[QECCode::gens]; $Failed
]

QECCode[vecs : {{__Integer} ..}] := fromVectors[vecs]

QECCode[ps_Wolfram`QuantumFramework`PauliStabilizer] := QECCode[ps["Stabilizers"]]


(* ---- properties ---- *)

$codeDirectProperties = {
    "CheckMatrix", "Phases", "Qubits", "Generators", "GeneratorVectors", "Signs",
    "StabilizerCount", "LogicalQubits"
};

$codeDerivedProperties = {
    "Parameters", "Distance", "MinimumWeightLogical", "LogicalOperators", "LogicalX", "LogicalZ",
    "StandardForm", "CompletedGenerators", "SyndromeTable", "Decoder", "PerfectQ", "CSSQ", "SyndromeCircuit",
    "EncodingGates", "EncodingCircuit", "EncodingCircuitValidQ", "PauliStabilizer", "State"
};

$codeParametrizedProperties = {
    "Syndrome", "Decode", "Decoder", "LogicalErrorRate", "CorrectionCycle", "PhysicalCorrectionCycle",
    "LogicalPauliQ", "StabilizerMemberQ"
};

QECCode::noprop = "`1` is not a property of QECCode. Use code[\"Properties\"] for the list.";

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
QECCode[a_Association]["Generators"] := QECPauliString /@ codeVectors[a]
QECCode[a_Association]["Signs"] := Replace[a["Phases"], {0 -> 1, 2 -> -1, e_ :> I^e}, {1}]

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
QECCode[a_Association]["Decoder", QECNoiseModel[noise_Association]] := If[
    exactEnumerableQ[a],
    QECPauliString /@ codeCosetRepresentatives[a, noise],
    Message[QECLogicalErrorRate::mlneedsall, a["Qubits"], 4^a["Qubits"], $QECExactEnumerationLimit];
    $Failed
]

QECCode[a_Association]["LogicalErrorRate", noise_QECNoiseModel, opts___] :=
    QECLogicalErrorRate[QECCode[a], noise, opts]

QECCode[a_Association]["LogicalErrorRate", noise_QECNoiseModel, count_Integer, opts___] :=
    QECLogicalErrorRate[QECCode[a], noise, count, opts]
QECCode[a_Association]["Decoder", maxWeight_Integer] := codeDecoderToWeight[a, maxWeight]
QECCode[a_Association]["PerfectQ"] := codePerfectQ[a]
QECCode[a_Association]["CSSQ"] := codeCSSQ[a]

QECCode[a_Association]["SyndromeCircuit"] := QECSyndromeCircuit[QECCode[a]]
QECCode[a_Association]["SyndromeCircuit", rounds_] := QECSyndromeCircuit[QECCode[a], rounds]

QECCode[a_Association]["EncodingGates"] := codeEncodingGates[a]
QECCode[a_Association]["EncodingCircuit"] := Wolfram`QuantumFramework`QuantumCircuitOperator[codeEncodingGates[a]]
QECCode[a_Association]["EncodingCircuitValidQ"] := codeEncodingValidQ[a]

QECCode[a_Association]["PauliStabilizer"] := Wolfram`QuantumFramework`PauliStabilizer[QECPauliString /@ codeVectors[a]]
QECCode[a_Association]["State"] := QECCode[a]["PauliStabilizer"]["State"]

QECCode[a_Association]["Syndrome", err_] := codeSyndrome[a, err]
QECCode[a_Association]["Decode", syn_] := codeDecode[a, syn]
QECCode[a_Association]["CorrectionCycle", err_] := codeCorrectionCycle[a, err]
QECCode[a_Association]["PhysicalCorrectionCycle", err_] := codePhysicalCorrectionCycle[a, err]
QECCode[a_Association]["LogicalPauliQ", p_] := codeLogicalPauliQ[a, p]
QECCode[a_Association]["StabilizerMemberQ", p_] := codeStabilizerMemberQ[a, p]

QECCode[a_Association][prop_String] := (Message[QECCode::noprop, prop]; Missing["NotFound", prop])


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
            BoxForm`SummaryItem[{"Stabilizers: ", Row[QECPauliString /@ codeVectors[a], ", "]}],
            BoxForm`SummaryItem[{"Parameters: ", Dynamic[Row[{"[[", Row[obj["Parameters"], ","], "]]"}]]}],
            BoxForm`SummaryItem[{"Logical X: ", Dynamic[obj["LogicalX"]]}],
            BoxForm`SummaryItem[{"Logical Z: ", Dynamic[obj["LogicalZ"]]}],
            BoxForm`SummaryItem[{"CSS: ", Dynamic[obj["CSSQ"]]}]
        },
        form,
        "Interpretable" -> False
    ]
