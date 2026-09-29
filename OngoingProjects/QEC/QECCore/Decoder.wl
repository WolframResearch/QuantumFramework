(* ::Package:: *)

Package["Wolfram`QuantumFramework`QEC`"]

PackageExport[QECDecoder]

PackageScope[decoderData]


(* ============================================================================ *)
(* The decoder as an object.                                                    *)
(*                                                                              *)
(* Three decoders already lived inside the layer, each as a private association *)
(* reached through an option string: the minimum-weight lookup table            *)
(* (codeDecoderToWeight), maximum likelihood over cosets                        *)
(* (codeMaximumLikelihoodDecoder), and the detector model's lightest-fault-set  *)
(* rule (demDecoderTable).  None of them could be held, inspected, compared or  *)
(* replaced.  QECDecoder gives them one head, the three become its built-ins,   *)
(* and the place an external decoder plugs in is a constructor rather than a    *)
(* rewrite:                                                                     *)
(*                                                                              *)
(*     QECDecoder[code]                 minimum weight, at the code's reach     *)
(*     QECDecoder[code, noise]          maximum likelihood over cosets          *)
(*                                      (code capacity); at a circuit or        *)
(*                                      phenomenological level, the detector    *)
(*                                      model's decoder for that noise          *)
(*     QECDecoder[dem]                  the detector model's lightest-set rule  *)
(*     QECDecoder[code, f]              f[syndrome] -> correction, any function *)
(*                                                                              *)
(* THE SEAM.  The last form is how PyMatching, BP-OSD or a neural decoder come  *)
(* in: f takes the syndrome as a 0/1 list in the order code["Syndrome", e]      *)
(* reports it, and returns the correction as a Pauli (string or row) on the     *)
(* code's qubits, or a Missing to say it gives up.  The layer does not know     *)
(* what f does; it only asks it for corrections, so f can call out through      *)
(* ExternalEvaluate as easily as it can be a lookup.                            *)
(*                                                                              *)
(* NOTHING IS RE-IMPLEMENTED.  Every built-in delegates to the engine function  *)
(* that already computed it, with the same memoisation, so the object and the   *)
(* old option strings answer decision for decision identically -- including the *)
(* detector model's herald filter (rows a verification check rejects are not    *)
(* offered as explanations), which is what moved the measured exponent of the   *)
(* verified-cat extraction from 1.73 to 1.965 and is the easiest thing to lose  *)
(* in a refactor.  The tests compare tables, not rates, for that reason.        *)
(*                                                                              *)
(* A decoder is passed back to the rate as "Decoder" -> dec at code capacity,   *)
(* so a custom decoder is scored by exactly the machinery the built-ins are.    *)
(* ============================================================================ *)

QECDecoder::usage = "QECDecoder[code] is the minimum-weight lookup decoder of a stabilizer code.\nQECDecoder[code, noise] is the maximum-likelihood decoder for code-capacity noise, and the detector-model decoder for circuit-level or phenomenological noise.\nQECDecoder[dem] is the lightest-fault-set decoder of a detector error model.\nQECDecoder[code, f] wraps a decoding function f, which maps a syndrome (a 0/1 list) to a correction Pauli or Missing[]; this is how an external decoder is plugged in.\ndec[\"Decode\", syndrome] gives the correction. dec[prop] gives a property; dec[\"Properties\"] lists them.";

QECDecoder::badsyn = "The syndrome must be a 0/1 list of length `1`.";
QECDecoder::badcorr = "The decoding function returned `1` for syndrome `2`, which is neither a Pauli on `3` qubits nor a Missing.";
QECDecoder::noprop = "`1` is not a property of QECDecoder. Use dec[\"Properties\"] for the list.";
QECDecoder::othercode = "The decoder was built for a different code; it can only score the code it decodes.";
QECDecoder::detector = "A detector-model decoder reads detector patterns, not code-capacity syndromes. Its rate is QECDetectorModel[code, noise][\"LogicalErrorRate\"].";

Options[QECDecoder] = {"DecoderReach" -> Automatic, "Rounds" -> Automatic, "Extraction" -> "BareAncilla"};

decoderData[QECDecoder[d_Association]] := d


(* ---- construction ---- *)

QECDecoder[QECCode[a_Association], opts : OptionsPattern[]] := With[
    {reach = decoderReachFor[a, OptionValue[QECDecoder, {opts}, "DecoderReach"]]},
    QECDecoder[<|"Method" -> "MinimumWeight", "Code" -> a, "Reach" -> reach|>]
]

QECDecoder[QECCode[a_Association], QECNoiseModel[noise_Association], opts : OptionsPattern[]] := Which[
    noiseLevel[noise] =!= "CodeCapacity",
        With[{r = OptionValue[QECDecoder, {opts}, "Rounds"], e = OptionValue[QECDecoder, {opts}, "Extraction"]},
            QECDecoder[
                If[ r === Automatic,
                    QECDetectorModel[QECCode[a], QECNoiseModel[noise], "Extraction" -> e],
                    QECDetectorModel[QECCode[a], QECNoiseModel[noise], r, "Extraction" -> e]
                ],
                "DecoderReach" -> OptionValue[QECDecoder, {opts}, "DecoderReach"]
            ]
        ],
    ! exactEnumerableQ[a],
        Message[QECLogicalErrorRate::mlneedsall, a["Qubits"], 4^a["Qubits"], $QECExactEnumerationLimit];
        $Failed,
    True,
        QECDecoder[<|"Method" -> "MaximumLikelihood", "Code" -> a, "Noise" -> noise|>]
]

QECDecoder[code_QECCode, qc_Wolfram`QuantumFramework`QuantumChannel, opts : OptionsPattern[]] :=
    Replace[QECNoiseModel[qc], {noise_QECNoiseModel :> QECDecoder[code, noise, opts], _ -> $Failed}]

QECDecoder[QECDetectorModel[dem_Association], opts : OptionsPattern[]] := With[
    {reach = decoderReachFor[dem["Code"], OptionValue[QECDecoder, {opts}, "DecoderReach"]]},
    QECDecoder[<|"Method" -> "DetectorModel", "Code" -> dem["Code"], "Model" -> dem, "Reach" -> reach|>]
]

QECDecoder[QECCode[a_Association], f : Except[_QECNoiseModel | _Wolfram`QuantumFramework`QuantumChannel | _Rule | _RuleDelayed | _String | _Integer | _List]] :=
    QECDecoder[<|"Method" -> "Function", "Code" -> a, "Function" -> f|>]


(* A constructor whose input already failed has said why; it passes the failure on. *)
QECDecoder[$Failed, ___] := $Failed


(* ---- decoding ---- *)

decoderSyndromeLength[d_Association] := If[
    d["Method"] === "DetectorModel",
    d["Model"]["Detectors"],
    codeStabilizerCount[d["Code"]]
]

validSyndromeQ[d_Association, syn_] :=
    VectorQ[syn, MatchQ[0 | 1]] && Length[syn] === decoderSyndromeLength[d]

(* Syndrome in, correction out: a Pauli string on the code's qubits, or Missing when
   the decoder has no answer.  The detector-model decoder answers in its own terms, the
   logical observables to flip (a 0/1 list), since a circuit decoder never names a
   data-qubit correction -- it only needs the class. *)
decoderDecode[d_Association, syn_] := Which[
    ! validSyndromeQ[d, syn],
        Message[QECDecoder::badsyn, decoderSyndromeLength[d]]; $Failed,
    True,
        decodeWith[d["Method"], d, syn]
]

decodeWith["MinimumWeight", d_, syn_] := Replace[
    Lookup[codeDecoderToWeight[d["Code"], d["Reach"]], Key[syn], Missing["UndecodableSyndrome", syn]],
    v : {__Integer} :> pauliString[v]
]

decodeWith["MaximumLikelihood", d_, syn_] := Replace[
    Lookup[codeCosetRepresentatives[d["Code"], d["Noise"]], FromDigits[syn, 2], Missing["UndecodableSyndrome", syn]],
    v : {__Integer} :> pauliString[v]
]

decodeWith["DetectorModel", d_, syn_] := With[
    {obs = Lookup[demDecoderTable[d["Model"], d["Reach"]], FromDigits[syn, 2], Missing["UndecodableSyndrome", syn]]},
    If[MissingQ[obs], obs, IntegerDigits[obs, 2, d["Model"]["Observables"]]]
]

decodeWith["Function", d_, syn_] := With[{n = d["Code"]["Qubits"], out = d["Function"][syn]},
    Which[
        MissingQ[out], out,
        QECPauliQ[out] && pauliQubits[pauliVector[out]] === n, pauliString[out],
        True, Message[QECDecoder::badcorr, out, syn, n]; $Failed
    ]
]


(* ---- the table ---- *)

(* Every syndrome the decoder answers, with its answer.  For a function decoder this
   asks the function once per syndrome, all 2^m of them, which is what scoring it
   exactly needs anyway. *)
decoderTable[d_Association] := Switch[d["Method"],
    "MinimumWeight",
        pauliString /@ codeDecoderToWeight[d["Code"], d["Reach"]],
    "MaximumLikelihood",
        With[{m = codeStabilizerCount[d["Code"]]},
            KeyMap[IntegerDigits[#, 2, m] &, pauliString /@ codeCosetRepresentatives[d["Code"], d["Noise"]]]
        ],
    "DetectorModel",
        With[{nd = d["Model"]["Detectors"], no = d["Model"]["Observables"]},
            KeyMap[IntegerDigits[#, 2, nd] &, IntegerDigits[#, 2, no] & /@ demDecoderTable[d["Model"], d["Reach"]]]
        ],
    "Function",
        DeleteMissing @ AssociationMap[
            decodeWith["Function", d, #] &,
            Tuples[{0, 1}, codeStabilizerCount[d["Code"]]]
        ]
]

(* syndrome key -> class key, the form the exact and sampled rates consume.  A
   syndrome the decoder does not answer is absent, which the rate counts as a failure,
   exactly as it does for the built-in table. *)
decoderClassTable[d_Association] := Switch[d["Method"],
    "MaximumLikelihood", codeMaximumLikelihoodDecoder[d["Code"], d["Noise"]],
    "MinimumWeight", codeMinimumWeightDecoder[d["Code"], d["Reach"]],
    _, With[{table = decoderTable[d]},
        If[ table === <||>, <||>,
            With[{keys = codeLabelKeys[d["Code"], pauliVector /@ Values[table]]},
                AssociationThread[FromDigits[#, 2] & /@ Keys[table] -> Last[keys]]
            ]
        ]
    ]
]

(* The rate's "Decoder" option accepts the object.  A detector-model decoder has no
   code-capacity class table -- it decodes detectors, not syndromes -- and a decoder
   built for one code cannot score another. *)
decoderFor[a_Association, noise_Association, QECDecoder[d_Association], reach_] := Which[
    d["Method"] === "DetectorModel",
        Message[QECDecoder::detector]; $Failed,
    d["Code"] =!= a,
        Message[QECDecoder::othercode]; $Failed,
    True,
        decoderClassTable[d]
]


(* ---- properties ---- *)

$decoderProperties = {"Method", "Code", "Reach", "Noise", "DetectorModel", "Table", "Decode", "Properties"};

QECDecoder[_Association]["Properties"] := $decoderProperties
QECDecoder[d_Association]["Method"] := d["Method"]
QECDecoder[d_Association]["Code"] := QECCode[d["Code"]]
QECDecoder[d_Association]["Reach"] := Lookup[d, "Reach", Missing["NotApplicable"]]
QECDecoder[d_Association]["Noise"] := If[KeyExistsQ[d, "Noise"], QECNoiseModel[d["Noise"]], Missing["NotApplicable"]]
QECDecoder[d_Association]["DetectorModel"] := If[KeyExistsQ[d, "Model"], QECDetectorModel[d["Model"]], Missing["NotApplicable"]]
QECDecoder[d_Association]["Table"] := decoderTable[d]
QECDecoder[d_Association]["Decode", syn_] := decoderDecode[d, syn]
QECDecoder[d_Association][prop_String] := (Message[QECDecoder::noprop, prop]; Missing["NotFound", prop])


(* ---- formatting ---- *)

QECDecoder /: MakeBoxes[obj : QECDecoder[d_Association] /; KeyExistsQ[d, "Method"], form : (StandardForm | TraditionalForm)] :=
    BoxForm`ArrangeSummaryBox[
        QECDecoder,
        obj,
        None,
        {
            BoxForm`SummaryItem[{"Method: ", d["Method"]}],
            BoxForm`SummaryItem[{"Code: ", {d["Code"]["Qubits"], codeLogicalQubits[d["Code"]]}}]
        },
        {
            BoxForm`SummaryItem[{"Reach: ", Lookup[d, "Reach", "-"]}],
            BoxForm`SummaryItem[{"Syndrome bits: ", decoderSyndromeLength[d]}]
        },
        form,
        "Interpretable" -> False
    ]
