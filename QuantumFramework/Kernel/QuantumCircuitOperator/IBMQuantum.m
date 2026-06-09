Package["Wolfram`QuantumFramework`"]

PackageExport["IBMJobSubmit"]
PackageExport["IBMJob"]

PackageScope["iMethodIBMJob"]



(* ::Section:: *)
(* IBM Quantum job submission + a lossless result object.

   IBMJobSubmit[qco, backend] submits a circuit to an IBM QPU through the
   IBMQuantumPlatform service connection and returns an IBMJob handle.
   IBMJob[<|...|>] carries the raw service responses verbatim under "Raw" and
   exposes curated, WL-native accessors on top, so no field IBM returns is lost.

   An IBMJob may also be built by the blocking Method path (qiskitApply), which
   stores a precomputed "Measurement" and "Status" and, when a ServiceConnect
   connection is active, a "ServiceObject" for the REST-backed accessors. *)


$ibmService = "IBMQuantumPlatform"


(* ---- active-connection check (never creates or prompts) ---- *)

ibmActiveConnection[] := With[{so = ServiceFramework`GetDefaultServiceObject[$ibmService]},
    If[MatchQ[so, _ServiceObject], so, $Failed]
]


(* ---- option handling: CamelCase -> snake_case, recursively ----
   the snake_cased tree is applied onto qiskit's own primitive options object, so qiskit owns
   the option schema; this only fixes up the key spelling *)

ibmSnakeKey[s_] := ToLowerCase @ StringReplace[ToString[s], RegularExpression["(?<=[a-z0-9])([A-Z])"] -> "_$1"]

ibmSnakeKeys[a_Association] := Association @ KeyValueMap[ibmSnakeKey[#1] -> ibmSnakeKeys[#2] &, a]
ibmSnakeKeys[x_] := x


(* ---- submission ---- *)

Options[IBMJobSubmit] = {
    "Primitive" -> "sampler",
    "Shots" -> 4096,
    "Observable" -> "Z",
    "Wait" -> False,
    "PrimitiveOptions" -> <||>
}

IBMJobSubmit[qco_QuantumCircuitOperator, backend : _String | Automatic : Automatic, opts : OptionsPattern[]] := Enclose @ Module[{
    so, primitive, shots, waitQ, observable, userOpts, bk, obsTerms, sub, id, measured, job
},
    so = ibmActiveConnection[];
    If[ so === $Failed,
        Return @ Failure["IBMJobSubmit", <|
            "MessageTemplate" -> "No active IBM Quantum connection. Run ServiceConnect[\"`1`\"] first.",
            "MessageParameters" -> {$ibmService}
        |>]
    ];

    primitive = ToLowerCase @ ToString @ OptionValue["Primitive"];
    If[! MemberQ[{"sampler", "estimator"}, primitive],
        Return @ Failure["IBMJobSubmit", <|
            "MessageTemplate" -> "Primitive must be \"sampler\" or \"estimator\", not `1`.",
            "MessageParameters" -> {primitive}
        |>]
    ];
    shots = OptionValue["Shots"];
    waitQ = TrueQ @ OptionValue["Wait"];
    observable = OptionValue["Observable"];
    userOpts = Replace[OptionValue["PrimitiveOptions"], Except[_Association] :> <||>];

    bk = Replace[backend, Automatic :> First[so["Backends"], $Failed]];
    If[! StringQ[bk],
        Return @ Failure["IBMJobSubmit", <|"MessageTemplate" -> "Could not determine a default backend."|>]
    ];

    (* reject non-qubit systems up front with the same clear message QuantumQASM gives *)
    If[! qasmQubitsQ[qco], Return @ qasmNonQubitFailure[qco]];

    (* estimator observable -> Pauli triples qiskit turns into a layout-mapped SparsePauliOp *)
    obsTerms = If[primitive === "estimator", Confirm @ ibmObsTerms[observable], Null];

    (* qiskit builds, transpiles, lays out the observable, serializes, and submits via its own
       SamplerV2 / EstimatorV2; run() is not awaited, so a job handle returns immediately *)
    sub = qiskitPrimitiveSubmit[
        qco["Qiskit"], primitive, obsTerms, shots, ibmSnakeKeys[userOpts],
        "Provider" -> "IBMProvider", "Backend" -> bk
    ];
    If[! AssociationQ[sub], Return @ sub];   (* surface the qiskit Failure verbatim *)
    id = Lookup[sub, "JobID", $Failed];
    If[! StringQ[id],
        Return @ Failure["IBMJobSubmit", <|
            "MessageTemplate" -> "Job submission did not return a job id.",
            "Response" -> sub
        |>]
    ];
    measured = Lookup[sub, "MeasuredQubits", Automatic];

    job = IBMJob[<|"ID" -> id, "Backend" -> bk, "ServiceObject" -> so, "MeasuredQubits" -> measured, "Primitive" -> primitive|>];
    (* one initial status fetch so the handle / summary box has a status without re-querying *)
    If[waitQ, ibmWaitFor[job], job["Refresh"]]
]


(* Normalize the estimator observable into {pauliString, reCoeff, imCoeff} triples that qiskit
   turns into a SparsePauliOp. Accepts a single Pauli string, a list of Pauli strings, or a list
   of {pauliString, coefficient} pairs. The Pauli width must match the circuit's qubit count;
   qiskit's apply_layout then widens it to the device register. *)
ibmObsTerms[s_String] := {{s, 1., 0.}}
ibmObsTerms[l : {__String}] := {#, 1., 0.} & /@ l
ibmObsTerms[l : {{_String, _} ..}] := {#[[1]], N @ Re @ #[[2]], N @ Im @ #[[2]]} & /@ l
ibmObsTerms[other_] := Failure["IBMJobSubmit", <|
    "MessageTemplate" -> "Observable `1` must be a Pauli string (e.g. \"ZZZ\"), a list of Pauli strings, or a list of {pauliString, coefficient} pairs.",
    "MessageParameters" -> {other}
|>]


(* Estimator results come back as RuntimeEncoder ndarrays; let qiskit deserialize them rather than
   parsing the ndarray JSON in WL. Reconstructs the runtime service from the saved account (the
   same API key the WL connection uses) and reads the expectation values + standard errors. *)
ibmPrimitiveResult[id_String] := Enclose @ Block[{
    $jid = id,
    $token = Replace[SystemCredential["IBMQuantumPlatform_APIKEY"], Except[_String] -> Null]
},
    ConfirmBy[
        PythonEvaluate[Context[$jid], "
from qiskit_ibm_runtime import QiskitRuntimeService
from wolframclient.language import wl
import numpy as np
try:
    service = QiskitRuntimeService()
except Exception:
    token = <* $token *>
    if token is not None:
        QiskitRuntimeService.save_account(channel='ibm_quantum_platform', token=token, overwrite=True)
    service = QiskitRuntimeService()
pr = service.job(<* $jid *>).result()[0]
evs = np.atleast_1d(np.asarray(pr.data.evs, dtype=float)).tolist()
stds = np.atleast_1d(np.asarray(pr.data.stds, dtype=float)).tolist()
wl.Association(wl.Rule('EVs', evs), wl.Rule('Stds', stds))
", "ibm"],
        AssociationQ
    ]
]


(* constructor used by the blocking Method path: measurement + exact bit-keyed counts
   already in hand, job id known *)
iMethodIBMJob[qm_, counts_Association, jobId_String, backend_] := IBMJob[<|
    "ID" -> jobId,
    "Backend" -> Replace[backend, Except[_String] -> Missing["NotAvailable"]],
    "ServiceObject" -> ibmActiveConnection[],
    "Measurement" -> qm,
    "Counts" -> counts,
    "Status" -> "Completed"
|>]


(* ---- raw fetch (prefer stored snapshot, else query the service object) ---- *)

ibmSO[d_Association] := Replace[Lookup[d, "ServiceObject", $Failed], Except[_ServiceObject] -> $Failed]

ibmRequestName = <|"Details" -> "RawJobDetails", "Metrics" -> "RawJobMetrics", "Results" -> "RawJobResults"|>

ibmFetch[d_Association, sect_String] := With[{so = ibmSO[d]},
    If[so === $Failed, Missing["NoConnection"], so[ibmRequestName[sect], "JobID" -> d["ID"]]]
]

ibmRaw[d_Association, sect_String] := With[{stored = Lookup[Lookup[d, "Raw", <||>], sect, Missing[]]},
    If[MissingQ[stored], ibmFetch[d, sect], stored]
]

(* always an Association, so Lookup chains never hit a Missing/non-Association *)
ibmRawAssoc[d_Association, sect_String] := Replace[ibmRaw[d, sect], Except[_Association] -> <||>]

ibmStatus[det_Association] := Lookup[det, "status", Lookup[Lookup[det, "state", <||>], "status", "Unknown"]]
ibmStatus[_] := "Unknown"

ibmIsEstimator[d_Association] := Lookup[d, "Primitive", "sampler"] === "estimator"

$ibmTerminal = {"Completed", "Cancelled", "Failed"}


(* ---- refresh / wait ---- *)

IBMJob[d_Association]["Refresh"] := With[{so = ibmSO[d]},
    If[ so === $Failed, IBMJob[d],
        Module[{id = d["ID"], det, raw, extra = <||>},
            det = so["RawJobDetails", "JobID" -> id];
            raw = <|"Details" -> det|>;
            If[ ibmStatus[det] === "Completed",
                raw["Results"] = so["RawJobResults", "JobID" -> id];
                raw["Metrics"] = so["RawJobMetrics", "JobID" -> id];
                (* cache the qiskit-decoded expectation values once, so estimator accessors and
                   the summary box never re-fetch from IBM on every call / render *)
                If[ ibmIsEstimator[d],
                    With[{er = ibmPrimitiveResult[id]}, If[AssociationQ[er], extra["EstimatorResult"] = er]]
                ]
            ];
            IBMJob[<|d, "Raw" -> raw, extra|>]
        ]
    ]
]

ibmWaitFor[job0_] := NestWhile[
    (Pause[2]; #["Refresh"]) &,
    job0["Refresh"],
    ! MemberQ[$ibmTerminal, #["Status"]] &,
    1, 1800
]


(* ---- accessors ---- *)

IBMJob[d_Association]["ID"] := d["ID"]
IBMJob[d_Association]["Backend"] := Lookup[d, "Backend", Missing["NotAvailable"]]
IBMJob[d_Association]["ServiceObject"] := ibmSO[d]

IBMJob[d_Association]["Status"] := If[KeyExistsQ[d, "Status"], d["Status"], ibmStatus @ ibmRawAssoc[d, "Details"]]

IBMJob[d_Association]["UserID"]    := Lookup[ibmRawAssoc[d, "Details"], "user_id", Missing["NotAvailable"]]
IBMJob[d_Association]["ProgramID"] := Lookup[Lookup[ibmRawAssoc[d, "Details"], "program", <||>], "id", Missing["NotAvailable"]]

ibmSeconds[x_] := If[NumberQ[x], Quantity[x, "Seconds"], Missing["NotAvailable"]]

IBMJob[d_Association]["Cost"] := ibmSeconds @ Lookup[ibmRawAssoc[d, "Details"], "cost", Missing[]]
IBMJob[d_Association]["EstimatedRunTime"] := ibmSeconds @ Lookup[ibmRawAssoc[d, "Details"], "estimated_running_time_seconds", Missing[]]
IBMJob[d_Association]["QuantumSeconds"] := ibmSeconds @ Lookup[Lookup[ibmRawAssoc[d, "Metrics"], "usage", <||>], "quantum_seconds", Missing[]]

IBMJob[d_Association]["SubmittedCircuit"] := Enclose @ Confirm[
    Lookup[ibmRawAssoc[d, "Details"], "params", <||>][["pubs", 1, 1]], Missing["NotAvailable"]]

IBMJob[d_Association]["Options"] := Lookup[Lookup[ibmRawAssoc[d, "Details"], "params", <||>], "options", <||>]

ibmDate[s_String] := DateObject[s, TimeZone -> 0]
ibmDate[_] := Missing["NotAvailable"]

IBMJob[d_Association]["Timestamps"] := ibmDate /@ Lookup[ibmRawAssoc[d, "Metrics"], "timestamps", <||>]

IBMJob[d_Association]["Duration"] := With[{ts = Lookup[ibmRawAssoc[d, "Metrics"], "timestamps", <||>]},
    Module[{a = ibmDate @ Lookup[ts, "running", Null], b = ibmDate @ Lookup[ts, "finished", Null]},
        If[DateObjectQ[a] && DateObjectQ[b], b - a, Missing["NotAvailable"]]]]

(* classical-register samples -> {hex-list, num_bits} *)
ibmSamples[res_Association] := Enclose @ Module[{data, creg},
    data = ConfirmBy[res[["results", 1, "data"]], AssociationQ];
    creg = First @ Keys @ data;
    {data[creg]["samples"], data[creg]["num_bits"]}
]
ibmSamples[_] := $Failed

IBMJob[d_Association]["Samples"] := If[ibmIsEstimator[d], Missing["NotApplicable", "estimator"],
    With[{s = ibmSamples @ ibmRawAssoc[d, "Results"]},
        If[ListQ[s], First[s], Missing["NotAvailable"]]]]

(* exact integer counts keyed by the WL bit-list outcome. Method path: stored verbatim.
   Raw-samples path: built from the per-shot hex once the job is Completed. IBM is
   little-endian (c[0] is the LSB), so Reverse[IntegerDigits[...]] gives the bits in
   classical-bit order; qiskitReorderByQubit then sorts them into ascending-qubit order
   using the "MeasuredQubits" map (classical bit -> original qubit) captured from the
   submitted, transpiled circuit, so a permuted measurement or any non-identity transpile
   layout decodes correctly. *)
ibmCounts[d_Association] := If[
    ibmIsEstimator[d],
    Missing["NotApplicable", "estimator"],
    If[
    KeyExistsQ[d, "Counts"],
    d["Counts"],
    Module[{status, sn, samples, size, measuredQubits},
        status = If[KeyExistsQ[d, "Status"], d["Status"], ibmStatus @ ibmRawAssoc[d, "Details"]];
        If[status =!= "Completed", Return[Missing["JobNotComplete", status], Module]];
        sn = ibmSamples @ ibmRawAssoc[d, "Results"];
        If[! ListQ[sn], Return[Missing["NoSamples"], Module]];
        {samples, size} = sn;
        measuredQubits = Lookup[d, "MeasuredQubits", Automatic];
        Counts[qiskitReorderByQubit[Reverse[IntegerDigits[Interpreter["HexInteger"][#], 2, size]], measuredQubits] & /@ samples]
    ]
]]

ibmMeasurement[d_Association] := If[
    KeyExistsQ[d, "Measurement"],
    d["Measurement"],
    With[{c = ibmCounts[d]},
        If[ AssociationQ[c],
            With[{size = Length @ First @ Keys @ c},
                QuantumMeasurement @ If[size <= 10, Join[AssociationMap[0 &, Tuples[{0, 1}, size]], c], c]],
            c
        ]
    ]
]

(* classical bit -> original-qubit map used to decode samples into ascending-qubit order *)
IBMJob[d_Association]["MeasuredQubits"] := Lookup[d, "MeasuredQubits", Missing["NotAvailable"]]

IBMJob[d_Association]["NumBits"] := With[{c = ibmCounts[d]},
    If[AssociationQ[c], Length @ First @ Keys @ c, Missing["NotAvailable"]]]
IBMJob[d_Association]["Qubits"] := IBMJob[d]["NumBits"]

IBMJob[d_Association]["Shots"] := With[{c = ibmCounts[d]},
    If[AssociationQ[c], Total[c], Missing["NotAvailable"]]]

IBMJob[d_Association]["Measurement"]   := ibmMeasurement[d]
IBMJob[d_Association]["Counts"]        := ibmCounts[d]
IBMJob[d_Association]["Probabilities"] := With[{m = ibmMeasurement[d]}, If[MatchQ[m, _QuantumMeasurement], m["Probabilities"], m]]

(* ---- estimator results: the qiskit-decoded expectation values + standard errors ---- *)

IBMJob[d_Association]["Primitive"] := Lookup[d, "Primitive", "sampler"]

(* prefer the value cached at Refresh; fall back to a live qiskit fetch for a handle that has
   not been refreshed since completion *)
ibmEstimatorResult[d_Association] := If[
    KeyExistsQ[d, "EstimatorResult"],
    d["EstimatorResult"],
    With[{status = IBMJob[d]["Status"]},
        If[status =!= "Completed", Missing["JobNotComplete", status], ibmPrimitiveResult[d["ID"]]]
    ]
]

IBMJob[d_Association]["ExpectationValues"] := If[! ibmIsEstimator[d], Missing["NotApplicable", "sampler"],
    With[{r = ibmEstimatorResult[d]}, If[AssociationQ[r], r["EVs"], r]]]

IBMJob[d_Association]["StandardErrors"] := If[! ibmIsEstimator[d], Missing["NotApplicable", "sampler"],
    With[{r = ibmEstimatorResult[d]}, If[AssociationQ[r], r["Stds"], r]]]

IBMJob[d_Association]["ExpectationValue"] := With[{evs = IBMJob[d]["ExpectationValues"]},
    If[ListQ[evs] && Length[evs] === 1, First[evs], evs]]

(* the server-side-transpiled circuit IBM ran; presigned URL -> QPY blob. Two
   constraints shape the failure handling:
   - IBM serves the blob from Cloud Object Storage's direct (VPC-private) endpoint
     (s3.direct.<region>...), reachable only from inside IBM Cloud; the AWS-SigV4
     presigned URL signs the host header, so it cannot be retargeted to the public
     endpoint. From elsewhere the download times out (URLRead returns a Failure).
   - The transpiled-circuit object exists only when the Runtime transpiles
     server-side. Jobs that submit already-ISA circuits to the V2 primitives carry
     out no server-side transpilation, so no object is written and the (still
     presigned) URL returns HTTP 404. *)
IBMJob[d_Association]["ExecutedCircuit"] := With[{so = ibmSO[d]},
    If[ so === $Failed, Missing["NoConnection"],
        Enclose @ Module[{url, resp, status, bytes},
            url = ConfirmBy[so["RawJobTranspiledCircuits", "JobID" -> d["ID"]]["url"], StringQ];
            resp = URLRead[url, TimeConstraint -> 30];
            If[ ! MatchQ[resp, _HTTPResponse],
                Confirm @ Failure["IBMExecutedCircuit", <|
                    "MessageTemplate" -> "Could not reach `1` to download the transpiled circuit. A Cloud Object Storage direct endpoint (host starting s3.direct) is reachable only from inside IBM Cloud.",
                    "MessageParameters" -> {Replace[URLParse[url, "Domain"], Except[_String] -> url]}
                |>]
            ];
            status = resp["StatusCode"];
            If[ status === 404,
                Confirm @ Failure["IBMExecutedCircuit", <|
                    "MessageTemplate" -> "IBM stored no transpiled circuit for job `1`. The object exists only when the Runtime transpiles server-side; already-ISA circuits submitted to the V2 primitives store none.",
                    "MessageParameters" -> {d["ID"]}
                |>]
            ];
            If[ status =!= 200,
                Confirm @ Failure["IBMExecutedCircuit", <|
                    "MessageTemplate" -> "Downloading the transpiled circuit for job `1` returned HTTP `2`.",
                    "MessageParameters" -> {d["ID"], status}
                |>]
            ];
            bytes = ConfirmBy[resp["BodyByteArray"], ByteArrayQ];
            Confirm @ ibmQPYToCircuit[bytes]
        ]
    ]
]

(* handles both zlib-compressed and raw QPY blobs *)
ibmQPYToCircuit[bytes_ByteArray] := Enclose @ Block[{$bytes = bytes}, QuantumCircuitOperator @ First @ ConfirmMatch[
    PythonEvaluate[Context[$bytes], "
import io, zlib, pickle
from qiskit import qpy
from wolframclient.language import wl
data = bytes(<* $bytes *>)
try:
    buf = io.BytesIO(zlib.decompress(data))
except Exception:
    buf = io.BytesIO(data)
[wl.Wolfram.QuantumFramework.QiskitCircuit(pickle.dumps(qc)) for qc in qpy.load(buf)]
", "default"],
    {__QiskitCircuit}
]]

IBMJob[d_Association]["ExecutionSpans"] := Lookup[
    Lookup[Lookup[ibmRawAssoc[d, "Results"], "metadata", <||>], "execution", <||>], "execution_spans", {}]

IBMJob[d_Association]["Cancel"] := With[{so = ibmSO[d]},
    If[so === $Failed, Missing["NoConnection"], so["RawJobCancel", "JobID" -> d["ID"]]]]

IBMJob[d_Association]["Raw"] := Lookup[d, "Raw", <||>]
IBMJob[d_Association]["Raw", sect_String] := ibmRaw[d, sect]

$ibmProperties = {
    "ID", "Backend", "ServiceObject", "Status", "Primitive", "UserID", "ProgramID",
    "Cost", "EstimatedRunTime", "SubmittedCircuit", "Options", "QuantumSeconds",
    "Timestamps", "Duration", "Samples", "MeasuredQubits", "NumBits", "Qubits", "Shots",
    "Measurement", "Counts", "Probabilities", "ExpectationValue", "ExpectationValues",
    "StandardErrors", "ExecutedCircuit", "ExecutionSpans",
    "Refresh", "Cancel", "Raw", "Properties"
}
IBMJob[_Association]["Properties"] := $ibmProperties

(* default application -> the primary result for the job's primitive *)
IBMJob[d_Association][] := If[ibmIsEstimator[d], IBMJob[d]["ExpectationValue"], ibmMeasurement[d]]

(* catch-all for unknown string properties *)
IBMJob[d_Association][prop_String] /; ! MemberQ[$ibmProperties, prop] :=
    Failure["IBMJob", <|
        "MessageTemplate" -> "Unknown property `1`. Use one of: `2`.",
        "MessageParameters" -> {prop, $ibmProperties}
    |>]


(* ---- formatting: status + the quantum-service summary-box icon ---- *)

$quantumServiceIcon := $quantumServiceIcon = Replace[
    Enclose @ With[{path = PacletObject[$ibmService]["AssetLocation", "icon"]},
        ImageTrim[
            ConfirmBy[If[StringQ[path] && FileExistsQ[path], Import[path], $Failed], ImageQ],
            ImageValuePositions[#, "Max"] &
        ]
    ],
    Except[_Image] -> None
]

$ibmStatusColor = <|
    "Completed" -> Darker[Green], "Running" -> RGBColor[0.05, 0.4, 0.85],
    "Queued" -> Orange, "Initializing" -> Orange, "Validating" -> Orange,
    "Cancelled" -> Gray, "Failed" -> Red
|>

ibmDisplayStatus[d_Association] := If[KeyExistsQ[d, "Status"], d["Status"],
    ibmStatus @ Lookup[Lookup[d, "Raw", <||>], "Details", <||>]]

IBMJob /: MakeBoxes[job : IBMJob[d_Association], fmt_] := Block[{
    st = ibmDisplayStatus[d],
    icon = Replace[$quantumServiceIcon, Except[_Image] -> QuantumCircuitOperator[{"H"}]["Icon"]],
    extra
},
    extra = Which[
        st =!= "Completed",
        {{BoxForm`SummaryItem[{"", Style["Results not available yet.", Gray]}]}},
        ibmIsEstimator[d],
        {
            {BoxForm`SummaryItem[{"Expectation: ", job["ExpectationValue"]}]},
            {BoxForm`SummaryItem[{"Quantum time: ", job["QuantumSeconds"]}]}
        },
        True,
        {
            {BoxForm`SummaryItem[{"Qubits: ", job["NumBits"]}]},
            {BoxForm`SummaryItem[{"Shots: ", job["Shots"]}]},
            {BoxForm`SummaryItem[{"Quantum time: ", job["QuantumSeconds"]}]}
        }
    ];
    BoxForm`ArrangeSummaryBox[
        "IBMJob",
        job,
        Show[icon, ImageSize -> 12],
        {
            {BoxForm`SummaryItem[{"Status: ",
                Row[{Style["\[FilledCircle] ", Lookup[$ibmStatusColor, st, Gray]], st}]}]},
            {BoxForm`SummaryItem[{"Backend: ", Lookup[d, "Backend", Missing[]]}]},
            {BoxForm`SummaryItem[{"Job ID: ", d["ID"]}]}
        },
        extra,
        fmt,
        "Interpretable" -> Automatic
    ]
]
