Needs["ServiceFramework`" -> "SF`"]

BeginPackage["IBMQuantumPlatform`"]

Begin["`Private`"]

(* device error model + ErrorMap graphic (iIBM* helpers, pure built-in graphics).
   Loaded here so the ProcessedRequests below can reference them. *)
Get[FileNameJoin[{DirectoryName[$InputFileName], "ErrorMap.m"}]]


IBMToken[key_] := Enclose[
    Block[{
        data, token
    },
        data = LocalSymbol["ibm_token"];
        If[ ! AssociationQ[data] || FromUnixTime[Lookup[data, "expiration", 0]] < Now,
            SystemCredential["IBMQuantumPlatform_APIKEY"] = key;
            data = LocalSymbol["ibm_token"] =
                Association @ Confirm @ URLExecute @ HTTPRequest[
                    "https://iam.cloud.ibm.com/identity/token",
                    <|
                        Method -> "POST",
                        "Content-Type" -> "application/x-www-form-urlencoded",
                        "Body" -> <|
                            "grant_type" -> "urn:ibm:params:oauth:grant-type:apikey", 
                            "apikey" -> key
                        |>
                    |>
                ]
        ];
        token = ConfirmBy[Lookup[data, "access_token"], StringQ];
        If[ TrueQ[Lookup[data, "token_type"] === "Bearer"],
            "Bearer " <> token,
            Confirm @ Failure["IBMTokenFailure", <|
                "MessageTemplate" -> "Unknown IBM Quantum Platform token type",
                "MessageParameters" -> {token}
            |>]
        ]
    ],
    Throw
]

IBMQuantumPlatform::missingCRN = "CRN is required for IBM Quantum Platform service connection. Please set crn to LocalSymbol[\"ibm_crn\"] before using the service connection.";

IBMCrn[] := With[{crn = LocalSymbol["ibm_crn"]},
    If[ MissingQ[crn],
        Message[IBMQuantumPlatform::missingCRN];
        Throw[$Failed]
    ];
    crn
]

params = <|
    "ServiceFrameworkVersion" -> "0.1.0",
    "ServiceName" -> "IBMQuantumPlatform",
    "Information" -> "IBM Quantum Platform service connection for Wolfram Language",
    "Icon" -> Show[ImageTrim[Import[PacletObject["IBMQuantumPlatform"]["AssetLocation", "icon"]], ImageValuePositions[#, "Max"] &], ImageSize -> 12],
    "AuthenticationMethod" -> "APIKey",
    "AuthenticationFunction" -> Function[key, <|"Headers" -> {"Authorization" -> IBMToken[key], "Service-CRN" -> IBMCrn[]}|>],
    "ClientDialog" -> Function[
        SF`MultipleKeyDialog[
            "IBM Quantum Platform",
            {"APIKey" -> {"apikey", FieldMasked -> True}},
            ##,
            "https://quantum.ibm.com/",
            "https://www.ibm.com/us-en/privacy"
        ]
    ],
    "RawRequests" -> <||>
|>


(* Raw Requests *)

$commonHeaders = {"Content-Type" -> "application/json"};

$serviceEndpoint = "https://quantum.cloud.ibm.com/api"

$version = "v1"

path[elements___] := URLBuild[{$serviceEndpoint, $version, elements}]

params["RawRequests", "RawBackends"] = {
    "URL"				-> path["backends"],
    "HTTPSMethod"		-> "GET",
    "Headers"			-> $commonHeaders
}

AppendTo[
    params["RawRequests"]
    ,
    Map[endpoint |->
        "RawBackend" <> Capitalize[endpoint] -> <|
            "URL"				-> Function[URLBuild[{path["backends"], #BackendName, endpoint}]],
            "HTTPSMethod"		-> "GET",
            "Headers"			-> $commonHeaders,
            "PathParameters"    -> {"BackendName"},
            "RequiredParameters" -> {"BackendName"}
        |>
        , 
        {"configuration", "defaults", "properties", "status"}
    ]
]

params["RawRequests", "RawInstances"] = {
    "URL"				-> path["instances"],
    "HTTPSMethod"		-> "GET",
    "Headers"			-> $commonHeaders
}

params["RawRequests", "RawInstances"] = {
    "URL"				-> path["instances"],
    "HTTPSMethod"		-> "GET",
    "Headers"			-> $commonHeaders
}

params["RawRequests", "RawInstancesUsage"] = {
    "URL"				-> path["instances", "usage"],
    "HTTPSMethod"		-> "GET",
    "Headers"			-> $commonHeaders
}

params["RawRequests", "RawInstancesConfiguration"] = {
    "URL"				-> path["instances", "configuration"],
    "HTTPSMethod"		-> "GET",
    "Headers"			-> $commonHeaders
}

params["RawRequests", "RawJobs"] = {
    "URL"				-> path["jobs"],
    "HTTPSMethod"		-> "GET",
    "Headers"			-> $commonHeaders
}

params["RawRequests", "RawWorkloads"] = {
    "URL"                -> path["workloads"],
    "HTTPSMethod"        -> "GET",
    "Headers"            -> $commonHeaders
}

params["RawRequests", "RawAccountDetails"] = {
    "URL"                -> Function[URLBuild[{path["accounts"], #AccountID}]],
    "HTTPSMethod"        -> "GET",
    "Headers"            -> $commonHeaders,
    "PathParameters"     -> {"AccountID"},
    "RequiredParameters" -> {"AccountID"}
}

params["RawRequests", "RawJobRun"] = {
    "URL"               -> path["jobs"],
    "HTTPSMethod"       -> "POST",
    "Headers"           -> $commonHeaders,
    "BodyParameters"    -> {"program_id", "backend", "runtime", "tags", "log_level", "cost", "session_id", "params", "private"},
    "RequiredParameters" -> {"program_id", "backend", "params"},
    "HiddenParameters"  -> {},
    "HTTPResponseProcessing" -> SF`ImportResponse
}

params["RawRequests", "RawJobDetails"] = {
    "URL"                -> Function[URLBuild[{path["jobs"], #JobID}]],
    "HTTPSMethod"        -> "GET",
    "Headers"            -> $commonHeaders,
    "PathParameters"     -> {"JobID"},
    "RequiredParameters" -> {"JobID"},
    "HiddenParameters"  -> {},
    "HTTPResponseProcessing" -> SF`ImportResponse
}

params["RawRequests", "RawJobResults"] = {
    "URL"                -> Function[URLBuild[{path["jobs"], #JobID, "results"}]],
    "HTTPSMethod"        -> "GET",
    "Headers"            -> $commonHeaders,
    "PathParameters"     -> {"JobID"},
    "RequiredParameters" -> {"JobID"},
    "HiddenParameters"  -> {},
    "HTTPResponseProcessing" -> Function[Enclose[ImportString[SF`ImportResponse[Confirm[#]], "RawJSON"]]]
}

params["RawRequests", "RawJobInterimResults"] = {
    "URL"                -> Function[URLBuild[{path["jobs"], #JobID, "interim_results"}]],
    "HTTPSMethod"        -> "GET",
    "Headers"            -> $commonHeaders,
    "PathParameters"     -> {"JobID"},
    "RequiredParameters" -> {"JobID"}
}

params["RawRequests", "RawJobTranspiledCircuits"] = {
    "URL"                -> Function[URLBuild[{path["jobs"], #JobID, "transpiled_circuits"}]],
    "HTTPSMethod"        -> "GET",
    "Headers"            -> $commonHeaders,
    "PathParameters"     -> {"JobID"},
    "RequiredParameters" -> {"JobID"}
}

params["RawRequests", "RawJobLogs"] = {
    "URL"                -> Function[URLBuild[{path["jobs"], #JobID, "logs"}]],
    "HTTPSMethod"        -> "GET",
    "Headers"            -> $commonHeaders,
    "PathParameters"     -> {"JobID"},
    "RequiredParameters" -> {"JobID"}
}

params["RawRequests", "RawJobCancel"] = {
    "URL"                -> Function[URLBuild[{path["jobs"], #JobID, "cancel"}]],
    "HTTPSMethod"        -> "POST",
    "Headers"            -> $commonHeaders,
    "PathParameters"     -> {"JobID"},
    "RequiredParameters" -> {"JobID"}
}

params["RawRequests", "RawJobMetrics"] = {
    "URL"                -> Function[URLBuild[{path["jobs"], #JobID, "metrics"}]],
    "HTTPSMethod"        -> "GET",
    "Headers"            -> $commonHeaders,
    "PathParameters"     -> {"JobID"},
    "RequiredParameters" -> {"JobID"}
}

params["RawRequests", "RawJobTags"] = {
    "URL"                -> Function[URLBuild[{path["jobs"], #JobID, "tags"}]],
    "HTTPSMethod"        -> "PUT",
    "Headers"            -> $commonHeaders,
    "PathParameters"     -> {"JobID"},
    "RequiredParameters" -> {"JobID"}
}

params["RawRequests", "RawBackends"] = {
    "URL"                -> path["backends"],
    "HTTPSMethod"        -> "GET",
    "Headers"            -> $commonHeaders
}

params["RawRequests", "RawBackendProperties"] = {
    "URL"                -> Function[URLBuild[{path["backends"], #BackendID, "properties"}]],
    "HTTPSMethod"        -> "GET",
    "Headers"            -> $commonHeaders,
    "PathParameters"     -> {"BackendID"},
    "RequiredParameters" -> {"BackendID"},
    "HTTPResponseProcessing" -> Function[Enclose[SF`ImportResponse[Confirm[#]]]]
}

params["RawRequests", "RawBackendConfiguration"] = {
    "URL"                -> Function[URLBuild[{path["backends"], #BackendID, "configuration"}]],
    "HTTPSMethod"        -> "GET",
    "Headers"            -> $commonHeaders,
    "PathParameters"     -> {"BackendID"},
    "RequiredParameters" -> {"BackendID"},
    "HTTPResponseProcessing" -> Function[Enclose[SF`ImportResponse[Confirm[#]]]]
}

params["RawRequests", "RawBackendStatus"] = {
    "URL"                -> Function[URLBuild[{path["backends"], #BackendID, "status"}]],
    "HTTPSMethod"        -> "GET",
    "Headers"            -> $commonHeaders,
    "PathParameters"     -> {"BackendID"},
    "RequiredParameters" -> {"BackendID"}
}

params["RawRequests", "RawBackendDefaults"] = {
    "URL"                -> Function[URLBuild[{path["backends"], #BackendID, "defaults"}]],
    "HTTPSMethod"        -> "GET",
    "Headers"            -> $commonHeaders,
    "PathParameters"     -> {"BackendID"},
    "RequiredParameters" -> {"BackendID"}
}

params["RawRequests", "RawSessions"] = {
    "URL"                -> path["sessions"],
    "HTTPSMethod"        -> "GET",
    "Headers"            -> $commonHeaders
}

params["RawRequests", "RawSessionDetails"] = {
    "URL"                -> Function[URLBuild[{path["sessions"], #SessionID}]],
    "HTTPSMethod"        -> "GET",
    "Headers"            -> $commonHeaders,
    "PathParameters"     -> {"SessionID"},
    "RequiredParameters" -> {"SessionID"}
}

params["RawRequests", "RawSessionClose"] = {
    "URL"                -> Function[URLBuild[{path["sessions"], #SessionID, "close"}]],
    "HTTPSMethod"        -> "DELETE",
    "Headers"            -> $commonHeaders,
    "PathParameters"     -> {"SessionID"},
    "RequiredParameters" -> {"SessionID"}
}

params["RawRequests", "RawInstance"] = {
    "URL"                -> path["instance"],
    "HTTPSMethod"        -> "GET",
    "Headers"            -> $commonHeaders
}

params["RawRequests", "RawInstanceConfiguration"] = {
    "URL"                -> path["instances", "configuration"],
    "HTTPSMethod"        -> "GET",
    "Headers"            -> $commonHeaders
}

params["RawRequests", "RawInstanceUsage"] = {
    "URL"                -> path["instances", "usage"],
    "HTTPSMethod"        -> "GET",
    "Headers"            -> $commonHeaders
}

params["RawRequests", "RawTags"] = {
    "URL"                -> path["tags"],
    "HTTPSMethod"        -> "GET",
    "Headers"            -> $commonHeaders
}



params["ParameterMap"] = <|
    "ProgramID" -> "program_id",
    "Backend" -> "backend",
    "Runtime" -> "runtime",
    "Tags" -> "tags",
    "LogLevel" -> "log_level",
    "Cost" -> "cost",
    "SessionID" -> "session_id",
    "Parameters" -> "params",
    "Private" -> "private"
|>

(* Processed Requests *)

(* These use immediate assignment (=), not SetDelayed (:=): := on an
   Association part stores a RuleDelayed entry, and ServiceFramework's
   addDefaultParameters cannot Part-assign default keys into a delayed-held
   entry, which surfaces as Set::pspec1 and Function::slota at load time. *)
params["ProcessedRequests"] = <||>;

params["ProcessedRequests", "Backends"] = <|
    "ExecuteFunction" -> SF`RequestExecute["RawBackends"],
    "SubmitFunction" -> SF`RequestExecute["RawBackends"],
    "Parameters" -> {},
    "RequiredParameters" -> {},
    "HiddenParameters" -> {},
    "PreprocessingFunction" -> Identity,
    "ExecuteResultProcessing" -> Function[Enclose[Lookup[ConfirmBy[#, AssociationQ], "devices"]]]
|>
params["ProcessedRequests", "JobRun"] = <|
		"ExecuteFunction" -> SF`RequestExecute["RawJobRun"],
		"SubmitFunction" -> SF`RequestSubmit["RawJobRun"],
		"Parameters" -> {
			"QASM" -> None, "Observable" -> "Z",
            "Backend" -> Automatic, "ProgramID" -> Automatic
		},
        "RequiredParameters" -> {"QASM"},
        "HiddenParameters" -> {},
        "PreprocessingFunction" -> Function[Enclose @ With[{program = Lookup[#, "ProgramID"]},
            <|
                "program_id" -> Replace[program, Automatic -> "sampler"],
                "backend" -> Replace[Lookup[#, "Backend"], _Missing | Automatic :> First[SF`GetDefaultServiceObject["IBMQuantumPlatform"]["Backends"]]],
                "params" -> <|
                    "pubs" -> {{ConfirmBy[Lookup[#, "QASM"], StringQ], Switch[program, "estimator", Lookup[#, "Observable"], "sampler" | _, Nothing]}},
                    "options" -> <|"dynamical_decoupling" -> <|"enable" -> True|>|>,
                    "version" -> 2
                |>
            |>]
        ],
        "ExecuteResultProcessing" -> Identity
	|>

params["ProcessedRequests", "JobResults"] = <|
    "ExecuteFunction" -> SF`RequestExecute["RawJobResults"],
    "SubmitFunction" -> SF`RequestExecute["RawJobResults"],
    "Parameters" -> {"JobID" -> Automatic},
    "RequiredParameters" -> {"JobID"},
    "HiddenParameters" -> {},
    "PreprocessingFunction" -> Function[Enclose[<|
        "JobID" -> ConfirmBy[Lookup[#, "JobID"], StringQ]
    |>]],
    "ExecuteResultProcessing" -> Function[Enclose[Lookup[#, "data"] & /@ ConfirmMatch[Lookup[ConfirmBy[#, AssociationQ], "results"], {__ ? AssociationQ}]]]
|>

(* ---- backend calibration: error map + device model + projection accessors ----
   Properties of the connection object, e.g. conn["ErrorMap", "Backend" -> "ibm_fez"].
   "Backend" defaults to the first available device, so conn["ErrorMap"] works bare.
   Each reads the two free, read-only metadata endpoints (RawBackendConfiguration /
   RawBackendProperties) through the connection's own auth: no circuit, no quantum
   time. Drawing + parse live in ErrorMap.m (iIBM* helpers). For many projections
   at once, conn["DeviceModel"] fetches once and carries every projection as a key.

   Each request uses Identity preprocessing (like "Backends") and resolves the
   backend inside the ExecuteFunction: the framework hands ExecuteFunction the raw
   parameters (a list of rules) as its first argument, which iIBMBackendFromParams /
   iIBMStyleRules read via Association. Every styling parameter is a STRING key. *)

params["ProcessedRequests", "ErrorMap"] = <|
    "ExecuteFunction" -> Function[iIBMErrorMap[iIBMFetchModel[iIBMBackendFromParams[#]], iIBMStyleRules[#]]],
    "SubmitFunction" -> "ExecuteFunction",
    "Parameters" -> {
        "Backend"             -> Automatic,
        "EdgeColorScheme"     -> "TemperatureMap",
        "VertexColorScheme"   -> "AvocadoColors",
        "VertexColorReversed" -> True,
        "VertexSizeReversed"  -> True,
        "EdgeThicknessMetric" -> "ZZ",
        "EdgeThicknessRange"  -> {4., 12.},
        "EdgeArrows"          -> True,
        "ShowQubitLabels"     -> True,
        "BackgroundColor"     -> Black
    },
    "RequiredParameters" -> {},
    "HiddenParameters" -> {},
    "PreprocessingFunction" -> Identity,
    "ExecuteResultProcessing" -> Identity
|>

params["ProcessedRequests", "DeviceModel"] = <|
    "ExecuteFunction" -> Function[iIBMFetchModel[iIBMBackendFromParams[#]]],
    "SubmitFunction" -> "ExecuteFunction",
    "Parameters" -> {"Backend" -> Automatic},
    "RequiredParameters" -> {},
    "HiddenParameters" -> {},
    "PreprocessingFunction" -> Identity,
    "ExecuteResultProcessing" -> Identity
|>

params["ProcessedRequests", "CouplingMap"] = <|
    "ExecuteFunction" -> Function[iIBMFetchModel[iIBMBackendFromParams[#]]["CouplingMap"]],
    "SubmitFunction" -> "ExecuteFunction",
    "Parameters" -> {"Backend" -> Automatic},
    "RequiredParameters" -> {},
    "HiddenParameters" -> {},
    "PreprocessingFunction" -> Identity,
    "ExecuteResultProcessing" -> Identity
|>

params["ProcessedRequests", "CZErrors"] = <|
    "ExecuteFunction" -> Function[iIBMCZErrors[iIBMFetchModel[iIBMBackendFromParams[#]]]],
    "SubmitFunction" -> "ExecuteFunction",
    "Parameters" -> {"Backend" -> Automatic},
    "RequiredParameters" -> {},
    "HiddenParameters" -> {},
    "PreprocessingFunction" -> Identity,
    "ExecuteResultProcessing" -> Identity
|>

params["ProcessedRequests", "ZZ"] = <|
    "ExecuteFunction" -> Function[iIBMZZ[iIBMFetchModel[iIBMBackendFromParams[#]]]],
    "SubmitFunction" -> "ExecuteFunction",
    "Parameters" -> {"Backend" -> Automatic},
    "RequiredParameters" -> {},
    "HiddenParameters" -> {},
    "PreprocessingFunction" -> Identity,
    "ExecuteResultProcessing" -> Identity
|>

params["ProcessedRequests", "GateErrors"] = <|
    "ExecuteFunction" -> Function[iIBMGateErrors[iIBMFetchModel[iIBMBackendFromParams[#]], Lookup[Association[#], "GateName", "cz"]]],
    "SubmitFunction" -> "ExecuteFunction",
    "Parameters" -> {"Backend" -> Automatic, "GateName" -> "cz"},
    "RequiredParameters" -> {},
    "HiddenParameters" -> {},
    "PreprocessingFunction" -> Identity,
    "ExecuteResultProcessing" -> Identity
|>

params["ProcessedRequests", "ReadoutErrors"] = <|
    "ExecuteFunction" -> Function[iIBMReadoutErrors[iIBMFetchModel[iIBMBackendFromParams[#]]]],
    "SubmitFunction" -> "ExecuteFunction",
    "Parameters" -> {"Backend" -> Automatic},
    "RequiredParameters" -> {},
    "HiddenParameters" -> {},
    "PreprocessingFunction" -> Identity,
    "ExecuteResultProcessing" -> Identity
|>

params["ProcessedRequests", "Coherence"] = <|
    "ExecuteFunction" -> Function[iIBMCoherence[iIBMFetchModel[iIBMBackendFromParams[#]]]],
    "SubmitFunction" -> "ExecuteFunction",
    "Parameters" -> {"Backend" -> Automatic},
    "RequiredParameters" -> {},
    "HiddenParameters" -> {},
    "PreprocessingFunction" -> Identity,
    "ExecuteResultProcessing" -> Identity
|>

params["ProcessedRequests", "Qubits"] = <|
    "ExecuteFunction" -> Function[iIBMFetchModel[iIBMBackendFromParams[#]]["Qubits"]],
    "SubmitFunction" -> "ExecuteFunction",
    "Parameters" -> {"Backend" -> Automatic},
    "RequiredParameters" -> {},
    "HiddenParameters" -> {},
    "PreprocessingFunction" -> Identity,
    "ExecuteResultProcessing" -> Identity
|>

End[]

EndPackage[]

SF`DefineServiceConnection[IBMQuantumPlatform`Private`params]