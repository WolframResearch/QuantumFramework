Needs["ServiceFramework`" -> "SF`"]

BeginPackage["IBMQuantumPlatform`"]

Begin["`Private`"]


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
    "Icon" -> ImageTrim[Import[PacletObject["IBMQuantumPlatform"]["AssetLocation", "logo"]], ImageValuePositions[#, "Max"] &],
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

params["RawRequests", "RawExperiments"] = {
    "URL"                -> path["experiments"],
    "HTTPSMethod"        -> "GET",
    "Headers"            -> $commonHeaders
}

params["RawRequests", "RawExperimentResults"] = {
    "URL"                -> path["experiments", "results"],
    "HTTPSMethod"        -> "GET",
    "Headers"            -> $commonHeaders
}

params["RawRequests", "RawReservations"] = {
    "URL"                -> path["reservations"],
    "HTTPSMethod"        -> "GET",
    "Headers"            -> $commonHeaders
}

params["RawRequests", "RawReservationDetails"] = {
    "URL"                -> Function[URLBuild[{path["reservations"], #ReservationID}]],
    "HTTPSMethod"        -> "GET",
    "Headers"            -> $commonHeaders,
    "PathParameters"     -> {"ReservationID"},
    "RequiredParameters" -> {"ReservationID"}
}

params["RawRequests", "RawUsers"] = {
    "URL"                -> path["users"],
    "HTTPSMethod"        -> "GET",
    "Headers"            -> $commonHeaders
}

params["RawRequests", "RawUserDetails"] = {
    "URL"                -> Function[URLBuild[{path["users"], #UserID}]],
    "HTTPSMethod"        -> "GET",
    "Headers"            -> $commonHeaders,
    "PathParameters"     -> {"UserID"},
    "RequiredParameters" -> {"UserID"}
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
    "HTTPResponseProcessing" -> Function[ImportString[SF`ImportResponse[#], "RawJSON"]]
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
    "RequiredParameters" -> {"BackendID"}
}

params["RawRequests", "RawBackendConfiguration"] = {
    "URL"                -> Function[URLBuild[{path["backends"], #BackendID, "configuration"}]],
    "HTTPSMethod"        -> "GET",
    "Headers"            -> $commonHeaders,
    "PathParameters"     -> {"BackendID"},
    "RequiredParameters" -> {"BackendID"}
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

params["ProcessedRequests"] = <||>;

params["ProcessedRequests", "JobRun"] := <|
		"ExecuteFunction" -> SF`RequestExecute["RawJobRun"],
		"SubmitFunction" -> SF`RequestSubmit["RawJobRun"],
		"Parameters" -> {
			"QASM" -> Automatic, "Backend" -> Automatic, "ProgramID" -> Automatic
		},
        "RequiredParameters" -> {"QASM", "Backend"},
        "HiddenParameters" -> {},
        "PreprocessingFunction" -> Function[Enclose[<|
            "program_id" -> Replace[Lookup[#, "ProgramID"], Automatic -> "sampler"],
            "backend" -> Replace[Lookup[#, "Backend"], Automatic -> "ibmq_qasm_simulator"],
            "params" -> <|
                "pubs" -> {{Confirm @ Lookup[#, "QASM"], "Z"}},
                "options" -> <|"dynamical_decoupling" -> <|"enable" -> True|>|>,
                "version" -> 2,
                "resilience_level" -> 1
            |>
        |>]],
        "ExecuteResultProcessing" -> Identity
	|>

params["ProcessedRequests", "JobResults"] := <|
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

End[]

EndPackage[]

SF`DefineServiceConnection[IBMQuantumPlatform`Private`params]