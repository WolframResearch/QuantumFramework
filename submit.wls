#!/usr/bin/env wolframscript
(* Submit the Wolfram/QuantumFramework paclet to the Wolfram Paclet Repository.

   The submission artifact is QuantumFramework/ResourceDefinition.nb.  The
   documentation pages beside it are generated from their literate markdown by
   build_docs.wls - rebuild those before submitting, never the notebooks.

   Run with:

       wolframscript -f submit.wls
*)

repoDir = DirectoryName[$InputFileName]

nbFile = FileNameJoin[{repoDir, "QuantumFramework", "ResourceDefinition.nb"}]

If[ ! FileExistsQ[nbFile],
    Print["FATAL: ", nbFile, " does not exist; run build_docs.wls first."];
    Exit[1]
]

Needs["DefinitionNotebookClient`"]

(* The account to submit from is WOLFRAM_CLOUD_USER and WOLFRAM_CLOUD_PASSWORD in
   the .env beside this script, which git ignores. Without them the kernel's own
   cloud session is used. *)

envFile = FileNameJoin[{repoDir, ".env"}]

env = If[ FileExistsQ[envFile],
    Association @ Cases[
        StringSplit[Select[ReadList[envFile, String], ! StringStartsQ[StringTrim[#], "#"] &], "=", 2],
        {key_, value_} /; StringMatchQ[StringTrim[key], RegularExpression["[A-Za-z_][A-Za-z0-9_]*"]] :>
            StringTrim[key] -> StringTrim[StringTrim[value], "\"" | "'"]
    ],
    <||>
]

(* Submitting needs the Internet, which a kernel can be set to refuse: it is
   allowed for the rest of this script alone, and the preference is left as it
   is. *)

Block[{$AllowInternet = True},

    If[ StringQ[env["WOLFRAM_CLOUD_USER"]] && StringQ[env["WOLFRAM_CLOUD_PASSWORD"]],
        If[TrueQ[$CloudConnected] && $WolframID =!= env["WOLFRAM_CLOUD_USER"], CloudDisconnect[]];
        If[! TrueQ[$CloudConnected], CloudConnect[env["WOLFRAM_CLOUD_USER"], env["WOLFRAM_CLOUD_PASSWORD"]]],
        If[! TrueQ[$CloudConnected], CloudConnect[]]
    ];

    If[ ! TrueQ[$CloudConnected],
        Print["FATAL: not cloud connected."];
        Exit[1]
    ];

    Print["submitting as ", $WolframID];

    (* The pre-submission check first: a submission with outstanding errors is
       rejected by the repository anyway, so failing here is strictly cheaper.
       Suggestions (such as PacletFilesChanged after a rebuild) do not block. *)

    check = UsingFrontEnd @ Block[{nb, r},
        nb = NotebookOpen[nbFile, Visible -> False];
        r = DefinitionNotebookClient`CheckDefinitionNotebook[nb];
        NotebookClose[nb];
        r
    ];

    issues = Normal[check];

    Scan[Print["  [", #Level, "] ", #Tag] &, issues];

    If[ MemberQ[issues[[All, "Level"]], "Error"],
        Print["FATAL: the definition notebook has errors; fix them before submitting."];
        Exit[1]
    ];

    (* Warnings, unlike errors, do not block a submission - the repository asks the
       submitter to acknowledge them, through a dialog that cannot open in a
       headless kernel (DialogInput::ninit, and the submission then reports
       $Failed against an expected True | $Canceled).  The acknowledgement is
       answered here instead, which is the same answer a submitter clicking through
       would give: the errors have already been checked for above and are fatal, so
       what reaches this point is warnings and suggestions only. *)

    result = UsingFrontEnd @ Block[
        {nb, r, DefinitionNotebookClient`ConfirmIssuesDialog = (True &)},
        nb = NotebookOpen[nbFile, Visible -> False];
        r = DefinitionNotebookClient`SubmitRepository[nb];
        NotebookClose[nb];
        r
    ];

    submission = FirstCase[result, _Success, Missing["NotFound"], Infinity];

    If[ MatchQ[submission, _Success],
        Print["submitted:"];
        KeyValueMap[Print["  ", #1, ": ", #2] &, KeyDrop[Last[submission], "Message"]],
        Print["SUBMISSION DID NOT SUCCEED:"];
        Print[result];
        Exit[1]
    ]
]
