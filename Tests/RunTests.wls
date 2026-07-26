#!/usr/bin/env wolframscript

(* Test runner for QuantumFramework.

   Usage: wolframscript -file Tests/RunTests.wls [pattern]
   Examples:
     wolframscript -file Tests/RunTests.wls            run all .wlt files
     wolframscript -file Tests/RunTests.wls Distance   run only files whose
                                                       name matches "Distance"

   Output carries two interleaved streams. Lines that begin in column 1 with an
   uppercase tag are machine readable, hold no dates or timings, and are emitted
   in sorted file order, so two runs diff cleanly:

     PACLET  <tab> location
     TEST    <tab> file <tab> TestID <tab> outcome
     FILE    <tab> file <tab> succeeded <tab> total
     OUTCOME <tab> outcome <tab> count
     TOTAL   <tab> succeeded <tab> total
     ALERT   <tab> reason <tab> count
     STATUS  <tab> GREEN or RED

   STATUS carries the verdict, and it is not implied by the rest: a load-time
   message, a message raised outside the tests, or a paclet resolved outside
   this repository each leave every TEST and TOTAL line looking perfect while
   the run is not to be believed. Diff TEST lines to see what moved; read
   STATUS to see whether the run counts at all.

   Everything else is indented prose for a human reader. To compare a run
   against another commit:

     wolframscript -file Tests/RunTests.wls | grep '^TEST' | sort > after.tsv
     diff before.tsv after.tsv

   A completed run always emits a TOTAL line and a STATUS line. That is the
   contract, and a caller that cares should check it:

     out=$(wolframscript -file Tests/RunTests.wls) || exit 1
     grep -q '^TOTAL' <<<"$out" || exit 1

   An abort inside a test file is caught and charged to that file, but Exit[]
   inside one cannot be: it takes the kernel down where it stands, before any
   total is written, and the shell then sees status 0. The absent TOTAL line is
   the only evidence left. Running one kernel per test file is the robust cure,
   at the price of paying a multi-second paclet load 32 times over.

   Helper names here (emit, row, flat, clip, shown) share the Global` context
   with every .wlt the runner reads, so a test file that assigns one of them
   would rewrite the report rather than fail. None of the current 32 do.

   TestIDs are written in InputForm, which is load bearing rather than
   decorative. An ID built as "Dim-JX-" <> ToString[j] with half-integer j
   embeds the newlines of a two-dimensional fraction: ToString[1/2] is the five
   characters "1\n-\n2". Written raw, one such test spans three lines and every
   line-oriented reader downstream loses sync for the rest of the file. *)

(* ------------------------------------------------------------------------ *)
(* Output.  wolframscript -file prints nothing by itself, so every line the   *)
(* runner produces has to be written to $Output explicitly.                   *)
(* ------------------------------------------------------------------------ *)

(* TestReport's progress monitor writes elapsed times into this same stream,
   which would put a non-reproducible line between two test results. *)
$ProgressReporting = False;

emit[] := WriteString[$Output, "\n"];
emit[line_String] := WriteString[$Output, line, "\n"];

(* One line, whatever the expression: InputForm escapes embedded newlines. *)
flat[expr_] := ToString[expr, InputForm, PageWidth -> Infinity];

(* Detail lines can carry a whole matrix; keep them readable but say when
   something was cut rather than dropping characters silently. *)
clip[s_String] := If[StringLength[s] <= 800, s,
    StringTake[s, 800] <> "... (" <> ToString[StringLength[s]] <> " chars total)"
];

(* Fields are strings today.  Converting anything else rather than pattern
   matching on String keeps a stray non-string from dropping a whole row in
   silence, which would leave the TEST lines disagreeing with TOTAL. *)
field[s_String] := s;
field[e_] := flat[e];
row[tag_String, fields___] := emit[StringRiffle[field /@ {tag, fields}, "\t"]];

$rule = StringRepeat["-", 71];
$heavyRule = StringRepeat["=", 71];

(* $MessageList entries are held; strip the wrapper textually.  Releasing the
   HoldForm would re-issue the message. *)
messageTag[m_] := StringTrim[flat[m], "HoldForm[" | "]"];

(* Test outputs arrive held, and a mismatch is unreadable as HoldForm[1] against
   HoldForm[999].  Safe here because these are already evaluated values; the same
   move on a message would re-issue it, so ActualMessages keeps its wrapper. *)
shown[expr_] := clip[flat[Replace[expr, HoldForm[x_] :> x]]];

(* ------------------------------------------------------------------------ *)
(* Paclet                                                                     *)
(* ------------------------------------------------------------------------ *)

$rootDir = AbsoluteFileName[DirectoryName[$InputFileName]];
$pacletCandidate = FileNameJoin[{ParentDirectory[$rootDir], "QuantumFramework"}];

(* Asked before AbsoluteFileName rather than after, so the path that is missing
   can be named.  AbsoluteFileName answers $Failed and raises its own message,
   and $Failed is not a String, so emit and row would then decline to print the
   header and PACLET lines at all instead of reporting anything. *)
If[! DirectoryQ[$pacletCandidate],
    emit["  ERROR: no QuantumFramework directory beside " <> $rootDir];
    emit["    looked for: " <> $pacletCandidate];
    row["STATUS", "RED"];
    Exit[1]
];

$pacletDir = AbsoluteFileName[$pacletCandidate];

(* A directory added with PacletDirectoryLoad outranks an installed paclet of
   the same name, so the working tree is tested without touching whatever the
   user has installed. *)
PacletDirectoryLoad[$pacletDir];
Needs["Wolfram`QuantumFramework`"];

(* Anything the load itself said.  A paclet that failed to load cleanly must
   not be able to produce an all-green run. *)
$loadMessages = $MessageList;

(* On a path boundary, so that a sibling directory named QuantumFramework-old is
   not read as living inside QuantumFramework. *)
insideRepo[p_String] := p === $pacletDir || StringStartsQ[p, $pacletDir <> $PathnameSeparator];
insideRepo[_] := False;

(* Ground truth for which copy won, rather than an assumption about priority. *)
$loadedFile = FindFile["Wolfram`QuantumFramework`"];
$wrongCopy = ! insideRepo[$loadedFile];

$shadowing = Select[
    Cases[PacletFind["Wolfram/QuantumFramework"],
        p_ /; StringQ[p["Location"]] :> AbsoluteFileName[p["Location"]]],
    ! insideRepo[#] &
];

(* ------------------------------------------------------------------------ *)
(* Test files                                                                 *)
(* ------------------------------------------------------------------------ *)

$pattern = If[Length[$ScriptCommandLine] > 1, $ScriptCommandLine[[2]], ""];

relative[f_String] :=
    StringTrim[StringReplace[f, StartOfString ~~ $rootDir -> ""], $PathnameSeparator];

(* Sorted, so that file order is a property of the repository rather than of
   the filesystem, and two runs line up. *)
$testFiles = Select[
    Sort @ FileNames["*.wlt", $rootDir, Infinity],
    $pattern === "" || StringContainsQ[FileNameTake[#], $pattern, IgnoreCase -> True] &
];

(* ------------------------------------------------------------------------ *)
(* Header                                                                     *)
(* ------------------------------------------------------------------------ *)

emit[$heavyRule];
emit["  QuantumFramework Test Coverage Report"];
emit[$heavyRule];
emit["  Date:       " <> DateString[]];
emit["  Paclet:     " <> $pacletDir];
emit["  Loaded:     " <> If[StringQ[$loadedFile], $loadedFile, flat[$loadedFile]]];
emit["  Test files: " <> ToString[Length[$testFiles]] <>
    If[$pattern === "", "", " matching \"" <> $pattern <> "\""]];

row["PACLET", $pacletDir];

If[$loadMessages =!= {},
    emit[];
    emit["  MESSAGES DURING LOAD (this run cannot be considered green):"];
    Scan[emit["    " <> messageTag[#]] &, $loadMessages];
];

If[$shadowing =!= {},
    emit[];
    emit["  NOTE: an installed Wolfram/QuantumFramework paclet is also visible:"];
    Scan[emit["    " <> #] &, $shadowing];
    emit["    Nothing has been uninstalled.  PacletDirectoryLoad gives the repo"];
    emit["    copy priority, so this is usually harmless.  To remove it yourself:"];
    emit["      PacletUninstall /@ PacletFind[\"Wolfram/QuantumFramework\"]"];
];

If[$wrongCopy,
    emit[];
    emit["  ERROR: Wolfram`QuantumFramework` did not resolve inside this repository."];
    emit["    resolved to: " <> flat[$loadedFile]];
    emit["    expected under: " <> $pacletDir];
    emit["    The results below would describe some other copy of the paclet."];
];

If[Length[$testFiles] == 0,
    emit[];
    emit["  No test files found" <>
        If[$pattern === "", "", " matching \"" <> $pattern <> "\""] <> "."];
    Exit[1]
];

(* ------------------------------------------------------------------------ *)
(* Run                                                                        *)
(* ------------------------------------------------------------------------ *)

emitResult[name_String, r_] := (
    row["TEST", name, flat[r["TestID"]], r["Outcome"]];
    If[r["Outcome"] =!= "Success",
        emit["      expected: " <> shown[r["ExpectedOutput"]]];
        emit["      actual:   " <> shown[r["ActualOutput"]]];
        (* ExpectedMessages comes back as HoldForm[{}] rather than {}, so testing
           it against {} alone would tack an empty messages line onto every
           ordinary wrong-result failure. *)
        If[r["ActualMessages"] =!= {} || ! MatchQ[r["ExpectedMessages"], {} | HoldForm[{}]],
            emit["      messages: " <> clip[flat[r["ActualMessages"]]]]
        ];
    ]
);

(* Every outcome TestReport reports is counted, including any this runner does
   not anticipate: the census is built from the outcomes themselves rather than
   from a fixed list of buckets that a new outcome could fall through. *)
runFile[file_] := With[{name = relative[file]},
    emit[];
    emit[$rule];
    emit["  " <> name];
    emit[$rule];
    Enclose[
        (* CheckAbort, because an Abort[] inside a test file stops TestReport from
           ever returning: without it the run ends mid-file with no summary, no
           totals and a status of 0.  Enclose alone does not catch an abort. *)
        With[{report = ConfirmMatch[CheckAbort[TestReport[file], $Aborted],
                _TestReportObject, "aborted, unreadable, or missing"]},
            With[{results = Values[report["TestResults"]]},
                Scan[emitResult[name, #] &, results];
                With[{census = Counts[#["Outcome"] & /@ results],
                      aborted = TrueQ[report["Aborted"]],
                      skipped = Length[report["TestsNotEvaluatedKeys"]]},
                    With[{total = Length[results], success = Lookup[census, "Success", 0]},
                        (* TestReport hands back a report even for a file that is
                           empty, unparseable, or absent, so a silent 0/0 is the
                           shape a broken file takes.  It is not a pass. *)
                        ConfirmAssert[total > 0, "file produced no tests"];
                        row["FILE", name, ToString[success], ToString[total]];
                        emit["  " <> ToString[success] <> "/" <> ToString[total] <> " succeeded"];
                        KeyValueMap[
                            emit["  " <> #1 <> ": " <> ToString[#2]] &,
                            KeyDrop[census, "Success"]
                        ];
                        If[aborted, emit["  ABORTED before the end of the file"]];
                        If[skipped > 0, emit["  not evaluated: " <> ToString[skipped]]];
                        <|"File" -> name, "Census" -> census, "Aborted" -> aborted,
                          "Skipped" -> skipped|>
                    ]
                ]
            ]
        ],
        (* A file that cannot even be turned into a report counts as one
           failure of its own, so it can never pass by being unreadable. *)
        Function[failure,
            emit["  RUNNER ERROR: " <>
                clip[Replace[failure["MessageTemplate"], Except[_String] :> flat[failure]]]];
            row["FILE", name, "0", "1"];
            (* Not flagged as aborted: an abort is only one of the ways to land
               here, alongside an empty or unparseable file, and the RunnerError
               in the census is already what makes the run red. *)
            <|"File" -> name, "Census" -> <|"RunnerError" -> 1|>, "Aborted" -> False,
              "Skipped" -> 0|>
        ]
    ]
];

$summaries = runFile /@ $testFiles;

(* ------------------------------------------------------------------------ *)
(* Summary                                                                    *)
(* ------------------------------------------------------------------------ *)

$census = Merge[#["Census"] & /@ $summaries, Total];
$totalTests = Total[Values[$census]];
$totalSuccess = Lookup[$census, "Success", 0];
$totalOther = $totalTests - $totalSuccess;

(* Messages raised after the load, by the runner or by anything a test file did
   outside VerificationTest, which absorbs its own. *)
$lateMessages = Drop[$MessageList, UpTo[Length[$loadMessages]]];

emit[];
emit[$heavyRule];
emit["  Summary"];
emit[$heavyRule];

Scan[
    With[{n = Total[Values[#["Census"]]], s = Lookup[#["Census"], "Success", 0]},
        emit["  " <> If[n === s && ! #["Aborted"], "PASS", "FAIL"] <> " " <>
            #["File"] <> ": " <> ToString[s] <> "/" <> ToString[n]]
    ] &,
    $summaries
];

emit[];
(* Floor, not Round: 2328 of 2329 is 99.957%, and a run with a failure in it
   must not be able to display 100%. *)
emit["  Total: " <> ToString[$totalSuccess] <> "/" <> ToString[$totalTests] <>
    " succeeded" <> If[$totalTests == 0, "",
        " (" <> ToString[NumberForm[Floor[100. $totalSuccess / $totalTests, 0.1],
            {Infinity, 1}]] <> "%)"]];
KeyValueMap[emit["  " <> #1 <> ": " <> ToString[#2]] &, KeyDrop[$census, "Success"]];

If[$lateMessages =!= {},
    emit[];
    emit["  MESSAGES OUTSIDE THE TESTS:"];
    Scan[emit["    " <> messageTag[#]] &, $lateMessages];
];

(* A green run means every test succeeded AND every test actually ran AND the
   paclet under test is the one in this repository AND nothing said anything on
   the way.  A test that never evaluated is not a test that passed. *)
$alerts = Join[
    If[$totalOther > 0, {{"NonSuccess", $totalOther}}, {}],
    If[$wrongCopy, {{"PacletOutsideRepo", 1}}, {}],
    If[$loadMessages =!= {}, {{"LoadMessages", Length[$loadMessages]}}, {}],
    If[$lateMessages =!= {}, {{"MessagesOutsideTests", Length[$lateMessages]}}, {}],
    With[{n = Count[$summaries, s_ /; s["Aborted"]]}, If[n > 0, {{"AbortedFiles", n}}, {}]],
    With[{n = Total[#["Skipped"] & /@ $summaries]}, If[n > 0, {{"NotEvaluated", n}}, {}]]
];

(* The verdict has to reach the machine-readable stream too.  Every reason a run
   is not to be believed other than a failing test leaves the TEST and TOTAL
   lines looking perfect, so a diff of those alone would call a red run green. *)
KeyValueMap[row["OUTCOME", #1, ToString[#2]] &, $census];
row["TOTAL", ToString[$totalSuccess], ToString[$totalTests]];
Scan[row["ALERT", First[#], ToString[Last[#]]] &, $alerts];
row["STATUS", If[$alerts === {}, "GREEN", "RED"]];
emit[$heavyRule];

Exit[If[$alerts === {}, 0, 1]]
