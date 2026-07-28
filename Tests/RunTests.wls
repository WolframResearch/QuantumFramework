#!/usr/bin/env wolframscript

(* Test runner for QuantumFramework.

   Usage: wolframscript -file Tests/RunTests.wls [pattern]
   Examples:
     wolframscript -file Tests/RunTests.wls            run all .wlt files
     wolframscript -file Tests/RunTests.wls Distance   run only files whose
                                                       name matches "Distance"

   Output carries two interleaved streams. The machine-readable one is the seven
   tags below, each starting in column 1 and followed by a tab; it holds no dates
   or timings and comes in sorted file order, so two runs diff cleanly. Match the
   tag and the tab, not merely column 1: kernel messages also land on stdout under
   wolframscript, at column 1, and Name::tag reads like a tag if you let it.

     PACLET  <tab> location
     TEST    <tab> file <tab> TestID <tab> outcome
     FILE    <tab> file <tab> succeeded <tab> total
     OUTCOME <tab> outcome <tab> count
     TOTAL   <tab> succeeded <tab> total
     ALERT   <tab> reason <tab> count
     STATUS  <tab> GREEN or RED

   Diff TEST lines to see what moved; read STATUS to see whether the run counts
   at all. STATUS is not implied by the other lines, for the reason given where
   it is computed.

   Everything else is indented prose for a human reader. To compare a run
   against another commit:

     wolframscript -file Tests/RunTests.wls | grep -E '^TEST\t' | sort > after.tsv
     diff before.tsv after.tsv

   Every test gets one TEST line, including a file charged as unreadable, so the
   TEST count and TOTAL always agree.

   A completed run always emits a TOTAL line and a STATUS line. That is the
   contract, and a caller that cares should check it:

     out=$(wolframscript -file Tests/RunTests.wls) || exit 1
     grep -q $'^TOTAL\t' <<<"$out" || exit 1

   An abort inside a test file is caught and charged to that file, but Exit[]
   inside one cannot be: it takes the kernel down where it stands, before any
   total is written, and the shell then sees status 0. The absent TOTAL line is
   the only evidence left. Running one kernel per test file is the robust cure,
   at the price of paying a multi-second paclet load 32 times over.

   Every name defined here, the helpers and the $-globals alike, shares the
   Global` context with every .wlt the runner reads, because that is the context
   TestReport reads a test file in. A test file assigning one of them, say
   emit = "x", installs an OwnValue ahead of the DownValue, and the effect ranges
   from the report going quiet to something worse: clobbering emitResult keeps the
   report looking normal while deleting every per-test line from it. None of the
   current 32 do. The TestLineMismatch alert catches a test line that was counted
   but never written, which covers the accidental cases; it is not armour against
   a file that redefines the writer itself, and the caller's TOTAL check is the
   backstop there. Defining the helpers inside a private context would close this
   properly, and does not disturb the context each .wlt is read in.

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

(* Each line is recorded as it is written, not as it is requested, so the tally
   at the end counts lines that really went out rather than calls that intended
   to write one. *)
emit[line_String] := (WriteString[$Output, line, "\n"]; Sow[line, "line"]);

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

(* A $MessageList entry is a held MessageName, so releasing it would yield the
   message text where the name is wanted.  Stripping the wrapper textually keeps
   the name.  (Releasing is harmless here, a bare MessageName cannot fire; the
   thing that must not be released is an ActualMessages entry, which holds a
   Message call and does fire.) *)
messageTag[HoldForm[MessageName[s_, t_String]]] := SymbolName[Unevaluated[s]] <> "::" <> t;
messageTag[m_] := StringTrim[flat[m], "HoldForm[" | "]"];

(* Reported by name: the Message call inside carries every argument, and a whole
   HoldCompleteForm of the offending expression buries the one useful word.  Named
   the same way as the load-time and late messages above, so one object does not
   get two renderings in one report. *)
messageNames[{}] := "none";
messageNames[msgs_List] :=
    StringRiffle[messageTag /@ Replace[msgs, HoldForm[Message[m_, ___]] :> HoldForm[m], {1}], ", "];
messageNames[other_] := flat[other];

(* A broken file's report says why in RuntimeFailures.  Named by tag as well as
   parameters: TestReport::rnterr puts the whole sentence in its parameters, but
   TestReport::fnfnd puts only the path there and keeps "not found" in the
   template, so parameters alone would report a bare path and say nothing. *)
breakReason[failures : {__}] := StringRiffle[
    Function[f,
        f["Tag"] <> ": " <> StringRiffle[
            field /@ Flatten[{Replace[f["MessageParameters"], Except[_List] :> flat[f]]}], "; "]
    ] /@ failures,
    "; "];
breakReason[_] := "file produced no tests";

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

(* The paclet's loader registers PackageScope as a package, so it can be put on
   the context path here and the tests can name internal symbols directly rather
   than spelling the full context at each use.  Wolfram`Arrays` joins it: the
   array-container contract is part of what the tests assert about. *)
Needs["Wolfram`QuantumFramework`PackageScope`"];
Needs["Wolfram`Arrays`"];

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
    row["STATUS", "RED"];
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
        (* Both sides.  A MessagesFailure whose declared message never arrived has
           an empty actual list, and printing that side alone leaves every visible
           field agreeing with a pass while the outcome says otherwise. *)
        If[r["ActualMessages"] =!= {} || ! MatchQ[r["ExpectedMessages"], {} | HoldForm[{}]],
            emit["      messages: expected " <>
                clip[messageNames[Replace[r["ExpectedMessages"], HoldForm[l_] :> l]]] <>
                ", got " <> clip[messageNames[r["ActualMessages"]]]]
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
             {results = report["Results"]},
             (* Aborted and TestsNotEvaluatedKeys are internal keys, not documented
                properties, and MUnit only ever writes them from the progress
                panel's Stop button.  With progress reporting off they are inert
                here, kept as belt and braces rather than as a live signal: an
                interrupted file arrives through the RunnerError branch instead. *)
             {census = Counts[#["Outcome"] & /@ results],
              aborted = TrueQ[report["Aborted"]],
              skipped = Length[report["TestsNotEvaluatedKeys"]],
              broke = report["RuntimeFailures"]},
             {total = Length[results], success = Lookup[census, "Success", 0]},
            Scan[emitResult[name, #] &, results];
            (* TestReport hands back a report even for a file that is empty,
               unparseable, or absent, so a silent 0/0 is the shape a broken file
               takes.  It is not a pass, and the report knows why: RuntimeFailures
               carries the syntax error or the missing path. *)
            ConfirmAssert[total > 0, breakReason[broke]];
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
        ],
        (* A file that cannot even be turned into a report counts as one failure
           of its own, so it can never pass by being unreadable.  It gets a TEST
           line as well as a FILE line: the header invites a reader to check the
           TEST count against TOTAL, and a file charged one test with no test line
           would break that identity. *)
        Function[failure,
            emit["  RUNNER ERROR: " <>
                clip[Replace[failure["MessageTemplate"], Except[_String] :> flat[failure]]]];
            row["TEST", name, flat[name], "RunnerError"];
            row["FILE", name, "0", "1"];
            (* Not flagged as aborted: an abort is only one of the ways to land
               here, alongside an empty or unparseable file, and the RunnerError
               in the census is already what makes the run red. *)
            <|"File" -> name, "Census" -> <|"RunnerError" -> 1|>, "Aborted" -> False,
              "Skipped" -> 0|>
        ]
    ]
];

{$summaries, {$writtenLines}} = Reap[runFile /@ $testFiles, "line"];
$testLines = Count[$writtenLines, _String?(StringStartsQ[#, "TEST\t"] &)];

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
    With[{census = #["Census"], aborted = #["Aborted"], file = #["File"]},
         {n = Total[Values[census]], s = Lookup[census, "Success", 0]},
        emit["  " <> If[n === s && ! aborted, "PASS", "FAIL"] <> " " <>
            file <> ": " <> ToString[s] <> "/" <> ToString[n]]
    ] &,
    $summaries
];

emit[];
(* Floor, not Round: 2328 of 2329 is 99.957%, and a run with a failure in it must
   not be able to display 100%.  Floored on the exact rational and only then made
   numeric, because flooring a machine float lands one tenth low whenever the true
   percentage sits exactly on a tenth: 23 of 40 is 57.5%, not 57.4%. *)
emit["  Total: " <> ToString[$totalSuccess] <> "/" <> ToString[$totalTests] <>
    " succeeded" <> If[$totalTests == 0, "",
        " (" <> ToString[NumberForm[N[Floor[100 $totalSuccess / $totalTests, 1/10]],
            {Infinity, 1}]] <> "%)"]];
KeyValueMap[emit["  " <> #1 <> ": " <> ToString[#2]] &, KeyDrop[$census, "Success"]];

If[$lateMessages =!= {},
    emit[];
    emit["  MESSAGES OUTSIDE THE TESTS:"];
    Scan[emit["    " <> messageTag[#]] &, $lateMessages];
];

(* A green run means every test succeeded AND every test actually ran AND every
   test left a line behind AND the paclet under test is the one in this
   repository AND nothing said anything on the way.  A test that never evaluated
   is not a test that passed. *)
alert[name_String, n_] := If[TrueQ[n > 0], <|name -> n|>, <||>];

$alerts = Join[
    alert["NonSuccess", $totalOther],
    alert["PacletOutsideRepo", Boole[$wrongCopy]],
    alert["LoadMessages", Length[$loadMessages]],
    alert["MessagesOutsideTests", Length[$lateMessages]],
    alert["AbortedFiles", Count[$summaries, s_ /; s["Aborted"]]],
    alert["NotEvaluated", Total[#["Skipped"] & /@ $summaries]],
    (* The header promises one TEST line per counted test.  Counting the lines
       actually written is what keeps that a fact rather than an intention: any
       path that drops a row, including a test file that assigns one of the
       helper names out from under this script, turns the run red here instead
       of reporting a full green tally with the per-test evidence missing. *)
    alert["TestLineMismatch", Boole[$testLines =!= $totalTests]]
];

(* The verdict has to reach the machine-readable stream too. *)
KeyValueMap[row["OUTCOME", #1, ToString[#2]] &, $census];
row["TOTAL", ToString[$totalSuccess], ToString[$totalTests]];
KeyValueMap[row["ALERT", #1, ToString[#2]] &, $alerts];
row["STATUS", If[$alerts === <||>, "GREEN", "RED"]];
emit[$heavyRule];

Exit[If[$alerts === <||>, 0, 1]]
