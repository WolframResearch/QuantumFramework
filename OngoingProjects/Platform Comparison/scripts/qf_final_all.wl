(* QF stabilizer end-to-end harness: times the USER path (constructor + ApplyCircuit +
   forced readout) and the interpreted per-gate route, plus measurement/display/compose rows.
   Emits RESULT|<label>|<ms> and a naive total/k per-gate figure.

   The naive per-gate figure is NOT the engine's throughput: the total is dominated by
   fixed per-call work (encode, pack, canonical materialization, constructor), and the
   tableau kernel itself is a small fraction of it. bench_phases_qf.wls in this directory
   times those phases separately (with Stim and QC.jl mirrors sharing its output schema);
   quote cross-engine ratios from the phase rows, never from the naive figure alone.

   Op streams are read from this script's own directory (no /tmp staging).
   Run: wolframscript -file qf_final_all.wl   (exit 1 on any failed validity gate)
   Expect run-to-run drift of ~10% with machine load. *)

$here = DirectoryName[$InputFileName];
$sizes = {100, 500, 1000};
$measSizes = {100, 500};   (* measurement rows stop at 500: the Pauli-string branch is the cost *)
$composeSize = 100;        (* composition is O(n^3) (Compose.m:24); 926 ms already at n=100 *)
PacletDirectoryLoad[FileNameDrop[$here, -3]];   (* the tree this script lives in, not a fixed path *)
Needs["Wolfram`QuantumFramework`"];

pr[label_, ms_] := WriteString["stdout", "RESULT|", label, "|", Round[ms, 0.001], "\n"];
prRate[label_, ms_, gates_] :=
  WriteString["stdout", "RESULT|", label, "_ns_per_gate|", Round[1.*^6 ms/gates, 0.1], "\n"];
encode = Symbol["Wolfram`QuantumFramework`PackageScope`encodeStabilizerGates"];

SetAttributes[tm, HoldFirst];
tm[e_, reps_ : 3] := 1000. Min@Table[ClearSystemCache[]; First@AbsoluteTiming[e], {reps}];

(* Total over the engine's own catalog ($stabilizerCliffordNames: H S X Y Z V T CNOT CX CZ SWAP),
   not a three-way Switch. An unrecognized spec aborts loudly; left to fall through it would
   produce an inert expression that still times and still prints a RESULT line. *)
$gate1 = {"H", "S", "X", "Y", "Z", "V", "T"};
$gate2 = {"CNOT", "CX", "CZ", "SWAP"};
toSpec[ops_] := Map[Replace[{
    {g_String /; MemberQ[$gate1, g], q_Integer} :> g -> q + 1,
    {g_String /; MemberQ[$gate2, g], a_Integer, b_Integer} :> g -> {a + 1, b + 1},
    bad_ :> Throw[bad, "badGateSpec"]}], ops];

(* Two-qubit targets must be splatted: the per-gate property form is ps["CNOT", a, b],
   not ps["CNOT", {a, b}]. *)
applyFold[n_, specs_] :=
  Fold[#1[First[#2], Sequence @@ Flatten[{Last[#2]}]] &, PauliStabilizer[n], specs];

stream[n_] := Import[FileNameJoin[{$here, "stab_ops_" <> ToString[n] <> ".json"}]];

(* Deeply scrambled reference state: ~20n gates, enough that measurement lands on the
   random-outcome branch rather than the deterministic fast path. A shallow state makes the
   measurement rows below silently compare two different code paths across n. *)
SeedRandom[3];
mk[n_] := PauliStabilizer[n]["ApplyCircuit", toSpec@Table[
   Switch[Mod[i, 3],
     0, {"H", RandomInteger[{0, n - 1}]},
     1, {"S", RandomInteger[{0, n - 1}]},
     2, With[{a = RandomInteger[{0, n - 1}]},
          {"CNOT", a, Mod[a + RandomInteger[{1, n - 1}], n]}]],
   {i, 20 n}]];

PauliStabilizer[4]["ApplyCircuit", {"H" -> 1}];   (* warm the compiled kernel *)

Catch[
  (* Correctness gate before any timing: the bulk route and the per-gate fold are independent
     implementations of the same evolution, so they must agree. A benchmark with no assertion
     can time a wrong computation indefinitely. *)
  With[{specs = toSpec@Take[stream[100], 200]},
    With[{bulk = PauliStabilizer[100]["ApplyCircuit", specs], fold = applyFold[100, specs]},
      If[bulk["Tableau"] =!= fold["Tableau"] || bulk["Signs"] =!= fold["Signs"],
        WriteString["stdout", "FAIL|route_disagreement\n"]; Exit[1]];
      If[! StabilizerStateQ[bulk], WriteString["stdout", "FAIL|not_stabilizer\n"]; Exit[1]]]];
  (* Exactly solvable checks: the fresh register (constructor path) and a GHZ circuit
     through the COMPILED path itself; route agreement alone would pass a bug shared by
     both routes. *)
  If[PauliStabilizer[6]["Stabilizers"] =!= Table[
       StringJoin@ReplacePart[ConstantArray["I", 6], i -> "Z"], {i, 6}],
    WriteString["stdout", "FAIL|closed_form_Z_generators\n"]; Exit[1]];
  If[PauliStabilizer[3]["ApplyCircuit",
       {"H" -> 1, "CNOT" -> {1, 2}, "CNOT" -> {1, 3}}]["Stabilizers"] =!=
     {"XXX", "ZZI", "ZIZ"},
    WriteString["stdout", "FAIL|ghz_closed_form\n"]; Exit[1]];
  WriteString["stdout", "CHECK|routes_agree_and_state_valid|True\n"];

  Scan[
    Function[n,
      With[{ops = stream[n]},
        With[{specs = toSpec[ops], m = Length[ops]},
          (* toSpec has validated every op in THIS stream; additionally assert the compiled
             encoder accepts it, since ApplyCircuit silently falls back to the interpreted
             per-gate fold on ANY encoder rejection (Compiled.m:244) and the row would then
             time the wrong route under a compiled label. *)
          If[! ListQ[encode[specs]],
            WriteString["stdout", "FAIL|encoder_rejected_stream_n", n, "\n"]; Exit[1]];
          With[{res = PauliStabilizer[n]["ApplyCircuit", specs]},
            If[! StabilizerStateQ[res],
              WriteString["stdout", "FAIL|timed_object_not_stabilizer_n", n, "\n"]; Exit[1]]];
          With[{ac = tm[PauliStabilizer[n]["ApplyCircuit", specs]["Signs"]]},
            pr["gates_n" <> ToString[n] <> "_AC_ms", ac];
            prRate["gates_n" <> ToString[n] <> "_AC", ac, m]];
          With[{pg = tm[applyFold[n, specs]["Signs"]]},
            pr["gates_n" <> ToString[n] <> "_pergate_ms", pg];
            prRate["gates_n" <> ToString[n] <> "_pergate", pg, m]]]]],
    $sizes];

  Scan[pr["ctor_n" <> ToString[#] <> "_ms", tm[PauliStabilizer[#]]] &, $sizes];

  Scan[
    Function[n,
      With[{ps = mk[n], zstr = StringJoin@ConstantArray["Z", n]},
        (* Report the measurement branch: 1 outcome key = deterministic O(n) fast path,
           2 keys = random-outcome rowsum branch. Timings across n are only comparable
           within the same branch. *)
        With[{keys = Length@Keys@ps["M", 1]},
          WriteString["stdout", "CHECK|measZ_n", n, "_branch_keys|", keys, "\n"];
          If[keys =!= 2,
            WriteString["stdout", "FAIL|meas_branch_not_random_n", n, "\n"]; Exit[1]]];
        pr["measZ_n" <> ToString[n] <> "_ms", tm[ps["M", 1]]];
        pr["measPauli_n" <> ToString[n] <> "_ms", tm[ps["M", zstr]]];
        pr["stabs_n" <> ToString[n] <> "_ms", tm[ps["Stabilizers"]]]]],
    $measSizes];

  With[{a = mk[$composeSize], b = mk[$composeSize]},
    pr["compose_n" <> ToString[$composeSize] <> "_ms", tm[a[b]]]];
  WriteString["stdout", "DONE\n"],

  "badGateSpec",
  (WriteString["stdout", "FAIL|unrecognized_gate_spec|", ToString[#1], "\n"]; Exit[1]) &]
