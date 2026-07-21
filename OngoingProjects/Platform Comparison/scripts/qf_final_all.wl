(* QF stabilizer end-to-end harness: times the USER path (constructor + ApplyCircuit +
   forced readout) and the interpreted per-gate route, plus measurement/display/compose rows.
   Emits RESULT|<label>|<ms> and a naive total/k per-gate figure.

   The naive per-gate figure is NOT the engine's throughput: the total is dominated by
   fixed per-call work (encode, pack, canonical materialization, constructor), and the
   tableau kernel itself is a small fraction of it. bench_phases_qf.wls in this directory
   times those phases separately (with Stim and QC.jl mirrors sharing its output schema);
   quote cross-engine ratios from the phase rows, never from the naive figure alone.

   Role: correctness-gated timings of the USER path only. The scaling-law gates
   (simulate ~ k, materialize ~ n^2), the n=1024 sparse-fallback regression
   guard, and the V-rejection edge probe live in bench_phases_qf.wls; this file
   deliberately does not duplicate them.

   Op streams are read from this script's own directory (no /tmp staging).
   Run: wolframscript -file qf_final_all.wl   (exit 1 on any failed validity gate)
   Expect run-to-run drift of ~10% with machine load. *)

$here = DirectoryName[$InputFileName];
$sizes = {100, 500, 1000};
$measSizes = {100, 500};   (* measurement rows stop at 500: the Pauli-string branch is the cost *)
$composeSize = 100;        (* composition is O(n^3) (Compose.m:24); 926 ms already at n=100 *)
PacletDirectoryLoad[FileNameDrop[$here, -3]];   (* the tree this script lives in, not a fixed path *)
Needs["Wolfram`QuantumFramework`"];
Get[FileNameJoin[{$here, "bench_qf_common.wl"}]];  (* fail, tm (min over 7), stream, gram, perGateFold *)

pr[label_, ms_] := WriteString["stdout", "RESULT|", label, "|", Round[ms, 0.001], "\n"];
prRate[label_, ms_, gates_] :=
  WriteString["stdout", "RESULT|", label, "_ns_per_gate|", Round[1.*^6 ms/gates, 0.1], "\n"];
encode = Symbol["Wolfram`QuantumFramework`PackageScope`encodeStabilizerGates"];

(* Total over the engine's own catalog ($stabilizerCliffordNames: H S X Y Z V T CNOT CX CZ SWAP),
   not a three-way Switch. An unrecognized spec aborts loudly; left to fall through it would
   produce an inert expression that still times and still prints a RESULT line. *)
$gate1 = {"H", "S", "X", "Y", "Z", "V", "T"};
$gate2 = {"CNOT", "CX", "CZ", "SWAP"};
toSpec[ops_] := Map[Replace[{
    {g_String /; MemberQ[$gate1, g], q_Integer} :> g -> q + 1,
    {g_String /; MemberQ[$gate2, g], a_Integer, b_Integer} :> g -> {a + 1, b + 1},
    bad_ :> fail["unrecognized_gate_spec: " <> ToString[bad]]}], ops];

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

(* Correctness gates before any timing (ghzGate/routeGate from
   bench_qf_common.wl): a benchmark with no assertion can time a wrong
   computation indefinitely. *)
With[{specs = toSpec@Take[stream[100], 200]},
  routeGate[100, specs, "route_disagreement"];
  If[! StabilizerStateQ[PauliStabilizer[100]["ApplyCircuit", specs]],
    fail["not_stabilizer"]]];
(* Exactly solvable checks: the fresh register (constructor path) and a GHZ circuit
   through the COMPILED path itself; route agreement alone would pass a bug shared by
   both routes. *)
If[PauliStabilizer[6]["Stabilizers"] =!= Table[
     StringJoin@ReplacePart[ConstantArray["I", 6], i -> "Z"], {i, 6}],
  fail["closed_form_Z_generators"]];
ghzGate[];
WriteString["stdout", "CHECK|routes_agree_and_state_valid|True\n"];

Scan[
  Function[n,
    With[{ops = stream[n]},
      With[{specs = toSpec[ops], m = Length[ops]},
        (* toSpec has validated every op in THIS stream; additionally assert the compiled
           encoder accepts it, since ApplyCircuit silently falls back to the interpreted
           per-gate fold on ANY encoder rejection (the Fold branch of ApplyCircuit in
           Compiled.m) and the row would then time the wrong route under a compiled label. *)
        If[! ListQ[encode[specs]],
          fail["encoder_rejected_stream_n" <> ToString[n]]];
        (* Per-size gates on the FULL timed stream: StabilizerStateQ alone is a
           shape check that accepts a bit-flipped tableau, so the bulk result must
           also match the independent per-gate fold bit-exactly at every timed n,
           and conserve the symplectic form at n=100 (the Gram Dot is O(n^3),
           seconds at n=1000). *)
        With[{res = PauliStabilizer[n]["ApplyCircuit", specs]},
          If[! StabilizerStateQ[res],
            fail["timed_object_not_stabilizer_n" <> ToString[n]]];
          If[n === 100, symplecticGate[res, 100]]];
        routeGate[n, specs, "route_disagreement_n" <> ToString[n]];
        With[{ac = tm[PauliStabilizer[n]["ApplyCircuit", specs]["Signs"]]},
          pr["gates_n" <> ToString[n] <> "_AC_ms", ac];
          prRate["gates_n" <> ToString[n] <> "_AC", ac, m]];
        With[{pg = tm[perGateFold[PauliStabilizer[n], specs]["Signs"]]},
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
        If[keys =!= 2, fail["meas_branch_not_random_n" <> ToString[n]]]];
      pr["measZ_n" <> ToString[n] <> "_ms", tm[ps["M", 1]]];
      pr["measPauli_n" <> ToString[n] <> "_ms", tm[ps["M", zstr]]];
      pr["stabs_n" <> ToString[n] <> "_ms", tm[ps["Stabilizers"]]]]],
  $measSizes];

With[{a = mk[$composeSize], b = mk[$composeSize]},
  pr["compose_n" <> ToString[$composeSize] <> "_ms", tm[a[b]]]];
WriteString["stdout", "DONE\n"];
