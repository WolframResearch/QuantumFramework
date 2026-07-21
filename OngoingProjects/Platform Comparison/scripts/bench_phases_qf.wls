(* Phase-resolved QF stabilizer benchmark. Companion to bench_phases_stim.py and
   bench_phases_qc.jl; the three scripts share one output schema

       RESULT|qf|<phase>|n=<n>|k=<k>|<ms>

   so cross-engine ratios can be composed phase by phase instead of quoted as one
   number whose boundaries are implicit.

   PHASES (each timed separately, none nested inside another):
     ingest      rule-form circuit -> packed integer gate matrix
                 (encodeStabilizerGates; WL Replace/Dispatch pass, n-independent)
     ctor        PauliStabilizer[n] register construction
     pack        canonical tableau + signs -> generator-major packed words
     simulate    compiledCliffordFold on pre-packed state and pre-encoded gates:
                 the ONLY phase that is tableau algebra
     materialize packed kernel output -> canonical {2,n,2n} arrays + object assembly
     e2e         ctor + ps["ApplyCircuit", specs]["Signs"]: what a user pays
     residual    DERIVED, not separately timed: e2e minus the sum of the measured
                 phases (property dispatch, validation, assembly)

   The neutral JSON op list is parsed OUTSIDE all timed regions on every engine;
   ingest starts from each engine's native in-memory circuit description (here the
   rule list, for Stim the Python objects or circuit text).

   Gate catalog: Keys[$stabilizerGateCodes] ONLY (H, S, S-dagger, CNOT, SWAP, X, Y, Z).
   $stabilizerCliffordNames also lists V and T, but the encoder rejects those and
   ApplyCircuit then falls back SILENTLY to the interpreted per-gate fold
   (Compiled.m:244), so a stream containing them would time the wrong engine under
   this label. The harness asserts the encoder accepted every stream it times.

   Validity gates run before any timing:
     - GHZ closed form through the compiled path itself,
     - symplectic-form (Gram) conservation, which a pure shape check like
       StabilizerStateQ cannot see (it accepts a bit-flipped tableau),
     - bulk-vs-per-gate route agreement on tableau and signs.
   Edge probes run after: the n=1024 packed-tableau regression guard (before kernel
   commit 253afdfb the fresh tableau's member matrices went SparseArray above
   n=1000, packGenRows fed the kernel an unpackable argument, cfta fired, and the
   fold silently ran interpreted at 163-459 us/gate, ~10^3 x the compiled kernel;
   the probe must now report compiled_path_available|True) and the V-fallback.

   Run: wolframscript -file bench_phases_qf.wls   (exit 1 on any failed gate) *)

$here = DirectoryName[$InputFileName];
PacletDirectoryLoad[FileNameDrop[$here, -3]];
Needs["Wolfram`QuantumFramework`"];

(* Kernel internals, grounded in Kernel/Stabilizer/Compiled.m (worktree copy):
   encodeStabilizerGates :196, packGenRows :79, unpackGenRows :88,
   stabilizerValidMasks :97, compiledCliffordFold :111, $bitsPerWord = 62 :60. *)
ctx = "Wolfram`QuantumFramework`PackageScope`";
encode = Symbol[ctx <> "encodeStabilizerGates"];
pack = Symbol[ctx <> "packGenRows"];
unpack = Symbol[ctx <> "unpackGenRows"];
masks = Symbol[ctx <> "stabilizerValidMasks"];
fold = Symbol[ctx <> "compiledCliffordFold"];
codes = Symbol[ctx <> "$stabilizerGateCodes"];
bits = Symbol[ctx <> "$bitsPerWord"];

fail[tag_] := (WriteString["stdout", "FAIL|", tag, "\n"]; Exit[1]);

pr[engine_, phase_, n_, k_, ms_] := WriteString["stdout",
  "RESULT|", engine, "|", phase, "|n=", n, "|k=", k, "|", Round[ms, 0.001], "\n"];

SetAttributes[tm, HoldFirst];
tm[e_, reps_ : 5] := 1000. Min@Table[ClearSystemCache[]; First@AbsoluteTiming[e], {reps}];

toSpec[ops_] := Map[Replace[{
    {g_String, q_Integer} :> g -> q + 1,
    {g_String, a_Integer, b_Integer} :> g -> {a + 1, b + 1},
    bad_ :> fail["bad_gate_spec: " <> ToString[bad]]}], ops];

stream[n_] := Import[FileNameJoin[{$here, "stab_ops_" <> ToString[n] <> ".json"}]];

gram[t_] := Mod[Transpose[t[[1]]].t[[2]] + Transpose[t[[2]]].t[[1]], 2];

PauliStabilizer[4]["ApplyCircuit", {"H" -> 1}];  (* build the compiled kernel once *)

(* ---- validity gates, all on the benchmarked (compiled) path ---- *)

(* GHZ closed form: |000> under H1, CNOT12, CNOT13 is stabilized by XXX, ZZI, ZIZ
   (generator i is the image of Z_i, so the order is fixed, all signs +). *)
If[PauliStabilizer[3]["ApplyCircuit",
     {"H" -> 1, "CNOT" -> {1, 2}, "CNOT" -> {1, 3}}]["Stabilizers"] =!=
   {"XXX", "ZZI", "ZIZ"}, fail["ghz_closed_form"]];

(* Symplectic form conservation at n=100 (the Gram Dot is O(n^3), cheap here,
   seconds at n=1000): Clifford conjugation preserves all commutation relations,
   so the evolved Gram matrix must equal the fresh register's exactly. *)
With[{specs = toSpec@stream[100]},
  With[{evolved = PauliStabilizer[100]["ApplyCircuit", specs]},
    If[gram[evolved["Tableau"]] =!= gram[PauliStabilizer[100]["Tableau"]],
      fail["symplectic_form_broken"]];
    (* route agreement: property route is a genuinely different code path *)
    With[{fld = Fold[#1[First[#2], Sequence @@ Flatten[{Last[#2]}]] &,
        PauliStabilizer[100], Take[specs, 500]],
        blk = PauliStabilizer[100]["ApplyCircuit", Take[specs, 500]]},
      If[blk["Tableau"] =!= fld["Tableau"] || blk["Signs"] =!= fld["Signs"],
        fail["route_disagreement"]]]]];
WriteString["stdout", "CHECK|qf|validity_gates_passed|True\n"];

(* ---- phase timings ---- *)
Scan[Function[n,
   With[{ops = stream[n]},
    With[{specs = toSpec[ops], k = Length[ops], nw = Ceiling[2 n/bits]},
     With[{enc = encode[specs]},
      If[! ListQ[enc], fail["encoder_rejected_stream_n" <> ToString[n]]];
      With[{reg = PauliStabilizer[n]},
       With[{t = reg["Tableau"], phase = Developer`ToPackedArray[(1 - reg["Signs"])/2],
             vm = masks[n, nw]},
        If[Head[t] === SparseArray || ! Developer`PackedArrayQ[t[[1]]],
          fail["tableau_not_packed_n" <> ToString[n]]];
        With[{x0 = pack[t[[1]], nw], z0 = pack[t[[2]], nw],
              s0 = First@pack[{phase}, nw]},
         With[{res = fold[x0, z0, s0, vm, enc]},
          (* the simulate phase must feed the compiled kernel packed input *)
          If[! Developer`PackedArrayQ[x0] || ! Developer`PackedArrayQ[enc],
            fail["kernel_input_unpacked_n" <> ToString[n]]];
          (* the hand-reconstructed pack -> fold -> unpack pipeline must equal the
             public ApplyCircuit path bit-exactly, or the simulate / materialize
             rows time a phantom (e.g. if the kernel's packing convention moves
             under these PackageScope signatures) *)
          With[{manual = PauliStabilizer[<|
               "Signs" -> 1 - 2 First@unpack[{res[[-1]]}, 2 n],
               "Tableau" -> {unpack[res[[1 ;; n]], 2 n],
                 unpack[res[[n + 1 ;; 2 n]], 2 n]}|>],
              pub = PauliStabilizer[n]["ApplyCircuit", specs]},
           If[manual["Tableau"] =!= pub["Tableau"] || manual["Signs"] =!= pub["Signs"],
             fail["reconstruction_drift_n" <> ToString[n]]]];
          With[{tIngest = tm[encode[specs]],
                tCtor = tm[PauliStabilizer[n]],
                tPack = tm[{pack[t[[1]], nw], pack[t[[2]], nw], First@pack[{phase}, nw]}],
                tSim = tm[fold[x0, z0, s0, vm, enc]],
                tMat = tm[PauliStabilizer[<|
                  "Signs" -> 1 - 2 First@unpack[{res[[-1]]}, 2 n],
                  "Tableau" -> {unpack[res[[1 ;; n]], 2 n],
                    unpack[res[[n + 1 ;; 2 n]], 2 n]}|>]],
                tE2e = tm[PauliStabilizer[n]["ApplyCircuit", specs]["Signs"]]},
           pr["qf", "ingest", n, k, tIngest];
           pr["qf", "ctor", n, k, tCtor];
           pr["qf", "pack", n, k, tPack];
           pr["qf", "simulate", n, k, tSim];
           pr["qf", "materialize", n, k, tMat];
           pr["qf", "e2e", n, k, tE2e];
           pr["qf", "residual", n, k, tE2e - (tIngest + tCtor + tPack + tSim + tMat)]];
          (* single Z measurement on the evolved state, for the workload Stim is
             actually built around; QF returns BOTH outcome branches as states,
             which is more work than a sampling collapse: see the caveat printed.
             The 2-key assert keeps this row on the random-outcome branch: a
             shallower stream would silently switch algorithm across n. *)
          With[{ev = PauliStabilizer[n]["ApplyCircuit", specs]},
           With[{keys = Length@Keys@ev["M", 1]},
            WriteString["stdout", "CHECK|qf|measure_branches|n=", n, "|", keys, "\n"];
            If[keys =!= 2, fail["meas_branch_not_random_n" <> ToString[n]]]];
           pr["qf", "measure1", n, k, tm[ev["M", 1]]]];
          ]]]]]]]],
  {100, 500, 1000}];
WriteString["stdout",
  "CAVEAT|qf|measure1 returns an association of BOTH outcome states, not one sampled collapse\n"];

(* ---- edge probes ---- *)
(* Regression guard for the n>1000 sparse-tableau fallback fixed at 253afdfb: the
   fresh register's X/Z tableau MATRICES (members of the length-2 Tableau list, not
   its head) must be packed dense, or packGenRows feeds the compiled kernel a
   sparse word matrix, cfta fires, and the fold runs interpreted at ~10^3 x the
   kernel cost. Probe the member matrix; probing Head[t] alone passes wrongly. *)
With[{t = PauliStabilizer[1024]["Tableau"]},
  WriteString["stdout", "EDGE|qf|n1024_xmatrix_head|", ToString[Head[t[[1]]]],
    "|compiled_path_available|",
    ToString[Developer`PackedArrayQ[Developer`ToPackedArray[t[[1]]]]], "\n"]];
(* V is in $stabilizerCliffordNames but not in $stabilizerGateCodes: the encoder
   must reject it (else ApplyCircuit silently times the interpreted route). *)
If[ListQ[encode[{"V" -> 1}]], fail["encoder_unexpectedly_accepts_V"]];
WriteString["stdout", "EDGE|qf|V_rejected_by_encoder|True\n"];
WriteString["stdout", "DONE|qf|", ToString[codes // Length], " encoded gate codes\n"];
