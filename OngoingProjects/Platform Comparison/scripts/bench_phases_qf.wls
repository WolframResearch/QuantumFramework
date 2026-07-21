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
                 phases (property dispatch, validation, assembly); e2e times the
                 public path while the phases time the manual kernel pipeline,
                 each min-over-reps, so the difference upper-bounds the
                 unattributed overhead and can go negative under noise

   The neutral JSON op list is parsed OUTSIDE all timed regions on every engine;
   ingest starts from each engine's native in-memory circuit description (here the
   rule list, for Stim the Python objects or circuit text).

   Estimator: every phase row is min over 7 evaluations (tm in
   bench_qf_common.wl), a budget pinned across all three engines: the min
   estimator is downward-biased by sample count, so asymmetric rep budgets
   would bias sub-ms cross-engine ratios for reasons unrelated to the engines.

   Gate catalog: Keys[$stabilizerGateCodes] ONLY (H, S, S-dagger, CNOT, SWAP, X, Y, Z).
   $stabilizerCliffordNames also lists V and T, but the encoder rejects those and
   ApplyCircuit then falls back SILENTLY to the interpreted per-gate fold
   (the Fold branch of ApplyCircuit in Compiled.m), so a stream containing them
   would time the wrong engine under this label. The harness asserts the encoder
   accepted every stream it times.

   Validity gates run before any timing:
     - GHZ closed form through the compiled path itself,
     - symplectic-form (Gram) conservation, which a pure shape check like
       StabilizerStateQ cannot see (it accepts a bit-flipped tableau),
     - bulk-vs-per-gate route agreement on tableau and signs,
     - the full encoded catalog (all eight gate codes) compiled-vs-per-gate,
       since the timed H/S/CNOT streams never exercise the other five codes.
   Asymptotic gates run after the sweep: simulate must track the gate count k,
   gated by a three-point log-log slope along the k=20n diagonal and at fixed
   n=1000 via an off-diagonal k=2000 fold that isolates the k-exponent (fixed
   per-call overhead keeps it mildly sublinear; the interpreted fallback breaks
   the diagonal by ~10^2); materialize must grow ~n^2 (three-point log-log
   slope in n; it unpacks the full {2, n, 2n} tableau).
   Edge probes run after and gate the exit code: the n=1024 packed-tableau
   regression guard (before kernel
   commit 253afdfb the fresh tableau's member matrices went SparseArray above
   n=1000, packGenRows fed the kernel an unpackable argument, cfta fired, and the
   fold silently ran interpreted at 163-459 us/gate, ~10^3 x the compiled kernel;
   the run now fails unless compiled_path_available|True) and the V-fallback.

   Run: wolframscript -file bench_phases_qf.wls   (exit 1 on any failed gate) *)

$here = DirectoryName[$InputFileName];
PacletDirectoryLoad[FileNameDrop[$here, -3]];
Needs["Wolfram`QuantumFramework`"];
Get[FileNameJoin[{$here, "bench_qf_common.wl"}]];  (* fail, tm (min over 7), stream, gram, perGateFold *)

(* Kernel internals, grounded in Kernel/Stabilizer/Compiled.m (worktree copy):
   encodeStabilizerGates, packGenRows, unpackGenRows, stabilizerValidMasks,
   compiledCliffordFold, $bitsPerWord = 62. *)
ctx = "Wolfram`QuantumFramework`PackageScope`";
encode = Symbol[ctx <> "encodeStabilizerGates"];
pack = Symbol[ctx <> "packGenRows"];
unpack = Symbol[ctx <> "unpackGenRows"];
masks = Symbol[ctx <> "stabilizerValidMasks"];
fold = Symbol[ctx <> "compiledCliffordFold"];
codes = Symbol[ctx <> "$stabilizerGateCodes"];
bits = Symbol[ctx <> "$bitsPerWord"];

pr[engine_, phase_, n_, k_, ms_] := WriteString["stdout",
  "RESULT|", engine, "|", phase, "|n=", n, "|k=", k, "|", Round[ms, 0.001], "\n"];

toSpec[ops_] := Map[Replace[{
    {g_String, q_Integer} :> g -> q + 1,
    {g_String, a_Integer, b_Integer} :> g -> {a + 1, b + 1},
    bad_ :> fail["bad_gate_spec: " <> ToString[bad]]}], ops];

(* One definition each for the pack and materialize pipelines: the
   reconstruction gate certifies exactly the expressions the timed rows
   measure, so both must reference a single source. *)
packTriple[t_, phase_, nw_] :=
  {pack[t[[1]], nw], pack[t[[2]], nw], First@pack[{phase}, nw]};
materialize[res_, n_] := PauliStabilizer[<|
   "Signs" -> 1 - 2 First@unpack[{res[[-1]]}, 2 n],
   "Tableau" -> {unpack[res[[1 ;; n]], 2 n], unpack[res[[n + 1 ;; 2 n]], 2 n]}|>];

PauliStabilizer[4]["ApplyCircuit", {"H" -> 1}];  (* build the compiled kernel once *)

(* ---- validity gates (ghzGate/symplecticGate/routeGate from
        bench_qf_common.wl), all on the benchmarked (compiled) path ---- *)

ghzGate[];
With[{specs = toSpec@stream[100]},
  symplecticGate[PauliStabilizer[100]["ApplyCircuit", specs], 100];
  routeGate[100, Take[specs, 500], "route_disagreement"]];

(* Full encoded catalog: the timed streams carry only the cross-engine H/S/CNOT
   contract, so the remaining codes (S-dagger, SWAP, X, Y, Z) and the encoder's
   CZ (H.CNOT.H expansion) and CX (CNOT alias) branches are exercised here,
   compiled fold vs per-gate route, bit-exact on tableau and signs, on a state
   entangled enough that the X/Y/Z sign updates are not vacuous. *)
With[{cat = {"H" -> 1, "CNOT" -> {1, 4}, "S" -> 2, SuperDagger["S"] -> 2,
     "SWAP" -> {2, 5}, "X" -> 3, "Y" -> 4, "Z" -> 5, "CZ" -> {4, 5},
     "CX" -> {2, 6}, "CNOT" -> {3, 6}, "H" -> 6}},
  If[! ListQ[encode[cat]] || Union[encode[cat][[All, 1]]] =!= Sort[Values[codes]],
    fail["catalog_not_fully_encoded"]];
  With[{blk = PauliStabilizer[6]["ApplyCircuit", cat],
        fld = perGateFold[PauliStabilizer[6], cat]},
    If[blk["Tableau"] =!= fld["Tableau"] || blk["Signs"] =!= fld["Signs"] ||
       blk["Signs"] === PauliStabilizer[6]["Signs"],
      fail["catalog_route_disagreement"]]]];
WriteString["stdout", "CHECK|qf|validity_gates_passed|True\n"];

(* ---- phase timings ---- *)

(* One timed size: gates first (packed kernel input, reconstruction drift),
   then the phase rows, then the measurement row. Each With level is a genuine
   data dependency; the sown timings feed the asymptotic gates below. *)
benchOne[n_] := With[{ops = stream[n]},
  With[{specs = toSpec[ops], k = Length[ops], nw = Ceiling[2 n/bits]},
   With[{enc = encode[specs], reg = PauliStabilizer[n]},
    If[! ListQ[enc], fail["encoder_rejected_stream_n" <> ToString[n]]];
    With[{t = reg["Tableau"], phase = Developer`ToPackedArray[(1 - reg["Signs"])/2],
          vm = masks[n, nw]},
     If[Head[t] === SparseArray || ! Developer`PackedArrayQ[t[[1]]],
       fail["tableau_not_packed_n" <> ToString[n]]];
     With[{p0 = packTriple[t, phase, nw]},
      With[{res = fold[Sequence @@ p0, vm, enc]},
       (* the simulate phase must feed the compiled kernel packed input *)
       If[! Developer`PackedArrayQ[First[p0]] || ! Developer`PackedArrayQ[enc],
         fail["kernel_input_unpacked_n" <> ToString[n]]];
       (* the hand-reconstructed pack -> fold -> unpack pipeline must equal the
          public ApplyCircuit path bit-exactly, or the simulate / materialize
          rows time a phantom (e.g. if the kernel's packing convention moves
          under these PackageScope signatures) *)
       With[{manual = materialize[res, n],
             pub = PauliStabilizer[n]["ApplyCircuit", specs]},
        If[manual["Tableau"] =!= pub["Tableau"] || manual["Signs"] =!= pub["Signs"],
          fail["reconstruction_drift_n" <> ToString[n]]]];
       With[{tIngest = tm[encode[specs]],
             tCtor = tm[PauliStabilizer[n]],
             tPack = tm[packTriple[t, phase, nw]],
             tSim = tm[fold[Sequence @@ p0, vm, enc]],
             tMat = tm[materialize[res, n]],
             tE2e = tm[PauliStabilizer[n]["ApplyCircuit", specs]["Signs"]]},
        MapThread[pr["qf", #1, n, k, #2] &,
          {{"ingest", "ctor", "pack", "simulate", "materialize", "e2e"},
           {tIngest, tCtor, tPack, tSim, tMat, tE2e}}];
        pr["qf", "residual", n, k, tE2e - (tIngest + tCtor + tPack + tSim + tMat)];
        Sow[<|"n" -> n, "k" -> k, "simulate" -> tSim, "materialize" -> tMat,
          "p0" -> p0, "vm" -> vm, "enc" -> enc|>]];
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
       ]]]]]];

sweep = First@Last@Reap@Scan[benchOne, {100, 500, 1000}];
WriteString["stdout",
  "CAVEAT|qf|measure1 returns an association of BOTH outcome states, not one sampled collapse\n"];

WriteString["stdout",
  "CAVEAT|qf|residual = e2e minus the sum of min-timed phases; e2e times the public path and the phases the manual kernel pipeline, so residual upper-bounds unattributed overhead and can go negative under noise\n"];

(* ---- asymptotic gates across the sweep ---- *)
(* simulate is one compiled pass over the gate list, so its cost must track the
   gate count k within noise (the interpreted fallback breaks this by ~10^2);
   materialize unpacks and validates the full {2, n, 2n} tableau, so it must
   grow ~quadratically in n. The sweep pins k = 20 n, so the diagonal slope
   alone cannot separate k-linearity from n-dependence: the diagonal is gated
   by a log-log slope fitted over all three sweep points, and an off-diagonal
   fold at fixed n = 1000 over the first k = 2000 encoded gates isolates the
   k-exponent (mildly sublinear, from fixed per-call overhead). materialize is
   gated the same way: a three-point log-log slope in n against its ~n^2 law.
   The bands absorb machine drift while still catching a phase that silently
   changes complexity class. *)
With[{lo = First[sweep], hi = Last[sweep]},
  With[{tKLo = tm[fold[Sequence @@ hi["p0"], hi["vm"], Take[hi["enc"], lo["k"]]]]},
    pr["qf", "simulate", hi["n"], lo["k"], tKLo];
    With[{slopeIn = Function[{xs, ys}, Coefficient[
        Fit[Transpose[{Log[xs], Log[ys]}], {1, \[FormalX]}, \[FormalX]],
        \[FormalX]]]},
     With[{simSlope = slopeIn[#["k"] & /@ sweep, #["simulate"] & /@ sweep],
           matSlope = slopeIn[#["n"] & /@ sweep, #["materialize"] & /@ sweep],
           kR = (hi["simulate"]/tKLo)/(hi["k"]/lo["k"])},
      If[! (1/2 <= simSlope <= 3/2), fail["simulate_not_linear_in_gate_count"]];
      If[! (1/3 <= kR <= 3), fail["simulate_k_exponent_off_at_fixed_n"]];
      If[! (3/2 <= matSlope <= 5/2), fail["materialize_not_quadratic_in_n"]];
      WriteString["stdout", "CHECK|qf|simulate_scaling_linear_in_k|True\n"];
      WriteString["stdout", "CHECK|qf|simulate_k_exponent_fixed_n|True\n"];
      WriteString["stdout", "CHECK|qf|materialize_scaling_quadratic_in_n|True\n"]]]]];

(* ---- edge probes ---- *)
(* Regression guard for the n>1000 sparse-tableau fallback fixed at 253afdfb: the
   fresh register's X/Z tableau MATRICES (members of the length-2 Tableau list, not
   its head) must be packed dense, or packGenRows feeds the compiled kernel a
   sparse word matrix, cfta fires, and the fold runs interpreted at ~10^3 x the
   kernel cost. Probe the member matrix; probing Head[t] alone passes wrongly. *)
With[{t = PauliStabilizer[1024]["Tableau"]},
  If[! Developer`PackedArrayQ[Developer`ToPackedArray[t[[1]]]],
    fail["n1024_sparse_fallback"]];
  WriteString["stdout", "EDGE|qf|n1024_xmatrix_head|", ToString[Head[t[[1]]]],
    "|compiled_path_available|True\n"]];
(* V is in $stabilizerCliffordNames but not in $stabilizerGateCodes: the encoder
   must reject it (else ApplyCircuit silently times the interpreted route). *)
If[ListQ[encode[{"V" -> 1}]], fail["encoder_unexpectedly_accepts_V"]];
WriteString["stdout", "EDGE|qf|V_rejected_by_encoder|True\n"];
WriteString["stdout", "DONE|qf|", ToString[codes // Length], " encoded gate codes\n"];
