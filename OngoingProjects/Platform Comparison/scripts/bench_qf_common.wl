(* Shared helpers for the two QF stabilizer benchmark harnesses in this
   directory (bench_phases_qf.wls, qf_final_all.wl). Get this file after
   setting $here (the scripts' own directory): stream[] resolves the op
   streams relative to it.

   Estimator: tm[] is min over 7 evaluations with ClearSystemCache[] before
   each. The rep budget is pinned to 7 across ALL engines in the mirror set
   (bench_phases_stim.py, bench_phases_qc.jl): the min estimator is
   downward-biased by sample count, so asymmetric budgets would bias sub-ms
   cross-engine ratios for reasons unrelated to the engines. Change it in
   every script of the set together or not at all. *)

fail[tag_] := (WriteString["stdout", "FAIL|", tag, "\n"]; Exit[1]);

SetAttributes[tm, HoldFirst];
tm[e_, reps_ : 7] := 1000. Min@Table[ClearSystemCache[]; First@AbsoluteTiming[e], {reps}];

stream[n_] := Import[FileNameJoin[{$here, "stab_ops_" <> ToString[n] <> ".json"}]];

(* Symplectic Gram of the tableau's two generator-row blocks (destabilizers,
   stabilizers), each row a 2n-wide symplectic vector: Clifford conjugation
   preserves all commutation relations, so this is conserved exactly; a shape
   check like StabilizerStateQ cannot see a violation (it accepts a bit-flipped
   tableau). *)
gram[t_] := Mod[Transpose[t[[1]]].t[[2]] + Transpose[t[[2]]].t[[1]], 2];

(* Per-gate property route; two-qubit targets are splatted: ps["CNOT", a, b],
   not ps["CNOT", {a, b}]. *)
perGateFold[ps_, specs_] :=
  Fold[#1[First[#2], Sequence @@ Flatten[{Last[#2]}]] &, ps, specs];

(* GHZ closed form through the compiled path: |000> under H1, CNOT12, CNOT13 is
   stabilized by XXX, ZZI, ZIZ (generator i is the image of Z_i, so the order
   is fixed, all signs +). *)
ghzGate[] := If[PauliStabilizer[3]["ApplyCircuit",
     {"H" -> 1, "CNOT" -> {1, 2}, "CNOT" -> {1, 3}}]["Stabilizers"] =!=
   {"XXX", "ZZI", "ZIZ"}, fail["ghz_closed_form"]];

(* Symplectic-form conservation: Clifford conjugation preserves all commutation
   relations, so the evolved Gram matrix must equal the fresh register's
   exactly. The Gram Dot is O(n^3): cheap at n=100, seconds at n=1000. *)
symplecticGate[evolved_, n_] :=
  If[gram[evolved["Tableau"]] =!= gram[PauliStabilizer[n]["Tableau"]],
    fail["symplectic_form_broken"]];

(* Bulk-vs-per-gate route agreement, bit-exact on tableau and signs: the two
   routes are independent implementations of the same evolution, so a benchmark
   that fails this would time a wrong computation. *)
routeGate[n_, specs_, tag_] :=
  With[{blk = PauliStabilizer[n]["ApplyCircuit", specs],
        fld = perGateFold[PauliStabilizer[n], specs]},
    If[blk["Tableau"] =!= fld["Tableau"] || blk["Signs"] =!= fld["Signs"],
      fail[tag]]];
