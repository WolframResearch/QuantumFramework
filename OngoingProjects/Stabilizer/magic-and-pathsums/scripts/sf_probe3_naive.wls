#!/usr/bin/env wolframscript
(* Probe 3: "naive, no rank reduction" demonstrations + P[theta] generic-angle blow-up.
   (a) P[theta] (generic numeric angle) blows up as 2^k, same as T.
   (b1) k T's on qubit 1 of |0> (T fixes |0>): true state stays |0> (rank 1), yet QF
        keeps 2^k components, ALL the identical tableau (distinct tableaux = 1).
   (b2) T^k on |+> with period 8 (T^8 = I): true state returns to |+> at k=8, yet 256
        components.
   No Quiet, no Print; final expression is the result Association.
   View with:  wolframscript -code 'Print[Get["sf_probe3_naive.wls"]]' *)

PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"];

theta = 0.7;
n = 3;
backbone = Join[Table["H" -> i, {i, n}], Table["CNOT" -> {i, i + 1}, {i, n - 1}]];

countP[k_] := Module[{res = QuantumCircuitOperator[Join[backbone, Table["P"[theta] -> (Mod[j-1,n]+1), {j, k}]]][Method -> "Stabilizer"]},
    <|"k" -> k, "components" -> If[Head[res] === StabilizerFrame, res["Length"], 1], "2^k" -> 2^k|>];

(* distinctTableaux (the point of b1) needs no materialization; the optional state==|0>
   correctness check materializes the StateVector, so we only run it for small k (large-k
   materialization of a 2^k-component frame underflows the coefficient products). *)
zeroDemo[k_] := Module[{res, sv},
    res = QuantumCircuitOperator[Table["T" -> 1, {k}]][Method -> "Stabilizer"];  (* on |0> *)
    <|"k" -> k,
      "components" -> If[Head[res] === StabilizerFrame, res["Length"], 1],
      "distinctTableaux" -> If[Head[res] === StabilizerFrame, Length @ DeleteDuplicates[#["Tableau"] & /@ res["Stabilizers"]], 1],
      "state==|0>" -> If[k <= 4,
         sv = Normal @ If[Head[res] === StabilizerFrame, res["StateVector"], res["State"]["StateVector"]];
         Chop[N[Norm[sv - UnitVector[Length[sv], 1]]]] < 10.^-12,   (* exact sv; N the scalar only *)
         "(not materialized)"]|>];

(* Keep the StateVector EXACT (Normal, no N) and N only the final real overlap, so the
   2^k-term coefficient products do not underflow as machine numbers. *)
plusDemo[k_] := Module[{res, sv, plus = {1, 1}/Sqrt[2]},
    res = QuantumCircuitOperator[Join[{"H" -> 1}, Table["T" -> 1, {k}]]][Method -> "Stabilizer"];
    sv = Normal @ If[Head[res] === StabilizerFrame, res["StateVector"], res["State"]["StateVector"]];
    <|"k" -> k,
      "components" -> If[Head[res] === StabilizerFrame, res["Length"], 1],
      "|<+|psi>|^2" -> Round[N @ Abs[Conjugate[plus] . sv]^2, 10.^-6]|>];

<|
    "(a) P[theta] blow-up"           -> (countP /@ Range[0, 8]),
    "(b1) T^k on |0> (rank-1 truth)" -> (zeroDemo /@ Range[1, 8]),
    "(b2) T^k on |+> (period 8)"     -> (plusDemo /@ {1, 2, 4, 7, 8})
|>
