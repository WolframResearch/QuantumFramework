#!/usr/bin/env wolframscript
(* Verification demos for Magic-Exponent-Impossibility-Audit.md.                    *)
(*                                                                                  *)
(* Audited sentence: "None removes the exponent in magic, because universality      *)
(* forbids that."  The audit's central computational claim is that the exponent is  *)
(* governed by a STRUCTURAL parameter, not by the T-count / magic: there exist      *)
(* circuit families with UNBOUNDED T-count whose exact amplitude is poly-time.      *)
(* These demos exhibit such families and check them against brute force, and        *)
(* cross-check the QF kernel paying 2^t where the intrinsic cost is O(1).            *)
(*                                                                                  *)
(* A Clifford+T amplitude is, up to a power-of-sqrt2 normalization, the Gauss sum   *)
(*   gap(f) = sum_{x in F_2^n} (-1)^{f(x)} = #{f=0} - #{f=1} = 2 #{f=0} - 2^n,       *)
(* of a degree-<=3 polynomial f over F_2 (real / Toffoli-Hadamard fragment).        *)
(* No Quiet, no Print; each block ends with a labeled Echo of the asserted result.  *)

(* ============================================================================ *)
(* Brute-force F_2 primitives (the ground truth every claim is checked against). *)
(* ============================================================================ *)
gapBrute[f_, vars_] := With[{n = Length[vars]},
   Total[(-1)^Table[Mod[f /. Thread[vars -> a], 2], {a, Tuples[{0, 1}, n]}]]];
countZeroBrute[f_, vars_] := Count[Tuples[{0, 1}, Length[vars]],
   a_ /; Mod[f /. Thread[vars -> a], 2] === 0];

(* ---------------------------------------------------------------------------- *)
Echo["==== (1) the gap identity and the degree-2 / degree-3 line ====", ""];
SeedRandom[2026];
vs6 = Array[x, 6];
(* a random cubic over F_2 in 6 variables *)
randCubic = Mod[Total[RandomInteger[1, Length[#]] #] &[
     Times @@@ Join[{1}, vs6, Subsets[vs6, {2}], Subsets[vs6, {3}]]], 2];
Echo[gapBrute[randCubic, vs6] === 2 countZeroBrute[randCubic, vs6] - 2^6,
   "(1) gap(f) == 2 #{f=0} - 2^n on a random cubic? (expect True):"];
(* degree-2 (quadratic) gap is poly-time via an F_2 rank normal form; here we just *)
(* confirm a single quadratic AND term has the textbook gap 2 over its 2 vars.      *)
Echo[{gapBrute[x[1] x[2], {x[1], x[2]}], 2},
   "(1) gap(x1 x2) {brute, expected 2}  (degree 2 = Gottesman-Knill-easy):"];
Echo[{gapBrute[x[1] x[2] x[3], {x[1], x[2], x[3]}], 6},
   "(1) gap(x1 x2 x3) {brute, expected 6}  (one CCZ = degree 3):"];

(* ============================================================================ *)
(* (2) THE WITNESS: a STAR of cubic terms sharing one hub variable.              *)
(*   f_m(x0,...,x_{2m}) = sum_{i=1}^m  x0 * x_{2i-1} * x_{2i}                     *)
(* T-count proxy  = number of cubic (CCZ) terms = m, UNBOUNDED.                  *)
(* Hitting set of the cubic monomials = {x0}, size 1 (Montanaro 2^{|S|}).        *)
(* Conditioning on x0 splits the sum into two QUADRATIC Gauss sums:              *)
(*   x0=0 : f = 0           -> 2^{2m}                                            *)
(*   x0=1 : f = sum x_{2i-1}x_{2i} (m independent ANDs) -> prod 2 = 2^m          *)
(*   => gap(f_m) = 2^{2m} + 2^m, computed in O(m) time, magic unbounded.         *)
(* ============================================================================ *)
Echo["==== (2) star family: unbounded T-count, gap is poly-time exact ====", ""];
starF[m_] := With[{vars = Array[x, 2 m + 1, 0]},
   Total[Table[x[0] x[2 i - 1] x[2 i], {i, m}]]];
starVars[m_] := Array[x, 2 m + 1, 0];
starHittingSetSplit[m_] := 2^(2 m) + 2^m;   (* the O(m) hitting-set-1 formula *)
starTable = Table[
   {m, m (* cubic-term / CCZ count *), gapBrute[starF[m], starVars[m]],
    starHittingSetSplit[m]}, {m, 1, 5}];
Echo[Grid[Prepend[starTable, {"m", "T-count (CCZ)", "gap (brute)", "gap (hitting-set split)"}]],
   "(2) star f_m: {m, magic=m, brute gap, O(m) split} (cols 3,4 must agree):"];
Echo[AllTrue[starTable, #[[3]] === #[[4]] &],
   "(2) gap(brute) === 2^{2m}+2^m for m=1..5, magic unbounded? (expect True):"];

(* ============================================================================ *)
(* (3) SECOND WITNESS: a MATCHING of disjoint cubic terms (no shared variable). *)
(*   g_m(x1,...,x_{3m}) = sum_{i=1}^m x_{3i-2} x_{3i-1} x_{3i}                   *)
(* Hitting set = m (no shared var), BUT the path-variable interaction graph is   *)
(* m disjoint triangles: bounded tree-width / rank-width. The gap FACTORIZES:    *)
(*   gap(g_m) = prod_i gap(one CCZ term) = 6^m, poly-time, magic = m unbounded.  *)
(* This is the de Colnet rank-width axis: small width, T-count irrelevant.       *)
(* ============================================================================ *)
Echo["==== (3) matching family: bounded width, gap factorizes, magic unbounded ====", ""];
matchG[m_] := With[{},
   Total[Table[x[3 i - 2] x[3 i - 1] x[3 i], {i, m}]]];
matchVars[m_] := Array[x, 3 m];
matchTable = Table[
   {m, m, gapBrute[matchG[m], matchVars[m]], 6^m}, {m, 1, 3}];
Echo[Grid[Prepend[matchTable, {"m", "T-count (CCZ)", "gap (brute)", "6^m (factorized)"}]],
   "(3) matching g_m: {m, magic=m, brute gap, 6^m} (cols 3,4 must agree):"];
Echo[AllTrue[matchTable, #[[3]] === #[[4]] &],
   "(3) gap(brute) === 6^m for m=1..3, magic unbounded? (expect True):"];

(* ============================================================================ *)
(* (4) LIVE KERNEL: QF's StabilizerFrame pays 2^t where the intrinsic cost is    *)
(* O(1).  T^t|+> = (|0> + e^{i t pi/4}|1>)/sqrt2 is ONE eighth-root phase: its    *)
(* exact amplitude is a closed form computed in O(1), independent of t (magic).   *)
(* QF nonetheless stores 2^t frame components for it.  So the exponent QF carries *)
(* (2^t, "in magic") is NOT the intrinsic strong-simulation cost of the family.   *)
(* ============================================================================ *)
Echo["==== (4) live QF cross-check: 2^t frame components for an O(1) amplitude ====", ""];
PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"];
stabRank[x_] := If[Head[x] === StabilizerFrame, x["Length"], 1];
TtPlus[t_] := PauliStabilizer[1][Join[{"H"}, Table["T", t]]];
trueTtPlus[t_] := {1, Exp[I t Pi/4]}/Sqrt[2];        (* the O(1) closed form *)
denseTtPlus[t_] := Normal[
   QuantumCircuitOperator[Join[{"H" -> 1}, Table["T" -> 1, t]]][
     QuantumState["0"], Method -> "Schrodinger"]["StateVector"]];
qfTable = Table[{t, stabRank[TtPlus[t]], 2^t,
    QuantumSimilarity[TtPlus[t]["State"],
     QuantumState[denseTtPlus[t]]]}, {t, 1, 6}];
Echo[Grid[Prepend[qfTable, {"t (T-count)", "frame components", "2^t", "similarity to dense"}]],
   "(4) T^t|+>: {t, QF components, 2^t, fidelity} (components=2^t, fidelity=1):"];
Echo[AllTrue[qfTable, #[[2]] === #[[3]] &],
   "(4) QF frame component count === 2^t (the blow-up 'in magic')? (expect True):"];
Echo[AllTrue[Range[6], Chop[Abs[trueTtPlus[#] . Conjugate[denseTtPlus[#]]] - 1] === 0 &],
   "(4) true O(1) closed form {1,e^{i t pi/4}}/sqrt2 == dense up to phase? (expect True):"];

(* ============================================================================ *)
(* (5) THE WALL IS REAL (no #P break, no tractable-class enlargement). The same  *)
(* gap, on a GENERIC cubic with no small hitting set and no bounded width, is the *)
(* #P-hard object; and the eighth-root phase cannot be hidden in a stabilizer     *)
(* sign (the involution (pm P)^2 = I breaks), so no Gottesman-Knill extension.    *)
(* ============================================================================ *)
Echo["==== (5) the wall: characteristic-2 obstruction (the impossibility that IS real) ====", ""];
Echo[FullSimplify[(Exp[I Pi/4] PauliMatrix[1]) . (Exp[I Pi/4] PauliMatrix[1])] === I IdentityMatrix[2],
   "(5) (e^{i pi/4} X)^2 = i I != I (no +1 eigenspace; stabilizes nothing)? (expect True):"];
Echo[Factor[Cyclotomic[8, x], Modulus -> 2],
   "(5) 1+x^4 over F_2 (expect (1+x)^4: omega's cyclotomic structure collapses mod 2):"];

Echo["ALL DEMOS COMPLETE", "status"];
