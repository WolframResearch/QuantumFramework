(* Experiment 2: where the algebra DOES transfer (diagonal Clifford+T as a   *)
(* Z_8 phase polynomial, rank 1 in principle) and where it CANNOT (H between  *)
(* T's: 2^k branching, exponential relocated into a Gauss sum).               *)
(* Run: wolframscript -file exp2-phasepoly-and-blowup.wls                      *)

PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"];

w = Exp[I Pi/4];  (* primitive 8th root, omega = e^{2 pi i / 8} *)
vEqual[a_, b_] := Max[Abs[N[a - b]]] < 10^-12;  (* numeric equality of complex vectors *)

(* ======================================================================== *)
(* D. Hypothesis 3: a DIAGONAL Clifford+T circuit is a phase polynomial.     *)
(*    U = T_1 . T_2 . CZ_{12} on |++>.                                       *)
(*    amplitude(x1,x2) = (1/2) w^{phi(x)}, phi = x1 + x2 + 4 x1 x2 (mod 8).  *)
(*    qubit 1 is most significant, so basis index 1+2 x1 + x2.              *)
(* ======================================================================== *)

denseDiag = QuantumCircuitOperator[{"H" -> 1, "H" -> 2, "T" -> 1, "T" -> 2, "CZ" -> {1, 2}}][
    QuantumState["00"], Method -> "Schrodinger"];
vDense = Chop[denseDiag["StateVector"] // Normal];

phi[x1_, x2_] := Mod[x1 + x2 + 4 x1 x2, 8];
vPoly = Flatten @ Table[(1/2) w^phi[x1, x2], {x1, 0, 1}, {x2, 0, 1}];
Echo[vEqual[vDense, vPoly], "D1 phase-polynomial reproduces dense T1.T2.CZ|++>:"];
Echo[Simplify[vPoly], "D1 the rank-1 phase-poly statevector (one support + one Z_8 poly):"];

(* QF's StabilizerFrame for the SAME circuit does NOT compress: it branches  *)
(* once per T gate, to 2^(#T) components, although the true object is rank 1 *)
(* as a phase polynomial.                                                    *)
qfFrame = PauliStabilizer[2][{"H" -> 1, "H" -> 2, "T" -> 1, "T" -> 2, "CZ" -> {1, 2}}];
Echo[Head[qfFrame], "D2 QF result head:"];
Echo[qfFrame["Length"], "D2 QF StabilizerFrame component count (true phase-poly rank = 1):"];
Echo[QuantumSimilarity[qfFrame["State"], denseDiag], "D2 QF frame still EXACT (=1) just not compressed:"];

(* Degree-3 term: a CCZ on |+++> contributes a cubic phase x1 x2 x3.        *)
(* Build the diagonal CCZ directly to avoid depending on a gate name.       *)
ccz = DiagonalMatrix[{1, 1, 1, 1, 1, 1, 1, -1}];
plus3 = (1/Sqrt[8]) ConstantArray[1, 8];
vCCZdense = Chop[ccz . plus3];
phi3[x1_, x2_, x3_] := Mod[4 x1 x2 x3, 8];  (* (-1)^{x1 x2 x3} = w^{4 x1 x2 x3} *)
vCCZpoly = Flatten @ Table[(1/Sqrt[8]) w^phi3[x1, x2, x3], {x1, 0, 1}, {x2, 0, 1}, {x3, 0, 1}];
Echo[vEqual[vCCZdense, vCCZpoly], "D3 cubic (degree-3) phase poly reproduces CCZ|+++>:"];

(* ======================================================================== *)
(* E. Hypothesis 4: interleaving H between T's forbids the diagonal form.    *)
(*    Frame component count = 2^(#T): the exponential is RELOCATED, not      *)
(*    removed. Materializing the state vector sums 2^k components.           *)
(* ======================================================================== *)

stabRank[x_] := If[MatchQ[x, _StabilizerFrame], x["Length"], 1];

(* Interleaved H/T on one qubit: H T H T H T ... -- each T doubles the rank. *)
interleaved[k_] := PauliStabilizer[1][Prepend[Flatten[Table[{"T", "H"}, k]], "H"]];
Echo[Table[{k, stabRank[interleaved[k]]}, {k, 0, 8}], "E1 (T-count k, frame rank) for interleaved H/T -> 2^k:"];

(* The amplitude after an interior H becomes a Gauss/exponential sum.        *)
(* For U = H T H: <x|U|0> = (1/2) Sum_{a in F_2} w^{a} (-1)^{a x}, an        *)
(* internal F_2 sum. Each interior H multiplies the number of summed terms,  *)
(* so a circuit with d interior H-layers carries a 2^d-term sum per amplitude.*)
ampHTH[x_] := (1/2) Sum[w^a (-1)^(a x), {a, 0, 1}];
denseHTH = QuantumCircuitOperator[{"H" -> 1, "T" -> 1, "H" -> 1}][QuantumState["0"], Method -> "Schrodinger"];
vHTHdense = Chop[denseHTH["StateVector"] // Normal];
vHTHsum = {ampHTH[0], ampHTH[1]};
Echo[vEqual[vHTHdense, vHTHsum], "E2 amplitude as an internal Gauss sum matches dense (H T H)|0>:"];

(* Benchmark: frame build + dense materialization cost grows as 2^k.        *)
ClearSystemCache[];
benchRow[k_] := Module[{f, t1, t2, bc},
    {t1, f} = AbsoluteTiming[interleaved[k]];
    bc = ByteCount[f];
    {t2} = AbsoluteTiming[f["StateVector"];][[{1}]];
    {k, stabRank[f], bc, Round[1000 t1], Round[1000 t2]}
];
Echo[Prepend[Table[ClearSystemCache[]; benchRow[k], {k, 2, 12, 2}],
    {"k", "rank=2^k", "ByteCount", "build ms", "materialize ms"}],
  "E3 cost vs T-count (ClearSystemCache each row):"];

(* ======================================================================== *)
(* F. Partial wins.                                                          *)
(* ======================================================================== *)

(* F1: a single expectation < psi | P | psi > on a frame is a sum of chi^2   *)
(* stabilizer inner products, each O(n): polynomial in n for fixed rank.     *)
fr = PauliStabilizer[1][{"H", "T"}];           (* T|+>, rank 2 *)
chi = fr["Length"];
(* <Y> on T|+> should be +1/sqrt2. Compute from the dense frame state.       *)
yExp = N[Re[Conjugate[#] . (PauliMatrix[2] . #)]] &[fr["StateVector"] // Normal];
Echo[{chi, chi^2, yExp}, "F1 {rank chi, chi^2 inner-product terms, <Y> on T|+>=+0.7071}:"];

(* F2: bounded T-count. QF keeps 2^t terms; the TRUE exact stabilizer rank   *)
(* of |T>^{otimes t} is far smaller (Bravyi et al: t+1 for t<=5, =t for      *)
(* 2<=t<=4). Tabulate QF's component count against the known exact rank.     *)
exactRankT = <|0 -> 1, 1 -> 1, 2 -> 2, 3 -> 3, 4 -> 4, 5 -> 5|>;  (* chi(T^t), BravyiEtAl18 *)
tcount[k_] := stabRank[PauliStabilizer[1][Join[{"H"}, Table["T", k]]]];
Echo[Prepend[Table[{k, 2^k, tcount[k], exactRankT[k]}, {k, 0, 5}],
    {"t", "2^t naive", "QF frame Length", "true exact rank chi"}],
  "F2 QF keeps 2^t, no compression to the true rank:"];

Echo["exp2 done", "STATUS:"];
