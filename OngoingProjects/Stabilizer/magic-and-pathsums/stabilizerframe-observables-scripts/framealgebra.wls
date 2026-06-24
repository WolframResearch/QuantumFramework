(* ============================================================================ *)
(* framealgebra.wls : out-of-tree StabilizerFrame observable prototypes.        *)
(*                                                                              *)
(* Loaded via Get by the verification scripts. Defines, in Global`:             *)
(*   symplectic Pauli algebra: sympTimes/sympDagger/sympId/sympFromRelating/    *)
(*     sympFromHermString/sympToHerm/refExpect                                  *)
(*   shared-reference core: frameSandwich -> frameExpectation/frameNorm2/        *)
(*     frameInnerProduct (exact, phase-correct, chi^2 poly(n))                  *)
(*   distinct-reference: pauliApplyVec/pairAmpDense (exact oracle, 2^n/pair),    *)
(*     pairMagClosed (poly magnitude), frameSandwichGeneralExact/                *)
(*     frameInnerProductGeneralExact, applyHermString                           *)
(*   Born: framePauliBorn/frameZmarginal                                        *)
(*   amplitudes: stabGenerators/stabXRank/stabSupportPoint/gaussF2WithPauli/     *)
(*     stabGroupElementWithX/stabRatio/frameAmplitudeRatio/frameCompProb/        *)
(*     frameAmplitude                                                           *)
(*                                                                              *)
(* Reuses two kernel PackageScope primitives, `stabilizerExpectation` (exact     *)
(* O(n^3) single-state Pauli expectation) and `stabilizerInnerProduct` (the      *)
(* closed form via Method -> "ClosedForm"), plus PauliStabilizer accessors. No   *)
(* kernel .m file is modified.                                                   *)
(* ============================================================================ *)

se = Wolfram`QuantumFramework`PackageScope`stabilizerExpectation;

(* ---- symplectic Pauli algebra: operator  s * X^x Z^z  ---------------------- *)
(* x,z are bit lists (qubit 1 = most significant, matching QF's IntegerDigits). *)
sympTimes[a_, b_] := <|
    "x" -> BitXor[a["x"], b["x"]],
    "z" -> BitXor[a["z"], b["z"]],
    "s" -> a["s"] b["s"] (-1)^(a["z"] . b["x"])
|>;
sympDagger[a_] := <|"x" -> a["x"], "z" -> a["z"], "s" -> Conjugate[a["s"]] (-1)^(a["x"] . a["z"])|>;
sympId[n_] := <|"x" -> ConstantArray[0, n], "z" -> ConstantArray[0, n], "s" -> 1|>;

(* relating Pauli {xb,zb,coeff} -> operator coeff*I^{x.z} X^x Z^z (kernel's     *)
(* sfApplyPauli convention).                                                     *)
sympFromRelating[{xb_, zb_, coeff_}] := <|"x" -> xb, "z" -> zb, "s" -> coeff I^(xb . zb)|>;

(* Hermitian Pauli string "-?[IXYZ]+" -> symp (Hermitian tensor i^{x.z} X^x Z^z)*)
sympFromHermString[str_String] := Module[{sign, body, chars, xb, zb},
    {sign, body} = If[StringStartsQ[str, "-"], {-1, StringDrop[str, 1]}, {1, str}];
    chars = Characters[body];
    xb = Replace[chars, {"I" -> 0, "X" -> 1, "Y" -> 1, "Z" -> 0}, {1}];
    zb = Replace[chars, {"I" -> 0, "X" -> 0, "Y" -> 1, "Z" -> 1}, {1}];
    <|"x" -> xb, "z" -> zb, "s" -> sign I^(xb . zb)|>
];

(* symp (x,z,s) -> {scalar, HermitianString} with op = scalar * HermTensor.     *)
sympToHerm[a_] := Module[{xb = a["x"], zb = a["z"], chars},
    chars = MapThread[
        Replace[{#1, #2}, {{0, 0} -> "I", {1, 0} -> "X", {1, 1} -> "Y", {0, 1} -> "Z"}] &,
        {xb, zb}];
    {a["s"] I^(- (xb . zb)), StringJoin[chars]}
];

(* <ref| (symp op) |ref> via the kernel's exact single-state expectation.       *)
refExpect[ref_, a_] := Module[{scalar, hstr},
    {scalar, hstr} = sympToHerm[a];
    scalar se[ref, hstr]
];

(* ---- core sandwich: requires fA, fB to share Components[[1,2]] (reference) -- *)
frameSandwich[fA_, Pobs_, fB_] := Module[{ref, aC, aP, bC, bP},
    ref = First[fA]["Components"][[1, 2]];
    aC = fA["Coefficients"]; aP = sympFromRelating /@ First[fA]["Paulis"];
    bC = fB["Coefficients"]; bP = sympFromRelating /@ First[fB]["Paulis"];
    Total @ Flatten @ Table[
        Conjugate[aC[[i]]] bC[[j]] *
            refExpect[ref, sympTimes[sympTimes[sympDagger[aP[[i]]], Pobs], bP[[j]]]],
        {i, Length[aC]}, {j, Length[bC]}
    ]
];

frameRefQubits[f_] := First[f]["Components"][[1, 2]]["Qubits"];
frameExpectation[f_, str_String] := frameSandwich[f, sympFromHermString[str], f];
frameNorm2[f_] := frameSandwich[f, sympId[frameRefQubits[f]], f];
frameInnerProduct[fA_, fB_] := frameSandwich[fA, sympId[frameRefQubits[fA]], fB];

(* ---- general two-reference sandwich (DISTINCT references) ------------------- *)
(* <psiA|Pobs|psiB> = Sum_ij conj(a_i) b_j <refA| PA_i^d Pobs PB_j |refB>, each   *)
(* term a phased Pauli transition amplitude <refA|Q|refB>.                        *)
(*                                                                              *)
(* IMPORTANT: a Pauli GATE on a PauliStabilizer is phase-blind (Z|1> -> |1>),    *)
(* so <refA|Q|refB> cannot be read off a gate-applied tableau. Two routes:       *)
(*   pairAmpDense   -- apply Q to refB's vector with explicit phase, then dot    *)
(*                     (2^n/pair, phase-correct: a correctness oracle).          *)
(*   pairMagClosed  -- stabilizerInnerProductClosedForm on the (phase-blind)     *)
(*                     gate-applied tableau (poly/pair, MAGNITUDE only).         *)
sip = Wolfram`QuantumFramework`PackageScope`stabilizerInnerProduct;
(* the closed form is reached through the public Method -> "ClosedForm" route     *)
(* (stabilizerInnerProductClosedForm itself is package-private, not PackageScope). *)
sipCF[a_, b_] := sip[a, b, Method -> "ClosedForm"];

applyHermString[ps_, str_] := Fold[
    Function[{s, qc}, If[qc[[2]] === "I", s, s[qc[[2]], qc[[1]]]]],
    ps, MapIndexed[{First[#2], #1} &, Characters[str]]];

(* operator (symp s X^x Z^z) applied to a dense vector, phase-correct. *)
pauliApplyVec[a_, v_, n_] := Module[{dim = 2^n, xmask = FromDigits[a["x"], 2], targets, signs},
    targets = 1 + BitXor[Range[0, dim - 1], xmask];
    signs = Table[(-1)^Mod[a["z"] . IntegerDigits[b, 2, n], 2], {b, 0, dim - 1}];
    SparseArray[Thread[targets -> a["s"] signs Normal[v]], dim]
];
pairAmpDense[refA_, Q_, refB_, n_] := Module[{vA = refA["State"]["StateVector"], vB = refB["State"]["StateVector"]},
    Conjugate[Normal @ vA] . Normal @ pauliApplyVec[Q, vB, n]
];
pairMagClosed[refA_, Q_, refB_] := Module[{scalar, qstr},
    {scalar, qstr} = sympToHerm[Q];
    Abs[scalar] sipCF[refA, applyHermString[refB, qstr]]
];

(* exact distinct-reference sandwich (oracle: phase-correct, 2^n/pair). *)
frameSandwichGeneralExact[fA_, Pobs_, fB_] := Module[{refA, refB, n, aC, aP, bC, bP},
    refA = First[fA]["Components"][[1, 2]]; refB = First[fB]["Components"][[1, 2]];
    n = refA["Qubits"];
    aC = fA["Coefficients"]; aP = sympFromRelating /@ First[fA]["Paulis"];
    bC = fB["Coefficients"]; bP = sympFromRelating /@ First[fB]["Paulis"];
    Total @ Flatten @ Table[
        Conjugate[aC[[i]]] bC[[j]] pairAmpDense[refA,
            sympTimes[sympTimes[sympDagger[aP[[i]]], Pobs], bP[[j]]], refB, n],
        {i, Length[aC]}, {j, Length[bC]}]
];
frameInnerProductGeneralExact[fA_, fB_] := frameSandwichGeneralExact[fA, sympId[frameRefQubits[fA]], fB];

(* ---- Born from the expectation primitive ---------------------------------- *)
(* P(P-eigenvalue = +1) for Hermitian Pauli string str, normalized.            *)
framePauliBorn[f_, str_String] := Module[{n2 = frameNorm2[f], ev = frameExpectation[f, str]},
    {(n2 + ev) / (2 n2), (n2 - ev) / (2 n2)}
];
(* Single-qubit Z marginal P(x_q = 0), P(x_q = 1). *)
frameZmarginal[f_, q_Integer] := Module[{n = frameRefQubits[f], str},
    str = StringJoin @ ReplacePart[ConstantArray["I", n], q -> "Z"];
    framePauliBorn[f, str]
];

(* ============================================================================ *)
(* Reference-state amplitudes (phase-correct, gauge-invariant ratios).          *)
(* gens: signed Hermitian generators as symp records {x,z,s}, s = sign*i^{x.z}. *)
(* ============================================================================ *)
stabGenerators[ps_] := Module[{gc = ps["GeneratorCount"], mat, signs, n},
    n = ps["Qubits"];
    mat = ps["Matrix"][[gc + 1 ;; 2 gc]];                 (* n rows, 2n symplectic *)
    signs = ps["StabilizerSigns"];
    MapThread[
        With[{x = #1[[1 ;; n]], z = #1[[n + 1 ;; 2 n]]},
            <|"x" -> x, "z" -> z, "s" -> #2 I^(x . z)|>] &,
        {mat, signs}]
];

(* X-part rank over F2 = log2(support size) of the stabilizer state.            *)
stabXRank[ps_] := MatrixRank[stabGenerators[ps][[All, "x"]], Modulus -> 2];

(* A computational support point x0 (poly): solve the pure-Z constraints.       *)
(* Reduce generator X-parts to row-echelon (tracking products) so the rows with *)
(* zero X-part are pure-Z stabilizers z.x = parity; solve over F2.              *)
stabSupportPoint[ps_] := Module[{n = ps["Qubits"], gens, xMat, redrows, zrows, A, bvec, sol},
    gens = stabGenerators[ps];
    (* Row-reduce [X | Z | parity] over F2 by Gaussian elimination on X.        *)
    {xMat, redrows} = gaussF2WithPauli[gens, n];
    (* pure-Z rows: X-part all zero *)
    zrows = Select[redrows, Total[#["x"]] == 0 &];
    If[zrows === {}, Return[ConstantArray[0, n]]];
    A = zrows[[All, "z"]];
    (* parity: pure-Z gen s*Z^z stabilizes x0 -> s*(-1)^{z.x0}=1.               *)
    (* s in {1,-1} for a pure-Z Hermitian gen; need z.x0 = (s==-1).             *)
    bvec = Map[If[# === -1, 1, 0] &, zrows[[All, "s"]]];
    (* A valid stabilizer state's pure-Z constraints are always consistent.     *)
    sol = LinearSolve[A, bvec, Modulus -> 2];
    Mod[sol, 2]
];

(* Gaussian elimination on generator X-parts over F2, multiplying generators    *)
(* (symplectic, phase-tracked) so the returned rows stay genuine stabilizers.   *)
(* Fold over columns carrying {generators, nextPivotRow}.                         *)
gaussF2WithPauli[gens0_, n_] := Module[{final},
    final = Fold[
        Function[{state, col},
            Module[{gens = state[[1]], piv = state[[2]], hit, pivoted},
                hit = SelectFirst[Range[piv, Length[gens]], gens[[#]]["x"][[col]] == 1 &];
                If[MissingQ[hit], state,
                    pivoted = ReplacePart[gens, {piv -> gens[[hit]], hit -> gens[[piv]]}];
                    {MapIndexed[
                        If[First[#2] =!= piv && #1["x"][[col]] == 1,
                            sympTimes[#1, pivoted[[piv]]], #1] &,
                        pivoted],
                     piv + 1}]]],
        {gens0, 1}, Range[n]];
    {final[[1]][[All, "x"]], final[[1]]}
];

(* Stabilizer-group element with X-part = a (a in column space of generator X). *)
(* Returns the symp record; selects generators via F2 LinearSolve on X-parts.   *)
stabGroupElementWithX[ps_, a_] := Module[{gens = stabGenerators[ps], Xt, coeffs, n = ps["Qubits"]},
    Xt = Transpose[gens[[All, "x"]]];                      (* n x (#gens) *)
    (* a in column space iff augmenting the target column does not raise rank.   *)
    If[MatrixRank[Xt, Modulus -> 2] != MatrixRank[Transpose @ Append[Transpose[Xt], a], Modulus -> 2],
        Return[$Failed]];
    coeffs = LinearSolve[Xt, a, Modulus -> 2];
    Fold[sympTimes, sympId[n], Pick[gens, coeffs, 1]]
];

(* <y|ref>/<x0|ref> for y,x0 in support. Choose g in S with X-part a = y(+)x0;  *)
(* then <y|ref> = <y|g|ref> = s_g (-1)^{b.x0} <x0|ref>, b = g's Z-part.          *)
stabRatio[ps_, y_, x0_] := Module[{a = BitXor[y, x0], g},
    g = stabGroupElementWithX[ps, a];
    If[g === $Failed, Return[0]];                          (* y not in support *)
    g["s"] (-1)^(g["z"] . x0)
];

(* Frame amplitude RATIO  <x|psi> / <x0|ref>, poly + gauge-invariant.           *)
(* <x|psi> = Sum_i c_i s_i (-1)^{z_i.(x(+)x_i)} <x(+)x_i|ref>, each ref         *)
(* amplitude expressed relative to the anchor x0.                               *)
frameAmplitudeRatio[f_, x_, x0_] := Module[{ref, n, coeffs, paulis},
    ref = First[f]["Components"][[1, 2]];
    n = ref["Qubits"];
    coeffs = f["Coefficients"];
    paulis = sympFromRelating /@ First[f]["Paulis"];
    Total @ MapThread[
        Function[{c, P},
            With[{y = BitXor[x, P["x"]]},
                c P["s"] (-1)^(P["z"] . y) stabRatio[ref, y, x0]]],
        {coeffs, paulis}]
];

(* Exact computational Born probability P(x) (poly, no global phase needed).     *)
(* |<x0|ref>|^2 = 2^{-k}, k = X-rank of the reference state.                     *)
frameCompProb[f_, x_] := Module[{ref = First[f]["Components"][[1, 2]], x0, k},
    x0 = stabSupportPoint[ref];
    k = stabXRank[ref];
    Abs[frameAmplitudeRatio[f, x, x0]]^2 * 2^(-k)
];

(* Full complex amplitude <x|psi>, pinned to the kernel's global-phase          *)
(* convention by reading ONE anchor amplitude <x0|ref> (the only non-poly step; *)
(* the kernel impl would supply it from the canonical CH-form).                  *)
frameAmplitude[f_, x_, anchorAmp_, x0_] := frameAmplitudeRatio[f, x, x0] anchorAmp;
