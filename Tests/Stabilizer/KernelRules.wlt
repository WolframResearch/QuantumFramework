BeginTestSection["KernelRules"];

(* What each Clifford gate does to one row of the tableau, written as polynomials mod 2 in that
   row's bits: r is the row's sign bit, xa, za the bits of the gate's first qubit, xb, zb of its
   second, xq, zq of any other qubit. The engine computes a gate by three routes: the canonical
   array code, which symbolic tableaux take; the packed code, which a single gate on a concrete
   tableau takes (62 qubits per machine word); and the compiled gate loop that "ApplyCircuit" runs.
   Each route is checked against these polynomials below. *)
krGateRules = <|
    "H" -> <|"r" -> "r" + "xa" "za", "xa" -> "za", "za" -> "xa", "xq" -> "xq", "zq" -> "zq"|>,
    "S" -> <|"r" -> "r" + "xa" "za", "xa" -> "xa", "za" -> "xa" + "za", "xq" -> "xq", "zq" -> "zq"|>,
    SuperDagger["S"] -> <|"r" -> "r" + "xa" + "xa" "za", "xa" -> "xa", "za" -> "xa" + "za", "xq" -> "xq", "zq" -> "zq"|>,
    "X" -> <|"r" -> "r" + "za", "xa" -> "xa", "za" -> "za", "xq" -> "xq", "zq" -> "zq"|>,
    "Y" -> <|"r" -> "r" + "xa" + "za", "xa" -> "xa", "za" -> "za", "xq" -> "xq", "zq" -> "zq"|>,
    "Z" -> <|"r" -> "r" + "xa", "xa" -> "xa", "za" -> "za", "xq" -> "xq", "zq" -> "zq"|>,
    "CNOT" -> <|"r" -> "r" + "xa" "zb" + "xa" "xb" "zb" + "xa" "za" "zb", "xa" -> "xa", "za" -> "za" + "zb",
        "xb" -> "xa" + "xb", "zb" -> "zb", "xq" -> "xq", "zq" -> "zq"|>,
    "SWAP" -> <|"r" -> "r", "xa" -> "xb", "za" -> "zb", "xb" -> "xa", "zb" -> "za", "xq" -> "xq", "zq" -> "zq"|>
|>;
krArity[g_] := If[MemberQ[{"CNOT", "SWAP"}, g], 2, 1];

krN = 3;
krPlacements[1] = List /@ Range[krN];
krPlacements[2] = Permutations[Range[krN], {2}];
krRoleName[c_, q_, qs_] := c <> Switch[Position[qs, q], {{1}}, "a", {{2}}, "b", _, "q"];

(* One row {r, x_1..x_k, z_1..z_k} after gate g on qubits qs, as the polynomials say. *)
krEvalRow[rule_, qs_, row_] := With[{k = (Length[row] - 1) / 2},
    With[{xs = row[[2 ;; k + 1]], zs = row[[k + 2 ;;]]},
        With[{base = Join[{"r" -> row[[1]], "xa" -> xs[[qs[[1]]]], "za" -> zs[[qs[[1]]]]},
                If[Length[qs] == 2, {"xb" -> xs[[qs[[2]]]], "zb" -> zs[[qs[[2]]]]}, {}]]},
            Mod[Join[{rule["r"] /. base},
                Table[rule[krRoleName["x", q, qs]] /. Join[base, {"xq" -> xs[[q]]}], {q, k}],
                Table[rule[krRoleName["z", q, qs]] /. Join[base, {"zq" -> zs[[q]]}], {q, k}]], 2]]]];

krTableauOf[rows_] := With[{k = (Length[First[rows]] - 1) / 2},
    PauliStabilizer[<|"Signs" -> 1 - 2 rows[[All, 1]],
        "Tableau" -> {Transpose[rows[[All, 2 ;; k + 1]]], Transpose[rows[[All, k + 2 ;;]]]}|>]];
krRowsOf[ps_] := Transpose[Join[{(1 - ps["Signs"]) / 2}, Normal[ps["Tableau"]][[1]], Normal[ps["Tableau"]][[2]]]];

(* The input check accepts only the bits 0, 1 and -1; inside these blocks it accepts any array of
   the right shape, so the kernel code runs unchanged on symbols. *)
SetAttributes[krSymbolic, HoldAll];
krSymbolic[body_] := Internal`InheritedBlock[{Wolfram`QuantumFramework`PackageScope`PauliTableauQ},
    Unprotect[Wolfram`QuantumFramework`PackageScope`PauliTableauQ];
    ClearAll[Wolfram`QuantumFramework`PackageScope`PauliTableauQ];
    Wolfram`QuantumFramework`PackageScope`PauliTableauQ[t_] := ArrayQ[t, 3] && MatchQ[Dimensions[t], {2, _, _}];
    body];

(* The canonical route on a tableau of symbolic bits: each output bit of row i must depend only on
   row i's bits and equal the polynomial on all 2^7 values of them. *)
krSymTab = {Array[krx, {krN, 2 krN}], Array[krz, {krN, 2 krN}]};
krSymSigns = Table[1 - 2 krr[i], {i, 2 krN}];
krRowVars[i_] := Join[{krr[i]}, Table[krx[q, i], {q, krN}], Table[krz[q, i], {q, krN}]];
krCanonicalOK[g_, qs_] := With[{res = krSymbolic[With[{out = First[PauliStabilizer[<|"Signs" -> krSymSigns, "Tableau" -> krSymTab|>][g, Sequence @@ qs]]},
        {Normal[out["Tableau"]], (1 - out["Signs"]) / 2}]]},
    With[{t = res[[1]], sbits = res[[2]]},
        AllTrue[Range[2 krN], Function[i,
            AllTrue[Tuples[{0, 1}, 2 krN + 1], Function[vals,
                (Join[{sbits[[i]]}, t[[1, All, i]], t[[2, All, i]]] /. Thread[krRowVars[i] -> vals]) ===
                    krEvalRow[krGateRules[g], qs, vals]]]]]]];

KeyValueMap[Function[{g, rule}, VerificationTest[
    AllTrue[krPlacements[krArity[g]], krCanonicalOK[g, #] &], True,
    TestID -> "KernelRules-Canonical-" <> ToString[g, InputForm]]], krGateRules];

(* The compiled gate loop falls back to its uncompiled source on any runtime error. A copy with
   that fallback switched off throws instead, so a test passing through it was computed by the
   compiled code itself. *)
krCF = Wolfram`QuantumFramework`PackageScope`compiledCliffordFold;
krLoud = ReplacePart[krCF, {{1, 3} -> BitAnd[krCF[[1, 3]], BitNot[1024]], 8 -> Function[Throw["fallback", "krLoud"]]}];
VerificationTest[
    Catch[krLoud @@ {{{a}}, {{b}}, {c}, {1}, {{0, 0, 0}}}, "krLoud"],
    "fallback",
    {CompiledFunction::cfta},
    TestID -> "KernelRules-CompiledFallbackSwitchedOff"];

krRun["Packed", g_, qs_, ps_] := ps[g, Sequence @@ qs];
krRun["Compiled", g_, qs_, ps_] := With[{loud = krLoud},
    Block[{Wolfram`QuantumFramework`PackageScope`compiledCliffordFold = loud}, ps["ApplyCircuit", {g -> qs}]]];
krRouteKey = <|"Packed" -> "PackedX", "Compiled" -> "GeneratorX"|>;
krRouteOK[route_, g_, qs_, rows_] := With[{out = Catch[krRun[route, g, qs, krTableauOf[rows]], "krLoud"]},
    out =!= "fallback" && KeyExistsQ[First[out], krRouteKey[route]] &&
        krRowsOf[out] === (krEvalRow[krGateRules[g], qs, #] & /@ rows)];

(* n = 3: every one of the 2^7 values a row can take, six rows per tableau. *)
krRowBatches = PadRight[#, {2 krN, 2 krN + 1}] & /@ Partition[Tuples[{0, 1}, 2 krN + 1], 2 krN, 2 krN, 1, {}];
Scan[Function[route, KeyValueMap[Function[{g, rule}, VerificationTest[
    AllTrue[Tuples[{krPlacements[krArity[g]], krRowBatches}], krRouteOK[route, g, #[[1]], #[[2]]] &], True,
    TestID -> "KernelRules-" <> route <> "-AllRows-" <> ToString[g, InputForm]]], krGateRules]],
    {"Packed", "Compiled"}];

(* n = 64: qubits 62 and 63 sit in different packed words and the 128 generators span three
   compiled words. The gate's qubits and the sign run through all their values; other bits are random. *)
SeedRandom[20261108];
krEdge = {1, 62, 63, 64};
krBigRows[qs_] := With[{cols = Join[{1}, 1 + qs, 1 + 64 + qs]},
    With[{rnd = RandomInteger[1, {128, 129}], pats = Tuples[{0, 1}, Length[cols]]},
        ReplacePart[rnd, Flatten[MapIndexed[Function[{row, i},
            Thread[({i[[1]], #} & /@ cols) -> pats[[Mod[i[[1]] - 1, Length[pats]] + 1]]]], rnd], 1]]]];
Scan[Function[route, KeyValueMap[Function[{g, rule}, VerificationTest[
    AllTrue[If[krArity[g] == 1, List /@ krEdge, Permutations[krEdge, {2}]], krRouteOK[route, g, #, krBigRows[#]] &], True,
    TestID -> "KernelRules-" <> route <> "-WordEdges-" <> ToString[g, InputForm]]], krGateRules]],
    {"Packed", "Compiled"}];


(* ---- Measurement ---- *)

(* agPhase is the power of i in the product of two one-qubit Paulis, Aaronson-Gottesman's g. *)
krG[x1_, z1_, x2_, z2_] := Which[x1 == 0 && z1 == 0, 0, x1 == 1 && z1 == 1, z2 - x2,
    x1 == 1, z2 (2 x2 - 1), True, x2 (1 - 2 z2)];
VerificationTest[
    AllTrue[Tuples[{0, 1}, 4], Wolfram`QuantumFramework`PackageScope`agPhase @@ # === krG @@ # &],
    True,
    TestID -> "KernelRules-AgPhase"];

(* rowsum[t, h, i] multiplies row i into row h: the bits add, and the sign bit becomes 1 exactly when
   2 r_h + 2 r_i + (sum of agPhase over qubits) is 2 mod 4. Checked on symbolic bits with agPhase
   left unevaluated, for every pair of rows at n = 3. *)
krRowsumTemplate[h_, i_] := {
    ReplacePart[krSymSigns, h -> 1 - 2 Boole[Mod[2 krr[h] + 2 krr[i] +
        Sum[Wolfram`QuantumFramework`PackageScope`agPhase[krx[q, i], krz[q, i], krx[q, h], krz[q, h]], {q, krN}], 4] == 2]],
    {ReplacePart[krSymTab[[1]], {q_, h} :> BitXor[krx[q, h], krx[q, i]]],
     ReplacePart[krSymTab[[2]], {q_, h} :> BitXor[krz[q, h], krz[q, i]]]}};
VerificationTest[
    Block[{Wolfram`QuantumFramework`PackageScope`agPhase},
        AllTrue[Select[Tuples[Range[2 krN], 2], Apply[UnsameQ]],
            Wolfram`QuantumFramework`PackageScope`rowsum[{krSymSigns, krSymTab}, Sequence @@ #] === krRowsumTemplate @@ # &]],
    True,
    TestID -> "KernelRules-Rowsum"];

(* Measuring qubit a branches only on the X bits of column a. With those fixed and every other bit
   symbolic, the canonical measurement makes these row products, in this order:
   - if some stabilizer row p has an X on a (the first one), row p is multiplied into every other row
     with an X on a, in increasing order; then row p - n takes row p, and row p becomes +-Z_a;
   - otherwise each destabilizer row i with an X on a has its stabilizer row n + i multiplied into
     a scratch row 2n + 1, whose sign is the outcome.
   The output must be exactly what these products and final steps give. *)
krTabFor[a_, col_] := ReplacePart[krSymTab, Thread[({1, a, #} & /@ Range[2 krN]) -> col]];
krCallsRule[col_] := With[{xs = Flatten[Position[col, 1]]},
    With[{p = SelectFirst[xs, GreaterThan[krN]]},
        If[MissingQ[p], {2 krN + 1, # + krN} & /@ Select[xs, # <= krN &], {#, p} & /@ DeleteCases[xs, p]]]];
krTraced[a_, col_] := krSymbolic[Internal`InheritedBlock[{Wolfram`QuantumFramework`PackageScope`rowsum},
    DownValues[Wolfram`QuantumFramework`PackageScope`rowsum] = Prepend[DownValues[Wolfram`QuantumFramework`PackageScope`rowsum],
        HoldPattern[Wolfram`QuantumFramework`PackageScope`rowsum[s_, h_Integer, i_Integer] /; ! TrueQ[$krTracing]] :>
            Block[{$krTracing = True}, Sow[{h, i}]; Wolfram`QuantumFramework`PackageScope`rowsum[s, h, i]]];
    {#1, Flatten[#2, 1]} & @@ Reap[PauliStabilizer[<|"Signs" -> krSymSigns, "Tableau" -> krTabFor[a, col]|>]["Measure", a]]]];
krRebuilt[a_, col_, calls_] := With[{t = krTabFor[a, col], rs = Wolfram`QuantumFramework`PackageScope`rowsum,
        mk = PauliStabilizer[<|"Signs" -> #1, "Tableau" -> #2|>] &},
    krSymbolic[If[calls =!= {} && calls[[1, 1]] === 2 krN + 1,
        <|(1 - Last[#[[1]]]) / 2 -> mk[krSymSigns, t]|> & @ Fold[rs[#1, Sequence @@ #2] &,
            {Append[krSymSigns, 1], PadRight[t, Dimensions[t] + {0, 0, 1}]}, calls],
        With[{p = SelectFirst[Flatten[Position[col, 1]], GreaterThan[krN]]},
            If[MissingQ[p], <|0 -> mk[krSymSigns, t]|>,
                Block[{r2, t2},
                    {r2, t2} = Fold[rs[#1, Sequence @@ #2] &, {krSymSigns, t}, calls];
                    t2[[All, All, p - krN]] = t2[[All, All, p]];
                    t2[[All, All, p]] = 0;
                    t2[[2, a, p]] = 1;
                    AssociationMap[mk[ReplacePart[r2, p -> 1 - 2 #], t2] &, {0, 1}]]]]]]];
krCols = Tuples[{0, 1}, 2 krN];
Scan[Function[a, VerificationTest[
    AllTrue[krCols, Function[col, With[{tr = krTraced[a, col]},
        tr[[2]] === krCallsRule[col] && tr[[1]] === krRebuilt[a, col, tr[[2]]]]]],
    True,
    TestID -> "KernelRules-MeasureRowProducts-Qubit" <> ToString[a]]], Range[krN]];

(* The packed measurement adds up the powers of i with a two-word counter: each bit position of the
   two words is a separate counter mod 4. Over one to three one-bit chunks, every step it can take
   is fed in; then random full words over three chunks. *)
krPRS = Wolfram`QuantumFramework`PackageScope`packedRowSumPhase;
krAg = Wolfram`QuantumFramework`PackageScope`agPhase;
VerificationTest[
    AllTrue[Range[3], Function[k, AllTrue[Tuples[{0, 1}, {4 k}], Function[bits,
        With[{xi = bits[[1 ;; k]], zi = bits[[k + 1 ;; 2 k]], xh = bits[[2 k + 1 ;; 3 k]], zh = bits[[3 k + 1 ;;]]},
            krPRS[Transpose[{xi, xh}], Transpose[{zi, zh}], 1, 2] === Mod[Total[MapThread[krAg, {xi, zi, xh, zh}]], 4]]]]]],
    True,
    TestID -> "KernelRules-PackedPhaseCounter-EveryStep"];
SeedRandom[20261108];
VerificationTest[
    And @@ Table[With[{b = RandomInteger[1, {4, 186}]},
        With[{w = Map[FromDigits[Reverse[#], 2] &, Partition[#, 62]] & /@ b},
            krPRS[Transpose[w[[{1, 3}]]], Transpose[w[[{2, 4}]]], 1, 2] === Mod[Total[MapThread[krAg, b]], 4]]], 200],
    True,
    TestID -> "KernelRules-PackedPhaseCounter-FullWords"];

(* The packed measurement procedure, with one qubit per word so it runs on symbolic bits (every step
   acts on whole words by BitXor or by copying, so the word size only groups qubits) and with the
   counter replaced by the agPhase sum it was just checked to equal. It must make the same row
   products as the canonical code and give the same output, for every column pattern and qubit. *)
krAsBits[res_Association, oneBitWords_ : False] := MapAll[Expand, KeyValueMap[
    #1 -> If[oneBitWords && KeyExistsQ[First[#2], "PackedX"],
        {(1 - First[#2]["Signs"]) / 2, {First[#2]["PackedX"], First[#2]["PackedZ"]}},
        {(1 - #2["Signs"]) / 2, Normal[#2["Tableau"]]}] &, res]];
krPackedTraced[a_, col_] := krSymbolic[Internal`InheritedBlock[{Wolfram`QuantumFramework`PackageScope`packedRowSumPhase},
    DownValues[Wolfram`QuantumFramework`PackageScope`packedRowSumPhase] = Prepend[
        DownValues[Wolfram`QuantumFramework`PackageScope`packedRowSumPhase],
        HoldPattern[Wolfram`QuantumFramework`PackageScope`packedRowSumPhase[x_, z_, i_Integer, h_Integer]] :>
            (Sow[{h, i}]; Total[MapThread[Wolfram`QuantumFramework`PackageScope`agPhase, {x[[All, i]], z[[All, i]], x[[All, h]], z[[All, h]]}]])];
    Block[{Wolfram`QuantumFramework`PackageScope`agPhase, Wolfram`QuantumFramework`PackageScope`$chunkWidth = 1},
        {krAsBits[#1, True], Flatten[#2, 1]} & @@ Reap[Wolfram`QuantumFramework`PackageScope`packedMeasureZ[
            PauliStabilizer[<|"Signs" -> krSymSigns, "Tableau" -> krTabFor[a, col]|>], a]]]]];
krCanonicalInert[a_, col_] := krSymbolic[Block[{Wolfram`QuantumFramework`PackageScope`agPhase},
    krAsBits[PauliStabilizer[<|"Signs" -> krSymSigns, "Tableau" -> krTabFor[a, col]|>]["Measure", a]]]];
VerificationTest[
    AllTrue[krCols, Function[col, AllTrue[Range[krN], Function[a,
        With[{pt = krPackedTraced[a, col]}, pt[[1]] === krCanonicalInert[a, col] && pt[[2]] === krCallsRule[col]]]]]],
    True,
    TestID -> "KernelRules-PackedMeasureProcedure"];

(* The packed route as it runs, 62 qubits per word, against the canonical route on random tableaux,
   measuring qubits on both sides of the word edge at n = 64. A random outcome comes back packed; a
   determined one returns the input. *)
krCanonicalRun[ps_, a_] := Internal`InheritedBlock[{Wolfram`QuantumFramework`PackageScope`psConcreteFastQ},
    DownValues[Wolfram`QuantumFramework`PackageScope`psConcreteFastQ] =
        {HoldPattern[Wolfram`QuantumFramework`PackageScope`psConcreteFastQ[_]] :> False};
    ps["Measure", a]];
krRandomPS[k_] := PauliStabilizer[<|"Signs" -> RandomChoice[{-1, 1}, 2 k], "Tableau" -> RandomInteger[1, {2, k, 2 k}]|>];
SeedRandom[20261108];
Scan[Function[c, VerificationTest[
    AllTrue[Tuples[{c[[2]], Range[c[[3]]]}], Function[t, With[{ps = krRandomPS[c[[1]]]}, With[{fast = ps["Measure", t[[1]]]},
        (Length[fast] == 1 || AllTrue[Values[fast], KeyExistsQ[First[#], "PackedX"] &]) &&
            krAsBits[fast] === krAsBits[krCanonicalRun[ps, t[[1]]]]]]]],
    True,
    TestID -> "KernelRules-PackedMeasureLive-n" <> ToString[c[[1]]]]],
    {{3, Range[3], 100}, {64, {1, 62, 63, 64}, 10}}];


(* ---- Expectation value and inner product ---- *)

(* Reference values from matrices: the state is the product of (1 + S)/2 over its signed stabilizer
   rows, with the bits (x, z) = (1, 1) standing for Y. *)
krPauli[{0, 0}] = IdentityMatrix[2]; krPauli[{1, 0}] = PauliMatrix[1];
krPauli[{1, 1}] = PauliMatrix[2]; krPauli[{0, 1}] = PauliMatrix[3];
krRowMatrix[x_, z_] := KroneckerProduct @@ (krPauli /@ Transpose[{x, z}]);
krRho[ps_] := With[{k = ps["Qubits"], t = Normal[ps["Tableau"]], s = ps["Signs"]},
    Dot @@ Table[(IdentityMatrix[2^k] + s[[i]] krRowMatrix[t[[1, All, i]], t[[2, All, i]]]) / 2, {i, k + 1, 2 k}]];
krLetters = <|{0, 0} -> "I", {1, 0} -> "X", {1, 1} -> "Y", {0, 1} -> "Z"|>;
krString[sign_, x_, z_] := If[sign == -1, "-", ""] <> StringJoin[krLetters /@ Transpose[{x, z}]];

(* Every valid tableau at n = 2: rows (x | z) whose symplectic pairings are those of the standard
   destabilizer/stabilizer pairs. *)
krSymp[u_, v_, k_] := Mod[u[[;; k]] . v[[k + 1 ;;]] + u[[k + 1 ;;]] . v[[;; k]], 2];
krJ[k_] := ArrayFlatten[{{0, IdentityMatrix[k]}, {IdentityMatrix[k], 0}}];
krValid2 = Select[Tuples[Tuples[{0, 1}, 4], 4], Outer[krSymp[#1, #2, 2] &, #, #, 1] === krJ[2] &];
krFromRows[rows_, signs_] := PauliStabilizer[<|"Signs" -> signs,
    "Tableau" -> {Transpose[rows[[All, ;; Length[First[rows]] / 2]]], Transpose[rows[[All, Length[First[rows]] / 2 + 1 ;;]]]}|>];
VerificationTest[Length[krValid2], 720, TestID -> "KernelRules-ValidTableaux-n2"];

SeedRandom[20261108];
VerificationTest[
    AllTrue[krValid2, Function[rows, With[{ps = krFromRows[rows, RandomChoice[{-1, 1}, 4]]},
        With[{rho = krRho[ps]}, AllTrue[Tuples[{{1, -1}, Tuples[{0, 1}, 4]}], Function[p,
            ps["Expectation", krString[p[[1]], p[[2, ;; 2]], p[[2, 3 ;;]]]] ===
                p[[1]] Tr[rho . krRowMatrix[p[[2, ;; 2]], p[[2, 3 ;;]]]]]]]]]],
    True,
    TestID -> "KernelRules-Expectation-AllTableaux-n2"];

(* The coefficients LinearSolve returns for a Pauli string in the stabilizer span are its pairings
   with the destabilizers: the unique solution on a valid tableau. *)
krSolveRecorded[ps_, str_] := Internal`InheritedBlock[{LinearSolve},
    Unprotect[LinearSolve];
    LinearSolve[m_, b_, Modulus -> 2] /; ! TrueQ[$krInSolve] := Block[{$krInSolve = True},
        With[{sol = LinearSolve[m, b, Modulus -> 2]}, Sow[sol]; sol]];
    Reap[ps["Expectation", str]][[2]]];
VerificationTest[
    AllTrue[krValid2, Function[rows, With[{ps = krFromRows[rows, {1, 1, 1, 1}]},
        AllTrue[Select[Tuples[{0, 1}, 4], Function[p, AllTrue[rows[[3 ;;]], krSymp[#, p, 2] == 0 &]]], Function[p,
            krSolveRecorded[ps, krString[1, p[[;; 2]], p[[3 ;;]]]] === {{krSymp[#, p, 2] & /@ rows[[;; 2]]}}]]]]],
    True,
    TestID -> "KernelRules-ExpectationCoefficients-n2"];

(* Random stabilizer states at n = 3, reached by random circuits from |000>. *)
krRandomCircuit[k_, len_] := Table[With[{g = RandomChoice[Keys[krGateRules]]},
    g -> If[krArity[g] == 2, RandomSample[Range[k], 2], RandomInteger[{1, k}]]], {len}];
krRandomState[k_] := PauliStabilizer[k]["ApplyCircuit", krRandomCircuit[k, 30]];
SeedRandom[20261108];
VerificationTest[
    And @@ Table[With[{ps = krRandomState[3], p = RandomInteger[1, 6], sign = RandomChoice[{1, -1}]},
        ps["Expectation", krString[sign, p[[;; 3]], p[[4 ;;]]]] === sign Tr[krRho[ps] . krRowMatrix[p[[;; 3]], p[[4 ;;]]]]], 300],
    True,
    TestID -> "KernelRules-Expectation-Random-n3"];

(* The closed-form inner product is |<psi|phi>| = Sqrt[Tr[rho_psi . rho_phi]]. *)
SeedRandom[20261108];
Scan[Function[k, VerificationTest[
    And @@ Table[With[{a = krRandomState[k], b = krRandomState[k]},
        a["InnerProduct", b, Method -> "ClosedForm"] === Sqrt[Tr[krRho[a] . krRho[b]]]], 200],
    True,
    TestID -> "KernelRules-InnerProductClosedForm-n" <> ToString[k]]], {2, 3}];

(* The same states with the tableau stored as a SparseArray give the same values. *)
krSparse[ps_] := PauliStabilizer[<|"Signs" -> ps["Signs"], "Tableau" -> SparseArray[Normal[ps["Tableau"]]]|>];
SeedRandom[20261108];
VerificationTest[
    And @@ Table[With[{ps = krRandomState[3]}, With[{sp = krSparse[ps]},
        Head[sp["Tableau"]] === SparseArray &&
            AllTrue[StringJoin /@ Tuples[{"I", "X", "Y", "Z"}, 3], sp["Expectation", #] === ps["Expectation", #] &]]], 50],
    True,
    TestID -> "KernelRules-Expectation-SparseTableau"];
VerificationTest[
    And @@ Table[With[{a = krRandomState[3], b = krRandomState[3]},
        krSparse[a]["InnerProduct", krSparse[b], Method -> "ClosedForm"] === a["InnerProduct", b, Method -> "ClosedForm"]], 100],
    True,
    TestID -> "KernelRules-InnerProductClosedForm-SparseTableau"];

EndTestSection[];
