(* ::Package:: *)

(* elimination.wl : the path-sum rewrite engine (Stage 3).                         *)
(*                                                                                 *)
(* A single internal (un-pinned) Hadamard variable v is summed out using the exact  *)
(* identity  Sum_{v in {0,1}} omega^{phi_0 + v C} = omega^{phi_0} (1 + omega^C),     *)
(* C = Coefficient[phi, v] the cofactor (phi is linear in the idempotent v). The     *)
(* sum factors cleanly only in these classes (each an exact, degree-checked rule):   *)
(*                                                                                  *)
(*   CONST  C is constant c:        1 + omega^c   (c=0 -> x2; c=4 -> 0 / dead;        *)
(*                                                 odd/other c -> the scalar 1+w^c).  *)
(*   HH     C = 4 Q, Q Boolean LINEAR (nonconstant):  1+(-1)^Q = 2 [Q=0]              *)
(*                                                 -> x2, drop v, impose Q=0.         *)
(*   S      C = 2 + 4 L, L Boolean (C ident. 2 mod 4): 1 + i (-1)^L = Sqrt2 w^{1-2L}  *)
(*                                                 -> xSqrt2, drop v, phase += 1-2L   *)
(*                                                 (applied only if degree stays <=3).*)
(*                                                                                  *)
(* Everything else (odd nonconstant cofactor = magic; mixed C mod 4; nonlinear Q;     *)
(* a degree-violating S step) is SKIPPED: the variable stays and the brute-force      *)
(* Gauss sum PathSum`amp completes it EXACTLY. Correctness therefore never depends on  *)
(* the rewriter's completeness; the rewriter only shrinks the residual. Termination   *)
(* is guaranteed (every rule removes >=1 variable). The exact prefactor accumulates   *)
(* in "scale" (a Z[omega,1/Sqrt2] element); a dead branch sets scale = 0.             *)
(*                                                                                  *)
(* Assumes Wolfram`QuantumFramework`, pathsum.wl loaded.                            *)

(* total degree of a (multilinear, idempotent) phase polynomial *)
PathSum`pdeg[poly_] := Max[Append[Length[Variables[#]] & /@ MonomialList[PathSum`idem[poly]], 0]];

(* classify variable v in a state: returns {ruleName, data...} or {"skip"}. *)
PathSum`classify[ps_, v_] := Module[{C, m4, Q, L},
    If[! PathSum`internalQ[ps, v], Return[{"skip"}]];     (* output-pinned: not summed here *)
    C = PathSum`cofactor[ps["phase"], v];
    If[FreeQ[C, _x], Return[{"const", Mod[C, 8]}]];        (* constant cofactor *)
    m4 = PolynomialMod[C, 4];
    Which[
        m4 === 0,                                          (* C = 4 Q *)
            Q = PathSum`f2[C/4];
            If[Length[Variables[Q]] >= 1 && AllTrue[MonomialList[Q], Length[Variables[#]] <= 1 &],
                {"hh", Q}, {"skip"}],
        m4 === 2,                                          (* C = 2 + 4 L *)
            L = PathSum`f2[(C - 2)/4];
            {"s", L},
        True, {"skip"}                                     (* odd / mixed: magic, residual *)
    ]
];

(* apply the rule that fires on v (assumes classify returned non-"skip"). *)
PathSum`fire[ps_, v_, rule_] := Module[{phase = ps["phase"], val, u, rhs, newphase},
    Switch[First[rule],
        "const", With[{c = rule[[2]]},
            <|ps, "phase" -> PathSum`z8[phase /. v -> 0], "vars" -> DeleteCases[ps["vars"], v],
              "scale" -> PathSum`z8expand[ps["scale"] (1 + PathSum`omega^c)]|>],
        "hh", With[{Q = rule[[2]]},
            u = First[Variables[Q]]; rhs = PathSum`f2[u + Q];          (* impose Q=0: u = rhs (F_2) *)
            (* val is an F_2 polynomial: substitute the F_2 form. The phase treats u as an   *)
            (* INTEGER 0/1 value, so it must receive the integer value boolToInt[rhs] (a pure *)
            (* F_2 substitution is only valid when every phase coefficient is 0 mod 4).       *)
            val = PathSum`f2[#] & /@ (ps["val"] /. (u -> rhs));
            <|ps, "phase" -> PathSum`z8[(phase /. v -> 0) /. (u -> PathSum`boolToInt[rhs])], "val" -> val,
              "vars" -> DeleteCases[ps["vars"], v | u], "scale" -> ps["scale"] 2|>],
        "s", With[{L = rule[[2]]},
            newphase = PathSum`z8[(phase /. v -> 0) + 1 - 2 PathSum`boolToInt[L]];
            If[PathSum`pdeg[newphase] > 3, ps,             (* guard: keep degree <= 3 *)
                <|ps, "phase" -> newphase, "vars" -> DeleteCases[ps["vars"], v],
                  "scale" -> PathSum`z8expand[ps["scale"] Sqrt[2] PathSum`omega^0]|>]]
    ]
];
PathSum`z8expand[s_] := Simplify[s];   (* keep the exact scalar tidy in Z[omega,1/Sqrt2] *)

(* one reduction step: fire the first eligible internal variable, else fixed point. *)
PathSum`step[ps_] := If[ps["scale"] === 0, ps,
    Module[{v, rule},
        v = SelectFirst[ps["vars"], (First[PathSum`classify[ps, #]] =!= "skip") &, Missing[]];
        If[MissingQ[v], Return[ps]];
        rule = PathSum`classify[ps, v];
        PathSum`fire[ps, v, rule]
    ]];

(* reduce to a fixed point, with a degree-<=3 invariant assertion and a step ceiling   *)
(* (a safety net; termination is already guaranteed since each fired rule drops a var).*)
PathSum`engineReduce[ps0_, opts : OptionsPattern[]] := Module[
    {ps = ps0, ceiling = OptionValue["StepCeiling"], steps = 0, next},
    While[steps < ceiling,
        next = PathSum`step[ps];
        If[next === ps, Break[]];                          (* fixed point *)
        If[next["scale"] =!= 0 && PathSum`pdeg[next["phase"]] > 3,
            (* a rule violated the degree bound: refuse it, treat as fixed point *)
            Break[]];
        ps = next; steps++;
    ];
    Append[ps, "steps" -> steps]
];
Options[PathSum`engineReduce] = {"StepCeiling" -> 100000};

(* number of summation variables left after reduction (the residual the brute force    *)
(* still has to enumerate; internal ones are the meaningful magic, output ones pinned). *)
PathSum`residualVars[ps_] := ps["vars"];
PathSum`internalResidual[ps_, n_] := Select[ps["vars"], PathSum`internalQ[ps, #] &];

(* engine amplitude = reduce, then complete with the exact brute-force Gauss sum. *)
PathSum`engineAmp[n_Integer, gates_List, y_List] :=
    PathSum`amp[PathSum`engineReduce[PathSum`build[n, gates]], n, y];

(* differential matchers: engine === brute-force build === QF Schrodinger, exactly. *)
PathSum`engineMatchesQ[n_Integer, gates_List] := With[{r = PathSum`engineReduce[PathSum`build[n, gates]],
                                                       b = PathSum`build[n, gates]},
    AllTrue[Tuples[{0, 1}, n], PathSum`exactEqualQ[PathSum`amp[r, n, #],
        PathSum`qfAmp[n, gates, #]] && PathSum`exactEqualQ[PathSum`amp[r, n, #], PathSum`amp[b, n, #]] &]];

(* confluence replay: reduce under a SHUFFLED variable-selection order; the residual    *)
(* set may differ but the completed amplitude must be invariant.                        *)
PathSum`stepShuffled[ps_, perm_] := If[ps["scale"] === 0, ps,
    Module[{cand, v, rule},
        cand = Select[Permute[ps["vars"], perm], (First[PathSum`classify[ps, #]] =!= "skip") &];
        If[cand === {}, Return[ps]];
        v = First[cand]; rule = PathSum`classify[ps, v]; PathSum`fire[ps, v, rule]]];
PathSum`engineReduceShuffled[ps0_, seed_] := Module[{ps = ps0, next, k = 0},
    While[k < 100000, SeedRandom[seed + k];
        next = PathSum`stepShuffled[ps, RandomPermutation[Length[ps["vars"]]]];
        If[next === ps || (next["scale"] =!= 0 && PathSum`pdeg[next["phase"]] > 3), Break[]];
        ps = next; k++]; ps];
PathSum`confluentQ[n_Integer, gates_List, seeds_List] := With[{ref = PathSum`engineReduce[PathSum`build[n, gates]]},
    AllTrue[seeds, Function[s, With[{r = PathSum`engineReduceShuffled[PathSum`build[n, gates], s]},
        AllTrue[Tuples[{0, 1}, n], PathSum`exactEqualQ[PathSum`amp[r, n, #], PathSum`amp[ref, n, #]] &]]]]];

(* ============================================================================== *)
(* Formal Inactive[Sum] carrier (the plan's structural representation). A path-sum *)
(* head PS holds its arguments (HoldAll) and does not thread (NonThreadable); the   *)
(* residual is rendered as an inert Sum that Activate collapses to the amplitude.   *)
(* ============================================================================== *)
SetAttributes[PathSum`PS, {HoldAll, NonThreadable}];

PathSum`formalSum[ps_, n_Integer, y_List] := With[
    {vars = ps["vars"], ph = ps["phase"], sc = ps["scale"], h = ps["h"], vv = ps["val"]},
    With[{body = PathSum`omega^Mod[ph, 8] Boole[And @@ Table[Mod[vv[q], 2] == y[[q]], {q, n}]]},
        (sc / Sqrt[2]^h) If[vars === {}, body,        (* fully reduced: no inert Sum needed *)
            Inactive[Sum][body, Sequence @@ (With[{u = #}, {u, 0, 1}] & /@ vars)]]]];

PathSum`activate[fs_] := Activate[fs];   (* Activate -> the kernel evaluates the inert Sum *)
