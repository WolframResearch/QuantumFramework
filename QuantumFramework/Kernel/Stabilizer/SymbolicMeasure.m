Package["Wolfram`QuantumFramework`"]

PackageScope[symbolicMeasure]
PackageScope[substituteOutcomes]
PackageScope[sampleOutcomes]



(* ============================================================================ *)
(* Symbolic measurement outcomes (FangYing23 SymPhase).                         *)
(*                                                                              *)
(* The standard `ps["M", q]` returns an Association of conditional outcomes,    *)
(* requiring the user to RandomChoice for an actual sample and re-traverse the  *)
(* circuit per sample. SymPhase amortizes this: a single circuit traversal      *)
(* produces a `PauliStabilizer` whose phases contain fresh F_2 symbols          *)
(* representing un-resolved random outcomes. Sampling then substitutes those    *)
(* symbols with random 0/1 values (one matrix-multiply over many samples).      *)
(*                                                                              *)
(* Public surface (method-grade on PauliStabilizer):                            *)
(*   ps["SymbolicMeasure", q]                                                   *)
(*   ps["SubstituteOutcomes", rules]                                            *)
(*   ps["SampleOutcomes", n]                                                    *)
(* The internal helpers below are PackageScope only; method dispatch lives in   *)
(* Stabilizer/Properties.m.                                                     *)
(*                                                                              *)
(* Reference: FangYing23, arxiv:2311.03906, Section 3.                          *)
(* ============================================================================ *)


(* ============================================================================ *)
(* Fresh-symbol allocator                                                       *)
(* ============================================================================ *)

(* Each call to symbolicMeasure[..., "Symbolic" -> True] allocates a fresh   *)
(* symbol in the form `\[FormalS][k]`. The counter is local to the user's      *)
(* session via $StabilizerSymbolCounter; reset by the user to keep symbol      *)
(* names short across long-running sessions.                                    *)

PackageScope[$StabilizerSymbolCounter]
$StabilizerSymbolCounter = 0

freshOutcome[] := \[FormalS][++$StabilizerSymbolCounter]


(* ============================================================================ *)
(* symbolicMeasure[ps, q]: symbolic Z-basis measurement.                      *)
(*                                                                              *)
(* Like `ps["M", q]` but instead of branching on outcome, allocates a fresh     *)
(* symbol \[FormalS][k] for the random outcome and returns ONE PauliStabilizer  *)
(* with that symbol embedded in the appropriate sign position.                  *)
(*                                                                              *)
(* For deterministic measurements (no anticommuting stabilizer), returns the    *)
(* PauliStabilizer with a concrete 0/1 outcome (no fresh symbol).               *)
(*                                                                              *)
(* Deterministic outcome after a prior symbolic measurement:                    *)
(* A deterministic measurement (no anticommuting stabilizer) that follows a     *)
(* symbolic one can have an outcome that is a polynomial in earlier fresh       *)
(* symbols (e.g. m_2 = m_1 for the Bell ZZ correlation). That correlation is    *)
(* not stamped into the post-state's signs, which carry the bare conserved      *)
(* generator (e.g. ZZZ), so reading the signs alone does not reveal it.         *)
(* The outcome polynomial is instead recorded in the "Outcomes" map (below),    *)
(* and SubstituteOutcomes / SampleOutcomes apply that map after substitution,   *)
(* so the correlation is resolved correctly in every sample.                    *)
(* ============================================================================ *)

(* Deterministic-outcome polynomials are tracked via an "Outcomes" key in the  *)
(* PauliStabilizer association: {fresh_symbol -> polynomial_in_prior_symbols}  *)
(* for each deterministic measurement whose outcome involves prior symbols.   *)
(* substituteOutcomes / sampleOutcomes apply this map iteratively after user  *)
(* / random rules so derived outcomes (e.g. Bell ZZ correlation m_2 = m_1)   *)
(* resolve correctly.                                                          *)

symbolicMeasure[ps_PauliStabilizer ? PauliStabilizerQ, a_Integer] := Module[{
    result, prevOutcomes
},
    result = ps["M", a];
    prevOutcomes = Lookup[First[ps], "Outcomes", <||>];
    Switch[Length[result],
        1,
            Module[{outcomeBit = First @ Keys @ result, postPS = First @ Values @ result, sym, newAssoc},
                newAssoc = First[postPS];
                If[FreeQ[outcomeBit, _\[FormalS]],
                    (* Outcome is concrete; just forward Outcomes map. *)
                    If[Length[prevOutcomes] > 0,
                        PauliStabilizer[Append[newAssoc, "Outcomes" -> prevOutcomes]],
                        postPS
                    ],
                    (* Outcome is a polynomial in earlier symbols. Allocate a  *)
                    (* fresh symbol and record fresh -> polynomial.            *)
                    sym = freshOutcome[];
                    PauliStabilizer[Append[newAssoc, "Outcomes" -> Append[prevOutcomes, sym -> outcomeBit]]]
                ]
            ],
        2,
            With[{ps0 = result[0]},
                Module[{sym = freshOutcome[], phase0, phase1, p, newAssoc},
                    phase0 = ps0["Phase"];
                    phase1 = result[1]["Phase"];
                    p = First @ Position[Mod[phase0 - phase1, 2], 1, {1}, 1];
                    newAssoc = <|
                        "Phase" -> ReplacePart[phase0, p -> sym],
                        "Tableau" -> ps0["Tableau"]
                    |>;
                    If[Length[prevOutcomes] > 0,
                        PauliStabilizer[Append[newAssoc, "Outcomes" -> prevOutcomes]],
                        PauliStabilizer[newAssoc]
                    ]
                ]
            ],
        _, $Failed
    ]
]


(* Multi-qubit symbolic measurement: fold over the qubit list *)
symbolicMeasure[ps_PauliStabilizer ? PauliStabilizerQ, qudits : {___Integer}] := Fold[
    symbolicMeasure[#1, #2] &, ps, qudits
]


(* ============================================================================ *)
(* substituteOutcomes[ps, rules]: replace measurement-outcome symbols with      *)
(* concrete 0/1 values, then reduce signs back to {-1, 1} via Mod 2.            *)
(* ============================================================================ *)

(* Apply user rules first, then iteratively close under the "Outcomes" map    *)
(* (deterministic-outcome polynomials) until a fixed point so derived         *)
(* outcomes resolve correctly (e.g. Bell m_2 = m_1).                           *)
substituteOutcomes[ps_PauliStabilizer ? PauliStabilizerQ, rules_] := Module[{
    phase, outcomesMap, allRules, fixed
},
    phase = ps["Phase"];
    outcomesMap = Normal @ Lookup[First[ps], "Outcomes", <||>];
    (* Combine user rules and outcomes map; iterate to fixed point on the     *)
    (* phase (closes any chain of derived outcomes).                           *)
    allRules = Join[Flatten[{rules}], outcomesMap];
    fixed = FixedPoint[(# /. allRules) &, phase, 10];
    PauliStabilizer[<|"Phase" -> Mod[fixed, 2], "Tableau" -> ps["Tableau"]|>]
]


(* ============================================================================ *)
(* sampleOutcomes[ps, n]: draw n independent samples by substituting each       *)
(* outcome symbol with a uniformly-random 0 or 1.                               *)
(* Returns a list of n PauliStabilizer objects, each with concrete signs.       *)
(* ============================================================================ *)

(* Only sample the FREE symbols (those appearing in the phase that are not   *)
(* derived via the Outcomes map). substituteOutcomes propagates the derived  *)
(* outcomes after random substitution.                                         *)
sampleOutcomes[ps_PauliStabilizer ? PauliStabilizerQ, n_Integer ? Positive] := Module[{
    outcomesMap, derivedSyms, allSyms, freeSyms
},
    outcomesMap = Lookup[First[ps], "Outcomes", <||>];
    derivedSyms = Keys[outcomesMap];
    allSyms = DeleteDuplicates @ Cases[
        Join[ps["Phase"], Flatten[Values[outcomesMap]]],
        _\[FormalS],
        Infinity
    ];
    freeSyms = Complement[allSyms, derivedSyms];
    Table[
        substituteOutcomes[ps, Thread[freeSyms -> RandomInteger[1, Length[freeSyms]]]],
        {n}
    ]
]


sampleOutcomes[ps_PauliStabilizer ? PauliStabilizerQ] := First @ sampleOutcomes[ps, 1]
