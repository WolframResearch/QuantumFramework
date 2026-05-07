Package["Wolfram`QuantumFramework`"]

PackageScope[symbolicMeasure]
PackageScope[substituteOutcomes]
PackageScope[sampleOutcomes]



(* ============================================================================ *)
(* Phase 3 \[Dash] symbolic measurement outcomes (FangYing23 SymPhase).         *)
(*                                                                              *)
(* The standard `ps["M", q]` returns an Association of conditional outcomes,    *)
(* requiring the user to RandomChoice for an actual sample and re-traverse the  *)
(* circuit per sample. SymPhase amortizes this: a single circuit traversal      *)
(* produces a `PauliStabilizer` whose phases contain fresh F_2 symbols          *)
(* representing un-resolved random outcomes. Sampling then substitutes those    *)
(* symbols with random 0/1 values (one matrix-multiply over many samples).      *)
(*                                                                              *)
(* Phase 6 (2026-05-06): demoted from top-level public symbols to method-grade  *)
(* operations on PauliStabilizer:                                               *)
(*   ps["SymbolicMeasure", q]   <-- was symbolicMeasure[ps, q]                *)
(*   ps["SubstituteOutcomes", rules] <-- was substituteOutcomes[ps, rules]      *)
(*   ps["SampleOutcomes", n]    <-- was sampleOutcomes[ps, n]                   *)
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
(* PHASE 3 LIMITATION (TODO: address in Phase 4 via StabilizerFrame):           *)
(* When a deterministic measurement is performed AFTER a prior symbolic         *)
(* measurement, the deterministic outcome polynomial (e.g. m2 = s_1 for Bell    *)
(* ZZ correlation) is computed correctly by the AG algorithm but is NOT        *)
(* explicitly stamped into the post-state's signs. The physics is preserved     *)
(* (the post-state is the right quantum state), but `SampleOutcomes` cannot     *)
(* directly recover the deterministic outcome from the stabilizer signs alone   *)
(* -- the correlation is "implicit" in the unrewritten stabilizer rows.         *)
(* For workflows that need explicit outcome correlation, use `ps["M", q]`       *)
(* (Association branching) until Phase 4 introduces StabilizerFrame.            *)
(* ============================================================================ *)

(* A.4 (2026-05-07): track deterministic-outcome polynomials.                 *)
(*                                                                              *)
(* Previously a deterministic symbolic measurement returned the post-state     *)
(* with no record of the outcome polynomial. Bell ZZ correlation: after        *)
(* SymbolicMeasure on qubit 1 (which allocates \[FormalS][1] for the random   *)
(* outcome), SymbolicMeasure on qubit 2 was deterministic with outcome        *)
(* polynomial \[FormalS][1] (i.e. m_2 = m_1) -- but this polynomial was       *)
(* dropped, so SampleOutcomes drew m_2 independently rather than mirroring   *)
(* m_1.                                                                         *)
(*                                                                              *)
(* Fix: an "Outcomes" key in the PauliStabilizer association records          *)
(* {fresh_symbol -> polynomial_in_prior_symbols} for each deterministic       *)
(* measurement whose outcome involves prior symbols. The non-deterministic    *)
(* case forwards the existing "Outcomes" map. substituteOutcomes and          *)
(* sampleOutcomes apply the Outcomes map iteratively after user / random     *)
(* rules so that derived outcomes resolve correctly.                           *)

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

(* A.4 (2026-05-07): apply user rules first, then iteratively close under the *)
(* "Outcomes" map (deterministic-outcome polynomials) until a fixed point so  *)
(* derived outcomes resolve correctly (e.g. Bell m_2 = m_1).                   *)
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

(* A.4 (2026-05-07): only sample the FREE symbols (those appearing in the    *)
(* phase that are not derived via the Outcomes map). substituteOutcomes      *)
(* propagates the derived outcomes after random substitution.                 *)
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
