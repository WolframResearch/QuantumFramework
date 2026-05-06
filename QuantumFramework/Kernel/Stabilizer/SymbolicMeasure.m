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

symbolicMeasure[ps_PauliStabilizer ? PauliStabilizerQ, a_Integer] := Module[{result},
    result = ps["M", a];
    Switch[Length[result],
        1, (* deterministic *)
        First @ Values @ result,
        2, (* non-deterministic: allocate fresh symbol *)
        Module[{sym = freshOutcome[]},
            (* Pick the post-state for outcome 0 and "stamp" it with the symbol.
               The post-state for outcome b has its `p` row sign flipped from outcome 0.
               The convention: phase[p] = sym, so substituteOutcomes[ps, sym -> 0] gives
               the outcome-0 branch and substituteOutcomes[ps, sym -> 1] gives outcome-1. *)
            With[{ps0 = result[0]},
                Module[{phase0, phase1, p},
                    phase0 = ps0["Phase"];
                    phase1 = result[1]["Phase"];
                    (* p = first index where phase0 differs from phase1 *)
                    p = First @ Position[Mod[phase0 - phase1, 2], 1, {1}, 1];
                    (* substitute symbolic phase at position p *)
                    PauliStabilizer[<|
                        "Phase" -> ReplacePart[phase0, p -> sym],
                        "Tableau" -> ps0["Tableau"]
                    |>]
                ]
            ]
        ],
        _,
        $Failed
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

substituteOutcomes[ps_PauliStabilizer ? PauliStabilizerQ, rules_] := Module[{phase},
    phase = ps["Phase"] /. rules;
    phase = Mod[phase, 2];
    PauliStabilizer[<|"Phase" -> phase, "Tableau" -> ps["Tableau"]|>]
]


(* ============================================================================ *)
(* sampleOutcomes[ps, n]: draw n independent samples by substituting each       *)
(* outcome symbol with a uniformly-random 0 or 1.                               *)
(* Returns a list of n PauliStabilizer objects, each with concrete signs.       *)
(* ============================================================================ *)

sampleOutcomes[ps_PauliStabilizer ? PauliStabilizerQ, n_Integer ? Positive] := Module[{symbols, draws},
    symbols = DeleteDuplicates @ Cases[ps["Phase"], _\[FormalS], Infinity];
    Table[
        Module[{rules},
            rules = Thread[symbols -> RandomInteger[1, Length[symbols]]];
            substituteOutcomes[ps, rules]
        ],
        {n}
    ]
]


sampleOutcomes[ps_PauliStabilizer ? PauliStabilizerQ] := First @ sampleOutcomes[ps, 1]
