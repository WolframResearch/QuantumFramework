(* ::Package:: *)

Package["Wolfram`QuantumFramework`QEC`"]

PackageExport[QECClearCache]

PackageScope[$QECMemoisedFunctions]


(* ============================================================================ *)
(* Memoisation, and how to undo it.                                             *)
(*                                                                              *)
(* Derived properties are memoised on the code's data association: the standard  *)
(* form, the logical operators, the distance, the decoder table and the encoder  *)
(* are each computed at most once per code.  Since a QECCode is an immutable      *)
(* expression, that cache lives in the DownValues of the internal functions,      *)
(* keyed by the association.                                                     *)
(*                                                                              *)
(* Which is fine in a session and a trap while developing the package: reloading  *)
(* a file adds the new general rule but leaves the old cached answers in place,   *)
(* and a literal-argument rule always outranks a pattern rule.  The stale answer  *)
(* then wins silently and the reload looks like it did nothing.  QECClearCache[]  *)
(* drops the cached answers and keeps the definitions, so a reload takes effect.  *)
(* ============================================================================ *)

QECClearCache::usage = "QECClearCache[] clears the memoised results of the QEC core (standard forms, logical operators, distances, decoder tables, encoders). Definitions are left intact. Call it after reloading the package during development, or to release memory held by cached results.";

(* Every function that caches its result on the code's data association. *)
$QECMemoisedFunctions := {
    codeBasis, codeSpanData, codeCompletedGenerators, codeStandardForm,
    codeLogicalVectors, codeMinimumLogical, codeSwapped, codeSyndromeTable,
    codeDecoderToWeight, codeDecoder, codeEncodingGates, codeCircuitInstructions,
    codeLabelMatrix, codeDetectorModel, demRowKeys, demGroups, demDecoderTable, codeErrorTally, codeMaximumLikelihoodDecoder, codeCosetRepresentatives, cosetProbabilities, codeMinimumWeightDecoder
};

(* A cached answer is a rule whose left-hand side holds no pattern: its arguments
   are literal.  Every defining rule carries one, so this keeps them all. *)
cachedRuleQ[rule_] := FreeQ[First[rule], Pattern | Blank | BlankSequence | BlankNullSequence]

QECClearCache[] := Total[
    Function[f,
        With[{kept = DeleteCases[DownValues[Evaluate[f]], _ ? cachedRuleQ]},
            With[{dropped = Length[DownValues[Evaluate[f]]] - Length[kept]},
                DownValues[Evaluate[f]] = kept;
                dropped
            ]
        ]
    ] /@ $QECMemoisedFunctions
]
