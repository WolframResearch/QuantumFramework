(* ::Package:: *)

Package["Wolfram`QuantumFramework`QEC`"]

PackageExport[QECPauliVector]
PackageExport[QECPauliString]
PackageExport[QECPauliWeight]
PackageExport[QECPauliCommuteQ]
PackageExport[QECPauliProduct]
PackageExport[QECPauliPhase]
PackageExport[QECFaultTolerant]
PackageExport[QECStimCircuit]
PackageExport[QECCodeCatalog]
PackageExport[QECClassicalHammingMatrix]


(* ============================================================================ *)
(* The names step 5 of the API redesign retired, kept for one release.          *)
(*                                                                              *)
(* Each still works and returns exactly what it returned before -- it calls     *)
(* the function that now does the job -- and says, once per session, what to   *)
(* write instead.  Once, because a loop over a thousand Paulis should not print *)
(* a thousand warnings; per session, because reloading the package (which every *)
(* test file does) is a new session for this purpose.  Remove the file, and the *)
(* ten exports, one release after the rename.                                   *)
(* ============================================================================ *)

General::qecdeprecated = "`1` is deprecated and will be removed in the next release; use `2` instead.";

$deprecationWarned = <||>;

SetAttributes[deprecated, HoldFirst];
deprecated[old_Symbol, new_String] := If[
    ! KeyExistsQ[$deprecationWarned, SymbolName[Unevaluated[old]]],
    $deprecationWarned[SymbolName[Unevaluated[old]]] = True;
    Message[MessageName[old, "qecdeprecated"], HoldForm[old], new]
]

$deprecatedNames = {
    {QECPauliVector, "QECPauli[p][\"Vector\"]"},
    {QECPauliString, "QECPauli[p][\"String\"]"},
    {QECPauliWeight, "QECPauli[p][\"Weight\"]"},
    {QECPauliCommuteQ, "QECPauli[p][\"CommuteQ\", q]"},
    {QECPauliProduct, "QECPauli[p][\"Product\", q] or QECPauli[p] ** QECPauli[q]"},
    {QECPauliPhase, "QECPauli[p][\"Phase\"]"},
    {QECFaultTolerant, "QECFaultTolerantCircuit"},
    {QECStimCircuit, "QECStim"},
    {QECCodeCatalog, "QECCode[\"Catalog\"]"},
    {QECClassicalHammingMatrix, "QECCode[\"Hamming\", r] for the quantum code, or the matrix itself"}
};

QECPauliVector::usage = "QECPauliVector is deprecated; use QECPauli[p][\"Vector\"].";
QECPauliString::usage = "QECPauliString is deprecated; use QECPauli[p][\"String\"].";
QECPauliWeight::usage = "QECPauliWeight is deprecated; use QECPauli[p][\"Weight\"].";
QECPauliCommuteQ::usage = "QECPauliCommuteQ is deprecated; use QECPauli[p][\"CommuteQ\", q].";
QECPauliProduct::usage = "QECPauliProduct is deprecated; use QECPauli[p][\"Product\", q] or QECPauli[p] ** QECPauli[q].";
QECPauliPhase::usage = "QECPauliPhase is deprecated; use QECPauli[p][\"Phase\"].";
QECFaultTolerant::usage = "QECFaultTolerant is deprecated; use QECFaultTolerantCircuit.";
QECStimCircuit::usage = "QECStimCircuit is deprecated; use QECStim.";
QECCodeCatalog::usage = "QECCodeCatalog is deprecated; use QECCode[\"Catalog\"].";
QECClassicalHammingMatrix::usage = "QECClassicalHammingMatrix is deprecated; QECCode[\"Hamming\", r] builds the quantum code directly.";

QECPauliVector[args___] := (deprecated[QECPauliVector, $deprecatedNames[[1, 2]]]; pauliVector[args])
QECPauliString[args___] := (deprecated[QECPauliString, $deprecatedNames[[2, 2]]]; pauliString[args])
QECPauliWeight[args___] := (deprecated[QECPauliWeight, $deprecatedNames[[3, 2]]]; pauliWeight[args])
QECPauliCommuteQ[args___] := (deprecated[QECPauliCommuteQ, $deprecatedNames[[4, 2]]]; pauliCommuteQ[args])
QECPauliProduct[args___] := (deprecated[QECPauliProduct, $deprecatedNames[[5, 2]]]; pauliProduct[args])
QECPauliPhase[args___] := (deprecated[QECPauliPhase, $deprecatedNames[[6, 2]]]; pauliPhase[args])
QECFaultTolerant[args___] := (deprecated[QECFaultTolerant, $deprecatedNames[[7, 2]]]; QECFaultTolerantCircuit[args])
QECStimCircuit[args___] := (deprecated[QECStimCircuit, $deprecatedNames[[8, 2]]]; QECStim[args])
QECCodeCatalog[args___] := (deprecated[QECCodeCatalog, $deprecatedNames[[9, 2]]]; codeCatalog[args])
QECClassicalHammingMatrix[args___] := (deprecated[QECClassicalHammingMatrix, $deprecatedNames[[10, 2]]]; classicalHammingMatrix[args])
