#!/usr/bin/env wolframscript
(* Probe 1: dichotomy + gate-convention check.
   Clifford circuit -> single PauliStabilizer; one non-Clifford gate -> ?
   No Quiet, no Print; final expression is the result Association.
   View with:  wolframscript -code 'Print[Get["sf_probe1_dichotomy.wls"]]'
   NOTE: the RZ/RY cases legitimately emit PauliStabilizer::nonclifford (that IS the
   behavior under test); they are not suppressed.
   PackageScope predicates (StabilizerFrameQ/PauliStabilizerQ) are NOT on the context
   path after Needs, so we classify via Head[...] and the public StabilizerStateQ. *)

PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"];

theta = 0.7;
cliffSpecs = {"H" -> 1, "CNOT" -> {1, 2}, "S" -> 2, "H" -> 3, "CNOT" -> {2, 3}};

classify[res_] := <|
    "Head" -> Head[res],
    "components" -> Which[Head[res] === StabilizerFrame, res["Length"], Head[res] === PauliStabilizer, 1, True, Missing[]],
    "failed" -> (res === $Failed)
|>;

apply[extra_] := classify @ QuantumCircuitOperator[Append[cliffSpecs, extra]][Method -> "Stabilizer"];

rT = QuantumCircuitOperator[Append[cliffSpecs, "T" -> 1]][Method -> "Stabilizer"];

<|
    (* gate-convention check: gate-level P[theta] vs the frame P-handler convention *)
    "P[theta] gate [2,2]"        -> QuantumOperator["P"[theta]]["Matrix"][[2, 2]],
    "e^{i theta}"                -> Exp[I theta],
    "e^{i theta/2} (frame c)"    -> Exp[I theta/2],
    "T gate [2,2]"               -> N @ QuantumOperator["T"]["Matrix"][[2, 2]],

    "Clifford backbone"          -> <|"Head" -> Head[#], "StabilizerStateQ" -> StabilizerStateQ[#]|> &[
                                       QuantumCircuitOperator[cliffSpecs][Method -> "Stabilizer"]],
    "+ T  on q1"                 -> apply["T" -> 1],
    "+ P[theta] on q1"           -> apply["P"[theta] -> 1],
    "+ RZ[theta] on q1"          -> apply["RZ"[theta] -> 1],
    "+ RY[theta] on q1"          -> apply["RY"[theta] -> 1],

    "T-frame component heads"    -> (Head[#[[2]]] & /@ rT["Components"]),
    "T-frame coeffs"             -> rT["Coefficients"],
    "T-frame coeffs == (1+-e^{iPi/4})/2" ->
        (Chop[rT["Coefficients"] - {(1 + Exp[I Pi/4])/2, (1 - Exp[I Pi/4])/2}] // Simplify)
|>
