#!/usr/bin/env wolframscript
(* Probe 5: corollaries of the materialization bug. Frame single-state observables route
   through fr["StateVector"], so they inherit the per-component global-phase incoherence.
   - <Y> on T|+> from the frame state is sign-flipped vs the dense value (WRONG).
   - The frame inner product (InnerProduct.m:111-116) also uses fr["StateVector"], but
     CAN coincidentally agree with dense when the two frames share the same component
     stabilizer states (hence the same spurious gauge, which cancels in the overlap), as
     for T|+> and S|+> here. So the inner product is not RELIABLY wrong -- it is just
     not trustworthy in general.
   No Quiet, no Print; final expression is the result Association.
   View with:  wolframscript -code 'Print[Get["sf_probe5_corollaries.wls"]]' *)

PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"];

svF = N @ Normal @ PauliStabilizer[1]["H", 1]["T", 1]["StateVector"];                       (* frame T|+> (wrong) *)
svD = N @ Normal @ QuantumCircuitOperator[{"H" -> 1, "T" -> 1}][Method -> "Schrodinger"]["StateVector"];
Ymat = {{0, -I}, {I, 0}};
expY[v_] := Re[Conjugate[v] . Ymat . v];

frA = PauliStabilizer[1]["H", 1]["T", 1];        (* T|+> *)
frB = PauliStabilizer[1]["H", 1]["P"[Pi], 1];    (* P[Pi] = diag(1,e^{iPi/2}) = S, so S|+> *)
ipFrame = Wolfram`QuantumFramework`PackageScope`stabilizerInnerProduct[frA, frB];
denseA = N @ Normal @ QuantumCircuitOperator[{"H" -> 1, "T" -> 1}][Method -> "Schrodinger"]["StateVector"];
denseB = N @ Normal @ QuantumCircuitOperator[{"H" -> 1, "S" -> 1}][Method -> "Schrodinger"]["StateVector"];

<|
    "<Y> dense (correct)"        -> Round[expY[svD], 10.^-6],
    "<Y> from frame state"       -> Round[expY[svF], 10.^-6],
    "<Y> true sin(Pi/4)"         -> Round[Sin[Pi/4], 10.^-6],
    "frame innerproduct <T+|S+>" -> Chop @ N @ ipFrame,
    "dense innerproduct <T+|S+>" -> Chop[Conjugate[denseA] . denseB],
    "innerproduct agree?"        -> (Chop[N @ ipFrame - Conjugate[denseA] . denseB] == 0)
|>
