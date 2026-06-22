(* Experiment 1: SymPhase Facts 1&2, the T-gate disanalogy, the Z_8 impossibility. *)
(* Run: wolframscript -file exp1-facts-and-disanalogy.wls                           *)

PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"];

(* ======================================================================== *)
(* A. Fang-Ying Fact 1 & 2: Pauli gates and measurement act on the SIGN     *)
(*    column of a FIXED X/Z tableau; the control flow ignores the signs.    *)
(* ======================================================================== *)

base = PauliStabilizer[3][{"H" -> 1, "CNOT" -> {1, 2}, "H" -> 3}];

(* Fact 1: X, Y, Z leave the X/Z tableau identical, change only the phase.  *)
Echo[base["Tableau"] === base["X", 1]["Tableau"], "A1 Xgate tableau unchanged:"];
Echo[base["Tableau"] === base["Z", 2]["Tableau"], "A1 Zgate tableau unchanged:"];
Echo[base["Tableau"] === base["Y", 1]["Tableau"], "A1 Ygate tableau unchanged:"];
Echo[base["Phase"] =!= base["X", 1]["Phase"], "A1 Xgate phase DID change:"];

(* A2: symbolicMeasure writes ONE fresh F_2 symbol into ONE sign of a       *)
(*     SINGLE tableau. Both outcome substitutions share that tableau.       *)
(* (Fresh kernel: the running symbol counter starts at 0, so this is s_1.)  *)
sm = PauliStabilizer[2][{"H", "CNOT"}]["SymbolicMeasure", 1];
Echo[Cases[sm["Phase"], _\[FormalS], Infinity], "A2 fresh symbols in phase:"];
b0 = sm["SubstituteOutcomes", {\[FormalS][1] -> 0}];
b1 = sm["SubstituteOutcomes", {\[FormalS][1] -> 1}];
Echo[sm["Tableau"] === b0["Tableau"] === b1["Tableau"], "A2 both branches SHARE one tableau:"];
Echo[b0["Signs"] =!= b1["Signs"], "A2 branches differ only in SIGNS:"];
(* The two branches are the two physical post-measurement states.           *)
Echo[Chop[b0["State"]["StateVector"] // Normal], "A2 branch s1=0 state:"];
Echo[Chop[b1["State"]["StateVector"] // Normal], "A2 branch s1=1 state:"];

(* A3 (Fact 2): control flow is sign-independent. Two states with the SAME  *)
(* tableau but different signs measure to the SAME tableau update.          *)
psA = PauliStabilizer[2][{"H", "CNOT"}];
psB = psA["X", 1];                          (* same tableau, flipped signs *)
Echo[psA["Tableau"] === psB["Tableau"], "A3 inputs share tableau:"];
mA = psA["SymbolicMeasure", 1]; mB = psB["SymbolicMeasure", 1];
Echo[mA["Tableau"] === mB["Tableau"], "A3 post-measure tableaux equal (flow ignores signs):"];

(* ======================================================================== *)
(* B. The T-gate disanalogy. T|+> is a StabilizerFrame of TWO components.   *)
(* ======================================================================== *)

frame = PauliStabilizer[1][{"H", "T"}];      (* T|+> *)
comps = frame["Stabilizers"];
coeffs = frame["Coefficients"];
Echo[frame["Length"], "B1 component count:"];
Echo[coeffs, "B1 component coefficients (COMPLEX, not +/-1):"];

(* SURPRISE vs the naive expectation: the two components SHARE the X/Z      *)
(* tableau and differ only by a relating Pauli (a sign pattern), exactly    *)
(* like a measurement -- because QF decomposes T|s> = a|s> + b Z|s>, and Z  *)
(* leaves the tableau fixed.                                                *)
Echo[comps[[1]]["Tableau"] === comps[[2]]["Tableau"], "B1 the two T-components SHARE one tableau:"];
Echo[comps[[1]]["Signs"] =!= comps[[2]]["Signs"], "B1 ...differing only in SIGNS:"];
Echo[{comps[[1]]["Stabilizers"], comps[[2]]["Stabilizers"]}, "B1 the two stabilizer states (+X vs -X):"];

(* Decisive difference: substituting a sign value yields a COMPONENT, never *)
(* the physical state. |+> and |-> are the components; T|+> is NEITHER.     *)
vPlus  = comps[[1]]["State"]["StateVector"] // Normal // Chop;
vMinus = comps[[2]]["State"]["StateVector"] // Normal // Chop;
vT     = frame["StateVector"] // Normal // Chop;
Echo[vPlus,  "B1 component 1 = |+>:"];
Echo[vMinus, "B1 component 2 = |->:"];
Echo[vT,     "B1 coherent sum c1|+>+c2|-> = T|+>:"];
dense = QuantumCircuitOperator[{"H" -> 1, "T" -> 1}][QuantumState["0"], Method -> "Schrodinger"];
Echo[QuantumSimilarity[frame["State"], dense], "B1 frame === dense T|+> up to global phase:"];
(* T|+> is NOT |+> and NOT |->: the symbolic-sign move cannot reach it.     *)
Echo[{Chop[vT - vPlus] =!= ConstantArray[0, 2], Chop[vT - vMinus] =!= ConstantArray[0, 2]},
  "B1 T|+> differs from BOTH components (no sign assignment reaches it):"];

(* T|+> has stabilizer RANK 2: no single Pauli (with any sign) stabilizes   *)
(* it, so no single tableau can represent it.                               *)
paulis1 = {PauliMatrix[1], PauliMatrix[2], PauliMatrix[3]};
stabilizedBy[v_, P_] := Chop[P . v - v] === ConstantArray[0, Length[v]];
Echo[Or @@ Flatten @ Table[stabilizedBy[vT, s P], {P, paulis1}, {s, {1, -1, I, -I}}],
  "B2 any signed/phased Pauli stabilizes T|+>? (expect False => rank>=2):"];

(* ======================================================================== *)
(* C. Hypothesis 2: a single tableau with a Z_8 sign cannot hold T|+>.      *)
(*    The sign trick works because +/-P is an involution ((+/-P)^2 = I).    *)
(*    A primitive 8th-root phase zeta*P is NOT an involution: it has no     *)
(*    +1 eigenstate at all, so it defines NO stabilizer state.              *)
(* ======================================================================== *)
zeta = Exp[I Pi/4];
genZ8 = zeta PauliMatrix[1];                 (* "e^{i pi/4} X" as a would-be generator *)
Echo[Chop[genZ8 . genZ8], "C zeta*X squared (NOT identity => not an involution):"];
Echo[Chop[Eigenvalues[genZ8]], "C eigenvalues of zeta*X (neither is +1 => no stabilized state):"];
(* Contrast: an ordinary +/- sign keeps the involution and a +1 eigenspace. *)
Echo[Chop[(-PauliMatrix[1]) . (-PauliMatrix[1])], "C (-X)^2 = I (sign stays an involution):"];
Echo[Chop[Eigenvalues[-PauliMatrix[1]]], "C eigenvalues of -X (a +1 eigenstate exists):"];

Echo["exp1 done", "STATUS:"];
