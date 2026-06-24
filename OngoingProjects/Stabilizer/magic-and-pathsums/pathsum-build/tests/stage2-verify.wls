(* Stage 2 verification: diagonal (CNOT-dihedral) phase-polynomial compressor +     *)
(* isolated single-variable Gauss-sum elimination lemmas (Z_2 / Z_4 / odd + quad).   *)
(* Run: wolframscript -file stage2-verify.wls                                        *)

PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"];
dir = DirectoryName[$InputFileName] // ParentDirectory;
Get[FileNameJoin[{dir, "pathsum.wl"}]];
Get[FileNameJoin[{dir, "harness.wl"}]];
Get[FileNameJoin[{dir, "diagonal.wl"}]];

Echo["==== A. diagonal compressor: closed-form amplitude === QF Schrodinger (exact) ====", ""];
(* diagonal-fragment gates (NO Hadamard); compressor prepends the |+...+> H-layer.   *)
d1 = {{"T", 1}, {"S", 2}, {"CZ", 1, 2}, {"CNOT", 1, 2}, {"T", 2}};                  (* S + CNOT + CZ + T *)
d2 = {{"T", 1}, {"S", 2}, {"CCZ", 1, 2, 3}, {"CZ", 1, 3}, {"Z", 2}, {"T", 3}, {"CNOT", 1, 2}}; (* CCZ + S *)
d3 = {{"T", 1}, {"T", 2}, {"T", 3}, {"S", 1}, {"CZ", 1, 2}, {"CZ", 2, 3}};          (* pure diagonal *)
Echo[{PathSum`diagMatchesQ[2, d1], PathSum`diagPolyMatchesQ[2, d1]}, "d1 (S,CNOT,CZ,T) {diagAmp===QF, psi===diagAmp}:"];
Echo[{PathSum`diagMatchesQ[3, d2], PathSum`diagPolyMatchesQ[3, d2]}, "d2 (CCZ,S,CZ,Z,T,CNOT) {diagAmp===QF, psi===diagAmp}:"];
Echo[{PathSum`diagMatchesQ[3, d3], PathSum`diagPolyMatchesQ[3, d3]}, "d3 (pure diagonal) {diagAmp===QF, psi===diagAmp}:"];

Echo["==== A2. rank-1 phase polynomial vs eager StabilizerFrame 2^{#T} ====", ""];
(* d3 has 3 T-gates and no CCZ, so QF builds a StabilizerFrame; its Length = 2^3.    *)
qfFrame = PauliStabilizer[3][{"H" -> 1, "H" -> 2, "H" -> 3, "T" -> 1, "T" -> 2, "T" -> 3,
                              "S" -> 1, "CZ" -> {1, 2}, "CZ" -> {2, 3}}];
c3 = PathSum`diagonalCompress[3, d3];
Echo[{qfFrame["Length"], c3["rank"], c3["bijective"]}, "d3 {eager frame components 2^3, compressor rank, bijective}:"];
Echo[Expand[c3["psi"]], "d3 closed-form phase polynomial psi(y) (one degree-<=3 Z_8 poly):"];

Echo["==== A3. randomized diagonal circuits n<=6 (diagAmp === QF, exact) ====", ""];
SeedRandom[20260622];
diagTrials = Table[Module[{n = RandomInteger[{2, 6}], d}, d = RandomInteger[{4, 14}];
    {n, d, PathSum`diagMatchesQ[n, PathSum`randomDiagCircuit[n, d]]}], {10}];
Echo[diagTrials, "per-circuit {n, depth, diagAmp===QF}:"];
Echo[{Count[diagTrials[[All, 3]], True], Length[diagTrials]}, "passed / total (expect all):"];

Echo["==== B. isolated single-variable Gauss-sum elimination lemmas ====", ""];
(* B1. general closed form reduceVar1 === brute-force single-var sum, all classes.    *)
SeedRandom[7];
randPhase[vars_] := PathSum`z8[RandomInteger[7] + Total[
    RandomInteger[7] Times @@ # & /@ Subsets[vars, {1, 3}]]];
b1 = Table[Module[{ph = randPhase[{v, u, w}], others = {u, w}},
    (* compare reduceVar1 and bruteVar1 at all assignments of the other variables *)
    AllTrue[Tuples[{0, 1}, 2], PathSum`exactEqualQ[
        PathSum`reduceVar1[ph, v] /. Thread[others -> #],
        PathSum`bruteVar1[ph, v] /. Thread[others -> #]] &]], {12}];
Echo[{Count[b1, True], Length[b1]}, "B1 reduceVar1 === bruteVar1 over random phases (expect all):"];

(* B2. Z_2 sign case: phi = 4 v Q (+const), Q linear -> sum over v = 2[Q=0] (a sign). *)
phZ2 = PathSum`z8[4 v (u + w) + 3];
Echo[PathSum`cofactorClass[phZ2, v], "B2 cofactor class of 4 v (u+w)+3 (expect Z2-sign):"];
Echo[Table[{u, w} -> Simplify[PathSum`bruteVar1[phZ2, v] / PathSum`omega^3], {u, 0, 1}, {w, 0, 1}] // Flatten,
    "B2 Sum_v omega^{4v(u+w)} (expect 2 when u(+)w=0, 0 when =1):"];

(* B3. Z_4 / S case: phi = 2 v Q -> 1 + i^Q. This is the A mod 4 != 0 branch elimQ    *)
(* REJECTS, but reduceVar1 handles it (a genuine Gauss factor, not a constraint).     *)
phZ4 = PathSum`z8[2 v u];
Echo[{PathSum`cofactorClass[phZ4, v], PathSum`elimQ[phZ4, v]}, "B3 {class of 2 v u, elimQ (expect Z4-S, False)}:"];
Echo[Table[u -> Simplify[PathSum`reduceVar1[phZ4, v]], {u, 0, 1}],
    "B3 reduceVar1 of 2 v u (expect u=0->2, u=1->1+i):"];
Echo[Table[u -> PathSum`exactEqualQ[PathSum`reduceVar1[phZ4, v], PathSum`bruteVar1[phZ4, v]] /. u -> u, {u, 0, 1}],
    "B3 reduceVar1 === brute (exact) at u=0,1:"];

(* B4. odd-T case: phi = v -> 1 + omega; the magic variable survives (no constraint). *)
Echo[{PathSum`cofactorClass[v, v], Simplify[PathSum`reduceVar1[v, v]]}, "B4 {class of v (a T phase), reduceVar1 = 1+omega}:"];

(* B5. quadratic completion / sequential elimination is SOUND: eliminate v (closed    *)
(* form), then sum the remaining variable u by brute force, === the double sum.        *)
phQ = PathSum`z8[4 u v + 2 v + u + 5];
lhs = (PathSum`reduceVar1[phQ, v] /. u -> 0) + (PathSum`reduceVar1[phQ, v] /. u -> 1);
rhs = Total[PathSum`omega^PathSum`z8[phQ /. Thread[{u, v} -> #]] & /@ Tuples[{0, 1}, 2]];
Echo[PathSum`exactEqualQ[lhs, rhs], "B5 sequential (reduce v, sum u) === double brute sum (expect True):"];

Echo["stage2 done", "STATUS:"];
