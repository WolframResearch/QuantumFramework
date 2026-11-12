# Inventory notes (raw, per file)

Key: S = acts on the user's symbolic data, answer depends on symbol facts; F = already reads facts (simplifyOnReals or Element[..., Reals]); U = runs Simplify-type work the user asked for (should receive facts); N = numeric only (after N, Chop, machine eigen-solvers, export); C = constant or structural (fixed matrices, basis names, dimensions, comments, property lists).

## QuantumOperator/Properties.m (111)
- S 531 "HermitianQ": HermitianMatrixQ on the raw matrix; symbolic entries with real angles test False.
- S 533 "UnitaryQ": UnitaryMatrixQ[N[matrix]]; any symbol makes it False.
- S 744 "Dagger" | "ConjugateTranspose" | "Inverse": the central conjugation.
- S 747 "Conjugate" | "Dual": elementwise conjugate of the state.
- S 752 property built from qo["Conjugate"] (doubling).
- S 793 "Normalize": Abs[qo["Norm"]].
- S 702-709 exact eigenvector orthonormalization: Conjugate[u] . v and Sqrt of it; exactReduce falls back to Simplify.
- S 668, 680 eigenvector orthonormality tests use Conjugate; 680 is the exact (=!=) path.
- S 890-899 UnitaryEulerAngles / with phase: Simplify of Arg, Abs, ArcSin of symbolic entries; Arg of a symbol is meaningless without facts.
- S 1005, 1007, 1141 PossibleZeroQ on symbolic entries and phases (zero tests that facts can decide).
- F 1016-1037 symbolic KAK (simp, real); 1016 also uses PowerExpand, which assumes facts nobody stated.
- U 816 "Simplify" | "FullSimplify" | "ComplexExpand" properties.
- N 572, 621-633, 945-996 (numeric KAK), 1063-1073, 1116, 1118, 1149-1150 (numeric probe), 1195-1258, 1354, 1371, 1408-1419 (QASM export).
- C 26-37, 421, 469, 778-779, 928, 932, 978, 1024, 1054-1058, 1081-1089 (ConjugateTranspose of a fixed basis), 1202, 1284-1311 (computational register states).

## QuantumState/Properties.m (93)
- S 631 "DensityMatrix" of a vector: KroneckerProduct[v, Conjugate[v]]. Every mixed-state quantity of a parametric ket starts here.
- S 1119 "Dual" | "Conjugate": ArrayConjugate of the amplitudes; 1121 "ConjugateTranspose" | "Dagger" built on it. Every bra.
- S 586 "NormalizedQ": Norm == 1 on a symbolic norm.
- S 590, 596 "Normalized" and the association form: Normalize / Norm of symbolic amplitudes (Abs).
- S 290 pure state read off a density matrix: Abs of the diagonal and Sqrt[m[[k, k]]].
- S 292-294 eigen route for a pure state: Chop of symbolic eigenvalues.
- S 653, 676, 684-685 "Eigenvalues", "VonNeumannEntropy": symbolic eigenvalues, Log, Chop @ Simplify.
- S 708 "Purity": Simplify @ Abs of the trace sum.
- S 720, 722 "PureStateQ" / "PhysicalQ": PositiveSemidefiniteMatrixQ[N[...]], False for any symbol.
- S 770-776, 793 "SchmidtBasis" on a ket: SingularValueDecomposition of symbolic coefficients, Sqrt of the norm.
- S 878 "Bipartition": silent Simplify of a whole symbolic state after a change of frame (unasked Simplify).
- S 962, 967 "DecomposeWithProbabilities": Normalize[..., Total] of symbolic weights.
- S 988-993 "Decompose": SingularValueDecomposition of symbolic coefficients.
- S 1072 "Purify": Sqrt of a symbolic state (matrix square root).
- S 1274-1321 "BlochVector", "BlochSphericalCoordinates", global phase: Simplify then Abs, Arg, Re, Im of entries.
- F 312, 324, 326, 542-562 (step 1 and 2 code: Probabilities, Weights, Weight, Norm through simplifyOnReals). To be replaced by the facts function without Simplify.
- U 528-535 "Simplify" | "FullSimplify" | "Chop" | "ComplexExpand" properties.
- N 853 numeric unitarity test; 888-900 exact-numeric simplification; 1131-1148 machine Hermitian part and eigen-pruning.
- C 21, 33, 41, 57, 343-404, 522, 621, 645-649, 693 (MatrixPower, needs no facts), 750, 912, 920-921, 970, 1050, 1100-1101 (basis names).

## Utilities.m (187)
- S 153 projector[v] = KroneckerProduct[v, Conjugate[v]] (pure-state density matrices across the package).
- S 197, 206-207 eigenvector normalization: Normalize (Abs) of symbolic vectors.
- S 219-245 eigensystem wrapper: Simplify @ Enclose[...] on a symbolic matrix (unasked Simplify), Arg and Abs of symbolic eigenvalues in the sorting key.
- S 440-499 matrixFunction: a symbolic matrix goes through its minimal polynomial (generic in parameters); Sqrt, Log, Power of a symbolic matrix pick a branch with no sign facts. Conjugate is passed straight through (463).
- S 984-1000 zeroBaseLimit and 1548-1558 0^m: Piecewise on Re[t] > 0 for symbolic eigenvalues; facts would settle the condition.
- S 1719-1721 QuditAdjacencyMatrix (state graph): Abs and Arg of amplitudes.
- F 1731-1735 simplifyOnReals.
- N 228-292, 268-331, 359-417 (fixed spin and Gell-Mann matrices), 430-432 (exact numeric unitary inverse), 469-1126 numeric matrix-function machinery (Schur, roundoff, branch-cut placement, power bounds), 1197-1210 exactConjugate (numeric entries only), 1231-1483 exact-number zero and sign tests, 1498-1528, 1594.
- C 84, 187, comments.

## QuantumEntanglement.m (125)
Contract today: every symbol of the input is taken real (QuantumEntangledQ::realparameters), passed as ps.
- F 88 renyiOrderQ, 702-733 symbolicContext (Hermiticity and positivity of a symbolic state under the reals), 1107 real-parameter density matrix.
- S 177-180 zero-state test: Simplify[Tr[a]], Conjugate[a] . a on symbolic amplitudes.
- S 271-273 exact Hermiticity test m =!= ConjugateTranspose[m].
- S 338 PossibleZeroQ on symbolic singular points.
- S 584, 593 signs of symbolic purity and partial-transpose eigenvalues (quantitySign).
- S 887-894 assumedSimplify: silent Simplify of every symbolic monotone (2 s cap) under the Assumptions option; 810-821 that option defaults to $Assumptions.
- S 955-975 concurrence: spin flip with Conjugate[mat], Eigenvalues, Sqrt, Re of a symbolic product. Conjugate of a symbolic density matrix is the core use.
- S 1017, 1022, 1038 purity test Simplify[p] == 1; Sqrt and Re of symbolic entries.
- S 1067-1069, 1139, 1173, 1259-1274 symbolic reduced spectra, purity, Renyi entropies: m . ConjugateTranspose[m], Eigenvalues, PossibleZeroQ.
- S 1311 measurement entanglement: Simplify of every post-measurement reduced state (unasked).
- N 163, 214-220, 295-302, 382, 407-434, 480, 578-591, 603 (hermitianPart, used on numbers), 625, 641, 662, 747-776 (exact and numeric sign), 851-871, 946, 1004-1011, 1116-1123, 1182-1248, 1279, 1290.

## QuantumDistance.m (81)
- S 45, 107-160 trace checks of a symbolic state: Re[t], Abs[t - 1], Element[t, Algebraics]; a symbolic trace that is 1 under the facts is reported unnormalized.
- S 54 zeroVectorQ: Conjugate[v] . v.
- S 203-204 physical-state test: HermitianMatrixQ and Re @ Eigenvalues of symbolic matrices.
- S 223-227 "RelativePurity" (Chop @ Tr) and "HilbertSchmidt" (Norm Frobenius = Abs of entries) on symbolic matrices.
- S 313 symbolic fidelity: Re[Total[Sqrt[Eigenvalues[r . s]]]].
- S 485 symbolic trace distance: Re @ Total @ SingularValueList.
- S 502-512 relative entropy: Re, Log, Conjugate[v] . s . v of symbolic eigenvectors.
- N 74, 111, 138-139 (numeric parts), 181-199, 231-480 (machine and exact-numeric fidelity and trace distance), 537-541.

## QuantumOperator/NamedOperators.m (50)
- S 170-238 "RX" | "RY" | "RZ" | "R" built through FullSimplify at construction (skipped only for an inexact angle or a lone symbol with no $Assumptions); records angleRealSymbols (step 1). The source of the gate-angle fact; the construction-time FullSimplify is a cost to remove once facts are applied directly.
- S 78-79, 1116 SuperDagger[...] shorthand goes to "Dagger".
- S 576 operator from a weighted decomposition: Sqrt of symbolic weights.
- S 1012-1022 Lindblad superoperator: Ldg["Dagger"] @ L and Ldg["Conjugate"]; a jump operator Sqrt[gamma] sigma needs gamma >= 0 for Conjugate[Sqrt[gamma]] to be Sqrt[gamma].
- N 746-762 random Hermitian u + u^dagger.
- C 57, 152, 165, 257 (U3 matrix), 272, 307-331, 377, 618-624, 639, 701-706, 892-915, 964-966, 991 (fixed operators), 998-1007, 1096.

## Stabilizer/GateMatrix.m (39)
- S 58, 69 zero tests of symbolic gate entries: Simplify[x] === 0, Chop[Simplify[#]] == 0 (unasked Simplify).
- S 214, 223 phase of a symbolic diagonal gate: Arg[TrigToExp[...]], - I PowerExpand[Log[...]] (PowerExpand assumes a branch nobody stated).
- F 90, 239, 344, 453 stabilizerUnitModulusQ (simplifyOnReals over every symbol of the entry).
- N 37-38, 59, 71, 96-162, 212, 361-383, 451.
- C 199, 517-521, 621 fixed Clifford matrices and phases.

## SecondQuantization/Interferometer.m (131)
The model QF already has: a domain built from the role of each element (opticalDomain: angles real, transmissivity and loss in [0, 1], photon counts nonnegative integers), applied with Refine, plus a cheap conjugation (realConjugate) that flips only explicit Complex numbers when every element parameter is real.
- F 614-681 "Amplitude", "Probability", "Probabilities" refine on opticalDomain and the query's Assumptions; 640 contradictory-assumption check.
- F 859 "Dagger": refineArray[ConjugateTranspose[U], opticalDomain].
- F 1129-1130, 1172, 1234-1385 loss channels and Kraus operators refined on 0 <= eta <= 1; 1244-1280 domainConjugate / realConjugate.
- S 100-109, 273-281, 324-349 decomposition of a symbolic unitary into beam splitters and phases: Abs, Arg, Conjugate of entries (complex by role, correct).
- S 155-169 squeezer: Abs[xi], Arg[xi] for a complex squeezing parameter (correct by role; collapses if xi is known real).
- S 1008-1029 coherent amplitudes: Abs, Conjugate of complex alpha (correct by role).
- S 127-138 identity and unitarity tests: PossibleZeroQ of w . w^dagger - 1.
- N 143-145, 159-199 (Gaussian kernel), 589-593 (plots), 712-713, 727-815 (sampling), 871-878, 939-1000, 1033, 1100, 1148-1151, 1212, 1284, 1319, 1349-1352, 1437.
- C 444 (asimp: Simplify on an exact unitary's phases), 548.

## SecondQuantization/PhaseSpaceRepresentations.m (48)
- S 171, 216, 321 Wigner, Husimi, s-ordered functions of real quadratures: 1/2 ComplexExpand[...] takes every symbol real, the state's own complex parameters (a coherent alpha inside rho) included. A role-correct version expands only x and p.
- S 113-133, 187-193, 328 kernels in alpha and Conjugate[alpha] (complex by role, correct); 149-152 Wirtinger form keeps alpha and alphaStar apart (the model for complex variables).
- S 320, 327, 360 order check Re[s] < 1 on a symbolic s.
- N 52, 69, 92-99, 273-289 (numeric grids).
- C 29, 159, 237, 298, 336-372.

## SecondQuantization/BosonicExp.m (40)
- S 365-388, 413, 537 user Assumptions (default $Assumptions) through Assuming and Simplify in the ordering identities.
- S 427, 434, 437, 447 unasked Simplify inside the algebra (inverse-exponent test, commutators, rotation).
- S 84 PossibleZeroQ[4 p q - w^2] decides the branch of a closed form; facts on p, q, w could settle it.
- C 69-215, 551 closed forms (Sqrt, Log of symbolic coefficients, generic branches); 452-516 expConjugate is conjugation by an exponential, not complex conjugation.

## SecondQuantization/Operators.m (35)
- S 132-178 squeeze operator: Abs[xi], Arg[xi], Conjugate[xi] (complex by role, correct; a real xi needs the fact to collapse).
- S 76-104 displacement: Conjugate[alpha] (complex by role, correct).
- C 29, 194, 241-243 (beam splitter with real theta, phi through SuperDagger of operators), 262-263 quadratures.

## SecondQuantization/Utils.m (26)
- S 175-213 g2 and field correlations: psi["Dagger"] @ ... @ psi of a symbolic state.
- S 241-249 CovarianceMatrix: Re[Conjugate[v] . v], Re[Conjugate[A] . v], norm warning Abs[norm - 1] > 10^-8 on a symbolic norm.
- C 62-118 ladder variable bookkeeping (SuperDagger as a name), 151, 154, 220.

## QuantumCircuitOperator/NamedCircuits.m (25)
- S 225 controlled-phase construction: QuantumOperator[gate]["Dagger"] of a possibly parametric gate.
- S 652-657 state preparation from a symbolic vector: AbsArg[Simplify[...]], Simplify of the angles and of a == b == 0.
- S 839-841 "CHSH"-type circuit: Simplify @ of states and operators carrying theta, at construction (unasked).
- S 535-540 "PowerExpand" option on a controlled-power circuit (user opt-in, assumes branches).
- N 739-807 numeric state preparation (Chop, Normalize, Orthogonalize).
- C 447-456 (Toffoli T^dagger), 529 ("InverseFourier" = Fourier dagger), 813, 842.

## Stabilizer/StabilizerFrame.m (24)
- S 255 Pauli dagger: Conjugate[s] of a coefficient that can be symbolic.
- S 285, 313 expectation values: Conjugate[a[[i]]] b[[j]], Conjugate[v] . P . v.
- S 363, 374 PossibleZeroQ pivots and Simplify @ LinearSolve on symbolic coefficients (unasked).
- C 103-146, 389-438 Clifford conjugation of Paulis (conjugation by a gate, not complex conjugation), fixed matrices.

## QuantumStateEstimation.m (24)
- N all: tomography and Bayesian sampling on numeric data (65 Conjugate of a sampled numeric state).

## QuantumOptimization.m (24)
- F 138, 154, 175 adiabatic gap: Simplify / FullSimplify with s real and the user's assumptions.
- S 156 Abs of a symbolic matrix element without facts (fallback path).
- S 402 quantum natural gradient metric: Re[ConjugateTranspose[d_i] . d_j] - ..., on symbolic derivatives.
- S 410 ComplexExpand // Simplify of a whole parametric circuit matrix: ComplexExpand takes every symbol real, declared or not.
- S 511, 565 expectation Re[ConjugateTranspose[v] . M . v] of symbolic vectors.
- S 736 ancilla["Conjugate"].
- N 179, 432, 792-794, 885-912.
- C 93.

## QuantumMeasurementOperator/NamedMeasurementOperators.m (24)
- S 25, 59 Simplify @ QuantumMeasurementOperator[...] at construction (unasked; fixed POVMs, so cost only).
- S 151-179 eigenvalue grouping of a symbolic observable: Simplify of a partial trace, Simplify[#1 - #2] === 0.
- S 222-238 normality and unitarity tests: m . m^dagger - m^dagger . m, Simplify then PossibleZeroQ.
- N 17 random Hermitian, 177, 193-196.
- C 29, 61-65 fixed SIC vectors.

## SecondQuantization/BosonicAlgebra.m (22)
- S 111, 131 Assumptions option (default $Assumptions) for BosonicVEV, BosonicMatrixElement.
- S 139-166, 201, 300-301 coherent and squeeze matrix elements: Abs, Arg, Conjugate of alpha, xi (complex by role, correct).
- S 230 Simplify[p >= 0, k >= 0] inside a rule (decides 0^p).
- C 135, 141, 153, 163, 240, 297, 336 closed forms.

## QuantumOperator/QuantumOperator.m (21)
- S 870 SuperDagger[qo] and 872 SuperStar[qo]: entry points to "Dagger" and "Conjugate".
- S 653 scalarBasePower: PossibleZeroQ of a symbolic base.
- S 796 exponentialAction: u . rho . ConjugateTranspose[u] of a symbolic exponential.
- C 139-461 basis "Dual" bookkeeping, 267-268, 633-671, 781-783, 893.

## QuantumMeasurement.m (21)
- S 46-51 branch weights: PossibleZeroQ[Im[w]], ComplexExpand[Re[w]] of a symbolic weight (ComplexExpand takes every symbol real).
- S 63-75 branch normalization: Norm, Sqrt of symbolic weights.
- S 248, 256, 291, 316 unasked Simplify of "ProbabilitiesList", distributions, entropy and amplitudes.
- U 338 "Simplify" | "FullSimplify" | "Chop" | "ComplexExpand".
- C 29, 180-226, 359.

## QuantumEvolution.m (20)
- S 40 Heisenberg picture: SuperDagger[evolution] of a symbolic propagator.
- S 108 Lindblad operators scaled by Sqrt[gammas]: rates need gamma >= 0 (a role fact, never recorded).
- S 140 Lindblad right-hand side: ConjugateTranspose of symbolic jump operators.
- F 338 Refine of time-dependent arrays on t0 < t < t1.
- S 391 discontinuity search: ComplexExpand[{Re[#], Im[#]}] takes every symbol real.
- N 118, 233, 413-443, 574, 628-630.
- C 81, 148, 273.

## QuantumChannel/NamedChannels.m (18)
- S 45-88 Kraus operators Sqrt[p], Sqrt[1 - p], Sqrt[gamma], Sqrt[1 - gamma], Sqrt[lambda]: where the [0, 1] facts come from by role; nothing records them, so Conjugate[Sqrt[p]] and Abs[Sqrt[1 - p]]^2 stay unsimplified.

## QuantumBasis/Properties.m (18)
- S 229-318 basis "Dual", "Conjugate", "Dagger", "Bend": conjugate the basis elements; fixed for named bases, symbolic only for a user basis with symbols.
- U 345 "Simplify" | "FullSimplify" | "Chop" | "ComplexExpand" on a basis.
- C 30.

## MIC.m (17)
- C 32-86, 117 Simplify @ of fixed MIC and Wigner bases at construction (cost only, no user symbols); 144-152 fixed projectors.
- N 120-126, 157 random POVMs.
- C 108, 134, 139.

## QuditBasis/NamedBases.m (15)
- C all: fixed bases; 199, 225, 272, 281, 288 Simplify @ GramDual of fixed POVMs (cost only).

## QuantumCircuitOperator/Draw.m (15)
- C all: drawing geometry and labels.

## QuditBasis/Properties.m (13)
- S 197, 210-219 "Conjugate" of basis elements (Conjugate[y], KroneckerProduct[y, Conjugate[y]]); symbolic only for a user basis with symbols.
- U 351-353 "Simplify" | "FullSimplify" | "Chop" | "ComplexExpand".
- C 16, 201-204, 282, 338 "Dual" bookkeeping on names.

## QuditBasis/QuditName.m (12)
- C all: "Dual" toggles ket and bra names.

## QuantumMeasurementOperator/Properties.m (10)
- S 275-280 eigenbasis of a symbolic observable: Simplify @ partial trace, Simplify @ Eigensystem (unasked).
- S 330 Simplify @ QuantumOperator[...] (unasked).
- S 359, 396 "Conjugate" | "Dual" | "Dagger" of a measurement operator.
- N 255. C 12, 136, 369, 393.

## QuantumState/NamedStates.m (9)
- S 238 random/Ginibre-type state m . m^dagger (numeric in practice).
- S 267 "Werner"[p, d]: Simplify of p sym/(...) + (1 - p) as/(...) at construction; p is a probability by role, never recorded.
- S 320-321 QuantumState[SuperDagger[...]] goes to "Dagger".
- C 246, 257, 299-300 ("BlochVector"[r]: r real by role, never recorded), 336.

## Labels.m (9)
Finding: simplifyLabel already applies the named-gate inverse rule to labels, with no check that the angles are real: conj[Subscript["R", ...][a]] -> R[-a], conj["P"[a]] -> "P"[-a], "PhaseShift", SuperDagger["U2"[a, b]] -> "U2"[Pi - b, Pi - a], SuperDagger["U"[a, b, c]] -> "U"[-a, -c, -b], SuperStar["U"[a, b, c]] -> "U"[a, -b, -c], spiders. So RZ[theta]^dagger is labeled RZ[-theta] while its matrix carries Conjugate[theta]: the label is right only for real theta, the matrix only stays general. The same unchecked rule as Cirq and Yao, at the label level.
- S 22-50 those rules; 118 shortcut labels with SuperDagger | SuperStar.

## Experimental.m (9)
- S 19-22 Abs[...]^2 of symbolic amplitudes in the experimental beam-splitter model.
- N 40-42, 237-258 (MPS by SVD, numeric).

## Stabilizer/Constructors.m (8)
- S 107-114 Pauli expectations Conjugate[v] . P . v with Chop and Re thresholds (numeric thresholds on possibly symbolic values).
- S 237 Collect with Simplify.
- N 170, 213, 250.

## QuantumCircuitOperator/Properties.m (8)
- S 388-404 circuit "Dagger", "Conjugate", "Dual", "Double", "Transpose": gate by gate, so each gate's dagger and label rule.
- U 394. N 438, 468.

## Stabilizer/Compiled.m (6)
- C all: gate codes for S^dagger, V^dagger.

## QuantumWignerTransform.m (6)
- C 18-23, 73-75 fixed matrices; 30, 51 Chop @ Simplify @ of fixed transforms and of the transformed state (the latter unasked when the state is symbolic: S).

## QuantumState/QuantumState.m (6)
- S 247 f[qs] for f = Times | Conjugate acts on the amplitudes; any other f goes through matrixFunction (Sqrt, Log, Power of a symbolic density matrix).
- S 304 SuperDagger[qs], 306 SuperStar[qs].
- S 383 ConjugateTranspose in a state product.
- C 93, 317.

## QuantumMeasurementOperator/QuantumMeasurementOperator.m (6)
- S 307-308 Conjugate[qmo]; 114 MatrixPower[#, 1/2] of POVM elements (symbolic square roots).
- C 10, 313, 333.

## QuantumCircuitOperator/Qiskit.m (6)
- N 53-62 export after N. C 65 SuperDagger name passed as a Dagger gate.

## QuantumChannel/Properties.m (6)
- S 68, 84-92 channel "Dagger", "Conjugate", "Dual", "Double" pass to the operator.
- U 85. C 116.

## Multiway.m (6)
- S 20-24 unasked Simplify of the norm and of every branch state; 164 Chop @ Simplify of amplitudes.
- C 159, 182.

## Stabilizer/GateUpdates.m (5), Stabilizer/PauliStabilizer.m (4), Stabilizer/Compiled.m, Stabilizer/PhasePolynomial.m (3)
- C all: S^dagger, V^dagger, T^dagger as Clifford and phase-polynomial updates; tableau "Dagger" | "Inverse" (159) is exact over integers.

## QSVT.m (5)
- S 12-16 block encoding: ConjugateTranspose[A] and MatrixPower[1 - A . A^dagger, 1/2] of a possibly symbolic A.
- C 38, 47.

## Stabilizer/Conversions.m (4)
- S 68-69, 140 operator to Pauli rows: qo @ P @ qo["Dagger"] (exact Clifford, symbolic only if the operator is).
- N 12.

## QuantumPhaseSpaceTransform.m (4)
- S 42 Chop @ Simplify @ of the transformed state (unasked).
- C 25-31.

## QuantumCircuitOperator/QuantumCircuitOperator.m (4)
- S 102 Method -> "Schrodinger": FullSimplify @ after every gate, on any state, symbolic or not. The single largest unasked simplification in the package.
- S 191, 193 SuperDagger[qco], SuperStar[qco].
- C 242.

## QuantumChannel/QuantumChannel.m (4)
- S 99 SuperDagger[qc] -> "Adjoint", 101 Transpose[qc].
- C 50, 108.

## ZX.m (3), Visualization.m (3)
- C all (ZX "Simplify" option is the external pyzx simplifier; plots).

## Stabilizer/InnerProduct.m (3)
- S 105, 121 Conjugate[v] . w of state vectors (symbolic only for symbolic states).
- C 216.

## QuantumBasis/QuantumBasis.m (3)
- C all: "Dual" of the input basis, property names.

## SecondQuantization/States.m (2)
- C 56, 97 Fock coefficients of coherent and cat states (alpha complex by role, phi real by role, neither recorded).

## QuditBasis/QuditBasis.m (2)
- S 177 Conjugate[qb]. C 170.

## QuantumOperator/Formatting.m (2), QuantumMeasurementOperator/Formatting.m (2), QuantumHamiltonianOperator/Formatting.m (2)
- S summary boxes show "HermitianQ" and "UnitaryQ" (1 s cap): a symbolic RZ[theta] is shown as not unitary.

## QuantumHamiltonianOperator/Properties.m (2)
- S 66, 73 gap Abs[Differences[eigenvalues /. parameter -> value]] (numeric after substitution; N in practice).

## QuantumCircuitOperator/OpenQASMImport.m (2)
- S 513, 517 inverse modifiers: op["Dagger"] of an imported parametric gate.

## QuantumHamiltonianOperator/NamedOperators.m (1)
- C 39 // Simplify of a fixed Hamiltonian.

## QuantumCircuitOperator/TensorNetwork.m (1)
- C 271 basis "Conjugate" when the result is a density matrix in a bent basis.

## Gates.m (1)
- C 86 SX^dagger gate.
