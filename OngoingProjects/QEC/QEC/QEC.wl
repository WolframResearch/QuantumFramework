(* ::Package:: *)

Package["Wolfram`QuantumFramework`QEC`"]

PackageExport["PauliToSymplectic"]
PackageExport["PauliStringSign"]
PackageExport["SymplecticToPauli"]
PackageExport["SymplecticInnerProduct"]
PackageExport["PauliCommuteQ"]
PackageExport["PauliWeight"]
PackageExport["QECCode"]
PackageExport["QECCodeCatalog"]
PackageExport["Syndrome"]
PackageExport["SyndromeTable"]
PackageExport["WeightOneErrors"]
PackageExport["WeightKPaulis"]
PackageExport["PerfectCodeQ"]
PackageExport["LookupDecoder"]
PackageExport["DecodeSyndrome"]
PackageExport["InStabilizerGroupQ"]
PackageExport["LogicalPauliQ"]
PackageExport["CorrectionCycle"]
PackageExport["PhysicalCorrectionCycle"]
PackageExport["ApplyPauliString"]
PackageExport["CompleteStabilizerGenerators"]
PackageExport["CodeStandardForm"]
PackageExport["LogicalOperators"]
PackageExport["CodeDistance"]
PackageExport["MinimumWeightLogical"]
PackageExport["EncodingGates"]
PackageExport["EncodingCircuit"]
PackageExport["EncodingCircuitValidQ"]
PackageExport["EncodingApplyGates"]
PackageExport["CSSCode"]
PackageExport["CSSCodeQ"]
PackageExport["ConcatenateCodes"]
PackageExport["RemoveQubit"]
PackageExport["PasteCodes"]
PackageExport["RepetitionCode"]
PackageExport["PhaseRepetitionCode"]
PackageExport["DistanceTwoCode"]
PackageExport["ClassicalHammingMatrix"]
PackageExport["HammingCSSCode"]

PackageScope["rowToPauli"]
PackageScope["pauliStringProduct"]
PackageScope["placeOnBlock"]
PackageScope["translateOuterGenerator"]
PackageScope["EncodingReductionGates"]
PackageScope["removeQubitAt"]
PackageScope["invertGate"]


PauliToSymplectic::usage =
"PauliToSymplectic[str] converts a Pauli string such as \"XZZXI\" into its binary symplectic vector (x1..xn|z1..zn).";

PauliStringSign::usage =
"PauliStringSign[str] returns -1 if the Pauli string carries a leading minus sign and +1 otherwise.";

SymplecticToPauli::usage =
"SymplecticToPauli[vec] converts a binary symplectic vector of length 2n back into a Pauli string.";

SymplecticInnerProduct::usage =
"SymplecticInnerProduct[u, v] gives the symplectic product x.z' + z.x' mod 2 of two symplectic vectors.";

PauliCommuteQ::usage =
"PauliCommuteQ[p, q] gives True if the two Paulis commute. Arguments may be Pauli strings or symplectic vectors.";

PauliWeight::usage =
"PauliWeight[str] gives the number of qubits on which the Pauli string acts nontrivially.";

QECCode::usage =
"QECCode[{gen1, gen2, ...}] represents a stabilizer code with the given Pauli string generators.\nQECCode[name] builds a named code, such as \"SteaneCode\" or \"5QubitCode\".\nQECCode[ps] converts a Wolfram`QuantumFramework`PauliStabilizer into a code.";

QECCodeCatalog::usage =
"QECCodeCatalog[] gives a Dataset of the available named codes and code families with their parameters.";

Syndrome::usage =
"Syndrome[code, error] gives the error syndrome as a vector of bits, one per stabilizer generator.";

SyndromeTable::usage =
"SyndromeTable[code] gives an Association from every weight-one Pauli error to its syndrome.";

WeightOneErrors::usage =
"WeightOneErrors[n] gives all 3n Pauli strings of weight one on n qubits.";

WeightKPaulis::usage =
"WeightKPaulis[n, k] gives all Pauli strings of weight exactly k on n qubits.";

PerfectCodeQ::usage =
"PerfectCodeQ[code] gives True if the code saturates the quantum Hamming bound for single-error correction.";

LookupDecoder::usage =
"LookupDecoder[code] gives an Association mapping each syndrome to a minimum-weight correction.";

DecodeSyndrome::usage =
"DecodeSyndrome[code, syndrome] gives the correction for the given syndrome, or Missing if no weight-one error produces it.";

InStabilizerGroupQ::usage =
"InStabilizerGroupQ[code, vec] gives True if the symplectic vector lies in the stabilizer group of the code.";

LogicalPauliQ::usage =
"LogicalPauliQ[code, pauli] gives True if the Pauli is a logical operator, i.e. lies in N(S) but not in S.";

CorrectionCycle::usage =
"CorrectionCycle[code, error] runs the symplectic error-correction cycle and reports syndrome, correction, residual and success.\nCorrectionCycle[code] uses a random weight-one error.";

PhysicalCorrectionCycle::usage =
"PhysicalCorrectionCycle[code, error] runs the cycle on a fiducial encoded state through Wolfram`QuantumFramework`PauliStabilizer, comparing stabilizers with signs.";

ApplyPauliString::usage =
"ApplyPauliString[ps, str] applies a Pauli string to a Wolfram`QuantumFramework`PauliStabilizer state as a sequence of gates.";

CompleteStabilizerGenerators::usage =
"CompleteStabilizerGenerators[code] extends the code generators to n commuting independent generators, pinning a fiducial encoded state.";

CodeStandardForm::usage =
"CodeStandardForm[code] gives the standard form of the code matrix together with the qubit permutation used and the rank of the X block.";

LogicalOperators::usage =
"LogicalOperators[code] gives an Association with the logical X and Z operators of the code as Pauli strings.";

CodeDistance::usage =
"CodeDistance[code] gives the code distance, the minimum weight of a logical operator.";

MinimumWeightLogical::usage =
"MinimumWeightLogical[code] gives a logical operator of minimum weight, the witness of the code distance.";

EncodingGates::usage =
"EncodingGates[code] gives the list of Clifford gates preparing a codeword of the code.";

EncodingCircuit::usage =
"EncodingCircuit[code] gives a Wolfram`QuantumFramework`QuantumCircuitOperator preparing a codeword of the code.";

EncodingCircuitValidQ::usage =
"EncodingCircuitValidQ[code] gives True if the encoding circuit prepares a state stabilized by every generator.";

EncodingApplyGates::usage =
"EncodingApplyGates[ps, gates] applies a list of gate rules to a Wolfram`QuantumFramework`PauliStabilizer state.";

CSSCode::usage =
"CSSCode[hx, hz] builds a CSS code from two classical parity-check matrices.\nCSSCode[h] uses the same matrix for both X-type and Z-type stabilizers.";

CSSCodeQ::usage =
"CSSCodeQ[code] gives True if every generator is purely X-type or purely Z-type.";

ConcatenateCodes::usage =
"ConcatenateCodes[outer, inner] concatenates two codes, re-encoding each qubit of the outer code with the inner code.";

RemoveQubit::usage =
"RemoveQubit[code] removes one qubit, turning an [[n,k,d]] code into an [[n-1,k+1,d-1]] code.\nRemoveQubit[code, q] removes the specified qubit.";

PasteCodes::usage =
"PasteCodes[code1, code2, solo1, solo2] pastes two codes together, keeping the first solo1 and solo2 generators separate and pairing the rest.";

RepetitionCode::usage =
"RepetitionCode[n] gives the n-qubit repetition code protecting against bit-flip errors.";

PhaseRepetitionCode::usage =
"PhaseRepetitionCode[n] gives the n-qubit repetition code protecting against phase-flip errors.";

DistanceTwoCode::usage =
"DistanceTwoCode[n] gives the [[n, n-2, 2]] code for even n.";

ClassicalHammingMatrix::usage =
"ClassicalHammingMatrix[r] gives the parity-check matrix of the classical Hamming code with r checks.";

HammingCSSCode::usage =
"HammingCSSCode[r] gives the CSS code built from the classical Hamming code, reducing to the Steane code for r = 3.";


$pauliLetterToXZ = <|"I" -> {0, 0}, "X" -> {1, 0}, "Y" -> {1, 1}, "Z" -> {0, 1}|>;
$xzToPauliLetter = <|{0, 0} -> "I", {1, 0} -> "X", {1, 1} -> "Y", {0, 1} -> "Z"|>;

PauliToSymplectic::invchar = "String `1` contains invalid characters: only I, X, Y, Z and an optional leading + or - are allowed.";

PauliToSymplectic[s_String] := Module[{str, chars, pairs},
  str = StringDelete[s, StartOfString ~~ ("+" | "-")];
  chars = Characters[str];
  If[chars === {} || ! SubsetQ[{"I", "X", "Y", "Z"}, DeleteDuplicates[chars]],
    Message[PauliToSymplectic::invchar, s]; Return[$Failed]
  ];
  pairs = $pauliLetterToXZ /@ chars;
  Join[pairs[[All, 1]], pairs[[All, 2]]]
]

PauliStringSign[s_String] := If[StringStartsQ[s, "-"], -1, 1]

SymplecticToPauli::badvec = "Argument must be a binary vector of even length.";

SymplecticToPauli[v_List] := Module[{n, x, z},
  If[OddQ[Length[v]] || ! SubsetQ[{0, 1}, DeleteDuplicates[v]],
    Message[SymplecticToPauli::badvec]; Return[$Failed]
  ];
  n = Length[v]/2;
  x = Take[v, n];
  z = Drop[v, n];
  StringJoin[$xzToPauliLetter /@ Transpose[{x, z}]]
]

SymplecticInnerProduct::len = "Arguments must be binary vectors of the same even length.";

SymplecticInnerProduct[u_List, v_List] := Module[{n},
  If[Length[u] =!= Length[v] || OddQ[Length[u]],
    Message[SymplecticInnerProduct::len]; Return[$Failed]
  ];
  n = Length[u]/2;
  Mod[Take[u, n] . Drop[v, n] + Drop[u, n] . Take[v, n], 2]
]

PauliCommuteQ[p_String, q_String] := SymplecticInnerProduct[PauliToSymplectic[p], PauliToSymplectic[q]] === 0

PauliCommuteQ[u_List, v_List] := SymplecticInnerProduct[u, v] === 0

QECCode::invgen = "Generators must be valid Pauli strings of equal length.";
QECCode::noncomm = "Generators `1` and `2` do not commute.";
QECCode::dep = "Generators are not independent: the symplectic matrix has rank `1` but there are `2` generators.";

QECCode[gens : {__String}] := Module[{vecs, n, m, rank, pairs, badPair},
  vecs = Quiet[PauliToSymplectic /@ gens];
  If[MemberQ[vecs, $Failed] || Length[DeleteDuplicates[Length /@ vecs]] =!= 1,
    Message[QECCode::invgen]; Return[$Failed]
  ];
  m = Length[gens];
  n = Length[First[vecs]]/2;
  pairs = Subsets[Range[m], {2}];
  badPair = SelectFirst[pairs, ! PauliCommuteQ[vecs[[#[[1]]]], vecs[[#[[2]]]]] &];
  If[! MissingQ[badPair],
    Message[QECCode::noncomm, gens[[badPair[[1]]]], gens[[badPair[[2]]]]]; Return[$Failed]
  ];
  rank = MatrixRank[vecs, Modulus -> 2];
  If[rank =!= m,
    Message[QECCode::dep, rank, m]; Return[$Failed]
  ];
  QECCode[<|
    "Generators" -> gens,
    "Signs" -> (PauliStringSign /@ gens),
    "Matrix" -> vecs,
    "Qubits" -> n,
    "StabilizerCount" -> m,
    "LogicalQubits" -> n - m
  |>]
]

QECCode[a_Association]["Parameters"] := {a["Qubits"], a["LogicalQubits"]}

QECCode[a_Association]["PauliStabilizer"] := Wolfram`QuantumFramework`PauliStabilizer[a["Generators"]]

QECCode[a_Association]["State"] := Wolfram`QuantumFramework`PauliStabilizer[a["Generators"]]["State"]

QECCode[a_Association]["Properties"] := Join[Keys[a], {"Parameters", "PauliStabilizer", "State"}]

QECCode[a_Association][key_String] := Lookup[a, key, Missing["NotFound", key]]

(* Summary box, like Wolfram`QuantumFramework`QuantumState / Wolfram`QuantumFramework`PauliStabilizer. Collapsed rows show cheap data (n, k,
   generator count); expanded rows compute the heavy properties (distance, logical operators)
   lazily via Dynamic, so they run only when the user opens the box, not on every display. *)
QECCode /: MakeBoxes[obj : QECCode[a_Association] /; KeyExistsQ[a, "Qubits"], form : (StandardForm | TraditionalForm)] :=
  BoxForm`ArrangeSummaryBox[
    QECCode,
    obj,
    ArrayPlot[a["Matrix"], Mesh -> All, MeshStyle -> GrayLevel[0.8],
      ColorRules -> {0 -> White, 1 -> RGBColor[0.15, 0.5, 0.65]},
      ImageSize -> {Automatic, 34}, Frame -> False],
    {
      BoxForm`SummaryItem[{"Qubits: ", a["Qubits"]}],
      BoxForm`SummaryItem[{"Logical qubits: ", a["LogicalQubits"]}],
      BoxForm`SummaryItem[{"Generators: ", a["StabilizerCount"]}]
    },
    {
      BoxForm`SummaryItem[{"Stabilizers: ", Row[a["Generators"], ", "]}],
      BoxForm`SummaryItem[{"Signs: ", a["Signs"]}],
      BoxForm`SummaryItem[{"Distance: ", Dynamic[CodeDistance[obj]]}],
      BoxForm`SummaryItem[{"Parameters: ", Dynamic[Row[{"[[", Row[obj["FullParameters"], ","], "]]"}]]}],
      BoxForm`SummaryItem[{"Logical X: ", Dynamic[LogicalOperators[obj]["X"]]}],
      BoxForm`SummaryItem[{"Logical Z: ", Dynamic[LogicalOperators[obj]["Z"]]}]
    },
    form,
    "Interpretable" -> False
  ]

Syndrome::len = "Error must act on `1` qubits.";

Syndrome[QECCode[a_Association], err_String] := Module[{v},
  v = PauliToSymplectic[err];
  If[v === $Failed, Return[$Failed]];
  Syndrome[QECCode[a], v]
]

Syndrome[QECCode[a_Association], v_List] := Module[{},
  If[Length[v] =!= 2 a["Qubits"],
    Message[Syndrome::len, a["Qubits"]]; Return[$Failed]
  ];
  SymplecticInnerProduct[v, #] & /@ a["Matrix"]
]

WeightOneErrors[n_Integer?Positive] := Flatten[Table[StringJoin[ReplacePart[ConstantArray["I", n], q -> p]], {q, n}, {p, {"X", "Y", "Z"}}]]

SyndromeTable[QECCode[a_Association]] := AssociationMap[Syndrome[QECCode[a], #] &, WeightOneErrors[a["Qubits"]]]

(* PerfectCodeQ tests perfection for single-error-correcting codes only (t = 1, Gottesman thesis sec. 8.4):
   it checks the quantum Hamming equality 1 + 3n == 2^(n-k) and that all weight-1 syndromes are distinct and nonzero.
   A general version parametrized by t would count errors up to weight t; add it once CodeDistance exists. *)
PerfectCodeQ[QECCode[a_Association]] := Module[{n, m, syn},
  n = a["Qubits"];
  m = a["StabilizerCount"];
  If[1 + 3 n =!= 2^m, Return[False]];
  syn = Values[SyndromeTable[QECCode[a]]];
  Length[DeleteDuplicates[syn]] === 3 n && FreeQ[syn, ConstantArray[0, m]]
]

LookupDecoder[QECCode[a_Association]] := Module[{n, decoder},
  n = a["Qubits"];
  decoder = <|ConstantArray[0, a["StabilizerCount"]] -> StringRepeat["I", n]|>;
  Do[
    With[{syn = Syndrome[QECCode[a], err]},
      If[! KeyExistsQ[decoder, syn], decoder[syn] = err]
    ],
    {err, WeightOneErrors[n]}
  ];
  decoder
]

DecodeSyndrome::badsyn = "Syndrome must be a binary vector of length `1`.";

DecodeSyndrome[QECCode[a_Association], syn_List] := Module[{},
  If[Length[syn] =!= a["StabilizerCount"] || ! SubsetQ[{0, 1}, DeleteDuplicates[syn]],
    Message[DecodeSyndrome::badsyn, a["StabilizerCount"]]; Return[$Failed]
  ];
  Lookup[LookupDecoder[QECCode[a]], Key[syn], Missing["UnknownSyndrome", syn]]
]

InStabilizerGroupQ[QECCode[a_Association], v_List] := MatrixRank[Append[a["Matrix"], v], Modulus -> 2] === a["StabilizerCount"]

(* CorrectionCycle verifies success up to phases: the residual error*correction is accepted if its
   symplectic vector lies in the row space of the stabilizer matrix over GF(2). Sign bookkeeping
   for the residual will be added together with a phase-tracking PauliMultiply. *)
CorrectionCycle[QECCode[a_Association], err_String] := Module[{code, errv, syn, corr, corrv, residual, success},
  code = QECCode[a];
  errv = PauliToSymplectic[err];
  If[errv === $Failed, Return[$Failed]];
  syn = Syndrome[code, errv];
  corr = DecodeSyndrome[code, syn];
  If[MissingQ[corr],
    Return[<|"Error" -> err, "Syndrome" -> syn, "Correction" -> corr, "Residual" -> Missing["UnknownSyndrome"], "Success" -> False|>]
  ];
  corrv = PauliToSymplectic[corr];
  residual = Mod[errv + corrv, 2];
  success = InStabilizerGroupQ[code, residual];
  <|"Error" -> err, "Syndrome" -> syn, "Correction" -> corr, "Residual" -> SymplecticToPauli[residual], "Success" -> success|>
]

CorrectionCycle[QECCode[a_Association]] := CorrectionCycle[QECCode[a], RandomChoice[WeightOneErrors[a["Qubits"]]]]

QECCode::unknown = "Unknown named code `1`. Available names: `2`.";

$QECNamedCodeRules := <|
  "BitFlipCode" -> {"ZZI", "IZZ"},
  "PhaseFlipCode" -> {"XXI", "IXX"},
  "ShorCode" -> {"ZZIIIIIII", "IZZIIIIII", "IIIZZIIII", "IIIIZZIII", "IIIIIIZZI", "IIIIIIIZZ", "XXXXXXIII", "IIIXXXXXX"},
  "5QubitCode" -> Most[Wolfram`QuantumFramework`PauliStabilizer["5QubitCode"]["Stabilizers"]],
  "SteaneCode" -> Most[Wolfram`QuantumFramework`PauliStabilizer["SteaneCode"]["Stabilizers"]]
|>

QECCode[name_String] := Module[{rules = $QECNamedCodeRules},
  If[! KeyExistsQ[rules, name],
    Message[QECCode::unknown, name, StringRiffle[Keys[rules], ", "]]; Return[$Failed]
  ];
  QECCode[rules[name]]
]

QECCode[ps_Wolfram`QuantumFramework`PauliStabilizer] := QECCode[ps["Stabilizers"]]

ApplyPauliString[ps_Wolfram`QuantumFramework`PauliStabilizer, str_String] := Module[{s},
  s = StringDelete[str, StartOfString ~~ ("+" | "-")];
  Fold[If[#2[[1]] === "I", #1, #1[#2[[1]], #2[[2]]]] &, ps, Transpose[{Characters[s], Range[StringLength[s]]}]]
]

(* CompleteStabilizerGenerators extends the m code generators to n commuting independent ones,
   pinning a fiducial encoded state (a specific logical basis state). Extra rows are picked greedily
   from the symplectic null space (the commutant of the code). Section 12's LogicalOperators will
   provide a canonical choice; any completion works as a reference state. *)
CompleteStabilizerGenerators::incomplete = "Could not complete the generator set to a full stabilizer state.";

CompleteStabilizerGenerators[QECCode[a_Association]] := Module[{mat, n, zero, J, comm, rows, extras},
  mat = a["Matrix"];
  n = a["Qubits"];
  zero = ConstantArray[0, {n, n}];
  J = ArrayFlatten[{{zero, IdentityMatrix[n]}, {IdentityMatrix[n], zero}}];
  comm = NullSpace[Mod[mat . J, 2], Modulus -> 2];
  rows = mat;
  extras = {};
  Do[
    If[Length[rows] < n &&
       MatrixRank[Append[rows, c], Modulus -> 2] > Length[rows] &&
       AllTrue[extras, SymplecticInnerProduct[#, c] === 0 &],
      rows = Append[rows, c]; extras = Append[extras, c]
    ],
    {c, comm}
  ];
  If[Length[rows] =!= n,
    Message[CompleteStabilizerGenerators::incomplete]; Return[$Failed]
  ];
  Join[a["Generators"], SymplecticToPauli /@ extras]
]

(* PhysicalCorrectionCycle runs the cycle on the fiducial encoded state (completed generators),
   applying error and correction as actual gates and comparing final stabilizers, signs included.
   Semantics differ from the symplectic CorrectionCycle: this checks restoration of THAT state, so a
   residual logical Z (which fixes the fiducial state) counts as success here, not as failure. *)
PhysicalCorrectionCycle[QECCode[a_Association], err_String] := Module[{code, ref, damaged, syn, corr, final, success},
  code = QECCode[a];
  ref = Wolfram`QuantumFramework`PauliStabilizer[CompleteStabilizerGenerators[code]];
  damaged = ApplyPauliString[ref, err];
  syn = (1 - (damaged["Expectation", #] & /@ a["Generators"]))/2;
  corr = DecodeSyndrome[code, syn];
  If[MissingQ[corr],
    Return[<|"Error" -> err, "Syndrome" -> syn, "Correction" -> corr, "Success" -> False|>]
  ];
  final = ApplyPauliString[damaged, corr];
  success = final["Stabilizers"] === ref["Stabilizers"];
  <|"Error" -> err, "Syndrome" -> syn, "Correction" -> corr, "Success" -> success|>
]

(* CodeStandardForm implements Gottesman thesis sec. 4.1: Gaussian elimination over GF(2) with row
   operations (generator products) and qubit column swaps, bringing the code matrix to
   [I A1 A2 | B C1 C2 ; 0 0 0 | D I E]. Returns the reduced matrix, the qubit permutation used,
   and r, the rank of the X block. *)
CodeStandardForm[QECCode[a_Association]] := Module[{mat, n, m, perm, r, s2, piv, i, c, j},
  mat = a["Matrix"];
  n = a["Qubits"];
  m = a["StabilizerCount"];
  perm = Range[n];
  r = 0;
  While[r < m,
    piv = Catch[
      Do[If[mat[[ii, cc]] === 1, Throw[{ii, cc}]], {cc, r + 1, n}, {ii, r + 1, m}];
      Throw[Missing["NoPivot"]]
    ];
    If[MissingQ[piv], Break[]];
    {i, c} = piv;
    If[c =!= r + 1,
      mat[[All, {r + 1, c}]] = mat[[All, {c, r + 1}]];
      mat[[All, {n + r + 1, n + c}]] = mat[[All, {n + c, n + r + 1}]];
      perm[[{r + 1, c}]] = perm[[{c, r + 1}]]
    ];
    If[i =!= r + 1, mat[[{i, r + 1}]] = mat[[{r + 1, i}]]];
    Do[If[j =!= r + 1 && mat[[j, r + 1]] === 1, mat[[j]] = Mod[mat[[j]] + mat[[r + 1]], 2]], {j, m}];
    r++
  ];
  s2 = 0;
  While[r + s2 < m,
    piv = Catch[
      Do[If[mat[[ii, n + cc]] === 1, Throw[{ii, cc}]], {cc, r + s2 + 1, n}, {ii, r + s2 + 1, m}];
      Throw[Missing["NoPivot"]]
    ];
    If[MissingQ[piv], Break[]];
    {i, c} = piv;
    If[c =!= r + s2 + 1,
      mat[[All, {r + s2 + 1, c}]] = mat[[All, {c, r + s2 + 1}]];
      mat[[All, {n + r + s2 + 1, n + c}]] = mat[[All, {n + c, n + r + s2 + 1}]];
      perm[[{r + s2 + 1, c}]] = perm[[{c, r + s2 + 1}]]
    ];
    If[i =!= r + s2 + 1, mat[[{i, r + s2 + 1}]] = mat[[{r + s2 + 1, i}]]];
    Do[If[j =!= r + s2 + 1 && j > r && mat[[j, n + r + s2 + 1]] === 1, mat[[j]] = Mod[mat[[j]] + mat[[r + s2 + 1]], 2]], {j, m}];
    s2++
  ];
  <|"Matrix" -> mat, "QubitPermutation" -> perm, "XRank" -> r|>
]

(* LogicalOperators implements the closed formulas of thesis sec. 4.1 on the standard form:
   Xbar = (0 E^T I | E^T C1^T + C2^T 0 0), Zbar = (0 0 0 | A2^T 0 I), then maps the result
   back to the original qubit ordering. *)
LogicalOperators[QECCode[a_Association]] := Module[{n, m, k, sf, mat, perm, r, mid, a2, c1, c2, eMat, unpermute, xRows, zRows},
  n = a["Qubits"];
  m = a["StabilizerCount"];
  k = n - m;
  If[k === 0, Return[<|"X" -> {}, "Z" -> {}|>]];
  sf = CodeStandardForm[QECCode[a]];
  mat = sf["Matrix"];
  perm = sf["QubitPermutation"];
  r = sf["XRank"];
  mid = m - r;
  a2 = mat[[1 ;; r, r + mid + 1 ;; n]];
  c1 = mat[[1 ;; r, n + r + 1 ;; n + r + mid]];
  c2 = mat[[1 ;; r, n + r + mid + 1 ;; 2 n]];
  eMat = mat[[r + 1 ;; m, n + r + mid + 1 ;; 2 n]];
  unpermute[v_] := Module[{x = Take[v, n], z = Drop[v, n], xo, zo},
    xo = ConstantArray[0, n];
    zo = ConstantArray[0, n];
    Do[xo[[perm[[s]]]] = x[[s]]; zo[[perm[[s]]]] = z[[s]], {s, n}];
    Join[xo, zo]
  ];
  xRows = Table[
    Join[
      ConstantArray[0, r],
      Table[eMat[[l, i]], {l, mid}],
      IdentityMatrix[k][[i]],
      Table[Mod[Sum[eMat[[l, i]] c1[[j, l]], {l, mid}] + c2[[j, i]], 2], {j, r}],
      ConstantArray[0, mid],
      ConstantArray[0, k]
    ], {i, k}];
  zRows = Table[
    Join[
      ConstantArray[0, n],
      Table[a2[[j, i]], {j, r}],
      ConstantArray[0, mid],
      IdentityMatrix[k][[i]]
    ], {i, k}];
  <|"X" -> (SymplecticToPauli[unpermute[#]] & /@ xRows),
    "Z" -> (SymplecticToPauli[unpermute[#]] & /@ zRows)|>
]

PauliWeight[s_String] := StringLength[StringDelete[s, StartOfString ~~ ("+" | "-")]] - StringCount[s, "I"]

WeightKPaulis[n_Integer?Positive, k_Integer] := If[k == 0, {StringRepeat["I", n]},
  Flatten[Table[
    StringJoin[ReplacePart[ConstantArray["I", n], Thread[pos -> #]]] & /@ Tuples[{"X", "Y", "Z"}, k],
    {pos, Subsets[Range[n], {k}]}
  ]]
]

LogicalPauliQ[QECCode[a_Association], pauli_String] := Module[{v},
  v = PauliToSymplectic[pauli];
  If[v === $Failed || Length[v] =!= 2 a["Qubits"], Return[False]];
  Syndrome[QECCode[a], v] === ConstantArray[0, a["StabilizerCount"]] && ! InStabilizerGroupQ[QECCode[a], v]
]

(* CodeDistance returns the minimum weight of an element of N(S) - S (a logical operator),
   Gottesman thesis sec. 3.2. Brute-force search over increasing weight; feasible for small n.
   Returns Infinity for codes with no logical qubits (k = 0). *)
CodeDistance[QECCode[a_Association]] := Module[{n = a["Qubits"], found},
  If[a["LogicalQubits"] === 0, Return[Infinity]];
  Do[
    found = SelectFirst[WeightKPaulis[n, w], LogicalPauliQ[QECCode[a], #] &];
    If[! MissingQ[found], Return[w, Module]],
    {w, 1, n}
  ];
  Infinity
]

MinimumWeightLogical[QECCode[a_Association]] := Module[{n = a["Qubits"], found},
  If[a["LogicalQubits"] === 0, Return[Missing["NoLogicalOperators"]]];
  Do[
    found = SelectFirst[WeightKPaulis[n, w], LogicalPauliQ[QECCode[a], #] &];
    If[! MissingQ[found], Return[found, Module]],
    {w, 1, n}
  ];
  Missing["NoLogicalOperators"]
]

QECCode[a_Association]["Distance"] := CodeDistance[QECCode[a]]

QECCode[a_Association]["FullParameters"] := {a["Qubits"], a["LogicalQubits"], CodeDistance[QECCode[a]]}

(* EncodingReductionGates implements Gottesman QECCbook sec. 6.4.1 (Procedure 6.6) with Table 6.1:
   the m generators are placed as columns of a 2n x m symplectic matrix G (rows 1..n = x-part = A,
   rows n+1..2n = z-part = C). Left-multiplication gates act as row operations (H swaps A/C rows,
   S adds an A row to its C row, CNOT/CZ combine rows, SWAP swaps qubits); generator relabeling is a
   free column operation. Reducing G to the initial Z-stabilizer form records the reduction gate list;
   the encoder is its inverse. *)
EncodingReductionGates[QECCode[a_Association]] := Module[
  {n, m, g, gates, hh, ss, cnot, cz, swap, coladd, piv, i},
  n = a["Qubits"];
  m = a["StabilizerCount"];
  g = Transpose[a["Matrix"]];
  gates = {};
  hh[q_] := (g[[{q, n + q}]] = g[[{n + q, q}]]; AppendTo[gates, "H" -> q]);
  ss[q_] := (g[[n + q]] = Mod[g[[n + q]] + g[[q]], 2]; AppendTo[gates, "S" -> q]);
  cnot[c_, t_] := (g[[t]] = Mod[g[[t]] + g[[c]], 2]; g[[n + c]] = Mod[g[[n + c]] + g[[n + t]], 2]; AppendTo[gates, "CNOT" -> {c, t}]);
  cz[c_, t_] := (g[[n + t]] = Mod[g[[n + t]] + g[[c]], 2]; g[[n + c]] = Mod[g[[n + c]] + g[[t]], 2]; AppendTo[gates, "CZ" -> {c, t}]);
  swap[u_, v_] := (g[[{u, v}]] = g[[{v, u}]]; g[[{n + u, n + v}]] = g[[{n + v, n + u}]]; AppendTo[gates, "SWAP" -> {u, v}]);
  coladd[src_, dst_] := (g[[All, dst]] = Mod[g[[All, dst]] + g[[All, src]], 2]);
  Do[
    piv = SelectFirst[Range[p, n], g[[#, p]] === 1 &];
    If[MissingQ[piv],
      piv = SelectFirst[Range[p, n], g[[n + #, p]] === 1 &];
      hh[piv]
    ];
    If[piv =!= p, swap[p, piv]];
    Do[If[i =!= p && g[[i, p]] === 1, cnot[p, i]], {i, n}];
    If[g[[n + p, p]] === 1, ss[p]];
    Do[If[i =!= p && g[[n + i, p]] === 1, cz[p, i]], {i, n}];
    Do[If[q =!= p && g[[p, q]] === 1, coladd[p, q]], {q, m}],
    {p, m}
  ];
  Do[hh[q], {q, m}];
  gates
]

(* invertGate reverses a reduction gate; H/CNOT/CZ/SWAP are self-inverse, S^-1 = S^3. *)
invertGate["H" -> q_] := {"H" -> q}
invertGate["CNOT" -> qs_] := {"CNOT" -> qs}
invertGate["CZ" -> qs_] := {"CZ" -> qs}
invertGate["SWAP" -> qs_] := {"SWAP" -> qs}
invertGate["S" -> q_] := {"S" -> q, "S" -> q, "S" -> q}

EncodingApplyGates[ps_Wolfram`QuantumFramework`PauliStabilizer, gates_List] := Fold[
  Which[
    MatchQ[#2, ("H" | "S" | "X" | "Y" | "Z") -> _], #1[#2[[1]], #2[[2]]],
    True, #1[#2[[1]], Sequence @@ #2[[2]]]
  ] &,
  ps, gates
]

(* EncodingGates returns the encoder: the reversed, inverted reduction sequence, with a Pauli sign
   fixup prepended. The fixup solves, over GF(2), for a Pauli whose symplectic products with the n
   completed generators match the pattern of wrong (-1) signs, then applies it as X/Z gates so the
   prepared state becomes the +1 codeword of every generator. *)
EncodingGates[QECCode[a_Association]] := Module[
  {red, enc, n = a["Qubits"], ps, full, fmat, sys, desired, p, fixGates},
  red = EncodingReductionGates[QECCode[a]];
  enc = Flatten[invertGate /@ Reverse[red], 1];
  ps = EncodingApplyGates[Wolfram`QuantumFramework`PauliStabilizer[n], enc];
  full = CompleteStabilizerGenerators[QECCode[a]];
  fmat = PauliToSymplectic /@ full;
  desired = Table[If[ps["Expectation", full[[i]]] === -1, 1, 0], {i, Length[full]}];
  If[Total[desired] === 0, Return[enc]];
  sys = Join[#[[n + 1 ;; 2 n]], #[[1 ;; n]]] & /@ fmat;
  p = LinearSolve[sys, desired, Modulus -> 2];
  fixGates = Flatten[Table[
    {If[p[[q]] === 1, "X" -> q, Nothing], If[p[[n + q]] === 1, "Z" -> q, Nothing]},
    {q, n}
  ]];
  Join[enc, fixGates]
]

EncodingCircuit[QECCode[a_Association]] := Wolfram`QuantumFramework`QuantumCircuitOperator[EncodingGates[QECCode[a]]]

EncodingCircuitValidQ[QECCode[a_Association]] := Module[{ps},
  ps = EncodingApplyGates[Wolfram`QuantumFramework`PauliStabilizer[a["Qubits"]], EncodingGates[QECCode[a]]];
  AllTrue[a["Generators"], ps["Expectation", #] === 1 &]
]

(* CSSCode builds a Calderbank-Shor-Steane code (Gottesman thesis sec. 8.6, QECCbook ch. 5) from two
   classical binary parity-check matrices: rows of hx become X-type stabilizers, rows of hz become
   Z-type stabilizers. They commute iff hx.hz^T = 0 (mod 2). CSSCode[h] uses hx = hz = h (weakly
   self-dual case, e.g. the [7,4] Hamming code gives the Steane code). *)
CSSCode::orthogonal = "HX and HZ must satisfy HX.Transpose[HZ] = 0 (mod 2) so X-type and Z-type stabilizers commute.";
CSSCode::dims = "HX and HZ must have the same number of columns (qubits).";

rowToPauli[row_List, letter_String] := StringJoin[If[# === 1, letter, "I"] & /@ row]

CSSCode[hx_?MatrixQ, hz_?MatrixQ] := Module[{xstabs, zstabs},
  If[Dimensions[hx][[2]] =!= Dimensions[hz][[2]],
    Message[CSSCode::dims]; Return[$Failed]
  ];
  If[! AllTrue[Flatten[Mod[hx . Transpose[hz], 2]], # === 0 &],
    Message[CSSCode::orthogonal]; Return[$Failed]
  ];
  xstabs = rowToPauli[#, "X"] & /@ hx;
  zstabs = rowToPauli[#, "Z"] & /@ hz;
  QECCode[Join[xstabs, zstabs]]
]

CSSCode[h_?MatrixQ] := CSSCode[h, h]

CSSCodeQ[QECCode[a_Association]] := AllTrue[a["Generators"],
  StringFreeQ[StringDelete[#, StartOfString ~~ ("+" | "-")], "Y"] &&
  (StringFreeQ[#, "X"] || StringFreeQ[#, "Z"]) &]

pauliStringProduct[p_String, q_String] := SymplecticToPauli[Mod[PauliToSymplectic[p] + PauliToSymplectic[q], 2]]

placeOnBlock[check_String, b_Integer, n2_Integer, bigN_Integer] :=
  StringJoin[StringRepeat["I", (b - 1) n2], check, StringRepeat["I", bigN - b n2]]

translateOuterGenerator[g_String, xbar_String, zbar_String, n2_Integer] :=
  StringJoin[Switch[#,
    "I", StringRepeat["I", n2],
    "X", xbar,
    "Z", zbar,
    "Y", pauliStringProduct[xbar, zbar]
  ] & /@ Characters[StringDelete[g, StartOfString ~~ ("+" | "-")]]]

(* ConcatenateCodes builds the concatenated code (Gottesman thesis sec. 3.5): each physical qubit of
   the outer [[n1,k,d1]] code is re-encoded with the inner [[n2,1,d2]] code, giving [[n1 n2, k, d1 d2]].
   Inner-floor generators are the inner checks copied on each block; outer-floor generators are the
   outer checks with each X/Z replaced by the inner logical Xbar/Zbar on the corresponding block. *)
ConcatenateCodes::inner = "The inner code must encode exactly one logical qubit (k = 1).";

ConcatenateCodes[QECCode[outer_Association], QECCode[inner_Association]] := Module[
  {n1, n2, bigN, lo, xbar, zbar, innerGens, outerGens},
  If[inner["LogicalQubits"] =!= 1, Message[ConcatenateCodes::inner]; Return[$Failed]];
  n1 = outer["Qubits"];
  n2 = inner["Qubits"];
  bigN = n1 n2;
  lo = LogicalOperators[QECCode[inner]];
  xbar = lo["X"][[1]];
  zbar = lo["Z"][[1]];
  innerGens = Flatten[Table[placeOnBlock[c, b, n2, bigN], {b, n1}, {c, inner["Generators"]}]];
  outerGens = translateOuterGenerator[#, xbar, zbar, n2] & /@ outer["Generators"];
  QECCode[Join[innerGens, outerGens]]
]

(* RemoveQubit implements the qubit-removal surgery of Gottesman thesis sec. 3.5, turning an
   [[n,k,d]] code into an [[n-1,k+1,d-1]] one. Generators are recombined (row operations) so that
   exactly one ends in X and one ends in Z at the removed qubit, with all others ending in I; the
   two are dropped (they become the new logical pair) and the remaining rows are truncated. *)
RemoveQubit::nopivot = "No qubit admits the removal surgery: some qubit must carry both X-type and Z-type generator tails.";

removeQubitAt[a_Association, q_Integer] := Module[{mat, n, m, xrow, zrow, keep, cols},
  n = a["Qubits"];
  m = a["StabilizerCount"];
  mat = a["Matrix"];
  xrow = SelectFirst[Range[m], mat[[#, q]] === 1 &];
  If[MissingQ[xrow], Return[$Failed]];
  Do[If[i =!= xrow && mat[[i, q]] === 1, mat[[i]] = Mod[mat[[i]] + mat[[xrow]], 2]], {i, m}];
  zrow = SelectFirst[Range[m], # =!= xrow && mat[[#, n + q]] === 1 &];
  If[MissingQ[zrow], Return[$Failed]];
  Do[If[i =!= zrow && mat[[i, n + q]] === 1, mat[[i]] = Mod[mat[[i]] + mat[[zrow]], 2]], {i, m}];
  keep = Complement[Range[m], {xrow, zrow}];
  cols = Join[Complement[Range[n], {q}], Complement[Range[n + 1, 2 n], {n + q}]];
  QECCode[SymplecticToPauli /@ mat[[keep, cols]]]
]

RemoveQubit[QECCode[a_Association], q_Integer] := Module[{res},
  res = removeQubitAt[a, q];
  If[res === $Failed, Message[RemoveQubit::nopivot]; Return[$Failed]];
  res
]

RemoveQubit[QECCode[a_Association]] := Module[{res},
  res = SelectFirst[Reverse[Range[a["Qubits"]]], removeQubitAt[a, #] =!= $Failed &];
  If[MissingQ[res], Message[RemoveQubit::nopivot]; Return[$Failed]];
  removeQubitAt[a, res]
]

(* PasteCodes implements the pasting construction of Gottesman thesis sec. 3.5. The first solo1
   generators of code1 and solo2 of code2 act alone (padded with identity on the other block); the
   remaining generators are paired across blocks, one fused generator per pair, which is what buys
   the extra logical qubits. Requires the same number of leftover generators on both sides. *)
PasteCodes::pairing = "Both codes must leave the same number of generators to pair: `1` vs `2`.";
PasteCodes::solo = "Solo counts must be between 0 and the number of generators of each code.";

PasteCodes[QECCode[a1_Association], QECCode[a2_Association], solo1_Integer, solo2_Integer] := Module[
  {m1, m2, n1, n2, pad1, pad2, soloGens, pairedGens},
  m1 = a1["StabilizerCount"]; m2 = a2["StabilizerCount"];
  n1 = a1["Qubits"]; n2 = a2["Qubits"];
  If[! (0 <= solo1 <= m1 && 0 <= solo2 <= m2), Message[PasteCodes::solo]; Return[$Failed]];
  If[m1 - solo1 =!= m2 - solo2, Message[PasteCodes::pairing, m1 - solo1, m2 - solo2]; Return[$Failed]];
  pad1 = StringRepeat["I", n2];
  pad2 = StringRepeat["I", n1];
  soloGens = Join[
    (# <> pad1) & /@ Take[a1["Generators"], solo1],
    (pad2 <> #) & /@ Take[a2["Generators"], solo2]
  ];
  pairedGens = MapThread[StringJoin, {Drop[a1["Generators"], solo1], Drop[a2["Generators"], solo2]}];
  QECCode[Join[soloGens, pairedGens]]
]

(* Parametric code families. RepetitionCode / PhaseRepetitionCode are the n-qubit repetition codes
   (Gottesman thesis sec. 2.2); DistanceTwoCode is the [[n,n-2,2]] family of thesis sec. 8.1, which
   saturates the quantum Singleton bound; HammingCSSCode[r] is the CSS code built from the classical
   [2^r-1, 2^r-1-r, 3] Hamming code (QECCbook sec. 5.1), giving [[2^r-1, 2^r-1-2r, 3]] and reducing
   to the Steane code at r = 3. *)
RepetitionCode::small = "The repetition code needs at least 2 qubits.";
DistanceTwoCode::even = "The [[n,n-2,2]] family needs an even number of qubits.";
HammingCSSCode::small = "HammingCSSCode needs r >= 3 for the classical Hamming code to be weakly self-dual.";

RepetitionCode[n_Integer] := Module[{},
  If[n < 2, Message[RepetitionCode::small]; Return[$Failed]];
  QECCode[Table[StringJoin[ReplacePart[ConstantArray["I", n], {i -> "Z", i + 1 -> "Z"}]], {i, n - 1}]]
]

PhaseRepetitionCode[n_Integer] := Module[{},
  If[n < 2, Message[RepetitionCode::small]; Return[$Failed]];
  QECCode[Table[StringJoin[ReplacePart[ConstantArray["I", n], {i -> "X", i + 1 -> "X"}]], {i, n - 1}]]
]

DistanceTwoCode[n_Integer] := Module[{},
  If[! EvenQ[n] || n < 4, Message[DistanceTwoCode::even]; Return[$Failed]];
  QECCode[{StringRepeat["X", n], StringRepeat["Z", n]}]
]

ClassicalHammingMatrix[r_Integer] := Transpose[Table[IntegerDigits[j, 2, r], {j, 2^r - 1}]]

HammingCSSCode[r_Integer] := Module[{},
  If[r < 3, Message[HammingCSSCode::small]; Return[$Failed]];
  CSSCode[ClassicalHammingMatrix[r]]
]

QECCodeCatalog[] := Dataset[Association @@ Join[
  Table[name -> Association[Thread[{"n", "k"} -> QECCode[name]["Parameters"]]], {name, Keys[$QECNamedCodeRules]}],
  {
    "RepetitionCode[5]" -> Association[Thread[{"n", "k"} -> RepetitionCode[5]["Parameters"]]],
    "DistanceTwoCode[6]" -> Association[Thread[{"n", "k"} -> DistanceTwoCode[6]["Parameters"]]],
    "HammingCSSCode[4]" -> Association[Thread[{"n", "k"} -> HammingCSSCode[4]["Parameters"]]]
  }
]]
