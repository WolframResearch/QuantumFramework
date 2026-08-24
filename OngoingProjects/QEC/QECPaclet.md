# QEC Paclet Development

### Setup

```wolfram
PacletInstall["https://www.wolfr.am/DevWQCF", ForceVersionInstall -> True]
```

```wolfram
Needs["Wolfram`QuantumFramework`"]
```

```wolfram
Get[FileNameJoin[{NotebookDirectory[], "QEC", "QEC.wl"}]]
```

```wolfram
Get[FileNameJoin[{NotebookDirectory[], "QEC", "QECTests.wl"}]]
```

```wolfram
QECRunTests[]
```

### 1. PauliToSymplectic

```wolfram
PauliToSymplectic /@ {"XI", "IX", "ZI", "IZ"} // MatrixForm
```

```wolfram
PauliToSymplectic["XZZXI"]
```

```wolfram
ps = PauliStabilizer[2];
PauliToSymplectic /@ Join[ps["Destabilizers"], ps["Stabilizers"]] === Normal[ps["Matrix"]]
```

### 2. PauliStringSign

```wolfram
withError = PauliStabilizer[2]["X", 1];
{withError["Stabilizers"], PauliStringSign /@ withError["Stabilizers"]}
```

### 3. SymplecticToPauli

```wolfram
SymplecticToPauli[{1, 0, 0, 1, 0, 0, 1, 1, 0, 0}]
```

```wolfram
SymplecticToPauli[{1, 0, 0, 1, 0, 1, 1, 0, 0, 0}]
```

### 4. SymplecticInnerProduct

```wolfram
SymplecticInnerProduct[PauliToSymplectic["X"], PauliToSymplectic["Z"]]
```

```wolfram
SymplecticInnerProduct[PauliToSymplectic["XX"], PauliToSymplectic["ZZ"]]
```

### 5. PauliCommuteQ

```wolfram
{PauliCommuteQ["X", "Z"], PauliCommuteQ["X", "Y"], PauliCommuteQ["X", "X"], PauliCommuteQ["X", "I"]}
```

```wolfram
Table[PauliCommuteQ[StringRepeat["X", n], StringRepeat["Z", n]], {n, 2, 7}]
```

```wolfram
fiveQubitChecks = {"XZZXI", "IXZZX", "XIXZZ", "ZXIXZ"};
Outer[PauliCommuteQ, fiveQubitChecks, fiveQubitChecks] // MatrixForm
```

```wolfram
PauliCommuteQ["IXIII", #] & /@ fiveQubitChecks
```

### 6. QECCode

```wolfram
bitFlip = QECCode[{"ZZI", "IZZ"}]
```

```wolfram
bitFlip["Parameters"]
```

```wolfram
bitFlip["Matrix"] // MatrixForm
```

```wolfram
fiveQubit = QECCode[fiveQubitChecks];
steane = QECCode[{"IIIXXXX", "XIXIXIX", "IXXIIXX", "IIIZZZZ", "ZIZIZIZ", "IZZIIZZ"}];
distanceTwo = QECCode[{StringRepeat["X", 6], StringRepeat["Z", 6]}];
{fiveQubit["Parameters"], steane["Parameters"], distanceTwo["Parameters"]}
```

```wolfram
QECCode[{"XXX", "ZZZ"}]
```

```wolfram
QECCode[{"ZZI", "IZZ", "ZIZ"}]
```

```wolfram
bitFlip["PauliStabilizer"]
```

```wolfram
QECCode[{"ZI", "IZ"}]["State"]["StateVector"] // Normal
```

### 7. Syndrome

```wolfram
Syndrome[bitFlip, #] & /@ {"III", "XII", "IXI", "IIX"}
```

```wolfram
Syndrome[bitFlip, #] & /@ {"ZII", "IZI", "IIZ"}
```

```wolfram
singleErrors = Flatten[Table[StringJoin[ReplacePart[ConstantArray["I", 5], q -> p]], {p, {"X", "Y", "Z"}}, {q, 1, 5}]];
Grid[Prepend[{#, Syndrome[fiveQubit, #]} & /@ singleErrors, {error, syndrome}], Frame -> All]
```

```wolfram
DeleteDuplicates[Syndrome[fiveQubit, #] & /@ singleErrors] // Length
```

```wolfram
Syndrome[fiveQubit, PauliToSymplectic["IXIII"]]
```

### 8. SyndromeTable and PerfectCodeQ

```wolfram
SyndromeTable[bitFlip]
```

```wolfram
SyndromeTable[fiveQubit] // Dataset
```

```wolfram
Values[SyndromeTable[fiveQubit]] // DeleteDuplicates // Length
```

```wolfram
{PerfectCodeQ[fiveQubit], PerfectCodeQ[bitFlip], PerfectCodeQ[steane], PerfectCodeQ[distanceTwo]}
```

### 9. LookupDecoder and DecodeSyndrome

```wolfram
LookupDecoder[bitFlip]
```

```wolfram
LookupDecoder[fiveQubit] // Dataset
```

```wolfram
DecodeSyndrome[fiveQubit, {1, 0, 0, 0}]
```

```wolfram
DecodeSyndrome[fiveQubit, Syndrome[fiveQubit, "IIYII"]]
```

```wolfram
AllTrue[WeightOneErrors[5], DecodeSyndrome[fiveQubit, Syndrome[fiveQubit, #]] === # &]
```

```wolfram
DecodeSyndrome[steane, {1, 1, 0, 0, 0, 1}]
```

### 10. CorrectionCycle

```wolfram
CorrectionCycle[bitFlip, "IXI"]
```

```wolfram
CorrectionCycle[bitFlip, "ZII"]
```

```wolfram
CorrectionCycle[fiveQubit, "IIYII"]
```

```wolfram
AllTrue[WeightOneErrors[5], CorrectionCycle[fiveQubit, #]["Success"] &]
```

```wolfram
CorrectionCycle[fiveQubit, "XXIII"]
```

```wolfram
CorrectionCycle[steane, "IXIIZII"]
```

```wolfram
SeedRandom[7];
Table[CorrectionCycle[fiveQubit]["Success"], {10}]
```

### 11. Framework Bridge: Named Codes and Physical Cycle

```wolfram
{#, QECCode[#]["Parameters"]} & /@ {"BitFlipCode", "PhaseFlipCode", "ShorCode", "5QubitCode", "SteaneCode"} // Grid
```

```wolfram
QECCode["5QubitCode"]["Generators"]
```

```wolfram
QECCode[PauliStabilizer[2]["X", 1]]
```

```wolfram
QECCode[PauliStabilizer[2]["X", 1]]["Signs"]
```

```wolfram
PhysicalCorrectionCycle[QECCode["BitFlipCode"], "IXI"]
```

```wolfram
{CorrectionCycle[QECCode["BitFlipCode"], "ZII"]["Success"], PhysicalCorrectionCycle[QECCode["BitFlipCode"], "ZII"]["Success"]}
```

```wolfram
AllTrue[WeightOneErrors[5], PhysicalCorrectionCycle[QECCode["5QubitCode"], #]["Success"] &]
```

### 12. CodeStandardForm and LogicalOperators

```wolfram
sf = CodeStandardForm[QECCode["5QubitCode"]];
{sf["XRank"], sf["QubitPermutation"]}
```

```wolfram
sf["Matrix"] // MatrixForm
```

```wolfram
LogicalOperators[QECCode["BitFlipCode"]]
```

```wolfram
LogicalOperators[QECCode["5QubitCode"]]
```

```wolfram
LogicalOperators[QECCode["SteaneCode"]]
```

```wolfram
LogicalOperators[QECCode[{StringRepeat["X", 8], StringRepeat["Z", 8]}]]
```

```wolfram
Module[{code = QECCode["5QubitCode"], zbar},
 zbar = LogicalOperators[code]["Z"][[1]];
 InStabilizerGroupQ[code, Mod[PauliToSymplectic[zbar] + PauliToSymplectic["ZZZZZ"], 2]]
]
```

```wolfram
xbar = LogicalOperators[QECCode["5QubitCode"]]["X"][[1]];
{PauliCommuteQ[xbar, #] & /@ QECCode["5QubitCode"]["Generators"], PauliCommuteQ[xbar, "ZZZZZ"]}
```

### 13. CodeDistance

```wolfram
Grid[Prepend[{#, QECCode[#]["FullParameters"]} & /@ {"BitFlipCode", "PhaseFlipCode", "ShorCode", "5QubitCode", "SteaneCode"}, {code, "[[n,k,d]]"}], Frame -> All]
```

```wolfram
CodeDistance[QECCode[{StringRepeat["X", 6], StringRepeat["Z", 6]}]]
```

```wolfram
MinimumWeightLogical[QECCode["5QubitCode"]]
```

```wolfram
{LogicalPauliQ[QECCode["BitFlipCode"], "IIZ"], LogicalPauliQ[QECCode["BitFlipCode"], "ZZI"]}
```

```wolfram
CodeDistance[QECCode[CompleteStabilizerGenerators[QECCode["5QubitCode"]]]]
```

### 14. EncodingCircuit

```wolfram
EncodingCircuit[QECCode["5QubitCode"]]["Diagram"]
```

```wolfram
EncodingGates[QECCode["5QubitCode"]]
```

```wolfram
Grid[Prepend[{#, EncodingCircuitValidQ[QECCode[#]]} & /@ {"BitFlipCode", "PhaseFlipCode", "ShorCode", "5QubitCode", "SteaneCode"}, {code, valid}], Frame -> All]
```

```wolfram
ps = EncodingApplyGates[PauliStabilizer[5], EncodingGates[QECCode["5QubitCode"]]];
ps["Expectation", #] & /@ QECCode["5QubitCode"]["Generators"]
```

### 15. CSSCode

```wolfram
hamming = {{0, 0, 0, 1, 1, 1, 1}, {0, 1, 1, 0, 0, 1, 1}, {1, 0, 1, 0, 1, 0, 1}};
CSSCode[hamming]
```

```wolfram
CSSCode[hamming]["Generators"]
```

```wolfram
Sort[CSSCode[hamming]["Generators"]] === Sort[QECCode["SteaneCode"]["Generators"]]
```

```wolfram
CSSCode[{{1, 1, 1, 1}}]
```

```wolfram
{CSSCodeQ[QECCode["SteaneCode"]], CSSCodeQ[QECCode["5QubitCode"]], CSSCodeQ[QECCode["ShorCode"]]}
```

```wolfram
CSSCode[{{1, 1, 0}}, {{1, 0, 0}}]
```

### 16. ConcatenateCodes

```wolfram
shor = ConcatenateCodes[QECCode["PhaseFlipCode"], QECCode["BitFlipCode"]]
```

```wolfram
shor["Generators"]
```

```wolfram
Sort[shor["Generators"]] === Sort[QECCode["ShorCode"]["Generators"]]
```

```wolfram
shor["FullParameters"]
```

```wolfram
ConcatenateCodes[QECCode[{StringRepeat["X", 4], StringRepeat["Z", 4]}], QECCode["BitFlipCode"]]["Parameters"]
```

### 17. RemoveQubit and PasteCodes

```wolfram
RemoveQubit[QECCode["5QubitCode"]]
```

```wolfram
RemoveQubit[QECCode["5QubitCode"]]["Generators"]
```

```wolfram
Grid[Prepend[{#, QECCode[#]["FullParameters"], RemoveQubit[QECCode[#]]["FullParameters"]} & /@ {"5QubitCode", "SteaneCode", "ShorCode"}, {code, before, "after removal"}], Frame -> All]
```

```wolfram
RemoveQubit[QECCode["BitFlipCode"]]
```

```wolfram
d2 = QECCode[{StringRepeat["X", 4], StringRepeat["Z", 4]}];
pasted = PasteCodes[d2, d2, 1, 1]
```

```wolfram
pasted["Generators"]
```

```wolfram
pasted["FullParameters"]
```

```wolfram
PasteCodes[QECCode["BitFlipCode"], QECCode["SteaneCode"], 1, 1]
```

### 18. Parametric Families and Catalog

```wolfram
RepetitionCode[5]
```

```wolfram
{RepetitionCode[3]["Generators"], PhaseRepetitionCode[3]["Generators"]}
```

```wolfram
Grid[Prepend[{#, DistanceTwoCode[#]["FullParameters"]} & /@ {4, 6, 8}, {n, "[[n,k,d]]"}], Frame -> All]
```

```wolfram
ClassicalHammingMatrix[3] // MatrixForm
```

```wolfram
HammingCSSCode[3]
```

```wolfram
Sort[HammingCSSCode[3]["Generators"]] === Sort[QECCode["SteaneCode"]["Generators"]]
```

```wolfram
Grid[Prepend[{#, HammingCSSCode[#]["Parameters"]} & /@ {3, 4, 5}, {r, "[[n,k]]"}], Frame -> All]
```

```wolfram
QECCodeCatalog[]
```
