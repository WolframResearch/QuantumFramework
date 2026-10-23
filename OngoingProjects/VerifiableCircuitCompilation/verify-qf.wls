(* Regime                   | Executable evidence
   General symbolic         | PhaseHalfTurnsConvention, NewtonExactError
   Exactly solvable         | PhaseBranchesPreserveCoherence, AdditionAllBranches
   Limiting                 | NewtonExactInitialGuess, NewtonQuadraticLimit
   Independent reference    | analytic phase diagonal and independently built adder matrix
   Failure / edge           | fixup, outcome, syntax, alias, wire and size checks below
   Independent exact matrices fulfill the reference row without floating-point error. *)
root = DirectoryName[$InputFileName];
repo = Nest[ParentDirectory, root, 2];
dependencyRoot = FileNameJoin[{$UserBaseDirectory, "Paclets", "Repository"}];
Scan[PacletDirectoryLoad, Select[
    FileNames[{"IBMQuantumPlatform-*", "Wolfram__Arrays-*", "Wolfram__TensorNetworks-*", "ImageMetadataTools-*"}, dependencyRoot], DirectoryQ]];
PacletDirectoryLoad[FileNameJoin[{repo, "QuantumFramework"}]];
Needs["Wolfram`QuantumFramework`"];
Get[FileNameJoin[{root, "KickmixBridge.wl"}]];
fragments = Import[FileNameJoin[{root, "results", "fragments.json"}], "RawJSON"];

phase = fragments["AndPhase"];
phaseBranches = KickmixBranchOperator[phase, {#}] & /@ {0, 1};
phaseMatrices = Normal[#["Matrix"]][[{1, 3, 5, 7}, {1, 3, 5, 7}]] & /@ phaseBranches;
idealPhase = DiagonalMatrix[{1, 1, 1, Exp[3 Pi I/16]}];
exactZeroQ[expr_] := TrueQ[FullSimplify[expr] === ConstantArray[0, Dimensions[expr]]];

addition = fragments["Addition"];
additionBranches = KickmixBranchOperator[addition, #] & /@ Tuples[{0, 1}, addition["Measurements"]];
dataIndices = 1 + Range[0, 15] 2^(addition["NumQubits"] - 4);
additionMatrices = Normal[#["Matrix"]][[dataIndices, dataIndices]] & /@ additionBranches;
outputIndex[bits_] := With[{a = bits[[1]] + 2 bits[[2]], b = bits[[3]] + 2 bits[[4]]},
    1 + FromDigits[Join[Reverse[IntegerDigits[Mod[a + b, 4], 2, 2]], bits[[3 ;; 4]]], 2]];
idealAddition = Normal @ SparseArray[MapIndexed[{outputIndex[#1], First[#2]} -> 1 &, Tuples[{0, 1}, 4]], {16, 16}];

Clear[x, e, theta];
newtonRelativeError = FullSimplify[Sqrt[x] ((1 + e)/Sqrt[x]) (3 - x ((1 + e)/Sqrt[x])^2)/2 - 1, x > 0];
symbolicPhase = Normal[QuantumOperator["Phase"[Pi theta]]["Matrix"]];
brokenPhase = Append[phase, "Lines" -> DeleteCases[phase["Lines"], s_String /; StringStartsQ[s, "CZ "]]];
brokenMatrix = Normal[KickmixBranchOperator[brokenPhase, {1}]["Matrix"]][[{1, 3, 5, 7}, {1, 3, 5, 7}]];
invalidInstructions = {"CX q0 q0", "CX q0 q00", "CCX q0 q1 q0", "HMR q0 nonsense",
    "CX nonsense if b0", "HMR q0 if b0", "Z_POW q0 hey if b0"};
aliasBranch = KickmixBranchOperator[<|"NumQubits" -> 1, "Lines" -> {"HMR q0 b00", "X q0 if b0"}|>, {1}];

checks = <|
    "PhaseHalfTurnsConvention" -> exactZeroQ[symbolicPhase - DiagonalMatrix[{1, Exp[I Pi theta]}]],
    "PhaseBranchesPreserveCoherence" -> AllTrue[phaseMatrices, exactZeroQ[Sqrt[2] # - idealPhase] &],
    "PhaseChannelTracePreserving" -> exactZeroQ[Total[ConjugateTranspose[#] . # & /@ phaseMatrices] - IdentityMatrix[4]],
    "PhaseScratchClean" -> AllTrue[phaseBranches, exactZeroQ[Normal[#["Matrix"]][[{2, 4, 6, 8}, {1, 3, 5, 7}]]] &],
    "MissingFixupDetected" -> !exactZeroQ[Sqrt[2] brokenMatrix - idealPhase],
    "AdditionAllBranches" -> AllTrue[additionMatrices, exactZeroQ[2^(addition["Measurements"]/2) # - idealAddition] &],
    "AdditionTracePreserving" -> exactZeroQ[Total[ConjugateTranspose[#] . # & /@ additionMatrices] - IdentityMatrix[16]],
    "NewtonExactError" -> TrueQ[Expand[newtonRelativeError] === -3 e^2/2 - e^3/2],
    "NewtonExactInitialGuess" -> TrueQ[(newtonRelativeError /. e -> 0) === 0],
    "NewtonQuadraticLimit" -> TrueQ[Limit[newtonRelativeError/e^2, e -> 0] === -3/2],
    "MissingMeasurementRejected" -> FailureQ[KickmixBranchOperator[phase, {}]],
    "UnsupportedInstructionRejected" -> FailureQ[KickmixBranchOperator[<|"NumQubits" -> 1, "Lines" -> {"H q0"}|>]],
    "OversizedFragmentRejected" -> FailureQ[KickmixBranchOperator[<|"NumQubits" -> 9, "Lines" -> {}|>]],
    "MalformedInstructionsRejected" -> AllTrue[invalidInstructions,
        FailureQ[KickmixBranchOperator[<|"NumQubits" -> 2, "Lines" -> {#}|>, {1}]] &],
    "ClassicalAliasesCanonicalized" -> exactZeroQ[Normal[aliasBranch["Matrix"]] - {{0, 0}, {1, -1}}/Sqrt[2]],
    "UnicodeOperandsRejected" -> AllTrue[{"Z_POW q0 ٠.٥", "HMR q0 b٢", "X q٢", "X\[NonBreakingSpace]q0"},
        FailureQ[KickmixBranchOperator[<|"NumQubits" -> 2, "Lines" -> {#}|>, {1}]] &],
    "InvalidArgumentShapesRejected" -> And[
        FailureQ[KickmixBranchOperator[phase, 1]],
        FailureQ[KickmixBranchOperator[phase, "0"]],
        FailureQ[KickmixBranchOperator[{}, {}]],
        FailureQ[KickmixBranchOperator[]]]
|>;
report = <|"Kernel" -> $Version, "QFVersion" -> PacletObject["Wolfram/QuantumFramework"]["Version"],
    "Checks" -> checks, "Passed" -> AllTrue[Values[checks], TrueQ]|>;
Export[FileNameJoin[{root, "results", "qf-verification.json"}], report, "RawJSON"];
WriteString[$Output, ExportString[report, "RawJSON"], "\n"];
Exit[If[TrueQ[report["Passed"]], 0, 1]];
