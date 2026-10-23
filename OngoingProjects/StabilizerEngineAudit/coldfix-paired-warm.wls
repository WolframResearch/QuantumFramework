PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"];
Needs["Wolfram`QuantumFramework`PackageScope`"];
beforeFold[ps_PauliStabilizer, gates_ /; MatchQ[gates, {{_Integer, _Integer, _Integer} ..}]] :=
    Module[{n = ps["Qubits"], phase, nw, xz, res, codes, controls, targets},
        (* RuntimeOptions -> "Speed" needs a checked boundary before native indexing. *)
        codes = gates[[All, 1]]; controls = gates[[All, 2]]; targets = gates[[All, 3]];
        If[Min[controls, targets] < 0 || Max[controls, targets] >= n ||
            MemberQ[Pick[controls - targets, Unitize[(codes - $stabilizerGateCodes["CNOT"])
                (codes - $stabilizerGateCodes["SWAP"])] , 0], 0],
            Message[PauliStabilizer::invalidorder, n]; Return[$Failed]];
        nw = Ceiling[2 n / $bitsPerWord];
        phase = Developer`ToPackedArray[(1 - ps["Signs"]) / 2];
        xz = psGetGeneratorPacked[ps];
        res = compiledCliffordFold[
            xz[[1]], xz[[2]], First @ packGenRows[{phase}, nw],
            stabilizerValidMasks[n, nw], Developer`ToPackedArray[gates]
        ];
        System`Private`HoldSetValid[PauliStabilizer[<|
            "Signs" -> 1 - 2 First @ unpackGenRows[{res[[-1]]}, 2 n],
            "GeneratorX" -> res[[1 ;; n]], "GeneratorZ" -> res[[n + 1 ;; 2 n]], "Qubits" -> n
        |>]]
    ]


afterFold[ps_PauliStabilizer, gates_ /; ListQ[gates] &&
    MatrixQ[gates, IntegerQ] && MatchQ[Dimensions[gates], {_Integer, 3}]] :=
    Module[{n = ps["Qubits"], phase, nw, xz, res, codes, controls, targets},
        (* RuntimeOptions -> "Speed" needs a checked boundary before native indexing. *)
        codes = gates[[All, 1]]; controls = gates[[All, 2]]; targets = gates[[All, 3]];
        If[Min[controls, targets] < 0 || Max[controls, targets] >= n ||
            MemberQ[Pick[controls - targets, Unitize[(codes - $stabilizerGateCodes["CNOT"])
                (codes - $stabilizerGateCodes["SWAP"])] , 0], 0],
            Message[PauliStabilizer::invalidorder, n]; Return[$Failed]];
        nw = Ceiling[2 n / $bitsPerWord];
        phase = Developer`ToPackedArray[(1 - ps["Signs"]) / 2];
        xz = psGetGeneratorPacked[ps];
        res = compiledCliffordFold[
            xz[[1]], xz[[2]], First @ packGenRows[{phase}, nw],
            stabilizerValidMasks[n, nw], Developer`ToPackedArray[gates]
        ];
        System`Private`HoldSetValid[PauliStabilizer[<|
            "Signs" -> 1 - 2 First @ unpackGenRows[{res[[-1]]}, 2 n],
            "GeneratorX" -> res[[1 ;; n]], "GeneratorZ" -> res[[n + 1 ;; 2 n]], "Qubits" -> n
        |>]]
    ]


records=<||>;
paired[name_,before_,after_]:=Module[{trials,row},
 before[];after[];
 trials=Table[If[OddQ[i],
   With[{a=(ClearSystemCache[];First@AbsoluteTiming[before[]]),b=(ClearSystemCache[];First@AbsoluteTiming[after[]])},{a,b}],
   With[{b=(ClearSystemCache[];First@AbsoluteTiming[after[]]),a=(ClearSystemCache[];First@AbsoluteTiming[before[]])},{a,b}]],{i,21}];
 row=<|"BeforeSeconds"->trials[[All,1]],"AfterSeconds"->trials[[All,2]],
 "BeforeMedian"->Median[trials[[All,1]]],"AfterMedian"->Median[trials[[All,2]]]|>;
 records[name]=row;Echo[name->row];
];
makeSpecs[n_,m_]:=Table[With[{j=Mod[37 i,n]+1},Switch[Mod[i,7],0,"H"->j,1,"S"->j,2,"CNOT"->{j,Mod[j,n]+1},3,"CZ"->{j,Mod[j,n]+1},4,"X"->j,5,"SWAP"->{j,Mod[j,n]+1},6,SuperDagger["S"]->j]],{i,m}];
makeCircuit[s_]:=QuantumCircuitOperator[<|"Elements"->(FromCircuitOperatorShorthand/@s),"Label"->None|>];
PauliStabilizer[2]["ApplyCircuit",{"H"->1}];
Scan[Function[n,With[{p=PauliStabilizer[n],specs=makeSpecs[n,1000]},With[{qc=makeCircuit[specs],long=Flatten[ConstantArray[specs,10],1]},
 paired["PublicWarm1000/n"<>ToString[n],Function[Block[{applyCompiledFold=beforeFold},qc[p]]],Function[Block[{applyCompiledFold=afterFold},qc[p]]]];
 paired["PublicTenCalls/n"<>ToString[n],Function[Block[{applyCompiledFold=beforeFold},Nest[qc[#]&,p,10]]],Function[Block[{applyCompiledFold=afterFold},Nest[qc[#]&,p,10]]]];
 paired["Direct1000/n"<>ToString[n],Function[Block[{applyCompiledFold=beforeFold},p["ApplyCircuit",specs]]],Function[Block[{applyCompiledFold=afterFold},p["ApplyCircuit",specs]]]];
 paired["Direct10000/n"<>ToString[n],Function[Block[{applyCompiledFold=beforeFold},p["ApplyCircuit",long]]],Function[Block[{applyCompiledFold=afterFold},p["ApplyCircuit",long]]]];
]]],{128,1024}];
Export["/private/tmp/stabilizer-coldfix-paired-warm.json",<|"Version"->$Version,"Rows"->records,"Messages"->ToString[$MessageList,InputForm]|>,"RawJSON"];
Echo[$MessageList];Exit[];
