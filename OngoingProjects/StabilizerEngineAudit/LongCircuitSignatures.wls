repo = Replace[Environment["STABILIZER_AUDIT_REPO"], $Failed -> "/Users/mohammadb/Documents/GitHub/QuantumFramework"];
phase = Environment["STABILIZER_AUDIT_PHASE"];
PacletDirectoryLoad[FileNameJoin[{repo,"QuantumFramework"}]];
Needs["Wolfram`QuantumFramework`"];
Needs["Wolfram`QuantumFramework`PackageScope`"];
makeSpecs[n_, m_] := Table[With[{j = Mod[37 i, n] + 1},
    Switch[Mod[i, 7], 0, "H" -> j, 1, "S" -> j,
        2, "CNOT" -> {j, Mod[j, n] + 1}, 3, "CZ" -> {j, Mod[j, n] + 1},
        4, "X" -> j, 5, "SWAP" -> {j, Mod[j, n] + 1}, 6, SuperDagger["S"] -> j]], {i, m}];

signature[ps_] := IntegerString[Hash[{ps["Signs"],ps["Tableau"]},"SHA256"],16,64];
rows=Map[Function[n,With[{ps=PauliStabilizer[n],base=makeSpecs[n,1000]},
    With[{specs=Flatten[ConstantArray[base,10],1]},
        With[{bulk=ps["ApplyCircuit",specs],split=Fold[#1["ApplyCircuit",#2]&,ps,ConstantArray[base,10]]},
            <|"Qubits"->n,"Gates"->Length[specs],"BulkHash"->signature[bulk],"SplitHash"->signature[split],
                "Equal"->(bulk["Signs"]===split["Signs"] && bulk["Tableau"]===split["Tableau"])|>]]]],{128,1024}];
rows=Append[rows, With[{ps=PauliStabilizer[128],specs=Insert[makeSpecs[128,10000],"V"->1,5000]},
    <|"Qubits"->128,"Gates"->10001,"OneVHash"->signature[ps["ApplyCircuit",specs]]|>]];
result=<|"Phase"->phase,"Version"->$Version,"LoadedFile"->FindFile["Wolfram`QuantumFramework`"],
    "Rows"->rows,"Messages"->ToString[$MessageList,InputForm]|>;
Echo[result];
Export["/private/tmp/stabilizer-long-signatures-"<>phase<>".json",result,"RawJSON"];
Exit[];
