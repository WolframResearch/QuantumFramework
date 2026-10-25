here = DirectoryName[$InputFileName];
Get[FileNameJoin[{here,"LoaderIdentities.wl"}]];
build = Import[FileNameJoin[{here,"build.json"}],"RawJSON"];
wordFile = FileNameJoin[{here,"local","words.u16"}];
If[IntegerString[FileHash[wordFile,"SHA256"],16,64] =!= build["words_sha256"], Exit[1]];
words = Partition[BinaryReadList[wordFile,"UnsignedInteger16",ByteOrdering->-1],75];
If[Dimensions[words] =!= {3976,75}, Exit[1]];
checks = DFTHCC30`SymbolicChecks[];
batches = Transpose[Partition[Transpose[words],3],{1,3,2}];
deltas = Prepend[BitXor[Rest[batches],Most[batches]],First[batches]];
forward = Rest[FoldList[BitXor,ConstantArray[0,{3976,3}],deltas]];
dataChecks = <|"all_four_address_rows" -> (Catenate[DFTHCC30`FourAddress /@ Partition[words,4]] === words),
  "all_delta_signed_bits" -> (BitGet[deltas,9] === BitGet[deltas,8]),
  "all_forward_batches" -> (forward === batches),
  "all_reverse_batches_clear" -> (Fold[BitXor,Last[batches],Reverse[deltas]] === ConstantArray[0,{3976,3}]),
  "all_angle_widths_128" -> AllTrue[build["profiles"], #["num_qubits"] === 128 &],
  "fixed_factors" -> TrueQ[build["weights_and_offsets_exactly_released"]]|>;
If[!AllTrue[Join[Values[checks],Values[dataChecks]],TrueQ], Echo[{checks,dataChecks}]; Exit[1]];
report = <|"status" -> "PASS", "kernel" -> $Version, "symbolic_checks" -> checks,
  "released_data_checks" -> dataChecks, "word_count" -> Times@@Dimensions[words],
  "source_sha256" -> AssociationMap[IntegerString[FileHash[FileNameJoin[{here,#}],"SHA256"],16,64]&,
    {"LoaderIdentities.wl","run.wls","build.json"}]|>;
Export[FileNameJoin[{here,"wl-checks.json"}],report,"RawJSON"];
