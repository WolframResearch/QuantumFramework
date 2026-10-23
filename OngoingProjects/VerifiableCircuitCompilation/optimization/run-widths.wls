root=DirectoryName[DirectoryName[$InputFileName]];
Get[FileNameJoin[{root,"optimization","WidthOptimization.wl"}]];
source=FileNameJoin[{root,"feasibility","B-quadratic-192.json"}];
If[FileHash[source,"SHA256","HexString"]!="2b768e13c7054d0e58bc4edf7afc0ff64559685aa309ad250b275fd9fd6a4c8d",Exit[1]];
rows=Import[source,"RawJSON"]["table"];
(* A bounded precision comparison, including deliberately inadequate widths. *)
widths={{29,25,19},{30,22,16},{30,23,14},{30,23,15},{30,23,16},{30,23,17},{30,24,16},{31,23,16}};
certificates=InvSqrtWidths`CertifyWidths[rows,#,15/(16 2^24)]& /@ widths;
Export[FileNameJoin[{root,"optimization","width-certificates.json"}],
  <|"kernel"->$Version,"source_sha256"->FileHash[source,"SHA256","HexString"],
    "certificate_source_sha256"->FileHash[FileNameJoin[{root,"optimization","WidthOptimization.wl"}],"SHA256","HexString"],
    "certificates"->certificates|>,"RawJSON"];
Echo[Map[KeyTake[#,{"fractional_bits","register_bits","certified"}]&,certificates]];
