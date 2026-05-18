(* ::Package:: *)

(* ::Package:: *)
(**)


(* SecondQuantization Loader *)

Get[FileNameJoin[{DirectoryName[$InputFileName], "SecondQuantization", "Utils.m"}]];
Get[FileNameJoin[{DirectoryName[$InputFileName], "SecondQuantization", "States.m"}]];
Get[FileNameJoin[{DirectoryName[$InputFileName], "SecondQuantization", "Operators.m"}]];
Get[FileNameJoin[{DirectoryName[$InputFileName], "SecondQuantization", "BosonicAlgebra.m"}]];
Get[FileNameJoin[{DirectoryName[$InputFileName], "SecondQuantization", "BlasiakFormula.m"}]];
Get[FileNameJoin[{DirectoryName[$InputFileName], "SecondQuantization", "PhaseSpaceRepresentations.m"}]];


SetFockSpaceSize[];
