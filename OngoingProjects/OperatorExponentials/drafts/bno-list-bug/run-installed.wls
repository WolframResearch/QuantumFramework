(* Runs the probe on whatever QF Needs finds first (the installed paclet). *)
Get[FileNameJoin[{DirectoryName[$InputFileName], "load.wl"}]];
Get[FileNameJoin[{DirectoryName[$InputFileName], "probe.wl"}]];
