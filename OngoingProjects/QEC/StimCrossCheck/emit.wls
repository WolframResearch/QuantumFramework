#!/usr/bin/env wolframscript
(* Emits the memory-experiment circuits in Stim syntax, together with this
   package's own EXACT detector firing rates and observable flip rates, so that
   check.py can compare the two.

   Usage:  wolframscript -file emit.wls [outputDirectory]
   Default output directory is ./out beside this script. *)

$root = ParentDirectory[DirectoryName[$InputFileName], 3];
PacletDirectoryLoad[FileNameJoin[{$root, "QuantumFramework"}]];
Needs["Wolfram`QuantumFramework`"];
Get[FileNameJoin[{$root, "OngoingProjects", "QEC", "QECCore", "QECCore.wl"}]];
Needs["Wolfram`QuantumFramework`QEC`"];
QECClearCache[];

out = Replace[$ScriptCommandLine[[2 ;;]], {{d_String, ___} :> d, _ -> FileNameJoin[{DirectoryName[$InputFileName], "out"}]}];
If[! DirectoryQ[out], CreateDirectory[out, CreateIntermediateDirectories -> True]];

cases = {
    {"bf-circ",     QECCode["BitFlipCode"], QECNoiseModel["Circuit", 1/200], 2},
    {"bf-phen",     QECCode["BitFlipCode"], QECNoiseModel["Depolarizing", 1/100, "MeasurementError" -> 1/50], 2},
    {"steane-circ", QECCode["SteaneCode"],  QECNoiseModel["Circuit", 1/500], 2},
    {"shor-phen",   QECCode["ShorCode"],    QECNoiseModel["Depolarizing", 1/50, "MeasurementError" -> 1/100], 1},
    {"5q-circ",     QECCode["5QubitCode"],  QECNoiseModel["Circuit", 1/500], 2}
};

Do[
    Module[{name = c[[1]], code = c[[2]], noise = c[[3]], r = c[[4]], dem},
        Export[FileNameJoin[{out, name <> ".stim"}], QECStimCircuit[code, noise, r], "Text"];
        dem = QECDetectorModel[code, noise, r];
        Export[FileNameJoin[{out, name <> ".json"}], <|
            "name" -> name, "rounds" -> r,
            "detectors" -> dem["Detectors"], "faults" -> dem["Faults"],
            "detectorRates" -> N[dem["DetectorRates"]],
            "observableRates" -> N[dem["ObservableRates"]]
        |>, "JSON"];
        Print[name, ": ", dem["Faults"], " faults, ", dem["Detectors"], " detectors"];
    ],
    {c, cases}];
