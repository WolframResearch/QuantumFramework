
Begin["Wolfram`QuantumFramework`Loader`"]

pacletInstalledQ[paclet_, version_] := AnyTrue[Through[PacletFind[paclet]["Version"]], ResourceFunction["VersionOrder"][#, version] <= 0 &]

If[ ! pacletInstalledQ["IBMQuantumPlatform", "0.0.4"],
    PacletInstall[PacletObject["Wolfram/QuantumFramework"]["AssetLocation", "IBMQuantumPlatform.paclet"]]
]

(* Register the IBM Quantum Platform service connection as part of loading
   QuantumFramework, so ServiceConnect["IBMQuantumPlatform"] works with no
   separate Needs["IBMQuantumPlatform`"] from the user. *)
Needs["IBMQuantumPlatform`"]

(* A dependency below its floor installs from the Paclet Repository: through the
   resource system first, which has a version as soon as it is published, then
   through the repository's paclet site, whose index lags behind. An install that
   fails reports its own messages. When neither route has a version that meets the
   floor, QuantumFramework loads with what is installed and says so.

   The floors are not cosmetic. Below TensorNetworks 1.1.0 the 20-qubit Trotter
   quench of the dense-circuit probe runs out of 4 GB, and a network whose diagonal
   gates share their wires' indices, the form the apply path builds from 1.1.0 on,
   runs out of memory even for the 14-qubit QFT; below 1.0.10 the "NetGraph"
   contraction method also fails on a phase-space contraction with
   FindPermutation::norel. Below Arrays 1.4.1 an evolved state read before its time
   is bound fails with Interpolation::inddp when the time grid repeats a point, as
   NDSolve's does at a discontinuity of the Hamiltonian. *)

requirePaclet::unmet = "QuantumFramework needs `1` `2` or later, which could not be installed from the Paclet Repository. Installed: `3`."

requirePaclet[paclet_String, version_String] := If[ ! pacletInstalledQ[paclet, version],
    PacletInstall[ResourceObject[paclet]];
    If[ ! pacletInstalledQ[paclet, version], PacletInstall[paclet]];
    If[ ! pacletInstalledQ[paclet, version],
        Message[requirePaclet::unmet, paclet, version, Replace[Through[PacletFind[paclet]["Version"]], {{v_, ___} :> v, {} -> "none"}]]
    ]
]

requirePaclet["Wolfram/TensorNetworks", "1.1.0"]

requirePaclet["Wolfram/Arrays", "1.4.1"]

$ContextAliases["H`"] = "WolframInstitute`Hypergraph`"

ClearAll["Wolfram`QuantumFramework`*", "Wolfram`QuantumFramework`**`*"]

PacletManager`Package`loadWolframLanguageCode[
    "Wolfram/QuantumFramework",
    "Wolfram`QuantumFramework`",
    ParentDirectory[DirectoryName[$InputFileName]],
    "Kernel/QuantumFramework.m",
    "AutoUpdate" -> False,
    "AutoloadSymbols" -> {},
    "HiddenImports" -> {},
    "SymbolsToProtect" -> {}
]

End[]

(* this turns PackageScope into a valid package usable with PackageImport *)
Block[{$ContextPath},
    BeginPackage["Wolfram`QuantumFramework`PackageScope`"];
    EndPackage[];

    Get["Wolfram`QuantumFramework`Init`"]
]

