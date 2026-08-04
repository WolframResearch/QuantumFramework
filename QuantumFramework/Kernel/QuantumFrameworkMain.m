
Begin["Wolfram`QuantumFramework`Loader`"]

pacletInstalledQ[paclet_, version_] := AnyTrue[Through[PacletFind[paclet]["Version"]], ResourceFunction["VersionOrder"][#, version] <= 0 &]

If[ ! pacletInstalledQ["IBMQuantumPlatform", "0.0.4"],
    PacletInstall[PacletObject["Wolfram/QuantumFramework"]["AssetLocation", "IBMQuantumPlatform.paclet"]]
]

(* Register the IBM Quantum Platform service connection as part of loading
   QuantumFramework, so ServiceConnect["IBMQuantumPlatform"] works with no
   separate Needs["IBMQuantumPlatform`"] from the user. *)
Needs["IBMQuantumPlatform`"]

(* Both are published, so they install from the repository by name.  The
   version floors are not cosmetic: the "NetGraph" contraction method needs
   TensorNetworks 1.0.10, and below it a phase-space contraction fails with
   FindPermutation::norel rather than producing a net. *)

If[ ! pacletInstalledQ["Wolfram/TensorNetworks", "1.0.10"],
    PacletInstall["Wolfram/TensorNetworks"]
]

If[ ! pacletInstalledQ["Wolfram/Arrays", "1.3.2"],
    PacletInstall["Wolfram/Arrays"]
]

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

