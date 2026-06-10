
Begin["Wolfram`QuantumFramework`Loader`"]

pacletInstalledQ[paclet_, version_] := AnyTrue[Through[PacletFind[paclet]["Version"]], ResourceFunction["VersionOrder"][#, version] <= 0 &]

If[ ! pacletInstalledQ["IBMQuantumPlatform", "0.0.4"],
    PacletInstall[PacletObject["Wolfram/QuantumFramework"]["AssetLocation", "IBMQuantumPlatform.paclet"]]
]

(* Register the IBM Quantum Platform service connection as part of loading
   QuantumFramework, so ServiceConnect["IBMQuantumPlatform"] works with no
   separate Needs["IBMQuantumPlatform`"] from the user. *)
Needs["IBMQuantumPlatform`"]

If[ ! pacletInstalledQ["Wolfram/TensorNetworks", "1.0.5"],
    PacletInstall["https://www.wolframcloud.com/obj/wolframquantumframework/TensorNetworks.paclet"]
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

