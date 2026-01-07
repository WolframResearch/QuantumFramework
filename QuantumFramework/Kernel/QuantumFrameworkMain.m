
Begin["Wolfram`QuantumFramework`Loader`"]

pacletInstalledQ[paclet_, version_] := AnyTrue[Through[PacletFind[paclet]["Version"]], ResourceFunction["VersionOrder"][#, version] <= 0 &]

If[ ! pacletInstalledQ["IBMQuantumPlatform", "0.0.2"],
    PacletInstall[PacletObject["Wolfram/QuantumFramework"]["AssetLocation", "IBMQuantumPlatform.paclet"]]
]

If[ ! pacletInstalledQ["Wolfram/TensorNetworks", "1.0.4"],
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

