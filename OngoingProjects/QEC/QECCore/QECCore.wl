(* ::Package:: *)

(* Loader for the rebuilt QEC core.

   During the rebuild these files live in OngoingProjects/QEC/QECCore/ and are
   read with a plain Get, so that a test cycle costs nothing.  They are written to
   move verbatim into QuantumFramework/Kernel/QEC/ once the core is clean, at
   which point this loader is replaced by the paclet's own and the contexts are
   registered in PacletInfo.wl. *)

Module[{dir = DirectoryName[$InputFileName]},
    Get[FileNameJoin[{dir, #}]] & /@ {
        "GF2.wl",
        "Pauli.wl",
        "Code.wl",
        "Structure.wl",
        "Syndrome.wl",
        "Encoder.wl",
        "Circuit.wl",
        "Constructions.wl",
        "Families.wl",
        "Noise.wl",
        "ErrorRate.wl",
        "DetectorModel.wl",
        "Memory.wl",
        "Stim.wl",
        "Cache.wl"
    }
]
