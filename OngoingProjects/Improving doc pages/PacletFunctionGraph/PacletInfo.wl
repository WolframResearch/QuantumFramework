(* ::Package:: *)

PacletObject[
    <|
        "Name"           -> "Wolfram/PacletFunctionGraph",
        "Description"    -> "Build and render a directed function-relation graph for any Wolfram Language paclet.",
        "Creator"        -> "Mohammad (Mads) Bahrami",
        "License"        -> "MIT",
        "PublisherID"    -> "Wolfram",
        "Version"        -> "0.1.0",
        "WolframVersion" -> "13.3+",
        "PrimaryContext" -> "Wolfram`PacletFunctionGraph`",
        "Extensions"     -> {
            {
                "Kernel",
                "Root"    -> "Kernel",
                "Context" -> {"Wolfram`PacletFunctionGraph`"},
                "Symbols" -> {
                    "Wolfram`PacletFunctionGraph`BuildFunctionGraph",
                    "Wolfram`PacletFunctionGraph`LoadFunctionGraph",
                    "Wolfram`PacletFunctionGraph`PacletFlowEdges",
                    "Wolfram`PacletFunctionGraph`SymbolNeighborhoodGraph",
                    "Wolfram`PacletFunctionGraph`PacletOverviewGraph",
                    "Wolfram`PacletFunctionGraph`$FunctionGraphKindColors"
                }
            }
        }
    |>
]
