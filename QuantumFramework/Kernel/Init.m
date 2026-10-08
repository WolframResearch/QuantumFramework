Package["Wolfram`QuantumFramework`Init`"]

PackageImport["Wolfram`QuantumFramework`"]

PackageImport["Wolfram`QuantumFramework`PackageScope`"]

(* FromOperatorShorthand is memoized by the Function Repository's Memoize. Until the repository
   publishes it, the name resolves only where the resource is registered locally; anywhere else the
   lookup fails and the resource's cloud deployment is used. *)
With[{memoize = ResourceFunction["Memoize"]},
    If[ FailureQ[memoize],
        ResourceFunction[ResourceObject["https://www.wolframcloud.com/obj/nikm/DeployedResources/Function/Memoize/"]],
        memoize
    ][FromOperatorShorthand]
]
