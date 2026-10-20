Package["Wolfram`QuantumFramework`Init`"]

PackageImport["Wolfram`QuantumFramework`"]

PackageImport["Wolfram`QuantumFramework`PackageScope`"]

(* Keep shorthand caching local instead of looking up Memoize during loading. Held keys
   distinguish pattern-bearing calls, and the saved rules preserve dispatch.
   Calls with messages, nested failures or no matching rule are not cached. *)
With[{cache = Unique["shorthandCache$"], unanswered = Unique["shorthandUnanswered$"]},
    With[{rules = Dispatch[Append[DownValues[FromOperatorShorthand], HoldPattern[_] :> unanswered]]},
        cache = <||>;
        DownValues[FromOperatorShorthand] = {
            HoldPattern[FromOperatorShorthand[args___]] :> With[{held = HoldComplete[args]},
                With[{key = Hash[held]},
                    With[{result = With[{hit = Lookup[cache, key, Missing[]]},
                        If[ ! MissingQ[hit] && First[hit] === held,
                            Last[hit],
                            Module[{computed, messageFree = True},
                                Check[computed = Replace[Unevaluated[FromOperatorShorthand[args]], rules], messageFree = False];
                                If[ messageFree && computed =!= unanswered && FreeQ[computed, $Failed | $Aborted | _Failure],
                                    cache[key] = {held, computed}
                                ];
                                computed
                            ]
                        ]
                    ]},
                        result /; result =!= unanswered
                    ]
                ]
            ]
        }
    ]
]
