(* ::Package:: *)

Package["Wolfram`QuantumFramework`SecondQuantization`"]

PackageImport["Wolfram`QuantumFramework`"]

PackageExport["CoherentState"]

PackageExport["FockState"]

PackageExport["ThermalState"]

PackageExport["CatState"]


FockState::usage = 
"\!\(FockState[n]\) Creates a Fock (number) state |n\[RightAngleBracket] in the Fock space.
\!\(FockState[{n1, n2, ...}]\) Creates a multi-mode Fock state |n1, n2, ...\[RightAngleBracket].
\!\(FockState[\[Ellipsis],size]\) Specifies the Fock space size of all the defined modes (default: $FockSize).";

FockState::clip = "Index `1` was outside the valid range {0, `2`} and has been clipped.";

FockState[n_Integer, size_ : $FockSize] := Block[{nEff = n},

  If[n < 0 || n >= size,
  
    Message[FockState::clip, n, size - 1];
    
    nEff = Clip[n, {0, size - 1}];
  ];
  
  QuantumState[SparseArray[{nEff + 1 -> 1}, size], size]
];



FockVals::len="Values of `1` must be non negative integers less that the desired size of the space: `2`";

FockState[vals_List, size_:$FockSize]:= If[!AllTrue[vals,IntegerQ[#]&&(0<=#<size)&], 

	Message[FockVals::len,vals,size],

	With[{idx = FromDigits[vals, size] + 1, dim = size^Length[vals]},
      QuantumState[SparseArray[{idx -> 1}, dim], ConstantArray[size, Length[vals]]]
    ]
  ]



CoherentState::usage = 
"\!\(CoherentState[]\) Returns a parametric coherent state |\[Alpha]\[RightAngleBracket] with formal parameter \[FormalAlpha].
\!\(CoherentState[size]\) Specifies the Fock space size of the parametric state (default: $FockSize)
\!\(CoherentState[size, 'Normalized'-> ]\) Specifies if the coherent state is normalized or not, the default is True";


Options[CoherentState] = {"Normalized" -> True};

CoherentState[size_Integer: $FockSize, OptionsPattern[]] := Block[{n = 0},
 
  If[OptionValue["Normalized"], #["Normalize"] &, Identity] @
  
   QuantumState[
    
    NestList[(n++; # \[FormalAlpha] / Sqrt[n]) &, 1, size - 1], 
    
    size, 
    
    "Parameters" -> \[FormalAlpha]
   ]
 ]



ThermalState::usage =
"\!\(ThermalState[nbar]\) Creates a thermal (Bose-Einstein) mixed state with mean photon number nbar.
\!\(ThermalState[nbar, size]\) Specifies the Fock space size (default: $FockSize).";

ThermalState[nbar_, size_:$FockSize] :=

    QuantumOperator[
    
        DiagonalMatrix[1 / (1 + nbar)
        
            Table[(nbar / (1 + nbar)) ^ n,
            
           {n, 0, size-1}
          ]
        ], size]["MatrixQuantumState"]["Normalize"]



CatState::usage = 
"\!\(CatState[]\)Returns a parametric Schr\[ODoubleDot]dinger cat state \!\(\*SubscriptBox[\(N\), \(\[Alpha]\)]\)(|\[Alpha]\[RightAngleBracket] + exp(i\[Phi])|-\[Alpha]\[RightAngleBracket]), with formal parameters \[FormalAlpha] and \[FormalPhi].
\!\(CatState[size]\) Specifies the Fock space size (default: $FockSize).";
 
CatState[size_: $FockSize] := With[{ns = Range[0, size - 1]},
  
    QuantumState[
    
      (\[FormalAlpha]^ns / Sqrt[ns!]) (1 + E^(I \[FormalPhi]) (-1)^ns),
      
      size,
      
      "Parameters" -> {\[FormalAlpha], \[FormalPhi]}
      
    ]["Normalize"]
  ]
