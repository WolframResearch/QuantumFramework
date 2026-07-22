(* ::Package:: *)

Package["Wolfram`QuantumFramework`SecondQuantization`"]

PackageImport["Wolfram`QuantumFramework`"]

PackageExport["CoherentState"]

PackageExport["FockState"]

PackageExport["ThermalState"]

PackageExport["CatState"]


FockState::usage =
"\!\(\*RowBox[{\"FockState\", \"[\", RowBox[{StyleBox[\"n\", \"TI\"]}], \"]\"}]\) gives the Fock (number) state |\!\(\*StyleBox[\"n\", \"TI\"]\)\[RightAngleBracket].\n\!\(\*RowBox[{\"FockState\", \"[\", RowBox[{RowBox[{\"{\", RowBox[{StyleBox[\"n1\", \"TI\"], \",\", StyleBox[\"n2\", \"TI\"], \",\", \"\[Ellipsis]\"}], \"}\"}]}], \"]\"}]\) gives the multi-mode Fock state |\!\(\*StyleBox[\"n1\", \"TI\"]\), \!\(\*StyleBox[\"n2\", \"TI\"]\), ...\[RightAngleBracket].\n\!\(\*RowBox[{\"FockState\", \"[\", RowBox[{\"\[Ellipsis]\", \",\", StyleBox[\"size\", \"TI\"]}], \"]\"}]\) specifies the Fock space \!\(\*StyleBox[\"size\", \"TI\"]\) of all the defined modes (default: \!\(\*StyleBox[\"$FockSize\", \"TI\"]\)).";
FockState::clip = "Index `1` was outside the valid range {0, `2`} and has been clipped.";

FockState[n_Integer, size_ : $FockSize] := Block[{nEff = n},

  If[n < 0 || n >= size,
  
    Message[FockState::clip, n, size - 1];
    
    nEff = Clip[n, {0, size - 1}];
  ];
  
  QuantumState[SparseArray[{nEff + 1 -> 1}, size], size, "Label"-> Ket[{n}]]
];



FockVals::len="Values of `1` must be non negative integers less that the desired size of the space: `2`";

FockState[vals_List, size_:$FockSize]:= If[!AllTrue[vals,IntegerQ[#]&&(0<=#<size)&], 

	Message[FockVals::len,vals,size],

	With[{idx = FromDigits[vals, size] + 1, dim = size^Length[vals]},
      QuantumState[SparseArray[{idx -> 1}, dim], ConstantArray[size, Length[vals]]
      ,"Label"->Ket[vals]]
    ]
  ]



CoherentState::usage =
"\!\(\*RowBox[{\"CoherentState\", \"[\", \"]\"}]\) returns a parametric coherent state |\[Alpha]\[RightAngleBracket] with formal parameter \[FormalAlpha].\n\!\(\*RowBox[{\"CoherentState\", \"[\", RowBox[{StyleBox[\"size\", \"TI\"]}], \"]\"}]\) specifies the Fock space \!\(\*StyleBox[\"size\", \"TI\"]\) of the parametric state (default: $FockSize)\n\!\(\*RowBox[{\"CoherentState\", \"[\", RowBox[{StyleBox[\"size\", \"TI\"], \",\", \"\\\"Normalized\\\"->\", StyleBox[\"bool\", \"TI\"]}], \"]\"}]\) option for defining if the coherent state is normalized or not(default: \!\(\*StyleBox[\"True\", \"TI\"]\))";


Options[CoherentState] = {"Normalized" -> True};

CoherentState[size_Integer: $FockSize, OptionsPattern[]] := Block[{n = 0},
 
  If[OptionValue["Normalized"], #["Normalize"] &, Identity] @
  
   QuantumState[
    
    NestList[(n++; # \[FormalAlpha] / Sqrt[n]) &, 1, size - 1], 
    
    size, 
    
    "Parameters" -> \[FormalAlpha],
    
    "Label"-> HoldForm[StringForm["CoherentState[``]",\[FormalAlpha]]]
   ]
 ]



ThermalState::usage =
"\!\(\*RowBox[{\"ThermalState\", \"[\", RowBox[{StyleBox[\"nbar\", \"TI\"]}], \"]\"}]\) gives a thermal mixed state with mean photon number \!\(\*StyleBox[\"nbar\", \"TI\"]\).\n\!\(\*RowBox[{\"ThermalState\", \"[\", RowBox[{StyleBox[\"nbar\", \"TI\"], \",\", StyleBox[\"size\", \"TI\"]}], \"]\"}]\) specifies the Fock space \!\(\*StyleBox[\"size\", \"TI\"]\) (default: \!\(\*StyleBox[\"$FockSize\", \"TI\"]\)).";

ThermalState[nbar_, size_:$FockSize] :=

    QuantumState[
    
        DiagonalMatrix[1 / (1 + nbar)
        
            Table[(nbar / (1 + nbar)) ^ n,
            
           {n, 0, size-1}
          ]
        ], size, "Label"->StringForm["ThermalState[``]",nbar]]["Normalize"]



CatState::usage =
"\!\(\*RowBox[{\"CatState\", \"[\", \"]\"}]\) returns a parametric Schr\[ODoubleDot]dinger cat state \!\(\*SubscriptBox[\(N\), \(\[Alpha]\)]\)(|\[Alpha]\[RightAngleBracket] + exp(i\[Phi])|-\[Alpha]\[RightAngleBracket]), with formal parameters \[FormalAlpha] and \[FormalPhi].\n\!\(\*RowBox[{\"CatState\", \"[\", RowBox[{StyleBox[\"size\", \"TI\"]}], \"]\"}]\) specifies the Fock space \!\(\*StyleBox[\"size\", \"TI\"]\) (default: \!\(\*StyleBox[\"$FockSize\", \"TI\"]\)).";
 
CatState[size_: $FockSize] := With[{ns = Range[0, size - 1]},
  
    QuantumState[
    
      (\[FormalAlpha]^ns / Sqrt[ns!]) (1 + E^(I \[FormalPhi]) (-1)^ns),
      
      size,
      
      "Parameters" -> {\[FormalAlpha], \[FormalPhi]},
      
      "Label"-> HoldForm[StringForm["CatState[``,``]",\[FormalAlpha], \[FormalPhi]]]
      
    ]["Normalize"]
  ]
