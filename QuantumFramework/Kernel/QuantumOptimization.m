(* ::Package:: *)

Package["Wolfram`QuantumFramework`QuantumOptimization`"]

PackageImport["Wolfram`QuantumFramework`"]

PackageExport["GradientDescent"]

PackageExport["QuantumNaturalGradientDescent"]

PackageExport["SPSRGradientValues"]

PackageExport["ASPSRGradientValues"]

PackageExport["FubiniStudyMetricTensor"]

PackageExport["QuantumLinearSolve"]

PackageExport["ParametrizedLayer"]

PackageExport["EntanglementLayer"]

PackageExport["GenerateParameters"]

PackageExport["ClassiqSetup"]

PackageExport["FubiniStudyMetricTensorLayers"]

PackageExport["QuantumLockingMechanism"]

PackageExport["QuantumUnlockingMechanism"]

PackageExport["QuantumAdiabaticEvolve"]


(* ::Section:: *)
(*Quantum Optimization Functionalities*)


ParametrizedLayer::dimensions = "Number of qubits and parameters of unequal length"

Options[ParametrizedLayer]={"Symbol"->"\[Theta]"};
ParametrizedLayer[op_,range_List,OptionsPattern[]]:=Sequence@@Table[op[Symbol[OptionValue["Symbol"]<>ToString[i]]]->Flatten[Position[range,i]],{i,range}]

ParametrizedLayer[op_,qubits_List,index_List,OptionsPattern[]]:=Module[{range},
	
	If[!MatchQ[Length[qubits],Length[index]],
		Message[ParametrizedLayer::dimensions];
		Return[$Failed]
	];
	
	
	Sequence@@Table[op[Symbol[OptionValue["Symbol"]<>ToString[i]]]->qubits[[Flatten[Position[index,i]]]],{i,index}]
	
	]


Options[EntanglementLayer] = {"Entanglement" -> Automatic};
EntanglementLayer[cop_, range_List, OptionsPattern[]] := Module[{opt, pairs},
  opt = OptionValue["Entanglement"];
  pairs = 
  Which[
	    MatchQ[opt, "Linear"], Partition[range, 2, 1],
	    MatchQ[opt, "ReverseLinear"], Reverse /@ Partition[range, 2, 1],
	    MatchQ[opt, "Pairwise"], Join @@ Reverse[Values[GroupBy[Drop[Sort@range,-1], EvenQ]]] /. x_Integer :> {x, x + 1},
	    MatchQ[opt, "Circular"], Append[#, {#[[-1, -1]], #[[1, 1]]}] &@Partition[range, 2, 1],
	    True, Subsets[range, {2}]
    ];
  
	Sequence @@ Thread[cop -> pairs]
   ]

Options[GenerateParameters]={"Symbol"->"\[Theta]"};
GenerateParameters[NQubits_,NLayers_,OptionsPattern[]]:=Table[Symbol[OptionValue["Symbol"]<>ToString[i] ],{i,1,NLayers *NQubits}]


ClearAll[QuantumAdiabaticEvolution];
QuantumAdiabaticEvolution/:MakeBoxes[obj:QuantumAdiabaticEvolution[asc_?myObjectAscQ],form:(StandardForm|TraditionalForm)]:=Module[
{s,energies,icon,above,below,props},

props={asc["MinimalEnergyGap"],Flatten@asc["MaxAdiabaticCoupling"],asc["AdiabaticTimeEstimate"]}/.x_Real:>SetPrecision[x,3];
		above={
			{BoxForm`SummaryItem[{"","Adiabatic Evolution"}],SpanFromLeft},
			{BoxForm`SummaryItem[{"MinimalEnergyGap: ",props[[1]]}],SpanFromLeft},
			{BoxForm`SummaryItem[{"MaxAdiabaticCoupling: ",Flatten@props[[2]]}],BoxForm`SummaryItem[{"AdiabaticTimeEstimate: ",props[[3]]}]}
	};
	below={
			{BoxForm`SummaryItem[{"Parameters: ",ToString@asc["Parameters"]}],BoxForm`SummaryItem[{"Parameters Domain: ","[0,1]"}]}
		};

	icon=adiabticPlot[{{1/8 (4-Sqrt[16-63s+63s^2]),1/8 (4+Sqrt[16-63s+63s^2])},s->1/2},"Icon",s];

	BoxForm`ArrangeSummaryBox[QuantumAdiabaticEvolution,obj,icon,above,below,form,"Interpretable"->Automatic]

];

myObjectAscQ[asc_?AssociationQ]:=AllTrue[{
"Hamiltonian","Energies","Eigenvectors","Eigensystem","LowestEnergies","EnergyGap","MinimalEnergyGap","LowestEnergyEigenvector","AdiabaticCoupling","MaxAdiabaticCoupling","AdiabaticTimeEstimate","Parameters","Properties"},KeyExistsQ[asc,#]&]

myObjectAscQ[_]=False;

QuantumAdiabaticEvolution[asc_?AssociationQ][prop_]:=Lookup[asc,prop]

QuantumAdiabaticEvolution[asc_?AssociationQ][plots_List/;SubsetQ[{"EnergySpectrumPlot","SpectralGapPlot","AdiabaticCouplingPlot","AdiabaticPathPlot"},plots]]:=QuantumAdiabaticEvolution[asc][#]&/@plots

QuantumAdiabaticEvolution[asc_?AssociationQ][plot:Alternatives@@{"EnergySpectrumPlot","SpectralGapPlot","AdiabaticCouplingPlot","AdiabaticPathPlot"}]:=Module[{ass},ass=Association@Thread[{"EnergySpectrumPlot","SpectralGapPlot","AdiabaticCouplingPlot","AdiabaticPathPlot"}->{"Energies","EnergyGap","AdiabaticCoupling","LowestEnergyEigenvector"}];
	adiabticPlot[asc[ass[plot]],plot,asc["Parameters"]]
]



ClearAll[QuantumAdiabaticEvolve];
QuantumAdiabaticEvolve::gap ="The minimum spectral gap is zero. The adiabatic condition is violated and the evolution is not guaranteed to reach the ground state.";
QuantumAdiabaticEvolve::num="The eigensystem contains Root objects. Symbolic expressions may remain unresolved and the computation will proceed using numerical approximations.";
Options[QuantumAdiabaticEvolve]={"Parameter"->\[FormalS],"TimeScaling"->10};

QuantumAdiabaticEvolve[Hb_QuantumOperator,Hp_QuantumOperator,opts:OptionsPattern[]]:=
Module[{s,H,output,plotoutput,cachedResults,reporter,energies,gap,gmin,\[Psi]E0,\[Psi]E1,lowestenergies,eigensystem,sortedstates,\[Xi]function,\[Xi],smin,smax,sortedeigensystem,plotnames,functions,plots,simplifier,mingap,prop},

prop={"Hamiltonian","Energies","Eigenvectors","Eigensystem","LowestEnergies","EnergyGap","MinimalEnergyGap","LowestEnergyEigenvector","AdiabaticCoupling","MaxAdiabaticCoupling","AdiabaticTimeEstimate","EnergySpectrumPlot","SpectralGapPlot","AdiabaticCouplingPlot","AdiabaticPathPlot"};

s=OptionValue["Parameter"];
H=QuantumOperator[(1-s)*Hb+s*Hp,"Parameters"->{s}];

eigensystem=Simplify[H["Eigensystem"],0<=s<=1];

{sortedeigensystem,gap,mingap}=adiabaticTest[Thread@eigensystem,s];

energies=sortedeigensystem[[All,1]];

lowestenergies=sortedeigensystem[[;;2]][[All,1]];

If[MatchQ[mingap,$Failed],Message[QuantumAdiabaticEvolve::gap];Return[$Failed,Module],{gmin,smin}=mingap];

	If[gmin<=0,Message[QuantumAdiabaticEvolve::gap];Return[$Failed,Module],
		{\[Psi]E0,\[Psi]E1}=QuantumState[#,"Parameters"->{s}]&/@sortedeigensystem[[;;2]][[All,2]];
		\[Xi]function=If[FreeQ[sortedeigensystem[[;;2]][[All,2]][[1]],_Root,Infinity],
			FullSimplify[Abs["Scalar"//SuperDagger[\[Psi]E1]@D[H,s]@\[Psi]E0],s\[Element]Reals&&0<=s<=1],
			Message[QuantumAdiabaticEvolve::num];
			Abs["Scalar"//SuperDagger[\[Psi]E1]@D[H,s]@\[Psi]E0]
		];
		{\[Xi],smax}=Maximize[{\[Xi]function,0<=s<=1},s]
	];


QuantumAdiabaticEvolution[<|"Hamiltonian"->H,"Eigensystem"->sortedeigensystem,"Energies"->energies,"Eigenvectors"->sortedeigensystem[[All,2]],"LowestEnergies"->lowestenergies,"EnergyGap"->gap,"MinimalEnergyGap"->{gmin,smin},"LowestEnergyEigenvector"->\[Psi]E0, "MaxAdiabaticCoupling"->{\[Xi],smax},"AdiabaticCoupling"->\[Xi]function,"AdiabaticTimeEstimate"->OptionValue["TimeScaling"]*\[Xi]/gmin^2,"Parameters"->s,"Properties"->prop|>]

]

adiabaticTest[eigensystem_,s_]:=Module[{sorted,energies,solutions,cases,gap,mingap},
sorted=SortBy[eigensystem,ReplaceAll[#[[1]],s->10^-10]&]; 
If[!DuplicateFreeQ[ReplaceAll[sorted[[All,1]],s->10^-10]],sorted=SortBy[eigensystem,ReplaceAll[#[[1]],s->0.1]&]];
	energies=sorted[[All,1]];
	gap=Simplify[Abs[Subtract@@Reverse@energies[[;;2]]],s\[Element]Reals&&0<=s<=1];
	solutions=Quiet@Flatten@(NSolve[{energies[[1]]==# && 0<=s<=1},s,Reals])&/@energies[[2;;]];
	cases=Cases[solutions,List[_Rule]];
	If[MatchQ[cases,{}],
		mingap=First@SortBy[Flatten@Minimize[{Abs[energies[[1]]-#] , 0<=s<=1},s]&/@energies[[2;;]],#[[1]]&]; 
		{sorted,gap,mingap} 
		,
		{sorted,gap,$Failed}
	]
]

adiabticPlot[f_,name_String,s_]:=Switch[name,
	"EnergySpectrumPlot", Plot[Evaluate@Values[#],{s,0,1},Sequence[PlotRange -> All, AxesLabel -> {"s", "Eigenvalues"}, PlotLegends -> Keys[#], PlotLabel -> "Energy Spectrum", LabelStyle -> 13]]&@energyPlotLabels[f],
	"SpectralGapPlot", Plot[f,{s,0,1},Sequence[AxesLabel -> {"s", "\!\(\*g(s)\)"}, PlotRange -> {{0, 1}, {0, Automatic}}, PlotLabel -> "Spectral Gap", LabelStyle -> 13]],
	"AdiabaticCouplingPlot",Plot[f,{s,0,1},Sequence[PlotRange -> All, AxesLabel -> {"s", "\[Xi](s)"}, PlotLabel -> "Adiabatic Coupling", LabelStyle -> 13]],
	"AdiabaticPathPlot",Plot[Evaluate@Keys[#],{s,0,1},Sequence[AxesLabel -> {"\!\(\*s = t/T\)", "Probability"}, LabelStyle -> 13, PlotLegends -> Map[adiabticPlotLabels[#]& , Values[#]], PlotLabel -> "Adiabatic Path", PlotRange -> {{0, 1}, {0, 1}}]]&@GroupBy[Normal[f[s]["Probability"]],Values[#]&,Keys],
"Icon",Rasterize[Plot[Evaluate@#[[1]],{s,0,1},Sequence[Ticks -> None, Axes -> False, PlotRange -> {Values[Part[#, 2]] + {-0.1, 0.1}, {ReplaceAll[First[Part[#, 1]], Keys[Part[#, 2]] -> Values[Part[#, 2]] - 0.1], ReplaceAll[Last[Part[#, 1]], Keys[Part[#, 2]] -> Values[Part[#, 2]] + 0.1]}}, Epilog -> {Gray, Arrowheads[{{0.1, 1}, {-0.1, 0}}], Arrow[{{Values[Part[#, 2]], ReplaceAll[First[Part[#, 1]], Keys[Part[#, 2]] -> Values[Part[#, 2]]]}, {Values[Part[#, 2]], ReplaceAll[Last[Part[#, 1]], Keys[Part[#, 2]] -> Values[Part[#, 2]]]}}]}]]&@f,Background->None,ImageSize->10^1.9],
	_,None
]

adiabticPlotLabels[label_]:=If[Length@#>1,Tooltip[ToString[Ket[{StringJoin[ToString/@First@First[#]]}], FormatType -> TraditionalForm]<>" \[Ellipsis]",#],First@#]&@label
energyPlotLabels[energies_]:=Module[{grouping,rangeformat},
	grouping=Replace[Normal[GroupBy[Thread[List[Range[0,Length@#-1],#]&@(energies)],#[[2]]&]],{i_,E_}:>i,{3}];
	rangeformat=Replace[grouping,Rule[e_,x_List]/;Length[x]>=2:>e->StringRiffle[ToString/@x[[{1,-1}]],"-"],{1}];
	Replace[rangeformat,Rule[e_,x_]:>Replace[Rule[Subscript["E", x],e],s_List:>Sequence@@s,{2}],{1}]
]



(* ::Section::Closed:: *)
(*Classiq Integration*)


ClassiqSetup::PythonEvaluators="WARNING: Non Python evaluator found";

ClassiqSetup::"ClassiqInstallation"="WARNING: Classiq Installation failed";

ClassiqSetup::"UpdatingClassiq"="Installing latest Classiq version";

ClassiqSetup::"PythonVersion"="WARNING: Incompatible Python version detected. Please install a supported Python version (3.10\[Dash]3.12). Make sure a compatible Python installation is properly \!\(\*TemplateBox[{\"configured\", {URL[\"https://reference.wolfram.com/language/workflow/ConfigurePythonForExternalEvaluate.html\"], None}, \"https://reference.wolfram.com/language/workflow/ConfigurePythonForExternalEvaluate.html\", \"HyperlinkActionRecycled\", {\"HyperlinkActive\"}, BaseStyle -> {\"Hyperlink\"}, HyperlinkAction -> \"Recycled\"},\n\"HyperlinkTemplate\"]\) and \!\(\*TemplateBox[{\"registered\", {URL[\"https://reference.wolfram.com/language/ref/RegisterExternalEvaluator.html\"], None}, \"https://reference.wolfram.com/language/ref/RegisterExternalEvaluator.html\", \"HyperlinkActionRecycled\", {\"HyperlinkActive\"}, BaseStyle -> {\"Hyperlink\"}, HyperlinkAction -> \"Recycled\"},\n\"HyperlinkTemplate\"]\) for ExternalEvaluate";

Options[ClassiqSetup]={"CheckDependencies"->{},"InstallPackages"->{},"LatestClassiqVersion"->"0.86.0"};

ClassiqSetup[]:=ClassiqSetup[{"ClassiqVersion","Session"}];

ClassiqSetup[prop : _String | {__String} | All ,opts:OptionsPattern[]]:=Module[
	
	{evaluators,cachedResults,reporter,output,classiq,ver,dependencies,dep,allpackages,packages,pak,evaluator,session,list,versions,compatible,path},
		
		cachedResults=<||>;
		
		dependencies=OptionValue["CheckDependencies"];
		
		packages=OptionValue["InstallPackages"];
		
		output = Replace[prop, 
							All -> {"ClassiqVersion","Session","Evaluators","CheckDependencies","InstallPackages"}
						];
		
		If[!MatchQ[prop,All],
								
			If[!MatchQ[dependencies,{}],
				output=Append[output,"CheckDependencies"]
			];
			
			If[!MatchQ[packages,{}],
				output=Append[output,"InstallPackages"]
			];
		];				
	
		reporter[data_Association] := (
				PrependTo[cachedResults, data];  
				If[ContainsAll[Keys[cachedResults], Flatten[{output}]],
					Return[Dataset[cachedResults[[output]]],Module]
				];
			);		
	
		evaluators=FindExternalEvaluators["Python"];
		
		
		
		If[Length@Normal[evaluators[All,"Evaluator"]]<1, Message[ClassiqSetup::PythonEvaluators]; Return[$Failed],
			
			If[$VersionNumber > 14.0,
			
			versions=Quiet[(#->ExternalEvaluate[#,"import sys; sys.version_info"]/._Failure->$Failed/.x_List:>x[[;;3]])&/@ExternalEvaluators["Python"]];
			compatible=DeleteDuplicates@ReverseSortBy[Select[DeleteCases[versions,_->$Failed],#[[2,2]]>=10&&#[[2,2]]<=12&],#[[2]]&];
			If[MatchQ[compatible,{}],Return[ClassiqSetup::"PythonVersion",Module]];
			path=First[First[compatible]]["Evaluator"];
			session=StartExternalSession[<|"System" -> "Python", "Evaluator"->path ,"ID" -> "default-python-session"|>],
			
			versions=FindExternalEvaluators["Python"][All,"Executable"]//Normal//Values;
			versions=Table[py->Quiet@ExternalEvaluate[<|"System"->"Python","Executable"->py|>,"import sys; sys.version_info"]/._Failure->$Failed/.x_List:>x[[;;3]],{py,versions}];
			compatible=DeleteDuplicates@ReverseSortBy[Select[DeleteCases[versions,_->$Failed],#[[2,2]]>=10&&#[[2,2]]<=12&],#[[2]]&];
			If[MatchQ[compatible,{}],Return[ClassiqSetup::"PythonVersion",Module]];
			path=First[First[compatible]];
			session=StartExternalSession[<|"System"->"Python","Evaluator"->path,"ID"->"DefaultPythonSession"|>]
			];
		
			reporter[<|"Evaluators"->evaluators|>];
			
			reporter[<|"Session"->session|>];
			
			evaluator=session["Evaluator"];
		
			classiq=RunProcess[{evaluator,"-m","pip","--quiet","install","classiq"}];
			
					
			If[!installedQ["classiq",evaluator],
				Print[ClassiqSetup::ClassiqInstallation]; Return[$Failed],
				ver=StringTrim@RunProcess[{evaluator,"-c","import classiq; print(classiq.__version__)"},"StandardOutput"];
				If[!MatchQ[ver,OptionValue["LatestClassiqVersion"]],Print[ClassiqSetup::UpdatingClassiq];RunProcess[{evaluator,"-m","pip","--quiet","install","classiq","--upgrade"}]]
			];
			
			ver=StringTrim@RunProcess[{evaluator,"-c","import classiq; print(classiq.__version__)"},"StandardOutput"];
			
			reporter[<|"ClassiqVersion"->ver|>];
			
			dep=If[Length@dependencies>=1,
					<|(#->installedQ[#,evaluator])&/@dependencies|>,
					None
				];
			
			reporter[<|"CheckDependencies"->dep|>];
			
			
			pak=If[Length@packages>=1,
					RunProcess[{evaluator,"-m","pip","--quiet","install",#}]&/@packages;
					<|(#->installedQ[#,evaluator])&/@packages|>,
					None
			];
			
			reporter[<|"InstallPackages"->pak|>]
	
		];
]


installedQ[name_String,eval_]:=!MatchQ[If[StringQ[#],If[!MatchQ[#,""],#,$Failed],$Failed]&@RunProcess[{eval,"-m","pip","show",name},"StandardOutput"],$Failed]


(* ::Section::Closed:: *)
(*Example specific functions *)


(* ::Subsection:: *)
(*Quantum Natural Gradient Descent*)


FubiniStudyMetricTensor::param = "QuantumState parameters ill\[Dash]defined.";

Options[FubiniStudyMetricTensor]={"Parameters"->Automatic};

FubiniStudyMetricTensor[qstate:_QuantumState| (_List ? VectorQ),prop : _String | {__String} | All : "Result", opts:OptionsPattern[]]:=
Module[

	{var,stateVector,result,derivatives,output,cachedResults,reporter,matrix},

	output = Replace[prop, 
					All -> {"Result","Matrix","MatrixForm","Parameters","SparseArray"}
				];	

	cachedResults=<||>;	
	
	reporter[data_Association] := (
		PrependTo[cachedResults, data];  
		If[ContainsAll[Keys[cachedResults], Flatten[{output}]],
			Return[cachedResults[[output]],Module]
		];
		);		
									
	
	If[MatchQ[OptionValue["Parameters"],Automatic],
		var=qstate["Parameters"],
		var=OptionValue["Parameters"]
	];
	
	reporter[<|"Parameters"->var|>];
	
	If[!MatchQ[var,{_Symbol..}],
			Message[FubiniStudyMetricTensor::param];
			Return[$Failed]
	];
		
	stateVector=If[MatchQ[qstate,_QuantumState],
					Normal[qstate["StateVector"]],
					qstate
				];
	
	derivatives=Table[D[stateVector,i],{i,var}];
	
	result=Table[
			Re[ConjugateTranspose[derivatives[[i]]] . derivatives[[j]]]-(ConjugateTranspose[derivatives[[i]]] . stateVector)(ConjugateTranspose[stateVector] . derivatives[[j]]),
			{i,Length@derivatives},{j,Length@derivatives}
			];
	
	result=QuantumOperator[result,"Parameters"->var];
	
	reporter[<|"Result"->result,"SparseArray"->result["Matrix"]|>];
	
	matrix=result["Matrix"]//Normal//ComplexExpand//Simplify;

	reporter[<|"Matrix"->matrix,"MatrixForm"->MatrixForm[matrix]|>];


]


FubiniStudyMetricTensor[layers_List, parameters_List, initParameters_?(VectorQ[#,NumericQ]&)]:=Module[{probabilities,covariance,variance1,variance2},

	probabilities=Values[

					N[#[[2,2]][#[[2,1]]][#[[1]][AssociationThread[Rule[parameters,initParameters]]]]["Probabilities"]]&/@layers

					];

	covariance=(#[[4]]-(#[[3]]+#[[4]])*(#[[2]]+#[[4]]))&/@probabilities;

	variance1=(#*(1-#))&@Total[#[[{3,4}]]]&/@probabilities;

	variance2=(#*(1-#))&@Total[#[[{2,4}]]]&/@probabilities;

	BlockDiagonalMatrix[{{#[[1]],#[[2]]},{#[[2]],#[[3]]}}&/@Thread[Chop[{variance1,covariance,variance2},10^-8]]]
]


FubiniStudyMetricTensorLayers[qc_QuantumCircuitOperator, parameters_List]:=Module[{elements,input,layers},

	(*Re-generating Original QuantumCircuitOperator Inputs*)

		elements=DeleteCases[qc["Elements"],"Barrier"];

		input=#["Label"]->#["InputOrder"]&/@elements;

	(*Cleaning inputs (subscripts and unnecessary brackets) to correctly build subcircuits*)

		input=input/.Subscript[x_String,y_String]:>x<>y/."CNOT"[___]:>"CNOT";

	(*Generating layers by splitting when parametric variables are found, dropping useless gates*)

		input=Most[SplitBy[input,FreeQ[#,Alternatives@@parameters]&]];

	(*Building layers*)

		layers=Table[input[[;;i]],{i,Range[2,Length@input,2]}];

	(*Changing last parametric Pauli gates from each layer to their corresponding QuantumMeasurementOperator*)

		layers=MapAt[#/.Rule[(pauli:("RX"|"RY"|"RZ"))[_],n_]:>QuantumMeasurementOperator[StringDelete[pauli,"R"],n]&,layers,{All,-1}];

	(*Output: QuantumState of first layer (parametrization is possible for QuantumState from a Parametrized QuantumCircuit) and Measurements*)

		{QuantumCircuitOperator[Flatten[#[[;;-2]]],"Parameters"->parameters][],#[[-1]]}&/@layers
]


(* ::Subsection::Closed:: *)
(*Stochastic Parameter Shift-Rule*)


(* ::Input::Initialization::Plain:: *)
Options[SPSRGradientValues]={
"Shift"->\[Pi]/4.,

"ParameterValues"->Subdivide[0,2\[Pi],50],

"RandomNumberCount"->10,

"MeasurementOperator"->QuantumOperator[{"PauliZ"->{1},"I"->{2}}]

};


(* ::Input::Initialization::Plain:: *)
SPSRGradientValues[generatorFunction_,pauli_,OptionsPattern[]]:=Module[{result,state,vector,rlist,\[Phi]value,rndlen,measurement,\[Theta]values},
	
	\[Theta]values=OptionValue["ParameterValues"];
		
	measurement=OptionValue["MeasurementOperator"]["Matrix"]//Normal;
	
		\[Phi]value=OptionValue["Shift"];
	
		rndlen=OptionValue["RandomNumberCount"];
	
		rlist=RandomReal[1,rndlen];

	result=Table[

				state=QuantumCircuitOperator[{
				
		Exp[I*(1.-s)*QuantumOperator[generatorFunction[\[Theta]]]],
		
				Exp[I*\[Phi]*pauli],
		
				Exp[I*s*QuantumOperator[generatorFunction[\[Theta]]]]
		
				},"Parameters"->{s,\[Phi]}][];
	
				Table[
					vector=state[<|s->sval,\[Phi]->shift|>]["StateVector"];
		
					Re[ConjugateTranspose[vector] . measurement . vector]
		
					,
					{shift,{\[Phi]value,-\[Phi]value}},{sval,rlist}
			],
		
			{\[Theta],\[Theta]values}
		];

		Thread[{\[Theta]values,Mean[(#[[1]]-#[[2]])]&/@result}]
]


(* ::Input::Initialization::Plain:: *)
Options[ASPSRGradientValues]={
"Shift"->\[Pi]/4.,

"ParameterValues"->Subdivide[0,2\[Pi],50],

"RandomNumberCount"->10,

"MeasurementOperator"->QuantumOperator[{"PauliZ"->{1},"I"->{2}}]

};


(* ::Input::Initialization::Plain:: *)
ASPSRGradientValues[generatorFunction_,pauli_,H_,OptionsPattern[]]:=Module[
{result,state,vector,rlist,\[Phi]value,rndlen,measurement,\[Theta]values},

	\[Theta]values=OptionValue["ParameterValues"];

	measurement=OptionValue["MeasurementOperator"]["Matrix"]//Normal;

	\[Phi]value=OptionValue["Shift"];

	rndlen=OptionValue["RandomNumberCount"];

	rlist=RandomReal[1,rndlen];

	result=Table[
				state=QuantumCircuitOperator[{
		
				Exp[I*(1.-s)*QuantumOperator[generatorFunction[\[Theta]]]],
		
				Exp[I*\[Pi]/4*(QuantumOperator[H]+shift*QuantumOperator[pauli])],
		
				Exp[I*s*QuantumOperator[generatorFunction[\[Theta]]]]
		
				},"Parameters"->{s}][];
			
				Table[
				vector=state[<|s->sval|>]["StateVector"];
	
				Re[ConjugateTranspose[vector] . measurement . vector],
	
			{sval,rlist}],

			{\[Theta],\[Theta]values},{shift,{\[Phi]value,-\[Phi]value}}
		];

			Thread[{\[Theta]values,Mean[(#[[1]]-#[[2]])]&/@result}]
]


(* ::Subsection::Closed:: *)
(*Quantum Locking Mechanism*)


QuantumLockingMechanism[key_String]:=Module[{pass,flipsigns, circuit},

	pass=IntegerString[ToCharacterCode[key],4,4];
	
	flipsigns=QuantumOperator[{
								"C"["FlipSign"[#,4]->(Range[StringLength@#]+1),{1}]
							}
				]&/@pass;

	circuit=QuantumCircuitOperator[{
								"H"->{1},
								#,
								"H"->{1},
								QuantumMeasurementOperator["Computational",{1}]
				}]&/@flipsigns

]


QuantumUnlockingMechanism[lock_List,key_String]:=Module[{result,pass,comb},
		
	pass=IntegerString[ToCharacterCode[key],4,4];

	If[!MatchQ[Length@lock,Length@pass],Return["Incorrect key, try again"]];

	result=#[[1]][
				QuantumState[
							<|Prepend[ToExpression[Characters[#[[2]]]],0]->1|>,
							Prepend[ConstantArray[4,StringLength@#[[2]]],2]
							]
				]["Mean"]&/@Transpose[{lock,pass}];

	If[
		MatchQ[result,{1..}],
			"Correct key",
			"Incorrect key, try again"
	]

]


(* ::Section::Closed:: *)
(*Gradient descent functions*)


QuantumNaturalGradientDescent::param = "Metric Tensor parameters ill\[Dash]defined.";

Options[QuantumNaturalGradientDescent]={
	
	"InitialPoint"->Automatic,
	
	"Gradient"->None,
	
	"MaxIterations"->50,
	
	"LearningRate"->0.8
	
};

QuantumNaturalGradientDescent[f_, g_QuantumOperator, OptionsPattern[]]:=Module[{p,gradf,init,steps,\[Eta],MetricTensorFunction},
	
	If[!MatchQ[g["Parameters"],{_Symbol..}],
		Message[QuantumNaturalGradientDescent::param];
		Return[$Failed]
	];
	
	If[VectorQ[OptionValue["InitialPoint"],NumericQ], 
		init=OptionValue["InitialPoint"],
		init=0.01*RandomVariate[NormalDistribution[],Length@g["Parameters"]]
	];
																																																		
	MetricTensorFunction[vect_ /; VectorQ[vect, NumericQ]]:=Normal[N[g[Sequence@@vect]["Matrix"]]];
	
	gradf=OptionValue["Gradient"];
	
	steps=OptionValue["MaxIterations"];
	
	\[Eta]=OptionValue["LearningRate"];
	
	If[
		MatchQ[gradf,None],
			
			NestList[(#-\[Eta] LinearSolve[MetricTensorFunction[#],CheapGradient[f,Table[Symbol["\[Theta]"<>ToString[i]],{i,Length@init}],#]])&,N@init,steps],
			
			NestList[(#-\[Eta] LinearSolve[MetricTensorFunction[#],gradf@@#])&,N@init,steps]
			
	]
]


(* ::Subsection:: *)
(*Auxiliary functions*)


CheapGradient[f_, vars_List, values_ ? VectorQ]:=Module[{permutedVars,nd},

		permutedVars=TakeDrop[#,{1}]&/@NestList[RotateLeft,Thread[vars->values],Length[vars]-1];
		
		centralFiniteDifference[f@@(vars/.#[[2]]),Sequence@@First@#[[1]]]&/@permutedVars
]


centralFiniteDifference[f_[vars__],var_,val_]:=With[{h=10.^-3},

			1/(2 h) ((f[vars]/.var->val+h)-(f[vars]/.var->val-h))
			
		]


Options[GradientDescent]={
	"Gradient"->None,
	
	"MaxIterations"->50,
	
	"LearningRate"->0.8
};
GradientDescent[f_, init_ ? VectorQ, OptionsPattern[]]:=Module[{gradf,steps,\[Eta]},

	gradf=OptionValue["Gradient"];
	
	steps=OptionValue["MaxIterations"];
	
	\[Eta]=OptionValue["LearningRate"];

	If[
		MatchQ[gradf,None],
		
			NestList[(#-\[Eta]*CheapGradient[f,Table[Symbol["\[Theta]"<>ToString[i]],{i,Length@init}],#])&,N@init,steps],
			
			NestList[(#-\[Eta]*gradf@@#)&,N@init,steps](*f->grad[f]*)
			
	]
]


(* ::Section::Closed:: *)
(*QuantumLinearSolver*)


QuantumLinearSolverCircuit[A_?MatrixQ,Ansatz_]:=Module[{pauliDecompose, multiplexer, ancillary,ansatzqubits,qc,parameters,state},

	parameters=Ansatz["Parameters"];
	
	pauliDecompose=QuantumOperator[A]["PauliDecompose"];
	
	multiplexer=QuantumCircuitOperator["Multiplexer"[Sequence@@Keys[pauliDecompose]]];
	
	ancillary=QuantumState[Sqrt@Values[pauliDecompose],"Label"->"\!\(\*SqrtBox[\(m\)]\)"];
	
	ansatzqubits=Range[ancillary["Qudits"]+1,multiplexer["Range"]];
	
	qc=N@QuantumCircuitOperator[
		{
		ancillary,
		Ansatz->ansatzqubits,
		multiplexer,
		SuperDagger[(ancillary["Conjugate"])]
		},
		"Parameters"->parameters
		];
	
	qc
	
	]


QuantumLinearSolve::error = "Global phase deviation exceeds threshold; possible errors."

Options[QuantumLinearSolve]=Join[{"Ansatz"->Automatic,"GlobalPhaseAccuracy"->10^-5},Options[NMinimize]]

QuantumLinearSolve[matrix_?MatrixQ, vector_?VectorQ/;!Mod[Log2[Length[vector]],1]==0, prop : _String | {__String} | All : "Result", opts:OptionsPattern[]]:=Module[
	{result,min,A,b},
	
	
	
	If[Dimensions[matrix]!={#,#}&@Length[vector],Message[QuantumLinearSolve::dim];Return[$Failed]];
	
	Progress`EvaluateWithProgress[
	
		min=Min@Select[Table[2^n,{n,1,10}],#>=Length[vector]&];
		
		A=ReplacePart[PadRight[matrix,{min,min},0],Table[{i,i}->1,{i,Length@vector+1,min}]];
		
		b=Flatten@Append[vector,ConstantArray[0,min-Length@vector]];
	
		,
				
		<|"Text" -> "Translating to \!\(\*SuperscriptBox[\(2\), \(n\)]\)\!\(\*SuperscriptBox[\(x2\), \(\(n\)\(\\\ \)\)]\)problem"|>
	
	];
	
	result=QuantumLinearSolve[A,b,prop,opts];
		
	If[MatchQ[prop,"Result"],result=Take[result,Length@vector]];
	
	If[MatchQ[result,_Association]&&KeyMemberQ[result,"Result"],result["Result"]=Take[result["Result"],Length@vector]];

	result				
]


QuantumLinearSolve::dim = "Vector and matrix dimension not compatible."

QuantumLinearSolve[matrix_?MatrixQ, vector_?VectorQ, prop : _String | {__String} | All : "Result", opts:OptionsPattern[]]:=Module[

	{A,b,cost,bstate,complexQ,\[Omega],v,var,output,plength,cachedResults,
	Ansatz,circuit,parameters,state,globalphase,QuantumDistanceCostFunction,optparameters,result, reporter},
	
	If[Dimensions[matrix]!={#,#}&@Length[vector],Message[QuantumLinearSolve::dim];Return[$Failed]];	
	
	$ModuleNumber=1;
	
	A = matrix/Norm[matrix];
	
	b = Normalize[vector];
																																																																																																																																																																																																									
	bstate = QuantumState[b];
	
	complexQ = !FreeQ[b,_Complex];

	output = Replace[prop, 
					All -> {"Result","Ansatz","CircuitOperator","GlobalPhase","OptimizedParameters","Parameters"}
				];
	
	cachedResults=<||>;
	
	reporter[data_Association] := (
			PrependTo[cachedResults, data];  
			If[ContainsAll[Keys[cachedResults], Flatten[{output}]],
				Return[cachedResults[[output]],Module]
			];
		);
			
 Progress`EvaluateWithProgress[		
	If[
		MatchQ[OptionValue["Ansatz"],Automatic], 
	
			plength = Length@vector;
			$ModuleNumber=1;
			
			If[complexQ,
				\[Omega] = Table[Unique[\[FormalOmega]],{i,2*plength}];
				(*If b is complex, then the QuantumState ansatz is the form of \[Sum]Subscript[\[Alpha], j]+Subscript[I\[Beta], j]|0\[Ellipsis]\[RightAngleBracket])*)		
				Ansatz = QuantumState[(#[[1]]+#[[2]]I)&@Partition[\[Omega],Length[\[Omega]]/2],"Label"->"Ansatz","Parameters"->\[Omega]]
				,
				\[Omega] = Table[Unique[\[FormalOmega]],{i,plength}];
				(*If b is real, then the QuantumState ansatz is the form of \[Sum]Subscript[\[Alpha], j]|0\[Ellipsis]\[RightAngleBracket])*)			
				Ansatz = QuantumState[\[Omega],"Label"->"Ansatz","Parameters"->\[Omega]]
			];
	
			, 
			
			Ansatz = OptionValue["Ansatz"];
			If[MatchQ[Ansatz,_QuantumCircuitOperator],Ansatz = Ansatz[]];
			\[Omega] = Ansatz["Parameters"];
			
		];
		,		
		<|"Text" -> "Ansatz initialization"|>
	];

	parameters=Ansatz["Parameters"];
	
	reporter[<|"Ansatz"->Ansatz,"Parameters"->parameters|>];


 Progress`EvaluateWithProgress[
 		 
		circuit = QuantumLinearSolverCircuit[A,Ansatz];
			
		state = circuit[];
		
	,
	
	<|"Text" -> "Variational circuit initialization", "ElapsedTime" -> Automatic|>
	
	];

	reporter[<|"CircuitOperator"->circuit|>];
	
	QuantumDistanceCostFunction[vect_ /; VectorQ[vect, RealValuedNumberQ]]:=QuantumDistance[bstate["Normalized"],state[AssociationThread[parameters->vect]]["Normalized"],"Fidelity"];


	
DynamicModule[{points},
	points={};
	Progress`EvaluateWithProgress[
		optparameters =
			Last[
			NMinimize[
				QuantumDistanceCostFunction[\[Omega]],
				\[Omega],
				FilterRules[{opts}, Options[NMinimize]],
				StepMonitor:>(AppendTo[points,cost=QuantumDistanceCostFunction[\[Omega]]])
				]
			],
			
			If[
				MatchQ[OptionValue[Method],Automatic],
				
				<|"Text" -> "Hybrid optimization procedure","ElapsedTime"->Automatic|>,
								
				<|"Text" -> Dynamic[
							ListLogPlot[points,
									PlotLabel->"Hybrid optimization procedure",
									FrameLabel->{"Step","Log[Cost function]"},
									PlotRange->{0,Log@First[points]},
									GridLines->Automatic,
									Frame->True,
									ImageSize->Small, 
									AspectRatio->1,
									PlotLegends->StringForm["Cost function\nminimization value:\n``",ScientificForm[Last[points]]],
									Joined->True
									]
								]								
				|>		
				
			]
		]
	];


	reporter[<|"OptimizedParameters"->optparameters|>];

 Progress`EvaluateWithProgress[

		globalphase=Divide@@(Extract[#,SparseArray[Chop[b]]["ExplicitPositions"]]&/@{state[<|optparameters|>]["AmplitudesList"],b});
	
		If[StandardDeviation[globalphase] > OptionValue["GlobalPhaseAccuracy"], Message[QuantumLinearSolve::error]];
	
		result=Ansatz[<|optparameters|>]["AmplitudesList"];
	
		result=1/Mean[globalphase]*result*Norm[vector]/Norm[matrix];

		,
		
		<|"Text" -> "Final normalization"|>
	
	];

	reporter[<|"Result"->result,"GlobalPhase"->Around[Mean[globalphase],StandardDeviation[globalphase]]|>]

]
