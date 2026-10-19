(* Execute the unchanged repository suite. "candidate" changes only a private
   definition in this fresh kernel; it never edits the implementation file. *)
$AllowInternet=False;
$PermanentAuditDirectory=DirectoryName[$InputFileName];
$PermanentAuditRoot=Nest[DirectoryName,$PermanentAuditDirectory,3];
$PermanentAuditCandidate=MemberQ[$ScriptCommandLine,"candidate"];
PacletDirectoryLoad/@{
 "/Users/mohammadb/Library/Wolfram/Paclets/Repository/Wolfram__Arrays-1.4.1",
 "/Users/mohammadb/Library/Wolfram/Paclets/Repository/Wolfram__TensorNetworks-1.2.0",
 "/Users/mohammadb/Library/Wolfram/Paclets/Repository/IBMQuantumPlatform-0.0.4"
};
If[$PermanentAuditCandidate,
 PacletDirectoryLoad[FileNameJoin[{$PermanentAuditRoot,"QuantumFramework"}]];
 Needs["Wolfram`QuantumFramework`"];Needs["Wolfram`QuantumFramework`SecondQuantization`"];
 Get[FileNameJoin[{$PermanentAuditDirectory,"Proposal.wl"}]];
 PermanentSelectionAudit`InstallProposal[];
 $Epilog:=Echo[With[{f=ToExpression[First@Names["Wolfram`QuantumFramework`SecondQuantization`*`permanentAmplitude"]]},!FreeQ[DownValues[f],PermanentSelectionAudit`ProposedAmplitude]],"PROPOSAL_STILL_INSTALLED"]
];
Block[{$ScriptCommandLine={"Tests/RunTests.wls","QuantumInterferometer"}},
 Get[FileNameJoin[{$PermanentAuditRoot,"Tests","RunTests.wls"}]]
];
