$ProgressReporting=False;
With[{report=TestReport[FileNameJoin[{DirectoryName[$InputFileName],"ProposalTests.wlt"}]]},
 Echo[{report["TestsSucceededCount"],report["TestsFailedCount"]},"PROPOSAL_TEST_COUNTS"];
 Echo[Map[<|"TestID"->#["TestID"],"Outcome"->#["Outcome"]|>&,report["Results"]],"PROPOSAL_TEST_RESULTS"];
 Export[FileNameJoin[{DirectoryName[$InputFileName],"proposal-test-results.json"}],Map[<|"TestID"->#["TestID"],"Outcome"->#["Outcome"]|>&,report["Results"]],"RawJSON"];
 Exit[If[report["TestsFailedCount"]===0,0,1]]
];
