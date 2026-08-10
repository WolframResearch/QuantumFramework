(* Cost profile of the 2D/3D conditional step on this machine.

   Run:  wolframscript -file benchmarks.wls
   ClearSystemCache[] precedes every timing (without it, cached FFT plans and results make
   second runs look far faster than reality). "cold step" is one full conditional step
   (outcome draw + Kraus update + Strang drift) timed once from a cold cache; "warm step"
   is the per-step average over a short Nest trajectory, the number a long run actually
   pays; the component lines are min-of-3, each from a cold cache. *)

Get[FileNameJoin[{DirectoryName[$InputFileName], "engine-2d3d.wl"}]];
w[s_] := WriteString[$Output, s, "\n"];
rep[label_, val_] := w[label <> " :: " <> ToString[val, InputForm]];

SetAttributes[minOf3, HoldFirst];
minOf3[expr_] := Min[Table[ClearSystemCache[]; First[AbsoluteTiming[expr]], {3}]];

bench2[{nx_, ny_}, ntW_] := Block[{g, ch, psi, dt = 0.002, tCold, tWarm},
   g = grid2[{nx, ny}, {nx/8., ny/8.}];
   ch = axisChannels2[g, {1., 1.}];
   psi = gaussian2[g, {0., 0.}, {0., 0.}, {1., 1.}];
   ClearSystemCache[];
   tCold = First@AbsoluteTiming[BlockRandom[SeedRandom[1];
       step2[g, 0 &, ch, dt, drawOutcomes2[g, ch, dt][psi]][psi]]];
   ClearSystemCache[];
   tWarm = First[AbsoluteTiming[BlockRandom[SeedRandom[2];
        Nest[Function[p, step2[g, 0 &, ch, dt, drawOutcomes2[g, ch, dt][p]][p]], psi, ntW]]]]/ntW;
   <|"grid" -> {nx, ny}, "stateMB" -> N[ByteCount[psi]/2^20],
     "coldStep" -> tCold, "warmStep" -> tWarm,
     "ftPair" -> minOf3[ift[ft[psi]]],
     "measUpdate" -> minOf3[measUpdate2[g, ch, dt, {0.1, -0.1}][psi]],
     "jointDraw" -> minOf3[BlockRandom[SeedRandom[3]; jointSample2[g][psi]]]|>];

bench3[{nx_, ny_, nz_}, ntW_] := Block[{g, ch, psi, dt = 0.002, tCold, tWarm},
   g = grid3[{nx, ny, nz}, {nx/8., ny/8., nz/8.}];
   ch = axisChannels3[g, {1., 1., 1.}];
   psi = gaussian3[g, {0., 0., 0.}, {0., 0., 0.}, {1., 1., 1.}];
   ClearSystemCache[];
   tCold = First@AbsoluteTiming[BlockRandom[SeedRandom[1];
       step3[g, 0 &, ch, dt, drawOutcomes3[g, ch, dt][psi]][psi]]];
   ClearSystemCache[];
   tWarm = First[AbsoluteTiming[BlockRandom[SeedRandom[2];
        Nest[Function[p, step3[g, 0 &, ch, dt, drawOutcomes3[g, ch, dt][p]][p]], psi, ntW]]]]/ntW;
   <|"grid" -> {nx, ny, nz}, "stateMB" -> N[ByteCount[psi]/2^20],
     "coldStep" -> tCold, "warmStep" -> tWarm,
     "ftPair" -> minOf3[ift[ft[psi]]],
     "measUpdate" -> minOf3[measUpdate3[g, ch, dt, {0.1, -0.1, 0.05}][psi]],
     "jointDraw" -> minOf3[BlockRandom[SeedRandom[3]; jointSample3[g][psi]]]|>];

rep["2D 128^2", bench2[{128, 128}, 40]]
rep["2D 256^2", bench2[{256, 256}, 20]]
rep["3D 64^3", bench3[{64, 64, 64}, 10]]
rep["3D 128^3", bench3[{128, 128, 128}, 6]]

(* The density-matrix route scales as N^2d in memory and N^3d per dense step: the numbers,
   not attempted beyond the 32^2 Lindblad check in validation.wls. *)
rhoLine[n_List] := With[{nn = Times @@ n},
   {n, "rho GB" -> N[16 nn^2/2^30], "dense step Tflop" -> N[8 nn^3/10^12]}];
rep["rho arithmetic", Map[rhoLine, {{128, 128}, {256, 256}, {64, 64, 64}, {128, 128, 128}}]]

rep["MaxMemoryUsed MB (this whole benchmark run)", N[MaxMemoryUsed[]/2^20]]
