(* Shared SDE-reproduction toolkit for the Cartesian-Bloch continuously-monitored-
   qubit trajectory SDE (dispersive-readout stochastic master equation, Gambetta07).
   Rates are {gci, gba, gd, gphi, g1, wx}; dep = g1/2 + gphi + gd = 1/T2eff. *)

depOf[g1_, gphi_, gd_] := g1/2 + gphi + gd;

driftV[{x_, y_, z_}, {gci_, gba_, gd_, gphi_, g1_, wx_}] :=
  With[{dep = depOf[g1, gphi, gd]}, {-dep x, -dep y - wx z, g1 (1 - z) + wx y}];

diffV[{x_, y_, z_}, {gci_, gba_, gd_, gphi_, g1_, wx_}] :=
  {-Sqrt[gci] x z - Sqrt[gba] y, -Sqrt[gci] y z + Sqrt[gba] x, Sqrt[gci] (1 - z^2)};

(* -- 2x2 Lindblad generator, factored once:  H = (wx/2) sigma_x, relaxation jump
      sqrt(g1)|g><e|, sigma_z dephasing rate (gphi+gd)/2.  Returns dρ/dt. -- *)
lindblad2x2Gen[rho_, {gci_, gba_, gd_, gphi_, g1_, wx_}] :=
  With[{sx = PauliMatrix[1], sz = PauliMatrix[3], sge = {{0, 1}, {0, 0}},
        dop = Function[{L, rr}, L . rr . ConjugateTranspose[L] -
           (ConjugateTranspose[L] . L . rr + rr . ConjugateTranspose[L] . L)/2]},
    (* (wx/2) multiplies the parenthesized commutator: Dot binds tighter than Times *)
    -I (wx/2) (sx . rho - rho . sx) + g1 dop[sge, rho] + ((gphi + gd)/2) dop[sz, rho]];

(* Unconditional Bloch rates {d<sx>,d<sy>,d<sz>}/dt from the 2x2 generator. *)
lindblad2x2Bloch[{x_, y_, z_}, rates_] :=
  With[{drho = lindblad2x2Gen[(IdentityMatrix[2] + x PauliMatrix[1] +
        y PauliMatrix[2] + z PauliMatrix[3])/2, rates]},
    Table[Tr[drho . PauliMatrix[k]], {k, 3}]];

(* -- Native ItoProcess + RandomFunction (the TransmonLab engine) -- *)
itoEnsemble[rates_, {x0_, y0_, z0_}, tf_, dt_, ntraj_, seed_] :=
  Block[{proc, td, vals, z, r},
    proc = ItoProcess[
      {driftV[{\[FormalX][t], \[FormalY][t], \[FormalZ][t]}, rates],
       List /@ diffV[{\[FormalX][t], \[FormalY][t], \[FormalZ][t]}, rates],
       {\[FormalX][t], \[FormalY][t], \[FormalZ][t]}},
      {{\[FormalX], \[FormalY], \[FormalZ]}, {x0, y0, z0}}, t];
    td = BlockRandom[SeedRandom[seed];
      RandomFunction[proc, {0., tf, dt}, ntraj,
        Method -> "StochasticRungeKuttaScalarNoise"]];
    vals = td["ValueList"];                (* ntraj x npts x 3 *)
    z = vals[[All, All, 3]];
    r = Map[Norm, vals, {2}];              (* Bloch radius per traj per step *)
    <|"meanZ" -> Mean[z], "times" -> td["Times"], "maxr" -> Max[r],
      "finalR" -> r[[All, -1]], "zTrajs" -> z|>
  ];

(* -- Hand-rolled Euler-Maruyama inside NDSolve via WhenEvent, one path.
      Module localizes x,y,z: they are NDSolve dependent variables and must be
      fresh symbols (the symbol-renaming case Module exists for).
      Frozen ODE (x'=0); increment ADDED at each dt; optional r<=1 projection.
      Midpoint phase Mod[t+dt/2,dt]==0 => exactly N fires, clean grid sampling. -- *)
emPath[rates_, {x0_, y0_, z0_}, tf_, dt_, project_] :=
  Module[{x, y, z, drift1, diff1, sol, grid, states},
    drift1 = driftV[{x[t], y[t], z[t]}, rates];
    diff1  = diffV[{x[t], y[t], z[t]}, rates];
    sol = NDSolve[{
       x'[t] == 0, y'[t] == 0, z'[t] == 0,
       WhenEvent[Mod[t + dt/2, dt] == 0,
        With[{dw = RandomVariate[NormalDistribution[0, Sqrt[dt]]]},
         With[{new = {x[t], y[t], z[t]} + drift1 dt + diff1 dw},
          With[{rr = Norm[new]},
           {x[t], y[t], z[t]} -> If[project && rr > 1, new/rr, new]]]]],
       x[0] == x0, y[0] == y0, z[0] == z0},
      {x, y, z}, {t, 0, tf}, MaxStepSize -> dt][[1]];
    grid = Range[0., tf, dt];
    states = Map[{x[#], y[#], z[#]} /. sol &, grid];
    {states[[All, 3]], Norm /@ states}     (* {zgrid, rgrid} *)
  ];

emEnsemble[rates_, init_, tf_, dt_, ntraj_, seed_, project_] :=
  Block[{paths, zs, rs},
    paths = BlockRandom[SeedRandom[seed];
      Table[emPath[rates, init, tf, dt, project], {ntraj}]];
    zs = paths[[All, 1]]; rs = paths[[All, 2]];
    <|"meanZ" -> Mean[zs], "maxr" -> Max[rs],
      "times" -> Range[0., tf, dt], "zTrajs" -> zs|>
  ];

(* -- Vectorized EM: ONE NDSolve for the whole ensemble.  Module localizes the
      dependent-variable vectors xs,ys,zs (symbol renaming, as Module is for);
      the scalar rates come in through a With that substitutes them into the
      held NDSolve.  One WhenEvent draws ntraj increments and (optionally)
      projects each path onto r<=1 via Clip[r,{1,Inf}] as a per-path divisor. -- *)
emEnsembleVec[rates_, {x0_, y0_, z0_}, tf_, dt_, ntraj_, seed_, project_] :=
  Module[{xs, ys, zs, sol, grid, xyzg, zg, rmat, tmax, cap = 10.^8},
    sol = BlockRandom[SeedRandom[seed];
      NDSolve[{
        xs'[t] == ConstantArray[0., ntraj],
        ys'[t] == ConstantArray[0., ntraj],
        zs'[t] == ConstantArray[0., ntraj],
        WhenEvent[Mod[t + dt/2, dt] == 0,
         Block[{dw = RandomVariate[NormalDistribution[0, Sqrt[dt]], ntraj], xyz, new, capr},
          xyz = {xs[t], ys[t], zs[t]};
          (* reuse the canonical driftV/diffV -- one source for the SDE, vectorized
             over the ntraj-length component vectors (dw is the shared noise) *)
          new = xyz + MapThread[#1 dt + #2 dw &, {driftV[xyz, rates], diffV[xyz, rates]}];
          capr = If[project, Clip[Sqrt[new[[1]]^2 + new[[2]]^2 + new[[3]]^2], {1., Infinity}],
                    ConstantArray[1., ntraj]];
          {xs[t], ys[t], zs[t]} -> (#/capr & /@ new)]],   (* divide each component by capr *)
        xs[0] == ConstantArray[N@x0, ntraj],
        ys[0] == ConstantArray[N@y0, ntraj],
        zs[0] == ConstantArray[N@z0, ntraj]},
       {xs, ys, zs}, {t, 0, tf}, MaxStepSize -> dt][[1]]];
    (* Clip sampling to the solution's valid domain (noProj may truncate early). *)
    tmax = (zs /. sol)["Domain"][[1, 2]];
    grid = Select[Range[0., tf, dt], # <= tmax + dt/4 &];
    xyzg = ({xs[#], ys[#], zs[#]} /. sol) & /@ grid;   (* npts x 3 x ntraj *)
    zg = xyzg[[All, 3]];
    rmat = Replace[
      Flatten[Sqrt[xyzg[[All, 1]]^2 + xyzg[[All, 2]]^2 + xyzg[[All, 3]]^2]],
      v_ /; (! NumberQ[v] || Abs[v] > cap) :> cap, {1}];
    <|"meanZ" -> Mean /@ zg, "maxr" -> Max[rmat],
      "diverged" -> (tmax < tf - dt/2 || Max[rmat] >= cap),
      "times" -> grid, "zTrajs" -> Transpose[zg]|>
  ];

(* -- Lindblad reference: the affine Bloch ODE CLOSES, so DSolve (exact), not
      NDSolve.  Returns the closed-form {<x>,<y>,<z>}(t) as pure functions. -- *)
lindbladBloch[rates_, {x0_, y0_, z0_}] :=
  Block[{gci, gba, gd, gphi, g1, wx, dep, sol},
    {gci, gba, gd, gphi, g1, wx} = rates; dep = depOf[g1, gphi, gd];
    sol = DSolve[{
       \[FormalX]'[t] == -dep \[FormalX][t],
       \[FormalY]'[t] == -dep \[FormalY][t] - wx \[FormalZ][t],
       \[FormalZ]'[t] == g1 (1 - \[FormalZ][t]) + wx \[FormalY][t],
       \[FormalX][0] == x0, \[FormalY][0] == y0, \[FormalZ][0] == z0},
      {\[FormalX], \[FormalY], \[FormalZ]}, t] // First;
    {\[FormalX], \[FormalY], \[FormalZ]} /. sol
  ];

lindbladZfun[rates_, init_, tf_ : Automatic] :=
  With[{fz = lindbladBloch[rates, init][[3]]}, Function[tt, Re[fz[tt]]]];

(* -- Independent machinery: 2x2 density-matrix Lindblad -> <sigma_z>(t) via the
      EXACT Liouvillian propagator rho(t) = exp(L t) rho(0) (no ODE integration).
      L is the 4x4 matrix of the factored generator in the vec(rho) basis, so the
      2x2 model lives in exactly one place, and the Hilbert-space evolution is a
      genuinely different route from the affine Bloch DSolve. -- *)
lindblad2x2Zfun[rates_, {x0_, y0_, z0_}, tf_ : Automatic] :=
  With[{basis = {{{1, 0}, {0, 0}}, {{0, 1}, {0, 0}}, {{0, 0}, {1, 0}}, {{0, 0}, {0, 1}}}},
    With[{Lmat = Transpose[Flatten[lindblad2x2Gen[#, rates]] & /@ basis],
          vec0 = {(1 + z0)/2, (x0 - I y0)/2, (x0 + I y0)/2, (1 - z0)/2}},
      Function[tt, With[{v = MatrixExp[Lmat tt] . vec0}, Re[v[[1]] - v[[4]]]]]]];  (* rho_gg - rho_ee *)
