<!-- #| style: Abstract -->
In this simple computational essay, I investigate the dynamics of a Bose-Einstein condensate (BEC) governed by the Gross-Pitaevskii equation under the influence of a time-dependent perturbation. The perturbation consists of an attractive localized potential traversing the condensate along various trajectories. Numerical simulations of the dynamical equation reveal the emergence of some sort of [wake](https://en.wikipedia.org/wiki/Wake_(physics)) patterns in response to the perturbation. While this is a simple study, I believe it hints at deeper concepts and applications yet to be explored.

This exploration was inspired by discussions with Sir Michael Berry during the Wolfram Innovator Award 2025 meetings. Also, Michael Trott suggested to examine trajectories beyond straight lines.

## Line

Trajectory of perturbation:

```wl
With[{initial = 
   DensityPlot[Exp[-(x^2 + y^2)], {x, -15, 15}, {y, -10, 10}, 
    PlotRange -> All, PlotPoints -> 50, 
    ColorFunction -> "SunsetColors"]},
 Manipulate[
  Show[initial, 
   Quiet@ContourPlot[
     Exp[-(x - 7 \[Tau] + 7)^2 - (y)^2], {x, -15, 15}, {y, -10, 10}, 
     PlotRange -> All, ContourShading -> None, PlotPoints -> 50, 
     ContourStyle -> Red], 
   ParametricPlot[{7 u - 7, 0}, {u, 0, 3}, 
    PlotStyle -> Directive[Dotted, Green]], 
   PlotLabel -> "Perturbation wrt the initial wavefunction"], {\[Tau],
    0, 3}, SaveDefinitions -> True]]
```

Hamiltonian:

```wl
\[ScriptCapitalH] = 
  SchrodingerPDEComponent[{\[CapitalPsi][t, x, y], t, {x, y}}, <|
    "ReducedPlanckConstant" -> 1, 
    "SchrodingerPotential" -> 
     3 \[CapitalPsi][t, x, 
        y] Conjugate[\[CapitalPsi][t, x, y]] + .5 (x^2 + y^2) - 
      29 Exp[-(x - 7 t + 7)^2 - (y)^2]|>];
```

Initial state:

```wl
ic = \[CapitalPsi][0, x, y] == Exp[-(x^2 + y^2)];
```

Solve dynamical equation numerically:

```wl
sol = NDSolveValue[{\[ScriptCapitalH] == 0, ic}, \[CapitalPsi], {t, 0,
      3}, {x, -20, 20}, {y, -10, 10}, 
    Method -> {"PDEDiscretization" -> {"MethodOfLines", 
        "SpatialDiscretization" -> {"FiniteElement", 
          "MeshOptions" -> 
           "MaxCellMeasure" -> .2}}}]; // AbsoluteTiming
```

Visualize results:

```wl
AnimationVideo[
 Show[DensityPlot[Abs[sol[t, x, y]], {x, -20, 20}, {y, -10, 10}, 
   PlotPoints -> 150, 
   ColorFunction -> (ColorData["SunsetColors"][Sqrt[#]] &), 
   PlotRange -> All], 
  Graphics[{Green, PointSize[Medium], Point[{7 t - 7, 0}]}]], {t, 0, 
  3}]
```

Get a snapshots of evolution:

```wl
imgs = Show[
     DensityPlot[Abs[sol[#, x, y]], {x, -20, 20}, {y, -10, 10}, 
      PlotPoints -> 150, 
      ColorFunction -> (ColorData["SunsetColors"][Sqrt[#]] &), 
      PlotRange -> All], 
     Graphics[{Red, PointSize[Medium], Point[{7 # - 7, 0}]}]] & /@ 
   Range[0, 3, .05];
```

Export as GIF:

```wl
Export[NotebookDirectory[] <> "Line.gif", imgs, 
 "DisplayDurations" -> .3, AnimationRepetitions -> 100]
```

## Circle

Radius:

```wl
R = 1;
```

Trajectory of perturbation:

```wl
With[{initial = 
   DensityPlot[Exp[-(x^2 + y^2)], {x, -3, 3}, {y, -3, 3}, 
    PlotRange -> All, PlotPoints -> 50, 
    ColorFunction -> "SunsetColors"]}, 
 Manipulate[
  Show[initial, 
   Quiet@ContourPlot[
     Exp[-(x - R Cos[2 \[Pi] \[Tau]])^2 - (y - 
          R Sin[2 \[Pi] \[Tau]])^2], {x, -3, 3}, {y, -3, 3}, 
     PlotRange -> All, ContourShading -> None, PlotPoints -> 100, 
     ContourStyle -> Red], 
   ParametricPlot[R {Cos[2 \[Pi] u], Sin[2 \[Pi] u]}, {u, 0, 1}, 
    PlotStyle -> Directive[Dotted, Green]], 
   PlotLabel -> "Perturbation wrt the initial wavefunction"], {\[Tau],
    0, 3}, SaveDefinitions -> True]]
```

Hamiltonian:

```wl
\[ScriptCapitalH] = 
  SchrodingerPDEComponent[{\[CapitalPsi][t, x, y], t, {x, y}}, <|
    "ReducedPlanckConstant" -> 1, 
    "SchrodingerPotential" -> 
     3 \[CapitalPsi][t, x, 
        y] Conjugate[\[CapitalPsi][t, x, y]] + .5 (x^2 + y^2) - 
      29 Exp[-(x - R Cos[2 \[Pi] t])^2 - (y - R Sin[2 \[Pi] t])^2]|>];
```

Initial state:

```wl
ic = \[CapitalPsi][0, x, y] == Exp[-(x^2 + y^2)];
```

Solve the dynamical equation numerically:

```wl
sol = NDSolveValue[{\[ScriptCapitalH] == 0, ic}, \[CapitalPsi], {t, 0,
      2}, {x, -15, 15}, {y, -15, 15}, 
    Method -> {"PDEDiscretization" -> {"MethodOfLines", 
        "SpatialDiscretization" -> {"FiniteElement", 
          "MeshOptions" -> 
           "MaxCellMeasure" -> .1}}}]; // AbsoluteTiming
```

Visualize results:

```wl
AnimationVideo[
 Show[DensityPlot[Abs[sol[t, x, y]], {x, -10, 10}, {y, -10, 10}, 
   PlotPoints -> 150, 
   ColorFunction -> (ColorData["SunsetColors"][Sqrt[#]] &), 
   PlotRange -> All],
  Graphics[{Green, Point[R {Cos[2 \[Pi] t], Sin[2 \[Pi] t]}]}]], {t, 
  0, 2}]
```

A series of snapshots:

```wl
imgs = Table[
   Show[DensityPlot[Abs[sol[t, x, y]], {x, -10, 10}, {y, -10, 10}, 
     PlotPoints -> 150, 
     ColorFunction -> (ColorData["SunsetColors"][Sqrt[#]] &), 
     PlotRange -> All],
    Graphics[{Green, Point[R {Cos[2 \[Pi] t], Sin[2 \[Pi] t]}]}]], {t,
     0, 2, .05}];
```

Export as GIF:

```wl
Export[NotebookDirectory[] <> "Circle.gif", imgs, 
 "DisplayDurations" -> .3, AnimationRepetitions -> 100]
```

## Spiral

Trajectory of perturbation

```wl
With[{initial = 
   DensityPlot[Exp[-(x^2 + y^2)], {x, -3, 3}, {y, -3, 3}, 
    PlotRange -> All, PlotPoints -> 50, 
    ColorFunction -> "SunsetColors"]}, 
 Manipulate[
  Show[initial, 
   Quiet@ContourPlot[
     Exp[-(x - \[Tau] Cos[2 \[Pi] \[Tau]])^2 - (y - \[Tau] Sin[
            2 \[Pi] \[Tau]])^2], {x, -3, 3}, {y, -3, 3}, 
     PlotRange -> All, ContourShading -> None, PlotPoints -> 100, 
     ContourStyle -> Red], 
   ParametricPlot[{u Cos[2 \[Pi] u], u Sin[2 \[Pi] u]}, {u, 0, 2}, 
    PlotStyle -> Directive[Dotted, Green]], 
   PlotLabel -> "Perturbation wrt the initial wavefunction"], {\[Tau],
    0, 2}, SaveDefinitions -> True]]
```

Hamiltonian:

```wl
\[ScriptCapitalH] = 
  SchrodingerPDEComponent[{\[CapitalPsi][t, x, y], t, {x, y}}, <|
    "ReducedPlanckConstant" -> 1, 
    "SchrodingerPotential" -> 
     3 \[CapitalPsi][t, x, 
        y] Conjugate[\[CapitalPsi][t, x, y]] + .5 (x^2 + y^2) - 
      29 Exp[-(x - t Cos[2 \[Pi] t])^2 - (y - t Sin[2 \[Pi] t])^2]|>];
```

Initial state:

```wl
ic = \[CapitalPsi][0, x, y] == Exp[-(x^2 + y^2)];
```

Solve the dynamical equation numerically:

```wl
sol = NDSolveValue[{\[ScriptCapitalH] == 0, ic}, \[CapitalPsi], {t, 0,
      2}, {x, -15, 15}, {y, -15, 15}, 
    Method -> {"PDEDiscretization" -> {"MethodOfLines", 
        "SpatialDiscretization" -> {"FiniteElement", 
          "MeshOptions" -> 
           "MaxCellMeasure" -> .1}}}]; // AbsoluteTiming
```

Visualize results:

```wl
AnimationVideo[
 Show[DensityPlot[Abs[sol[t, x, y]], {x, -10, 10}, {y, -10, 10}, 
   PlotPoints -> 150, 
   ColorFunction -> (ColorData["SunsetColors"][Sqrt[#]] &), 
   PlotRange -> All],
  Graphics[{Green, Point[{t Cos[2 \[Pi] t], t Sin[2 \[Pi] t]}]}]], {t,
   0, 2}]
```

A series of snapshots:

```wl
imgs = Table[
   Show[DensityPlot[Abs[sol[t, x, y]], {x, -10, 10}, {y, -10, 10}, 
     PlotPoints -> 150, 
     ColorFunction -> (ColorData["SunsetColors"][Sqrt[#]] &), 
     PlotRange -> All],
    Graphics[{Green, 
      Point[{t Cos[2 \[Pi] t], t Sin[2 \[Pi] t]}]}]], {t, 0, 2, .05}];
```

Export as GIF:

```wl
Export[NotebookDirectory[] <> "Spiral.gif", imgs, 
 "DisplayDurations" -> .3, AnimationRepetitions -> 100]
```

## Bowtie

Trajectory of perturbation

```wl
With[{initial = 
   DensityPlot[Exp[-(x^2 + y^2)], {x, -3, 3}, {y, -3, 3}, 
    PlotRange -> All, PlotPoints -> 50, 
    ColorFunction -> "SunsetColors"]}, 
 Manipulate[
  Show[initial, 
   Quiet@ContourPlot[
     Exp[-(x - 1.5 Sin[2 \[Pi] \[Tau]])^2 - (y - 
          1.5 Sin[4 \[Pi] \[Tau]])^2], {x, -3, 3}, {y, -3, 3}, 
     PlotRange -> All, ContourShading -> None, PlotPoints -> 100, 
     ContourStyle -> Red], 
   ParametricPlot[1.5 { Sin[2 \[Pi] u], Sin[4 \[Pi] u]}, {u, 0, 2}, 
    PlotStyle -> Directive[Dotted, Green]], 
   PlotLabel -> "Perturbation wrt the initial wavefunction"], {\[Tau],
    0, 2}, SaveDefinitions -> True]]
```

Hamiltonian:

```wl
\[ScriptCapitalH] = 
  SchrodingerPDEComponent[{\[CapitalPsi][t, x, y], t, {x, y}}, <|
    "ReducedPlanckConstant" -> 1, 
    "SchrodingerPotential" -> 
     3 \[CapitalPsi][t, x, 
        y] Conjugate[\[CapitalPsi][t, x, y]] + .5 (x^2 + y^2) - 
      29 Exp[-(x - 1.5 Sin[2 \[Pi] t])^2 - (y - 1.5 Sin[4 \[Pi] t])^2]|>];
```

Initial state:

```wl
ic = \[CapitalPsi][0, x, y] == Exp[-(x^2 + y^2)];
```

Solve the dynamical equation numerically:

```wl
sol = NDSolveValue[{\[ScriptCapitalH] == 0, ic}, \[CapitalPsi], {t, 0,
      2}, {x, -15, 15}, {y, -15, 15}, 
    Method -> {"PDEDiscretization" -> {"MethodOfLines", 
        "SpatialDiscretization" -> {"FiniteElement", 
          "MeshOptions" -> 
           "MaxCellMeasure" -> .1}}}]; // AbsoluteTiming
```

Visualize results:

```wl
AnimationVideo[
 Show[DensityPlot[Abs[sol[t, x, y]], {x, -10, 10}, {y, -10, 10}, 
   PlotPoints -> 150, 
   ColorFunction -> (ColorData["SunsetColors"][Sqrt[#]] &), 
   PlotRange -> All],
  Graphics[{Green, Point[{Sin[2 \[Pi] t], Sin[2 \[Pi] 2 t]}]}]], {t, 
  0, 2}]
```

A series of snapshots:

```wl
imgs = Table[
   Show[DensityPlot[Abs[sol[t, x, y]], {x, -10, 10}, {y, -10, 10}, 
     PlotPoints -> 150, 
     ColorFunction -> (ColorData["SunsetColors"][Sqrt[#]] &), 
     PlotRange -> All],
    Graphics[{Green, Point[{Sin[2 \[Pi] t], Sin[2 \[Pi] 2 t]}]}]], {t,
     0, 2, .05}];
```

Export as GIF:

```wl
Export[NotebookDirectory[] <> "BowTie.gif", imgs, 
 "DisplayDurations" -> .3, AnimationRepetitions -> 100]
```

## No movement, but time-dependent amplitude of perturbation

Hamiltonian:

```wl
\[ScriptCapitalH] = 
  SchrodingerPDEComponent[{\[CapitalPsi][t, x, y], t, {x, y}}, <|
    "ReducedPlanckConstant" -> 1, 
    "SchrodingerPotential" -> 
     3 \[CapitalPsi][t, x, 
        y] Conjugate[\[CapitalPsi][t, x, y]] + .5 (x^2 + y^2) - 
      29 Sin[2 \[Pi] t] Exp[-(x)^2 - (y)^2]|>];
```

Initial state:

```wl
ic = \[CapitalPsi][0, x, y] == Exp[-(x^2 + y^2)];
```

Solve dynamical equation numerically:

```wl
sol = NDSolveValue[{\[ScriptCapitalH] == 0, ic}, \[CapitalPsi], {t, 0,
      2}, {x, -15, 15}, {y, -15, 15}, 
    Method -> {"PDEDiscretization" -> {"MethodOfLines", 
        "SpatialDiscretization" -> {"FiniteElement", 
          "MeshOptions" -> 
           "MaxCellMeasure" -> .1}}}]; // AbsoluteTiming
```

Visualize results:

```wl
ListAnimate[
 Table[DensityPlot[Abs[sol[t, x, y]], {x, -12, 12}, {y, -12, 12}, 
   PlotPoints -> 150, 
   ColorFunction -> (ColorData["SunsetColors"][Sqrt[#]] &), 
   PlotRange -> All], {t, 0, 2, .1}], SaveDefinitions -> True]
```

Get a snapshots of evolution:

```wl
imgs = DensityPlot[Abs[sol[#, x, y]], {x, -12, 12}, {y, -12, 12}, 
     PlotPoints -> 150, 
     ColorFunction -> (ColorData["SunsetColors"][Sqrt[#]] &), 
     PlotRange -> All] & /@ Range[0, 2, .1];
```

Export as GIF:

```wl
Export[NotebookDirectory[] <> "TimeVaryingAmplitude.gif", imgs, 
 "DisplayDurations" -> .3, AnimationRepetitions -> 100]
```
