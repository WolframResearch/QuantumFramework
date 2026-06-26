# NDSolveValue | [SpanFromLeft]

> [NDSolveValue](https://reference.wolfram.com/language/ref/NDSolveValue.html)[*eqns*,*expr*,{*x*,*x*_*min*,*x*_*max*}] — gives the value of `*expr*` with functions determined by a numerical solution to the ordinary differential equations `*eqns*` with the independent variable `*x*` in the range `*x*_*min*` to `*x*_*max*`.
> [NDSolveValue](https://reference.wolfram.com/language/ref/NDSolveValue.html)[*eqns*,*expr*,{*x*,*x*_*min*,*x*_*max*},{*y*,*y*_*min*,*y*_*max*}] — solves the partial differential equations `*eqns*` over a rectangular region.
> [NDSolveValue](https://reference.wolfram.com/language/ref/NDSolveValue.html)[*eqns*,*expr*,{*x*,*y*}∈Ω] — solves the partial differential equations `*eqns*` over the region `Ω`.
> [NDSolveValue](https://reference.wolfram.com/language/ref/NDSolveValue.html)[*eqns*,*u*,{*t*,*t*_*min*,*t*_*max*},{*x*,*y*}∈Ω] — solves the time-dependent partial differential equations `*eqns*` over the region `Ω`.

## Details and Options

[NDSolveValue](https://reference.wolfram.com/language/ref/NDSolveValue.html)[*eqns*,*y*[*x*],{*x*,*x*_*min*,*x*_*max*}] gives solutions for `*y*[*x*]` rather than for the function `*y*` itself.

Differential equations must be stated in terms of derivatives such as `*y*'[*x*]`, obtained with [D](https://reference.wolfram.com/language/ref/D.html), not total derivatives obtained with [Dt](https://reference.wolfram.com/language/ref/Dt.html).

Partial differential equations may also be specified using the differential operators [Grad](https://reference.wolfram.com/language/ref/Grad.html) (`∇`), [Div](https://reference.wolfram.com/language/ref/Div.html) (`∇.`), [Laplacian](https://reference.wolfram.com/language/ref/Laplacian.html) (`∇^(2)`), and [Curl](https://reference.wolfram.com/language/ref/Curl.html) (`∇⨯`). Typically these operators are used as in [Inactive](https://reference.wolfram.com/language/ref/Inactive.html)[*op*] to keep the operator form from evaluating.

[NDSolveValue](https://reference.wolfram.com/language/ref/NDSolveValue.html) solves a wide range of ordinary differential equations as well as many partial differential equations.

[NDSolveValue](https://reference.wolfram.com/language/ref/NDSolveValue.html) can also solve many delay differential equations.

In ordinary differential equations, the functions `*y*_*i*` must depend only on the single variable `*x*`. In partial differential equations, they may depend on more than one variable.

The differential equations must contain enough initial or boundary conditions to determine the solutions for the `*y*_*i*` completely.

Initial and boundary conditions are typically stated in the form `*y*[*x*_0]==*c*_0, `*y*'[*x*_0]==*dc*_0, etc., but may consist of more complicated equations.

The `*c*_0, `*dc*_0, etc. can be lists, specifying that `*y*[*x*]` is a function with vector or general list values.

Periodic boundary conditions can be specified using `*y*[*x*_0]==*y*[*x*_1]`.

The point `*x*_0 that appears in the initial or boundary conditions need not lie in the range `*x*_*min*` to `*x*_*max*` over which the solution is sought.

In delay differential equations, initial history functions are given in the form `*y*[*x*/;*x*<*x*_0]==*c*_0, where `*c*_0 is in general a function of `*x*`.

[WhenEvent](https://reference.wolfram.com/language/ref/WhenEvent.html)[*event*,*action**]* may be included in the equations `*eqns*` to specify an `*action*` that occurs when `*event*` becomes [True](https://reference.wolfram.com/language/ref/True.html).

Boundary values may also be specified using [DirichletCondition](https://reference.wolfram.com/language/ref/DirichletCondition.html) and [NeumannValue](https://reference.wolfram.com/language/ref/NeumannValue.html).

The differential equations in [NDSolveValue](https://reference.wolfram.com/language/ref/NDSolveValue.html) can involve complex numbers.

[NDSolveValue](https://reference.wolfram.com/language/ref/NDSolveValue.html) can solve many differential[Hyphen]algebraic equations, in which some of the `*eqns*` are purely algebraic, or some of the variables are implicitly algebraic.

The `*y*_*i*` can be functions of the dependent variables, and need not include all such variables.

The following options can be given:

| [AccuracyGoal](https://reference.wolfram.com/language/ref/AccuracyGoal.html) | [Automatic](https://reference.wolfram.com/language/ref/Automatic.html) | digits of absolute accuracy sought |
| --- | --- | --- |
| [Compiled](https://reference.wolfram.com/language/ref/Compiled.html) | [Automatic](https://reference.wolfram.com/language/ref/Automatic.html) | whether expressions should be compiled automatically |
| [DependentVariables](https://reference.wolfram.com/language/ref/DependentVariables.html) | [Automatic](https://reference.wolfram.com/language/ref/Automatic.html) | the list of all dependent variables |
| [EvaluationMonitor](https://reference.wolfram.com/language/ref/EvaluationMonitor.html) | [None](https://reference.wolfram.com/language/ref/None.html) | expression to evaluate whenever the function is evaluated |
| [InitialSeeding](https://reference.wolfram.com/language/ref/InitialSeeding.html) | {} | seeding equations for some algorithms |
| [InterpolationOrder](https://reference.wolfram.com/language/ref/InterpolationOrder.html) | [Automatic](https://reference.wolfram.com/language/ref/Automatic.html) | the continuity degree of the final output |
| [MaxStepFraction](https://reference.wolfram.com/language/ref/MaxStepFraction.html) | 1/10 | maximum fraction of range to cover in each step |
| [MaxSteps](https://reference.wolfram.com/language/ref/MaxSteps.html) | [Automatic](https://reference.wolfram.com/language/ref/Automatic.html) | maximum number of steps to take |
| [MaxStepSize](https://reference.wolfram.com/language/ref/MaxStepSize.html) | [Automatic](https://reference.wolfram.com/language/ref/Automatic.html) | maximum size of each step |
| [Method](https://reference.wolfram.com/language/ref/Method.html) | [Automatic](https://reference.wolfram.com/language/ref/Automatic.html) | method to use |
| [NormFunction](https://reference.wolfram.com/language/ref/NormFunction.html) | [Automatic](https://reference.wolfram.com/language/ref/Automatic.html) | the norm to use for error estimation |
| [PrecisionGoal](https://reference.wolfram.com/language/ref/PrecisionGoal.html) | [Automatic](https://reference.wolfram.com/language/ref/Automatic.html) | digits of precision sought |
| [StartingStepSize](https://reference.wolfram.com/language/ref/StartingStepSize.html) | [Automatic](https://reference.wolfram.com/language/ref/Automatic.html) | initial step size used |
| [StepMonitor](https://reference.wolfram.com/language/ref/StepMonitor.html) | [None](https://reference.wolfram.com/language/ref/None.html) | expression to evaluate when a step is taken |
| [WorkingPrecision](https://reference.wolfram.com/language/ref/WorkingPrecision.html) | [MachinePrecision](https://reference.wolfram.com/language/ref/MachinePrecision.html) | precision to use in internal computations |

[NDSolveValue](https://reference.wolfram.com/language/ref/NDSolveValue.html) adapts its step size so that the estimated error in the solution is just within the tolerances specified by [PrecisionGoal](https://reference.wolfram.com/language/ref/PrecisionGoal.html) and [AccuracyGoal](https://reference.wolfram.com/language/ref/AccuracyGoal.html).

The option [NormFunction](https://reference.wolfram.com/language/ref/NormFunction.html)->*f* specifies that the estimated errors for each of the `*y*_*i*` should be combined using `*f*[{*e*_1,*e*_2,…}]`.

[AccuracyGoal](https://reference.wolfram.com/language/ref/AccuracyGoal.html) effectively specifies the absolute local error allowed at each step in finding a solution, while [PrecisionGoal](https://reference.wolfram.com/language/ref/PrecisionGoal.html) specifies the relative local error.

If solutions must be followed accurately when their values are close to zero, [AccuracyGoal](https://reference.wolfram.com/language/ref/AccuracyGoal.html) should be set larger, or to [Infinity](https://reference.wolfram.com/language/ref/Infinity.html).

The default setting of [Automatic](https://reference.wolfram.com/language/ref/Automatic.html) for [AccuracyGoal](https://reference.wolfram.com/language/ref/AccuracyGoal.html) and [PrecisionGoal](https://reference.wolfram.com/language/ref/PrecisionGoal.html) is equivalent to [WorkingPrecision](https://reference.wolfram.com/language/ref/WorkingPrecision.html)/2.

The default setting of [Automatic](https://reference.wolfram.com/language/ref/Automatic.html) for [MaxSteps](https://reference.wolfram.com/language/ref/MaxSteps.html) estimates the maximum number of steps to be taken by [NDSolve](https://reference.wolfram.com/language/ref/NDSolve.html), depending on start and stop time and an estimate of the step size. Should this not be possible, a fixed number of steps is taken.

The setting for [MaxStepFraction](https://reference.wolfram.com/language/ref/MaxStepFraction.html) specifies the maximum step to be taken by [NDSolve](https://reference.wolfram.com/language/ref/NDSolve.html) as a fraction of the range of values for each independent variable.

With [DependentVariables](https://reference.wolfram.com/language/ref/DependentVariables.html)->[Automatic](https://reference.wolfram.com/language/ref/Automatic.html), [NDSolve](https://reference.wolfram.com/language/ref/NDSolve.html) attempts to determine the dependent variables by analyzing the equations given.

[NDSolveValue](https://reference.wolfram.com/language/ref/NDSolveValue.html) typically solves differential equations by going through several different stages depending on the type of equations. With [Method](https://reference.wolfram.com/language/ref/Method.html)->{*s*_1->*m*_1,*s*_2->*m*_2,…}, stage `*s*_*i*` is handled by method `*m*_*i*`. The actual stages used and their order are determined by [NDSolve](https://reference.wolfram.com/language/ref/NDSolve.html), based on the problem to solve.

Possible solution stages are:

"TimeIntegration" | time integration for systems of differential equations
"BoundaryValues" | ordinary differential equation boundary value solutions
"DiscontinuityProcessing" | symbolic processing for handling of discontinuous differential equations
"EquationSimplification" | simplification of equation form for numerical evaluation
"IndexReduction" | symbolic index reduction for differential algebraic equations
"DAEInitialization" | consistent initialization for differential algebraic equations
"PDEDiscretization" | discretization for partial differential equations

With [Method](https://reference.wolfram.com/language/ref/Method.html)->*m*_1 or [Method](https://reference.wolfram.com/language/ref/Method.html)->{*m*_1,*s*_2->*m*_2,…}, the method `*m*_1 is assumed to be for time integration, so [Method](https://reference.wolfram.com/language/ref/Method.html)->*m*_1 is equivalent to [Method](https://reference.wolfram.com/language/ref/Method.html)->{"TimeIntegration"->*m*_1}.

Possible explicit time integration settings for the [Method](https://reference.wolfram.com/language/ref/Method.html) option include:

"Adams" | predictor[Hyphen]corrector Adams method with orders 1 through 12
"BDF" | implicit backward differentiation formulas with orders 1 through 5
"ExplicitRungeKutta" | adaptive embedded pairs of 2(1) through 9(8) Runge-Kutta methods
"ImplicitRungeKutta" | families of arbitrary[Hyphen]order implicit Runge-Kutta methods
"SymplecticPartitionedRungeKutta" | interleaved Runge-Kutta methods for separable Hamiltonian systems

With [Method](https://reference.wolfram.com/language/ref/Method.html)->{"StyleBox[controller, "TI"]",[Method](https://reference.wolfram.com/language/ref/Method.html)->"StyleBox[submethod, "TI"]"} or [Method](https://reference.wolfram.com/language/ref/Method.html)->{"StyleBox[controller, "TI"]",[Method](https://reference.wolfram.com/language/ref/Method.html)->{*m*_1,*m*_2,…}}, possible controller methods include:

"Composition" | compose a list of submethods
"DoubleStep" | adapt step size by the double[Hyphen]step method
"EventLocator" | respond to specified events
"Extrapolation" | adapt order and step size using polynomial extrapolation
"FixedStep" | use a constant step size
"OrthogonalProjection" | project solutions to fulfill orthogonal constraints
"Projection" | project solutions to fulfill general constraints
"Splitting" | split equations and use different submethods
"StiffnessSwitching" | switch from explicit to implicit methods if stiffness is detected

Methods used mainly as submethods include:

"ExplicitEuler" | forward Euler method
"ExplicitMidpoint" | midpoint rule method
"ExplicitModifiedMidpoint" | midpoint rule method with Gragg smoothing
"LinearlyImplicitEuler" | linearly implicit Euler method
"LinearlyImplicitMidpoint" | linearly implicit midpoint rule method
"LinearlyImplicitModifiedMidpoint" | linearly implicit Bader[Hyphen]smoothed midpoint rule method
"LocallyExact" | numerical approximation to locally exact symbolic solution

The setting [InterpolationOrder](https://reference.wolfram.com/language/ref/InterpolationOrder.html)->[All](https://reference.wolfram.com/language/ref/All.html) specifies that [NDSolve](https://reference.wolfram.com/language/ref/NDSolve.html) should generate solutions that use interpolation of the same order as the underlying method used.

## Examples

### Basic Examples

Solve a first-order ordinary differential equation:

```wolfram
ysol=NDSolveValue[{y'[x]==y[x]Cos[x+y[x]],y[0]==1},y,{x,0,30}]
(* Output *)
InterpolatingFunction[...]
```

Use the solution in a plot:

```wolfram
Plot[ysol[x],{x,0,30},PlotRange->All]
```

*([Graphics])*

Use the function and its derivative in a plot:

```wolfram
ParametricPlot[{ysol[x],ysol'[x]},{x,0,20}]
```

*([Graphics])*

Find specific values:

```wolfram
{ysol[10.5],ysol'[12.5]}
(* Output *)
{0.046930201363369874,0.1137977923447328}
```

[NDSolveValue](https://reference.wolfram.com/language/ref/NDSolveValue.html) can also directly return the above values:

```wolfram
NDSolveValue[{y'[x]==y[x]Cos[x+y[x]],y[0]==1},{y[10.5],y'[12.5]},{x,0,30}]
(* Output *)
{0.046930201363369874,0.1137977923447328}
```

Second-order nonlinear ordinary differential equation:

```wolfram
ysol=NDSolveValue[{y''[x]+Sin[y[x]]y[x]==0,y[0]==1,y'[0]==0},y,{x,0,30}]
(* Output *)
InterpolatingFunction[...]
```

Plot the function and its first two derivatives:

```wolfram
Plot[{ysol[x],ysol'[x],ysol''[x]},{x,0,30}]
```

*([Graphics])*

Alternatively, substitute the functions to plot directly into [NDSolveValue](https://reference.wolfram.com/language/ref/NDSolveValue.html):

```wolfram
Plot[Evaluate[NDSolveValue[{y''[x]+Sin[y[x]]y[x]==0,y[0]==1,y'[0]==0},{y[x],y'[x],y''[x]},{x,0,30}]],{x,0,30}]
```

*([Graphics])*

System of ordinary differential equations:

```wolfram
{xsol,ysol}=NDSolveValue[{x'[t]==-y[t]-x[t]^2,y'[t]==2x[t]-y[t]^3,x[0]==y[0]==1},{x,y},{t,20}]
(* Output *)
{InterpolatingFunction[...],InterpolatingFunction[...]}
```

```wolfram
ParametricPlot[{xsol[t],ysol[t]},{t,0,20}]
```

*([Graphics])*

This solves the heat equation in one dimension:

```wolfram
usol=NDSolveValue[{D[u[t,x],t]==D[u[t,x],x,x],u[0,x]==0,u[t,0]==Sin[t],u[t,5]==0},u,{t,0,10},{x,0,5}]
(* Output *)
InterpolatingFunction[...]
```

```wolfram
Plot3D[usol[t,x],{t,0,10},{x,0,5},PlotRange->All]
```

*([Graphics3D])*

Alternative form of equation:

```wolfram
NDSolveValue[{∂_tu[t,x]==∂_x,xu[t,x],u[0,x]==0,u[t,0]==Sin[t],u[t,5]==0},u,{t,0,10},{x,0,5}]
(* Output *)
InterpolatingFunction[...]
```

Solve the Poisson equation over a [Disk](https://reference.wolfram.com/language/ref/Disk.html):

```wolfram
NDSolveValue[{-Laplacian[u[x,y],{x,y}]==1,DirichletCondition[u[x,y]==0,True]},u[x,y],{x,y}∈Disk[]]
(* Output *)
InterpolatingFunction[...][x,y]
```

```wolfram
Plot3D[%,{x,y}∈Disk[]]
```

*([Graphics3D])*

Find a minimal surface over a [Disk](https://reference.wolfram.com/language/ref/Disk.html) with a sinusoidal boundary condition:

```wolfram
NDSolveValue[{-((1)/(Sqrt[1+∇_{x,y}u[x,y].∇_{x,y}u[x,y]]))u[x,y]==0,DirichletCondition[u[x,y]==Sin[2π*(x+y)],True]},u[x,y],{x,y}∈Disk[]]
(* Output *)
InterpolatingFunction[...][x,y]
```

```wolfram
Plot3D[%,{x,y}∈Disk[]]
```

*([Graphics3D])*

Solve a coupled nonlinear sine-Gordon equation over a region:

```wolfram
Ω=Polygon[2{{-Pi,-E},{0,-1},{0,0},{Pi,0},{Pi,E},{-Pi,E}}];
NDSolveValue[{Laplacian[u[x,y],{x,y}]+λ Cos[v[x,y]]Sin[u[x,y]]==0,Laplacian[v[x,y],{x,y}]+λ Sin[u[x,y]]Sin[v[x,y]]==0,DirichletCondition[{u[x,y]==1/10,v[x,y]==1/20},True]}/.λ->5,{u[x,y],v[x,y]},{x,y}∈Ω]
(* Output *)
{InterpolatingFunction[...][x,y],InterpolatingFunction[...][x,y]}
```

```wolfram
Plot3D[%,{x,y}∈Ω]
```

*([Graphics3D])*

### Scope

#### Ordinary Differential Equations

Specify any order equation; reduction to normal form is done automatically:

```wolfram
xsol=NDSolveValue[{x''[t]+1/10x'[t]+Sin[x[t]]==1/2Cos[t],x[0]==x'[0]==0},x,{t,0,100}]
(* Output *)
InterpolatingFunction[...]
```

Directly differentiate the solution to make a phase plot:

```wolfram
ParametricPlot[{xsol[t],xsol'[t]},{t,0,100},ColorFunction->Hue]
```

*([Graphics])*

Directly specify a system of equations:

```wolfram
{xsol,ysol,zsol}=NDSolveValue[{x^′[t]==-3(x[t]-y[t]),y^′[t]==-x[t]z[t]+26.5x[t]-y[t],z^′[t]==x[t]y[t]-z[t],x[0]==z[0]==0,y[0]==1},{x,y,z},{t,0,5}]
(* Output *)
{InterpolatingFunction[...],InterpolatingFunction[...],InterpolatingFunction[...]}
```

```wolfram
Plot[{xsol[t],ysol[t],zsol[t]},{t,0,5}]
```

*([Graphics])*

Solve for a vector-valued function:

```wolfram
ysol=NDSolveValue[{y^′′[x]+({{1, 1, 1, 1}, {1, 2, 1, 2}, {1, 1, 3, 1}, {1, 2, 1, 4}}).y[x]==0,y[0]==y^′[0]=={1,1,1,1}},y,{x,0,8}]
(* Output *)
InterpolatingFunction[...]
```

```wolfram
ysol[5]
(* Output *)
{0.36310359036080914,0.6584610925358916,1.215740221077966,1.4692304088955441}
```

Plot the four components of the solution:

```wolfram
Plot[ysol[x],{x,0,8}]
```

*([Graphics])*

Different equivalent ways of specifying a harmonic oscillator as a second-order equation:

```wolfram
xsol=NDSolveValue[{x''[t]+x[t]==0,x[0]==1,x'[0]==0},x,{t,10}]
(* Output *)
InterpolatingFunction[...]
```

```wolfram
Plot[Evaluate[{xsol[t],xsol'[t]}],{t,0,10}]
```

*([Graphics])*

As a system of first-order equations:

```wolfram
{xsol,ysol}=NDSolveValue[{x'[t]==y[t],y'[t]==-x[t],x[0]==1,y[0]==0},{x,y},{t,0,10}]
(* Output *)
{InterpolatingFunction[...],InterpolatingFunction[...]}
```

```wolfram
Plot[Evaluate[{xsol[t],ysol[t]}],{t,0,10}]
```

*([Graphics])*

Using a vector variable with the dimension deduced from the initial condition:

```wolfram
zsol=NDSolveValue[{z'[t]=={{0,1},{-1,0}}.z[t],z[0]=={1,0}},z,{t,0,10}]
(* Output *)
InterpolatingFunction[...]
```

```wolfram
Plot[zsol[t],{t,0,10}]
```

*([Graphics])*

Use matrix-valued variables to compute the fundamental matrix solution:

```wolfram
A=RandomReal[{0,1},{5,5}];
```

```wolfram
fs=NDSolveValue[{x'[t]==A.x[t],x[0]==IdentityMatrix[5]},x,{t,0,1}]
(* Output *)
InterpolatingFunction[...]
```

Compare to the exact solution:

```wolfram
Plot[Norm[Flatten[fs[t]-MatrixExp[A t]]],{t,0,1}]
```

*([Graphics])*

Define a Van der Pol equation:

```wolfram
vdp={x'[t]==y[t],y'[t]==-x[t]+1000(1-x[t]^2)y[t],x[0]==2,y[0]==0};
```

The solution's stiff behavior that the default solver automatically handles:

```wolfram
{xsol,ysol}=NDSolveValue[vdp,{x,y},{t,2000}]
(* Output *)
{InterpolatingFunction[...],InterpolatingFunction[...]}
```

```wolfram
Plot[Evaluate[xsol[t]],{t,0,2000},PlotRange->All]
```

*([Graphics])*

Other methods may not be able to solve this system:

```wolfram
NDSolveValue[{x^′[t]==y[t],y^′[t]==-x[t]+1000(1-x[t]^2)y[t],x[0]==2,y[0]==0},{x,y},{t,2000},Method->"ExplicitRungeKutta"]
(* Output *)
NDSolveValue
(* Output *)
{InterpolatingFunction[...],InterpolatingFunction[...]}
```

The solution `*ysol*[*x*]` is continuous, as it integrates the piecewise function once:

```wolfram
ysol=NDSolveValue[{y'[x]+Cos[y[x]]==Floor[x],y[0]==1},y,{x,0,3}]
(* Output *)
InterpolatingFunction[...]
```

```wolfram
Plot[ysol[x],{x,0,3}]
```

*([Graphics])*

The solution `*ysol*[*x*]` is differentiable, whereas `*ysol*'[*x*]` is continuous only:

```wolfram
ysol=NDSolveValue[{y''[x]+y[x]==Floor[x],y[0]==y'[0]==1},y,{x,0,5}]
(* Output *)
InterpolatingFunction[...]
```

```wolfram
Plot[{ysol[x],ysol'[x]},{x,0,5}]
```

*([Graphics])*

#### Partial Differential Equations

Nonlinear advection-diffusion equation in one dimension:

```wolfram
usol=NDSolveValue[{D[u[t,x],t]==0.5D[u[t,x],x,x]+u[t,x]D[u[t,x],x],u[t,-Pi]==u[t,Pi]==0,u[0,x]==Sin[x]},u,{t,0,2},{x,-Pi,Pi}]
(* Output *)
InterpolatingFunction[...]
```

```wolfram
Plot3D[usol[t,x],{t,0,2},{x,-Pi,Pi},PlotRange->All]
```

*([Graphics3D])*

Define a system of PDEs of mixed parabolic-hyperbolic type:

```wolfram
pde={∂_tu[t,x]==∂_x((v[t,x]-1)∂_xu[t,x])+(16x t-2t-16(v[t,x]-1))(u[t,x]-1)+10x ℯ^-4x,∂_tv[t,x]==∂_{x,2}v[t,x]+∂_xu[t,x]+4u[t,x]-4+x^2-2t-10t ℯ^-4x};
```

```wolfram
bc={u[0,x]==1,v[0,x]==1,u[t,0]==1,v[t,0]==1,3u[t,1]+u^(0,1)[t,1]==3,5v^(0,1)[t,1]==ℯ^4(u[t,1]-1)};
```

```wolfram
{usol,vsol}=NDSolveValue[{pde,bc},{u,v},{x,0,1},{t,0,2}]
(* Output *)
{InterpolatingFunction[...],InterpolatingFunction[...]}
```

```wolfram
Plot3D[{usol[t,x],vsol[t,x]},{x,0,1},{t,0,2},PlotStyle->{Directive[Opacity[0.7],Red],Directive[Opacity[0.7],Blue]}]
(* Output *)
![image](img/image_001.png)
```

Nonlinear sine-Gordon equation in two spatial dimensions with periodic boundary conditions:

```wolfram
L=4;
usol=NDSolveValue[{D[u[t,x,y],t,t]==D[u[t,x,y],x,x]+D[u[t,x,y],y,y]+Sin[u[t,x,y]],u[t,-L,y]==u[t,L,y],u[t,x,-L]==u[t,x,L],u[0,x,y]==Exp[-(x^2+y^2)],Derivative[1,0,0][u][0,x,y]==0},u,{t,0,L/2},{x,-L,L},{y,-L,L}]
(* Output *)
InterpolatingFunction[]
```

Plot the solution at the final time:

```wolfram
Plot3D[usol[L/2,x,y],{x,-L,L},{y,-L,L}]
```

*([Graphics3D])*

Plot the time evolution of a radial cross section of the solution:

```wolfram
Plot3D[usol[t,x,0],{t,0,L/2},{x,0,L}]
```

*([Graphics3D])*

Solve a wave equation over a region with a slit:

```wolfram
Ω=RegionDifference[RegionDifference[Rectangle[{0,0},{2,1}],Rectangle[{9/10,0},{11/10,4/10}]],Rectangle[{9/10,6/10},{11/10,1}]];
sol=NDSolveValue[{D[u[t,x,y],{t,2}]-Laplacian[u[t,x,y],{x,y}]==0,DirichletCondition[u[t,x,y]==0,True],u[0,x,y]==2*Exp[-125((x-0.25)^2+(y-0.5)^2)],Derivative[1,0,0][u][0,x,y]==0},u,{t,0,2},{x,y}∈Ω]
(* Output *)
InterpolatingFunction[]
```

```wolfram
ListAnimate[Table[Rasterize[Plot3D[sol[t,x,y],{x,y}∈Ω,PlotRange->{-0.75,2},AspectRatio->Automatic]],{t,0,2,1/25}],SaveDefinitions->True]
```

Solve a Poisson equation with periodic boundary conditions on curved boundaries:

```wolfram
Ω=RegionDifference[RegionUnion[Disk[],Rectangle[{0,-1},{2,1}]],Disk[{2,0}]];
```

```wolfram
ufun=NDSolveValue[{-∇_{x,y}^2u[x,y]==1.,PeriodicBoundaryCondition[u[x,y],(x-2)^2+y^2==1,Function[x,x-{2,0}]],DirichletCondition[u[x,y]==0,(0<x<(2-10^-6)&&(y<=-1||y>=1))]},u,{x,y}∈Ω];
```

Visualize the solution:

```wolfram
ContourPlot[ufun[x,y],{x,y}∈Ω,ColorFunction->"TemperatureMap",AspectRatio->Automatic]
(* Output *)
![image](img/image_003.png)
```

#### Boundary Value Problems

A nonlinear multipoint boundary value problem:

```wolfram
{xsol,ysol}=NDSolveValue[{x''[t]==y[t]x[t],y'[t]==2-x[t],x[0]==x[4]==1,y[1]==1},{x,y},t]
(* Output *)
{InterpolatingFunction[...],InterpolatingFunction[...]}
```

```wolfram
Plot[{xsol[t],ysol[t]},{t,0,4},Epilog->{Point[{0,1}],Point[{4,1}],Red,Point[{1,1}]}]
```

*([Graphics])*

Solve a nonlinear diffusion equation $\nabla \cdot(-u^{2} \nabla u)=4$ with Dirichlet and Neumann boundary conditions starting from an initial seed of $u=1$:

```wolfram
sol=NDSolveValue[{(-u[x]^2u[x])==4+NeumannValue[2.,x==1],DirichletCondition[u[x]==1.,x==0]},u,{x}∈Line[{{0},{1}}],InitialSeeding->{u[x]==1}]
(* Output *)
InterpolatingFunction[...]
```

Visualize the result:

```wolfram
Plot[sol[x],{x,0,1}]
```

*([Graphics])*

Solve a nonlinear equation $\nabla^{2}u+4 u \nabla u=2$ with Dirichlet boundary conditions starting from an initial seed of $u=1$:

```wolfram
sol=NDSolveValue[{∇_{x}^2u[x]+4u[x]D[u[x],x]==2,DirichletCondition[u[x]==1.,True]},u,{x}∈Line[{{0},{1}}],InitialSeeding->{u[x]==1}]
(* Output *)
InterpolatingFunction[...]
```

Visualize the result:

```wolfram
Plot[sol[x],{x,0,1}]
```

*([Graphics])*

Solve a complex-valued nonlinear reaction equation $-\nabla^{2}u+i u^{2}=10$ with Dirichlet boundary conditions;

```wolfram
sol=NDSolveValue[{-∇_{x}^2u[x]+ⅈ u[x]^2==10,DirichletCondition[u[x]==-6ⅈ,x==0],DirichletCondition[u[x]==-Sqrt[2]ⅈ,x==1]},u,{x}∈Line[{{0},{1}}]]
(* Output *)
InterpolatingFunction[...]
```

Visualize the result:

```wolfram
ReImPlot[sol[x],{x,0,1}]
```

*([Graphics])*

Solve a boundary value problem with a nonlinear load term $\nabla^{2}u=sin(u)$:

```wolfram
sol=NDSolveValue[{∇_{x}^2u[x]==Sin[u[x]],DirichletCondition[u[x]==5.,x==0],DirichletCondition[u[x]==5.1,x==π/2]},u,{x}∈Line[{{0},{π/2}}]]
(* Output *)
InterpolatingFunction[...]
```

Visualize the result:

```wolfram
Plot[sol[x],{x,0,π/2}]
```

*([Graphics])*

#### Delay Differential Equations

Solve a delay differential with two constant delays and initial history function $cos(t)$:

```wolfram
xsol=NDSolveValue[{x'[t]==x[t](x[t-Pi]-x'[t-1]),x[t/;t<=0]==Cos[t]},x,{t,0,8}]
(* Output *)
InterpolatingFunction[...]
```

Discontinuities are propagated from $t=0$ at intervals equal to the delays:

```wolfram
Plot[Evaluate[{xsol[t],xsol'[t]}],{t,0,8},PlotRange->All]
```

*([Graphics])*

Investigate stability for a linear delay differential equation:

```wolfram
Manipulate[
Module[{sol,y,t},
sol=NDSolveValue[{y'[t]==λ y[t]+μ y[t-1],y[t/;t<=0]==1-t},y,
{t,0,10}];
If[pp,ParametricPlot[{sol[t],sol[t-1]},{t,1,10},PlotRange->{{-3,3},{-3,3}}],
Plot[sol[t],{t,0,10},PlotRange->{{0,10},{-3,3}}]]],{{pp,False,"Plot in phase plane"},{True,False}},{{λ,-1},-5,5},{{μ,1},-5,5},FrameLabel->ToString[y'[t]==λ y[t]+μ y[t-1],]]
```

#### Hybrid and Discontinuous Systems

A differential equation with a discontinuous right-hand side:

```wolfram
ifun=NDSolveValue[{y'[t]==-Sign[y[t]],y[0]==1},y[t],{t,0,2}];
```

```wolfram
Plot[ifun,{t,0,2},PlotRange->All]
```

*([Graphics])*

A differential equation with multiple right-hand sides:

```wolfram
ifun=NDSolveValue[{y'[t]==Piecewise[{{.5,0<=y[t]<=1},{1,1<y[t]<=2}},0],y[0]==0.5},y[t],{t,0,3}];
```

```wolfram
Plot[ifun,{t,0,3}]
```

*([Graphics])*

A differential equation whose right-hand side changes at regular time intervals:

```wolfram
ifun=NDSolveValue[{y'[t]==a[t],y[0]==0,a[0]==1,WhenEvent[Mod[t,1]==0,a[t]->-a[t]]},y[t],{t,0,5},DiscreteVariables->a[t]];
```

```wolfram
Plot[ifun,{t,0,5}]
```

*([Graphics])*

Reflect a solution across the $y$ axis each time it crosses the negative $x$ axis:

```wolfram
de={x'[t]==y[t],y'[t]==-x[t]+.2y[t]+1};
ic={x[0]==.5,y[0]==.5};
```

```wolfram
ifun=NDSolveValue[{de,ic,WhenEvent[And[y[t]==0,x[t]<0],x[t]->-x[t]]},{x[t],y[t]},{t,0,1000}];
```

```wolfram
ParametricPlot[ifun,{t,0,200},PlotPoints->200]
```

*([Graphics])*

Periodic solution with sliding mode:

```wolfram
Manipulate[Module[{x0,y0,xsol,ysol},
{x0,y0}=p;{xsol,ysol}=NDSolveValue[{x'[t]==y[t],y'[t]==If[y[t]-1>0,-1,.1y[t]-x[t]],x[0]==x0,y[0]==y0},{x,y},{t,0,20}];
ParametricPlot[{xsol[t],ysol[t]},{t,0,20},PlotRange->{{-1.5,1.5},{-1.5,1.5}}]],
{{p,{-1.5,1.5}},Locator}]
```

#### Differential-Algebraic Equations

A differential equation with an algebraic constraint:

```wolfram
{xsol,ysol}=NDSolveValue[{x'[t]==y[t]^2+x[t]y[t],2x[t]^2+y[t]^2==1,x[0]==0},{x,y},{t,0,10}]
(* Output *)
{InterpolatingFunction[...],InterpolatingFunction[...]}
```

```wolfram
ParametricPlot[{xsol[t],ysol[t]},{t,0,10},PlotRange->All,AspectRatio->1]
```

*([Graphics])*

### Generalizations & Extensions

The names of functions need not be symbols:

```wolfram
eqns={Table[y[i]^′[x]==y[i-1][x]-y[i][x],{i,2,4}],{y[1]^′[x]==-y[1][x],y[5]^′[x]==y[4][x],y[1][0]==1},Table[y[i][0]==0,{i,2,5}]}
(* Output *)
{{y[2]^′[x]==y[1][x]-y[2][x],y[3]^′[x]==y[2][x]-y[3][x],y[4]^′[x]==y[3][x]-y[4][x]},{y[1]^′[x]==-y[1][x],y[5]^′[x]==y[4][x],y[1][0]==1},{y[2][0]==0,y[3][0]==0,y[4][0]==0,y[5][0]==0}}
```

```wolfram
sol=NDSolveValue[eqns,Table[y[i],{i,5}],{x,10}]
(* Output *)
{InterpolatingFunction[...],InterpolatingFunction[...],InterpolatingFunction[...],InterpolatingFunction[...],InterpolatingFunction[...]}
```

```wolfram
Plot[Evaluate[Table[sol[[i]][x],{i,5}]],{x,0,10}]
```

*([Graphics])*

### Options

#### AccuracyGoal and PrecisionGoal

Use defaults to solve a celestial mechanics equation with sensitive dependence on initial conditions:

```wolfram
eqn=z''[t]==-z[t]/(z[t]^2+((1+Sin[2π t]/2)/2)^2)^(3/2);
```

```wolfram
Plot[Evaluate[NDSolveValue[{eqn,z[0]==1,z'[0]==0},z[t],{t,0,40}]],{t,0,40}]
```

*([Graphics])*

Higher accuracy and precision goals give a different result:

```wolfram
Plot[Evaluate[NDSolveValue[{eqn,z[0]==1,z'[0]==0},z[t],{t,0,40},AccuracyGoal->10,PrecisionGoal->10]],{t,0,40}]
```

*([Graphics])*

Increasing the goals extends the correct solution further:

```wolfram
Plot[Evaluate[NDSolveValue[{eqn,z[0]==1,z'[0]==0},z[t],{t,0,40},AccuracyGoal->20,PrecisionGoal->20,WorkingPrecision->35,MaxSteps->Infinity]],{t,0,40}]
```

*([Graphics])*

#### DependentVariables

Set up a very large system of equations:

```wolfram
n=1000;
vars=Table[x_i[t],{i,n}];
eqns=Table[j=Mod[i,n]+1;{x_i'[t]==1/(x_i[t]+x_j[t])^2,x_i[0]==1/i},{i,n}];
```

```wolfram
Short[eqns,3]
(* Output *)
{{x_1^′[t]==(1)/((x_1[t]+x_2[t])^2),x_1[0]==1},{x_2^′[t]==(1)/((<<1>>[t]+<<1>>)^2),x_2[0]==(1)/(2)},<<996>>,{x_999^′[t]==(1)/((<<1>>)^2),x_999[0]==(1)/(999)},{x_1000^′[t]==(1)/((x_1[t]+x_<<4>>[t])^2),x_1000[0]==(1)/(1000)}}
```

Solve for all the dependent variables, but save only the solution for `x_1:

```wolfram
x1sol=NDSolveValue[eqns,x_1,{t,0,100},DependentVariables->vars]
(* Output *)
InterpolatingFunction[...]
```

```wolfram
Plot[x1sol[t],{t,0,100},PlotRange->All]
```

*([Graphics])*

#### EvaluationMonitor

Total number of evaluations:

```wolfram
Module[{c=0},NDSolveValue[{y'[x]==y[x]Cos[x+y[x]],y[0]==1},y,{x,0,30},EvaluationMonitor:>c++];c]
(* Output *)
722
```

The distance between successive evaluations; negative distance means a rejected step:

```wolfram
ListLinePlot[Differences[Reap[NDSolveValue[{y'[x]==y[x]Cos[x+y[x]],y[0]==1},y,{x,0,30},EvaluationMonitor:>Sow[x]]][[2,1]]]]
```

*([Graphics])*

#### InitialSeedings

Specify an initial seeding of 0 for a boundary value problem:

```wolfram
NDSolveValue[{-2u[x]u^′[x]-u^′′[x]==1,DirichletCondition[u[x]==0,x==0],DirichletCondition[u[x]==0,x==1]},u,{x}∈Line[{{0},{1}}],InitialSeeding->{u[x]==0}]
(* Output *)
InterpolatingFunction[...]
```

Specify an initial seeding that depends on a spatial coordinate:

```wolfram
NDSolveValue[{-2u[x]u^′[x]-u^′′[x]==0,DirichletCondition[u[x]==0,x==0],DirichletCondition[u[x]==1,x==1]},u,{x}∈Line[{{0},{1}}],InitialSeeding->{u[x]==x}]
(* Output *)
InterpolatingFunction[...]
```

#### InterpolationOrder

Use [InterpolationOrder](https://reference.wolfram.com/language/ref/InterpolationOrder.html)->[All](https://reference.wolfram.com/language/ref/All.html) to get interpolation the same order as the method:

```wolfram
Timing[{xa,ya}=NDSolveValue[{x'[t]==y[t],y'[t]==-x[t],x[0]==1,y[0]==0},{x,y},{t,0,1000},InterpolationOrder->All]]
(* Output *)
{0.019703,{InterpolatingFunction[],InterpolatingFunction[]}}
```

This is more time consuming than the default interpolation order used:

```wolfram
Timing[{xd,yd}=NDSolveValue[{x'[t]==y[t],y'[t]==-x[t],x[0]==1,y[0]==0},{x,y},{t,0,1000}]]
(* Output *)
{0.007657,{InterpolatingFunction[...],InterpolatingFunction[...]}}
```

It is much better in between steps:

```wolfram
Plot[Evaluate[{xa[t],xd[t]}-Cos[t]],{t,0,10},PlotRange->All,PlotStyle->{Directive[Thick,Blue],Red}]
```

*([Graphics])*

#### MaxStepFraction

Features with small relative size in the integration interval can be missed:

```wolfram
Table[Plot[Evaluate[NDSolveValue[{y''[x]+y[x]==Exp[-(1000(x/L-1/3))^2],y[0]==y'[0]==0},y[x],{x,L}]],{x,L/4,L/2},PlotRange->All],{L,{1,10,100,1000}}]
(* Output *)
{[Graphics],[Graphics],[Graphics],[Graphics]}
```

Use [MaxStepFraction](https://reference.wolfram.com/language/ref/MaxStepFraction.html) to ensure features are not missed, independent of interval size:

```wolfram
Table[Plot[Evaluate[NDSolveValue[{y''[x]+y[x]==Exp[-(1000(x/L-1/3))^2],y[0]==y'[0]==0},y[x],{x,L},MaxStepFraction->0.001]],{x,L/4,L/2},PlotRange->All],{L,{1,10,100,1000}}]
(* Output *)
{[Graphics],[Graphics],[Graphics],[Graphics]}
```

#### MaxSteps

Integration stops short of the requested interval:

```wolfram
NDSolveValue[{y''[x]+x y[x]==0,y[0]==1,y'[0]==0},y,{x,0,200}]
(* Output *)
InterpolatingFunction[...]
```

More steps are needed to resolve the solution:

```wolfram
ysol=NDSolveValue[{y''[x]+x y[x]==0,y[0]==1,y'[0]==0},y,{x,0,200},MaxSteps->10^6]
(* Output *)
InterpolatingFunction[...]
```

Plot the solution in the phase plane:

```wolfram
ParametricPlot[Evaluate[{ysol[x],ysol'[x]}],{x,0,200},ColorFunction->Function[Hue[#3]],AspectRatio->1]
```

*([Graphics])*

#### MaxStepSize

The default step control may miss a suddenly varying feature:

```wolfram
NDSolveValue[{y''[x]+(1+Sech[1000(x-Pi)])y[x]==0,y[0]==1,y'[0]==0},y,{x,0,10}][x]
(* Output *)
InterpolatingFunction[...][x]
```

```wolfram
Plot[Evaluate[NDSolveValue[{y''[x]+(1+Sech[1000(x-Pi)])y[x]==0,y[0]==1,y'[0]==0},y,{x,0,10}][x]-Cos[x]],{x,0,10}]
```

*([Graphics])*

A smaller [MaxStepSize](https://reference.wolfram.com/language/ref/MaxStepSize.html) setting ensures that [NDSolveValue](https://reference.wolfram.com/language/ref/NDSolveValue.html) catches the feature:

```wolfram
Plot[Evaluate[NDSolveValue[{y''[x]+(1+Sech[1000(x-Pi)])y[x]==0,y[0]==1,y'[0]==0},y,{x,0,10}, MaxStepSize->0.01][x]-Cos[x]],{x,0,10}]
```

*([Graphics])*

Attempting to compute the number of positive integers less than $e^{5}$ misses several events:

```wolfram
Block[{n=0},NDSolveValue[{y'[t]==y[t],y[-1]==ℯ^-1,WhenEvent[Sin[π y[t]]==0,n++]},y,{t,5}];n]
(* Output *)
18
```

Setting a small enough [MaxStepSize](https://reference.wolfram.com/language/ref/MaxStepSize.html) ensures that none of the events are missed:

```wolfram
Block[{n=0},NDSolveValue[{y'[t] ==y[t],y[-1]==ℯ^-1,WhenEvent[Sin[π y[t]]==0,n++]},y,{t,5},MaxStepSize->0.001];n]
(* Output *)
148
```

#### Method

TimeIntegration  (5)

Differences between values of $x$ at successive steps with the default solution method:

```wolfram
ListLinePlot[Differences[Reap[NDSolveValue[{y'[x]==y[x]Cos[x+y[x]],y[0]==1},y,{x,0,30},StepMonitor:>Sow[x]]][[2,1]]]]
```

*([Graphics])*

With an explicit Runge-Kutta method, the step size is changed more often:

```wolfram
ListLinePlot[Differences[Reap[NDSolveValue[{y'[x]==y[x]Cos[x+y[x]],y[0]==1},y,{x,0,30},StepMonitor:>Sow[x],Method->"ExplicitRungeKutta"]][[2,1]]]]
```

*([Graphics])*

Difference order of 8:

```wolfram
ListLinePlot[Differences[Reap[NDSolveValue[{y'[x]==y[x]Cos[x+y[x]],y[0]==1},y,{x,0,30},StepMonitor:>Sow[x],Method->{"ExplicitRungeKutta","DifferenceOrder"->8}]][[2,1]]]]
```

*([Graphics])*

With a difference order of 3, the steps are much smaller:

```wolfram
ListLinePlot[Differences[Reap[NDSolveValue[{y'[x]==y[x]Cos[x+y[x]],y[0]==1},y,{x,0,30},StepMonitor:>Sow[x],Method->{"ExplicitRungeKutta","DifferenceOrder"->3}]][[2,1]]]]
```

*([Graphics])*

Extrapolation tends to take very large steps:

```wolfram
ListLinePlot[Differences[Reap[NDSolveValue[{y'[x]==y[x]Cos[x+y[x]],y[0]==1},y,{x,0,30},StepMonitor:>Sow[x],Method->"Extrapolation"]][[2,1]]]]
```

*([Graphics])*

BoundaryValues  (1)

Solve a boundary value problem:

```wolfram
xsol0=NDSolveValue[{x''[t]+Sin[x[t]]==0,x[0]==x[10]==0},x,t]
(* Output *)
InterpolatingFunction[...]
```

With the default option, the method finds the trivial solution:

```wolfram
Table[xsol0[t],{t,0,10}]
(* Output *)
{0.,0.7832451461046588,0.9161379502703044,0.302630716501597,-0.5726632750658163,-0.9614507521954894,-0.5726634011878152,0.30263050291269733,0.9161378790954307,0.7832451752551658,1.2447046624073997×10^-7}
```

Specify different starting conditions for the `"Shooting"` method to find different solutions:

```wolfram
xsols=Map[NDSolveValue[{x''[t]+Sin[x[t]]==0,x[0]==x[10]==0},x,t,Method->{"Shooting","StartingInitialConditions"->{x[0]==0,x'[0]==#}}]&,
{1.5,1.75,2}];
Plot[Evaluate[Through[xsols[t]]],{t,0,10},PlotStyle->{StandardGray,Blue,Green}]
```

*([Graphics])*

DiscontinuityProcessing  (1)

[NDSolveValue](https://reference.wolfram.com/language/ref/NDSolveValue.html) automatically does processing for discontinuous functions like [Sign](https://reference.wolfram.com/language/ref/Sign.html):

```wolfram
xsol1=NDSolveValue[{x'[t]==Sign[1-x[t]],x[0]==0},x,{t,0,2}]
(* Output *)
InterpolatingFunction[...]
```

If the processing is turned off, [NDSolveValue](https://reference.wolfram.com/language/ref/NDSolveValue.html) may fail at the discontinuity point:

```wolfram
NDSolveValue[{x'[t]==Sign[1-x[t]],x[0]==0},x,{t,0,2},Method->{"DiscontinuityProcessing"->False}]
(* Output *)
NDSolveValue
(* Output *)
InterpolatingFunction[]
```

With some time integration methods, the solution may be very inaccurate:

```wolfram
xsol2=NDSolveValue[{x'[t]==Sign[1-x[t]],x[0]==0},x,{t,0,2},Method->{"DiscontinuityProcessing"->False,"TimeIntegration"->"Extrapolation"}]
(* Output *)
InterpolatingFunction[...]
```

```wolfram
Plot[{xsol1[t],xsol2[t]},{t,0,2}]
```

*([Graphics])*

An equivalent way to find the solution is to use `"DiscontinuitySignature"`:

```wolfram
xsol3=NDSolveValue[{x^′[t]==v[t],WhenEvent[1==x[t],v[t]->"DiscontinuitySignature"],x[0]==0,v[0]==1},x,{t,0,2},DiscreteVariables->{Element[v,{-1,0,1}]}]
(* Output *)
InterpolatingFunction[...]
```

The solutions are effectively identical:

```wolfram
Table[(xsol1[t])-(xsol3[t]),{t,0,2,.25}]
(* Output *)
{0.,0.,0.,0.,0.,0.,0.,0.,0.}
```

EquationSimplification  (1)

The solution cannot be completed because the square root function is not sufficiently smooth:

```wolfram
NDSolveValue[{x^′[t]^2+x[t]^2==1,x[0]==1/2},x,{t,0,10Pi}]
(* Output *)
NDSolveValue
(* Output *)
NDSolveValue
(* Output *)
InterpolatingFunction[]
```

One solution can be found by forming a residual and solving as a DAE system:

```wolfram
xsol1=NDSolveValue[{x^′[t]^2+x[t]^2==1,x[0]==1/2},x,{t,0,10Pi},Method->{"EquationSimplification"->"Residual"}]
(* Output *)
InterpolatingFunction[...]
```

The other solution branch can be given by specifying a consistent value of $x^{'}$:

```wolfram
xsol2=NDSolveValue[{x^′[t]^2+x[t]^2==1,x[0]==1/2,x'[0]==-Sqrt[3/4]},x,{t,0,10Pi},Method->{"EquationSimplification"->"Residual"}]
(* Output *)
InterpolatingFunction[...]
```

```wolfram
Plot[{xsol1[t],xsol2[t]},{t,0,10Pi}]
```

*([Graphics])*

IndexReduction  (1)

An index 3 formulation of a constrained pendulum using index reduction:

```wolfram
eqns={x''[t]==T[t]x[t],y''[t]==T[t]y[t]-1,x[t]^2+y[t]^2==1};
init={x[0]==1,y[0]==0,y'[0]==-14/10};
```

The default method can only solve index 1 problems:

```wolfram
NDSolveValue[{eqns,init},{T,x,y},{t,0,10}]
(* Output *)
NDSolveValue
(* Output *)
{InterpolatingFunction[...],InterpolatingFunction[...],InterpolatingFunction[...]}
```

The problem resulting from symbolic index reduction can be solved:

```wolfram
s1=NDSolveValue[{x''[t]==T[t]x[t],y''[t]==T[t]y[t]-1,x[t]^2+y[t]^2==1,x[0]==1,y[0]==0,y'[0]==-14/10},{T,x,y},{t,0,10},Method->{"IndexReduction"->True}]
(* Output *)
{InterpolatingFunction[...],InterpolatingFunction[...],InterpolatingFunction[...]}
```

Solve using reduction to index 0 and a projection method to maintain the constraints:

```wolfram
s2=NDSolveValue[{x''[t]==T[t]x[t],y''[t]==T[t]y[t]-1,x[t]^2+y[t]^2==1,x[0]==1,y[0]==0,y'[0]==-14/10},{T,x,y},{t,0,10},Method->{"IndexReduction"->{True,"ConstraintMethod"->"Projection"}}]
(* Output *)
{InterpolatingFunction[...],InterpolatingFunction[...],InterpolatingFunction[...]}
```

Plot implicit energy constraint for the two solutions at the time steps:

```wolfram
energy[sol_]:=Module[{T,x,y,times,es},
{T,x,y}=sol;times=Part[x,3,1];es[t_]=1+y[t]+(1)/(2)(-y[t]x^′[t]+x[t]y^′[t])^2;Transpose[{times,es[times]-es[0]}]];
```

```wolfram
ListPlot[{energy[s1],energy[s2]}]
```

*([Graphics])*

DAEInitialization  (1)

Use forward collocation for initialization to avoid problems with the [Abs](https://reference.wolfram.com/language/ref/Abs.html) term at 0:

```wolfram
NDSolveValue[{x'[t]==Abs[u[t]],u[t]==Sin[t],y[t]==x[t],y[0]==10},{x,y,u},{t,0,10},Method->{"DAEInitialization"->{"Collocation","CollocationDirection"->"Forward"}}]
(* Output *)
{InterpolatingFunction[...],InterpolatingFunction[...],InterpolatingFunction[...]}
```

PDEDiscretization  (2)

Solutions of Burgers' equation may steepen, leading to numerical instability:

```wolfram
usol1=NDSolveValue[{∂_tu[t,x]==(1)/(1000)∂_x,xu[t,x]-u[t,x]∂_xu[t,x],u[0,x]==Sin[2π x],u[t,0]==u[t,1]},u,{t,0,2},{x,0,1}]
(* Output *)
NDSolveValue
(* Output *)
NDSolveValue
(* Output *)
InterpolatingFunction[...]
```

```wolfram
end=Part[usol1,1,1,-1];
Plot[usol1[end,x],{x,0,1}]
```

*([Graphics])*

Specify a spatial discretization sufficiently fine to resolve the front:

```wolfram
usol2=NDSolveValue[{∂_tu[t,x]==(1)/(1000)∂_x,xu[t,x]-u[t,x]∂_xu[t,x],u[0,x]==Sin[2π x],u[t,0]==u[t,1]},u,{t,0,2},{x,0,1},Method->{"PDEDiscretization"->{"MethodOfLines","SpatialDiscretization"->{"TensorProductGrid","MinPoints"->1000}}}]
(* Output *)
InterpolatingFunction[]
```

```wolfram
end=Part[usol2,1,1,-1];
Plot[usol2[end,x],{x,0,1}]
```

*([Graphics])*

After the front forms, the solution decays relatively rapidly:

```wolfram
X=Part[usol2,3,2];
Plot[Max[usol2[t,X]],{t,0,2}]
```

*([Graphics])*

Solve a time-dependent Schrödinger equation on a refined mesh:

```wolfram
solution=NDSolveValue[{({{-(1)/(2)}}.u[t,x])==ⅈ u^(1,0)[t,x],u[0,x]==ℯ^((-3ⅈ  x-(1)/(20) x^2)),DirichletCondition[u[t,x]==0,True]},u,{t,0,5},{x,-10,10},Method->{"PDEDiscretization"->{"MethodOfLines","SpatialDiscretization"->{"FiniteElement","MeshOptions"->{"MaxCellMeasure"->0.1}}}}]
(* Output *)
InterpolatingFunction[...]
```

Plot the higher-accuracy solution:

```wolfram
Plot3D[Evaluate[Abs[solution[t,x]]],{x,-10,10},{t,0,5},PlotPoints->50,AxesLabel->{x,t}]
(* Output *)
![image](img/image_005.png)
```

#### NormFunction

Plot the actual solution error when using different error estimation norms:

```wolfram
L=5;Table[usol[p]=NDSolveValue[{D[u[t,x],t,t]==D[u[t,x],x,x],u[0,x]==Exp[-x^2],Derivative[1,0][u][0,x]==-Exp[-x^2],u[t,-L]==u[t,L]},u,{t,0,4L},{x,-L,L},NormFunction->(Norm[#,p]&)];
Plot[Evaluate[(usol[p][4L,x])-(Exp[-x^2]-2Sqrt[π])],{x,-L,L},PlotRange->All],{p,{1,2,3,∞}}]
(* Output *)
![image](img/image_007.png)
```

A plot of the best solution:

```wolfram
Plot3D[usol[1][t,x],{t,0,4L},{x,-10,10},Mesh->False,PlotPoints->20]
```

*([Graphics3D])*

#### StartingStepSize

For a very large interval, a short-lived feature near the start may be missed:

```wolfram
feature[t_?NumberQ]:=If[t<1,t Sin[10Pi t],0];xsol=NDSolveValue[{x''[t]+10^-4x[t]==feature[t],x[0]==x'[0]==0},x,{t,0,10000}]
(* Output *)
InterpolatingFunction[...]
```

```wolfram
Max[Abs[xsol[RandomReal[{0,10000},10000]]]]
(* Output *)
0.
```

Setting a sufficiently small step size to start with ensures that the input is not missed:

```wolfram
xsol=NDSolveValue[{x''[t]+10^-4x[t]==feature[t],x[0]==x'[0]==0},x,{t,0,10000},StartingStepSize->0.1]
(* Output *)
InterpolatingFunction[...]
```

```wolfram
{Plot[xsol[t],{t,0,2}],Plot[xsol[t],{t,0,10000}]}
(* Output *)
{[Graphics],[Graphics]}
```

#### StepMonitor

Plot the solution at each point where a step is taken in the solution process:

```wolfram
ListLinePlot[Reap[NDSolveValue[{y'[x]==y[x]Cos[x+y[x]],y[0]==1},y,{x,0,30},StepMonitor:>Sow[{x,y[x]}]]][[2,1]]]
```

*([Graphics])*

Total number of steps involved in finding the solution:

```wolfram
Module[{c=0},NDSolveValue[{y'[x]==y[x]Cos[x+y[x]],y[0]==1},y,{x,0,30},StepMonitor:>c++];c]
(* Output *)
337
```

Differences between values of $x$ at successive steps:

```wolfram
ListLinePlot[Differences[Reap[NDSolveValue[{y'[x]==y[x]Cos[x+y[x]],y[0]==1},y,{x,0,30},StepMonitor:>Sow[x]]][[2,1]]]]
```

*([Graphics])*

#### WorkingPrecision

Error in the solution to a harmonic oscillator over 100 periods:

```wolfram
T=200Pi;
Timing[NDSolveValue[{y''[x]+y[x]==0,y[0]==0,y'[0]==1},y,{x,0,T},MaxSteps->Infinity][T]]
(* Output *)
{0.015625,-3.25315926311991×10^-6}
```

When the working precision is increased, the local tolerances are correspondingly increased:

```wolfram
Timing[NDSolveValue[{y''[x]+y[x]==0,y[0]==0, y'[0]==1},y,{x,0,T},WorkingPrecision->32,MaxSteps->10^6][T]]
(* Output *)
{1.953125,-9.2263155697439299826536658177097896465727132×10^-14}
```

With a large working precision, sometimes the `"Extrapolation"` method is quite effective:

```wolfram
Timing[NDSolveValue[{y''[x]+y[x]==0,y[0]==0,y'[0]==1},y,{x,0,T},WorkingPrecision->32,Method->"Extrapolation"][T]]
(* Output *)
{0.484375,-4.6019234683096408493642766888056357614777×10^-16}
```

### Applications

#### Duffing Equation

Simulate Duffing's equation for a particle in a double potential well:

```wolfram
xsol=NDSolveValue[{x''[t]+0.15x'[t]-x[t]+x[t]^3==0.3Cos[t],x[0]==-1,x'[0]==1},x,{t,0,50}]
(* Output *)
InterpolatingFunction[...]
```

```wolfram
Plot[xsol[t],{t,0,50}]
```

*([Graphics])*

The solution depends strongly on initial conditions:

```wolfram
xsol=NDSolveValue[{x''[t]+0.15x'[t]-x[t]+x[t]^3==0.3Cos[t],
x[0]==-1,x'[0]==1.001},x,{t,0,50}]
(* Output *)
InterpolatingFunction[...]
```

```wolfram
Plot[xsol[t],{t,0,50}]
```

*([Graphics])*

#### Lotka-Volterra Equations

The Lotka-Volterra predator-prey equations [[more info](http://mathworld.wolfram.com/Lotka-VolterraEquations.html)]:

```wolfram
{xsol,ysol}=NDSolveValue[{y'[t]==y[t](x[t]-1),x'[t]==x[t](2-y[t]),x[0]==1,y[0]==2.7},{x,y},{t,0,10}]
(* Output *)
{InterpolatingFunction[...],InterpolatingFunction[...]}
```

Phase plane plot:

```wolfram
ParametricPlot[{xsol[t],ysol[t]},{t,0,10}]
```

*([Graphics])*

#### Lorenz Equations

The Lorenz equations [[more info](http://mathworld.wolfram.com/LorenzEquations.html)]:

```wolfram
{xsol,ysol,zsol}=NDSolveValue[{x^′[t]==-3(x[t]-y[t]),y^′[t]==-x[t]z[t]+26.5x[t]-y[t],z^′[t]==x[t]y[t]-z[t],x[0]==z[0]==0,y[0]==1},{x,y,z},{t,0,200},MaxSteps->∞]
(* Output *)
{InterpolatingFunction[...],InterpolatingFunction[...],InterpolatingFunction[...]}
```

```wolfram
ParametricPlot3D[{xsol[t],ysol[t],zsol[t]},{t,0,200},PlotPoints->1000,ColorFunction->(Hue[#4]&)]
(* Output *)
![image](img/image_009.png)
```

#### Gavrilov-Shilnikov Model

Look at the appearance of the blue sky catastrophe orbit:

```wolfram
eqns={
x'[t]==x[t](2+μ-10(x[t]^2+y[t]^2))+z[t]^2+y[t]^2+2y[t],
y'[t]==-z[t]^3-(1+y[t])(z[t]^2+y[t]^2+2y[t])-4x[t]+μ y[t],
z'[t]==(1+y[t])z[t]^2+x[t]^2-ε};
```

```wolfram
Table[Block[{ε=0.0357,T=1000},
{xsol,ysol,zsol}=NDSolveValue[{eqns,x[0]==y[0]==z[0]==1},{x,y,z},{t,0,T},MaxSteps->∞];
ParametricPlot3D[{xsol[t],ysol[t],zsol[t]},{t,T/4,T},BoxRatios->1,PlotPoints->1000,PlotRange->All]],{μ,.455,.456,.001}]
(* Output *)
![image](img/image_011.png)
```

#### Heat Equation

Simple model for soil temperature at depth $x$ with periodic heating at the surface:

```wolfram
usol=NDSolveValue[{D[u[t,x],t]==D[u[t,x],x,x],u[0,x]==0,u[t,0]==Sin[t],u[t,5]==0},u,{t,0,10},{x,0,5}]
(* Output *)
InterpolatingFunction[...]
```

```wolfram
Plot3D[Evaluate[usol[t,x]],{t,0,10},{x,0,5},PlotRange->All]
```

*([Graphics3D])*

Model a temperature field with a heat source in a rod:

```wolfram
vars={Θ[x],{x}};
pars=<|"ThermalConductivity"->0.026,"HeatSource"->1|>;
Ω=Line[{{0},{1}}];
eqn=HeatTransferPDEComponent[vars,pars]==0
(* Output *)
-1.+({{-0.026}}.Θ[x])==0
```

Solve the PDE:

```wolfram
Tfun=NDSolveValue[{eqn,HeatTemperatureCondition[x==0,vars,pars,<|"SurfaceTemperature"->0|>]},Θ,x∈Ω];
```

Visualize the solution:

```wolfram
Plot[Tfun[x],{x}∈Ω]
```

*([Graphics])*

More information about the heat equation and applicable boundary conditions can be found in the [Heat Transfer tutorial](https://reference.wolfram.com/language/PDEModels/tutorial/HeatTransfer/HeatTransfer.html) and the [Heat Transfer PDE models guide page](https://reference.wolfram.com/language/guide/HeatTransferPDEModels.html).

#### Wave Equation

Simple wave evolution with periodic boundary conditions:

```wolfram
usol=NDSolveValue[{∂_t,tu[t,x]==∂_x,xu[t,x],u[0,x]==ℯ^(-x^2),u^(1,0)[0,x]==0,u[t,-10]==u[t,10]},u,{t,0,40},{x,-10,10}]
(* Output *)
InterpolatingFunction[]
```

Plot the solution:

```wolfram
DensityPlot[usol[t,x],{t,0,40},{x,-10,10},PlotPoints->30]
```

*([Graphics])*

#### Wolfram's Nonlinear Wave Equation

Wolfram's nonlinear wave equation [[more info](http://www.wolframscience.com/nksonline/page-923a-text)]:

```wolfram
usol=NDSolveValue[{∂_t,tu[t,x]==∂_x,xu[t,x]+(1-u[t,x]^2)(1+2u[t,x]),u[0,x]==ℯ^(-x^2),u^(1,0)[0,x]==0,u[t,-10]==u[t,10]},u,{t,0,10},{x,-10,10}]
(* Output *)
InterpolatingFunction[...]
```

```wolfram
Plot3D[usol[t,x],{t,0,10},{x,-10,10}]
```

*([Graphics3D])*

```wolfram
DensityPlot[usol[-t,x],{x,-10,10},{t,0,-10}]
```

*([Graphics])*

Model with Wolfram's nonlinear wave equation [[more info](http://www.wolframscience.com/nksonline/page-923a-text)] with a Dirichlet condition.

Define initial and boundary conditions:

```wolfram
ics={u[0,x]==Exp[-x^2],Derivative[1,0][u][0,x]==0};
```

```wolfram
bcs=DirichletCondition[u[t,x]==0,True];
```

Solve the equation:

```wolfram
sol=NDSolveValue[{(1+2 u[t,x]) (-1+u[t,x]^2)+(-u[t,x])+u^(2,0)[t,x]==0,ics,bcs},u,{t,0,10},{x,-10,10}]
(* Output *)
InterpolatingFunction[...]
```

Visualize the result:

```wolfram
Plot3D[sol[t,x],{x,-10,10},{t,0,10},ColorFunction -> SouthwestColors, Axes -> None, Boxed -> False, ImageSize -> Medium]
```

*([Graphics3D])*

This solves a `2+1-dimensional version of the equation:

```wolfram
usol=NDSolveValue[{D[u[t,x,y],t,t]==D[u[t,x,y],x,x]+D[u[t,x,y],y,y]/2+(1-u[t,x,y]^2)(1+2u[t,x,y]),u[0,x,y]==ℯ^(-(x^2+y^2)),u[t,-5,y]==u[t,5,y],u[t,x,-5]==u[t,x,5],u^(1,0,0)[0,x,y]==0},u,{t,0,4},{x,-5,5},{y,-5,5}]
(* Output *)
InterpolatingFunction[]
```

```wolfram
Table[Plot3D[usol[t,x,y],{x,-5,5},{y,-5,5},PlotRange->All],{t,1,4}]
(* Output *)
{[Graphics3D],[Graphics3D],[Graphics3D],[Graphics3D]}
```

#### *n*-Body Problems

Reduced 3-body problem [[more info]](http://www.wolframscience.com/nksonline/page-973a-text):

```wolfram
eqn=z''[t]==-(z[t])/((z[t]^2+((1)/(2)(1+0.1Sin[2Pi t]))^2)^(3/2));
```

```wolfram
zsol=NDSolveValue[{eqn,z[0]==1,z'[0]==0},z,{t,0,30}]
(* Output *)
InterpolatingFunction[...]
```

```wolfram
Plot[zsol[t],{t,0,30}]
```

*([Graphics])*

A formulation suitable for a number of different initial conditions:

```wolfram
eqn=D[z[t,s],t,t]==-(z[t,s])/((z[t,s]^2+((1)/(2)(1+0.1Sin[2Pi t]))^2)^(3/2));
```

```wolfram
zsol=NDSolveValue[{eqn,z[0,s]==s,Derivative[1,0][z][0,s]==0},z,{t,0,10},{s,1,2}]
(* Output *)
InterpolatingFunction[...]
```

```wolfram
Plot3D[zsol[tp,z0],{tp,0,10},{z0,1,2}]
```

*([Graphics3D])*

#### Nonlinear Schrödinger Equation

A moving soliton driven by a periodic source term:

```wolfram
L=50;
usol=NDSolveValue[{I D[u[t,x],t]+D[u[t,x],x,x]+2Abs[u[t,x]]^2u[t,x]==0.1(1-Cos[Pi x/L]),u[0,x]==Sech[x]Exp[I x],u[t,-L]==u[t,L]},u,{t,0,10},{x,-L,L}]
(* Output *)
InterpolatingFunction[]
```

```wolfram
Plot3D[Abs[usol[t,x]],{t,0,10},{x,-L,L},PlotPoints->60,MaxRecursion->3,Mesh->None]
(* Output *)
![image](img/image_013.png)
```

A moving soliton under periodic potential:

```wolfram
L=20;
usol=NDSolveValue[{I D[u[t,x],t]+D[u[t,x],x,x]+2 Abs[u[t,x]]^2 u[t,x]+0.1 (1-Cos[Pi x/L]) u[t,x]==0,u[0,x]==Sech[x]Exp[I x],u[t,-L]==u[t,L]},u,{t,0,15},{x,-L,L}]
(* Output *)
InterpolatingFunction[...]
```

```wolfram
Plot3D[Abs[usol[t,x]],{t,0,15},{x,-L,L},PlotPoints->60,PlotRange->All,MaxRecursion->3]
(* Output *)
![image](img/image_015.png)
```

#### Korteweg-de Vries (KdV) Equation

Define a KdV equation in canonical form:

```wolfram
KdV = {D[u[x,t],t]+D[u[x,t],{x, 3}]+ 6u[x,t]D[u[x,t],x]==0}
(* Output *)
{u^(0,1)[x,t]+6 u[x,t] u^(1,0)[x,t]+u^(3,0)[x,t]==0}
```

Find a traveling soliton solution with periodic boundary conditions:

```wolfram
sol = NDSolveValue[{KdV,u[x,0]==2-2 Tanh[(1)/(4)+x]^2,u[-50,t]== u[50,t]},u,{x,-50,50},{t,0,10}]
(* Output *)
InterpolatingFunction[]
```

```wolfram
Plot[Evaluate[Table[sol[x,t],{t,1,12,3}]], {x,-5,50},PlotRange->All,Filling->Axis,Epilog->{Arrow[{{5,1.6},{9,1.6}}],Arrow[{{17,1.6},{21,1.6}}],Arrow[{{29,1.6},{33,1.6}}],Arrow[{{41,1.6},{45,1.6}}]}]
```

*([Graphics])*

Define an alternative form of the KdV equation:

```wolfram
KdV = {D[u[x,t],t]+D[u[x,t],{x, 3}]- 6u[x,t]D[u[x,t],x]==0}
(* Output *)
{u^(0,1)[x,t]-6 u[x,t] u^(1,0)[x,t]+u^(3,0)[x,t]==0}
```

Fission into two solitons from a sech-squared initial condition:

```wolfram
sol=NDSolveValue[{KdV,u[x,0]==-6 Sech[x]^2,u[-50,t]== u[50,t]},u,{x,-50,50},{t,0,10},Method->{"MethodOfLines","SpatialDiscretization"->{"TensorProductGrid","MaxPoints"->5000,"MinPoints"->5000},Method->"StiffnessSwitching"}]
(* Output *)
InterpolatingFunction[]
```

```wolfram
ListAnimate[Plot[#,{x,-50,50},PlotRange->All]&/@Table[sol[x,t],{t,0,10,0.1}],SaveDefinitions->True]
```

#### Mackey-Glass Equation

View solutions of the Mackey-Glass delay differential equation:

```wolfram
Manipulate[
Module[{xsol,x,t},
xsol=NDSolveValue[{x'[t]==a x[t-τ]/(1+x[t-τ]^10)-b x[t],x[t/;t<=0]==1/2},x,{t,0,500}];
ParametricPlot[{xsol[t],xsol[t-τ]},{t,300,500},PlotRange->{{0,2},{0,2}}]],{{a,.2},0,1},{{b,.1},0,1},
{{τ,17},1,20}]
```

#### Acoustics Equation

Define model variables `*vars*` for a transient acoustic pressure field with model parameters `*pars*`:

```wolfram
vars={p[t,x],t,{x}};
pars=<|"SoundSpeed"->343,"MassDensity"->1.2|>;
```

Define initial conditions `*ics*` of a right-going sound wave $p_{0}$:

```wolfram
p0=D[0.125 Erf[(x-0.5)/0.15],x];
ics={p[0,x]==p0,Derivative[1,0][p][0,x]==-343*D[p0,x]};
```

Set up the equation with a sound hard boundary at the right end:

```wolfram
eqn=AcousticPDEComponent[vars,pars]==AcousticSoundHardValue[x==1,vars,pars]
(* Output *)
({{-0.8333333333333334}}.p[t,x])+7.083216460261739×10^-6 p^(2,0)[t,x]==NeumannValue[0,x==1]
```

Solve the PDE:

```wolfram
pfun=NDSolveValue[{eqn,ics},p,{t,0,0.003},x∈Line[{{0},{1}}]];
```

Visualize the sound field in the time domain:

```wolfram
Manipulate[Plot[pfun[t,x],{x,0,1},PlotRange -> 02, AxesLabel -> xP],{{t,0},0,0.003,10^-4},SaveDefinitions -> True]
```

More information about the wave equation and acoustics boundary conditions can be found in the [Acoustics in the Time Domain tutorial](https://reference.wolfram.com/language/PDEModels/tutorial/Acoustics/AcousticsTimeDomain.html#481530965) and the [Acoustics PDE models guide page](https://reference.wolfram.com/language/guide/AcousticPDEModels.html).

#### Mass Transport

Model a 1D chemical species transport through different material with a reaction rate in one. The right side and left side are subjected to a mass concentration and inflow condition, respectively:

$\overset{\overset{           mass transport model              }{\text{OverBrace}}}{\nabla \cdot(-d \nabla c(x))+a c(x)} =\overset{\overset{ mass flux value  }{\text{OverBrace}}}{|_{\Gamma_{x=0}}q(x)}$

Set up the stationary mass transport model variables `*vars*`:

```wolfram
vars={c[x],{x}};
```

Set up a region $\Omega$:

```wolfram
Ω=Line[{{0},{1}}];
```

Specify the mass transport model parameters species diffusivity $d$ and a reaction rate $a$ active in the region $x>1/2$:

```wolfram
pars=<|"DiffusionCoefficient"->0.01,"MassReactionRate"->If[x>1/2,1,0]|>;
```

Specify a species flux boundary condition:

```wolfram
Γ_flux=MassFluxValue[x==0,vars,pars,<|"MassFlux"->1|>]
(* Output *)
NeumannValue[1+c[x] NeumannBoundaryUnitNormal[{x}].{0},x==0]
```

Specify a mass concentration boundary condition:

```wolfram
Γ_C=MassConcentrationCondition[x==1,vars,pars,<|"MassConcentration"->0|>]
(* Output *)
DirichletCondition[c[x]==0,x==1]
```

Set up the equation:

```wolfram
eqn=MassTransportPDEComponent[vars,pars]==Γ_flux
(* Output *)
c[x] If[x>0.5,1.,0.]+({{-0.01}}.c[x])==NeumannValue[1+c[x] NeumannBoundaryUnitNormal[{x}].{0},x==0]
```

Solve the PDE:

```wolfram
cfun=NDSolveValue[{eqn,Γ_C},c,x∈Ω];
```

Visualize the solution:

```wolfram
Show[Plot[cfun[x],x∈Ω,AxesLabel -> xc],Graphics[DashedGrayLine[0.5-0.10.560]Text[Material 1, 0.2560]Text[Material 2, 0.7560]]]
```

*([Graphics])*

More information about the mass transport equation and boundary conditions can be found in the [Mass Transport tutorial](https://reference.wolfram.com/language/PDEModels/tutorial/MassTransport/MassTransport.html) and the [MassTransport PDE models guide page](https://reference.wolfram.com/language/guide/MassTransportPDEModels.html).

#### Transient Gray-Scott

Set up the equations and specify coefficients:

```wolfram
eqn = {u[t,x,y] (f+v[t,x,y]^2)+(-c1 u[t,x,y])+u^(1,0,0)[t,x,y]==f,v[t,x,y] (f+k-u[t,x,y] v[t,x,y])+(-c2 v[t,x,y])+v^(1,0,0)[t,x,y]==0}//. {c1 -> 2*10^-5, c2 -> c1/4, f -> 1/25, k -> 3/50};
```

As initial conditions, set $u=1/2$ and $v=1$ if $x^{2}+y^{2}<=1/40$, else 0:

```wolfram
ics={u[0,x,y]==1/2,v[0,x,y]==If[x^2+y^2<=1/40,1,0]};
```

Set both $u$ and $v$ to zero on the boundary:

```wolfram
bcs=DirichletCondition[{u[t,x,y]==0,v[t,x,y]==0},True];
```

Solve the system of equations on a refined finite element mesh:

```wolfram
{ufun,vfun}=NDSolveValue[{eqn,bcs,ics},{u,v},{x,y}∈Disk[],{t,0,3000},Method->{"PDEDiscretization"->{"MethodOfLines","SpatialDiscretization"->{"FiniteElement","MeshOptions"->{"MaxCellMeasure"->0.003}}}}];
```

Visualize the result:

```wolfram
frames=Table[ContourPlot[vfun[t,x,y],{x,-1,1},{y,-1,1},Contours -> Range[Part[#, 1], Part[#, 2], Part[#, 2], -, Part[#, 1], ), /, 4], &, ), [, MinMax[vfun[ValuesOnGrid]]SelectWithContents -> TrueSelectable -> False]],{t,100,2000,25}];
ListAnimate[frames,SaveDefinitions -> True]
```

#### Transient Landau-Ginzburg

Set up the equations:

```wolfram
eqn={u^(1,0)[t,x]+(-u[t,x])+(-v[t,x])-u[t,x]+(u[t,x]-2 v[t,x]) (u[t,x]^2+v[t,x]^2)==0,v^(1,0)[t,x]+u[t,x]+(-v[t,x])-v[t,x]+(2 u[t,x]+v[t,x]) (u[t,x]^2+v[t,x]^2)==0};
```

Set up initial conditions such that machine underflow is avoided:

```wolfram
exp[x_]:=If[x<Log[$MinMachineNumber],0.,Exp[x]]
ics={u[0,x]==1-exp[-x^2],v[0,x]==exp[-(10-x)^2]};
```

Specify boundary conditions:

```wolfram
bcs={u[t,-50]==1,u[t,50]==1,v[t,-50]==0,v[t,50]==0};
```

Solve the equations on a mesh with specified spacing:

```wolfram
{ufun,vfun}=NDSolveValue[{eqn,ics,bcs},{u,v},{t,0,40},{x,-50,50},Method->{"MethodOfLines","SpatialDiscretization"->{"FiniteElement","MeshOptions"->{"MaxCellMeasure"->0.4}}}];
```

Visualize the result:

```wolfram
DensityPlot[Evaluate[Sqrt[ufun[t,x]^2+vfun[t,x]^2]],{x,-50,50},{t,0,40},PlotRange -> All, MaxRecursion -> 8, PlotPoints -> 100, Frame -> None, Axes -> None, ImageSize -> Medium]
(* Output *)
![image](img/image_017.png)
```

#### Black-Scholes Equation

The goal is to determine the price of a European vanilla call option using the Black-Scholes model, given that both the underlying asset price $S$ and the strike price $K$ are $$100$, the risk-free interest rate $r$ is $5 %$, the asset's volatility $\sigma$ is $20 %$, and the option has a maturity period $t_{end}$ of one year.

Define the parameters:

```wolfram
pars={tEnd->1,k->100,σ->0.2,r-> 0.05}
(* Output *)
{tEnd->1,k->100,σ->0.2,r->0.05}
```

Define the bounds for the asset's price $S$:

```wolfram
{lower,upper}={0,3k}/.pars
(* Output *)
{0,300}
```

Define the PDE and replace parameters:

```wolfram
BSPDE=-r V[t,s]+r s V^(0,1)[t,s]+(1)/(2) s^2 σ^2 V^(0,2)[t,s]+V^(1,0)[t,s]==0/.pars
(* Output *)
-0.05 V[t,s]+0.05 s V^(0,1)[t,s]+0.020000000000000004 s^2 V^(0,2)[t,s]+V^(1,0)[t,s]==0
```

Set up the terminal condition:

```wolfram
terminalBC=V[tEnd,s]==Max[0,s-k]/.pars
(* Output *)
V[1,s]==Max[0,-100+s]
```

The boundary conditions are chosen so that at the lower bound, where the asset price is very low, the option is worth zero.The condition at the upper bound comes from the put-call parity.

Specify boundary conditions:

```wolfram
bcs={DirichletCondition[V[t,s]==upper-k ℯ^(-r(tEnd-t)),s==upper],DirichletCondition[V[t,s]==0,s==lower]}/.pars
(* Output *)
{DirichletCondition[V[t,s]==300-100 ℯ^(-0.05 (1-t)),s==300],DirichletCondition[V[t,s]==0,s==0]}
```

Solve the PDE on a refined mesh:

```wolfram
sol=NDSolveValue[{BSPDE,terminalBC,bcs},V,{t,0,tEnd}/.pars,{s,lower,upper},Method->{"MethodOfLines","SpatialDiscretization"->{"FiniteElement","MeshOptions"->{"MaxCellMeasure"->0.1}}}]
(* Output *)
InterpolatingFunction[]
```

Check the value of the option for an asset price of $100$ at time $t=0$:

```wolfram
sol[0,100]
(* Output *)
10.450516851337541
```

Plot the obtained numerical solution for the whole time range:

```wolfram
Plot3D[sol[t,s],{s,50,150},{t,0,1},ColorFunction -> BlueGreenYellow, AxesLabel -> stock price, stimeoption price]
(* Output *)
![image](img/image_019.png)
```

Define a function to get the price of the option using the [FinancialDerivative](https://reference.wolfram.com/language/ref/FinancialDerivative.html) function:

```wolfram
optionPrice[s0_]:=FinancialDerivative@@ReplaceAll[{{"European","Call"}, {"StrikePrice"-> k, "Expiration"->tEnd}, {"InterestRate"-> r, "Volatility" -> σ ,"CurrentPrice"->s0}},pars]
```

Inspect the conventionally true value given by [FinancialDerivative](https://reference.wolfram.com/language/ref/FinancialDerivative.html):

```wolfram
optionPrice[100]
(* Output *)
10.450583572185565
```

Calculate the error between the two values:

```wolfram
trueVal=optionPrice[100];
numericalVal=sol[0,100];
PercentForm[(Abs[trueVal-numericalVal])/(trueVal)]
(* Output *)
"0.0006384%"
```

Plot the difference between the numerical solution at time $0$ and the conventionally true solution:

```wolfram
Plot[optionPrice[s0]-sol[0,s0],{s0,k-100,k+100}/.pars,AxesLabel -> stock priceoption price, PlotLabel -> Difference, PlotRange -> All]
```

*([Graphics])*

### Properties & Relations

[NDSolveValue](https://reference.wolfram.com/language/ref/NDSolveValue.html) gives the solution function, while [NDSolve](https://reference.wolfram.com/language/ref/NDSolve.html) returns the result as a rule:

```wolfram
ysol=NDSolveValue[{y'[x]==y[x]Cos[x+y[x]],y[0]==1},y,{x,0,30}]
(* Output *)
InterpolatingFunction[...]
```

```wolfram
sol=NDSolve[{y'[x]==y[x]Cos[x+y[x]],y[0]==1},y,{x,0,30}]
(* Output *)
{{y->InterpolatingFunction[...]}}
```

The following two are equivalent:

```wolfram
ysol===First[y/.sol]
(* Output *)
True
```

```wolfram
{Plot[ysol[x],{x,0,30}],Plot[y[x]/.sol,{x,0,30}]}
(* Output *)
{[Graphics],[Graphics]}
```

[NDSolveValue](https://reference.wolfram.com/language/ref/NDSolveValue.html) can directly return the result of the solution function at a specified point:

```wolfram
NDSolveValue[{y'[x]==y[x]Cos[x+y[x]],y[0]==1},y[1.],{x,0,30}]
(* Output *)
0.9913873390555912
```

With [NDSolve](https://reference.wolfram.com/language/ref/NDSolve.html):

```wolfram
y[1.]/.First@NDSolve[{y'[x]==y[x]Cos[x+y[x]],y[0]==1},y,{x,0,30}]
(* Output *)
0.9913873390555912
```

[NDSolveValue](https://reference.wolfram.com/language/ref/NDSolveValue.html) can directly return the solution evaluated at several points:

```wolfram
NDSolveValue[{y'[x]==y[x]Cos[x+y[x]],y[0]==1},Table[y[i],{i,10,30,10}],{x,0,30}]
(* Output *)
{0.06434901313534486,0.19102584812031645,0.019301764000191374}
```

With [NDSolve](https://reference.wolfram.com/language/ref/NDSolve.html):

```wolfram
Table[y[i],{i,10,30,10}]/.First@NDSolve[{y'[x]==y[x]Cos[x+y[x]],y[0]==1},y,{x,0,30}]
(* Output *)
{0.06434901313534486,0.19102584812031645,0.019301764000191374}
```

[NDSolveValue](https://reference.wolfram.com/language/ref/NDSolveValue.html) can directly return the result of any expression involving the solution:

```wolfram
f[x_]=NDSolveValue[{y'[x]==y[x]Cos[x+y[x]],y[0]==1},y''[x]Cos[x+y[x]],{x,0,30}]
(* Output *)
Cos[x+InterpolatingFunction[...][x]] InterpolatingFunction[...][x]
```

```wolfram
Plot[f[x],{x,0,30}]
```

*([Graphics])*

```wolfram
NDSolveValue[{y'[x]==y[x]Cos[x+y[x]],y[0]==1},y''[25]Cos[25+y[25]],{x,0,30}]
(* Output *)
0.06139060296818053
```

With [NDSolve](https://reference.wolfram.com/language/ref/NDSolve.html):

```wolfram
sol=First@NDSolve[{y'[x]==y[x]Cos[x+y[x]],y[0]==1},y,{x,0,30}];
```

```wolfram
y''[x]Cos[x+y[x]]/.sol
(* Output *)
Cos[x+InterpolatingFunction[...][x]] InterpolatingFunction[...][x]
```

```wolfram
y''[25]Cos[25+y[25]]/.sol
(* Output *)
0.06139060296818053
```

If an equation has multiple solutions, [NDSolveValue](https://reference.wolfram.com/language/ref/NDSolveValue.html) will pick one:

```wolfram
ysol=NDSolveValue[{y'[x]^2==(1-y[x]^2)(1-y[x]^2/2),y[0]==0},y,{x,0,1.5}]
(* Output *)
NDSolveValue
(* Output *)
InterpolatingFunction[...]
```

```wolfram
Plot[ysol[x],{x,0,1.5}]
```

*([Graphics])*

[NDSolve](https://reference.wolfram.com/language/ref/NDSolve.html) will get both solution branches:

```wolfram
sol=NDSolve[{y'[x]^2==(1-y[x]^2)(1-y[x]^2/2),y[0]==0},y[x],{x,0,1.5}];
Plot[Evaluate[y[x]/.sol],{x,0,1.5}]
```

*([Graphics])*

Numerically compute values of an integral at different points in an interval:

```wolfram
iv=Table[{x,NIntegrate[Sin[Exp[s]],{s,0,x}]},{x,0,3,.25}]
(* Output *)
{{0.,0.},{0.25,0.2259882299920528},{0.5,0.4730479976194906},{0.75,0.709540439767007},{1.,0.8749571987803856},{1.25,0.8879990172574581},{1.5,0.7120350249749674},{1.75,0.49508378585658475},{2.,0.5509351737392814},{2.25,0.7284694489767151},{2.5,0.5519902414453264},{2.75,0.6877455856287154},{3.,0.6061244734187705}}
```

For functions of the independent variable, [NDSolve](https://reference.wolfram.com/language/ref/NDSolve.html) effectively gives an indefinite integral:

```wolfram
ysol=NDSolveValue[{y'[x]==Sin[Exp[x]],y[0]==0},y,{x,0,3}]
(* Output *)
InterpolatingFunction[...]
```

```wolfram
Plot[ysol[x],{x,0,3},Epilog->Map[Point,iv]]
```

*([Graphics])*

Finding an event is related to finding a root of a function of the solution:

```wolfram
ev[t_?NumberQ]:=NDSolveValue[{y''[x]+y[x]==0,y[0]==1,y'[0]==0},y,{x,t}][t]
```

```wolfram
Plot[ev[t],{t,0,3}]
```

*([Graphics])*

```wolfram
FindRoot[ev[t],{t,1}]
(* Output *)
{t->1.570796364662529}
```

Event location finds the root accurately and efficiently:

```wolfram
NDSolveValue[{y''[x]+y[x]==0,y[0]==1,y'[0]==0,WhenEvent[y[x]==0,"StopIntegration"]},y,{x,∞}]
(* Output *)
InterpolatingFunction[...]
```

This gives $y(10)$ as a function of $y^{'}(0)$ for a differential equation:

```wolfram
y_10[yp0_?NumberQ]:=NDSolveValue[{y''[x]+Sin[y[x]]==0,y[0]==1,y'[0]==yp0},y,{x,10}][10]
```

```wolfram
Plot[y_10[y_0^′],{y_0^′,-2,2}]
```

*([Graphics])*

Find a root of $y_{10}(y_{0}^{'})=0$:

```wolfram
FindRoot[y_10[y_0^′]==0,{y_0^′,-1.5}]
(* Output *)
{y->-1.5608322791316405}
```

Solve the equivalent boundary value problem:

```wolfram
ysol=NDSolveValue[{y''[x]+Sin[y[x]]==0,y[0]==1,y[10]==0},y,x];Plot[{ysol[x],ysol'[x]},{x,0,10}]
```

*([Graphics])*

Use [NDSolve](https://reference.wolfram.com/language/ref/NDSolve.html) as a solver for a [SystemModel](https://reference.wolfram.com/language/ref/SystemModel.html):

```wolfram
sim=SystemModelSimulate[[Graphics],Method->"NDSolve"]
(* Output *)
SystemModelSimulationData[...]
```

Plot variables from the simulation result:

```wolfram
SystemModelPlot[sim,{y,z}]
(* Output *)
![image](img/image_021.png)
```

Use [SystemModel](https://reference.wolfram.com/language/ref/SystemModel.html) to model larger hierarchical models:

```wolfram
sim=SystemModelSimulate[[Graphics]]
(* Output *)
SystemModelSimulationData[...]
```

Plot the tank levels in the tank system over time:

```wolfram
SystemModelPlot[sim,{"tank1.h","tank2.h","tank3.h"},PlotRange->All]
(* Output *)
![image](img/image_023.png)
```

### Possible Issues

Many [NDSolveValue](https://reference.wolfram.com/language/ref/NDSolveValue.html) messages have specific message reference pages implemented. See how to access them in the [Understand Error Messages](https://reference.wolfram.com/language/workflow/UnderstandErrorMessages.html) workflow.

#### Multiple Solutions

For equations with multiple solutions, [NDSolveValue](https://reference.wolfram.com/language/ref/NDSolveValue.html) will return only one solution:

```wolfram
NDSolveValue[{y^′[x]^2-y[x]^2==0,y[0]^2==4},y[x],{x,1}]
(* Output *)
NDSolveValue
(* Output *)
InterpolatingFunction[...][x]
```

Use [NDSolve](https://reference.wolfram.com/language/ref/NDSolve.html) to get all the solutions:

```wolfram
NDSolve[{y^′[x]^2-y[x]^2==0,y[0]^2==4},y[x],{x,1}]
(* Output *)
{{y[x]->InterpolatingFunction[...][x]},{y[x]->InterpolatingFunction[...][x]},{y[x]->InterpolatingFunction[...][x]},{y[x]->InterpolatingFunction[...][x]}}
```

```wolfram
Plot[Evaluate[y[x]/.%],{x,0,1}]
```

*([Graphics])*

#### Numerical Error

The error tends to grow as one goes further from the initial conditions:

```wolfram
ysol=NDSolveValue[{y''[x]+y[x]==0,y[0]==1,y'[0]==0},y,{x,0,100}]
(* Output *)
InterpolatingFunction[...]
```

Find the difference between numerical and exact solutions:

```wolfram
Plot[ysol[x]-Cos[x],{x,0,100}]
```

*([Graphics])*

Error for a nonlinear equation:

```wolfram
ysol=NDSolveValue[{y''[x]==y[x](y[x]^2-3/2),y[0]==0,y'[0]==1},y,{x,0,10}]
(* Output *)
InterpolatingFunction[...]
```

```wolfram
Plot[ysol[x]-JacobiSN[x,1/2],{x,0,10}]
```

*([Graphics])*

For high-order methods, the default interpolation may have large errors between steps:

```wolfram
ysol=NDSolveValue[{y''[x]==y[x](y[x]^2-3/2),y[0]==0,y'[0]==1},y,{x,0,10},Method->{"ExplicitRungeKutta","DifferenceOrder"->8}]
(* Output *)
InterpolatingFunction[...]
```

```wolfram
Plot[(ysol''[x]-ysol[x](ysol[x]^2-3/2)),{x,0,10},PlotRange->All]
```

*([Graphics])*

Interpolation with the order corresponding to the method reduces the error between steps:

```wolfram
ysol2=NDSolveValue[{y''[x]==y[x](y[x]^2-3/2),y[0]==0,y'[0]==1},y,{x,0,10},Method->{"ExplicitRungeKutta","DifferenceOrder"->8},InterpolationOrder->All]
(* Output *)
InterpolatingFunction[...]
```

```wolfram
Plot[(ysol2''[x]-ysol2[x](ysol2[x]^2-3/2)),{x,0,10},PlotRange->All]
```

*([Graphics])*

#### Differential Algebraic Equations

Here is a system of differential-algebraic equations:

```wolfram
DAE={x_1^′[t]==x_3[t],x_2[t](1-x_2[t])==0,x_1[t]x_2[t]+x_3[t](1-x_2[t])==t};
```

Find the solution with $x_{1}^{'}(0)=1$:

```wolfram
NDSolveValue[{DAE,x_1'[0]==1},{x_1,x_2,x_3},{t,0,1}]
(* Output *)
{InterpolatingFunction[...],InterpolatingFunction[...],InterpolatingFunction[...]}
```

[NDSolveValue](https://reference.wolfram.com/language/ref/NDSolveValue.html) may change the specified initial conditions if it cannot find the solution with $x_{1}^{'}(0)=0$:

```wolfram
NDSolveValue[{DAE,x_1'[0]==0},{x_1,x_2,x_3},{t,0,1}]
(* Output *)
NDSolveValue
(* Output *)
{InterpolatingFunction[...],InterpolatingFunction[...],InterpolatingFunction[...]}
```

[NDSolveValue](https://reference.wolfram.com/language/ref/NDSolveValue.html) is limited to index 1, but the solution with $x_{1}^{'}(0)=0$ has index 2:

```wolfram
NDSolveValue[{DAE,x_1[0] == 0, x_2[0]== 1, x_3[0] == 0},{x_1,x_2,x_3}, {t,0,1}]
(* Output *)
NDSolveValue
(* Output *)
{InterpolatingFunction[...],InterpolatingFunction[...],InterpolatingFunction[...]}
```

The default method may not be able to converge to the default tolerances:

```wolfram
s=100;
dae={x'[t]==-(xs[t]-ys[t]),y'[t]==(xs[t]-ys[t]),xs[t]==s x[t],ys[t]==s y[t],x[0]==1,y[0]==0,xs[0]==s,ys[0]==0};
```

```wolfram
NDSolveValue[dae,{x,y},{t,0,1}]
(* Output *)
NDSolveValue
(* Output *)
{InterpolatingFunction[...],InterpolatingFunction[...]}
```

With lower [AccuracyGoal](https://reference.wolfram.com/language/ref/AccuracyGoal.html) and [PrecisionGoal](https://reference.wolfram.com/language/ref/PrecisionGoal.html) settings, a solution is found:

```wolfram
NDSolveValue[dae,{x,y},{t,0,1},AccuracyGoal->6,PrecisionGoal->6]
(* Output *)
{InterpolatingFunction[...],InterpolatingFunction[...]}
```

The `"StateSpace"` time integration method can solve this with default tolerances:

```wolfram
NDSolveValue[dae,{x,y},{t,0,1},Method->{"TimeIntegration"->"StateSpace"}]
(* Output *)
{InterpolatingFunction[...],InterpolatingFunction[...]}
```

#### Partial Differential Equations

A large collection of PDE models from various fields with extensive explanation can found in the [PDE models overview](https://reference.wolfram.com/language/PDEModels/tutorial/PDEModelsOverview.html).

Define a nonlinear PDE:

```wolfram
L=25;
wpde={∂_t,tu[t,x]==∂_x,xu[t,x]+(1-u[t,x]^2)(1+2u[t,x]),u[0,x]==ℯ^(-x^2),u^(1,0)[0,x]==0,u[t,-L]==u[t,L]};
```

The spatial discretization is based on the initial value, which varies less than the final value:

```wolfram
usol1=NDSolveValue[wpde,u,{t,0,L},{x,-L,L}]
(* Output *)
NDSolveValue
(* Output *)
InterpolatingFunction[...]
```

By increasing the minimal number of spatial grid points, you can accurately compute the final value:

```wolfram
usol2=NDSolveValue[wpde,u,{t,0,L},{x,-L,L},Method->{"MethodOfLines", "SpatialDiscretization"->{"TensorProductGrid", "MinPoints"->1000}}]
(* Output *)
InterpolatingFunction[]
```

The plot demonstrates the onset of a spatially more complex solution:

```wolfram
DensityPlot[usol2[-t,x],{x,-L,L},{t,0,-L},PlotPoints->50]
(* Output *)
![image](img/image_025.png)
```

Define a heat equation with an initial value that is a step function:

```wolfram
hpde={D[u[t,x],t]==D[u[t,x],x,x],u[0,x]==UnitStep[1/2-x],u[t,0]==1,u[t,1]==0};
```

Discontinuities in the initial value may result in too many spatial grid points:

```wolfram
Plot3D[Evaluate[NDSolveValue[hpde,u[t,x],{t,0,.01},{x,0,1}]],{x,0,1},{t,0,.01}]
(* Output *)
NDSolveValue
(* Output *)
![image](img/image_027.png)
```

Setting the number of spatial grid points smaller results in essentially as good a solution:

```wolfram
Plot3D[Evaluate[NDSolveValue[hpde,u[t,x],{t,0,.01},{x,0,1}]],{x,0,1},{t,0,.01},Method->{"MethodOfLines","SpatialDiscretization"->{"TensorProductGrid","MaxPoints"->100}}]
(* Output *)
NDSolveValue
(* Output *)
![image](img/image_029.png)
```

Define a Laplace equation with initial values:

```wolfram
lpde={D[u[x,y],x,x]+D[u[x,y],y,y]==0,u[x,0]==Sin[Pi x],Derivative[0,1][u][x,0]==0,u[0,y]==u[1,y]==0};
```

The solver only works for equations well posed as initial value (Cauchy) problems:

```wolfram
usol=NDSolveValue[lpde,u,{x,0,1},{y,0,1}]
(* Output *)
NDSolveValue
(* Output *)
InterpolatingFunction[...]
```

The ill-posedness shows up as numerical instability:

```wolfram
Plot3D[Evaluate[usol[x,y]],{x,0,1},{y,0,1},PlotRange->All]
(* Output *)
![image](img/image_031.png)
```

#### Boundary Value Problems

This finds a trivial solution of a boundary value problem:

```wolfram
ysol=NDSolveValue[{y''[x]+Sin[y[x]]==0,y[0]==0,y[10]==0},y,x];Plot[{ysol[x],ysol'[x]},{x,0,10}]
```

*([Graphics])*

You can find other solutions by giving starting conditions for the solution search:

```wolfram
Table[
ysol=NDSolveValue[{y''[x]+Sin[y[x]]==0,y[0]==0,y[10]==0},y,x,Method->{"Shooting","StartingInitialConditions"->{y[0]==0,y'[0]==i}}];
Plot[{ysol[x],ysol'[x]},{x,0,10}],{i,1,2}]
(* Output *)
{[Graphics],[Graphics]}
```

### Neat Examples

Solve a maze by solving a Laplace equation and computing the gradient of the solution:

```wolfram
mesh=ImageMesh[[Graphics]];
if=NDSolveValue[{Laplacian[u[x,y],{x,y}]==0,DirichletCondition[u[x,y]==1,x<=6],DirichletCondition[u[x,y]==0,x>=595]},u,{x,y}∈mesh];
ContourPlot[Evaluate[Sqrt[Total[Grad[if[x,y],{x,y}]^2]]],{x,y}∈mesh,Contours->2,PlotPoints->30]
(* Output *)
![image](img/image_033.png)
```

## Tech Notes ▪Numerical Differential Equations ▪Numerical Solution of Differential Equations ▪Advanced Numerical Differential Equation Solving ▪NDSolve Options for the Finite Element Method ▪Numerical Mathematics ▪Implementation Notes: Numerical and Related Functions

## Related Guides ▪Solvers over Regions ▪Partial Differential Equations ▪Differential Equations ▪Differential Operators ▪Symbolic Vectors, Matrices and Arrays

## History Introduced in 2012 (9.0) | Updated in 2014 (10.0) ▪ 2019 (12.0)
