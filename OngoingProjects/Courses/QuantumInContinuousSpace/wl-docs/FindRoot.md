# FindRoot | [SpanFromLeft]

> [FindRoot](https://reference.wolfram.com/language/ref/FindRoot.html)[*f*,{*x*,*x*_0}]  — searches for a numerical root of `*f*`, starting from the point `*x*=*x*_0
> [FindRoot](https://reference.wolfram.com/language/ref/FindRoot.html)[*lhs*==*rhs*,{*x*,*x*_0}]  — searches for a numerical solution to the equation `*lhs*==*rhs*`.
> [FindRoot](https://reference.wolfram.com/language/ref/FindRoot.html)[{*f*_1,*f*_2,…},{{*x*,*x*_0},{*y*,*y*_0},…}]  — searches for a simultaneous numerical root of all the `*f*_*i*`.
> [FindRoot](https://reference.wolfram.com/language/ref/FindRoot.html)[{*eqn*_1,*eqn*_2,…},{{*x*,*x*_0},{*y*,*y*_0},…}] — searches for a numerical solution to the simultaneous equations `*eqn*_*i*`.

## Details and Options

If the starting point for a variable is given as a list, the values of the variable are taken to be lists with the same dimensions.

[FindRoot](https://reference.wolfram.com/language/ref/FindRoot.html) returns a list of replacements for `*x*`, `*y*`, … , in the same form as obtained from [Solve](https://reference.wolfram.com/language/ref/Solve.html).

[FindRoot](https://reference.wolfram.com/language/ref/FindRoot.html) first localizes the values of all variables, then evaluates `*f*` with the variables being symbolic, and then repeatedly evaluates the result numerically.

[FindRoot](https://reference.wolfram.com/language/ref/FindRoot.html) has attribute [HoldAll](https://reference.wolfram.com/language/ref/HoldAll.html), and effectively uses [Block](https://reference.wolfram.com/language/ref/Block.html) to localize variables.

[FindRoot](https://reference.wolfram.com/language/ref/FindRoot.html)[*lhs*==*rhs*,{*x*,*x*_0,*x*_1}] searches for a solution using `*x*_0 and `*x*_1 as the first two values of `*x*`, avoiding the use of derivatives.

[FindRoot](https://reference.wolfram.com/language/ref/FindRoot.html)[*lhs*==*rhs*,{*x*,*x*_*start*,*x*_*min*,*x*_*max*}] searches for a solution, stopping the search if `*x*` ever gets outside the range `*x*_*min*` to `*x*_*max*`.

If you specify only one starting value of `*x*`, [FindRoot](https://reference.wolfram.com/language/ref/FindRoot.html) searches for a solution using Newton methods. If you specify two starting values, [FindRoot](https://reference.wolfram.com/language/ref/FindRoot.html) uses a variant of the secant method.

If all equations and starting values are real, then [FindRoot](https://reference.wolfram.com/language/ref/FindRoot.html) will search only for real roots. If any are complex, it will also search for complex roots.

You can always tell [FindRoot](https://reference.wolfram.com/language/ref/FindRoot.html) to search for complex roots by adding `0.[I](https://reference.wolfram.com/language/ref/I.html)` to the starting value.

The following options can be given:

| [AccuracyGoal](https://reference.wolfram.com/language/ref/AccuracyGoal.html) | [Automatic](https://reference.wolfram.com/language/ref/Automatic.html) | the accuracy sought |
| --- | --- | --- |
| [EvaluationMonitor](https://reference.wolfram.com/language/ref/EvaluationMonitor.html) | [None](https://reference.wolfram.com/language/ref/None.html) | expression to evaluate whenever equations are evaluated |
| Jacobian | [Automatic](https://reference.wolfram.com/language/ref/Automatic.html) | the Jacobian of the system |
| [MaxIterations](https://reference.wolfram.com/language/ref/MaxIterations.html) | 100 | maximum number of iterations to use |
| $Method$ | [Automatic](https://reference.wolfram.com/language/ref/Automatic.html) | method to be used |
| [PrecisionGoal](https://reference.wolfram.com/language/ref/PrecisionGoal.html) | [Automatic](https://reference.wolfram.com/language/ref/Automatic.html) | the precision sought |
| [StepMonitor](https://reference.wolfram.com/language/ref/StepMonitor.html) | [None](https://reference.wolfram.com/language/ref/None.html) | expression to evaluate whenever a step is taken |
| [WorkingPrecision](https://reference.wolfram.com/language/ref/WorkingPrecision.html) | [MachinePrecision](https://reference.wolfram.com/language/ref/MachinePrecision.html) | the precision to use in internal computations |

The default settings for [AccuracyGoal](https://reference.wolfram.com/language/ref/AccuracyGoal.html) and [PrecisionGoal](https://reference.wolfram.com/language/ref/PrecisionGoal.html) are [WorkingPrecision](https://reference.wolfram.com/language/ref/WorkingPrecision.html)/2.

The setting for [AccuracyGoal](https://reference.wolfram.com/language/ref/AccuracyGoal.html) specifies the number of digits of accuracy to seek in both the value of the position of the root, and the value of the function at the root.

The setting for [PrecisionGoal](https://reference.wolfram.com/language/ref/PrecisionGoal.html) specifies the number of digits of precision to seek in the value of the position of the root.

[FindRoot](https://reference.wolfram.com/language/ref/FindRoot.html) continues until either of the goals specified by [AccuracyGoal](https://reference.wolfram.com/language/ref/AccuracyGoal.html) or [PrecisionGoal](https://reference.wolfram.com/language/ref/PrecisionGoal.html) is achieved.

If [FindRoot](https://reference.wolfram.com/language/ref/FindRoot.html) does not succeed in finding a solution to the accuracy you specify within [MaxIterations](https://reference.wolfram.com/language/ref/MaxIterations.html) steps, it returns the most recent approximation to a solution that it found. You can then apply [FindRoot](https://reference.wolfram.com/language/ref/FindRoot.html) again, with this approximation as a starting point.

## Examples

### Basic Examples

Find a root of $e^{x}+sin(x)$ near $x=0$:

```wolfram
FindRoot[Sin[x]+Exp[x],{x,0}]
(* Output *)
{x->-0.5885327439818611}
```

Find a solution to $cos(x)=x$ near $x=0$:

```wolfram
FindRoot[Cos[x]==x,{x,0}]
(* Output *)
{x->0.7390851332151607}
```

Solve a nonlinear system of equations:

```wolfram
FindRoot[{Exp[x-2]==y,y^2==x},{{x,1},{y,1}}]
(* Output *)
{x->0.019026016103714054,y->0.13793482556524314}
```

### Scope

Find the solution of a system of two nonlinear equations:

```wolfram
FindRoot[{y==Exp[x],x+y==2},{{x,1},{y,1}}]
(* Output *)
{x->0.4428544010023886,y->1.5571455989976115}
```

Find a root for a three-component function of three variables:

```wolfram
FindRoot[{Sin[x+y],Cos[x-y],x^2+y^2-z},{{x,1},{y,0},{z,0}}]
(* Output *)
{x->0.7853981633974483,y->-0.7853981633974483,z->1.2337005501361697}
```

You can cause the search to use complex values by giving a complex starting value:

```wolfram
FindRoot[Zeta[z],{z,1/2+14 I}]
(* Output *)
{z->0.5000000000000016+14.134725141734695 ⅈ}
```

When the function is complex for real input, a real starting value may give a complex result:

```wolfram
FindRoot[(Cos[z+I]-2)(z+2),{z,1}]
(* Output *)
{z->-5.9471781882505156×10^-24+0.31695789692481674 ⅈ}
```

```wolfram
FindRoot[(Cos[z+I]-2)(z+2),{z,-1}]
(* Output *)
{z->-1.9999999999999993+6.700856866951783×10^-16 ⅈ}
```

### Generalizations & Extensions

A variable may be considered as vector valued if specified in the starting value:

```wolfram
A={{1,2,3},{4,5,6},{6,7,8}};
```

```wolfram
FindRoot[{A.x==λ x,x.x==1},{{x,RandomReal[1,3]},{λ,1}}]
(* Output *)
{x->{0.832050293744123,1.2861385633584948×10^-9,-0.5547001971158102},λ->-0.9999999983825539}
```

### Options

#### AccuracyGoal

Change tolerances for error estimates:

```wolfram
10-x/.FindRoot[Sin[x-10]-x+10,{x,0}]
(* Output *)
1.9716335408759278×10^-7
```

Relax error tolerances for stopping:

```wolfram
10-x/.FindRoot[Sin[x-10]-x+10,{x,0},AccuracyGoal->4,PrecisionGoal->4]
(* Output *)
0.0014761758665553515
```

Make estimated relative distance to the root the main criterion for stopping:

```wolfram
10-x/.FindRoot[Sin[x-10]-x+10,{x,0},AccuracyGoal->Infinity,PrecisionGoal->8]
(* Output *)
8.615632651753913×10^-9
```

#### DampingFactor

`DampingFactor` can be used to help speed convergence to higher-order roots:

```wolfram
Block[{s=0},{FindRoot[Sin[x]-1,{x,1.5},StepMonitor:>s++],s}]
(* Output *)
{{x->1.5707963102749631},22}
```

```wolfram
Block[{s=0},{FindRoot[Sin[x]-1,{x,1.5},DampingFactor->2,StepMonitor:>s++],s}]
(* Output *)
{{x->1.5707963267931118},2}
```

#### EvaluationMonitor

[EvaluationMonitor](https://reference.wolfram.com/language/ref/EvaluationMonitor.html) can be used to keep track of function evaluations used:

```wolfram
{res,{evx}}=Reap[FindRoot[x^2==Exp[x],{x,0},EvaluationMonitor:>Sow[x]]]
(* Output *)
{{x->-0.7034674224983924},{{0.,-1.,-0.7330436052454454,-0.703807786324133,-0.7034674683317975,-0.7034674224983924}}}
```

```wolfram
Show[Plot[{x^2,Exp[x]},{x,-1,0}],ListPlot[{Transpose[{evx,evx^2}],Transpose[{evx,Exp[evx]}]}]]
```

*([Graphics])*

#### Jacobian

Specify the Jacobian for a "black-box" function:

```wolfram
n=1000;
A=SparseArray[{{i_,i_}->-2.,{i_,j_}/;Abs[i-j]==1->1.},{n,n},0.];
blackbox[x_?VectorQ]:=n^2 A.x+x-x^3-1;
```

```wolfram
J[x_?VectorQ]:=Module[{d=1-3 x^2},n^2 A+SparseArray[{i_,i_}:>d[[i]],{n,n}]]
```

```wolfram
sc=ConstantArray[0.,n];
Block[{e=0},Timing[FindRoot[blackbox[x],{x,sc},Jacobian->J[x],EvaluationMonitor:>e++];e]]
(* Output *)
{0.023126,4}
```

Without a specified Jacobian, extra evaluations are used to compute finite differences:

```wolfram
Block[{e=0},Timing[FindRoot[blackbox[x],{x,sc},EvaluationMonitor:>e++];e]]
(* Output *)
{0.538854,3007}
```

If you just know the sparse form, specifying the sparse pattern template saves evaluations:

```wolfram
sp=SparseArray[{{i_,j_}/;Abs[i-j]<=1->_},{n,n}]
(* Output *)
SparseArray[...]
```

```wolfram
Block[{e=0},Timing[FindRoot[blackbox[x],{x,sc},Jacobian->{"FiniteDifference",Sparse->sp},EvaluationMonitor:>e++];e]]
(* Output *)
{0.013116,16}
```

Inspect the number of Jacobian evaluations needed by different methods:

```wolfram
Block[{j=0},FindRoot[blackbox[x],{x,sc},Jacobian->{J[x],EvaluationMonitor:>j++}];"Jacobian Evaluations"->j]
(* Output *)
"Jacobian Evaluations"->3
```

```wolfram
Block[{j=0},FindRoot[blackbox[x],{x,sc},
Method->"AffineCovariantNewton",Jacobian->{J[x],EvaluationMonitor:>j++}];"Jacobian Evaluations"->j]
(* Output *)
"Jacobian Evaluations"->1
```

#### MaxIterations

Limit or increase the number of steps taken:

```wolfram
FindRoot[Exp[1/(x^2-1)]-10^-8,{x,.1},MaxIterations->10]
(* Output *)
FindRoot
(* Output *)
{x->0.9665826038411666}
```

The default number of iterations is 100:

```wolfram
FindRoot[Exp[1/(x^2-1)],{x,.1}]
(* Output *)
FindRoot
(* Output *)
{x->0.9953265163164798}
```

Eventually the algorithm stalls out since this mollifier function has all derivatives 0 at $x=1$:

```wolfram
FindRoot[Exp[1/(x^2-1)],{x,.1},MaxIterations->1000]
(* Output *)
{x->0.9993287975228579}
```

#### Method

Method options are also explained in [Unconstrained Optimization](https://reference.wolfram.com/language/tutorial/UnconstrainedOptimizationOverview.html).

Find a root for$f(x)=tan^{-1}(1000 cos(x))$ using different methods:

```wolfram
f[x_]=ArcTan[1000 Cos[x]];
Plot[f[x],{x,0,12},Exclusions->None]
```

*([Graphics])*

Define a function that monitors the steps and evaluations used by [FindRoot](https://reference.wolfram.com/language/ref/FindRoot.html):

```wolfram
monitoredFindRoot[args__] := Module[{s = 0, e = 0,j = 0},
{FindRoot[args, StepMonitor:>s++,EvaluationMonitor:>e++,Jacobian->{Automatic,EvaluationMonitor:>j++}],"Steps"->s,"Evaluations"->e,"Jacobian Evaluations"->j}]
```

The default (Newton's) method:

```wolfram
monitoredFindRoot[ArcTan[1000Cos[x]],{x,1}]
(* Output *)
{{x->10.995574287564276},"Steps"->9,"Evaluations"->13,"Jacobian Evaluations"->9}
```

Brent's root-bracketing method requiring two initial conditions bracketing the root:

```wolfram
monitoredFindRoot[ArcTan[1000Cos[x]],{x,1,2},Method->"Brent"]
(* Output *)
{{x->1.5707963267948966},"Steps"->14,"Evaluations"->16,"Jacobian Evaluations"->0}
```

Secant method, starting with two initial conditions:

```wolfram
monitoredFindRoot[ArcTan[1000Cos[x]],{x,1,2},Method->"Secant"]
(* Output *)
{{x->1.5707963267948966},"Steps"->12,"Evaluations"->27,"Jacobian Evaluations"->0}
```

Select the affine covariant Newton method:

```wolfram
monitoredFindRoot[ArcTan[1000Cos[x]],{x,1},Method->"AffineCovariantNewton"]
(* Output *)
{{x->1.5707963267964609},"Steps"->14,"Evaluations"->20,"Jacobian Evaluations"->13}
```

#### PrecisionGoal

Change tolerances for error estimates:

```wolfram
10-x/.FindRoot[Sin[x-10]-x+10,{x,0}]
(* Output *)
1.9716335408759278×10^-7
```

Relax error tolerances for stopping:

```wolfram
10-x/.FindRoot[Sin[x-10]-x+10,{x,0},AccuracyGoal->4,PrecisionGoal->4]
(* Output *)
0.0014761758665553515
```

Make estimated relative distance to the root the main criterion for stopping:

```wolfram
10-x/.FindRoot[Sin[x-10]-x+10,{x,0},AccuracyGoal->Infinity,PrecisionGoal->8]
(* Output *)
8.615632651753913×10^-9
```

#### StepMonitor

Monitor when iterative steps have been taken:

```wolfram
f[x_,y_]={x-1,10(y-Cos[2x]+1)};
```

```wolfram
{res,{stxy}}=Reap[FindRoot[f[x,y],{{x,-1},{y,-1}},StepMonitor:>Sow[{x,y}]]]
(* Output *)
{{x->1.,y->-1.4161468365471424},{{{-0.8,-0.6778957129244416},{-0.6200000000000001,-0.3531795967671845},{-0.4580000000000001,-0.07894799282493775},{-0.31220000000000003,0.12113552694504906},{-0.18098000000000006,0.243578580619981},{-0.06288200000000004,0.2963800521257964},{0.04340619999999998,0.29261629897930774},{0.13906558000000002,0.24639009924822125},{0.22515902200000004,0.17063263443660714},{0.3499057621544118,0.018518901242232877},{0.5885565636880948,-0.3819790474514864},{1.,-1.376345614207783},{1.,-1.4161468365471424}}}}
```

Show the steps on a contour plot of $‖f‖^{2}$:

```wolfram
ContourPlot[f[x,y].f[x,y],{x,-1,1},{y,-2,1},Epilog->{Red,Map[Point,stxy]}]
```

*([Graphics])*

Show steps (red) and evaluations (green). A step may require several evaluations:

```wolfram
{res,{exy}}=Reap[FindRoot[f[x,y],{{x,-1},{y,-1}},EvaluationMonitor:>Sow[{x,y}]]];
```

```wolfram
ContourPlot[f[x,y].f[x,y],{x,-1,1},{y,-2,1},Epilog->{{PointSize[0.04],Red,Map[Point,stxy]},{Green,Map[Point,exy]}}]
```

*([Graphics])*

#### WorkingPrecision

Find a root using 100-digit precision arithmetic:

```wolfram
FindRoot[Cos[x^2]-x,{x,1},WorkingPrecision->100]
(* Output *)
{x->0.8010707652092183662168678540865655855421529449709737551110566847790409145678013209735696971621054444230752988149312}
```

Find the root starting with machine precision and adaptively working up to precision 100:

```wolfram
Block[{prec=MachinePrecision},FindRoot[Cos[x^2]-x,{x,1.0},WorkingPrecision->100,StepMonitor:>If[Precision[x]≠prec,prec=Precision[x];Print["Increased precision to ",prec," at x = ",x];]]]
(* Output *)
"Increased precision to "25.075257498915995" at x = "0.80107076520930349497172655477694735993
(* Output *)
"Increased precision to "50.15051499783199" at x = "0.801070765209218366216867860104169246449814460829029839602555358317
(* Output *)
"Increased precision to "100." at x = "0.8010707652092183662168678540865655855421529449709737851800008407376871607365151432799043919751751341927025380985114
(* Output *)
{x->0.8010707652092183662168678540865655855421529449709737551110566847790409145678013209735696971621054444230752979310096}
```

### Applications

#### Computing Inverse Functions

For an isomorphism $f$, the inverse $s=f ^{-1}(t)$ is the root of $f (s)-t=0$:

```wolfram
inv[f_,s_]:=Function[{t},s/.FindRoot[f-t,{s,1}]]
```

An approximate inverse for the exponential function:

```wolfram
einv=inv[Exp[x],x]
(* Output *)
Function[{t$},x/.FindRoot[ℯ^x-t$,{x,1}]]
```

It is very close to the built-in [Log](https://reference.wolfram.com/language/ref/Log.html) function:

```wolfram
Plot[einv[x]-Log[x],{x,0,1},PlotRange->All]
```

*([Graphics])*

A "black-box" function giving the period of an oscillation:

```wolfram
f[a_?NumberQ]:=Module[{p=0},NDSolve[{x''[t]+x[t]+x[t]^3==0,x[0]==a,x'[0]==0},x,{t,∞},Method->{"EventLocator","Event"->x[t],"EventAction":>Throw[p=t,"StopIntegration"]}];
4p];
```

```wolfram
Plot[f[a],{a,0,4}]
```

*([Graphics])*

Plot its inverse:

```wolfram
finv=inv[f[a],a]
(* Output *)
Function[{t$},a/.FindRoot[f[a]-t$,{a,1}]]
```

```wolfram
Plot[finv[s],{s,3,6}]
```

*([Graphics])*

#### Solving Boundary Value Problems

Solve a boundary value problem $x^{''}+x^{'}+x^{3}=sin(4 t)$, $x(0)=x(1)=0$ using a shooting method:

```wolfram
x1[s_?NumberQ]:=First[x[1]/.NDSolve[{x''[t]+x'[t]+x[t]^3==Sin[4 t],x[0]==0,x'[0]==s},x,{t,0,1}]]
```

```wolfram
Plot[x1[s],{s,-10,10}]
```

*([Graphics])*

Use points on either side of the root to give bracketing starting values:

```wolfram
FindRoot[x1[s],{s,-1,0}]
(* Output *)
{s->-0.34432343230768403}
```

Plot the solution:

```wolfram
sol=NDSolve[{x''[t]+x'[t]+x[t]^3==Sin[4t],x[0]==0,x'[0]==s}/.%,x,{t,0,1}]
(* Output *)
{{x->InterpolatingFunction[...]}}
```

```wolfram
Plot[x[t]/.sol,{t,0,1}]
```

*([Graphics])*

Solve the boundary-value problem $\epsilon x'' + x + x^{3}=1$, $x(0)=x(1)=0$ with `*n*` collocation points:

```wolfram
n=100;
```

Consider as a first-order system $\{u,v \}' = f(\{u,v \})$:

```wolfram
f[{u_,v_}]:={v,(1-u-u^3)/ε};
```

Equations for collocation using the trapezoidal rule:

```wolfram
eqns=Flatten[Join[{u_0==0,u_n==0},Table[Thread[{u_i,v_i}=={u_i-1,v_i-1}+(1)/(2 n)(f[{u_i-1,v_i-1}]+f[{u_i,v_i}])],{i,1,n}]]];
```

Use 0 as a starting value:

```wolfram
sv=Flatten[Table[{{u_i,0},{v_i,0}},{i,0,n}],1];
```

Find a solution for a particular value of `ε`:

```wolfram
froot=FindRoot[eqns/.ε->0.01,sv];
sol=Table[u_i,{i,0,n}]/.froot;
ListLinePlot[sol]
```

*([Graphics])*

### Properties & Relations

For a polynomial system of equations, [NSolve](https://reference.wolfram.com/language/ref/NSolve.html) finds all solutions and [FindRoot](https://reference.wolfram.com/language/ref/FindRoot.html) finds one:

```wolfram
polysys=Table[Total[RandomInteger[{-9,9},5]RandomChoice[{x,y,z},5]^RandomInteger[{0,3},5]]==0,{3}]
(* Output *)
{5 x-6 z==0,-7-3 x^2+5 y^2-7 z^2==0,10+7 y^2+9 z-9 z^3==0}
```

[FindRoot](https://reference.wolfram.com/language/ref/FindRoot.html) will find a single solution using an iterative method:

```wolfram
FindRoot[polysys,Transpose[{{x,y,z},RandomComplex[1+I,3]}]]
(* Output *)
{x->-0.4488388150294671+1.0301560673611714 ⅈ,y->-0.8668817950625549+0.8385857527938984 ⅈ,z->-0.37403234585788925+0.8584633894676427 ⅈ}
```

[NSolve](https://reference.wolfram.com/language/ref/NSolve.html) will find all solutions using a direct method:

```wolfram
NSolve[polysys,{x,y,z}]
(* Output *)
{{x->3.0107442967255955,y->-3.9562022230467404,z->2.508953580604663},{x->3.0107442967255693,y->3.9562022230467293,z->2.5089535806046412},{x->-0.44883881502946704-1.030156067361171 ⅈ,y->0.8668817950625548+0.8385857527938978 ⅈ,z->-0.3740323458578892-0.8584633894676424 ⅈ},{x->-0.44883881502946704+1.030156067361171 ⅈ,y->0.8668817950625548-0.8385857527938978 ⅈ,z->-0.3740323458578892+0.8584633894676424 ⅈ},{x->-0.4488388150294671+1.0301560673611692 ⅈ,y->-0.8668817950625541+0.8385857527938962 ⅈ,z->-0.37403234585788925+0.858463389467641 ⅈ},{x->-0.4488388150294671-1.0301560673611692 ⅈ,y->-0.8668817950625541-0.8385857527938962 ⅈ,z->-0.37403234585788925-0.858463389467641 ⅈ}}
```

For equations involving parameters or exact solutions use [Solve](https://reference.wolfram.com/language/ref/Solve.html), [Reduce](https://reference.wolfram.com/language/ref/Reduce.html), or [FindInstance](https://reference.wolfram.com/language/ref/FindInstance.html):

```wolfram
eq=x Exp[2x+a]+3==0;
```

[Solve](https://reference.wolfram.com/language/ref/Solve.html) will return some solutions:

```wolfram
Solve[eq,x]
(* Output *)
Solve
(* Output *)
{{x->(1)/(2) ProductLog[-6 ℯ^-a]}}
```

[Reduce](https://reference.wolfram.com/language/ref/Reduce.html) will enumerate all solutions:

```wolfram
Reduce[eq,x]
(* Output *)
1∈Integers&&x==(1)/(2) ProductLog[1,-6 ℯ^-a]
```

[FindInstance](https://reference.wolfram.com/language/ref/FindInstance.html) will find particular instances:

```wolfram
FindInstance[eq,{x,a},2]
(* Output *)
{{x->(92)/(5)-(191 ⅈ)/(10),a->(1)/(5) ⅈ ((191+184 ⅈ)+1780 π-5 ⅈ Log[-(5520)/(70337)-(5730 ⅈ)/(70337)])},{x->43-(411 ⅈ)/(10),a->-(1)/(5) ⅈ ((-411-430 ⅈ)+1200 π+5 ⅈ Log[-(12900)/(353821)-(12330 ⅈ)/(353821)])}}
```

### Possible Issues

If a function is complex, variables are allowed to have complex values:

```wolfram
FindRoot[Zeta[1/2+I t],{t,1}]
(* Output *)
{t->9.920598770752625×10^-17+2.5 ⅈ}
```

If the function is kept real, variables are also taken to be real:

```wolfram
FindRoot[Abs[Zeta[1/2+I t]],{t,1}]
(* Output *)
{t->14.134725141734691}
```

It can be time-consuming to compute functions symbolically:

```wolfram
f[x_]:=Nest[Sin[#+Sin[2#]]&,x,20]
```

```wolfram
FindRoot[f[x],{x,0}]//Timing
(* Output *)
{2.093208999999999,{x->0.}}
```

Restricting the function definition avoids symbolic evaluation:

```wolfram
g[x_?NumericQ]:=Nest[Sin[#+Sin[2#]]&,x,20]
```

```wolfram
FindRoot[g[x],{x,0}]//Timing
(* Output *)
{0.00025400000000086465,{x->0.}}
```

Restricting the search domain too much can be problematic. Find the solution of a system of two nonlinear equations starting from `1/2 in the region from `[0,1]`:

```wolfram
FindRoot[{y==Exp[x],x+y==2},{{x,1/2,0,1},{y,1/2,0,1}}]
(* Output *)
FindRoot
(* Output *)
{x->0.47341835113694847,y->1.}
```

## Tech Notes ▪Numerical Mathematics: Basic Operations ▪Solving Equations ▪Numerical Equation Solving ▪Numerical Root Finding ▪Constrained Optimization ▪Unconstrained Optimization ▪Symbolic Evaluation ▪Implementation notes: Numerical and Related Functions

## Related Guides ▪Equation Solving ▪Inverse Functions ▪Code Compilation ▪Polynomial Equations ▪Symbolic Vectors, Matrices and Arrays

## Related Links [NKS|Online](http://www.wolframscience.com/nks/search/?q=FindRoot) ([A New Kind of Science](http://www.wolframscience.com/nks/))

## History Introduced in 1988 (1.0) | Updated in 2003 (5.0)
