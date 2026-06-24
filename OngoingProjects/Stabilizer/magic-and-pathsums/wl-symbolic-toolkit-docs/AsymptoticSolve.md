# AsymptoticSolve | [SpanFromLeft]

> [AsymptoticSolve](https://reference.wolfram.com/language/ref/AsymptoticSolve.html)[*eqn*,*y->b*,*x->a*]  — computes asymptotic approximations of solutions `*y*[*x*]` of the equation `*eqn*` passing through `{*a*,*b*}`.
> [AsymptoticSolve](https://reference.wolfram.com/language/ref/AsymptoticSolve.html)[*eqn*,{*y*},*x->a*] — computes asymptotic approximations of solutions `*y*[*x*]` of the equation `*eqn*` for `*x*` near `*a*`.
> [AsymptoticSolve](https://reference.wolfram.com/language/ref/AsymptoticSolve.html)[*eqns*,{*y*_1,*y*_2,…}->{*b*_1,*b*_2,…},{*x*_1,*x*_2,…}->{*a*_1,*a*_2,…}] — computes asymptotic approximations of solutions `{*y*_1[*x*_1,*x*_2,…],*y*_2[*x*_1,*x*_2,…],…}` of the system of equations `*eqns*`.
> [AsymptoticSolve](https://reference.wolfram.com/language/ref/AsymptoticSolve.html)[*eqns*,…,{{*x*_1,*x*_2,…},{*a*_1,*a*_2,…},*n*}] — computes the asymptotic approximation to order `*n*`.
> [AsymptoticSolve](https://reference.wolfram.com/language/ref/AsymptoticSolve.html)[…,[Reals](https://reference.wolfram.com/language/ref/Reals.html)] — computes only solutions that are real valued for real argument values.

## Details and Options

Asymptotic approximations are typically used to solve problems for which no exact solution can be found or to get simpler answers for computation, comparison and interpretation.

[AsymptoticSolve](https://reference.wolfram.com/language/ref/AsymptoticSolve.html)[*eqn*,…,*x*->*a*] computes the leading term in an asymptotic expansion for `*eqn*`. Use [SeriesTermGoal](https://reference.wolfram.com/language/ref/SeriesTermGoal.html) to specify more terms.

The asymptotic approximation `*y*_*n*[*x*]` is often given as a sum `*y*_*n*[*x*]==∑_*k*=1^*n*α_*k*φ_*k*[*x*]`, where `{φ_1[*x*],…,φ_*n*[*x*]}` is an asymptotic scale `φ_1[*x*][Succeeds]φ_2[*x*][Succeeds]⋯>φ_*n*[*x*]` as `*x*->*a*`. Then the result satisfies [AsymptoticLess](https://reference.wolfram.com/language/ref/AsymptoticLess.html)[*y*[*x*]-*y*_*n*[*x*],φ_*n*[*x*],*x*->*a*] or `*y*[*x*]-*y*_*n*[*x*]∈*o*[φ_*n*[*x*]]` as `*x*->*a*`.

Common asymptotic scales include:

`$(x-a)^{0}\text{Succeeds}(x-a)^{1}\text{Succeeds}(x-a)^{2}\text{Succeeds}\cdots$` | Taylor scale when `*x*->*a*`
`$(x-a)^{-3}\text{Succeeds}(x-a)^{-2}\text{Succeeds}(x-a)^{-1}\text{Succeeds}\cdots$` | Laurent scale when `*x*->*a*`
`$x^{-1}\text{Succeeds}x^{-2}\text{Succeeds}x^{-3}\text{Succeeds}\cdots$` | Laurent scale when `*x*->±∞`
`$(x-a)^{1/p}\text{Succeeds}(x-a)^{2/p}\text{Succeeds}(x-a)^{3/p}\text{Succeeds}\ldots$` | Puiseux scale when `*x*->*a*`

The scales used to express the asymptotic approximation are automatically inferred from the problem and can often include more exotic scales.

The center coordinates `*a*` and `*b*` can be any finite or infinite real or complex numbers.

The order `*n*` must be a positive integer and specifies order of approximation for the asymptotic solution. It may not be related to polynomial degree.

The system of equations `*eqns*` can be any logical combination of equations.

The following options can be given:

| [Assumptions](https://reference.wolfram.com/language/ref/Assumptions.html) | [$Assumptions](https://reference.wolfram.com/language/ref/$Assumptions.html) | assumptions to make about parameters |
| --- | --- | --- |
| [Direction](https://reference.wolfram.com/language/ref/Direction.html) | [Automatic](https://reference.wolfram.com/language/ref/Automatic.html) | direction in which `*x*` approaches `*a*` |
| [GenerateConditions](https://reference.wolfram.com/language/ref/GenerateConditions.html) | [Automatic](https://reference.wolfram.com/language/ref/Automatic.html) | whether to generate answers that involve conditions on parameters  |
| [Method](https://reference.wolfram.com/language/ref/Method.html) | [Automatic](https://reference.wolfram.com/language/ref/Automatic.html) | method to use |
| [PerformanceGoal](https://reference.wolfram.com/language/ref/PerformanceGoal.html) | [$PerformanceGoal](https://reference.wolfram.com/language/ref/$PerformanceGoal.html) | aspects of performance to optimize |
| [SeriesTermGoal](https://reference.wolfram.com/language/ref/SeriesTermGoal.html) | [Automatic](https://reference.wolfram.com/language/ref/Automatic.html) | number of terms in the approximation |

Possible settings for [Direction](https://reference.wolfram.com/language/ref/Direction.html) include:

[Reals](https://reference.wolfram.com/language/ref/Reals.html) or "TwoSided" | from both real directions
"FromAbove" or  `-1 | from above or larger values
"FromBelow" or `+1 | from below or smaller values
[Complexes](https://reference.wolfram.com/language/ref/Complexes.html) | from all complex directions
[Exp](https://reference.wolfram.com/language/ref/Exp.html)[ⅈ θ] | in the direction $\theta$
{*dir*_1,…,*dir*_*n*} | use direction `*dir*_*i*` for variable `*x*_*i*` independently

[Direction](https://reference.wolfram.com/language/ref/Direction.html)->[Exp](https://reference.wolfram.com/language/ref/Exp.html)[ⅈ θ] at `*x*^(*)` indicates the direction tangent of a curve approaching the limit point `*x*^(*)`.

$$
\text{[Graphics]}
$$

For finite values of `*a*`, the [Automatic](https://reference.wolfram.com/language/ref/Automatic.html) setting means from above.

When domain [Reals](https://reference.wolfram.com/language/ref/Reals.html) is specified, the solutions are real valued when `*x*` approaches `*a*` in the indicated [Direction](https://reference.wolfram.com/language/ref/Direction.html).

Possible settings for [GenerateConditions](https://reference.wolfram.com/language/ref/GenerateConditions.html) include:

[Automatic](https://reference.wolfram.com/language/ref/Automatic.html) | nongeneric conditions only
[True](https://reference.wolfram.com/language/ref/True.html) | all conditions
[False](https://reference.wolfram.com/language/ref/False.html) | no conditions
[None](https://reference.wolfram.com/language/ref/None.html) | return unevaluated if conditions are needed

Possible settings for [PerformanceGoal](https://reference.wolfram.com/language/ref/PerformanceGoal.html) include [$PerformanceGoal](https://reference.wolfram.com/language/ref/$PerformanceGoal.html), `"Quality"` and `"Speed"`. With the `"Quality"` setting, [AsymptoticSolve](https://reference.wolfram.com/language/ref/AsymptoticSolve.html) typically solves more problems or produces simpler results, but it potentially uses more time and memory.

## Examples

### Basic Examples

Find asymptotic approximations of solutions passing through the point `{0,0}`:

```wolfram
AsymptoticSolve[ℯ^y-Cos[x y]+x==0,{y,0},{x,0,7}]
(* Output *)
{{y->-x-(x^2)/(2)-(x^3)/(3)-(3 x^4)/(4)-(6 x^5)/(5)-(13 x^6)/(8)-(141 x^7)/(56)}}
```

Find asymptotic approximations of solutions for `*x*` near `0:

```wolfram
AsymptoticSolve[x y^4-(x+1) y^2+x==1,{y},{x,0,3}]
(* Output *)
{{y->-ⅈ+(3 ⅈ x)/(2)-(27 ⅈ x^2)/(8)+(207 ⅈ x^3)/(16)},{y->ⅈ-(3 ⅈ x)/(2)+(27 ⅈ x^2)/(8)-(207 ⅈ x^3)/(16)},{y->-(1)/(Sqrt[x])-Sqrt[x]+2 x^(3/2)},{y->(1)/(Sqrt[x])+Sqrt[x]-2 x^(3/2)}}
```

Find the leading terms of asymptotic approximations of solutions as $x \to \infty$:

```wolfram
AsymptoticSolve[y^3+E^x y+x Log[x]==0,y,x->Infinity]
(* Output *)
{{y->-ℯ^-x x Log[x]},{y->-ⅈ ℯ^(x/2)},{y->ⅈ ℯ^(x/2)}}
```

Find only the solutions that are real valued when `*x*` approaches `0 from above:

```wolfram
AsymptoticSolve[x y^4-(x+1) y^2+x==1,{y},{x,0,3},Reals]
(* Output *)
{{y->-(1)/(Sqrt[x])-Sqrt[x]+2 x^(3/2)},{y->(1)/(Sqrt[x])+Sqrt[x]-2 x^(3/2)}}
```

Find asymptotic approximations of solutions of a system of equations:

```wolfram
AsymptoticSolve[x^2-u (y+1)==y^2-v (x+1)&&u^2-u (x+1)==v^2+v (y+1),{{u,v},{0,0}},{{x,y},{0,0},3}]
(* Output *)
{{u->(x^2)/(2)-(x^3)/(2)-(y^2)/(2)+(x y^2)/(2),v->-(x^2)/(2)+(x^2 y)/(2)+(y^2)/(2)-(y^3)/(2)}}
```

### Scope

#### One-Dimensional Solutions in 2D

Power series solutions of polynomial equations passing through a specified point:

```wolfram
AsymptoticSolve[y^3-x y==24,{y,3},{x,1,3}]
(* Output *)
{{y->3+(3)/(26) (-1+x)-(3 (-1+x)^2)/(17576)-(309 (-1+x)^3)/(5940688)}}
```

Plot the asymptotic solution and the solution it approximates:

```wolfram
asympt=First[y/.%];
exact[x_]:=Last[y/.NSolve[y^3-x y==24,y,Reals]]
```

```wolfram
Plot[{asympt,exact[x]},{x,-50,50},PlotStyle->{Dashed,DotDashed}]
```

*([Graphics])*

Power series solutions of analytic equations passing through a specified point:

```wolfram
f=Erf[x+y]-Log[1+Sin[x y]];
```

```wolfram
AsymptoticSolve[f==0,{y,0},{x,0,3}]
(* Output *)
{{y->-x-(Sqrt[π] x^2)/(2)-(π x^3)/(4)}}
```

Plot the asymptotic solution and the solution it approximates:

```wolfram
asympt=First[y/.%];
exact[x_]:=y/.FindRoot[f,{y,asympt}]
```

```wolfram
Plot[{asympt,exact[x]},{x,-0.4,0.4},PlotStyle->{Dashed,DotDashed}]
```

*([Graphics])*

Puiseux series solutions of polynomial equations passing through a specified point:

```wolfram
AsymptoticSolve[y^5-x y+x^2==0,{y,0},{x,0,5}]
(* Output *)
{{y->x+x^4},{y->-x^(1/4)-(x)/(4)+(5 x^(7/4))/(32)-(5 x^(5/2))/(32)+(385 x^(13/4))/(2048)-(x^4)/(4)+(23205 x^(19/4))/(65536)},{y->-ⅈ x^(1/4)-(x)/(4)-(5)/(32) ⅈ x^(7/4)+(5 x^(5/2))/(32)+(385 ⅈ x^(13/4))/(2048)-(x^4)/(4)-(23205 ⅈ x^(19/4))/(65536)},{y->ⅈ x^(1/4)-(x)/(4)+(5)/(32) ⅈ x^(7/4)+(5 x^(5/2))/(32)-(385 ⅈ x^(13/4))/(2048)-(x^4)/(4)+(23205 ⅈ x^(19/4))/(65536)},{y->x^(1/4)-(x)/(4)-(5 x^(7/4))/(32)-(5 x^(5/2))/(32)-(385 x^(13/4))/(2048)-(x^4)/(4)-(23205 x^(19/4))/(65536)}}
```

Solutions that are real valued when `*x*` approaches `0 from above:

```wolfram
AsymptoticSolve[y^5-x y+x^2==0,{y,0},{x,0,5},Reals]
(* Output *)
{{y->x+x^4},{y->-x^(1/4)-(x)/(4)+(5 x^(7/4))/(32)-(5 x^(5/2))/(32)+(385 x^(13/4))/(2048)-(x^4)/(4)+(23205 x^(19/4))/(65536)},{y->x^(1/4)-(x)/(4)-(5 x^(7/4))/(32)-(5 x^(5/2))/(32)-(385 x^(13/4))/(2048)-(x^4)/(4)-(23205 x^(19/4))/(65536)}}
```

Solutions that are real valued when `*x*` approaches `0 from below:

```wolfram
AsymptoticSolve[y^5-x y+x^2==0,{y,0},{x,0,5},Reals,Direction->"FromBelow"]
(* Output *)
{{y->x+x^4}}
```

Puiseux series solutions of analytic equations passing through a specified point:

```wolfram
AsymptoticSolve[Sin[x y]-y^2+x==0,{y,0},{x,0,5}]
(* Output *)
{{y->-Sqrt[x]+(x)/(2)-(x^(3/2))/(8)+(x^(5/2))/(128)-(x^(7/2))/(1024)-(x^4)/(12)+(4101 x^(9/2))/(32768)-(x^5)/(12)},{y->Sqrt[x]+(x)/(2)+(x^(3/2))/(8)-(x^(5/2))/(128)+(x^(7/2))/(1024)-(x^4)/(12)-(4101 x^(9/2))/(32768)-(x^5)/(12)}}
```

Asymptotic series solutions passing through a specified point:

```wolfram
AsymptoticSolve[Sin[y]-y^2+Cos[y]/Log[x]==0,{y,0},{x,0,5}]
(* Output *)
{{y->-(51)/(5 Log[x]^5)+(23)/(6 Log[x]^4)-(5)/(3 Log[x]^3)+(1)/(Log[x]^2)-(1)/(Log[x])}}
```

Solutions of polynomial equations near a specified value of the independent variable:

```wolfram
AsymptoticSolve[x y^4+y^2-x^2 y==1,{y},{x,0,3}]
(* Output *)
{{y->-1+(x)/(2)-(3 x^2)/(8)+(17 x^3)/(16)},{y->1-(x)/(2)+(11 x^2)/(8)-(49 x^3)/(16)},{y->(ⅈ)/(Sqrt[x])+(ⅈ Sqrt[x])/(2)-(5)/(8) ⅈ x^(3/2)-(x^2)/(2)},{y->-(ⅈ)/(Sqrt[x])-(ⅈ Sqrt[x])/(2)+(5)/(8) ⅈ x^(3/2)-(x^2)/(2)}}
```

Real asymptotic series solutions at infinity:

```wolfram
f=y^5-y/Sqrt[x]+x Log[x];
```

```wolfram
AsymptoticSolve[f==0,{y},{x,Infinity,3},Reals]
(* Output *)
{{y->-(1)/(25) ((1)/(x))^(21/10) (5 x+(25)/(((1)/(x))^(23/10) ((1)/(Log[x]))^(4/5))-((1)/(x))^(3/10) ((1)/(Log[x]))^(4/5)) ((1)/(Log[x]))^(3/5)}}
```

Plot the asymptotic solution and the solution it approximates:

```wolfram
asympt=First[y/.%];
exact[x_]:=First[y/.NSolve[f==0,y,Reals]]
```

```wolfram
Plot[{asympt,exact[x]},{x,1,5},PlotStyle->{Dashed,DotDashed}]
```

*([Graphics])*

Equations with symbolic parameters:

```wolfram
AsymptoticSolve[Sin[x y]-y Log[a+y]+x==0,{y,0},{x,0,3}]
(* Output *)
{{y->(x)/(Log[a])+(x^2 (-1+a Log[a]))/(a Log[a]^3)+(x^3 (4+Log[a]-6 a Log[a]+2 a^2 Log[a]^2))/(2 a^2 Log[a]^5)}}
```

Conditions on parameters may be generated:

```wolfram
AsymptoticSolve[y^4+a y^3+a x y-x^2+a x^2==0,{y,Root[#^4+a #^3&,1]},{x,0,3}]
(* Output *)
{{y->-a-(x)/(a)+((1+a) x^2)/(a^3)+((-2-5 a) x^3)/(a^5)}}
```

#### One-Dimensional Solutions in nD

Power series solutions of polynomial systems passing through a specified point:

```wolfram
eqns=y^3-t x y==21&&x^2-t^2 y==1;
```

```wolfram
AsymptoticSolve[eqns,{{x,y},{2,3}},{t,1,3}]
(* Output *)
{{x->2+(156)/(97) (-1+t)+(345195 (-1+t)^2)/(912673)-(618244515 (-1+t)^3)/(8587340257),y->3+(42)/(97) (-1+t)+(212997 (-1+t)^2)/(912673)+(247596072 (-1+t)^3)/(8587340257)}}
```

Plot the asymptotic solution and the solution it approximates:

```wolfram
asympt=First[{t,x,y}/.%];
exact={t,x,y}/.Solve[eqns,{x,y},Reals][[2]];
```

```wolfram
ParametricPlot3D[{asympt,exact},{t,0,2},PlotStyle->{Dashed,DotDashed}]
```

*([Graphics3D])*

Power series solutions of analytic systems passing through a specified point:

```wolfram
eqns=Sin[t+ x]+Cos[y]==t+1&&Log[1+x+y]-Sin[t]==0;
```

```wolfram
AsymptoticSolve[eqns,{{x,y},{0,0}},{t,0,5}]
(* Output *)
{{x->(t^2)/(2)+(t^3)/(6)+(t^4)/(24)+(t^5)/(30),y->t-(t^3)/(6)-(t^4)/(6)-(t^5)/(10)}}
```

Plot the asymptotic solution and the solution it approximates:

```wolfram
asympt=First[{t,x,y}/.%];
exact[t_]:={t,x,y}/.FindRoot[eqns,{{x,asympt[[2]]},{y,asympt[[3]]}}]
```

```wolfram
ParametricPlot3D[{asympt,exact[t]},{t,-0.9,0.9},PlotStyle->{Dashed,DotDashed}]
```

*([Graphics3D])*

Puiseux series solutions of polynomial systems passing through a specified point:

```wolfram
AsymptoticSolve[y^2-t  x+t==0&&x^2+t y==0,{{x,y},{0,0}},{t,0,2}]
(* Output *)
{{x->-(-1)^(1/4) t^(3/4)-(1)/(4) ⅈ t^(3/2),y->-ⅈ Sqrt[t]-(1)/(2) (-1)^(3/4) t^(5/4)},{x->(-1)^(1/4) t^(3/4)-(1)/(4) ⅈ t^(3/2),y->-ⅈ Sqrt[t]+(1)/(2) (-1)^(3/4) t^(5/4)},{x->-(-1)^(3/4) t^(3/4)+(1)/(4) ⅈ t^(3/2),y->ⅈ Sqrt[t]-(1)/(2) (-1)^(1/4) t^(5/4)},{x->(-1)^(3/4) t^(3/4)+(1)/(4) ⅈ t^(3/2),y->ⅈ Sqrt[t]+(1)/(2) (-1)^(1/4) t^(5/4)}}
```

None of the solutions are real valued when `*t*` approaches `0 from above:

```wolfram
AsymptoticSolve[y^2-t  x+t==0&&x^2+t y==0,{{x,y},{0,0}},{t,0,2},Reals]
(* Output *)
{}
```

Two of the solutions are real valued when `*t*` approaches `0 from below:

```wolfram
AsymptoticSolve[y^2-t  x+t==0&&x^2+t y==0,{{x,y},{0,0}},{t,0,2},Reals,Direction->"FromBelow"]
(* Output *)
{{x->-(-t)^(3/4)-(1)/(4) (-t)^(3/2),y->Sqrt[-t]+(1)/(2) (-t)^(5/4)},{x->(-t)^(3/4)-(1)/(4) (-t)^(3/2),y->Sqrt[-t]-(1)/(2) (-t)^(5/4)}}
```

Solutions of polynomial systems near a specified value of the independent variable:

```wolfram
AsymptoticSolve[y^3-t x y==0&&x^2-t^2 y==1,{{x,y}},{t,0,3}]
(* Output *)
{{x->-1,y->0},{x->-1+(1)/(2) ⅈ t^(5/2),y->-ⅈ Sqrt[t]-(t^3)/(4)},{x->-1-(1)/(2) ⅈ t^(5/2),y->ⅈ Sqrt[t]-(t^3)/(4)},{x->1,y->0},{x->1-(t^(5/2))/(2),y->-Sqrt[t]+(t^3)/(4)},{x->1+(t^(5/2))/(2),y->Sqrt[t]+(t^3)/(4)}}
```

Equations with symbolic parameters:

```wolfram
AsymptoticSolve[{Sin[t+a x]+Cos[y]==t+1,Log[1+x+y]-Sin[t]==0},{{x,y},{0,0}},{t,0,5}]
(* Output *)
{{x->(t^2)/(2 a)+((-3+4 a) t^3)/(6 a^2)+((15-22 a+8 a^2) t^4)/(24 a^3)+((-105+180 a-100 a^2+29 a^3) t^5)/(120 a^4),y->t+((-1+a) t^2)/(2 a)+((3-4 a) t^3)/(6 a^2)+((-15+22 a-8 a^2-3 a^3) t^4)/(24 a^3)+((105-180 a+100 a^2-29 a^3-8 a^4) t^5)/(120 a^4)}}
```

#### Higher-Dimensional Solutions in nD

Power series solutions of polynomial equations passing through a specified point:

```wolfram
AsymptoticSolve[z^5-2z^3+2z+x^2-y==0,{z,0},{{x,y},{0,0},5}]
(* Output *)
{{z->-(x^2)/(2)+(y)/(2)+(3 x^4 y)/(8)-(3 x^2 y^2)/(8)+(y^3)/(8)+(5 y^5)/(64)}}
```

Plot the asymptotic solution and the solution it approximates:

```wolfram
asympt=First[z/.%];
exact[x_,y_]:=First[z/.NSolve[z^5-2z^3+2z+x^2-y==0,z,Reals]];
```

```wolfram
Plot3D[{asympt,exact[x,y]},{x,-1,1},{y,-1,1},PlotLegends->{"Asymptotic","Exact"},PlotStyle->Opacity[0.7],BoxRatios->1]
```

*([Graphics3D])*

Power series solutions of analytic equations passing through a specified point:

```wolfram
AsymptoticSolve[Exp[z+x-y]-Cos[x y]+x-y==0,{z,0},{{x,y},{0,0},3}]
(* Output *)
{{z->-2 x-(x^2)/(2)-(x^3)/(3)+2 y+x y+x^2 y-(y^2)/(2)-x y^2+(y^3)/(3)}}
```

Plot the asymptotic solution and the solution it approximates:

```wolfram
asympt=First[z/.%];
exact[x_,y_]:=z/.FindRoot[Exp[z+x-y]-Cos[x y]+x-y==0,{z,asympt}]
```

```wolfram
Plot3D[{asympt,exact[x,y]},{x,-1,1},{y,-1,1},PlotLegends->{"Asymptotic","Exact"},PlotStyle->Opacity[0.7],BoxRatios->1]
(* Output *)
![image](img/image_001.png)
```

Power series solutions of polynomial systems passing through a specified point:

```wolfram
AsymptoticSolve[x^2-u (y+1)==y^2-v (x+1)+w&&u^2-u (x+1)==v^2+v (y+1)-w&&x+y==u+v+w,{{u,v,w},{0,0,0}},{{x,y},{0,0},3}]
(* Output *)
{{u->(5 x^2)/(8)-(5 x^3)/(16)-(x y)/(4)+(9 x^2 y)/(16)-(7 y^2)/(8)+(13 x y^2)/(16)-(y^3)/(16),v->(x)/(2)-(3 x^2)/(4)+(3 x^3)/(8)+(y)/(2)-(x y)/(4)+(7 x^2 y)/(16)+(y^2)/(2)-(3 x y^2)/(8)-(7 y^3)/(16),w->(x)/(2)+(x^2)/(8)-(x^3)/(16)+(y)/(2)+(x y)/(2)-x^2 y+(3 y^2)/(8)-(7 x y^2)/(16)+(y^3)/(2)}}
```

Power series solutions of analytic systems passing through a specified point:

```wolfram
AsymptoticSolve[{Sin[u+x]+Cos[y]==v+1,Log[1+x+y]-u-Sin[u v]==0},{{u,v},{0,0}},{{x,y},{0,0},3}]
(* Output *)
{{u->x-(5 x^2)/(2)+(47 x^3)/(6)+y-4 x y+18 x^2 y-(3 y^2)/(2)+14 x y^2+(23 y^3)/(6),v->2 x-(5 x^2)/(2)+(13 x^3)/(2)+y-4 x y+16 x^2 y-2 y^2+13 x y^2+(11 y^3)/(3)}}
```

Power series solutions of polynomial systems near specified values of independent variables:

```wolfram
AsymptoticSolve[x^2-u==y^2-v(x+1)-u^2&&u^2-u(x+1)==v^2+v,{{u,v}},{{x,y},{0,0},2}]
(* Output *)
{{u->-1-(7 x)/(6)+(83 x^2)/(216)-(y^2)/(2),v->-2-(3 x)/(2)+(7 x^2)/(24)-(y^2)/(2)},{u->(x^2)/(2)-(y^2)/(2),v->-(x^2)/(2)+(y^2)/(2)},{u->1+(x)/(2)-(x^2)/(8)+(y^2)/(2),v->-(x)/(2)-(5 x^2)/(8)+(y^2)/(2)},{u->2+(2 x)/(3)-(41 x^2)/(54)+(y^2)/(2),v->-2+(5 x^2)/(6)-(y^2)/(2)}}
```

### Options

#### Assumptions

Specify conditions on parameters using [Assumptions](https://reference.wolfram.com/language/ref/Assumptions.html):

```wolfram
AsymptoticSolve[y^3-a y+x^2y+x==0,{y},{x,0,3},Reals,Assumptions->a>0]
(* Output *)
{{y->(x)/(a)+((1+a^2) x^3)/(a^4)},{y->-Sqrt[a]-(x)/(2 a)+((3+4 a^2) x^2)/(8 a^(5/2))+((-1-a^2) x^3)/(2 a^4)},{y->Sqrt[a]-(x)/(2 a)+((-3-4 a^2) x^2)/(8 a^(5/2))+((-1-a^2) x^3)/(2 a^4)}}
```

Different assumptions can produce different results:

```wolfram
AsymptoticSolve[y^3-a y+x^2y+x==0,{y},{x,0,3},Reals,Assumptions->a<0]
(* Output *)
{{y->(x)/(a)+((1+a^2) x^3)/(a^4)}}
```

#### Direction

By default, [AsymptoticSolve](https://reference.wolfram.com/language/ref/AsymptoticSolve.html) gives solutions valid when `*x*` approaches `0 from above:

```wolfram
AsymptoticSolve[y^3+x^2y+x==0,{y},{x,0,3},Reals]
(* Output *)
{{y->-x^(1/3)+(x^(5/3))/(3)}}
```

This finds the solutions valid when `*x*` approaches `0 from below:

```wolfram
AsymptoticSolve[y^3+x^2y+x==0,{y},{x,0,3},Reals,Direction->"FromBelow"]
(* Output *)
{{y->(-x)^(1/3)-(1)/(3) (-x)^(5/3)}}
```

Complex solutions may also depend on the direction:

```wolfram
AsymptoticSolve[y^2+E^(-2/x) y+E^(-1/x)==0,{y,0},{x,0,5},Direction->"FromAbove"]
(* Output *)
{{y->-(1)/(2) ℯ^(-2/x)-ⅈ Sqrt[ℯ^(-1/x)]+(1)/(8) ⅈ (ℯ^(-1/x))^(7/2)},{y->-(1)/(2) ℯ^(-2/x)+ⅈ Sqrt[ℯ^(-1/x)]-(1)/(8) ⅈ (ℯ^(-1/x))^(7/2)}}
```

```wolfram
AsymptoticSolve[y^2+E^(-2/x) y+E^(-1/x)==0,{y,0},{x,0,5},Direction->"FromBelow"]
(* Output *)
{{y->-ℯ^((1)/(x))-ℯ^(4/x)}}
```

This gives solutions that are real when `*x*` approaches `0 from a complex direction:

```wolfram
AsymptoticSolve[x y^7+2 x^4 y^6+3 x^4 y^3-2 x^7 y^2-x^10==0,{y,0},{x,0,3},Reals,Direction->E^(2 I Pi/3)]
(* Output *)
{{y->(ℯ^((2 ⅈ π)/(3)) x^2)/(3^(1/3))+(2 x^3)/(9)},{y->-(2 x^3)/(3)-3^(1/4) (-ℯ^(-(2 ⅈ π)/(3)) x)^(3/4)},{y->-(2 x^3)/(3)+3^(1/4) (-ℯ^(-(2 ⅈ π)/(3)) x)^(3/4)}}
```

#### GenerateConditions

By default, [AsymptoticSolve](https://reference.wolfram.com/language/ref/AsymptoticSolve.html) gives conditions it assumed to obtain the result:

```wolfram
AsymptoticSolve[y^3+a y^2-x==0,{y},{x,0,3},Reals]
(* Output *)
{{y->-a+(x)/(a^2)+(2 x^2)/(a^5)+(7 x^3)/(a^8)}}
```

This gives the result without the assumed conditions:

```wolfram
AsymptoticSolve[y^3+a y^2-x==0,{y},{x,0,3},Reals,GenerateConditions->False]
(* Output *)
{{y->-Sqrt[(1)/(a)] Sqrt[x]-(x)/(2 a^2)-(5)/(8) ((1)/(a))^(7/2) x^(3/2)-(x^2)/(a^5)-(231)/(128) ((1)/(a))^(13/2) x^(5/2)-(7 x^3)/(2 a^8)},{y->Sqrt[(1)/(a)] Sqrt[x]-(x)/(2 a^2)+(5)/(8) ((1)/(a))^(7/2) x^(3/2)-(x^2)/(a^5)+(231)/(128) ((1)/(a))^(13/2) x^(5/2)-(7 x^3)/(2 a^8)},{y->-a+(x)/(a^2)+(2 x^2)/(a^5)+(7 x^3)/(a^8)}}
```

By default, assumed conditions that are generically true are not reported:

```wolfram
AsymptoticSolve[y^3+a y^2-x==0,{y},{x,0,2}]
(* Output *)
{{y->-Sqrt[(1)/(a)] Sqrt[x]-(x)/(2 a^2)-(5)/(8) ((1)/(a))^(7/2) x^(3/2)-(x^2)/(a^5)},{y->Sqrt[(1)/(a)] Sqrt[x]-(x)/(2 a^2)+(5)/(8) ((1)/(a))^(7/2) x^(3/2)-(x^2)/(a^5)},{y->-a+(x)/(a^2)+(2 x^2)/(a^5)}}
```

With [GenerateConditions](https://reference.wolfram.com/language/ref/GenerateConditions.html)->[True](https://reference.wolfram.com/language/ref/True.html), all conditions are reported:

```wolfram
AsymptoticSolve[y^3+a y^2-x==0,{y},{x,0,2},GenerateConditions->True]
(* Output *)
{{y->-Sqrt[(1)/(a)] Sqrt[x]-(x)/(2 a^2)-(5)/(8) ((1)/(a))^(7/2) x^(3/2)-(x^2)/(a^5)},{y->Sqrt[(1)/(a)] Sqrt[x]-(x)/(2 a^2)+(5)/(8) ((1)/(a))^(7/2) x^(3/2)-(x^2)/(a^5)},{y->-a+(x)/(a^2)+(2 x^2)/(a^5)}}
```

With [GenerateConditions](https://reference.wolfram.com/language/ref/GenerateConditions.html)->[None](https://reference.wolfram.com/language/ref/None.html), [AsymptoticSolve](https://reference.wolfram.com/language/ref/AsymptoticSolve.html) returns only generically valid results:

```wolfram
AsymptoticSolve[y^3+a y^2-x==0,{y},{x,0,2},GenerateConditions->None]
(* Output *)
{{y->-Sqrt[(1)/(a)] Sqrt[x]-(x)/(2 a^2)-(5)/(8) ((1)/(a))^(7/2) x^(3/2)-(x^2)/(a^5)},{y->Sqrt[(1)/(a)] Sqrt[x]-(x)/(2 a^2)+(5)/(8) ((1)/(a))^(7/2) x^(3/2)-(x^2)/(a^5)},{y->-a+(x)/(a^2)+(2 x^2)/(a^5)}}
```

If nongeneric conditions are needed, [AsymptoticSolve](https://reference.wolfram.com/language/ref/AsymptoticSolve.html) returns unevaluated:

```wolfram
AsymptoticSolve[y^3+a y^2-x==0,{y},{x,0,2},Reals,GenerateConditions->None]
(* Output *)
AsymptoticSolve[-x+a y^2+y^3==0,{y},{x,0,2},Reals,GenerateConditions->None]
```

#### Method

Return a series whenever the result is a power series or a Puiseux series:

```wolfram
AsymptoticSolve[Tan[x+y]^3-y^2+x==0,{y,0},{x,0,3},Method->{"SeriesOutput"->Automatic}]
(* Output *)
{{y->-Sqrt[x]+(x)/(2)-(17 x^(3/2))/(8)+6 x^2-(2399 x^(5/2))/(128)+(958 x^3)/(15)+O[x]^(7/2)},{y->Sqrt[x]+(x)/(2)+(17 x^(3/2))/(8)+6 x^2+(2399 x^(5/2))/(128)+(958 x^3)/(15)+O[x]^(7/2)}}
```

Check that the series solutions satisfy the equation:

```wolfram
Tan[x+y]^3-y^2+x/.%
(* Output *)
{O[x]^4,O[x]^4}
```

#### SeriesTermGoal

By default, [AsymptoticSolve](https://reference.wolfram.com/language/ref/AsymptoticSolve.html)[*eqn*,…,*x*->*a*] computes the leading terms of the solutions:

```wolfram
AsymptoticSolve[Cos[x y]-Exp[y^2-2x^2]==0,y->0,x->0]
(* Output *)
{{y->-Sqrt[2] x},{y->Sqrt[2] x}}
```

Use [SeriesTermGoal](https://reference.wolfram.com/language/ref/SeriesTermGoal.html) to obtain more terms:

```wolfram
AsymptoticSolve[Cos[x y]-Exp[y^2-2x^2]==0,y->0,x->0,SeriesTermGoal->5]
(* Output *)
{{y->-Sqrt[2] x+(x^3)/(2 Sqrt[2])-(3 x^5)/(16 Sqrt[2])},{y->Sqrt[2] x-(x^3)/(2 Sqrt[2])+(3 x^5)/(16 Sqrt[2])}}
```

### Applications

#### Implicit Functions

The equation $y^{2}=x$ implicitly defines two different functions $y(x)$ near each $x>0$. Compute third-order asymptotic approximations for these two functions near $x=4$:

```wolfram
AsymptoticSolve[y^2==x,y,{x,4,3}]
(* Output *)
{{y->-2+(4-x)/(4)+(1)/(64) (-4+x)^2-(1)/(512) (-4+x)^3},{y->2+(1)/(4) (-4+x)-(1)/(64) (-4+x)^2+(1)/(512) (-4+x)^3}}
```

Define functions based on these expansions:

```wolfram
{yA1[x_],yA2[x_]}= y/.%;
```

At the point $x=4$, these two functions have different values:

```wolfram
{yA1[4],yA2[4]}
(* Output *)
{-2,2}
```

However, both exactly satisfy the equation $y^{2}=x$ at $x=4$:

```wolfram
{yA1[4]^2==4,yA2[4]^2==4}
(* Output *)
{True,True}
```

Visualize the equation and the approximations to its two branches:

```wolfram
ContourPlot[{y^2==x,y==yA1[x],y==yA2[x]},{x,-1,6},{y,-3,3},Epilog->{PointSize[Large],Point[{{4,-2},{4,2}}]}]
```

*([Graphics])*

In this case, it is easy to solve exactly for the two implicitly defined functions:

```wolfram
{y1[x_],y2[x_]}=y/.Solve[y^2==x,y]
(* Output *)
{-Sqrt[x],Sqrt[x]}
```

The two expressions returned by [AsymptoticSolve](https://reference.wolfram.com/language/ref/AsymptoticSolve.html) are the series of the exact solutions:

```wolfram
Series[yA1[x]-y1[x],{x,4,3}]//Normal
(* Output *)
0
```

```wolfram
Series[yA2[x]-y2[x],{x,4,3}]//Normal
(* Output *)
0
```

Compute the second-order asymptotic approximations to the unit circle $x^{2}+y^{2}=1$ at $x=0$:

```wolfram
{yR1[x_],yR2[x_]}=y/.AsymptoticSolve[y^2+x^2==1,y,{x,0,2}]
(* Output *)
{-1+(x^2)/(2),1-(x^2)/(2)}
```

The approximations closely track the circle of both larger and smaller values of $x$ at these regular points:

```wolfram
ContourPlot[{y^2+x^2==1,y==yR1[x],y==yR2[x]},{x,-3/2,3/2},{y,-3/2,3/2},Epilog->{PointSize[Large],Point[{{0,1},{0,-1}}]}]
```

*([Graphics])*

At the singular point $x=1$, the approximation uses fractional powers:

```wolfram
{yS1[x_],yS2[x_]}=y/.AsymptoticSolve[y^2==1-x^2,y,{x,1,2}]
(* Output *)
{-ⅈ Sqrt[2] Sqrt[-1+x]-(ⅈ (-1+x)^(3/2))/(2 Sqrt[2]),ⅈ Sqrt[2] Sqrt[-1+x]+(ⅈ (-1+x)^(3/2))/(2 Sqrt[2])}
```

Visually, the approximations seem to only be defined for values of $x<1$:

```wolfram
ContourPlot[{y^2+x^2==1,y==yS1[x],y==yS2[x]},{x,-3/2,3/2},{y,-3/2,3/2},Epilog->{PointSize[Large],Point[{1,0}]}]
```

*([Graphics])*

This is because at $x=1$ the functions switch from being real to purely imaginary:

```wolfram
Plot[{Im[yS1[x]],Im[yS2[x]]},{x,0,2},PlotTheme->{"Detailed","DashedLines"}]
```

*([Graphics])*

Trying to find solutions over the reals will therefore fail:

```wolfram
AsymptoticSolve[y^2+x^2==1,y,{x,1,2},Reals]
(* Output *)
{}
```

However, it is possible to find purely real expressions if restricting to smaller values of $x$:

```wolfram
AsymptoticSolve[y^2==1-x^2,y,{x,1,2},Reals,Direction->"FromBelow"]
(* Output *)
{{y->-Sqrt[2] Sqrt[1-x]+((1-x)^(3/2))/(2 Sqrt[2])},{y->Sqrt[2] Sqrt[1-x]-((1-x)^(3/2))/(2 Sqrt[2])}}
```

The curve $sin(y)=x$ crosses the line $x=0$ infinitely many times. On any section that passes the vertical line test--any vertical line intersects the curve only once; no vertical line intersects the section more than once--a function is implicitly defined:

```wolfram
ContourPlot[{Sin[y]==x,x==0},{x,-1,1},{y,-3π,3π},PlotTheme->"DashedLines"]
```

*([Graphics])*

Compute an approximation for the section that goes through the origin:

```wolfram
AsymptoticSolve[Sin[y]==x,{y,0},{x,0,3}]
(* Output *)
{{y->x+(x^3)/(6)}}
```

Note that this matches the Taylor series of $sin^{-1}(x)$, the inverse function of $sin(y)$:

```wolfram
Series[ArcSin[x],{x,0,3}]
(* Output *)
x+(x^3)/(6)+O[x]^4
```

Compute an approximation for the section that goes through the point $(0,\pi)$:

```wolfram
AsymptoticSolve[Sin[y]==x,{y,Pi},{x,0,3}]
(* Output *)
{{y->π-x-(x^3)/(6)}}
```

Visualize the curve and the two approximations:

```wolfram
ContourPlot[{Sin[y]==x,x==0,y==x+(x^3)/(6),y==π-x-(x^3)/(6)},{x,-1.5,1.5},{y,-π,2π},PlotTheme -> DashedLines, PlotLegends -> Expressions, Epilog -> RGBColor[0, 0, 1]PointSize[Large]Point[000Pi]]
```

*([Graphics])*

#### Perturbed Equations

Find solutions of a perturbed polynomial equation:

```wolfram
f=x^3-(6+ε)x^2+(11+ε)x-6+ε;
```

```wolfram
AsymptoticSolve[f==0,{x},{ε,0,3}]
(* Output *)
{{x->1-(ε)/(2)+(ε^2)/(8)+(ε^3)/(16)},{x->2-ε+3 ε^2-11 ε^3},{x->3+(5 ε)/(2)-(25 ε^2)/(8)+(175 ε^3)/(16)}}
```

Plot the asymptotic solutions and the solutions they approximate:

```wolfram
asympt=x/.%;
exact=x/.Solve[f==0,x];
```

```wolfram
Plot[{asympt,exact},{ε,0,0.5},PlotStyle->{Dashed,DotDashed}]
```

*([Graphics])*

Investigate the behavior of solutions of an analytic equation under a small perturbation:

```wolfram
f=Cos[x]-1+ε E^x;
```

```wolfram
AsymptoticSolve[f==0,{x,0},{ε,0,3}]
(* Output *)
{{x->-Sqrt[2] Sqrt[ε]+ε-(5 ε^(3/2))/(3 Sqrt[2])+(5 ε^2)/(3)-(221 ε^(5/2))/(60 Sqrt[2])+(13 ε^3)/(3)},{x->Sqrt[2] Sqrt[ε]+ε+(5 ε^(3/2))/(3 Sqrt[2])+(5 ε^2)/(3)+(221 ε^(5/2))/(60 Sqrt[2])+(13 ε^3)/(3)}}
```

Plot the asymptotic solutions and the solutions they approximate:

```wolfram
asympt=x/.%;
exact[s_]:=x/.FindRoot[f,{x,s}]
```

```wolfram
Show[Plot[{#,exact[#]},{ε,0,0.2},PlotStyle->{Dashed,DotDashed}]&/@asympt,PlotRange->All]
```

*([Graphics])*

#### Series Solutions of Equations

Find a series solution of an equation at a nonsingular point:

```wolfram
f=Cos[x+y]-x+y-1;
```

The derivative of $f$ with respect to $y$ does not vanish at $0$:

```wolfram
D[f,y]/.{x->0,y->0}
(* Output *)
1
```

Find the series solution up to order five:

```wolfram
AsymptoticSolve[f==0,{y,0},{x,0,5}]
(* Output *)
{{y->x+2 x^2+4 x^3+(28 x^4)/(3)+24 x^5}}
```

The result satisfies the equation:

```wolfram
Series[f/.%[[1]],{x,0,5}]
(* Output *)
O[x]^6
```

Find a multivariate series solution of a system of equations at a nonsingular point:

```wolfram
f=Sin[x+y+u+v]-u x+v y;
g=Exp[x-u]+Exp[y-v]-2 Cos[x-v]+y-u;
```

The Jacobian of $\{f,g \}$ with respect to $\{x,y \}$ does not vanish at $0$:

```wolfram
Det[D[{f,g},{{x,y}}]]/.{x->0,y->0,u->0,v->0}
(* Output *)
1
```

Find the series solution up to order three:

```wolfram
AsymptoticSolve[f==0&&g==0,{{x,y},{0,0}},{{u,v},{0,0},3}]
(* Output *)
{{x->-4 u+25 u^2-(1135 u^3)/(3)-3 v+38 u v-829 u^2 v+17 v^2-639 u v^2-(517 v^3)/(3),y->3 u-29 u^2+(1210 u^3)/(3)+2 v-44 u v+896 u^2 v-19 v^2+700 u v^2+(574 v^3)/(3)}}
```

The result satisfies the equations:

```wolfram
Series[{f,g}/.%[[1]]/.{u->t u,v->t v},{t,0,3}]
(* Output *)
{O[t]^4,O[t]^4}
```

#### Asymptotic Approximations of Curves

Find Puiseux series solutions in a neighborhood of a singular point of an algebraic plane curve:

```wolfram
f=y^6-5x y^4+3x^2 y^4+10x^3 y^2+3x^4 y^2-x^5+x^6;
```

```wolfram
ls=AsymptoticSolve[f==0,{y,0},{x,0,3},Reals,Direction->"FromBelow"]
(* Output *)
{{y->Sqrt[(1)/(5) (5-2 Sqrt[5])] x+(8 (-5+2 Sqrt[5]) x^2)/(25 Sqrt[5-2 Sqrt[5]])-(16 (-5+2 Sqrt[5]) (-5+3 Sqrt[5]) x^3)/(625 Sqrt[5-2 Sqrt[5]])},{y->-Sqrt[(1)/(5) (5-2 Sqrt[5])] x-(8 (-5+2 Sqrt[5]) x^2)/(25 Sqrt[5-2 Sqrt[5]])+(16 (-5+2 Sqrt[5]) (-5+3 Sqrt[5]) x^3)/(625 Sqrt[5-2 Sqrt[5]])},{y->Sqrt[(1)/(5) (5+2 Sqrt[5])] x+(8)/(25) Sqrt[5+2 Sqrt[5]] x^2+(16)/(625) Sqrt[5+2 Sqrt[5]] (5+3 Sqrt[5]) x^3},{y->-Sqrt[(1)/(5) (5+2 Sqrt[5])] x-(8)/(25) Sqrt[5+2 Sqrt[5]] x^2-(16)/(625) Sqrt[5+2 Sqrt[5]] (5+3 Sqrt[5]) x^3}}
```

```wolfram
rs=AsymptoticSolve[f==0,{y,0},{x,0,3},Reals,Direction->"FromAbove"]
(* Output *)
{{y->-Sqrt[(1)/(5) (5-2 Sqrt[5])] x-(8 (-5+2 Sqrt[5]) x^2)/(25 Sqrt[5-2 Sqrt[5]])+(16 (-5+2 Sqrt[5]) (-5+3 Sqrt[5]) x^3)/(625 Sqrt[5-2 Sqrt[5]])},{y->Sqrt[(1)/(5) (5-2 Sqrt[5])] x+(8 (-5+2 Sqrt[5]) x^2)/(25 Sqrt[5-2 Sqrt[5]])-(16 (-5+2 Sqrt[5]) (-5+3 Sqrt[5]) x^3)/(625 Sqrt[5-2 Sqrt[5]])},{y->-Sqrt[(1)/(5) (5+2 Sqrt[5])] x-(8)/(25) Sqrt[5+2 Sqrt[5]] x^2-(16)/(625) Sqrt[5+2 Sqrt[5]] (5+3 Sqrt[5]) x^3},{y->Sqrt[(1)/(5) (5+2 Sqrt[5])] x+(8)/(25) Sqrt[5+2 Sqrt[5]] x^2+(16)/(625) Sqrt[5+2 Sqrt[5]] (5+3 Sqrt[5]) x^3},{y->-Sqrt[5] Sqrt[x]+(1)/(2) Sqrt[5] x^(3/2)+(381 x^(5/2))/(200 Sqrt[5])},{y->Sqrt[5] Sqrt[x]-(1)/(2) Sqrt[5] x^(3/2)-(381 x^(5/2))/(200 Sqrt[5])}}
```

Plot the asymptotic solutions and the curve they approximate near 0:

```wolfram
lp=ParametricPlot[Evaluate[{x,y}/.ls],{x,-1,0},PlotStyle->{{Red,Dashed}}];
```

```wolfram
rp=ParametricPlot[Evaluate[{x,y}/.rs],{x,0,1},PlotStyle->{{Red,Dashed}}];
```

```wolfram
exact=ContourPlot[f==0,{x,-1,1},{y,-1,1},PlotPoints->100];
```

```wolfram
Show[{exact,lp,rp}]
```

*([Graphics])*

Find Puiseux series solutions in a neighborhood of a singular point of an algebraic space curve:

```wolfram
eqns=z^2-y^2+x y z+x^2 z+x==0&&x^2 z^2+x y^2+y z+x y+x^2==0;
```

```wolfram
ls=AsymptoticSolve[eqns,{{y, z},{0,0}},{x,0,3},Reals,Direction->"FromBelow"]
(* Output *)
{{y->(-x)^(3/2)+2 (-x)^(5/2)-x^2+(5 x^3)/(2),z->-Sqrt[-x]-(x^2)/(2)-(x^3)/(2)},{y->-(-x)^(3/2)-2 (-x)^(5/2)-x^2+(5 x^3)/(2),z->Sqrt[-x]-(x^2)/(2)-(x^3)/(2)}}
```

```wolfram
rs=AsymptoticSolve[eqns,{{y, z},{0,0}},{x,0,3},Reals,Direction->"FromAbove"]
(* Output *)
{{y->-Sqrt[x]-(x^(3/2))/(2)+(3 x^2)/(2)-(3 x^(5/2))/(8)-2 x^3,z->-x+2 x^(3/2)},{y->Sqrt[x]+(x^(3/2))/(2)+(3 x^2)/(2)+(3 x^(5/2))/(8)-2 x^3,z->-x-2 x^(3/2)}}
```

Plot the asymptotic solutions:

```wolfram
lasym={x,y,z}/.ls;
rasym={x,y,z}/.rs;
lp=ParametricPlot3D[lasym,{x,-0.3,0},PlotStyle->{{Red,Dashed}}];
rp=ParametricPlot3D[rasym,{x,0,0.3},PlotStyle->{{Red,Dashed}}];
```

Find numeric solutions using asymptotic solution values as starting points:

```wolfram
nsol[s_]:={x,y,z}/.FindRoot[eqns,{{y,s[[2]]},{z,s[[3]]}}]
```

```wolfram
lpn=ParametricPlot3D[nsol/@lasym,{x,-0.3,0},PlotStyle->Blue];
rpn=ParametricPlot3D[nsol/@rasym,{x,0,0.3},PlotStyle->Blue];
```

Compare the numeric solutions and the asymptotic solutions:

```wolfram
Show[{lpn,rpn,lp,rp},PlotRange->All,BoxRatios->1]
```

*([Graphics3D])*

Show the curve as an intersection of two surfaces:

```wolfram
surf=ContourPlot3D[Evaluate[List@@eqns],{x,-0.3,0.3},{y,-0.6,0.6},{z,-0.6,0.6},Mesh->None,ContourStyle->{Directive[Orange,Opacity[0.3]],Directive[Yellow,Opacity[0.3]]}];
```

```wolfram
Show[{surf,lpn,rpn,lp,rp},PlotRange->All,BoxRatios->1]
```

*([Graphics3D])*

Approximate Fermat's spiral near 0:

```wolfram
AsymptoticSolve[x==r Cos[r^2]&&y==r Sin[r^2],{{y,r},{0,0}},{x,0,12}]
(* Output *)
{{y->x^3+(4 x^7)/(3)+(19 x^11)/(5),r->x+(x^5)/(2)+(29 x^9)/(24)}}
```

Compare plots:

```wolfram
asympt=First[y/.%];
ap=Plot[asympt,{x,0,1},PlotStyle->{Red,DotDashed}];
spiral=ParametricPlot[{r Cos[r^2],r Sin[r^2]},{r,0,5},PlotStyle->Dashed];
```

```wolfram
Show[{ap,spiral},PlotRange->All,AspectRatio->Automatic]
```

*([Graphics])*

#### Asymptotic Solutions of Physics Problems

Solve Kepler's equation for the eccentric anomaly $E$ in terms of the mean anomaly $M$:

```wolfram
AsymptoticSolve[M==Ε-e Sin[Ε],{Ε,0},{M,0,7}]
(* Output *)
{{Ε->(M)/(1-e)-(e M^3)/(6 (-1+e)^4)-(e (1+9 e) M^5)/(120 (-1+e)^7)+((-e-54 e^2-225 e^3) M^7)/(5040 (-1+e)^10)}}
```

Compare with the exact solution for eccentricity $e=\frac{1}{2}$:

```wolfram
asympt=First[Ε/.%]/.e->1/2;
exact[M_]:=First[Ε/.Solve[M==Ε-1/2Sin[Ε],Ε,Reals]//Quiet]
```

```wolfram
Plot[{asympt,exact[M]},{M,0,0.7},PlotStyle->{Dashed,DotDashed}]
```

*([Graphics])*

Study the energy levels of a particle of mass $m$ in a one-dimensional box of width $L$ and depth $V$. Solutions $\psi_{1}$, $\psi_{2}$ and $\psi_{3}$ of the time-independent Schrödinger equation to the left of the box, inside the box, and to the right of the box are given by:

```wolfram
α=(Sqrt[2 m (V-Ε)])/(ℏ);k=(Sqrt[2 m Ε])/(ℏ);
ψ_1=g Exp[α x];ψ_2=a Sin[k x]+b Cos[k x];ψ_3=h Exp[-α x];
```

The solution must be continuously differentiable on the boundary of the box:

```wolfram
bconds={(ψ_1/.x->-L/2)==(ψ_2/.x->-L/2),(ψ_2/.x->L/2)==(ψ_3/.x->L/2),(∂_xψ_1/.x->-L/2)==(∂_xψ_2/.x->-L/2),(∂_xψ_2/.x->L/2)==(∂_xψ_3/.x->L/2)}
(* Output *)
{ℯ^(-(L Sqrt[m (V-Ε)])/(Sqrt[2] ℏ)) g==b Cos[(L Sqrt[m Ε])/(Sqrt[2] ℏ)]-a Sin[(L Sqrt[m Ε])/(Sqrt[2] ℏ)],b Cos[(L Sqrt[m Ε])/(Sqrt[2] ℏ)]+a Sin[(L Sqrt[m Ε])/(Sqrt[2] ℏ)]==ℯ^(-(L Sqrt[m (V-Ε)])/(Sqrt[2] ℏ)) h,(Sqrt[2] ℯ^(-(L Sqrt[m (V-Ε)])/(Sqrt[2] ℏ)) g Sqrt[m (V-Ε)])/(ℏ)==(Sqrt[2] a Sqrt[m Ε] Cos[(L Sqrt[m Ε])/(Sqrt[2] ℏ)])/(ℏ)+(Sqrt[2] b Sqrt[m Ε] Sin[(L Sqrt[m Ε])/(Sqrt[2] ℏ)])/(ℏ),(Sqrt[2] a Sqrt[m Ε] Cos[(L Sqrt[m Ε])/(Sqrt[2] ℏ)])/(ℏ)-(Sqrt[2] b Sqrt[m Ε] Sin[(L Sqrt[m Ε])/(Sqrt[2] ℏ)])/(ℏ)==-(Sqrt[2] ℯ^(-(L Sqrt[m (V-Ε)])/(Sqrt[2] ℏ)) h Sqrt[m (V-Ε)])/(ℏ)}
```

The homogenous linear equations admit nonzero solutions if their coefficient matrix is singular:

```wolfram
eqn=Det[CoefficientArrays[bconds,{a,b,g,h}][[2]]]==0
(* Output *)
(4 ℯ^(-(Sqrt[2] L Sqrt[m (V-Ε)])/(ℏ)) Sqrt[m (V-Ε)] Sqrt[m Ε] Cos[(L Sqrt[m Ε])/(Sqrt[2] ℏ)]^2)/(ℏ^2)+(4 ℯ^(-(Sqrt[2] L Sqrt[m (V-Ε)])/(ℏ)) m V Cos[(L Sqrt[m Ε])/(Sqrt[2] ℏ)] Sin[(L Sqrt[m Ε])/(Sqrt[2] ℏ)])/(ℏ^2)-(8 ℯ^(-(Sqrt[2] L Sqrt[m (V-Ε)])/(ℏ)) m Ε Cos[(L Sqrt[m Ε])/(Sqrt[2] ℏ)] Sin[(L Sqrt[m Ε])/(Sqrt[2] ℏ)])/(ℏ^2)-(4 ℯ^(-(Sqrt[2] L Sqrt[m (V-Ε)])/(ℏ)) Sqrt[m (V-Ε)] Sqrt[m Ε] Sin[(L Sqrt[m Ε])/(Sqrt[2] ℏ)]^2)/(ℏ^2)==0
```

Assume that $m$ and $L$ are 1 and the units are chosen so that $\hbar=1$:

```wolfram
m=1;L=1; ℏ=1;
```

Find the possible energy levels for $V=1$:

```wolfram
e1=Ε/.Solve[(eqn/.V->1)&&Ε>0,Ε,Reals]
(* Output *)
{Root}
```

Compute the asymptotic solution near $V=1$:

```wolfram
AsymptoticSolve[eqn,{Ε,N[e1[[1]],20]},{V,1,5}]
(* Output *)
{{Ε->0.69207840729638742868601466337493310476+0.49704765254350857004019065453328707396 (-1+V)-0.1272320183888753406318151940606459238 (-1+V)^2+0.04167687353587762619072553198505408427 (-1+V)^3-0.01541795805439759373304626420886660118 (-1+V)^4+0.00613293609311328310773052385856846242 (-1+V)^5}}
```

Compare the asymptotic solution and the minimum exact solution:

```wolfram
asympt=First[Ε/.%];
exact=Table[{V,Min[Ε/.Solve[eqn&&Ε>=1/1000,Ε,Reals]]},{V,1/4,4,1/4}];
```

```wolfram
Show[{Plot[asympt,{V,0.25,4}],ListPlot[exact,PlotStyle->Red]}]
```

*([Graphics])*

## Related Guides ▪Asymptotics ▪Equation Solving ▪Polynomial Algebra

## History Introduced in 2019 (12.0) | Updated in 2020 (12.1)
