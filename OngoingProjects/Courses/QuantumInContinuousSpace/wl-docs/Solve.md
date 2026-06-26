# Solve | [SpanFromLeft]

> [Solve](https://reference.wolfram.com/language/ref/Solve.html)[*expr*,*vars*] — attempts to solve the system `*expr*` of equations or inequalities for the variables `*vars*`.
> [Solve](https://reference.wolfram.com/language/ref/Solve.html)[*expr*,*vars*,*dom*] — solves over the domain `*dom*`. Common choices of `*dom*` are [Reals](https://reference.wolfram.com/language/ref/Reals.html), [Integers](https://reference.wolfram.com/language/ref/Integers.html), and [Complexes](https://reference.wolfram.com/language/ref/Complexes.html).

## Details and Options

The system `*expr*` can be any logical combination of:

*lhs*==*rhs* | equations
*lhs*!=*rhs* | inequations
`*lhs*>*rhs*` or `*lhs*>=*rhs*` | inequalities
*expr*∈*dom* | domain specifications
{*x*,*y*,…}∈*reg* | region specification
[ForAll](https://reference.wolfram.com/language/ref/ForAll.html)[*x*,*cond*,*expr*] | universal quantifiers
[Exists](https://reference.wolfram.com/language/ref/Exists.html)[*x*,*cond*,*expr*] | existential quantifiers

[Solve](https://reference.wolfram.com/language/ref/Solve.html)[{*expr*_1,*expr*_2,…},*vars*] is equivalent to [Solve](https://reference.wolfram.com/language/ref/Solve.html)[*expr*_1&&*expr*_2&&…,*vars*].

A single variable or a list of variables can be specified.

[Solve](https://reference.wolfram.com/language/ref/Solve.html) gives solutions in terms of rules of the form:

{} | no solutions
{{*x*->*sol*_*x*,*y*->*sol*_*y*,…},…} | several solutions
{{}} | solution set is full dimensional

When a single variable is specified and a particular root of an equation has multiplicity greater than one, [Solve](https://reference.wolfram.com/language/ref/Solve.html) gives several copies of the corresponding solution.

[Solve](https://reference.wolfram.com/language/ref/Solve.html)[*expr*,*vars*] assumes by default that quantities appearing algebraically in inequalities are real, while all other quantities are complex.

[Solve](https://reference.wolfram.com/language/ref/Solve.html)[*expr*,*vars*,*dom*] restricts all variables and parameters to belong to the domain `*dom*`.

If `*dom*` is [Reals](https://reference.wolfram.com/language/ref/Reals.html), or a subset such as [Integers](https://reference.wolfram.com/language/ref/Integers.html) or [Rationals](https://reference.wolfram.com/language/ref/Rationals.html), then all constants and function values are also restricted to be real.

[Solve](https://reference.wolfram.com/language/ref/Solve.html)[*expr*&&*vars*∈[Reals](https://reference.wolfram.com/language/ref/Reals.html),*vars*,[Complexes](https://reference.wolfram.com/language/ref/Complexes.html)] solves for real values of variables, but function values are allowed to be complex.

[Solve](https://reference.wolfram.com/language/ref/Solve.html)[*expr*,*vars*,[Integers](https://reference.wolfram.com/language/ref/Integers.html)] solves Diophantine equations over the integers.

[Solve](https://reference.wolfram.com/language/ref/Solve.html)[…,*x*∈*reg*,[Reals](https://reference.wolfram.com/language/ref/Reals.html)] constrains `*x*` to be in the region `*reg*`. The different coordinates for `*x*` can be referred to using [Indexed](https://reference.wolfram.com/language/ref/Indexed.html)[*x*,*i*].

Algebraic variables in `*expr*` free of `*vars*` and of each other are treated as independent parameters.

[Solve](https://reference.wolfram.com/language/ref/Solve.html) deals primarily with linear and polynomial equations.

When `*expr*` involves only polynomial equations and inequalities over real or complex domains, then [Solve](https://reference.wolfram.com/language/ref/Solve.html) can always in principle solve directly for `*vars*`.

When `*expr*` involves transcendental conditions or integer domains, [Solve](https://reference.wolfram.com/language/ref/Solve.html) will often introduce additional parameters in its results.

[Solve](https://reference.wolfram.com/language/ref/Solve.html) can give explicit representations for solutions to all linear equations and inequalities over the integers and can solve a large fraction of Diophantine equations described in the literature.

When `*expr*` involves only polynomial conditions over real or complex domains, [Solve](https://reference.wolfram.com/language/ref/Solve.html)[*expr*,*vars*] will always be able to eliminate quantifiers.

[Solve](https://reference.wolfram.com/language/ref/Solve.html) gives generic solutions only. Solutions that are valid only when continuous parameters satisfy equations are removed. Other solutions that are only conditionally valid are expressed as [ConditionalExpression](https://reference.wolfram.com/language/ref/ConditionalExpression.html) objects.

Conditions included in [ConditionalExpression](https://reference.wolfram.com/language/ref/ConditionalExpression.html) solutions may involve inequalities, [Element](https://reference.wolfram.com/language/ref/Element.html) statements, equations and inequations on non-continuous parameters, and equations with full-dimensional solutions. Inequations and [NotElement](https://reference.wolfram.com/language/ref/NotElement.html) conditions on continuous parameters and variables are dropped.

[Solve](https://reference.wolfram.com/language/ref/Solve.html) may use non-equivalent transformations to find solutions of transcendental equations and hence it may not find some solutions and may not establish exact conditions on the validity of the solutions found. If this happens, an error message is issued.

[Solve](https://reference.wolfram.com/language/ref/Solve.html) uses special efficient techniques for handling sparse systems of linear equations with approximate numerical coefficients.

The following options can be given:

| [Assumptions](https://reference.wolfram.com/language/ref/Assumptions.html) | [$Assumptions](https://reference.wolfram.com/language/ref/$Assumptions.html) | assumptions on parameters |
| --- | --- | --- |
| [Cubics](https://reference.wolfram.com/language/ref/Cubics.html) | [Automatic](https://reference.wolfram.com/language/ref/Automatic.html) | whether to use explicit radicals to solve all cubics |
| [GeneratedParameters](https://reference.wolfram.com/language/ref/GeneratedParameters.html) | [C](https://reference.wolfram.com/language/ref/C.html) | how to name parameters that are generated |
| [InverseFunctions](https://reference.wolfram.com/language/ref/InverseFunctions.html) | [Automatic](https://reference.wolfram.com/language/ref/Automatic.html) | whether to use symbolic inverse functions |
| [MaxExtraConditions](https://reference.wolfram.com/language/ref/MaxExtraConditions.html) | 0 | how many extra equational conditions on continuous parameters to allow |
| [MaxRoots](https://reference.wolfram.com/language/ref/MaxRoots.html) | [Infinity](https://reference.wolfram.com/language/ref/Infinity.html) | maximum number of roots returned |
| [Method](https://reference.wolfram.com/language/ref/Method.html) | [Automatic](https://reference.wolfram.com/language/ref/Automatic.html) | what method should be used |
| [Modulus](https://reference.wolfram.com/language/ref/Modulus.html) | 0 | modulus to assume for integers |
| [Quartics](https://reference.wolfram.com/language/ref/Quartics.html) | [Automatic](https://reference.wolfram.com/language/ref/Automatic.html) | whether to use explicit radicals to solve all quartics |
| [VerifySolutions](https://reference.wolfram.com/language/ref/VerifySolutions.html) | [Automatic](https://reference.wolfram.com/language/ref/Automatic.html) | whether to verify solutions obtained using non-equivalent transformations |
| [WorkingPrecision](https://reference.wolfram.com/language/ref/WorkingPrecision.html) | [Infinity](https://reference.wolfram.com/language/ref/Infinity.html) | precision to be used in computations |

With [MaxExtraConditions](https://reference.wolfram.com/language/ref/MaxExtraConditions.html)->[Automatic](https://reference.wolfram.com/language/ref/Automatic.html), only solutions that require the minimal number of equational conditions on continuous parameters are included.

With [MaxExtraConditions](https://reference.wolfram.com/language/ref/MaxExtraConditions.html)->[All](https://reference.wolfram.com/language/ref/All.html), solutions that require arbitrary conditions on parameters are given and all conditions are included.

With [MaxExtraConditions](https://reference.wolfram.com/language/ref/MaxExtraConditions.html)->*k*, only solutions that require at most `*k*` equational conditions on continuous parameters are included.

With [Method](https://reference.wolfram.com/language/ref/Method.html)->[Reduce](https://reference.wolfram.com/language/ref/Reduce.html), [Solve](https://reference.wolfram.com/language/ref/Solve.html) uses only equivalent transformations and finds all solutions.

[Solve](https://reference.wolfram.com/language/ref/Solve.html)[*eqns*,…,[Modulus](https://reference.wolfram.com/language/ref/Modulus.html)->*m*] solves equations over the integers modulo `*m*`. With [Modulus](https://reference.wolfram.com/language/ref/Modulus.html)->[Automatic](https://reference.wolfram.com/language/ref/Automatic.html), [Solve](https://reference.wolfram.com/language/ref/Solve.html) will attempt to find the largest modulus for which the equations have solutions.

## Examples

### Basic Examples

Solve a quadratic equation:

```wolfram
Solve[x^2+a x+1==0,x]
(* Output *)
{{x->(1)/(2) (-a-Sqrt[-4+a^2])},{x->(1)/(2) (-a+Sqrt[-4+a^2])}}
```

Solve simultaneous equations in $x$ and $y$:

```wolfram
Solve[a x+y==7&&b x-y==1,{x,y}]
(* Output *)
{{x->(8)/(a+b),y->-(a-7 b)/(a+b)}}
```

Solve an equation over the reals:

```wolfram
Solve[(x^2+2)(x^2-2)==0,x,Reals]
(* Output *)
{{x->-Sqrt[2]},{x->Sqrt[2]}}
```

Solve an equation over the positive integers:

```wolfram
Solve[x^2+2y^3==3681&&x>0 &&y>0,{x,y},Integers]
(* Output *)
{{x->15,y->12},{x->41,y->10},{x->57,y->6}}
```

Solve equations in a geometric region:

```wolfram
Solve[{x,y}∈InfiniteLine[{{0,0},{2,1}}]&&{x,y}∈Circle[],{x,y}]
(* Output *)
{{x->-(2)/(Sqrt[5]),y->-(1)/(Sqrt[5])},{x->(2)/(Sqrt[5]),y->(1)/(Sqrt[5])}}
```

```wolfram
Graphics[{{Blue,InfiniteLine[{{0,0},{2,1}}],Circle[]},{PointSize[Large],Red,Point[{x,y}]/.%}}]
```

*([Graphics])*

### Scope

#### Basic Uses

Solutions are given as lists of replacements:

```wolfram
sol=Solve[x^2+y^2==2&&x-y==1,{x,y}]
(* Output *)
{{x->(1)/(2) (1-Sqrt[3]),y->(1)/(2) (-1-Sqrt[3])},{x->(1)/(2) (1+Sqrt[3]),y->(1)/(2) (-1+Sqrt[3])}}
```

Use [ReplaceAll](https://reference.wolfram.com/language/ref/ReplaceAll.html) (`/.`) to replace $x$ by solutions:

```wolfram
x/.sol
(* Output *)
{(1)/(2) (1-Sqrt[3]),(1)/(2) (1+Sqrt[3])}
```

Check that solutions satisfy the equations:

```wolfram
x^2+y^2==2&&x-y==1/.sol
(* Output *)
{True,True}
```

[Solve](https://reference.wolfram.com/language/ref/Solve.html) uses `{}` to represent the empty solution or no solution:

```wolfram
Solve[x==1&&x==2,x]
(* Output *)
{}
```

[Solve](https://reference.wolfram.com/language/ref/Solve.html) uses `{{}}` to represent the universal solution or all points satisfying the equations:

```wolfram
Solve[x==x,x]
(* Output *)
{{}}
```

Find solutions over specified domains:

```wolfram
Solve[(x^4-1)(x^4-4)==0,x,Complexes]
(* Output *)
{{x->-1},{x->-ⅈ},{x->ⅈ},{x->1},{x->-Sqrt[2]},{x->-ⅈ Sqrt[2]},{x->ⅈ Sqrt[2]},{x->Sqrt[2]}}
```

```wolfram
Solve[(x^4-1)(x^4-4)==0,x,Reals]
(* Output *)
{{x->-1},{x->1},{x->-Sqrt[2]},{x->Sqrt[2]}}
```

```wolfram
Solve[(x^4-1)(x^4-4)==0,x,Integers]
(* Output *)
{{x->-1},{x->1}}
```

Solve equations with coefficients involving a symbolic parameter:

```wolfram
sol=Solve[x^2+y^2==1&&x+y==a,{x,y}]
(* Output *)
{{x->(1)/(2) (a-Sqrt[2-a^2]),y->(1)/(2) (a+Sqrt[2-a^2])},{x->(1)/(2) (a+Sqrt[2-a^2]),y->(1)/(2) (a-Sqrt[2-a^2])}}
```

Plot the real parts of the solutions for `y` as a function of the parameter `a`:

```wolfram
Plot[Evaluate[Re[y/.sol]],{a,0,2}]
```

*([Graphics])*

Solution of this equation over the reals requires conditions on the parameters:

```wolfram
sol=Solve[x^2+a x+b==0,x,Reals]
(* Output *)
{{x->-(a)/(2)-(1)/(2) Sqrt[a^2-4 b]},{x->-(a)/(2)+(1)/(2) Sqrt[a^2-4 b]}}
```

Replace `x` by solutions and simplify the results:

```wolfram
x^2+a x/.sol//Simplify
(* Output *)
{-b,-b}
```

Use [Normal](https://reference.wolfram.com/language/ref/Normal.html) to remove the conditions:

```wolfram
Normal[sol]
(* Output *)
{{x->-(a)/(2)-(1)/(2) Sqrt[a^2-4 b]},{x->-(a)/(2)+(1)/(2) Sqrt[a^2-4 b]}}
```

Solution of this equation over the positive integers requires introduction of a new parameter:

```wolfram
sol=Solve[x^2-2y^2==1&&x>0&&y>0,{x,y},Integers]
(* Output *)
{{x->(1)/(2) ((3-2 Sqrt[2])^1+(3+2 Sqrt[2])^1),y->-((3-2 Sqrt[2])^1-(3+2 Sqrt[2])^1)/(2 Sqrt[2])}}
```

List the first 10 solutions:

```wolfram
{x,y}/.First[sol]/.Table[{1->i},{i,10}]//Simplify
(* Output *)
{{3,2},{17,12},{99,70},{577,408},{3363,2378},{19601,13860},{114243,80782},{665857,470832},{3880899,2744210},{22619537,15994428}}
```

#### Complex Equations in One Variable

Polynomial equations solvable in radicals:

```wolfram
Solve[x^4-x^2-5==0,x]
(* Output *)
{{x->-ⅈ Sqrt[(1)/(2) (-1+Sqrt[21])]},{x->ⅈ Sqrt[(1)/(2) (-1+Sqrt[21])]},{x->-Sqrt[(1)/(2) (1+Sqrt[21])]},{x->Sqrt[(1)/(2) (1+Sqrt[21])]}}
```

To use general formulas for solving cubic equations, set [Cubics](https://reference.wolfram.com/language/ref/Cubics.html)->[True](https://reference.wolfram.com/language/ref/True.html):

```wolfram
Solve[x^3-2x^2+5x+7==0,x,Cubics->True]
(* Output *)
{{x->(1)/(3) (2-11 ((2)/(-263+3 Sqrt[8277]))^(1/3)+((1)/(2) (-263+3 Sqrt[8277]))^(1/3))},{x->(2)/(3)-(1)/(6) (1+ⅈ Sqrt[3]) ((1)/(2) (-263+3 Sqrt[8277]))^(1/3)+(11 (1-ⅈ Sqrt[3]))/(3 2^(2/3) (-263+3 Sqrt[8277])^(1/3))},{x->(2)/(3)-(1)/(6) (1-ⅈ Sqrt[3]) ((1)/(2) (-263+3 Sqrt[8277]))^(1/3)+(11 (1+ⅈ Sqrt[3]))/(3 2^(2/3) (-263+3 Sqrt[8277])^(1/3))}}
```

By default, [Solve](https://reference.wolfram.com/language/ref/Solve.html) uses [Root](https://reference.wolfram.com/language/ref/Root.html) objects to represent solutions of general cubic equations with numeric coefficients:

```wolfram
Solve[x^3-2x^2+5x+7==0,x]
(* Output *)
{{x->Root},{x->Root},{x->Root}}
```

General polynomial equations:

```wolfram
Solve[x^5-2x+17==0,x]
(* Output *)
{{x->Root},{x->Root},{x->Root},{x->Root},{x->Root}}
```

Polynomial equations with multiple roots:

```wolfram
Solve[(x^2-1)(x^4-1)==0,x]
(* Output *)
{{x->-1},{x->-1},{x->-ⅈ},{x->ⅈ},{x->1},{x->1}}
```

Find five roots of a polynomial of a high degree:

```wolfram
Solve[x^1234567+9x^2+7x-1==0,x,MaxRoots->5]
(* Output *)
{{x->Root},{x->Root},{x->Root},{x->Root},{x->Root}}
```

Polynomial equations with symbolic coefficients:

```wolfram
Solve[a x^2+b x+c==0,x]
(* Output *)
{{x->(-b-Sqrt[b^2-4 a c])/(2 a)},{x->(-b+Sqrt[b^2-4 a c])/(2 a)}}
```

```wolfram
Solve[a x^5+b x+c==0,x]
(* Output *)
{{x->Root[c+b #1+a #1^5&,1]},{x->Root[c+b #1+a #1^5&,2]},{x->Root[c+b #1+a #1^5&,3]},{x->Root[c+b #1+a #1^5&,4]},{x->Root[c+b #1+a #1^5&,5]}}
```

Algebraic equations:

```wolfram
Solve[Sqrt[x]+3x^(1/3)==5,x]
(* Output *)
{{x->-2-957 ((2)/(1217+5125 Sqrt[5]))^(1/3)+3 ((1)/(2) (1217+5125 Sqrt[5]))^(1/3)}}
```

Transcendental equations:

```wolfram
Solve[Sin[x]==1/3,x]
(* Output *)
{{x->π-ArcSin[(1)/(3)]+2 π 1},{x->ArcSin[(1)/(3)]+2 π 1}}
```

```wolfram
Solve[x^(2a)+2x^a+1==0,x]
(* Output *)
Solve
(* Output *)
{{x->(-1)^((1)/(a))}}
```

```wolfram
Solve[(x^5-1)^x==0,x]
(* Output *)
{{x->1},{x->(-1)^(2/5)},{x->-(-1)^(3/5)}}
```

```wolfram
Solve[ArcTan[1/Sqrt[1+x^2],x/Sqrt[1+x^2]]==a,x]
(* Output *)
Solve
(* Output *)
{{x->-(ⅈ (-1+Cos[2 a]+ⅈ Sin[2 a]))/(1+Cos[2 a]+ⅈ Sin[2 a])}}
```

[Solve](https://reference.wolfram.com/language/ref/Solve.html) cannot find all solutions here:

```wolfram
Solve[5==x*2^(x^2),x]
(* Output *)
Solve
(* Output *)
{{x->Sqrt[(ProductLog[50 Log[2]])/(2 Log[2])]}}
```

Find three solutions:

```wolfram
Solve[5==x*2^(x^2),x,MaxRoots->3]
(* Output *)
{{x->Root},{x->Root},{x->Root}}
```

Univariate elementary function equations over bounded regions:

```wolfram
Solve[Sin[E^x]-Cos[2 x]==1&&-1<=Re[x]<=1&&-1<=Im[x]<=1,x]
(* Output *)
{{x->Root},{x->Root}}
```

Univariate holomorphic function equations over bounded regions:

```wolfram
Solve[Gamma[x]-Log[x]==I/2&&Abs[x-2]<3/2,x]
(* Output *)
{{x->Root},{x->Root}}
```

Here [Solve](https://reference.wolfram.com/language/ref/Solve.html) finds some solutions but is not able to prove there are no other solutions:

```wolfram
Solve[x==E^(1/x)&&Abs[x]<5,x]
(* Output *)
Solve
(* Output *)
{{x->Root},{x->Root},{x->Root},{x->Root},{x->Root}}
```

Equation with a purely imaginary period over a vertical stripe in the complex plane:

```wolfram
Solve[Cos[Exp[x]]==3 Exp[-x]+1&&0<=Re[x]<=1,x]
(* Output *)
{{x->2 ⅈ π 1+Root},{x->2 ⅈ π 1+Root},{x->2 ⅈ π 1+Root}}
```

Find a specified number of roots of an unrestricted complex equation:

```wolfram
Solve[Sin[FresnelS[x]+BesselJ[3,x^2-1]]==2^Cos[x]-3,x,MaxRoots->5]
(* Output *)
{{x->Root},{x->Root},{x->Root},{x->Root},{x->Root}}
```

Symbolic functions:

```wolfram
Solve[f[x]^3+2f[x]-3==0,x]
(* Output *)
Solve
(* Output *)
{{x->f^((-1))[1]},{x->f^((-1))[(1)/(2) (-1-ⅈ Sqrt[11])]},{x->f^((-1))[(1)/(2) (-1+ⅈ Sqrt[11])]}}
```

Nonanalytic complex equations:

```wolfram
Solve[{Re[x]+Im[x]^2==1,Re[x]-Im[x]==-1},x]
(* Output *)
{{x->-3-2 ⅈ},{x->ⅈ}}
```

```wolfram
Solve[{Abs[x]==1,Arg[x]==π/4},x]
(* Output *)
{{x->(1+ⅈ)/(Sqrt[2])}}
```

#### Systems of Complex Equations in Several Variables

Systems of linear equations:

```wolfram
Solve[{x+2y+3z==4,3x+4y+5z==6,7x+7y+8z==9},{x,y,z}]
(* Output *)
{{x->0,y->-1,z->2}}
```

Linear equations with symbolic coefficients:

```wolfram
Solve[a x+b y+c z==d &&3x+4y+5z==6&&7x+7y+8z==9,{x,y,z}]
(* Output *)
{{x->(3 (b-2 c+d))/(3 a-11 b+7 c),y->-1-(11 (b-2 c+d))/(3 a-11 b+7 c),z->2+(7 (b-2 c+d))/(3 a-11 b+7 c)}}
```

Underdetermined systems of linear equations:

```wolfram
Solve[x+2y+3z==4&&3x+4y+5z==6&&6x+7y+8z==9,{x,y,z}]
(* Output *)
Solve
(* Output *)
{{y->-1-2 x,z->2+x}}
```

Linear equations with no solutions:

```wolfram
Solve[x+2y+3z==4&&3x+4y+5z==6&&6x+7y+8z==0,{x,y,z}]
(* Output *)
{}
```

Systems of polynomial equations:

```wolfram
Solve[x^2+y^2==1&&x+2y==3,{x,y}]
(* Output *)
{{x->(3)/(5)-(4 ⅈ)/(5),y->(6)/(5)+(2 ⅈ)/(5)},{x->(3)/(5)+(4 ⅈ)/(5),y->(6)/(5)-(2 ⅈ)/(5)}}
```

Find five out of a trillion roots of a polynomial system:

```wolfram
Solve[x^10000==y^2+3y+2&&y^10000==z^2+3z+2&& z^10000==x^2+3x+2,{x,y,z},MaxRoots->5]
(* Output *)
{{x->Root,y->Root,z->Root},{x->Root,y->Root,z->Root},{x->Root,y->Root,z->Root},{x->Root,y->Root,z->Root},{x->Root,y->Root,z->Root}}
```

Polynomial equations with symbolic coefficients:

```wolfram
Solve[a x^2+b y^2==c&&x+2y==3,{x,y}]
(* Output *)
{{x->(3 b-2 Sqrt[-9 a b+4 a c+b c])/(4 a+b),y->(1)/(2) (3-(3 b)/(4 a+b)+(2 Sqrt[-9 a b+4 a c+b c])/(4 a+b))},{x->(3 b+2 Sqrt[-9 a b+4 a c+b c])/(4 a+b),y->(1)/(2) (3-(3 b)/(4 a+b)-(2 Sqrt[-9 a b+4 a c+b c])/(4 a+b))}}
```

Algebraic equations:

```wolfram
Solve[Sqrt[x y]==Sqrt[x+y]&&Sqrt[x]-y^(1/3)==1,{x,y}]
(* Output *)
Solve
(* Output *)
{{x->Root,y->1-7 Root-4 Root^3+Root^4}}
```

Transcendental equations:

```wolfram
Solve[Sin[x+y]==1&&Log[x-y]==1/2,{x,y}]
(* Output *)
{{x->(1)/(2) (Sqrt[ℯ]+(π)/(2)+2 π 1),y->(1)/(2) (-Sqrt[ℯ]+(π)/(2)+2 π 1)}}
```

```wolfram
Solve[3^x^2==7&&x^2-y^3==4,{x,y}]
(* Output *)
Solve
(* Output *)
{{x->-Sqrt[(Log[7])/(Log[3])],y->-(-(1)/(Log[3]))^(1/3) (-4 Log[3]+Log[7])^(1/3)},{x->Sqrt[(Log[7])/(Log[3])],y->-(-(1)/(Log[3]))^(1/3) (-4 Log[3]+Log[7])^(1/3)},{x->-Sqrt[(Log[7])/(Log[3])],y->((-4 Log[3]+Log[7])/(Log[3]))^(1/3)},{x->Sqrt[(Log[7])/(Log[3])],y->((-4 Log[3]+Log[7])/(Log[3]))^(1/3)},{x->-Sqrt[(Log[7])/(Log[3])],y->(-1)^(2/3) ((-4 Log[3]+Log[7])/(Log[3]))^(1/3)},{x->Sqrt[(Log[7])/(Log[3])],y->(-1)^(2/3) ((-4 Log[3]+Log[7])/(Log[3]))^(1/3)}}
```

```wolfram
Solve[x^a+y^b==1&&(x y)^(a b)==2,{x,y}]
(* Output *)
Solve
(* Output *)
Solve
(* Output *)
{{y->2^((1)/(a b)) (x^a)^(-1/a)}}
```

Find a specified number of solutions of transcendental equations:

```wolfram
Solve[Sin[x+y]==x y+1&&Cos[x-y]==AiryAi[x y]+2,{x,y},MaxRoots->3]
(* Output *)
{{x->Root,y->Root},{x->Root,y->Root},{x->Root,y->Root}}
```

Square analytic systems over bounded boxes:

```wolfram
Solve[FresnelS[x-y]-AiryAi[x]+y==(3+I)/2&&CosIntegral[x] y-Sin[x y]==-5(1-I)/4&&1<=Re[x]<=2&&0<=Im[x]<=1&&1<=Re[y]<=2&&0<=Im[y]<=1,{x,y}]
(* Output *)
{{x->Root,y->Root}}
```

Nonanalytic equations:

```wolfram
Solve[{Re[x^2]+Im[y]==1,x^2-y^2==0,Im[x^3]-Im[y]^2==2},{x,y}]
(* Output *)
{{x->ⅈ AlgebraicNumber+AlgebraicNumber,y->AlgebraicNumber+ⅈ AlgebraicNumber},{x->ⅈ AlgebraicNumber+AlgebraicNumber,y->AlgebraicNumber+ⅈ AlgebraicNumber},{x->AlgebraicNumber+ⅈ AlgebraicNumber,y->AlgebraicNumber+ⅈ AlgebraicNumber},{x->AlgebraicNumber+ⅈ AlgebraicNumber,y->AlgebraicNumber+ⅈ AlgebraicNumber}}
```

#### Real Equations in One Variable

Polynomial equations:

```wolfram
Solve[x^5-2x+1==0,x,Reals]
(* Output *)
{{x->1},{x->Root},{x->Root}}
```

Polynomial equations with multiple roots:

```wolfram
Solve[(x^2-1)(x^4-1)==0,x,Reals]
(* Output *)
{{x->-1},{x->-1},{x->1},{x->1}}
```

Polynomial equations with symbolic coefficients:

```wolfram
Solve[x^3+a x+a^2==0,x,Reals]
(* Output *)
{{x->Root[a^2+a #1+#1^3&,1]},{x->Root[a^2+a #1+#1^3&,2]},{x->Root[a^2+a #1+#1^3&,3]}}
```

Algebraic equations:

```wolfram
Solve[Sqrt[x]+3x^(1/3)==5,x,Reals]
(* Output *)
{{x->Root^6}}
```

Piecewise equations:

```wolfram
Solve[Abs[x]^2-x+UnitStep[x]==9,x,Reals]
(* Output *)
{{x->(1)/(2) (1+Sqrt[33])},{x->(1)/(2) (1-Sqrt[37])}}
```

Transcendental equations, solvable using inverse functions:

```wolfram
Solve[E^x-x==7,x,Reals]
(* Output *)
{{x->-7-ProductLog[-(1)/(ℯ^7)]},{x->-7-ProductLog[-1,-(1)/(ℯ^7)]}}
```

```wolfram
Solve[ (27^(2x-1))^(1/x)==Sqrt[9^(2x-1)], x , Reals]
(* Output *)
{{x->(1)/(2)},{x->3}}
```

```wolfram
Solve[JacobiSN[x,y]==1,x,Reals]
(* Output *)
{{x->EllipticK[y]+4 1 EllipticK[y]}}
```

Transcendental equations, solvable using special function zeros:

```wolfram
Solve[AiryBi[1-x^2]==0&&2<x<3,x,Reals]
(* Output *)
{{x->Sqrt[1-AiryBiZero[2]]},{x->Sqrt[1-AiryBiZero[3]]},{x->Sqrt[1-AiryBiZero[4]]},{x->Sqrt[1-AiryBiZero[5]]}}
```

Transcendental inequalities, solvable using special function zeros:

```wolfram
Solve[900<AiryAiZero[2t+1]^2<1000,t,Reals]
(* Output *)
{{t->(35)/(2)},{t->18}}
```

Exp-log equations:

```wolfram
Solve[E^(2E^x)-Log[x^2+1]-20x==11,x,Reals]
(* Output *)
{{x->Root},{x->Root}}
```

High-degree sparse polynomial equations:

```wolfram
Solve[x^1000000-2x^777777+3x^12345+9x^67-10==0,x,Reals]
(* Output *)
{{x->Root},{x->Root},{x->Root},{x->Root}}
```

Algebraic equations involving high-degree radicals:

```wolfram
Solve[2x^(123451/67890)-x^2+4Sqrt[x]-4x-9/8==0,x,Reals]
(* Output *)
{{x->Root^67890},{x->Root^67890},{x->Root^67890},{x->Root^67890}}
```

Equations involving non-rational real powers:

```wolfram
Solve[x^Pi-x^x^Sqrt[2]-Sqrt[3]x+2^(1/3)==0,x,Reals]
(* Output *)
{{x->Root},{x->Root},{x->Root}}
```

Equation with a double root:

```wolfram
Solve[E^(2x)+x^4+4(x^2+1)==(2x^2+4)E^x,x,Reals]
(* Output *)
{{x->Root},{x->Root}}
```

Tame elementary function equations:

```wolfram
Solve[10Sin[Tan[E^-x^2]]-x==3,x,Reals]
(* Output *)
{{x->Root},{x->Root},{x->Root}}
```

Elementary function equations in bounded intervals:

```wolfram
Solve[2Sin[Exp[x]]-Cos[Pi x]==3/2&&-1<x<1,x,Reals]
(* Output *)
{{x->Root},{x->Root}}
```

Holomorphic function equations in bounded intervals:

```wolfram
Solve[Cos[x]-BesselJ[5,x]==1/2&&0<=x<=10,x]
(* Output *)
{{x->Root},{x->Root},{x->Root}}
```

Periodic elementary function equations over the reals:

```wolfram
Solve[Exp[Sin[x]]-Sin[3 Cos[x]]==0,x,Reals]
(* Output *)
{{x->2 π 1+Root},{x->2 π 1+Root}}
```

#### Systems of Real Equations and Inequalities in Several Variables

Linear systems:

```wolfram
Solve[2 x+3y-5z==1&&3x-4y+7z==3,{y,z},Reals]
(* Output *)
{{y->22-29 x,z->13-17 x}}
```

Polynomial systems:

```wolfram
Solve[x y==z^2-x&&x y z==2&&x^2+y^2+z^2<=5,{y,z},Reals]
(* Output *)
{{y->(-x+Root[-2-x #1+#1^3&,1]^2)/(x),z->Root[-2-x #1+#1^3&,1]}}
```

Quantified polynomial systems:

```wolfram
Solve[Exists[x, x^2+a x+b==0&&2 x+a==0],a,Reals]
(* Output *)
{{a->-2 Sqrt[b]},{a->2 Sqrt[b]}}
```

```wolfram
Solve[ForAll[x,Exists[y,a x^2+b y^2-3y==1&&y<0&&a y-y==b+1]],{a,b},Reals]
(* Output *)
{{a->0,b->(1)/(3) (-2-(8)/((1+3 Sqrt[57])^(1/3))+(1+3 Sqrt[57])^(1/3))},{a->1,b->-1}}
```

Algebraic systems:

```wolfram
Solve[Sqrt[x+2y]-3x+4y>=5&&x+y^(1/3)==1,y,Reals]
(* Output *)
{{y->1-3 x+3 x^2-x^3}}
```

Piecewise systems:

```wolfram
Solve[Max[x,y]==Min[y^2-x,x]+y&&z+UnitStep[x-y]==1,{y,z},Reals]
(* Output *)
{{y->Sqrt[x],z->1},{y->-(1)/(2)+(1)/(2) Sqrt[1+8 x],z->0}}
```

Transcendental systems, solvable using inverse functions:

```wolfram
Solve[Sin[x+y]==1/2&&E^x-y<=1,y,Reals]
(* Output *)
{{y->(1)/(6) (π-6 x+12 π 1)},{y->(1)/(6) (5 π-6 x+12 π 1)}}
```

```wolfram
Solve[ 3^x-2^2y==77 && Sqrt[3^x]-2^y==7, {x, y}, Reals]
(* Output *)
{{x->4,y->1}}
```

Systems exp-log in the first variable and polynomial in the other variables:

```wolfram
Solve[E^x y^3+Log[x]y==1&&x y+E^x/x>=2,y,Reals]
(* Output *)
{{y->Root[-1+Log[x] #1+ℯ^x #1^3&,1]},{y->Root[-1+Log[x] #1+ℯ^x #1^3&,2]},{y->Root[-1+Log[x] #1+ℯ^x #1^3&,3]}}
```

Quantified system:

```wolfram
Solve[Exists[a,a x^2+Sinh[x^2+1]a^2==1&&x^2-a^2==1],x,Reals]
(* Output *)
{{x->Root},{x->Root},{x->Root},{x->Root}}
```

Systems elementary and bounded in the first variable and polynomial in the other variables:

```wolfram
Solve[Sin[x-Cos[x]] y^3-x==1&&x^2+y^2<=1,y,Reals]
(* Output *)
{{y->Root[-1-x+Sin[x-Cos[x]] #1^3&,1]}}
```

Quantified system:

```wolfram
Solve[Exists[y,y^3-Cos[x] y+2 x^2 Sin[x^2-1]==0&&x^2+y^2==2],x,Reals]
(* Output *)
{{x->Root},{x->Root}}
```

Systems analytic and bounded in the first variable and polynomial in the other variables:

```wolfram
Solve[y^3-BesselJ[2,x+2] y-y-3 x==-2&&y<0 &&x^2<2,y,Reals]
(* Output *)
{{y->Root[2-3 x+(-1-BesselJ[2,2+x]) #1+#1^3&,1]},{y->Root[2-3 x+(-1-BesselJ[2,2+x]) #1+#1^3&,2]}}
```

Quantified system:

```wolfram
Solve[Exists[y,y^4-Gamma[x+2] y-y-3 ArcSin[x/3]==1&&x^2+y^3==1&&0<x<2],x,Reals]
(* Output *)
{{x->Root}}
```

Square systems of analytic equations over bounded regions:

```wolfram
Solve[Gamma[x+y+1]-Sin[x y]==1&&Erf[x^2-y]-E^y-x+4==0&&0<x<3&&0<y<3,{x,y}, Reals]
(* Output *)
{{x->Root,y->Root},{x->Root,y->Root}}
```

#### Diophantine Equations

Linear systems of equations:

```wolfram
Solve[2 x+3y-5z==1&&3x-4y+7z==3,{x,y,z},Integers]
(* Output *)
{{x->1,y->22-29 1,z->13-17 1}}
```

Linear systems of equations and inequalities:

```wolfram
Solve[2 x+3y==4&&3x-4y<=5&&x-2y>-21,{x,y},Integers]
(* Output *)
{{x->-7,y->6},{x->-4,y->4},{x->-1,y->2}}
```

Univariate polynomial equations:

```wolfram
Solve[x^12345-2x^777+1==0,x,Integers]
(* Output *)
{{x->1}}
```

Binary quadratic equations:

```wolfram
Solve[x^2+x y+y^2==109,{x,y},Integers]
(* Output *)
{{x->-12,y->5},{x->-12,y->7},{x->-7,y->-5},{x->-7,y->12},{x->-5,y->-7},{x->-5,y->12},{x->5,y->-12},{x->5,y->7},{x->7,y->-12},{x->7,y->5},{x->12,y->-7},{x->12,y->-5}}
```

```wolfram
Solve[x^2-3y^2==22&&x>0&&y>0,{x,y},Integers]
(* Output *)
{{x->5,y->1},{x->(1)/(2) (5 (2-Sqrt[3])^1+Sqrt[3] (2-Sqrt[3])^1+5 (2+Sqrt[3])^1-Sqrt[3] (2+Sqrt[3])^1),y->(1)/(6) (-3 (2-Sqrt[3])^1-5 Sqrt[3] (2-Sqrt[3])^1-3 (2+Sqrt[3])^1+5 Sqrt[3] (2+Sqrt[3])^1)},{x->(1)/(2) (5 (2-Sqrt[3])^1-Sqrt[3] (2-Sqrt[3])^1+5 (2+Sqrt[3])^1+Sqrt[3] (2+Sqrt[3])^1),y->(1)/(6) (3 (2-Sqrt[3])^1-5 Sqrt[3] (2-Sqrt[3])^1+3 (2+Sqrt[3])^1+5 Sqrt[3] (2+Sqrt[3])^1)}}
```

```wolfram
Solve[x^2-6 x y+9y^2-x+2y==1,{x,y},Integers]
(* Output *)
{{x->-3-2 1+3 1^2,y->-1-1+1^2}}
```

Thue equations:

```wolfram
Solve[x^3-2x^2 y+y^3==2,{x,y},Integers]
(* Output *)
{{x->1,y->-1},{x->5,y->3}}
```

Sum of squares equations:

```wolfram
Solve[x^2+4y^2+9z^2+16t^2==354&&x>0&&y>0&&z>0&&t>0,{x,y,z,t},Integers]
(* Output *)
{{x->1,y->2,z->3,t->4},{x->1,y->4,z->5,t->2},{x->1,y->8,z->3,t->1},{x->5,y->4,z->1,t->4},{x->5,y->8,z->1,t->2},{x->7,y->2,z->5,t->2},{x->7,y->4,z->5,t->1}}
```

The Pythagorean equation:

```wolfram
Solve[x^2+y^2==z^2,{x,y,z},Integers]
(* Output *)
{{x->2 1 2 3,y->1 (2^2-3^2),z->1 (2^2+3^2)},{x->1 (2^2-3^2),y->2 1 2 3,z->1 (2^2+3^2)}}
```

Bounded systems of equations and inequalities:

```wolfram
Solve[x^4+y^4+z^4<=500&&x+y^2+z^3==32,{x,y,z},Integers]
(* Output *)
{{x->-4,y->-3,z->3},{x->-4,y->3,z->3},{x->1,y->-2,z->3},{x->1,y->2,z->3},{x->4,y->-1,z->3},{x->4,y->1,z->3}}
```

High[Hyphen]degree systems with no solutions:

```wolfram
Solve[2x^7+8y^15+14 x y z==3,{x,y,z},Integers]
(* Output *)
{}
```

Transcendental Diophantine systems:

```wolfram
Solve[Exp[y^2]<x&&Abs[x]<5&&Abs[y]<5,{x,y},Integers]
(* Output *)
{{x->2,y->0},{x->3,y->-1},{x->3,y->0},{x->3,y->1},{x->4,y->-1},{x->4,y->0},{x->4,y->1}}
```

```wolfram
Solve[Exp[x^2-5y^2+1]+x^2-5y^2==0&&x>0&&y>0,{x,y},Integers]
(* Output *)
{{x->2,y->1},{x->(1)/(2) (-2 (9-4 Sqrt[5])^1-Sqrt[5] (9-4 Sqrt[5])^1-2 (9+4 Sqrt[5])^1+Sqrt[5] (9+4 Sqrt[5])^1),y->(1)/(10) (5 (9-4 Sqrt[5])^1+2 Sqrt[5] (9-4 Sqrt[5])^1+5 (9+4 Sqrt[5])^1-2 Sqrt[5] (9+4 Sqrt[5])^1)}}
```

Polynomial systems of congruences:

```wolfram
Solve[Mod[x^2+y^2,2]==1&&Mod[x-2y,3]==2,{x,y},Integers]
(* Output *)
{{x->6 1,y->5+6 2},{x->1+6 1,y->4+6 2},{x->2+6 1,y->3+6 2},{x->3+6 1,y->2+6 2},{x->4+6 1,y->1+6 2},{x->5+6 1,y->6 2}}
```

#### Modular Equations

Linear systems:

```wolfram
Solve[2 x+3y-5z==1&&3x-4y+7z==3&&2x-10y+3z==5,{x,y,z},Modulus->12]
(* Output *)
{{z->7,y->4,x->6}}
```

```wolfram
Solve[2 x+3y-5z==1&&3x-4y+7z==3,{x,y,z},Modulus->12]
(* Output *)
Solve
(* Output *)
{{z->1+7 x,y->10+7 x}}
```

Univariate polynomial equations:

```wolfram
Solve[x^3-2x+1==0,x,Modulus->5]
(* Output *)
{{x->1},{x->2},{x->2}}
```

Systems of polynomial equations and inequations:

```wolfram
Solve[x^2+y^3==1&&x+2y^2==4,{x,y},Modulus->7]
(* Output *)
{{x->0,y->4},{x->3,y->5}}
```

```wolfram
Solve[x^2+y^3==z&&x+2y==3z+1&&x y z≠0,{x,y,z},Modulus->7]
(* Output *)
{{x->5,y->2,z->5},{x->5,y->6,z->3},{x->6,y->4,z->2}}
```

Quantified polynomial systems:

```wolfram
Solve[ForAll[x,Exists[y,a x^2+b y^2-3y==1&&y≠0]],{a,b},Modulus->3]
(* Output *)
{{a->0,b->1}}
```

#### Equations over Finite Fields

Univariate equations:

```wolfram
ℱ=FiniteField[53,4];
Solve[x^5+ℱ[123]x==ℱ[234],x]
(* Output *)
{{x-><|interpretation -> FiniteFieldElement[FiniteField[53, 2, +, 38, #, +, 9, #, ^, 2, +, #, ^, 4, &, Polynomial], 10134132], index -> 4879932, shortIndex -> 4879932, indexShortened -> True, characteristic -> 53, shortCharacteristic -> 53, extensionDegree -> 4, field -> FiniteField[...], fieldDisplayed -> False|>}}
```

```wolfram
Solve[x^7+2 x+3==0,x,ℱ]
(* Output *)
![image](img/image_001.png)
```

Systems of linear equations:

```wolfram
ℱ=FiniteField[71,2];
Solve[ℱ[123]x+ℱ[234]y==ℱ[345]&&ℱ[321]x+ℱ[432]y==ℱ[543],{x,y}]
(* Output *)
{{x-><|interpretation -> FiniteFieldElement[FiniteField[71, 7, +, 69, #, +, #, ^, 2, &, Polynomial], 1114], index -> 1005, shortIndex -> 1005, indexShortened -> True, characteristic -> 71, shortCharacteristic -> 71, extensionDegree -> 2, field -> FiniteField[...], fieldDisplayed -> False|>,y-><|interpretation -> FiniteFieldElement[FiniteField[71, 7, +, 69, #, +, #, ^, 2, &, Polynomial], 6157], index -> 4108, shortIndex -> 4108, indexShortened -> True, characteristic -> 71, shortCharacteristic -> 71, extensionDegree -> 2, field -> FiniteField[...], fieldDisplayed -> False|>}}
```

```wolfram
Solve[ℱ[1234]x+ℱ[2345]y+ℱ[3456]z==ℱ[4567]&&ℱ[1]x+ℱ[2]y+ℱ[3]z==ℱ[4],{y,z}]
(* Output *)
![image](img/image_003.png)
```

Systems of polynomial equations:

```wolfram
ℱ=FiniteField[7,5];
Solve[x^2+y^2==3&&x^5+y^5==5,{x,y},ℱ]
(* Output *)
![image](img/image_005.png)
```

```wolfram
Solve[ℱ[123]x^2+ℱ[234]y^3+ℱ[345]z^4==ℱ[456]&&ℱ[21]x+ℱ[32]y^2+ℱ[43]z^3==ℱ[54]&&x y z==ℱ[1],{x,y,z}]
(* Output *)
![image](img/image_007.png)
```

#### Systems with Mixed Variable Domains

Mixed real and complex variables:

```wolfram
Solve[x^2+y^2==1&&Element[x,Reals],y]
(* Output *)
{{y->-Sqrt[1-x^2]},{y->Sqrt[1-x^2]},{y->-ⅈ Sqrt[-1+x^2]},{y->ⅈ Sqrt[-1+x^2]}}
```

Mixed real and integer variables:

```wolfram
Solve[x^2+y^2==7 &&Element[x,Integers]&&Element[y,Reals],{x,y}]
(* Output *)
{{x->-2,y->-Sqrt[3]},{x->-2,y->Sqrt[3]},{x->-1,y->-Sqrt[6]},{x->-1,y->Sqrt[6]},{x->0,y->-Sqrt[7]},{x->0,y->Sqrt[7]},{x->1,y->-Sqrt[6]},{x->1,y->Sqrt[6]},{x->2,y->-Sqrt[3]},{x->2,y->Sqrt[3]}}
```

#### Systems with Geometric Region Constraints

Solve over special regions in 2D:

```wolfram
ℛ_1=Circle[];
ℛ_2=Line[{{-2,1},{1,-2}}];
```

```wolfram
Solve[{x,y}∈ℛ_1&&{x,y}∈ℛ_2,{x,y}]
(* Output *)
{{x->-1,y->0},{x->0,y->-1}}
```

Plot it:

```wolfram
Graphics[{{Blue,ℛ_1,ℛ_2},{Red,Point[{x,y}]/.%}}]
```

*([Graphics])*

Solve over special regions in 3D:

```wolfram
ℛ_1=Sphere[];
ℛ_2=InfinitePlane[{{0,0,0},{0,1,0},{1,0,1}}];
```

```wolfram
Solve[2 x y==z^2&&{x,y,z}∈ℛ_1&&{x,y,z}∈ℛ_2,{x,y,z},Reals]
(* Output *)
{{x->-(2)/(3),y->-(1)/(3),z->-(2)/(3)},{x->0,y->-1,z->0},{x->0,y->1,z->0},{x->(2)/(3),y->(1)/(3),z->(2)/(3)}}
```

Plot it:

```wolfram
Show[{ContourPlot3D[2 x y==z^2,{x,-1.2,1.2},{y,-1.2,1.2},{z,-1.2,1.2},Mesh->None,ContourStyle->Opacity[0.5]],Graphics3D[{{Opacity[0.5],Green,ℛ_1},{Opacity[0.5],Yellow,ℛ_2},{PointSize[Large],Red,Point[{x,y,z}/.%]}}]}]
```

*([Graphics3D])*

A quantified system:

```wolfram
ℛ=Cone[{{0,0,0},{1,2,3}},4];
```

```wolfram
Solve[∃_z(x^2+y^2==z^2&&x+y==z-1&&{x,y,z}∈ℛ),y,Reals]
(* Output *)
{{y->(-1-2 x)/(2+2 x)}}
```

An implicitly defined region:

```wolfram
ℛ=ImplicitRegion[a+2 b-3 c>=1&&a b c==7,{a,b,c}];
```

```wolfram
Solve[y^2+x z==1&&{x,y,z}∈ℛ,{y,z},Reals]
(* Output *)
{{y->(7)/(x Root[49-x^2 #1^2+x^3 #1^3&,1]),z->Root[49-x^2 #1^2+x^3 #1^3&,1]}}
```

A parametrically defined region:

```wolfram
ℛ=ParametricRegion[{s+t,s-t,s t},{s,t}];
```

```wolfram
Solve[x y==z&&x+2 y+3 z==1&&{x,y,z}∈ℛ,{x,y,z},Reals]
(* Output *)
{{x->(2 Sqrt[5])/((1-Sqrt[5]) (3+(2)/(1-Sqrt[5]))),y->(1)/(2) (1-Sqrt[5]),z->(Sqrt[5])/(3+(2)/(1-Sqrt[5]))},{x->(6 (1+(1)/(3) (-3+Sqrt[5])))/((3-Sqrt[5]) (3+(6)/(3-Sqrt[5]))),y->(1)/(6) (3-Sqrt[5]),z->(1+(1)/(3) (-3+Sqrt[5]))/(3+(6)/(3-Sqrt[5]))},{x->-(2 Sqrt[5])/((1+Sqrt[5]) (3+(2)/(1+Sqrt[5]))),y->(1)/(2) (1+Sqrt[5]),z->-(Sqrt[5])/(3+(2)/(1+Sqrt[5]))},{x->(6 (1+(1)/(3) (-3-Sqrt[5])))/((3+Sqrt[5]) (3+(6)/(3+Sqrt[5]))),y->(1)/(6) (3+Sqrt[5]),z->(1+(1)/(3) (-3-Sqrt[5]))/(3+(6)/(3+Sqrt[5]))}}
```

Derived regions:

```wolfram
ℛ_1=Disk[{0,0},2];
ℛ_2=Circle[{1,1},2];
ℛ_3=RegionIntersection[ℛ_1,ℛ_2];
```

```wolfram
Solve[x^2==x y+1&&{x,y}∈ℛ_3,{x,y},Reals]
(* Output *)
{{x->Root,y->1-Sqrt[3+2 Root-Root^2]},{x->Root,y->1-Sqrt[3+2 Root-Root^2]}}
```

Plot it:

```wolfram
Show[{ContourPlot[x^2==x y+1,{x,-2,3},{y,-2,3}],Graphics[{{Opacity[0.5],Yellow,ℛ_1},{Green,ℛ_2},{Red,Point[{x,y}/.%]}}]}]
```

*([Graphics])*

Eliminate quantifiers over a Cartesian product of regions:

```wolfram
ℛ=RegionProduct[Circle[],Circle[]];
```

```wolfram
Solve[∃_{a,b}(4 x y a b==1&&{x,a,y,b}∈ℛ),{x,y},Reals]
(* Output *)
{{x->-(1)/(Sqrt[2]),y->-(1)/(Sqrt[2])},{x->-(1)/(Sqrt[2]),y->(1)/(Sqrt[2])},{x->(1)/(Sqrt[2]),y->-(1)/(Sqrt[2])},{x->(1)/(Sqrt[2]),y->(1)/(Sqrt[2])}}
```

Regions dependent on parameters:

```wolfram
ℛ_1=InfiniteLine[{{2,0},{0,t}}];
ℛ_2=Circle[];
```

The answer depends on the parameter value $t$:

```wolfram
Solve[{x,y}∈ℛ_1&&{x,y}∈ℛ_2,{x,y},Reals]
(* Output *)
{{x->(2 t-2 ((4 t)/(4+t^2)-Sqrt[(4 t^2-3 t^4)/((4+t^2)^2)]))/(t),y->(4 t)/(4+t^2)-Sqrt[(4 t^2-3 t^4)/((4+t^2)^2)]},{x->(2 t-2 ((4 t)/(4+t^2)+Sqrt[(4 t^2-3 t^4)/((4+t^2)^2)]))/(t),y->(4 t)/(4+t^2)+Sqrt[(4 t^2-3 t^4)/((4+t^2)^2)]}}
```

Use $x \in \mathcal{R}$ to specify that $x$ is a vector in $\mathbb{R}^{2}$:

```wolfram
ℛ=RegionIntersection[Circle[],Line[{{-2,-1},{1,2}}]];
```

```wolfram
Solve[x∈ℛ,x]
(* Output *)
{{x->{-1,0}},{x->{0,1}}}
```

In this case $x$ is a vector in $\mathbb{R}^{3}$:

```wolfram
ℛ=Sphere[];
```

```wolfram
Solve[x.{1,2,3}==0&&x.{-3,-2,-1}==0&&x∈ℛ,x]
(* Output *)
{{x->{2 Sqrt[(2)/(3)]-Sqrt[(3)/(2)],-Sqrt[(2)/(3)],(1)/(Sqrt[6])}},{x->{-2 Sqrt[(2)/(3)]+Sqrt[(3)/(2)],Sqrt[(2)/(3)],-(1)/(Sqrt[6])}}}
```

### Generalizations & Extensions

All variables are solved for:

```wolfram
Solve[x^2+3x+1==0]
(* Output *)
{{x->(1)/(2) (-3-Sqrt[5])},{x->(1)/(2) (-3+Sqrt[5])}}
```

```wolfram
Solve[{x+y==2,y-x==1}]
(* Output *)
{{x->(1)/(2),y->(3)/(2)}}
```

### Options

#### Assumptions

Specify conditions on parameters using [Assumptions](https://reference.wolfram.com/language/ref/Assumptions.html):

```wolfram
Solve[x^2==a,x,Reals,Assumptions->a>0]
(* Output *)
{{x->-Sqrt[a]},{x->Sqrt[a]}}
```

```wolfram
Solve[x^2==a,x,Reals,Assumptions->a<0]
(* Output *)
{}
```

By default, no solutions that require parameters to satisfy equations are produced:

```wolfram
Solve[s==(2t)/(1+t^2)&&c==(1-t^2)/(1+t^2),t]
(* Output *)
{}
```

With an equation on parameters given as an assumption, a solution is returned:

```wolfram
Solve[s==(2t)/(1+t^2)&&c==(1-t^2)/(1+t^2),t,Assumptions->s^2+c^2==1]
(* Output *)
{{t->(s)/(1+c)}}
```

Assumptions that contain solve variables are considered to be a part of the system to solve:

```wolfram
Solve[2^x==8,x,Assumptions->x>0]
(* Output *)
{{x->3}}
```

Equivalent statement without using [Assumptions](https://reference.wolfram.com/language/ref/Assumptions.html):

```wolfram
Solve[2^x==8&&x>0,x]
(* Output *)
{{x->3}}
```

With parameters assumed to belong to a discrete set, solutions involving arbitrary conditions are returned:

```wolfram
Solve[a x^2+b x+c==0,x,Assumptions->Element[a,Integers]]
(* Output *)
Solve
(* Output *)
{{x->-(c)/(b)},{x->(-b-Sqrt[b^2-4 a c])/(2 a)},{x->(-b+Sqrt[b^2-4 a c])/(2 a)}}
```

#### Cubics

By default, [Solve](https://reference.wolfram.com/language/ref/Solve.html) uses general formulas for solving cubics in radicals only when symbolic parameters are present:

```wolfram
Solve[x^3+a x^2+2 x+3==0,x]
(* Output *)
{{x->-(a)/(3)-(2^(1/3) (6-a^2))/(3 (-81+18 a-2 a^3+3 Sqrt[3] Sqrt[275-108 a-4 a^2+12 a^3])^(1/3))+((-81+18 a-2 a^3+3 Sqrt[3] Sqrt[275-108 a-4 a^2+12 a^3])^(1/3))/(3 2^(1/3))},{x->-(a)/(3)+((1+ⅈ Sqrt[3]) (6-a^2))/(3 2^(2/3) (-81+18 a-2 a^3+3 Sqrt[3] Sqrt[275-108 a-4 a^2+12 a^3])^(1/3))-((1-ⅈ Sqrt[3]) (-81+18 a-2 a^3+3 Sqrt[3] Sqrt[275-108 a-4 a^2+12 a^3])^(1/3))/(6 2^(1/3))},{x->-(a)/(3)+((1-ⅈ Sqrt[3]) (6-a^2))/(3 2^(2/3) (-81+18 a-2 a^3+3 Sqrt[3] Sqrt[275-108 a-4 a^2+12 a^3])^(1/3))-((1+ⅈ Sqrt[3]) (-81+18 a-2 a^3+3 Sqrt[3] Sqrt[275-108 a-4 a^2+12 a^3])^(1/3))/(6 2^(1/3))}}
```

For polynomials with numeric coefficients, [Solve](https://reference.wolfram.com/language/ref/Solve.html) does not use the formulas:

```wolfram
Solve[x^3+2 x^2+3 x+4==0,x]
(* Output *)
{{x->Root},{x->Root},{x->Root}}
```

With [Cubics](https://reference.wolfram.com/language/ref/Cubics.html)->[False](https://reference.wolfram.com/language/ref/False.html), [Solve](https://reference.wolfram.com/language/ref/Solve.html) never uses the formulas:

```wolfram
Solve[x^3+a x^2+2 x+3==0,x,Cubics->False]
(* Output *)
{{x->Root[3+2 #1+a #1^2+#1^3&,1]},{x->Root[3+2 #1+a #1^2+#1^3&,2]},{x->Root[3+2 #1+a #1^2+#1^3&,3]}}
```

With [Cubics](https://reference.wolfram.com/language/ref/Cubics.html)->[True](https://reference.wolfram.com/language/ref/True.html), [Solve](https://reference.wolfram.com/language/ref/Solve.html) always uses the formulas:

```wolfram
Solve[x^3+2 x^2+3 x+4==0,x,Cubics->True]
(* Output *)
{{x->(1)/(3) (-2-(5^(2/3))/((-7+3 Sqrt[6])^(1/3))+(5 (-7+3 Sqrt[6]))^(1/3))},{x->-(2)/(3)+(5^(2/3) (1+ⅈ Sqrt[3]))/(6 (-7+3 Sqrt[6])^(1/3))-(1)/(6) (1-ⅈ Sqrt[3]) (5 (-7+3 Sqrt[6]))^(1/3)},{x->-(2)/(3)+(5^(2/3) (1-ⅈ Sqrt[3]))/(6 (-7+3 Sqrt[6])^(1/3))-(1)/(6) (1+ⅈ Sqrt[3]) (5 (-7+3 Sqrt[6]))^(1/3)}}
```

Real roots of irreducible cubics still contain [I](https://reference.wolfram.com/language/ref/I.html) in their algebraic forms by [casus irreducibilis](https://mathworld.wolfram.com/CasusIrreducibilis.html):

```wolfram
Solve[-61+110x-60 x^2+10 x^3==0,x,Cubics->True]
(* Output *)
{{x->2+(((1)/(5) (9+ⅈ Sqrt[1119]))^(1/3))/(6^(2/3))+(2^(2/3))/(((3)/(5) (9+ⅈ Sqrt[1119]))^(1/3))},{x->2-((1+ⅈ Sqrt[3]) ((1)/(5) (9+ⅈ Sqrt[1119]))^(1/3))/(2 6^(2/3))-(1-ⅈ Sqrt[3])/(((6)/(5) (9+ⅈ Sqrt[1119]))^(1/3))},{x->2-((1-ⅈ Sqrt[3]) ((1)/(5) (9+ⅈ Sqrt[1119]))^(1/3))/(2 6^(2/3))-(1+ⅈ Sqrt[3])/(((6)/(5) (9+ⅈ Sqrt[1119]))^(1/3))}}
```

Machine-precision numerical evaluation gives a spurious imaginary part:

```wolfram
N[%]
(* Output *)
{{x->3.0466805318046024+0. ⅈ},{x->1.8989687421189894+1.1102230246251565×10^-16 ⅈ},{x->1.0543507260764087-5.551115123125783×10^-17 ⅈ}}
```

Arbitrary-precision evaluation still leaves an imaginary part:

```wolfram
N[%%,20]
(* Output *)
{{x->3.04668053180460226115865631199604352874+0`19.666688080382638 ⅈ},{x->1.89896874211898918231164565878620696388+0`19.871997181718967 ⅈ},{x->1.05435072607640855652969802921774950737+0`20.1275298963836 ⅈ}}
```

With the default setting [Cubics](https://reference.wolfram.com/language/ref/Cubics.html)->[Automatic](https://reference.wolfram.com/language/ref/Automatic.html), explicitly real results are obtained:

```wolfram
Solve[-61+110x-60 x^2+10 x^3==0,x]
(* Output *)
{{x->Root},{x->Root},{x->Root}}
```

```wolfram
N[%]
(* Output *)
{{x->1.0543507260764085},{x->1.898968742118988},{x->3.0466805318046015}}
```

#### GeneratedParameters

[Solve](https://reference.wolfram.com/language/ref/Solve.html) may introduce new parameters to represent the solution:

```wolfram
Solve[x^2-3y^2==1&&x>0&&y>0,{x,y},Integers]
(* Output *)
{{x->(1)/(2) ((2-Sqrt[3])^1+(2+Sqrt[3])^1),y->-((2-Sqrt[3])^1-(2+Sqrt[3])^1)/(2 Sqrt[3])}}
```

Use [GeneratedParameters](https://reference.wolfram.com/language/ref/GeneratedParameters.html) to control how the parameters are generated:

```wolfram
Solve[x^2-3y^2==1&&x>0&&y>0,{x,y},Integers,GeneratedParameters->(k_#&)]
(* Output *)
{{x->(1)/(2) ((2-Sqrt[3])^(k_1)+(2+Sqrt[3])^(k_1)),y->-((2-Sqrt[3])^(k_1)-(2+Sqrt[3])^(k_1))/(2 Sqrt[3])}}
```

#### InverseFunctions

By default, [Solve](https://reference.wolfram.com/language/ref/Solve.html) uses inverse functions but prints warning messages:

```wolfram
Solve[f[x]==0,x]
(* Output *)
Solve
(* Output *)
{{x->f^((-1))[0]}}
```

```wolfram
Solve[x+E^x==1,x]
(* Output *)
Solve
(* Output *)
{{x->0}}
```

For symbols with the [NumericFunction](https://reference.wolfram.com/language/ref/NumericFunction.html) attribute, symbolic inverses are not used:

```wolfram
Solve[Gamma[x]==a,x]
(* Output *)
Solve
(* Output *)
Solve[Gamma[x]==a,x]
```

With [InverseFunctions](https://reference.wolfram.com/language/ref/InverseFunctions.html)->[True](https://reference.wolfram.com/language/ref/True.html), [Solve](https://reference.wolfram.com/language/ref/Solve.html) does not print inverse function warning messages:

```wolfram
Solve[f[x]==0,x,InverseFunctions->True]
(* Output *)
{{x->f^((-1))[0]}}
```

```wolfram
Solve[x+E^x==1,x,InverseFunctions->True]
(* Output *)
{{x->0}}
```

Symbolic inverses are used for all symbols:

```wolfram
Solve[Gamma[x]==a,x,InverseFunctions->True]
(* Output *)
{{x->Gamma^((-1))[a]}}
```

With [InverseFunctions](https://reference.wolfram.com/language/ref/InverseFunctions.html)->[False](https://reference.wolfram.com/language/ref/False.html), [Solve](https://reference.wolfram.com/language/ref/Solve.html) does not use inverse functions:

```wolfram
Solve[f[x]==0,x,InverseFunctions->False]
(* Output *)
Solve
(* Output *)
Solve[f[x]==0,x,InverseFunctions->False]
```

Solving algebraic equations does not require using inverse functions:

```wolfram
Solve[x^2+1==0,x,InverseFunctions->False]
(* Output *)
{{x->-ⅈ},{x->ⅈ}}
```

Here, a method based on [Reduce](https://reference.wolfram.com/language/ref/Reduce.html) is used, as it does not require using inverse functions:

```wolfram
Solve[x+E^x==1,x,InverseFunctions->False]
(* Output *)
{{x->0},{x->1-ProductLog[1,ℯ]}}
```

#### MaxExtraConditions

By default, no solutions requiring extra conditions are produced:

```wolfram
Solve[a==0&&x==0,x]
(* Output *)
{}
```

Unless the parameters are discrete:

```wolfram
Solve[a==0&&x==0,x,Integers]
(* Output *)
{{x->0}}
```

The default setting, [MaxExtraConditions](https://reference.wolfram.com/language/ref/MaxExtraConditions.html)->0, gives no solutions requiring conditions:

```wolfram
Solve[(x-a)(x-b)==0&&(x-c)(x-d)==0&& (x-a)(x-e)==0,x]
(* Output *)
{}
```

```wolfram
Solve[(x-a)(x-b)==0&&(x-c)(x-d)==0&& (x-a)(x-e)==0,x,MaxExtraConditions->0]
(* Output *)
{}
```

[MaxExtraConditions](https://reference.wolfram.com/language/ref/MaxExtraConditions.html)->1 gives solutions requiring up to one equation on parameters:

```wolfram
Solve[(x-a)(x-b)==0&&(x-c)(x-d)==0&& (x-a)(x-e)==0,x,MaxExtraConditions->1]
(* Output *)
{{x->c},{x->d}}
```

[MaxExtraConditions](https://reference.wolfram.com/language/ref/MaxExtraConditions.html)->2 gives solutions requiring up to two equations on parameters:

```wolfram
Solve[(x-a)(x-b)==0&&(x-c)(x-d)==0&& (x-a)(x-e)==0,x,MaxExtraConditions->2]
(* Output *)
{{x->c},{x->d},{x->e},{x->e}}
```

Give solutions requiring the minimal number of parameter equations:

```wolfram
Solve[(x-a)(x-b)==0&&(x-c)(x-d)==0&& (x-a)(x-e)==0,x,MaxExtraConditions->Automatic]
(* Output *)
{{x->c},{x->d}}
```

Give all solutions:

```wolfram
Solve[(x-a)(x-b)==0&&(x-c)(x-d)==0&& (x-a)(x-e)==0,x,MaxExtraConditions->All]
(* Output *)
{{x->c},{x->d},{x->e}}
```

By default, [Solve](https://reference.wolfram.com/language/ref/Solve.html) drops inequation conditions on continuous parameters:

```wolfram
Solve[a x==1,x]
(* Output *)
{{x->(1)/(a)}}
```

With [MaxExtraConditions](https://reference.wolfram.com/language/ref/MaxExtraConditions.html)->[All](https://reference.wolfram.com/language/ref/All.html), [Solve](https://reference.wolfram.com/language/ref/Solve.html) includes all conditions:

```wolfram
Solve[a x==1,x,MaxExtraConditions->All]
(* Output *)
{{x->(1)/(a)}}
```

#### MaxRoots

Find $3$ out of $12345$ roots of a polynomial:

```wolfram
Solve[x^12345+x+1==0,x,MaxRoots->3]
(* Output *)
{{x->Root},{x->Root},{x->Root}}
```

Find $3$ out of $1000000000$ roots of a polynomial system:

```wolfram
Solve[x^1000==y^3+1&&y^1000==z^3+1&&z^1000==x^3+1,{x,y,z},MaxRoots->3]
(* Output *)
{{x->Root,y->Root,z->Root},{x->Root,y->Root,z->Root},{x->Root,y->Root,z->Root}}
```

Find $5$ solutions of a transcendental system:

```wolfram
Solve[Sin[x+y]==2x+3y+1&&AiryAi[x^2-y]==y^2-1,{x,y},MaxRoots->5]
(* Output *)
{{x->Root,y->Root},{x->Root,y->Root},{x->Root,y->Root},{x->Root,y->Root},{x->Root,y->Root}}
```

When the system contains symbolic parameters, the option value is ignored:

```wolfram
Solve[x^4==a,x,MaxRoots->3]
(* Output *)
Solve
(* Output *)
{{x->-a^(1/4)},{x->-ⅈ a^(1/4)},{x->ⅈ a^(1/4)},{x->a^(1/4)}}
```

#### Method

By default, [Solve](https://reference.wolfram.com/language/ref/Solve.html) uses inverse functions to solve non-polynomial complex equations:

```wolfram
Solve[x E^x==1/2,x]
(* Output *)
Solve
(* Output *)
{{x->ProductLog[(1)/(2)]}}
```

With [Method](https://reference.wolfram.com/language/ref/Method.html)->[Reduce](https://reference.wolfram.com/language/ref/Reduce.html), [Solve](https://reference.wolfram.com/language/ref/Solve.html) uses [Reduce](https://reference.wolfram.com/language/ref/Reduce.html) to find the complete solution set:

```wolfram
Solve[x E^x==1/2,x,Method->Reduce]
(* Output *)
{{x->ProductLog[1,(1)/(2)]}}
```

#### Modulus

Solve equations over the integers modulo 9:

```wolfram
Solve[x^2+3y^2==4&&3x^3-4y^2+x y==1,{x,y},Modulus->9]
(* Output *)
{{x->8,y->1},{x->8,y->4},{x->8,y->7}}
```

Find a modulus for which a system of equations has a solution:

```wolfram
Solve[x^2+y^2==1&&2x+3y==5&&x^3-x y==7,{x,y},Modulus->Automatic]
(* Output *)
{{Modulus->191721,y->96231,x->47377}}
```

#### Quartics

By default, [Solve](https://reference.wolfram.com/language/ref/Solve.html) uses the general formulas for solving quartics in radicals only when symbolic parameters are present:

```wolfram
Solve[x^4+a x^2+2 x+3==0,x][[1]]
(* Output *)
{x->(1)/(2) Sqrt(-(2 a)/(3)+(36+a^2)/(3 (54-108 a+a^3+6 Sqrt[3] Sqrt[-405-108 a+72 a^2+a^3-3 a^4])^(1/3))+(1)/(3) (54-108 a+a^3+6 Sqrt[3] Sqrt[-405-108 a+72 a^2+a^3-3 a^4])^(1/3))-(1)/(2) Sqrt(-(4 a)/(3)-(36+a^2)/(3 (54-108 a+a^3+6 Sqrt[3] Sqrt[-405-108 a+72 a^2+a^3-3 a^4])^(1/3))-(1)/(3) (54-108 a+a^3+6 Sqrt[3] Sqrt[-405-108 a+72 a^2+a^3-3 a^4])^(1/3)-4/(Sqrt(-(2 a)/(3)+(36+a^2)/(3 (54-108 a+a^3+6 Sqrt[3] Sqrt[-405-108 a+72 a^2+a^3-3 a^4])^(1/3))+(1)/(3) (54-108 a+a^3+6 Sqrt[3] Sqrt[-405-108 a+72 a^2+a^3-3 a^4])^(1/3))))}
```

For polynomials with numeric coefficients, [Solve](https://reference.wolfram.com/language/ref/Solve.html) does not use the formulas:

```wolfram
Solve[x^4+2 x^2+3 x+4==0,x]
(* Output *)
{{x->Root},{x->Root},{x->Root},{x->Root}}
```

With [Quartics](https://reference.wolfram.com/language/ref/Quartics.html)->[False](https://reference.wolfram.com/language/ref/False.html), [Solve](https://reference.wolfram.com/language/ref/Solve.html) never uses the formulas:

```wolfram
Solve[x^4+a x^2+2 x+3==0,x,Quartics->False]
(* Output *)
{{x->Root[3+2 #1+a #1^2+#1^4&,1]},{x->Root[3+2 #1+a #1^2+#1^4&,2]},{x->Root[3+2 #1+a #1^2+#1^4&,3]},{x->Root[3+2 #1+a #1^2+#1^4&,4]}}
```

With [Quartics](https://reference.wolfram.com/language/ref/Quartics.html)->[True](https://reference.wolfram.com/language/ref/True.html), [Solve](https://reference.wolfram.com/language/ref/Solve.html) always uses the formulas:

```wolfram
Solve[x^4+2 x^2+3 x+4==0,x,Quartics->True][[1]]
(* Output *)
{x->(1)/(2) Sqrt[(1)/(3) (-4+(52)/(((1)/(2) (-317+9 ⅈ Sqrt[5703]))^(1/3))+((1)/(2) (-317+9 ⅈ Sqrt[5703]))^(1/3))]-(1)/(2) Sqrt(-(8)/(3)-(52)/(3 ((1)/(2) (-317+9 ⅈ Sqrt[5703]))^(1/3))-(1)/(3) ((1)/(2) (-317+9 ⅈ Sqrt[5703]))^(1/3)-(6)/(Sqrt[(1)/(3) (-4+(52)/(((1)/(2) (-317+9 ⅈ Sqrt[5703]))^(1/3))+((1)/(2) (-317+9 ⅈ Sqrt[5703]))^(1/3))]))}
```

#### VerifySolutions

[Solve](https://reference.wolfram.com/language/ref/Solve.html) verifies solutions obtained using non-equivalent transformations:

```wolfram
eqns=Sqrt[x+Sqrt[x]]==2&&y^9-y-2 x^(1/7)==3;
(sol1=Solve[eqns,{x,y}]);//Timing
(* Output *)
{5.296875,Null}
```

With [VerifySolutions](https://reference.wolfram.com/language/ref/VerifySolutions.html)->[False](https://reference.wolfram.com/language/ref/False.html), [Solve](https://reference.wolfram.com/language/ref/Solve.html) does not verify the solutions:

```wolfram
(sol2=Solve[eqns,{x,y},VerifySolutions->False]);//Timing
(* Output *)
{3.34375,Null}
```

Some of the solutions returned with [VerifySolutions](https://reference.wolfram.com/language/ref/VerifySolutions.html)->[False](https://reference.wolfram.com/language/ref/False.html) are not correct:

```wolfram
Length/@{sol1, sol2}
(* Output *)
{9,126}
```

This uses a fast numeric test in an attempt to select correct solutions:

```wolfram
(sol3=Select[sol2,TrueQ[eqns/.N[#,20]]&]);//Timing
(* Output *)
{0.078125,Null}
```

In this case numeric verification gives the correct solution set:

```wolfram
sol3===sol1
(* Output *)
True
```

#### WorkingPrecision

By default, [Solve](https://reference.wolfram.com/language/ref/Solve.html) finds exact solutions of equations:

```wolfram
SeedRandom[0];
mat=Table[RandomInteger[{-10^100,10^100}], {200},{200}];
b=Table[RandomInteger[{-10^100,10^100}],{200}];
vars=x/@Range[200];
(sol=Solve[Thread[mat.vars==b],vars]);//Timing
(* Output *)
{7.34375,Null}
```

Computing the solution using 100-digit numbers is faster:

```wolfram
(sol1=Solve[Thread[mat.vars==b],vars,WorkingPrecision->100]);//Timing
(* Output *)
{2.328125,Null}
```

The result agrees with the exact solution in the first 100 digits:

```wolfram
Max[Abs[(vars/.sol)-(vars/.sol1)]]
(* Output *)
0`96.5798549905083
```

Computing the solution using machine numbers is much faster:

```wolfram
(sol2=Solve[Thread[mat.vars==b],vars,WorkingPrecision->MachinePrecision]);//Timing
(* Output *)
{0.09375,Null}
```

The result is still quite close to the exact solution:

```wolfram
Max[Abs[(vars/.sol)-(vars/.sol2)]]
(* Output *)
4.1300296516055823×10^-13
```

### Applications

Solve a quadratic equation:

```wolfram
Solve[a x^2+b x+c==0,x]
(* Output *)
{{x->(-b-Sqrt[b^2-4 a c])/(2 a)},{x->(-b+Sqrt[b^2-4 a c])/(2 a)}}
```

Find intersection points of a circle and a parabola:

```wolfram
pts=Solve[x^2+y^2==1&&y-2x^2+3/2==0,{x,y}]
(* Output *)
{{x->-(1)/(2) Sqrt[(1)/(2) (5-Sqrt[5])],y->(1)/(4) (-1-Sqrt[5])},{x->(1)/(2) Sqrt[(1)/(2) (5-Sqrt[5])],y->(1)/(4) (-1-Sqrt[5])},{x->-(1)/(2) Sqrt[(1)/(2) (5+Sqrt[5])],y->(1)/(4) (-1+Sqrt[5])},{x->(1)/(2) Sqrt[(1)/(2) (5+Sqrt[5])],y->(1)/(4) (-1+Sqrt[5])}}
```

```wolfram
Show[{ContourPlot[{x^2+y^2==1,y-2x^2+3/2==0},{x,-1.5,1.5},{y,-1.5,1.5}],Graphics[{Red,PointSize[Medium],Point[{x,y}/.pts]}]}]
```

*([Graphics])*

Find conditions for a quartic to have all roots equal:

```wolfram
f[x_]:= x^4 + a x^3 + b x^2 + c x + d
```

A method using [Subresultants](https://reference.wolfram.com/language/ref/Subresultants.html):

```wolfram
Solve[Thread[Drop[Subresultants[f[x],D[f[x],x],x],-1]==0],{b,c,d}]
(* Output *)
{{b->(3 a^2)/(8),c->(a^3)/(16),d->(a^4)/(256)}}
```

A method using quantifier elimination:

```wolfram
Solve[Exists[x,f[x]==0,ForAll[y,f[y]==0, x==y]],{b,c,d}]
(* Output *)
{{b->(3 a^2)/(8),c->(a^3)/(16),d->(a^4)/(256)}}
```

Plot a space curve given by an implicit description:

```wolfram
curve=x^2+y^2+z^2==1&&x^3+x y^2==z^2;
```

```wolfram
sol=Solve[curve,{y,z},Reals]
(* Output *)
{{y->-Sqrt[(1-x^2-x^3)/(1+x)],z->-Sqrt[x^3+(x (1-x^2-x^3))/(1+x)]},{y->-Sqrt[(1-x^2-x^3)/(1+x)],z->Sqrt[x^3+(x (1-x^2-x^3))/(1+x)]},{y->Sqrt[(1-x^2-x^3)/(1+x)],z->-Sqrt[x^3+(x (1-x^2-x^3))/(1+x)]},{y->Sqrt[(1-x^2-x^3)/(1+x)],z->Sqrt[x^3+(x (1-x^2-x^3))/(1+x)]}}
```

```wolfram
cond=Union[#[[2]]&/@Cases[sol, _ConditionalExpression, Infinity]]
(* Output *)
{0<x<Root}
```

```wolfram
{l,u}=N[{cond[[1,1]],cond[[1,5]]}]
(* Output *)
{0.,0.7548776662466927}
```

```wolfram
ParametricPlot3D[{x,y,z}/.sol,{x,l,u}]
```

*([Graphics3D])*

Plot the projection of the space curve on the `{*x*,*y*}` plane:

```wolfram
proj=Solve[Exists[z,curve],y,Reals]
(* Output *)
{{y->-Sqrt[(1-x^2-x^3)/(1+x)]},{y->Sqrt[(1-x^2-x^3)/(1+x)]}}
```

```wolfram
ParametricPlot[{x,y}/.proj,{x,l,u}]
```

*([Graphics])*

Find a Pythagorean triple:

```wolfram
Solve[x^2+y^2==5^2&&y>x>0,{x,y},Integers]
(* Output *)
{{x->3,y->4}}
```

Find a sequence of Pythagorean triples:

```wolfram
Table[Solve[x^2+y^2==z^2&&y>x>0,{x,y},Integers],{z,30}]
(* Output *)
{{},{},{},{},{{x->3,y->4}},{},{},{},{},{{x->6,y->8}},{},{},{{x->5,y->12}},{},{{x->9,y->12}},{},{{x->8,y->15}},{},{},{{x->12,y->16}},{},{},{},{},{{x->7,y->24},{x->15,y->20}},{{x->10,y->24}},{},{},{{x->20,y->21}},{{x->18,y->24}}}
```

Find how to pay $2.27 postage with 10-, 23-, and 37-cent stamps:

```wolfram
Solve[a 10 + b 23 + c 37 == 227 &&a>=0&&b>=0&&c>=0,{a,b,c},Integers]
(* Output *)
{{a->1,b->3,c->4},{a->2,b->9,c->0},{a->7,b->2,c->3},{a->13,b->1,c->2},{a->19,b->0,c->1}}
```

The same task can be accomplished with [IntegerPartitions](https://reference.wolfram.com/language/ref/IntegerPartitions.html):

```wolfram
IntegerPartitions[227,All,{37,23,10}]
(* Output *)
{{10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,37},{10,10,10,10,10,10,10,10,10,10,10,10,10,23,37,37},{10,10,10,10,10,10,10,23,23,37,37,37},{10,10,23,23,23,23,23,23,23,23,23},{10,23,23,23,37,37,37,37}}
```

Find 200 roots of a complex analytic function:

```wolfram
rts=Solve[Sin[x^3]-x^5==1,x,MaxRoots->200];
```

Show the roots on the complex plot for the function:

```wolfram
ComplexPlot[Sin[x^3]-x^5-1,{x,-7.5-7.5I,7.5+7.5I},Epilog -> LightYellowPoint[ReplaceAll[ReIm[x], rts]], ColorFunction -> None]
```

![image](img/image_009.png)

### Properties & Relations

Solutions satisfy the equations:

```wolfram
eqns={x^2+y^2==1&&2x+3y==4};
Solve[eqns,{x,y}]
(* Output *)
{{x->(1)/(13) (8-3 ⅈ Sqrt[3]),y->(2)/(13) (6+ⅈ Sqrt[3])},{x->(1)/(13) (8+3 ⅈ Sqrt[3]),y->(2)/(13) (6-ⅈ Sqrt[3])}}
```

Solutions are given as replacement rules and can be directly used for substitution:

```wolfram
eqns/.%
(* Output *)
{{True},{True}}
```

[Solve](https://reference.wolfram.com/language/ref/Solve.html) uses `{}` to represent the empty solution or no solution:

```wolfram
Solve[x==1&&x==2,x]
(* Output *)
{}
```

[Solve](https://reference.wolfram.com/language/ref/Solve.html) uses `{{}}` to represent the universal solution or all points satisfying the equations:

```wolfram
Solve[x==x,x]
(* Output *)
{{}}
```

For univariate equations, [Solve](https://reference.wolfram.com/language/ref/Solve.html) repeats solutions according to their multiplicity:

```wolfram
Solve[(x-1)^2(x-2)^3==0,x]
(* Output *)
{{x->1},{x->1},{x->2},{x->2},{x->2}}
```

Solutions of algebraic equations are often given in terms of [Root](https://reference.wolfram.com/language/ref/Root.html) objects:

```wolfram
Solve[x^5-2x+7==0,x]
(* Output *)
{{x->Root},{x->Root},{x->Root},{x->Root},{x->Root}}
```

Use [N](https://reference.wolfram.com/language/ref/N.html) to compute numeric approximations of [Root](https://reference.wolfram.com/language/ref/Root.html) objects:

```wolfram
N[%,20]
(* Output *)
{{x->-1.590595373196720625937165890567850784546490647021529838258},{x->-0.364106049434746317144012958059943726025419824823305155234-1.482002193950815920567877380381905672897040661675151447088 ⅈ},{x->-0.364106049434746317144012958059943726025419824823305155234+1.482002193950815920567877380381905672897040661675151447088 ⅈ},{x->1.159403736033106630112595903343869654869510169795992678404-0.738550307091452775151024140797594919566349478309649321681 ⅈ},{x->1.159403736033106630112595903343869654869510169795992678404+0.738550307091452775151024140797594919566349478309649321681 ⅈ}}
```

[Root](https://reference.wolfram.com/language/ref/Root.html) objects may involve parameters:

```wolfram
Solve[x^5+a x+1==0,x]
(* Output *)
{{x->Root[1+a #1+#1^5&,1]},{x->Root[1+a #1+#1^5&,2]},{x->Root[1+a #1+#1^5&,3]},{x->Root[1+a #1+#1^5&,4]},{x->Root[1+a #1+#1^5&,5]}}
```

Use [Series](https://reference.wolfram.com/language/ref/Series.html) to compute series expansions of [Root](https://reference.wolfram.com/language/ref/Root.html) objects:

```wolfram
Series[x/.%[[1]],{a,0,10}]
(* Output *)
-1+(a)/(5)+(a^2)/(25)+(a^3)/(125)-(21 a^5)/(15625)-(78 a^6)/(78125)-(187 a^7)/(390625)-(286 a^8)/(1953125)+(9367 a^10)/(244140625)+O[a]^11
```

The series satisfies the equation up to order 11:

```wolfram
x^5+a x+1/.x->%
(* Output *)
O[a]^11
```

[Solve](https://reference.wolfram.com/language/ref/Solve.html) represents solutions in terms of replacement rules:

```wolfram
Solve[x^4==4,x]
(* Output *)
{{x->-Sqrt[2]},{x->-ⅈ Sqrt[2]},{x->ⅈ Sqrt[2]},{x->Sqrt[2]}}
```

[Reduce](https://reference.wolfram.com/language/ref/Reduce.html) represents solutions in terms of Boolean combinations of equations and inequalities:

```wolfram
Reduce[x^4==4,x]
(* Output *)
x==-Sqrt[2]||x==-ⅈ Sqrt[2]||x==ⅈ Sqrt[2]||x==Sqrt[2]
```

[Solve](https://reference.wolfram.com/language/ref/Solve.html) uses fast heuristics to solve transcendental equations, but may give incomplete solutions:

```wolfram
Solve[x-E^x==1, x]
(* Output *)
Solve
(* Output *)
{{x->1-ProductLog[-ℯ]}}
```

```wolfram
Solve[x Log[x]==a,x]
(* Output *)
Solve
(* Output *)
{{x->(a)/(ProductLog[a])}}
```

[Reduce](https://reference.wolfram.com/language/ref/Reduce.html) uses methods that are often slower, but finds all solutions and gives all necessary conditions:

```wolfram
Reduce[x-E^x==1,x]
(* Output *)
1∈Integers&&x==1-ProductLog[1,-ℯ]
```

```wolfram
Reduce[x Log[x]==a,x]
(* Output *)
(a≠0&&((Im[ProductLog[-1,a]]>-π&&x==ℯ^(ProductLog[-1,a]))||x==ℯ^ProductLog[a]||(Im[ProductLog[1,a]]<=π&&x==ℯ^ProductLog[1,a])))||(a==0&&x==1)
```

Use [FindInstance](https://reference.wolfram.com/language/ref/FindInstance.html) to find solution instances:

```wolfram
FindInstance[x^2+y^2+z^2==1&&2x y==z^3,{x,y,z}]
(* Output *)
{{x->-2,y->AlgebraicNumber,z->AlgebraicNumber}}
```

Like [Reduce](https://reference.wolfram.com/language/ref/Reduce.html), [FindInstance](https://reference.wolfram.com/language/ref/FindInstance.html) can be given inequalities and domain specifications:

```wolfram
FindInstance[x^2-3y^2==1&&0<x<10^10,{x,y},Integers,3]
(* Output *)
{{x->1,y->0},{x->18817,y->-10864},{x->2,y->-1}}
```

Use [DSolve](https://reference.wolfram.com/language/ref/DSolve.html) to solve differential equations:

```wolfram
DSolve[y''[x]==-y[x],y[x],x]
(* Output *)
{{y[x]->1 Cos[x]+2 Sin[x]}}
```

```wolfram
DSolve[{y''[x]==-y[x],y[0]==0,y'[0]==1},y[x],x]
(* Output *)
{{y[x]->Sin[x]}}
```

Use [RSolve](https://reference.wolfram.com/language/ref/RSolve.html) to solve recurrence equations:

```wolfram
RSolve[f[n+1]==n f[n],f[n],n]
(* Output *)
{{f[n]->1 Pochhammer[1,-1+n]}}
```

```wolfram
RSolve[{f[n+1]==n+f[n],f[0]==0},f[n],n]
(* Output *)
{{f[n]->(1)/(2) (-1+n) n}}
```

[SolveAlways](https://reference.wolfram.com/language/ref/SolveAlways.html) gives the values of parameters for which complex equations are always true:

```wolfram
SolveAlways[(a -2b+1)x^2+(a-b^2-c)x ==a^2-b+3c-1,x]
(* Output *)
{{a->-2-Sqrt[13],b->(1)/(2) (-1-Sqrt[13]),c->(1)/(2) (-11-3 Sqrt[13])},{a->-2+Sqrt[13],b->(1)/(2) (-1+Sqrt[13]),c->-(11)/(2)+(3 Sqrt[13])/(2)}}
```

The same problem can be expressed using [ForAll](https://reference.wolfram.com/language/ref/ForAll.html) and solved with [Solve](https://reference.wolfram.com/language/ref/Solve.html) or [Reduce](https://reference.wolfram.com/language/ref/Reduce.html):

```wolfram
Solve[ForAll[x,(a -2b+1)x^2+(a-b^2-c)x ==a^2-b+3c-1],{a,b,c}]
(* Output *)
{{a->-2-Sqrt[13],b->(1)/(2) (-1-Sqrt[13]),c->(1)/(2) (-11-3 Sqrt[13])},{a->-2+Sqrt[13],b->(1)/(2) (-1+Sqrt[13]),c->(1)/(2) (-11+3 Sqrt[13])}}
```

```wolfram
Reduce[ForAll[x,(a -2b+1)x^2+(a-b^2-c)x ==a^2-b+3c-1],{a,b,c}]
(* Output *)
(a==-2-Sqrt[13]||a==-2+Sqrt[13])&&b==(1+a)/(2)&&c==(1)/(2) (-5+3 a)
```

[Resolve](https://reference.wolfram.com/language/ref/Resolve.html) eliminates quantifiers, possibly without solving the resulting quantifier-free system:

```wolfram
Resolve[Exists[x,x^2+y^2+z^2==1&&x y>z^3],Reals]
(* Output *)
(z<0&&y^2+z^2<=1)||(y^2+z^2<=1&&-y^2+y^4+y^2 z^2+z^6<0)
```

```wolfram
Resolve[ForAll[{x,y},a x^2+b x+c==0&&a y^2+b y+c==0,x==y]]
(* Output *)
(a==0&&b≠0)||(a==0&&c≠0)||(a≠0&&b^2-4 a c==0)
```

[Eliminate](https://reference.wolfram.com/language/ref/Eliminate.html) eliminates variables from systems of complex equations:

```wolfram
Eliminate[x^2+y^2+z^2==1&&x z==y^3,x]
(* Output *)
(-1+y^2) z^2+z^4==-y^6
```

This solves the same problem using [Resolve](https://reference.wolfram.com/language/ref/Resolve.html):

```wolfram
Resolve[Exists[x,x^2+y^2+z^2==1&&x z==y^3]]
(* Output *)
(y==0&&z==0)||(z≠0&&y^6-z^2+y^2 z^2+z^4==0)
```

[Reduce](https://reference.wolfram.com/language/ref/Reduce.html) and [Solve](https://reference.wolfram.com/language/ref/Solve.html) additionally solve the resulting equations:

```wolfram
Reduce[Exists[x,x^2+y^2+z^2==1&&x z==y^3],{y,z}]
(* Output *)
(y==0&&z==0)||((z==-(Sqrt[1-y^2-Sqrt[1-2 y^2+y^4-4 y^6]])/(Sqrt[2])||z==(Sqrt[1-y^2-Sqrt[1-2 y^2+y^4-4 y^6]])/(Sqrt[2])||z==-(Sqrt[1-y^2+Sqrt[1-2 y^2+y^4-4 y^6]])/(Sqrt[2])||z==(Sqrt[1-y^2+Sqrt[1-2 y^2+y^4-4 y^6]])/(Sqrt[2]))&&z≠0)
```

```wolfram
Solve[Exists[x,x^2+y^2+z^2==1&&x z==y^3],z]
(* Output *)
{{z->-(Sqrt[1-y^2-Sqrt[1-2 y^2+y^4-4 y^6]])/(Sqrt[2])},{z->(Sqrt[1-y^2-Sqrt[1-2 y^2+y^4-4 y^6]])/(Sqrt[2])},{z->-(Sqrt[1-y^2+Sqrt[1-2 y^2+y^4-4 y^6]])/(Sqrt[2])},{z->(Sqrt[1-y^2+Sqrt[1-2 y^2+y^4-4 y^6]])/(Sqrt[2])}}
```

$f(x)$ is bijective iff the equation $f(x)=y$ has exactly one solution for each $y$:

```wolfram
Solve[x^3+x+1==y,x,Reals]
(* Output *)
{{x->Root[1-y+#1+#1^3&,1]}}
```

```wolfram
Solve[x^2+x+1==y,x,Reals]
(* Output *)
{{x->-(1)/(2)-(1)/(2) Sqrt[-3+4 y]},{x->-(1)/(2)+(1)/(2) Sqrt[-3+4 y]}}
```

Use [FunctionBijective](https://reference.wolfram.com/language/ref/FunctionBijective.html) to test whether a function is bijective:

```wolfram
FunctionBijective[x^3+x+1,x]
(* Output *)
True
```

```wolfram
FunctionBijective[x^2+x+1,x]
(* Output *)
False
```

Use [FunctionAnalytic](https://reference.wolfram.com/language/ref/FunctionAnalytic.html) to test whether a function is analytic:

```wolfram
f=Sin[3x^3]-2x+1;
```

```wolfram
FunctionAnalytic[f,x,Complexes]
(* Output *)
True
```

An analytic function can have only finitely many zeros in a closed and bounded region:

```wolfram
Solve[f==0&&Abs[x]<=1,x]
(* Output *)
{{x->Root},{x->Root},{x->Root},{x->Root},{x->Root}}
```

```wolfram
ComplexPlot[f,{x,-1-I,1+I},RegionFunction->(Abs[#]<=1&),Epilog->{PointSize[Medium],Blue,Point[{Re[x],Im[x]}]/.%}]
```

![image](img/image_010.png)

### Possible Issues

[Solve](https://reference.wolfram.com/language/ref/Solve.html) gives generic solutions; solutions involving equations on parameters are not given:

```wolfram
Solve[a x^2+x==1,x]
(* Output *)
{{x->(-1-Sqrt[1+4 a])/(2 a)},{x->(-1+Sqrt[1+4 a])/(2 a)}}
```

[Reduce](https://reference.wolfram.com/language/ref/Reduce.html) gives all solutions, including those that require equations on parameters:

```wolfram
Reduce[a x^2+x==1,x]
(* Output *)
(a==0&&x==1)||(a≠0&&(x==(-1-Sqrt[1+4 a])/(2 a)||x==(-1+Sqrt[1+4 a])/(2 a)))
```

With [MaxExtraConditions](https://reference.wolfram.com/language/ref/MaxExtraConditions.html)->[All](https://reference.wolfram.com/language/ref/All.html), [Solve](https://reference.wolfram.com/language/ref/Solve.html) also gives non-generic solutions:

```wolfram
Solve[a x^2+x==1,x,MaxExtraConditions->All]
(* Output *)
{{x->1},{x->(-1-Sqrt[1+4 a])/(2 a)},{x->(-1+Sqrt[1+4 a])/(2 a)}}
```

[Solve](https://reference.wolfram.com/language/ref/Solve.html) results do not depend on whether some of the input equations contain only parameters. The following two systems are equivalent and have no generic solutions:

```wolfram
Solve[x==1&& a==2,x]
(* Output *)
{}
```

```wolfram
Solve[x==1&&x+a==x+2,x]
(* Output *)
{}
```

Use [MaxExtraConditions](https://reference.wolfram.com/language/ref/MaxExtraConditions.html) to specify the number of parameter conditions allowed:

```wolfram
Solve[x==1&& a==2,x,MaxExtraConditions->1]
(* Output *)
{{x->1}}
```

```wolfram
Solve[x==1&& x+a==x+2,x,MaxExtraConditions->1]
(* Output *)
{{x->1}}
```

Use the [Exists](https://reference.wolfram.com/language/ref/Exists.html) quantifier to find solutions that are valid for some value of parameter $a$:

```wolfram
Solve[Exists[a,x==1&& a==2],x]
(* Output *)
{{x->1}}
```

```wolfram
Solve[Exists[a,x==1&& x+a==x+2],x]
(* Output *)
{{x->1}}
```

[Solve](https://reference.wolfram.com/language/ref/Solve.html) does not eliminate solutions that are neither generically correct nor generically incorrect:

```wolfram
Solve[Sqrt[((1-z)/(1+z))^2]==(1-a)/(1+a),z]
(* Output *)
{{z->(1)/(a)},{z->a}}
```

The solutions are correct for $0<|a|<1$ and incorrect for $|a|>1$:

```wolfram
Sqrt[((1-z)/(1+z))^2]-(1-a)/(1+a)/.%/.{{a->(2+3I)/4},{a->2+3I}}
(* Output *)
{{0,0},{(4)/(3)+(2 ⅈ)/(3),(4)/(3)+(2 ⅈ)/(3)}}
```

For transcendental equations, [Solve](https://reference.wolfram.com/language/ref/Solve.html) may not give all solutions:

```wolfram
Solve[x+E^x==1/2,x]
(* Output *)
Solve
(* Output *)
{{x->(1)/(2) (1-2 ProductLog[Sqrt[ℯ]])}}
```

Use [Reduce](https://reference.wolfram.com/language/ref/Reduce.html) to get all solutions:

```wolfram
Reduce[x+E^x==1/2,x]
(* Output *)
1∈Integers&&x==(1)/(2)-ProductLog[1,Sqrt[ℯ]]
```

[Solve](https://reference.wolfram.com/language/ref/Solve.html) with [Method](https://reference.wolfram.com/language/ref/Method.html)->"Reduce" uses [Reduce](https://reference.wolfram.com/language/ref/Reduce.html) to find solutions, but returns replacement rules:

```wolfram
Solve[x+E^x==1/2,x,Method->"Reduce"]
(* Output *)
{{x->(1)/(2)-ProductLog[1,Sqrt[ℯ]]}}
```

Using inverse functions allows [Solve](https://reference.wolfram.com/language/ref/Solve.html) to find some solutions fast:

```wolfram
Solve[x^n==1,x]//Timing
(* Output *)
Solve
(* Output *)
{0.03125,{{x->1}}}
```

Finding the complete solution may take much longer, and the solution may be large:

```wolfram
(red=Reduce[x^n==1,x])//LeafCount//Timing
(* Output *)
{4.5,611}
```

This finds the values of `n` for which `x==2 is a solution:

```wolfram
Reduce[red/.x->2,n]
(* Output *)
(1∈Integers&&n==0)||(1∈Integers&&1>=1&&n==(2 ⅈ π 1)/(Log[2]))||(1∈Integers&&1<=-1&&n==(2 ⅈ π 1)/(Log[2]))
```

```wolfram
Simplify[x^n/.{x->2,n->2I Pi 1/Log[2]},Element[1,Integers]]
(* Output *)
1
```

Interpretation of assumptions depends on their syntactic properties. Here the solution is generic in the parameter space restricted by the assumptions:

```wolfram
Solve[x==a&&x==b,x, Assumptions->a==b]
(* Output *)
{{x->b}}
```

This mathematically equivalent assumption contains the solve variable, and hence is treated as a part of the system to solve:

```wolfram
Solve[x==a&&x==b,x, Assumptions->a+x==b+x]
(* Output *)
{}
```

There are no generic solutions, because the input is interpreted as:

```wolfram
Solve[x==a &&x==b&&a+x==b+x,x]
(* Output *)
{}
```

The solution is non-generic, since it requires the parameters to satisfy an equation:

```wolfram
Solve[x==a &&x==b&&a+x==b+x,x, MaxExtraConditions->All]
(* Output *)
{{x->b}}
```

When parameters are restricted to a discrete set, the notion of genericity is not well defined, and all solutions are returned:

```wolfram
Solve[a x^2+b x==1,x,Assumptions->Element[a|b,Integers]]
(* Output *)
{{x->(1)/(b)},{x->(-b-Sqrt[4 a+b^2])/(2 a)},{x->(-b+Sqrt[4 a+b^2])/(2 a)}}
```

Removable singularities of input equations are generally not considered valid solutions:

```wolfram
Solve[(x^2-2Sqrt[2]x+2)/(x-Sqrt[2])==0,x]
(* Output *)
{}
```

```wolfram
Limit[(x^2-2Sqrt[2]x+2)/(x-Sqrt[2]),x->Sqrt[2]]
(* Output *)
0
```

However, solutions may include removable singularities that are cancelled by automatic simplification:

```wolfram
Solve[x^2/x==0,x]
(* Output *)
{{x->0}}
```

The removable singularity at $x=0$ is cancelled by evaluation:

```wolfram
x^2/x==0
(* Output *)
x==0
```

Here the removable singularity at $x=1$ is cancelled by [Together](https://reference.wolfram.com/language/ref/Together.html), which is used to preprocess the equation:

```wolfram
Solve[(x^2-2x+1)/(x-1)==0,x]
(* Output *)
{{x->1}}
```

```wolfram
Together[(x^2-2x+1)/(x-1)==0]
(* Output *)
-1+x==0
```

The value of [MaxRoots](https://reference.wolfram.com/language/ref/MaxRoots.html) is used only for systems with numeric coefficients:

```wolfram
Solve[x^4==y,{x,y},MaxRoots->3]
(* Output *)
{{x->0,y->0},{x->(-(3)/(13))^(1/3) 2^(2/3),y->-(12)/(13) (-(3)/(13))^(1/3) 2^(2/3)},{x->-((3)/(13))^(1/3) 2^(2/3),y->(12)/(13) ((3)/(13))^(1/3) 2^(2/3)}}
```

When symbolic parameters are present, the option value is ignored:

```wolfram
Solve[x^4==y,{x},MaxRoots->3]
(* Output *)
Solve
(* Output *)
{{x->-y^(1/4)},{x->-ⅈ y^(1/4)},{x->ⅈ y^(1/4)},{x->y^(1/4)}}
```

Expressions given as variables are treated as atomic objects and not as functions of their subexpressions:

```wolfram
Solve[x+Sin[x]==1,Sin[x]]
(* Output *)
{{Sin[x]->1-x}}
```

Effectively, variables are replaced with new symbols before the equations are solved:

```wolfram
Solve[f[x]+Integrate[f[x],{x,0,1}]==2x,f[x]]
(* Output *)
{{f[x]->x}}
```

The result comes from:

```wolfram
f[x]+Integrate[f[x],{x,0,1}]==2x/.f[x]->z
(* Output *)
2 z==2 x
```

## Tech Notes ▪Symbolic Mathematics: Basic Operations ▪Solving Equations ▪Simultaneous Equations ▪Solving Equations Involving Power Series ▪Solving Linear Systems ▪Solving Logical Combinations of Equations ▪Generic and Non[Hyphen]Generic Solutions ▪Eliminating Variables

## Related Guides ▪Equation Solving ▪Polynomial Algebra ▪Polynomial Equations ▪Solvers over Regions ▪Polynomial Systems ▪Plane Geometry ▪Solid Geometry ▪Polygons ▪Precollege Education ▪Symbolic Vectors, Matrices and Arrays ▪Finite Mathematics ▪Finite Fields ▪Polyhedra ▪Geometric Computation ▪Additive Number Theory

## History Introduced in 1988 (1.0) | Updated in 1996 (3.0) ▪ 2014 (10.0) ▪ 2020 (12.2) ▪ 2024 (14.0)
