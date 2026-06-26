# Reduce | [SpanFromLeft]

> [Reduce](https://reference.wolfram.com/language/ref/Reduce.html)[*expr*,*vars*] — reduces the statement `*expr*` by solving equations or inequalities for `*vars*` and eliminating quantifiers.
> [Reduce](https://reference.wolfram.com/language/ref/Reduce.html)[*expr*,*vars*,*dom*] — does the reduction over the domain `*dom*`. Common choices of `*dom*` are [Reals](https://reference.wolfram.com/language/ref/Reals.html), [Integers](https://reference.wolfram.com/language/ref/Integers.html), and [Complexes](https://reference.wolfram.com/language/ref/Complexes.html).

## Details and Options

The statement `*expr*` can be any logical combination of:

*lhs*==*rhs* | equations
*lhs*!=*rhs* | inequations
`*lhs*>*rhs*` or `*lhs*>=*rhs*` | inequalities
*expr*∈*dom* | domain specifications
{*x*,*y*,…}∈*reg* | region specification
[ForAll](https://reference.wolfram.com/language/ref/ForAll.html)[*x*,*cond*,*expr*] | universal quantifiers
[Exists](https://reference.wolfram.com/language/ref/Exists.html)[*x*,*cond*,*expr*] | existential quantifiers

The result of [Reduce](https://reference.wolfram.com/language/ref/Reduce.html)[*expr*,*vars*] always describes exactly the same mathematical set as `*expr*`.

[Reduce](https://reference.wolfram.com/language/ref/Reduce.html)[{*expr*_1,*expr*_2,…},*vars*] is equivalent to [Reduce](https://reference.wolfram.com/language/ref/Reduce.html)[*expr*_1&&*expr*_2&&…,*vars*].

[Reduce](https://reference.wolfram.com/language/ref/Reduce.html)[*expr*,*vars*] assumes by default that quantities appearing algebraically in inequalities
are real, while all other quantities are complex.

[Reduce](https://reference.wolfram.com/language/ref/Reduce.html)[*expr*,*vars*,*dom*] restricts all variables and parameters to belong to the domain `*dom*`.

If `*dom*` is [Reals](https://reference.wolfram.com/language/ref/Reals.html), or a subset such as [Integers](https://reference.wolfram.com/language/ref/Integers.html) or [Rationals](https://reference.wolfram.com/language/ref/Rationals.html), then all constants and function values are also restricted to be real.

[Reduce](https://reference.wolfram.com/language/ref/Reduce.html)[*expr*&&*vars*∈[Reals](https://reference.wolfram.com/language/ref/Reals.html),*vars*,[Complexes](https://reference.wolfram.com/language/ref/Complexes.html)] performs reductions with variables assumed real, but function values allowed to be complex.

[Reduce](https://reference.wolfram.com/language/ref/Reduce.html)[*expr*,*vars*,[Integers](https://reference.wolfram.com/language/ref/Integers.html)] reduces Diophantine equations over the integers.

[Reduce](https://reference.wolfram.com/language/ref/Reduce.html)[…,*x*∈*reg*,[Reals](https://reference.wolfram.com/language/ref/Reals.html)] constrains `*x*` to be in the region `*reg*`. The different coordinates for `*x*` can be referred to using [Indexed](https://reference.wolfram.com/language/ref/Indexed.html)[*x*,*i*].

[Reduce](https://reference.wolfram.com/language/ref/Reduce.html)[*expr*,{*x*_1,*x*_2,…},…] effectively writes `*expr*` as a combination of conditions on `*x*_1, `*x*_2, … where each condition involves only the earlier $x_{i}$.

Algebraic variables in `*expr*` free of the $x_{i}$ and of each other are treated as independent parameters.

Applying [LogicalExpand](https://reference.wolfram.com/language/ref/LogicalExpand.html) to the results of [Reduce](https://reference.wolfram.com/language/ref/Reduce.html)[*expr*,…] yields an expression of the form $e_{1}||e_{2}||\ldots$, where each of the $e_{i}$ can be thought of as representing a separate component in the set defined by `*expr*`.

The $e_{i}$ may not be disjoint and may have different dimensions. After [LogicalExpand](https://reference.wolfram.com/language/ref/LogicalExpand.html), each of the $e_{i}$ has the form $e&&e&&\ldots$.

Without [LogicalExpand](https://reference.wolfram.com/language/ref/LogicalExpand.html), [Reduce](https://reference.wolfram.com/language/ref/Reduce.html) by default returns a nested collection of conditions on the $x_{i}$, combined alternately by [Or](https://reference.wolfram.com/language/ref/Or.html) and [And](https://reference.wolfram.com/language/ref/And.html) on successive levels.

When `*expr*` involves only polynomial equations and inequalities over real or complex domains, then [Reduce](https://reference.wolfram.com/language/ref/Reduce.html) can always in principle solve directly for all the $x_{i}$.

When `*expr*` involves transcendental conditions or integer domains, [Reduce](https://reference.wolfram.com/language/ref/Reduce.html) will often introduce additional parameters in its results.

When `*expr*` involves only polynomial conditions, [Reduce](https://reference.wolfram.com/language/ref/Reduce.html)[*expr*,*vars*,[Reals](https://reference.wolfram.com/language/ref/Reals.html)] gives a cylindrical algebraic decomposition of `*expr*`.

[Reduce](https://reference.wolfram.com/language/ref/Reduce.html) can give explicit representations for solutions to all linear equations and inequalities over the integers and can solve a large fraction of Diophantine equations described in the literature.

When `*expr*` involves only polynomial conditions over real or complex domains, [Reduce](https://reference.wolfram.com/language/ref/Reduce.html)[*expr*,*vars*] will always eliminate quantifiers, so that quantified variables do not appear in the result.

The following options can be given:

| Backsubstitution | [False](https://reference.wolfram.com/language/ref/False.html) | whether to give results unwound by backsubstitution  » |
| --- | --- | --- |
| [Cubics](https://reference.wolfram.com/language/ref/Cubics.html) | [False](https://reference.wolfram.com/language/ref/False.html) | whether to use explicit radicals to solve all cubics  » |
| [GeneratedParameters](https://reference.wolfram.com/language/ref/GeneratedParameters.html) | [C](https://reference.wolfram.com/language/ref/C.html) | how to name parameters that are generated  » |
| [Modulus](https://reference.wolfram.com/language/ref/Modulus.html) | 0 | modulus to assume for integers  » |
| [Quartics](https://reference.wolfram.com/language/ref/Quartics.html) | [False](https://reference.wolfram.com/language/ref/False.html) | whether to use explicit radicals to solve all quartics  » |

[Reduce](https://reference.wolfram.com/language/ref/Reduce.html)[*expr*,{*x*_1,*x*_2,…},Backsubstitution->[True](https://reference.wolfram.com/language/ref/True.html)] yields a form in which values from equations generated for earlier $x_{i}$ are backsubstituted so that the conditions for a particular $x_{i}$ have only minimal dependence on earlier $x_{i}$.

## Examples

### Basic Examples

Reduce equations and inequalities:

```wolfram
Reduce[x^2-y^3==1,{x,y}]
(* Output *)
y==(-1+x^2)^(1/3)||y==-(-1)^(1/3) (-1+x^2)^(1/3)||y==(-1)^(2/3) (-1+x^2)^(1/3)
```

```wolfram
Reduce[x^2+y^2<1,{x,y}]
(* Output *)
-1<x<1&&-Sqrt[1-x^2]<y<Sqrt[1-x^2]
```

Use specific domains:

```wolfram
Reduce[x^2-7y^2==1&&x>0&&y>0,{x,y},Integers]
(* Output *)
1∈Integers&&1>=1&&x==(1)/(2) ((8-3 Sqrt[7])^1+(8+3 Sqrt[7])^1)&&y==-((8-3 Sqrt[7])^1-(8+3 Sqrt[7])^1)/(2 Sqrt[7])
```

```wolfram
Reduce[x^2-7y^2==1&&x>0&&y>0,{x,y},Reals]
(* Output *)
x>1&&y==(Sqrt[-1+x^2])/(Sqrt[7])
```

Reduce a quantified expression:

```wolfram
Reduce[Exists[{x,y},x^2+a y^2<=1 &&x-y>=2],a]
(* Output *)
a<=(1)/(3)
```

Reduce with geometric region constraints:

```wolfram
Reduce[{x,y}∈InfiniteLine[{{0,0},{2,2}}]&&{x,y}∈Circle[],{x,y}]
(* Output *)
(x==-(1)/(Sqrt[2])||x==(1)/(Sqrt[2]))&&y==x
```

```wolfram
Graphics[{Blue,Circle[],InfiniteLine[{{0,0},{2,2}}],{Red,Point[{x,y}]//.{ToRules[%]}}}]
```

*([Graphics])*

### Scope

#### Basic Uses

Find an explicit description of the solution set of a system of equations:

```wolfram
Reduce[x^2+y^2==2&&x-y==1,{x,y}]
(* Output *)
(x==(1)/(2) (1-Sqrt[3])||x==(1)/(2) (1+Sqrt[3]))&&y==-1+x
```

Use [ToRules](https://reference.wolfram.com/language/ref/ToRules.html) and [ReplaceRepeated](https://reference.wolfram.com/language/ref/ReplaceRepeated.html) (`//.`) to list the solutions:

```wolfram
{x,y}//.{ToRules[%]}
(* Output *)
{{(1)/(2) (1-Sqrt[3]),-1+(1)/(2) (1-Sqrt[3])},{(1)/(2) (1+Sqrt[3]),-1+(1)/(2) (1+Sqrt[3])}}
```

Find an explicit description of the solution set of a system of inequalities:

```wolfram
Reduce[x^2+y^2<2&&x-y>=1,{x,y}]
(* Output *)
((1)/(2) (1-Sqrt[3])<x<(1)/(2) (1+Sqrt[3])&&-Sqrt[2-x^2]<y<=-1+x)||(x==(1)/(2) (1+Sqrt[3])&&-Sqrt[2-(1)/(4) (1+Sqrt[3])^2]<y<-1+(1)/(2) (1+Sqrt[3]))||((1)/(2) (1+Sqrt[3])<x<Sqrt[2]&&-Sqrt[2-x^2]<y<Sqrt[2-x^2])
```

Find solutions over specified domains:

```wolfram
Reduce[(x^4-1)(x^4-4)==0,x,Complexes]
(* Output *)
x==-1||x==-ⅈ||x==ⅈ||x==1||x==-Sqrt[2]||x==-ⅈ Sqrt[2]||x==ⅈ Sqrt[2]||x==Sqrt[2]
```

```wolfram
Reduce[(x^4-1)(x^4-4)==0,x,Reals]
(* Output *)
x==-1||x==1||x==-Sqrt[2]||x==Sqrt[2]
```

```wolfram
Reduce[(x^4-1)(x^4-4)==0,x,Integers]
(* Output *)
x==-1||x==1
```

The solution set may depend on symbolic parameters:

```wolfram
Reduce[a x^2+b x+c==0,x]
(* Output *)
(a≠0&&(x==(-b-Sqrt[b^2-4 a c])/(2 a)||x==(-b+Sqrt[b^2-4 a c])/(2 a)))||(a==0&&b≠0&&x==-(c)/(b))||(c==0&&b==0&&a==0)
```

Representing solutions may require introduction of new parameters:

```wolfram
Reduce[Sin[x]==1/2,x]
(* Output *)
1∈Integers&&(x==(π)/(6)+2 π 1||x==(5 π)/(6)+2 π 1)
```

```wolfram
Reduce[x^2-2y^2==1&&x>0&&y>0,{x,y},Integers]
(* Output *)
1∈Integers&&1>=1&&x==(1)/(2) ((3-2 Sqrt[2])^1+(3+2 Sqrt[2])^1)&&y==-((3-2 Sqrt[2])^1-(3+2 Sqrt[2])^1)/(2 Sqrt[2])
```

List the first 10 solutions:

```wolfram
%/.Table[{1->i},{i,10}]//Simplify
(* Output *)
{x==3&&y==2,x==17&&y==12,x==99&&y==70,x==577&&y==408,x==3363&&y==2378,x==19601&&y==13860,x==114243&&y==80782,x==665857&&y==470832,x==3880899&&y==2744210,x==22619537&&y==15994428}
```

#### Complex Domain

A linear system:

```wolfram
Reduce[2 x+3y-5z==1&&3x-4y+7z==3,{x,y,z}]
(* Output *)
y==22-29 x&&z==13-17 x
```

A univariate polynomial equation:

```wolfram
Reduce[x^3-2x+1==0,x]
(* Output *)
x==1||x==(1)/(2) (-1-Sqrt[5])||x==(1)/(2) (-1+Sqrt[5])
```

A multivariate polynomial equation:

```wolfram
Reduce[x^2-y z==1,{x,y,z}]
(* Output *)
((x==-1||x==1)&&y==0)||(y≠0&&z==(-1+x^2)/(y))
```

Systems of polynomial equations and inequations can always be reduced:

```wolfram
Reduce[x^2+y^3==z&&x+2y==3z+1&&x y z≠0,{x,y,z}]
(* Output *)
(y==Root[1-x+3 x^2-2 #1+3 #1^3&,1]||y==Root[1-x+3 x^2-2 #1+3 #1^3&,2]||y==Root[1-x+3 x^2-2 #1+3 #1^3&,3])&&z==(1)/(3) (-1+x+2 y)&&-x y+x^2 y+2 x y^2≠0
```

A quantified polynomial system:

```wolfram
Reduce[ForAll[x,Exists[y,a x^2+b y^2-3y==1&&y≠0]],{a,b}]
(* Output *)
a==0||b≠0
```

An algebraic system:

```wolfram
Reduce[Sqrt[x+2y]-3x+4y==5&&x+y^(1/3)==1,{x,y}]
(* Output *)
x==Root&&y==-(-1+Root)^3
```

Transcendental equations solvable in terms of inverse functions:

```wolfram
Reduce[Sin[x]==1/3,x]
(* Output *)
1∈Integers&&(x==π-ArcSin[(1)/(3)]+2 π 1||x==ArcSin[(1)/(3)]+2 π 1)
```

```wolfram
Reduce[ 4^(x^2)2^x==8,x]
(* Output *)
1∈Integers&&(x==(1)/(4) (-1-Sqrt[1+8 (3+(2 ⅈ π 1)/(Log[2]))])||x==(1)/(4) (-1+Sqrt[1+8 (3+(2 ⅈ π 1)/(Log[2]))]))
```

In this case there is no solution:

```wolfram
Reduce[Log[x]==75/11 I Pi+17,x]
(* Output *)
False
```

Equations involving elliptic functions:

```wolfram
Reduce[JacobiSN[x,y]==1,x]
(* Output *)
(1|2)∈Integers&&x==2 ⅈ 2 EllipticK[1-y]+EllipticK[y]+4 1 EllipticK[y]
```

Equations solvable using special function zeros:

```wolfram
Reduce[Zeta[x]==0,x]
(* Output *)
Reduce
(* Output *)
1∈Integers&&((1>=1&&x==-2 1)||(1≠0&&x==ZetaZero[1]))
```

Solving this system does not require the Riemann hypothesis:

```wolfram
Reduce[Zeta[x]==0&&Re[x]==1/2&&Im[x]^2<500,x]
(* Output *)
x==ZetaZero[-2]||x==ZetaZero[-1]||x==ZetaZero[1]||x==ZetaZero[2]
```

Elementary function equation in a bounded region:

```wolfram
Reduce[Sin[E^x]-Cos[2 x]==1&&-1<=Re[x]<=1&&-1<=Im[x]<=1,x]
(* Output *)
x==Root||x==Root
```

Holomorphic function equation in a bounded region:

```wolfram
Reduce[Gamma[x]-Log[x]==I/2&&Abs[x-2]<3/2,x]
(* Output *)
x==Root||x==Root
```

Here [Reduce](https://reference.wolfram.com/language/ref/Reduce.html) finds some solutions but is not able to prove there are no other solutions:

```wolfram
Reduce[x==E^(1/x)&&Abs[x]<5,x]
(* Output *)
Reduce
(* Output *)
x==Root||x==Root||x==Root||x==Root||x==Root
```

Equation with a purely imaginary period over a vertical stripe in the complex plane:

```wolfram
Reduce[Cos[Exp[x]]==3 Exp[-x]+1&&0<=Re[x]<=1,x]
(* Output *)
1∈Integers&&(x==2 ⅈ π 1+Root||x==2 ⅈ π 1+Root||x==2 ⅈ π 1+Root)
```

Doubly periodic transcendental equation:

```wolfram
Reduce[JacobiCS[x,3]==1+I,x]
(* Output *)
(1|2)∈Integers&&(x==2 ⅈ EllipticK[-2]+4 ⅈ 2 EllipticK[-2]+2 1 EllipticK[3]-InverseJacobiCS[1+ⅈ,3]||x==4 ⅈ 2 EllipticK[-2]+2 1 EllipticK[3]+InverseJacobiCS[1+ⅈ,3])
```

A system of transcendental equations solvable using inverse functions:

```wolfram
Reduce[Sin[x+y]==1/2&&E^x-y==1,{x,y}]
(* Output *)
(1|2)∈Integers&&(x==1-(7 π)/(6)-2 π 2-ProductLog[1,ℯ^(1-(7 π)/(6)-2 π 2)]||x==1+(π)/(6)-2 π 2-ProductLog[1,ℯ^(1+(π)/(6)-2 π 2)])&&y==-1+ℯ^x
```

A square system of analytic equations over a bounded box:

```wolfram
Reduce[FresnelS[x-y]-AiryAi[x]+y==(3+I)/2&&CosIntegral[x] y-Sin[x y]==-5(1-I)/4&&1<=Re[x]<=2&&0<=Im[x]<=1&&1<=Re[y]<=2&&0<=Im[y]<=1,{x,y}]
(* Output *)
x==Root&&y==Root
```

#### Real Domain

A linear system:

```wolfram
Reduce[2 x+3y-5z==1&&3x-4y+7z==3,{x,y,z},Reals]
(* Output *)
y==22-29 x&&z==13-17 x
```

A univariate polynomial equation:

```wolfram
Reduce[x^5-2x+1==0,x,Reals]
(* Output *)
x==1||x==Root||x==Root
```

A univariate polynomial inequality:

```wolfram
Reduce[x^5-2x+1<0,x,Reals]
(* Output *)
x<Root||Root<x<1
```

A multivariate polynomial equation:

```wolfram
Reduce[x^2-y z==1,{x,y,z},Reals]
(* Output *)
(x==-1&&y==0)||(x==1&&y==0)||((y<0||y>0)&&z==-(1-x^2)/(y))
```

A multivariate polynomial inequality:

```wolfram
Reduce[x^2-2y +z^2<=1,{x,y,z},Reals]
(* Output *)
(y==(1)/(2) (-1+x^2)&&z==0)||(y>(1)/(2) (-1+x^2)&&-Sqrt[1-x^2+2 y]<=z<=Sqrt[1-x^2+2 y])
```

Systems of polynomial equations and inequalities can always be reduced:

```wolfram
Reduce[x^2+y z==1&&x+2y<=3z+1&&x y z>7,{x,y,z},Reals]
(* Output *)
x<Root&&y<0&&z==(1-x^2)/(y)
```

A quantified polynomial system:

```wolfram
Reduce[ForAll[x,Exists[y,a x^2+b y^2-3y==1&&y<0]],{a,b},Reals]
(* Output *)
(a<0&&b>=0)||(a==0&&b>=-(9)/(4))||(a>0&&-(9)/(4)<=b<0)
```

An algebraic system:

```wolfram
Reduce[Sqrt[x+2y]-3x+4y>=5&&x+y^(1/3)==1,{x,y},Reals]
(* Output *)
x<=Root&&y==1-3 x+3 x^2-x^3
```

Piecewise equations:

```wolfram
Reduce[Abs[(x+Abs[x+2])^2-1]^2==9,x,Reals]
(* Output *)
x<=-2||x==0
```

```wolfram
Reduce[Max[x,y]==Min[y^2-x,x],{x,y}]
(* Output *)
(x<=0&&y<=x)||(0<x<2&&y<=-Sqrt[2] Sqrt[x])||(x==2&&(y<=-2||y==2))||(x>2&&(y<=-Sqrt[2] Sqrt[x]||Sqrt[2] Sqrt[x]<=y<=x))
```

Piecewise inequalities:

```wolfram
Reduce[Abs[3x^2-7x-6]<Abs[x^2+x],x,Reals]
(* Output *)
(1)/(4) (3-Sqrt[33])<x<2-Sqrt[7]||(1)/(4) (3+Sqrt[33])<x<2+Sqrt[7]
```

```wolfram
Reduce[Floor[x^2+Ceiling[x^2]]<10,x,Reals]
(* Output *)
-Sqrt[5]<x<Sqrt[5]
```

Transcendental equations, solvable using inverse functions:

```wolfram
Reduce[E^x-x==7,x,Reals]
(* Output *)
x==-7-ProductLog[-(1)/(ℯ^7)]||x==-7-ProductLog[-1,-(1)/(ℯ^7)]
```

```wolfram
Reduce[ (27^(2x-1))^(1/x)==Sqrt[9^(2x-1)], x , Reals]
(* Output *)
x==(1)/(2)||x==3
```

Transcendental inequalities, solvable using inverse functions:

```wolfram
Reduce[Sin[x]<1/3,x,Reals]
(* Output *)
1∈Integers&&-π-ArcSin[(1)/(3)]+2 π 1<x<ArcSin[(1)/(3)]+2 π 1
```

```wolfram
Reduce[ (1)/(2^x-1)>(1)/(1-2^(x-1)), x ,Reals]
(* Output *)
0<x<(2 Log[2]-Log[3])/(Log[2])||x>1
```

Inequalities involving elliptic functions:

```wolfram
Reduce[1<JacobiNC[x,3]<=2,x,Reals]
(* Output *)
1∈Integers&&(2 1 EllipticK[(1)/(3)])/(Sqrt[3])<x<(1)/(3) (2 Sqrt[3] EllipticK[(1)/(3)]+2 Sqrt[3] 1 EllipticK[(1)/(3)])
```

Transcendental equation, solvable using special function zeros:

```wolfram
Reduce[AiryBi[1-x^2]==0&&99<x<100,x,Reals]
(* Output *)
1∈Integers&&205874<=1<=212175&&x==Sqrt[1-AiryBiZero[1]]
```

Transcendental inequality, solvable using special function zeros:

```wolfram
Reduce[900<AiryAiZero[2t+1]^2<1000,t,Reals]
(* Output *)
t==(35)/(2)||t==18
```

Exp-log equations:

```wolfram
Reduce[E^(2E^x)-Log[x^2+1]-20x==11,x,Reals]
(* Output *)
x==Root||x==Root
```

High-degree sparse polynomial equation:

```wolfram
Reduce[x^1000000-2x^777777+3x^12345+9x^67-10==0,x,Reals]
(* Output *)
x==Root||x==Root||x==Root||x==Root
```

Algebraic equation involving high-degree radicals:

```wolfram
Reduce[2x^(123451/67890)-x^2+4Sqrt[x]-4x-9/8==0,x,Reals]
(* Output *)
x==Root^67890||x==Root^67890||x==Root^67890||x==Root^67890
```

Equation involving irrational real powers:

```wolfram
Reduce[x^Pi-x^x^Sqrt[2]-Sqrt[3]x+2^(1/3)==0,x,Reals]
(* Output *)
x==Root||x==Root||x==Root
```

Exp-log inequality:

```wolfram
Reduce[E^(2E^x)-Log[x^2+1]-20x<11,x,Reals]
(* Output *)
Root<x<Root
```

Elementary function equation in a bounded interval:

```wolfram
Reduce[2Sin[Exp[x]]-Cos[Pi x]==3/2&&-1<x<1,x,Reals]
(* Output *)
x==Root||x==Root
```

Holomorphic function equation in a bounded interval:

```wolfram
Reduce[Cos[x]-BesselJ[5,x]==1/2&&0<=x<=10,x]
(* Output *)
x==Root||x==Root||x==Root
```

Meromorphic function inequality in a bounded interval:

```wolfram
Reduce[Gamma[x]<=2x&&-2<=x<=2,x,Reals]
(* Output *)
Reduce
(* Output *)
Gamma[x]∈Reals&&(-1<x<0||Root<=x<=2)
```

Periodic elementary function equation over the reals:

```wolfram
Reduce[Exp[Sin[x]]-Sin[3 Cos[x]]==0,x,Reals]
(* Output *)
1∈Integers&&(x==2 π 1+Root||x==2 π 1+Root)
```

Transcendental systems solvable using inverse functions:

```wolfram
Reduce[Sin[x+y]==1/2&&E^x-y<=1,{x,y},Reals]
(* Output *)
1∈Integers&&((1>=0&&x==0&&(y==(1)/(6) (π+12 π 1)||y==(1)/(6) (5 π+12 π 1)))||(x<=(1)/(6) (6+π+12 π 1-6 ProductLog[ℯ^(1+(π)/(6)+2 π 1)])&&y==(1)/(6) (π-6 x+12 π 1))||(x<=(1)/(6) (6+5 π+12 π 1-6 ProductLog[ℯ^(1+(5 π)/(6)+2 π 1)])&&y==(1)/(6) (5 π-6 x+12 π 1)))
```

```wolfram
Reduce[ 3^x-2^2y==77 && Sqrt[3^x]-2^y==7, {x, y}, Reals]
(* Output *)
x==4&&y==1
```

Systems exp-log in the first variable and polynomial in the other variables:

```wolfram
Reduce[E^x y^3+Log[x]y==1&&x y+E^x/x>=2,{x,y},Reals]
(* Output *)
(0<x<Root&&(y==Root[-1+Log[x] #1+ℯ^x #1^3&,1]||y==Root[-1+Log[x] #1+ℯ^x #1^3&,2]||y==Root[-1+Log[x] #1+ℯ^x #1^3&,3]))||(x==Root&&(y==Root[-1+Log[x] #1+ℯ^x #1^3&,1]||y==Root[-1+Log[x] #1+ℯ^x #1^3&,3]))||(x>Root&&y==Root[-1+Log[x] #1+ℯ^x #1^3&,1])
```

```wolfram
Reduce[E^x y^3+Log[x] y>3x, {x, y}, Reals]
(* Output *)
(0<x<Root&&(Root[-3 x+Log[x] #1+ℯ^x #1^3&,1]<y<Root[-3 x+Log[x] #1+ℯ^x #1^3&,2]||y>Root[-3 x+Log[x] #1+ℯ^x #1^3&,3]))||(x==Root&&y>Root[-3 x+Log[x] #1+ℯ^x #1^3&,3])||(x>Root&&y>Root[-3 x+Log[x] #1+ℯ^x #1^3&,1])
```

Quantified system:

```wolfram
Reduce[Exists[a,a x^2+Sinh[x^2+1]a^2>=1&&x^2+a^2<=1],x,Reals]
(* Output *)
Root<=x<=Root
```

Systems elementary and bounded in the first variable and polynomial in the other variables:

```wolfram
Reduce[Sin[x-Cos[x]] y^3-x==1&&x^2+y^2<=1,{x,y},Reals]
(* Output *)
(x==-1&&y==0)||(Root<=x<=Root&&y==Root[-1-x+Sin[x-Cos[x]] #1^3&,1])
```

Quantified system:

```wolfram
Reduce[Exists[y,y^3-Cos[x] y+2 x^2 Sin[x^2-1]>0&&x^2+Cos[x]^2<=2],x,Reals]
(* Output *)
Root<=x<=Root
```

Systems analytic and bounded in the first variable and polynomial in the other variables:

```wolfram
Reduce[y^3-BesselJ[2,x+2] y-y-3 x==-2&&y<0 &&x^2<2,{x,y},Reals]
(* Output *)
(-Sqrt[2]<x<(2)/(3)&&y==Root[2-3 x+(-1-BesselJ[2,2+x]) #1+#1^3&,1])||(x==(2)/(3)&&y==-Sqrt[1+BesselJ[2,(8)/(3)]])||((2)/(3)<x<Root&&(y==Root[2-3 x+(-1-BesselJ[2,2+x]) #1+#1^3&,1]||y==Root[2-3 x+(-1-BesselJ[2,2+x]) #1+#1^3&,2]))||(x==Root&&y==Root[2-3 x+(-1-BesselJ[2,2+x]) #1+#1^3&,1])
```

Quantified system:

```wolfram
Reduce[Exists[y,y^4-Gamma[x+2] y-y-3 ArcSin[x/3]==1&&y<-1/2&&0<x<2],x,Reals]
(* Output *)
Gamma[2+x]∈Reals&&Root<x<Root
```

Square systems of analytic equations over bounded regions:

```wolfram
Reduce[Gamma[x+y+1]-Sin[x y]==1&&Erf[x^2-y]-E^y-x+4==0&&0<x<3&&0<y<3,{x,y}, Reals]
(* Output *)
(x==Root&&y==Root)||(x==Root&&y==Root)
```

#### Integer Domain

Linear system of equations:

```wolfram
Reduce[2 x+3y-5z==1&&3x-4y+7z==3,{x,y,z},Integers]
(* Output *)
1∈Integers&&x==1&&y==22-29 1&&z==13-17 1
```

A linear system of equations and inequalities:

```wolfram
Reduce[2 x+3y==4&&3x-4y<=5&&x-2y>-21,{x,y,z},Integers]
(* Output *)
z∈Integers&&((x==-7&&y==6)||(x==-4&&y==4)||(x==-1&&y==2))
```

A univariate polynomial equation:

```wolfram
Reduce[x^1000-2x^777+1==0,x,Integers]
(* Output *)
x==1
```

A univariate polynomial inequality:

```wolfram
Reduce[x^5-2x+1<0,x,Integers]
(* Output *)
x∈Integers&&x<=-2
```

Binary quadratic equations:

```wolfram
Reduce[x^2+x y+y^2==109,{x,y},Integers]
(* Output *)
(x==-12&&y==5)||(x==-12&&y==7)||(x==-7&&y==-5)||(x==-7&&y==12)||(x==-5&&y==-7)||(x==-5&&y==12)||(x==5&&y==-12)||(x==5&&y==7)||(x==7&&y==-12)||(x==7&&y==5)||(x==12&&y==-7)||(x==12&&y==-5)
```

```wolfram
Reduce[x^2-3y^2==22&&x>0&&y>0,{x,y},Integers]
(* Output *)
(x==5&&y==1)||(1∈Integers&&1>=0&&x==(1)/(2) (5 (2-Sqrt[3])^1-Sqrt[3] (2-Sqrt[3])^1+5 (2+Sqrt[3])^1+Sqrt[3] (2+Sqrt[3])^1)&&y==(1)/(6) (3 (2-Sqrt[3])^1-5 Sqrt[3] (2-Sqrt[3])^1+3 (2+Sqrt[3])^1+5 Sqrt[3] (2+Sqrt[3])^1))||(1∈Integers&&1>=1&&x==(1)/(2) (5 (2-Sqrt[3])^1+Sqrt[3] (2-Sqrt[3])^1+5 (2+Sqrt[3])^1-Sqrt[3] (2+Sqrt[3])^1)&&y==(1)/(6) (-3 (2-Sqrt[3])^1-5 Sqrt[3] (2-Sqrt[3])^1-3 (2+Sqrt[3])^1+5 Sqrt[3] (2+Sqrt[3])^1))
```

```wolfram
Reduce[x^2-6 x y+9y^2-x+2y==1,{x,y},Integers]
(* Output *)
1∈Integers&&x==-3-2 1+3 1^2&&y==-1-1+1^2
```

A Thue equation:

```wolfram
Reduce[x^3-2x^2 y+y^3==2,{x,y},Integers]
(* Output *)
(x==1&&y==-1)||(x==5&&y==3)
```

A sum of squares equation:

```wolfram
Reduce[x^2+4y^2+9z^2+16t^2==354&&x>0&&y>0&&z>0&&t>0,{x,y,z,t},Integers]
(* Output *)
(x==1&&y==2&&z==3&&t==4)||(x==1&&y==4&&z==5&&t==2)||(x==1&&y==8&&z==3&&t==1)||(x==5&&y==4&&z==1&&t==4)||(x==5&&y==8&&z==1&&t==2)||(x==7&&y==2&&z==5&&t==2)||(x==7&&y==4&&z==5&&t==1)
```

The Pythagorean equation:

```wolfram
Reduce[x^2+y^2==z^2,{x,y,z},Integers]
(* Output *)
(1|2|3)∈Integers&&3>=0&&((x==1 (2^2-3^2)&&y==2 1 2 3&&z==1 (2^2+3^2))||(x==2 1 2 3&&y==1 (2^2-3^2)&&z==1 (2^2+3^2)))
```

A bounded system of equations and inequalities:

```wolfram
Reduce[x^4+y^4+z^4<=500&&x+y^2+z^3==32,{x,y,z},Integers]
(* Output *)
(x==-4&&y==-3&&z==3)||(x==-4&&y==3&&z==3)||(x==1&&y==-2&&z==3)||(x==1&&y==2&&z==3)||(x==4&&y==-1&&z==3)||(x==4&&y==1&&z==3)
```

A high-degree system with no solution:

```wolfram
Reduce[2x^7+8y^15+14 x y z==3,{x,y,z},Integers]
(* Output *)
False
```

Transcendental Diophantine systems:

```wolfram
Reduce[Exp[y^2]<x&&Abs[x]<5&&Abs[y]<5,{x,y},Integers]
(* Output *)
(x==2&&y==0)||(x==3&&(y==-1||y==0||y==1))||(x==4&&(y==-1||y==0||y==1))
```

```wolfram
Reduce[Exp[x^2-5y^2+1]+x^2-5y^2==0&&x>0&&y>0,{x,y},Integers]
(* Output *)
(x==2&&y==1)||(1∈Integers&&1>=1&&x==(1)/(2) (-2 (9-4 Sqrt[5])^1-Sqrt[5] (9-4 Sqrt[5])^1-2 (9+4 Sqrt[5])^1+Sqrt[5] (9+4 Sqrt[5])^1)&&y==(1)/(10) (5 (9-4 Sqrt[5])^1+2 Sqrt[5] (9-4 Sqrt[5])^1+5 (9+4 Sqrt[5])^1-2 Sqrt[5] (9+4 Sqrt[5])^1))
```

A polynomial system of congruences:

```wolfram
Reduce[Mod[x^2+y^2,2]==1&&Mod[x-2y,3]==2,{x,y},Integers]
(* Output *)
(1|2)∈Integers&&((x==6 1&&y==5+6 2)||(x==4+6 1&&y==1+6 2)||(x==2+6 1&&y==3+6 2)||(x==3+6 1&&y==2+6 2)||(x==1+6 1&&y==4+6 2)||(x==5+6 1&&y==6 2))
```

Diophantine equations with irrational coefficients:

```wolfram
Reduce[GoldenRatio x^2-Sqrt[5] y^2-z==Sqrt[5],{x,y,z},PositiveIntegers]
(* Output *)
(x==2&&y==1&&z==2)||(1∈Integers&&1>=1&&x==(1)/(2) (2 (3-2 Sqrt[2])^1+Sqrt[2] (3-2 Sqrt[2])^1+2 (3+2 Sqrt[2])^1-Sqrt[2] (3+2 Sqrt[2])^1)&&y==(1)/(2) (-(3-2 Sqrt[2])^1-Sqrt[2] (3-2 Sqrt[2])^1-(3+2 Sqrt[2])^1+Sqrt[2] (3+2 Sqrt[2])^1)&&z==(1)/(8) (2 (3-2 Sqrt[2])^1+Sqrt[2] (3-2 Sqrt[2])^1+2 (3+2 Sqrt[2])^1-Sqrt[2] (3+2 Sqrt[2])^1)^2)
```

```wolfram
Reduce[Pi x^2+Sqrt[6]x y-Sqrt[3] y^3 +Pi^2 Sqrt[3] y^2z^2  -3Sqrt[3] y z==Sqrt[2(109+12Sqrt[3])]+4Pi+144Sqrt[3]Pi^2-Sqrt[2]- 63 Sqrt[3],{x,y,z},Integers]
(* Output *)
x==2&&y==3&&z==4
```

#### Modular Domains

A linear system:

```wolfram
Reduce[2 x+3y-5z==1&&3x-4y+7z==3,{x,y,z},Modulus->12]
(* Output *)
x==1&&y==10+7 1&&z==1+7 1
```

A univariate polynomial equation:

```wolfram
Reduce[x^3-2x+1==0,x,Modulus->5]
(* Output *)
x==1||x==2
```

A multivariate polynomial equation:

```wolfram
Reduce[x^2-y z==1,{x,y,z},Modulus->4]
(* Output *)
(x==0&&y==1&&z==3)||(x==0&&y==3&&z==1)||(x==1&&y==0&&z==0)||(x==1&&y==0&&z==1)||(x==1&&y==0&&z==2)||(x==1&&y==0&&z==3)||(x==1&&y==1&&z==0)||(x==1&&y==2&&z==0)||(x==1&&y==2&&z==2)||(x==1&&y==3&&z==0)||(x==2&&y==1&&z==3)||(x==2&&y==3&&z==1)||(x==3&&y==0&&z==0)||(x==3&&y==0&&z==1)||(x==3&&y==0&&z==2)||(x==3&&y==0&&z==3)||(x==3&&y==1&&z==0)||(x==3&&y==2&&z==0)||(x==3&&y==2&&z==2)||(x==3&&y==3&&z==0)
```

A system of polynomial equations and inequations:

```wolfram
Reduce[x^2+y^3==z&&x+2y==3z+1&&x y z≠0,{x,y,z},Modulus->7]
(* Output *)
(x==5&&y==2&&z==5)||(x==5&&y==6&&z==3)||(x==6&&y==4&&z==2)
```

Reduce a quantified polynomial system:

```wolfram
Reduce[ForAll[x,Exists[y,a x^2+b y^2-3y==1&&y≠0]],{a,b},Modulus->3]
(* Output *)
a==0&&b==1
```

#### Finite Field Domains

Univariate equations:

```wolfram
ℱ=FiniteField[53,4];
Reduce[x^5+ℱ[123]x==ℱ[234],x]
(* Output *)
x==<|interpretation -> FiniteFieldElement[FiniteField[53, 2, +, 38, #, +, 9, #, ^, 2, +, #, ^, 4, &, Polynomial], 10134132], index -> 4879932, shortIndex -> 4879932, indexShortened -> True, characteristic -> 53, shortCharacteristic -> 53, extensionDegree -> 4, field -> FiniteField[...], fieldDisplayed -> False|>
```

```wolfram
Reduce[x^7+2 x+3==0,x,ℱ]
(* Output *)
![image](img/image_001.png)
```

Systems of linear equations:

```wolfram
ℱ=FiniteField[71,2];
Reduce[ℱ[123]x+ℱ[234]y==ℱ[345]&&ℱ[321]x+ℱ[432]y==ℱ[543],{x,y}]
(* Output *)
x==<|interpretation -> FiniteFieldElement[FiniteField[71, 7, +, 69, #, +, #, ^, 2, &, Polynomial], 1114], index -> 1005, shortIndex -> 1005, indexShortened -> True, characteristic -> 71, shortCharacteristic -> 71, extensionDegree -> 2, field -> FiniteField[...], fieldDisplayed -> False|>&&y==<|interpretation -> FiniteFieldElement[FiniteField[71, 7, +, 69, #, +, #, ^, 2, &, Polynomial], 6157], index -> 4108, shortIndex -> 4108, indexShortened -> True, characteristic -> 71, shortCharacteristic -> 71, extensionDegree -> 2, field -> FiniteField[...], fieldDisplayed -> False|>
```

```wolfram
Reduce[ℱ[1234]x+ℱ[2345]y+ℱ[3456]z==ℱ[4567]&&ℱ[1]x+ℱ[2]y+ℱ[3]z==ℱ[4],{x,y,z}]
(* Output *)
![image](img/image_003.png)
```

Systems of polynomial equations:

```wolfram
ℱ=FiniteField[7,5];
Reduce[x^2+y^2==3&&x^5+y^5==5,{x,y},ℱ]
(* Output *)
![image](img/image_005.png)
```

```wolfram
Reduce[ℱ[123]x^2+ℱ[234]y^3+ℱ[345]z^4==ℱ[456]&&ℱ[21]x+ℱ[32]y^2+ℱ[43]z^3==ℱ[54]&&x y z==ℱ[1],{x,y,z}]
(* Output *)
![image](img/image_007.png)
```

Systems involving quantifiers:

```wolfram
ℱ=FiniteField[2,5];
Reduce[Exists[z,ℱ[1]x+ℱ[3]y+ℱ[5]z==ℱ[7]&&ℱ[21]x+ℱ[23]y+ℱ[25]z==ℱ[27]],{x,y}]
(* Output *)
![image](img/image_009.png)
```

```wolfram
Reduce[Exists[{y,z},ℱ[1]x^2+ℱ[2]y^3+ℱ[3]z^4==ℱ[4]&&ℱ[5]x^4+ℱ[6]y^3+ℱ[7]z^2==ℱ[8]&&x y z!=ℱ[0]],x]
(* Output *)
![image](img/image_011.png)
```

#### Mixed Domains

Mixed real and complex variables:

```wolfram
Reduce[x^2+y^2==1&&Element[x,Reals],{x,y}]
(* Output *)
(x<-1&&(y==-ⅈ Sqrt[-1+x^2]||y==ⅈ Sqrt[-1+x^2]))||(x==-1&&y==0)||(-1<x<1&&(y==-Sqrt[1-x^2]||y==Sqrt[1-x^2]))||(x==1&&y==0)||(x>1&&(y==-ⅈ Sqrt[-1+x^2]||y==ⅈ Sqrt[-1+x^2]))
```

Find real values of $x$ and complex values of $y$ for which $x^{2}+ y^{2}$ is real and less than $1$:

```wolfram
Reduce[x^2+y^2<1 &&Element[x,Reals],{x,y},Complexes]
(* Output *)
(x<-1&&Re[y]==0&&(Im[y]<-Sqrt[-1+x^2]||Im[y]>Sqrt[-1+x^2]))||(x==-1&&Re[y]==0&&(Im[y]<0||Im[y]>0))||(-1<x<1&&((-Sqrt[1-x^2]<Re[y]<0&&Im[y]==0)||Re[y]==0||(0<Re[y]<Sqrt[1-x^2]&&Im[y]==0)))||(x==1&&Re[y]==0&&(Im[y]<0||Im[y]>0))||(x>1&&Re[y]==0&&(Im[y]<-Sqrt[-1+x^2]||Im[y]>Sqrt[-1+x^2]))
```

Reduce an inequality involving [Abs](https://reference.wolfram.com/language/ref/Abs.html)[*x*]:

```wolfram
Reduce[1<Abs[ (z-2)/(2z-1)]<2,z]
(* Output *)
(-1<Re[z]<0&&-Sqrt[1-Re[z]^2]<Im[z]<Sqrt[1-Re[z]^2])||(0<=Re[z]<=(4)/(5)&&(-Sqrt[1-Re[z]^2]<Im[z]<-(Sqrt[4 Re[z]-5 Re[z]^2])/(Sqrt[5])||(Sqrt[4 Re[z]-5 Re[z]^2])/(Sqrt[5])<Im[z]<Sqrt[1-Re[z]^2]))||((4)/(5)<Re[z]<1&&-Sqrt[1-Re[z]^2]<Im[z]<Sqrt[1-Re[z]^2])
```

Plot the solution set:

```wolfram
Block[{z=u+I v},RegionPlot[1<Abs[ (z-2)/(2z-1)]<2,{u,-1,1},{v,-1,1}]]
```

*([Graphics])*

Mixed integer and real variables:

```wolfram
Reduce[x^2+y^2<5&&Element[x,Integers],{x,y},Reals]
(* Output *)
(x==-2&&-1<y<1)||(x==-1&&-2<y<2)||(x==0&&-Sqrt[5]<y<Sqrt[5])||(x==1&&-2<y<2)||(x==2&&-1<y<1)
```

#### Geometric Regions

Constrain variables to basic geometric regions in 2D:

```wolfram
ℛ_1=Circle[];
ℛ_2=Line[{{-2,1},{1,-2}}];
```

```wolfram
Reduce[{x,y}∈ℛ_1&&{x,y}∈ℛ_2,{x,y}]
(* Output *)
(x==-1&&y==0)||(x==0&&y==-1)
```

Plot the solution:

```wolfram
Graphics[{{Blue,ℛ_1,ℛ_2},{Red,Point[{x,y}]//.{ToRules[%]}}}]
```

*([Graphics])*

Constrain variables to basic geometric regions in 3D:

```wolfram
Reduce[2 x+3 y-5 z==1&&2 x y==z^2&&{x,y,z}∈Sphere[],{x,y,z},Reals]
(* Output *)
((x==-(2)/(3)&&y==-(1)/(3))||(x==-(8)/(17)&&y==-(9)/(17))||(x==(1)/(51) (27-5 Sqrt[21])&&y==(5189)/(8131)+(1690 (27-5 Sqrt[21]))/(24393)+(41 (27-5 Sqrt[21])^2)/(829362)-(14 (27-5 Sqrt[21])^3)/(414681))||(x==(1)/(51) (27+5 Sqrt[21])&&y==(5189)/(8131)+(1690 (27+5 Sqrt[21]))/(24393)+(41 (27+5 Sqrt[21])^2)/(829362)-(14 (27+5 Sqrt[21])^3)/(414681)))&&z==(1)/(5) (-1+2 x+3 y)
```

Plot the solution:

```wolfram
Show[{ContourPlot3D[{2 x+3 y-5 z==1,2 x y==z^2},{x,-1.2,1.2},{y,-1.2,1.2},{z,-1.2,1.2},Mesh->None,ContourStyle->Opacity[0.5]],Graphics3D[{{Opacity[0.5],Green,Sphere[]},{PointSize[Large],Red,Point[{x,y,z}//.{ToRules[%]}]}}]}]
```

*([Graphics3D])*

Project a 3D region onto the $x$-$y$ plane:

```wolfram
ℛ=Cone[{{0,0,0},{1,1,1}},2];
```

```wolfram
Reduce[∃_z{x,y,z}∈ℛ,{x,y},Reals]
(* Output *)
(x==-2 Sqrt[(2)/(3)]&&y==-(x)/(2)-(1)/(2) Sqrt[8-3 x^2])||(-2 Sqrt[(2)/(3)]<x<=(1)/(3) (2-Sqrt[6])&&(y==-(x)/(2)-(1)/(2) Sqrt[8-3 x^2]||-(x)/(2)-(1)/(2) Sqrt[8-3 x^2]<y<-(x)/(2)+(1)/(2) Sqrt[8-3 x^2]||y==-(x)/(2)+(1)/(2) Sqrt[8-3 x^2]))||((1)/(3) (2-Sqrt[6])<x<(1)/(3) (2+Sqrt[6])&&(y==-(x)/(2)-(1)/(2) Sqrt[8-3 x^2]||-(x)/(2)-(1)/(2) Sqrt[8-3 x^2]<y<=-(x)/(2)+(1)/(2) Sqrt[8-3 x^2]||-(x)/(2)+(1)/(2) Sqrt[8-3 x^2]<y<(1)/(5) (12-7 x)-(2)/(5) Sqrt[6] Sqrt[1-2 x+x^2]||y==(1)/(5) (12-7 x)-(2)/(5) Sqrt[6] Sqrt[1-2 x+x^2]))||((1)/(3) (2+Sqrt[6])<=x<2 Sqrt[(2)/(3)]&&(y==-(x)/(2)-(1)/(2) Sqrt[8-3 x^2]||-(x)/(2)-(1)/(2) Sqrt[8-3 x^2]<y<-(x)/(2)+(1)/(2) Sqrt[8-3 x^2]||y==-(x)/(2)+(1)/(2) Sqrt[8-3 x^2]))||(x==2 Sqrt[(2)/(3)]&&y==-(x)/(2)-(1)/(2) Sqrt[8-3 x^2])
```

Plot the projection:

```wolfram
RegionPlot[%,{x,-2,2},{y,-2,2}]
```

*([Graphics])*

An implicitly defined region:

```wolfram
ℛ=ImplicitRegion[a+2 b-3 c>=1&&a b c==7,{a,b,c}];
```

```wolfram
Reduce[x^2+y z==1&&x+2 y<3 z+7&&{x,y,z}∈ℛ,{x,y,z},Reals]
(* Output *)
x==Root&&(7-x)/(4)-(1)/(4) Sqrt[(168+49 x-14 x^2+x^3)/(x)]<y<(7-x)/(4)+(1)/(4) Sqrt[(168+49 x-14 x^2+x^3)/(x)]&&z==(7)/(x y)
```

A parametrically defined region:

```wolfram
ℛ=ParametricRegion[{s^2,t^2,s t},{s,t}];
```

```wolfram
Reduce[x^2+y^2+z^2==1&&{x,y,z}∈ℛ,{x,y,z},Reals]
(* Output *)
(x==0&&y==1&&z==0)||(0<x<1&&y==-(x)/(2)+(1)/(2) Sqrt[4-3 x^2]&&(z==-Sqrt[1-x^2-y^2]||z==Sqrt[1-x^2-y^2]))||(x==1&&y==0&&z==0)
```

Derived regions:

```wolfram
ℛ=RegionIntersection[Disk[{0,0},2],Disk[{1,1},2]];
```

```wolfram
Reduce[x^2==x y+1&&{x,y}∈ℛ,{x,y},Reals]
(* Output *)
(Root<=x<=Root||Root<=x<=Root)&&y==-(1-x^2)/(x)
```

The solution of $x^{2}=x y+1$ restricted to the intersection:

```wolfram
Show[{DiscretizeRegion[ℛ],ContourPlot[x^2==x y+1,{x,-2,3},{y,-2,3}],
Plot@@{Last[%[[2]]],{x,First[#],Last[#]},PlotStyle->Red}&/@(List@@%[[1]])}]
```

*([Graphics])*

Eliminate quantifiers over a Cartesian product of regions:

```wolfram
ℛ=RegionProduct[Circle[],Circle[]];
```

```wolfram
Reduce[∃_{a,b}(4 x y a b==1&&{x,a,y,b}∈ℛ),{x,y}]
(* Output *)
(x==-(1)/(Sqrt[2])&&(y==-(1)/(Sqrt[2])||y==(1)/(Sqrt[2])))||(x==(1)/(Sqrt[2])&&(y==-(1)/(Sqrt[2])||y==(1)/(Sqrt[2])))
```

Regions dependent on parameters:

```wolfram
ℛ_1=InfiniteLine[{{2,0},{0,t}}];
ℛ_2=Circle[];
```

```wolfram
Reduce[{x,y}∈ℛ_1&&{x,y}∈ℛ_2,{t,x,y},Reals]
(* Output *)
((t==-(2)/(Sqrt[3])&&x==(1)/(2))||(-(2)/(Sqrt[3])<t<(2)/(Sqrt[3])&&(x==(2 t^2)/(4+t^2)-2 Sqrt[-(-4+3 t^2)/((4+t^2)^2)]||x==(2 t^2)/(4+t^2)+2 Sqrt[-(-4+3 t^2)/((4+t^2)^2)]))||(t==(2)/(Sqrt[3])&&x==(1)/(2)))&&y==(1)/(2) (2 t-t x)
```

A condition for $\mathcal{R}_{1}\subseteq \mathcal{R}_{2}$:

```wolfram
ℛ_1=Disk[{a,b}];
ℛ_2=Triangle[{{0,0},{0,5},{7,0}}];
Reduce[Subscript[∀, {x,y},{x,y}∈ℛ_1]{x,y}∈ℛ_2,{a,b},Reals]
(* Output *)
(1<=a<(1)/(5) (28-Sqrt[74])&&1<=b<=-(Sqrt[74])/(7)-(5)/(7) (-7+a))||(a==(1)/(5) (28-Sqrt[74])&&b==1)
```

Use $x \in \mathcal{R}$ to specify that $x$ is a vector in $\mathbb{R}^{2}$:

```wolfram
ℛ=RegionIntersection[Circle[],Line[{{-2,-1},{1,2}}]];
```

```wolfram
Reduce[x∈ℛ,x]
(* Output *)
(x==-1&&x==0)||(x==0&&x==1)
```

In this case $x$ is a vector in $\mathbb{R}^{3}$:

```wolfram
ℛ=Sphere[];
```

```wolfram
Reduce[x.{1,2,3}==0&&x∈ℛ,x]
(* Output *)
((x==-Sqrt[(13)/(14)]&&x==Sqrt[(2)/(91)])||(-Sqrt[(13)/(14)]<x<Sqrt[(13)/(14)]&&(x==-(2 x)/(13)-(3)/(13) Sqrt[13-14 x^2]||x==-(2 x)/(13)+(3)/(13) Sqrt[13-14 x^2]))||(x==Sqrt[(13)/(14)]&&x==-Sqrt[(2)/(91)]))&&x==(1)/(3) (-x-2 x)
```

### Options

#### Backsubstitution

Since `y` appears after `x` in the variable list, [Reduce](https://reference.wolfram.com/language/ref/Reduce.html) may use `x` to express the solution for `y`:

```wolfram
Reduce[x^2+y^2==1&&x^2-y==2,{x,y}]
(* Output *)
(x==-Sqrt[(3)/(2)-(ⅈ Sqrt[3])/(2)]||x==Sqrt[(3)/(2)-(ⅈ Sqrt[3])/(2)]||x==-Sqrt[(3)/(2)+(ⅈ Sqrt[3])/(2)]||x==Sqrt[(3)/(2)+(ⅈ Sqrt[3])/(2)])&&y==-2+x^2
```

With `Backsubstitution->[True](https://reference.wolfram.com/language/ref/True.html)`, [Reduce](https://reference.wolfram.com/language/ref/Reduce.html) gives explicit numeric values for `y`:

```wolfram
Reduce[x^2+y^2==1&&x^2-y==2,{x,y},Backsubstitution->True]
(* Output *)
(x==-Sqrt[(1)/(2) (3-ⅈ Sqrt[3])]&&y==-(1)/(2) ⅈ (-ⅈ+Sqrt[3]))||(x==Sqrt[(1)/(2) (3-ⅈ Sqrt[3])]&&y==-(1)/(2) ⅈ (-ⅈ+Sqrt[3]))||(x==-Sqrt[(1)/(2) (3+ⅈ Sqrt[3])]&&y==(1)/(2) ⅈ (ⅈ+Sqrt[3]))||(x==Sqrt[(1)/(2) (3+ⅈ Sqrt[3])]&&y==(1)/(2) ⅈ (ⅈ+Sqrt[3]))
```

#### Cubics

By default, [Reduce](https://reference.wolfram.com/language/ref/Reduce.html) does not use general formulas for solving cubics in radicals:

```wolfram
Reduce[x^3+2 x^2+3 x+4==0,x]
(* Output *)
x==Root||x==Root||x==Root
```

With [Cubics](https://reference.wolfram.com/language/ref/Cubics.html)->[True](https://reference.wolfram.com/language/ref/True.html), [Reduce](https://reference.wolfram.com/language/ref/Reduce.html) solves all cubics in terms of radicals:

```wolfram
Reduce[x^3+2 x^2+3 x+4==0,x,Cubics->True]
(* Output *)
x==(1)/(3) (-2-(5^(2/3))/((-7+3 Sqrt[6])^(1/3))+(5 (-7+3 Sqrt[6]))^(1/3))||x==-(2)/(3)+(5^(2/3) (1+ⅈ Sqrt[3]))/(6 (-7+3 Sqrt[6])^(1/3))-(1)/(6) (1-ⅈ Sqrt[3]) (5 (-7+3 Sqrt[6]))^(1/3)||x==-(2)/(3)+(5^(2/3) (1-ⅈ Sqrt[3]))/(6 (-7+3 Sqrt[6])^(1/3))-(1)/(6) (1+ⅈ Sqrt[3]) (5 (-7+3 Sqrt[6]))^(1/3)
```

#### GeneratedParameters

[Reduce](https://reference.wolfram.com/language/ref/Reduce.html) may introduce new parameters to represent the solution:

```wolfram
Reduce[Cos[x]==1/2,x]
(* Output *)
1∈Integers&&(x==-(π)/(3)+2 π 1||x==(π)/(3)+2 π 1)
```

Use [GeneratedParameters](https://reference.wolfram.com/language/ref/GeneratedParameters.html) to control how the parameters are generated:

```wolfram
Reduce[Cos[x]==1/2,x,GeneratedParameters->(k_#&)]
(* Output *)
k_1∈Integers&&(x==-(π)/(3)+2 π k_1||x==(π)/(3)+2 π k_1)
```

#### Method

This locally sets system options in ["InequalitySolvingOptions"](https://reference.wolfram.com/language/tutorial/RealPolynomialSystems.html#343227045) and ["ReduceOptions"](https://reference.wolfram.com/language/tutorial/RealPolynomialSystems.html#227757976) groups:

```wolfram
Reduce[y^5-3 y^2+x^5+4==0&&y^2+y-x^5-3 x-11==0,{x,y},Reals,Backsubstitution->True,Method->{"AlgebraicNumberOutput"->False}]
(* Output *)
x==<|icon -> Root, small -> "-1.30", approx -> -1.297917436600305, interp -> Root[150524, +, 210474, #, +, 117189, #, ^, 2, +, 32427, #, ^, 3, +, 4455, #, ^, 4, +, 71123, #, ^, 5, +, 78489, #, ^, 6, +, 32472, #, ^, 7, +, 5940, #, ^, 8, +, 405, #, ^, 9, +, 13141, #, ^, 10, +, 10839, #, ^, 11, +, 2970, #, ^, 12, +, 270, #, ^, 13, +, 1206, #, ^, 15, +, 660, #, ^, 16, +, 90, #, ^, 17, +, 55, #, ^, 20, +, 15, #, ^, 21, +, #, ^, 25, &, 1, 0], head -> Root, big -> 150524+210474 #1+117189 #1^2+32427 #1^3+4455 #1^4+71123 #1^5+78489 #1^6+32472 #1^7+5940 #1^8+405 #1^9+13141 #1^10+10839 #1^11+2970 #1^12+270 #1^13+1206 #1^15+660 #1^16+90 #1^17+55 #1^20+15 #1^21+#1^25&, degree -> 25, shortDegree -> 25, number -> 1|>&&y==(1)/(25907338)(9674965580-9375045027 <|icon -> Root, small -> "-1.30", approx -> -1.297917436600305, interp -> Root[150524, +, 210474, #, +, 117189, #, ^, 2, +, 32427, #, ^, 3, +, 4455, #, ^, 4, +, 71123, #, ^, 5, +, 78489, #, ^, 6, +, 32472, #, ^, 7, +, 5940, #, ^, 8, +, 405, #, ^, 9, +, 13141, #, ^, 10, +, 10839, #, ^, 11, +, 2970, #, ^, 12, +, 270, #, ^, 13, +, 1206, #, ^, 15, +, 660, #, ^, 16, +, 90, #, ^, 17, +, 55, #, ^, 20, +, 15, #, ^, 21, +, #, ^, 25, &, 1, 0], head -> Root, big -> 150524+210474 #1+117189 #1^2+32427 #1^3+4455 #1^4+71123 #1^5+78489 #1^6+32472 #1^7+5940 #1^8+405 #1^9+13141 #1^10+10839 #1^11+2970 #1^12+270 #1^13+1206 #1^15+660 #1^16+90 #1^17+55 #1^20+15 #1^21+#1^25&, degree -> 25, shortDegree -> 25, number -> 1|>-9342236349 <|icon -> Root, small -> "-1.30", approx -> -1.297917436600305, interp -> Root[150524, +, 210474, #, +, 117189, #, ^, 2, +, 32427, #, ^, 3, +, 4455, #, ^, 4, +, 71123, #, ^, 5, +, 78489, #, ^, 6, +, 32472, #, ^, 7, +, 5940, #, ^, 8, +, 405, #, ^, 9, +, 13141, #, ^, 10, +, 10839, #, ^, 11, +, 2970, #, ^, 12, +, 270, #, ^, 13, +, 1206, #, ^, 15, +, 660, #, ^, 16, +, 90, #, ^, 17, +, 55, #, ^, 20, +, 15, #, ^, 21, +, #, ^, 25, &, 1, 0], head -> Root, big -> 150524+210474 #1+117189 #1^2+32427 #1^3+4455 #1^4+71123 #1^5+78489 #1^6+32472 #1^7+5940 #1^8+405 #1^9+13141 #1^10+10839 #1^11+2970 #1^12+270 #1^13+1206 #1^15+660 #1^16+90 #1^17+55 #1^20+15 #1^21+#1^25&, degree -> 25, shortDegree -> 25, number -> 1|>^2+1291893057 <|icon -> Root, small -> "-1.30", approx -> -1.297917436600305, interp -> Root[150524, +, 210474, #, +, 117189, #, ^, 2, +, 32427, #, ^, 3, +, 4455, #, ^, 4, +, 71123, #, ^, 5, +, 78489, #, ^, 6, +, 32472, #, ^, 7, +, 5940, #, ^, 8, +, 405, #, ^, 9, +, 13141, #, ^, 10, +, 10839, #, ^, 11, +, 2970, #, ^, 12, +, 270, #, ^, 13, +, 1206, #, ^, 15, +, 660, #, ^, 16, +, 90, #, ^, 17, +, 55, #, ^, 20, +, 15, #, ^, 21, +, #, ^, 25, &, 1, 0], head -> Root, big -> 150524+210474 #1+117189 #1^2+32427 #1^3+4455 #1^4+71123 #1^5+78489 #1^6+32472 #1^7+5940 #1^8+405 #1^9+13141 #1^10+10839 #1^11+2970 #1^12+270 #1^13+1206 #1^15+660 #1^16+90 #1^17+55 #1^20+15 #1^21+#1^25&, degree -> 25, shortDegree -> 25, number -> 1|>^3-398647575 <|icon -> Root, small -> "-1.30", approx -> -1.297917436600305, interp -> Root[150524, +, 210474, #, +, 117189, #, ^, 2, +, 32427, #, ^, 3, +, 4455, #, ^, 4, +, 71123, #, ^, 5, +, 78489, #, ^, 6, +, 32472, #, ^, 7, +, 5940, #, ^, 8, +, 405, #, ^, 9, +, 13141, #, ^, 10, +, 10839, #, ^, 11, +, 2970, #, ^, 12, +, 270, #, ^, 13, +, 1206, #, ^, 15, +, 660, #, ^, 16, +, 90, #, ^, 17, +, 55, #, ^, 20, +, 15, #, ^, 21, +, #, ^, 25, &, 1, 0], head -> Root, big -> 150524+210474 #1+117189 #1^2+32427 #1^3+4455 #1^4+71123 #1^5+78489 #1^6+32472 #1^7+5940 #1^8+405 #1^9+13141 #1^10+10839 #1^11+2970 #1^12+270 #1^13+1206 #1^15+660 #1^16+90 #1^17+55 #1^20+15 #1^21+#1^25&, degree -> 25, shortDegree -> 25, number -> 1|>^4+1033763855 <|icon -> Root, small -> "-1.30", approx -> -1.297917436600305, interp -> Root[150524, +, 210474, #, +, 117189, #, ^, 2, +, 32427, #, ^, 3, +, 4455, #, ^, 4, +, 71123, #, ^, 5, +, 78489, #, ^, 6, +, 32472, #, ^, 7, +, 5940, #, ^, 8, +, 405, #, ^, 9, +, 13141, #, ^, 10, +, 10839, #, ^, 11, +, 2970, #, ^, 12, +, 270, #, ^, 13, +, 1206, #, ^, 15, +, 660, #, ^, 16, +, 90, #, ^, 17, +, 55, #, ^, 20, +, 15, #, ^, 21, +, #, ^, 25, &, 1, 0], head -> Root, big -> 150524+210474 #1+117189 #1^2+32427 #1^3+4455 #1^4+71123 #1^5+78489 #1^6+32472 #1^7+5940 #1^8+405 #1^9+13141 #1^10+10839 #1^11+2970 #1^12+270 #1^13+1206 #1^15+660 #1^16+90 #1^17+55 #1^20+15 #1^21+#1^25&, degree -> 25, shortDegree -> 25, number -> 1|>^5-6166452360 <|icon -> Root, small -> "-1.30", approx -> -1.297917436600305, interp -> Root[150524, +, 210474, #, +, 117189, #, ^, 2, +, 32427, #, ^, 3, +, 4455, #, ^, 4, +, 71123, #, ^, 5, +, 78489, #, ^, 6, +, 32472, #, ^, 7, +, 5940, #, ^, 8, +, 405, #, ^, 9, +, 13141, #, ^, 10, +, 10839, #, ^, 11, +, 2970, #, ^, 12, +, 270, #, ^, 13, +, 1206, #, ^, 15, +, 660, #, ^, 16, +, 90, #, ^, 17, +, 55, #, ^, 20, +, 15, #, ^, 21, +, #, ^, 25, &, 1, 0], head -> Root, big -> 150524+210474 #1+117189 #1^2+32427 #1^3+4455 #1^4+71123 #1^5+78489 #1^6+32472 #1^7+5940 #1^8+405 #1^9+13141 #1^10+10839 #1^11+2970 #1^12+270 #1^13+1206 #1^15+660 #1^16+90 #1^17+55 #1^20+15 #1^21+#1^25&, degree -> 25, shortDegree -> 25, number -> 1|>^6-2368888182 <|icon -> Root, small -> "-1.30", approx -> -1.297917436600305, interp -> Root[150524, +, 210474, #, +, 117189, #, ^, 2, +, 32427, #, ^, 3, +, 4455, #, ^, 4, +, 71123, #, ^, 5, +, 78489, #, ^, 6, +, 32472, #, ^, 7, +, 5940, #, ^, 8, +, 405, #, ^, 9, +, 13141, #, ^, 10, +, 10839, #, ^, 11, +, 2970, #, ^, 12, +, 270, #, ^, 13, +, 1206, #, ^, 15, +, 660, #, ^, 16, +, 90, #, ^, 17, +, 55, #, ^, 20, +, 15, #, ^, 21, +, #, ^, 25, &, 1, 0], head -> Root, big -> 150524+210474 #1+117189 #1^2+32427 #1^3+4455 #1^4+71123 #1^5+78489 #1^6+32472 #1^7+5940 #1^8+405 #1^9+13141 #1^10+10839 #1^11+2970 #1^12+270 #1^13+1206 #1^15+660 #1^16+90 #1^17+55 #1^20+15 #1^21+#1^25&, degree -> 25, shortDegree -> 25, number -> 1|>^7+1136838537 <|icon -> Root, small -> "-1.30", approx -> -1.297917436600305, interp -> Root[150524, +, 210474, #, +, 117189, #, ^, 2, +, 32427, #, ^, 3, +, 4455, #, ^, 4, +, 71123, #, ^, 5, +, 78489, #, ^, 6, +, 32472, #, ^, 7, +, 5940, #, ^, 8, +, 405, #, ^, 9, +, 13141, #, ^, 10, +, 10839, #, ^, 11, +, 2970, #, ^, 12, +, 270, #, ^, 13, +, 1206, #, ^, 15, +, 660, #, ^, 16, +, 90, #, ^, 17, +, 55, #, ^, 20, +, 15, #, ^, 21, +, #, ^, 25, &, 1, 0], head -> Root, big -> 150524+210474 #1+117189 #1^2+32427 #1^3+4455 #1^4+71123 #1^5+78489 #1^6+32472 #1^7+5940 #1^8+405 #1^9+13141 #1^10+10839 #1^11+2970 #1^12+270 #1^13+1206 #1^15+660 #1^16+90 #1^17+55 #1^20+15 #1^21+#1^25&, degree -> 25, shortDegree -> 25, number -> 1|>^8-519316353 <|icon -> Root, small -> "-1.30", approx -> -1.297917436600305, interp -> Root[150524, +, 210474, #, +, 117189, #, ^, 2, +, 32427, #, ^, 3, +, 4455, #, ^, 4, +, 71123, #, ^, 5, +, 78489, #, ^, 6, +, 32472, #, ^, 7, +, 5940, #, ^, 8, +, 405, #, ^, 9, +, 13141, #, ^, 10, +, 10839, #, ^, 11, +, 2970, #, ^, 12, +, 270, #, ^, 13, +, 1206, #, ^, 15, +, 660, #, ^, 16, +, 90, #, ^, 17, +, 55, #, ^, 20, +, 15, #, ^, 21, +, #, ^, 25, &, 1, 0], head -> Root, big -> 150524+210474 #1+117189 #1^2+32427 #1^3+4455 #1^4+71123 #1^5+78489 #1^6+32472 #1^7+5940 #1^8+405 #1^9+13141 #1^10+10839 #1^11+2970 #1^12+270 #1^13+1206 #1^15+660 #1^16+90 #1^17+55 #1^20+15 #1^21+#1^25&, degree -> 25, shortDegree -> 25, number -> 1|>^9-349540011 <|icon -> Root, small -> "-1.30", approx -> -1.297917436600305, interp -> Root[150524, +, 210474, #, +, 117189, #, ^, 2, +, 32427, #, ^, 3, +, 4455, #, ^, 4, +, 71123, #, ^, 5, +, 78489, #, ^, 6, +, 32472, #, ^, 7, +, 5940, #, ^, 8, +, 405, #, ^, 9, +, 13141, #, ^, 10, +, 10839, #, ^, 11, +, 2970, #, ^, 12, +, 270, #, ^, 13, +, 1206, #, ^, 15, +, 660, #, ^, 16, +, 90, #, ^, 17, +, 55, #, ^, 20, +, 15, #, ^, 21, +, #, ^, 25, &, 1, 0], head -> Root, big -> 150524+210474 #1+117189 #1^2+32427 #1^3+4455 #1^4+71123 #1^5+78489 #1^6+32472 #1^7+5940 #1^8+405 #1^9+13141 #1^10+10839 #1^11+2970 #1^12+270 #1^13+1206 #1^15+660 #1^16+90 #1^17+55 #1^20+15 #1^21+#1^25&, degree -> 25, shortDegree -> 25, number -> 1|>^10-1095160986 <|icon -> Root, small -> "-1.30", approx -> -1.297917436600305, interp -> Root[150524, +, 210474, #, +, 117189, #, ^, 2, +, 32427, #, ^, 3, +, 4455, #, ^, 4, +, 71123, #, ^, 5, +, 78489, #, ^, 6, +, 32472, #, ^, 7, +, 5940, #, ^, 8, +, 405, #, ^, 9, +, 13141, #, ^, 10, +, 10839, #, ^, 11, +, 2970, #, ^, 12, +, 270, #, ^, 13, +, 1206, #, ^, 15, +, 660, #, ^, 16, +, 90, #, ^, 17, +, 55, #, ^, 20, +, 15, #, ^, 21, +, #, ^, 25, &, 1, 0], head -> Root, big -> 150524+210474 #1+117189 #1^2+32427 #1^3+4455 #1^4+71123 #1^5+78489 #1^6+32472 #1^7+5940 #1^8+405 #1^9+13141 #1^10+10839 #1^11+2970 #1^12+270 #1^13+1206 #1^15+660 #1^16+90 #1^17+55 #1^20+15 #1^21+#1^25&, degree -> 25, shortDegree -> 25, number -> 1|>^11-65551410 <|icon -> Root, small -> "-1.30", approx -> -1.297917436600305, interp -> Root[150524, +, 210474, #, +, 117189, #, ^, 2, +, 32427, #, ^, 3, +, 4455, #, ^, 4, +, 71123, #, ^, 5, +, 78489, #, ^, 6, +, 32472, #, ^, 7, +, 5940, #, ^, 8, +, 405, #, ^, 9, +, 13141, #, ^, 10, +, 10839, #, ^, 11, +, 2970, #, ^, 12, +, 270, #, ^, 13, +, 1206, #, ^, 15, +, 660, #, ^, 16, +, 90, #, ^, 17, +, 55, #, ^, 20, +, 15, #, ^, 21, +, #, ^, 25, &, 1, 0], head -> Root, big -> 150524+210474 #1+117189 #1^2+32427 #1^3+4455 #1^4+71123 #1^5+78489 #1^6+32472 #1^7+5940 #1^8+405 #1^9+13141 #1^10+10839 #1^11+2970 #1^12+270 #1^13+1206 #1^15+660 #1^16+90 #1^17+55 #1^20+15 #1^21+#1^25&, degree -> 25, shortDegree -> 25, number -> 1|>^12+175232808 <|icon -> Root, small -> "-1.30", approx -> -1.297917436600305, interp -> Root[150524, +, 210474, #, +, 117189, #, ^, 2, +, 32427, #, ^, 3, +, 4455, #, ^, 4, +, 71123, #, ^, 5, +, 78489, #, ^, 6, +, 32472, #, ^, 7, +, 5940, #, ^, 8, +, 405, #, ^, 9, +, 13141, #, ^, 10, +, 10839, #, ^, 11, +, 2970, #, ^, 12, +, 270, #, ^, 13, +, 1206, #, ^, 15, +, 660, #, ^, 16, +, 90, #, ^, 17, +, 55, #, ^, 20, +, 15, #, ^, 21, +, #, ^, 25, &, 1, 0], head -> Root, big -> 150524+210474 #1+117189 #1^2+32427 #1^3+4455 #1^4+71123 #1^5+78489 #1^6+32472 #1^7+5940 #1^8+405 #1^9+13141 #1^10+10839 #1^11+2970 #1^12+270 #1^13+1206 #1^15+660 #1^16+90 #1^17+55 #1^20+15 #1^21+#1^25&, degree -> 25, shortDegree -> 25, number -> 1|>^13-126571896 <|icon -> Root, small -> "-1.30", approx -> -1.297917436600305, interp -> Root[150524, +, 210474, #, +, 117189, #, ^, 2, +, 32427, #, ^, 3, +, 4455, #, ^, 4, +, 71123, #, ^, 5, +, 78489, #, ^, 6, +, 32472, #, ^, 7, +, 5940, #, ^, 8, +, 405, #, ^, 9, +, 13141, #, ^, 10, +, 10839, #, ^, 11, +, 2970, #, ^, 12, +, 270, #, ^, 13, +, 1206, #, ^, 15, +, 660, #, ^, 16, +, 90, #, ^, 17, +, 55, #, ^, 20, +, 15, #, ^, 21, +, #, ^, 25, &, 1, 0], head -> Root, big -> 150524+210474 #1+117189 #1^2+32427 #1^3+4455 #1^4+71123 #1^5+78489 #1^6+32472 #1^7+5940 #1^8+405 #1^9+13141 #1^10+10839 #1^11+2970 #1^12+270 #1^13+1206 #1^15+660 #1^16+90 #1^17+55 #1^20+15 #1^21+#1^25&, degree -> 25, shortDegree -> 25, number -> 1|>^14-52601296 <|icon -> Root, small -> "-1.30", approx -> -1.297917436600305, interp -> Root[150524, +, 210474, #, +, 117189, #, ^, 2, +, 32427, #, ^, 3, +, 4455, #, ^, 4, +, 71123, #, ^, 5, +, 78489, #, ^, 6, +, 32472, #, ^, 7, +, 5940, #, ^, 8, +, 405, #, ^, 9, +, 13141, #, ^, 10, +, 10839, #, ^, 11, +, 2970, #, ^, 12, +, 270, #, ^, 13, +, 1206, #, ^, 15, +, 660, #, ^, 16, +, 90, #, ^, 17, +, 55, #, ^, 20, +, 15, #, ^, 21, +, #, ^, 25, &, 1, 0], head -> Root, big -> 150524+210474 #1+117189 #1^2+32427 #1^3+4455 #1^4+71123 #1^5+78489 #1^6+32472 #1^7+5940 #1^8+405 #1^9+13141 #1^10+10839 #1^11+2970 #1^12+270 #1^13+1206 #1^15+660 #1^16+90 #1^17+55 #1^20+15 #1^21+#1^25&, degree -> 25, shortDegree -> 25, number -> 1|>^15-72464502 <|icon -> Root, small -> "-1.30", approx -> -1.297917436600305, interp -> Root[150524, +, 210474, #, +, 117189, #, ^, 2, +, 32427, #, ^, 3, +, 4455, #, ^, 4, +, 71123, #, ^, 5, +, 78489, #, ^, 6, +, 32472, #, ^, 7, +, 5940, #, ^, 8, +, 405, #, ^, 9, +, 13141, #, ^, 10, +, 10839, #, ^, 11, +, 2970, #, ^, 12, +, 270, #, ^, 13, +, 1206, #, ^, 15, +, 660, #, ^, 16, +, 90, #, ^, 17, +, 55, #, ^, 20, +, 15, #, ^, 21, +, #, ^, 25, &, 1, 0], head -> Root, big -> 150524+210474 #1+117189 #1^2+32427 #1^3+4455 #1^4+71123 #1^5+78489 #1^6+32472 #1^7+5940 #1^8+405 #1^9+13141 #1^10+10839 #1^11+2970 #1^12+270 #1^13+1206 #1^15+660 #1^16+90 #1^17+55 #1^20+15 #1^21+#1^25&, degree -> 25, shortDegree -> 25, number -> 1|>^16+15761592 <|icon -> Root, small -> "-1.30", approx -> -1.297917436600305, interp -> Root[150524, +, 210474, #, +, 117189, #, ^, 2, +, 32427, #, ^, 3, +, 4455, #, ^, 4, +, 71123, #, ^, 5, +, 78489, #, ^, 6, +, 32472, #, ^, 7, +, 5940, #, ^, 8, +, 405, #, ^, 9, +, 13141, #, ^, 10, +, 10839, #, ^, 11, +, 2970, #, ^, 12, +, 270, #, ^, 13, +, 1206, #, ^, 15, +, 660, #, ^, 16, +, 90, #, ^, 17, +, 55, #, ^, 20, +, 15, #, ^, 21, +, #, ^, 25, &, 1, 0], head -> Root, big -> 150524+210474 #1+117189 #1^2+32427 #1^3+4455 #1^4+71123 #1^5+78489 #1^6+32472 #1^7+5940 #1^8+405 #1^9+13141 #1^10+10839 #1^11+2970 #1^12+270 #1^13+1206 #1^15+660 #1^16+90 #1^17+55 #1^20+15 #1^21+#1^25&, degree -> 25, shortDegree -> 25, number -> 1|>^17+6065172 <|icon -> Root, small -> "-1.30", approx -> -1.297917436600305, interp -> Root[150524, +, 210474, #, +, 117189, #, ^, 2, +, 32427, #, ^, 3, +, 4455, #, ^, 4, +, 71123, #, ^, 5, +, 78489, #, ^, 6, +, 32472, #, ^, 7, +, 5940, #, ^, 8, +, 405, #, ^, 9, +, 13141, #, ^, 10, +, 10839, #, ^, 11, +, 2970, #, ^, 12, +, 270, #, ^, 13, +, 1206, #, ^, 15, +, 660, #, ^, 16, +, 90, #, ^, 17, +, 55, #, ^, 20, +, 15, #, ^, 21, +, #, ^, 25, &, 1, 0], head -> Root, big -> 150524+210474 #1+117189 #1^2+32427 #1^3+4455 #1^4+71123 #1^5+78489 #1^6+32472 #1^7+5940 #1^8+405 #1^9+13141 #1^10+10839 #1^11+2970 #1^12+270 #1^13+1206 #1^15+660 #1^16+90 #1^17+55 #1^20+15 #1^21+#1^25&, degree -> 25, shortDegree -> 25, number -> 1|>^18-9615105 <|icon -> Root, small -> "-1.30", approx -> -1.297917436600305, interp -> Root[150524, +, 210474, #, +, 117189, #, ^, 2, +, 32427, #, ^, 3, +, 4455, #, ^, 4, +, 71123, #, ^, 5, +, 78489, #, ^, 6, +, 32472, #, ^, 7, +, 5940, #, ^, 8, +, 405, #, ^, 9, +, 13141, #, ^, 10, +, 10839, #, ^, 11, +, 2970, #, ^, 12, +, 270, #, ^, 13, +, 1206, #, ^, 15, +, 660, #, ^, 16, +, 90, #, ^, 17, +, 55, #, ^, 20, +, 15, #, ^, 21, +, #, ^, 25, &, 1, 0], head -> Root, big -> 150524+210474 #1+117189 #1^2+32427 #1^3+4455 #1^4+71123 #1^5+78489 #1^6+32472 #1^7+5940 #1^8+405 #1^9+13141 #1^10+10839 #1^11+2970 #1^12+270 #1^13+1206 #1^15+660 #1^16+90 #1^17+55 #1^20+15 #1^21+#1^25&, degree -> 25, shortDegree -> 25, number -> 1|>^19-1649643 <|icon -> Root, small -> "-1.30", approx -> -1.297917436600305, interp -> Root[150524, +, 210474, #, +, 117189, #, ^, 2, +, 32427, #, ^, 3, +, 4455, #, ^, 4, +, 71123, #, ^, 5, +, 78489, #, ^, 6, +, 32472, #, ^, 7, +, 5940, #, ^, 8, +, 405, #, ^, 9, +, 13141, #, ^, 10, +, 10839, #, ^, 11, +, 2970, #, ^, 12, +, 270, #, ^, 13, +, 1206, #, ^, 15, +, 660, #, ^, 16, +, 90, #, ^, 17, +, 55, #, ^, 20, +, 15, #, ^, 21, +, #, ^, 25, &, 1, 0], head -> Root, big -> 150524+210474 #1+117189 #1^2+32427 #1^3+4455 #1^4+71123 #1^5+78489 #1^6+32472 #1^7+5940 #1^8+405 #1^9+13141 #1^10+10839 #1^11+2970 #1^12+270 #1^13+1206 #1^15+660 #1^16+90 #1^17+55 #1^20+15 #1^21+#1^25&, degree -> 25, shortDegree -> 25, number -> 1|>^20-1617714 <|icon -> Root, small -> "-1.30", approx -> -1.297917436600305, interp -> Root[150524, +, 210474, #, +, 117189, #, ^, 2, +, 32427, #, ^, 3, +, 4455, #, ^, 4, +, 71123, #, ^, 5, +, 78489, #, ^, 6, +, 32472, #, ^, 7, +, 5940, #, ^, 8, +, 405, #, ^, 9, +, 13141, #, ^, 10, +, 10839, #, ^, 11, +, 2970, #, ^, 12, +, 270, #, ^, 13, +, 1206, #, ^, 15, +, 660, #, ^, 16, +, 90, #, ^, 17, +, 55, #, ^, 20, +, 15, #, ^, 21, +, #, ^, 25, &, 1, 0], head -> Root, big -> 150524+210474 #1+117189 #1^2+32427 #1^3+4455 #1^4+71123 #1^5+78489 #1^6+32472 #1^7+5940 #1^8+405 #1^9+13141 #1^10+10839 #1^11+2970 #1^12+270 #1^13+1206 #1^15+660 #1^16+90 #1^17+55 #1^20+15 #1^21+#1^25&, degree -> 25, shortDegree -> 25, number -> 1|>^21+842562 <|icon -> Root, small -> "-1.30", approx -> -1.297917436600305, interp -> Root[150524, +, 210474, #, +, 117189, #, ^, 2, +, 32427, #, ^, 3, +, 4455, #, ^, 4, +, 71123, #, ^, 5, +, 78489, #, ^, 6, +, 32472, #, ^, 7, +, 5940, #, ^, 8, +, 405, #, ^, 9, +, 13141, #, ^, 10, +, 10839, #, ^, 11, +, 2970, #, ^, 12, +, 270, #, ^, 13, +, 1206, #, ^, 15, +, 660, #, ^, 16, +, 90, #, ^, 17, +, 55, #, ^, 20, +, 15, #, ^, 21, +, #, ^, 25, &, 1, 0], head -> Root, big -> 150524+210474 #1+117189 #1^2+32427 #1^3+4455 #1^4+71123 #1^5+78489 #1^6+32472 #1^7+5940 #1^8+405 #1^9+13141 #1^10+10839 #1^11+2970 #1^12+270 #1^13+1206 #1^15+660 #1^16+90 #1^17+55 #1^20+15 #1^21+#1^25&, degree -> 25, shortDegree -> 25, number -> 1|>^22-111942 <|icon -> Root, small -> "-1.30", approx -> -1.297917436600305, interp -> Root[150524, +, 210474, #, +, 117189, #, ^, 2, +, 32427, #, ^, 3, +, 4455, #, ^, 4, +, 71123, #, ^, 5, +, 78489, #, ^, 6, +, 32472, #, ^, 7, +, 5940, #, ^, 8, +, 405, #, ^, 9, +, 13141, #, ^, 10, +, 10839, #, ^, 11, +, 2970, #, ^, 12, +, 270, #, ^, 13, +, 1206, #, ^, 15, +, 660, #, ^, 16, +, 90, #, ^, 17, +, 55, #, ^, 20, +, 15, #, ^, 21, +, #, ^, 25, &, 1, 0], head -> Root, big -> 150524+210474 #1+117189 #1^2+32427 #1^3+4455 #1^4+71123 #1^5+78489 #1^6+32472 #1^7+5940 #1^8+405 #1^9+13141 #1^10+10839 #1^11+2970 #1^12+270 #1^13+1206 #1^15+660 #1^16+90 #1^17+55 #1^20+15 #1^21+#1^25&, degree -> 25, shortDegree -> 25, number -> 1|>^23-186867 <|icon -> Root, small -> "-1.30", approx -> -1.297917436600305, interp -> Root[150524, +, 210474, #, +, 117189, #, ^, 2, +, 32427, #, ^, 3, +, 4455, #, ^, 4, +, 71123, #, ^, 5, +, 78489, #, ^, 6, +, 32472, #, ^, 7, +, 5940, #, ^, 8, +, 405, #, ^, 9, +, 13141, #, ^, 10, +, 10839, #, ^, 11, +, 2970, #, ^, 12, +, 270, #, ^, 13, +, 1206, #, ^, 15, +, 660, #, ^, 16, +, 90, #, ^, 17, +, 55, #, ^, 20, +, 15, #, ^, 21, +, #, ^, 25, &, 1, 0], head -> Root, big -> 150524+210474 #1+117189 #1^2+32427 #1^3+4455 #1^4+71123 #1^5+78489 #1^6+32472 #1^7+5940 #1^8+405 #1^9+13141 #1^10+10839 #1^11+2970 #1^12+270 #1^13+1206 #1^15+660 #1^16+90 #1^17+55 #1^20+15 #1^21+#1^25&, degree -> 25, shortDegree -> 25, number -> 1|>^24)
```

By default, [Reduce](https://reference.wolfram.com/language/ref/Reduce.html) uses [AlgebraicNumber](https://reference.wolfram.com/language/ref/AlgebraicNumber.html) objects to represent the solution here:

```wolfram
Reduce[y^5-3 y^2+x^5+4==0&&y^2+y-x^5-3 x-11==0,{x,y},Reals,Backsubstitution->True]
(* Output *)
x==<|icon -> AlgebraicNumber, small -> "-1.30", approx -> -1.297917436600305, interp -> AlgebraicNumber[Root[150524, +, 210474, #, +, 117189, #, ^, 2, +, 32427, #, ^, 3, +, 4455, #, ^, 4, +, 71123, #, ^, 5, +, 78489, #, ^, 6, +, 32472, #, ^, 7, +, 5940, #, ^, 8, +, 405, #, ^, 9, +, 13141, #, ^, 10, +, 10839, #, ^, 11, +, 2970, #, ^, 12, +, 270, #, ^, 13, +, 1206, #, ^, 15, +, 660, #, ^, 16, +, 90, #, ^, 17, +, 55, #, ^, 20, +, 15, #, ^, 21, +, #, ^, 25, &, 1, 0], 0100000000000000000000000], head -> AlgebraicNumber, big -> {0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, generator -> <|icon -> Root, small -> "-1.30", approx -> -1.297917436600305, interp -> Root[150524, +, 210474, #, +, 117189, #, ^, 2, +, 32427, #, ^, 3, +, 4455, #, ^, 4, +, 71123, #, ^, 5, +, 78489, #, ^, 6, +, 32472, #, ^, 7, +, 5940, #, ^, 8, +, 405, #, ^, 9, +, 13141, #, ^, 10, +, 10839, #, ^, 11, +, 2970, #, ^, 12, +, 270, #, ^, 13, +, 1206, #, ^, 15, +, 660, #, ^, 16, +, 90, #, ^, 17, +, 55, #, ^, 20, +, 15, #, ^, 21, +, #, ^, 25, &, 1, 0], head -> Root, big -> 150524+210474 #1+117189 #1^2+32427 #1^3+4455 #1^4+71123 #1^5+78489 #1^6+32472 #1^7+5940 #1^8+405 #1^9+13141 #1^10+10839 #1^11+2970 #1^12+270 #1^13+1206 #1^15+660 #1^16+90 #1^17+55 #1^20+15 #1^21+#1^25&, degree -> 25, shortDegree -> 25, number -> 1|>, degree -> 25|>&&y==<|icon -> AlgebraicNumber, small -> "1.42", approx -> 1.4164974834097863, interp -> AlgebraicNumber[...], head -> AlgebraicNumber, big -> {(4837482790)/(12953669),-(9375045027)/(25907338),-(9342236349)/(25907338),(1291893057)/(25907338),-(398647575)/(25907338),(1033763855)/(25907338),-(3083226180)/(12953669),-(1184444091)/(12953669),(1136838537)/(25907338),-(519316353)/(25907338),-(349540011)/(25907338),-(547580493)/(12953669),-(32775705)/(12953669),(87616404)/(12953669),-(63285948)/(12953669),-(26300648)/(12953669),-(36232251)/(12953669),(7880796)/(12953669),(3032586)/(12953669),-(9615105)/(25907338),-(1649643)/(25907338),-(808857)/(12953669),(421281)/(12953669),-(55971)/(12953669),-(186867)/(25907338)}, generator -> <|icon -> Root, small -> "-1.30", approx -> -1.297917436600305, interp -> Root[150524, +, 210474, #, +, 117189, #, ^, 2, +, 32427, #, ^, 3, +, 4455, #, ^, 4, +, 71123, #, ^, 5, +, 78489, #, ^, 6, +, 32472, #, ^, 7, +, 5940, #, ^, 8, +, 405, #, ^, 9, +, 13141, #, ^, 10, +, 10839, #, ^, 11, +, 2970, #, ^, 12, +, 270, #, ^, 13, +, 1206, #, ^, 15, +, 660, #, ^, 16, +, 90, #, ^, 17, +, 55, #, ^, 20, +, 15, #, ^, 21, +, #, ^, 25, &, 1, 0], head -> Root, big -> 150524+210474 #1+117189 #1^2+32427 #1^3+4455 #1^4+71123 #1^5+78489 #1^6+32472 #1^7+5940 #1^8+405 #1^9+13141 #1^10+10839 #1^11+2970 #1^12+270 #1^13+1206 #1^15+660 #1^16+90 #1^17+55 #1^20+15 #1^21+#1^25&, degree -> 25, shortDegree -> 25, number -> 1|>, degree -> 25|>
```

#### Modulus

Solve equations over the integers modulo 9:

```wolfram
Reduce[x^2+3y^2==4&&3x^3-4y^2+x y==1,{x,y},Modulus->9]
(* Output *)
(x==8&&y==1)||(x==8&&y==4)||(x==8&&y==7)
```

#### Quartics

By default, [Reduce](https://reference.wolfram.com/language/ref/Reduce.html) does not use general formulas for solving quartics in radicals:

```wolfram
Reduce[x^4+2 x^2+3 x+4==0,x]
(* Output *)
x==Root||x==Root||x==Root||x==Root
```

With [Quartics](https://reference.wolfram.com/language/ref/Quartics.html)->[True](https://reference.wolfram.com/language/ref/True.html), [Reduce](https://reference.wolfram.com/language/ref/Reduce.html) solves all quartics in terms of radicals:

```wolfram
Reduce[x^4+2 x^2+3 x+4==0,x,Quartics->True]
(* Output *)
x==(1)/(2) Sqrt[-(4)/(3)+(52)/(3 ((1)/(2) (-317+9 ⅈ Sqrt[5703]))^(1/3))+(1)/(3) ((1)/(2) (-317+9 ⅈ Sqrt[5703]))^(1/3)]-(1)/(2) Sqrt(-(8)/(3)-(52)/(3 ((1)/(2) (-317+9 ⅈ Sqrt[5703]))^(1/3))-(1)/(3) ((1)/(2) (-317+9 ⅈ Sqrt[5703]))^(1/3)-(6)/(Sqrt[-(4)/(3)+(52)/(3 ((1)/(2) (-317+9 ⅈ Sqrt[5703]))^(1/3))+(1)/(3) ((1)/(2) (-317+9 ⅈ Sqrt[5703]))^(1/3)]))||x==(1)/(2) Sqrt[-(4)/(3)+(52)/(3 ((1)/(2) (-317+9 ⅈ Sqrt[5703]))^(1/3))+(1)/(3) ((1)/(2) (-317+9 ⅈ Sqrt[5703]))^(1/3)]+(1)/(2) Sqrt(-(8)/(3)-(52)/(3 ((1)/(2) (-317+9 ⅈ Sqrt[5703]))^(1/3))-(1)/(3) ((1)/(2) (-317+9 ⅈ Sqrt[5703]))^(1/3)-(6)/(Sqrt[-(4)/(3)+(52)/(3 ((1)/(2) (-317+9 ⅈ Sqrt[5703]))^(1/3))+(1)/(3) ((1)/(2) (-317+9 ⅈ Sqrt[5703]))^(1/3)]))||x==-(1)/(2) Sqrt[-(4)/(3)+(52)/(3 ((1)/(2) (-317+9 ⅈ Sqrt[5703]))^(1/3))+(1)/(3) ((1)/(2) (-317+9 ⅈ Sqrt[5703]))^(1/3)]-(1)/(2) Sqrt(-(8)/(3)-(52)/(3 ((1)/(2) (-317+9 ⅈ Sqrt[5703]))^(1/3))-(1)/(3) ((1)/(2) (-317+9 ⅈ Sqrt[5703]))^(1/3)+(6)/(Sqrt[-(4)/(3)+(52)/(3 ((1)/(2) (-317+9 ⅈ Sqrt[5703]))^(1/3))+(1)/(3) ((1)/(2) (-317+9 ⅈ Sqrt[5703]))^(1/3)]))||x==-(1)/(2) Sqrt[-(4)/(3)+(52)/(3 ((1)/(2) (-317+9 ⅈ Sqrt[5703]))^(1/3))+(1)/(3) ((1)/(2) (-317+9 ⅈ Sqrt[5703]))^(1/3)]+(1)/(2) Sqrt(-(8)/(3)-(52)/(3 ((1)/(2) (-317+9 ⅈ Sqrt[5703]))^(1/3))-(1)/(3) ((1)/(2) (-317+9 ⅈ Sqrt[5703]))^(1/3)+(6)/(Sqrt[-(4)/(3)+(52)/(3 ((1)/(2) (-317+9 ⅈ Sqrt[5703]))^(1/3))+(1)/(3) ((1)/(2) (-317+9 ⅈ Sqrt[5703]))^(1/3)]))
```

#### WorkingPrecision

Finding the solution with exact computations takes a long time:

```wolfram
TimeConstrained[Reduce[x^77+3 x-11==ℯ&&y^5-x^2 y+21==π&&x^2 y^3+z^4==ℯ^π,{x,y,z},Reals],60]
(* Output *)
$Aborted
```

With [WorkingPrecision](https://reference.wolfram.com/language/ref/WorkingPrecision.html)->100, [Reduce](https://reference.wolfram.com/language/ref/Reduce.html) finds a solution fast, but it may be incorrect:

```wolfram
Timing[Reduce[x^77+3 x-11==ℯ&&y^5-x^2 y+21==π&&x^2 y^3+z^4==ℯ^π,{x,y,z},Reals,WorkingPrecision->100]]
(* Output *)
{0.015625,x==Root&&y==Root[21-π-x^2 #1+#1^5&,1]&&(z==-(ℯ^π-x^2 y^3)^(1/4)||z==(ℯ^π-x^2 y^3)^(1/4))}
```

### Applications

#### Basic Applications

Prove geometric inequalities for $a$, $b$, and $c$ sides of a triangle:

```wolfram
triangle=a>0&&b>0&&c>0&&a+b>c&&a+c>b&&b+c>a;
acute=a^2+b^2>c^2&&a^2+c^2>b^2&&b^2+c^2>a^2;
s=1/2(a+b+c);
F=Sqrt[s(s-a)(s-b)(s-c)];
```

Prove an inequality for triangles:

```wolfram
Reduce[ForAll[{a,b,c},triangle,a b c(a^2/b^2+b^2/c^2+c^2/a^2)>=a^3+b^3+c^3+a b(b-a)+a c(a-c)+b c(c-b)]]
(* Output *)
True
```

Prove an inequality for acute triangles:

```wolfram
Reduce[ForAll[{a,b,c},triangle&&acute,27(b^2+c^2-a^2)^2(a^2+c^2-b^2)^2(a^2+b^2-c^2)^2<=(4F)^6]]
(* Output *)
True
```

#### Polynomial Root Problems

Find conditions for a quartic to have all roots equal:

```wolfram
f[x_]:= x^4 + a x^3 + b x^2 + c x + d
```

Using quantifier elimination:

```wolfram
Reduce[Exists[x,f[x]==0,ForAll[y,f[y]==0, x==y]],{a,b,c,d}]
(* Output *)
b==(3 a^2)/(8)&&c==(a^3)/(16)&&d==(a^4)/(256)
```

Using [Subresultants](https://reference.wolfram.com/language/ref/Subresultants.html):

```wolfram
Reduce[Thread[Drop[Subresultants[f[x],D[f[x],x],x],-1]==0],{a,b,c,d},Backsubstitution->True]
(* Output *)
b==(3 a^2)/(8)&&c==(a^3)/(16)&&d==(a^4)/(256)
```

#### Parametrization Problems

Plot a space curve given by an implicit description:

```wolfram
curve=x^2+y^2+z^2==1&&x^3+x y^2==z^2;
```

```wolfram
red=Reduce[curve,{x,y,z},Reals]
(* Output *)
0<=x<=Root&&((y==-Sqrt[(1-x^2-x^3)/(1+x)]&&(z==-Sqrt[x^3+x y^2]||z==Sqrt[x^3+x y^2]))||(y==Sqrt[(1-x^2-x^3)/(1+x)]&&(z==-Sqrt[x^3+x y^2]||z==Sqrt[x^3+x y^2])))
```

```wolfram
bound=N[First[red]]
(* Output *)
0.<=x<=0.7548776662466927
```

```wolfram
pieces={x,y,z}//.{ToRules[Rest[red]]}
(* Output *)
{{x,-Sqrt[(1-x^2-x^3)/(1+x)],-Sqrt[x^3+(x (1-x^2-x^3))/(1+x)]},{x,-Sqrt[(1-x^2-x^3)/(1+x)],Sqrt[x^3+(x (1-x^2-x^3))/(1+x)]},{x,Sqrt[(1-x^2-x^3)/(1+x)],-Sqrt[x^3+(x (1-x^2-x^3))/(1+x)]},{x,Sqrt[(1-x^2-x^3)/(1+x)],Sqrt[x^3+(x (1-x^2-x^3))/(1+x)]}}
```

```wolfram
ParametricPlot3D[pieces,{x,First[bound],Last[bound]}]
```

*([Graphics3D])*

Plot the projection of the space curve on the $x$-$y$ plane:

```wolfram
proj=Reduce[Exists[z,curve],{x,y},Reals]
(* Output *)
0<=x<=Root&&(y==-Sqrt[(1-x^2-x^3)/(1+x)]||y==Sqrt[(1-x^2-x^3)/(1+x)])
```

```wolfram
bounds=N[First[proj]]
(* Output *)
0.<=x<=0.7548776662466927
```

```wolfram
pieces2d={x,y}//.{ToRules[Rest[proj]]}
(* Output *)
{{x,-Sqrt[(1-x^2-x^3)/(1+x)]},{x,Sqrt[(1-x^2-x^3)/(1+x)]}}
```

```wolfram
ParametricPlot[pieces2d,{x,First[bounds],Last[bounds]}]
```

*([Graphics])*

#### Integer Problems

Find a Pythagorean triple:

```wolfram
Reduce[x^2+y^2==5^2&&y>x>0,{x,y},Integers]
(* Output *)
x==3&&y==4
```

Find a sequence of Pythagorean triples:

```wolfram
Table[Reduce[x^2+y^2==z^2&&y>x>0,{x,y},Integers],{z,30}]
(* Output *)
{False,False,False,False,x==3&&y==4,False,False,False,False,x==6&&y==8,False,False,x==5&&y==12,False,x==9&&y==12,False,x==8&&y==15,False,False,x==12&&y==16,False,False,False,False,(x==7&&y==24)||(x==15&&y==20),x==10&&y==24,False,False,x==20&&y==21,x==18&&y==24}
```

Find how to pay $2.27 postage with 10-, 23- and 37-cent stamps:

```wolfram
Reduce[a 10 + b 23 + c 37 == 227 &&a>=0&&b>=0&&c>=0,{a,b,c},Integers]
(* Output *)
(a==1&&b==3&&c==4)||(a==2&&b==9&&c==0)||(a==7&&b==2&&c==3)||(a==13&&b==1&&c==2)||(a==19&&b==0&&c==1)
```

The same task can be accomplished with [IntegerPartitions](https://reference.wolfram.com/language/ref/IntegerPartitions.html):

```wolfram
IntegerPartitions[227,All,{37,23,10}]
(* Output *)
{{10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,37},{10,10,10,10,10,10,10,10,10,10,10,10,10,23,37,37},{10,10,10,10,10,10,10,23,23,37,37,37},{10,10,23,23,23,23,23,23,23,23,23},{10,23,23,23,37,37,37,37}}
```

Show that there are only five regular polyhedrons:

```wolfram
Show[PolyhedronData[#],ImageSize->Tiny]&/@PolyhedronData["Platonic"]
(* Output *)
{[Graphics3D],[Graphics3D],[Graphics3D],[Graphics3D],[Graphics3D]}
```

Each face for a regular $n$-gon contributes $n$ edges, but they are shared, so they are counted twice:

```wolfram
e=(1)/(2)(n f);
```

Each face for a regular $n$-gon contributes $n$ vertices, but they are shared, so they are counted $m$ times:

```wolfram
v=(1)/(m)(n f);
```

Using Euler's formula $f-e+v=2$, find the number of faces:

```wolfram
Reduce[f-e+v==2&&m>=3&&n>=3,f,Reals]
(* Output *)
((3<=n<=6&&(3<=m<(2 n)/(-2+n)||m>(2 n)/(-2+n)))||(n>6&&m>=3))&&f==-(4 m)/(-2 m-2 n+m n)
```

For this last formula to be well defined, the denominator needs to be positive and an integer:

```wolfram
Reduce[2 m+2 n-m n>0&&m>=3&&n>=3,{m,n},Integers]
(* Output *)
(m==3&&n==3)||(m==3&&n==4)||(m==3&&n==5)||(m==4&&n==3)||(m==5&&n==3)
```

Hence the following five cases:

```wolfram
{f,e,v}/.f->-(4 m)/(-2 m-2 n+m n)/.FindInstance[2 m+2 n-m n>0&&m>=3&&n>=3,{m,n},Integers,5]//Sort
(* Output *)
{{4,6,4},{6,12,8},{8,12,6},{12,30,20},{20,30,12}}
```

Compare this to the actual counts in [PolyhedronData](https://reference.wolfram.com/language/ref/PolyhedronData.html):

```wolfram
Table[{PolyhedronData[p,"FaceCount"],PolyhedronData[p,"EdgeCount"],PolyhedronData[p,"VertexCount"]},{p,PolyhedronData["Platonic"]}]//Sort
(* Output *)
{{4,6,4},{6,12,8},{8,12,6},{12,30,20},{20,30,12}}
```

#### Geometry Problems

The region `ℛ` is a subset of `𝒮` if $\forall_{\{x,y,\ldots \}}\{x,y,\ldots \}\in \mathcal{R}\Rightarrow \{x,y,\ldots \}\in \mathcal{S}$ is true. Show that [Disk](https://reference.wolfram.com/language/ref/Disk.html)[{0,0},{2,1}] is a subset of [Rectangle](https://reference.wolfram.com/language/ref/Rectangle.html)[{-2,-1},{2,1}]:

```wolfram
ℛ=Disk[{0,0},{2,1}];
𝒮=Rectangle[{-2,-1},{2,1}];
```

```wolfram
Reduce[∀_{x,y}{x,y}∈ℛ⟹{x,y}∈𝒮,Reals]
(* Output *)
True
```

Plot it:

```wolfram
Graphics[{{LightRed,EdgeForm[Gray],𝒮},{LightBlue,EdgeForm[Gray],ℛ}}]
```

*([Graphics])*

Show that [Cylinder](https://reference.wolfram.com/language/ref/Cylinder.html)[]⊆[Ball](https://reference.wolfram.com/language/ref/Ball.html)[{0,0,0},2]:

```wolfram
ℛ=Cylinder[];
𝒮=Ball[{0,0,0},2];
```

```wolfram
Reduce[∀_{x,y,z}{x,y,z}∈ℛ⟹{x,y,z}∈𝒮,Reals]
(* Output *)
True
```

Plot it:

```wolfram
Graphics3D[{{Opacity[0.3],𝒮},{LightBlue,EdgeForm[Gray],ℛ}}]
```

*([Graphics3D])*

For a finite point set $\{p_{1},\ldots,p_{n}\}$, the Voronoi cell for a point $p_{i}$ can be defined by $\{p|\wedge_{1 \leq j \leq n \wedge j \neq i}d(p_{i},p)\leq d(p_{j},p)\}$, which corresponds to all points closer to $p_{i}$ than any other point $p_{j}$ for $j \neq i$. Find a simple formula for a Voronoi cell, using [Reduce](https://reference.wolfram.com/language/ref/Reduce.html):

```wolfram
pts={{0,0},{1,0},{1,1},{2,0}};
```

The Voronoi cell associated with `*pts*[[1]]` is given by:

```wolfram
voronoi1=And@@Table[EuclideanDistance[pts[[1]],{x,y}]<=EuclideanDistance[pts[[j]],{x,y}],{j,2,4}]
(* Output *)
Sqrt[Abs[x]^2+Abs[y]^2]<=Sqrt[Abs[1-x]^2+Abs[y]^2]&&Sqrt[Abs[x]^2+Abs[y]^2]<=Sqrt[Abs[1-x]^2+Abs[1-y]^2]&&Sqrt[Abs[x]^2+Abs[y]^2]<=Sqrt[Abs[2-x]^2+Abs[y]^2]
```

The resulting cell is given by an intersection of half-spaces:

```wolfram
Reduce[voronoi1,{x,y},Reals]
(* Output *)
x<=(1)/(2)&&y<=1-x
```

```wolfram
RegionPlot[%,{x,-1,1},{y,-1,1}]
```

*([Graphics])*

Find simple formulas for all Voronoi cells:

```wolfram
vcells=And@@@Table[EuclideanDistance[pts[[i]],{x,y}]<=EuclideanDistance[pts[[j]],{x,y}],{i,4},{j,Complement[Range[4],{i}]}];
```

```wolfram
vscells=Reduce[#,{x,y},Reals]&/@vcells
(* Output *)
{x<=(1)/(2)&&y<=1-x,(1)/(2)<=x<=(3)/(2)&&y<=(1)/(2),(x<=(1)/(2)&&y>=1-x)||((1)/(2)<x<=(3)/(2)&&y>=(1)/(2))||(x>(3)/(2)&&y>=-1+x),x>=(3)/(2)&&y<=-1+x}
```

Plot them:

```wolfram
RegionPlot[Evaluate@Table[vscells[[i]],{i,4}],{x,-1,3},{y,-1,3},PlotLegends->SwatchLegend[Map[Style[#,Small]&,vscells]],Epilog->Point[pts]]
```

*([Graphics])*

### Properties & Relations

The result of reduction is equivalent to the original system:

```wolfram
syst=x^2+y^2+z^2<=1&&2x y>z^2;
```

```wolfram
red=Reduce[syst,{x,y,z},Reals]
(* Output *)
(-1<x<0&&((y==-Sqrt[1-x^2]&&z==0)||(-Sqrt[1-x^2]<y<-1-x&&-Sqrt[1-x^2-y^2]<=z<=Sqrt[1-x^2-y^2])||(y==-1-x&&-Sqrt[1-x^2-y^2]<z<Sqrt[1-x^2-y^2])||(-1-x<y<0&&-Sqrt[2] Sqrt[x y]<z<Sqrt[2] Sqrt[x y])))||(0<x<1&&((0<y<1-x&&-Sqrt[2] Sqrt[x y]<z<Sqrt[2] Sqrt[x y])||(y==1-x&&-Sqrt[1-x^2-y^2]<z<Sqrt[1-x^2-y^2])||(1-x<y<Sqrt[1-x^2]&&-Sqrt[1-x^2-y^2]<=z<=Sqrt[1-x^2-y^2])||(y==Sqrt[1-x^2]&&z==0)))
```

```wolfram
SeedRandom[0];Table[{syst,red}/.{x->RandomReal[{-1,1}],y->RandomReal[{-1,1}],z->RandomReal[{-1,1}]},{10}]
(* Output *)
{{True,True},{False,False},{False,False},{False,False},{True,True},{False,False},{False,False},{False,False},{True,True},{False,False}}
```

[ToRules](https://reference.wolfram.com/language/ref/ToRules.html) and [ReplaceRepeated](https://reference.wolfram.com/language/ref/ReplaceRepeated.html) can be used to backsubstitute finite solution sets:

```wolfram
Reduce[x^2+y==1&&x^2-y^2==2,{x,y}]
(* Output *)
(x==-Sqrt[(3)/(2)-(ⅈ Sqrt[3])/(2)]||x==Sqrt[(3)/(2)-(ⅈ Sqrt[3])/(2)]||x==-Sqrt[(3)/(2)+(ⅈ Sqrt[3])/(2)]||x==Sqrt[(3)/(2)+(ⅈ Sqrt[3])/(2)])&&y==1-x^2
```

```wolfram
x^2+y==1&&x^2-y^2==2//.{ToRules[%]}
(* Output *)
{(3)/(2)-(ⅈ Sqrt[3])/(2)-(-(1)/(2)+(ⅈ Sqrt[3])/(2))^2==2,(3)/(2)-(ⅈ Sqrt[3])/(2)-(-(1)/(2)+(ⅈ Sqrt[3])/(2))^2==2,(3)/(2)+(ⅈ Sqrt[3])/(2)-(-(1)/(2)-(ⅈ Sqrt[3])/(2))^2==2,(3)/(2)+(ⅈ Sqrt[3])/(2)-(-(1)/(2)-(ⅈ Sqrt[3])/(2))^2==2}
```

Use [Expand](https://reference.wolfram.com/language/ref/Expand.html) to simplify a result of substitution involving simple radicals:

```wolfram
Expand[%]
(* Output *)
{True,True,True,True}
```

To simplify expressions involving algebraic numbers, use [RootReduce](https://reference.wolfram.com/language/ref/RootReduce.html):

```wolfram
Reduce[x^2+y^2==1&&x^3-2y^3==3,{x,y}]
(* Output *)
(x==Root||x==Root||x==Root||x==Root||x==Root||x==Root)&&y==(1)/(16) (-9+5 x+27 x^2+7 x^3-15 x^4-5 x^5)
```

```wolfram
RootReduce[x^2+y^2==1&&x^3-2y^3==3//.{ToRules[%]}]
(* Output *)
{True,True,True,True,True,True}
```

To find solution instances, use [FindInstance](https://reference.wolfram.com/language/ref/FindInstance.html):

```wolfram
FindInstance[x^2+y^2+z^2<=1&&2x y>z^2,{x,y,z},Reals]
(* Output *)
{{x->-(1)/(2),y->-(1)/(2),z->0}}
```

```wolfram
FindInstance[x^2-3y^2==1&&0<x<10^10,{x,y},Integers,3]
(* Output *)
{{x->708158977,y->408855776},{x->2642885282,y->-1525870529},{x->978122,y->564719}}
```

[Solve](https://reference.wolfram.com/language/ref/Solve.html) represents solutions of complex equations in terms of replacement rules:

```wolfram
Solve[x^2+y^2==1&&2x+3y==1,{x,y}]
(* Output *)
{{x->(2)/(13) (1-3 Sqrt[3]),y->(1)/(13) (3+4 Sqrt[3])},{x->(2)/(13) (1+3 Sqrt[3]),y->(1)/(13) (3-4 Sqrt[3])}}
```

[Solve](https://reference.wolfram.com/language/ref/Solve.html) omits solutions involving equations on parameters:

```wolfram
Reduce[a x^2+x==1,x]
(* Output *)
(a==0&&x==1)||(a≠0&&(x==(-1-Sqrt[1+4 a])/(2 a)||x==(-1+Sqrt[1+4 a])/(2 a)))
```

```wolfram
Solve[a x^2+x==1,x]
(* Output *)
{{x->(-1-Sqrt[1+4 a])/(2 a)},{x->(-1+Sqrt[1+4 a])/(2 a)}}
```

For transcendental equations, [Solve](https://reference.wolfram.com/language/ref/Solve.html) may not give all solutions:

```wolfram
Solve[x+E^x==1/2,x]
(* Output *)
Solve
(* Output *)
{{x->(1)/(2) (1-2 ProductLog[Sqrt[ℯ]])}}
```

```wolfram
Reduce[x+E^x==1/2,x]
(* Output *)
1∈Integers&&x==(1)/(2)-ProductLog[1,Sqrt[ℯ]]
```

Using inverse functions allows [Solve](https://reference.wolfram.com/language/ref/Solve.html) to find some solutions fast:

```wolfram
Solve[x^n==1,x]//Timing
(* Output *)
Solve
(* Output *)
{0.,{{x->1}}}
```

Finding the complete solution may take much longer, and the solution may be large:

```wolfram
(red=Reduce[x^n==1,x])//LeafCount//Timing
(* Output *)
{4.484375,611}
```

This finds the values of $n$ for which `*x*==2 is a solution:

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

[SolveAlways](https://reference.wolfram.com/language/ref/SolveAlways.html) gives the values of parameters for which complex equations are always true:

```wolfram
SolveAlways[(a -2b+1)x^2+(a-b^2-c)x ==a^2-b+3c-1,x]
(* Output *)
{{a->-2-Sqrt[13],b->(1)/(2) (-1-Sqrt[13]),c->(1)/(2) (-11-3 Sqrt[13])},{a->-2+Sqrt[13],b->(1)/(2) (-1+Sqrt[13]),c->-(11)/(2)+(3 Sqrt[13])/(2)}}
```

This solves the same problem using [Reduce](https://reference.wolfram.com/language/ref/Reduce.html):

```wolfram
Reduce[ForAll[x,(a -2b+1)x^2+(a-b^2-c)x ==a^2-b+3c-1],{a,b,c}]
(* Output *)
(a==-2-Sqrt[13]||a==-2+Sqrt[13])&&b==(1+a)/(2)&&c==(1)/(2) (-5+3 a)
```

[Resolve](https://reference.wolfram.com/language/ref/Resolve.html) eliminates quantifiers, possibly without solving the resulting quantifier[Hyphen]free system:

```wolfram
Resolve[Exists[x,x^2+y^2+z^2==1&&x y>z^3],Reals]
(* Output *)
(z^3<0&&y^2+z^2<=1)||(y^2+z^2<=1&&-y^2+y^4+y^2 z^2+z^6<0)
```

```wolfram
Resolve[ForAll[{x,y},a x^2+b x+c==0&&a y^2+b y+c==0,x==y]]
(* Output *)
(a==0&&b≠0)||(a==0&&c≠0)||(a≠0&&b^2-4 a c==0)
```

[Eliminate](https://reference.wolfram.com/language/ref/Eliminate.html) eliminates variables from systems of complex equations:

```wolfram
Eliminate[x^2+y^2+z^2==1&&x y==z^3,x]
(* Output *)
y^4+y^2 (-1+z^2)==-z^6
```

This solves the same problem using [Resolve](https://reference.wolfram.com/language/ref/Resolve.html):

```wolfram
Resolve[Exists[x,x^2+y^2+z^2==1&&x y==z^3]]
(* Output *)
(y==0&&z==0)||(y≠0&&-y^2+y^4+y^2 z^2+z^6==0)
```

[Reduce](https://reference.wolfram.com/language/ref/Reduce.html) additionally solves the resulting equations:

```wolfram
Reduce[Exists[x,x^2+y^2+z^2==1&&x y==z^3],{y,z}]
(* Output *)
(y==0&&z==0)||((z==-Sqrt[Root[-y^2+y^4+y^2 #1+#1^3&,1]]||z==Sqrt[Root[-y^2+y^4+y^2 #1+#1^3&,1]]||z==-Sqrt[Root[-y^2+y^4+y^2 #1+#1^3&,2]]||z==Sqrt[Root[-y^2+y^4+y^2 #1+#1^3&,2]]||z==-Sqrt[Root[-y^2+y^4+y^2 #1+#1^3&,3]]||z==Sqrt[Root[-y^2+y^4+y^2 #1+#1^3&,3]])&&y≠0)
```

### Possible Issues

Because $x$ appears in an inequality, it is assumed to be real; $y$ is allowed to be complex:

```wolfram
Reduce[x^2<1 && y-Sqrt[x]==0, {x, y}]
(* Output *)
-1<x<1&&y==Sqrt[x]
```

When domain [Reals](https://reference.wolfram.com/language/ref/Reals.html) is specified, $x$, $y$, and [Sqrt](https://reference.wolfram.com/language/ref/Sqrt.html)[*x*] are required to be real:

```wolfram
Reduce[x^2<1 && y-Sqrt[x]==0, {x, y},Reals]
(* Output *)
0<=x<1&&y==Sqrt[x]
```

This allows complex values of $x$ for which both sides of the inequality are real:

```wolfram
Reduce[x^2<1 && y-Sqrt[x]==0, {x, y},Complexes]
(* Output *)
((-1<Re[x]<0&&Im[x]==0)||Re[x]==0||(0<Re[x]<1&&Im[x]==0))&&y==Sqrt[x]
```

The order of variables may affect solvability of the problem. This is easy to solve for $y$ in terms of $x$:

```wolfram
Reduce[y==x-Sin[x],{x,y}]
(* Output *)
y==x-Sin[x]
```

Solving for $x$ in terms of $y$ would be much harder:

```wolfram
Reduce[y==x-Sin[x],{y,x}]
(* Output *)
Reduce
(* Output *)
Reduce[y==x-Sin[x],{y,x}]
```

[Reduce](https://reference.wolfram.com/language/ref/Reduce.html) does not solve equations that depend on branch cuts of Wolfram Language functions:

```wolfram
Reduce[Sqrt[x+2y]-x^2+y-1==0,{x,y}]
(* Output *)
Reduce
(* Output *)
(0==-1+Sqrt[3+x+2 x^2]-Sqrt[4+x+2 x^2-2 Sqrt[3+x+2 x^2]]&&y==2+x^2-Sqrt[3+x+2 x^2])||(0==-1-Sqrt[3+x+2 x^2]-Sqrt[4+x+2 x^2+2 Sqrt[3+x+2 x^2]]&&y==2+x^2+Sqrt[3+x+2 x^2])
```

Plot the region where the first condition is nonzero:

```wolfram
Block[{x=u+I v},RegionPlot[Abs[-1+Sqrt[3+x+2 x^2]-Sqrt[4+x+2 x^2-2 Sqrt[3+x+2 x^2]]]>0.001,{u,-10,10},{v,-10,10}]]
```

*([Graphics])*

[Reduce](https://reference.wolfram.com/language/ref/Reduce.html) will not find real solutions that fall outside the real domains of subexpressions:

```wolfram
Reduce[z>0&&Im[HankelH1[0,z]]==0,z,Reals]
(* Output *)
False
```

[Reduce](https://reference.wolfram.com/language/ref/Reduce.html) uses [FunctionDomain](https://reference.wolfram.com/language/ref/FunctionDomain.html) to determine real domains of mathematical functions:

```wolfram
FunctionDomain[HankelH1[0,z],z]
(* Output *)
False
```

[FunctionDomain](https://reference.wolfram.com/language/ref/FunctionDomain.html) has real domain information that is accurate up to lower-dimensional sets:

```wolfram
Plot[Im[HankelH1[0,z]],{z,-3,3}]
```

*([Graphics])*

There exist isolated points at which [HankelH1](https://reference.wolfram.com/language/ref/HankelH1.html)[0,*z*] is real valued:

```wolfram
N[HankelH1[0,BesselYZero[0,1]],20]
(* Output *)
0.81012385935356425413537361505687207685+0`39.106962020856386 ⅈ
```

Removable singularities of input equations are generally not considered valid solutions:

```wolfram
Reduce[(x^2-2Sqrt[2]x+2)/(x-Sqrt[2])==0,x]
(* Output *)
False
```

```wolfram
Limit[(x^2-2Sqrt[2]x+2)/(x-Sqrt[2]),x->Sqrt[2]]
(* Output *)
0
```

However, solutions may include removable singularities that are cancelled by automatic simplification:

```wolfram
Reduce[x^2/x==0,x]
(* Output *)
x==0
```

The removable singularity at $x=0$ is cancelled by evaluation:

```wolfram
x^2/x==0
(* Output *)
x==0
```

Here the removable singularity at $x=1$ is cancelled by [Together](https://reference.wolfram.com/language/ref/Together.html), which is used to preprocess the equation:

```wolfram
Reduce[(x^2-2x+1)/(x-1)==0,x]
(* Output *)
x==1
```

```wolfram
Together[(x^2-2x+1)/(x-1)==0]
(* Output *)
-1+x==0
```

### Neat Examples

Find the vertical asymptotes of $\frac{1}{y^{2}-2}$ by directly using the definition of limit:

```wolfram
Reduce[ForAll[M, M > 0, Exists[δ, δ > 0,
    ForAll[y, Element[x | y, Reals] && 0 < Abs[y - x] < δ,
     Abs[1/(y^2 - 2)] > M]]], x]
(* Output *)
x==-Sqrt[2]||x==Sqrt[2]
```

## Tech Notes ▪Solving Equations ▪Inequalities ▪Generic and Non[Hyphen]Generic Solutions ▪Equations and Inequalities over Domains ▪Solving Logical Combinations of Equations ▪The Representation of Solution Sets ▪Complex Polynomial Systems ▪Real Polynomial Systems ▪Diophantine Polynomial Systems ▪Implementation notes: Algebra and Calculus

## Related Guides ▪Manipulating Equations ▪Formula Manipulation ▪Number Theory ▪Computational Geometry ▪Polynomial Algebra ▪Polynomial Equations ▪Diophantine Equations ▪Polynomial Systems ▪Equation Solving ▪Theorem Proving ▪Inequalities ▪Solvers over Regions ▪Solid Geometry ▪Plane Geometry ▪Finite Fields ▪Finite Mathematics ▪Polygons ▪Polyhedra ▪Symbolic Vectors, Matrices and Arrays

## History Introduced in 1988 (1.0) | Updated in 2003 (5.0) ▪ 2014 (10.0) ▪ 2024 (14.0)
