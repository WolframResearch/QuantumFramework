# CylindricalDecomposition | [SpanFromLeft]

> [CylindricalDecomposition](https://reference.wolfram.com/language/ref/CylindricalDecomposition.html)[*expr*,{*x*_1,*x*_2,…}] — finds a decomposition of the region represented by the statement `*expr*` into cylindrical parts whose directions correspond to the successive `*x*_*i*`.
> [CylindricalDecomposition](https://reference.wolfram.com/language/ref/CylindricalDecomposition.html)[*expr*,{*x*_1,*x*_2,…},*op*] — finds a decomposition of the result of applying the topological operation `*op*` to the region represented by the statement `*expr*`.
> [CylindricalDecomposition](https://reference.wolfram.com/language/ref/CylindricalDecomposition.html)[*expr*,{*x*_1,*x*_2,…},"Function"] — represents the result as [CylindricalDecompositionFunction](https://reference.wolfram.com/language/ref/CylindricalDecompositionFunction.html)[…][*x*_1,*x*_2,…] that can be efficiently used in further computation.

## Details and Options

The statement `*expr*` can be any logical combination of:

*lhs*==*rhs* | equations
*lhs*!=*rhs* | inequations
`*lhs*>*rhs*` or `*lhs*>=*rhs*` | inequalities
[CylindricalDecompositionFunction](https://reference.wolfram.com/language/ref/CylindricalDecompositionFunction.html)[…][*x*,*y*,…] | cylindrical algebraic formulas
[ForAll](https://reference.wolfram.com/language/ref/ForAll.html)[*x*,*cond*,*expr*] | universal quantifiers
[Exists](https://reference.wolfram.com/language/ref/Exists.html)[*x*,*cond*,*expr*] | existential quantifiers

Equations and inequalities may involve polynomials, rational functions or real algebraic functions.

[CylindricalDecomposition](https://reference.wolfram.com/language/ref/CylindricalDecomposition.html) assumes that all variables are real.

[CylindricalDecomposition](https://reference.wolfram.com/language/ref/CylindricalDecomposition.html) returns inequalities whose bounds in general involve algebraic functions.

The topological operation `*op*` can be one of:

`"Boundary"` | boundary of the solution set
`"Closure"` | closure of the solution set
`"Interior"` | interior of the solution set
`"Exterior"` | exterior of the solution set
`"ClosureOfInterior"` | closure of the interior of the solution set
`"InteriorOfClosure"` | interior of the closure of the solution set
`"Components"` | connected components of the solution set

[CylindricalDecompositionFunction](https://reference.wolfram.com/language/ref/CylindricalDecompositionFunction.html) objects provide an explicit compact representation of semialgebraic sets that can be efficiently used in further computation.

[CylindricalDecompositionFunction](https://reference.wolfram.com/language/ref/CylindricalDecompositionFunction.html) objects are typically used for repeated computations with semialgebraic sets, including computing Boolean combinations of sets, restricting sets with additional conditions, eliminating some variables or optimizing over sets.

A cylindrical algebraic formula in `*x*_1,…,*x*_*n*` has the form $F=C_{1}\vee \ldots \vee C_{k}$, where $C_{i}=B_{i,1}(x_{1})\wedge B_{i,2}(x_{1},x_{2})\ldots \wedge B_{i,n}(x_{1},\ldots,x_{n})$. Each $B_{i,j}$ has the form either $x_{j}=r_{i,j}(x_{1},\ldots,x_{j-1})$ or $r_{i,j}(x_{1},\ldots,x_{j-1})<x_{j}<s_{i,j}(x_{1},\ldots,x_{j-1})$, where $r_{i,j}$ and $s_{i,j}$ are algebraic functions that are defined and continuous on the solution set of $B_{i,1}\wedge \ldots \wedge B_{i,j-1}$. The solution set $Z(C_{i})$ of $C_{i}$ in $\mathbb{R}^{n}$ is called a cell. Projections of cells $Z(C_{i_{1}})$ and $Z(C_{i_{2}})$ on $\mathbb{R}^{m}$, for any $1 \leq m \leq n$, are either disjoint or identical.

Without the `"Function"` specification, [CylindricalDecomposition](https://reference.wolfram.com/language/ref/CylindricalDecomposition.html) returns cylindrical algebraic formulas written explicitly as Boolean combinations of equations and inequalities.

[CylindricalDecompositionFunction](https://reference.wolfram.com/language/ref/CylindricalDecompositionFunction.html) provides an encapsulated representation of cylindrical algebraic formulas, which is often more efficient when used in the input to [CylindricalDecomposition](https://reference.wolfram.com/language/ref/CylindricalDecomposition.html) or to solvers such as [Reduce](https://reference.wolfram.com/language/ref/Reduce.html), [Resolve](https://reference.wolfram.com/language/ref/Resolve.html), [FindInstance](https://reference.wolfram.com/language/ref/FindInstance.html), [Solve](https://reference.wolfram.com/language/ref/Solve.html) or [Minimize](https://reference.wolfram.com/language/ref/Minimize.html).

[Normal](https://reference.wolfram.com/language/ref/Normal.html) converts [CylindricalDecompositionFunction](https://reference.wolfram.com/language/ref/CylindricalDecompositionFunction.html) objects to explicit Boolean combinations of equations and inequalities.

Specifications of topological operations and output format can be combined, e.g. [CylindricalDecomposition](https://reference.wolfram.com/language/ref/CylindricalDecomposition.html)[*ineqs*,{*x*_1,*x*_2,…},"BoundaryFunction"] gives the boundary of the solution set of `*ineqs*` represented as a [CylindricalDecompositionFunction](https://reference.wolfram.com/language/ref/CylindricalDecompositionFunction.html) object.

## Examples

### Basic Examples

Find a cylindrical decomposition of the unit disc:

```wolfram
CylindricalDecomposition[x^2+y^2<1,{x,y}]
(* Output *)
-1<x<1&&-Sqrt[1-x^2]<y<Sqrt[1-x^2]
```

Find a cylindrical decomposition of the boundary of the unit disc:

```wolfram
CylindricalDecomposition[x^2+y^2<1,{x,y}, "Boundary"]
(* Output *)
(x==-1&&y==Sqrt[1-x^2])||(-1<x<1&&(y==-Sqrt[1-x^2]||y==Sqrt[1-x^2]))||(x==1&&y==Sqrt[1-x^2])
```

### Scope

#### Basic Uses

For univariate polynomials, the result consists of intervals:

```wolfram
CylindricalDecomposition[(x-1)(x-2)^2(x-3)(x-4)(x-5)>0,x]
(* Output *)
x<1||3<x<4||x>5
```

```wolfram
Show[Plot[(x-1)(x-2)^2(x-3)(x-4)(x-5),{x,0.5,5.5}],Graphics[{Red,Point[{{1,0},{2,0},{3,0},{4,0},{5,0}}]}]]
```

*([Graphics])*

```wolfram
NumberLinePlot[%%,{x,0.5,5.5}]
```

*([Graphics])*

In general, individual points can occur:

```wolfram
CylindricalDecomposition[(x-1)(x-2)^2(x-3)(x-4)(x-5)>=0,x]
(* Output *)
x<=1||x==2||3<=x<=4||x>=5
```

```wolfram
NumberLinePlot[%,{x,0.5,5.5}]
```

*([Graphics])*

This is the form for any logical combination as well:

```wolfram
CylindricalDecomposition[(x-1)(x-2)>=0&&(x-1/2)(x-3/2)>=0,x]
(* Output *)
x<=(1)/(2)||x>=2
```

```wolfram
NumberLinePlot[%,{x,0,3}]
```

*([Graphics])*

```wolfram
CylindricalDecomposition[Xor[(x-1)(x-2)>=0,(x-1/2)(x-3/2)>=0],x]
(* Output *)
(1)/(2)<x<=1||(3)/(2)<=x<2
```

```wolfram
NumberLinePlot[%,{x,0,3}]
```

*([Graphics])*

For multivariate polynomials, the result is in cylinder form $l_{1}<x<u_{1}\wedge l_{2}[x]<y<u_{2}[x]\wedge \cdots$:

```wolfram
CylindricalDecomposition[x^2+y^2<1,{x,y}]
(* Output *)
-1<x<1&&-Sqrt[1-x^2]<y<Sqrt[1-x^2]
```

In general, several cylinders will result:

```wolfram
CylindricalDecomposition[x^2+y^2<1&&y>x,{x,y}]
(* Output *)
(-1<x<=-(1)/(Sqrt[2])&&-Sqrt[1-x^2]<y<Sqrt[1-x^2])||(-(1)/(Sqrt[2])<x<(1)/(Sqrt[2])&&x<y<Sqrt[1-x^2])
```

Plot the individual cylinders using [RegionPlot](https://reference.wolfram.com/language/ref/RegionPlot.html):

```wolfram
{cyl1,cyl2}=List@@%;
```

```wolfram
RegionPlot[{cyl1,cyl2},{x,-1,1},{y,-1,1},PlotLegends->"Expressions"]
```

*([Graphics])*

By changing the order of variables, the cylinders take the form $l_{1}<y<u_{1}\wedge l_{2}[y]<x<u_{2}[y]$:

```wolfram
CylindricalDecomposition[x^2+y^2<1&&y>x,{y,x}]
(* Output *)
(-(1)/(Sqrt[2])<y<=(1)/(Sqrt[2])&&-Sqrt[1-y^2]<x<y)||((1)/(Sqrt[2])<y<1&&-Sqrt[1-y^2]<x<Sqrt[1-y^2])
```

Plot the individual cylinders:

```wolfram
{cyl1,cyl2}=List@@%;
```

```wolfram
RegionPlot[{cyl1,cyl2},{x,-1,1},{y,-1,1},PlotLegends->"Expressions"]
```

*([Graphics])*

Here cylinders of dimensions 0, 2, and 1 occur in the result:

```wolfram
CylindricalDecomposition[x^2+y^2<=1||x==1,{x,y}]
(* Output *)
(x==-1&&y==-Sqrt[1-x^2])||(-1<x<1&&-Sqrt[1-x^2]<=y<=Sqrt[1-x^2])||x==1
```

Three- and four-dimensional decompositions:

```wolfram
CylindricalDecomposition[x^2+y^2+z^2<1,{x,y,z}]
(* Output *)
-1<x<1&&-Sqrt[1-x^2]<y<Sqrt[1-x^2]&&-Sqrt[1-x^2-y^2]<z<Sqrt[1-x^2-y^2]
```

```wolfram
CylindricalDecomposition[x^2+y^2+z^2+w^2<1,{x,y,z,w}]
(* Output *)
-1<x<1&&-Sqrt[1-x^2]<y<Sqrt[1-x^2]&&-Sqrt[1-x^2-y^2]<z<Sqrt[1-x^2-y^2]&&-Sqrt[1-x^2-y^2-z^2]<w<Sqrt[1-x^2-y^2-z^2]
```

[CylindricalDecomposition](https://reference.wolfram.com/language/ref/CylindricalDecomposition.html) allows quantified formulas:

```wolfram
CylindricalDecomposition[Exists[y,x^2+y^2<1&&y>x],x]
(* Output *)
-1<x<(1)/(Sqrt[2])
```

Coefficients can include real algebraic numbers:

```wolfram
CylindricalDecomposition[Sqrt[2]x^2+Root[#^5-2#-11&,1]y^2<1,{x,y}]
(* Output *)
-(1)/(2^(1/4))<x<(1)/(2^(1/4))&&-(Sqrt[1-Sqrt[2] x^2])/(Sqrt[Root])<y<Sqrt[(1-Sqrt[2] x^2)/(Root)]
```

Coefficients can include real exact transcendental numbers:

```wolfram
CylindricalDecomposition[E^Pi x^2+Pi^E y^2<=77&&Log[2]x+Log[3]y==2,{x,y}]
(* Output *)
(2 π^ℯ Log[2])/(π^ℯ Log[2]^2+ℯ^π Log[3]^2)-(Sqrt[-4 ℯ^π π^ℯ Log[3]^2+77 π^ℯ Log[2]^2 Log[3]^2+77 ℯ^π Log[3]^4])/(π^ℯ Log[2]^2+ℯ^π Log[3]^2)<=x<=(2 π^ℯ Log[2])/(π^ℯ Log[2]^2+ℯ^π Log[3]^2)+(Sqrt[-4 ℯ^π π^ℯ Log[3]^2+77 π^ℯ Log[2]^2 Log[3]^2+77 ℯ^π Log[3]^4])/(π^ℯ Log[2]^2+ℯ^π Log[3]^2)&&y==(2-x Log[2])/(Log[3])
```

Functions can be real algebraic:

```wolfram
CylindricalDecomposition[Sqrt[x-y]<=Root[#^3-y #+x&,1],{x,y}]
(* Output *)
(x<0&&Root[x^2-x^3+5 x^2 #1-8 x #1^2+4 #1^3&,1]<=y<=x)||(x==0&&y==0)
```

#### Topological Operations

Find decompositions of the solution region and its boundary:

[Graphics]

```wolfram
CylindricalDecomposition[x^4>=x^2+y^2,{x,y}]
(* Output *)
(x<-1&&-Sqrt[-x^2+x^4]<=y<=Sqrt[-x^2+x^4])||(x==-1&&y==0)||(x==0&&y==0)||(x==1&&y==0)||(x>1&&-Sqrt[-x^2+x^4]<=y<=Sqrt[-x^2+x^4])
```

```wolfram
CylindricalDecomposition[x^4>=x^2+y^2,{x,y},"Boundary"]
(* Output *)
(x<-1&&(y==-Sqrt[-x^2+x^4]||y==Sqrt[-x^2+x^4]))||(x==-1&&y==-Sqrt[-x^2+x^4])||(x==0&&y==-Sqrt[-x^2+x^4])||(x==1&&y==-Sqrt[-x^2+x^4])||(x>1&&(y==-Sqrt[-x^2+x^4]||y==Sqrt[-x^2+x^4]))
```

Find decompositions of the solution region and its closure:

[Graphics]

```wolfram
CylindricalDecomposition[x^4>x^2+y^2,{x,y}]
(* Output *)
(x<-1&&-Sqrt[-x^2+x^4]<y<Sqrt[-x^2+x^4])||(x>1&&-Sqrt[-x^2+x^4]<y<Sqrt[-x^2+x^4])
```

```wolfram
CylindricalDecomposition[x^4>x^2+y^2,{x,y},"Closure"]
(* Output *)
(x<-1&&(y==-Sqrt[-x^2+x^4]||-Sqrt[-x^2+x^4]<y<Sqrt[-x^2+x^4]||y==Sqrt[-x^2+x^4]))||(x==-1&&y==Sqrt[-x^2+x^4])||(x==1&&y==Sqrt[-x^2+x^4])||(x>1&&(y==-Sqrt[-x^2+x^4]||-Sqrt[-x^2+x^4]<y<Sqrt[-x^2+x^4]||y==Sqrt[-x^2+x^4]))
```

Find decompositions of the solution region and its interior:

[Graphics]

```wolfram
CylindricalDecomposition[x^4>=x^2+y^2,{x,y}]
(* Output *)
(x<-1&&-Sqrt[-x^2+x^4]<=y<=Sqrt[-x^2+x^4])||(x==-1&&y==0)||(x==0&&y==0)||(x==1&&y==0)||(x>1&&-Sqrt[-x^2+x^4]<=y<=Sqrt[-x^2+x^4])
```

```wolfram
CylindricalDecomposition[x^4>=x^2+y^2,{x,y},"Interior"]
(* Output *)
(x<-1&&-Sqrt[-x^2+x^4]<y<Sqrt[-x^2+x^4])||(x>1&&-Sqrt[-x^2+x^4]<y<Sqrt[-x^2+x^4])
```

Find decompositions of the solution region and its exterior:

[Graphics]

```wolfram
CylindricalDecomposition[x^4>=x^2+y^2,{x,y}]
(* Output *)
(x<-1&&-Sqrt[-x^2+x^4]<=y<=Sqrt[-x^2+x^4])||(x==-1&&y==0)||(x==0&&y==0)||(x==1&&y==0)||(x>1&&-Sqrt[-x^2+x^4]<=y<=Sqrt[-x^2+x^4])
```

```wolfram
CylindricalDecomposition[x^4>=x^2+y^2,{x,y},"Exterior"]
(* Output *)
(x<-1&&(y<-Sqrt[-x^2+x^4]||y>Sqrt[-x^2+x^4]))||(x==-1&&(y<-Sqrt[-x^2+x^4]||y>-Sqrt[-x^2+x^4]))||-1<x<0||(x==0&&(y<-Sqrt[-x^2+x^4]||y>-Sqrt[-x^2+x^4]))||0<x<1||(x==1&&(y<-Sqrt[-x^2+x^4]||y>-Sqrt[-x^2+x^4]))||(x>1&&(y<-Sqrt[-x^2+x^4]||y>Sqrt[-x^2+x^4]))
```

Find decompositions of the solution region and its connected components:

[Graphics]

```wolfram
CylindricalDecomposition[x^4>=x^2+y^2,{x,y}]
(* Output *)
(x<-1&&-Sqrt[-x^2+x^4]<=y<=Sqrt[-x^2+x^4])||(x==-1&&y==0)||(x==0&&y==0)||(x==1&&y==0)||(x>1&&-Sqrt[-x^2+x^4]<=y<=Sqrt[-x^2+x^4])
```

```wolfram
CylindricalDecomposition[x^4>=x^2+y^2,{x,y},"Components"]
(* Output *)
{(x<-1&&(y==-Sqrt[-x^2+x^4]||-Sqrt[-x^2+x^4]<y<Sqrt[-x^2+x^4]||y==Sqrt[-x^2+x^4]))||(x==-1&&y==-Sqrt[-x^2+x^4]),x==0&&y==-Sqrt[-x^2+x^4],(x==1&&y==-Sqrt[-x^2+x^4])||(x>1&&(y==-Sqrt[-x^2+x^4]||-Sqrt[-x^2+x^4]<y<Sqrt[-x^2+x^4]||y==Sqrt[-x^2+x^4]))}
```

### Options

#### WorkingPrecision

This computation takes a long time due to high degrees of numeric coefficients:

```wolfram
TimeConstrained[CylindricalDecomposition[E x^2 + Pi y^2 + E^77 z^2 <= Pi^33 && x^3 - 9y^3 + 7z^3 >= 1,{x,y,z},"Function"],60]
(* Output *)
$Aborted
```

This finds a decomposition using [WorkingPrecision](https://reference.wolfram.com/language/ref/WorkingPrecision.html)->1000, but the result may be incorrect:

```wolfram
CylindricalDecomposition[E x^2 + Pi y^2 + E^77 z^2 <= Pi^33 && x^3 - 9y^3 + 7z^3 >= 1,{x,y,z},"Function",WorkingPrecision->1000]//Timing
(* Output *)
{1.140625,CylindricalDecompositionFunction[...][x,y,z]}
```

### Applications

Find the connected components of an algebraic curve:

[Graphics]

```wolfram
CylindricalDecomposition[y^2==x (-1+x^4),{x,y}]
(* Output *)
(x==-1&&y==0)||(-1<x<0&&(y==-Sqrt[-x+x^5]||y==Sqrt[-x+x^5]))||(x==0&&y==0)||(x==1&&y==0)||(x>1&&(y==-Sqrt[-x+x^5]||y==Sqrt[-x+x^5]))
```

```wolfram
CylindricalDecomposition[y^2==x (-1+x^4),{x,y},"Components"]
(* Output *)
{(x==-1&&y==-Sqrt[-x+x^5])||(-1<x<0&&(y==-Sqrt[-x+x^5]||y==Sqrt[-x+x^5]))||(x==0&&y==-Sqrt[-x+x^5]),(x==1&&y==-Sqrt[-x+x^5])||(x>1&&(y==-Sqrt[-x+x^5]||y==Sqrt[-x+x^5]))}
```

Find the connected components of a region:

[Graphics]

```wolfram
CylindricalDecomposition[(x^2+y^2) (x (2+x)+y^2)<8 x y^2,{x,y}]
(* Output *)
(-2<x<0&&Root[2 x^3+x^4+(-6 x+2 x^2) #1^2+#1^4&,1]<y<Root[2 x^3+x^4+(-6 x+2 x^2) #1^2+#1^4&,2])||(0<x<(9)/(8)&&(Root[2 x^3+x^4+(-6 x+2 x^2) #1^2+#1^4&,1]<y<Root[2 x^3+x^4+(-6 x+2 x^2) #1^2+#1^4&,2]||Root[2 x^3+x^4+(-6 x+2 x^2) #1^2+#1^4&,3]<y<Root[2 x^3+x^4+(-6 x+2 x^2) #1^2+#1^4&,4]))
```

```wolfram
CylindricalDecomposition[(x^2+y^2) (x (2+x)+y^2)<8 x y^2,{x,y},"Components"]
(* Output *)
{-2<x<0&&Root[2 x^3+x^4+(-6 x+2 x^2) #1^2+#1^4&,1]<y<Root[2 x^3+x^4+(-6 x+2 x^2) #1^2+#1^4&,2],0<x<(9)/(8)&&Root[2 x^3+x^4+(-6 x+2 x^2) #1^2+#1^4&,1]<y<Root[2 x^3+x^4+(-6 x+2 x^2) #1^2+#1^4&,2],0<x<(9)/(8)&&Root[2 x^3+x^4+(-6 x+2 x^2) #1^2+#1^4&,3]<y<Root[2 x^3+x^4+(-6 x+2 x^2) #1^2+#1^4&,4]}
```

### Properties & Relations

Use [RegionPlot](https://reference.wolfram.com/language/ref/RegionPlot.html) to visualize 2D semialgebraic sets:

```wolfram
RegionPlot[x^2+y^2<1&&y>x,{x,-1,1},{y,-1,1}]
```

*([Graphics])*

Use [RegionPlot3D](https://reference.wolfram.com/language/ref/RegionPlot3D.html) to visualize 3D semialgebraic sets:

```wolfram
RegionPlot3D[x^2+y^2+z^2<1,{x,-1,1},{y,-1,1},{z,-1,1}]
```

*([Graphics3D])*

[Resolve](https://reference.wolfram.com/language/ref/Resolve.html) performs quantifier elimination and may avoid computing cylindrical decomposition:

```wolfram
CylindricalDecomposition[Exists[z,x^2+y^2+z^2<1&&y>x],{x,y}]
(* Output *)
(-1<x<=-(1)/(Sqrt[2])&&-Sqrt[1-x^2]<y<Sqrt[1-x^2])||(-(1)/(Sqrt[2])<x<(1)/(Sqrt[2])&&x<y<Sqrt[1-x^2])
```

```wolfram
Resolve[Exists[z,x^2+y^2+z^2<1&&y>x]]
(* Output *)
(x|y)∈Reals&&x-y<0&&x^2+y^2<1
```

[Reduce](https://reference.wolfram.com/language/ref/Reduce.html) in addition deals with different domains and transcendental functions:

```wolfram
Reduce[x^2-2y^2==1,{x,y},Complexes]
(* Output *)
y==-(Sqrt[-1+x^2])/(Sqrt[2])||y==(Sqrt[-1+x^2])/(Sqrt[2])
```

```wolfram
Reduce[x^2-2y^2==1&&x>0&&y>0,{x,y},Reals]
(* Output *)
x>1&&y==(Sqrt[-1+x^2])/(Sqrt[2])
```

```wolfram
Reduce[x^2-2y^2==1&&x>0&&y>0,{x,y},Integers]
(* Output *)
1∈Integers&&1>=1&&x==(1)/(2) ((3-2 Sqrt[2])^1+(3+2 Sqrt[2])^1)&&y==-((3-2 Sqrt[2])^1-(3+2 Sqrt[2])^1)/(2 Sqrt[2])
```

```wolfram
Reduce[Exp[x]-2y==1&&y==x+1,{x,y}]
(* Output *)
1∈Integers&&x==-(3)/(2)-ProductLog[1,-(1)/(2 ℯ^(3/2))]&&y==(1)/(2) (-1+ℯ^x)
```

Use [FindInstance](https://reference.wolfram.com/language/ref/FindInstance.html) to find points that satisfy equations and inequalities:

```wolfram
FindInstance[x^2+y^2<1,{x,y},5]
(* Output *)
{{x->-(220)/(251),y->(47)/(105)},{x->-(218)/(251),y->-(23)/(51)},{x->(16)/(251),y->-(22)/(51)},{x->(84)/(251),y->(31)/(54)},{x->(149)/(251),y->(47)/(63)}}
```

```wolfram
x^2+y^2<1/.%
(* Output *)
{True,True,True,True,True}
```

[SemialgebraicComponentInstances](https://reference.wolfram.com/language/ref/SemialgebraicComponentInstances.html) will give sample points in each cylinder:

```wolfram
SemialgebraicComponentInstances[x^2+y^2<1,{x,y}]
(* Output *)
{{x->-(1)/(2),y->-(1)/(Sqrt[2])},{x->-(1)/(2),y->(1)/(Sqrt[2])},{x->0,y->0},{x->0,y->-(1)/(Sqrt[2])},{x->0,y->(1)/(Sqrt[2])},{x->(1)/(2),y->-(1)/(Sqrt[2])},{x->(1)/(2),y->(1)/(Sqrt[2])},{x->-(1)/(Sqrt[2]),y->0},{x->(1)/(Sqrt[2]),y->0}}
```

[CylindricalDecomposition](https://reference.wolfram.com/language/ref/CylindricalDecomposition.html) merges several cylinders to get a more compact representation:

```wolfram
CylindricalDecomposition[x^2+y^2<1,{x,y}]
(* Output *)
-1<x<1&&-Sqrt[1-x^2]<y<Sqrt[1-x^2]
```

[GenericCylindricalDecomposition](https://reference.wolfram.com/language/ref/GenericCylindricalDecomposition.html) will compute the full-dimensional part only:

```wolfram
CylindricalDecomposition[x^2+y^2<=1||x==1,{x,y}]
(* Output *)
(x==-1&&y==0)||(-1<x<1&&-Sqrt[1-x^2]<=y<=Sqrt[1-x^2])||x==1
```

```wolfram
GenericCylindricalDecomposition[x^2+y^2<=1||x==1,{x,y}]
(* Output *)
{-1<x<1&&-Sqrt[1-x^2]<y<Sqrt[1-x^2],-1+x==0||-1+x^2+y^2==0}
```

The output and input are equal as sets:

```wolfram
set1=x^2+y^2<1&&y>x
(* Output *)
x^2+y^2<1&&y>x
```

```wolfram
set2=CylindricalDecomposition[set1,{x,y}]
(* Output *)
(-1<x<=-(1)/(Sqrt[2])&&-Sqrt[1-x^2]<y<Sqrt[1-x^2])||(-(1)/(Sqrt[2])<x<(1)/(Sqrt[2])&&x<y<Sqrt[1-x^2])
```

```wolfram
set3=CylindricalDecomposition[set1,{y,x}]
(* Output *)
(-(1)/(Sqrt[2])<y<=(1)/(Sqrt[2])&&-Sqrt[1-y^2]<x<y)||((1)/(Sqrt[2])<y<1&&-Sqrt[1-y^2]<x<Sqrt[1-y^2])
```

Points are simultaneously inside or outside of the sets:

```wolfram
{set1,set2,set3}/.{{x->1,y->1},{x->0,y->1/2}}
(* Output *)
{{False,False,False},{True,True,True}}
```

### Possible Issues

[CylindricalDecomposition](https://reference.wolfram.com/language/ref/CylindricalDecomposition.html) requires exact, infinite[Hyphen]precision input:

```wolfram
CylindricalDecomposition[1.2x^2+3.2y^2<1.0,{x,y}]
(* Output *)
CylindricalDecomposition[1.2 x^2+3.2 y^2<1.,{x,y}]
```

[Rationalize](https://reference.wolfram.com/language/ref/Rationalize.html) will convert inexact numbers to exact ones:

```wolfram
Rationalize[1.2x^2+3.2y^2<1.0,0]
(* Output *)
(6 x^2)/(5)+(16 y^2)/(5)<1
```

```wolfram
CylindricalDecomposition[%,{x,y}]
(* Output *)
-Sqrt[(5)/(6)]<x<Sqrt[(5)/(6)]&&-(1)/(4) Sqrt[5-6 x^2]<y<(1)/(4) Sqrt[5-6 x^2]
```

In general, the output can be in a nested and more compact form:

```wolfram
cd=CylindricalDecomposition[x^2+y^2+y^2<=1&&x^2+2y^2+3z^2>=4x y z+1,{x,y,z}]
(* Output *)
(x==-1&&y==0)||(-1<x<0&&-(Sqrt[1-x^2])/(Sqrt[2])<=y<=(Sqrt[1-x^2])/(Sqrt[2])&&(z<=(2 x y)/(3)-(1)/(3) Sqrt[3-3 x^2-6 y^2+4 x^2 y^2]||z>=(2 x y)/(3)+(1)/(3) Sqrt[3-3 x^2-6 y^2+4 x^2 y^2]))||(x==0&&(y==-(1)/(Sqrt[2])||(-(1)/(Sqrt[2])<y<(1)/(Sqrt[2])&&(z<=-(1)/(3) Sqrt[3-6 y^2]||z>=(1)/(3) Sqrt[3-6 y^2]))||y==(1)/(Sqrt[2])))||(0<x<1&&-(Sqrt[1-x^2])/(Sqrt[2])<=y<=(Sqrt[1-x^2])/(Sqrt[2])&&(z<=(2 x y)/(3)-(1)/(3) Sqrt[3-3 x^2-6 y^2+4 x^2 y^2]||z>=(2 x y)/(3)+(1)/(3) Sqrt[3-3 x^2-6 y^2+4 x^2 y^2]))||(x==1&&y==0)
```

Flatten the result into disjunctive normal form without splitting the inequalities:

```wolfram
ToDNF[e_Or]:=ToDNF/@e;
ToDNF[And[a___,b_Or,c___]]:=ToDNF[And[a,#,c]]&/@b;
ToDNF[e_]:=e;
```

```wolfram
ToDNF[cd]
(* Output *)
(x==-1&&y==0)||(-1<x<0&&-(Sqrt[1-x^2])/(Sqrt[2])<=y<=(Sqrt[1-x^2])/(Sqrt[2])&&z<=(2 x y)/(3)-(1)/(3) Sqrt[3-3 x^2-6 y^2+4 x^2 y^2])||(-1<x<0&&-(Sqrt[1-x^2])/(Sqrt[2])<=y<=(Sqrt[1-x^2])/(Sqrt[2])&&z>=(2 x y)/(3)+(1)/(3) Sqrt[3-3 x^2-6 y^2+4 x^2 y^2])||(x==0&&y==-(1)/(Sqrt[2]))||(x==0&&-(1)/(Sqrt[2])<y<(1)/(Sqrt[2])&&z<=-(1)/(3) Sqrt[3-6 y^2])||(x==0&&-(1)/(Sqrt[2])<y<(1)/(Sqrt[2])&&z>=(1)/(3) Sqrt[3-6 y^2])||(x==0&&y==(1)/(Sqrt[2]))||(0<x<1&&-(Sqrt[1-x^2])/(Sqrt[2])<=y<=(Sqrt[1-x^2])/(Sqrt[2])&&z<=(2 x y)/(3)-(1)/(3) Sqrt[3-3 x^2-6 y^2+4 x^2 y^2])||(0<x<1&&-(Sqrt[1-x^2])/(Sqrt[2])<=y<=(Sqrt[1-x^2])/(Sqrt[2])&&z>=(2 x y)/(3)+(1)/(3) Sqrt[3-3 x^2-6 y^2+4 x^2 y^2])||(x==1&&y==0)
```

### Neat Examples

Semialgebraic sets are quite general:

```wolfram
disk[{x0_,y0_},r_]:=(x-x0)^2+(y-y0)^2<r^2
```

```wolfram
disk[{-7/4,9/4},1/3]||disk[{9/4,9/4},1/3]||(disk[{0,0},5]&&!disk[{2,2},1]&&!disk[{-2,2},1]&&(!disk[{0,-1},Sqrt[5]]||disk[{0,1},3]))
(* Output *)
((7)/(4)+x)^2+(-(9)/(4)+y)^2<(1)/(9)||(-(9)/(4)+x)^2+(-(9)/(4)+y)^2<(1)/(9)||(x^2+y^2<25&&(-2+x)^2+(-2+y)^2>=1&&(2+x)^2+(-2+y)^2>=1&&(x^2+(1+y)^2>=5||x^2+(-1+y)^2<9))
```

```wolfram
RegionPlot[%,{x,-5,5},{y,-5,5}]
```

*([Graphics])*

## Tech Notes ▪The Representation of Solution Sets ▪Implementation notes: Algebra and Calculus

## Related Guides ▪Polynomial Systems ▪Computational Geometry ▪Theorem Proving ▪Polynomial Algebra

## History Introduced in 2003 (5.0) | Updated in 2017 (11.2) ▪ 2020 (12.2)
