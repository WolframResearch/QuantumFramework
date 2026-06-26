# DEigensystem | [SpanFromLeft]

> [DEigensystem](https://reference.wolfram.com/language/ref/DEigensystem.html)[*ℒ*[*u*[*x,y,…*]],*u*,*{x,y,…}*∈Ω,*n*]  — gives the `*n*` smallest magnitude eigenvalues and eigenfunctions for the linear differential operator `*ℒ*` over the region `Ω`.
> [DEigensystem](https://reference.wolfram.com/language/ref/DEigensystem.html)[*eqns*,*u*,*t*,*{x,y,…}*∈Ω,*n*] — gives the eigenvalues and eigenfunctions for solutions `*u*` of the time-dependent differential equations `*eqns*`.

## Details and Options

[DEigensystem](https://reference.wolfram.com/language/ref/DEigensystem.html) can compute eigenvalues and eigenfunctions for ordinary and partial differential operators with given boundary conditions.

[DEigensystem](https://reference.wolfram.com/language/ref/DEigensystem.html) gives lists `{{λ_1,…,λ_*n*},{*u*_1,…,*u*_*n*}}` of eigenvalues `λ_*i*` and eigenfunctions `*u*_*i*`.

An eigenvalue and eigenfunction pair `{λ_*i*,*u*_*i*}` for the differential operator `*ℒ*` satisfy `*ℒ*[*u*_*i*[*x,y,…*]]==λ_*i* *u*_*i*[*x,y,…*]`.

Homogeneous [DirichletCondition](https://reference.wolfram.com/language/ref/DirichletCondition.html) or [NeumannValue](https://reference.wolfram.com/language/ref/NeumannValue.html) boundary conditions may be included. Inhomogeneous boundary conditions will be replaced with corresponding homogeneous boundary conditions.

When no boundary condition is specified on the boundary `∂Ω`, then this is equivalent to specifying a Neumann 0 condition.

The equations `*eqns*` are specified as in [DSolve](https://reference.wolfram.com/language/ref/DSolve.html).

[N](https://reference.wolfram.com/language/ref/N.html)[DEigensystem[…]] calls [NDEigensystem](https://reference.wolfram.com/language/ref/NDEigensystem.html) for eigensystems that cannot be computed symbolically.

The following options can be given:

| [Assumptions](https://reference.wolfram.com/language/ref/Assumptions.html) | [$Assumptions](https://reference.wolfram.com/language/ref/$Assumptions.html) | assumptions on parameters |
| --- | --- | --- |
| [Method](https://reference.wolfram.com/language/ref/Method.html) | [Automatic](https://reference.wolfram.com/language/ref/Automatic.html) | what method to use |

Eigenfunctions are not automatically normalized. The setting [Method](https://reference.wolfram.com/language/ref/Method.html)->"Normalize" can be used to give normalized eigenfunctions.

## Examples

### Basic Examples

Find the 4 smallest eigenvalues and eigenfunctions of the Laplacian operator on `[0,π]`:

```wolfram
DEigensystem[{-Laplacian[u[x],{x}],DirichletCondition[u[x]==0,True]},u[x],{x,0,π},4]
(* Output *)
{{1,4,9,16},{Sin[x],Sin[2 x],Sin[3 x],Sin[4 x]}}
```

Visualize the eigenfunctions:

```wolfram
Plot[Evaluate[%[[2]]],{x,0,π}]
```

*([Graphics])*

Compute the first 6 eigenfunctions for a circular membrane with the edges clamped:

```wolfram
{vals,funs}=DEigensystem[{-Laplacian[u[x,y],{x,y}],DirichletCondition[u[x,y]==0,True]},u[x,y],{x,y}∈Disk[],6];
```

```wolfram
vals//N
(* Output *)
{5.783185962946785,14.681970642123892,14.681970642123892,26.374616427163392,26.374616427163392,30.471262343662094}
```

Visualize the eigenfunctions:

```wolfram
Table[Plot3D[funs[[i]]//N//Evaluate,{x,y}∈Disk[],PlotRange->All,PlotLabel->vals[[i]],PlotTheme->"Minimal"],{i,Length[vals]}]
(* Output *)
![image](img/image_001.png)
```

### Scope

#### 1D

Specify a Laplacian operator:

```wolfram
ℒ=-Laplacian[u[x],{x}];
```

Specify homogeneous Dirichlet boundary conditions:

```wolfram
ℬ=DirichletCondition[u[x]==0,True];
```

Find the 5 smallest eigenvalues and eigenfunctions in an interval:

```wolfram
{vals,funs}=DEigensystem[{ℒ,ℬ},u[x],{x,0,π},5];
```

```wolfram
vals
(* Output *)
{1,4,9,16,25}
```

Visualize the eigenfunctions:

```wolfram
Plot[Evaluate[funs],{x,0,π}]
```

*([Graphics])*

Specify a Laplacian operator:

```wolfram
ℒ=-Laplacian[u[x],{x}];
```

Specify homogeneous Neumann boundary conditions:

```wolfram
ℬ=NeumannValue[0,True];
```

Find the 5 smallest eigenvalues and eigenfunctions in an interval:

```wolfram
{vals,funs}=DEigensystem[ℒ+ℬ,u[x],{x,0,π},5];
```

```wolfram
vals
(* Output *)
{0,1,4,9,16}
```

This is equivalent:

```wolfram
{vals,funs}=DEigensystem[ℒ,u[x],{x,0,π},5];
```

```wolfram
vals
(* Output *)
{0,1,4,9,16}
```

Specify a Laplacian operator:

```wolfram
ℒ=-Laplacian[u[x],{x}];
```

Specify a homogeneous Dirichlet boundary condition:

```wolfram
ℬ1=DirichletCondition[u[x]==0,x==0];
```

Specify a homogeneous Neumann boundary condition:

```wolfram
ℬ2=NeumannValue[0,x==π];
```

Find the 5 smallest eigenvalues and eigenfunctions in an interval:

```wolfram
{vals,funs}=DEigensystem[{ℒ+ℬ2,ℬ1},u[x],{x,0,π},5];
```

```wolfram
vals
(* Output *)
{(1)/(4),(9)/(4),(25)/(4),(49)/(4),(81)/(4)}
```

Visualize the eigenfunctions:

```wolfram
Plot[Evaluate[funs],{x,0,π}]
```

*([Graphics])*

Specify a Laplacian operator:

```wolfram
ℒ=-Laplacian[u[x],{x}];
```

Specify a homogeneous Dirichlet boundary condition:

```wolfram
ℬ1=DirichletCondition[u[x]==0,x==π];
```

Specify a homogeneous non-zero Neumann boundary condition:

```wolfram
ℬ2=NeumannValue[-u[x]/3,x==0];
```

Find the 5 smallest eigenvalues and eigenfunctions in an interval:

```wolfram
{vals,funs}=DEigensystem[{ℒ+ℬ2,ℬ1},u[x],{x,0,π},5];
```

The eigenvalues are roots of a transcendental equation:

```wolfram
vals[[1]]//
(* Output *)
Root
```

Visualize the eigenfunctions:

```wolfram
Plot[Evaluate[funs],{x,0,π}]
```

*([Graphics])*

Specify an Airy operator:

```wolfram
ℒ=-Laplacian[u[x],{x}]+x u[x];
```

Specify a homogeneous Dirichlet boundary condition:

```wolfram
ℬ=DirichletCondition[u[x]==0,True];
```

Find the 5 smallest eigenvalues and eigenfunctions in an interval:

```wolfram
{vals,funs}=DEigensystem[{ℒ,ℬ},u[x],{x,0,1},5];
```

The eigenvalues are roots of a transcendental equation:

```wolfram
vals[[1]]//
(* Output *)
Root
```

Visualize the eigenfunctions:

```wolfram
Plot[Evaluate[funs],{x,0,1}]
```

*([Graphics])*

Specify an Airy operator:

```wolfram
ℒ=-Laplacian[u[x],{x}]+x u[x];
```

Specify a homogeneous Neumann boundary condition:

```wolfram
ℬ=NeumannValue[0,True];
```

Find the 5 smallest eigenvalues and eigenfunctions in an interval:

```wolfram
{vals,funs}=DEigensystem[{ℒ+ℬ},u[x],{x,0,1},5];
```

The eigenvalues are roots of a transcendental equation:

```wolfram
vals[[1]]//
(* Output *)
Root
```

Visualize the eigenfunctions:

```wolfram
Plot[Evaluate[funs],{x,0,1}]
```

*([Graphics])*

Find symbolic expressions for the eigenvalues and eigenfunctions of a Laplace operator:

```wolfram
{vals,funs}=DEigensystem[{-Laplacian[u[x],{x}],DirichletCondition[u[x]==0,True]},u[x],{x,a,b},3];
```

Symbolic eigenvalues:

```wolfram
vals
(* Output *)
{(π^2)/((-a+b)^2),(4 π^2)/((-a+b)^2),(9 π^2)/((-a+b)^2)}
```

Symbolic eigenfunctions:

```wolfram
funs
(* Output *)
{Sin[(π (-a+x))/(-a+b)],Sin[(2 π (-a+x))/(-a+b)],Sin[(3 π (-a+x))/(-a+b)]}
```

Enter the quantum harmonic operator:

```wolfram
qho=-(ℏ^2)/(2m)Laplacian[u[x],{x}]+(m ω^2)/(2)x^2 u[x];
```

Find symbolic expressions for the eigenvalues and eigenfunctions on the real line:

```wolfram
DEigensystem[qho,u[x],{x,-∞,∞},4,Assumptions->ℏ>0 && m>0&&ω>0]
(* Output *)
{{(ω ℏ)/(2),(3 ω ℏ)/(2),(5 ω ℏ)/(2),(7 ω ℏ)/(2)},{ℯ^(-(m x^2 ω)/(2 ℏ)),2 ℯ^(-(m x^2 ω)/(2 ℏ)) x Sqrt[(m ω)/(ℏ)],ℯ^(-(m x^2 ω)/(2 ℏ)) (-2+(4 m x^2 ω)/(ℏ)),ℯ^(-(m x^2 ω)/(2 ℏ)) (-12 x Sqrt[(m ω)/(ℏ)]+8 x^3 ((m ω)/(ℏ))^(3/2))}}
```

Specify a heat equation with homogeneous Dirichlet boundary conditions:

```wolfram
{ℒ,ℬ}={D[u[t,x],t]==Laplacian[u[t,x],{x}],DirichletCondition[u[t,x]==0,True]};
```

Find the 4 smallest eigenvalues and eigenfunctions:

```wolfram
{vals,funs}=DEigensystem[{ℒ,ℬ},u[t,x],t,{x,0,π},4]
(* Output *)
{{-1,-4,-9,-16},{ℯ^-t Sin[x],ℯ^(-4 t) Sin[2 x],ℯ^(-9 t) Sin[3 x],ℯ^(-16 t) Sin[4 x]}}
```

Visualize the eigenfunctions:

```wolfram
Table[Plot3D[funs[[i]]//Evaluate,{x,-3,3},{t,0,1/3},PlotRange->All,Ticks->False,Mesh->False],{i,4}]
(* Output *)
{[Graphics3D],[Graphics3D],[Graphics3D],[Graphics3D]}
```

#### 2D

Specify a Laplacian operator with homogeneous Dirichlet boundary conditions:

```wolfram
{ℒ,ℬ}={-Laplacian[u[x,y],{x,y}],DirichletCondition[u[x,y]==0,True]};
```

Find the 9 smallest eigenvalues and eigenfunctions in a rectangle:

```wolfram
{vals,funs}=DEigensystem[{ℒ,ℬ},u[x,y],{x,0,π},{y,0,π},9];
```

```wolfram
vals
(* Output *)
{2,5,5,8,10,10,13,13,17}
```

Visualize the eigenfunctions:

```wolfram
Plot3D[#,{x,0,π},{y,0,π}]&/@funs
(* Output *)
![image](img/image_003.png)
```

Specify a Laplacian operator with homogeneous Neumann boundary conditions:

```wolfram
ℒ=-Laplacian[u[x,y],{x,y}]+NeumannValue[0,True];
```

Find the 4 smallest eigenvalues and eigenfunctions in a rectangle:

```wolfram
{vals,funs}=DEigensystem[ℒ,u[x,y],{x,0,π},{y,0,π},4];
```

```wolfram
vals
(* Output *)
{0,1,1,2}
```

```wolfram
funs
(* Output *)
{1,Cos[x],Cos[y],Cos[x] Cos[y]}
```

Visualize the eigenfunctions:

```wolfram
ContourPlot[#,{x,0,π},{y,0,π}]&/@funs
(* Output *)
{[Graphics],[Graphics],[Graphics],[Graphics]}
```

Specify a Laplacian operator:

```wolfram
ℒ=-Laplacian[u[x,y],{x,y}];
```

Specify homogeneous Dirichlet boundary conditions:

```wolfram
ℬ=DirichletCondition[u[x,y]==0,True];
```

Find the 4 smallest eigenvalues and eigenfunctions of the operator in a unit disk:

```wolfram
{vals,funs}=DEigensystem[{ℒ,ℬ},u[x,y],{x,y}∈Disk[],4];
```

```wolfram
vals//
(* Output *)
{(0)^(2),(1)^(2),(1)^(2),(2)^(2)}
```

Visualize the eigenfunctions:

```wolfram
ContourPlot[#,{x,y}∈Disk[],Contours->10]&/@N[funs]
(* Output *)
![image](img/image_005.png)
```

Specify a quantum harmonic oscillator operator:

```wolfram
ℒ=-Laplacian[u[x,y],{x,y}]+(x^2+y^2)u[x,y];
```

Specify homogeneous Dirichlet boundary conditions:

```wolfram
ℬ=DirichletCondition[u[x,y]==0,True];
```

Find the 6 smallest eigenvalues and eigenfunctions of the operator on the plane:

```wolfram
{vals,funs}=DEigensystem[{ℒ,ℬ},u[x,y],{x,y}∈FullRegion[2],6];
```

```wolfram
vals
(* Output *)
{1,3,3,5,5,5}
```

Visualize the eigenfunctions:

```wolfram
Table[ContourPlot[funs[[i]],{x,y}∈Rectangle[{-3,-3},{3,3}]],{i,6}]
(* Output *)
{[Graphics],[Graphics],[Graphics],[Graphics],[Graphics],[Graphics]}
```

Specify a Laplacian operator:

```wolfram
ℒ=-Laplacian[u[x,y],{x,y}];
```

Specify homogeneous Dirichlet boundary conditions:

```wolfram
ℬ=DirichletCondition[u[x,y]==0,True];
```

Find the 6 smallest eigenvalues and eigenfunctions of the operator in a triangle:

```wolfram
{vals,funs}=DEigensystem[{ℒ,ℬ},u[x,y],{x,y}∈Triangle[],6];
```

```wolfram
vals
(* Output *)
{5 π^2,10 π^2,13 π^2,17 π^2,20 π^2,25 π^2}
```

Visualize the eigenfunctions:

```wolfram
Table[Plot3D[funs[[i]],{x,y}∈Triangle[],Boxed->False,Axes->False],{i,6}]
(* Output *)
![image](img/image_007.png)
```

Specify a Laplacian operator:

```wolfram
ℒ=-Laplacian[u[x,y],{x,y}];
```

Specify homogeneous Dirichlet boundary conditions:

```wolfram
ℬ=DirichletCondition[u[x,y]==0,True];
```

Find the 4 smallest eigenvalues and eigenfunctions of the operator in a sector of a disk:

```wolfram
{vals,funs}=DEigensystem[{ℒ,ℬ},u[x,y],{x,y}∈Disk[{0,0},1,{0,Pi/5}],4];
```

```wolfram
vals//
(* Output *)
{(5)^(2),(5)^(2),(10)^(2),(5)^(2)}
```

Visualize the eigenfunctions:

```wolfram
ContourPlot[#,{x,y}∈Disk[{0,0},1,{0,Pi/5}]]&/@N[funs]//Quiet
(* Output *)
![image](img/image_009.png)
```

#### 3D

Specify a Laplacian operator:

```wolfram
ℒ=-Laplacian[u[x,y,z],{x,y,z}];
```

Specify homogeneous Dirichlet boundary conditions:

```wolfram
ℬ=DirichletCondition[u[x,y,z]==0,True];
```

Find the 7 smallest eigenvalues and eigenfunctions in a cuboid:

```wolfram
{vals,funs}=DEigensystem[{ℒ,ℬ},u[x,y,z],{x,y,z}∈Cuboid[{2,1,1},{4,2,3}],7];
```

```wolfram
vals
(* Output *)
{(3 π^2)/(2),(9 π^2)/(4),(9 π^2)/(4),3 π^2,(7 π^2)/(2),(7 π^2)/(2),(17 π^2)/(4)}
```

Visualize an eigenfunction:

```wolfram
ContourPlot3D[Evaluate[funs[[7]]],{x,2,4},{y,1,2},{z,1,3},Boxed->False,Axes->False]
(* Output *)
![image](img/image_011.png)
```

Specify a Laplacian operator:

```wolfram
ℒ=-Laplacian[u[x,y,z],{x,y,z}];
```

Specify homogeneous Dirichlet boundary conditions:

```wolfram
ℬ=DirichletCondition[u[x,y,z]==0,True];
```

Find the 7 smallest eigenvalues and eigenfunctions in a cylinder:

```wolfram
{vals,funs}=DEigensystem[{ℒ,ℬ},u[x,y,z],{x,y,z}∈Cylinder[{{0,0,0},{0,0,5}},2],7];
```

```wolfram
vals[[1;;5]]//
(* Output *)
{(π^(2))/(25)+((0)^(2))/(4),(4 π^(2))/(25)+((0)^(2))/(4),(π^(2))/(25)+((1)^(2))/(4),(π^(2))/(25)+((1)^(2))/(4),(9 π^(2))/(25)+((0)^(2))/(4)}
```

Visualize an eigenfunction:

```wolfram
DensityPlot3D[funs[[7]]//N//Evaluate,{x,y,z}∈Cylinder[{{0,0,0},{0,0,5}},2],Boxed->False,Axes->False]
(* Output *)
![image](img/image_013.png)
```

Specify a Laplacian operator:

```wolfram
ℒ=-Laplacian[u[x,y,z],{x,y,z}];
```

Specify homogeneous Dirichlet boundary conditions:

```wolfram
ℬ=DirichletCondition[u[x,y,z]==0,True];
```

Find the 7 smallest eigenvalues and eigenfunctions in a ball:

```wolfram
{vals,funs}=DEigensystem[{ℒ,ℬ},u[x,y,z],{x,y,z}∈Ball[{0,0,0},2],7];
```

Visualize an eigenfunction:

```wolfram
DensityPlot3D[funs[[7]]//N//Evaluate,{x,y,z}∈Ball[{0,0,0},2],Boxed->False,Axes->False,ColorFunction->Hue]
```

*([Graphics3D])*

Specify a Laplacian operator:

```wolfram
ℒ=-Laplacian[u[x,y,z],{x,y,z}];
```

Specify homogeneous Dirichlet boundary conditions:

```wolfram
ℬ=DirichletCondition[u[x,y,z]==0,True];
```

Find the 7 smallest eigenvalues and eigenfunctions in a prism:

```wolfram
{vals,funs}=DEigensystem[{ℒ,ℬ},u[x,y,z],{x,y,z}∈Prism[{{0,0,0},{1,0,0},{0,1,0},{0,0,2},{1,0,2},{0,1,2}}],7];
```

```wolfram
vals
(* Output *)
{(21 π^2)/(4),6 π^2,(29 π^2)/(4),9 π^2,(41 π^2)/(4),11 π^2,(45 π^2)/(4)}
```

Visualize an eigenfunction:

```wolfram
DensityPlot3D[funs[[7]]//N//Evaluate,{x,y,z}∈Prism[{{0,0,0},{1,0,0},{0,1,0},{0,0,2},{1,0,2},{0,1,2}}],Boxed->False,ColorFunction->Hue,PlotPoints-> 50,Axes->False]
```

*([Graphics3D])*

Specify a quantum harmonic oscillator operator:

```wolfram
ℒ=-Laplacian[u[x,y,z],{x,y,z}]+2(x^2+y^2+z^2)u[x,y,z];
```

Specify homogeneous Dirichlet boundary conditions:

```wolfram
ℬ=DirichletCondition[u[x,y,z]==0,True];
```

Find the 8 smallest eigenvalues and eigenfunctions throughout space:

```wolfram
{vals,funs}=DEigensystem[{ℒ,ℬ},u[x,y,z],{x,y,z}∈FullRegion[3],8];
```

```wolfram
vals
(* Output *)
{Sqrt[2],3 Sqrt[2],3 Sqrt[2],3 Sqrt[2],5 Sqrt[2],5 Sqrt[2],5 Sqrt[2],5 Sqrt[2]}
```

Visualize the eigenfunctions:

```wolfram
Table[DensityPlot3D[funs[[i]]//N//Evaluate,{x,y,z}∈Ball[{0,0,0},3],Boxed->False,ColorFunction->Hue,PlotPoints-> 50,Axes->False],{i,8}]
(* Output *)
![image](img/image_015.png)
```

### Options

#### Assumptions

Use [Assumptions](https://reference.wolfram.com/language/ref/Assumptions.html) to simplify a result:

```wolfram
DEigensystem[-Laplacian[u[x,y,z],{x,y,z}]+2ω^2(x^2+y^2+z^2)u[x,y,z],u[x,y,z],{x,y,z}∈FullRegion[3],4,Assumptions->ω>0]
(* Output *)
{{Sqrt[2] ω,3 Sqrt[2] ω,3 Sqrt[2] ω,3 Sqrt[2] ω},{ℯ^(-(x^2 ω)/(Sqrt[2])-(y^2 ω)/(Sqrt[2])-(z^2 ω)/(Sqrt[2])),2 2^(1/4) ℯ^(-(x^2 ω)/(Sqrt[2])-(y^2 ω)/(Sqrt[2])-(z^2 ω)/(Sqrt[2])) x Sqrt[ω],2 2^(1/4) ℯ^(-(x^2 ω)/(Sqrt[2])-(y^2 ω)/(Sqrt[2])-(z^2 ω)/(Sqrt[2])) y Sqrt[ω],2 2^(1/4) ℯ^(-(x^2 ω)/(Sqrt[2])-(y^2 ω)/(Sqrt[2])-(z^2 ω)/(Sqrt[2])) z Sqrt[ω]}}
```

Without the option, an equivalent but more complicated answer would be returned:

```wolfram
DEigensystem[-Laplacian[u[x,y,z],{x,y,z}]+2ω^2(x^2+y^2+z^2)u[x,y,z],u[x,y,z],{x,y,z}∈FullRegion[3],4]
(* Output *)
{{Sqrt[2] Sqrt[ω^2],3 Sqrt[2] Sqrt[ω^2],3 Sqrt[2] Sqrt[ω^2],3 Sqrt[2] Sqrt[ω^2]},{ℯ^(-(x^2 Sqrt[ω^2])/(Sqrt[2])-(y^2 Sqrt[ω^2])/(Sqrt[2])-(z^2 Sqrt[ω^2])/(Sqrt[2])),2 2^(1/4) ℯ^(-(x^2 Sqrt[ω^2])/(Sqrt[2])-(y^2 Sqrt[ω^2])/(Sqrt[2])-(z^2 Sqrt[ω^2])/(Sqrt[2])) x (ω^2)^(1/4),2 2^(1/4) ℯ^(-(x^2 Sqrt[ω^2])/(Sqrt[2])-(y^2 Sqrt[ω^2])/(Sqrt[2])-(z^2 Sqrt[ω^2])/(Sqrt[2])) y (ω^2)^(1/4),2 2^(1/4) ℯ^(-(x^2 Sqrt[ω^2])/(Sqrt[2])-(y^2 Sqrt[ω^2])/(Sqrt[2])-(z^2 Sqrt[ω^2])/(Sqrt[2])) z (ω^2)^(1/4)}}
```

#### Method

Normalize the eigenfunctions for a differential operator:

```wolfram
DEigensystem[{-Laplacian[u[x],{x}],DirichletCondition[u[x]==0,True]},u[x],{x,0,Pi},3,Method->"Normalize"]
(* Output *)
{{1,4,9},{Sqrt[(2)/(π)] Sin[x],Sqrt[(2)/(π)] Sin[2 x],Sqrt[(2)/(π)] Sin[3 x]}}
```

Verify that the eigenfunctions are normalized to unity:

```wolfram
Integrate[%[[2]]^2,{x,0,Pi}]
(* Output *)
{1,1,1}
```

### Applications

Compute the first three terms in the eigenfunction expansion of the function $f(x)=x^{2} (\pi-x)^{3}$ with respect to the basis provided by a 1D Laplacian with a Dirichlet condition on the interval $\{0,\pi \}$:

```wolfram
basis=DEigensystem[{-Laplacian[u[x],{x}]+u[x],DirichletCondition[u[x]==0,True]},u[x],{x,0,π},3,Method->"Normalize"][[2]]
(* Output *)
{Sqrt[(2)/(π)] Sin[x],Sqrt[(2)/(π)] Sin[2 x],Sqrt[(2)/(π)] Sin[3 x]}
```

Compute the Fourier coefficients:

```wolfram
f[x_]:=x^2(π-x)^3
```

```wolfram
coeffs=Table[Integrate[f[x] basis[[i]],{x,0,Pi}],{i,3}]
(* Output *)
{-2 Sqrt[2 π] (-12+π^2),-(1)/(2) Sqrt[(π)/(2)] (-15+π^2),(2)/(81) Sqrt[2 π] (4-3 π^2)}
```

The required eigenfunction expansion is:

```wolfram
eigexp[x_]=Sum[coeffs[[i]]basis[[i]],{i,3}]
(* Output *)
-4 (-12+π^2) Sin[x]-(1)/(2) (-15+π^2) Sin[2 x]+(4)/(81) (4-3 π^2) Sin[3 x]
```

Compare the function with its eigenfunction expansion:

```wolfram
Plot[{f[x],eigexp[x]}//Evaluate,{x,0,Pi}]
```

*([Graphics])*

Build a solution of the heat equation by using a linear combination of eigenfunctions for the heat equation with a Dirichlet condition:

```wolfram
eqns={D[u[t,x],t]==Laplacian[u[t,x],{x}],DirichletCondition[u[t,x]==0,True]};
```

```wolfram
eigfuns=DEigensystem[eqns,u[t,x],t,{x,0,Pi},5][[2]]
(* Output *)
{ℯ^-t Sin[x],ℯ^(-4 t) Sin[2 x],ℯ^(-9 t) Sin[3 x],ℯ^(-16 t) Sin[4 x],ℯ^(-25 t) Sin[5 x]}
```

Form a linear combination of eigenfunctions:

```wolfram
sol[t_,x_]=Sum[RandomInteger[{1,10}] eigfuns[[i]],{i,5}]
(* Output *)
6 ℯ^-t Sin[x]+4 ℯ^(-4 t) Sin[2 x]+7 ℯ^(-9 t) Sin[3 x]+7 ℯ^(-16 t) Sin[4 x]+5 ℯ^(-25 t) Sin[5 x]
```

Verify that this is indeed a solution of the heat equation:

```wolfram
eqns[[1]]/.{u-> sol}
(* Output *)
True
```

The solution satisfies the homogeneous Dirichlet condition:

```wolfram
{sol[t,0],sol[t,π]}
(* Output *)
{0,0}
```

Visualize the solution:

```wolfram
Plot3D[sol[t,x],{x,0,Pi},{t,0,0.7},PlotRange->All]
```

*([Graphics3D])*

Experimentally, a CO molecule oscillates about its equilibrium length with an effective spring constant of $k=kN/m$. The oscillations will be governed by the quantum harmonic oscillator equation. In the following, $m$ is the reduced mass of the molecule, $\omega=\sqrt{k/m}$ is the natural frequency, $x$ is the displacement from the equilibrium position, and $\hbar$ is the reduced Planck's constant:

```wolfram
qho=-(ℏ^2)/(2m)Laplacian[u[x],{x}]+(m ω^2)/(2)x^2 u[x];
```

Compute the eigenvalues--the energies of the respective states--and normalized eigenfunctions:

```wolfram
{ℰs,efuns}=DEigensystem[qho,u[x],{x,-∞,∞},4,Assumptions->ℏ>0 && m>0&&ω>0,Method->"Normalize"]
(* Output *)
{{(ω ℏ)/(2),(3 ω ℏ)/(2),(5 ω ℏ)/(2),(7 ω ℏ)/(2)},{(ℯ^(-(m x^2 ω)/(2 ℏ)) (m ω ℏ)^(1/4))/(π^(1/4) Sqrt[ℏ]),(Sqrt[2] ℯ^(-(m x^2 ω)/(2 ℏ)) x ((m ω)/(ℏ))^(3/4))/(π^(1/4)),(ℯ^(-(m x^2 ω)/(2 ℏ)) (-2+(4 m x^2 ω)/(ℏ)) (m ω ℏ)^(1/4))/(2 Sqrt[2] π^(1/4) Sqrt[ℏ]),(ℯ^(-(m x^2 ω)/(2 ℏ)) (-12 x Sqrt[(m ω)/(ℏ)]+8 x^3 ((m ω)/(ℏ))^(3/2)) (m ω ℏ)^(1/4))/(4 Sqrt[3] π^(1/4) Sqrt[ℏ])}}
```

If the particle is in an equal superposition of the four states, the wavefunction has the following form:

```wolfram
ψ[x_,t_]=Total[MapThread[(1)/(2)Exp[I  t #1/ℏ]#2&,{ℰs,efuns}]]
(* Output *)
(ℯ^((3 ⅈ t ω)/(2)-(m x^2 ω)/(2 ℏ)) x ((m ω)/(ℏ))^(3/4))/(Sqrt[2] π^(1/4))+(ℯ^((ⅈ t ω)/(2)-(m x^2 ω)/(2 ℏ)) (m ω ℏ)^(1/4))/(2 π^(1/4) Sqrt[ℏ])+(ℯ^((7 ⅈ t ω)/(2)-(m x^2 ω)/(2 ℏ)) (-12 x Sqrt[(m ω)/(ℏ)]+8 x^3 ((m ω)/(ℏ))^(3/2)) (m ω ℏ)^(1/4))/(8 Sqrt[3] π^(1/4) Sqrt[ℏ])+(ℯ^((5 ⅈ t ω)/(2)-(m x^2 ω)/(2 ℏ)) (-2+(4 m x^2 ω)/(ℏ)) (m ω ℏ)^(1/4))/(4 Sqrt[2] π^(1/4) Sqrt[ℏ])
```

Compute $m$, $\omega$, and $\hbar$ using base units of atomic mass units, femtoseconds, and picometers, which give values close to order unity:

```wolfram
m=QuantityMagnitude[(carbon["atomic mass"]oxygen["atomic mass"])/(carbon["atomic mass"]+oxygen["atomic mass"]),"AtomicMassUnits"]
(* Output *)
6.8605494109246697575
```

```wolfram
ω=Sqrt[QuantityMagnitude[Quantity[1.86, "Kilonewtons"/"Meters"],"AtomicMassUnit"/"Femtoseconds"^2]/m]
(* Output *)
0.40406615832190024
```

```wolfram
ℏ=QuantityMagnitude[Quantity[1.,"ReducedPlanckConstant"],"AtomicMassUnit"*"Picometers"^2/"Femtoseconds"]
(* Output *)
63.507799250537346
```

The response of the eigenfunctions to the potential energy $\frac{1}{2} m \omega^{2}x^{2}$ can be visualized by rescaling them to fit in the band $\mathcal{E}_{n}\pm \frac{\hbar \omega}{2}$:

```wolfram
Show[Plot[Evaluate[ℏ ω efuns+ℰs],{x,-18,18}],Plot[Evaluate[Append[ℰs,(1)/(2)m ω^2x^2]],{x,-18,18},PlotStyle->Dashed]]
```

*([Graphics])*

The probability density function of the displacement from equilibrium is given by $\psi^{*}(x,t) \psi(x,t)$:

```wolfram
ρ[x_, t_] =FullSimplify[ComplexExpand[Conjugate[ψ[x, t]] ψ[x, t]]]
(* Output *)
ℯ^(-0.043650006411243295 x^2) (0.04420266839268609+0.003858893517469619 x^2-0.00011229381785190271 x^4+3.2677505794523595×10^-6 x^6+x (0.0051003936277304335+0.0010749634235769716 x^2) Cos[0.4040661583219002 t]-0.04167467542267814 Cos[0.8081323166438004 t]+x ((0.015080821839470845-0.0017554079199459203 x^2+0.00003831178347999334 x^4) Cos[0.4040661583219003 t]+x (-0.002663347027583451+0.0001833750366707988 x^2) Cos[0.8081323166438004 t]+(-0.021327502777112036+0.0006206304219711663 x^2) Cos[1.2121984749657007 t]+8.673617379884035×10^-19 Cos[2.4243969499314013 t]))
```

As a probability distribution, the integral of $\rho$ over the reals is 1 for all $t$:

```wolfram
Chop[Integrate[ρ[x,t],{x,-∞,∞}]]
(* Output *)
1.
```

Visualize the probability density over time:

```wolfram
Animate[Plot[ρ[x,t],{x,-25,25},PlotRange->{0,.16},PlotTheme->"Detailed",FrameLabel->Style[x,FontSize->Larger],PlotLegends->Placed[{HoldForm[ρ][x,NumberForm[t,{3,2}]]},Above]],{t,0.,5.7},AnimationRate->1,SaveDefinitions->True,Alignment->Center]
```

### Properties & Relations

Use [NDEigensystem](https://reference.wolfram.com/language/ref/NDEigensystem.html) to find numerical eigenvalues and eigenvectors:

```wolfram
{ℒ,ℬ}={-Laplacian[u[x],{x}],DirichletCondition[u[x]==0,True]};
```

Find exact eigenvalues and eigenvectors:

```wolfram
exacteigsys=DEigensystem[{ℒ,ℬ},u[x],{x,0,Pi},3,Method->"Normalize"]
(* Output *)
{{1,4,9},{Sqrt[(2)/(π)] Sin[x],Sqrt[(2)/(π)] Sin[2 x],Sqrt[(2)/(π)] Sin[3 x]}}
```

```wolfram
exacteigsys/.{x-> 0.3}
(* Output *)
{{1,4,9},{0.23579101030035493,0.4505195118954414,0.6250044472531904}}
```

Find numerical eigenvalues and eigenvectors:

```wolfram
neigsys=NDEigensystem[{ℒ,ℬ},u[x],{x,0,Pi},3]/.{x-> 0.3}
(* Output *)
{{1.0000008444742254,4.000053838423359,9.000609362487232},{0.2357744171180486,0.45040373369760645,0.6247026399235258}}
```

Use [DEigenvalues](https://reference.wolfram.com/language/ref/DEigenvalues.html) to find the eigenvalues for a differential operator:

```wolfram
{ℒ,ℬ}={-Laplacian[u[x],{x}],DirichletCondition[u[x]==0,True]};
```

Find eigenvalues and eigenfunctions:

```wolfram
{vals,funs}=DEigensystem[{ℒ,ℬ},u[x],{x,0,Pi},3]
(* Output *)
{{1,4,9},{Sin[x],Sin[2 x],Sin[3 x]}}
```

```wolfram
vals
(* Output *)
{1,4,9}
```

Find eigenvalues only:

```wolfram
DEigenvalues[{ℒ,ℬ},u[x],{x,0,Pi},3]
(* Output *)
{1,4,9}
```

Use [DSolve](https://reference.wolfram.com/language/ref/DSolve.html) to solve an eigenvalue problem:

```wolfram
{ℒ,ℬ}={-Laplacian[u[x],{x}],DirichletCondition[u[x]==0,True]};
```

Find eigenvalues and eigenfunctions:

```wolfram
{vals,funs}=DEigensystem[{ℒ,ℬ},u[x],{x,0,Pi},3]
(* Output *)
{{1,4,9},{Sin[x],Sin[2 x],Sin[3 x]}}
```

Find the complete eigensystem:

```wolfram
DSolve[{u''[x]+λ u[x]==0,u[0]==0,u[π]==0},u[x],x]
(* Output *)
{{u[x]->{, {{1 Sin[x Sqrt[λ]], n∈Integers&&n>=1&&λ==n^2}, {0, True}}}}}
```

The eigenfunctions given by [DEigensystem](https://reference.wolfram.com/language/ref/DEigensystem.html) are orthogonal:

```wolfram
{ℒ,ℬ}={-Laplacian[u[x],{x}],DirichletCondition[u[x]==0,True]};
```

Find eigenvalues and eigenfunctions:

```wolfram
{vals,funs}=DEigensystem[{ℒ,ℬ},u[x],{x,0,Pi},4]
(* Output *)
{{1,4,9,16},{Sin[x],Sin[2 x],Sin[3 x],Sin[4 x]}}
```

Verify that the eigenfunctions are orthogonal:

```wolfram
Table[Integrate[funs[[i]]*funs[[j]],{x,0,Pi}],{i,1,4},{j,1,i-1}]//Flatten
(* Output *)
{0,0,0,0,0,0}
```

The system of eigenfunctions given by [DEigensystem](https://reference.wolfram.com/language/ref/DEigensystem.html) is not orthonormal by default:

```wolfram
{ℒ,ℬ}={-Laplacian[u[x],{x}],DirichletCondition[u[x]==0,True]};
```

Find eigenvalues and eigenfunctions:

```wolfram
{vals,funs}=DEigensystem[{ℒ,ℬ},u[x],{x,0,Pi},3]
(* Output *)
{{1,4,9},{Sin[x],Sin[2 x],Sin[3 x]}}
```

The eigenfunctions are not orthonormalized by default:

```wolfram
Table[Integrate[funs[[i]] funs[[j]],{x,0,Pi}],{i,1,3},{j,1,3}]//MatrixForm
(* Output *)
({{(π)/(2), 0, 0}, {0, (π)/(2), 0}, {0, 0, (π)/(2)}})
```

Use [Method](https://reference.wolfram.com/language/ref/Method.html)->"Normalize" to obtain an orthonormal system:

```wolfram
{vals,funs}=DEigensystem[{ℒ,ℬ},u[x],{x,0,Pi},3,Method->"Normalize"]
(* Output *)
{{1,4,9},{Sqrt[(2)/(π)] Sin[x],Sqrt[(2)/(π)] Sin[2 x],Sqrt[(2)/(π)] Sin[3 x]}}
```

```wolfram
Table[Integrate[funs[[i]]*funs[[j]],{x,0,Pi}],{i,1,3},{j,1,3}]//MatrixForm
(* Output *)
({{1, 0, 0}, {0, 1, 0}, {0, 0, 1}})
```

Apply [N](https://reference.wolfram.com/language/ref/N.html)[DEigensystem[...]] to invoke [NDEigensystem](https://reference.wolfram.com/language/ref/NDEigensystem.html) if symbolic evaluation fails:

```wolfram
DEigensystem[{-Laplacian[u[x,y],{x,y}]+E^(-x^3+x) u[x,y],DirichletCondition[u[x,y]==0,True]},u[x,y],{x,y}∈Rectangle[],2]
(* Output *)
DEigensystem[{ℯ^(x-x^3) u[x,y]-u^(0,2)[x,y]-u^(2,0)[x,y],DirichletCondition[u[x,y]==0,True]},u[x,y],{x,y}∈Rectangle[{0,0}],2]
```

```wolfram
N[%]
(* Output *)
{{21.127458154737347,50.66140983559559},{InterpolatingFunction[...][x,y],InterpolatingFunction[...][x,y]}}
```

### Possible Issues

Inhomogeneous Dirichlet conditions are replaced with homogeneous ones:

```wolfram
First[DEigensystem[{-Laplacian[u[x],{x}],DirichletCondition[u[x]==1,True]},u[x],{x,0,1},2]]
(* Output *)
DEigensystem
(* Output *)
{π^2,4 π^2}
```

The same result:

```wolfram
First[DEigensystem[{-Laplacian[u[x],{x}],DirichletCondition[u[x]==0,True]},u[x],{x,0,1},2]]
(* Output *)
{π^2,4 π^2}
```

Inhomogeneous Neumann values are replaced with homogeneous ones:

```wolfram
First[DEigensystem[{-Laplacian[u[x],{x}]+NeumannValue[1,True]},u[x],{x,0,1},2]]
(* Output *)
DEigensystem
(* Output *)
{0,π^2}
```

The same result:

```wolfram
First[DEigensystem[-Laplacian[u[x],{x}],u[x],{x,0,1},2]]
(* Output *)
{0,π^2}
```

## Tech Notes ▪Symbolic Solutions of PDEs

## Related Guides ▪Partial Differential Equations ▪Differential Operators ▪Differential Equations ▪Solvers over Regions

## History Introduced in 2015 (10.3)
