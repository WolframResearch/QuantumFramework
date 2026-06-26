# DEigenvalues | [SpanFromLeft]

> [DEigenvalues](https://reference.wolfram.com/language/ref/DEigenvalues.html)[*ℒ*[*u*[*x*,*y*,…]],*u*,{*x*,*y*,…}∈Ω,*n*]  — gives the `*n*` smallest magnitude eigenvalues for the linear differential operator `*ℒ*` over the region `Ω`.
> [DEigenvalues](https://reference.wolfram.com/language/ref/DEigenvalues.html)[*eqns*,*u*,*t*,{*x*,*y*,…}∈Ω,*n*] — gives the eigenvalues for solutions `*u*` of the time-dependent differential equations `*eqns*`.

## Details and Options

[DEigenvalues](https://reference.wolfram.com/language/ref/DEigenvalues.html) can compute eigenvalues for ordinary and partial differential operators with given boundary conditions.

[DEigenvalues](https://reference.wolfram.com/language/ref/DEigenvalues.html) gives a list `{λ_1,…,λ_*n*}` of the `*n*` smallest magnitude eigenvalues `λ_*i*`.

An eigenvalue and eigenfunction pair `{λ_*i*,*u*_*i*}` for the differential operator `*ℒ*` satisfy `*ℒ*[*u*_*i*[*x*,*y*,…]]==λ_*i* *u*_*i*[*x*,*y*,…]`.

Homogeneous [DirichletCondition](https://reference.wolfram.com/language/ref/DirichletCondition.html) or [NeumannValue](https://reference.wolfram.com/language/ref/NeumannValue.html) boundary conditions may be included. Inhomogeneous boundary conditions will be replaced with corresponding homogeneous boundary conditions.

When no boundary condition is specified on the boundary `∂Ω`, then this is equivalent to specifying a Neumann 0 condition.

The equations `*eqns*` are specified as in [DSolve](https://reference.wolfram.com/language/ref/DSolve.html).

[N](https://reference.wolfram.com/language/ref/N.html)[DEigenvalues[…]] calls [NDEigenvalues](https://reference.wolfram.com/language/ref/NDEigenvalues.html) for eigenvalues that cannot be computed symbolically.

The [Assumptions](https://reference.wolfram.com/language/ref/Assumptions.html) option can be used to specify assumptions on parameters.

## Examples

### Basic Examples

Find the 4 smallest eigenvalues of the Laplacian operator on `[0,π]`:

```wolfram
DEigenvalues[{-Laplacian[u[x],{x}],DirichletCondition[u[x]==0,True]},u[x],{x,0,π},4]
(* Output *)
{1,4,9,16}
```

Compute the first 6 eigenvalues for a circular membrane with the edges clamped:

```wolfram
DEigenvalues[{-Laplacian[u[x,y],{x,y}],DirichletCondition[u[x,y]==0,True]},u[x,y],{x,y}∈Disk[],6]//
(* Output *)
{(0)^(2),(1)^(2),(1)^(2),(2)^(2),(2)^(2),(0)^(2)}
```

```wolfram
N[%]
(* Output *)
{5.783185962946785,14.681970642123892,14.681970642123892,26.374616427163392,26.374616427163392,30.471262343662094}
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

Find the 5 smallest eigenvalues in an interval:

```wolfram
DEigenvalues[{ℒ,ℬ},u[x],{x,0,π},5]
(* Output *)
{1,4,9,16,25}
```

Specify a Laplacian operator:

```wolfram
ℒ=-Laplacian[u[x],{x}];
```

Specify homogeneous Neumann boundary conditions:

```wolfram
ℬ=NeumannValue[0,True];
```

Find the 5 smallest eigenvalues in an interval:

```wolfram
DEigenvalues[ℒ+ℬ,u[x],{x,0,π},5]
(* Output *)
{0,1,4,9,16}
```

This is equivalent:

```wolfram
DEigenvalues[ℒ,u[x],{x,0,π},5]
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

Find the 5 smallest eigenvalues in an interval:

```wolfram
DEigenvalues[{ℒ+ℬ2,ℬ1},u[x],{x,0,π},5]
(* Output *)
{(1)/(4),(9)/(4),(25)/(4),(49)/(4),(81)/(4)}
```

Find the $n$$^{th}$ eigenvalue:

```wolfram
FindSequenceFunction[%,n]//Factor
(* Output *)
(1)/(4) (-1+2 n)^2
```

Find symbolic expressions for the eigenvalues of a Laplace operator:

```wolfram
DEigenvalues[{-Laplacian[u[x],{x}],DirichletCondition[u[x]==0,True]},u[x],{x,a,b},3]
(* Output *)
{(π^2)/((-a+b)^2),(4 π^2)/((-a+b)^2),(9 π^2)/((-a+b)^2)}
```

Specify a Laplacian operator:

```wolfram
ℒ=-Laplacian[u[x],{x}];
```

Specify a homogeneous Dirichlet boundary condition:

```wolfram
ℬ1=DirichletCondition[u[x]==0,x==π];
```

Specify a homogeneous nonzero Neumann boundary condition:

```wolfram
ℬ2=NeumannValue[-u[x]/3,x==0];
```

Find the 5 smallest eigenvalues in an interval:

```wolfram
vals=DEigenvalues[{ℒ+ℬ2,ℬ1},u[x],{x,0,π},5];
```

```wolfram
vals[[1]]//
(* Output *)
Root
```

Specify an Airy operator:

```wolfram
ℒ=-Laplacian[u[x],{x}]+x u[x];
```

Specify a homogeneous Dirichlet boundary condition:

```wolfram
ℬ=DirichletCondition[u[x]==0,True];
```

Find the 3 smallest eigenvalues in an interval:

```wolfram
vals=DEigenvalues[{ℒ,ℬ},u[x],{x,0,1},3];
```

The eigenvalues are roots of a transcendental equation:

```wolfram
vals//
(* Output *)
{Root,Root,Root}
```

Numerical approximations for the eigenvalues:

```wolfram
N[vals,20]
(* Output *)
{10.36850716183633712655089182337070031722,39.97874478988335432509624230367350900171,89.32663454247874607960470633957266540657}
```

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
vals=DEigenvalues[{ℒ+ℬ},u[x],{x,0,1},5];
```

The eigenvalues are roots of a transcendental equation:

```wolfram
vals[[1]]//
(* Output *)
Root
```

Specify a heat equation with homogeneous Dirichlet boundary conditions:

```wolfram
{ℒ,ℬ}={D[u[t,x],t]==Laplacian[u[t,x],{x}],DirichletCondition[u[t,x]==0,True]};
```

Find the 4 smallest eigenvalues:

```wolfram
DEigenvalues[{ℒ,ℬ},u[t,x],t,{x,0,π},4]
(* Output *)
{-1,-4,-9,-16}
```

#### 2D

Specify a Laplacian operator with homogeneous Dirichlet boundary conditions:

```wolfram
{ℒ,ℬ}={-Laplacian[u[x,y],{x,y}],DirichletCondition[u[x,y]==0,True]};
```

Find the 9 smallest eigenvalues in a rectangle:

```wolfram
DEigenvalues[{ℒ,ℬ},u[x,y],{x,0,π},{y,0,π},9]
(* Output *)
{2,5,5,8,10,10,13,13,17}
```

Specify a Laplacian operator with homogeneous Neumann boundary conditions:

```wolfram
ℒ=-Laplacian[u[x,y],{x,y}]+NeumannValue[0,True];
```

Find the 4 smallest eigenvalues in a rectangle:

```wolfram
DEigenvalues[ℒ,u[x,y],{x,0,π},{y,0,π},4]
(* Output *)
{0,1,1,2}
```

Specify a Laplacian operator:

```wolfram
ℒ=-Laplacian[u[x,y],{x,y}];
```

Specify homogeneous Dirichlet boundary conditions:

```wolfram
ℬ=DirichletCondition[u[x,y]==0,True];
```

Find the 4 smallest eigenvalues of the operator in a unit disk:

```wolfram
DEigenvalues[{ℒ,ℬ},u[x,y],{x,y}∈Disk[],4]//
(* Output *)
{(0)^(2),(1)^(2),(1)^(2),(2)^(2)}
```

Specify a Laplacian operator:

```wolfram
ℒ=-Laplacian[u[x,y],{x,y}];
```

Specify homogeneous Dirichlet boundary conditions:

```wolfram
ℬ=DirichletCondition[u[x,y]==0,True];
```

Find the 6 smallest eigenvalues of the operator in a triangle:

```wolfram
DEigenvalues[{ℒ,ℬ},u[x,y],{x,y}∈Triangle[],6]
(* Output *)
{5 π^2,10 π^2,13 π^2,17 π^2,20 π^2,25 π^2}
```

Specify a Laplacian operator:

```wolfram
ℒ=-Laplacian[u[x,y],{x,y}];
```

Specify homogeneous Dirichlet boundary conditions:

```wolfram
ℬ=DirichletCondition[u[x,y]==0,True];
```

Find the 4 smallest eigenvalues in a sector of a disk:

```wolfram
DEigenvalues[{ℒ,ℬ},u[x,y],{x,y}∈Disk[{0,0},1,{0,Pi/5}],4]//
(* Output *)
{(5)^(2),(5)^(2),(10)^(2),(5)^(2)}
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

Find the 7 smallest eigenvalues in a cuboid:

```wolfram
DEigenvalues[{ℒ,ℬ},u[x,y,z],{x,y,z}∈Cuboid[{2,1,1},{4,2,3}],7]
(* Output *)
{(3 π^2)/(2),(9 π^2)/(4),(9 π^2)/(4),3 π^2,(7 π^2)/(2),(7 π^2)/(2),(17 π^2)/(4)}
```

Specify a Laplacian operator:

```wolfram
ℒ=-Laplacian[u[x,y,z],{x,y,z}];
```

Specify homogeneous Dirichlet boundary conditions:

```wolfram
ℬ=DirichletCondition[u[x,y,z]==0,True];
```

Find the 5 smallest eigenvalues in a cylinder:

```wolfram
DEigenvalues[{ℒ,ℬ},u[x,y,z],{x,y,z}∈Cylinder[{{0,0,0},{0,0,5}},2],5]//
(* Output *)
{(π^(2))/(25)+((0)^(2))/(4),(4 π^(2))/(25)+((0)^(2))/(4),(π^(2))/(25)+((1)^(2))/(4),(π^(2))/(25)+((1)^(2))/(4),(9 π^(2))/(25)+((0)^(2))/(4)}
```

Specify a Laplacian operator:

```wolfram
ℒ=-Laplacian[u[x,y,z],{x,y,z}];
```

Specify homogeneous Dirichlet boundary conditions:

```wolfram
ℬ=DirichletCondition[u[x,y,z]==0,True];
```

Find the 7 smallest eigenvalues in a ball:

```wolfram
DEigenvalues[{ℒ,ℬ},u[x,y,z],{x,y,z}∈Ball[{0,0,0},2],7]//
(* Output *)
{(1)/(4) ((1)/(2))^(2),(1)/(4) ((3)/(2))^(2),(1)/(4) ((3)/(2))^(2),(1)/(4) ((3)/(2))^(2),(1)/(4) ((5)/(2))^(2),(1)/(4) ((5)/(2))^(2),(1)/(4) ((5)/(2))^(2)}
```

Specify a Laplacian operator:

```wolfram
ℒ=-Laplacian[u[x,y,z],{x,y,z}];
```

Specify homogeneous Dirichlet boundary conditions:

```wolfram
ℬ=DirichletCondition[u[x,y,z]==0,True];
```

Find the 7 smallest eigenvalues in a prism:

```wolfram
DEigenvalues[{ℒ,ℬ},u[x,y,z],{x,y,z}∈Prism[{{0,0,0},{1,0,0},{0,1,0},{0,0,2},{1,0,2},{0,1,2}}],7]
(* Output *)
{(21 π^2)/(4),6 π^2,(29 π^2)/(4),9 π^2,(41 π^2)/(4),11 π^2,(45 π^2)/(4)}
```

### Properties & Relations

Use [NDEigenvalues](https://reference.wolfram.com/language/ref/NDEigenvalues.html) to find numerical eigenvalues and eigenvectors:

```wolfram
{ℒ,ℬ}={-Laplacian[u[x],{x}],DirichletCondition[u[x]==0,True]};
```

Exact eigenvalues:

```wolfram
exacteigvals=DEigenvalues[{ℒ,ℬ},u[x],{x,0,Pi},3]
(* Output *)
{1,4,9}
```

Numerical eigenvalues:

```wolfram
numeigvals=NDEigenvalues[{ℒ,ℬ},u[x],{x,0,Pi},3]
(* Output *)
{1.0000008444742254,4.000053838423359,9.000609362487232}
```

Use [DEigensystem](https://reference.wolfram.com/language/ref/DEigensystem.html) to find the eigensystem for a differential operator:

```wolfram
{ℒ,ℬ}={-Laplacian[u[x],{x}],DirichletCondition[u[x]==0,True]};
```

Find the eigenvalues:

```wolfram
DEigenvalues[{ℒ,ℬ},u[x],{x,0,Pi},3]
(* Output *)
{1,4,9}
```

Find the eigensystem:

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

Apply [N](https://reference.wolfram.com/language/ref/N.html)[DEigenvalues[…]] to invoke [NDEigenvalues](https://reference.wolfram.com/language/ref/NDEigenvalues.html) if symbolic evaluation fails:

```wolfram
DEigenvalues[{-Laplacian[u[x,y],{x,y}]+E^(-x^3+x) u[x,y],DirichletCondition[u[x,y]==0,True]},u[x,y],{x,y}∈Rectangle[],4]
(* Output *)
DEigenvalues[{ℯ^(x-x^3) u[x,y]-u^(0,2)[x,y]-u^(2,0)[x,y],DirichletCondition[u[x,y]==0,True]},u[x,y],{x,y}∈Rectangle[{0,0}],4]
```

```wolfram
N[%]
(* Output *)
{21.127458154737383,50.66140983559569,50.73679491487303,80.27075462597699}
```

### Possible Issues

Inhomogeneous Dirichlet conditions are replaced with homogeneous ones:

```wolfram
DEigenvalues[{-Laplacian[u[x],{x}],DirichletCondition[u[x]==1,True]},u[x],{x,0,1},2]
(* Output *)
DEigenvalues
(* Output *)
{π^2,4 π^2}
```

The same result:

```wolfram
DEigenvalues[{-Laplacian[u[x],{x}],DirichletCondition[u[x]==0,True]},u[x],{x,0,1},2]
(* Output *)
{π^2,4 π^2}
```

Inhomogeneous Neumann values are replaced with homogeneous ones:

```wolfram
DEigenvalues[{-Laplacian[u[x],{x}]+NeumannValue[1,True]},u[x],{x,0,1},2]
(* Output *)
DEigensystem
(* Output *)
{0,π^2}
```

The same result:

```wolfram
DEigenvalues[-Laplacian[u[x],{x}],u[x],{x,0,1},2]
(* Output *)
{0,π^2}
```

## Tech Notes ▪Symbolic Solutions of PDEs

## Related Guides ▪Partial Differential Equations ▪Differential Equations ▪Differential Operators ▪Solvers over Regions

## History Introduced in 2015 (10.3)
