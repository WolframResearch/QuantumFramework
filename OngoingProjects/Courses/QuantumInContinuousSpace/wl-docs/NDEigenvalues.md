# NDEigenvalues | [SpanFromLeft]

> [NDEigenvalues](https://reference.wolfram.com/language/ref/NDEigenvalues.html)[*ℒ*[*u*[*x*,*y*,…]],*u*,{*x*,*y*,…}∈Ω,*n*]  — gives the `*n*` smallest magnitude eigenvalues for the linear differential operator `*ℒ*` over the region `Ω`.
> [NDEigenvalues](https://reference.wolfram.com/language/ref/NDEigenvalues.html)[{*ℒ*_1[*u*[*x*,*y*,…],*v*[*x*,*y*,…],…],*ℒ*_2[*u*[*x*,*y*,…],*v*[*x*,*y*,…],…],…},{*u*,*v*,…},{*x*,*y*,…}∈Ω,*n*] — gives eigenvalues for the coupled differential operators `{*op*_1,*op*_2,…}` over the region `Ω`.
> [NDEigenvalues](https://reference.wolfram.com/language/ref/NDEigenvalues.html)[*eqns*,{*u*,…},*t*,{*x*,*y*,…}∈Ω,*n*] — gives the eigenvalues in the spatial variables `{*x*,*y*,…}` for solutions `*u*,…` of the coupled time-dependent differential equations `*eqns*`.

## Details and Options

[NDEigenvalues](https://reference.wolfram.com/language/ref/NDEigenvalues.html), also known as an eigenmode solver, is a numerical eigen solver that finds eigenvalues of differential equations over regions.

[NDEigenvalues](https://reference.wolfram.com/language/ref/NDEigenvalues.html) gives a list `{λ_1,…,λ_*n*}` of the `*n*` smallest magnitude eigenvalues `λ_*i*`.

The equations `*eqns*` are specified as in [NDSolve](https://reference.wolfram.com/language/ref/NDSolve.html).

Eigenvalues are sorted in order of increasing absolute value.

Homogeneous [DirichletCondition](https://reference.wolfram.com/language/ref/DirichletCondition.html), [NeumannValue](https://reference.wolfram.com/language/ref/NeumannValue.html) or generalized Robin boundary conditions may be included.

[PeriodicBoundaryCondition](https://reference.wolfram.com/language/ref/PeriodicBoundaryCondition.html) may be included.

When no boundary condition is specified on part of the boundary `∂Ω`, then this is equivalent to specifying a Neumann 0 condition.

For a system of first-order time-dependent equations, the time derivatives [D](https://reference.wolfram.com/language/ref/D.html)[*u*[*t*,*x*,*y*,*…*],*t*],[D](https://reference.wolfram.com/language/ref/D.html)[*v*[*t*,*x*,*y*,*…*],*t*],… are effectively replaced with `*λ****u*[*x*,*y*,…],*λ****v*[*x*,*y*,…],…`.

Systems of time-dependent equations that are higher than first order are reduced to a coupled first-order system with intermediate variables `*u*_*t*=*u*^(*),**u*_t^**=…`, `*v*_*t*=*v*^(*),*v_t^**=…`, `…`. Only the functions `*u*`, `*v*`, `…` are returned.

[NDEigenvalues](https://reference.wolfram.com/language/ref/NDEigenvalues.html) accepts a [Method](https://reference.wolfram.com/language/ref/Method.html) option that may be used to control different stages of the solution. With [Method](https://reference.wolfram.com/language/ref/Method.html)->{*s*_1->*m*_1,*s*_2->*m*_2,…}, stage `*s*_*i*` is handled by method `*m*_*i*`. When stages are not given explicitly, [NDEigenvalues](https://reference.wolfram.com/language/ref/NDEigenvalues.html) tries to automatically determine what stage to apply a given method to.

Possible solution stages are:

"PDEDiscretization" | discretization of spatial operators.
"Eigensystem" | computation of the eigensystem from the discretized system.
"VectorNormalization" | normalization of the eigenvectors that are used to construct the eigenfunctions.

## Examples

### Basic Examples

Find the 4 smallest eigenvalues of the Laplacian operator on `[0,π]`:

```wolfram
NDEigenvalues[-Laplacian[u[x],{x}],u[x],{x,0,π},4]
(* Output *)
{4.080266017155255×10^-14,1.0000008444742063,4.000053838423353,9.000609362487243}
```

Compute the first 6 eigenvalues for a circular membrane with the edges clamped:

```wolfram
NDEigenvalues[{-Laplacian[u[x,y],{x,y}],DirichletCondition[u[x,y]==0,True]},u[x,y],{x,y}∈Disk[],6]
(* Output *)
{5.783226344335697,14.68266911817372,14.682674600316588,26.378393529310333,26.379103039158753,30.47873170940837}
```

Specify a Schrödinger operator with parameter $h$ and potential $V$:

```wolfram
h=1/2;V[x_]:=x^2
ℒ=-h^2*u''[x]+V[x]*u[x];
```

Find the 5 smallest eigenvalues:

```wolfram
NDEigenvalues[ℒ,u[x],{x,-3,3},5]
(* Output *)
{0.5000411901298615,1.5002826947466013,2.5009631168612256,3.5020737612401676,4.502060720755467}
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

Find the 4 smallest eigenvalues:

```wolfram
NDEigenvalues[{ℒ,ℬ},u[x],{x,0,π},4]
(* Output *)
{1.0000008444742223,4.000053838423355,9.000609362487232,16.00339383619117}
```

Specify a Laplacian operator:

```wolfram
ℒ=-Laplacian[u[x],{x}];
```

Specify homogeneous Neumann boundary conditions:

```wolfram
ℬ=NeumannValue[0,True];
```

Find the four smallest eigenvalues:

```wolfram
NDEigenvalues[ℒ+ℬ,u[x],{x,0,π},4]
(* Output *)
{4.080266017155255×10^-14,1.0000008444742063,4.000053838423353,9.000609362487243}
```

This is equivalent:

```wolfram
NDEigenvalues[ℒ,u[x],{x,0,π},4]
(* Output *)
{4.080266017155255×10^-14,1.0000008444742063,4.000053838423353,9.000609362487243}
```

Specify a transient equation with homogeneous Dirichlet boundary conditions:

```wolfram
{ℒ,ℬ}={D[u[t,x],t]==-Laplacian[u[t,x],{x}],DirichletCondition[u[t,x]==0,True]};
```

Find the 4 smallest eigenvalues:

```wolfram
NDEigenvalues[{ℒ,ℬ},u[t,x],t,{x,0,π},4]
(* Output *)
{1.0000008444742223,4.000053838423355,9.000609362487232,16.00339383619117}
```

Specify a wave equation with homogeneous Dirichlet boundary conditions:

```wolfram
{ℒ,ℬ}={D[u[t,x],t,t]==Laplacian[u[t,x],{x}],DirichletCondition[u[t,x]==0,True]};
```

Find the 4 smallest eigenvalues:

```wolfram
NDEigenvalues[{ℒ,ℬ},u[t,x],t,{x,0,π},4]
(* Output *)
{-1.1104114731174438×10^-16-1.0000004222370051 ⅈ,-1.1104114731174438×10^-16+1.0000004222370051 ⅈ,-9.162451867121918×10^-16-2.00001345956054 ⅈ,-9.162451867121918×10^-16+2.00001345956054 ⅈ}
```

Compute the eigenvalues of a generalized wave equation $\frac{\partial^{2}u}{\partial t^{2}}+\gamma \frac{\partial u}{\partial t} - c^{2}\frac{\partial^{2}u}{\partial x^{2}}+\gamma u=0$:

```wolfram
γ=1.3;c=1.1;
NDEigenvalues[D[u[t,x],{t,2}]+γ D[u[t,x],{t,1}]-c^2 D[u[t,x],{x,2}]+γ u[t,x]==0,u,t,{x,0,π},4]
(* Output *)
{-0.6500000000000006-0.93674969975975 ⅈ,-0.6500000000000006+0.93674969975975 ⅈ,-0.6499999999999995-1.4448186812931718 ⅈ,-0.6499999999999995+1.4448186812931718 ⅈ}
```

Compare to the exact solution of an equivalent first-order system of ordinary differential equations:

```wolfram
Transpose[Eigenvalues[({{0, 1}, {-ω^2c^2-γ, -γ}})]/.ω->{0,1}]
(* Output *)
{{-0.65-0.9367496997597597 ⅈ,-0.65+0.9367496997597597 ⅈ},{-0.65-1.444818327679989 ⅈ,-0.65+1.444818327679989 ⅈ}}
```

Specify a Liouville operator:

```wolfram
ℒ=-Inactive[Div][{{1-x^2}}.Inactive[Grad][u[x],{x}],{x}]
(* Output *)
-(({{1-x^2}}.u[x]))
```

Compute the 4 smallest eigenvalues:

```wolfram
NDEigenvalues[ℒ,u[x],{x,-1,1},4]
(* Output *)
{-1.2345680033821396×10^-15,2.000000000000004,5.999999999999998,12.00014396503262}
```

Compare the eigenvalues with the analytical eigenvalues:

```wolfram
(n+1)n/.n->{0,1,2,3}
(* Output *)
{0,2,6,12}
```

Write a function to compute parametric, complex-valued periodic eigenvalues of the Laplace operator:

```wolfram
fun[alpha_]:=NDEigenvalues[{-u''[x],PeriodicBoundaryCondition[Exp[alpha I] u[x],x==2 π,Function[x,x-2 π]]},u[x],{x,0,2 π},5];
```

Find the eigenvalues:

```wolfram
res=Transpose[Table[fun[α],{α,0,4π,π/10}]];
```

Visualize the eigenvalues over the range from 0 to 4$\pi$:

```wolfram
ListLinePlot[res,DataRange->{0,4π}]
```

*([Graphics])*

#### 2D

Specify a Laplacian operator:

```wolfram
ℒ=-Laplacian[u[x,y],{x,y}];
```

Find the 4 smallest eigenvalues:

```wolfram
NDEigenvalues[ℒ,u[x,y],{x,0,π},{y,0,π},4]
(* Output *)
{2.2478597284837017×10^-14,1.0000015382281866,1.0000017695290675,2.0000125935329556}
```

Specify a Laplacian operator:

```wolfram
ℒ=-Laplacian[u[x,y],{x,y}];
```

Find the 4 smallest eigenvalues of the operator in a unit disk:

```wolfram
NDEigenvalues[ℒ,u[x,y],{x,y}∈Disk[],4]
(* Output *)
{1.5295779233167×10^-14,3.3899670948303995,3.3899672003966015,9.328561151575382}
```

Specify a Laplacian operator with homogeneous Dirichlet boundary conditions:

```wolfram
{ℒ,ℬ}={-Laplacian[u[x,y],{x,y}],DirichletCondition[u[x,y]==0,True]};
```

Find the 9 smallest eigenvalues in a rectangle:

```wolfram
NDEigenvalues[{ℒ,ℬ},u[x,y],{x,0,π},{y,0,π},9]
(* Output *)
{2.0000128656842495,5.000200550012288,5.0002113917591755,8.000798754375737,10.001479784631819,10.001664582371912,13.003262795810858,13.003586436347206,17.007275353139207}
```

Specify a wave equation with homogeneous Dirichlet boundary conditions:

```wolfram
{ℒ,ℬ}={D[u[t,x,y],t,t]==-Laplacian[u[t,x,y],{x,y}],DirichletCondition[u[t,x,y]==0.,True]};
```

Find the 4 smallest eigenvalues in a disk:

```wolfram
NDEigenvalues[{ℒ,ℬ},u[t,x,y],t,{x,y}∈Disk[],4]
(* Output *)
{-2.404833953589233,2.4048339535892334,3.831797113388663,-3.8317971133886655}
```

Solve a partially constrained eigenvalue problem:

```wolfram
NDEigenvalues[{-Laplacian[u[x,y],{x,y}],DirichletCondition[u[x,y]==0.,x==0]},u[x,y],{x,y}∈Rectangle[{0,0},{1,1}],9]
(* Output *)
{2.4674012306269595,12.337013968499967,22.206704686385642,32.0763175985178,41.9463502339355,61.68565651704908,61.687048683700525,71.55666280624403,91.2998551885231}
```

### Options

#### Method

"Eigensystem"  (4)

Specify a method to use for finding the eigenvalues:

```wolfram
AbsoluteTiming[NDEigenvalues[{-Laplacian[u[x,y],{x,y}]},u,{x,0,1},{y,0,1},4,Method->{"Eigensystem"->"Direct"}]]
(* Output *)
{1.203082,{0.,9.869621049921076,9.869621616269733,19.73934710540828}}
```

In this case, the default method is faster:

```wolfram
AbsoluteTiming[NDEigenvalues[{-Laplacian[u[x,y],{x,y}]},u,{x,0,1},{y,0,1},4]]
(* Output *)
{0.087045,{-3.6652834904648225×10^-14,9.869621049920202,9.869621616268848,19.739347105408093}}
```

Arnoldi is used as the default method:

```wolfram
AbsoluteTiming[NDEigenvalues[{-Laplacian[u[x,y],{x,y}]},u,{x,0,1},{y,0,1},4,Method->{"Eigensystem"->"Arnoldi"}]]
(* Output *)
{0.082279,{-3.6652834904648225×10^-14,9.869621049920202,9.869621616268848,19.739347105408093}}
```

Specify a maximum number of iterations for the Arnoldi method:

```wolfram
AbsoluteTiming[NDEigenvalues[{-Laplacian[u[x,y],{x,y}]},u,{x,0,1},{y,0,1},4,Method->{"Eigensystem"->{"Arnoldi","MaxIterations"->100}}]]
(* Output *)
{0.105803,{-3.6652834904648225×10^-14,9.869621049920202,9.869621616268848,19.739347105408093}}
```

Find two eigenvalues of a Sturm-Liouville operator within the band of $(5,6)$ with the FEAST method for [Eigenvalues](https://reference.wolfram.com/language/ref/Eigenvalues.html):

```wolfram
vals=NDEigenvalues[{u''[x]+RealAbs[x] u[x],DirichletCondition[u[x]==0,True]},u[x],{x,-8,8},2,Method->{"SpatialDiscretization"->{"FiniteElement",{"MeshOptions"->{"MaxCellMeasure"->0.05}}},"Eigensystem"->{"FEAST","Interval"->{5,6}}}]
(* Output *)
{5.661892563705698,5.661892573969469}
```

According to the Sturm-Liouville theory, the eigenvalues must be distinct, but for this example they are close to degenerate:

```wolfram
Subtract@@vals
(* Output *)
-1.0263771166307833×10^-8
```

The interval end points $(5,6)$ are not included in the interval FEAST finds eigenvalues in; for more information, please refer to the [Eigenvalues](https://reference.wolfram.com/language/ref/Eigenvalues#2045383922.html) reference page.

The usage of the `"Shift"` option is [explained in an example below](https://reference.wolfram.com/language/ref/NDEigenvalues.html#42461913).

"PDEDiscretization"  (1)

Change the [MaxCellMeasure](https://reference.wolfram.com/language/ref/MaxCellMeasure.html) for the underlying computation:

```wolfram
NDEigenvalues[{-Laplacian[u[x],{x}]},u,{x,0,π},4,Method->{"PDEDiscretization"->{"FiniteElement",{"MeshOptions"->{"MaxCellMeasure"->0.5}}}}]
(* Output *)
{0.,1.0000557617465335,4.003462494088496,9.037620966219041}
```

The exact eigenvalues are `0,1,4,9,…`, so the eigenvalue error is:

```wolfram
err=Abs[% - {0,1,4,9}]
(* Output *)
{0.,0.00005576174653354471,0.003462494088496193,0.03762096621904121}
```

A finer mesh results in a decreased discretization error:

```wolfram
NDEigenvalues[-Laplacian[u[x],{x}],u[x],{x,0,π},4,Method->{"PDEDiscretization"->{"FiniteElement","MeshOptions"->{"MaxCellMeasure"->0.01}}}]
(* Output *)
{-3.143239121948122×10^-12,1.0000000000107128,4.0000000008770495,9.00000001001247}
```

```wolfram
err=Abs[%-{0,1,4,9}]
(* Output *)
{3.143239121948122×10^-12,1.071276400921306×10^-11,8.770495441012827×10^-10,1.00124708524163×10^-8}
```

"VectorNormalization"  (1)

Compute without any normalization of the computed eigenvalues:

```wolfram
NDEigenvalues[-Laplacian[u[x],{x}],u[x],{x,0,π},4,Method->{"VectorNormalization"->None}]
(* Output *)
{4.080266017155255×10^-14,1.0000008444742063,4.000053838423353,9.000609362487243}
```

The normalization does not have an effect on the eigenvalues:

```wolfram
NDEigenvalues[-Laplacian[u[x],{x}],u[x],{x,0,π},4]
(* Output *)
{4.080266017155255×10^-14,1.0000008444742063,4.000053838423353,9.000609362487243}
```

### Applications

#### Acoustics

Compute the acoustic eigenvalues and eigenfunctions for an approximation of a cross section through a [Mini](https://en.wikipedia.org/wiki/Mini). Import an image of the cross section:

```wolfram
img=Import["http://upload.wikimedia.org/wikipedia/commons/d/d3/Mini_cross_section.jpg"]
(* Output *)
![image](img/image_001.png)
```

Use the mask tool to create a boundary graphic:

```wolfram
boundary=[Graphics];
```

Discretize the graphic:

```wolfram
bdr=BoundaryDiscretizeGraphics[boundary]
```

*([Graphics])*

Compute 6 eigenvalues of the cross section:

```wolfram
NDEigenvalues[{-Laplacian[u[x,y],{x,y}]},u[x,y],{x,y}∈bdr,6]
(* Output *)
{1.1576383174794323×10^-20,1.270729393334167×10^-6,3.956142113696531×10^-6,6.1204073770419415×10^-6,9.976858517360195×10^-6,0.000011841218019956593}
```

#### Structural Mechanics

Specify a plane stress PDE:

```wolfram
ps=SolidMechanicsPDEComponent[{{u[x,y],v[x,y]},{x,y}},<|"YoungModulus"->10^3,"PoissonRatio"->33/100,"Thickness"->1|>];
```

Compute constraint eigenvalues:

```wolfram
NDEigenvalues[{ps,DirichletCondition[{u[x,y]==0.,v[x,y]==0.},x==0]}, {u[x,y],v[x,y]},{x,y}∈Rectangle[{0,0},{5,1}],9]
(* Output *)
{1.558451891907176,45.29890210583204,99.21385118911145,257.60383626643574,713.0087913127235,886.4698406606385,1446.6716761256528,2416.979263774888,2439.7851066438293}
```

#### Interval of Eigenvalues and Eigenfunctions

Find an eigenvalue in an interval:

```wolfram
NDEigenvalues[{-Laplacian[u[x,y],{x,y}],DirichletCondition[u[x,y]==0,True]},u,{x,y}∈ImplicitRegion[(x^2+y^2+2 y)^2<4 (x^2+y^2),{x,y}],1,Method->{"Eigensystem"->{"FEAST","Interval"->{400,405}},"PDEDiscretization"->{"FiniteElement","MeshOptions"->{"MaxCellMeasure"->0.001}}}]
(* Output *)
{400.79135867884787}
```

#### Quantum Mechanics

Specify a Schrödinger operator with parameter $h$ and potential $V$:

```wolfram
h=1/10;V[x_]:=x^2
ℒ=-h^2*u''[x]+V[x]*u[x];
```

Find the 10 smallest eigenvalues on a refined mesh:

```wolfram
NDEigenvalues[ℒ,u[x],{x,-3,3},10,Method->{"PDEDiscretization"->{"FiniteElement",{"MeshOptions"->{MaxCellMeasure->0.05}}}}]
(* Output *)
{0.10000016225988044,0.30000113385045535,0.5000040417329503,0.7000101643270643,0.9000207683332292,1.1000371088172076,1.3000604292934634,1.50009196180729,1.7001329270156675,1.9001845342669057}
```

### Properties & Relations

Specify a transient equation with homogeneous Dirichlet boundary conditions:

```wolfram
{ℒ,ℬ}={D[u[t,x],t]==-Laplacian[u[t,x],{x}],DirichletCondition[u[t,x]==0,True]};
```

Find the 4 smallest eigenvalues:

```wolfram
NDEigenvalues[{ℒ,ℬ},u[t,x],t,{x,0,π},4]
(* Output *)
{1.0000008444742172,4.000053838423344,9.000609362487225,16.003393836191144}
```

This is equivalent:

```wolfram
{ℒ,ℬ}={-Laplacian[u[x],{x}],DirichletCondition[u[x]==0,True]};
NDEigenvalues[{ℒ,ℬ},u[x],{x,0,π},4]
(* Output *)
{1.0000008444742172,4.000053838423344,9.000609362487225,16.003393836191144}
```

Compare an analytical solution to a higher-order time-dependent PDE.

Find the six smallest eigenvalues of a wave equation $u_{\text{t t}}=\mathit{u}_{\text{x x}}$ for $\mathit{x}$ between 0 and $\pi$:

```wolfram
NDEigenvalues[{D[u[t,x],{t,2}]==Laplacian[u[t,x],{x}]},{u[t,x]},t,{x,0,π},6]
(* Output *)
{1.1406479921142865×10^-17-1.122427229034216×10^-7 ⅈ,1.1406479921142865×10^-17+1.122427229034216×10^-7 ⅈ,9.714272158138526×10^-16-1.0000004222369994 ⅈ,9.714272158138526×10^-16+1.0000004222369994 ⅈ,1.4152545456455693×10^-15-2.000013459560544 ⅈ,1.4152545456455693×10^-15+2.000013459560544 ⅈ}
```

Compare the eigenvalues with the exact solution of $u''(t)=-\omega^{2}u(t)$ transformed into a system of two first-order equations:

```wolfram
Eigenvalues[({{0, 1}, {-ω^2, 0}})]
(* Output *)
{-ⅈ ω,ⅈ ω}
```

```wolfram
Transpose[%/.ω->Range[0,3]]
(* Output *)
{{0,0},{-ⅈ,ⅈ},{-2 ⅈ,2 ⅈ},{-3 ⅈ,3 ⅈ}}
```

Show the relation between higher-order time-dependent PDEs and systems of first-order PDEs.

Find the six smallest eigenvalues of a wave equation $u_{tt}=\mathit{u}_{xx}$ within 0 and $\pi$:

```wolfram
NDEigenvalues[{D[u[t,x],{t,2}]==Laplacian[u[t,x],{x}]},{u[t,x]},t,{x,0,π},6]
(* Output *)
{1.1406479921142865×10^-17-1.122427229034216×10^-7 ⅈ,1.1406479921142865×10^-17+1.122427229034216×10^-7 ⅈ,9.714272158138526×10^-16-1.0000004222369994 ⅈ,9.714272158138526×10^-16+1.0000004222369994 ⅈ,1.4152545456455693×10^-15-2.000013459560544 ⅈ,1.4152545456455693×10^-15+2.000013459560544 ⅈ}
```

Find the six smallest eigenvalues of a wave equation given as a system of first-order PDEs $(\begin{pmatrix}
u_{t} & = & v \\
v_{t} & = & u_{xx}
\end{pmatrix})$:

```wolfram
NDEigenvalues[{D[u[t,x],t]==v[t,x],D[v[t,x],t]==Laplacian[u[t,x],{x}]},{u[t,x],v[t,x]},t,{x,0,π},6]
(* Output *)
{1.1406479921142865×10^-17-1.122427229034216×10^-7 ⅈ,1.1406479921142865×10^-17+1.122427229034216×10^-7 ⅈ,9.714272158138526×10^-16-1.0000004222369994 ⅈ,9.714272158138526×10^-16+1.0000004222369994 ⅈ,1.4152545456455693×10^-15-2.000013459560544 ⅈ,1.4152545456455693×10^-15+2.000013459560544 ⅈ}
```

The eigenvalues for the second-order system and the system of first-order equations are the same:

```wolfram
%%-%
(* Output *)
{0.+0. ⅈ,0.+0. ⅈ,0.+0. ⅈ,0.+0. ⅈ,0.+0. ⅈ,0.+0. ⅈ}
```

### Possible Issues

The computed eigenvalues depend on the granularity of the discretization:

```wolfram
NDEigenvalues[-Laplacian[u[x],{x}],u[x],{x,0,π},4]
(* Output *)
{4.080266017155255×10^-14,1.0000008444742063,4.000053838423353,9.000609362487243}
```

The exact eigenvalues are `0,1,4,9,…`, so the eigenvalue error is:

```wolfram
err=Abs[% - {0,1,4,9}]
(* Output *)
{4.080266017155255×10^-14,8.444742063407062×10^-7,0.000053838423353269604,0.0006093624872427483}
```

A finer mesh results in a decreased discretization error:

```wolfram
NDEigenvalues[-Laplacian[u[x],{x}],u[x],{x,0,π},4,Method->{"SpatialDiscretization"->{"FiniteElement","MeshOptions"->{"MaxCellMeasure"->0.01}}}]
(* Output *)
{-3.0865446100242853×10^-12,1.0000000000113762,4.000000000877267,9.0000000100125}
```

```wolfram
err=Abs[% - {0,1,4,9}]
(* Output *)
{3.0865446100242853×10^-12,1.1376233288729054×10^-11,8.772671478141092×10^-10,1.001249927412573×10^-8}
```

The eigenvalues of the wave equation will be the square root of the angular frequencies:

```wolfram
NDEigenvalues[{D[u[t,x],{t,2}]==Laplacian[u[t,x],{x}]},{u[t,x]},t,{x,0,π},6]
(* Output *)
{1.1406479921142865×10^-17-1.122427229034216×10^-7 ⅈ,1.1406479921142865×10^-17+1.122427229034216×10^-7 ⅈ,9.714272158138526×10^-16-1.0000004222369994 ⅈ,9.714272158138526×10^-16+1.0000004222369994 ⅈ,1.4152545456455693×10^-15-2.000013459560544 ⅈ,1.4152545456455693×10^-15+2.000013459560544 ⅈ}
```

Compare to the exact solution of an equivalent first-order system of ordinary differential equations:

```wolfram
Transpose[Eigenvalues[({{0, 1}, {-ω^2, 0}})]/.ω->{0,1,2}]
(* Output *)
{{0,0},{-ⅈ,ⅈ},{-2 ⅈ,2 ⅈ}}
```

Eigenvalues with inhomogeneous Dirichlet conditions cannot be solved for:

```wolfram
NDEigenvalues[{-Laplacian[u[x],{x}],DirichletCondition[u[x]==1,True]},u[x],{x,0,1},2]
(* Output *)
NDEigenvalues
(* Output *)
NDEigenvalues[{-u^′′[x],DirichletCondition[u[x]==1,True]},u[x],{x,0,1},2]
```

Eigenvalues with homogeneous Dirichlet conditions can be solved for:

```wolfram
NDEigenvalues[{-Laplacian[u[x],{x}],DirichletCondition[u[x]==0,True]},u[x],{x,0,1},2]
(* Output *)
{9.869612735715425,39.47894896829713}
```

Eigenvalues with inhomogeneous Neumann values cannot be solved for:

```wolfram
NDEigenvalues[-Laplacian[u[x],{x}]+NeumannValue[1,x==1],u[x],{x,0,1},2]
(* Output *)
NDEigenvalues
(* Output *)
NDEigenvalues[NeumannValue[1,x==1]-u^′′[x],u[x],{x,0,1},2]
```

Eigenvalues with homogeneous Neumann values can be solved for:

```wolfram
NDEigenvalues[-Laplacian[u[x],{x}]+NeumannValue[0,x==1],u[x],{x,0,1},2]
(* Output *)
{6.365278674517519×10^-14,9.86961273571546}
```

The same result:

```wolfram
NDEigenvalues[-Laplacian[u[x],{x}],u[x],{x,0,1},2]
(* Output *)
{6.365278674517519×10^-14,9.86961273571546}
```

Eigenvalues with inhomogeneous generalized Neumann values cannot be solved for:

```wolfram
NDEigenvalues[-Laplacian[u[x],{x}]+NeumannValue[u[x]+1,x==1],u[x],{x,0,1},2]
(* Output *)
NDEigenvalues
(* Output *)
NDEigenvalues[NeumannValue[1+u[x],x==1]-u^′′[x],u[x],{x,0,1},2]
```

The operator and possible boundary conditions need to be stationary and linear:

```wolfram
NDEigenvalues[{D[u[t,x,y],t]==(D[u[t,x,y],{x,2}]+D[u[t,x,y],{y,2}])+(.01 t)(D[u[t,x,y],{x,2}]+D[u[t,x,y],{y,2}]),DirichletCondition[u[t,x,y]==0.,True]},u,t,{x,y}∈Disk[],10]
(* Output *)
NDEigenvalues
(* Output *)
NDEigenvalues[{u^(1,0,0)[t,x,y]==u^(0,0,2)[t,x,y]+u^(0,2,0)[t,x,y]+0.01 t (u^(0,0,2)[t,x,y]+u^(0,2,0)[t,x,y]),DirichletCondition[u[t,x,y]==0.,True]},u,t,{x,y}∈Disk[{0,0}],10]
```

Initial conditions will be set to zero and ignored:

```wolfram
γ=1.3;c=1.1;
f[x_]=D[0.125 Erf[(x-0.5)/0.125],x];
vInit[x_]=-c*D[f[x],x]-γ f[x]/2;
NDEigenvalues[{D[u[t,x],{t,2}]+γ D[u[t,x],{t,1}]-c^2 D[u[t,x],{x,2}]+γ u[t,x]==0,u[0,x]==f[x],
Derivative[1,0][u][0,x]==vInit[x],
DirichletCondition[u[t,x]==0,x==0],DirichletCondition[u[t,x]==0,x==1]},u,t,{x,0,1},4]
(* Output *)
NDEigenvalues
(* Output *)
NDEigenvalues
(* Output *)
{-0.6499999999999928-3.58046525052479 ⅈ,-0.6499999999999928+3.58046525052479 ⅈ,-0.6500000000000045-6.97474216381078 ⅈ,-0.6500000000000045+6.97474216381078 ⅈ}
```

The same result:

```wolfram
NDEigenvalues[{D[u[t,x],{t,2}]+γ D[u[t,x],{t,1}]-c^2 D[u[t,x],{x,2}]+γ u[t,x]==0,DirichletCondition[u[t,x]==0,x==0],DirichletCondition[u[t,x]==0,x==1]},u,t,{x,0,1},4]
(* Output *)
{-0.6499999999999928-3.58046525052479 ⅈ,-0.6499999999999928+3.58046525052479 ⅈ,-0.6500000000000045-6.97474216381078 ⅈ,-0.6500000000000045+6.97474216381078 ⅈ}
```

`NDEigenvalue` converts PDEs to time-dependent PDEs. This transformation is not unique and may lead to what seem to be unexpected results for coupled PDEs:

```wolfram
NDEigenvalues[{
-u[x]-Laplacian[u[x],{x}],
-v[x]-Laplacian[v[x],{x}]},{v,u},{x}∈Line[{{0},{1}}],1]
(* Output *)
Eigensystem
(* Output *)
{-0.027280306458415445}
```

Internally, the given equations are rewritten as a system of time-dependent PDEs. In the previous case from the given dependent variables `{*v*[*x*],*u*[*x*]}`, the following temporal system is generated: `{[D](https://reference.wolfram.com/language/ref/D.html)[*v*[*t*, *x*], *t*] == - *u*[*t*, *x*] - [Laplacian](https://reference.wolfram.com/language/ref/Laplacian.html)[*u*[*t*, *x*], {*x*}],[D](https://reference.wolfram.com/language/ref/D.html)[*u*[*t*, *x*], *t*] == -*v*[*t*, *x*] - [Laplacian](https://reference.wolfram.com/language/ref/Laplacian.html)[*v*[*t*, *x*], {*x*}]}`

To uniquely specify the system of equations, it is best to use the temporal description:

```wolfram
NDEigenvalues[{
D[u[t,x],t]==-u[t,x]-Laplacian[u[t,x],{x}],
D[v[t,x],t]==-v[t,x]-Laplacian[v[t,x],{x}]},{v,u},t,{x}∈Line[{{0},{1}}],1]
(* Output *)
{-0.9999999999998526}
```

More information on this topic can be found in [Finite Element Method Usage Tips](https://reference.wolfram.com/language/FEMDocumentation/tutorial/FiniteElementBestPractice.html#898113620).

In some cases, [NDEigenvalues](https://reference.wolfram.com/language/ref/NDEigenvalues.html) may return what seem to be unexpected results for coupled PDEs:

```wolfram
eqTest={
-u[x]-Derivative[2][u][x],
-d[x]-Derivative[2][d][x]};
NDEigenvalues[eqTest,{u[x],d[x]},{x}∈Line[{{0},{700}}],1]
(* Output *)
Eigensystem
(* Output *)
{0.031085077593515247}
```

One way to avoid this issue is to specify an ordering of the dependent variables via the `"InterpolationOrder"` option:

```wolfram
NDEigenvalues[eqTest,{u[x],d[x]},{x}∈Line[{{0},{700}}],1,Method->{"PDEDiscretization"->{"FiniteElement","InterpolationOrder"->{u->2,d->2}}}]
(* Output *)
{-0.9510204081632653}
```

Alternatively, the `"Direct"` method can be used:

```wolfram
NDEigenvalues[eqTest,{u[x],d[x]},{x}∈Line[{{0},{700}}],1,Method->"Direct"]
(* Output *)
{-0.9510204081632647}
```

More information on this topic can be found in [Finite Element Method Usage Tips](https://reference.wolfram.com/language/FEMDocumentation/tutorial/FiniteElementBestPractice.html#898113620).

[NDEigenvalues](https://reference.wolfram.com/language/ref/NDEigenvalues.html) finds the $n$ smallest magnitude eigenvalues for a given linear differential operator. Particularly in cases where one is interested in finding the most negative eigenvalues, e.g. many quantum mechanical problems, the $n$ most negative eigenvalues might not correspond with the $n$ eigenvalues [NDEigenvalues](https://reference.wolfram.com/language/ref/NDEigenvalues.html) returns by default. To see this, consider the following example.

Take, for instance, the dimensionless radial Schrödinger equation for the hydrogen atom, where the energy units are Rydbergs and length is measured in Bohr radii.

Define the radial Schrödinger equation:

```wolfram
H[l_]:=-ψ''[r]+((l(l+1))/(r^2)-(2)/(r))ψ[r]
```

Solve the eigenvalue problem for $l=0$. A refined mesh is used for a good approximation quality:

```wolfram
en=NDEigenvalues[{H[0],DirichletCondition[ψ[r]==0,True]},ψ[r],{r,0,100},4,Method->{"SpatialDiscretization"->{"FiniteElement",{"MeshOptions"->{"MaxCellMeasure"->0.001}}}}];
```

Look at the eigenvalues obtained in electronvolts:

```wolfram
en UnitConvert[,"Electronvolts"]
(* Output *)
{0.04506395633158937,-0.12688225963282415,0.25216093851455273,-0.2611306700368917}
```

These correspond to the eigenvalues closest to $0$. But they are not in the desired order. For instance, take a look at the analytical eigenvalues for this problem.

The analytical energies for this equation follow the following relation:

```wolfram
enHydrogen[n_]:=-(1)/(n^2)
UnitConvert[Table[N[enHydrogen[n],5],{n,4}],"Electronvolts"]
(* Output *)
{-13.6056931229942343775,-3.4014232807485585944,-1.5117436803326927088,-0.8503558201871396486}
```

If you know in advance a lower bound for the eigenvalues, say $\lambda_{0}$, you can redefine the eigenvalue problem $\mathcal{H} \psi=E \psi$ as $(\mathcal{H}-\lambda_{0})\psi=(E-\lambda_{0})\psi$. In other words, a new eigenvalue problem: $\mathcal{L} \psi=\mathcal{E} \psi$, where $\mathcal{L}=\mathcal{H}-\lambda_{0}$ and $\mathcal{E}=E-\lambda_{0}$. Given that $\lambda_{0}$ is a negative lower bound for $E$, the shift ensures that $\mathcal{E}$ will only assume positive values. Therefore, the $n$ smallest magnitude values of $\mathcal{E}$ will correspond to the $n$ most negative values of $E$. That way, one can calculate for $\mathcal{E}$ and then get $E$. The easiest way to do this is to use the `"Shift"` option.

Set a lower bound:

```wolfram
λ0=QuantityMagnitude[-2 ]
(* Output *)
-2
```

Define a shift and solve the eigenvalue problem:

```wolfram
en=NDEigenvalues[{H[0],DirichletCondition[ψ[r]==0,True]},ψ[r],{r,0,100},4,Method->{"Eigensystem"->{"Arnoldi","Shift"->λ0},"SpatialDiscretization"->{"FiniteElement",{"MeshOptions"->{"MaxCellMeasure"->0.001}}}}];
```

Look at the eigenvalues obtained in electronvolts:

```wolfram
en UnitConvert[,"Electronvolts"]
(* Output *)
{-0.8503558211282636,-1.5117436835664173,-3.401423284831847,-13.605693122766406}
```

These correspond to the most negative eigenvalues, but in reverse order, and then it is necessary to sort the eigenstates.

Sort the eigenstates with [Sort](https://reference.wolfram.com/language/ref/Sort.html):

```wolfram
en=Sort[en];
```

Take a look at the difference between the analytical eigenvalues and the obtained results:

```wolfram
Table[QuantityMagnitude[enHydrogen[n]],{n,4}]-en
(* Output *)
{-1.6745049791211386×10^-11,3.001161541504871×10^-10,2.376743801768555×10^-10,6.917133532624575×10^-11}
```

## Tech Notes ▪Advanced Numerical Differential Equation Solving in the Wolfram Language ▪Numerical Mathematics in the Wolfram Language ▪PDE Models Overview

## Related Guides ▪Partial Differential Equations ▪Differential Equations ▪Solvers over Regions ▪Differential Operators

## History Introduced in 2015 (10.2)
