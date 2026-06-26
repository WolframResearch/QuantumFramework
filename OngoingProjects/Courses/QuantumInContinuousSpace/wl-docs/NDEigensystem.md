# NDEigensystem | [SpanFromLeft]

> [NDEigensystem](https://reference.wolfram.com/language/ref/NDEigensystem.html)[*ℒ*[*u*[*x*,*y*,…]],*u*,{*x*,*y*,…}∈Ω,*n*] — gives the `*n*` smallest magnitude eigenvalues and eigenfunctions for the linear differential operator `*ℒ*` over the region `Ω`.
> [NDEigensystem](https://reference.wolfram.com/language/ref/NDEigensystem.html)[{**ℒ*_1*[*u*[*x*,*y*,…],*v*[*x*,*y*,…],…],**ℒ*_2*[*u*[*x*,*y*,…],*v*[*x*,*y*,…],…],…},{*u*,*v*,…},{*x*,*y*,…}∈Ω,*n*] — gives eigenvalues and eigenfunctions for the coupled differential operators `{*op*_1,*op*_2,…}` over the region `Ω`.
> [NDEigensystem](https://reference.wolfram.com/language/ref/NDEigensystem.html)[*eqns*,{*u*,…},*t*,{*x*,*y*,…}∈Ω,*n*] — gives the eigenvalues and eigenfunctions in the spatial variables `{*x*,*y*,…}` for solutions `*u*,…` of the coupled time-dependent differential equations `*eqns*`.

## Details and Options

[NDEigensystem](https://reference.wolfram.com/language/ref/NDEigensystem.html), also known as an eigenmode solver, is a numerical eigen solver that finds eigenvectors and eigenvalues of differential equations over regions.

[NDEigensystem](https://reference.wolfram.com/language/ref/NDEigensystem.html) gives lists `{{λ_1,…,λ_*n*},{*u*_1,…,*u*_*n*}}` of eigenvalues `λ_*i*` and eigenfunctions `*u*_*i*` or `{{λ_1,…,λ_*n*},{{*u*_1,*v*_1,…},…,{*u*_*n*,*v*_*n*,…}}}` in case of coupled systems.

The equations `*eqns*` are specified as in [NDSolve](https://reference.wolfram.com/language/ref/NDSolve.html).

An eigenvalue and eigenfunction pair `{λ_*i*,*u*_*i*}` for the differential operator `*ℒ*` satisfy `*ℒ*[*u*[*x*,*y*,…]]==λ_*i* *u*_*i*[*x*,*y*,…]`.

An eigenvalue and eigenfunctions pair `{λ_*i*,{*u*_*i*,*v*_*i*,…}}` for coupled differential operators satisfy:

``*ℒ*_1[*u*_*i*[*x*,*y*,…],*v*_*i*[*x*,*y*,…],…]==λ_*i* *u*_*i*[*x*,*y*,…]``
``*ℒ*_2[*u*_*i*[*x*,*y*,…],*v*_*i*[*x*,*y*,…],…]==λ_*i* *v*_*i*[*x*,*y*,…]``
`⋮`

Eigenvalues are sorted in order of increasing absolute value.

With the default normalization, the eigenfunctions `*u*_*i*` computed by [NDEigensystem](https://reference.wolfram.com/language/ref/NDEigensystem.html)[*ℒ*[*u*[*x*,*y*,…]],*u*,{*x*,*y*,…}∈Ω,*n*] approximately satisfy $\int_{\Omega}u_{i}\mathit{u_{i}}^{*}=1$.

With the default normalization, the eigenfunctions `{*u*_*i*,*v*_*i*,…}` for coupled differential operators approximately satisfy $\int_{\Omega}\{u_{i},v_{i},\ldots \}.\{\mathit{u_{i}}^{*},\mathit{v_{i}}^{*},\ldots \}=1$.

Homogeneous [DirichletCondition](https://reference.wolfram.com/language/ref/DirichletCondition.html), [NeumannValue](https://reference.wolfram.com/language/ref/NeumannValue.html) or generalized Robin boundary conditions may be included.

[PeriodicBoundaryCondition](https://reference.wolfram.com/language/ref/PeriodicBoundaryCondition.html) may be included.

When no boundary condition is specified on part of the boundary `∂Ω`, then this is equivalent to specifying a Neumann 0 condition.

For a system of first-order time-dependent equations, the time derivatives [D](https://reference.wolfram.com/language/ref/D.html)[*u*[*t*,*x*,*y*,…],*t*], [D](https://reference.wolfram.com/language/ref/D.html)[*v*[*t*,*x*,*y*,…],*t*],… are effectively replaced with `λ *u*[*x*,*y*,…], λ *v*[*x*,*y*,…],…`.

Systems of time-dependent equations that are higher than first order are reduced to a coupled first-order system with intermediate variables `*u*_*t*=*u*^(*),**u*_t^**=…`, `*v*_*t*=*v*^(*),*v_t^**=…`, `…` Only the functions `*u*`, `*v*`, `…` are returned.

[NDEigensystem](https://reference.wolfram.com/language/ref/NDEigensystem.html) accepts a [Method](https://reference.wolfram.com/language/ref/Method.html) option that may be used to control different stages of the solution. With [Method](https://reference.wolfram.com/language/ref/Method.html)->{*s*_1->*m*_1,*s*_2->*m*_2,…}, stage `*s*_*i*` is handled by method `*m*_*i*`. When stages are not given explicitly, [NDEigensystem](https://reference.wolfram.com/language/ref/NDEigensystem.html) tries to automatically determine what stage to apply a given method to.

Possible solution stages are:

"PDEDiscretization" | discretization of spatial operators
"Eigensystem" | computation of the eigensystem from the discretized system
"Interpolation" | creation of interpolating functions
"VectorNormalization" | normalization of the eigenvectors that are used to construct the eigenfunctions

## Examples

### Basic Examples

Find the 4 smallest eigenvalues and eigenfunctions of the Laplacian operator on `[0,π]`:

```wolfram
NDEigensystem[-Laplacian[u[x],{x}],u[x],{x,0,π},4]
(* Output *)
{{1.9718986802676216×10^-14,1.0000008444741935,4.000053838423326,9.000609362487214},{InterpolatingFunction[...][x],InterpolatingFunction[...][x],InterpolatingFunction[...][x],InterpolatingFunction[...][x]}}
```

Visualize the eigenfunctions:

```wolfram
Plot[Evaluate[%[[2]]],{x,0,π}]
```

*([Graphics])*

Compute the first 6 eigenfunctions for a circular membrane with the edges clamped:

```wolfram
{vals,funs}=NDEigensystem[{-Laplacian[u[x,y],{x,y}],DirichletCondition[u[x,y]==0,True]},u[x,y],{x,y}∈Disk[],6];
```

```wolfram
vals
(* Output *)
{5.783226344335697,14.68266911817372,14.682674600316588,26.378393529310333,26.379103039158753,30.47873170940837}
```

Visualize the eigenfunctions:

```wolfram
Table[Plot3D[funs[[i]],{x,y}∈Disk[],PlotRange->All,PlotLabel->vals[[i]],PlotTheme->"Minimal"],{i,Length[vals]}]
(* Output *)
![image](img/image_001.png)
```

Specify a Schrödinger operator with parameter $h$ and potential $V$:

```wolfram
h=1/2;V[x_]:=x^2
ℒ=-h^2*u''[x]+V[x]*u[x];
```

Find the 5 smallest eigenvalues and eigenfunctions:

```wolfram
{vals,funs}=NDEigensystem[ℒ,u[x],{x,-3,3},5];
```

```wolfram
vals
(* Output *)
{0.5000411901298615,1.5002826947466013,2.5009631168612256,3.5020737612401676,4.502060720755467}
```

Visualize the eigenfunctions scaled by $h$ and offset by the respective eigenvalue:

```wolfram
Show[Plot[Evaluate[h*funs+vals],{x,-3,3}],
Plot[V[x],{x,-3,3}],PlotRange->{{-3,3},{0,5}},AxesOrigin->{-3,0}]
(* Output *)
![image](img/image_003.png)
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

Find the 4 smallest eigenvalues and eigenfunctions:

```wolfram
{vals,funs}=NDEigensystem[{ℒ,ℬ},u[x],{x,0,π},4];
```

```wolfram
vals
(* Output *)
{1.0000008444742223,4.000053838423355,9.000609362487232,16.00339383619117}
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

Find the 4 smallest eigenvalues and eigenfunctions:

```wolfram
{vals,funs}=NDEigensystem[ℒ+ℬ,u[x],{x,0,π},4];
```

```wolfram
vals
(* Output *)
{4.080266017155255×10^-14,1.0000008444742063,4.000053838423353,9.000609362487243}
```

This is equivalent:

```wolfram
{vals,funs}=NDEigensystem[ℒ,u[x],{x,0,π},4];
```

```wolfram
vals
(* Output *)
{4.080266017155255×10^-14,1.0000008444742063,4.000053838423353,9.000609362487243}
```

Specify a transient equation with homogeneous Dirichlet boundary conditions:

```wolfram
{ℒ,ℬ}={D[u[t,x],t]==-Laplacian[u[t,x],{x}],DirichletCondition[u[t,x]==0,True]};
```

Find the 4 smallest eigenvalues and eigenfunctions:

```wolfram
{vals,funs}=NDEigensystem[{ℒ,ℬ},u[t,x],t,{x,0,π},4];
```

```wolfram
vals
(* Output *)
{1.0000008444742243,4.000053838423356,9.000609362487232,16.003393836191155}
```

This is equivalent:

```wolfram
{ℒ,ℬ}={-Laplacian[u[x],{x}],DirichletCondition[u[x]==0,True]};
{vals,funs}=NDEigensystem[{ℒ,ℬ},u[x],{x,0,π},4];
```

```wolfram
vals
(* Output *)
{1.0000008444742243,4.000053838423356,9.000609362487232,16.003393836191155}
```

Find 5 eigenvalues and eigenvectors of a Laplacian with periodic boundary conditions:

```wolfram
{vals,funs}=NDEigensystem[{-u''[x],u[0]==u[2π]},u[x],{x,0,2π},5];
```

Compare the eigenvalues with the expected analytical eigenvalues:

```wolfram
vals-Flatten[Join[{0},Table[{n^2,n^2},{n,1,2}]]]
(* Output *)
{0.,0.00001345960579635097,0.000013459605862076174,0.0008484590477646492,0.0008484590477717546}
```

Visualize the eigenfunctions:

```wolfram
Plot[funs,{x,0,2π}]
```

*([Graphics])*

Inspect to see that the bounds are periodic:

```wolfram
(funs/.x->0)-(funs/.x->2π)
(* Output *)
{5.551115123125783×10^-17,0.,-1.1102230246251565×10^-16,-5.551115123125783×10^-17,1.1102230246251565×10^-16}
```

An equivalent formulation using [PeriodicBoundaryCondition](https://reference.wolfram.com/language/ref/PeriodicBoundaryCondition.html):

```wolfram
{vals,funs}=NDEigensystem[{-u''[x],PeriodicBoundaryCondition[u[x],x==2π,Function[x,x-2π]]},u[x],{x,0,2π},5];
```

Specify a wave equation with homogeneous Dirichlet boundary conditions:

```wolfram
{ℒ,ℬ}={D[u[t,x],t,t]==Laplacian[u[t,x],{x}],DirichletCondition[u[t,x]==0,True]};
```

Find the 4 smallest eigenvalues and eigenfunctions:

```wolfram
{vals,funs}=NDEigensystem[{ℒ,ℬ},u[t,x],t,{x,0,π},4];
```

```wolfram
vals
(* Output *)
{-1.1104114731174438×10^-16-1.0000004222370051 ⅈ,-1.1104114731174438×10^-16+1.0000004222370051 ⅈ,-9.162451867121918×10^-16-2.00001345956054 ⅈ,-9.162451867121918×10^-16+2.00001345956054 ⅈ}
```

Compute the eigensystem of a generalized wave equation $\frac{\partial^{2}u}{\partial t^{2}}+\gamma \frac{\partial u}{\partial t} - c^{2}\frac{\partial^{2}u}{\partial x^{2}}+\gamma u=0$:

```wolfram
γ=1.3;c=1.1;
{vals,funs}=NDEigensystem[D[u[t,x],{t,2}]+γ D[u[t,x],{t,1}]-c^2 D[u[t,x],{x,2}]+γ u[t,x]==0,u,t,{x,0,π},4];
```

```wolfram
vals
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

Compute the 4 smallest eigenvalues and eigenfunctions:

```wolfram
{vals,funs}=NDEigensystem[ℒ,u[x],{x,-1,1},4];
```

```wolfram
vals
(* Output *)
{-1.2345680033821396×10^-15,2.000000000000004,5.999999999999998,12.00014396503262}
```

Compare the eigenvalues with the analytical eigenvalues:

```wolfram
(n+1)n/.n->{0,1,2,3}
(* Output *)
{0,2,6,12}
```

Visualize the eigenfunctions:

```wolfram
Plot[funs,{x,-1,1}]
```

*([Graphics])*

Visualize the analytical eigenfunctions:

```wolfram
Plot[Evaluate[LegendreP[#,x]&/@{0,1,2,3}],{x,-1,1}]
```

*([Graphics])*

#### 2D

Specify a Laplacian operator:

```wolfram
ℒ=-Laplacian[u[x,y],{x,y}];
```

Find the 4 smallest eigenvalues and eigenfunctions:

```wolfram
{vals,funs}=NDEigensystem[ℒ,u[x,y],{x,0,π},{y,0,π},4];
```

```wolfram
vals
(* Output *)
{2.2478597284837017×10^-14,1.0000015382281866,1.0000017695290675,2.0000125935329556}
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

Find the 4 smallest eigenvalues and eigenfunctions of the operator in a unit disk:

```wolfram
{vals,funs}=NDEigensystem[ℒ,u[x,y],{x,y}∈Disk[],4];
```

```wolfram
vals
(* Output *)
{1.5295779233167×10^-14,3.3899670948303995,3.3899672003966015,9.328561151575382}
```

Visualize the eigenfunctions:

```wolfram
ContourPlot[#,{x,y}∈Disk[]]&/@funs
(* Output *)
{[Graphics],[Graphics],[Graphics],[Graphics]}
```

Specify a Laplacian operator with homogeneous Dirichlet boundary conditions:

```wolfram
{ℒ,ℬ}={-Laplacian[u[x,y],{x,y}],DirichletCondition[u[x,y]==0,True]};
```

Find the 9 smallest eigenvalues and eigenfunctions in a rectangle:

```wolfram
{vals,funs}=NDEigensystem[{ℒ,ℬ},u[x,y],{x,0,π},{y,0,π},9];
```

```wolfram
vals
(* Output *)
{2.0000128656842495,5.000200550012288,5.0002113917591755,8.000798754375737,10.001479784631819,10.001664582371912,13.003262795810858,13.003586436347206,17.007275353139207}
```

Visualize the eigenfunctions:

```wolfram
Plot3D[#,{x,0,π},{y,0,π}]&/@funs
(* Output *)
![image](img/image_005.png)
```

Specify a wave equation with homogeneous Dirichlet boundary conditions:

```wolfram
{ℒ,ℬ}={D[u[t,x,y],t,t]==-Laplacian[u[t,x,y],{x,y}],DirichletCondition[u[t,x,y]==0.,True]};
```

Find the 4 smallest eigenvalues and eigenfunctions in a disk:

```wolfram
{vals,funs}=NDEigensystem[{ℒ,ℬ},u[t,x,y],t,{x,y}∈Disk[],4];
```

```wolfram
vals
(* Output *)
{-2.404833953589233,2.4048339535892334,3.831797113388663,-3.8317971133886655}
```

Visualize the eigenfunctions:

```wolfram
Plot3D[#,{x,y}∈Disk[]]&/@funs
(* Output *)
![image](img/image_007.png)
```

The eigenvalues of the wave equation will be the square root of the angular frequencies:

```wolfram
{vals,funs}=NDEigensystem[{D[u[t,x],{t,2}]==Laplacian[u[t,x],{x}]},{u[t,x]},t,{x,0,π},6];
```

```wolfram
vals
(* Output *)
{-3.355111767626083×10^-18-1.1488257254232319×10^-7 ⅈ,-3.355111767626083×10^-18+1.1488257254232319×10^-7 ⅈ,1.1100364512453329×10^-16-1.0000004222369983 ⅈ,1.1100364512453329×10^-16+1.0000004222369983 ⅈ,-8.607332883226083×10^-16-2.000013459560535 ⅈ,-8.607332883226083×10^-16+2.000013459560535 ⅈ}
```

Compare to the exact solution of an equivalent first-order system of ordinary differential equations:

```wolfram
Transpose[Eigenvalues[({{0, 1}, {-ω^2, 0}})]/.ω->{0,1,2}]
(* Output *)
{{0,0},{-ⅈ,ⅈ},{-2 ⅈ,2 ⅈ}}
```

Solve a partially constrained eigensystem problem:

```wolfram
{vals,funs}=NDEigensystem[{-Laplacian[u[x,y],{x,y}],DirichletCondition[u[x,y]==0.,x==0]},u[x,y],{x,y}∈Rectangle[{0,0},{1,1}],9];
```

```wolfram
vals
(* Output *)
{2.4674012306269595,12.337013968499967,22.206704686385642,32.0763175985178,41.9463502339355,61.68565651704908,61.687048683700525,71.55666280624403,91.2998551885231}
```

Visualize the eigenfunctions:

```wolfram
Show[
Graphics3D[Line[{{0,0,0},{1,0,0},{1,1,0},{0,1,0},{0,0,0}}]],
Plot3D[#,{x,y}∈Rectangle[{0,0},{1,1}],ColorFunction->"TemperatureMap",Axes->False],Boxed->False,AspectRatio->1]&/@funs
(* Output *)
![image](img/image_009.png)
```

### Options

#### Method

"Eigensystem"  (4)

The discretized PDE is solved with [Eigensystem](https://reference.wolfram.com/language/ref/Eigensystem.html) or [Eigenvalues](https://reference.wolfram.com/language/ref/Eigenvalues.html). All [Eigensystem](https://reference.wolfram.com/language/ref/Eigensystem.html) options available in [NDEigensystem](https://reference.wolfram.com/language/ref/NDEigensystem.html) are discussed in the following.

Specify a method to use for solving the eigensystem:

```wolfram
AbsoluteTiming[First[NDEigensystem[{-Laplacian[u[x,y],{x,y}]},u,{x,0,1},{y,0,1},4,Method->{"Eigensystem"->"Direct"}]]]
(* Output *)
{1.315071,{0.,9.869621049921076,9.869621616269733,19.73934710540828}}
```

In this case, the default method is faster:

```wolfram
AbsoluteTiming[First[NDEigensystem[{-Laplacian[u[x,y],{x,y}]},u,{x,0,1},{y,0,1},4]]]
(* Output *)
{0.115699,{-3.6652834904648225×10^-14,9.869621049920202,9.869621616268848,19.739347105408093}}
```

Arnoldi is used as the default method:

```wolfram
AbsoluteTiming[First[NDEigensystem[{-Laplacian[u[x,y],{x,y}]},u,{x,0,1},{y,0,1},4,Method->{"Eigensystem"->"Arnoldi"}]]]
(* Output *)
{0.105259,{-3.6652834904648225×10^-14,9.869621049920202,9.869621616268848,19.739347105408093}}
```

The Arnoldi method is numerically more efficient but less stable than the direct method, which is usable only for small-scale problems.

Specify a maximum number of iterations for the Arnoldi method:

```wolfram
AbsoluteTiming[First[NDEigensystem[{-Laplacian[u[x,y],{x,y}]},u,{x,0,1},{y,0,1},4,Method->{"Eigensystem"->{"Arnoldi","MaxIterations"->100}}]]]
(* Output *)
{0.087271,{-3.6652834904648225×10^-14,9.869621049920202,9.869621616268848,19.739347105408093}}
```

The FEAST method is suitable for finding eigenvalues and eigenvectors within a band. Find two eigenvalues and eigenfunctions of a Sturm-Liouville operator within the band of $(5,6)$ with the FEAST method for [Eigensystem](https://reference.wolfram.com/language/ref/Eigensystem.html):

```wolfram
{vals,funs}=NDEigensystem[{u''[x]+RealAbs[x] u[x],DirichletCondition[u[x]==0,True]},u[x],{x,-8,8},2,Method->{"SpatialDiscretization"->{"FiniteElement",{"MeshOptions"->{"MaxCellMeasure"->0.05}}},"Eigensystem"->{"FEAST","Interval"->{5,6}}}]
(* Output *)
{{5.661892563705724,5.661892573969445},{InterpolatingFunction[...][x],InterpolatingFunction[...][x]}}
```

According to the Sturm-Liouville theory, the eigenvalues must be distinct, but for this example they are close to degenerate:

```wolfram
Subtract@@vals
(* Output *)
-1.026372142831633×10^-8
```

The almost degenerate eigenfunctions are:

```wolfram
Plot[funs,{x,-8,8}]
```

*([Graphics])*

The interval end points $(5,6)$ are not included in the interval FEAST finds eigenvalues in; for more information, please refer to the [Eigensystem](https://reference.wolfram.com/language/ref/Eigensystem#155953046.html) reference page.

The usage of the `"Shift"` option is [explained in an example below](https://reference.wolfram.com/language/ref/NDEigensystem.html#42461913).

"Interpolation"  (1)

By default, the eigenfunctions give [Indeterminate](https://reference.wolfram.com/language/ref/Indeterminate.html) for points outside the region $\Omega$:

```wolfram
Ω=Disk[];{vals,funs}=NDEigensystem[{-Laplacian[u[x,y],{x,y}]},u,{x,y}∈Ω,2];
```

```wolfram
funs[[2]][1,1]
(* Output *)
Indeterminate
```

Modify the default behavior to eliminate the warning messages:

```wolfram
{vals,funs}=NDEigensystem[{-Laplacian[u[x,y],{x,y}]},u,{x,y}∈Ω,2,Method->{"Interpolation"->{"ExtrapolationHandler"->{(Indeterminate&),"WarningMessage"->False}}}];
```

Evaluate the second eigenfunction outside of the region $\Omega$:

```wolfram
funs[[2]][1,1]
(* Output *)
Indeterminate
```

Visualize part of the solution:

```wolfram
ContourPlot[funs[[2]][x,y],{x,0,1},{y,0,1}]
```

*([Graphics])*

Change the eigenfunctions to return an extrapolated value outside of the region $\Omega$:

```wolfram
{vals,funs}=NDEigensystem[{-Laplacian[u[x,y],{x,y}]},u,{x,y}∈Ω,2,Method->{"Interpolation"->{"ExtrapolationHandler"->{Automatic,"WarningMessage"->False}}}];
```

Evaluate the second eigenfunction outside of the region $\Omega$:

```wolfram
funs[[2]][1,1]
(* Output *)
0.6065013082536627
```

"PDEDiscretization"  (1)

Change the [MaxCellMeasure](https://reference.wolfram.com/language/ref/MaxCellMeasure.html) for the underlying computation:

```wolfram
{vals,funs}=NDEigensystem[{-Laplacian[u[x],{x}]},u,{x,0,π},4,Method->{"PDEDiscretization"->{"FiniteElement",{"MeshOptions"->{"MaxCellMeasure"->0.5}}}}];
```

The exact eigenvalues are `0,1,4,9,…`, so the eigenvalue error is as follows:

```wolfram
err=Abs[vals - {0,1,4,9}]
(* Output *)
{0.,0.00005576174653354471,0.003462494088496193,0.03762096621904121}
```

A finer mesh results in a decreased discretization error:

```wolfram
{vals,funs}=NDEigensystem[-Laplacian[u[x],{x}],u[x],{x,0,π},4,Method->{"PDEDiscretization"->{"FiniteElement","MeshOptions"->{"MaxCellMeasure"->0.01}}}];
```

```wolfram
err=Abs[vals-{0,1,4,9}]
(* Output *)
{3.143239121948122×10^-12,1.071276400921306×10^-11,8.770495441012827×10^-10,1.00124708524163×10^-8}
```

"VectorNormalization"  (3)

Compute without any normalization of the computed eigenvalues:

```wolfram
{vals,funs}=NDEigensystem[-Laplacian[u[x],{x}],u[x],{x,0,π},4,Method->{"VectorNormalization"->None}];
```

```wolfram
NIntegrate[#^2,{x,0,π}]&/@funs
(* Output *)
{0.07662421106316575,0.07479974990831866,0.07479782758346251,0.07479153511462865}
```

The default is to normalize with respect to the system matrices:

```wolfram
{vals,funs}=NDEigensystem[-Laplacian[u[x],{x}],u,{x,0,π},4,Method->{"VectorNormalization"->Function[{values,vectors,stiffness,damping},
s=stiffness;d=damping;
norm=vectors/Diagonal[vectors.damping.ConjugateTranspose[vectors]]^(1/2)]}];
```

The functions are constructed by interpolating the eigenvectors:

```wolfram
Max[Abs[#["ValuesOnGrid"]&/@funs-norm]]
(* Output *)
0.
```

[NDEigensystem](https://reference.wolfram.com/language/ref/NDEigensystem.html) constructs the damping (mass) and stiffness system matrices to be solved as a generalized eigensystem. The normalization of the eigenvectors is constructed such that applying the normalized eigenvector to the damping (mass) matrix forms an identity matrix:

```wolfram
Chop[norm.d.ConjugateTranspose[norm]]
(* Output *)
{{1.,0,0,0},{0,1.,0,0},{0,0,0.9999999999999997,0},{0,0,0,0.9999999999999999}}
```

Applying the normalized eigenvector to the stiffness matrix produces eigenvalues on the diagonal:

```wolfram
Chop[Diagonal[norm.s.ConjugateTranspose[norm]]-vals]
(* Output *)
{0,0,0,0}
```

With the default normalization, the L2 norm of the eigenfunctions will be approximately 1:

```wolfram
NIntegrate[#[x]^2,{x,0,π}]&/@funs
(* Output *)
{1.0000000000000009,1.0000003815928011,0.9999953617073144,0.9999998810301433}
```

Degrees of freedom associated with Dirichlet conditions are removed before the eigensystem is computed, and 0 is reinserted after the solution. This results in different dimensions of the normalization vector and the eigenvector for constraint problems:

```wolfram
First[{vals,funs}=NDEigensystem[{-Laplacian[u[x],{x}],DirichletCondition[u[x]==0,True]},u,{x,0,π},4,Method->{"VectorNormalization"->Function[{values,vectors,stiffness,damping},
s=stiffness;d=damping;
norm=vectors/Diagonal[(vectors.damping.ConjugateTranspose[vectors])]^(1/2)]}]];
```

During the eigensystem computation, the two Dirichlet conditions are removed:

```wolfram
Dimensions[norm]
(* Output *)
{4,39}
```

After the eigensystem computation, zero values are reintroduced into the functions:

```wolfram
Dimensions[#["ValuesOnGrid"]&/@funs]
(* Output *)
{4,41}
```

### Applications

#### Acoustics

Compute the acoustic eigenvalues and eigenfunctions for an approximation of a cross section through a [Mini](https://en.wikipedia.org/wiki/Mini). Import an image of the cross section:

```wolfram
img=Import["http://upload.wikimedia.org/wikipedia/commons/d/d3/Mini_cross_section.jpg"]
```

![image](img/image_011.png)

Use the mask tool to create a boundary graphic:

```wolfram
boundary=[Graphics];
```

Discretize the graphic:

```wolfram
bdr=BoundaryDiscretizeGraphics[boundary]
```

*([Graphics])*

Compute 6 eigenvalues and eigenfunctions of the cross section:

```wolfram
{vals,funs}=NDEigensystem[{-Laplacian[u[x,y],{x,y}]},u[x,y],{x,y}∈bdr,6];
```

```wolfram
vals
(* Output *)
{1.508326094631877×10^-20,1.2707293933341714×10^-6,3.956142113696531×10^-6,6.120407377041952×10^-6,9.976858517360204×10^-6,0.000011841218019956618}
```

Visualize the second eigenfunction in the cross section of the car:

```wolfram
Show[img,ContourPlot[funs[[2]],{x,y}∈bdr,Axes->None,Frame->None,AspectRatio->Automatic,ColorFunction->Function[f,{Opacity[0.75],ColorData["TemperatureMap"][f]}]],ImageSize->Automatic]
```

![image](img/image_012.png)

The visualization shows qualitatively a standing wave assuming acoustically rigid walls. A large amplitude corresponds to a large change in pressure and thus a large acoustic volume.

#### Structural Mechanics

Specify a plane stress PDE:

```wolfram
ps=SolidMechanicsPDEComponent[{{u[x,y],v[x,y]},{x,y}},<|"YoungModulus"->10^3,"PoissonRatio"->33/100,"Thickness"->1|>];
```

Compute constraint eigenvalues and eigenfunctions:

```wolfram
{vals,funs}=NDEigensystem[{ps,DirichletCondition[{u[x,y]==0.,v[x,y]==0.},x==0]}, {u[x,y],v[x,y]},{x,y}∈Rectangle[{0,0},{5,1}],9];
```

```wolfram
vals
(* Output *)
{1.5584518919016852,45.29890210584211,99.21385118911388,257.6038362664441,713.0087913127295,886.4698406606493,1446.6716761256569,2416.979263774896,2439.785106643841}
```

Visualize the scaled norm of the eigenfunctions:

```wolfram
Show[{Graphics3D[
{Gray,GraphicsComplex[{{0,0,0},{5,0,0},{5,1,0},{0,1,0}},Line[{{1,2,3,4,1}}]]}],
Plot3D[Sqrt[Total[#^2]],{x,y}∈Rectangle[{0,0},{5,1}],ColorFunction->"TemperatureMap",Axes->False,Mesh->False]},Boxed->False]&/@funs
(* Output *)
![image](img/image_013.png)
```

The visualization shows the principal deformation of the region.

Specify variables and parameters for a plane stress model for steel:

```wolfram
vars={{u[t,x,y],v[t,x,y]},t,{x,y}};
pars=<|"YoungModulus"->210 10^9,"PoissonRatio"->3/10,"MassDensity"->7850,"Thickness"->1,"AnalysisType"->"Eigenmode"|>;
```

Compute the first eigenvalue and eigenfunction of a beam with length 10 meters and height 0.7 meters:

```wolfram
length=10;height=0.7;
{vals,funs}=NDEigensystem[{SolidMechanicsPDEComponent[vars,pars]=={0,0},SolidFixedCondition[x==0,vars,pars]}, {u,v},t,{x,y}∈Rectangle[{0,0},{length,height}],1];
```

Compute resonant frequency from the eigenvalue:

```wolfram
nfreq=Sqrt[vals[[1]]]/(2 π)
(* Output *)
5.831448270426073
```

Compute the first resonance frequency analytically utilizing Euler-Bernoulli beam theory:

```wolfram
intertia=(width*(height^3))/12.;
freq=0.56/(length^2)*Sqrt[(intertia*pars["YoungModulus"])/(height*width*pars["MassDensity"])]
(* Output *)
5.852888665649132
```

Compare the difference in results; about half the difference is attributable to the different theoretical approaches:

```wolfram
freq-nfreq
(* Output *)
0.02144039522305974
```

#### Eigenfunction Expansions

Construct a solution for a plucked clamped string using an eigenfunction expansion:

```wolfram
ρ[x_] = 1+x-x^2;{vals,funs}=NDEigensystem[{D[-ρ[x] D[u[x],x],x],DirichletCondition[u[x]==0,True]},u[x],{x,0,1},10];
```

```wolfram
vals
(* Output *)
{10.970577444088587,45.396103372762276,102.78978978602873,183.16544140683675,286.57239649199755,413.11230101079644,562.9632314374769,736.406341070613,933.8531379050805,1155.8714141990072}
```

Approximate the initial condition $1-\left|2 x-1\right|$ as a linear combination of eigenfunctions:

```wolfram
f[x_]:= 1-Abs[2x-1];
coeffs = Chop[Quiet[Map[NIntegrate[# f[x],{x,0,1}]&, funs]]]
(* Output *)
{0.5717871741189454,0,0.07467357026055196,0,-0.023194336374750928,0,-0.012845461783552888,0,-0.007278243657362491,0}
```

Construct the time-dependent solution using the coefficients and the eigenvalues as frequencies, assuming that the time derivative at $t=0$ is 0:

```wolfram
tap[t_,x_] = (Cos[Sqrt[vals] t] coeffs).funs;
```

Make an animation of the approximation, including each of the functions with nonzero coefficients:

```wolfram
cnz = Flatten[Position[coeffs,c_ /; c ≠ 0]];
```

```wolfram
plots = Table[Show[Plot[tap[t,x],{x,0,1},PlotStyle->{Black,Thickness[0.01]},PlotRange->{-1,1}],Plot[Evaluate[Cos[Sqrt[vals[[cnz]]] t]coeffs[[cnz]] funs[[cnz]]],{x,0,1},PlotRange->{-1,1}]],{t,0,1,.025}];
```

```wolfram
ListAnimate[plots]
```

Compare to the time-dependent solution found by [NDSolve](https://reference.wolfram.com/language/ref/NDSolve.html):

```wolfram
f[x_]:= 1-Abs[2x-1];
ρ[x_] = 1+x-x^2;tu = NDSolveValue[{D[u[t,x],t,t]==D[ρ[x] D[u[t,x],x],x],u[0,x]==f[x],Derivative[1,0][u][0,x]==0,DirichletCondition[u[t,x]==0,True]},u,{t,0,1},{x,0,1},Method->{"PDEDiscretization"->"MethodOfLines"}];
```

Make a plot of the difference:

```wolfram
Table[Plot[tu[t,x]-tap[t,x],{x,0,1}],{t,{0,.25,.5,1}}]
(* Output *)
![image](img/image_015.png)
```

Compare different harmonic frequencies on a circular drum:

```wolfram
Ω=Disk[];
{vals, funs} = NDEigensystem[{-Laplacian[u[x,y],{x,y}],DirichletCondition[u[x,y]==0,True]},u,{x,y}∈Ω,6];
```

```wolfram
vals
(* Output *)
{5.783226344331008,14.682669114667387,14.6826745954688,26.37839354203904,26.379103017216053,30.478731698301047}
```

Show a time evolution for each of the eigenfunctions over one cycle of the first:

```wolfram
ranges = Map[Max[Abs[Last[PlotRange /. Options[Plot3D[#[x,y],{x,y}∈Ω]]]]]&,funs];
row[t_] := Grid[{MapThread[Plot3D[Sin[Sqrt[#1] t] #2[x,y] ,{x,y}∈Ω,PlotRange->#3{-1,1},Axes->False,Mesh->False]&,{vals,funs,ranges}]}]
```

```wolfram
rows = Table[row[t],{t,0,2 Pi/Sqrt[vals[[1]]],.1}];
```

Reduce the size of the graphics:

```wolfram
rows=Map[Rasterize[#,"Image",ImageResolution->80]&,rows];
```

Note that the eigenfunctions with larger eigenvalues oscillate more quickly:

```wolfram
ListAnimate[rows]
```

![image](img/image_017.png)

Show the time evolution of a linear combination of eigenfunctions:

```wolfram
lct[t_,x_,y_]:=4(Cos[Sqrt[vals] t]/vals).Map[#[x,y]&,funs]
```

```wolfram
plots = Table[Plot3D[lct[t,x,y],{x,y}∈Ω,PlotRange->{-1.2,1.2},Mesh->False],{t,0,10,.25}];
ListAnimate[plots]
```

This is the solution of the wave equation with consistent initial conditions:

```wolfram
w = NDSolveValue[{D[u[t,x,y],t,t]== Laplacian[u[t,x,y],{x,y}],DirichletCondition[u[t,x,y]==0,True],u[0,x,y]==lct[0,x,y],(D[u[t,x,y],t] /.t->0)== (D[lct[t,x,y],t] /. t->0)},u,{t,0,10},{x,y}∈Ω];
```

Show the difference between a finite eigenexpansion and full solution:

```wolfram
Block[{t = 10},Plot3D[w[t,x,y] - lct[t,x,y],{x,y}∈Ω]]
```

*([Graphics3D])*

Use an eigenfunction expansion to get a solution of the heat equation on a region:

```wolfram
Ω=ImplicitRegion[x^2+y^2>=1/4,{{x,-1,1},{y,-1,1}}];
```

Compute the first 20 eigenvalues and eigenfunctions for the Laplacian on $\Omega$:

```wolfram
{vals, funs} = NDEigensystem[{Laplacian[u[x,y],{x,y}],DirichletCondition[u[x,y]==0,x==1||x==-1||y==1||y==-1]},u[x,y],{x,y}∈Ω,20];
```

```wolfram
vals
(* Output *)
{-9.611839984124455,-11.796051473941674,-11.796095933976925,-16.54077515108097,-20.976623561140293,-29.482980990806936,-29.483176852841982,-39.9563909667683,-44.765551730092454,-51.81234263867105,-51.81397714940477,-53.87582589670583,-61.54658132588645,-67.23451449861948,-67.2381114382582,-74.84000645916699,-86.06196684456357,-86.08077415536471,-86.08662442981034,-93.19577053845879}
```

```wolfram
Map[Plot3D[#,{x,y}∈Ω,Mesh->False]&,funs]
(* Output *)
![image](img/image_018.png)
```

Let $u_{0}$ be an initial condition:

```wolfram
u_0=(x^2-1)(y^2-1);
```

For computing integrals, it is faster and more consistent to use the discretized domain `Ω_d`:

```wolfram
Ω_d=Head[funs[[1]]]["ElementMesh"];
```

Since $u_{0}$ does not satisfy the homogeneous boundary conditions, find the steady state solution $u_{s}$ that does on the discretized domain `Ω_d`, so that $u_{0}-u_{s}$ does satisfy them:

```wolfram
u_s=NDSolveValue[{Laplacian[u[x,y],{x,y}] == NeumannValue[Grad[u_0,{x,y}].{2x,2y}, x^2+y^2==1/4],DirichletCondition[u[x,y]==u_0,x==1||x==-1||y==1||y==-1]},u[x,y],{x,y}∈Ω_d];
```

```wolfram
Map[Plot3D[#,{x,y}∈Ω_d]&,{u_0,u_s,u_0-u_s}]
(* Output *)
![image](img/image_020.png)
```

Compute the coefficients of the eigenfunction expansion for the difference $u_{0}-u_{s}$:

```wolfram
coeffs=Map[NIntegrate[(u_0-u_s) #,{x,y}∈Ω_d]&,funs]
(* Output *)
{-0.44586288504179283,4.677292822665986×10^-6,-6.868962223915444×10^-6,-4.7859285256320223×10^-8,1.0385328058225458×10^-7,-6.513972558134812×10^-7,1.134359344803374×10^-6,-0.00423297087144628,-2.280675416285792×10^-6,4.067048772199739×10^-6,2.3672881595776605×10^-6,-4.620393796229451×10^-6,-0.013134309160629605,-6.122550925709465×10^-6,2.054963726763246×10^-6,3.1384819957689227×10^-6,-4.177940701115485×10^-6,-7.661170816313629×10^-6,2.7150424702555333×10^-6,0.003354570211716862}
```

The approximation to the solution of the heat equation is as follows:

```wolfram
approx[t_,x_,y_]=u_s+(coeffs Exp[vals t]).funs;
```

Compute the time-dependent solution:

```wolfram
ut=NDSolveValue[{D[u[t,x,y],t]==Laplacian[u[t,x,y],{x,y}]-NeumannValue[Grad[u_0,{x,y}].{2x,2y}, x^2+y^2==1/4],u[0,x,y]==u_0,DirichletCondition[u[t,x,y]==0,x==1||x==-1||y==1||y==-1]},u,{t,0,1},{x,y}∈Ω];
```

Make a plot of the difference:

```wolfram
Table[Plot3D[{ut[t,x,y]-approx[t,x,y]},{x,y}∈Ω,PlotRange->All],{t,0,0.05,0.025}]
(* Output *)
![image](img/image_022.png)
```

Eigenfunctions in an L-shaped region:

```wolfram
L=Polygon[{{1,0},{2,0},{2,2},{0,2},{0,1},{1,1}}];
```

```wolfram
{vals,funs}=NDEigensystem[{Laplacian[u[x,y],{x,y}],DirichletCondition[u[x,y]==0.,True]},u[x,y],{x,y}∈L,6];
```

The first six eigenfunctions:

```wolfram
ContourPlot[#,{x,y}∈L,PlotPoints->35]&/@funs
(* Output *)
![image](img/image_024.png)
```

The first one:

```wolfram
Plot3D[funs[[1]],{x,y}∈L,PlotPoints->75,Mesh->None,PlotStyle->Directive[Orange,Specularity[White,30]],BoxRatios->{1,1,0.8}]
(* Output *)
![image](img/image_026.png)
```

#### Interval of Eigenvalues and Eigenfunctions

Find an eigenvalue and a eigenfunction in an interval:

```wolfram
nds=NDEigensystem[{-Laplacian[u[x,y],{x,y}],DirichletCondition[u[x,y]==0,True]},u,{x,y}∈ImplicitRegion[(x^2+y^2+2 y)^2<4 (x^2+y^2),{x,y}],1,Method->{"Eigensystem"->{"FEAST","Interval"->{400,405}},"PDEDiscretization"->{"FiniteElement","MeshOptions"->{"MaxCellMeasure"->0.001}}}]
(* Output *)
{{400.7913586788476},{InterpolatingFunction[]}}
```

Visualize the eigenfunction found:

```wolfram
ContourPlot[Evaluate[Abs[nds[[2,1]][x,y]]^(1/2)],{x,y}∈nds[[2,1]]["ElementMesh"],PlotPoints->200,ColorFunction->(GrayLevel[1-#]&),AspectRatio->Automatic,ContourStyle->None,PlotRange->All,MaxRecursion->0]
(* Output *)
![image](img/image_028.png)
```

#### Quantum Mechanics

Specify a Schrödinger operator with parameter $h$ and potential $V$:

```wolfram
h=1/10;V[x_]:=x^2
ℒ=-h^2*u''[x]+V[x]*u[x];
```

Find the 10 smallest eigenvalues and eigenfunctions on a refined mesh:

```wolfram
{vals,funs}=NDEigensystem[ℒ,u[x],{x,-3,3},10,Method->{"PDEDiscretization"->{"FiniteElement",{"MeshOptions"->{MaxCellMeasure->0.05}}}}];
```

```wolfram
vals
(* Output *)
{0.10000016225988075,0.30000113385045624,0.5000040417329513,0.7000101643270646,0.9000207683332282,1.1000371088172063,1.3000604292934648,1.5000919618072899,1.700132927015669,1.9001845342669068}
```

Visualize the eigenfunctions scaled by $h$ and offset by the respective eigenvalue:

```wolfram
Show[Plot[Evaluate[h*funs+vals],{x,-3,3}],
Plot[V[x],{x,-3,3}],PlotRange->{{-3,3},{0,2}},AxesOrigin->{-3,0}]
(* Output *)
![image](img/image_030.png)
```

This example shows the avoided-crossing phenomenon of eigenvalues. For this, a region difference between two disks is computed The outer disk has a radius of $R$ and its center is at $\{0,0 \}$. The inner disk has a radius of $R/2$. The center of the inner disk starts from the central $\{0,0 \}$ position and slowly moves to the right by $\{xPos,0 \}$. Write a helper function:

```wolfram
fun[xPos_]:=Module[{Ω,radius=1},
Ω=RegionDifference[Disk[{0,0},radius],Disk[{xPos,0},radius/2]];
NDEigensystem[{-Laplacian[u[x,y],{x,y}],DirichletCondition[u[x,y]==0,True]},u[x,y],{x,y}∈Ω,20,Method->{"PDEDiscretization"->{"FiniteElement",{"MeshOptions"->{"MaxCellMeasure"->0.001}}}}]
]
```

Compute the eigenvalues and vectors:

```wolfram
Monitor[results=Table[{{xr,#}&/@#[[1]],#[[2]]}&[fun[xr]],{xr, 0, 0.49, 0.49/200}], xr];
```

Visualize the eigenvalues of the Laplacian of the geometry at the various positions of $xPos$:

```wolfram
ListLinePlot[Transpose[results[[All,1]]]]
```

*([Graphics])*

Note that some curves cross, while others show an avoided-crossing phenomenon.

Visualize the eigenvectors for a particular solution:

```wolfram
solutionNumber=Length[results];
DensityPlot[Abs[#]^2,{x,y}∈#[[0]]["ElementMesh"],ColorFunction->"TemperatureMap", Frame->False,
PlotRange->All]&/@results[[All,2]][[solutionNumber]]
(* Output *)
![image](img/image_032.png)
```

### Properties & Relations

With the default normalization, the eigenfunctions `*u*_*i*` approximately satisfy $\int_{\Omega}u_{i}\mathit{u_{i}}^{*}=1$:

```wolfram
{vals,funs}=NDEigensystem[-Laplacian[u[x],{x}],u[x],{x,0,π},4];
```

```wolfram
vals
(* Output *)
{4.080266017155255×10^-14,1.0000008444742063,4.000053838423353,9.000609362487243}
```

```wolfram
NIntegrate[#^2,{x,0,π}]&/@funs
(* Output *)
{1.0000000000000009,1.0000003815928011,0.9999953617073144,0.9999998810301433}
```

When the eigenvalues are complex, it is necessary to use the conjugate:

```wolfram
{vals,funs}=NDEigensystem[D[u[t,x],t,t]==Laplacian[u[t,x],{x}],u[t,x],t,{x,0,π},4];
```

```wolfram
vals
(* Output *)
{-8.237588075998181×10^-18-1.1225319332930581×10^-7 ⅈ,-8.237588075998181×10^-18+1.1225319332930581×10^-7 ⅈ,1.026938413922829×10^-15-1.0000004222369996 ⅈ,1.026938413922829×10^-15+1.0000004222369996 ⅈ}
```

```wolfram
NIntegrate[# Conjugate[#],{x,0,π}]&/@funs
(* Output *)
{0.9999999999999898,0.9999999999999898,0.4999999796778648,0.4999999796778648}
```

Compare the approximated eigensystem to one computed analytically for the wave equation.

Find the 6 smallest eigenvalues and eigenfunctions of a wave equation $u_{tt}=\mathit{u}_{xx}$ within 0 and $\pi$:

```wolfram
{vals,funs}=NDEigensystem[{D[u[t,x],{t,2}]==Laplacian[u[t,x],{x}]},{u[t,x]},t,{x,0,π},6];
```

```wolfram
vals
(* Output *)
{1.1406479921142865×10^-17-1.122427229034216×10^-7 ⅈ,1.1406479921142865×10^-17+1.122427229034216×10^-7 ⅈ,9.714272158138526×10^-16-1.0000004222369994 ⅈ,9.714272158138526×10^-16+1.0000004222369994 ⅈ,1.4152545456455693×10^-15-2.000013459560544 ⅈ,1.4152545456455693×10^-15+2.000013459560544 ⅈ}
```

The eigenvalues of the PDE system converted to first order in time correspond to eigenvalues of the matrix $(\begin{pmatrix}
0 & 1 \\
-\omega^{2} & 0
\end{pmatrix})$, where $\omega^{2}$ is an eigenvalue of $u_{xx}=-\omega^{2}u$ with $u_{x}=0$ at the boundaries:

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

Show the relationship between higher-order time-dependent PDEs and systems of first-order PDEs.

Find the 6 smallest eigenvalues and eigenfunctions of the wave equation $u_{tt}=u_{xx}$ within 0 and $\pi$:

```wolfram
{vals,funs}=NDEigensystem[{D[u[t,x],{t,2}]==Laplacian[u[t,x],{x}]},{u[t,x]},t,{x,0,π},6];
```

```wolfram
vals
(* Output *)
{1.1406479921142865×10^-17-1.122427229034216×10^-7 ⅈ,1.1406479921142865×10^-17+1.122427229034216×10^-7 ⅈ,9.714272158138526×10^-16-1.0000004222369994 ⅈ,9.714272158138526×10^-16+1.0000004222369994 ⅈ,1.4152545456455693×10^-15-2.000013459560544 ⅈ,1.4152545456455693×10^-15+2.000013459560544 ⅈ}
```

Plot the modulus and the real and imaginary parts of the eigenfunctions:

```wolfram
Map[Plot[{Abs[#],Re[#],Im[#]},{x,0,π},PlotRange->All]&, funs]
(* Output *)
![image](img/image_034.png)
```

Find the 6 smallest eigenvalues and eigenfunctions of a wave equation given as a system of first-order PDEs $(\begin{pmatrix}
u_{t} & = & v \\
v_{t} & = & u_{xx}
\end{pmatrix})$:

```wolfram
{svals,sfuns}=NDEigensystem[{D[u[t,x],t]==v[t,x],D[v[t,x],t]==Laplacian[u[t,x],{x}]},{u[t,x],v[t,x]},t,{x,0,π},6];
```

```wolfram
vals
(* Output *)
{1.1406479921142865×10^-17-1.122427229034216×10^-7 ⅈ,1.1406479921142865×10^-17+1.122427229034216×10^-7 ⅈ,9.714272158138526×10^-16-1.0000004222369994 ⅈ,9.714272158138526×10^-16+1.0000004222369994 ⅈ,1.4152545456455693×10^-15-2.000013459560544 ⅈ,1.4152545456455693×10^-15+2.000013459560544 ⅈ}
```

Plot the modulus and the real and imaginary parts of the eigenfunctions:

```wolfram
Map[Plot[{Abs[#],Re[#],Im[#]},{x,0,π},PlotRange->All]&, sfuns[[All,1]]]
(* Output *)
![image](img/image_036.png)
```

The eigenvalues for the second-order system and the system of first-order equations are the same:

```wolfram
vals-svals
(* Output *)
{0.+0. ⅈ,0.+0. ⅈ,0.+0. ⅈ,0.+0. ⅈ,0.+0. ⅈ,0.+0. ⅈ}
```

For a given eigenvalue, the eigenfunction returned for the second-order system is the same as the first eigenfunction returned for the system:

```wolfram
Table[SameQ[funs[[i]]- sfuns[[i,1]]],{i,6}]
(* Output *)
{True,True,True,True,True,True}
```

For a given eigenvalue, $\lambda$ times the eigenfunction returned for the second-order system is the same as the second eigenfunction returned for the system:

```wolfram
Plot[Evaluate[Table[RealExponent[vals[[i]] funs[[i]]- sfuns[[i,2]]],{i,6}]],{x,0,π},PlotRange->All]
(* Output *)
![image](img/image_038.png)
```

The eigenfunctions for the system are normalized so that $\int_{\Omega}\{u_{i},v_{i},\ldots \}.\{\mathit{u_{i}}^{*},\mathit{v_{i}}^{*},\ldots \}\simeq 1$:

```wolfram
Chop[NIntegrate[#.Conjugate[#],{x,0,π}]&/@sfuns]
(* Output *)
{1.0000000000000022,1.0000000000000022,1.000000381592801,1.000000381592801,0.9999953618280051,0.9999953618280051}
```

The eigenfunctions for the second-order time-dependent system are normalized in the same way as the system, and since $v_{i}=\lambda_{i}u_{i}$, $\int_{\Omega}\{u_{i},v_{i},\ldots \}.\{\mathit{u_{i}}^{*},\mathit{v_{i}}^{*},\ldots \}=(1+\lambda_{i}^{2})\int_{\Omega}\left|u_{i}\right|^{2}\simeq 1$ so $\int_{\Omega}\left|u_{i}\right|^{2}\simeq \frac{1}{(1+\lambda_{i}^{2})}$:

```wolfram
Map[NIntegrate[# Conjugate[#],{x,0,π}]&, funs]
(* Output *)
{0.9999999999999896,0.9999999999999896,0.499999979677865,0.499999979677865,0.19999691883641213,0.19999691883641213}
```

```wolfram
Map[1/(1 + # Conjugate[#])&, vals]
(* Output *)
{0.9999999999999873+0. ⅈ,0.9999999999999873+0. ⅈ,0.49999978888154495+0. ⅈ,0.49999978888154495+0. ⅈ,0.19999784648625493+0. ⅈ,0.19999784648625493+0. ⅈ}
```

### Possible Issues

The computed eigensystem depends on the granularity of the discretization:

```wolfram
{vals,funs}=NDEigensystem[-Laplacian[u[x],{x}],u[x],{x,0,π},4];
```

```wolfram
vals
(* Output *)
{4.080266017155256×10^-14,1.0000008444742072,4.000053838423357,9.000609362487246}
```

The exact eigenvalues are `0,1,4,9,…`, so the eigenvalue error is as follows:

```wolfram
err=Abs[vals - {0,1,4,9}]
(* Output *)
{4.080266017155256×10^-14,8.444742072288847×10^-7,0.00005383842335682232,0.000609362487246301}
```

Visualize the discretization error in $-(u x)=\lambda (u x)$:

```wolfram
Plot[Evaluate[-Laplacian[funs[[#]],{x}]-vals[[#]]*funs[[#]]],{x,0,π}]&/@Range[Length[vals]]
(* Output *)
{[Graphics],[Graphics],[Graphics],[Graphics]}
```

A finer mesh results in a decreased discretization error:

```wolfram
{vals,funs}=NDEigensystem[-Laplacian[u[x],{x}],u[x],{x,0,π},4,Method->{"PDEDiscretization"->{"FiniteElement","MeshOptions"->{"MaxCellMeasure"->0.01}}}];
```

```wolfram
vals
(* Output *)
{-3.143239121948122×10^-12,1.0000000000107128,4.0000000008770495,9.00000001001247}
```

```wolfram
err=Abs[vals - {0,1,4,9}]
(* Output *)
{3.143239121948122×10^-12,1.071276400921306×10^-11,8.770495441012827×10^-10,1.00124708524163×10^-8}
```

Visualize the discretization error in $\mathit{-}\mathit{\nabla_{\{x \}}^{2}}\mathit{u}[\mathit{x}]=\lambda u[x]$:

```wolfram
Plot[Evaluate[-Laplacian[funs[[#]],{x}]-vals[[#]]*funs[[#]]],{x,0,π}]&/@Range[Length[vals]]
(* Output *)
![image](img/image_040.png)
```

The granularity of the discretization of the region influences the number of eigenvalues and eigenfunctions that can be found:

```wolfram
{vals1,funs1}=NDEigensystem[-Laplacian[u[x],{x}],u[x],{x,0,π},100];
Length[vals1]
(* Output *)
NDEigensystem
(* Output *)
41
```

```wolfram
{vals2,funs2}=NDEigensystem[-Laplacian[u[x],{x}],u[x],{x,0,π},100,Method->{"PDEDiscretization"->{"FiniteElement",{"MeshOptions"->{"MaxCellMeasure"->0.01}}},"Eigensystem"->"Direct"}];
Length[vals2]
(* Output *)
100
```

The differences between the eigenvalues and eigenfunctions in a coarser versus a finer discretization tend to increase with the mode number:

```wolfram
vals1[[Range[10]]]-vals2[[Range[10]]]
(* Output *)
{0.,8.444769618032311×10^-7,0.00005383753154841742,0.0006093524572285958,0.0033937799180492334,0.012803298903151017,0.0377259127966596,0.09368372566930105,0.20517901070080313,0.4081235501053726}
```

Eigensystems with inhomogeneous Dirichlet conditions cannot be solved:

```wolfram
NDEigensystem[{-Laplacian[u[x],{x}],DirichletCondition[u[x]==1,True]},u[x],{x,0,1},2]
(* Output *)
NDEigensystem
(* Output *)
NDEigensystem[{-u^′′[x],DirichletCondition[u[x]==1,True]},u[x],{x,0,1},2]
```

Eigensystems with homogeneous Dirichlet conditions can be solved:

```wolfram
First[NDEigensystem[{-Laplacian[u[x],{x}],DirichletCondition[u[x]==0,True]},u[x],{x,0,1},2]]
(* Output *)
{9.869612735715425,39.47894896829713}
```

Eigensystems with inhomogeneous Neumann values cannot be solved:

```wolfram
NDEigensystem[-Laplacian[u[x],{x}]+NeumannValue[1,x==1],u[x],{x,0,1},2]
(* Output *)
NDEigensystem
(* Output *)
NDEigensystem[NeumannValue[1,x==1]-u^′′[x],u[x],{x,0,1},2]
```

Eigensystems with homogeneous Neumann values can be solved:

```wolfram
First[NDEigensystem[-Laplacian[u[x],{x}]+NeumannValue[0,x==1],u[x],{x,0,1},2]]
(* Output *)
{6.365278674517519×10^-14,9.86961273571546}
```

The same result:

```wolfram
First[NDEigensystem[-Laplacian[u[x],{x}],u[x],{x,0,1},2]]
(* Output *)
{6.365278674517519×10^-14,9.86961273571546}
```

Eigensystems with inhomogeneous generalized Neumann values cannot be solved:

```wolfram
NDEigensystem[-Laplacian[u[x],{x}]+NeumannValue[u[x]+1,x==1],u[x],{x,0,1},2]
(* Output *)
NDEigensystem
(* Output *)
NDEigensystem[NeumannValue[1+u[x],x==1]-u^′′[x],u[x],{x,0,1},2]
```

Eigensystems with homogeneous generalized Neumann values can be solved:

```wolfram
NDEigensystem[-Laplacian[u[x],{x}]+NeumannValue[u[x],x==1],u[x],{x,0,1},2]
(* Output *)
{{0.7401738853453642,11.734873791907017},{InterpolatingFunction[...][x],InterpolatingFunction[...][x]}}
```

The operator and possible boundary conditions need to be stationary and linear:

```wolfram
NDEigensystem[{D[u[t,x,y],t]==(D[u[t,x,y],{x,2}]+D[u[t,x,y],{y,2}])+(.01 t)(D[u[t,x,y],{x,2}]+D[u[t,x,y],{y,2}]),DirichletCondition[u[t,x,y]==0.,True]},u,t,{x,y}∈Disk[],10]
(* Output *)
NDEigensystem
(* Output *)
NDEigensystem[{u^(1,0,0)[t,x,y]==u^(0,0,2)[t,x,y]+u^(0,2,0)[t,x,y]+0.01 t (u^(0,0,2)[t,x,y]+u^(0,2,0)[t,x,y]),DirichletCondition[u[t,x,y]==0.,True]},u,t,{x,y}∈Disk[{0,0}],10]
```

Initial conditions will be set to zero and ignored:

```wolfram
γ=1.3;c=1.1;
f[x_]=D[0.125 Erf[(x-0.5)/0.125],x];
vInit[x_]=-c*D[f[x],x]-γ f[x]/2;
First[NDEigensystem[{D[u[t,x],{t,2}]+γ D[u[t,x],{t,1}]-c^2 D[u[t,x],{x,2}]+γ u[t,x]==0,u[0,x]==f[x],
Derivative[1,0][u][0,x]==vInit[x],
DirichletCondition[u[t,x]==0,x==0],DirichletCondition[u[t,x]==0,x==1]},u,t,{x,0,1},4]]
(* Output *)
NDEigensystem
(* Output *)
NDEigensystem
(* Output *)
{-0.6499999999999928-3.58046525052479 ⅈ,-0.6499999999999928+3.58046525052479 ⅈ,-0.6500000000000045-6.97474216381078 ⅈ,-0.6500000000000045+6.97474216381078 ⅈ}
```

The same result:

```wolfram
First[NDEigensystem[{D[u[t,x],{t,2}]+γ D[u[t,x],{t,1}]-c^2 D[u[t,x],{x,2}]+γ u[t,x]==0,DirichletCondition[u[t,x]==0,x==0],DirichletCondition[u[t,x]==0,x==1]},u,t,{x,0,1},4]]
(* Output *)
{-0.6499999999999928-3.58046525052479 ⅈ,-0.6499999999999928+3.58046525052479 ⅈ,-0.6500000000000045-6.97474216381078 ⅈ,-0.6500000000000045+6.97474216381078 ⅈ}
```

The eigenvalues and eigenfunctions will be sorted by magnitude. Define a potential function:

```wolfram
dressingTransformationPotentials[evals_,xMax_]:=Module[{ε,F,f,evals1=Sort[evals]-Max[evals],W=0&},Do[
ε=evals1[[-k]];
F=NDSolveValue[{f'[x]-f[x]^2+W[x]==ε,f[0]==0},f,{x,0,xMax}];
W=Evaluate[2 ε+2 F[#]^2-W[#]]&
,{k,2,Length[evals1]}];
W]
V=dressingTransformationPotentials[{1,2,3,4,5},20];
```

Note the negative eigenvalues:

```wolfram
NDEigensystem[-ψ''[x]+V[Abs[x]] ψ[x],ψ,{x,-20,20},20][[1]]
(* Output *)
{0.006862666671842475,0.00887663973509894,0.06268891708457003,0.07952713522865552,0.177715676496226,0.21934123747403309,0.35488905779444196,0.42696029768605265,0.5958237070641151,0.7029403568046996,-0.8164637065930467,0.9023150861443376,1.0508535375797672,1.2775588960612732,1.4776408755725714,1.7248954726096934,-1.7424464845247147,1.9920679039591416,2.2252452386335912,2.8228106980406147}
```

Obtain the eigenvalues by the real part:

```wolfram
NDEigenvalues[-ψ''[x]+V[Abs[x]] ψ[x],ψ,{x,-20,20},20,Method->{"Eigensystem"->{"Arnoldi","Criteria"->"RealPart"}}]
(* Output *)
{0.006862666671843817,0.008876639735100497,0.06268891708457044,0.07952713522865654,0.17771567649622558,0.21934123747403222,0.3548890577944428,0.42696029768605315,0.5958237070641158,0.7029403568046986,0.9023150861443389,1.0508535375797667,1.2775588960612752,1.4776408755725734,1.7248954726096957,1.9920679039591422,2.225245238633591,2.822810698040617,3.1520441309736453,3.3940007352620305}
```

[NDEigensystem](https://reference.wolfram.com/language/ref/NDEigensystem.html) finds the $n$ smallest magnitude eigenvalues and their correspondent eigenfunctions for a given differential operator. Particularly in cases where one is interested in finding the most negative eigenvalues, e.g. many quantum mechanical problems, the $n$ most negative eigenvalues might not correspond with the $n$ eigenvalues [NDEigensystem](https://reference.wolfram.com/language/ref/NDEigensystem.html) returns by default. To see this, consider the following example.

Take, for instance, the dimensionless radial Schrödinger equation for the hydrogen atom, where the energy units are Rydbergs and length is measured in Bohr radii.

Define the radial Schrödinger equation:

```wolfram
H[l_]:=-ψ''[r]+((l(l+1))/(r^2)-(2)/(r))ψ[r]
```

Solve the eigenvalue problem for $l=0$. A refined mesh is used for a good approximation quality:

```wolfram
{en,funs}=NDEigensystem[{H[0],DirichletCondition[ψ[r]==0,True]},ψ[r],{r,0,100},4,Method->{"SpatialDiscretization"->{"FiniteElement",{"MeshOptions"->{"MaxCellMeasure"->0.001}}}}];
```

Look at the eigenvalues obtained in electronvolts:

```wolfram
en UnitConvert[,"Electronvolts"]
(* Output *)
{0.04506395633158937,-0.12688225963282415,0.25216093851455273,-0.2611306700368917}
```

These correspond to the eigenvalues closest to $0$. But they are not in the desired order. For instance, take a look at the analytical eigenvalues for this problem.

The analytical energies for this equation follow this equation:

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
{en,funs}=NDEigensystem[{H[0],DirichletCondition[ψ[r]==0,True]},ψ[r],{r,0,100},4,Method->{"Eigensystem"->{"Arnoldi","Shift"->λ0},"SpatialDiscretization"->{"FiniteElement",{"MeshOptions"->{"MaxCellMeasure"->0.001}}}}];
```

Look at the eigenvalues obtained in electronvolts:

```wolfram
en UnitConvert[,"Electronvolts"]
(* Output *)
{-0.8503558211282636,-1.5117436835664173,-3.401423284831847,-13.605693122766406}
```

These correspond to the most negative eigenvalues, but in reverse order, and then it is necessary to sort the eigenstates.

Sort the eigenstates with [SortBy](https://reference.wolfram.com/language/ref/SortBy.html):

```wolfram
{en,funs}=SortBy[{en,funs}ᵀ,Re]ᵀ;
```

Plot the probability densities:

```wolfram
Plot[Evaluate[funs^2],{r,0,30},PlotLegends->en UnitConvert[,"Electronvolts"],PlotRange->All]
(* Output *)
![image](img/image_042.png)
```

Take a look at the difference between the analytical eigenvalues and the obtained results:

```wolfram
Table[QuantityMagnitude[enHydrogen[n]],{n,4}]-en
(* Output *)
{-1.6745049791211386×10^-11,3.001161541504871×10^-10,2.376743801768555×10^-10,6.917133532624575×10^-11}
```

[NDEigensystem](https://reference.wolfram.com/language/ref/NDEigensystem.html) converts PDEs to time-dependent PDEs. This transformation is not unique and may lead to what seem to be unexpected results for coupled PDEs:

```wolfram
NDEigensystem[{
-u[x]-Laplacian[u[x],{x}],
-v[x]-Laplacian[v[x],{x}]},{v,u},{x}∈Line[{{0},{1}}],1]
(* Output *)
Eigensystem
(* Output *)
{{-0.027280306458415445},{{InterpolatingFunction[...],InterpolatingFunction[...]}}}
```

Internally, the given equations are rewritten as a system of time-dependent PDEs. This is done by adding time-dependent derivative terms and a time variable $t$ to the operator and setting them to 0. In the previous case, from the given dependent variables `{*v*[*x*],*u*[*x*]}`, the following temporal system is generated: `{[D](https://reference.wolfram.com/language/ref/D.html)[*v*[*t*,*x*],*t*]==- *u*[*t*,*x*]-[Laplacian](https://reference.wolfram.com/language/ref/Laplacian.html)[*u*[*t*,*x*],{*x*}],[D](https://reference.wolfram.com/language/ref/D.html)[*u*[*t*,*x*],*t*]==-*v*[*t*,*x*]-[Laplacian](https://reference.wolfram.com/language/ref/Laplacian.html)[*v*[*t*,*x*],{*x*}]}`

Because of the ordering of the dependent variables `{*v*[*x*],*u*[*x*]}`, the first equation has a time derivative term on `*v*[*t*,*x*]` and the second equation has a time derivative term on `*u*[*t*,*x*]`. This is most likely not the intended equation.

To uniquely specify the system of equations, it is best to use the temporal description:

```wolfram
NDEigensystem[{
D[u[t,x],t]==-u[t,x]-Laplacian[u[t,x],{x}],
D[v[t,x],t]==-v[t,x]-Laplacian[v[t,x],{x}]},{v,u},t,{x}∈Line[{{0},{1}}],1]
(* Output *)
{{-0.9999999999998526},{{InterpolatingFunction[...],InterpolatingFunction[...]}}}
```

As an alternative, the dependent variables can be given in the order `{*u*[*x*],*v*[*x*]}`:

```wolfram
NDEigensystem[{
-u[x]-Laplacian[u[x],{x}],
-v[x]-Laplacian[v[x],{x}]},{u,v},{x}∈Line[{{0},{1}}],1]
(* Output *)
{{-0.9999999999998526},{{InterpolatingFunction[...],InterpolatingFunction[...]}}}
```

More information on this topic can be found in [Finite Element Method Usage Tips](https://reference.wolfram.com/language/FEMDocumentation/tutorial/FiniteElementBestPractice.html#898113620).

Caution needs to be exercised when converting equations. Consider the Sturm-Liouville problem $\mathit{\nabla_{\{x \}}^{2}}\mathit{u}[\mathit{x}]=\lambda x^{2} u[x]$ with homogeneous Dirichlet conditions. The analytical solution of the eigenvalues is given by:

```wolfram
N[4*BesselJZero[1/4,1]^2]
(* Output *)
30.933346133863868
```

To avoid numerical stability problems at $x=0$, it is best to use the temporal description:

```wolfram
{vals,funs}=NDEigensystem[{x^2*D[u[t,x],{t,1}]+Laplacian[u[t,x],{x}]==0,DirichletCondition[u[t,x]==0,x==0||x==1]},{u[t,x]},t,{x,0,1},1];
vals
(* Output *)
{30.933443250231416}
```

Compare this to the following formulation of the problem:

```wolfram
{vals,funs}=NDEigensystem[{-(u''[x]/x^2),DirichletCondition[u[x]==0,x==0||x==1]},u,{x,0,1},1];
vals
(* Output *)
{16.10354581943299}
```

In some cases, [NDEigensystem](https://reference.wolfram.com/language/ref/NDEigensystem.html) may return what seem to be unexpected results for coupled PDEs:

```wolfram
eqTest={
-u[x]-Derivative[2][u][x],
-d[x]-Derivative[2][d][x]};
NDEigensystem[eqTest,{u[x],d[x]},{x}∈Line[{{0},{700}}],1]
(* Output *)
{{-0.9510204081632653},{{InterpolatingFunction[...][x],InterpolatingFunction[...][x]}}}
```

One way to avoid this issue is to specify an ordering of the dependent variables via the `"InterpolationOrder"` option:

```wolfram
NDEigensystem[eqTest,{u[x],d[x]},{x}∈Line[{{0},{700}}],1,Method->{"PDEDiscretization"->{"FiniteElement","InterpolationOrder"->{u->2,d->2}}}]
(* Output *)
{{-0.9510204081632653},{{InterpolatingFunction[...][x],InterpolatingFunction[...][x]}}}
```

Alternatively, the `"Direct"` method can be used:

```wolfram
NDEigensystem[eqTest,{u[x],d[x]},{x}∈Line[{{0},{700}}],1,Method->"Direct"]
(* Output *)
{{-0.9510204081632656},{{InterpolatingFunction[...][x],InterpolatingFunction[...][x]}}}
```

More information on this topic can be found in [Finite Element Method Usage Tips](https://reference.wolfram.com/language/FEMDocumentation/tutorial/FiniteElementBestPractice.html#898113620).

Eigenfunctions can be negative multiples of an analytical solution.

Define parameters for an eigenvalue problem with an analytical solution:

```wolfram
h=1/2; m=1;xl=-10;xh=+10;k=1; ω=1;
λ=Sqrt[h/(m ω)];
φ1[x_]:= 1/2 k (x)^2;
ψ[n_,x_]:= (Pi*λ^(2 ))^(-1/4)*1/Sqrt[2^nn!]*HermiteH[n,x/λ]*Exp[-x^2/(2λ^(2 ))];
```

Solve the eigensystem:

```wolfram
{vals1,funs1}=NDEigensystem[-h^2/(2 m)*u''[x]+φ1[x]*u[x],u[x],{x,xl,xh},30,Method->{"SpatialDiscretization"->{"FiniteElement",{"MeshOptions"->{MaxCellMeasure->0.01}}}}];
```

Compare the analytical to the numerical eigenfunction, which differ by a multiplicative $-1$:

```wolfram
GraphicsPlot[Evaluate[funs1[[1]]],{x,xl,xh},PlotRange->{{xl,xh},{-2,2}},Frame->True,GridLines->Automatic,PlotLabel->"StyleBox["NDEigensystem",FontColor->RGBColor[1, 0, 0]]"],Plot[ψ[0,x],{x,xl,xh},PlotRange->{{xl,xh},{-2,2}},Frame->True,GridLines->Automatic,PlotLabel->"StyleBox["Analytic",FontColor->RGBColor[1, 0, 0]]"]
```

*([Graphics])*

There is no overall sign guaranteed for the eigenfunctions; ($-1$) times an eigenfunction is still an eigenfunction. In fact, any constant times an eigenfunction is still an eigenfunction. [NDEigensystem](https://reference.wolfram.com/language/ref/NDEigensystem.html) chooses eigenfunctions that have a unit norm, as explained in the [properties and relations](https://reference.wolfram.com/language/ref/NDEigensystem.html#1729336659) section.

Currently, [NDEigensystem](https://reference.wolfram.com/language/ref/NDEigensystem.html) cannot deal with indexed variables:

```wolfram
NDEigensystem[-Laplacian[Indexed[u,1][x],{x}],Indexed[u,1][x],{x,0,π},4]
(* Output *)
NDEigensystem
(* Output *)
NDEigensystem[-u^′′[x],u[x],{x,0,π},4]
```

Construct a list of $n$ dependent variables:

```wolfram
n=1;
Table[Symbol[StringJoin[u,StringPadLeft[ToString[i],2,"0"]]][x],{i,n}]
(* Output *)
{u01[x]}
```

```wolfram
NDEigensystem[-Laplacian[u01[x],{x}],u01[x],{x,0,π},4]
(* Output *)
{{3.84028924710301×10^-14,1.000000844474201,4.000053838423349,9.000609362487232},{InterpolatingFunction[...][x],InterpolatingFunction[...][x],InterpolatingFunction[...][x],InterpolatingFunction[...][x]}}
```

When integrating over the eigenfunctions, [NIntegrate](https://reference.wolfram.com/language/ref/NIntegrate.html) may have a hard time when using symbolic regions:

```wolfram
{vals,funs}=NDEigensystem[{-Laplacian[u[x,y],{x,y}],DirichletCondition[u[x,y]==0,True]},u[x,y],{x,y}∈Disk[],3]
(* Output *)
{{5.783226344330992,14.68266911466738,14.682674595468772},{InterpolatingFunction[...][x,y],InterpolatingFunction[...][x,y],InterpolatingFunction[...][x,y]}}
```

```wolfram
NIntegrate[# Conjugate[#],{x,y}∈Disk[]]&/@funs
(* Output *)
NIntegrate
(* Output *)
NIntegrate
(* Output *)
NIntegrate
(* Output *)
General
(* Output *)
{0.999999258477718,0.9999993625066292,0.9999993503329363}
```

While the result is correct, it may be more efficient to use the same mesh that was used for the numerical eigensystem computation:

```wolfram
NIntegrate[# Conjugate[#],{x,y}∈Head[#]["ElementMesh"]]&/@funs
(* Output *)
{0.9999992708092124,0.9999993787100159,0.999999378708732}
```

### Neat Examples

Explore how eigenfunctions $u$ in this quantum mechanics example tend to cluster when the number of wells increases and in the infinite space limit form allowed energy bands. At the same time, the influence of the boundary conditions decreases with increasing structure size. Also note that within a band the eigenfunctions $u$ are quite different (they must be orthogonal to each other), while the squared eigenfunctions are similar. This indicates the various possibilities to form superpositions from the approximate eigenfunctions of a single well:

```wolfram
potentials = <|SquareWave->1 + SquareWave[(-1/2 + x)/2],
 TriangleWave->1 + TriangleWave[(-1/2 + x)/2],
     Cos->1 + Cos[π( x+1)],Circle->Piecewise[{{Sqrt[1-(-1+2*Abs[Mod[x,2,-1]])^2],0<Abs[Mod[x,2,-1]]<1/2},{2-Sqrt[1-(-1+2*Abs[Mod[x,2,-1]])^2],1/2<Abs[Mod[x,2,-1]]<1}},0]|>;
Manipulate[Module[{V=α potentials[P] /.   If[shift,x->x+1,{}],nds},
nds=NDEigensystem[{-u''[x]+ V u[x],
Which[bc==="Dirichlet", DirichletCondition[u[x]==0,True],
bc==="Neumann",Nothing,
bc===periodic,PeriodicBoundaryCondition[u[x],x==+L+1/2,Function[x,x-(2L+1)]],
bc===antiperiodic,PeriodicBoundaryCondition[-u[x],x==+L+1/2,Function[x,x-(2L+1)]]]}, u[x],{x,-L-1/2, +L+1/2},n,
Method->{"PDEDiscretization"->{"FiniteElement",
{"MeshOptions"->{MaxCellMeasure->ControlActive[0.01,0.002]}}}}];
Show[{Plot[V,{x,-L-1/2, +L+1/2},Exclusions->None,Filling->-1],
Plot[Evaluate[MapIndexed[(nds[[1,#2[[1]]]]+σ If[wfs,#1^2,#1])&, nds[[2]]]],{x,-L-1/2, +L+1/2},
PlotStyle->If[hl===none,Directive[Gray,Thickness[0.002]],
MapAt[Directive[Darker[Red],Thickness[0.003]]&,
Table[Directive[Gray,Thickness[0.002]],{n}],hl]],Filling->If[wfs, MapIndexed[(#2[[1]]->#1)&,nds[[1]]],None]]},
PlotRange -> All,Axes -> False,Frame -> True]],
Style[@HoldForm[-u''[x]+V[x]u[x]==λ u[x]],Bold],
{{P, SquareWave,potential, Style["V",Italic],"(",Style[x,Italic],")"}, (#1->Plot[#2,{x,-2,2},PlotTheme->"Minimal",Exclusions->{},ImageSize->40])&@@@Normal[potentials],SetterBar},
{{shift,False,"potential position"},{True->"barrier in center",False -> "well in center"}},
{{bc,"Dirichlet","boundary condition" },
{"Dirichlet","Neumann",periodic,antiperiodic}},
{{L, 5, "potential range"}, 1, 10,1,Appearance->"Labeled"},
{{α, 8.75,"potential strength"},0,30,Appearance->"Labeled"},
{{n,10,"max eigenfunction count"}, 1, 20,1,Appearance->"Labeled"},
{{σ,2,"eigenfunction scaling"}, 0.1,20},
{{wfs,False,show},{False->eigenfunctions, True->"eigenfunctions squared"}},
{{hl,none,"highlight eigenfunction"},Join[{none},(#->Style[#,8])&/@ Range[n]],SetterBar, ImageSize->Tiny},
TrackedSymbols:>True, SaveDefinitions->True,ControlPlacement->Top]
```

## Tech Notes ▪Advanced Numerical Differential Equation Solving in the Wolfram Language ▪Numerical Mathematics in the Wolfram Language ▪PDE Models Overview

## Related Guides ▪Partial Differential Equations ▪Differential Equations ▪Solvers over Regions ▪Matrix Decompositions ▪Differential Operators

## History Introduced in 2015 (10.2) | Updated in 2023 (13.3)
