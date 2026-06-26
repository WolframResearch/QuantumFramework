# SchrodingerPDEComponent | [SpanFromLeft]

> [SchrodingerPDEComponent](https://reference.wolfram.com/language/ref/SchrodingerPDEComponent.html)[*vars*,*pars*]  — yields a Schrödinger PDE term with model variables `*vars*` and model parameters `*pars*`.

## Details

Generates the Schrödinger equation to be used for eigensystems and time-dependent analysis. Given the variables and parameters, a PDE operator is returned.

[SchrodingerPDEComponent](https://reference.wolfram.com/language/ref/SchrodingerPDEComponent.html) is the non-relativistic quantum mechanics analogous to Newton's second law and describes the time evolution of a wavefunction in quantum mechanics.

[SchrodingerPDEComponent](https://reference.wolfram.com/language/ref/SchrodingerPDEComponent.html) returns a sum of differential operators to be used as a part of partial differential equations:

$$
SchrodingerPDEComponent[\mathit{vars},\ldots]=\langle \mathit{boundary condition}\rangle
$$

[SchrodingerPDEComponent](https://reference.wolfram.com/language/ref/SchrodingerPDEComponent.html) can be used to model Schrödinger equations with independent variables $\mathit{x}_{\mathit{i}}\in \mathbb{R}^{n}$ in units of meter $m$, dependent variable $\Psi$ in units of $m^{-n/2}$ and time variable $t$ in units of $s$.

Stationary model variables `*vars*` are `*vars**=**{**Ψ**[**x_1**,*…*,***x*_*n***]**,**{***x*_1**,*…*,***x*_*n***}**}*`.

Time-dependent model variables `*vars*` are `*vars**=**{**Ψ**[**t**,***x*_1**,*…*,***x*_*n***]**,**t**,**{*x*_1,…,*x*_*n*}**}*`.

The [SchrodingerPDEComponent](https://reference.wolfram.com/language/ref/SchrodingerPDEComponent.html) is based on a kinetic term and a potential term:

$$
-i \hbar \frac{\partial \mathit{\Psi}(t,x)}{\partial t}+ \nabla \cdot(-\frac{\hbar^{2}}{2 m} \nabla \mathit{\Psi}(t,x)) +V(t,x) \mathit{\Psi}(t,x)
$$

$\hbar$ is the reduced Planck constant in units of $s\,J$, and $V$ is the Schrödinger potential in units of $J$.

The mass $m$ can be isotropic, orthotropic or anisotropic. $1/m$ represents a $n$ by $n$ matrix given by:

$$
\frac{1}{m}=(\begin{pmatrix}
\frac{1}{m_{11}} & \cdot s & \frac{1}{m_{1 n}} \\
\vdots & \ddots & \vdots \\
\frac{1}{m_{n 1}} & \cdot s & \frac{1}{m_{\text{n n}}}
\end{pmatrix})\mathit{  }
$$

For a magnetic field interaction, the [SchrodingerPDEComponent](https://reference.wolfram.com/language/ref/SchrodingerPDEComponent.html) is:

$$
-i \hbar \frac{\partial \mathit{\Psi}(t,x)}{\partial t}+\nabla \cdot(-\frac{\hbar^{2}}{2 m} \nabla \mathit{\Psi}(t,x)) + V(t,x) \mathit{\Psi}(t,x) +\nabla \cdot(\frac{i \hbar q }{2 m}\overset{⇀}{A}(t,x)\mathit{\Psi}(t,x))+\frac{i \hbar q }{2}\overset{⇀}{A}(t,x)\cdot \frac{1}{m}\nabla \mathit{\Psi}(t,x)+ \frac{q^{2}}{2 }\overset{⇀}{A}(t,x)\cdot(\frac{1}{m}\overset{⇀}{A}(t,x) )\mathit{\Psi}(t,x)
$$

$q$ represents the electric charge of the particle in units of coulombs $C$.

$\overset{⇀}{A}(t,x)$ is the magnetic vector potential in units of $s\,V/m$, defined such that $\overset{⇀}{B} =\nabla \times \overset{⇀}{A}$, where $\overset{⇀}{B}$ is the magnetic flux density in units of $T$.

The units of the Schrödinger PDE terms are joules $J$ times the units of the wavefunction in $J/\sqrt{m}^{n}$.

The following model parameters `*pars*` can be given:

| *parameter* | *default* | *symbol* |
| --- | --- | --- |
| `"AzimuthalQuantumNumber"` | [None](https://reference.wolfram.com/language/ref/None.html) | $l$ |
| `"MagneticVectorPotential"` | 0 | $\overset{⇀}{A}$, magnetic vector potential in $V s m^{-1}$ |
| `"Mass"` | 1 | $m$, mass in $kg$ |
| `"PlanckConstant"` | $h$ | $h$, Planck's constant in $J s$ |
| `"ReducedPlanckConstant"` | $\hbar$ | $\hbar$, reduced Planck's constant in $J s$ |
| `"SchrodingerPotential"` | 0 | $V, potential in J$ |
| `"RegionSymmetry"` | [None](https://reference.wolfram.com/language/ref/None.html) | $-$ |
| `"ParticleCharge"` | 1 | $q$, charge in $C$ |

All parameters may depend on any of $x$, $t$ and $\Psi$, as well as other dependent variables.

Parameters given with units are converted to SI base units.

Default values for Planck's and reduced Planck's constants are in SI base units.

A possible choice for the parameter `"RegionSymmetry"` is `"Axisymmetric"`.

`"Axisymmetric"` region symmetry represents a truncated cylindrical coordinate system where the cylindrical coordinates are reduced by removing the angle variable as follows:

| *dimension* | *reduction* | *general equation* |
| --- | --- | --- |
| 1D | $\{r,\theta \}\to r$ | $\nabla_{\{r \}}\cdot(-\frac{\hbar^{2}}{2 m }\nabla_{\{r \}}\Psi)+-\frac{\hbar^{2}}{2 m r}\nabla_{\{r \}}\mathit{\Psi}\mathit{ }\mathit{+}\mathit{V \Psi}$ |
| 2D | $\{r,\theta,z \}\to \{r,z \}$ | $\nabla_{\{r,z \}}\cdot(-\frac{\hbar^{2}}{2 m }\nabla_{\{r,z \}}\mathit{\Psi}\mathit{)}+(-\frac{\hbar^{2}}{2 m r},0)\cdot \nabla_{\{r,z \}}\mathit{\Psi}\mathit{+}\mathit{V \Psi}$ |

If `"RegionSymmetry"` is set to `"Axisymmetric"`, $m$ has to be isotropic.

With an `"AzimuthalQuantumNumber"` $l$, the [SchrodingerPDEComponent](https://reference.wolfram.com/language/ref/SchrodingerPDEComponent.html) is:

$$
-i \hbar \frac{\partial \Psi(t,r,z)}{\partial t}-\frac{\hbar^{2}}{2 m r  }\frac{\partial \Psi(t,r,z)}{\partial r} -\frac{\partial}{\partial r}(\frac{\hbar^{2}}{2 m  }\frac{\partial \Psi(t,r,z)}{\partial r}) -\frac{\partial}{\partial z}( \frac{\hbar^{2}}{2 m}\frac{\partial \Psi(t,r,z)}{\partial z}) +(V(t,r,z)+\frac{\hbar^{2} }{2 m r^{2}}l^{2})\Psi(t,r,z)
$$

The diffusion component affects the meaning of [NeumannValue](https://reference.wolfram.com/language/ref/NeumannValue.html).

If the [SchrodingerPDEComponent](https://reference.wolfram.com/language/ref/SchrodingerPDEComponent.html) depends on parameters $\mathit{p_{i}}$ that are specified in the association `*pars*` as `*<|*…*,**key**->**p_i*…,*p_i**->**v*_*i**,*…*|>*`, the parameters $\mathit{p_{i}}$ are replaced with $\mathit{v_{i}}$.

## Examples

### Basic Examples

Specify a PDE operator for a particle with mass $m=1$ and potential $V(x)=0$:

```wolfram
SchrodingerPDEComponent[{Ψ[x],{x}},<||>]
(* Output *)
({{-(17561922253088409/(320000000000000000000000000000000000000000000000000000000000000000000000000000000000 π^2))}}.Ψ[x])
```

The exact value is used for Planck's constant:

```wolfram
(ℏ^2)/(2)/.{ℏ->QuantityMagnitude[,"SIBase"]}
(* Output *)
17561922253088409/(320000000000000000000000000000000000000000000000000000000000000000000000000000000000 π^2)
```

Specify a PDE operator for a particle with mass $m$ and potential $V(x)$ with a reduced Planck constant $\hbar$:

```wolfram
SchrodingerPDEComponent[{Ψ[x],{x}},<|"Mass"->m,"SchrodingerPotential"->V[x],"ReducedPlanckConstant"->ℏ|>]
(* Output *)
V[x] Ψ[x]+({{-(ℏ^2)/(2 m)}}.Ψ[x])
```

Define a Hamiltonian with a reduced Planck constant of $1/10$, a harmonic potential $x^{2}$ and a mass of $1/2$:

```wolfram
ℋ=SchrodingerPDEComponent[{Ψ[x],{x}},<|"ReducedPlanckConstant"->1/10,"SchrodingerPotential"->x^2,"Mass"->(1)/(2)|>]
(* Output *)
x^2 Ψ[x]+({{-(1)/(100)}}.Ψ[x])
```

Find the 10 smallest eigenvalues and eigenfunctions on a refined mesh:

```wolfram
{vals,funs}=NDEigensystem[ℋ,Ψ[x],{x,-3,3},10,Method->{"PDEDiscretization"->{"FiniteElement",{"MeshOptions"->{MaxCellMeasure->0.05}}}}];
```

Visualize the eigenfunctions scaled by the reduced Planck constant and offset by the respective eigenvalue:

```wolfram
Show[Plot[Evaluate[1/10*funs+vals],{x,-3,3}],
Plot[x^2,{x,-3,3}],PlotRange -> -3302, AxesOrigin -> -30]
```

*([Graphics])*

### Scope

#### 1D

Define a Schrödinger PDE operator:

```wolfram
SchrodingerPDEComponent[{Ψ[x],{x}},<||>]
(* Output *)
({{-(17561922253088409/(320000000000000000000000000000000000000000000000000000000000000000000000000000000000 π^2))}}.Ψ[x])
```

Activate the PDE operator:

```wolfram
Activate[%]
(* Output *)
-((17561922253088409 Ψ^′′[x])/(320000000000000000000000000000000000000000000000000000000000000000000000000000000000 π^2))
```

Define a Schrödinger PDE operator with $h$ as the Planck constant:

```wolfram
SchrodingerPDEComponent[{Ψ[x],{x}},<|"PlanckConstant"->h|>]
(* Output *)
({{-(h^2)/(8 π^2)}}.Ψ[x])
```

Define a Schrödinger PDE operator with $\hbar$ as the reduced Planck constant:

```wolfram
SchrodingerPDEComponent[{Ψ[x],{x}},<|"ReducedPlanckConstant"->ℏ|>]
(* Output *)
({{-(ℏ^2)/(2)}}.Ψ[x])
```

Define a Schrödinger PDE operator with mass $m$:

```wolfram
SchrodingerPDEComponent[{Ψ[x],{x}},<|"Mass"->m|>]
(* Output *)
({{-(17561922253088409/(320000000000000000000000000000000000000000000000000000000000000000000000000000000000 m π^2))}}.Ψ[x])
```

Define a Schrödinger PDE operator with a potential $V$:

```wolfram
SchrodingerPDEComponent[{Ψ[x],{x}},<|"SchrodingerPotential"->V|>]
(* Output *)
V Ψ[x]+({{-(17561922253088409/(320000000000000000000000000000000000000000000000000000000000000000000000000000000000 π^2))}}.Ψ[x])
```

Define a Schrödinger PDE operator with $\hbar$ as the reduced Planck constant, mass $m$, a magnetic vector potential $\overset{⇀}{A}=(A_{x},A_{y},A_{z})$ and charge $q$:

```wolfram
SchrodingerPDEComponent[{Ψ[x],{x}},<|"ReducedPlanckConstant"->ℏ,"Mass"->m,"MagneticVectorPotential"->{A_x,A_y,A_z},"ParticleCharge"->q|>]
(* Output *)
{(ⅈ q ℏ A_x)/(2 m)}.Ψ[x]+(q^2 A_x^2 Ψ[x])/(2 m)+({{-(ℏ^2)/(2 m)}}.Ψ[x])+{(ⅈ q ℏ A_x Ψ[x])/(2 m)}
```

#### 2D

Specify a PDE operator for a particle with mass $m$, reduced Planck constant $\hbar$ and a potential $V(x,y)$ in two spatial dimensions:

```wolfram
SchrodingerPDEComponent[{Ψ[x,y],{x,y}},<|"Mass"->m,"ReducedPlanckConstant"->ℏ,"SchrodingerPotential"->V[x,y]|>]
(* Output *)
V[x,y] Ψ[x,y]+({{-(ℏ^2)/(2 m),0},{0,-(ℏ^2)/(2 m)}}.Ψ[x,y])
```

#### 3D

Specify a PDE operator for a particle with mass $m$, reduced Planck constant $\hbar$ and a potential $V(x,y,z)$ in three spatial dimensions:

```wolfram
SchrodingerPDEComponent[{Ψ[x,y,z],{x,y,z}},<|"Mass"->m,"ReducedPlanckConstant"->ℏ,"SchrodingerPotential"->V[x,y,z]|>]
(* Output *)
V[x,y,z] Ψ[x,y,z]+({{-(ℏ^2)/(2 m),0,0},{0,-(ℏ^2)/(2 m),0},{0,0,-(ℏ^2)/(2 m)}}.Ψ[x,y,z])
```

Define a PDE operator for a particle with electric charge $q$ and mass $m$ in an external magnetic field with a correspondent magnetic vector potential $\overset{⇀}{A}$:

```wolfram
SchrodingerPDEComponent[{Ψ[x,y,z],{x,y,z}},<|"ParticleCharge"->q,"Mass"->m,"ReducedPlanckConstant"->ℏ,"MagneticVectorPotential"->{A_x,A_y,A_z}|>]
(* Output *)
{(ⅈ q ℏ A_x)/(2 m),(ⅈ q ℏ A_y)/(2 m),(ⅈ q ℏ A_z)/(2 m)}.Ψ[x,y,z]+(1)/(2) q^2 ((A_x^2)/(m)+(A_y^2)/(m)+(A_z^2)/(m)) Ψ[x,y,z]+({{-(ℏ^2)/(2 m),0,0},{0,-(ℏ^2)/(2 m),0},{0,0,-(ℏ^2)/(2 m)}}.Ψ[x,y,z])+{(ⅈ q ℏ A_x Ψ[x,y,z])/(2 m),(ⅈ q ℏ A_y Ψ[x,y,z])/(2 m),(ⅈ q ℏ A_z Ψ[x,y,z])/(2 m)}
```

#### Axisymmetric

Specify an axisymmetric time-independent Schrödinger PDE term with the reduction $\{r,\theta \}\to \{r \}$, for a particle with mass $m$, reduced Planck constant $\hbar$ and a potential $V(x)=0$:

```wolfram
ℋ=SchrodingerPDEComponent[{Ψ[r],{r}},<|"RegionSymmetry"->"Axisymmetric","Mass"->m,"ReducedPlanckConstant"->ℏ|>]
(* Output *)
{-(ℏ^2)/(2 m r)}.Ψ[r]+({{-(ℏ^2)/(2 m)}}.Ψ[r])
```

Activate the PDE operator:

```wolfram
Activate[ℋ]
(* Output *)
-(ℏ^2 Ψ^′[r])/(2 m r)-(ℏ^2 Ψ^′′[r])/(2 m)
```

Specify an axisymmetric time-independent Schrödinger PDE term with the reduction $\{r,\theta,z \}\to \{r,z \}$, for a particle with mass $m$, reduced Planck constant $\hbar$ and a potential $V(x)=0$:

```wolfram
ℋ=SchrodingerPDEComponent[{Ψ[r,z],{r,z}},<|"RegionSymmetry"->"Axisymmetric","Mass"->m,"ReducedPlanckConstant"->ℏ|>]
(* Output *)
{-(ℏ^2)/(2 m r),0}.Ψ[r,z]+({{-(ℏ^2)/(2 m),0},{0,-(ℏ^2)/(2 m)}}.Ψ[r,z])
```

Activate the PDE operator:

```wolfram
Activate[ℋ]
(* Output *)
-(ℏ^2 Ψ^(0,2)[r,z])/(2 m)-(ℏ^2 Ψ^(1,0)[r,z])/(2 m r)-(ℏ^2 Ψ^(2,0)[r,z])/(2 m)
```

Define a PDE operator with a reduced Planck constant $\hbar$, a mass $m$ and an azimuthal quantum number $l$:

```wolfram
ℋ=SchrodingerPDEComponent[{Ψ[r,z],{r,z}},<|"ReducedPlanckConstant"->ℏ,"Mass"->m,"RegionSymmetry"->"Axisymmetric","AzimuthalQuantumNumber"->l|>]
(* Output *)
{-(ℏ^2)/(2 m r),0}.Ψ[r,z]+(l^2 ℏ^2 Ψ[r,z])/(2 m r^2)+({{-(ℏ^2)/(2 m),0},{0,-(ℏ^2)/(2 m)}}.Ψ[r,z])
```

Replace the parameters in the PDE with $\hbar=1$, $m=1/2$ and $l =1$:

```wolfram
ℋ/.{ℏ->1,m->(1)/(2),l->1}
(* Output *)
{-(1)/(r),0}.Ψ[r,z]+(Ψ[r,z])/(r^2)+({{-1,0},{0,-1}}.Ψ[r,z])
```

Activate the PDE operator:

```wolfram
Activate[ℋ]
(* Output *)
(l^2 ℏ^2 Ψ[r,z])/(2 m r^2)-(ℏ^2 Ψ^(0,2)[r,z])/(2 m)-(ℏ^2 Ψ^(1,0)[r,z])/(2 m r)-(ℏ^2 Ψ^(2,0)[r,z])/(2 m)
```

#### Time-Dependent

Specify a time-dependent Schrödinger PDE term for a particle with mass $m$, reduced Planck constant $\hbar$ and a potential $V(x)$:

```wolfram
SchrodingerPDEComponent[{Ψ[t,x],t,{x}},<|"Mass"->m,"ReducedPlanckConstant"->ℏ,"SchrodingerPotential"->V[x]|>]
(* Output *)
V[x] Ψ[t,x]+({{-(ℏ^2)/(2 m)}}.Ψ[t,x])-ⅈ ℏ Ψ^(1,0)[t,x]
```

#### Material Symmetry

Define a PDE operator with a potential $V$ and with an isotropic mass $m$:

```wolfram
SchrodingerPDEComponent[{Ψ[x,y,z],{x,y,z}},<|"SchrodingerPotential"->V,"Mass"->m,"PlanckConstant"->1|>]
(* Output *)
V Ψ[x,y,z]+({{-(1)/(8 m π^2),0,0},{0,-(1)/(8 m π^2),0},{0,0,-(1)/(8 m π^2)}}.Ψ[x,y,z])
```

Define a PDE operator with a potential $V$, an orthotropic mass $\{m_{x},m_{y},m_{z}\}$ and a Planck constant of 1:

```wolfram
SchrodingerPDEComponent[{Ψ[x,y,z],{x,y,z}},<|"SchrodingerPotential"->V,"Mass"->{m_x,m_y,m_z},"PlanckConstant"->1|>]
(* Output *)
V Ψ[x,y,z]+({{-(1)/(8 π^2 m_x),0,0},{0,-(1)/(8 π^2 m_y),0},{0,0,-(1)/(8 π^2 m_z)}}.Ψ[x,y,z])
```

Define a Hamiltonian with a potential $V$, an anisotropic mass $m_{\text{i j}}$ and a Planck constant of 1:

```wolfram
SchrodingerPDEComponent[{Ψ[x,y,z],{x,y,z}},<|"SchrodingerPotential"->V,"Mass"->{{m_xx,m_xy,m_xz},{m_yx,m_yy,m_yz},{m_zx,m_zy,m_zz}},"PlanckConstant"->1|>]
(* Output *)
V Ψ[x,y,z]+({{-(1)/(8 π^2 m_xx),-(1)/(8 π^2 m_xy),-(1)/(8 π^2 m_xz)},{-(1)/(8 π^2 m_yx),-(1)/(8 π^2 m_yy),-(1)/(8 π^2 m_yz)},{-(1)/(8 π^2 m_zx),-(1)/(8 π^2 m_zy),-(1)/(8 π^2 m_zz)}}.Ψ[x,y,z])
```

### Generalizations & Extensions

Numerically specify a particle with mass $m$ and potential $V(x)=0$:

```wolfram
N[SchrodingerPDEComponent[{Ψ[x],{x}},<|"Mass"->m|>]]
(* Output *)
({{-(5.560608592867591×10^-69)/(m)}}.Ψ[x])
```

### Applications

Compute the eigensystem of a time-independent Schrödinger PDE operator with a Planck constant of 1, a mass of 1 and a potential $1/(2 \pi) x^{2}$:

```wolfram
{en,funs}=NDEigensystem[{SchrodingerPDEComponent[{Ψ[x],{x}},<|"Mass"->1,"PlanckConstant"->1,"SchrodingerPotential"->1/(2π) x^2|>]},Ψ[x],{x,-3,3},4];
```

Visualize the eigenfunctions and label the eigenvalues:

```wolfram
Plot[funs,{x,-3,3},PlotRange -> All, AxesLabel -> xΨ, PlotLegends -> en]
(* Output *)
![image](img/image_001.png)
```

Find the 10 smallest eigenvalues and eigenfunctions with a reduced Planck constant $h$ and a piecewise potential that is $\frac{x^{2}}{2}$ for negative values of $x$ and $2 x^{2}$ for positive values of $x$.

Define the PDE operator:

```wolfram
h=1/10;
V[x_]:=UnitStep[x] (2 x^2) + UnitStep[-x] (x^2/2)
ℋ=SchrodingerPDEComponent[{Ψ[x],{x}},<|"SchrodingerPotential"->V[x],"ReducedPlanckConstant"->h,"Mass"->(1)/(2)|>]
(* Output *)
((1)/(2) x^2 UnitStep[-x]+2 x^2 UnitStep[x]) Ψ[x]+({{-(1)/(100)}}.Ψ[x])
```

Plot the potential:

```wolfram
Plot[V[x],{x,-3,3},AxesLabel -> xHoldForm[V[x]], PlotRange -> 05]
```

*([Graphics])*

Find the 10 smallest eigenvalues and eigenfunctions on a refined mesh:

```wolfram
{vals,funs}=NDEigensystem[ℋ,Ψ[x],{x,-3,3},10,Method->{"PDEDiscretization"->{"FiniteElement",{"MeshOptions"->{MaxCellMeasure->0.05}}}}];
```

Look at the energy eigenvalues:

```wolfram
vals
(* Output *)
{0.09597142005510859,0.2819043568351261,0.47161515234760987,0.6600901636273434,0.8484006034235657,1.037170842687354,1.2257386760345745,1.4142362738368157,1.6029085238208953,1.7915094934925397}
```

Visualize the eigenfunctions scaled by $h$ and offset by the respective eigenvalue:

```wolfram
Show[Plot[Evaluate[h*funs+vals],{x,-3,3}],
Plot[V[x],{x,-3,3}],PlotRange -> -3302, AxesOrigin -> -30]
```

*([Graphics])*

Describe a particle confined in a two-dimensional disk of radius $R$, assuming symmetry around the axis that passes through the center of the disk. Assume that the wavefunction for the particle has no polar coordinate dependence.

Set up an axisymmetric Schrödinger PDE operator, with a Planck constant equal to $1$:

```wolfram
ℋ=SchrodingerPDEComponent[{Ψ[r],{r}},<|"ReducedPlanckConstant"->1,"RegionSymmetry"->"Axisymmetric"|>]
(* Output *)
{-(1)/(2 r)}.Ψ[r]+({{-(1)/(2)}}.Ψ[r])
```

Set the radius of the disk:

```wolfram
R=1;
```

The potential energy outside the region is infinite, then the wavefunction vanishes at the boundary.

Define a Dirichlet boundary condition:

```wolfram
bc=DirichletCondition[Ψ[r]==0,r==R];
```

Solve the eigenvalue problem:

```wolfram
{en,funs}=NDEigensystem[{ℋ,bc},Ψ[r],{r,0,R},5,Method->{"SpatialDiscretization"->{"FiniteElement","MeshOptions"->{"MaxCellMeasure"->(R)/(150)}}}];
```

Visualize the probabiity densities:

```wolfram
Table[RevolutionPlot3D[funs[[i]]^2,{r,0,R},PlotLabel -> Part[en, i], Mesh -> None, PlotRange -> All, ImageSize -> Small],{i,Length[en]}]
(* Output *)
{[Graphics3D],[Graphics3D],[Graphics3D],[Graphics3D],[Graphics3D]}
```

Consider a model in which the wavefunction can be separated in the following way: $\Psi=\psi(r,z) \Phi(\phi)=\psi(r,z)e^{-i l \phi}$, where $r$, $\phi$ and $z$ are the radial, azimuthal and height coordinates in a cylindrical coordinate system, respectively. Here, $l$ is the azimuthal quantum number.

This particular example considers a particle confined in a torus with infinite potential walls. Given these considerations, it is only necessary to solve for $\psi(r,z)$ in a cross section of the torus.

Define the major radius of the torus:

```wolfram
R=4;
```

Define the region with a minor radius of $1$:

```wolfram
Ω=DiscretizeRegion[Disk[{R,0},1],"MaxCellMeasure"->0.001R];
```

Visualize the mesh as the torus cross section it represents:

```wolfram
Graphics3D[{{Opacity[0.1],Torus[{0,0,0},{R-1,R+1}]},{MeshRegion[ReplaceAll[MeshCoordinates[, Ω], PatternTest[Pattern[r, Blank], NumericQ]PatternTest[Pattern[z, Blank], NumericQ] -> 0rz], MeshCells[, Ω, All]]}},ViewPoint -> 5.3.52.]
```

*([Graphics3D])*

Define the boundary condition such that the wavefunction vanishes at the boundary of the torus:

```wolfram
bc=DirichletCondition[ψ[r,z]==0,True]
(* Output *)
DirichletCondition[ψ[r,z]==0,True]
```

Create a helper function to find six eigenstates for different values of the azimuthal quantum number with a PDE operator such that $\hbar=1$ and $m=1/2$ for a simplified model and $V=0$ for the particle inside the torus:

```wolfram
findEigenstates[AzimuthalQuantumNumber_Integer]:=Module[{en,funs,l=AzimuthalQuantumNumber},
ℋ=SchrodingerPDEComponent[{ψ[r,z],{r,z}},<|"ReducedPlanckConstant"->1,"Mass"->(1)/(2),"RegionSymmetry"->"Axisymmetric","AzimuthalQuantumNumber"->l|>];
{en,funs}=NDEigensystem[{ℋ,bc},ψ[r,z],{r,z}∈Ω,6];
{l,{en,funs}}
]
```

It is possible to explore how the energy levels vary for different values of the azimuthal quantum number. To exemplify this, calculate the energy for $l$ between 0 and 5.

Calculate the energy eigenvalues for each value of $l$ between 0 and 5:

```wolfram
energiesVsl=Transpose[Table[{l,#}&/@findEigenstates[l][[2,1]],{l,-5,5,1}]];
```

Plot the energy eigenvalues as a function of the azimuthal quantum number:

```wolfram
ListPlot[energiesVsl,PlotLabel -> En [meV] vs l, PlotTheme -> Scientific, PlotStyle -> PointSize[Medium], PlotRange -> 530]
```

*([Graphics])*

Observe how the energy levels grow nonlinearly as $|l|$ grows. Also, it can be noticed how the ground state is non-degenerate, while states two and three, as well as four and five, are two-fold degenerate.

To explore the effect of the azimuthal quantum number on the wavefunctions, plot the real and imaginary parts of the total wavefunction $\Psi(r,\phi,z)=\psi(r,z)e^{-i l \phi}$. For this purpose, you can choose any value of $l$, for instance, $l =5$.

Calculate the energy and the part of the wavefunction that depends on $r$ and $z$:

```wolfram
{enl5,funsl5}=findEigenstates[5][[2]];
```

Plot the real and imaginary parts of the wavefunction $\psi(r,z)e^{-i l \phi}$:

```wolfram
Table[DensityPlot3D[#/.{r->Sqrt[x^2+y^2],φ->ArcTan[x,y],l->5},{x,-6,6},{y,-6,6},{z,-6,6},PlotPoints -> 65, ColorFunction -> SunsetColors, PlotLabel -> Head[#], [, [_Ψ, i]SelectWithContents -> TrueSelectable -> False]]&/@{Re[ℯ^(-ⅈ l φ  ) funsl5[[i]]],Im[ℯ^(-ⅈ l φ  ) funsl5[[i]]]},{i,Length[enl5]}]
(* Output *)
{{[Graphics3D],[Graphics3D]},{[Graphics3D],[Graphics3D]},{[Graphics3D],[Graphics3D]},{[Graphics3D],[Graphics3D]},{[Graphics3D],[Graphics3D]},{[Graphics3D],[Graphics3D]}}
```

Now plot the probability densities for the same value of $l=5$:

```wolfram
Table[DensityPlot3D[Abs[ℯ^(-ⅈ l φ  ) funsl5[[i]]]^2/.{r->Sqrt[x^2+y^2],φ->ArcTan[x,y],l->5},{x,-6,6},{y,-6,6},{z,-6,6},PlotPoints -> 65, ColorFunction -> SunsetColors, PlotLabel -> [Abs[_Ψ, i], ^, 2]],{i,Length[enl5]}]
(* Output *)
{[Graphics3D],[Graphics3D],[Graphics3D],[Graphics3D],[Graphics3D],[Graphics3D]}
```

The probability density is independent of the coordinate $\phi$, as expected.

You can study a quantum billiard system, in particular, a particle confined in a 2D Bunimovich stadium.

Define the region and visualize it:

```wolfram
Ω = StadiumShape[{{-0.4,0},{0.4,0}},0.4];
Graphics[{LightBlue,Ω}]
```

*([Graphics])*

Discretize the region:

```wolfram
Ω=DiscretizeRegion[Ω,MaxCellMeasure->0.01];
```

Generate a Schrödinger operator with no potential:

```wolfram
ℋ=SchrodingerPDEComponent[{Ψ[x,y],{x,y}},<|"ReducedPlanckConstant"->1,"Mass"->(1)/(2)|>]
(* Output *)
({{-1,0},{0,-1}}.Ψ[x,y])
```

Given that the potential energy outside the region is infinite, the wavefunction vanishes at the boundary.

Define the boundary condition:

```wolfram
bc=DirichletCondition[Ψ[x,y]==0,True];
```

Find the first 100 eigenstates:

```wolfram
{en,funs}=NDEigensystem[{ℋ,bc},Ψ[x,y],{x,y}∈Ω,100,Method->{"Eigensystem"->{"Arnoldi","MaxIterations"->2000}}];
```

Visualize the energy levels:

```wolfram
ListPlot[en,AxesLabel->{n,"E_n"}]
```

*([Graphics])*

Plot the wavefunctions for different states:

```wolfram
ContourPlot[funs[[#]],{x,y}∈Ω,AspectRatio -> Automatic, ImageSize -> Small]&/@{1,7,50}
(* Output *)
![image](img/image_003.png)
```

Set up a time-dependent Schrödinger PDE operator with a reduced Planck constant:

```wolfram
ℋ=SchrodingerPDEComponent[{Ψ[t,x],t,{x}},<|"ReducedPlanckConstant"->1|>]
(* Output *)
({{-(1)/(2)}}.Ψ[t,x])-ⅈ Ψ^(1,0)[t,x]
```

Solve the equation with a Gaussian as an initial setting:

```wolfram
solution=NDSolveValue[{ℋ==0,Ψ[0,x]==Exp[-(x^2)],Ψ[t,5]==0,Ψ[t,-5]==0},Ψ,{t,0,20},{x,-5,5}]
(* Output *)
InterpolatingFunction[...]
```

Apply [Animate](https://reference.wolfram.com/language/ref/Animate.html) to the solution:

```wolfram
Animate[Plot[Abs[solution[t,x]]^2,{x,-5,5},PlotRange->{0,1}],{t,0,20,0.01},SaveDefinitions -> True, AnimationRunning -> False]
```

Set up a nonlinear, time-dependent Schrödinger PDE operator:

```wolfram
ℋ =SchrodingerPDEComponent[{Ψ[t,x],t,{x}},<|"ReducedPlanckConstant"->1,"SchrodingerPotential"->Ψ[t,x]Ψ[t,x]*|>]
(* Output *)
Conjugate[Ψ[t,x]] Ψ[t,x]^2+({{-(1)/(2)}}.Ψ[t,x])-ⅈ Ψ^(1,0)[t,x]
```

Define boundary and initial conditions:

```wolfram
bc={Ψ[t,-10]==0, Ψ[t,10]==0};
ic=Ψ[0,x]==ℯ^(-x^2);
```

Compute the solution:

```wolfram
sol=NDSolveValue[{ℋ==0,ic,bc},Ψ[t,x],{t,0,10},{x,-10,10}]
(* Output *)
InterpolatingFunction[...][t,x]
```

Plot the absolute value of the solution:

```wolfram
DensityPlot[Abs[sol],{x,-10,10},{t,0,10},PlotPoints -> 40, ColorFunction -> SunsetColors]
```

*([Graphics])*

Set up a time-dependent Schrödinger PDE operator with a harmonic potential:

```wolfram
ℋ=SchrodingerPDEComponent[{Ψ[t,x],t,{x}},<|"ReducedPlanckConstant"->1,"SchrodingerPotential"->(x^2)/(2)|>]
(* Output *)
(1)/(2) x^2 Ψ[t,x]+({{-(1)/(2)}}.Ψ[t,x])-ⅈ Ψ^(1,0)[t,x]
```

Solve the equation with a coherent state as an initial setting:

```wolfram
solution=NDSolveValue[{ℋ==0,Ψ[0,x]==((1)/(π))^((1)/(4))ℯ^(((-(x-1)^2)/(2)))},Ψ,{t,0,20},{x,-4,4}]
(* Output *)
InterpolatingFunction[...]
```

[Animate](https://reference.wolfram.com/language/ref/Animate.html) the probability density with the potential scaled down by a factor of $\frac{1}{4}$:

```wolfram
Animate[Plot[{(Norm@solution[t,x])^2,(1)/(4)(x^2)/(2)},{x,-4,4},PlotRange -> 00.6, PlotPoints -> 80, PlotLabel -> Row[Abs[, Ψ[N[t, 4], x]]^2[V[x]], , ]],{t,0,20},SaveDefinitions -> True, AnimationRunning -> False]
```

Set up a time-dependent Schrödinger PDE operator with a reduced Planck constant and a Gaussian potential $V(x)=6 e^{-x^{2}}$:

```wolfram
V[x_]:=6 ℯ^(-x^2)
ℋ=SchrodingerPDEComponent[{Ψ[t,x],t,{x}},<|"ReducedPlanckConstant"->1,"SchrodingerPotential"->V[x]|>]
(* Output *)
6 ℯ^(-x^2) Ψ[t,x]+({{-(1)/(2)}}.Ψ[t,x])-ⅈ Ψ^(1,0)[t,x]
```

Specify an initial condition of a wave packet with positive momentum:

```wolfram
ic=((1)/(2 π ))^(1/4)ℯ^(3ⅈ x-((7+x)/(2))^2);
```

Visualize it along with the potential that is being scaled down by a $1/6$ factor to fit in the plot:

```wolfram
Plot[{Re[ic],Im[ic],V[x]/6 },{x,-12,12},PlotRange -> All, AxesLabel -> xNone, PlotLegends -> Re[, Ψ[0, x]]Im[, Ψ[0, x]]V[x]]
(* Output *)
![image](img/image_005.png)
```

Solve the equation with a refined mesh:

```wolfram
solution=NDSolveValue[{ℋ== 0,Ψ[0,x]==ic,DirichletCondition[Ψ[t,x]==0,x==15||x==-15]},Ψ,{t,0,5},{x,-15,15},Method->{"MethodOfLines","SpatialDiscretization"->{"FiniteElement","MeshOptions"->{"MaxCellMeasure"->(1)/(30)}}}]
(* Output *)
InterpolatingFunction[]
```

Visualize the solution where the vertical and horizontal axes represent time, $t$, and position, $x$, respectively:

```wolfram
DensityPlot[Abs[solution[t,x]],{x,-15,15},{t,0,5},PlotPoints -> 100, Mesh -> False, ColorFunction -> SunsetColors, AxesLabel -> xt, Frame -> False, Axes -> True, AxesOrigin -> -150]
(* Output *)
![image](img/image_007.png)
```

Note that at about $t \approx 4.5$, the wave packed hits the wall and gets reflected.

Visualize the solution with an animation:

```wolfram
Animate[Plot[{Abs[solution[t,x]],(1)/(6) 6 ℯ^(-x^2)},{x,-15,15},PlotPoints -> 50, PlotLabel -> Row[Abs[, Ψ[t, x]][V[x]], , ], PlotRange -> All, AxesLabel -> xNone],{t,0,5,0.0001},SaveDefinitions -> True, AnimationRunning -> False]
```

Define an operator for a free particle:

```wolfram
ℋ=SchrodingerPDEComponent[{Ψ[t,x,y],t,{x,y}},<|"ReducedPlanckConstant"->1|>]
(* Output *)
({{-(1)/(2),0},{0,-(1)/(2)}}.Ψ[t,x,y])-ⅈ Ψ^(1,0,0)[t,x,y]
```

Define a region that defines the double slit and visualize it:

```wolfram
Ω=RegionDifference[Rectangle[{-18,-18},{18,18}],RegionUnion[Rectangle[{-18,0},{-slitWidth-(separation)/(2),thickness}],Rectangle[{-(separation)/(2),0},{(separation)/(2),thickness}], Rectangle[{(separation)/(2)+slitWidth,0},{18,thickness}]]]/.{slitWidth->2,
separation->2,
thickness->0.25};
RegionPlot[Ω]
```

*([Graphics])*

Set up the initial conditions as a wave packet:

```wolfram
ic=(ℯ^(2 ⅈ y+(1)/(18) (-x^2-(10+y)^2)))/(Sqrt[3] π^(1/4));
```

Visualize the initial condition:

```wolfram
GraphicsDensityPlot[Re[ic],x-1818, y-1818, PlotRange -> All, PlotPoints -> 60, PlotLabel -> Real part],DensityPlot[Im[ic],x-1818, y-1818, PlotRange -> All, PlotPoints -> 60, PlotLabel -> Imaginary part]
(* Output *)
![image](img/image_009.png)
```

Solve the PDE:

```wolfram
solution=NDSolveValue[{ℋ== 0,Ψ[0,x,y]==ic,DirichletCondition[Ψ[t,x,y]==0,True]},Ψ,{t,0,8},{x,y}∈Ω,Method->{"MethodOfLines","SpatialDiscretization"->{"FiniteElement","MeshOptions"->{"MaxCellMeasure"->0.1}}}]
(* Output *)
InterpolatingFunction[]
```

Visualize the solution:

```wolfram
ListAnimate[Table[DensityPlot[Abs[solution[t,x,y]],{x,y}∈ Ω,PlotPoints -> 85, ColorFunction -> SunsetColors, PlotRange -> All],{t,0,8,0.1} ],AnimationRunning -> False, SaveDefinitions -> True]
```

You can model a spinless charged particle in a magnetic field that is pointing in the $z$ direction. The particle is described by an initial wave packet, and the goal is to study the evolution of the probability density under magnetic field interaction. Since this is a charged particle in a magnetic field, one can expect for its probability density to reflect some kind of circular trajectory, similar to the behavior of classical charged particles. For that reason, a very high magnetic flux density has been chosen, such that the cyclotron radius is in the order of angstroms.

Set up the reduced Planck constant, particle's mass, magnetic flux density and the particle's electric charge, respectively:

```wolfram
ℏ=UnitConvert[,"SIBase"];
m=UnitConvert[electron[mass],"SIBase"];
B=UnitConvert[15000,"SIBase"];
e=UnitConvert[1,"SIBase"];
```

Define a constant $fs$ to measure the time in femtoseconds:

```wolfram
fs= UnitConvert[,"SIBase"]//QuantityMagnitude ;
```

The particle's initial condition will be described by a wave packet with a momentum $\hbar k$, meaning it would move in a straight line if there were no magnetic field interaction.

Define $k$ with units of $/Å$:

```wolfram
k=Quantity[1,"Angstroms"^-1];
```

A magnetic flux density is chosen to have a cyclotron radius in the order of angstroms.

Compute the cyclotron radius of the particle:

```wolfram
N[(ℏ k)/(e B)]
(* Output *)
4.388079713006044
```

The magnetic vector potential can be chosen such that $\overset{⇀}{A}=\frac{1}{2}\overset{⇀}{B}\times \overset{⇀}{r}$:

```wolfram
A=(1)/(2) {0,0,B[[1]]}⨯{x,y,z}
(* Output *)
{-7500 y,7500 x,0}
```

The work region can be defined as a rectangle of $50 \times 50$ $Å$.

Define the region:

```wolfram
Ω=Rectangle[{-25,-25},{25,25}];
```

Since the length units are Angstroms ($Å$), you need to take that into account when defining the PDE operator. For that reason, `"ScaleUnits" -> {"Meters" -> "Angstroms"}` are used in the parameters argument of the [SchrodingerPDEComponent](https://reference.wolfram.com/language/ref/SchrodingerPDEComponent.html).

Define a Schrödinger PDE operator that assumes Coulomb's gauge for the magnetic vector potential:

```wolfram
ℋ=SchrodingerPDEComponent[{Ψ[t,x,y],t,{x,y}},<|"ReducedPlanckConstant"->ℏ,"Mass"->m,"ParticleCharge"->-e,"MagneticVectorPotential"->A,"ScaleUnits"->{"Meters"->"Angstroms"}|>]
(* Output *)
{6.9555076118046073965 ⅈ y,-6.9555076118046073965 ⅈ x}.Ψ[t,x,y]+1.283484983267784978×10^-38 (6.1749512669689715726807527386225×10^37 x^2+6.1749512669689715726807527386225×10^37 y^2) Ψ[t,x,y]+({{-61.0426436900378294492,0},{0,-61.0426436900378294492}}.Ψ[t,x,y])+{6.9555076118046073965 ⅈ y Ψ[t,x,y],-6.9555076118046073965 ⅈ x Ψ[t,x,y]}-1.054571817646156391262428×10^-14 ⅈ Ψ^(1,0,0)[t,x,y]
```

Set up an initial condition:

```wolfram
ic=Ψ[0,x,y]==(ℯ^(-ⅈ k[[1]] x+(1)/(18) (-x^2-(-8+y)^2)))/(Sqrt[3] π^(1/4))
(* Output *)
Ψ[0,x,y]==(ℯ^(-ⅈ x+(1)/(18) (-x^2-(-8+y)^2)))/(Sqrt[3] π^(1/4))
```

Change the units of $\hbar$ to be used in the boundary condition:

```wolfram
ℏ=UnitConvert[ℏ,QuantityUnit[ℏ]/.{"Meters"->"Angstroms"}]
(* Output *)
(132521403)/(4000000000000000000000 π)
```

Solve the problem with a transparent boundary condition with a [NeumannValue](https://reference.wolfram.com/language/ref/NeumannValue.html):

```wolfram
Monitor[solution=NDSolveValue[{ℋ==NeumannValue[QuantityMagnitude[-ⅈ (ℏ )/(2 m)]Ψ[t,x,y],True],ic},Ψ,{t,0,18fs},{x,y}∈Ω,Method->{"MethodOfLines","SpatialDiscretization"->{"FiniteElement","MeshOptions"->{"MaxCellMeasure"->0.4}}},EvaluationMonitor:>(monitor="t = ",CForm[t])];,monitor]
```

Animate the solution:

```wolfram
animation=Table[DensityPlot[Abs[solution[t,x,y]]^2,{x,-25,25},{y,-25,25},ColorFunction->"SunsetColors", PlotPoints->50,PlotRange->{0,0.7}],{t,0,18fs,0.1 fs}];
```

```wolfram
ListAnimate[animation, AnimationRunning->False,SaveDefinitions->True]
```

As a demonstration of the Aharonov-Bohm effect, the scattering process of a wave packet that is passing through two slits under the influence of a magnetic vector potential is modeled. Notice that the magnetic flux density is zero for the entire solving domain, but the magnetic vector potential is not.

For instance, a domain is defined with two slits and a small hole that represents a long solenoid producing a magnetic flux density inside it, and in turn, producing a magnetic vector potential field throughout the rest of the domain. Note that the wavefunction cannot interact with a region in which the magnetic flux density is nonzero.

Set up the domain:

```wolfram
Ω=[Graphics];
```

Define variables and parameters such that the mass is non-isotropic, with values of the effective masses of heavy holes in $GaAs$, a common semiconductor material. Here, the `"MagneticVectorPotential"` parameter is for a solenoid with a magnetic flux density of $T$ and a radius of $nm$, transformed into Cartesian coordinates.

Define variables and parameters:

```wolfram
vars={Ψ[t,x,y],t,{x,y}};
pars=<|"Mass"->{0.34,0.7}electron[mass],"ReducedPlanckConstant"->,"ParticleCharge"->electron["electric charge"],"MagneticVectorPotential"->TransformedField[ "Polar"->"Cartesian",{0,(1000)/(2 r)}, {r,φ}->{x,y}],"ScaleUnits"->{"Meters"->"Nanometers"}|>;
fs=QuantityMagnitude[,"SIBase"];
```

Set up the PDE:

```wolfram
ℋ=SchrodingerPDEComponent[vars,pars] ;
```

Solve the PDE:

```wolfram
solution=NDSolveValue[{ℋ==0,Ψ[0,x,y]==ℯ^(-2 ⅈ y -(1)/(20)(x^2+(y-10)^2)),DirichletCondition[Ψ[t,x,y]==0,True]},Ψ,{t,0,70fs},{x,y}∈Ω,Method->{"PDEDiscretization"->{"MethodOfLines","SpatialDiscretization"->{"FiniteElement","MeshOptions"->{"MaxCellMeasure"->0.1}}}}]
(* Output *)
InterpolatingFunction[]
```

Animate the evolution of the probability density:

```wolfram
animation=Table[GraphicsDensityPlot[Abs[solution[t,x,y]]^2,{x,y}∈Ω,PlotPoints -> 60, PlotRange -> 04, Epilog -> DashedWhiteLine[-25-1025-10], ColorFunction -> SunsetColors, AspectRatio -> Automatic],Plot[Abs[solution[t,x,-10]]^2,{x,-25,25},PlotRange -> 01, AxesLabel -> xAbs[, Ψ[N[t, 3], x, -10]]^2, PlotPoints -> 80, PlotStyle -> RGBColor[1, 0.5, 0]],{t,0,70fs,2 fs}];
ListAnimate[animation, AnimationRunning->False,SaveDefinitions->True]
```

Notice how the interference pattern is shifted compared with the [case of no magnetic field](No magnetic field two slit scattering), and is non-symmetric.

### Properties & Relations

Specify an operator for a particle with mass $m$ and potential $V(x)$ with a reduced Planck constant $\hbar$:

```wolfram
SchrodingerPDEComponent[{Ψ[x],{x}},<|"Mass"->m,"ReducedPlanckConstant"->ℏ,"SchrodingerPotential"->V[x]|>]
(* Output *)
V[x] Ψ[x]+({{-(ℏ^2)/(2 m)}}.Ψ[x])
```

Activate the component:

```wolfram
Activate[%]
(* Output *)
V[x] Ψ[x]-(ℏ^2 Ψ^′′[x])/(2 m)
```

When parameters are given with units, they are automatically converted to SI base units.

Specify a PDE operator with the electron's mass and the reduced Planck constant $\hbar$:

```wolfram
SchrodingerPDEComponent[{Ψ[x],{x}},<|"Mass"->electron[mass],"ReducedPlanckConstant"->|>]
(* Output *)
({{-6.10426436900378294492744587242310891607×10^-39}}.Ψ[x])
```

Inspect the difference between each unit system in the diffusion coefficient. Compute the Planck constant divided by the electron mass:

```wolfram
-(^2)/(2 electron[mass])
(* Output *)
-0.00097847559924742474923858992902935013
```

Convert the result to SI base units:

```wolfram
UnitConvert[%,"SIBase"]
(* Output *)
-6.10426436900378294492147880806883147682×10^-39
```

This example shows that at an internal interface, the [BenDaniel-Duke](https://reference.wolfram.com/language/PDEModels/tutorial/Physics/ModelCollection/QuantumRing.html#735128210) boundary conditions are applied automatically.

Consider a quantum well with width $l_{w}$​, a finite potential $V_{0}$​ and wavefunction $\psi$. The particle's mass within the well is $m_{w}$​, while the particle's mass in the barrier region is $m_{b}$​. The mass thus is a function of the position $m(x)$. The quantities $\psi$ and $1/m(x) \partial \psi/\partial x$, should be continuous at interfaces between the barrier and the well region, and the BenDaniel-Duke boundary conditions should be applied at the interface. In fact, the BenDaniel-Duke boundary conditions are applied automatically at the interface.

Define the parameters:

```wolfram
mw=0.067 QuantityMagnitude[electron[mass],"SIBase"];
mb=0.092 QuantityMagnitude[electron[mass],"SIBase"];
lw=QuantityMagnitude[100];
V0=QuantityMagnitude[250.65,"SIBase"];
ℏ=QuantityMagnitude[,"SIBase"];
```

Note that the domain's length unit of the well, $l_{w}$ is in angstroms. To define the PDE operator, the parameter `"ScaleUnits"->{"Meters"->"Angstroms"}` is used. This converts the rest of the model quantities such that all length units are now in angstroms. This [improves the quality of the numerical solution](https://reference.wolfram.com/language/PDEModels/tutorial/PDEModelsBestPractice.html#294118184).

Define the model parameters:

```wolfram
pars=<|"ScaleUnits"->{"Meters"->"Angstroms"},"Mass"->Quantity[Piecewise[{{mb,Abs[x]>=(lw)/(2)}},mw],"Kilograms"],"SchrodingerPotential"->Quantity[Piecewise[{{V0,Abs[x]>=(lw)/(2)}},0],"Joules"],"ReducedPlanckConstant"->Quantity[ℏ ,"Joules""Seconds"]|>;
```

Define the PDE operator:

```wolfram
ℋ=SchrodingerPDEComponent[{ψ[x],{x}},pars];
```

Physically, the wavefunction should decay until it goes to $0$ for all bound states at infinity. A practical modeling approach is to use a Dirichlet boundary condition with $\psi=0$ at the external boundary, which is heuristically chosen to be placed at a long enough distance from the center of the quantum well. In this case, a distance of $4 l_{w}$ works well for this purpose.​

Define the Dirichlet condition:

```wolfram
dc=DirichletCondition[ψ[x]==0,x==-4lw||x==4lw];
```

Solve the eigenvalue problem with [NDEigensystem](https://reference.wolfram.com/language/ref/NDEigensystem.html) on a refined mesh:

```wolfram
{en,funs}=NDEigensystem[{ℋ, dc},ψ[x],{x,-4lw,4lw},3,Method->{"PDEDiscretization"->{"FiniteElement","MeshOptions"->{"MaxCellMeasure"->0.01}}}];
```

Plot the wavefunctions:

```wolfram
Show[Plot[Evaluate[funs],{x,-2lw,2lw},PlotRange -> All, PlotLegends -> Table[Ket[n], nLength[en]]]]
(* Output *)
![image](img/image_011.png)
```

The units of the eigenvalues are $kg\,Å^{2}/s^{2}$, given that the meters units were scaled in the parameters to angstroms. Knowing this,the eigenvalues can then be transformed to the desired energy unit.

Redefine the energy eigenvalues with the correct units:

```wolfram
en=UnitConvert[Quantity[en,"Kilograms"("Angstroms"^2)/("Seconds"^2)],"Millielectronvolts"]
(* Output *)
{30.5598203765541,120.10615028412593,244.36411574139}
```

Finally, the energy eigenvalues and probability densities can be compared between the numerical and analytical approaches.

Define the code for the analytical eigenvalues:

```wolfram
analyticEnergies=Module[...]
(* Output *)
{30.559820378431155,120.10615027083826,244.3613818059211}
```

Define the analytic wavefunctions code:

```wolfram
analyticFuns[x_,n_Integer,E0_]:=analyticFuns[...]
```

Plot the numerical and analytical probability densities together:

```wolfram
Column[{GraphicsRow[
Table[Plot[{funs[[i]]^2,analyticFuns[x ,i,analyticEnergies[[i]]]^2
},{x,-lw,lw},PlotRange -> All, AxesLabel -> x[, Å]ToString[Subscript[, ψ, ^, 2, ToString[i]], FormatType -> ], PlotStyle -> BlueOrangeDashed, ImageSize -> Small],{i,3}]
],LineLegend[BlueOrangeDashed, NumericalAnalytical, LegendLayout -> Row]},Alignment -> Center]
(* Output *)
{{[Graphics]}, {"Numerical"}}
```

Compare the difference between the analytical solution and the numerical one:

```wolfram
en-analyticEnergies
(* Output *)
{-1.877054955912172×10^-9,1.3287674960338336×10^-8,0.002733935468910431}
```

The above has shown that the energy eigenvalues provided by [NDEigensystem](https://reference.wolfram.com/language/ref/NDEigensystem.html) are numerically the same as the analytical solution to the Schrödinger equation when BenDaniel-Duke boundary conditions are implemented. Furthermore, the wavefunctions obtained by [NDEigensystem](https://reference.wolfram.com/language/ref/NDEigensystem.html) are equivalent to the analytical solution. This is the case because the numerical procedure ensures the continuity of both the wavefunction $\psi$ and the quantity $1/m(x) \partial \psi/\partial x$.

### Possible Issues

The mass has to be isotropic for the `"Axisymmetric"` case:

```wolfram
SchrodingerPDEComponent[{Ψ[r,z],{r,z}},<|"SchrodingerPotential"->V,"Mass"->{{m_rr,m_rz},{m_zr,m_zz}},"ReducedPlanckConstant"->ℏ,"RegionSymmetry"->"Axisymmetric","AzimuthalQuantumNumber"->l|>]
(* Output *)
SchrodingerPDEComponent[{Ψ[r,z],{r,z}},<|"SchrodingerPotential"->V,"Mass"->{{m_rr,m_rz},{m_zr,m_zz}},"ReducedPlanckConstant"->ℏ,"RegionSymmetry"->"Axisymmetric","AzimuthalQuantumNumber"->l|>]
```

### Neat Examples

Define a quantum well potential energy as a piecewise function:

```wolfram
V[z_]:=Piecewise[{{V00,Abs[z]>=d0/2}}]
```

Set up the variables and parameters, with the mass and the reduced Planck constant equal to $1$, for a simplified Schrödinger PDE operator:

```wolfram
ℋ=SchrodingerPDEComponent[{Ψ[z],{z}},<|"SchrodingerPotential"->V[z],"Mass"->1,"ReducedPlanckConstant"->1|>]
(* Output *)
({, {{V00, Abs[z]>=(d0)/(2)}, {0, True}}}) Ψ[z]+({{-(1)/(2)}}.Ψ[z])
```

Use [NDEigensystem](https://reference.wolfram.com/language/ref/NDEigensystem.html) to solve the time-independent Schrödinger equation in one dimension. The objective of this visualization is to see the effect that the barrier height and the well's length have on the eigenstates. One can see how increasing the well's length decreases the energy for each wavefunction. On the other hand, increasing the barrier height results in higher energy eigenvalues.

Define a function to calculate the eigenstates for any given value of potential barrier, well's width or number of eigenstates:

```wolfram
solver[V0_,d_,n_]:=NDEigensystem[{With[{V00=V0,d0=d},ℋ//Evaluate],DirichletCondition[Ψ[z]==0,True]},Ψ[z],{z,-4d,4d},n,Method->{"PDEDiscretization"->{"FiniteElement",{"MeshOptions"->{"MaxCellMeasure"->0.01 d}}}}];
```

Set up a [Manipulate](https://reference.wolfram.com/language/ref/Manipulate.html):

```wolfram
Manipulate[
Block[{en,funs},
{en,funs}=solver[V0,d,n];
ListPlot[en->en,Filling -> Bottom, PlotLabel -> Energy eigenstates, ImageSize -> Small, PlotRange -> 01.],Show[Plot[With[{V00=V0,d0=d},V[z]//Evaluate],{z,-2d,2d},PlotRange -> All, AxesLabel -> zV0,Ψ, Exclusions -> None],Plot[Evaluate[0.25 funs+en],{z,-2 d,2 d},PlotRange -> All],ImageSize -> Small]
],{{V0,0.75},0.5,1},{{d,10},5,25},{{n,4},{2,3,4,5}},SaveDefinitions -> True]
```

## Tech Notes ▪Quantum ring

## Related Guides ▪Partial Differential Equation Terms

## History Introduced in 2024 (14.0) | Updated in 2024 (14.1)
