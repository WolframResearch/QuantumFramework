# Det | [SpanFromLeft]

> [Det](https://reference.wolfram.com/language/ref/Det.html)[*m*] — gives the determinant of the square matrix `*m*`.

## Details and Options

[Det](https://reference.wolfram.com/language/ref/Det.html)[*m*,[Modulus](https://reference.wolfram.com/language/ref/Modulus.html)->*n*] computes the determinant modulo `*n*`.

## Examples

### Basic Examples

Find the determinant of a symbolic matrix:

```wolfram
Det[({{a_1,1, a_1,2}, {a_2,1, a_2,2}})]
(* Output *)
-a_1,2 a_2,1+a_1,1 a_2,2
```

The determinant of an exact matrix:

```wolfram
Det[{{1,2,3},{4,5,6},{7,8,9}}]
(* Output *)
0
```

### Scope

#### Basic Uses

Find the determinant of a [MachinePrecision](https://reference.wolfram.com/language/ref/MachinePrecision.html) matrix:

```wolfram
Det[{{1.7,7.1,-2.7},{2.2,8.7,3.2},{3.2,-9.2,1.2}}]
(* Output *)
251.572
```

Determinant of a complex matrix:

```wolfram
Det[{{1.+I,2,3-2 I},{0,4 π,5I},{3,0,6}}]
(* Output *)
-37.6991118430775+180.79644737231007 ⅈ
```

Determinant of an exact matrix:

```wolfram
Det[{{1,2,4},{5,4,5},{9,2,7}}]
(* Output *)
-66
```

Determinant of an arbitrary-precision matrix:

```wolfram
Det[RandomReal[2,{3, 3},WorkingPrecision->20]]
(* Output *)
4.58610424731056414521241979187962180887
```

Determinant of a symbolic matrix:

```wolfram
Det[{{a,b,c},{d,e,f},{g,h,i}}]
(* Output *)
-c e g+b f g+c d h-a f h-b d i+a e i
```

The determinant of a large numerical matrix is computed efficiently:

```wolfram
mat=BlockRandom[RandomReal[3,{1500,1500}],RandomSeeding->1234];
```

```wolfram
Det[mat]//Timing
(* Output *)
{0.19226,-2.049012094100383401426529096056327401×10^1963}
```

Note that the result may not be a machine number:

```wolfram
Divide[Abs[Last[%]],$MaxMachineNumber]
(* Output *)
1.139800811586965352140215033639564923×10^1655
```

Determinant of a matrix with finite field elements:

```wolfram
ℱ=FiniteField[17,3];
Det[{{ℱ[12],ℱ[23],ℱ[34]},{ℱ[45],ℱ[56],ℱ[67]},{ℱ[78],ℱ[89],ℱ[90]}}]
(* Output *)
<|interpretation -> FiniteFieldElement[FiniteField[17, 14, +, #, +, #, ^, 3, &, Polynomial], 3216], index -> 4661, shortIndex -> 4661, indexShortened -> True, characteristic -> 17, shortCharacteristic -> 17, extensionDegree -> 3, field -> FiniteField[...], fieldDisplayed -> False|>
```

Determinant of a [CenteredInterval](https://reference.wolfram.com/language/ref/CenteredInterval.html) matrix:

```wolfram
(m=Map[CenteredInterval,RandomReal[{-10,10},{3,3},WorkingPrecision->10],{2}])//MatrixForm
(* Output *)
({{<|Interpretation -> interpretation, Center -> 7.1564109169412404299, Radius -> 7.1564109239280471186361864966×10^-10, Type -> Real|>, <|Interpretation -> interpretation, Center -> -2.2167962766252458096, Radius -> 2.2167962777930316775609753677×10^-10, Type -> Real|>, <|Interpretation -> interpretation, Center -> -6.3833727594465017319, Radius -> 6.3833727626666503240926431317×10^-10, Type -> Real|>}, {<|Interpretation -> interpretation, Center -> -6.9442660396452993155, Radius -> 6.9442660476093376331618856057×10^-10, Type -> Real|>, <|Interpretation -> interpretation, Center -> -9.6399236260913312435, Radius -> 9.6399236336475180308980270638×10^-10, Type -> Real|>, <|Interpretation -> interpretation, Center -> 6.8045913032256066799, Radius -> 6.8045913088204956764570852101×10^-10, Type -> Real|>}, {<|Interpretation -> interpretation, Center -> 6.5127600240521132946, Radius -> 6.5127600274578201222652751312×10^-10, Type -> Real|>, <|Interpretation -> interpretation, Center -> 9.8257320572156459093, Radius -> 9.8257320682582083293254981982×10^-10, Type -> Real|>, <|Interpretation -> interpretation, Center -> 6.2193128990475088358, Radius -> 6.219312903218976451569233177×10^-10, Type -> Real|>}})
```

```wolfram
mdet=Det[m]
(* Output *)
<|Interpretation -> interpretation, Center -> -1066.7235770104452967644, Radius -> 5.01482577242×10^-7, Type -> Real|>
```

Find a random representative `*mrep*` of `*m*`:

```wolfram
ranrep[e_CenteredInterval]:=e["Center"]+RandomInteger[{-1000,1000}]/1000 e["Radius"]
(mrep=Map[ranrep,m,{2}])//MatrixForm
(* Output *)
({{(206269501055054616439)/(28823037615171174400), -(255579209855421424731)/(115292150460684697600), -(7359527725662597535723)/(1152921504606846976000)}, {-(4003096825192476941497)/(576460752303423488000), -(173657425808603988919)/(18014398509481984000), (1569031928573243750897)/(230584300921369395200)}, {(7508701085643238787543)/(1152921504606846976000), (283207444668182492129)/(28823037615171174400), (3585189792941641974197)/(576460752303423488000)}})
```

Verify that `*mdet*` contains the determinant of `*mrep*`:

```wolfram
IntervalMemberQ[mdet,Det[mrep]]
(* Output *)
True
```

#### Special Matrices

Determinants of sparse matrices:

```wolfram
SparseArray[{{1,3}->2,{2,2}->3,{3,1}->1,{4,2}->5},{4,4}]
(* Output *)
SparseArray[...]
```

```wolfram
Det[%]
(* Output *)
0
```

```wolfram
SparseArray[{{x_,y_}/;Abs[x-y]<3->1},{10,10}]
(* Output *)
SparseArray[...]
```

```wolfram
Det[%]
(* Output *)
1
```

Determinants of structured matrices:

```wolfram
SymmetrizedArray[{{1,1}->2,{1,2}->1},{2,2},Symmetric[All]]
(* Output *)
SymmetrizedArray[...]
```

```wolfram
Det[%]
(* Output *)
-1
```

```wolfram
QuantityArray[{{1,2},{3,4}},{"Meters","Seconds"}]
(* Output *)
QuantityArray[...]
```

```wolfram
Det[%]
(* Output *)
-2
```

[IdentityMatrix](https://reference.wolfram.com/language/ref/IdentityMatrix.html) always has unit determinant:

```wolfram
Det[IdentityMatrix[22]]
(* Output *)
1
```

Determinant of [HilbertMatrix](https://reference.wolfram.com/language/ref/HilbertMatrix.html):

```wolfram
HilbertMatrix[9]//Det
(* Output *)
(1)/(1028781784378569697887052962909388800000000)
```

Compute the determinant of a $10 \times 10$ matrix of univariate polynomials of degree $100$:

```wolfram
rpoly[n_]:=RandomInteger[{-2^10,2^10},{n+1}].x^Range[0,n]
SeedRandom[1234];
m=Table[rpoly[100],{10},{10}];
```

```wolfram
Det[m]//Short//AbsoluteTiming
(* Output *)
{0.2601462,3034638622673506403581743718005+<<1485>>+998816316987602849360464925610 x^1000}
```

### Options

#### Modulus

Compute a determinant using arithmetic modulo 47:

```wolfram
m = RandomInteger[46, {100, 100}];
```

```wolfram
Det[m, Modulus->47]//AbsoluteTiming
(* Output *)
{0.006647,29}
```

This is faster than computing [Mod](https://reference.wolfram.com/language/ref/Mod.html)[Det[*m*],47]:

```wolfram
Mod[Det[m], 47]//AbsoluteTiming
(* Output *)
{0.012585,29}
```

### Applications

#### Area and Volumes

Use [Det](https://reference.wolfram.com/language/ref/Det.html) to find area of a parallelogram spanned by $(1,4)$ and $(-5,2)$:

```wolfram
v={1,4};
w={-5,2};
```

Visualize the parallelogram when one vertex is at the origin:

```wolfram
p = Parallelogram[{0,0},{v,w}];
Graphics[p,Axes->True]
```

*([Graphics])*

The area is given by the absolute value of the determinant:

```wolfram
Abs[Det[{v,w}]]
(* Output *)
22
```

Compare with the result given by [Area](https://reference.wolfram.com/language/ref/Area.html):

```wolfram
Area[p]
(* Output *)
22
```

Use [Det](https://reference.wolfram.com/language/ref/Det.html) to find the volume of a parallelepiped spanned by $(1,4,3)$, $(-2,-5,2)$ and $(-1,2,-2)$:

```wolfram
{v_1,v_2,v_3}={{1,4,3},{-2,-5,2},{-1,2,-2}};
```

Visualize the parallelepiped when one vertex is at the origin:

```wolfram
p=Parallelepiped[{0,0,0},{v_1,v_2,v_3}];
Graphics3D[p,Axes->True]
```

*([Graphics3D])*

The volume is given by the absolute value of the determinant:

```wolfram
Abs[Det[{v_1,v_2,v_3}]]
(* Output *)
45
```

Compare with a direct computation using [Volume](https://reference.wolfram.com/language/ref/Volume.html):

```wolfram
Volume[p]
(* Output *)
45
```

Use [Det](https://reference.wolfram.com/language/ref/Det.html) to find hypervolume of a hyper-parallelepiped spanned by the following vectors:

```wolfram
{v1,v2,v3,v4}={{6,-6,-1,6},{10,-1,0,-7},{3,-2,-9,-3},{4,5,6,-3}}
(* Output *)
{{6,-6,-1,6},{10,-1,0,-7},{3,-2,-9,-3},{4,5,6,-3}}
```

The hypervolume is given by the absolute value of the determinant:

```wolfram
Abs[Det[{v1,v2,v3,v4}]]
(* Output *)
3580
```

Compare with the result given by [RegionMeasure](https://reference.wolfram.com/language/ref/RegionMeasure.html):

```wolfram
RegionMeasure[Parallelepiped[{0,0,0,0},{v1,v2,v3,v4}]]
(* Output *)
3580
```

The determinant itself is negative, so the $v_{i}$ are not right-handed:

```wolfram
Det[{v1,v2,v3,v4}]
(* Output *)
-3580
```

Simply reorder any two vectors, say the middle two, to produce a right-handed set:

```wolfram
Det[{v1,v3,v2,v4}]
(* Output *)
3580
```

Find the area of the image of the unit disk $\mathbb{D}$ under the linear transformation associated to the matrix $m$:

```wolfram
m={{7,-3},{5,7},{-10,4}};
m//MatrixForm
(* Output *)
({{7, -3}, {5, 7}, {-10, 4}})
```

The area of the image $f(\mathbb{D})$ is given by $\sqrt{m.m} Area[\mathbb{D}]=\pi \sqrt{m.m}$:

```wolfram
πSqrt[Det[Transpose[m].m]]
(* Output *)
10 Sqrt[122] π
```

Compare with a direct computation:

```wolfram
f𝔻=ParametricRegion[ {m.{x,y}, {x,y}∈Disk[]},{x,y}];
Area[f𝔻]
(* Output *)
10 Sqrt[122] π
```

Visualize the image $f(\mathbb{D})$:

```wolfram
Region[f𝔻,PlotTheme->"Scientific"]
```

*([Graphics3D])*

Find the volume factor $v(r,\theta)$ in the change of variables formula $\int \int f(x,y)\mathrm{d}x \mathrm{d}y=\int \int f(r,\theta) v(r,\theta)\mathrm{d}r \mathrm{d}\theta$ between Cartesian and polar coordinates. The mapping from polar to Cartesian coordinates is given by:

```wolfram
polar[r_,θ_]:={r Cos[θ],r Sin[θ]}
```

Compute the Jacobian of the mapping using [Grad](https://reference.wolfram.com/language/ref/Grad.html):

```wolfram
jac=Grad[polar[r,θ],{r,θ}]
(* Output *)
{{Cos[θ],-r Sin[θ]},{Sin[θ],r Cos[θ]}}
```

By the change of variables theorem, the volume is the determinant of the Jacobian:

```wolfram
Simplify[Abs[Det[jac]], r>0]
(* Output *)
r
```

Compare with the result given by [CoordinateChartData](https://reference.wolfram.com/language/ref/CoordinateChartData.html):

```wolfram
CoordinateChartData["Polar","VolumeFactor",{r,θ}]
(* Output *)
r
```

The same procedure will work with any coordinate system, for example, spherical coordinates:

```wolfram
spherical[r_,θ_,φ_]:={r Sin[θ]Cos[φ],r Sin[θ]Sin[φ],r Cos[θ]}
Simplify[Abs[Det[Grad[spherical[r,θ,φ],{r,θ,φ}]]], r>0&&0<=θ<=π]
(* Output *)
r^2 Sin[θ]
```

```wolfram
CoordinateChartData["Spherical","VolumeFactor",{r,θ,φ}]
(* Output *)
r^2 Sin[θ]
```

Use the change of variables theorem to compute $\mathbb{I}=\int_{D}(x^{2}+y^{2})$, where $\mathbb{D}$ is the following region:

```wolfram
𝔻=ImplicitRegion[x>0&&y>0&& 1<=x^2-y^2<=9 && 2<=x y<=4,{x,y}];
RegionPlot[𝔻]
(* Output *)
![image](img/image_001.png)
```

First, define hyperbolic coordinates $(u,v)$ as follows:

```wolfram
hyper[x_,y_]:={x^2-y^2,x y}
```

The region $\mathbb{D}$ clearly corresponds to $1 \leq u \leq 9$ and $2 \leq v \leq 4$. By the change of variables formula, $\int \int 1 \mathrm{d}u \mathrm{d}v=\int \int \{u,v \}\mathrm{d}x \mathrm{d}y$. The gradient is given by:

```wolfram
Grad[hyper[x,y],{x,y}]
(* Output *)
{{2 x,-2 y},{y,x}}
```

The determinant of the gradient is twice the function whose integral is $\mathbb{I}$:

```wolfram
Det[%]
(* Output *)
2 x^2+2 y^2
```

Hence, $\mathbb{I}$ is given by the trivial integral $\int_{2}^{4}\int_{1}^{9}\frac{1}{2}\mathrm{d}u \mathrm{d}v$:

```wolfram
∫_2^4∫_1^9(1)/(2)ⅆuⅆv
(* Output *)
8
```

Compare with a direct integration over the region:

```wolfram
∫_{x,y}∈𝔻(x^2+y^2)
(* Output *)
8
```

#### Orientation and Rotations

Determine whether the following basis for $\mathbb{R}^{3}$ is right-handed:

```wolfram
{b1,b2,b3}={{1,0,1},{0,1,1},{1,0,0}};
```

The determinant of the matrix formed by the basis is negative, so it is not right-handed:

```wolfram
Det[{b1,b2,b3}]
(* Output *)
-1
```

Determine if linear transformation corresponding to $m$ is orientation-preserving or orientation-reversing:

```wolfram
m={{-8,10,-7,-10},{2,6,-7,-9},{6,2,4,-10},{-7,-4,-3,-10}};
```

As $m>0$, the mapping is orientation-preserving:

```wolfram
Det[m]
(* Output *)
13890
```

Show that the following matrix is not a rotation matrix:

```wolfram
m={{0.969655,-0.170187,-0.0303448,-0.170187,-0.0303448},{0.170187,0.254483,0.170187,-0.0455171,0.170187},{-0.0303448,-0.170187,0.969655,-0.170187,-0.0303448},{0.170187,-0.0455171,0.170187,0.954483,0.170187},{-0.0303448,-0.170187,-0.0303448,-0.170187,0.969655}};
```

All rotation matrices have unit determinant; since $m \neq 1$, it cannot be a rotation matrix:

```wolfram
Det[m]
(* Output *)
0.3318622031923233
```

Show that the matrix $m$ is orthogonal and determine if it is a rotation matrix or includes a reflection:

```wolfram
m={{-0.969655,-0.170187,-0.0303448,-0.170187,-0.0303448},{-0.170187,0.954483,0.170187,-0.0455171,0.170187},{0.0303448,-0.170187,0.969655,-0.170187,-0.0303448},{-0.170187,-0.0455171,0.170187,0.954483,0.170187},{0.0303448,-0.170187,-0.0303448,-0.170187,0.969655}};
```

Up to the input precision, $m=m$, which shows that $m$ is orthogonal:

```wolfram
Chop[Transpose[m].m,10^-6]//MatrixForm
(* Output *)
({{0.9999996627370801, 0, 0, 0, 0}, {0, 1.00000044858841, 0, 0, 0}, {0, 0, 0.9999996627370801, 0, 0}, {0, 0, 0, 1.00000044858841, 0}, {0, 0, 0, 0, 0.9999996627370801}})
```

All orthogonal matrices have $m=\pm 1$, but rotations have $m=1$; as $m=-1$, $m$ includes a reflection:

```wolfram
Det[m]
(* Output *)
-0.9999999426937871
```

The generalization of a rotation matrix to complex vector spaces is a special unitary matrix that is unitary and has unit determinant. Show that the following matrix is a special unitary matrix:

```wolfram
u=(1)/(Sqrt[Cosh[2 Im[α]]])({{Cosh[Im[α]], ⅈ Sinh[Im[α]]}, {ⅈ Sinh[Im[α]], Cosh[Im[α]]}});
```

The matrix is unitary because $u^{\dagger}=u$:

```wolfram
ConjugateTranspose[u].u//Simplify
(* Output *)
{{1,0},{0,1}}
```

It also has unit determinant, so it is in fact an element of the special unitary group $SU(2)$:

```wolfram
Det[u]//Simplify
(* Output *)
1
```

#### Linear and Abstract Algebra

Determine the values of the parameter $s$ for which the system $2 s x+y=1$, $3 s x+6 s y=2$ has a unique solution and describe that solution. First, form the coefficient matrix $a$ and constant vector $b$:

```wolfram
a=({{2s, 1}, {3s, 6s}});b=({{1}, {2}});
```

```wolfram
b
(* Output *)
{1,2}
```

The solutions will be unique as $a \neq 0$:

```wolfram
Det[a]!=0
(* Output *)
-3 s+12 s^2≠0
```

Solving over the reals gives three open intervals separated at $s=0$ and $s=\frac{1}{4}$:

```wolfram
Reduce[%,s,Reals]
(* Output *)
s<0||0<s<(1)/(4)||s>(1)/(4)
```

Since the matrix is invertible for these values of $s$, the solution is simply $a.b$:

```wolfram
Inverse[a].b//Simplify
(* Output *)
{(2-6 s)/(3 s-12 s^2),-(1)/(3-12 s)}
```

Verify the solution in the original system of equations:

```wolfram
Simplify[2 s x+y==1&&3 s x+6 s y==2 /. Thread[{x,y}->%]]
(* Output *)
True
```

Use Cramer's rule to solve the system of equations $6 x+9 y=11$, $3 z-7 x=-12$, $5 y+9 z=-9$. First, form the coefficient matrix $a$ and constant vector $b$:

```wolfram
a = {{6,9,0},{-7,0,3},{0,5,9}};
b={11,-12,-9};
{a//MatrixForm,b//MatrixForm}
(* Output *)
{({{6, 9, 0}, {-7, 0, 3}, {0, 5, 9}}),({{11}, {-12}, {-9}})}
```

Form the three matrices $d_{j}$ where $b$ replaces the corresponding columns of $a$:

```wolfram
{dx,dy,dz}=Table[ReplacePart[a,{j_,i}:>b[[j]]],{i,3}];
{dx//MatrixForm,dy//MatrixForm, dz//MatrixForm}
(* Output *)
{({{11, 9, 0}, {-12, 0, 3}, {-9, 5, 9}}),({{6, 11, 0}, {-7, -12, 3}, {0, -9, 9}}),({{6, 9, 11}, {-7, 0, -12}, {0, 5, -9}})}
```

The entries of the solution are given by $d_{j}/a$:

```wolfram
({Det[dx],Det[dy],Det[dz]})/(Det[a])
(* Output *)
{(188)/(159),(23)/(53),-(592)/(477)}
```

Verify the result:

```wolfram
a.%==b
(* Output *)
True
```

Write a function implementing Cramer's rule for solving a linear system `*m*.*x*=*b*`:

```wolfram
crule[m_,b_]:=Module[{d=Det[m],a},
Table[a=m;a[[All,k]]=b;Det[a]/d,{k,Length[m]}]]
```

Use the function to solve a system for particular values of `*m*` and `*b*`:

```wolfram
m={{1,2,3},{1,4,9},{1,8,27}};
b={4,16,46};
```

```wolfram
x=crule[m,{4,16,46}]
(* Output *)
{-5,3,1}
```

Verify the solution:

```wolfram
m.x==b
(* Output *)
True
```

For numerical systems, `[LinearSolve](https://reference.wolfram.com/language/ref/LinearSolve.html)` is much faster and more accurate:

```wolfram
n=500;
m=RandomReal[1,{n,n}];
x=ConstantArray[1,n];
b=m.x;
```

```wolfram
AbsoluteTiming[Norm[crule[m,b]-x]]
(* Output *)
{0.962589,1.2989004755123241×10^-11}
```

```wolfram
AbsoluteTiming[Norm[LinearSolve[m,b]-x]]
(* Output *)
{0.00513,2.8405282793103536×10^-11}
```

Determine if the matrix $a$ has a nontrivial kernel (null space):

```wolfram
a={{-10,-7,-7,8,-1},{5,0,2,0,-8},{7,-4,-2,-10,1},{-5,7,10,0,3},{-8,-5,-4,1,-2}};
```

Since the determinant is nonzero, the kernel is trivial:

```wolfram
Det[a]
(* Output *)
-24415
```

Confirm the result using [NullSpace](https://reference.wolfram.com/language/ref/NullSpace.html):

```wolfram
NullSpace[a]
(* Output *)
{}
```

Determine if the mapping corresponding to the matrix $a$ is injective:

```wolfram
a= {{3,1,1}, {2,-5,3}, {1,-11,5}};
```

Since $a=0$, the mapping is not injective:

```wolfram
Det[a]
(* Output *)
0
```

Confirm the result using [FunctionInjective](https://reference.wolfram.com/language/ref/FunctionInjective.html):

```wolfram
f[x_,y_,z_]:=a.{x,y,z}
```

```wolfram
FunctionInjective[f[x,y,z],{x,y,z}]
(* Output *)
False
```

Since $a$ defines a linear function $f:\mathbb{R}^{3}\to \mathbb{R}^{3}$, the failure to be injective implies a failure to be surjective:

```wolfram
FunctionSurjective[f[x,y,z],{x,y,z}]
(* Output *)
False
```

Determine if the matrix $a$ defines an automorphism (a bijective linear map):

```wolfram
a= RandomInteger[{-10,10},{3,3}]
(* Output *)
{{4,-9,-3},{2,-6,-9},{-3,-7,-2}}
```

Since $a!=0$, the mapping is an automorphism:

```wolfram
Det[a]
(* Output *)
-387
```

Confirm the result using [FunctionBijective](https://reference.wolfram.com/language/ref/FunctionBijective.html):

```wolfram
f[x_,y_,z_]:=a.{x,y,z}
```

```wolfram
FunctionBijective[f[x,y,z],{x,y,z}]
(* Output *)
True
```

Compute the cofactor obtained from removing row `*i*` and column `*j*`:

```wolfram
cofactor[m_,{i_Integer,j_Integer}]:=(-1)^(i+j)Det[Drop[m,{i},{j}]]
```

```wolfram
cofactor[{{1,2,3,4},{5,6,7,8},{8,7,6,4},{5,3,2,1}},{3,2}]
(* Output *)
4
```

Check the result:

```wolfram
(-1)^(3+2)Det[{{1,3,4},{5,7,8},{5,2,1}}]
(* Output *)
4
```

Modular computation of a determinant:

```wolfram
a={{12,13},{14,15}};
```

Modular determinants:

```wolfram
d1=Det[a,Modulus->5]
(* Output *)
3
```

```wolfram
d2=Det[a,Modulus->7]
(* Output *)
5
```

Recover the result:

```wolfram
ChineseRemainder[{d1,d2},{5,7}]
(* Output *)
33
```

Shift the residue to be symmetric:

```wolfram
Mod[33,5 7,(-5 7+1)/2]
(* Output *)
-2
```

Confirm that the non-modular determinant was recovered:

```wolfram
Det[a]
(* Output *)
-2
```

### Properties & Relations

The determinant is the product of the eigenvalues:

```wolfram
m=RandomReal[1,{100,100}];
```

```wolfram
Det[m]
(* Output *)
2.614110724973401×10^26
```

```wolfram
Apply[Times,Eigenvalues[m]]
(* Output *)
2.614110724973328×10^26-5.7935921152×10^10 ⅈ
```

[Det](https://reference.wolfram.com/language/ref/Det.html) satisfies $a=\sum_{\sigma}^{S_{n}}sgn[\sigma]\prod_{i}^{n}a[[i,\sigma[[i]]]]$, where $S_{n}$ is all $n$-permutations and $sgn$ is [Signature](https://reference.wolfram.com/language/ref/Signature.html):

```wolfram
a = {{-4,8,-5,6,-3},{10,2,-5,-4,8},{8,8,-1,8,0},{-7,-6,-4,-5,2},{8,-7,-6,-7,10}};
```

```wolfram
s5=Permutations[Range[5]];
```

```wolfram
Det[a]==∑_{σ}^{s5}Signature[σ]∏_{i}^{5}a[[i,σ[[i]]]]
(* Output *)
True
```

[Det](https://reference.wolfram.com/language/ref/Det.html) can be computed recursively via cofactor expansion along any row:

```wolfram
m=RandomReal[1,{5,5}];
Block[{n=Length[m],i=RandomInteger[{1,Length[m]}]},
Sum[(-1)^(i+k)m[[i,k]]Det[Drop[m,{i},{k}]],{k,n}]==Det[m]]
(* Output *)
True
```

Or any column:

```wolfram
Block[{n=Length[m],j=RandomInteger[{1,Length[m]}]},
Sum[(-1)^(j+k)m[[k,j]]Det[Drop[m,{k},{j}]],{k,n}]==Det[m]]
(* Output *)
True
```

The determinant is the signed volume of the parallelepiped generated by its rows:

```wolfram
m={{1,2,3},{-1,1,0},{3,2,1}};
```

```wolfram
Det[m]
(* Output *)
-12
```

This equals the volume up to sign:

```wolfram
Volume[Parallelepiped[{0,0,0},m]]
(* Output *)
12
```

A square matrix has an inverse if and only if its determinant is nonzero:

```wolfram
m={{1,2,1},{1,0,2},{-1,2,-3}};
```

```wolfram
Det[m]
(* Output *)
0
```

```wolfram
Inverse[m]
(* Output *)
Inverse
(* Output *)
Inverse[{{1,2,1},{1,0,2},{-1,2,-3}}]
```

The determinant of a triangular matrix is the product of its diagonal elements:

```wolfram
MatrixForm[m=SparseArray[{i_,j_}/;i>=j:>RandomReal[],{5,5}]]
(* Output *)
({{0.5571369763906076, 0, 0, 0, 0}, {0.13867149076773755, 0.3528109527764214, 0, 0, 0}, {0.9472373141722348, 0.4923969645209132, 0.3156886414885063, 0, 0}, {0.6228632604607458, 0.3484836177422703, 0.4736633908080914, 0.5879115298148636, 0}, {0.11301300433309591, 0.4037711705303699, 0.13699394460188485, 0.8874752918551632, 0.3400282783339075}})
```

```wolfram
Det[m]
(* Output *)
0.012404807011684988
```

```wolfram
Apply[Times,Diagonal[m]]
(* Output *)
0.012404807011684988
```

The determinant of a matrix product is the product of the determinants:

```wolfram
{a,b}=RandomReal[1,{2,100,100}];
```

```wolfram
Det[a.b]
(* Output *)
-2.085948473695922×10^51
```

```wolfram
Det[a]Det[b]
(* Output *)
-2.0859484736910962×10^51
```

The determinant of the inverse is the reciprocal of the determinant:

```wolfram
m=({{a_1,1, a_1,2, a_1,3}, {a_2,1, a_2,2, a_2,3}, {a_3,1, a_3,2, a_3,3}});
```

```wolfram
Det[Inverse[m]]//Simplify
(* Output *)
1/(a_1,3 (-a_2,2 a_3,1+a_2,1 a_3,2)+a_1,2 (a_2,3 a_3,1-a_2,1 a_3,3)+a_1,1 (-a_2,3 a_3,2+a_2,2 a_3,3))
```

```wolfram
1/Det[m]
(* Output *)
1/(-a_1,3 a_2,2 a_3,1+a_1,2 a_2,3 a_3,1+a_1,3 a_2,1 a_3,2-a_1,1 a_2,3 a_3,2-a_1,2 a_2,1 a_3,3+a_1,1 a_2,2 a_3,3)
```

```wolfram
%%==%//Simplify
(* Output *)
True
```

A matrix and its transpose have equal determinants:

```wolfram
m=RandomReal[1,{5,5}];
Det[m]==Det[Transpose[m]]
(* Output *)
True
```

The determinant of the matrix exponential is the exponential of the trace:

```wolfram
m=RandomReal[1,{10,10},WorkingPrecision->$MachinePrecision];
```

```wolfram
Det[MatrixExp[m]]==Exp[Tr[m]]
(* Output *)
True
```

[CharacteristicPolynomial](https://reference.wolfram.com/language/ref/CharacteristicPolynomial.html)[*m*] is equal to $|m-\lambda I |$:

```wolfram
m=RandomInteger[9,{3,3}];
```

```wolfram
Det[m-λ IdentityMatrix[3]]
(* Output *)
137+3 λ+10 λ^2-λ^3
```

```wolfram
CharacteristicPolynomial[m,λ]
(* Output *)
137+3 λ+10 λ^2-λ^3
```

[Det](https://reference.wolfram.com/language/ref/Det.html)[*m*] can be computed from [LUDecomposition](https://reference.wolfram.com/language/ref/LUDecomposition.html)[*m*]:

```wolfram
m=RandomInteger[9,{5,5}];
{l,u,p,c}=LUDecomposition[m];
Det[m]==Signature[p["PermutationList"]]Times@@Diagonal[u]
(* Output *)
True
```

Consider two rectangular matrices $a$ and $b$ such that $a.b$ and $b.a$ are both square:

```wolfram
a=RandomReal[1,{5,3}];b=RandomReal[1,{3,5}];
```

Sylvester's determinant theorem states that $\text{DoubleStruckOne}+a.b=\text{DoubleStruckOne}+b.a$, where $\text{DoubleStruckOne}$ is the matching identity matrix:

```wolfram
Det[IdentityMatrix[5]+a.b]==Det[IdentityMatrix[3]+b.a]
(* Output *)
True
```

If a matrix $m$ is the [TensorProduct](https://reference.wolfram.com/language/ref/TensorProduct.html) of two vectors $u$ and $v$, then $\text{DoubleStruckOne}+m=1+u.v$:

```wolfram
{u,v}=RandomReal[1,{2,5}];
```

```wolfram
Det[IdentityMatrix[5]+u⊗v]==1+u.v
(* Output *)
True
```

This can be expressed equally in terms of [KroneckerProduct](https://reference.wolfram.com/language/ref/KroneckerProduct.html):

```wolfram
Det[IdentityMatrix[5]+KroneckerProduct[u,v]]==1+u.v
(* Output *)
True
```

This follows from Sylvester's determinant theorem for the corresponding row and column matrices:

```wolfram
c=List/@u;r={v};
```

```wolfram
{Det[IdentityMatrix[5]+c.r]==Det[IdentityMatrix[1]+r.c], r.c=={{u.v}}}
(* Output *)
{True,True}
```

### Neat Examples

Determinants of tridiagonal matrices:

```wolfram
tridiagonal[n_]:=SparseArray[{Band[{2,1}]->a,Band[{1,1}]->b,Band[{1,2}]->c},{n,n}]
```

```wolfram
tridiagonal[5]//MatrixForm
(* Output *)
({{b, c, 0, 0, 0}, {a, b, c, 0, 0}, {0, a, b, c, 0}, {0, 0, a, b, c}, {0, 0, 0, a, b}})
```

```wolfram
Table[Det[tridiagonal[n]],{n,2,12}]//TableForm
(* Output *)
{{b^2-a c}, {b^3-2 a b c}, {b^4-3 a b^2 c+a^2 c^2}, {b^5-4 a b^3 c+3 a^2 b c^2}, {b^6-5 a b^4 c+6 a^2 b^2 c^2-a^3 c^3}, {b^7-6 a b^5 c+10 a^2 b^3 c^2-4 a^3 b c^3}, {b^8-7 a b^6 c+15 a^2 b^4 c^2-10 a^3 b^2 c^3+a^4 c^4}, {b^9-8 a b^7 c+21 a^2 b^5 c^2-20 a^3 b^3 c^3+5 a^4 b c^4}, {b^10-9 a b^8 c+28 a^2 b^6 c^2-35 a^3 b^4 c^3+15 a^4 b^2 c^4-a^5 c^5}, {b^11-10 a b^9 c+36 a^2 b^7 c^2-56 a^3 b^5 c^3+35 a^4 b^3 c^4-6 a^5 b c^5}, {b^12-11 a b^10 c+45 a^2 b^8 c^2-84 a^3 b^6 c^3+70 a^4 b^4 c^4-21 a^5 b^2 c^5+a^6 c^6}}
```

A closed-form formula for these determinants is given by $(a c)^{n/2} n$:

```wolfram
%==Table[(a c)^(n/2) ChebyshevU[n,b/(2 Sqrt[a c])]//Simplify,{n,2,12}]
(* Output *)
True
```

## Tech Notes ▪Vectors and Matrices ▪Basic Matrix Operations ▪Implementation notes: Numerical and Related Functions

## Related Guides ▪Matrix Operations ▪Finite Fields ▪Linear Systems ▪Matrices and Linear Algebra ▪Finite Mathematics ▪Symbolic Vectors, Matrices and Arrays

## Related Links [NKS|Online](http://www.wolframscience.com/nks/search/?q=Det) ([A New Kind of Science](http://www.wolframscience.com/nks/))

## History Introduced in 1988 (1.0) | Updated in 2022 (13.2) ▪ 2024 (14.0)
