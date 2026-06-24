# LinearSolve

> [LinearSolve](https://reference.wolfram.com/language/ref/LinearSolve.html)[*m*,*b*] — finds an `*x*` that solves the matrix equation `*m*.*x*==*b*`.
> [LinearSolve](https://reference.wolfram.com/language/ref/LinearSolve.html)[*m*] — generates a [LinearSolveFunction](https://reference.wolfram.com/language/ref/LinearSolveFunction.html)[…] that can be applied repeatedly to different `*b*`.
> [LinearSolve](https://reference.wolfram.com/language/ref/LinearSolve.html)[*a,b*] — finds an `*x*` that solves the array equation `*a*.*x*==*b*`.

## Details and Options

[LinearSolve](https://reference.wolfram.com/language/ref/LinearSolve.html) works on both numerical and symbolic matrices, as well as [SparseArray](https://reference.wolfram.com/language/ref/SparseArray.html) objects.

The argument `*b*` can be either a vector or a matrix.

The matrix `*m*` can be square or rectangular.

[LinearSolve](https://reference.wolfram.com/language/ref/LinearSolve.html)[*m*] and [LinearSolveFunction](https://reference.wolfram.com/language/ref/LinearSolveFunction.html)[…] provide an efficient way to solve the same approximate numerical linear system many times.

[LinearSolve](https://reference.wolfram.com/language/ref/LinearSolve.html)[*m*,*b*] is equivalent to [LinearSolve](https://reference.wolfram.com/language/ref/LinearSolve.html)[*m*][*b*].

For underdetermined systems, [LinearSolve](https://reference.wolfram.com/language/ref/LinearSolve.html) will return one of the possible solutions; [Solve](https://reference.wolfram.com/language/ref/Solve.html) will return a general solution.

For an `*n*_1×…×`*n*_k`×`*m*` array `*a*` and an `*n*_1×…×`*n*_*k*`×`*d*_1×…×`*d*_*l*` array `*b*`, [LinearSolve](https://reference.wolfram.com/language/ref/LinearSolve.html)[*a*,*b*] gives an `*m*`×`*d*_1×…×`*d*_*l*` array `*x*`, such that `*a*.*x*==*b*`.

[LinearSolve](https://reference.wolfram.com/language/ref/LinearSolve.html) has the following options and settings:

| [Method](https://reference.wolfram.com/language/ref/Method.html) | [Automatic](https://reference.wolfram.com/language/ref/Automatic.html) | method to use |
| --- | --- | --- |
| [Modulus](https://reference.wolfram.com/language/ref/Modulus.html) | 0 | prime modulus to use |
| [ZeroTest](https://reference.wolfram.com/language/ref/ZeroTest.html) | [Automatic](https://reference.wolfram.com/language/ref/Automatic.html) | test to determine when expressions are zero |

[LinearSolve](https://reference.wolfram.com/language/ref/LinearSolve.html)[*m*,…[Modulus](https://reference.wolfram.com/language/ref/Modulus.html)->*p*] solves rational systems modulo the prime `*p*`. If `*p*` is zero, ordinary arithmetic is used.

The [ZeroTest](https://reference.wolfram.com/language/ref/ZeroTest.html) option only applies to exact and symbolic matrices.

With [Method](https://reference.wolfram.com/language/ref/Method.html)->[Automatic](https://reference.wolfram.com/language/ref/Automatic.html), the method is automatically selected depending upon input.

Explicit [Method](https://reference.wolfram.com/language/ref/Method.html) settings for exact and symbolic matrices include:

"CofactorExpansion" | Laplace cofactor expansion
"DivisionFreeRowReduction" | Bareiss method of division-free row reduction
"OneStepRowReduction" | standard row reduction

Explicit [Method](https://reference.wolfram.com/language/ref/Method.html) settings for approximate numeric matrices include:

|  | "Banded" | banded matrix solver |
| --- | --- | --- |
|  | "Cholesky" | Cholesky method for positive definite Hermitian matrices |
|  | "Krylov" | iterative Krylov sparse solver |
|  | "Multifrontal" | direct sparse LU decomposition |
| ▪ | "Mumps" | parallel direct sparse solver |
|  | "Pardiso" | parallel direct sparse solver |

## Examples

### Basic Examples

Solve the matrix-vector equation $m.x=b$ with $m=\begin{pmatrix}
r & s \\
t & u
\end{pmatrix}$ and $b=\begin{pmatrix}
y \\
z
\end{pmatrix}$:

```wolfram
LinearSolve[{{r,s},{t,u}},{y,z}]
(* Output *)
{(u y-s z)/(-s t+r u),(t y-r z)/(s t-r u)}
```

Verify the solution:

```wolfram
{{r,s},{t,u}}.%=={y,z}//Simplify
(* Output *)
True
```

Solve the matrix equation $m.x=b$ with $m=\begin{pmatrix}
1 & 2 \\
3 & 4
\end{pmatrix}$ and $b=\begin{pmatrix}
5 & 6 \\
7 & 8
\end{pmatrix}$:

```wolfram
m={{1,2},{3,4}};
b={{5,6},{7,8}};
LinearSolve[m,b]//MatrixForm
(* Output *)
({{-3, -4}, {4, 5}})
```

Verify the solution:

```wolfram
m.%==b
(* Output *)
True
```

Solve a rectangular matrix equation:

```wolfram
m=({{1, 5}, {2, 6}, {3, 7}, {4, 8}});b=({{9}, {10}, {11}, {12}});
```

```wolfram
LinearSolve[m,b]
(* Output *)
{-1,2}
```

Verify the solution:

```wolfram
m.%==b
(* Output *)
True
```

### Scope

#### Basic Uses

Solve $m.x=b$ at machine precision:

```wolfram
m=({{1.2, 3.2, 3.2}, {7.9, -1.4, 5.1}, {1.1, 2.5, -1.5}});
```

```wolfram
LinearSolve[m,{.5,1.3,-2.2}]
(* Output *)
{-0.30970149253731344,-0.36268656716417924,0.6350746268656717}
```

Solve a case where $b$ is a matrix:

```wolfram
LinearSolve[m,{{1.5,4.75,-3.2},{-1.7,6.7,-9.3},{4.9,-8.65,15.4}}]//MatrixForm
(* Output *)
({{0.5770994425463046, -1.301452076964575, 2.116166157165977}, {1.1609242941916924, -1.0649433555115992, 2.595468440927891}, {-0.9085865851465565, 3.0373628843733145, -4.389030749865132}})
```

Solve $m.x=b$ for a complex matrix:

```wolfram
LinearSolve[{{1+I,2,3-2 I},{0,4,5I},{1+I,6,3+3I}},{1,2,3}]
(* Output *)
{0,(1)/(2),0}
```

Find a solution for an exact, rectangular matrix:

```wolfram
(m={{1,2,3,4},{5,6,7,8},{9,10,11,12}})//MatrixForm
(* Output *)
({{1, 2, 3, 4}, {5, 6, 7, 8}, {9, 10, 11, 12}})
```

```wolfram
LinearSolve[m,{13,14,15}]
(* Output *)
{-(25)/(2),(51)/(4),0,0}
```

Compute a solution at arbitrary precision:

```wolfram
(m=RandomReal[4,{2, 3},WorkingPrecision->20])//MatrixForm
(* Output *)
({{3.68615328219987110582771608346952518787, 2.36790123934776455848823933925118012667, 1.38088439047715436089401423724876138976}, {1.48933102010743352006084082306269777973, 0.50491076107417303962311173892718585421, 0.1407954757441195365881125151075181634}})
```

```wolfram
LinearSolve[m,{Pi,E}]
(* Output *)
{2.91243481928634601801773938909437805672,-3.20709681152392225887691563182937270082,0`20.}
```

Solve $m.x=b$ for a symbolic matrix:

```wolfram
LinearSolve[{{a,b,c},{d,e,f},{g,h,i}},{3,2,1}]
(* Output *)
{(c e-b f-2 c h+3 f h+2 b i-3 e i)/(c e g-b f g-c d h+a f h+b d i-a e i),(c d-a f-2 c g+3 f g+2 a i-3 d i)/(-c e g+b f g+c d h-a f h-b d i+a e i),(b d-a e-2 b g+3 e g+2 a h-3 d h)/(c e g-b f g-c d h+a f h+b d i-a e i)}
```

Solve the system when $b$ is a matrix:

```wolfram
LinearSolve[{{a,b,c},{d,e,f},{g,h,i}},DiagonalMatrix[{1,0,-1}]]//MatrixForm
(* Output *)
({{(f h-e i)/(c e g-b f g-c d h+a f h+b d i-a e i), 0, (-c e+b f)/(c e g-b f g-c d h+a f h+b d i-a e i)}, {(f g-d i)/(-c e g+b f g+c d h-a f h-b d i+a e i), 0, (-c d+a f)/(-c e g+b f g+c d h-a f h-b d i+a e i)}, {(e g-d h)/(c e g-b f g-c d h+a f h+b d i-a e i), 0, (-b d+a e)/(c e g-b f g-c d h+a f h+b d i-a e i)}})
```

Solve $m.x=b$ over a finite field:

```wolfram
ℱ=FiniteField[29,4];
m={{ℱ[12],ℱ[23],ℱ[34]},{ℱ[45],ℱ[56],ℱ[67]},{ℱ[78],ℱ[89],ℱ[90]}};
b={{ℱ[123],ℱ[234]},{ℱ[345],ℱ[456]},{ℱ[567],ℱ[678]}};
LinearSolve[m,b]//MatrixForm
(* Output *)
![image](img/image_001.png)
```

```wolfram
m.%===b
(* Output *)
True
```

Solve $m.x=b$ for [CenteredInterval](https://reference.wolfram.com/language/ref/CenteredInterval.html) matrices:

```wolfram
m=Map[CenteredInterval,RandomReal[{-10,10},{3,3},WorkingPrecision->10],{2}];
b=Map[CenteredInterval,RandomReal[{-10,10},{3,2},WorkingPrecision->10],{2}];
(sol=LinearSolve[m,b])//MatrixForm
(* Output *)
({{<|Interpretation -> interpretation, Center -> -0.3776445852557515082, Radius -> 9.9581047574570646929714712314×10^-10, Type -> Real|>, <|Interpretation -> interpretation, Center -> 0.9203669134110441519, Radius -> 1.47437663119975859160604159115×10^-9, Type -> Real|>}, {<|Interpretation -> interpretation, Center -> -0.1662569332450232196, Radius -> 9.4154247246014399763680557953×10^-10, Type -> Real|>, <|Interpretation -> interpretation, Center -> -1.1146872568878052334, Radius -> 1.35284616842312743756338022649×10^-9, Type -> Real|>}, {<|Interpretation -> interpretation, Center -> -1.0721459805605491056, Radius -> 8.1106371098427221255633412511×10^-10, Type -> Real|>, <|Interpretation -> interpretation, Center -> -0.7937628841620494313, Radius -> 8.8020001703936756598523061257×10^-10, Type -> Real|>}})
```

Find random representatives `*mrep*` and `*brep*` of `*m*` and `*b*`:

```wolfram
ranrep[e_CenteredInterval]:=e["Center"]+RandomInteger[{-1000,1000}]/1000 e["Radius"]
(mrep=Map[ranrep,m,{2}])//MatrixForm
(* Output *)
({{-(672980458564084734751)/(72057594037927936000), -(2413497498866728391)/(281474976710656000), (670241489912381454361)/(2305843009213693952000)}, {-(3854015692704537686367)/(576460752303423488000), -(293422105772630194481)/(36028797018963968000), (2859375052400080104329)/(288230376151711744000)}, {(3002767587093399590663)/(576460752303423488000), (527630747089909040631)/(57646075230342348800), (5428595781884011830961)/(2305843009213693952000)}})
```

```wolfram
(brep=Map[ranrep,b,{2}])//MatrixForm
(* Output *)
({{(334414015954834426303)/(72057594037927936000), (42160918386492610281)/(57646075230342348800)}, {-(3895351429657908445429)/(576460752303423488000), -(178328624310696070181)/(36028797018963968000)}, {-(1733131702738263520209)/(288230376151711744000), -(4195039226539635706087)/(576460752303423488000)}})
```

Verify that `*sol*` contains [LinearSolve](https://reference.wolfram.com/language/ref/LinearSolve.html)[*mrep*,*brep*]:

```wolfram
MapThread[IntervalMemberQ,{sol,LinearSolve[mrep,brep]},2]//MatrixForm
(* Output *)
({{True, True}, {True, True}, {True, True}})
```

Solve $m.x=b$ for $x$ when $b$ is a matrix of different dimensions:

```wolfram
m={{1,1,1},{1,2,3},{1,4,9}};
b={{1,2},{3,4},{5,6}};
LinearSolve[m,b]
(* Output *)
{{-2,-1},{4,4},{-1,-1}}
```

When no right[Hyphen]hand side $b$ for $m.x=b$ is given, a [LinearSolveFunction](https://reference.wolfram.com/language/ref/LinearSolveFunction.html) is returned:

```wolfram
m={{1,1,1},{1,2,3},{1,4,9}};
ls=LinearSolve[m]
(* Output *)
LinearSolveFunction[...]
```

This contains data to solve the problem quickly for a few values of $b$:

```wolfram
Map[ls,{{1,2,3},{0,1,3},{1,0,2}}]
(* Output *)
{{-(1)/(2),2,-(1)/(2)},{-1,1,0},{4,-5,2}}
```

#### Special Matrices

Solve $m.x=b$  with sparse matrices:

```wolfram
m=SparseArray[{{1,2}->1,{3,1}->5,{2,3}->2},{3,3}]
(* Output *)
SparseArray[...]
```

```wolfram
b=SparseArray[{1,0,2}]
(* Output *)
SparseArray[...]
```

As the result is typically not sparse, the result is returned as an ordinary list:

```wolfram
LinearSolve[m,b]
(* Output *)
{(2)/(5),1,0}
```

Sparse methods are used to efficiently solve sparse matrices:

```wolfram
s=SparseArray[{{i_,i_}->-2.,{i_,j_}/;Abs[i-j]==1->1.},{999,999}]
(* Output *)
SparseArray[...]
```

```wolfram
b=SparseArray[{499->1},{999}]
(* Output *)
SparseArray[...]
```

```wolfram
AbsoluteTiming[x=LinearSolve[s,b];]
(* Output *)
{0.000895,Null}
```

Visualize the result:

```wolfram
ListPlot[x]
```

*([Graphics])*

Solve a system with structured matrices:

```wolfram
m=SymmetrizedArray[{{1,2}->3,{2,2}->5,{2,3}->7,{3,2}->13},{3,3},Symmetric[All]]
(* Output *)
SymmetrizedArray[...]
```

```wolfram
b=SymmetrizedArray[{{2}->1},{3}]
(* Output *)
SymmetrizedArray[...]
```

```wolfram
LinearSolve[m,b]
(* Output *)
{(1)/(3),0,0}
```

Use a different type of matrix structure:

```wolfram
m=QuantityArray[{{1,2},{3,4}},{"Meters","Seconds"}]
(* Output *)
QuantityArray[...]
```

```wolfram
b=QuantityArray[{5,6},"Meters"]
(* Output *)
QuantityArray[...]
```

```wolfram
LinearSolve[m,b]
(* Output *)
{-4,(9)/(2)}
```

An identity matrix always produces a trivial solution:

```wolfram
LinearSolve[IdentityMatrix[3],{a,b,c}]
(* Output *)
{a,b,c}
```

Solve a linear system whose coefficient matrix is a Hilbert matrix:

```wolfram
LinearSolve[HilbertMatrix[4],{a,b,c,d}]
(* Output *)
{4 (4 a-30 b+60 c-35 d),-60 (2 a-20 b+45 c-28 d),60 (4 a-45 b+108 c-70 d),-140 (a-12 b+30 c-20 d)}
```

Solve a $10 \times 10$ system whose coefficients are univariate polynomials of degree $100$:

```wolfram
rpoly[n_]:=RandomInteger[{-2^10,2^10},{n+1}].x^Range[0,n]
SeedRandom[1234];
m=Table[rpoly[100],{10},{10}];
b=Table[rpoly[100],{10}];
```

```wolfram
LinearSolve[m,b]//Short[#,3]&//AbsoluteTiming
(* Output *)
{0.7114798,{(5625322287210941851316952361428+<<1497>>+2578979863996329764873509807266 x^1000)/(3034638622673506403581743718005+<<1485>>+998816316987602849360464925610 x^1000),(<<1>>)/(<<1>>),<<6>>,(<<1>>)/(<<1>>),(<<1>>)/(<<1>>)}}
```

#### Arrays

Solve $a.x=b$ with a 2×3×6 array $a$ and a 2×3×4×5 array $b$:

```wolfram
a=RandomInteger[{-10,10},{2,3,6}];
```

```wolfram
b=RandomInteger[{-10,10},{2,3,4,5}];
```

The result is a 6×4×5 array:

```wolfram
(x=LinearSolve[a,b])//Dimensions
(* Output *)
{6,4,5}
```

Verify the solution:

```wolfram
a.x==b
(* Output *)
True
```

### Options

#### Method

"Banded"  (1)

Solve $m.x=b$ using a banded matrix method:

```wolfram
m=SparseArray[{{i_,i_}->-2.,{i_,j_}/;Abs[i-j]==1->1.},{19999,19999}];
b=Table[0.,{19999}];b[[9999]]=1.;
```

```wolfram
Timing[x=LinearSolve[m,b,Method->"Banded"];]
(* Output *)
{0.,Null}
```

Check a relative error of the computed solution:

```wolfram
Norm[m.x-b]/Norm[x]
(* Output *)
1.3669827509714134×10^-16
```

"Cholesky"  (1)

Solve $m.x=b$ using the Cholesky decomposition:

```wolfram
m=SparseArray[{{i_,i_}->2.,{i_,j_}/;Abs[i-j]==1->-1.},{19999,19999}];
b=Table[0.,{19999}];b[[9999]]=1.;
```

```wolfram
Timing[x=LinearSolve[m,b,Method->"Cholesky"];]
(* Output *)
{0.0625,Null}
```

Check a relative error of the computed solution:

```wolfram
Norm[m.x-b]/Norm[x]
(* Output *)
2.898012744861011×10^-16
```

"Krylov"  (2)

The following suboptions can be specified for the method `"Krylov"`:

"BasisSize" | the size of the Krylov basis (GMRES only)
"MaxIterations" | the maximum number of iterations
"Method" | methods to be used
"Preconditioner" | which preconditioner to apply
"PreconditionerSide" | how to apply a preconditioner (`"Left"` or `"Right"`)
"ResidualNormFunction" | A norm function that computes a norm of the residual of the solution
"StartingVector" | the initial vector to start iterations
"Tolerance" | the tolerance used to terminate iterations

Possible settings for `"Method"` include:

"BiCGSTAB" | iterative method for arbitrary square matrices
"ConjugateGradient" | iterative method for Hermitian positive definite matrices
"GMRES" | iterative method for arbitrary square matrices

Possible settings for `"Preconditioner"` include:

"ILU0" | a preconditioner based on an incomplete LU factorization of the original matrix without fill-in
"ILUT" | a variant of ILU0 with fill-in
"ILUTP" | a variant of ILUT with column permutation

Possible suboptions for `"Preconditioner"` include:

"FillIn" | upper bound on the number of additional nonzero elements in a row introduced by the ILUT preconditioner
"PermutationTolerance" | when to permute columns
"Tolerance" | drop tolerance (any element of magnitude smaller than this tolerance is treated as zero)

Solve $m.x=b$ using a Krylov method:

```wolfram
m=SparseArray[{{i_,i_}->-2.,{i_,j_}/;Abs[i-j]==1->1.},{19999,19999}];
b=Table[0.,{19999}];b[[9999]]=1.;
```

```wolfram
AbsoluteTiming[x=LinearSolve[m,b,Method->"Krylov"];]
(* Output *)
{0.016112,Null}
```

Check a relative error of the computed solution:

```wolfram
Norm[m.x-b]/Norm[x]
(* Output *)
1.366982750971415×10^-16
```

"Multifrontal"  (1)

Solve $m.x=b$ using a direct multifrontal method:

```wolfram
m=SparseArray[{{i_,i_}->-2.,{i_,j_}/;Abs[i-j]==1->1.},{19999,19999}];
b=Table[0.,{19999}];b[[9999]]=1.;
```

```wolfram
AbsoluteTiming[x=LinearSolve[m,b,Method->"Multifrontal"];]
(* Output *)
{0.059621,Null}
```

Check a relative error of the computed solution:

```wolfram
Norm[m.x-b]/Norm[x]
(* Output *)
1.3685767137209222×10^-16
```

"Pardiso"  (1)

Solve $m.x=b$ using Pardiso:

```wolfram
m=SparseArray[{{i_,i_}->-2.,{i_,j_}/;Abs[i-j]==1->1.},{19999,19999}];
b=Table[0.,{19999}];b[[9999]]=1.;
```

```wolfram
AbsoluteTiming[x=LinearSolve[m,b,Method->"Pardiso"];]
(* Output *)
{0.039811,Null}
```

Check a relative error of the computed solution:

```wolfram
Norm[m.x-b]/Norm[x]
(* Output *)
1.632657682200382×10^-16
```

#### Modulus

Find the solution `x` to `m.x==b` modulo 47:

```wolfram
m={{1,1,1},{1,2,3},{1,4,9}};
b={1,2,3};
LinearSolve[m,b,Modulus->47]
(* Output *)
{23,2,23}
```

Verify the solution:

```wolfram
Mod[m.%,47]
(* Output *)
{1,2,3}
```

### Applications

#### Spans and Linear Independence

The following three vectors are not linearly independent:

```wolfram
v1={1,2,3};v2={4,5,6};v3={7,8,9};
Solve[a v1 + b v2 == v3,{a,b}]
(* Output *)
{{a->-1,b->2}}
```

The equation with a generic right-hand side does not have a solution:

```wolfram
LinearSolve[{v1,v2,v3},{x,y,z}]
(* Output *)
LinearSolve
(* Output *)
LinearSolve[{{1,2,3},{4,5,6},{7,8,9}},{x,y,z}]
```

Equivalently, the equation with the identity matrix on the right-hand side has no solution:

```wolfram
LinearSolve[{v1,v2,v3},IdentityMatrix[3]]
(* Output *)
LinearSolve
(* Output *)
LinearSolve[{{1,2,3},{4,5,6},{7,8,9}},{{1,0,0},{0,1,0},{0,0,1}}]
```

The following three vectors are linearly independent:

```wolfram
v1={1,2,3};v2={4,5,6};v3={7,8,10};
Solve[a v1 + b v2 == v3,{a,b}]
(* Output *)
{}
```

The equation with a generic right-hand side has a solution:

```wolfram
LinearSolve[{v1,v2,v3},{x,y,z}]
(* Output *)
{(1)/(3) (-2 x-4 y+3 z),(1)/(3) (-2 x+11 y-6 z),x-2 y+z}
```

Equivalently, the equation with the identity matrix on the right-hand side has a solution:

```wolfram
LinearSolve[{v1,v2,v3},IdentityMatrix[3]]
(* Output *)
{{-(2)/(3),-(4)/(3),1},{-(2)/(3),(11)/(3),-2},{1,-2,1}}
```

The solution is the inverse of $\{v_{1},v_{2},v_{3}\}$:

```wolfram
%==Inverse[{v1,v2,v3}]
(* Output *)
True
```

Determine if the following vectors are linearly independent or not:

```wolfram
v_1={108,90,252,186};
v_2={240,260,520,420};
v_3={264,340,536,468};
v_4={522,705,1038,929};
```

As $\{v_{1},v_{2},v_{3},v_{4}\}.x$ does not have a solution for an arbitrary $b$, they are not linearly independent:

```wolfram
LinearSolve[{v_1,v_2,v_3,v_4},{w,x,y,z}]
(* Output *)
LinearSolve
(* Output *)
LinearSolve[{{108,90,252,186},{240,260,520,420},{264,340,536,468},{522,705,1038,929}},{w,x,y,z}]
```

#### Equation Solving and Invertibility

Solve the following system of equations:

```wolfram
system={2x+3y+5z==-2, -3x+4y-z==7,4x+y-6z==3};
```

Rewrite the system in matrix form:

```wolfram
a={{2,3,5},{-3,4,-1},{4,1,-6}};
b={-2,7,3};
system ==Thread[ a.{x,y,z}==b]
(* Output *)
True
```

Use [LinearSolve](https://reference.wolfram.com/language/ref/LinearSolve.html) to find a solution:

```wolfram
LinearSolve[a,b]
(* Output *)
{-(2)/(3),(73)/(69),-(53)/(69)}
```

Show that the solution is unique using [NullSpace](https://reference.wolfram.com/language/ref/NullSpace.html):

```wolfram
NullSpace[a]=={}
(* Output *)
True
```

Verify the result using [SolveValues](https://reference.wolfram.com/language/ref/SolveValues.html):

```wolfram
SolveValues[system,{x,y,z}]
(* Output *)
{{-(2)/(3),(73)/(69),-(53)/(69)}}
```

Find all solutions of the following system of equations:

```wolfram
system={10 x+12 y+14 z==1,22 x+27 y+32 z==0,34 x+42 y+50 z==-1}
(* Output *)
{10 x+12 y+14 z==1,22 x+27 y+32 z==0,34 x+42 y+50 z==-1}
```

First, write the coefficient matrix $m$, variable vector $v$ and constant vector $b$:

```wolfram
m={{10,12,14},{22,27,32},{34,42,50}};
b={1,0,-1};
v={x,y,z};
```

Verify the rewrite:

```wolfram
system==Thread[m.v==b]
(* Output *)
True
```

[LinearSolve](https://reference.wolfram.com/language/ref/LinearSolve.html) gives a particular solution:

```wolfram
x1=LinearSolve[m,b]
(* Output *)
{(9)/(2),-(11)/(3),0}
```

[NullSpace](https://reference.wolfram.com/language/ref/NullSpace.html) gives a basis for solutions to the homogeneous equation $m.b=0$:

```wolfram
ns=NullSpace[m]
(* Output *)
{{1,-2,1}}
```

Define $x_{0}$ to be an arbitrary linear combination of the elements of $ns$:

```wolfram
x0=Sum[C[i]ns[[i]],{i,Length[ns]}]
(* Output *)
{1,-2 1,1}
```

The general solution is the sum of $x_{1}$ and $x_{0}$:

```wolfram
x1+x0
(* Output *)
{(9)/(2)+1,-(11)/(3)-2 1,1}
```

```wolfram
m.%==b//Simplify
(* Output *)
True
```

Determine if the following matrix has an inverse:

```wolfram
a=({{108, 90, 252, 186}, {240, 260, 520, 420}, {264, 340, 536, 468}, {522, 705, 1038, 929}});
```

Since the system $a.x=I$ has no solution, $a$ does not have an inverse:

```wolfram
LinearSolve[a,IdentityMatrix[4]]
(* Output *)
LinearSolve
(* Output *)
LinearSolve[{{108,90,252,186},{240,260,520,420},{264,340,536,468},{522,705,1038,929}},{{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}}]
```

Verify the result using [Inverse](https://reference.wolfram.com/language/ref/Inverse.html):

```wolfram
Inverse[a]
(* Output *)
Inverse
(* Output *)
Inverse[{{108,90,252,186},{240,260,520,420},{264,340,536,468},{522,705,1038,929}}]
```

Determine if the following matrix has a nonzero determinant:

```wolfram
a=({{1, 4, 2, -9}, {4, 12, 2, 5}, {6, 7, -11, 9}, {5, 15, 10, 12}});
```

Since the system $a.x=I$ has a solution, $a$'s determinant must be nonzero:

```wolfram
LinearSolve[a,IdentityMatrix[4]]
(* Output *)
{{(1387)/(3395),-(2901)/(3395),(1052)/(3395),(292)/(679)},{-(367)/(3395),(1311)/(3395),-(342)/(3395),-(113)/(679)},{(253)/(3395),-(654)/(3395),(23)/(3395),(89)/(679)},{-(66)/(679),(23)/(679),-(6)/(679),(2)/(679)}}
```

Confirm the result using [Det](https://reference.wolfram.com/language/ref/Det.html):

```wolfram
Det[a]
(* Output *)
-3395
```

Find the inverse of the following matrix:

```wolfram
m={{a,b,c},{d,e,f},{g,h,i}};
```

To find the inverse, first solve the system $m.x=I$:

```wolfram
LinearSolve[m,IdentityMatrix[3]]//MatrixForm
(* Output *)
({{(f h-e i)/(c e g-b f g-c d h+a f h+b d i-a e i), (c h-b i)/(-c e g+b f g+c d h-a f h-b d i+a e i), (c e-b f)/(c e g-b f g-c d h+a f h+b d i-a e i)}, {(f g-d i)/(-c e g+b f g+c d h-a f h-b d i+a e i), (c g-a i)/(c e g-b f g-c d h+a f h+b d i-a e i), (c d-a f)/(-c e g+b f g+c d h-a f h-b d i+a e i)}, {(e g-d h)/(c e g-b f g-c d h+a f h+b d i-a e i), (b g-a h)/(-c e g+b f g+c d h-a f h-b d i+a e i), (b d-a e)/(c e g-b f g-c d h+a f h+b d i-a e i)}})
```

Verify the result using [Inverse](https://reference.wolfram.com/language/ref/Inverse.html):

```wolfram
Simplify[%==Inverse[m]]
(* Output *)
True
```

Solve the system $m.x=b$, with several different $b$ by means of computing a [LinearSolveFunction](https://reference.wolfram.com/language/ref/LinearSolveFunction.html):

```wolfram
m=RandomReal[1,{500,500}];
```

```wolfram
AbsoluteTiming[
lsf=LinearSolve[m];
x1=lsf[UnitVector[500,1]];
x250=lsf[UnitVector[500,250]];
x500=lsf[UnitVector[500,500]];]
(* Output *)
{0.005914,Null}
```

Perform the computation by inverting the matrix and multiplying by the inverse:

```wolfram
AbsoluteTiming[inv=Inverse[m];
y1=inv.UnitVector[500,1];
y250=inv.UnitVector[500,250];
y500=inv.UnitVector[500,500];]
(* Output *)
{0.026855,Null}
```

The results are practically identical, even though [LinearSolveFunction](https://reference.wolfram.com/language/ref/LinearSolveFunction.html) is multiple times faster:

```wolfram
Norm[{x1-y1,x250-y250,x500-y500}]
(* Output *)
1.4285339682116405×10^-13
```

#### Calculus

Newton's method for finding a root of a multivariate function:

```wolfram
f[{x_,y_,z_}]={x-Sin[y],z-Cos[y],Cos[x]Sin[z]};
J[{x_,y_,z_}]=Grad[f[{x,y,z}],{x,y,z}];
FixedPoint[(#-LinearSolve[J[#],f[#]])&,{.5,.5,.5}]
(* Output *)
{1.,1.5707963267948966,0.}
```

Compare with the answer found by [FindRoot](https://reference.wolfram.com/language/ref/FindRoot.html):

```wolfram
x/.FindRoot[f[x],{x, {.5,.5,.5}}]
(* Output *)
{1.,1.5707963267948966,-6.468659422812297×10^-29}
```

Approximately solve the boundary value problem $u_{\text{x x}} + u = f; u(0) = u(1) = 0$ using discrete differences:

```wolfram
f[x_]:=Exp[-100(x-1/2)^2];
n=100;h=1./n;grid=N[Range[n-1]]/n;
s=SparseArray[{{i_,i_}->1.-2/h^2,{i_,j_}/;Abs[i-j]==1->1/h^2},{n-1,n-1}];
b=f[grid];
ux=LinearSolve[s,b];
ListPlot[Transpose[{grid,ux}]]
```

*([Graphics])*

Show the error compared with the exact solution:

```wolfram
ex=DSolveValue[{u''[x]+u[x]==f[x],u[0]==u[1]==0},u,x][grid];
ListPlot[Transpose[{grid,Abs[ux-ex]}]]
```

*([Graphics])*

### Properties & Relations

For an invertible matrix $m$, [LinearSolve](https://reference.wolfram.com/language/ref/LinearSolve.html)[*m*,*b*] gives the same result as [SolveValues](https://reference.wolfram.com/language/ref/SolveValues.html) for the corresponding system of equations:

```wolfram
m={{1,1,1},{1,0,1},{0,1,1}};
b={6,4,5};
```

```wolfram
sol=LinearSolve[m,b]
(* Output *)
{1,2,3}
```

Create the corresponding system of linear equations:

```wolfram
eqns=Thread[m.{x,y,z}=={6,4,5}]
(* Output *)
{x+y+z==6,x+z==4,y+z==5}
```

Confirm that [SolveValues](https://reference.wolfram.com/language/ref/SolveValues.html) gives the same result:

```wolfram
sol==First[SolveValues[eqns,{x,y,z}]]
(* Output *)
True
```

[LinearSolve](https://reference.wolfram.com/language/ref/LinearSolve.html) always returns the trivial solution $x=0$ to the homogenous equation $m.x=0$:

```wolfram
m={{1,2,3},{4,5,6},{7,8,9}};
```

```wolfram
LinearSolve[m,{0,0,0}]
(* Output *)
{0,0,0}
```

Use [NullSpace](https://reference.wolfram.com/language/ref/NullSpace.html) to get the complete spanning set of solutions if $m$ is singular:

```wolfram
NullSpace[m]
(* Output *)
{{1,-2,1}}
```

Compare with the result of [SolveValues](https://reference.wolfram.com/language/ref/SolveValues.html):

```wolfram
SolveValues[{{1,2,3},{4,5,6},{7,8,9}}.{x,y,z}=={0,0,0},{x,y,z}]//Quiet
(* Output *)
{{x,-2 x,x}}
```

If $m$ is nonsingular, the solution $x$ of $m.x = b$ is the inverse of $m$ when $b$ is the identity matrix:

```wolfram
m={{1,1,1},{1,2,3},{1,4,9}};
LinearSolve[m,IdentityMatrix[3]]==Inverse[m]
(* Output *)
True
```

In this case there is no solution to $m.x=b$:

```wolfram
m={{1,2,3},{4,5,6},{7,8,9}};b={1,-2,1};
```

```wolfram
LinearSolve[m,b]
(* Output *)
LinearSolve
(* Output *)
LinearSolve[{{1,2,3},{4,5,6},{7,8,9}},{1,-2,1}]
```

Use [LeastSquares](https://reference.wolfram.com/language/ref/LeastSquares.html) to minimize $‖m.x-b‖$:

```wolfram
x=LeastSquares[m,b]
(* Output *)
{0,0,0}
```

```wolfram
Norm[m.x-b]
(* Output *)
Sqrt[6]
```

Compare to general minimization:

```wolfram
Minimize[Norm[m.{x1,x2,x3}-b],{x1,x2,x3}]
(* Output *)
{Sqrt[6],{x1->0,x2->0,x3->0}}
```

If $m.x=b$ can be solved, [LeastSquares](https://reference.wolfram.com/language/ref/LeastSquares.html) is equivalent to [LinearSolve](https://reference.wolfram.com/language/ref/LinearSolve.html):

```wolfram
m={{1,1},{1,2},{1,9}};
b={7,7,7};
```

```wolfram
LeastSquares[m,b]==LinearSolve[m,b]
(* Output *)
True
```

For a square matrix, [LinearSolve](https://reference.wolfram.com/language/ref/LinearSolve.html)[*m*,*b*] has a solution for a generic `*b*` iff [Det](https://reference.wolfram.com/language/ref/Det.html)[*m*]!=0:

```wolfram
m=RandomReal[{-1,1},{3,3}];
```

```wolfram
LinearSolve[m, {b_x,b_y,b_z}]
(* Output *)
{0.8939042337671677 (1. b_x+0.003946673292165344 b_y-0.5484044351143206 b_z),-0.14447340188117364 (1. b_x-9.650977139938437 b_y-10.094664203412352 b_z),-0.39020172326428443 (1. b_x+0.535567685443482 b_y+3.485093134849929 b_z)}
```

```wolfram
Det[m]!=0
(* Output *)
True
```

For a square matrix, [LinearSolve](https://reference.wolfram.com/language/ref/LinearSolve.html)[*m*,*b*] has a solution for a generic `*b*` iff `*m*` has full rank:

```wolfram
m=RandomReal[{-1,1},{3,3}];
```

```wolfram
LinearSolve[m, {b_x,b_y,b_z}]
(* Output *)
{-1.334186247193961 (1. b_x+1.9701590219535943 b_y+0.9738301270008943 b_z),1.1731362341036626 (1. b_x-0.05699167600311511 b_y-0.20892320500344064 b_z),-0.46481079927281316 (1. b_x-3.0853826208975788 b_y-4.462548713289445 b_z)}
```

```wolfram
MatrixRank[m]==Length[m]
(* Output *)
True
```

For a square matrix, [LinearSolve](https://reference.wolfram.com/language/ref/LinearSolve.html)[*m*,*b*] has a solution for a generic `*b*` iff `*m*` has an inverse:

```wolfram
m={{1,2,3},{4,5,6},{7,8,9}};
```

```wolfram
LinearSolve[m, {b_x,b_y,b_z}]
(* Output *)
LinearSolve
(* Output *)
LinearSolve[{{1,2,3},{4,5,6},{7,8,9}},{b_x,b_y,b_z}]
```

```wolfram
Inverse[m]
(* Output *)
Inverse
(* Output *)
Inverse[{{1,2,3},{4,5,6},{7,8,9}}]
```

For a square matrix, [LinearSolve](https://reference.wolfram.com/language/ref/LinearSolve.html)[*m*,*b*] has a solution for a generic `*b*` iff `*m*` has a trivial null space:

```wolfram
m=RandomReal[{-1,1},{3,3}];
```

```wolfram
LinearSolve[m, {b_x,b_y,b_z}]
(* Output *)
{0.4198611413357452 (1. b_x-1.3984516522506294 b_y+5.6961100541508785 b_z),1.4556905022511686 (1. b_x-0.09969897712488225 b_y+1.3467089659771374 b_z),1.3358543599095152 (1. b_x-0.5500529080005379 b_y-0.3880309480901083 b_z)}
```

```wolfram
NullSpace[m]=={}
(* Output *)
True
```

### Possible Issues

Solution found for an underdetermined system is not unique:

```wolfram
LinearSolve[{{1,2,3},{4,5,6}},{6,15}]
(* Output *)
{0,3,0}
```

All solutions are found by [Solve](https://reference.wolfram.com/language/ref/Solve.html):

```wolfram
{x,y,z}/.Solve[{{1,2,3},{4,5,6}}.{x,y,z}=={6,15},{x,y}]
(* Output *)
{{z,3-2 z,z}}
```

[LinearSolve](https://reference.wolfram.com/language/ref/LinearSolve.html) gave the solution corresponding to $z=0$:

```wolfram
%/.z->0
(* Output *)
{{0,3,0}}
```

With ill-conditioned matrices, numerical solutions may not be sufficiently accurate:

```wolfram
m=HilbertMatrix[20];
b=m.Table[1,{20}];
LinearSolve[N[m],N[b]]
(* Output *)
LinearSolve
(* Output *)
{1.0000000158380435,0.9999980941042604,1.0000589652040834,0.999168386331059,1.006727084192484,0.965595128379418,1.1094553509032565,0.8368492822354289,0.7731000071706434,2.7442639839074308,-3.060323085678936,5.7517830311152265,-0.24875546386849343,-5.663909758774916,18.527457157686417,-26.336178151690163,29.634006401717794,-18.121030384895864,8.282113709004278,-0.20037977940688634}
```

The solution is more accurate if sufficiently high precision is used:

```wolfram
LinearSolve[N[m,30],N[b,30]]
(* Output *)
{0.999999999999999999999999999999999678037405721605,1.000000000000000000000000000000120122697125283322,0.999999999999999999999999999988831506361773154847,1.000000000000000000000000000457417558760437646325,0.999999999999999999999999989616020129296707856559,1.000000000000000000000000147771319876341160341606,0.999999999999999999999998578703237010673636436332,1.000000000000000000000009708184690116879899300315,0.999999999999999999999951293349132265585415875626,1.000000000000000000000183674607134490596215090917,0.999999999999999999999471413184552282585735398904,1.000000000000000000001171136463761855887905923013,0.99999999999999999999799632635901980249092235518,1.000000000000000000002639451411603353653302440503,0.999999999999999999997348913931677623368777646628,1.000000000000000000001992278434342151519334552173,0.999999999999999999998915968517585512553088559254,1.000000000000000000000403183361701235390215875024,0.999999999999999999999908326534258958142277584221,1.000000000000000000000009609123594703870655611528}
```

Some of the linear solvers available are not deterministic. Set up a system of equations:

```wolfram
n=1000;m=SparseArray[{{i_,i_}->-2.,{i_,j_}/;Abs[i-j]==1->1.},{n,n}];
b=Table[0.,{n}];b[[99]]=1.;
```

Create a solver function:

```wolfram
solver[method_]:=LinearSolve[m,b,Method->method]
```

The `"Pardiso"` solver is not deterministic:

```wolfram
MinMax[solver["Pardiso"]-solver["Pardiso"]]
(* Output *)
{2.6645352591003757×10^-15,1.7905676941154525×10^-12}
```

The [Automatic](https://reference.wolfram.com/language/ref/Automatic.html) solver method is deterministic:

```wolfram
MinMax[solver[Automatic]-solver[Automatic]]
(* Output *)
{0.,0.}
```

[LinearSolve](https://reference.wolfram.com/language/ref/LinearSolve.html) may not give a result if a non-prime [Modulus](https://reference.wolfram.com/language/ref/Modulus.html) is provided:

```wolfram
LinearSolve[{{2,0},{0,2}},Modulus->4]
(* Output *)
LinearSolve
(* Output *)
LinearSolve[{{2,0},{0,2}},Modulus->4]
```

### Neat Examples

Solve 100,000 equations using a direct method:

```wolfram
s=SparseArray[{{i_,i_}->-2.,{i_,j_}/;Abs[i-j]==1->1.},{10^5,10^5}]
(* Output *)
SparseArray[]
```

```wolfram
Timing[LinearSolve[s,RandomReal[1,10^5]];]
(* Output *)
{0.312791,Null}
```

Solve a million equations using an iterative method:

```wolfram
s=SparseArray[{{i_,i_}->-2.,{i_,j_}/;Abs[i-j]==1->1.},{10^6,10^6}]
(* Output *)
SparseArray[]
```

```wolfram
Timing[x1=LinearSolve[s,b1=RandomReal[1,10^6],Method->"Krylov"];]
(* Output *)
{1.894551,Null}
```

Check a relative error of the solution:

```wolfram
Norm[s.x1-b1]/Norm[x1]
(* Output *)
1.5342499634813195×10^-16
```

Solve the same system of equations using a banded matrix method:

```wolfram
Timing[x2=LinearSolve[s,b2=RandomReal[1,10^6],Method->"Banded"];]
(* Output *)
{1.734355,Null}
```

Check a relative error of the solution:

```wolfram
Norm[s.x2-b2]/Norm[x2]
(* Output *)
1.3029978169639978×10^-16
```

## Tech Notes ▪Solving Linear Systems ▪Linear Algebra in the Wolfram Language ▪Implementation notes: Numerical and Related Functions

## Related Guides ▪Linear Systems ▪Equation Solving ▪Matrix Operations ▪Finite Mathematics ▪GPU Computing ▪Matrices and Linear Algebra ▪Matrix Decompositions ▪GPU Computing with NVIDIA ▪Systems Modeling ▪Finite Fields ▪GPU Programming ▪Symbolic Vectors, Matrices and Arrays ▪Structured Arrays

## Related Links [NKS|Online](http://www.wolframscience.com/nks/search/?q=LinearSolve) ([A New Kind of Science](http://www.wolframscience.com/nks/))

## History Introduced in 1988 (1.0) | Updated in 1996 (3.0) ▪ 2003 (5.0) ▪ 2014 (10.0) ▪ 2022 (13.2) ▪ 2023 (13.3) ▪ 2024 (14.0) ▪ 2024 (14.1) ▪ 2026 (15.0)
