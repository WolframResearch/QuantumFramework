# MatrixRank | [SpanFromLeft]

> [MatrixRank](https://reference.wolfram.com/language/ref/MatrixRank.html)[*m*] — gives the rank of the matrix `*m*`.

## Details and Options

[MatrixRank](https://reference.wolfram.com/language/ref/MatrixRank.html) works on both numerical and symbolic matrices.

The rank of a matrix is the number of linearly independent rows or columns.

[MatrixRank](https://reference.wolfram.com/language/ref/MatrixRank.html)[*m*,[Modulus](https://reference.wolfram.com/language/ref/Modulus.html)->*p*] finds the rank for rational matrices modulo the prime `*p*`. If `*p*` is zero, ordinary arithmetic is used.

[MatrixRank](https://reference.wolfram.com/language/ref/MatrixRank.html)[*m*,[ZeroTest](https://reference.wolfram.com/language/ref/ZeroTest.html)->*test*] evaluates `*test*[*m*[[*i*,*j*]]]` to determine whether matrix elements are zero. The default setting is [ZeroTest](https://reference.wolfram.com/language/ref/ZeroTest.html)->[Automatic](https://reference.wolfram.com/language/ref/Automatic.html).

[MatrixRank](https://reference.wolfram.com/language/ref/MatrixRank.html)[*m*,[Tolerance](https://reference.wolfram.com/language/ref/Tolerance.html)->*t*] gives the minimum rank with each element in a numerical matrix assumed to be correct only to within tolerance `*t*`.

[MatrixRank](https://reference.wolfram.com/language/ref/MatrixRank.html) works with sparse arrays and structured arrays.

## Examples

### Basic Examples

Find the number of linearly independent rows of a numerical matrix:

```wolfram
MatrixRank[{{1,2,3},{4,5,6},{7,8,9}}]
(* Output *)
2
```

Find the number of linearly independent rows of a symbolic matrix:

```wolfram
MatrixRank[{{a,b,c},{d,e,f},{g,h,i}}]
(* Output *)
3
```

Compute the rank of a rectangular matrix:

```wolfram
MatrixRank[{{0,5,2,4,4},{2,5,0,4,0},{5,1,5,4,5}}]
(* Output *)
3
```

### Scope

#### Basic Uses

Find the rank of a machine-precision matrix:

```wolfram
MatrixRank[{{1.25,3.2,3.2},{7.9,-1.4,5.1},{1.1,2.5,-1.5}}]
(* Output *)
3
```

Rank of a complex matrix:

```wolfram
MatrixRank[{{1.+I,2,3-2 I},{0,4,5I},{1+I,6,3+3I}}]
(* Output *)
2
```

Rank of an exact, rectangular matrix:

```wolfram
(m={{2,2,3,4},{3,2,1,3}})//MatrixForm
(* Output *)
({{2, 2, 3, 4}, {3, 2, 1, 3}})
```

```wolfram
MatrixRank[m]
(* Output *)
2
```

Find the rank of an arbitrary-precision matrix with more rows than columns:

```wolfram
(m=RandomReal[4,{5, 3},WorkingPrecision->20])//MatrixForm
(* Output *)
({{2.98920363445597803446138196470371894975, 2.85647670720339556446869949679623346128, 0.71497063111313092828428585046030008243}, {2.37267767251821478009278509890833674945, 1.08664396974596870438459911278883396335, 0.15612815521459883329538631535626791447}, {0.27810560063663226825696021815570446734, 0.29521172630037192409751395372197180222, 0.93873604867411056826374121664358085582}, {1.62032284128501151922630181545503802454, 2.64840354391735912435962112254639322373, 0.78935021470386928157535362526875388767}, {1.95096538059908812193102636833863527954, 2.84665784786715920452108800897672580277, 2.68391074407155711203671472064868908092}})
```

```wolfram
MatrixRank[m]
(* Output *)
3
```

Compute the rank of a symbolic matrix:

```wolfram
MatrixRank[{{a,b},{2a,2b}}]
(* Output *)
1
```

[MatrixRank](https://reference.wolfram.com/language/ref/MatrixRank.html) assumes all symbols to be independent:

```wolfram
MatrixRank[{{a,b},{c,d}}]
(* Output *)
2
```

Computing the rank of large machine-precision matrices is efficient:

```wolfram
mat=RandomReal[{0,10},{800,600}];
```

```wolfram
MatrixRank[mat] //Timing
(* Output *)
{0.322724,600}
```

Compute the rank of a matrix with finite field elements:

```wolfram
ℱ=FiniteField[29,4];
MatrixRank[{{ℱ[12],ℱ[23],ℱ[34]},{ℱ[45],ℱ[56],ℱ[67]},{ℱ[101],ℱ[98],ℱ[240]}}]
(* Output *)
2
```

#### Special Matrices

The rank of a sparse matrix:

```wolfram
SparseArray[{{1,3}->1,{2,2}->2,{3,1}->3},{3,4}]
(* Output *)
SparseArray[...]
```

```wolfram
MatrixRank[%]
(* Output *)
3
```

The rank of structured matrices:

```wolfram
QuantityArray[{{1,2,3},{4,5,6}},"Meters"]
(* Output *)
QuantityArray[...]
```

```wolfram
MatrixRank[%]
(* Output *)
2
```

```wolfram
SymmetrizedArray[{{1,1}->3,{2,2}->1,{3,1}->-5},{3,3},Symmetric[All]]
(* Output *)
SymmetrizedArray[...]
```

```wolfram
MatrixRank[%]
(* Output *)
3
```

[IdentityMatrix](https://reference.wolfram.com/language/ref/IdentityMatrix.html) always has full rank:

```wolfram
MatrixRank[IdentityMatrix[3]]
(* Output *)
3
```

```wolfram
MatrixRank[IdentityMatrix[{4,5}]]
(* Output *)
4
```

[HilbertMatrix](https://reference.wolfram.com/language/ref/HilbertMatrix.html) always has full rank:

```wolfram
MatrixRank[HilbertMatrix[3]]
(* Output *)
3
```

Compute the rank of a $15 \times 15$ matrix of univariate polynomials of degree $100$:

```wolfram
rpoly[n_]:=RandomInteger[{-2^10,2^10},{n+1}].x^Range[0,n]
SeedRandom[1234];
m=Table[rpoly[100],{15},{15}];
```

```wolfram
MatrixRank[m]//AbsoluteTiming
(* Output *)
{2.6068339,15}
```

### Options

#### Modulus

The rank of a matrix depends on the modulus used:

```wolfram
m = {{1,0,4},{2,0,3},{2,1,2}};
```

With ordinary arithmetic, `*m*` has the full rank of 3:

```wolfram
MatrixRank[m]
(* Output *)
3
```

With arithmetic modulo 5, the rank is only 2:

```wolfram
MatrixRank[m,Modulus->5]
(* Output *)
2
```

#### Tolerance

The setting of [Tolerance](https://reference.wolfram.com/language/ref/Tolerance.html) can affect the estimated rank for numerical ill-conditioned matrices:

```wolfram
m = {{1,1,1},{0,10^-10, 0},{0,0,10^-20}};
```

In exact arithmetic, `m` has full rank:

```wolfram
MatrixRank[m]
(* Output *)
3
```

With machine arithmetic, the default is to consider elements that are too small as zero:

```wolfram
MatrixRank[N[m]]
(* Output *)
2
```

With zero tolerance, even small terms may be taken into account:

```wolfram
MatrixRank[N[m], Tolerance->0]
(* Output *)
3
```

With a tolerance greater than the pivot in the middle row, the last two rows are considered zero:

```wolfram
MatrixRank[N[m], Tolerance->10^-8]
(* Output *)
1
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

Therefore the matrix rank of the matrix whose rows are the vectors is `2:

```wolfram
MatrixRank[{v1,v2,v3}]
(* Output *)
2
```

The following three vectors are linearly independent:

```wolfram
v1={1,2,3};v2={4,5,6};v3={7,8,10};
Solve[a v1 + b v2 == v3,{a,b}]
(* Output *)
{}
```

Therefore the matrix rank of the matrix whose rows are the vectors is `3:

```wolfram
MatrixRank[{v1,v2,v3}]
(* Output *)
3
```

Determine if the following vectors are linearly independent or not:

```wolfram
v_1={108,90,252,186};
v_2={240,260,520,420};
v_3={264,340,536,468};
v_4={522,705,1038,929};
```

The rank of the matrix formed from the vectors is less than four, so they are not linearly independent:

```wolfram
MatrixRank[{v_1,v_2,v_3,v_4}]
(* Output *)
2
```

Find the dimension of the column space of the following matrix:

```wolfram
a=({{1, 4, 2, -9}, {4, 12, 2, 5}, {6, 7, -11, 9}, {5, 15, 10, 12}});
```

The dimension of the space of all linear combinations of the columns equals the matrix rank:

```wolfram
MatrixRank[a]
(* Output *)
4
```

Find the dimension of the subspace spanned by the following vectors:

```wolfram
v_1={8,8,0,-4,6};v_2={2,-2,-11,1,0};v_3={-4,-8,-12,4,-4};v_4={8,12,14,-6,6};v_5={4,4,6,-2,0};
```

Since the matrix rank of the matrix formed by the vectors is three, that is the dimension of the subspace:

```wolfram
MatrixRank[{v_1,v_2,v_3,v_4,v_5}]
(* Output *)
3
```

#### Equation Solving and Invertibility

Determine if the following system of equations has a unique solution:

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

The coefficient matrix $a$ has full rank, so the system has a unique solution:

```wolfram
MatrixRank[a]
(* Output *)
3
```

Verify the result using [Solve](https://reference.wolfram.com/language/ref/Solve.html):

```wolfram
Solve[system,{x,y,z}]
(* Output *)
{{x->-(2)/(3),y->(73)/(69),z->-(53)/(69)}}
```

Determine if the following matrix has an inverse:

```wolfram
a=({{108, 90, 252, 186}, {240, 260, 520, 420}, {264, 340, 536, 468}, {522, 705, 1038, 929}});
```

The rank is less than the dimension of the matrix, so it is not invertible:

```wolfram
MatrixRank[a]
(* Output *)
2
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

Since it has full rank, its determinant must be nonzero:

```wolfram
MatrixRank[a]
(* Output *)
4
```

Confirm the result using [Det](https://reference.wolfram.com/language/ref/Det.html):

```wolfram
Det[a]
(* Output *)
-3395
```

$\lambda$ is an eigenvalue of $m$ if $m-\lambda I$ does not have full rank. Moreover, a matrix is deficient if it has an eigenvalue whose multiplicity is greater than the difference between the rank of $m-\lambda I$ and the number of columns. Show that $-4$ is an eigenvalue for the following matrix $m$:

```wolfram
m={{-8,-6,9,10},{-6,-4,0,6},{0,0,-4,0},{-10,-6,9,12}};
```

```wolfram
MatrixRank[m - (-4) IdentityMatrix[4]]
(* Output *)
2
```

Confirm the result using [Eigenvalues](https://reference.wolfram.com/language/ref/Eigenvalues.html):

```wolfram
Eigenvalues[m]
(* Output *)
{-4,-4,2,2}
```

The matrix $m$ is deficient because 2 appears twice, but the difference in rank is only one:

```wolfram
4-MatrixRank[m - 2 IdentityMatrix[4]]
(* Output *)
1
```

Confirm the result with [Eigensystem](https://reference.wolfram.com/language/ref/Eigensystem.html), which indicates deficiency by padding the eigenvector list with zeros:

```wolfram
Eigensystem[m]
(* Output *)
{{-4,-4,2,2},{{1,1,0,1},{0,3,2,0},{1,0,0,1},{0,0,0,0}}}
```

Most but not all random 10×10 0-1 matrices have full rank:

```wolfram
Table[MatrixRank[RandomInteger[1,{10,10}]],{20}]
(* Output *)
{8,10,10,10,10,10,10,10,10,10,10,10,10,10,10,9,10,8,9,9}
```

Estimate the average rank of a random 10×10 0-1 matrix:

```wolfram
Block[{n = 10000}, N[Sum[MatrixRank[RandomInteger[1,{10,10}]],{n}]/n]]
(* Output *)
9.686
```

Find the ranks of coprimality arrays:

```wolfram
m=Table[If[CoprimeQ[x,y],1,0],{x,20},{y,20}];
```

```wolfram
MatrixRank[m]
(* Output *)
13
```

```wolfram
ArrayPlot[m,PlotTheme->"Marketing"]
```

*([Graphics])*

Compute the first 50 such arrays. Only the first three have full rank:

```wolfram
Table[MatrixRank[Table[If[CoprimeQ[x,y],1,0],{x,n},{y,n}]],{n,50}]
(* Output *)
{1,2,3,3,4,5,6,6,6,7,8,8,9,10,11,11,12,12,13,13,14,15,16,16,16,17,17,17,18,19,20,20,21,22,23,23,24,25,26,26,27,28,29,29,29,30,31,31,31,31}
```

Visualize the growth in rank versus the dimension of the matrix:

```wolfram
ListLinePlot[%,Epilog->Line[{{0,0},{50,50}}],PlotRange->{0,50}]
```

*([Graphics])*

### Properties & Relations

By the rank-nullity theorem, [MatrixRank](https://reference.wolfram.com/language/ref/MatrixRank.html)[*m*] is the number of columns minus the dimension of the null space:

```wolfram
m=N[RandomInteger[1,{4,4}]];
```

```wolfram
{r,c}=Dimensions[m];
```

```wolfram
ns=NullSpace[m]
(* Output *)
{{0.37796447300922725,-0.37796447300922725,-0.3779644730092272,0.7559289460184543}}
```

```wolfram
MatrixRank[m]==c-Length[ns]
(* Output *)
True
```

The column and row rank of a matrix are equal:

```wolfram
m={{1,2,3},{4,5,6},{7,8,9}};
```

```wolfram
{MatrixRank[m],MatrixRank[mᵀ]}
(* Output *)
{2,2}
```

[MatrixRank](https://reference.wolfram.com/language/ref/MatrixRank.html)[*m*] equals the number of nonzero rows in [RowReduce](https://reference.wolfram.com/language/ref/RowReduce.html)[*m*]:

```wolfram
m={{1,4,6},{3,7,0},{9,2,3},{2,1,8},{7,8,10}};
```

```wolfram
MatrixRank[m]
(* Output *)
3
```

```wolfram
RowReduce[m]//MatrixForm
(* Output *)
({{1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {0, 0, 0}, {0, 0, 0}})
```

For a square matrix, `*m*` has full rank if and only if [Det](https://reference.wolfram.com/language/ref/Det.html)[*m*]!=0:

```wolfram
m=RandomReal[{-1,1},{3,3}];
```

```wolfram
MatrixRank[m]==Length[m]
(* Output *)
True
```

```wolfram
Det[m]!=0
(* Output *)
True
```

For a square matrix, `*m*` has full rank if and only if the null space is empty:

```wolfram
m=RandomReal[{-1,1},{3,3}];
```

```wolfram
MatrixRank[m]==Length[m]
(* Output *)
True
```

```wolfram
NullSpace[m]==={}
(* Output *)
True
```

For a square matrix, `*m*` has full rank if and only if `*m*` has an inverse:

```wolfram
m={{1,2,3},{4,5,6},{7,8,9}};
```

```wolfram
MatrixRank[m]==Length[m]
(* Output *)
False
```

```wolfram
Inverse[m]
(* Output *)
Inverse
(* Output *)
Inverse[{{1,2,3},{4,5,6},{7,8,9}}]
```

For a square matrix, `*m*` has full rank iff [LinearSolve](https://reference.wolfram.com/language/ref/LinearSolve.html)[*m*,*b*] has a solution for a generic `*b*`:

```wolfram
m=RandomReal[{-1,1},{3,3}];
```

```wolfram
MatrixRank[m]==Length[m]
(* Output *)
True
```

```wolfram
LinearSolve[m, {b_x,b_y,b_z}]
(* Output *)
{0.17494426075519107 (1. b_x+5.245603567145893 b_y+7.010786785255493 b_z),-0.38441261923041203 (1. b_x+1.446425989841002 b_y-3.2664729968477277 b_z),1.0384786892590354 (1. b_x-0.8026893013495976 b_y-0.46410297105469317 b_z)}
```

[MatrixRank](https://reference.wolfram.com/language/ref/MatrixRank.html)[*m*] is equal to [Length](https://reference.wolfram.com/language/ref/Length.html)[[SingularValueList](https://reference.wolfram.com/language/ref/SingularValueList.html)[*m*]]:

```wolfram
m = N[RandomInteger[1, {4,4}]];
```

```wolfram
SingularValueList[m]
(* Output *)
{2.052880840033716,1.2086401975729215,0.5699729199123005}
```

```wolfram
MatrixRank[m] == Length[%]
(* Output *)
True
```

The outer product of vectors has matrix rank 1:

```wolfram
v=RandomReal[1,1000];
```

```wolfram
MatrixRank[KroneckerProduct[v,v]]
(* Output *)
1
```

### Possible Issues

[MatrixRank](https://reference.wolfram.com/language/ref/MatrixRank.html) may depend on the precision of the given matrix:

```wolfram
m = {{1,2,3},{4,5,6},{7,8,9 + 10^-20}};
```

Use exact arithmetic to compute the matrix rank exactly:

```wolfram
MatrixRank[m]
(* Output *)
3
```

Use machine arithmetic. Machine numbers cannot distinguish between $9$ and $9+10^{-20}$:

```wolfram
MatrixRank[N[m]]
(* Output *)
2
```

Use 24[Hyphen]digit-precision arithmetic:

```wolfram
MatrixRank[N[m, 24]]
(* Output *)
3
```

[MatrixRank](https://reference.wolfram.com/language/ref/MatrixRank.html) assumes all symbols to be independent:

```wolfram
MatrixRank[{{1,0},{0,a}}]
(* Output *)
2
```

The special case $a=0$ gives a different result:

```wolfram
MatrixRank[{{1,0},{0,0}}]
(* Output *)
1
```

## Tech Notes ▪Solving Linear Systems ▪Implementation notes: Numerical and Related Functions ▪Implementation notes: Algebra and Calculus

## Related Guides ▪Linear Systems ▪Matrix Operations ▪Matrices and Linear Algebra ▪Finite Fields

## History Introduced in 2003 (5.0) | Updated in 2007 (6.0) ▪ 2014 (10.0) ▪ 2022 (13.2) ▪ 2024 (14.0)
