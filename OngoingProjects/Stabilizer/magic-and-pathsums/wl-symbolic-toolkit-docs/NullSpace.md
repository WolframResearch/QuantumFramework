# NullSpace | [SpanFromLeft]

> [NullSpace](https://reference.wolfram.com/language/ref/NullSpace.html)[*m*] — gives a list of vectors that forms a basis for the null space of the matrix `*m*`.

## Details and Options

[NullSpace](https://reference.wolfram.com/language/ref/NullSpace.html) works on both numerical and symbolic matrices.

The following options can be given:

| [Method](https://reference.wolfram.com/language/ref/Method.html) | [Automatic](https://reference.wolfram.com/language/ref/Automatic.html) | method to use |
| --- | --- | --- |
| [Modulus](https://reference.wolfram.com/language/ref/Modulus.html) | 0 | integer modulus to use |
| [Tolerance](https://reference.wolfram.com/language/ref/Tolerance.html) | [Automatic](https://reference.wolfram.com/language/ref/Automatic.html) | numerical tolerance to use |
| [ZeroTest](https://reference.wolfram.com/language/ref/ZeroTest.html) | [Automatic](https://reference.wolfram.com/language/ref/Automatic.html) | function to test whether matrix elements should be considered to be zero |

[NullSpace](https://reference.wolfram.com/language/ref/NullSpace.html)[*m*,[Modulus](https://reference.wolfram.com/language/ref/Modulus.html)->*n*] finds the null space of rational matrices modulo the integer `*n*`. If `*n*` is zero, ordinary arithmetic is used. If `*n*` is not prime, the computation may fail.

[NullSpace](https://reference.wolfram.com/language/ref/NullSpace.html)[*m*,[ZeroTest](https://reference.wolfram.com/language/ref/ZeroTest.html)->*test*] evaluates `*test*[*m*[[*i*,*j*]]]` to determine whether matrix elements are zero.

Possible settings for the [Method](https://reference.wolfram.com/language/ref/Method.html) option include `"CofactorExpansion"`, `"DivisionFreeRowReduction"`, and `"OneStepRowReduction"`. The default setting of [Automatic](https://reference.wolfram.com/language/ref/Automatic.html) switches among these methods depending on the matrix given.

## Examples

### Basic Examples

Find the null space of a 3×3 matrix:

```wolfram
m = {{1,2,3},{4,5,6},{7,8,9}};
```

```wolfram
NullSpace[m]
(* Output *)
{{1,-2,1}}
```

The action of `m` on the vector is the zero vector:

```wolfram
m.{1,-2,1}
(* Output *)
{0,0,0}
```

The null space of a symbolic matrix:

```wolfram
m=({{a, b}, {2a, 2b}});
```

```wolfram
NullSpace[m]
(* Output *)
{{-(b)/(a),1}}
```

Verify that multiplication by `m` gives the zero vector:

```wolfram
m.First[%]
(* Output *)
{0,0}
```

Compute the null space of a rectangular matrix:

```wolfram
m=({{0, 5, 2, 4, 4}, {2, 5, 0, 4, 0}, {5, 1, 5, 4, 5}});
```

```wolfram
ns=NullSpace[m]
(* Output *)
{{25,-10,-71,0,48},{-1,-2,-1,3,0}}
```

Verify that both vectors returned are in the null space:

```wolfram
m.First[ns]
(* Output *)
{0,0,0}
```

```wolfram
m.Last[ns]
(* Output *)
{0,0,0}
```

### Scope

#### Basic Uses

Null space of a machine-precision matrix:

```wolfram
NullSpace[{{1.5,4.75,-3.2},{-1.7,6.7,-9.3},{4.9,-8.65,15.4}}]
(* Output *)
{{-0.650536544993824,0.5548231188665165,0.5186265616016303}}
```

Null space of a complex matrix:

```wolfram
NullSpace[{{1.+I,2,3-2 I},{0,4,5I},{1+I,6,3+3I}}]
(* Output *)
{{-0.6817850187554061-0.6213519768293823 ⅈ,0.17869058441972133+0.24285544249373833 ⅈ,-0.19428435399499075+0.14295246753577706 ⅈ}}
```

Null space of an exact matrix:

```wolfram
NullSpace[{{1,2,3,10},{4,5,6,11},{7,8,9,12},{13,14,15,16}}]
(* Output *)
{{1,-2,1,0}}
```

Null space of an arbitrary-precision matrix:

```wolfram
m=N[{{π,Sqrt[2],ℯ},{(1)/(Sqrt[2]),1+Sqrt[2],2 π},{-(1)/(Sqrt[2])+2 π,-1+Sqrt[2],2 ℯ-2 π}},20]
(* Output *)
{{3.1415926535897932384626433832795028845,1.41421356237309504880168872420969807857,2.71828182845904523536028740430832960766},{0.70710678118654752440084436210484903928,2.41421356237309504880168872420969807857,6.28318530717958647692528676655900576901},{5.57607852599303895252444240445415673236,0.41421356237309504880168872420969807857,-0.84662165026149600620471182900334086863}}
```

```wolfram
NullSpace[m]
(* Output *)
{{0.12140494312399516235164166862737674301,-0.9310579514741113840586765129011873151,0.34408128513752495421377835322892111534}}
```

Find the null space symbolically:

```wolfram
NullSpace[{{a,b,c},{c,b,a},{0,0,0}}]
(* Output *)
{{1,-(a+c)/(b),1}}
```

Null space of non-square matrices:

```wolfram
{{3,2,2,4},{2,3,-2,7},{3,2,5,7}}//MatrixForm
(* Output *)
({{3, 2, 2, 4}, {2, 3, -2, 7}, {3, 2, 5, 7}})
```

```wolfram
NullSpace[%]
(* Output *)
{{12,-23,-5,5}}
```

```wolfram
{{a,2a,π a,(a)/(3),ⅈ a}}//MatrixForm
(* Output *)
({a, 2 a, a π, (a)/(3), ⅈ a})
```

```wolfram
NullSpace[%]//MatrixForm
(* Output *)
({{-ⅈ, 0, 0, 0, 1}, {-(1)/(3), 0, 0, 1, 0}, {-π, 0, 1, 0, 0}, {-2, 1, 0, 0, 0}})
```

The null space of a large numerical matrix is computed efficiently:

```wolfram
mat=RandomReal[{2,10},{200,200}];
```

```wolfram
NullSpace[mat]//Timing
(* Output *)
{0.087139,{}}
```

Null space of a matrix with finite field elements:

```wolfram
ℱ=FiniteField[29,4];
NullSpace[{{ℱ[12],ℱ[23],ℱ[34]},{ℱ[45],ℱ[56],ℱ[67]},{ℱ[101],ℱ[98],ℱ[240]}}]
(* Output *)
![image](img/image_001.png)
```

#### Special Matrices

Null space of a sparse matrix:

```wolfram
SparseArray[{{2,1}->12,{1,1}->4,{2,3}->13,{3,1}->7,{4,4}->2},{4,4}]
(* Output *)
SparseArray[...]
```

```wolfram
NullSpace[%]
(* Output *)
{{0,1,0,0}}
```

Null space of structured matrices:

```wolfram
SymmetrizedArray[{{1,2}->3,{2,2}->5,{2,3}->7,{3,2}->13},{3,3},Symmetric[All]]
(* Output *)
SymmetrizedArray[...]
```

```wolfram
NullSpace[%]
(* Output *)
{{-10,0,3}}
```

```wolfram
QuantityArray[{{1,2},{2,4}},{"Meters","Seconds"}]
(* Output *)
QuantityArray[...]
```

```wolfram
NullSpace[%]
(* Output *)
{{-2,1}}
```

[IdentityMatrix](https://reference.wolfram.com/language/ref/IdentityMatrix.html)[*n*] always has an empty null space:

```wolfram
NullSpace[IdentityMatrix[3]]
(* Output *)
{}
```

The null space of [IdentityMatrix](https://reference.wolfram.com/language/ref/IdentityMatrix.html)[{*m*,*n*}] is nonempty:

```wolfram
NullSpace[IdentityMatrix[{3,4}]]
(* Output *)
{{0,0,0,1}}
```

Compute the null space for [HilbertMatrix](https://reference.wolfram.com/language/ref/HilbertMatrix.html):

```wolfram
NullSpace[HilbertMatrix[{3,4}]]
(* Output *)
{{-(1)/(20),(3)/(5),-(3)/(2),1}}
```

Compute the null space of a $10 \times 12$ matrix of univariate polynomials of degree $100$:

```wolfram
rpoly[n_]:=RandomInteger[{-2^10,2^10},{n+1}].x^Range[0,n]
SeedRandom[1234];
m=Table[rpoly[100],{10},{12}];
```

```wolfram
NullSpace[m]//Short[#,3]&//AbsoluteTiming
(* Output *)
{0.5974312,{{3313658063082096615328956582013-31717721400268291637200422634499 x+<<1490>>+35355576464418195311795993707288 x^999-1447804039860325240039469411136 x^1000,<<1>>,<<1>>,<<1>>,<<4>>,<<1>>,<<1>>,<<1>>,0},{<<1>>}}}
```

### Options

#### Modulus

`m` is a 3×3 random matrix of integers between 0 and 4:

```wolfram
m = {{1,0,4},{2,0,3},{2,1,2}};
```

Use arithmetic modulo 5 to compute the null space:

```wolfram
NullSpace[m, Modulus->5]
(* Output *)
{{1,1,1}}
```

The vector is in the null space modulo 5:

```wolfram
m.Transpose[%]
(* Output *)
{{5},{5},{5}}
```

```wolfram
Mod[%,5]
(* Output *)
{{0},{0},{0}}
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

Therefore the null space of the matrix whose rows are the vectors is nonempty:

```wolfram
NullSpace[{v1,v2,v3}]
(* Output *)
{{1,-2,1}}
```

The following three vectors are linearly independent:

```wolfram
v1={1,2,3};v2={4,5,6};v3={7,8,10};
Solve[a v1 + b v2 == v3,{a,b}]
(* Output *)
{}
```

Therefore the null space of the matrix whose rows are the vectors is empty:

```wolfram
NullSpace[{v1,v2,v3}]
(* Output *)
{}
```

Determine if the following vectors are linearly independent or not:

```wolfram
v_1={108,90,252,186};
v_2={240,260,520,420};
v_3={264,340,536,468};
v_4={522,705,1038,929};
```

The matrix formed from the vectors has a nonempty null space, so they are not linearly independent:

```wolfram
NullSpace[{v_1,v_2,v_3,v_4}]
(* Output *)
{{-44,-3,0,27},{-26,6,9,0}}
```

Find the dimension of the column space of the following matrix:

```wolfram
a=({{1, 4, 2, -9}, {4, 12, 2, 5}, {6, 7, -11, 9}, {5, 15, 10, 12}});
```

Since the null space is empty, the dimension of the column space equals the number of columns:

```wolfram
NullSpace[a]
(* Output *)
{}
```

Find the dimension of the subspace spanned by the following vectors:

```wolfram
v_1={8,8,0,-4,6};v_2={2,-2,-11,1,0};v_3={-4,-8,-12,4,-4};v_4={8,12,14,-6,6};v_5={4,4,6,-2,0};
```

Since the matrix rank of the matrix formed by the vectors is three, that is the dimension of the subspace:

```wolfram
With[{a={v_1,v_2,v_3,v_4,v_5}},Length[a]-Length[NullSpace[a]]]
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

The coefficient matrix $a$ has an empty null space, so the system has a unique solution:

```wolfram
NullSpace[a]
(* Output *)
{}
```

Verify the result using [Solve](https://reference.wolfram.com/language/ref/Solve.html):

```wolfram
Solve[system,{x,y,z}]
(* Output *)
{{x->-(2)/(3),y->(73)/(69),z->-(53)/(69)}}
```

$m$ is a 3×3 singular matrix with a nonempty null space:

```wolfram
m = {{1,2,3},{4,5,6},{7,8,9}};
```

```wolfram
ns = NullSpace[m]
(* Output *)
{{1,-2,1}}
```

Find a solution $x_{1}$ for $m.x=b$:

```wolfram
b = {1,1,1};
x_1= LinearSolve[m,b]
(* Output *)
{-1,1,0}
```

All solutions are given by $x_{0}+x_{1}$, where $x_{0}$ is any vector in the null space:

```wolfram
x = x_1 + Sum[c[i] ns[[i]], {i, Length[ns]}]
(* Output *)
{-1+c[1],1-2 c[1],c[1]}
```

```wolfram
Simplify[m.x==b]
(* Output *)
True
```

Determine if the following matrix has an inverse:

```wolfram
a=({{108, 90, 252, 186}, {240, 260, 520, 420}, {264, 340, 536, 468}, {522, 705, 1038, 929}});
```

The null space is not trivial, so the matrix is not invertible:

```wolfram
NullSpace[a]
(* Output *)
{{-44,-3,0,27},{-26,6,9,0}}
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

Since the null space is empty, its determinant must be nonzero:

```wolfram
NullSpace[a]
(* Output *)
{}
```

Confirm the result using [Det](https://reference.wolfram.com/language/ref/Det.html):

```wolfram
Det[a]
(* Output *)
-3395
```

$\lambda$ is an eigenvalue of $m$ if the null space of $m-\lambda I$ is nontrivial. A matrix is deficient if it has an eigenvalue whose multiplicity is greater than the dimension of the null space of $m-\lambda I$. Show that $-4$ is an eigenvalue for the following matrix $m$:

```wolfram
m={{-8,-6,9,10},{-6,-4,0,6},{0,0,-4,0},{-10,-6,9,12}};
```

```wolfram
NullSpace[m - (-4) IdentityMatrix[4]]
(* Output *)
{{1,1,0,1},{0,3,2,0}}
```

Confirm the result using [Eigenvalues](https://reference.wolfram.com/language/ref/Eigenvalues.html):

```wolfram
Eigenvalues[m]
(* Output *)
{-4,-4,2,2}
```

The matrix $m$ is deficient because 2 appears twice, but eigenspace is one-dimensional:

```wolfram
NullSpace[m - 2 IdentityMatrix[4]]
(* Output *)
{{1,0,0,1}}
```

Confirm the result with [Eigensystem](https://reference.wolfram.com/language/ref/Eigensystem.html), which indicates deficiency by padding the eigenvector list with zeros:

```wolfram
Eigensystem[m]
(* Output *)
{{-4,-4,2,2},{{1,1,0,1},{0,3,2,0},{1,0,0,1},{0,0,0,0}}}
```

Find a basis for the eigenvectors of the following matrix:

```wolfram
m={{6,4,-6,-4},{4,6,-15,-4},{0,0,-4,0},{4,4,-6,-2}};
```

First, compute the eigenvalues:

```wolfram
λ=Eigenvalues[m]
(* Output *)
{6,-4,2,2}
```

For each unique eigenvalue, find the null space:

```wolfram
Table[NullSpace[m-k IdentityMatrix[4]],{k,DeleteDuplicates[λ]}]
(* Output *)
{{{1,1,0,1}},{{0,3,2,0}},{{1,0,0,1},{-1,1,0,0}}}
```

Combine the two one-dimensional spaces and one two-dimensional space using [Join](https://reference.wolfram.com/language/ref/Join.html) and [Apply](https://reference.wolfram.com/language/ref/Apply.html):

```wolfram
Join@@%
(* Output *)
{{1,1,0,1},{0,3,2,0},{1,0,0,1},{-1,1,0,0}}
```

Confirm the result using [Eigenvectors](https://reference.wolfram.com/language/ref/Eigenvectors.html):

```wolfram
%==Eigenvectors[m]
(* Output *)
True
```

Estimate the fraction of random 10×10 0-1 matrix that are invertible:

```wolfram
Block[{n = 10000}, N[Sum[If[NullSpace[RandomInteger[1,{10,10}]]==={},1,0],{n}]/n]]
(* Output *)
0.7136
```

### Properties & Relations

For any linear combination $x$ of elements of the null space of $m$, $m.x$ gives zero:

```wolfram
m = {{1,1,1,1,1},{1,0,0,0,0},{0,0,0,0,1},{0,1,1,1,0}, {1,0,0,0,1}};
```

```wolfram
ns = NullSpace[m];
x=Sum[c[i] ns[[i]], {i, Length[ns]}]
(* Output *)
{0,-c[1]-c[2],c[2],c[1],0}
```

```wolfram
m.x
(* Output *)
{0,0,0,0,0}
```

By the rank-nullity theorem, the dimension of the null space is the number of columns minus [MatrixRank](https://reference.wolfram.com/language/ref/MatrixRank.html)[*m*]:

```wolfram
m=N[RandomInteger[1,{4,4}]];
```

```wolfram
{r,c}=Dimensions[m];
```

```wolfram
ns=NullSpace[m]
(* Output *)
{{-0.5773502691896257,-0.5773502691896257,0.,0.5773502691896257}}
```

```wolfram
Length[ns]==c-MatrixRank[m]
(* Output *)
True
```

For a square matrix, `*m*` has a trivial null space if and only if [Det](https://reference.wolfram.com/language/ref/Det.html)[*m*]!=0:

```wolfram
m=RandomReal[{-1,1},{3,3}];
```

```wolfram
NullSpace[m]=={}
(* Output *)
True
```

```wolfram
Det[m]!=0
(* Output *)
True
```

For a square matrix, `*m*` has a trivial null space if and only if `*m*` has full rank:

```wolfram
m=RandomReal[{-1,1},{3,3}];
```

```wolfram
NullSpace[m]==={}
(* Output *)
True
```

```wolfram
MatrixRank[m]==Length[m]
(* Output *)
True
```

For a square matrix, `*m*` has a trivial null space if and only if `*m*` has an inverse:

```wolfram
m={{1,2,3},{4,5,6},{7,8,9}};
```

```wolfram
NullSpace[m]==={}
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

For a square matrix, `*m*` has a trivial null space iff [LinearSolve](https://reference.wolfram.com/language/ref/LinearSolve.html)[*m*,*b*] has a solution for a generic `*b*`:

```wolfram
m=RandomReal[{-1,1},{3,3}];
```

```wolfram
NullSpace[m]=={}
(* Output *)
True
```

```wolfram
LinearSolve[m, {b_x,b_y,b_z}]
(* Output *)
{-1.4897328857789365 (1. b_x-0.20269386815208917 b_y-0.2906547444200082 b_z),1.045121702643124 (1. b_x-2.1682170679361716 b_y-1.3806775809438867 b_z),2.3160999489410004 (1. b_x-0.6529454619297507 b_y-0.010861294209304373 b_z)}
```

The dimension of the null space of a matrix $m$ is known as its nullity $\nu(m)$:

```wolfram
ν[m_?MatrixQ]:=Length[NullSpace[m]]
```

The nullity of a product of two square matrices satisfies Sylvester's law of nullity $max(\nu(m_{1}),\nu(m_{2}))\leq \nu(m_{1}.m_{2})\leq \nu(m_{1})+\nu(m_{2})$:

```wolfram
m1={{1,-1,0},{1,-1,0},{-1,1,0}};
m2={{1,1,0},{1,1,0},{0,0,0}};
```

```wolfram
{ν[m1],ν[m2],ν[m1.m2]}
(* Output *)
{2,2,3}
```

```wolfram
Max[ν[m1],ν[m2]]<=ν[m1.m2]<=ν[m1]+ν[m2]
(* Output *)
True
```

In this case, the product in the opposite order has a different nullity:

```wolfram
ν[m2.m1]
(* Output *)
2
```

But it still satisfies the law:

```wolfram
Max[ν[m1],ν[m2]]<=ν[m2.m1]<=ν[m1]+ν[m2]
(* Output *)
True
```

The null space of a square matrix `*m*` can be computed using [RowReduce](https://reference.wolfram.com/language/ref/RowReduce.html):

```wolfram
m={{1,2,3,4},{5,6,7,8},{9,10,11,12},{13,14,15,16}};
```

Do row reduction on the matrix augmented with the identity matrix:

```wolfram
MatrixForm[r=RowReduce[Join[m,IdentityMatrix[Length[m]],2]]]
(* Output *)
({{1, 0, -1, -2, 0, 0, -(7)/(2), (5)/(2)}, {0, 1, 2, 3, 0, 0, (13)/(4), -(9)/(4)}, {0, 0, 0, 0, 1, 0, -3, 2}, {0, 0, 0, 0, 0, 1, -2, 1}})
```

The augmented half of a row is in the null space if the row has a leading 1 in the augmented half:

```wolfram
{v1,v2}=r[[3;;4,5;;8]]
(* Output *)
{{1,0,-3,2},{0,1,-2,1}}
```

Get null vectors using [NullSpace](https://reference.wolfram.com/language/ref/NullSpace.html):

```wolfram
{n1,n2}=NullSpace[m]
(* Output *)
{{2,-3,0,1},{1,-2,1,0}}
```

Even though the vectors are not the same, they are a basis for the same vector subspace:

```wolfram
RowReduce[{v1,v2,n1,n2}]//MatrixForm
(* Output *)
({{1, 0, -3, 2}, {0, 1, -2, 1}, {0, 0, 0, 0}, {0, 0, 0, 0}})
```

### Possible Issues

The computation may fail if a non-prime modulus is provided:

```wolfram
NullSpace[{{1,2},{3,2}},Modulus->6]
(* Output *)
RowReduce
(* Output *)
RangeSpace[{{1,2},{3,2}},Modulus->6]
```

## Tech Notes ▪Solving Linear Systems ▪Implementation notes: Numerical and Related Functions ▪Implementation notes: Algebra and Calculus

## Related Guides ▪Linear Systems ▪Matrix Operations ▪Matrices and Linear Algebra ▪Finite Fields

## Related Links [NKS|Online](http://www.wolframscience.com/nks/search/?q=NullSpace) ([A New Kind of Science](http://www.wolframscience.com/nks/))

## History Introduced in 1988 (1.0) | Updated in 1996 (3.0) ▪ 2022 (13.2) ▪ 2024 (14.0)
