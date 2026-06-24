# RowReduce | [SpanFromLeft]

> [RowReduce](https://reference.wolfram.com/language/ref/RowReduce.html)[*m*] — gives the row[Hyphen]reduced form of the matrix `*m*`.

## Details and Options

[RowReduce](https://reference.wolfram.com/language/ref/RowReduce.html) performs a version of Gaussian elimination, adding multiples of rows together so as to produce zero elements when possible. The final matrix is in reduced row echelon form.

If `*m*` is a non[Hyphen]degenerate square matrix, [RowReduce](https://reference.wolfram.com/language/ref/RowReduce.html)[*m*] is [IdentityMatrix](https://reference.wolfram.com/language/ref/IdentityMatrix.html)[[Length](https://reference.wolfram.com/language/ref/Length.html)[*m*]].

If `*m*` is a sufficiently non[Hyphen]degenerate rectangular matrix with $k$ rows and more than $k$ columns, then the first $k$ columns of [RowReduce](https://reference.wolfram.com/language/ref/RowReduce.html)[*m*] will form an identity matrix.

[RowReduce](https://reference.wolfram.com/language/ref/RowReduce.html) works on both numerical and symbolic matrices.

The following options can be given:

| [Method](https://reference.wolfram.com/language/ref/Method.html) | [Automatic](https://reference.wolfram.com/language/ref/Automatic.html) | method to use |
| --- | --- | --- |
| [Modulus](https://reference.wolfram.com/language/ref/Modulus.html) | 0 | integer modulus to use |
| [Tolerance](https://reference.wolfram.com/language/ref/Tolerance.html) | [Automatic](https://reference.wolfram.com/language/ref/Automatic.html) | numerical tolerance to use |
| [ZeroTest](https://reference.wolfram.com/language/ref/ZeroTest.html) | [Automatic](https://reference.wolfram.com/language/ref/Automatic.html) | function to test whether matrix elements should be considered to be zero |

[RowReduce](https://reference.wolfram.com/language/ref/RowReduce.html)[*m*,[Modulus](https://reference.wolfram.com/language/ref/Modulus.html)->*n*] finds the reduction of rational matrices modulo the integer `*n*`. If `*n*` is zero, ordinary arithmetic is used. If `*n*` is not prime, the computation may fail.

[RowReduce](https://reference.wolfram.com/language/ref/RowReduce.html)[*m*,[ZeroTest](https://reference.wolfram.com/language/ref/ZeroTest.html)->*test*] evaluates `*test*[*m*[[*i*,*j*]]]` to determine whether matrix elements are zero.

Possible settings for the [Method](https://reference.wolfram.com/language/ref/Method.html) option include `"CofactorExpansion"`, `"DivisionFreeRowReduction"`, and `"OneStepRowReduction"`. The default setting of [Automatic](https://reference.wolfram.com/language/ref/Automatic.html) switches among these methods depending on the matrix given.

## Examples

### Basic Examples

Do row reduction on a square matrix:

```wolfram
RowReduce[{{1,2,3},{4,5,6},{7,8,9}}]
(* Output *)
{{1,0,-1},{0,1,2},{0,0,0}}
```

```wolfram
MatrixForm[%]
(* Output *)
({{1, 0, -1}, {0, 1, 2}, {0, 0, 0}})
```

Do row reduction on a rectangular matrix:

```wolfram
RowReduce[{{1,2,3,1,0,0},{4,5,6,0,1,0},{7,8,9,0,0,1}}]
(* Output *)
{{1,0,-1,0,-(8)/(3),(5)/(3)},{0,1,2,0,(7)/(3),-(4)/(3)},{0,0,0,1,-2,1}}
```

```wolfram
MatrixForm[%]
(* Output *)
({{1, 0, -1, 0, -(8)/(3), (5)/(3)}, {0, 1, 2, 0, (7)/(3), -(4)/(3)}, {0, 0, 0, 1, -2, 1}})
```

Row reduce a matrix with symbolic entries:

```wolfram
RowReduce[{{3,1,a},{2,1,b}}]
(* Output *)
{{1,0,a-b},{0,1,-2 a+3 b}}
```

### Scope

#### Basic Uses

Find the row reduction of a real machine-number matrix:

```wolfram
RowReduce[{{1.5,4.75,-3.2},{-1.7,6.7,-9.3},{4.9,-8.65,15.4}}]//MatrixForm
(* Output *)
({{1, 0., 1.2543448275862072}, {0, 1, -1.0697931034482757}, {0, 0, 0}})
```

Row reduce a complex machine-precision matrix:

```wolfram
MatrixForm[RowReduce[RandomComplex[3+ⅈ,{5,5}]]]
(* Output *)
({{1, 0.+0. ⅈ, 0.+0. ⅈ, 0.+0. ⅈ, 0.+0. ⅈ}, {0, 1, 0.+0. ⅈ, 0.+0. ⅈ, 0.+0. ⅈ}, {0, 0, 1, 0.+0. ⅈ, 0.+0. ⅈ}, {0, 0, 0, 1, 0.+0. ⅈ}, {0, 0, 0, 0, 1}})
```

Row reduce an arbitrary-precision matrix:

```wolfram
m=N[{{π,Sqrt[2],ℯ},{(1)/(Sqrt[2]),1+Sqrt[2],2 π},{-(1)/(Sqrt[2])+2 π,-1+Sqrt[2],2 ℯ-2 π}},20]
(* Output *)
{{3.1415926535897932384626433832795028845,1.41421356237309504880168872420969807857,2.71828182845904523536028740430832960766},{0.70710678118654752440084436210484903928,2.41421356237309504880168872420969807857,6.28318530717958647692528676655900576901},{5.57607852599303895252444240445415673236,0.41421356237309504880168872420969807857,-0.84662165026149600620471182900334086863}}
```

```wolfram
RowReduce[m]//MatrixForm
(* Output *)
({{1, 0, -0.35283797279317628768556164351162160563}, {0, 1, 2.70592441870815279090198991045952890046}, {0, 0, 0}})
```

Row reduce an exact matrix:

```wolfram
RowReduce[{{0,0,0,π},{2,2,2,2},{3+ISqrt[2],0,0,0},{0,4,4,0}}]//MatrixForm
(* Output *)
({{1, 0, 0, 0}, {0, 1, 1, 0}, {0, 0, 0, 1}, {0, 0, 0, 0}})
```

Row reduction of a symbolic matrix:

```wolfram
RowReduce[{{a,b,c},{d,e,f},{a+d,b+e,c+f}}]
(* Output *)
{{1,0,(c e-b f)/(-b d+a e)},{0,1,(c d-a f)/(b d-a e)},{0,0,0}}
```

[RowReduce](https://reference.wolfram.com/language/ref/RowReduce.html) assumes all symbols to be independent:

```wolfram
RowReduce[{{a,b,c},{d,e,f},{g,h,i}}]
(* Output *)
{{1,0,0},{0,1,0},{0,0,1}}
```

Row reduction of rectangular matrices:

```wolfram
RowReduce[({{1.2, 3.4}, {5.6, 7.8}, {9.10, 11.12}})]//MatrixForm
(* Output *)
({{1, 0.}, {0, 1}, {0, 0}})
```

```wolfram
RowReduce[RandomReal[1,{4, 6},WorkingPrecision->20]]//MatrixForm
(* Output *)
({{1, 0, 0, 0, -1.24556656091145600112175921998584676548, -1.98726785432734826516524732276509443563}, {0, 1, 0, 0, 0.94770932280781210390107097032519134457, 0.53815936428448247780670681313356161673}, {0, 0, 1, 0, -0.04044866368993816175991716994035621277, 1.27046201363706873607829210055756224079}, {0, 0, 0, 1, 1.15964962016624241911148428331238726339, 1.70819587784857479776530011094248576872}})
```

```wolfram
RowReduce[RandomComplex[1,{3, 5}]]//MatrixForm
(* Output *)
({{1, 0.+0. ⅈ, 0.+0. ⅈ, -0.16183252356595518+0. ⅈ, 0.6106246260558967+0. ⅈ}, {0, 1, 0.+0. ⅈ, 0.05152202660516861+0. ⅈ, -0.39185652060507936+0. ⅈ}, {0, 0, 1, 0.666171640353511+0. ⅈ, 0.6592750358048632+0. ⅈ}})
```

The row reduction of a large numerical matrix is computed efficiently:

```wolfram
RandomReal[{0,9},{500,300}];
RowReduce[%]; //AbsoluteTiming
(* Output *)
{0.013067,Null}
```

#### Special Matrices

Row reduction of a sparse matrix:

```wolfram
SparseArray[{{1,2}->1,{3,1}->5,{2,3}->2},{3,3}]
(* Output *)
SparseArray[...]
```

```wolfram
RowReduce[%]
(* Output *)
{{1,0,0},{0,1,0},{0,0,1}}
```

Row reduction of a structured matrix:

```wolfram
SymmetrizedArray[{{1,2}->3,{2,2}->5,{2,3}->7,{3,2}->13},{3,3},Symmetric[All]]
(* Output *)
SymmetrizedArray[...]
```

```wolfram
RowReduce[%]
(* Output *)
{{1,0,(10)/(3)},{0,1,0},{0,0,0}}
```

```wolfram
QuantityArray[{{1,2},{2,4}},{"Meters","Seconds"}]
(* Output *)
QuantityArray[...]
```

```wolfram
RowReduce[%]
(* Output *)
{{1,2},{0,0}}
```

[RowReduce](https://reference.wolfram.com/language/ref/RowReduce.html) does not alter an identity matrix:

```wolfram
RowReduce[IdentityMatrix[3]]
(* Output *)
{{1,0,0},{0,1,0},{0,0,1}}
```

```wolfram
RowReduce[IdentityMatrix[{3,4}]]
(* Output *)
{{1,0,0,0},{0,1,0,0},{0,0,1,0}}
```

Row reduction of a Hilbert matrix:

```wolfram
RowReduce[HilbertMatrix[7]]//MatrixForm
(* Output *)
({{1, 0, 0, 0, 0, 0, 0}, {0, 1, 0, 0, 0, 0, 0}, {0, 0, 1, 0, 0, 0, 0}, {0, 0, 0, 1, 0, 0, 0}, {0, 0, 0, 0, 1, 0, 0}, {0, 0, 0, 0, 0, 1, 0}, {0, 0, 0, 0, 0, 0, 1}})
```

Row reduce a matrix with finite field elements:

```wolfram
ℱ=FiniteField[29,4];
RowReduce[{{ℱ[12],ℱ[23],ℱ[34]},{ℱ[45],ℱ[56],ℱ[67]},{ℱ[101],ℱ[98],ℱ[240]}}]//MatrixForm
(* Output *)
![image](img/image_001.png)
```

Compute the row reduction of a $10 \times 11$ matrix of univariate polynomials of degree $100$:

```wolfram
rpoly[n_]:=RandomInteger[{-2^10,2^10},{n+1}].x^Range[0,n]
SeedRandom[1234];
m=Table[rpoly[100],{10},{11}];
```

```wolfram
RowReduce[m]//Short[#,4]&//AbsoluteTiming
(* Output *)
{0.9994312,{{1,0,0,0,0,0,0,0,0,0,(<<1498>>+6851379239732580085987753872528 x^1000)/(<<1506>>+5879966341995275387562752769861 x^1000)},{0,1,0,0,0,0,0,0,0,0,(<<1>>)/(<<1>>)},{0,0,1,0,0,0,0,0,0,0,(<<1>>)/(<<1>>)},<<4>>,{0,0,0,0,0,0,0,1,0,0,(<<1>>)/(<<1>>)},{0,0,0,0,0,0,0,0,1,0,(<<1>>)/(<<1>>)},{0,0,0,0,0,0,0,0,0,1,(<<1>>)/(<<1>>)}}}
```

### Options

#### Modulus

`*m*` is a 3×3 matrix of integers between 0 and 4:

```wolfram
m={{1,0,4},{2,0,3},{2,1,2}};
```

Row reduction of `*m*` in ordinary arithmetic:

```wolfram
RowReduce[m]
(* Output *)
{{1,0,0},{0,1,0},{0,0,1}}
```

Row reduction of `*m*` in arithmetic modulo 5:

```wolfram
RowReduce[m,Modulus->5]
(* Output *)
{{1,0,4},{0,1,4},{0,0,0}}
```

#### Tolerance

`*m*` is an ill-conditioned matrix:

```wolfram
m={{1,0,1},{1,10^-20,0},{0,0,1}};
```

In exact arithmetic, `*m*` is clearly non-degenerate:

```wolfram
RowReduce[m]
(* Output *)
{{1,0,0},{0,1,0},{0,0,1}}
```

With machine arithmetic, the default is to consider elements that are too small as zero:

```wolfram
RowReduce[N[m]]
(* Output *)
{{1,0.,0.},{0,0,1},{0,0,0}}
```

With zero tolerance, even small terms may be taken into account:

```wolfram
RowReduce[N[m],Tolerance->0]
(* Output *)
RowReduce
(* Output *)
{{1,0.,0.},{0,1,0.},{0,0,1}}
```

With an augmented matrix, you can see how possible solution components are amplified:

```wolfram
MatrixForm[RowReduce[Transpose[Append[Transpose[m],{1,1,1}]]]]
(* Output *)
({{1, 0, 0, 0}, {0, 1, 0, 100000000000000000000}, {0, 0, 1, 1}})
```

#### ZeroTest

By default, symbolic expressions are considered nonzero:

```wolfram
m=({{1, -2, 3, 2}, {1, 1, 1, k}, {2, -1, 4, k^2}});
RowReduce[m]//MatrixForm
(* Output *)
({{1, 0, (5)/(3), 0}, {0, 1, -(2)/(3), 0}, {0, 0, 0, 1}})
```

In this case, [RowReduce](https://reference.wolfram.com/language/ref/RowReduce.html) is missing the special case $k=0$ that produces a singular matrix:

```wolfram
RowReduce[m /. k->2]//MatrixForm
(* Output *)
({{1, 0, (5)/(3), 2}, {0, 1, -(2)/(3), 0}, {0, 0, 0, 0}})
```

Write a function that tests if an expression might be zero:

```wolfram
test[k_][e_]:=Reduce[e==0,k,Reals]=!=False
{test[k][k^2+1],test[k][k^2]}
(* Output *)
{False,True}
```

Pass `test[k]` as the value of [ZeroTest](https://reference.wolfram.com/language/ref/ZeroTest.html) to get a more symbolic reduction:

```wolfram
MatrixForm[gen=RowReduce[{{1,-2,3,2},{1,1,1,k},{2,-1,4,k^2}},ZeroTest->test[k]]]
(* Output *)
({{1, 0, (5)/(3), (2 (1+k))/(3)}, {0, 1, -(2)/(3), (1)/(3) (-2+k)}, {0, 0, 0, -3 (-2+k)+3 (-4+k^2)}})
```

This catches the special case at $k=0$, but gives a non-reduced matrix for other values of $k$:

```wolfram
{MatrixForm[gen/.k->2], MatrixForm[gen /. k->1]}
(* Output *)
{({{1, 0, (5)/(3), 2}, {0, 1, -(2)/(3), 0}, {0, 0, 0, 0}}),({{1, 0, (5)/(3), (4)/(3)}, {0, 1, -(2)/(3), -(1)/(3)}, {0, 0, 0, -6}})}
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

The row-reduced form has a zero row:

```wolfram
RowReduce[{v1,v2,v3}]
(* Output *)
{{1,0,-1},{0,1,2},{0,0,0}}
```

The following three vectors are linearly independent:

```wolfram
v1={1,2,3};v2={4,5,6};v3={7,8,10};
Solve[a v1 + b v2 == v3,{a,b}]
(* Output *)
{}
```

Therefore the row-reduced form is the identity matrix:

```wolfram
RowReduce[{v1,v2,v3}]
(* Output *)
{{1,0,0},{0,1,0},{0,0,1}}
```

Determine if the following vectors are linearly independent or not:

```wolfram
v_1={108,90,252,186};
v_2={240,260,520,420};
v_3={264,340,536,468};
v_4={522,705,1038,929};
```

As $\{v_{1},v_{2},v_{3},v_{4}\}$ does not reduce to an identity matrix, they are not linearly independent:

```wolfram
RowReduce[{v_1,v_2,v_3,v_4}]
(* Output *)
{{1,0,(26)/(9),(44)/(27)},{0,1,-(2)/(3),(1)/(9)},{0,0,0,0},{0,0,0,0}}
```

Find the dimension of the column space of the following matrix:

```wolfram
a=({{1, 4, 2, -9}, {4, 12, 2, 5}, {6, 7, -11, 9}, {5, 15, 10, 12}});
```

Since the row-reduced form is an identity matrix, the dimension of the column space equals the number of columns:

```wolfram
RowReduce[a]
(* Output *)
{{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}}
```

Find the dimension of the subspace spanned by the following vectors:

```wolfram
v_1={8,8,0,-4,6};v_2={2,-2,-11,1,0};v_3={-4,-8,-12,4,-4};v_4={8,12,14,-6,6};v_5={4,4,6,-2,0};
```

Since the row-reduced form has three nonzero rows, that is the dimension of the subspace:

```wolfram
RowReduce[{v_1,v_2,v_3,v_4,v_5}]
(* Output *)
{{1,0,0,0,-1},{0,1,0,-(1)/(2),(7)/(4)},{0,0,1,0,-(1)/(2)},{0,0,0,0,0},{0,0,0,0,0}}
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

The coefficient matrix $a$ reduces to the identity matrix, so the system has a unique solution:

```wolfram
RowReduce[a]
(* Output *)
{{1,0,0},{0,1,0},{0,0,1}}
```

Verify the result using [Solve](https://reference.wolfram.com/language/ref/Solve.html):

```wolfram
Solve[system,{x,y,z}]
(* Output *)
{{x->-(2)/(3),y->(73)/(69),z->-(53)/(69)}}
```

Solve the system $m.x=b$, with $m$ a matrix and $b$ a vector, using row reduction:

```wolfram
m={{1,2,3},{5,6,7},{7,8,9}};
b={1,1,1};
```

Form an augmented matrix $m|b$:

```wolfram
MatrixForm[a=Join[m,Transpose[{b}],2]]
(* Output *)
({{1, 2, 3, 1}, {5, 6, 7, 1}, {7, 8, 9, 1}})
```

Do row reduction on the augmented matrix:

```wolfram
MatrixForm[r=RowReduce[a]]
(* Output *)
({{1, 0, -1, -1}, {0, 1, 2, 1}, {0, 0, 0, 0}})
```

The last column is the solution of $m.x=b$:

```wolfram
x=r[[All,-1]]
(* Output *)
{-1,1,0}
```

```wolfram
m.x==b
(* Output *)
True
```

Do it for another right[Hyphen]hand side:

```wolfram
b1={1,-2,1};
a1=Transpose[Join[Transpose[m],{b1}]];
MatrixForm[r1=RowReduce[a1]]
(* Output *)
({{1, 0, -1, 0}, {0, 1, 2, 0}, {0, 0, 0, 1}})
```

There is no solution to $m.x=b 1$ since there is a leading 1 in the last column; confirm using [Solve](https://reference.wolfram.com/language/ref/Solve.html):

```wolfram
Solve[m.{x1,x2,x3}==b1,{x1,x2,x3}]
(* Output *)
{}
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

Reduce the augmented matrix $m|b$ to find one solution:

```wolfram
RowReduce[Join[m,Transpose[{b}],2]]//MatrixForm
(* Output *)
({{1, 0, -1, (9)/(2)}, {0, 1, 2, -(11)/(3)}, {0, 0, 0, 0}})
```

Extract the last column, which forms one solution:

```wolfram
x1=%[[All,-1]]
(* Output *)
{(9)/(2),-(11)/(3),0}
```

Since there is a zero row, the null space is nonempty. Row reduce the augmented matrix $m|I$:

```wolfram
RowReduce[Join[m,IdentityMatrix[3],2]]//MatrixForm
(* Output *)
({{1, 0, -1, 0, 7, -(9)/(2)}, {0, 1, 2, 0, -(17)/(3), (11)/(3)}, {0, 0, 0, 1, -2, 1}})
```

The last half of the last row is an element of the null space:

```wolfram
x0=%[[-1,-3;;-1]]
(* Output *)
{1,-2,1}
```

The general solution is the sum of $x1$ and any multiple of $x0$:

```wolfram
x1+C[1]x0
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

It does not reduce to an identity matrix, so it is not invertible:

```wolfram
RowReduce[a]
(* Output *)
{{1,0,(26)/(9),(44)/(27)},{0,1,-(2)/(3),(1)/(9)},{0,0,0,0},{0,0,0,0}}
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

Since it reduces to an identity matrix, its determinant must be nonzero:

```wolfram
RowReduce[a]
(* Output *)
{{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}}
```

Confirm the result using [Det](https://reference.wolfram.com/language/ref/Det.html):

```wolfram
Det[a]
(* Output *)
-3395
```

$\lambda$ is an eigenvalue of $m$ if $m-\lambda I$ does not reduce to an identity matrix. A matrix is deficient if it has an eigenvalue whose multiplicity is greater than the number of zero rows in the row-reduced form. Show that $-4$ is an eigenvalue for the following matrix $m$:

```wolfram
m={{-8,-6,9,10},{-6,-4,0,6},{0,0,-4,0},{-10,-6,9,12}};
```

```wolfram
RowReduce[m - (-4) IdentityMatrix[4]]
(* Output *)
{{1,0,0,-1},{0,1,-(3)/(2),-1},{0,0,0,0},{0,0,0,0}}
```

Confirm the result using [Eigenvalues](https://reference.wolfram.com/language/ref/Eigenvalues.html):

```wolfram
Eigenvalues[m]
(* Output *)
{-4,-4,2,2}
```

The matrix $m$ is deficient because 2 appears twice, but the row-reduced form only has one zero row:

```wolfram
RowReduce[m - 2 IdentityMatrix[4]]
(* Output *)
{{1,0,0,-1},{0,1,0,0},{0,0,1,0},{0,0,0,0}}
```

Confirm the result with [Eigensystem](https://reference.wolfram.com/language/ref/Eigensystem.html), which indicates deficiency by padding the eigenvector list with zeros:

```wolfram
Eigensystem[m]
(* Output *)
{{-4,-4,2,2},{{1,1,0,1},{0,3,2,0},{1,0,0,1},{0,0,0,0}}}
```

Find the inverse of the following matrix using row reduction:

```wolfram
m={{a,b,c},{d,e,f},{g,h,i}};
```

Form the augmented matrix $m|I$:

```wolfram
MatrixForm[aug=Join[m,IdentityMatrix[3],2]]
(* Output *)
({{a, b, c, 1, 0, 0}, {d, e, f, 0, 1, 0}, {g, h, i, 0, 0, 1}})
```

Row-reduce the augmented matrix:

```wolfram
r=RowReduce[aug];
```

The first three columns of $r$ are the identity matrix:

```wolfram
r[[All,1;;3]]
(* Output *)
{{1,0,0},{0,1,0},{0,0,1}}
```

The last three columns are the inverse of $m$:

```wolfram
MatrixForm[minv= r[[All, -3;;-1]]]
(* Output *)
({{(f h-e i)/(c e g-b f g-c d h+a f h+b d i-a e i), (c h-b i)/(-c e g+b f g+c d h-a f h-b d i+a e i), (c e-b f)/(c e g-b f g-c d h+a f h+b d i-a e i)}, {(f g-d i)/(-c e g+b f g+c d h-a f h-b d i+a e i), (c g-a i)/(c e g-b f g-c d h+a f h+b d i-a e i), (c d-a f)/(-c e g+b f g+c d h-a f h-b d i+a e i)}, {(e g-d h)/(c e g-b f g-c d h+a f h+b d i-a e i), (b g-a h)/(-c e g+b f g+c d h-a f h-b d i+a e i), (b d-a e)/(c e g-b f g-c d h+a f h+b d i-a e i)}})
```

Verify the result using [Inverse](https://reference.wolfram.com/language/ref/Inverse.html):

```wolfram
Simplify[minv==Inverse[m]]
(* Output *)
True
```

Estimate the fraction of random 10×10 0-1 matrices that are invertible:

```wolfram
Block[{n = 10000}, N[Sum[If[RowReduce[RandomInteger[1,{10,10}]]===IdentityMatrix[10],1,0],{n}]/n]]
(* Output *)
0.7046
```

### Properties & Relations

[RowReduce](https://reference.wolfram.com/language/ref/RowReduce.html) returns the leading zeroes as exact integers, irrespective the precision of the input:

```wolfram
m={{0,0,0,1},{2,2,2,2},{3,0,0,0},{0,4,4,0}};
MatrixForm /@ {RowReduce[m],RowReduce[N[m]],RowReduce[N[m,10]]}
(* Output *)
{({{1, 0, 0, 0}, {0, 1, 1, 0}, {0, 0, 0, 1}, {0, 0, 0, 0}}),({{1, 0., 0., 0.}, {0, 1, 1., 0.}, {0, 0, 0, 1}, {0, 0, 0, 0}}),({{1, 0, 0, 0}, {0, 1, 1., 0}, {0, 0, 0, 1}, {0, 0, 0, 0}})}
```

A square matrix `*m*` is invertible iff [RowReduce](https://reference.wolfram.com/language/ref/RowReduce.html)[*m*] is [IdentityMatrix](https://reference.wolfram.com/language/ref/IdentityMatrix.html)[[Length](https://reference.wolfram.com/language/ref/Length.html)[*m*]]:

```wolfram
m={{1,0,4},{2,0,3},{2,1,2}};
```

```wolfram
RowReduce[m]==IdentityMatrix[Length[m]]
(* Output *)
True
```

```wolfram
Inverse[m]
(* Output *)
{{-(3)/(5),(4)/(5),0},{(2)/(5),-(6)/(5),1},{(2)/(5),-(1)/(5),0}}
```

Indeed, the inverse can be found by inverting the augmented matrix form by `*m*` and the identity matrix:

```wolfram
RowReduce[Join[m,IdentityMatrix[Length[m]],2]][[All,-Length[m];;-1]]
(* Output *)
{{-(3)/(5),(4)/(5),0},{(2)/(5),-(6)/(5),1},{(2)/(5),-(1)/(5),0}}
```

For a square matrix, `*m*` reduces to an identity matrix if and only if [Det](https://reference.wolfram.com/language/ref/Det.html)[*m*]!=0:

```wolfram
m=RandomReal[{-1,1},{3,3}];
```

```wolfram
RowReduce[m]==IdentityMatrix[Length[m]]
(* Output *)
True
```

```wolfram
Det[m]!=0
(* Output *)
True
```

For a square matrix, `*m*` reduces to an identity matrix if and only if the null space is empty:

```wolfram
m=RandomReal[{-1,1},{3,3}];
```

```wolfram
RowReduce[m]==IdentityMatrix[Length[m]]
(* Output *)
True
```

```wolfram
NullSpace[m]==={}
(* Output *)
True
```

For a square matrix, `*m*` reduces to an identity matrix iff [LinearSolve](https://reference.wolfram.com/language/ref/LinearSolve.html)[*m*,*b*] has a solution for a generic `*b*`:

```wolfram
m={{1,2,3},{4,5,6},{7,8,9}};
```

```wolfram
RowReduce[m]==IdentityMatrix[Length[m]]
(* Output *)
False
```

```wolfram
LinearSolve[m, {b_x,b_y,b_z}]
(* Output *)
LinearSolve
(* Output *)
LinearSolve[{{1,2,3},{4,5,6},{7,8,9}},{b_x,b_y,b_z}]
```

An $m \times n$ matrix $a$ with $n \geq m$ has rank $m$ iff the first $m$ columns of [RowReduce](https://reference.wolfram.com/language/ref/RowReduce.html)[*m*] form an identity matrix:

```wolfram
a={{3,7,1,2,7},{4,5,6,6,6},{0,5,5,1,4}};
m=Length[a];
r=RowReduce[a];
```

```wolfram
MatrixRank[a]==m
(* Output *)
True
```

```wolfram
r[[All,1;;Length[a]]]==IdentityMatrix[Length[a]]
(* Output *)
True
```

Format the row-reduced matrix:

```wolfram
MatrixForm[r]
(* Output *)
({{1, 0, 0, (17)/(15), (67)/(135)}, {0, 1, 0, -(4)/(15), (106)/(135)}, {0, 0, 1, (7)/(15), (2)/(135)}})
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

Even though the vectors are not the same, they are a basis for the same space:

```wolfram
RegionEqual[ParametricRegion[a v1+b v2,{a,b}],ParametricRegion[a n1+b n2,{a,b}]]
(* Output *)
True
```

### Possible Issues

The computation may fail if a non-prime modulus is provided:

```wolfram
RowReduce[{{1,2},{3,2}},Modulus->6]
(* Output *)
RowReduce
(* Output *)
RangeSpace[{{1,2},{3,2}},Modulus->6]
```

## Tech Notes ▪Solving Linear Systems ▪Implementation notes: Numerical and Related Functions

## Related Guides ▪Linear Systems ▪Matrix Operations ▪Finite Mathematics ▪Matrices and Linear Algebra ▪Finite Fields

## History Introduced in 1988 (1.0) | Updated in 1996 (3.0) ▪ 2022 (13.2) ▪ 2024 (14.0)
