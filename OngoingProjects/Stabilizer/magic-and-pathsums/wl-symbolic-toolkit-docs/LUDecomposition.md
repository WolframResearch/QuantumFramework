# LUDecomposition

> [LUDecomposition](https://reference.wolfram.com/language/ref/LUDecomposition.html)[*m*] — generates a representation of the LU decomposition of matrix `*m*` as a list `{*l*,*u*,*p*,*c*}`, where `*l*` is lower triangular, `*u*` is upper triangular, `*p*` is a permutation matrix, and `*c*` is an approximate condition number.

## Details and Options

LU decompositions are typically used for solving linear systems of equations and for other tasks such as inverting or computing the determinant of a matrix.

An LU decomposition of `*m*` is a pair of matrices `*l*` and `*u*`, where `*l*` is a lower triangular matrix and `*u*` is an upper triangular matrix, with the property that $l.u=m$.

$\begin{pmatrix}
1 & 0 & 0 & \cdots & 0 \\
\ell_{2 1} & 1 & 0 & \cdots & 0 \\
\ell_{3 1} & \ell_{3 2} & 1 & \cdots & 0 \\
\vdots & \vdots & \vdots & \ddots & \vdots \\
\ell_{n 1} & \ell_{n 2} & \ell_{n 3} & \cdots & 1
\end{pmatrix}$ $\begin{pmatrix}
u_{11} & u_{12} & u_{13} & \cdots & u_{1 n} \\
0 & u_{22} & u_{23} & \cdots & u_{2 n} \\
0 & 0 & u_{33} & \cdots & u_{3 n} \\
\vdots & \vdots & \vdots & \ddots & \vdots \\
0 & 0 & 0 & \cdots & u_{\text{n n}}
\end{pmatrix}$ =$\begin{pmatrix}
m_{11} & m_{12} & m_{13} & \cdots & m_{1 n} \\
m_{2 1} & m_{22} & m_{23} & \cdots & m_{2 n} \\
m_{3 1} & m_{3 2} & m_{33} & \cdots & m_{3 n} \\
\vdots & \vdots & \vdots & \ddots & \vdots \\
m_{n 1} & m_{n 2} & m_{n 3} & \cdots & m_{\text{n n}}
\end{pmatrix}$

The linear system $m.x=b$ can be solved in two steps using the above decomposition by first solving the system $u.y=b$ using "forward substitution" and then solving the system $l.x=y$ using "backward substitution".

By default, [LUDecomposition](https://reference.wolfram.com/language/ref/LUDecomposition.html)[*m*] returns a list $\{l,u,p,c \}$ where $l$ is a [LowerTriangularMatrix](https://reference.wolfram.com/language/ref/LowerTriangularMatrix.html) object with $1$s on the diagonal, $u$ is an [UpperTriangularMatrix](https://reference.wolfram.com/language/ref/UpperTriangularMatrix.html) object, $p$ is a [PermutationMatrix](https://reference.wolfram.com/language/ref/PermutationMatrix.html) object that reflects the row permutations or pivoting operations used during the decomposition process, and $c$ is a condition number. The result satisfies $l.u=p.m$.

$\begin{pmatrix}
1 & 0 & 0 & \cdots & 0 \\
\ell_{2 1} & 1 & 0 & \cdots & 0 \\
\ell_{3 1} & \ell_{3 2} & 1 & \cdots & 0 \\
\vdots & \vdots & \vdots & \ddots & \vdots \\
\ell_{n 1} & \ell_{n 2} & \ell_{n 3} & \cdots & 1
\end{pmatrix}$ $\begin{pmatrix}
u_{11} & u_{12} & u_{13} & \cdots & u_{1 n} \\
0 & u_{22} & u_{23} & \cdots & u_{2 n} \\
0 & 0 & u_{33} & \cdots & u_{3 n} \\
\vdots & \vdots & \vdots & \ddots & \vdots \\
0 & 0 & 0 & \cdots & u_{\text{n n}}
\end{pmatrix}$ =$(\begin{pmatrix}
0 & 1 & \cdots & \cdots & 0 \\
0 & \cdots & \cdots & 1 & 0 \\
1 & 0 & \cdots & \cdots & 0 \\
\vdots & \vdots & \vdots & \ddots & \vdots \\
0 & 0 & 1 & \cdots & 0
\end{pmatrix})$$\begin{pmatrix}
m_{11} & m_{12} & m_{13} & \cdots & m_{1 n} \\
m_{2 1} & m_{22} & m_{23} & \cdots & m_{2 n} \\
m_{3 1} & m_{3 2} & m_{33} & \cdots & m_{3 n} \\
\vdots & \vdots & \vdots & \ddots & \vdots \\
m_{n 1} & m_{n 2} & m_{n 3} & \cdots & m_{\text{n n}}
\end{pmatrix}$

For approximate numerical matrices $m$, `*c*` is an estimate of the $L^{\infty}$ condition number of $m$, $m m$. If the matrix is not numerical or the condition number is otherwise unavailable, it is reported as `0

The following options can be given:

|  | [Modulus](https://reference.wolfram.com/language/ref/Modulus.html) | `0 | prime modulus to assume for integers |
| --- | --- | --- | --- |
|  | [Pivoting](https://reference.wolfram.com/language/ref/Pivoting.html) | [True](https://reference.wolfram.com/language/ref/True.html) | whether to allow pivoting operations |
| ▪ | [TargetStructure](https://reference.wolfram.com/language/ref/TargetStructure.html) | [Automatic](https://reference.wolfram.com/language/ref/Automatic.html) | the structure of the returned matrices |

[LUDecomposition](https://reference.wolfram.com/language/ref/LUDecomposition.html)[*m*,[Modulus](https://reference.wolfram.com/language/ref/Modulus.html)->*p*] computes the decomposition of the integer matrix `*m*` with respect to the prime modulus `*p*`.

By default, [LUDecomposition](https://reference.wolfram.com/language/ref/LUDecomposition.html) will pivot (permute) rows to improve numerical stability. With the setting [Pivoting](https://reference.wolfram.com/language/ref/Pivoting.html)->[False](https://reference.wolfram.com/language/ref/False.html), pivoting is disabled and $p$ will be an identity matrix.

Possible settings for [TargetStructure](https://reference.wolfram.com/language/ref/TargetStructure.html) include:

[Automatic](https://reference.wolfram.com/language/ref/Automatic.html), `"Structured"` | return `*l*`, `*u*` and `*p*` as appropriate structured matrices
`"Dense"` | return `*l*`, `*u*` and `*p*` as normal lists of lists
`"Sparse"` | returns `*l*`, `*u*` and `*p*` as sparse matrices
`"Classic"` | returns the classic output `{*lu*,*p*,*c*}` with `*lu*` being a packed version of `*l*`and `*u*`.

With the setting `Target"Classic"`, [LUDecomposition](https://reference.wolfram.com/language/ref/LUDecomposition.html) returns a list `{*lu*,*p*,*c*}` of three elements. The first element, `*lu*`, is a combination of lower and upper triangular matrices `*l*` and `*u*`, respectively, the second element `*p*` is a permutation list specifying how rows are pivoted, and for approximate numerical matrices `*m*`, the third element `*c*` is an estimate of the `**L**^(∞)` condition number of `*m*`. The result satisfies `*l*.*u*==*m*[[*p*]]`, where `*l*` and `*u*` must be computed from `*lu*` using [LowerTriangularize](https://reference.wolfram.com/language/ref/LowerTriangularize.html) and [UpperTriangularize](https://reference.wolfram.com/language/ref/UpperTriangularize.html).

## Examples

### Basic Examples

Compute the LU decomposition of a matrix:

```wolfram
m=({{1, 1, 1}, {2, 4, 8}, {3, 9, 27}});
{l,u,p,c}=LUDecomposition[m]
(* Output *)
{LowerTriangularMatrix[...],UpperTriangularMatrix[...],PermutationMatrix[...],0}
```

Format the result:

```wolfram
MatrixForm/@%
(* Output *)
{({{1, 0, 0}, {2, 1, 0}, {3, 3, 1}}),({{1, 1, 1}, {0, 2, 6}, {0, 0, 6}}),({{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}),0}
```

Check the result:

```wolfram
l.u==p.m
(* Output *)
True
```

Find the LU decomposition of a symbolic matrix:

```wolfram
{l,u}=Take[LUDecomposition[({{a, b, c}, {d, e, f}, {g, h, i}})],2];
{MatrixForm[l],MatrixForm[u]}
(* Output *)
{({{1, 0, 0}, {(d)/(a), 1, 0}, {(g)/(a), (-(b g)/(a)+h)/(-(b d)/(a)+e), 1}}),({{a, b, c}, {0, -(b d)/(a)+e, -(c d)/(a)+f}, {0, 0, -(c g)/(a)-((-(c d)/(a)+f) (-(b g)/(a)+h))/(-(b d)/(a)+e)+i}})}
```

Verify that $l.u$ equals the original matrix:

```wolfram
l.u//Simplify//MatrixForm
(* Output *)
({{a, b, c}, {d, e, f}, {g, h, i}})
```

### Scope

#### Basic Uses

Find the LU decomposition of a machine-precision matrix:

```wolfram
LUDecomposition[{{1.6,2.7,3.6},{1.2,3.2,5.2},{3.3,3.4,6.5}}]
(* Output *)
{LowerTriangularMatrix[...],UpperTriangularMatrix[...],PermutationMatrix[...],20.839100346020764}
```

Format the result:

```wolfram
MatrixForm/@%
(* Output *)
{({{1., 0., 0.}, {0.36363636363636365, 1., 0.}, {0.48484848484848486, 0.5354938271604939, 1.}}),({{3.3, 3.4, 6.5}, {0., 1.9636363636363638, 2.8363636363636364}, {0., 0., -1.0703703703703706}}),({{0, 0, 1}, {0, 1, 0}, {1, 0, 0}}),20.839100346020764}
```

LU decomposition for a complex matrix:

```wolfram
MatrixForm/@LUDecomposition[({{2+4 ⅈ, 9+9 ⅈ, 9+2 ⅈ}, {2+9 ⅈ, 1+3 ⅈ, 4 ⅈ}, {3+8 ⅈ, 0, 7+4 ⅈ}})]
(* Output *)
{({{1, 0, 0}, {2+(ⅈ)/(2), 1, 0}, {(19)/(10)+(ⅈ)/(5), (5598)/(5365)-(621 ⅈ)/(5365), 1}}),({{2+4 ⅈ, 9+9 ⅈ, 9+2 ⅈ}, {0, -(25)/(2)-(39 ⅈ)/(2), -17-(9 ⅈ)/(2)}, {0, 0, (9184)/(1073)+(1210 ⅈ)/(1073)}}),({{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}),0}
```

Use [LUDecomposition](https://reference.wolfram.com/language/ref/LUDecomposition.html) for an exact matrix:

```wolfram
MatrixForm/@LUDecomposition[({{3, 1, 5}, {2, 4, 6}, {8, 7, 9}})]
(* Output *)
{({{1, 0, 0}, {(3)/(2), 1, 0}, {4, (9)/(5), 1}}),({{2, 4, 6}, {0, -5, -4}, {0, 0, -(39)/(5)}}),({{0, 1, 0}, {1, 0, 0}, {0, 0, 1}}),0}
```

LU decomposition for an arbitrary-precision matrix:

```wolfram
MatrixForm/@LUDecomposition[RandomReal[5,{3,3},WorkingPrecision->10]]
(* Output *)
{({{1., 0, 0}, {0.3436671263447392762208765178767952702, 1., 0}, {0.81219318397447619234588255427532831113, -0.49368519775532140640466140188835976374, 1.}}),({{2.581661394797265529632568359375, 2.934443809208460152149200439453125, 0.165688143461011350154876708984375}, {0, 3.71394295273310781206914030216914568427, 4.8377699895848875198}, {0, 0, 3.93085191281383815560014784341266713535}}),({{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}),6.5885670551332453983}
```

The precision of the result follows the precision of the input:

```wolfram
Precision[%]
(* Output *)
10.
```

Use [LUDecomposition](https://reference.wolfram.com/language/ref/LUDecomposition.html) with a symbolic matrix:

```wolfram
MatrixForm/@LUDecomposition[{{a,b},{c,d}}]
(* Output *)
{({{1, 0}, {(c)/(a), 1}}),({{a, b}, {0, -(b c)/(a)+d}}),({{1, 0}, {0, 1}}),0}
```

The LU decomposition for a large numerical matrix is computed efficiently:

```wolfram
mat=RandomReal[{0,9},{1000,1000}];
```

```wolfram
AbsoluteTiming[LUDecomposition[mat];]
(* Output *)
{0.119568,Null}
```

LU decomposition of a wide matrix:

```wolfram
w=({{3, 2, 2}, {2, 3, -2}});
{lw,uw,pw,cw}=LUDecomposition[w];
MatrixForm/@%
(* Output *)
{({{1, 0}, {(3)/(2), 1}}),({{2, 3, -2}, {0, -(5)/(2), 5}}),({{0, 1}, {1, 0}}),0}
```

LU decomposition of a tall matrix:

```wolfram
t=({{3, 2}, {2, 3}, {2, -2}});
{lt,ut,pt,ct}=LUDecomposition[t];
MatrixForm/@%
(* Output *)
{({{1, 0, 0}, {(3)/(2), 1, 0}, {1, 2, 1}}),({{2, 3}, {0, -(5)/(2)}, {0, 0}}),({{0, 1, 0}, {1, 0, 0}, {0, 0, 1}}),0}
```

The upper triangular matrix has the same shape as the input matrix:

```wolfram
Dimensions[ut]==Dimensions[t]&& Dimensions[uw]==Dimensions[w]
(* Output *)
True
```

The lower triangular matrix is square with the same number of rows as the input:

```wolfram
SquareMatrixQ[lw]&&Length[lw]==Length[w]&&SquareMatrixQ[lt]&&Length[lt]==Length[t]
(* Output *)
True
```

The product of the triangular matrices gives the permuted matrices:

```wolfram
lw.uw==pw.w&&lw.uw==pw.w
(* Output *)
True
```

#### Special Matrices

Find the LU decomposition for sparse matrices:

```wolfram
SparseArray[{{1,3}->1,{2,2}->2,{3,1}->3},{3,3}]
(* Output *)
SparseArray[...]
```

```wolfram
MatrixForm/@LUDecomposition[%]
(* Output *)
{({{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}),({{3, 0, 0}, {0, 2, 0}, {0, 0, 1}}),({{0, 0, 1}, {0, 1, 0}, {1, 0, 0}}),0}
```

```wolfram
SparseArray[{{x_,y_}/;Abs[x-y]<3->3},{5,5}]
(* Output *)
SparseArray[...]
```

```wolfram
MatrixForm/@LUDecomposition[%]
(* Output *)
{({{1, 0, 0, 0, 0}, {0, 1, 0, 0, 0}, {0, 0, 1, 0, 0}, {1, 0, 0, 1, 0}, {1, 0, 0, 1, 1}}),({{3, 3, 3, 0, 0}, {0, 3, 3, 3, 3}, {0, 0, 3, 3, 3}, {0, 0, 0, 3, 0}, {0, 0, 0, 0, 3}}),({{1, 0, 0, 0, 0}, {0, 0, 0, 1, 0}, {0, 0, 0, 0, 1}, {0, 1, 0, 0, 0}, {0, 0, 1, 0, 0}}),0}
```

LU decompositions of structured matrices:

```wolfram
SymmetrizedArray[{{1,1}->2,{1,2}->1},{2,2},Symmetric[All]]
(* Output *)
SymmetrizedArray[...]
```

```wolfram
MatrixForm/@LUDecomposition[%]
(* Output *)
{({{1, 0}, {2, 1}}),({{1, 0}, {0, 1}}),({{0, 1}, {1, 0}}),0}
```

Use with a [QuantityArray](https://reference.wolfram.com/language/ref/QuantityArray.html) structured matrix:

```wolfram
QuantityArray[{{1,4},{3,5}},{"Meters","Seconds"}]
(* Output *)
QuantityArray[...]
```

The units go in the $u$ matrix; the $l$ matrix is dimensionless:

```wolfram
MatrixForm/@LUDecomposition[%]
(* Output *)
{({{1, 0}, {3, 1}}),({{1, 4}, {0, -7}}),({{1, 0}, {0, 1}}),0}
```

The LU decomposition of identity matrices is trivial:

```wolfram
MatrixForm/@LUDecomposition[IdentityMatrix[3]]
(* Output *)
{({{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}),({{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}),({{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}),0}
```

```wolfram
MatrixForm/@LUDecomposition[IdentityMatrix[{3,4}]]
(* Output *)
{({{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}),({{1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}}),({{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}),0}
```

LU decomposition of [HilbertMatrix](https://reference.wolfram.com/language/ref/HilbertMatrix.html):

```wolfram
MatrixForm/@LUDecomposition[HilbertMatrix[5]]
(* Output *)
{({{1, 0, 0, 0, 0}, {(5)/(4), 1, 0, 0, 0}, {(5)/(3), (10)/(3), 1, 0, 0}, {(5)/(2), 10, (15)/(2), 1, 0}, {5, 40, 60, 20, 1}}),({{(1)/(5), (1)/(6), (1)/(7), (1)/(8), (1)/(9)}, {0, -(1)/(120), -(1)/(84), -(3)/(224), -(1)/(72)}, {0, 0, (1)/(630), (1)/(336), (1)/(252)}, {0, 0, 0, -(1)/(1120), -(1)/(504)}, {0, 0, 0, 0, (1)/(630)}}),({{0, 0, 0, 0, 1}, {0, 0, 0, 1, 0}, {0, 0, 1, 0, 0}, {0, 1, 0, 0, 0}, {1, 0, 0, 0, 0}}),0}
```

LU decomposition for a matrix with finite field elements:

```wolfram
ℱ=FiniteField[19,3];
m={{ℱ[12],ℱ[23],ℱ[34]},{ℱ[45],ℱ[56],ℱ[67]},{ℱ[78],ℱ[89],ℱ[90]}};
{l,u,p,c}=LUDecomposition[m]
(* Output *)
{LowerTriangularMatrix[...],UpperTriangularMatrix[...],PermutationMatrix[...],0}
```

Format the results; the nontrivial elements of $l$ and $u$ are in the finite field:

```wolfram
MatrixForm/@%
(* Output *)
![image](img/image_001.png)
```

Verify the decomposition:

```wolfram
l.u===p.m
(* Output *)
True
```

Compute LU decomposition for a matrix of [CenteredInterval](https://reference.wolfram.com/language/ref/CenteredInterval.html) objects:

```wolfram
(m=Map[CenteredInterval,RandomReal[{-10,10},{3,3},WorkingPrecision->10],{2}])//MatrixForm
(* Output *)
({{<|Interpretation -> interpretation, Center -> 6.933338511735200882, Radius -> 6.9333385126907343121160920418×10^-10, Type -> Real|>, <|Interpretation -> interpretation, Center -> -0.1182841812260448933, Radius -> 1.182841812778084866764061189×10^-11, Type -> Real|>, <|Interpretation -> interpretation, Center -> 7.3283616919070482254, Radius -> 7.3283616935232442912706574134×10^-10, Type -> Real|>}, {<|Interpretation -> interpretation, Center -> -8.4950379142537713051, Radius -> 8.4950379196369052436921265326×10^-10, Type -> Real|>, <|Interpretation -> interpretation, Center -> -1.2864069151692092419, Radius -> 1.2864069163197322520630905274×10^-10, Type -> Real|>, <|Interpretation -> interpretation, Center -> 7.0337377930991351604, Radius -> 7.033737793593619613830014714×10^-10, Type -> Real|>}, {<|Interpretation -> interpretation, Center -> 5.3904261335264891386, Radius -> 5.3904261349685400617204322771×10^-10, Type -> Real|>, <|Interpretation -> interpretation, Center -> 4.2584752652328461409, Radius -> 4.2584752694627137081795353879×10^-10, Type -> Real|>, <|Interpretation -> interpretation, Center -> 6.5186058566905558109, Radius -> 6.5186058634258969846086984035×10^-10, Type -> Real|>}})
```

The nontrivial elements are more centered intervals:

```wolfram
{l,u,p,c}=LUDecomposition[m];
MatrixForm/@%
(* Output *)
{({{1, 0, 0}, {<|Interpretation -> interpretation, Center -> -0.6345382078262105097, Radius -> 1.2696642218493203735363294982×10^-10, Type -> Real|>, 1, 0}, {<|Interpretation -> interpretation, Center -> -0.8161633393185638852, Radius -> 1.6329200274931010561374478129×10^-10, Type -> Real|>, <|Interpretation -> interpretation, Center -> -0.3393765703088149621, Radius -> 1.6118553583868711065463230625×10^-10, Type -> Real|>, 1}}),({{<|Interpretation -> interpretation, Center -> -8.4950379142537713051, Radius -> 8.4950379196369052436921265326×10^-10, Type -> Real|>, <|Interpretation -> interpretation, Center -> -1.2864069151692092419, Radius -> 1.2864069163197322520630905274×10^-10, Type -> Real|>, <|Interpretation -> interpretation, Center -> 7.0337377930991351604, Radius -> 7.033737793593619613830014714×10^-10, Type -> Real|>}, {0, <|Interpretation -> interpretation, Center -> 3.442200926745954348, Radius -> 6.7103281958802440954059420619×10^-10, Type -> Real|>, <|Interpretation -> interpretation, Center -> 10.9817812302426318638, Radius -> 1.99213614701654329053326364374×10^-9, Type -> Real|>}, {0, 0, <|Interpretation -> interpretation, Center -> 16.7959998668138723588, Radius -> 4.9043743233223580091362236999×10^-9, Type -> Real|>}}),({{0, 1, 0}, {0, 0, 1}, {1, 0, 0}}),0}
```

Confirm equality in the sense that each element of $l.u$ is a subinterval of $p.m$:

```wolfram
MapThread[IntervalMemberQ,{l.u,p.m},2]//MatrixForm
(* Output *)
({{True, True, True}, {True, True, True}, {True, True, True}})
```

The opposite inclusion does not hold:

```wolfram
MapThread[IntervalMemberQ,{p.m,l.u},2]//MatrixForm
(* Output *)
({{False, False, False}, {False, False, False}, {False, False, False}})
```

### Options

#### Modulus

Compute the LU decomposition of an integer matrix modulus 5:

```wolfram
m=({{1, 2}, {3, 4}});
```

```wolfram
{l,u,p,c}=LUDecomposition[m,Modulus->5]
(* Output *)
{LowerTriangularMatrix[...],UpperTriangularMatrix[...],PermutationMatrix[...],0}
```

The product of the matrices is not the original matrix in ordinary arithmetic:

```wolfram
l.u==m
(* Output *)
False
```

However, the two are equal modulus 5:

```wolfram
Mod[l.u,5]==Mod[m,5]
(* Output *)
True
```

#### Pivoting

By default, decompositions are computed with pivoting to reduce numerical error:

```wolfram
m=({{3, 1, 5}, {2, 4, 6}, {8, 7, 9}});
{lt,ut,pt,ct}=LUDecomposition[m]
(* Output *)
{LowerTriangularMatrix[...],UpperTriangularMatrix[...],PermutationMatrix[...],0}
```

Compute the decomposition with pivoting disabled:

```wolfram
{lf,uf,pf,cf}=LUDecomposition[m,Pivoting->False]
(* Output *)
{LowerTriangularMatrix[...],UpperTriangularMatrix[...],PermutationMatrix[...],0}
```

Now the permutation is trivial:

```wolfram
pf==IdentityMatrix[3]
(* Output *)
True
```

However, both are valid decompositions:

```wolfram
{lt.ut==pt.m, lf.uf==m}
(* Output *)
{True,True}
```

#### TargetStructure

By default, the matrices are returned as structured objects:

```wolfram
m={{1.6,2.7,3.6},{1.2,3.2,5.2},{3.3,3.4,6.5}};{l,u,p,c}=LUDecomposition[m]
(* Output *)
{LowerTriangularMatrix[...],UpperTriangularMatrix[...],PermutationMatrix[...],20.839100346020764}
```

Have them returned as normal lists of lists:

```wolfram
{ld,ud,pd,cd}=LUDecomposition[m,TargetStructure->"Dense"]
(* Output *)
{{{1.,0.,0.},{0.36363636363636365,1.,0.},{0.48484848484848486,0.5354938271604939,1.}},{{3.3,3.4,6.5},{0.,1.9636363636363638,2.8363636363636364},{0.,0.,-1.0703703703703706}},{{0,0,1},{0,1,0},{1,0,0}},20.839100346020764}
```

Although the matrices are different objects, the results are equal:

```wolfram
{l,u,p,c}=={ld,ud,pd,cd}
(* Output *)
True
```

Obtain the results a spare arrays:

```wolfram
m=RandomReal[1,{3,3}];
LUDecomposition[m,TargetStructure->"Sparse"]
(* Output *)
{SparseArray[...],SparseArray[...],SparseArray[...],18.19356983083102}
```

This result is equal to the default, structured output:

```wolfram
%==LUDecomposition[m]
(* Output *)
True
```

Obtain the `"Classic"` form with lower and upper parts merged:

```wolfram
{lu,pList,c}=LUDecomposition[m={{1,1,1},{2,4,8},{3,9,27}},TargetStructure->"Classic"]
(* Output *)
{{{1,1,1},{2,2,6},{3,3,6}},{1,2,3},0}
```

Extract the lower matrix:

```wolfram
l=LowerTriangularize[lu,-1]+IdentityMatrix[3]
(* Output *)
{{1,0,0},{2,1,0},{3,3,1}}
```

Extract the upper matrix:

```wolfram
u=UpperTriangularize[lu]
(* Output *)
{{1,1,1},{0,2,6},{0,0,6}}
```

Convert the permutation list to a permutation matrix:

```wolfram
p=PermutationMatrix[pList]
(* Output *)
PermutationMatrix[...]
```

Check that the results agree with the results from the default return form:

```wolfram
{l,u,p,c}==LUDecomposition[m]
(* Output *)
True
```

### Applications

The LU decomposition of a matrix $m$ decomposes a matrix into lower triangular ($l$) and upper triangular ($u$) parts that satisfy $l.u=p.m$, where $p$ is a column permutation of $m$:

```wolfram
d=100;
m=RandomReal[1,{d,d}];
{l,u,p,n}=LUDecomposition[m];
```

Verify the decomposition `*l*.*u*==*p*.*m*`:

```wolfram
Max[l.u-p.m]//Chop
(* Output *)
0
```

Illustrate the structure by using [MatrixPlot](https://reference.wolfram.com/language/ref/MatrixPlot.html):

```wolfram
MatrixPlot[l,PlotTheme->"Minimal",ImageSize->100],".",MatrixPlot[u,PlotTheme->"Minimal",ImageSize->100]," == ",MatrixPlot[p.m,PlotTheme->"Minimal",ImageSize->100]
(* Output *)
![image](img/image_003.png)
```

A triangular linear system is a system of linear equations in which the first equation has one variable and each subsequent equation introduces exactly one additional variable. Rewrite the following system in four variables as two triangular linear systems in eight variables:

```wolfram
system={w+8 x+y+2 z==3,6 w+8 x+3 y+z==5,7 w+2 x+3 y+8 z==5,2 w+8 x+3 y+4 z==6};
```

Rewrite the system in matrix form $m.v=b$ and compute the LU decomposition of $m$:

```wolfram
v={w,x,y,z};
m={{1,8,1,2},{6,8,3,1},{7,2,3,8},{2,8,3,4}};
b={3,5,5,6};
{l,u,p,c}=LUDecomposition[m]
(* Output *)
{LowerTriangularMatrix[...],UpperTriangularMatrix[...],PermutationMatrix[...],0}
```

Extract the $l$ and $u$ matrices; as $l.u=m[[p]]$, the original system can be reordered into $l.u.v=b[[p]]$:

```wolfram
p.system==Thread[l.u.v==p.b]
(* Output *)
True
```

Introduce new variables $\nu=\{\omega,\xi,\psi,\zeta \}$ and set $u.v=\nu$; the result is a triangular system in $\{w,x,y,z \}$:

```wolfram
ν={ω,ξ,ψ,ζ};
tri1=Reverse[Thread[u.v==ν]];
Column[tri1]
(* Output *)
{{(281 z)/(32)==ζ}, {-8 y-11 z==ψ}, {-8 x+y==ξ}, {w+8 x+y+2 z==ω}}
```

Substituting $u.v=\nu$ into $l.u.v=b[[p]]$ gives $l.v=b[[p]]$, a triangular system in $\{\omega,\xi,\psi,\zeta \}$:

```wolfram
tri2=Thread[l.ν==p.b];
Column[tri2]
(* Output *)
{{ω==3}, {ξ+2 ω==6}, {5 ξ+ψ+6 ω==5}, {ζ+(27 ξ)/(4)+(43 ψ)/(32)+7 ω==5}}
```

The eight triangular equations together give the same result for $v$ as the original system of equations:

```wolfram
v/.Solve[Join[tri1,tri2],Join[v,ν]]
(* Output *)
{{-(35)/(281),(49)/(281),(392)/(281),(47)/(281)}}
```

```wolfram
SolveValues[system,v]
(* Output *)
{{-(35)/(281),(49)/(281),(392)/(281),(47)/(281)}}
```

LU decompositions are mainly used to solve linear systems. Here is a 5×5 random matrix:

```wolfram
m=RandomReal[1,{5,5}];
```

[LinearSolve](https://reference.wolfram.com/language/ref/LinearSolve.html)[*m*] sets up an LU decomposition in a functional form convenient for solving:

```wolfram
luf=LinearSolve[m]
(* Output *)
LinearSolveFunction[...]
```

This solves the system `*m*.*x*=*b*` for `*x*`:

```wolfram
b=ConstantArray[1,5];
x=luf[b]
(* Output *)
{0.578870471209769,1.3915209727715059,0.8830264650458737,0.4617603999365414,-1.420059281884867}
```

Verify that `*x*` is indeed the solution:

```wolfram
m.x
(* Output *)
{0.9999999999999999,1.0000000000000002,1.0000000000000002,1.,1.}
```

This can be done manually with the output of [LUDecomposition](https://reference.wolfram.com/language/ref/LUDecomposition.html) as well:

```wolfram
{l,u,p,c}=LUDecomposition[m]
(* Output *)
{LowerTriangularMatrix[...],UpperTriangularMatrix[...],PermutationMatrix[...],33.0868300467525}
```

Solve the system with two backsolves:

```wolfram
y=LinearSolve[l,p.b];
x==LinearSolve[u,y]
(* Output *)
True
```

Up to sign, the determinant of $m$ is given by the product of the diagonal elements of $u$:

```wolfram
m=RandomReal[1,{100,100},WorkingPrecision->10];{l,u,p,c}=LUDecomposition[m];
```

Compute the product of the diagonal elements of $u$

```wolfram
Tr[u,Times]
(* Output *)
-7.8110149833624858887509858×10^25
```

This equals $m$ up to sign:

```wolfram
Det[m]
(* Output *)
7.8110149833624858886327909×10^25
```

The sign can be fixed using the signature of the permutation represented by $p$:

```wolfram
Det[m]==p["PermutationSignature"]Tr[u,Times]
(* Output *)
True
```

### Properties & Relations

Compute the LU decomposition of a matrix:

```wolfram
m=({{"0.737", "0.475", "0.737", "0.723", "0.931", "0.626"}, {"0.0789", "0.199", "0.0626", "0.475", "0.79", "0.288"}, {"0.973", "0.795", "0.659", "0.576", "0.224", "0.596"}, {"0.734", "0.902", "0.883", "0.494", "0.566", "0.0264"}, {"0.854", "0.268", "0.284", "0.101", "0.7", "0.0822"}, {"0.715", "0.612", "0.593", "0.039", "0.0389", "0.857"}});
{l,u,p,c}=LUDecomposition[m];
```

`*l*` is strictly lower triangular and has ones along the diagonal:

```wolfram
{LowerTriangularMatrixQ[l], Diagonal[l]}
(* Output *)
{True,{1.,1.,1.,1.,1.,1.}}
```

`*u*` is upper triangular:

```wolfram
UpperTriangularMatrixQ[u]
(* Output *)
True
```

`*l*.*u*` is equal to the permutation of the rows of `*m*` given by `*p*`:

```wolfram
l.u-p.m//Chop
(* Output *)
{{0,0,0,0,0,0},{0,0,0,0,0,0},{0,0,0,0,0,0},{0,0,0,0,0,0},{0,0,0,0,0,0},{0,0,0,0,0,0}}
```

If $m$ is a rectangular matrix, the $l$ matrix is square with the same number of rows as $m$:

```wolfram
m=({{"0.308", "0.208", "0.68", "0.736"}, {"0.029", "0.132", "0.427", "0.943"}, {"0.476", "0.248", "0.132", "0.561"}});
{l,u,p,c}=LUDecomposition[m];
l//MatrixForm
(* Output *)
({{1., 0., 0.}, {0.06084231044686835, 1., 0.}, {0.6475813372549272, 0.40676156770101546, 1.}})
```

The $u$ matrix has the same shape as $m$:

```wolfram
u//MatrixForm
(* Output *)
({{0.47606447358059123, 0.2480736415405338, 0.13202945134919064, 0.5612849334301935}, {0., 0.11698635971343478, 0.4188294632400183, 0.9088135511348062}, {0., 0., 0.42440923686131027, 0.002824144671109141}})
```

The basic identity $l.u=p.m$ still holds:

```wolfram
l.u==p.m
(* Output *)
True
```

If `*m*` decomposes to `{*l*,*u*,*p*,*c*}`, then [Det](https://reference.wolfram.com/language/ref/Det.html)[*m*] is the product of the diagonal entries of `*u*` and the signature of the permutation represented by `*p*`:

```wolfram
m=RandomInteger[9,{5,5}];
{l,u,p,c}=LUDecomposition[m];
Det[m]==p["PermutationSignature"]*Times@@Diagonal[u]
(* Output *)
True
```

If a matrix is singular, the $u$ matrix will have a zero along the diagonal:

```wolfram
m=({{1, 2, 3}, {4, 5, 6}, {7, 8, 9}});
LUDecomposition[m][[2]]//MatrixForm
(* Output *)
LUDecomposition
(* Output *)
({{1, 2, 3}, {0, -3, -6}, {0, 0, 0}})
```

The [CholeskyDecomposition](https://reference.wolfram.com/language/ref/CholeskyDecomposition.html) of a positive definite Hermitian matrix `*h*`:

```wolfram
h={{2,1},{1,2}};
PositiveDefiniteMatrixQ[h]&&HermitianMatrixQ[h]
(* Output *)
True
```

```wolfram
MatrixForm[u=CholeskyDecomposition[h]]
(* Output *)
({{Sqrt[2], (1)/(Sqrt[2])}, {0, Sqrt[(3)/(2)]}})
```

This gives a kind of LU decomposition via [ConjugateTranspose](https://reference.wolfram.com/language/ref/ConjugateTranspose.html):

```wolfram
MatrixForm[l=ConjugateTranspose[u]]
(* Output *)
({{Sqrt[2], 0}, {(1)/(Sqrt[2]), Sqrt[(3)/(2)]}})
```

```wolfram
l.u==h
(* Output *)
True
```

This is generally a different decomposition from the one given by [LUDecomposition](https://reference.wolfram.com/language/ref/LUDecomposition.html):

```wolfram
First[LUDecomposition[h]]//MatrixForm
(* Output *)
({{1, 0}, {2, 1}})
```

[LDLDecomposition](https://reference.wolfram.com/language/ref/LDLDecomposition.html)[*m*] can be computed directly from [LUDecomposition](https://reference.wolfram.com/language/ref/LUDecomposition.html)[*m*] when the latter returns a trivial permutation:

```wolfram
mat={{3,-4,5},{-4,8,-2},{5,-2,11}};
{l_1,d}=LDLDecomposition[mat];
{l_2,u,perm,c}=LUDecomposition[mat];
perm==IdentityMatrix[3]
(* Output *)
True
```

The two lower diagonal matrices are the same, and the diagonal entries of $u$ and $d$ are the same:

```wolfram
l_1==l_2&&Diagonal[d]==Diagonal[u]
(* Output *)
True
```

The condition number of an invertible numerical matrix is $m m$:

```wolfram
m=({{"0.964", "0.304", "0.481"}, {"0.359", "0.132", "0.397"}, {"0.771", "0.143", "0.237"}});
```

```wolfram
{Last[LUDecomposition[m]],Norm[m,Infinity]Norm[Inverse[m],Infinity]}
(* Output *)
{54.561271910806525,54.561271910806525}
```

The value returned by [LUDecomposition](https://reference.wolfram.com/language/ref/LUDecomposition.html) is only an estimate:

```wolfram
a=({{"0.252", "0.528", "0.862"}, {"0.738", "0.0933", "0.645"}, {"0.815", "0.999", "0.433"}});
```

```wolfram
{Last[LUDecomposition[a]],Norm[a,Infinity]Norm[Inverse[a],Infinity]}
(* Output *)
{5.14189680093829,6.379036973708673}
```

If the condition number is not available, it is reported as `0:

```wolfram
LUDecomposition[{{1.5,a},{3.2,4.1}}] //Last
(* Output *)
0
```

The condition number $c$ is the same condition number reported by [LinearSolve](https://reference.wolfram.com/language/ref/LinearSolve.html)[*m*]:

```wolfram
m=RandomReal[1,{100,100}];
```

```wolfram
LinearSolve[m]
(* Output *)
LinearSolveFunction[...]
```

```wolfram
Last[LUDecomposition[m]]
(* Output *)
4026.0196813735593
```

### Interactive Examples

Visualize LU decompositions, where $P=p=p$ and hence $m=P.l.u$:

```wolfram
Manipulate[
SeedRandom[sr];
With[{m=Table[Random[Integer,{0,max}],{n},{n}]},
{l,u,perm,c}=Quiet[LUDecomposition[m]];
p = Transpose[perm];
Text[Column[Row[Graphics[Text[Style[If[m, ==, Dot[p, l, u], StyleBox[m,FontSlant->"Italic"] = P l u, StyleBox[m,FontSlant->"Italic"] ≠ P l u; the decomposition failed.], Label, 14]], ImageSize -> 40020]]Graphics[Text[Row[If[m, ==, Dot[p, l, u], MatrixForm[m, TableSpacing -> 3] = MatrixForm[p, TableSpacing -> 3] MatrixForm[l, TableSpacing -> 3] MatrixForm[u, TableSpacing -> 3], MatrixForm[m, TableSpacing -> 3] ≠ MatrixForm[Dot[p, l, u], TableSpacing -> 3]]]], ImageSize -> 400180]]]],
{{n,4,"matrix rank"},2,5,1,Appearance-> "Labeled"},{{max,4,"maximum entry"},1,9,1,Appearance-> "Labeled"},{{sr,1,"seed random"},1,100,1}, SaveDefinitions-> True, TrackedSymbols:> {n,max,sr}]
```

## Tech Notes ▪Advanced Matrix Operations ▪Implementation notes: Numerical and Related Functions

## Related Guides ▪Matrix Decompositions ▪Matrices and Linear Algebra ▪Linear Systems ▪Finite Fields

## History Introduced in 1996 (3.0) | Updated in 2024 (14.0) ▪ 2025 (14.3) ▪ 2026 (15.0)
