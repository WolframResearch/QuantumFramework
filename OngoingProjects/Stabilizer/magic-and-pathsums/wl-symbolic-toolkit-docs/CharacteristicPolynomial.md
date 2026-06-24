# CharacteristicPolynomial | [SpanFromLeft]

> [CharacteristicPolynomial](https://reference.wolfram.com/language/ref/CharacteristicPolynomial.html)[*m*,*x*] — gives the characteristic polynomial for the matrix `*m*`.
> [CharacteristicPolynomial](https://reference.wolfram.com/language/ref/CharacteristicPolynomial.html)[{*m*,*a*},*x*] — gives the generalized characteristic polynomial with respect to `*a*`.

## Details

`*m*` must be a square matrix.

It can contain numeric or symbolic entries.

[CharacteristicPolynomial](https://reference.wolfram.com/language/ref/CharacteristicPolynomial.html)[*m*,*x*] is essentially equivalent to [Det](https://reference.wolfram.com/language/ref/Det.html)[*m*-*id* *x*] where `*id*` is the identity matrix of appropriate size.

[CharacteristicPolynomial](https://reference.wolfram.com/language/ref/CharacteristicPolynomial.html)[{*m*,*a*},*x*] is essentially [Det](https://reference.wolfram.com/language/ref/Det.html)[*m*-*a* *x*].

## Examples

### Basic Examples

Find the characteristic polynomial of a matrix with integer entries:

```wolfram
CharacteristicPolynomial[{{1,2},{3,4}},x]
(* Output *)
-2-5 x+x^2
```

Visualize the polynomial:

```wolfram
Plot[%, {x,-5,10}]
```

*([Graphics])*

Find the characteristic polynomial in $x$ of the symbolic matrix $m$:

```wolfram
m=({{a, b}, {c, d}});
```

```wolfram
CharacteristicPolynomial[m,x]
(* Output *)
-b c+a d-a x-d x+x^2
```

Compare with a direct computation:

```wolfram
Det[m-x IdentityMatrix[2]]
(* Output *)
-b c+a d-a x-d x+x^2
```

Compute the characteristic polynomials of the identity matrix and zero matrix:

```wolfram
CharacteristicPolynomial[IdentityMatrix[3],λ]
(* Output *)
1-3 λ+3 λ^2-λ^3
```

```wolfram
CharacteristicPolynomial[ConstantArray[0,{10,10}],λ]
(* Output *)
λ^10
```

### Scope

#### Basic Uses

Find the characteristic polynomial of a machine-precision matrix:

```wolfram
CharacteristicPolynomial[{{1.1,2.2,3.25},{0.76,4.6,5},{0.1,0.1,6.1}},x]
(* Output *)
19.968799999999977-37.332999999999984 x+11.799999999999999 x^2-x^3
```

Arbitrary-precision matrix:

```wolfram
CharacteristicPolynomial[Table[N[1/(i+j+1),20],{i,3},{j,3}],x]
(* Output *)
2.64550264550264550264550264550265×10^-6-0.01257936507936507936507936507936511278 x+0.67619047619047619047619047619047624131 x^2-x^3
```

Characteristic polynomial of a complex matrix:

```wolfram
CharacteristicPolynomial[{{1.2+I,3-2 I,3 π},{-0.2,5I,2},{1,2.3,ℯ}},x]
(* Output *)
(-15.815837907173742-40.501511564476296 ⅈ)+(15.162839766618532-21.909690970754298 ⅈ) x+(3.918281828459049+6.000000000000002 ⅈ) x^2-x^3
```

Characteristic polynomial of an exact matrix:

```wolfram
CharacteristicPolynomial[{{(1)/(3),(1)/(2),(3)/(5)},{(1)/(2),(4)/(5),1},{(3)/(5),1,(9)/(7)}},x]
(* Output *)
(1)/(10500)-(239 x)/(2100)+(254 x^2)/(105)-x^3
```

Visualize the result:

```wolfram
Plot[%, {x,-1.5,2.5}]
```

*([Graphics])*

The characteristic polynomials of large numerical matrices are computed efficiently:

```wolfram
m=RandomReal[{1,9},{100,100}];
```

```wolfram
CharacteristicPolynomial[m,x];//Timing
(* Output *)
{0.020851,Null}
```

Characteristic polynomial of a matrix with finite field elements:

```wolfram
ℱ=FiniteField[43,2];
CharacteristicPolynomial[{{ℱ[12],ℱ[23],ℱ[34]},{ℱ[45],ℱ[56],ℱ[67]},{ℱ[78],ℱ[89],ℱ[90]}},x]
(* Output *)
![image](img/image_001.png)
```

Characteristic polynomial of a matrix containing [CenteredInterval](https://reference.wolfram.com/language/ref/CenteredInterval.html) objects:

```wolfram
(m=Map[CenteredInterval,RandomReal[{-10,10},{3,3},WorkingPrecision->10],{2}])//MatrixForm
(* Output *)
({{<|Interpretation -> interpretation, Center -> -4.5751592656597495079, Radius -> 4.5751592661971707265422537603×10^-10, Type -> Real|>, <|Interpretation -> interpretation, Center -> 8.437393291387706995, Radius -> 8.4373932927178652008137760276×10^-10, Type -> Real|>, <|Interpretation -> interpretation, Center -> 2.226566969184204936, Radius -> 2.2265669711619973103466918474×10^-10, Type -> Real|>}, {<|Interpretation -> interpretation, Center -> -0.9152538678608834743, Radius -> 9.15253868787740954005016647×10^-11, Type -> Real|>, <|Interpretation -> interpretation, Center -> -2.2509080066811293364, Radius -> 2.2509080081604959655550146636×10^-10, Type -> Real|>, <|Interpretation -> interpretation, Center -> -0.6426929100416600704, Radius -> 6.426929106094492194500844562×10^-11, Type -> Real|>}, {<|Interpretation -> interpretation, Center -> 5.5024141946341842413, Radius -> 5.5024142013582189036924319225×10^-10, Type -> Real|>, <|Interpretation -> interpretation, Center -> 8.6770807742141187191, Radius -> 8.6770807813246753070757222304×10^-10, Type -> Real|>, <|Interpretation -> interpretation, Center -> 3.4514569758903235197, Radius -> 3.4514569761251445711991436838×10^-10, Type -> Real|>}})
```

```wolfram
p=CharacteristicPolynomial[m,x]
(* Output *)
<|Interpretation -> interpretation, Center -> 16.7396155849201022647, Radius -> 1.6893276622020891863940050825477×10^-7, Type -> Real|>+<|Interpretation -> interpretation, Center -> 12.2140534405943981255, Radius -> 1.629611331988023437133961124346×10^-8, Type -> Real|> x+<|Interpretation -> interpretation, Center -> -3.3746102964505553246, Radius -> 1.02775242656616416780934741837×10^-9, Type -> Real|> x^2+<|Interpretation -> interpretation, Center -> -1., Radius -> 0, Type -> Integer|> x^3
```

Find a random representative `*mrep*` of `*m*`:

```wolfram
ranrep[e_CenteredInterval]:=e["Center"]+RandomInteger[{-1000,1000}]/1000 e["Radius"]
(mrep=Map[ranrep,m,{2}])//MatrixForm
(* Output *)
({{-(10549599008221706530589)/(2305843009213693952000), (486382608439858608689)/(57646075230342348800), (2053645552326751460219)/(922337203685477580800)}, {-(263803966606585260959)/(288230376151711744000), -(2076096196815107441127)/(922337203685477580800), -(296389790752083792177)/(461168601842738790400)}, {(317192582610165443919)/(57646075230342348800), (10003993021260870377077)/(1152921504606846976000), (124351842812396123837)/(36028797018963968000)}})
```

Verify that the coefficients of `*p*` contain the coefficients of the characteristic polynomial of `*mrep*`:

```wolfram
MapThread[IntervalMemberQ,{CoefficientList[p,x],CoefficientList[CharacteristicPolynomial[mrep,x],x]}]
(* Output *)
{True,True,True,True}
```

#### Generalized Eigenvalues

The generalized characteristic polynomial of two matrices:

```wolfram
m1={{1,2},{5,4}};m2={{4,3},{6,4}};
```

```wolfram
CharacteristicPolynomial[{m1,m2},x]
(* Output *)
-6+7 x-2 x^2
```

This is equivalent to $|m_{1}-x m_{2}|$:

```wolfram
Det[m1-x m2]
(* Output *)
-6+7 x-2 x^2
```

A generalized characteristic polynomial of machine-precision matrices:

```wolfram
a={{1.,1.5,2.},{3.1,2.,2.9},{3.,2.,1.}};
```

```wolfram
b={{1.3,.5,1.1},{0.,1.5,2.3},{1.,0.,1.}};
```

```wolfram
CharacteristicPolynomial[{a,b},x]
(* Output *)
5.+5.9700000000000095 x-3.2800000000000047 x^2-1.4500000000000024 x^3
```

Find a generalized exact characteristic polynomial:

```wolfram
a={{1,1,1},{1,0,1},{0,0,1}};
```

```wolfram
b={{0,1,1},{0,1,1},{1,0,0}};
```

```wolfram
CharacteristicPolynomial[{a,b},x]
(* Output *)
-1-x+x^2
```

The absence of an $x^{3}$ term indicates an infinite generalized eigenvalue:

```wolfram
Eigenvalues[{a,b}]
(* Output *)
{∞,(1)/(2) (1+Sqrt[5]),(1)/(2) (1-Sqrt[5])}
```

Compute the result at finite precision:

```wolfram
CharacteristicPolynomial[N[{a,b},20],y]
(* Output *)
-1.-1. y+1. y^2+0`19.397940008672045 y^3
```

Find the generalized characteristic polynomial of symbolic matrices:

```wolfram
a={{x,1+x},{1-x,x}};b={{1,1},{1,2x}};
CharacteristicPolynomial[{a,b},y]
(* Output *)
-1+2 x^2+2 y-x y-2 x^2 y-y^2+2 x y^2
```

#### Special Matrices

Characteristic polynomials of sparse matrices:

```wolfram
SparseArray[{{1,3}->2,{2,2}->3,{3,1}->1,{4,2}->5},{4,4}]
(* Output *)
SparseArray[...]
```

```wolfram
CharacteristicPolynomial[%,x]
(* Output *)
6 x-2 x^2-3 x^3+x^4
```

```wolfram
SparseArray[{{x_,y_}/;Abs[x-y]<3->1},{10,10}]
(* Output *)
SparseArray[...]
```

```wolfram
CharacteristicPolynomial[%,x]
(* Output *)
1+2 x-15 x^2-20 x^3+61 x^4+28 x^5-76 x^6+28 x^8-10 x^9+x^10
```

Characteristic polynomials of structured matrices:

```wolfram
SymmetrizedArray[{{1,1}->2,{1,2}->1},{2,2},Symmetric[All]]
(* Output *)
SymmetrizedArray[...]
```

```wolfram
CharacteristicPolynomial[%,x]
(* Output *)
-1-2 x+x^2
```

```wolfram
QuantityArray[{{1,2},{3,4}},{"Meters","Meters"}]
(* Output *)
QuantityArray[...]
```

```wolfram
CharacteristicPolynomial[%,x]
(* Output *)
x^2+x (-4)+-2+x (-1)
```

The characteristic polynomial of an identity matrix is a binomial expansion:

```wolfram
CharacteristicPolynomial[IdentityMatrix[12],λ]
(* Output *)
1-12 λ+66 λ^2-220 λ^3+495 λ^4-792 λ^5+924 λ^6-792 λ^7+495 λ^8-220 λ^9+66 λ^10-12 λ^11+λ^12
```

```wolfram
Factor[%]
(* Output *)
(-1+λ)^12
```

Characteristic polynomial of [HilbertMatrix](https://reference.wolfram.com/language/ref/HilbertMatrix.html)[*n*]:

```wolfram
CharacteristicPolynomial[HilbertMatrix[5],λ]
(* Output *)
(1)/(266716800000)-(61501 λ)/(53343360000)+(852401 λ^2)/(222264000)-(735781 λ^3)/(2116800)+(563 λ^4)/(315)-λ^5
```

The characteristic polynomial of [JordanMatrix](https://reference.wolfram.com/language/ref/JordanMatrix.html)[λ,*n*] is $(\lambda-x)^{n}$:

```wolfram
CharacteristicPolynomial[JordanMatrix[λ,4],x]==(λ-x)^4//Simplify
(* Output *)
True
```

```wolfram
CharacteristicPolynomial[JordanMatrix[λ,5],x]==(λ-x)^5//Simplify
(* Output *)
True
```

The minimal polynomial of [CompanionMatrix](https://reference.wolfram.com/language/ref/CompanionMatrix.html)[{*c*_0,*c*_1,…,*c*_*n*}]:

```wolfram
CompanionMatrix[{c0,c1,c2,c3}]//MatrixForm
(* Output *)
({{0, 0, 0, -c0}, {1, 0, 0, -c1}, {0, 1, 0, -c2}, {0, 0, 1, -c3}})
```

```wolfram
CharacteristicPolynomial[%,x]
(* Output *)
c0+c1 x+c2 x^2+c3 x^3+x^4
```

### Applications

Find the characteristic polynomial of the matrix $m$ and compare the behavior for $a<32$, $a=32$ and $a>32$:

```wolfram
m=({{-6, 28, 21}, {4, -15, -12}, {-8, a, 25}})
(* Output *)
{{-6,28,21},{4,-15,-12},{-8,a,25}}
```

```wolfram
cp=CharacteristicPolynomial[m,x]
(* Output *)
-382+12 a+379 x-12 a x+4 x^2-x^3
```

Examining the roots, there is a root at $x=1$ independent of $a$:

```wolfram
roots=SolveValues[cp==0,x]
(* Output *)
{1,(1)/(2) (3-Sqrt[1537-48 a]),(1)/(2) (3+Sqrt[1537-48 a])}
```

For $a=32$ the root at $x=1$ is repeated:

```wolfram
roots/.a->32
(* Output *)
{1,1,2}
```

For $a<32$ there are three distinct real roots:

```wolfram
roots /. a->31
(* Output *)
{1,-2,5}
```

And for $a>32$, $x=1$ is the only real root, with the other two roots a complex conjugate pair:

```wolfram
roots /. a->33
(* Output *)
{1,(1)/(2) (3-ⅈ Sqrt[47]),(1)/(2) (3+ⅈ Sqrt[47])}
```

Visualize the three polynomials, zooming in on the bounce of the plot at the double root $x=1$:

```wolfram
{Plot[Table[cp,{a,{32,31,33}}]//Evaluate,{x,-3,6},PlotLegends->LineLegend[{a==32,a==31,a==33}],ImageSize->Small],
Plot[cp /. a->32,{x,0.5,2.5},PlotLegends->LineLegend[{a==32}],ImageSize->Small]}//Row
(* Output *)
![image](img/image_003.png)
```

Compute the determinant of a matrix as the constant term in its characteristic polynomial:

```wolfram
m={{6,5,8,1,4},{8,0,2,7,1},{7,4,0,8,3},{2,3,1,1,6},{0,4,8,2,3}};
```

```wolfram
cp=CharacteristicPolynomial[m,x]
(* Output *)
-8481-1710 x+638 x^2+148 x^3+10 x^4-x^5
```

Substitute in $x=0$:

```wolfram
cp /. x->0
(* Output *)
-8481
```

This result is also the product of the roots of the characteristic polynomial:

```wolfram
Times @@ SolveValues[cp==0,x] // FullSimplify
(* Output *)
-8481
```

Compare with a direct computation using [Det](https://reference.wolfram.com/language/ref/Det.html):

```wolfram
Det[m]
(* Output *)
-8481
```

Compute the trace of a matrix as the coefficient of the subleading power term in the characteristic polynomial:

```wolfram
m={{4,3,8,9},{8,1,0,7},{1,2,6,7},{2,5,1,2}};
```

```wolfram
cp=CharacteristicPolynomial[m,x]
(* Output *)
-568-30 x-36 x^2-13 x^3+x^4
```

Extract the coefficient of $(-x)^{n-1}$, where $n$ is the height or width of the matrix:

```wolfram
(-1)^(Length[m]-1)Coefficient[cp,x^(Length[m]-1)]
(* Output *)
13
```

This result is also the sum of the roots of the characteristic polynomial:

```wolfram
Total[SolveValues[cp==0,x]]//FullSimplify
(* Output *)
13
```

Compare with a direct computation using [Tr](https://reference.wolfram.com/language/ref/Tr.html):

```wolfram
Tr[m]
(* Output *)
13
```

Find the eigenvalues of a matrix as the roots of the characteristic polynomial:

```wolfram
m={{1,2,3},{4,5,6},{7,8,9}};
```

```wolfram
SolveValues[CharacteristicPolynomial[m,x]==0,x]
(* Output *)
{0,(3)/(2) (5-Sqrt[33]),(3)/(2) (5+Sqrt[33])}
```

Compare with a direct computation using [Eigenvalues](https://reference.wolfram.com/language/ref/Eigenvalues.html):

```wolfram
Eigenvalues[m]
(* Output *)
{(3)/(2) (5+Sqrt[33]),(3)/(2) (5-Sqrt[33]),0}
```

Use the characteristic polynomial to find the eigenvalues and eigenvectors of the matrices $a$ and $a$:

```wolfram
a={{3,9,4},{3,7,1},{8,1,8}};
```

```wolfram
cp=CharacteristicPolynomial[a,x]
(* Output *)
-191-41 x+18 x^2-x^3
```

The two matrices have the same characteristic polynomial:

```wolfram
CharacteristicPolynomial[Transpose[a],x]==cp
(* Output *)
True
```

Thus, they will both have the same eigenvalues, which are the roots of the polynomial:

```wolfram
λ=SolveValues[cp==0,x]
(* Output *)
{<|icon -> Root, small -> "-2.22", approx -> -2.2224141172805587, interp -> Root[191, +, 41, #, -, 18, #, ^, 2, +, #, ^, 3, &, 1, 0], head -> Root, big -> 191+41 #1-18 #1^2+#1^3&, degree -> 3, shortDegree -> 3, number -> 1|>,<|icon -> Root, small -> "6.07", approx -> 6.074633070347852, interp -> Root[191, +, 41, #, -, 18, #, ^, 2, +, #, ^, 3, &, 2, 0], head -> Root, big -> 191+41 #1-18 #1^2+#1^3&, degree -> 3, shortDegree -> 3, number -> 2|>,<|icon -> Root, small -> "14.1", approx -> 14.147781046932707, interp -> Root[191, +, 41, #, -, 18, #, ^, 2, +, #, ^, 3, &, 3, 0], head -> Root, big -> 191+41 #1-18 #1^2+#1^3&, degree -> 3, shortDegree -> 3, number -> 3|>}
```

The eigenvectors are given by the null space of $a-\lambda I$:

```wolfram
v_λ=Join[NullSpace[a-λ[[1]]IdentityMatrix[3]],NullSpace[a-λ[[2]]IdentityMatrix[3]],NullSpace[a-λ[[3]]IdentityMatrix[3]]]//FullSimplify
(* Output *)
{{<|icon -> Root, small -> "-1.32", approx -> -1.317833310267928, interp -> Root[-101, -, 504, #, +, 473, #, ^, 2, +, 605, #, ^, 3, &, 1, 0], head -> Root, big -> -101-504 #1+473 #1^2+605 #1^3&, degree -> 3, shortDegree -> 3, number -> 1|>,<|icon -> Root, small -> "0.320", approx -> 0.3202523648628665, interp -> Root[43, -, 147, #, -, 154, #, ^, 2, +, 605, #, ^, 3, &, 2, 0], head -> Root, big -> 43-147 #1-154 #1^2+605 #1^3&, degree -> 3, shortDegree -> 3, number -> 2|>,1},{<|icon -> Root, small -> "-0.178", approx -> -0.17753403275884397, interp -> Root[-101, -, 504, #, +, 473, #, ^, 2, +, 605, #, ^, 3, &, 2, 0], head -> Root, big -> -101-504 #1+473 #1^2+605 #1^3&, degree -> 3, shortDegree -> 3, number -> 2|>,<|icon -> Root, small -> "-0.505", approx -> -0.5050946675813957, interp -> Root[43, -, 147, #, -, 154, #, ^, 2, +, 605, #, ^, 3, &, 1, 0], head -> Root, big -> 43-147 #1-154 #1^2+605 #1^3&, degree -> 3, shortDegree -> 3, number -> 1|>,1},{<|icon -> Root, small -> "0.714", approx -> 0.7135491612085902, interp -> Root[-101, -, 504, #, +, 473, #, ^, 2, +, 605, #, ^, 3, &, 3, 0], head -> Root, big -> -101-504 #1+473 #1^2+605 #1^3&, degree -> 3, shortDegree -> 3, number -> 3|>,<|icon -> Root, small -> "0.439", approx -> 0.43938775726398366, interp -> Root[43, -, 147, #, -, 154, #, ^, 2, +, 605, #, ^, 3, &, 3, 0], head -> Root, big -> 43-147 #1-154 #1^2+605 #1^3&, degree -> 3, shortDegree -> 3, number -> 3|>,1}}
```

[Eigensystem](https://reference.wolfram.com/language/ref/Eigensystem.html) gives the same result, though it sorts eigenvalues by absolute value:

```wolfram
Eigensystem[a]
(* Output *)
{{<|icon -> Root, small -> "14.1", approx -> 14.147781046932707, interp -> Root[191, +, 41, #, -, 18, #, ^, 2, +, #, ^, 3, &, 3, 0], head -> Root, big -> 191+41 #1-18 #1^2+#1^3&, degree -> 3, shortDegree -> 3, number -> 3|>,<|icon -> Root, small -> "6.07", approx -> 6.074633070347852, interp -> Root[191, +, 41, #, -, 18, #, ^, 2, +, #, ^, 3, &, 2, 0], head -> Root, big -> 191+41 #1-18 #1^2+#1^3&, degree -> 3, shortDegree -> 3, number -> 2|>,<|icon -> Root, small -> "-2.22", approx -> -2.2224141172805587, interp -> Root[191, +, 41, #, -, 18, #, ^, 2, +, #, ^, 3, &, 1, 0], head -> Root, big -> 191+41 #1-18 #1^2+#1^3&, degree -> 3, shortDegree -> 3, number -> 1|>},{{<|icon -> Root, small -> "0.714", approx -> 0.7135491612085902, interp -> Root[-101, -, 504, #, +, 473, #, ^, 2, +, 605, #, ^, 3, &, 3, 0], head -> Root, big -> -101-504 #1+473 #1^2+605 #1^3&, degree -> 3, shortDegree -> 3, number -> 3|>,<|icon -> Root, small -> "0.439", approx -> 0.43938775726398366, interp -> Root[43, -, 147, #, -, 154, #, ^, 2, +, 605, #, ^, 3, &, 3, 0], head -> Root, big -> 43-147 #1-154 #1^2+605 #1^3&, degree -> 3, shortDegree -> 3, number -> 3|>,1},{<|icon -> Root, small -> "-0.178", approx -> -0.17753403275884397, interp -> Root[-101, -, 504, #, +, 473, #, ^, 2, +, 605, #, ^, 3, &, 2, 0], head -> Root, big -> -101-504 #1+473 #1^2+605 #1^3&, degree -> 3, shortDegree -> 3, number -> 2|>,<|icon -> Root, small -> "-0.505", approx -> -0.5050946675813957, interp -> Root[43, -, 147, #, -, 154, #, ^, 2, +, 605, #, ^, 3, &, 1, 0], head -> Root, big -> 43-147 #1-154 #1^2+605 #1^3&, degree -> 3, shortDegree -> 3, number -> 1|>,1},{<|icon -> Root, small -> "-1.32", approx -> -1.317833310267928, interp -> Root[-101, -, 504, #, +, 473, #, ^, 2, +, 605, #, ^, 3, &, 1, 0], head -> Root, big -> -101-504 #1+473 #1^2+605 #1^3&, degree -> 3, shortDegree -> 3, number -> 1|>,<|icon -> Root, small -> "0.320", approx -> 0.3202523648628665, interp -> Root[43, -, 147, #, -, 154, #, ^, 2, +, 605, #, ^, 3, &, 2, 0], head -> Root, big -> 43-147 #1-154 #1^2+605 #1^3&, degree -> 3, shortDegree -> 3, number -> 2|>,1}}}
```

While $a$ has the same eigenvalues as $a$, it has different eigenvectors:

```wolfram
v_λ^*=Join[NullSpace[aᵀ-λ[[1]]IdentityMatrix[3]], NullSpace[aᵀ-λ[[2]]IdentityMatrix[3]], NullSpace[aᵀ-λ[[3]]IdentityMatrix[3]]]//FullSimplify
(* Output *)
{{<|icon -> Root, small -> "-3.34", approx -> -3.3444419213684182, interp -> Root[31, -, 228, #, +, 113, #, ^, 2, +, 55, #, ^, 3, &, 1, 0], head -> Root, big -> 31-228 #1+113 #1^2+55 #1^3&, degree -> 3, shortDegree -> 3, number -> 1|>,<|icon -> Root, small -> "3.16", approx -> 3.1553535681931146, interp -> Root[689, -, 381, #, -, 122, #, ^, 2, +, 55, #, ^, 3, &, 3, 0], head -> Root, big -> 689-381 #1-122 #1^2+55 #1^3&, degree -> 3, shortDegree -> 3, number -> 3|>,1},{<|icon -> Root, small -> "0.148", approx -> 0.14752592654596802, interp -> Root[31, -, 228, #, +, 113, #, ^, 2, +, 55, #, ^, 3, &, 2, 0], head -> Root, big -> 31-228 #1+113 #1^2+55 #1^3&, degree -> 3, shortDegree -> 3, number -> 2|>,<|icon -> Root, small -> "-2.52", approx -> -2.5154706358360195, interp -> Root[689, -, 381, #, -, 122, #, ^, 2, +, 55, #, ^, 3, &, 1, 0], head -> Root, big -> 689-381 #1-122 #1^2+55 #1^3&, degree -> 3, shortDegree -> 3, number -> 1|>,1},{<|icon -> Root, small -> "1.14", approx -> 1.1423705402769957, interp -> Root[31, -, 228, #, +, 113, #, ^, 2, +, 55, #, ^, 3, &, 3, 0], head -> Root, big -> 31-228 #1+113 #1^2+55 #1^3&, degree -> 3, shortDegree -> 3, number -> 3|>,<|icon -> Root, small -> "1.58", approx -> 1.5782988858247229, interp -> Root[689, -, 381, #, -, 122, #, ^, 2, +, 55, #, ^, 3, &, 2, 0], head -> Root, big -> 689-381 #1-122 #1^2+55 #1^3&, degree -> 3, shortDegree -> 3, number -> 2|>,1}}
```

Visualize the two sets of eigenvectors:

```wolfram
{MatrixPlot[v_λ],MatrixPlot[v_λ^*]}
(* Output *)
{[Graphics],[Graphics]}
```

Find the generalized eigensystem of $a$ with respect to $b$ as the roots of the characteristic polynomial:

```wolfram
a={{1.,1.5,2.},{3.1,2.,2.9},{3.,2.,1.}};
```

```wolfram
b={{1.3,.5,1.1},{0.,1.5,2.3},{1.,0.,1.}};
```

```wolfram
CharacteristicPolynomial[{a,b},x]
(* Output *)
5.+5.9700000000000095 x-3.2800000000000047 x^2-1.4500000000000024 x^3
```

The roots of the generalized characteristic polynomial are the generalized eigenvalues:

```wolfram
λ=SolveValues[%==0,x]
(* Output *)
{-3.210040153039578,-0.6656977598648105,1.6136689473871475}
```

The generalized eigenvectors are given by the null space of $a-\lambda b$:

```wolfram
v_λ=Join[NullSpace[a-λ[[1]]b],NullSpace[a-λ[[2]]b],NullSpace[a-λ[[3]]b]]//FullSimplify
(* Output *)
{{-0.12925008940601218,-0.8081783372442966,0.5745800114845468},{0.30511503726232764,-0.8740704652144278,0.3780286178006138},{0.26560296832774916,0.10978275656721893,0.9578114687014305}}
```

Compare with a direct computation using [Eigensystem](https://reference.wolfram.com/language/ref/Eigensystem.html):

```wolfram
Eigensystem[{a,b}]
(* Output *)
{{-3.210040153039578,1.6136689473871475,-0.6656977598648105},{{-0.12925008940601218,-0.8081783372442966,0.5745800114845468},{-0.26560296832774916,-0.10978275656721843,-0.9578114687014306},{-0.3051150372623274,0.8740704652144277,-0.3780286178006139}}}
```

### Properties & Relations

The characteristic polynomial is equivalent to [Det](https://reference.wolfram.com/language/ref/Det.html)[*m*-*id* *x*]:

```wolfram
m=RandomInteger[9,{10,10}];
```

```wolfram
cp=CharacteristicPolynomial[m,x]
(* Output *)
307960090+150589872 x-30349953 x^2-7786386 x^3+2056940 x^4+6398 x^5-36248 x^6+2905 x^7+12 x^8-49 x^9+x^10
```

```wolfram
%==Det[m-IdentityMatrix[Length[m]] x]
(* Output *)
True
```

The generalized characteristic polynomial is equivalent to [Det](https://reference.wolfram.com/language/ref/Det.html)[*m*-*a* *x*]:

```wolfram
{m,a} = RandomInteger[9,{2,5,5}];
```

```wolfram
CharacteristicPolynomial[{m,a},x]
(* Output *)
-176+15012 x-20780 x^2+7604 x^3-4680 x^4-1880 x^5
```

```wolfram
%==Det[m - a x]
(* Output *)
-783+79 x-2671 x^2+16020 x^3-26327 x^4+9390 x^5
```

A matrix is a root of its characteristic polynomial (Cayley-Hamilton theorem [[more...](http://mathworld.wolfram.com/Cayley-HamiltonTheorem.html)]):

```wolfram
m=RandomInteger[10,{10,10}];
cp=CharacteristicPolynomial[m,x];
MatrixPolynomialValue[cp,m,x]==ConstantArray[0,{10,10}]
(* Output *)
True
```

```wolfram
MatrixPolynomialValue[cp,m,x]//AbsoluteTiming
(* Output *)
{0.000302,{{0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0}}}
```

$\prod_{\mathit{i}=1}^{\mathit{n}}(\mathit{v}_{\mathit{i}}-\mathit{x})$ where $v_{i}$ are the eigenvalues is equivalent to the characteristic polynomial:

```wolfram
m=RandomReal[1, {10,10}];
```

```wolfram
Product[(v-x),{v,Eigenvalues[m]}]-CharacteristicPolynomial[m,x]//Simplify//Chop
(* Output *)
0
```

The sum of the roots of the characteristic polynomial is the trace ([Tr](https://reference.wolfram.com/language/ref/Tr.html)) of the matrix:

```wolfram
m=RandomReal[1, {10,10}];
```

```wolfram
roots = SolveValues[CharacteristicPolynomial[m,x]==0,x];
```

```wolfram
Total[roots]==Tr[m]
(* Output *)
True
```

Similarly, the product of the roots is the determinant ([Det](https://reference.wolfram.com/language/ref/Det.html)):

```wolfram
Times @@ roots == Det[m]
(* Output *)
True
```

A matrix and its transpose have the same characteristic polynomial:

```wolfram
m=RandomReal[1, {10,10}];
```

```wolfram
CharacteristicPolynomial[m,x]==CharacteristicPolynomial[Transpose[m],x]
(* Output *)
True
```

All triangular matrices with a common diagonal have the same characteristic polynomial:

```wolfram
t1=({{1, 0, 0}, {0, 1, 0}, {0, 0, 1}});t2=({{1, 4, 0}, {0, 1, 3}, {0, 0, 1}});t3=({{1, 0, 0}, {5, 1, 0}, {6, 7, 1}});
```

```wolfram
CharacteristicPolynomial[t1,x]==CharacteristicPolynomial[t2,x]==CharacteristicPolynomial[t3,x]
(* Output *)
True
```

If $p(x)$ is a monic polynomial, then the characteristic polynomial of its companion matrix is $p(x)$:

```wolfram
cl=RandomInteger[9,10];
p=cl.x^Range[0,9]+x^10
(* Output *)
5+8 x+9 x^2+3 x^3+7 x^4+3 x^5+4 x^7+x^8+2 x^9+x^10
```

```wolfram
CharacteristicPolynomial[CompanionMatrix[cl],x]==p
(* Output *)
True
```

[MatrixMinimalPolynomial](https://reference.wolfram.com/language/ref/MatrixMinimalPolynomial.html)[*m*,λ] divides [CharacteristicPolynomial](https://reference.wolfram.com/language/ref/CharacteristicPolynomial.html)[*m*,λ] with a (possibly constant) polynomial quotient:

```wolfram
mat=({{0, -1, 1}, {1, 2, -1}, {1, 1, 0}});
PolynomialQuotientRemainder[CharacteristicPolynomial[mat,t],MatrixMinimalPolynomial[mat,t],t]
(* Output *)
{1-t,0}
```

The minimal polynomial and the characteristic polynomial have the same distinct roots:

```wolfram
a=({{-6, 4, 0, 9}, {-3, 0, 1, 6}, {-1, -2, 1, 0}, {-4, 4, 0, 7}});
χ=CharacteristicPolynomial[a,λ];
μ=MatrixMinimalPolynomial[a,λ];
Roots[μ==0,λ]==DeleteDuplicates[Roots[χ==0,λ]]
(* Output *)
True
```

In particular, this means $\chi$ divides $\mu$ raised to the $lcm(deg(\mu),deg(\chi))$ power:

```wolfram
PolynomialRemainder[μ^LCM@@(Max@Cases[#,λ^k_:>k,-1]&/@{μ,χ}),χ,λ]
(* Output *)
0
```

## Tech Notes ▪Basic Matrix Operations ▪Eigenvalues and Eigenvectors

## Related Guides ▪Matrix Operations ▪Finite Fields

## History Introduced in 2003 (5.0) | Updated in 2007 (6.0) ▪ 2024 (14.0)
