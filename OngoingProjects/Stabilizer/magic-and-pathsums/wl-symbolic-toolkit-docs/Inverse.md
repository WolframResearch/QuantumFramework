# Inverse | [SpanFromLeft]

> [Inverse](https://reference.wolfram.com/language/ref/Inverse.html)[*m*] — gives the inverse of a square matrix `*m*`.

## Details and Options

[Inverse](https://reference.wolfram.com/language/ref/Inverse.html) works on both symbolic and numerical matrices.

For matrices with approximate real or complex numbers, the inverse is generated to the maximum possible precision given the input. A warning is given for ill[Hyphen]conditioned matrices.

The following options can be given

| [Method](https://reference.wolfram.com/language/ref/Method.html) | [Automatic](https://reference.wolfram.com/language/ref/Automatic.html) | method to use |
| --- | --- | --- |
| [Modulus](https://reference.wolfram.com/language/ref/Modulus.html) | 0 | integer modulus to use |
| [ZeroTest](https://reference.wolfram.com/language/ref/ZeroTest.html) | [Automatic](https://reference.wolfram.com/language/ref/Automatic.html) | test for whether an element is zero |

[Inverse](https://reference.wolfram.com/language/ref/Inverse.html)[*m*,[Modulus](https://reference.wolfram.com/language/ref/Modulus.html)->*p*] evaluates the inverse of rational matrices modulo the integer `*n*`. If `*n*` is zero, ordinary arithmetic is used. If `*n*` is not prime, the computation may fail.

[Inverse](https://reference.wolfram.com/language/ref/Inverse.html)[*m*,[ZeroTest](https://reference.wolfram.com/language/ref/ZeroTest.html)->*test*] evaluates `*test*[*m*[[*i*,*j*]]]` to determine whether matrix elements are zero.

Settings of [Method](https://reference.wolfram.com/language/ref/Method.html) applicable to exact and symbolic matrices include `"CofactorExpansion"`, `"DivisionFreeRowReduction"`, and `"OneStepRowReduction"`. The default setting of [Automatic](https://reference.wolfram.com/language/ref/Automatic.html) switches among these methods depending on the matrix given.

[Inverse](https://reference.wolfram.com/language/ref/Inverse.html)[*m*] formats as `*m*` in [StandardForm](https://reference.wolfram.com/language/ref/StandardForm.html) and [](https://reference.wolfram.com/language/ref/.html).

## Examples

### Basic Examples

Inverse of a 2×2 matrix:

```wolfram
Inverse[{{1.4,2},{3,-6.7}}]
(* Output *)
{{0.43563068920676207,0.13003901170351104},{0.1950585175552666,-0.09102730819245773}}
```

Enter the matrix in a grid:

```wolfram
Inverse[({{1, 2, 3}, {4, 2, 2}, {5, 1, 7}})]
(* Output *)
{{-(2)/(7),(11)/(42),(1)/(21)},{(3)/(7),(4)/(21),-(5)/(21)},{(1)/(7),-(3)/(14),(1)/(7)}}
```

Inverse of a symbolic matrix:

```wolfram
Inverse[{{u,v},{v,u}}]
(* Output *)
{{(u)/(u^2-v^2),-(v)/(u^2-v^2)},{-(v)/(u^2-v^2),(u)/(u^2-v^2)}}
```

### Scope

#### Basic Uses

Find the inverse of a machine-precision matrix:

```wolfram
Inverse[{{1.2,2.5,-3.2},{0.7,-9.4,5.8},{-0.2,0.3,6.4}}]
(* Output *)
{{0.7454598005684282,0.20424875957416058,0.18762946192013105},{0.06792234693386,-0.08478250397417986,0.1107953176935305},{0.02011175875523869,0.010356953610482198,0.15691989016811983}}
```

Invert a complex matrix:

```wolfram
Inverse[{{1.+I,2,3-2I},{0,4π,5I},{ℯ,0,6}}]
(* Output *)
{{-0.06819300426548888-0.4303810422940135 ⅈ,0.01085325371313926+0.06849727029413433 ⅈ,0.23463790814252758+0.18341514163089442 ⅈ},{0.07758120214149704-0.01229258431170088 ⅈ,0.06723003973411688+0.001956425556581079 ⅈ,-0.033062718336364-0.024018340242081278 ⅈ},{0.030894634053818096+0.19498282776351344 ⅈ,-0.004917033724680351-0.031032480856598814 ⅈ,0.06036467300475411-0.08309567442658373 ⅈ}}
```

Inverse of an exact matrix:

```wolfram
Inverse[{{2,3,2},{4,9,2},{7,2,4}}]
(* Output *)
{{-(8)/(13),(2)/(13),(3)/(13)},{(1)/(26),(3)/(26),-(1)/(13)},{(55)/(52),-(17)/(52),-(3)/(26)}}
```

Inverse of an arbitrary-precision matrix:

```wolfram
Inverse[RandomReal[1,{2,2},WorkingPrecision->20]]
(* Output *)
{{-0.41246390063006539541877763713524680758,1.33765080599459824958158407750181734448},{1.71369835877580405609849322563797233781,-1.35333364739239005972567809739194934816}}
```

Inverse of a symbolic matrix:

```wolfram
Inverse[{{a,b},{c,d}}]
(* Output *)
{{(d)/(-b c+a d),-(b)/(-b c+a d)},{-(c)/(-b c+a d),(a)/(-b c+a d)}}
```

Verifying a symbolic inverse may require simplification:

```wolfram
%.{{a,b},{c,d}}
(* Output *)
{{-(b c)/(-b c+a d)+(a d)/(-b c+a d),0},{0,-(b c)/(-b c+a d)+(a d)/(-b c+a d)}}
```

```wolfram
Simplify[%]
(* Output *)
{{1,0},{0,1}}
```

The inversion of large machine-precision matrices is efficient:

```wolfram
a=RandomReal[{0,10},{800,800}];
```

```wolfram
Inverse[a];//Timing
(* Output *)
{0.365246,Null}
```

Inverse of a matrix over a finite field:

```wolfram
ℱ=FiniteField[23,4];
Inverse[{{ℱ[12],ℱ[23],ℱ[34]},{ℱ[45],ℱ[56],ℱ[67]},{ℱ[78],ℱ[89],ℱ[90]}}]//MatrixForm
(* Output *)
![image](img/image_001.png)
```

Inverse of a [CenteredInterval](https://reference.wolfram.com/language/ref/CenteredInterval.html) matrix:

```wolfram
(m=Map[CenteredInterval,RandomReal[{-10,10},{3,3},WorkingPrecision->10],{2}])//MatrixForm
(* Output *)
({{<|Interpretation -> interpretation, Center -> 4.7654618904925882816, Radius -> 4.7654618966219697284714129637×10^-10, Type -> Real|>, <|Interpretation -> interpretation, Center -> -3.3854539901949465275, Radius -> 3.3854539917971271378860365076×10^-10, Type -> Real|>, <|Interpretation -> interpretation, Center -> -0.4094120825175195932, Radius -> 4.094120826863764661673883438×10^-11, Type -> Real|>}, {<|Interpretation -> interpretation, Center -> 5.4861263034399598837, Radius -> 5.4861263101835255895366572076×10^-10, Type -> Real|>, <|Interpretation -> interpretation, Center -> -8.3527811791282147169, Radius -> 8.3527811837319498877718615404×10^-10, Type -> Real|>, <|Interpretation -> interpretation, Center -> 4.5083638816140592098, Radius -> 4.5083638818693705374585078971×10^-10, Type -> Real|>}, {<|Interpretation -> interpretation, Center -> 0.7107301882933825254, Radius -> 7.107301888947120671602419861×10^-11, Type -> Real|>, <|Interpretation -> interpretation, Center -> -1.3443382747936993837, Radius -> 1.3443382756821165013860763793×10^-10, Type -> Real|>, <|Interpretation -> interpretation, Center -> -3.7522752629593014717, Radius -> 3.7522752656740654408906721073×10^-10, Type -> Real|>}})
```

```wolfram
(minv=Inverse[m])//MatrixForm
(* Output *)
({{<|Interpretation -> interpretation, Center -> 0.3805293000586971175, Radius -> 4.2221103609824117697257861437×10^-10, Type -> Real|>, <|Interpretation -> interpretation, Center -> -0.1236403749939825047, Radius -> 4.353050272673519849320200592×10^-10, Type -> Real|>, <|Interpretation -> interpretation, Center -> -0.1900737670203369589, Radius -> 5.7354282781418852721344592283×10^-10, Type -> Real|>}, {<|Interpretation -> interpretation, Center -> 0.2420325968486025658, Radius -> 2.0599072534516205124077714572×10^-10, Type -> Real|>, <|Interpretation -> interpretation, Center -> -0.1789614371330543463, Radius -> 2.304508346304134924764639436×10^-10, Type -> Real|>, <|Interpretation -> interpretation, Center -> -0.2414306748302976757, Radius -> 2.9466502675379313558323701727×10^-10, Type -> Real|>}, {<|Interpretation -> interpretation, Center -> -0.0146364588844383192, Radius -> 3.937907823176739485759867421×10^-11, Type -> Real|>, <|Interpretation -> interpretation, Center -> 0.0406979104537938952, Radius -> 5.044200517348440548914823012×10^-11, Type -> Real|>, <|Interpretation -> interpretation, Center -> -0.2160093837841685627, Radius -> 6.243180972198558720265282318×10^-11, Type -> Real|>}})
```

Find a random representative `*mrep*` of `*m*`:

```wolfram
ranrep[e_CenteredInterval]:=e["Center"]+RandomInteger[{-1000,1000}]/1000e["Radius"]
(mrep=Map[ranrep,m,{2}])//MatrixForm
(* Output *)
({{(54942034925388205087)/(11529215046068469760), -(487895338557742750829)/(144115188075855872000), -(7552319906422183269517)/(18446744073709551616000)}, {(63250729920727638207)/(11529215046068469760), -(481505052205144626787)/(57646075230342348800), (5197789669818438115839)/(1152921504606846976000)}, {(327766447230420984161)/(461168601842738790400), -(1549916506611497486057)/(1152921504606846976000), -(8652157683249603138977)/(2305843009213693952000)}})
```

Verify that `*minv*` contains the inverse of `*mrep*`:

```wolfram
MapThread[IntervalMemberQ,{minv,Inverse[mrep]},2]//MatrixForm
(* Output *)
({{True, True, True}, {True, True, True}, {True, True, True}})
```

Output formatting:

```wolfram
Inverse[m]
(* Output *)
m
```

#### Special Matrices

The inverse of a sparse matrix is returned as a normal matrix:

```wolfram
SparseArray[{{1,3}->1,{2,2}->2,{3,1}->3},{3,3}]
(* Output *)
SparseArray[...]
```

```wolfram
Inverse[%]
(* Output *)
{{0,0,(1)/(3)},{0,(1)/(2),0},{1,0,0}}
```

Format the result:

```wolfram
%//MatrixForm
(* Output *)
({{0, 0, (1)/(3)}, {0, (1)/(2), 0}, {1, 0, 0}})
```

When possible, the inverse of a structured matrix is returned as another structured matrix:

```wolfram
QuantityArray[{{1,2},{3,4}},{"Meters","Seconds"}]
(* Output *)
QuantityArray[...]
```

```wolfram
Inverse[%]
(* Output *)
QuantityArray[...]
```

This is not always possible:

```wolfram
SymmetrizedArray[{{1,1}->3,{2,2}->1,{3,1}->-5},{3,3},Symmetric[All]]
(* Output *)
SymmetrizedArray[...]
```

```wolfram
Inverse[%]
(* Output *)
{{0,0,-(1)/(5)},{0,1,0},{-(1)/(5),0,-(3)/(25)}}
```

[IdentityMatrix](https://reference.wolfram.com/language/ref/IdentityMatrix.html) is its own inverse:

```wolfram
Inverse[IdentityMatrix[3]]
(* Output *)
{{1,0,0},{0,1,0},{0,0,1}}
```

Inverse of [HilbertMatrix](https://reference.wolfram.com/language/ref/HilbertMatrix.html):

```wolfram
Inverse[HilbertMatrix[3]]//MatrixForm
(* Output *)
({{9, -36, 30}, {-36, 192, -180}, {30, -180, 180}})
```

Visualize the inverses for several matrix sizes:

```wolfram
Table[ArrayPlot[Inverse[HilbertMatrix[n]],Mesh->True,ColorFunction->"Rainbow"],{n,{3,7,15,25}}]
(* Output *)
{[Graphics],[Graphics],[Graphics],[Graphics]}
```

Compute the inverse of a $10 \times 10$ matrix of univariate polynomials of degree $100$:

```wolfram
rpoly[n_]:=RandomInteger[{-2^10,2^10},{n+1}].x^Range[0,n]
SeedRandom[1234];
m=Table[rpoly[100],{10},{10}];
```

```wolfram
Inverse[m]//Short[#,4]&//AbsoluteTiming
(* Output *)
{2.4998049,{{(<<1335>>+822285305206069049790162945 x^900)/(3034638622673506403581743718005+<<1485>>+998816316987602849360464925610 x^1000),(<<1>>)/(<<1>>),(<<1>>)/(<<1>>),(<<1>>)/(<<1>>),(<<1>>)/(<<1>>),(<<1>>)/(<<1>>),(<<1>>)/(<<1>>),(<<1>>)/(<<1>>),(<<1>>)/(3034638622673506403581743718005+<<1485>>+<<1>>),(<<1>>)/(<<1>>)},<<9>>}}
```

### Options

#### Method

The [Automatic](https://reference.wolfram.com/language/ref/Automatic.html) method attempts to use a method appropriate for a given matrix:

```wolfram
h=N[HilbertMatrix[5],10];
```

Inverse[h,Method->Automatic]//MatrixForm

```
(* Output *)
({{25., -300., 1050., -1400., 630.}, {-300., 4800., -18900., 26880., -12600.}, {1050., -18900., 79380., -117600., 56700.}, {-1400., 26880., -117600., 179200., -88200.}, {630., -12600., 56700., -88200., 44100.}})
```

Here, the [Automatic](https://reference.wolfram.com/language/ref/Automatic.html) method has avoided numerical instability possible for some methods when applied to an ill-conditioned matrix:

```wolfram
Inverse[h,Method->"CofactorExpansion"]
(* Output *)
Inverse
(* Output *)
Inverse[{{1.,0.5,0.33333333333333333333333333333333333333,0.25,0.2},{0.5,0.33333333333333333333333333333333333333,0.25,0.2,0.16666666666666666666666666666666666667},{0.33333333333333333333333333333333333333,0.25,0.2,0.16666666666666666666666666666666666667,0.14285714285714285714285714285714285714},{0.25,0.2,0.16666666666666666666666666666666666667,0.14285714285714285714285714285714285714,0.125},{0.2,0.16666666666666666666666666666666666667,0.14285714285714285714285714285714285714,0.125,0.11111111111111111111111111111111111111}},Method->"CofactorExpansion"]
```

Generate a random numeric matrix with small off-diagonal elements and compute its inverse using various methods:

```wolfram
SeedRandom[1111];
matrix=RandomReal[{0,10^-10},{4,4}]+DiagonalMatrix[{1,1,1,1}]
(* Output *)
{{1.0000000000077796,5.562659362959752×10^-11,8.076816868335859×10^-11,5.2798284115799055×10^-11},{5.558767623869201×10^-11,1.000000000098109,9.075966782672407×10^-11,3.643901151484488×10^-11},{7.22450567447879×10^-11,4.201865534344489×10^-11,1.0000000000156009,1.300203320896276×10^-11},{2.6639760381147483×10^-11,8.798874749233894×10^-11,1.853749394363635×10^-11,1.0000000000982971}}
```

```wolfram
inverses=Table[Inverse[matrix,Method->meth],{meth,{Automatic,"CofactorExpansion","DivisionFreeRowReduction","OneStepRowReduction"}}];
```

The results are very close but can differ slightly in the least significant digits due to numerical precision handling:

```wolfram
MatrixForm[inverses[[1]]-#]&/@Rest[inverses]
(* Output *)
{({{0., 0., 0., 0.}, {0., 0., 0., 0.}, {-1.2924697071141057×10^-26, -6.462348535570529×10^-27, 0., 1.6155871338926322×10^-27}, {3.2311742677852644×10^-27, 0., 0., 0.}}),({{0., 0., 0., 0.}, {0., 0., 0., 0.}, {-1.2924697071141057×10^-26, 0., 0., 0.}, {3.2311742677852644×10^-27, 0., 0., 0.}}),({{0., 0., 0., 0.}, {0., 0., 0., 0.}, {-1.2924697071141057×10^-26, 0., 0., 0.}, {3.2311742677852644×10^-27, 0., 0., 0.}})}
```

Compare timings using different methods on a matrix with exact numeric entries:

```wolfram
matrix=HilbertMatrix[11];
```

```wolfram
inverses=Table[RepeatedTiming@Inverse[matrix,Method->meth],{meth,{Automatic,"CofactorExpansion","DivisionFreeRowReduction","OneStepRowReduction"}}];
```

```wolfram
inverses[[All,1]]
(* Output *)
{0.0005844140625,0.0843365,0.0004314990234375,0.0003071240234375}
```

Create a matrix with several symbolic elements:

```wolfram
matrix=HilbertMatrix[7];
matrix[[1,{1,2}]]=a;
matrix[[2,{3,4}]]={b,c};
matrix[[{3,4},2]]={c,b};
matrix[[4,4]]=d;
MatrixForm[matrix]
(* Output *)
({{a, a, (1)/(3), (1)/(4), (1)/(5), (1)/(6), (1)/(7)}, {(1)/(2), (1)/(3), b, c, (1)/(6), (1)/(7), (1)/(8)}, {(1)/(3), c, (1)/(5), (1)/(6), (1)/(7), (1)/(8), (1)/(9)}, {(1)/(4), b, (1)/(6), d, (1)/(8), (1)/(9), (1)/(10)}, {(1)/(5), (1)/(6), (1)/(7), (1)/(8), (1)/(9), (1)/(10), (1)/(11)}, {(1)/(6), (1)/(7), (1)/(8), (1)/(9), (1)/(10), (1)/(11), (1)/(12)}, {(1)/(7), (1)/(8), (1)/(9), (1)/(10), (1)/(11), (1)/(12), (1)/(13)}})
```

Compare timings and result sizes using the different methods on this matrix:

```wolfram
inverses=Table[Timing[Inverse[matrix,Method->meth]],{meth,{"CofactorExpansion","DivisionFreeRowReduction","OneStepRowReduction"}}];
Map[{#[[1]],LeafCount[#[[2]]]}&,inverses]
(* Output *)
{{0.014489,9740},{0.461269,22386},{0.15414,6671}}
```

#### Modulus

Invert a matrix using arithmetic modulo five:

```wolfram
m={{1,2,4},{2,0,3},{2,1,2}};
```

```wolfram
Inverse[m,Modulus->5]
(* Output *)
{{3,0,4},{3,1,0},{3,2,4}}
```

The inverse using normal arithmetic:

```wolfram
Inverse[m]
(* Output *)
{{-(1)/(3),0,(2)/(3)},{(2)/(9),-(2)/(3),(5)/(9)},{(2)/(9),(1)/(3),-(4)/(9)}}
```

Visualize the two results:

```wolfram
{MatrixPlot[%],MatrixPlot[%%]}
(* Output *)
{[Graphics],[Graphics]}
```

#### ZeroTest

The automatic zero test cannot detect that the following matrix is nonsingular:

```wolfram
m={{1,0},{1,π-Root}}
(* Output *)
{{1,0},{1,π-Root}}
```

```wolfram
Inverse[m]
(* Output *)
Inverse
(* Output *)
Inverse[{{1,0},{1,π-Root}}]
```

The problem is that machine-precision underflows for the bottom right entry:

```wolfram
N[π-Root]
(* Output *)
0.
```

The entry is very small but nonzero:

```wolfram
N[π-Root,10]
(* Output *)
5.836637342329589842368578323586993777×10^-10001
```

Use a zero test employing arbitrary-precision arithmetic to invert the matrix:

```wolfram
Inverse[m,ZeroTest->(N[#,$MachinePrecision]==0&)]
(* Output *)
![image](img/image_003.png)
```

```wolfram
%.m//Simplify
(* Output *)
{{1,0},{0,1}}
```

#### Modulus

Invert a matrix of integers modulo 5:

```wolfram
m={{12,-2},{1,3}};
minv=Inverse[m,Modulus->5]
(* Output *)
{{1,4},{3,4}}
```

Check it:

```wolfram
Mod[m.minv,5]
(* Output *)
{{1,0},{0,1}}
```

### Applications

#### Solving Equations

Solve the system of equations $6 x+9 y=11$, $3 z-7 x=-12$, $5 y+9 z=-9$. First, form the coefficient matrix $a$ and constant vector $b$:

```wolfram
a={{6,9,0},{-7,0,3},{0,5,9}};
b={11,-12,-9};
{a//MatrixForm,b//MatrixForm}
(* Output *)
{({{6, 9, 0}, {-7, 0, 3}, {0, 5, 9}}),({{11}, {-12}, {-9}})}
```

The solution $\{x,y,z \}$ is given by $a.b$:

```wolfram
Inverse[a].b
(* Output *)
{(188)/(159),(23)/(53),-(592)/(477)}
```

Substitute the solution into the original system of equations to verify the solution:

```wolfram
6x+9y==11&&3z-7x==-12&&5y+9z==-9/.Thread[{x,y,z}->%]
(* Output *)
True
```

Solve the matrix equation $m.x=b$:

```wolfram
m={{1,2},{3,4}};b={5,6};
```

Multiplying both sides of the equation on the left by $m$ shows $x=m.b$:

```wolfram
Inverse[m].b
(* Output *)
{-4,(9)/(2)}
```

Confirm the result using [LinearSolve](https://reference.wolfram.com/language/ref/LinearSolve.html):

```wolfram
LinearSolve[m,b]
(* Output *)
{-4,(9)/(2)}
```

For numerical and especially sparse systems, [LinearSolve](https://reference.wolfram.com/language/ref/LinearSolve.html) can be considerably faster:

```wolfram
n=3500;
s=SparseArray[{{i_,i_}->-2.,{i_,j_}/;Abs[i-j]==1->1.},{n,n}];
v=RandomReal[1,n];
```

```wolfram
AbsoluteTiming[LinearSolve[s,v];]
(* Output *)
{0.002592,Null}
```

```wolfram
AbsoluteTiming[Inverse[s].v;]
(* Output *)
{1.042624,Null}
```

Solve the matrix equation $m.x=y$:

```wolfram
m={{-6,7,9},{-2,1,3},{7,-10,4}};y={{-5,6,-4},{-7,7,8},{-2,1,12}};
```

Multiplying both sides of the equation on the left by $m$ shows $x=m.y$:

```wolfram
Inverse[m].y
(* Output *)
{{(158)/(29),-(305)/(58),-(234)/(29)},{4,-(15)/(4),-7},{-(1)/(29),(9)/(116),-(11)/(29)}}
```

Substitute the solution into the equation for verification:

```wolfram
m.x==y/.x->%
(* Output *)
True
```

Solve the system of ODEs $x^{'}=y$, $y^{'}=z$, $z^{'}=-2 x+y+z$. First, construct the coefficient matrix $a$ for the right-hand side:

```wolfram
a={{0,1,0},{0,0,1},{-2,1,2}};
```

Find the eigenvalues and eigenvectors:

```wolfram
{λ,v}=Eigensystem[a]
(* Output *)
{{2,-1,1},{{1,2,4},{1,-1,1},{1,1,1}}}
```

Construct a diagonal matrix whose entries are the exponential of $t \lambda$:

```wolfram
d=DiagonalMatrix[Exp[t λ]]
(* Output *)
{{ℯ^(2 t),0,0},{0,ℯ^-t,0},{0,0,ℯ^t}}
```

Construct the matrix whose columns are the corresponding eigenvectors:

```wolfram
p=Transpose[v]
(* Output *)
{{1,1,1},{2,-1,1},{4,1,1}}
```

The general solution is $p.d.p.\{1,2,3 \}$, for three arbitrary starting values:

```wolfram
p.d.Inverse[p].{C[1],C[2],C[3]}
(* Output *)
{((ℯ^-t)/(3)+ℯ^t-(ℯ^(2 t))/(3)) 1+(-(ℯ^-t)/(2)+(ℯ^t)/(2)) 2+((ℯ^-t)/(6)-(ℯ^t)/(2)+(ℯ^(2 t))/(3)) 3,(-(ℯ^-t)/(3)+ℯ^t-(2 ℯ^(2 t))/(3)) 1+((ℯ^-t)/(2)+(ℯ^t)/(2)) 2+(-(ℯ^-t)/(6)-(ℯ^t)/(2)+(2 ℯ^(2 t))/(3)) 3,((ℯ^-t)/(3)+ℯ^t-(4 ℯ^(2 t))/(3)) 1+(-(ℯ^-t)/(2)+(ℯ^t)/(2)) 2+((ℯ^-t)/(6)-(ℯ^t)/(2)+(4 ℯ^(2 t))/(3)) 3}
```

Verify the solution using [DSolveValue](https://reference.wolfram.com/language/ref/DSolveValue.html):

```wolfram
Simplify[%==DSolveValue[{x'[t]==y[t],y'[t]==z[t],z'[t]==-2x[t]+y[t]+2z[t]},{x[t],y[t],z[t]},t]]
(* Output *)
True
```

#### Change of Basis/Coordinates

Express a general vector in $\mathbb{R}^{3}$ as a linear combination of the vectors $\{1,0,1 \}$, $\{2,2,3 \}$ and $\{-1,-1,1 \}$. First, verify the vectors are linearly independent by checking that their null space is empty:

```wolfram
b_1={1,0,1};b_2={2,2,3};b_3={-1,-1,1};
NullSpace[{b_1,b_2,b_3}]
(* Output *)
{}
```

Form the matrix $p$ whose columns are the basis vectors:

```wolfram
p=Transpose[{b_1,b_2,b_3}];
```

The coefficients of a general vector $v$ are given by $p.v$:

```wolfram
{c_1,c_2,c_3}=Inverse[p].{x,y,z}
(* Output *)
{x-y,-(x)/(5)+(2 y)/(5)+(z)/(5),-(2 x)/(5)-(y)/(5)+(2 z)/(5)}
```

Verify that $v$ does indeed equal the linear combination $\sum_{i=1}^{3}c_{i}b_{i}$:

```wolfram
Simplify[∑_{i=1}^{3}c_ib_i]
(* Output *)
{x,y,z}
```

Find the change-of-basis matrix that transforms coordinates with respect to the basis $\{b_{i}\}$ to coordinates with respect to the basis $\{d_{i}\}$:

```wolfram
{b_1,b_2,b_3}={{2,5,-4},{1,0,3},{-3,3,-2}};
```

```wolfram
{d_1,d_2,d_3}={{-2,4,1},{3,-4,-1},{3,3,-4}};
```

The matrix $p$ whose columns are the $\{b_{i}\}$ transforms from $b$-coordinates to standard coordinates:

```wolfram
(p=Transpose[{b_1,b_2,b_3}])//MatrixForm
(* Output *)
({{2, 1, -3}, {5, 0, 3}, {-4, 3, -2}})
```

The matrix $q$ whose columns are the $\{d_{i}\}$ transforms from $d$-coordinates to standard coordinates:

```wolfram
(q=Transpose[{d_1,d_2,d_3}])//MatrixForm
(* Output *)
({{-2, 3, 3}, {4, -4, 3}, {1, -1, -4}})
```

Its inverse converts from standard coordinates back to $d$-coordinates:

```wolfram
Inverse[q]
(* Output *)
{{1,(9)/(19),(21)/(19)},{1,(5)/(19),(18)/(19)},{0,(1)/(19),-(4)/(19)}}
```

Therefore, $q.p$ converts from $b$-coordinates to $d$-coordinates:

```wolfram
Inverse[q].p
(* Output *)
{{-(1)/(19),(82)/(19),-(72)/(19)},{-(9)/(19),(73)/(19),-(78)/(19)},{(21)/(19),-(12)/(19),(11)/(19)}}
```

Express the linear operator $\mathcal{L}$ whose representation in the standard is given by $m$ in the basis $b_{1}=\{1,0,1 \}$, $b_{2}=\{2,2,3 \}$, $b_{3}=\{-1,-1,2 \}$

```wolfram
(m={{1,2,-3},{2,0,1},{-3,1,5}})//MatrixForm
(* Output *)
({{1, 2, -3}, {2, 0, 1}, {-3, 1, 5}})
```

```wolfram
b_1={1,0,1};b_2={2,2,3};b_3={-1,-1,2};
```

The matrix $p$ whose columns are the $b_{i}$ transforms from $b$-coordinates to standard coordinates:

```wolfram
(p=Transpose[{b_1,b_2,b_3}])//MatrixForm
(* Output *)
({{1, 2, -1}, {0, 2, -1}, {1, 3, 2}})
```

Its inverse converts from standard coordinates to $b$-coordinates:

```wolfram
Inverse[p]
(* Output *)
{{1,-1,0},{-(1)/(7),(3)/(7),(1)/(7)},{-(2)/(7),-(1)/(7),(2)/(7)}}
```

Therefore the representation of $\mathcal{L}$ in $b$-coordinates is:

```wolfram
p.m.Inverse[p]
(* Output *)
{{(67)/(7),-(47)/(7),-(11)/(7)},{8,-7,-1},{-(17)/(7),-(5)/(7),(24)/(7)}}
```

[Inverse](https://reference.wolfram.com/language/ref/Inverse.html) can be used to diagonalize a matrix as $m=p.d.p$. Compute $m$'s eigenvalues and eigenvectors:

```wolfram
m={{9,-7,3},{12,-10,3},{16,-16,1}};
```

```wolfram
Eigensystem[m]
(* Output *)
{{-3,2,1},{{-1,0,4},{1,1,0},{-3,-3,1}}}
```

Construct a diagonal matrix $d$ from the eigenvalues and a matrix $p$ whose columns are the eigenvectors:

```wolfram
{d=DiagonalMatrix[First[%]],p=Transpose[Last[%]]}
(* Output *)
{{{-3,0,0},{0,2,0},{0,0,1}},{{-1,1,-3},{0,1,-3},{4,0,1}}}
```

Confirm the identity $m=p.d.p$:

```wolfram
m==p.d.Inverse[p]
(* Output *)
True
```

Any function of the matrix can now be computed as $f(m)=p.f(d).p$. For example, [MatrixPower](https://reference.wolfram.com/language/ref/MatrixPower.html):

```wolfram
MatrixPower[m,k]==p.MatrixPower[d,k].Inverse[p]
(* Output *)
True
```

Similarly, [MatrixExp](https://reference.wolfram.com/language/ref/MatrixExp.html) becomes trivial, requiring only exponentiating the diagonal elements of $d$:

```wolfram
MatrixExp[m]==p.MatrixExp[d].Inverse[p]
(* Output *)
True
```

```wolfram
MatrixExp[d]
(* Output *)
{{(1)/(ℯ^3),0,0},{0,ℯ^2,0},{0,0,ℯ}}
```

[Inverse](https://reference.wolfram.com/language/ref/Inverse.html) of a transformation matrix gives the matrix for the reverse operation. For example, consider a translation by $\{\Deltax,\Deltat,\Deltaz \}$:

```wolfram
trans=TranslationTransform[{Δx,Δy,Δz}]
(* Output *)
TransformationFunction[({{1, 0, 0, Δx}, {0, 1, 0, Δy}, {0, 0, 1, Δz}, {0, 0, 0, 1}})]
```

The inverse of its transformation matrix gives a translation by the opposite motion:

```wolfram
Inverse[TransformationMatrix[trans]]//MatrixForm
(* Output *)
({{1, 0, 0, -Δx}, {0, 1, 0, -Δy}, {0, 0, 1, -Δz}, {0, 0, 0, 1}})
```

```wolfram
TransformationFunction[%]==TranslationTransform[-{Δx,Δy,Δz}]
(* Output *)
True
```

Consider a general affine transformation:

```wolfram
aff=LinearFractionalTransform[{Array[a_#1,#2&,{2,2}],Array[b_#&,{2}]}]
(* Output *)
TransformationFunction[({{a_1,1, a_1,2, b_1}, {a_2,1, a_2,2, b_2}, {0, 0, 1}})]
```

Construct the inverse transformation:

```wolfram
affInv=TransformationFunction[Inverse[TransformationMatrix[aff]]]
(* Output *)
TransformationFunction[({{(a_2,2)/(-a_1,2 a_2,1+a_1,1 a_2,2), (a_1,2)/(a_1,2 a_2,1-a_1,1 a_2,2), (-b_2 a_1,2+b_1 a_2,2)/(a_1,2 a_2,1-a_1,1 a_2,2)}, {(a_2,1)/(a_1,2 a_2,1-a_1,1 a_2,2), (a_1,1)/(-a_1,2 a_2,1+a_1,1 a_2,2), (-b_2 a_1,1+b_1 a_2,1)/(-a_1,2 a_2,1+a_1,1 a_2,2)}, {0, 0, 1}})]
```

Verify that the two transformations really do undo each other:

```wolfram
affInv[aff[{x,y}]]//Simplify
(* Output *)
{x,y}
```

```wolfram
aff[affInv[{x,y}]]//Simplify
(* Output *)
{x,y}
```

For a mapping $f:x \in \mathbb{R}^{n}\to y \in \mathbb{R}^{n}$, the Jacobian of the inverse mapping $f^{(-1)}$ is given by $(f)$. Consider the mapping from Cartesian to spherical coordinates:

```wolfram
f[x_,y_,z_]:={Sqrt[x^2+y^2+z^2],ArcTan[z,Sqrt[x^2+y^2]],ArcTan[x,y]}
```

Compute the Jacobian at the point $\{1,1,0 \}$:

```wolfram
fJac=Grad[f[x,y,z],{x,y,z}]/.{x->1,y->1,z->0}
(* Output *)
{{(1)/(Sqrt[2]),(1)/(Sqrt[2]),0},{0,0,-(1)/(Sqrt[2])},{-(1)/(2),(1)/(2),0}}
```

The inverse mapping is the transformation from spherical back to Cartesian coordinates:

```wolfram
fInv[r_,θ_,φ_]:={r Cos[φ]Sin[θ],r Sin[θ]Sin[φ],r Cos[θ]}
```

Verify the inverse function:

```wolfram
fInv@@f[x,y,z]
(* Output *)
{x,y,z}
```

Compute the Jacobian of the inverse mapping at the $\{r,\theta,\varphi \}$ coordinates corresponding to $\{1,1,0 \}$:

```wolfram
fInvJac=Grad[fInv[r,θ,φ],{r,θ,φ}]/.Thread[{r,θ,φ}->f[1,1,0]]
(* Output *)
{{(1)/(Sqrt[2]),0,-1},{(1)/(Sqrt[2]),0,1},{0,-Sqrt[2],0}}
```

Confirm the identity $f^{(-1)}=(f)$:

```wolfram
fInvJac==Inverse[fJac]
(* Output *)
True
```

### Properties & Relations

[Inverse](https://reference.wolfram.com/language/ref/Inverse.html) satisfies the relation $a.a=a.a=IdentityMatrix[n]$ for an $n \times n$ matrix $a$:

```wolfram
a={{1,2},{3,4}}
(* Output *)
{{1,2},{3,4}}
```

```wolfram
a.Inverse[a]==Inverse[a].a==IdentityMatrix[2]
(* Output *)
True
```

[Inverse](https://reference.wolfram.com/language/ref/Inverse.html) satisfies the relation $(a.b)=b.a$:

```wolfram
a={{1,1,1},{6,9,7},{8,1,9}};
```

```wolfram
b={{0,3,9},{7,9,7},{4,4,1}};
```

```wolfram
Inverse[a.b]==Inverse[b].Inverse[a]
(* Output *)
True
```

[Inverse](https://reference.wolfram.com/language/ref/Inverse.html) satisfies the relation $(a)=(a)$:

```wolfram
a={{-(1)/(25),(12)/(25),-(7)/(25)},{-(7)/(25),(9)/(25),(1)/(25)},{(13)/(75),-(31)/(75),(16)/(75)}};
```

```wolfram
Inverse[Transpose[a]]==Transpose[Inverse[a]]
(* Output *)
True
```

A square matrix has an inverse if and only if its determinant is nonzero:

```wolfram
m=({{1, 2, 3}, {4, 5, 6}, {7, 8, 9}});
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
Inverse[{{1,2,3},{4,5,6},{7,8,9}}]
```

Moreover, determinant of the inverse $a$ equals $\frac{1}{a}$:

```wolfram
a=RandomReal[1,{5,5}];
```

```wolfram
Det[Inverse[a]]==(1)/(Det[a])
(* Output *)
True
```

[MatrixPower](https://reference.wolfram.com/language/ref/MatrixPower.html)[*m*,-1] equals [Inverse](https://reference.wolfram.com/language/ref/Inverse.html)[*m*]:

```wolfram
m={{9,4,1,8},{4,8,3,4},{7,6,8,4},{6,9,0,9}};
```

```wolfram
MatrixPower[m,-1]==Inverse[m]
(* Output *)
True
```

For an invertible matrix $a$, [Inverse](https://reference.wolfram.com/language/ref/Inverse.html)[*a*] equals [Adjugate](https://reference.wolfram.com/language/ref/Adjugate.html)[*a*]/[Det](https://reference.wolfram.com/language/ref/Det.html)[*a*]:

```wolfram
a={{1,3,5},{7,1,-1},{8,4,1}};
```

```wolfram
Inverse[a]==Adjugate[a]/Det[a]
(* Output *)
True
```

[Inverse](https://reference.wolfram.com/language/ref/Inverse.html)[*m*] equals [LinearSolve](https://reference.wolfram.com/language/ref/LinearSolve.html)[*m*,[IdentityMatrix](https://reference.wolfram.com/language/ref/IdentityMatrix.html)[*n*]] for an invertible `*n*×*n*` matrix:

```wolfram
m=RandomReal[1,{5,5}];
```

```wolfram
Inverse[m]==LinearSolve[m,IdentityMatrix[5]]
(* Output *)
True
```

The inverse of an orthogonal matrix is given by [Transpose](https://reference.wolfram.com/language/ref/Transpose.html):

```wolfram
r ={{-(6)/(7),(2)/(7),(3)/(7)},{(2)/(7),-(3)/(7),(6)/(7)},{(3)/(7),(6)/(7),(2)/(7)}};
```

```wolfram
OrthogonalMatrixQ[r]
(* Output *)
True
```

```wolfram
Inverse[r]==Transpose[r]
(* Output *)
True
```

The inverse of a unitary matrix is given by [ConjugateTranspose](https://reference.wolfram.com/language/ref/ConjugateTranspose.html):

```wolfram
u={{-(8)/(9),-(2ⅈ)/(9),(2)/(9)-(ⅈ)/(3)},{(2ⅈ)/(9),-(5)/(9),(2)/(3)+(4ⅈ)/(9)},{(2)/(9)+(ⅈ)/(3),(2)/(3)-(4ⅈ)/(9),(4)/(9)}};
```

```wolfram
UnitaryMatrixQ[u]
(* Output *)
True
```

```wolfram
Inverse[u]==ConjugateTranspose[u]
(* Output *)
True
```

A matrix and its inverse have the same symmetry:

```wolfram
s={{3,0,-5},{0,1,0},{-5,0,0}};
```

```wolfram
{TensorSymmetry[s],TensorSymmetry[Inverse[s]]}
(* Output *)
{Symmetric[{1,2}],Symmetric[{1,2}]}
```

```wolfram
a={{0,-(3)/(2),(5)/(2),-6},{(3)/(2),0,4,3},{-(5)/(2),-4,0,-5},{6,-3,5,0}};
```

```wolfram
{TensorSymmetry[a],TensorSymmetry[Inverse[a]]}
(* Output *)
{Antisymmetric[{1,2}],Antisymmetric[{1,2}]}
```

A [QuantityArray](https://reference.wolfram.com/language/ref/QuantityArray.html) and its inverse have reciprocal units:

```wolfram
QuantityArray[{{1,2},{3,4}},{("Meters")/("Seconds"),("Meters")/("Seconds")}]
(* Output *)
QuantityArray[...]
```

```wolfram
Inverse[%]
(* Output *)
QuantityArray[...]
```

For an invertible matrix $a$, [Inverse](https://reference.wolfram.com/language/ref/Inverse.html)[*a*] and [PseudoInverse](https://reference.wolfram.com/language/ref/PseudoInverse.html)[*a*] coincide:

```wolfram
a=({{4, 3, 8}, {9, 5, 2}, {9, 6, 7}});
```

```wolfram
Inverse[a]===PseudoInverse[a]
(* Output *)
True
```

[PseudoInverse](https://reference.wolfram.com/language/ref/PseudoInverse.html) extends to singular as well as rectangular matrices:

```wolfram
PseudoInverse[({{1, 2}, {1, 2}})]//MatrixForm
(* Output *)
({{(1)/(10), (1)/(10)}, {(1)/(5), (1)/(5)}})
```

```wolfram
PseudoInverse[({{1, 2}, {3, 4}, {5, 6}})]//MatrixForm
(* Output *)
({{-(4)/(3), -(1)/(3), (2)/(3)}, {(13)/(12), (1)/(3), -(5)/(12)}})
```

For an invertible matrix $a$, [Inverse](https://reference.wolfram.com/language/ref/Inverse.html)[*a*] and [DrazinInverse](https://reference.wolfram.com/language/ref/DrazinInverse.html)[*a*] coincide:

```wolfram
a=({{4, 3, 8}, {9, 5, 2}, {9, 6, 7}});
```

```wolfram
Inverse[a]===DrazinInverse[a]
(* Output *)
True
```

[DrazinInverse](https://reference.wolfram.com/language/ref/DrazinInverse.html) extends to singular square matrices:

```wolfram
DrazinInverse[({{1, 2}, {1, 2}})]//MatrixForm
(* Output *)
({{(1)/(9), (2)/(9)}, {(1)/(9), (2)/(9)}})
```

### Possible Issues

The inverse may not exist:

```wolfram
Inverse[({{1, 2}, {1, 2}})]
(* Output *)
Inverse
(* Output *)
Inverse[{{1,2},{1,2}}]
```

Typically a pseudo inverse does:

```wolfram
PseudoInverse[({{1, 2}, {1, 2}})]
(* Output *)
{{(1)/(10),(1)/(10)},{(1)/(5),(1)/(5)}}
```

Full inverses do not exist for rectangular matrices:

```wolfram
Inverse[({{1, 2, 2}, {3, 1, 4}})]
(* Output *)
Inverse
(* Output *)
Inverse[{{1,2,2},{3,1,4}}]
```

Use [PseudoInverse](https://reference.wolfram.com/language/ref/PseudoInverse.html) instead:

```wolfram
PseudoInverse[({{1, 2, 2}, {3, 1, 4}})]
(* Output *)
{{-(1)/(5),(14)/(65)},{(3)/(5),-(17)/(65)},{0,(2)/(13)}}
```

Accurate inverses cannot be found for ill-conditioned machine-precision numerical matrices:

```wolfram
Inverse[N@HilbertMatrix[15]][[1,1]]
(* Output *)
Inverse
(* Output *)
160.68145087361336
```

Exact result:

```wolfram
Inverse[HilbertMatrix[15]][[1,1]]
(* Output *)
225
```

Arbitrary-precision result:

```wolfram
Inverse[N[HilbertMatrix[15],30]][[1,1]]
(* Output *)
225.0000000000000000000000000000000000003971413445837896359332
```

The computation may fail if a non-prime modulus is provided:

```wolfram
Inverse[{{1,2},{3,2}},Modulus->6]
(* Output *)
Inverse
(* Output *)
Inverse[{{1,2},{3,2}},Modulus->6]
```

## Tech Notes ▪Vectors and Matrices ▪Matrix Inversion ▪Implementation notes: Numerical and Related Functions

## Related Guides ▪Matrix Operations ▪Linear Systems ▪Matrices and Linear Algebra ▪Finite Mathematics ▪Finite Fields ▪GPU Computing ▪Symbolic Vectors, Matrices and Arrays ▪Structured Arrays ▪GPU Computing with NVIDIA

## History Introduced in 1988 (1.0) | Updated in 1996 (3.0) ▪ 2022 (13.2) ▪ 2024 (14.0) ▪ 2025 (14.3)
