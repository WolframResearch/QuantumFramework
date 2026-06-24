# Dot (.) | [SpanFromLeft]

> `*a*.*b*.*c*` or [Dot](https://reference.wolfram.com/language/ref/Dot.html)[*a*,*b*,*c*] — gives products of vectors, matrices, and tensors.

## Details

`*a.b*` gives an explicit result when `*a*` and `*b*` are lists with appropriate dimensions. It contracts the last index in `*a*` with the first index in `*b*`.

Various applications of [Dot](https://reference.wolfram.com/language/ref/Dot.html):

{*a*_1,*a*_2}.{*b*_1,*b*_2} | scalar product of vectors
{*a*_1,*a*_2}.{{*m*_11,*m*_12},{*m*_21,*m*_22}} | [SpanFromLeft]
product of a vector and a matrix
{{*m*_11,*m*_12},{*m*_21,*m*_22}}.{*a*_1,*a*_2} | [SpanFromLeft]
product of a matrix and a vector
{{*m*_11,*m*_12},{*m*_21,*m*_22}}.{{*n*_11,*n*_12},{*n*_21,*n*_22}} | [SpanFromLeft]
product of two matrices

The result of applying [Dot](https://reference.wolfram.com/language/ref/Dot.html) to two tensors $T_{i_{1} i_{2}... i_{n}}$ and $U_{j_{1} j_{2}... j_{m}}$ is the tensor $\sum_{k} T_{i_{1}i_{2}...i_{n-1}k} U_{k j_{2}...j_{m}}$. Applying [Dot](https://reference.wolfram.com/language/ref/Dot.html) to a rank $n$ tensor and a rank $m$ tensor gives a rank $m+n-2$ tensor.

[Dot](https://reference.wolfram.com/language/ref/Dot.html) can be used on [SparseArray](https://reference.wolfram.com/language/ref/SparseArray.html) and structured array objects. It will return an object of the same type as the input when possible.

[Dot](https://reference.wolfram.com/language/ref/Dot.html) is linear in all arguments. » It does not define a complex (Hermitian) inner product on vectors.

When its arguments are not lists or sparse arrays, [Dot](https://reference.wolfram.com/language/ref/Dot.html) remains unevaluated. It has the attribute [Flat](https://reference.wolfram.com/language/ref/Flat.html).

## Examples

### Basic Examples

Scalar product of vectors in three dimensions:

```wolfram
{a,b,c} . {x,y,z}
(* Output *)
a x+b y+c z
```

Scalar product of vectors in two dimensions:

```wolfram
u={1,1};
v={-1,1};
u.v
(* Output *)
0
```

Vectors are perpendicular if their inner product is zero:

```wolfram
θ=VectorAngle[u,v]
(* Output *)
(π)/(2)
```

Visualize the vectors:

```wolfram
GeometricScene[{u->u,v->v}, {PlanarAngle[{u,{0,0},v}]==θ}]["Graphics"]
```

*([Graphics])*

The product $M v$ of a matrix and a vector:

```wolfram
{{a,b},{c,d}} . {x, y}
(* Output *)
{a x+b y,c x+d y}
```

The product $v M$ of a vector and a matrix:

```wolfram
{x,y} . {{a,b},{c,d}}
(* Output *)
{a x+c y,b x+d y}
```

The product $v M w$ of a matrix and two vectors:

```wolfram
{x,y} . {{a,b},{c,d}} . {r,s}
(* Output *)
r (a x+c y)+s (b x+d y)
```

The product of two matrices:

```wolfram
{{a,b},{c,d}} . {{1,2},{3,4}} //MatrixForm
(* Output *)
({{a+3 b, 2 a+4 b}, {c+3 d, 2 c+4 d}})
```

Multiply in the other order:

```wolfram
{{1,2},{3,4}} .{{a,b},{c,d}} //MatrixForm
(* Output *)
({{a+2 c, b+2 d}, {3 a+4 c, 3 b+4 d}})
```

Use rectangular matrices:

```wolfram
{{1,2,3},{4,5,6}} .{{a,b},{c,d},{e,f}} //MatrixForm
(* Output *)
({{a+2 c+3 e, b+2 d+3 f}, {4 a+5 c+6 e, 4 b+5 d+6 f}})
```

### Scope

#### Dot Products of Vectors

Scalar product of machine-precision vectors:

```wolfram
Dot[{3.2,4.2,5.2},{0.75,1.1,0.0625}]
(* Output *)
7.3450000000000015
```

Dot product of exact vectors:

```wolfram
{1,2,3,4,5}.{1,8,9,0,-1}
(* Output *)
39
```

Inner product of symbolic vectors:

```wolfram
{a,b,c} . {0,b,d+e}
(* Output *)
b^2+c (d+e)
```

The dot product of arbitrary-precision vectors:

```wolfram
{u,v}=RandomReal[4,{2, 3},WorkingPrecision->20]
(* Output *)
{{3.66166430850333286841345216411358620689,2.71775849268176529241968462924372573752,2.96377971465521267290682029538206876396},{1.30193422338032984423216928404620773563,0.86344100369857373434082092877872582903,1.81290881141712464057333051292708603341}}
```

```wolfram
u.v
(* Output *)
12.48693255829999190020413690101423411852
```

[Dot](https://reference.wolfram.com/language/ref/Dot.html) allows complex inputs, but does not conjugate any of them:

```wolfram
u={1,2-I};
v={z,3};
u.v
(* Output *)
(6-3 ⅈ)+z
```

To compute the complex or Hermitian inner product, apply [Conjugate](https://reference.wolfram.com/language/ref/Conjugate.html) to one of the inputs:

```wolfram
cDot[a_,b_]:=Conjugate[a].b
```

```wolfram
cDot[u,v]
(* Output *)
(6+3 ⅈ)+z
```

Some sources, particularly in the mathematical literature, conjugate the second argument:

```wolfram
cDot2[a_,b_]:=a.Conjugate[b]
```

```wolfram
cDot2[u,v]
(* Output *)
(6-3 ⅈ)+Conjugate[z]
```

Compute the norm of `u` using the two inner products:

```wolfram
{Sqrt[cDot[u,u]],Sqrt[cDot2[u,u]]}
(* Output *)
{Sqrt[6],Sqrt[6]}
```

Verify the result using [Norm](https://reference.wolfram.com/language/ref/Norm.html):

```wolfram
Norm[u]
(* Output *)
Sqrt[6]
```

Dot product of sparse vectors:

```wolfram
v=SparseArray[{1->1, 50->3},{100}]
(* Output *)
SparseArray[...]
```

```wolfram
w=SparseArray[{50->7,100->5},{100}]
(* Output *)
SparseArray[...]
```

```wolfram
v.w
(* Output *)
21
```

Compute the scalar product of two [QuantityArray](https://reference.wolfram.com/language/ref/QuantityArray.html) vectors:

```wolfram
x=QuantityArray[{1,2,3},"Meters"]
(* Output *)
QuantityArray[...]
```

```wolfram
f=QuantityArray[{0,0,-9.8},"Newtons"]
(* Output *)
QuantityArray[...]
```

```wolfram
x.f
(* Output *)
-29.400000000000002
```

#### Matrix-Vector Multiplication

Define a rectangular matrix of dimensions $2 \times 3$:

```wolfram
r={{0.187902,0.498054,0.767621},{0.226789,0.852257,0.819982}};
r//MatrixForm
(* Output *)
({{0.187902, 0.498054, 0.767621}, {0.226789, 0.852257, 0.819982}})
```

Define a 2-vector and a 3-vector:

```wolfram
v2={0.618678,0.213605};
v3={0.804978,0.587651,0.2951};
{v2//MatrixForm,v3//MatrixForm}
(* Output *)
{({{0.618678}, {0.213605}}),({{0.804978}, {0.587651}, {0.2951}})}
```

The $2 \times 3$ matrix can be multiplied by the 2-vector only on the left:

```wolfram
v2.r
(* Output *)
{0.16469409790099998,0.490181409097,0.6500624801479999}
```

Multiplying in the opposite order produces an error message due to the incompatible shapes:

```wolfram
r.v2
(* Output *)
Dot
(* Output *)
{{0.187902,0.498054,0.767621},{0.226789,0.852257,0.819982}}.{0.618678,0.213605}
```

Similarly, the $2 \times 3$ matrix can be multiplied by the 3-vector only on the right:

```wolfram
r.v3
(* Output *)
{0.67046386441,0.9253665221490001}
```

Multiply the matrix on both sides at once:

```wolfram
v2.r.v3
(* Output *)
0.6124641586690871
```

Define a square matrix and a compatible vector:

```wolfram
m={{1,2},{3,4}};
v={5,6};
```

The products `m.v` and `v.m` return different vectors:

```wolfram
m . v
(* Output *)
{17,39}
```

```wolfram
v.m
(* Output *)
{23,34}
```

The product `v.m.v` is a scalar:

```wolfram
v.m.v
(* Output *)
319
```

Define a column and row matrices `c` and `r` with the same numerical entries as `v`:

```wolfram
c={{5},{6}};
r = {{5,6}};
```

Products involving `m`, `c` and `r` have the same entries as those involving `m` and `v`, but are all matrices:

```wolfram
r.m
(* Output *)
{{23,34}}
```

```wolfram
m.c
(* Output *)
{{17},{39}}
```

```wolfram
r.m.c
(* Output *)
{{319}}
```

Moreover, the products must be done in an order that respects the matrices' shapes:

```wolfram
c.m.r
(* Output *)
Dot
(* Output *)
Dot
(* Output *)
{{5},{6}}.{{1,2},{3,4}}.{{5,6}}
```

Define a matrix and two vectors:

```wolfram
m={{.5,.32},{.19, .73}};
u={1.5,.27};
v={-3.2,5.5};
```

Since $m.u$ is a vector, $m.u.v$ is an allowed product:

```wolfram
m.u.v
(* Output *)
-0.024929999999999897
```

Note that it is effectively multiplying $v$ on the left side of the matrix, not the right:

```wolfram
{v.m.u,u.m.v}
(* Output *)
{-0.02493000000000023,1.1598899999999999}
```

The product of a sparse matrix and sparse vector is a sparse vector:

```wolfram
SparseArray[{{1,1}->5, {10,10}->10}].SparseArray[{{10}->7}]
(* Output *)
SparseArray[...]
```

Format the result as a row matrix:

```wolfram
MatrixForm[{%}]
(* Output *)
({0, 0, 0, 0, 0, 0, 0, 0, 0, 70})
```

The product of a sparse matrix and an ordinary vector is a normal vector:

```wolfram
SparseArray[{{1,1}->5, {10,10}->10}].{1,2,3,4,5,6,7,8,9,10}
(* Output *)
{5,0,0,0,0,0,0,0,0,100}
```

The product of a structured matrix with a vector will retain the structure if possible:

```wolfram
mat=SymmetrizedArray[{{1,1}->2,{1,2}->3,{2,2}->5},{2,2},Symmetric[{1,2}]]
(* Output *)
SymmetrizedArray[...]
```

```wolfram
mat.{7,11}
(* Output *)
SymmetrizedArray[...]
```

```wolfram
mat.SparseArray[{1,2}]
(* Output *)
SymmetrizedArray[...]
```

The product of a normal matrix with a structured vector may have the structure of the vector:

```wolfram
{{2,3},{3,5}}.QuantityArray[{a,b},"Meters"]
(* Output *)
QuantityArray[...]
```

#### Matrix-Matrix Multiplication

Multiply real machine-number matrices:

```wolfram
m={{1.2,3.2,5.2},{2.2,4.2,-6.4},{3.1,5.1,7.3}};
n={{4.2,6.3,8.2},{2.5,-7.3,9.3},{6.3,8.3,-1.10}};
m.n // MatrixForm
(* Output *)
({{45.8, 27.360000000000003, 33.88}, {-20.58, -69.92, 64.14}, {71.75999999999999, 42.89000000000001, 64.82}})
```

Product of complex matrices:

```wolfram
m={{1+I,2,3-2 I},{0,4,5I},{0,0,6}};
n={{6+I,4,5-7 I},{5,3,2I},{5,2,7}};
m.n // MatrixForm
(* Output *)
({{30-3 ⅈ, 16, 33-12 ⅈ}, {20+25 ⅈ, 12+10 ⅈ, 43 ⅈ}, {30, 12, 42}})
```

Product of exact matrices:

```wolfram
m={{1,2},{3,4},{5,6}};
n={{6,5,4},{3,2,1}};
m.n // MatrixForm
(* Output *)
({{12, 9, 6}, {30, 23, 16}, {48, 37, 26}})
```

Multiply in the other order:

```wolfram
n.m // MatrixForm
(* Output *)
({{41, 56}, {14, 20}})
```

Visualize the input and output matrices:

```wolfram
{MatrixPlot[m],MatrixPlot[n],MatrixPlot[m.n],MatrixPlot[n.m]}
(* Output *)
{[Graphics],[Graphics],[Graphics],[Graphics]}
```

Multiply arbitrary-precision matrices:

```wolfram
m=RandomReal[1,{4, 3},WorkingPrecision->20];
n=RandomReal[1,{3, 2},WorkingPrecision->20];
m.n//MatrixForm
(* Output *)
({{0.84173182903817583766227029380659489026, 0.20736824783724487078347807927183769689}, {0.84558988485292398970689172291943892921, 0.21618558618735590812433736392960315309}, {0.48073040247608544567468780077817013134, 0.14786405264561511401467340094638619303}, {0.93340095260269193753731971236541820009, 0.24432890618735353846283378198669105176}})
```

Since $2 \neq 4$, these matrices cannot be multiplied in the opposite order:

```wolfram
n.m
(* Output *)
Dot
(* Output *)
{{0.48760262601522473449832924610736206716,0.16722367150868444635208225415579086359},{0.57111873569377479163884127333350448907,0.1527562867416397200434727450535365989},{0.54297160522291044448275031120143552243,0.06497936096616323168684379446569110428}}.{{0.29807446359421927058864056417353260997,0.87886522522462858804772126153359579348,0.35812838037121656357069472351173367031},{0.6308547488450924121584245142513314164,0.54863374293381510870919442329762816257,0.41373884176936144358631506495749263763},{0.76958362489562906655381820802963588335,0.07757936120504195174881010427236915916,0.11266220694209821092160971056639340304},{0.94862849167543497226753402076782073493,0.34768780469687570364406385692301881818,0.50145566894328597130821864938376997145}}
```

Product of symbolic matrices:

```wolfram
{{a,b},{c,d}} . {{r,s},{t,u}}
(* Output *)
{{a r+b t,a s+b u},{c r+d t,c s+d u}}
```

Product of finite field element matrices:

```wolfram
ℱ=FiniteField[29,4];
m={{ℱ[12],ℱ[23],ℱ[34]},{ℱ[45],ℱ[56],ℱ[67]},{ℱ[78],ℱ[89],ℱ[90]}};
n={{ℱ[123],ℱ[234]},{ℱ[345],ℱ[456]},{ℱ[567],ℱ[678]}};
m.n // MatrixForm
(* Output *)
![image](img/image_001.png)
```

Product of [CenteredInterval](https://reference.wolfram.com/language/ref/CenteredInterval.html) matrices:

```wolfram
m=Map[CenteredInterval,RandomReal[{-10,10},{3,3},WorkingPrecision->10],{2}];
n=Map[CenteredInterval,RandomReal[{-10,10},{3,2},WorkingPrecision->10],{2}];
(mn=m.n)//MatrixForm
(* Output *)
({{<|Interpretation -> interpretation, Center -> -6.6624473402880539652, Radius -> 1.33294422548224655855619857903×10^-9, Type -> Real|>, <|Interpretation -> interpretation, Center -> -7.2770320155568697373, Radius -> 1.45586115905815827176184029668×10^-9, Type -> Real|>}, {<|Interpretation -> interpretation, Center -> 21.5888116891455865698, Radius -> 4.31958136620780308589928608853×10^-9, Type -> Real|>, <|Interpretation -> interpretation, Center -> 39.1475113060878356919, Radius -> 7.83314028918180582650165888481×10^-9, Type -> Real|>}, {<|Interpretation -> interpretation, Center -> -69.2062364945959416218, Radius -> 1.384852334196384759934517205693×10^-8, Type -> Real|>, <|Interpretation -> interpretation, Center -> -115.0480946206080261618, Radius -> 2.301689497330450251411093631759×10^-8, Type -> Real|>}})
```

Find random representatives `*mrep*` and `*nrep*` of `*m*` and `*n*`:

```wolfram
ranrep[e_CenteredInterval]:=e["Center"]+RandomInteger[{-1000,1000}]/1000 e["Radius"]
(mrep=Map[ranrep,m,{2}])//MatrixForm
(* Output *)
({{(6551174605072648920581)/(36893488147419103232000), (1111534564596958329701)/(2305843009213693952000), -(8029344229444342909783)/(9223372036854775808000)}, {-(1346218409338961650963)/(1152921504606846976000), -(601542105980571238219)/(115292150460684697600), (116320476289456897721)/(46116860184273879040)}, {(184348245749988189461)/(23058430092136939520), (1413459733922495261001)/(144115188075855872000), -(9007626060915918993739)/(1152921504606846976000)}})
```

```wolfram
(nrep=Map[ranrep,n,{2}])//MatrixForm
(* Output *)
({{-(517735765950126532629)/(461168601842738790400), -(19636563731511719131)/(4503599627370496000)}, {-(374610193631715320723)/(922337203685477580800), -(367287916988837610003)/(92233720368547758080)}, {(83002382928117078027)/(11529215046068469760), (3034923626689095477493)/(576460752303423488000)}})
```

Verify that `*mn*` contains the product of `*mrep*` and `*nrep*`:

```wolfram
MapThread[IntervalMemberQ,{mn,mrep.nrep},2]//MatrixForm
(* Output *)
({{True, True}, {True, True}, {True, True}})
```

The product of sparse matrices is another sparse matrix:

```wolfram
Dot[SparseArray[...],SparseArray[...]]
(* Output *)
SparseArray[...]
```

Format the result:

```wolfram
%//MatrixForm
(* Output *)
({{27, 0, 0}, {0, 0, 22}})
```

The product of structured matrices preserves the structure if possible:

```wolfram
Dot[SymmetrizedArray[...],SymmetrizedArray[...]]
(* Output *)
SymmetrizedArray[...]
```

Format the result:

```wolfram
%//MatrixForm
(* Output *)
({{16, 0, 0}, {0, 175, 0}, {0, 0, 175}})
```

Cube a matrix:

```wolfram
m={{a,1},{0,b}};
```

```wolfram
m.m.m//MatrixForm
(* Output *)
({{a^3, a^2+b (a+b)}, {0, b^3}})
```

Compare with [MatrixPower](https://reference.wolfram.com/language/ref/MatrixPower.html):

```wolfram
MatrixPower[m,3]//MatrixForm
(* Output *)
({{a^3, b^2+a (a+b)}, {0, b^3}})
```

Raise a matrix to the tenth power using [Dot](https://reference.wolfram.com/language/ref/Dot.html) in combination with [Apply](https://reference.wolfram.com/language/ref/Apply.html) (`@@`) and [ConstantArray](https://reference.wolfram.com/language/ref/ConstantArray.html):

```wolfram
Dot @@ ConstantArray[m,{10}]//Expand
(* Output *)
{{a^10,a^9+a^8 b+a^7 b^2+a^6 b^3+a^5 b^4+a^4 b^5+a^3 b^6+a^2 b^7+a b^8+b^9},{0,b^10}}
```

Verify the result:

```wolfram
MatrixPower[m,10]//Expand
(* Output *)
{{a^10,a^9+a^8 b+a^7 b^2+a^6 b^3+a^5 b^4+a^4 b^5+a^3 b^6+a^2 b^7+a b^8+b^9},{0,b^10}}
```

Efficiently multiply large matrices:

```wolfram
mat=RandomReal[{0,9},{1000,1000}];
mat2=RandomComplex[1+ⅈ,{1000,300}];
Dot[mat,mat2];//AbsoluteTiming
(* Output *)
{0.033001,Null}
```

#### Higher-Rank Arrays

[Dot](https://reference.wolfram.com/language/ref/Dot.html) works for arrays of any rank:

```wolfram
a = RandomInteger[9, {2,3,4}];
b =RandomInteger[9, {4,5,2}];
c = RandomInteger[9, {2}];
```

```wolfram
a.b
(* Output *)
{{{{14,33},{53,52},{35,58},{58,58},{43,15}},{{14,38},{113,93},{55,62},{74,84},{80,34}},{{37,78},{160,143},{92,142},{127,130},{139,67}}},{{{26,49},{158,137},{74,100},{94,102},{115,61}},{{21,29},{120,106},{51,74},{58,62},{86,54}},{{39,69},{174,150},{90,138},{105,108},{156,90}}}}
```

The dimensions of the result are those of the input with the common dimension collapsed:

```wolfram
Dimensions[%]
(* Output *)
{2,3,5,2}
```

Any combination is allowed as long as products are done with a common dimension:

```wolfram
c.a.b.c
(* Output *)
{{2314,7717,5115,6122,4651},{1994,8781,4819,5682,5264},{4369,12740,9148,9603,9375}}
```

Create a rank-three array with three equal dimensions:

```wolfram
(a= {{{8,0},{4,2}}, {{0,7},{5,6}}})//MatrixForm
(* Output *)
({{({{8}, {0}}), ({{4}, {2}})}, {({{0}, {7}}), ({{5}, {6}})}})
```

Create three vectors of the same dimension:

```wolfram
x={1,2};
y={(1)/(2),(1)/(3)};
z={-1,1};
```

$a.x.y.z$ is the complete contraction $\sum_{ i j k} a_{\text{i j k}}x_{k}y_{j}z_{i}$ that pairs $x$ with $a$'s last level and $z$ with its first:

```wolfram
a.x.y.z
(* Output *)
6
```

$a.z.y.x$ is the different contraction $\sum_{ i j k} a_{\text{i j k}}x_{i}y_{j}z_{k}$ that pairs $x$ with $a$'s first level and $z$ with its last:

```wolfram
a.z.y.x
(* Output *)
3
```

Contract both levels of `m` with the second and third levels of `a`, respectively:

```wolfram
a=RandomReal[1,{2,3,4}];
m=RandomReal[1,{3,4}];
Flatten[m].Flatten[a,{2,3}]
(* Output *)
{2.278186162557501,3.3158749203975413}
```

[Dot](https://reference.wolfram.com/language/ref/Dot.html) of two sparse arrays is generally another sparse array:

```wolfram
a=SparseArray[...];
```

```wolfram
m = SparseArray[...];
```

```wolfram
a.m
(* Output *)
SparseArray[...]
```

[Dot](https://reference.wolfram.com/language/ref/Dot.html) of a sparse array and an ordinary list may be another sparse array or an ordinary list:

```wolfram
a.{1,2,3,4}
(* Output *)
SparseArray[...]
```

```wolfram
a.{{1},{2},{3},{4}}
(* Output *)
{{{{y},{0},{0}},{{0},{0},{4 x}}},{{{0},{0},{0}},{{0},{0},{0}}}}
```

Format the rank-three array:

```wolfram
%//MatrixForm
(* Output *)
({{({{y}, {0}, {0}}), ({{0}, {0}, {4 x}})}, {({{0}, {0}, {0}}), ({{0}, {0}, {0}})}})
```

The product of two [SymmetrizedArray](https://reference.wolfram.com/language/ref/SymmetrizedArray.html) objects is generally another symmetrized array:

```wolfram
s=SymmetrizedArray[...];
m=SymmetrizedArray[...];
s.m
(* Output *)
SymmetrizedArray[...]
```

The symmetry of the new array may be much more complicated than the symmetry of either input:

```wolfram
a=SymmetrizedArray[...];
s.a
(* Output *)
SymmetrizedArray[...]
```

```wolfram
TensorSymmetry[%]
(* Output *)
{{Cycles[{{1,2}}],1},{Cycles[{{2,3}}],1},{Cycles[{{4,5}}],-1}}
```

### Applications

#### Projections and Bases

Project the vector $v=\{-1,3 \}$ on the line spanned by the vector $l=\{1,1 \}$:

```wolfram
v={-1,3};l={1,1};
p=(v.l)/(l.l)l
(* Output *)
{1,1}
```

Visualize $v$ and its projection onto the line spanned by $l$:

```wolfram
Graphics[{InfiniteLine[{0,0},l],Arrow[{{0,0},v}], Directive[Red,Thick],Arrow[{{0,0},p}],Directive[Black,Dotted],Arrow[{p,v}]},PlotRangePadding->.5,Axes->True]
```

*([Graphics])*

Project the vector $v=\{-1,-4,-2 \}$ on the plane spanned by the vectors $\{1,1,-1 \}$ and $\{1,0,1 \}$:

```wolfram
v={1,2,1/2};
b1={2,4,-2};
b2={-3,3,0};
```

First, replace $b_{2}$ with a vector in the plane perpendicular to $b_{1}$:

```wolfram
b3=b2-(b2.b1)/(b1.b1)b1
(* Output *)
{-(7)/(2),2,(1)/(2)}
```

The projection in the plane is the sum of the projections onto $b_{1}$ and $b_{3}$:

```wolfram
p=(v.b1)/(b1.b1)b1+(v.b3)/(b3.b3)b3
(* Output *)
{(13)/(22),(35)/(22),-(8)/(11)}
```

Find the component perpendicular to the plane:

```wolfram
v-p
(* Output *)
{(9)/(22),(9)/(22),(27)/(22)}
```

Confirm the result by projecting $v$ onto the normal to the plane:

```wolfram
% == (v.(b1⨯b2))/((b1⨯b2).(b1⨯b2))(b1⨯b2)
(* Output *)
True
```

Visualize the plane, the vector and its parallel and perpendicular components:

```wolfram
Graphics3D[{FaceForm[LightGray],InfinitePlane[{0,0,0},{b1,b2}],Arrow[{{0,0,0},v}],Directive[Red,Thick],Arrow[Tube[{{0,0,0},p}]],Directive[Black,Dotted],Arrow[{p,v}]},PlotRangePadding->1,Axes->True]
```

*([Graphics3D])*

Apply the Gram-Schmidt process to construct an orthonormal basis from the following vectors:

```wolfram
{v_1,v_2,v_3,v_4} ={{-0.449,-0.028,-0.209,0.376},{0.547,-0.943,0.141,-0.522},{0.405,-0.078,-0.511,0.532},{-0.358,-0.452,0.651,-0.13}};
```

The first vector in the orthonormal basis, $e_{1}$, is merely the normalized multiple $v_{1}$:

```wolfram
e_1=Normalize[v_1]
(* Output *)
{-0.7213449438994951,-0.04498364906277475,-0.3357708090757115,0.6040661445572609}
```

For subsequent vectors, components parallel to earlier basis vectors are subtracted prior to normalization:

```wolfram
e_2=Normalize[v_2-v_2.e_1e_1]
(* Output *)
{0.03185032246042375,-0.9901957244980933,-0.10054365764665733,-0.09159124986533951}
```

```wolfram
e_3=Normalize[v_3-v_3.e_1e_1-v_3.e_2e_2]
(* Output *)
{0.6742979135528568,0.028309479224670437,-0.531503917202115,0.5118832710326633}
```

```wolfram
e_4=Normalize[v_4-v_4.e_1e_1-v_4.e_2e_2-v_4.e_3e_3]
(* Output *)
{0.15482038834748285,-0.12917999802106084,0.7711371620395957,0.6038962268343226}
```

Confirm the answers using [Orthogonalize](https://reference.wolfram.com/language/ref/Orthogonalize.html):

```wolfram
{e_1,e_2,e_3,e_4}==Orthogonalize[{v_1,v_2,v_3,v_4}]
(* Output *)
True
```

Define a basis for $\mathbb{R}^{4}$:

```wolfram
e_1={1,1,0,0}/Sqrt[2];
e_2={1,-1,0,0}/Sqrt[2];
e_3={0,0,1,-1}/Sqrt[2];
e_4={0,0,1,1}/Sqrt[2];
```

Verify that the basis is orthonormal:

```wolfram
Table[e_i.e_j,{i,4},{j,4}]//MatrixForm
(* Output *)
({{1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}})
```

Find the components of a general vector with respect to this new basis:

```wolfram
v={w,x,y,z};
```

```wolfram
{c_1,c_2,c_3,c_4}=Table[v.e_j,{j,4}]
(* Output *)
{(w)/(Sqrt[2])+(x)/(Sqrt[2]),(w)/(Sqrt[2])-(x)/(Sqrt[2]),(y)/(Sqrt[2])-(z)/(Sqrt[2]),(y)/(Sqrt[2])+(z)/(Sqrt[2])}
```

Verify the components with respect to the $e_{j}$:

```wolfram
v==Sum[c_j e_j,{j,4}]//Simplify
(* Output *)
True
```

Define a basis for $\mathbb{R}^{5}$:

```wolfram
{b_1,b_2,b_3,b_4,b_5}={{0,1,4,-2,-5},{-4,-1,0,3,-5},{-3,-3,0,4,-4},{0,-4,-2,-3,3},{1,-2,-1,3,1}};
```

Verify that it is a basis by showing that the matrix formed by the vectors has nonzero determinant:

```wolfram
Det[{b_1,b_2,b_3,b_4,b_5}]
(* Output *)
255
```

The change of basis matrix is the inverse of the matrix whose columns are the $b_{i}$:

```wolfram
c= Inverse[Transpose[{b_1,b_2,b_3,b_4,b_5}]]
(* Output *)
{{(94)/(255),-(14)/(255),-(62)/(255),-(31)/(255),-(91)/(255)},{(227)/(255),(53)/(255),-(421)/(255),-(83)/(255),-(293)/(255)},{-(59)/(51),-(14)/(51),(91)/(51),(20)/(51),(62)/(51)},{(49)/(255),-(29)/(255),-(92)/(255),-(46)/(255),-(61)/(255)},{(278)/(255),(2)/(255),-(319)/(255),-(32)/(255),-(242)/(255)}}
```

A vector whose coordinates are $x_{i}$ in the standard bases will have coordinates $c x$ with respect to $b_{i}$:

```wolfram
{y_1,y_2,y_3,y_4,y_5}=c.{x_1,x_2,x_3,x_4,x_5}
(* Output *)
{(94 x_1)/(255)-(14 x_2)/(255)-(62 x_3)/(255)-(31 x_4)/(255)-(91 x_5)/(255),(227 x_1)/(255)+(53 x_2)/(255)-(421 x_3)/(255)-(83 x_4)/(255)-(293 x_5)/(255),-(59 x_1)/(51)-(14 x_2)/(51)+(91 x_3)/(51)+(20 x_4)/(51)+(62 x_5)/(51),(49 x_1)/(255)-(29 x_2)/(255)-(92 x_3)/(255)-(46 x_4)/(255)-(61 x_5)/(255),(278 x_1)/(255)+(2 x_2)/(255)-(319 x_3)/(255)-(32 x_4)/(255)-(242 x_5)/(255)}
```

Verify that these coordinates give back the vector $x$:

```wolfram
Sum[y_ib_i,{i,5}]//Simplify
(* Output *)
{x_1,x_2,x_3,x_4,x_5}
```

The Frenet-Serret system encodes every space curve's properties in a vector basis and scalar functions. Consider the following curve:

```wolfram
c[t_] := {Cos[t],Sin[t],(t)/(2π)}
```

Construct an orthonormal basis from the first three derivatives by subtracting parallel projections:

```wolfram
{e_1,e_2,e_3}=Table[Normalize[c^(j)[t] - ∑_{k=1}^{j-1}(c^(j)[t].c^(k)[t])/(c^(k)[t].c^(k)[t])c^(k)[t]],{j,3}]~Simplify~(t∈Reals)
(* Output *)
{{-(2 π Sin[t])/(Sqrt[1+4 π^2]),(2 π Cos[t])/(Sqrt[1+4 π^2]),(1)/(Sqrt[1+4 π^2])},{-Cos[t],-Sin[t],0},{(Sin[t])/(Sqrt[1+4 π^2]),-(Cos[t])/(Sqrt[1+4 π^2]),(2 π)/(Sqrt[1+4 π^2])}}
```

Ensure that the basis is right-handed:

```wolfram
e_3=Simplify[Det[{e_1,e_2,e_3}]]e_3;
```

Compute the curvature, $\kappa$, and torsion, $\tau$, which quantify how the curve bends:

```wolfram
{κ,τ}={(D[e_1,t].e_2)/(Sqrt[c'[t].c'[t]]),-(D[e_3,t].e_2)/(Sqrt[c'[t].c'[t]])}//Simplify
(* Output *)
{(4 π^2)/(1+4 π^2),(2 π)/(1+4 π^2)}
```

Verify the answers using [FrenetSerretSystem](https://reference.wolfram.com/language/ref/FrenetSerretSystem.html):

```wolfram
Simplify[{{κ,τ},{e_1,e_2,e_3}}==FrenetSerretSystem[c[t],t],t∈Reals]
(* Output *)
True
```

Visualize the curve and the associated moving basis, also called a frame:

```wolfram
DynamicModule[{s},Labeled[Show[ParametricPlot3D[c[t],{t,0,6π},PlotRangePadding->{.5,.5,.9}],Graphics3D[{Blue,Arrow[{c[t],c[t]+e_1}],Red,Arrow[{c[t],c[t]+e_2}],Purple,Arrow[{c[t],c[t]+e_3}]}]/. t->],Animator[,{0,6π}],Top]]
```

#### Matrices and Linear Operators

A matrix is orthogonal of $A A=Id$. Show that a rotation matrix is orthogonal:

```wolfram
r=RotationMatrix[θ,{1,1,1}]
(* Output *)
{{(1)/(3) (1+2 Cos[θ]),(1)/(3) (1-Cos[θ]-Sqrt[3] Sin[θ]),(1)/(3) (1-Cos[θ]+Sqrt[3] Sin[θ])},{(1)/(3) (1-Cos[θ]+Sqrt[3] Sin[θ]),(1)/(3) (1+2 Cos[θ]),(1)/(3) (1-Cos[θ]-Sqrt[3] Sin[θ])},{(1)/(3) (1-Cos[θ]-Sqrt[3] Sin[θ]),(1)/(3) (1-Cos[θ]+Sqrt[3] Sin[θ]),(1)/(3) (1+2 Cos[θ])}}
```

```wolfram
r.Transpose[r]//Simplify
(* Output *)
{{1,0,0},{0,1,0},{0,0,1}}
```

Confirm using [OrthogonalMatrixQ](https://reference.wolfram.com/language/ref/OrthogonalMatrixQ.html):

```wolfram
OrthogonalMatrixQ[r]
(* Output *)
True
```

A matrix is unitary of $A A^{\dagger}=Id$. Show that Pauli matrices are unitary:

```wolfram
p=PauliMatrix[2]
(* Output *)
{{0,-ⅈ},{ⅈ,0}}
```

```wolfram
p.ConjugateTranspose[p]
(* Output *)
{{1,0},{0,1}}
```

Confirm with [UnitaryMatrixQ](https://reference.wolfram.com/language/ref/UnitaryMatrixQ.html):

```wolfram
UnitaryMatrixQ[p]
(* Output *)
True
```

A matrix is normal if $A A^{\dagger}=A^{\dagger}A$. Show that the following matrix is normal:

```wolfram
m={{1,2,-1},{-1,1,2},{2,-1,1}};
MatrixForm/@{m.Transpose[m],Transpose[m].m}
(* Output *)
{({{6, -1, -1}, {-1, 6, -1}, {-1, -1, 6}}),({{6, -1, -1}, {-1, 6, -1}, {-1, -1, 6}})}
```

Confirm using [NormalMatrixQ](https://reference.wolfram.com/language/ref/NormalMatrixQ.html):

```wolfram
NormalMatrixQ[m]
(* Output *)
True
```

Normal matrices include many other types of matrices as special cases. Unitary matrices are normal:

```wolfram
p={{0,-ⅈ},{ⅈ,0}};
p . ConjugateTranspose[p]==ConjugateTranspose[p].p=={{1,0},{0,1}}
(* Output *)
True
```

Hermitian or self-adjoint matrices for which $A=A^{\dagger}$ are also normal, as the matrix $p$ shows:

```wolfram
p==ConjugateTranspose[p]
(* Output *)
True
```

However, the matrix $m$ is not a named type of normal matrix such as unitary or Hermitian:

```wolfram
{UnitaryMatrixQ[m],HermitianMatrixQ[m]}
(* Output *)
{False,False}
```

In quantum mechanics, systems with finitely many states are represented by unit vectors and physical quantities by matrices that act on them. Consider a spin-1/2 particle such as an electron. It might be in a state such as the following:

```wolfram
s={(1)/(Sqrt[5]),(2I)/(Sqrt[5])};
```

The angular momentum in the $z$ direction is given by the following matrix:

```wolfram
jz = (ℏ)/(2)PauliMatrix[3]
(* Output *)
{{(ℏ)/(2),0},{0,-(ℏ)/(2)}}
```

The angular momentum of this state is $\langle J_{z}\rangle=s^{*} J_{z}s$:

```wolfram
Conjugate[s].jz.s
(* Output *)
-(3 ℏ)/(10)
```

The uncertainty in the angular momentum of this state is $\sigma_{z}=\sqrt{\langle J_{z}^{2}\rangle-\langle J_{z}\rangle^{2}}=s^{*} J_{z}^{2}s-(s^{*} J_{z}s)^{2}$:

```wolfram
σz=Simplify[Sqrt[Conjugate[s].jz.jz.s-(Conjugate[s].jz.s)^2],ℏ>0]
(* Output *)
(2 ℏ)/(5)
```

The uncertainty in the $y$ direction is computed analogously:

```wolfram
jy = (ℏ)/(2)PauliMatrix[2]
(* Output *)
{{0,-(ⅈ ℏ)/(2)},{(ⅈ ℏ)/(2),0}}
```

```wolfram
σy=Simplify[Sqrt[Conjugate[s].jy.jy.s-(Conjugate[s].jy.s)^2],ℏ>0]
(* Output *)
(3 ℏ)/(10)
```

The uncertainty principle gives a lower bound on the product of uncertainties, $\sigma_{y} \sigma_{z}\geq \frac{1}{2} \hbar \langle J_{y} J_{z}-J_{y} J_{z}\rangle$:

```wolfram
Simplify[σy σz>(ℏ)/(2)Abs[ Conjugate[s].(jy.jz-jz.jy).s] ,ℏ>0]
(* Output *)
True
```

Consider a linear mapping $l:R^{n}\to R^{n}$:

```wolfram
n = 10;
u = Array[x_#&, n];
l = ListConvolve[{1,-2,1} n^2, u,{2,2}]
(* Output *)
{-200 x_1+100 x_2+100 x_10,100 x_1-200 x_2+100 x_3,100 x_2-200 x_3+100 x_4,100 x_3-200 x_4+100 x_5,100 x_4-200 x_5+100 x_6,100 x_5-200 x_6+100 x_7,100 x_6-200 x_7+100 x_8,100 x_7-200 x_8+100 x_9,100 x_8-200 x_9+100 x_10,100 x_1+100 x_9-200 x_10}
```

Get the matrix representation $m$ for $l$:

```wolfram
{c, m} = N[CoefficientArrays[%, u]]
(* Output *)
{SparseArray[...],SparseArray[...]}
```

Create a vector to be acted upon:

```wolfram
v = Sin[2. Pi Range[n]/n];
```

Apply the linear mapping to the vector using different methods:

```wolfram
ListConvolve[{1,-2,1} n^2, v,{2,2}]
(* Output *)
{-22.451398828979308,-36.327126400268014,-36.32712640026806,-22.45139882897927,0.,22.451398828979265,36.32712640026804,36.32712640026806,22.451398828979254,3.4638958368304884×10^-14}
```

```wolfram
l /. Thread[u->v]
(* Output *)
{-22.451398828979304,-36.327126400268014,-36.32712640026806,-22.45139882897927,0.,22.451398828979265,36.32712640026804,36.32712640026806,22.451398828979254,4.1880444608293127×10^-14}
```

```wolfram
m.v
(* Output *)
{-22.451398828979304,-36.327126400268014,-36.32712640026806,-22.45139882897927,0.,22.451398828979265,36.32712640026804,36.32712640026806,22.451398828979254,4.1880444608293127×10^-14}
```

Using $m$ with [Dot](https://reference.wolfram.com/language/ref/Dot.html) is the faster method:

```wolfram
trials = 10000;
Map[First, {AbsoluteTiming[Do[ListConvolve[{1,-2,1} n^2, v,{2,2}],{trials}]],
AbsoluteTiming[Do[l /. Thread[u->v],{trials}]],AbsoluteTiming[Do[m.v, {trials}]]}]
(* Output *)
{0.087007,0.714883,0.027799}
```

The application of a single matrix to multiple vectors $\{v_{i}\}$ can be computed as $\{v_{i}\}.m$:

```wolfram
m=RandomReal[1,{3,3}];
```

```wolfram
vs=RandomReal[1,{10^6,3}];
```

```wolfram
Map[m.#&,vs]==vs.Transpose[m]
(* Output *)
True
```

The matrix method is significantly faster than repeated application:

```wolfram
vs.Transpose[m];//AbsoluteTiming
(* Output *)
{0.010558,Null}
```

```wolfram
Map[m.#&,vs];//AbsoluteTiming
(* Output *)
{0.122327,Null}
```

#### Matrices and Arrays with Symmetry

A real symmetric matrix $s$ gives a quadratic form $q:\mathbb{R}^{n}\to \mathbb{R}$ by the formula $q(v)=v.s.v$:

```wolfram
s={{12,7,13},{7,10,-8},{13,-8,4}};
s//MatrixForm
(* Output *)
({{12, 7, 13}, {7, 10, -8}, {13, -8, 4}})
```

```wolfram
q[v_] := v.s.v
```

Quadratic forms have the property that $q(\alpha v)=\alpha^{2} v$:

```wolfram
q[α{x,y,z}]==α ^2q[{x,y,z}]//Simplify
(* Output *)
True
```

Equivalently, they define a homogeneous quadratic polynomial in the variables of $\mathbb{R}^{n}$:

```wolfram
q[{x,y,z}]//Expand
(* Output *)
12 x^2+14 x y+10 y^2+26 x z-16 y z+4 z^2
```

The range of the polynomial can be $NonNegativeReals$, $NonPositiveReals$, $\mathbb{R}$ or $\{0 \}$.  In this case it is $\mathbb{R}$:

```wolfram
FunctionRange[%,{x,y,z},q]
(* Output *)
True
```

Visualize the polynomial:

```wolfram
DensityPlot3D[q[{x,y,z}], {x,y,z}∈Cuboid[5{-1,-1,-1},5{1,1,1}]]
(* Output *)
![image](img/image_003.png)
```

A positive-definite, real symmetric matrix or metric $g$ defines an inner product by $\langle u,v \rangle=u.g.v$:

```wolfram
(g={{7,2,0},{2,6,-2},{0,-2,5}})//MatrixForm
(* Output *)
({{7, 2, 0}, {2, 6, -2}, {0, -2, 5}})
```

```wolfram
⟨u_,v_⟩ := u.g.v
```

Being positive-definite means that the associated quadratic form $q(u)=\langle u,u \rangle$ is positive for $u \neq 0$:

```wolfram
FunctionSign[{⟨{x,y,z},{x,y,z}⟩,{x,y,z}!={0,0,0}},{x,y,z},StrictInequalities->True]
(* Output *)
1
```

Note that [Dot](https://reference.wolfram.com/language/ref/Dot.html) itself is the inner product associated with the identity matrix:

```wolfram
{u,v,w}.{x,y,z} == {u,v,w}.IdentityMatrix[3].{x,y,z}
(* Output *)
True
```

Apply the Gram-Schmidt process to the standard basis to obtain an orthonormal basis:

```wolfram
{b_1,b_2,b_3} = IdentityMatrix[3];
Do[ e_j=b_j-∑_{k=1}^{j-1}⟨b_j,e_k⟩e_k; e_j=(e_j)/(Sqrt[⟨e_j,e_j⟩]),{j,3}];
{e_1,e_2,e_3}
(* Output *)
{{(1)/(Sqrt[7]),0,0},{-Sqrt[(2)/(133)],Sqrt[(7)/(38)],0},{-(2)/(9 Sqrt[19]),(7)/(9 Sqrt[19]),(Sqrt[19])/(9)}}
```

Confirm that this basis is orthonormal with respect to the inner product $\langle•,•\rangle$:

```wolfram
Table[⟨e_i,e_j⟩,{i,3},{j,3}]//Simplify
(* Output *)
{{1,0,0},{0,1,0},{0,0,1}}
```

An antisymmetric matrix $\Omega$ for which $\Omega.\Omega=-Id$ defines a Hamiltonian 2-form $\omega(u,v)=u.\Omega.v$:

```wolfram
(Ω={{0,0,-1,0},{0,0,0,-1},{1,0,0,0},{0,1,0,0}})//MatrixForm
(* Output *)
({{0, 0, -1, 0}, {0, 0, 0, -1}, {1, 0, 0, 0}, {0, 1, 0, 0}})
```

```wolfram
ω[u_,v_] := u.Ω.v
```

Verify the conditions on $\Omega$:

```wolfram
{AntisymmetricMatrixQ[Ω],Ω.Ω==-IdentityMatrix[4]}
(* Output *)
{True,True}
```

$\omega(v,v)$ is identically zero:

```wolfram
ω[{q1,q2,p1,p2},{q1,q2,p1,p2}]
(* Output *)
0
```

However, the form is nondegenerate, meaning $\omega(u,v)=0 \forall v$ implies $u=0$:

```wolfram
Solve[∀_{q3,q4,p3,p4}ω[{q1,q2,p1,p2},{q3,q4,p3,p4}]==0,{q1,q2,p1,p2}]
(* Output *)
{{q1->0,q2->0,p1->0,p2->0}}
```

Construct six vectors in dimension six:

```wolfram
{a,b,c,d,e,f}=RandomReal[1,{6,6}]
(* Output *)
{{0.11669342872150978,0.3218712579364169,0.45714022082331196,0.6758930933086318,0.060555986525428596,0.8690207449945133},{0.4033263101221545,0.9739084244717564,0.805778036232107,0.11493594286068065,0.3976796438157717,0.3510855694564283},{0.9230193463597514,0.4546516471378037,0.48158394876995714,0.3418256282433647,0.7188263003559883,0.2974475864246624},{0.7341829422219861,0.6334300742705772,0.15214949427631264,0.8350344959037792,0.9681086411257702,0.9326329014919938},{0.9022124377833491,0.7770659029268778,0.9836995326872155,0.10076450879361931,0.9553415554989404,0.5288729126022071},{0.9954308988894269,0.9284957029938059,0.8092404664997201,0.6555504088879394,0.3427951059701284,0.1889266616448022}}
```

Construct the totally antisymmetric array in dimension six using [LeviCivitaTensor](https://reference.wolfram.com/language/ref/LeviCivitaTensor.html):

```wolfram
ε=LeviCivitaTensor[6,SymmetrizedArray]
(* Output *)
SymmetrizedArray[...]
```

Compute the complete contraction $\sum_{ i_{1}\ldots i_{6}} \epsilon_{i_{1},i_{2},\ldots,i_{6}}a b \ldots f$:

```wolfram
ε.f.e.d.c.b.a
(* Output *)
0.014178764920819525
```

This is equal to the determinant of the matrix formed by the vectors:

```wolfram
Det[{a,b,c,d,e,f}]
(* Output *)
0.014178764920819435
```

By the antisymmetry of $\epsilon$, the reversed contraction differs by $(-1)^{ n (n+1)/2}$ in dimension $n$:

```wolfram
Block[{n=6}, (-1)^((n(n+1))/(2))ε.a.b.c.d.e.f]
(* Output *)
0.014178764920819556
```

### Properties & Relations

[Dot](https://reference.wolfram.com/language/ref/Dot.html) is linear in each argument:

```wolfram
u_1={x_1,y_1,z_1};
u_2={x_2,y_2,z_2};
u_3={x_3,y_3,z_3};
```

```wolfram
Expand[Dot[a u_1+u_3,u_2]==a Dot[u_1,u_2]+Dot[u_3,u_2]]
(* Output *)
True
```

```wolfram
Expand[Dot[u_1,b u_2+u_3]==b Dot[u_1,u_2]+Dot[u_1,u_3]]
(* Output *)
True
```

For a vector $v$ with real entries, [Norm](https://reference.wolfram.com/language/ref/Norm.html)[*v*] equals $\sqrt{v.v}$:

```wolfram
vr = RandomReal[1 , {3}]
(* Output *)
{0.699144551288992,0.1533754248463779,0.47869590323248024}
```

```wolfram
Norm[vr] == Sqrt[vr.vr]
(* Output *)
True
```

For a vector with complex values, the norm is given by $\sqrt{v.\overline{v}}=\sqrt{\overline{v}.v}$:

```wolfram
vc=RandomComplex[1+I,{3}]
(* Output *)
{0.504729648243269+0.15686343638621736 ⅈ,0.13147502999191074+0.4411652956578127 ⅈ,0.7245116026807861+0.03598857631843044 ⅈ}
```

```wolfram
Norm[vc]==Sqrt[vc.Conjugate[vc]]==Sqrt[Conjugate[vc].vc]
(* Output *)
True
```

For two vectors with real entries, $u_{1}.u_{2}=\left\|u_{1}\right\| \left\|u_{2}\right\|cos(\theta)$, with $\theta$ the angle between $u_{1}$ and $u_{2}$:

```wolfram
{u1,u2}=RandomReal[1,{2,5}]
(* Output *)
{{0.6548277587702362,0.21221641150313086,0.1303554574603396,0.824516345174636,0.2063680867500357},{0.5979843970273768,0.9287287240066968,0.6134504505164753,0.7708729867782815,0.4974195456793804}}
```

```wolfram
u1.u2==Norm[u1]Norm[u2]Cos[VectorAngle[u1,u2]]
(* Output *)
True
```

The scalar product of vectors is invariant under rotations:

```wolfram
{u,v}=RandomReal[1,{2,3}];
rt=RotationTransform[RandomReal[{0,2Pi}],RandomReal[{-1,1},3]];
u.v ==rt[u].rt[v]
(* Output *)
True
```

For two matrices, the $i$, $j$$^{th}$ entry of $m.n$ is the dot product of the $i$$^{th}$ row of $m$ with the $j$$^{th}$ column of $n$:

```wolfram
m=Array[x,{3,4}];
n=Array[y,{4,5}];
m.n == Table[m[[i,All]].n[[All,j]],{i,3},{j,5}]
(* Output *)
True
```

Matrix multiplication is non-commutative, $m.n \neq n.m$:

```wolfram
m=RandomInteger[{-10,10},{3,3}];
n=RandomInteger[{-10,10},{3,3}];
m.n ==n.m
(* Output *)
False
```

Use [MatrixPower](https://reference.wolfram.com/language/ref/MatrixPower.html) to compute repeated matrix products:

```wolfram
a = {{1,1,0},{0,1,1},{0,0,1}};
```

```wolfram
b=MatrixPower[a,4]
(* Output *)
{{1,4,6},{0,1,4},{0,0,1}}
```

Compare with a direct computation:

```wolfram
b==Dot[a,a,a,a]
(* Output *)
True
```

The action of `b` on a vector is the same as acting four times with `a` on that vector:

```wolfram
b.{v1,v2,v3}==Nest[x↦a.x,{v1,v2,v3},4]
(* Output *)
True
```

For two tensors $T_{i_{1}i_{2}...i_{n}}$ and $U_{j_{1}j_{2}...j_{m}}$, $T.U$ is the tensor $\sum_{k} T_{i_{1}i_{2}...i_{n-1}k} U_{k j_{2}...j_{m}}$:

```wolfram
t = RandomInteger[9, {2,3,4}];
u =RandomInteger[9, {4,5}];
t.u==Table[Sum[t[[i1, i2, k]] u[[k, j2]], {k, 4}], {i1, 2},{i2, 3}, {j2, 5}]
(* Output *)
True
```

Applying [Dot](https://reference.wolfram.com/language/ref/Dot.html) to a rank-$n$ tensor and a rank-$m$ tensor gives a rank-$n+m-2$ tensor:

```wolfram
a = RandomInteger[9, {2,3,4}];
b =RandomInteger[9, {4,5}];
```

```wolfram
TensorRank[a.b] == TensorRank[a]+TensorRank[b]-2
(* Output *)
True
```

[Dot](https://reference.wolfram.com/language/ref/Dot.html) with two arrays is a special case of [Inner](https://reference.wolfram.com/language/ref/Inner.html):

```wolfram
Inner[Times,{a,b,c},{x,y,z},Plus]
(* Output *)
a x+b y+c z
```

```wolfram
Dot[{a,b,c},{x,y,z}]
(* Output *)
a x+b y+c z
```

[Dot](https://reference.wolfram.com/language/ref/Dot.html) implements the standard inner product of arrays:

```wolfram
{{a,b},{c,d}}.{{w,x},{y,z}}
(* Output *)
{{a w+b y,a x+b z},{c w+d y,c x+d z}}
```

Use [Times](https://reference.wolfram.com/language/ref/Times.html) to do elementwise multiplication:

```wolfram
{{a,b},{c,d}}*{{w,x},{y,z}}
(* Output *)
{{a w,b x},{c y,d z}}
```

[Dot](https://reference.wolfram.com/language/ref/Dot.html) can be implemented as a combination of [TensorProduct](https://reference.wolfram.com/language/ref/TensorProduct.html) and [TensorContract](https://reference.wolfram.com/language/ref/TensorContract.html):

```wolfram
v = Array[x,{3}];
a= Array[y, {3,4,5}];
m = Array[z,{5,6}];
```

```wolfram
v.a.m==TensorContract[v⊗a⊗m,{{1,2},{4,5}}]//Simplify
(* Output *)
True
```

Use [Dot](https://reference.wolfram.com/language/ref/Dot.html) in combination with [Flatten](https://reference.wolfram.com/language/ref/Flatten.html) to contract multiple levels of one array with those of another:

```wolfram
a=Array[x,{3,4,5}];
b=Array[y,{3,4,5}];
```

```wolfram
TensorContract[a⊗b,{{2,5},{3,6}}]==Flatten[a,{{1},{2,3}}].Flatten[b,{{2,3},{1}}]
(* Output *)
True
```

[TensorReduce](https://reference.wolfram.com/language/ref/TensorReduce.html) can simplify expressions involving [Dot](https://reference.wolfram.com/language/ref/Dot.html):

```wolfram
TensorReduce[v.m.v,Assumptions->{m∈Matrices[{n,n}, Antisymmetric[{1,2}]],v∈Vectors[n]}]
(* Output *)
0
```

[Outer](https://reference.wolfram.com/language/ref/Outer.html) of two vectors can be computed with [Dot](https://reference.wolfram.com/language/ref/Dot.html):

```wolfram
u={1,2,3};
v={x,y,z};
Outer[Times,u,v]
(* Output *)
{{x,y,z},{2 x,2 y,2 z},{3 x,3 y,3 z}}
```

Construct the column and row matrices corresponding to `u` and `v`:

```wolfram
c=List/@u
(* Output *)
{{1},{2},{3}}
```

```wolfram
r={v}
(* Output *)
{{x,y,z}}
```

The outer product equals `c.r`:

```wolfram
c.r
(* Output *)
{{x,y,z},{2 x,2 y,2 z},{3 x,3 y,3 z}}
```

[Dot](https://reference.wolfram.com/language/ref/Dot.html) of a row and column matrix equals the [KroneckerProduct](https://reference.wolfram.com/language/ref/KroneckerProduct.html) of the corresponding vectors:

```wolfram
u={{1},{2},{3}};
v={{4,5,6}};
u.v==KroneckerProduct[Flatten[u],Flatten[v]]
(* Output *)
True
```

### Possible Issues

[Dot](https://reference.wolfram.com/language/ref/Dot.html) effectively treats vectors multiplied from the right as column vectors:

```wolfram
a ={{1,2},{3,4},{5,6}};
```

```wolfram
a.{1,1}
(* Output *)
{3,7,11}
```

```wolfram
a.{{1},{1}}
(* Output *)
{{3},{7},{11}}
```

[Dot](https://reference.wolfram.com/language/ref/Dot.html) effectively treats vectors multiplied from the left as row vectors:

```wolfram
{1,1,1}.a
(* Output *)
{9,12}
```

```wolfram
{{1,1,1}}.a
(* Output *)
{{9,12}}
```

[Dot](https://reference.wolfram.com/language/ref/Dot.html) does not give the standard inner product on $\mathbb{C}^{n}$:

```wolfram
a={1+I,2-I,-1-2I};
```

```wolfram
a.a
(* Output *)
2 ⅈ
```

Use [Conjugate](https://reference.wolfram.com/language/ref/Conjugate.html) in one argument to get the Hermitian inner product:

```wolfram
Conjugate[a].a
(* Output *)
12
```

Check that the result coincides with the square of the norm of `a`:

```wolfram
Norm[a]^2
(* Output *)
12
```

## Tech Notes ▪Vectors and Matrices ▪Multiplying Vectors and Matrices

## Related Guides ▪Matrix Operations ▪GPU Computing ▪Tensors ▪Operations on Vectors ▪Matrices and Linear Algebra ▪Symbolic Vectors, Matrices and Arrays ▪Finite Mathematics ▪GPU Computing with NVIDIA ▪Math & Counting Operations on Lists ▪Finite Fields ▪Structured Arrays ▪Graph Programming

## Related Links [NKS|Online](http://www.wolframscience.com/nks/search/?q=Dot) ([A New Kind of Science](http://www.wolframscience.com/nks/))

## History Introduced in 1988 (1.0) | Updated in 2003 (5.0) ▪ 2012 (9.0) ▪ 2024 (14.0)
