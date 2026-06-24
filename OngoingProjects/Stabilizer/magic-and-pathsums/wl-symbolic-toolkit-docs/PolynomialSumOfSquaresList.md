# PolynomialSumOfSquaresList | [SpanFromLeft]

> [PolynomialSumOfSquaresList](https://reference.wolfram.com/language/ref/PolynomialSumOfSquaresList.html)[*f*,*vars*]  — attempts to find polynomials with real coefficients `{*f*_1,…,*f*_*n*}` such that `*f*==*f*_1^(2)+…+*f*_*n*^(2)`.

## Details

[PolynomialSumOfSquaresList](https://reference.wolfram.com/language/ref/PolynomialSumOfSquaresList.html) is typically used to show that a real polynomial is non-negative.

`*f*` should be a polynomial in `*vars*` with real numeric coefficients.

## Examples

### Basic Examples

Find a sum-of-squares representation of a polynomial:

```wolfram
PolynomialSumOfSquaresList[x^2+x y+y^2,{x,y}]
(* Output *)
{(x)/(2)+y,(Sqrt[3] x)/(2)}
```

```wolfram
x^2+x y+y^2-%.%//Expand
(* Output *)
0
```

A polynomial representable as a sum of squares is non-negative:

```wolfram
Plot3D[{x^2+x y+y^2,0}, x-11, y-11, Mesh -> None, Ticks -> None]
```

*([Graphics3D])*

### Scope

Find a sum-of-squares representation of a polynomial:

```wolfram
PolynomialSumOfSquaresList[x^4+x y+y^2+1/60,{x,y}]
(* Output *)
{(1)/(2 Sqrt[15])-(31 Sqrt[15] x^2)/(122),(x)/(2)+y,(x)/(2 Sqrt[61]),(Sqrt[469] x^2)/(122)}
```

```wolfram
x^4+x y+y^2+1/60-%.%//Expand
(* Output *)
0
```

The polynomial is non-negative:

```wolfram
Plot3D[{x^4+x y+y^2+1/60,0}, x-11, y-11, Ticks -> None, Mesh -> None]
```

*([Graphics3D])*

Find a sum-of-squares representation of a polynomial with exact real numeric coefficients:

```wolfram
PolynomialSumOfSquaresList[E x^2+Sqrt[2] x y+Pi y^4+Log[2],{x,y}]
(* Output *)
{-(8 y^2)/(9 Sqrt[Log[2]])+Sqrt[Log[2]],(3)/(4) ((29)/(41)+(1)/(2) (-(58)/(41)+Sqrt[2])) x+(4 y)/(3),y^2 Sqrt[π-(64)/(81 Log[2])],Sqrt[-(9)/(16) ((29)/(41)+(1)/(2) (-(58)/(41)+Sqrt[2]))^2+ℯ] x}
```

```wolfram
E x^2+Sqrt[2] x y+Pi y^4+Log[2]-%.%//Expand
(* Output *)
0
```

Find a sum-of-squares representation of a polynomial with machine-precision real coefficients:

```wolfram
f=1.23 x^2+4.56 x y+7.89 y^4+3.21x^2z^6+9.87;
```

```wolfram
PolynomialSumOfSquaresList[f,{x,y,z}]
(* Output *)
{3.1416556144810017-1.9886443831777683 y^2,0.+0.6450033897747808 x+3.534865143570981 y,0.+1.9837574239950593 y^2,0.+0.9022032072537994 x-1.4282225979546164 x z^2,0.+1.6053329925887703 x z-1.26748099829319 x z^3,0.+1.4246537610137362 x z^2,0.+1.2662906139451948 x z^3}
```

The difference of $f$ and the sum of squares has small coefficients:

```wolfram
f-%.%//Expand
(* Output *)
-1.7763568394002505×10^-15-1.7763568394002505×10^-15 y^2+4.440892098500626×10^-16 x^2 z^2
```

Find a sum-of-squares representation of a polynomial with arbitrary-precision real coefficients:

```wolfram
f=1.23 x^2+4.56 x y+7.89 y^4+3.21x^2z^6+9.87;
```

```wolfram
PolynomialSumOfSquaresList[f,{x,y,z}]
(* Output *)
{3.14165561448100163975653338790041009239-1.98864438317776851745411336165316570491 y^2,0.64500338977478087722357109088621087856 x+3.53486514357098065293584223583750650087 y,1.98375742399505907521424200796609014453 y^2,0.9022032072537993894828329587054709162 x-1.42822259795461636349745294476825265076 x z^2,1.60533299258877037944487826778008159199 x z-1.26748099829319000172321837526103924163 x z^3,1.42465376101373627772727657198944730599 x z^2,1.26629061394519484290179810396564289762 x z^3}
```

The difference of $f$ and the sum of squares has arbitrary-precision zero coefficients:

```wolfram
f-%.%//Expand
(* Output *)
0`18.30627849418264+0`19.21069053541288 x^2+0`18.641630804187844 x y+0`17.99962144493133 y^2+0`18.403518643642855 y^4+0`18.685236882139606 x^2 z^2+0`18.486830560585812 x^2 z^4+0`18.794090614447406 x^2 z^6
```

A polynomial may not be representable as a sum of squares:

```wolfram
PolynomialSumOfSquaresList[x^4+x y+y^2,{x,y}]
(* Output *)
PolynomialSumOfSquaresList[x^4+x y+y^2,{x,y}]
```

The polynomial attains negative values, hence it cannot be a sum of squares:

```wolfram
FindInstance[x^4+x y+y^2<0,{x,y}]
(* Output *)
{{x->-(1)/(4),y->(1)/(7)}}
```

```wolfram
Show[Plot3D[{x^4+x y+y^2,0},{x,-1,1},{y,-1,1},Ticks -> None, Mesh -> None, PlotPoints -> 40],Graphics3D[{Red,Arrow[{{x,y,1},{x,y,x^4+x y+y^2}}/.%[[1]]]}]]
```

*([Graphics3D])*

[PolynomialSumOfSquaresList](https://reference.wolfram.com/language/ref/PolynomialSumOfSquaresList.html) may fail even if a sum-of-squares representation exists:

```wolfram
PolynomialSumOfSquaresList[16 x^4+16 x y+16 y^4+2,{x,y}]
(* Output *)
PolynomialSumOfSquaresList[2+16 x^4+16 x y+16 y^4,{x,y}]
```

```wolfram
Expand[16 x^4+16 x y+16 y^4+2-((4 y^2-1)^2+(2 Sqrt[2] (x+y))^2+(4 x^2-1)^2)]
(* Output *)
0
```

### Applications

Prove non-negativity of a large polynomial:

```wolfram
f=717+6 u+77 u^2+87 u^4+14 v+10 u v-6 u^3 v+43 v^2+10 u v^2+205 u^2 v^2+20 u v^3+38 v^4+26 x+66 u x+36 u^2 x-90 u^3 x-14 v x+38 u v x-14 u^2 v x+8 v^2 x-2 x^2+239 u^2 x^2-44 u v x^2+27 v^2 x^2-12 v x^3+25 x^4+62 y-184 u y+2 u^2 y+14 u^3 y-42 v y-32 u v y-94 u^2 v y+14 u v^2 y+16 v^3 y-30 x y+16 u x y-106 u^2 x y-6 v x y-38 u v x y-10 v^2 x y+6 x^2 y+32 x^3 y+83 y^2+22 u y^2+501 u^2 y^2+10 v y^2+108 u v y^2+146 v^2 y^2-14 u x y^2+48 v x y^2+42 x^2 y^2-24 y^3-6 u y^3-4 v y^3+66 y^4-14 z+154 u z+74 u^2 z-4 u^3 z+110 v z-38 u v z+10 u^2 v z-22 v^2 z+18 u v^2 z+12 v^3 z-6 x z-52 u x z+64 u^2 x z-2 v x z+10 u v x z-12 v^2 x z-4 u x^2 z+158 y z+70 u y z-66 u^2 y z+8 v y z-106 u v y z+4 v^2 y z-8 x y z-42 u x y z+36 v x y z-6 x^2 y z+6 y^2 z+60 u y^2 z-22 v y^2 z+40 z^2-10 u z^2+251 u^2 z^2-12 v z^2+100 u v z^2+251 v^2 z^2-20 x z^2+6 u x z^2+4 v x z^2+56 x^2 z^2-30 y z^2-108 u y z^2+74 v y z^2-16 x y z^2+210 y^2 z^2+8 z^3-10 u z^3+22 v z^3-40 y z^3+121 z^4;
```

Compute a sum-of-squares representation:

```wolfram
(ff=PolynomialSumOfSquaresList[f,{x,y,z,u,v}];)//AbsoluteTiming
(* Output *)
{0.4591895,Null}
```

$f$ is a sum of squares, hence it is non-negative:

```wolfram
f-ff.ff//Expand
(* Output *)
0
```

[FindInstance](https://reference.wolfram.com/language/ref/FindInstance.html) uses a method based on sum-of-squares representation here:

```wolfram
FindInstance[f<0,{x,y,z,u,v}]//AbsoluteTiming
(* Output *)
{0.2041877,{}}
```

Deciding non-negativity of `*f*` using [CylindricalDecomposition](https://reference.wolfram.com/language/ref/CylindricalDecomposition.html) takes much more time:

```wolfram
TimeConstrained[CylindricalDecomposition[f<0,{x,y,z,u,v}],300]
(* Output *)
$Aborted
```

### Properties & Relations

A polynomial that has a sum-of-squares representation is non-negative:

```wolfram
f=x^2+2y^2+3 z^2-x y-x z-y z;
```

```wolfram
PolynomialSumOfSquaresList[f, {x,y,z}]
(* Output *)
{-(x)/(2 Sqrt[3])-(y)/(2 Sqrt[3])+Sqrt[3] z,-(7 x)/(2 Sqrt[69])+(1)/(2) Sqrt[(23)/(3)] y,Sqrt[(17)/(23)] x}
```

Use [FindInstance](https://reference.wolfram.com/language/ref/FindInstance.html) to show that $f$ is non-negative:

```wolfram
FindInstance[f<0,{x,y,z}]
(* Output *)
{}
```

The Motzkin polynomial is non-negative, but is not a sum of squares:

```wolfram
f=x^4 y^2+x^2 y^4-3 x^2 y^2+1;
```

```wolfram
PolynomialSumOfSquaresList[f, {x,y}]
(* Output *)
PolynomialSumOfSquaresList[1-3 x^2 y^2+x^4 y^2+x^2 y^4,{x,y}]
```

Use [Resolve](https://reference.wolfram.com/language/ref/Resolve.html) to show that $f$ is non-negative:

```wolfram
Resolve[ForAll[{x,y},f>=0],Reals]
(* Output *)
True
```

Use [MinValue](https://reference.wolfram.com/language/ref/MinValue.html) to find the infimum of values of a polynomial:

```wolfram
f=16 x^4+16 x y+16 y^4;
```

```wolfram
m=MinValue[f, {x,y}]
(* Output *)
-2
```

[PolynomialSumOfSquaresList](https://reference.wolfram.com/language/ref/PolynomialSumOfSquaresList.html) fails for the non-negative polynomial $f-m$:

```wolfram
PolynomialSumOfSquaresList[f-m, {x,y}]
(* Output *)
PolynomialSumOfSquaresList[2+16 x^4+16 x y+16 y^4,{x,y}]
```

[PolynomialSumOfSquaresList](https://reference.wolfram.com/language/ref/PolynomialSumOfSquaresList.html) succeeds for the strictly positive polynomial $f-m+\frac{1}{100}$:

```wolfram
PolynomialSumOfSquaresList[f-m+1/100, {x,y}]
(* Output *)
{(Sqrt[201])/(10)-(2570)/(317) Sqrt[(3)/(67)] x^2+(1080)/(103) Sqrt[(3)/(67)] x y-(2570)/(317) Sqrt[(3)/(67)] y^2,(250)/(103) Sqrt[(634)/(771)] x+Sqrt[(1542)/(317)] y,-(23895386209957 x^2)/(243576460 Sqrt[1472484259])+(4163400 x y)/(103 Sqrt[1472484259])+(2)/(317) Sqrt[(21977377)/(67)] y^2,(1)/(103) Sqrt[(51724138)/(244407)] x,(565235017007310)/(103) Sqrt[(5)/(23374357768798460397416389)] x^2+(1)/(103) Sqrt[(2851379330765809)/(40987808105)] x y,(Sqrt[(133250566104133777499433909)/(2851379330765809)] x^2)/(768380)}
```

```wolfram
%.%-f+m-1/100//Expand
(* Output *)
0
```

### Possible Issues

[PolynomialSumOfSquaresList](https://reference.wolfram.com/language/ref/PolynomialSumOfSquaresList.html) may fail even if a sum-of-squares representation exists:

```wolfram
PolynomialSumOfSquaresList[16 x^4+16 x y+16 y^4+2,{x,y}]
(* Output *)
PolynomialSumOfSquaresList[2+16 x^4+16 x y+16 y^4,{x,y}]
```

```wolfram
Expand[16 x^4+16 x y+16 y^4+2-((4 y^2-1)^2+(2 Sqrt[2] (x+y))^2+(4 x^2-1)^2)]
(* Output *)
0
```

## Related Guides ▪Polynomial Algebra

## History Introduced in 2021 (13.0)
