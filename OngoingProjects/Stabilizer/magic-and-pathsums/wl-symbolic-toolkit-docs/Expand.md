# Expand | [SpanFromLeft]

> [Expand](https://reference.wolfram.com/language/ref/Expand.html)[*expr*] — expands out products and positive integer powers in `*expr*`.
> [Expand](https://reference.wolfram.com/language/ref/Expand.html)[*expr*,*patt*] — leaves unexpanded any parts of `*expr*` that are free of the pattern `*patt*`.

## Details and Options

[Expand](https://reference.wolfram.com/language/ref/Expand.html) works only on positive integer powers.

[Expand](https://reference.wolfram.com/language/ref/Expand.html) applies only to the top level in `*expr*`.

[Expand](https://reference.wolfram.com/language/ref/Expand.html)[*expr*,[Modulus](https://reference.wolfram.com/language/ref/Modulus.html)->*p*] expands `*expr*` reducing the result modulo `*p*`.

[Expand](https://reference.wolfram.com/language/ref/Expand.html) automatically threads over lists in `*expr*`, as well as equations, inequalities and logic functions.

[Expand](https://reference.wolfram.com/language/ref/Expand.html) takes the following options:

| [Modulus](https://reference.wolfram.com/language/ref/Modulus.html) | 0 | modulus to assume for integers |
| --- | --- | --- |
| [Trig](https://reference.wolfram.com/language/ref/Trig.html) | [False](https://reference.wolfram.com/language/ref/False.html) | whether to do trigonometric as well as algebraic transformations |

## Examples

### Basic Examples

Expand a polynomial as a simple sum of terms:

```wolfram
Expand[(x+3)(x+2)]
(* Output *)
6+5 x+x^2
```

Expand a product of polynomials:

```wolfram
Expand[(1+x)^3(5x^2+1)(x^5+1)]
(* Output *)
1+3 x+8 x^2+16 x^3+15 x^4+6 x^5+3 x^6+8 x^7+16 x^8+15 x^9+5 x^10
```

Expand a polynomial of several variables:

```wolfram
Expand[(x+y)^5 -(x^5+y^5)]
(* Output *)
5 x^4 y+10 x^3 y^2+10 x^2 y^3+5 x y^4
```

### Scope

#### Basic Uses

Expand a polynomial:

```wolfram
Expand[(x+y)^2(x-y)^2]
(* Output *)
x^4-2 x^2 y^2+y^4
```

Expand a rational function:

```wolfram
Expand[(x+y)^2/(x-y)^2]
(* Output *)
(x^2)/((x-y)^2)+(2 x y)/((x-y)^2)+(y^2)/((x-y)^2)
```

Expand expressions involving trig functions:

```wolfram
Expand[(Sin[x]+Cos[x])^2]
(* Output *)
Cos[x]^2+2 Cos[x] Sin[x]+Sin[x]^2
```

Expand expressions involving transcendental functions:

```wolfram
Expand[(x^s+2E^-x)^3]
(* Output *)
8 ℯ^(-3 x)+12 ℯ^(-2 x) x^s+6 ℯ^-x x^(2 s)+x^(3 s)
```

Expand expressions involving arbitrary functions:

```wolfram
Expand[(f[x_1]+f[x_2])^3]
(* Output *)
f[x_1]^3+3 f[x_1]^2 f[x_2]+3 f[x_1] f[x_2]^2+f[x_2]^3
```

[Expand](https://reference.wolfram.com/language/ref/Expand.html) does not go into subexpressions:

```wolfram
Expand[Sqrt[(1+x)^2]]
(* Output *)
Sqrt[(1+x)^2]
```

In contrast, [ExpandAll](https://reference.wolfram.com/language/ref/ExpandAll.html) does:

```wolfram
ExpandAll[Sqrt[(1+x)^2]]
(* Output *)
Sqrt[1+2 x+x^2]
```

[Expand](https://reference.wolfram.com/language/ref/Expand.html) applies to the numerator only:

```wolfram
Expand[(x+y)^2/(z+y)^2]
(* Output *)
(x^2)/((y+z)^2)+(2 x y)/((y+z)^2)+(y^2)/((y+z)^2)
```

[ExpandAll](https://reference.wolfram.com/language/ref/ExpandAll.html) applies to both the numerator and denominator:

```wolfram
ExpandAll[(x+y)^2/(z+y)^2]
(* Output *)
(x^2)/(y^2+2 y z+z^2)+(2 x y)/(y^2+2 y z+z^2)+(y^2)/(y^2+2 y z+z^2)
```

Verify the equality of an expanded polynomial and its factored form:

```wolfram
1+5 x+10 x^2+10 x^3+5 x^4+x^5==Expand[(x+1)^5]
(* Output *)
True
```

[Expand](https://reference.wolfram.com/language/ref/Expand.html) accepts polynomials with real or complex coefficients as inputs:

```wolfram
Expand[(x+I)(x-I)]
(* Output *)
1+x^2
```

[Expand](https://reference.wolfram.com/language/ref/Expand.html) threads over lists:

```wolfram
Expand[Table[(x+1)^n,{n,3}]]
(* Output *)
{1+x,1+2 x+x^2,1+3 x+3 x^2+x^3}
```

[Expand](https://reference.wolfram.com/language/ref/Expand.html) threads over equations and inequalities:

```wolfram
Expand[1<(x+y)^2<2]
(* Output *)
1<x^2+2 x y+y^2<2
```

#### Advanced Uses

Apply [Expand](https://reference.wolfram.com/language/ref/Expand.html) to a polynomial over the integers modulo $4$:

```wolfram
Expand[(2y+3x)^6,Modulus->4]
(* Output *)
x^6
```

Expand a polynomial over a finite field:

```wolfram
ℱ=FiniteField[17,3];
```

```wolfram
Expand[(ℱ[123]x-ℱ[345])(ℱ[567]x-ℱ[789])]
(* Output *)
![image](img/image_001.png)
```

Apply [Expand](https://reference.wolfram.com/language/ref/Expand.html) to a trigonometric expression:

```wolfram
Expand[(Sin[2x]+x)^2,Trig->True]
(* Output *)
(1)/(2)+x^2-(Cos[x]^4)/(2)+4 x Cos[x] Sin[x]+3 Cos[x]^2 Sin[x]^2-(Sin[x]^4)/(2)
```

Polynomials with high-order powers are expanded efficiently:

```wolfram
Expand[(x+y+z)^1000]//Length
(* Output *)
501501
```

Leave parts free of `x` unexpanded:

```wolfram
Expand[(a+b)^2(1+x)^2,x]
(* Output *)
(a+b)^2+2 (a+b)^2 x+(a+b)^2 x^2
```

Leave parts free of `1+x` unexpanded:

```wolfram
Expand[(1+x)^2+(2+x)^2,1+x]
(* Output *)
1+2 x+x^2+(2+x)^2
```

Leave anything not matching `x[_]` unexpanded:

```wolfram
Expand[(a[1]+a[2])(x[1]+x[2])^2,x[_]]
(* Output *)
(a[1]+a[2]) x[1]^2+2 (a[1]+a[2]) x[1] x[2]+(a[1]+a[2]) x[2]^2
```

### Options

#### Modulus

Work in the field GF(2):

```wolfram
Expand[(1+x)^10,Modulus->2]
(* Output *)
1+x^2+x^8+x^10
```

The modulus does not have to be a prime:

```wolfram
Expand[(1+x)^10,Modulus->4]
(* Output *)
1+2 x+x^2+2 x^4+2 x^6+x^8+2 x^9+x^10
```

#### Trig

Expand a trigonometric expression:

```wolfram
Expand[Sin[x+y],Trig->True]
(* Output *)
Cos[y] Sin[x]+Cos[x] Sin[y]
```

### Applications

Calculate the characteristic polynomial of a diagonal matrix:

```wolfram
(matrix = DiagonalMatrix[{1,-1,1,-1}])//MatrixForm
(* Output *)
({{1, 0, 0, 0}, {0, -1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, -1}})
```

The characteristic polynomial is the determinant of the following matrix:

```wolfram
matrix -x*IdentityMatrix[4]//MatrixForm
(* Output *)
({{1-x, 0, 0, 0}, {0, -1-x, 0, 0}, {0, 0, 1-x, 0}, {0, 0, 0, -1-x}})
```

Expand the product of diagonal elements to get the characteristic polynomial:

```wolfram
Expand[(1-x)(-1-x)(1-x)(-1-x)]
(* Output *)
1-2 x^2+x^4
```

This can be directly computed using the [CharacteristicPolynomial](https://reference.wolfram.com/language/ref/CharacteristicPolynomial.html) function:

```wolfram
CharacteristicPolynomial[matrix ,x]
(* Output *)
1-2 x^2+x^4
```

[Cyclotomic](https://reference.wolfram.com/language/ref/Cyclotomic.html) polynomials are monic, with integer coefficients, and are irreducible over the rational numbers.

```wolfram
cyclo[n_]:=Times@@Table[If[GCD[k,n]==1,(x-(-1)^(2k/n)),1],{k,1,n}]
```

The 16th cyclotomic polynomial is given below:

```wolfram
cyclo[16]
(* Output *)
(-(-1)^(1/8)+x) ((-1)^(1/8)+x) (-(-1)^(3/8)+x) ((-1)^(3/8)+x) (-(-1)^(5/8)+x) ((-1)^(5/8)+x) (-(-1)^(7/8)+x) ((-1)^(7/8)+x)
```

Use [Expand](https://reference.wolfram.com/language/ref/Expand.html) to show that it has integer coefficients:

```wolfram
Expand[%]
(* Output *)
1+x^8
```

These polynomials can be directly computed using the [Cyclotomic](https://reference.wolfram.com/language/ref/Cyclotomic.html) function:

```wolfram
Cyclotomic[16,x]
(* Output *)
1+x^8
```

[Expand](https://reference.wolfram.com/language/ref/Expand.html) can be used to verify that two expressions are equal:

```wolfram
Cos[3x] == 4Cos[x]^3 -3Cos[x]
(* Output *)
Cos[3 x]==-3 Cos[x]+4 Cos[x]^3
```

```wolfram
Expand[%,Trig->True]
(* Output *)
True
```

### Properties & Relations

Many functions give results in unexpanded form:

```wolfram
prod=Product[x+i,{i,6}]
(* Output *)
(1+x) (2+x) (3+x) (4+x) (5+x) (6+x)
```

Apply [Expand](https://reference.wolfram.com/language/ref/Expand.html):

```wolfram
Expand[prod]
(* Output *)
720+1764 x+1624 x^2+735 x^3+175 x^4+21 x^5+x^6
```

Apply [Expand](https://reference.wolfram.com/language/ref/Expand.html) to a polynomial:

```wolfram
Expand[(1+x+y)(2-x)^3]
(* Output *)
8-4 x-6 x^2+5 x^3-x^4+8 y-12 x y+6 x^2 y-x^3 y
```

[Factor](https://reference.wolfram.com/language/ref/Factor.html) is essentially the inverse of [Expand](https://reference.wolfram.com/language/ref/Expand.html):

```wolfram
Factor[%]
(* Output *)
-(-2+x)^3 (1+x+y)
```

When no powers are involved, [Distribute](https://reference.wolfram.com/language/ref/Distribute.html) gives the same results as [Expand](https://reference.wolfram.com/language/ref/Expand.html):

```wolfram
Distribute[(1+x)(2+x)(3+x)]===Expand[(1+x)(2+x)(3+x)]
(* Output *)
True
```

Direct application of the distributive law often generates far more terms than are needed:

```wolfram
Distribute[Factor[x^6-1],Plus,Times,List,Times]
(* Output *)
{-1,-x,-x^2,x,x^2,x^3,-x^2,-x^3,-x^4,-x,-x^2,-x^3,x^2,x^3,x^4,-x^3,-x^4,-x^5,x,x^2,x^3,-x^2,-x^3,-x^4,x^3,x^4,x^5,x^2,x^3,x^4,-x^3,-x^4,-x^5,x^4,x^5,x^6}
```

Sum them:

```wolfram
Total[%]
(* Output *)
-1+x^6
```

Use [NonCommutativeExpand](https://reference.wolfram.com/language/ref/NonCommutativeExpand.html) to expand noncommutative polynomials:

```wolfram
NonCommutativeExpand[(a+2b)**(3a+4b)]
(* Output *)
3 NonCommutativeMultiply+8 NonCommutativeMultiply+4 a**b+6 b**a
```

Specify names for addition and multiplication operations:

```wolfram
alg=NonCommutativeAlgebra[<|"Multiplication"->mult,"Addition"->add|>];
```

```wolfram
NonCommutativeExpand[mult[add[a,2b],add[3a,4b]],alg]
(* Output *)
add[3 mult,4 mult[a,b],6 mult[b,a],8 mult]
```

### Neat Examples

Stylize a binomial expansion involving two emojis:

```wolfram
Style[Expand[([HappySmiley]+[SadSmiley])^5],20]
(* Output *)
[HappySmiley]^5+5 [HappySmiley]^4 [SadSmiley]+10 [HappySmiley]^3 [SadSmiley]^2+10 [HappySmiley]^2 [SadSmiley]^3+5 [HappySmiley] [SadSmiley]^4+[SadSmiley]^5
```

Expand a binomial raised to a large power:

```wolfram
Expand[(1+x)^100]
(* Output *)
1+100 x+4950 x^2+161700 x^3+3921225 x^4+75287520 x^5+1192052400 x^6+16007560800 x^7+186087894300 x^8+1902231808400 x^9+17310309456440 x^10+141629804643600 x^11+1050421051106700 x^12+7110542499799200 x^13+44186942677323600 x^14+253338471349988640 x^15+1345860629046814650 x^16+6650134872937201800 x^17+30664510802988208300 x^18+132341572939212267400 x^19+535983370403809682970 x^20+2041841411062132125600 x^21+7332066885177656269200 x^22+24865270306254660391200 x^23+79776075565900368755100 x^24+242519269720337121015504 x^25+699574816500972464467800 x^26+1917353200780443050763600 x^27+4998813702034726525205100 x^28+12410847811948286545336800 x^29+29372339821610944823963760 x^30+66324638306863423796047200 x^31+143012501349174257560226775 x^32+294692427022540894366527900 x^33+580717429720889409486981450 x^34+1095067153187962886461165020 x^35+1977204582144932989443770175 x^36+3420029547493938143902737600 x^37+5670048986634686922786117600 x^38+9013924030034630492634340800 x^39+13746234145802811501267369720 x^40+20116440213369968050635175200 x^41+28258808871162574166368460400 x^42+38116532895986727945334202400 x^43+49378235797073715747364762200 x^44+61448471214136179596720592960 x^45+73470998190814997343905056800 x^46+84413487283064039501507937600 x^47+93206558875049876949581681100 x^48+98913082887808032681188722800 x^49+100891344545564193334812497256 x^50+98913082887808032681188722800 x^51+93206558875049876949581681100 x^52+84413487283064039501507937600 x^53+73470998190814997343905056800 x^54+61448471214136179596720592960 x^55+49378235797073715747364762200 x^56+38116532895986727945334202400 x^57+28258808871162574166368460400 x^58+20116440213369968050635175200 x^59+13746234145802811501267369720 x^60+9013924030034630492634340800 x^61+5670048986634686922786117600 x^62+3420029547493938143902737600 x^63+1977204582144932989443770175 x^64+1095067153187962886461165020 x^65+580717429720889409486981450 x^66+294692427022540894366527900 x^67+143012501349174257560226775 x^68+66324638306863423796047200 x^69+29372339821610944823963760 x^70+12410847811948286545336800 x^71+4998813702034726525205100 x^72+1917353200780443050763600 x^73+699574816500972464467800 x^74+242519269720337121015504 x^75+79776075565900368755100 x^76+24865270306254660391200 x^77+7332066885177656269200 x^78+2041841411062132125600 x^79+535983370403809682970 x^80+132341572939212267400 x^81+30664510802988208300 x^82+6650134872937201800 x^83+1345860629046814650 x^84+253338471349988640 x^85+44186942677323600 x^86+7110542499799200 x^87+1050421051106700 x^88+141629804643600 x^89+17310309456440 x^90+1902231808400 x^91+186087894300 x^92+16007560800 x^93+1192052400 x^94+75287520 x^95+3921225 x^96+161700 x^97+4950 x^98+100 x^99+x^100
```

Create a nested pattern corresponding to an additive cellular automaton (rule 60):

```wolfram
Column[Table[CoefficientList[Expand[(1+x)^t,Modulus->2],x],{t,0,31}]]
(* Output *)
{{{1}}, {{1,1}}, {{1,0,1}}, {{1,1,1,1}}, {{1,0,0,0,1}}, {{1,1,0,0,1,1}}, {{1,0,1,0,1,0,1}}, {{1,1,1,1,1,1,1,1}}, {{1,0,0,0,0,0,0,0,1}}, {{1,1,0,0,0,0,0,0,1,1}}, {{1,0,1,0,0,0,0,0,1,0,1}}, {{1,1,1,1,0,0,0,0,1,1,1,1}}, {{1,0,0,0,1,0,0,0,1,0,0,0,1}}, {{1,1,0,0,1,1,0,0,1,1,0,0,1,1}}, {{1,0,1,0,1,0,1,0,1,0,1,0,1,0,1}}, {{1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1}}, {{1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1}}, {{1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1}}, {{1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1}}, {{1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1}}, {{1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1}}, {{1,1,0,0,1,1,0,0,0,0,0,0,0,0,0,0,1,1,0,0,1,1}}, {{1,0,1,0,1,0,1,0,0,0,0,0,0,0,0,0,1,0,1,0,1,0,1}}, {{1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1}}, {{1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1}}, {{1,1,0,0,0,0,0,0,1,1,0,0,0,0,0,0,1,1,0,0,0,0,0,0,1,1}}, {{1,0,1,0,0,0,0,0,1,0,1,0,0,0,0,0,1,0,1,0,0,0,0,0,1,0,1}}, {{1,1,1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1}}, {{1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1}}, {{1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1}}, {{1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1}}, {{1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1}}}
```

## Tech Notes ▪Transforming Algebraic Expressions ▪Putting Expressions into Different Forms ▪Structural Operations on Rational Expressions ▪Structural Operations on Polynomials ▪Polynomials Modulo Primes

## Related Guides ▪Algebraic Transformations ▪Formula Manipulation ▪Polynomial Factoring & Decomposition ▪Rational Functions ▪Polynomial Algebra ▪Theorem Proving ▪Finite Fields

## Related Links [NKS|Online](http://www.wolframscience.com/nks/search/?q=Expand) ([A New Kind of Science](http://www.wolframscience.com/nks/))

## History Introduced in 1988 (1.0) | Updated in 1996 (3.0) ▪ 2007 (6.0) ▪ 2023 (13.3)
