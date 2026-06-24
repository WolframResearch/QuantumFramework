# FactorList | [SpanFromLeft]

> [FactorList](https://reference.wolfram.com/language/ref/FactorList.html)[*poly*] — gives a list of the factors of a polynomial, together with their exponents.

## Details and Options

The first element of the list is always the overall numerical factor. It is `{1,1}` if there is no overall numerical factor.

[FactorList](https://reference.wolfram.com/language/ref/FactorList.html) takes the following options:

| [Extension](https://reference.wolfram.com/language/ref/Extension.html) | [None](https://reference.wolfram.com/language/ref/None.html) | coefficient field to be used |
| --- | --- | --- |
| [GaussianIntegers](https://reference.wolfram.com/language/ref/GaussianIntegers.html) | [False](https://reference.wolfram.com/language/ref/False.html) | whether to allow Gaussian integer coefficients |
| [Modulus](https://reference.wolfram.com/language/ref/Modulus.html) | 0 | modulus to assume for integers |
| [Trig](https://reference.wolfram.com/language/ref/Trig.html) | [False](https://reference.wolfram.com/language/ref/False.html) | whether to do trigonometric as well as algebraic transformations |

[FactorList](https://reference.wolfram.com/language/ref/FactorList.html)[*poly*,[Modulus](https://reference.wolfram.com/language/ref/Modulus.html)->*p*] requires that `*p*` be prime.

[FactorList](https://reference.wolfram.com/language/ref/FactorList.html)[*poly*,[Extension](https://reference.wolfram.com/language/ref/Extension.html)->{*a*_1,*a*_2,…}] allows coefficients that are arbitrary rational combinations of the `*a*_*i*`.

## Examples

### Basic Examples

List the irreducible factors of polynomials:

```wolfram
FactorList[x^2-1]
(* Output *)
{{1,1},{-1+x,1},{1+x,1}}
```

```wolfram
FactorList[2x^3+2x^2-2x-2]
(* Output *)
{{2,1},{-1+x,1},{1+x,2}}
```

List the irreducible factors of multivariate polynomials:

```wolfram
FactorList[x^3+2x^2y+x y^2+x^2z-y^2z-x z^2-2y z^2-z^3]
(* Output *)
{{1,1},{x-z,1},{x+y+z,2}}
```

### Scope

#### Basic Uses

A univariate polynomial:

```wolfram
FactorList[x^3-6x^2+11x-6]
(* Output *)
{{1,1},{-3+x,1},{-2+x,1},{-1+x,1}}
```

A multivariate polynomial:

```wolfram
FactorList[2x^3 y-2a^2x y-3a^2x^2+3a^4]
(* Output *)
{{1,1},{a-x,1},{a+x,1},{3 a^2-2 x y,1}}
```

A polynomial with multiple factors:

```wolfram
FactorList[x^8+11x^7+43x^6+59x^5-35x^4-151x^3-63x^2+81x+54]
(* Output *)
{{1,1},{-1+x,2},{1+x,2},{2+x,1},{3+x,3}}
```

A rational function:

```wolfram
FactorList[(x^3+2x^2)/(x^2-4y^2)-(x+2)/(x^2-4y^2)]
(* Output *)
{{1,1},{-1+x,1},{1+x,1},{2+x,1},{x-2 y,-1},{x+2 y,-1}}
```

#### Advanced Uses

The irreducible factors over the Gaussian integers:

```wolfram
FactorList[x^2+1, GaussianIntegers -> True]
(* Output *)
{{1,1},{-ⅈ+x,1},{ⅈ+x,1}}
```

The irreducible factors over an algebraic extension:

```wolfram
FactorList[x^4-2,Extension->Sqrt[2]]
(* Output *)
{{-1,1},{Sqrt[2]-x^2,1},{Sqrt[2]+x^2,1}}
```

The irreducible factors over the integers modulo 2:

```wolfram
FactorList[x^4+1,Modulus->2]
(* Output *)
{{1,1},{1+x,4}}
```

The irreducible factors over a finite field:

```wolfram
ℱ=FiniteField[17,3];
```

```wolfram
FactorList[ℱ[3]x^2+ℱ[1771]]
(* Output *)
![image](img/image_001.png)
```

```wolfram
FactorList[x^3+5x+19,Extension->ℱ]
(* Output *)
![image](img/image_003.png)
```

The irreducible factors over an extension of a finite field:

```wolfram
ℱ=FiniteField[29,3];
𝒢=FiniteField[29,6];
ℰ=FiniteFieldEmbedding[ℱ,𝒢];
```

A polynomial irreducible over $\mathcal{F}$ factors after embedding $\mathcal{F}$ in a larger field $\mathcal{G}$:

```wolfram
FactorList[ℱ[12]x^2+ℱ[34]x+ℱ[56]]
(* Output *)
![image](img/image_005.png)
```

```wolfram
FactorList[ℱ[12]x^2+ℱ[34]x+ℱ[56],Extension->ℰ]
(* Output *)
![image](img/image_007.png)
```

List factors of non-polynomial expressions:

```wolfram
FactorList[E^(2 x)-2E^x+1]
(* Output *)
{{1,1},{-1+ℯ^x,2}}
```

```wolfram
FactorList[x^2-Sqrt[x]]
(* Output *)
{{1,1},{-1+Sqrt[x],1},{Sqrt[x],1},{1+Sqrt[x]+x,1}}
```

List factors of a polynomial of degree $2000$:

```wolfram
rpoly[n_]:=RandomInteger[{-2^10,2^10},{n+1}].x^Range[0,n]
SeedRandom[1234];
p=rpoly[1000]; q=rpoly[1000];
r=Expand[p q];
```

```wolfram
(facs=FactorList[r]);//AbsoluteTiming
(* Output *)
{1.3238984,Null}
```

```wolfram
{Exponent[#[[1]],x],#[[2]]}&/@facs
(* Output *)
{{0,1},{1000,1},{1000,1}}
```

### Options

#### Extension

Factor over algebraic number fields:

```wolfram
FactorList[1+x^4, Extension -> Sqrt[2]]
(* Output *)
{{-1,1},{-1+Sqrt[2] x-x^2,1},{1+Sqrt[2] x+x^2,1}}
```

```wolfram
FactorList[1+x^4,Extension->{Sqrt[2],I}]
(* Output *)
{{1,1},{4,-1},{Sqrt[2]+(1+ⅈ) x,1},{Sqrt[2]+(1-ⅈ) x,1},{Sqrt[2]-(1-ⅈ) x,1},{Sqrt[2]-(1+ⅈ) x,1}}
```

[Extension](https://reference.wolfram.com/language/ref/Extension.html)->[Automatic](https://reference.wolfram.com/language/ref/Automatic.html) automatically extends to a field that covers the coefficients:

```wolfram
FactorList[2+2Sqrt[2]x+x^2]
(* Output *)
{{1,1},{2+2 Sqrt[2] x+x^2,1}}
```

```wolfram
FactorList[2+2Sqrt[2]x+x^2,Extension->Automatic]
(* Output *)
{{1,1},{Sqrt[2]+x,2}}
```

Factor a polynomial with integer coefficients over a finite field:

```wolfram
ℱ=FiniteField[73,3];
```

```wolfram
FactorList[x^3+21x+3,Extension->ℱ]
(* Output *)
![image](img/image_009.png)
```

Factor a polynomial with coefficients in a finite field:

```wolfram
ℱ=FiniteField[7,2];
```

```wolfram
FactorList[ℱ[1]x^4+ℱ[234]x+ℱ[567]]
(* Output *)
![image](img/image_011.png)
```

Embedding $\mathcal{F}$ in a larger field $\mathcal{G}$ allows further factorization:

```wolfram
𝒢=FiniteField[7,4];
ℰ=FiniteFieldEmbedding[ℱ,𝒢];
```

```wolfram
FactorList[ℱ[1]x^4+ℱ[234]x+ℱ[567],Extension->ℰ]
(* Output *)
![image](img/image_013.png)
```

#### GaussianIntegers

Factor over Gaussian integers:

```wolfram
FactorList[1 + x^2, GaussianIntegers -> True]
(* Output *)
{{1,1},{-ⅈ+x,1},{ⅈ+x,1}}
```

```wolfram
FactorList[1 + x^2, Extension ->I]
(* Output *)
{{1,1},{-ⅈ+x,1},{ⅈ+x,1}}
```

#### Modulus

Factor over finite fields:

```wolfram
FactorList[5+x^12]
(* Output *)
{{1,1},{5+x^12,1}}
```

```wolfram
FactorList[5+x^12,Modulus->7]
(* Output *)
{{1,1},{2+x^3,1},{5+x^3,1},{4+x^6,1}}
```

#### Trig

Factor a trigonometric expression:

```wolfram
FactorList[Sin[x]+Sin[y],Trig->True]
(* Output *)
{{2,1},{Cos[(x)/(2)-(y)/(2)],1},{Sin[(x)/(2)+(y)/(2)],1}}
```

### Applications

[FactorList](https://reference.wolfram.com/language/ref/FactorList.html) can be useful in determining the behavior of functions:

```wolfram
f=x^3-x^2-16x-20;
```

```wolfram
FactorList[f]
(* Output *)
{{1,1},{-5+x,1},{2+x,2}}
```

$f$ has roots at $-2$ and $5$, at $-2$ it does not cross the $x$ axis, and at $5$ it crosses the $x$ axis:

```wolfram
Plot[-20-16 x-x^2+x^3,{x,-4,6}]
```

*([Graphics])*

Show that a polynomial is a perfect square:

```wolfram
f=62500-112500 x+15625 x^2+47500 x^3-3500 x^4-9380 x^5-1106 x^6+628 x^7+208 x^8+24 x^9+x^10;
```

Compute the factors and note that the constant factor is positive:

```wolfram
fac = FactorList[f]
(* Output *)
{{1,1},{-2+x,2},{-1+x,2},{5+x,6}}
```

Extract the exponents of nonconstant factors:

```wolfram
exp=Transpose[Rest[fac]][[2]]
(* Output *)
{2,2,6}
```

Check that every factor has even exponent and thus $f$ is a square:

```wolfram
AllTrue[exp,EvenQ]
(* Output *)
True
```

### Properties & Relations

[FactorList](https://reference.wolfram.com/language/ref/FactorList.html) gives a list of irreducible factors:

```wolfram
FactorList[x^8+11x^7+43x^6+59x^5-35x^4-151x^3-63x^2+81x+54]
(* Output *)
{{1,1},{-1+x,2},{1+x,2},{2+x,1},{3+x,3}}
```

This multiplies the factors together:

```wolfram
Times@@Power@@@%
(* Output *)
(-1+x)^2 (1+x)^2 (2+x) (3+x)^3
```

[Factor](https://reference.wolfram.com/language/ref/Factor.html) gives a product of factors:

```wolfram
Factor[x^8+11x^7+43x^6+59x^5-35x^4-151x^3-63x^2+81x+54]
(* Output *)
(-1+x)^2 (1+x)^2 (2+x) (3+x)^3
```

[Expand](https://reference.wolfram.com/language/ref/Expand.html) combines all the factors back together:

```wolfram
Expand[%]
(* Output *)
54+81 x-63 x^2-151 x^3-35 x^4+59 x^5+43 x^6+11 x^7+x^8
```

[FactorSquareFreeList](https://reference.wolfram.com/language/ref/FactorSquareFreeList.html) gives a list of square-free factors:

```wolfram
FactorSquareFreeList[x^8+11x^7+43x^6+59x^5-35x^4-151x^3-63x^2+81x+54]
(* Output *)
{{1,1},{2+x,1},{3+x,3},{-1+x^2,2}}
```

[FactorInteger](https://reference.wolfram.com/language/ref/FactorInteger.html) gives a list of prime factors of an integer:

```wolfram
FactorInteger[1234567]
(* Output *)
{{127,1},{9721,1}}
```

```wolfram
Times@@Power@@@%
(* Output *)
1234567
```

## Tech Notes ▪Algebraic Operations on Polynomials

## Related Guides ▪Polynomial Factoring & Decomposition ▪Polynomial Algebra ▪Finite Fields

## History Introduced in 1988 (1.0) | Updated in 1996 (3.0) ▪ 2022 (13.2) ▪ 2023 (13.3)
