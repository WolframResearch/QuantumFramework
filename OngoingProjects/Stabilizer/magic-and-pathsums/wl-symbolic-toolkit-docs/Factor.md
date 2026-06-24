# Factor | [SpanFromLeft]

> [Factor](https://reference.wolfram.com/language/ref/Factor.html)[*poly*] — factors a polynomial over the integers.
> [Factor](https://reference.wolfram.com/language/ref/Factor.html)[*poly*,[Modulus](https://reference.wolfram.com/language/ref/Modulus.html)->*p*] — factors a polynomial modulo the prime `*p*`.
> [Factor](https://reference.wolfram.com/language/ref/Factor.html)[*poly*,[Extension](https://reference.wolfram.com/language/ref/Extension.html)->{*a*_1,*a*_2,…}] — factors a polynomial allowing coefficients that are rational combinations of the algebraic numbers `*a*_*i*`.

## Details and Options

[Factor](https://reference.wolfram.com/language/ref/Factor.html) applies only to the top algebraic level in an expression. You may have to use [Map](https://reference.wolfram.com/language/ref/Map.html), or apply [Factor](https://reference.wolfram.com/language/ref/Factor.html) again, to reach other levels.

[Factor](https://reference.wolfram.com/language/ref/Factor.html)[*poly*,[GaussianIntegers](https://reference.wolfram.com/language/ref/GaussianIntegers.html)->[True](https://reference.wolfram.com/language/ref/True.html)] factors allowing Gaussian integer coefficients.

If any coefficients in `*poly*` are complex numbers, factoring is done allowing Gaussian integer coefficients.

The exponents of variables need not be positive integers. [Factor](https://reference.wolfram.com/language/ref/Factor.html) can deal with exponents that are linear combinations of symbolic expressions.

When given a rational expression, [Factor](https://reference.wolfram.com/language/ref/Factor.html) effectively first calls [Together](https://reference.wolfram.com/language/ref/Together.html), then factors numerator and denominator.

[Factor](https://reference.wolfram.com/language/ref/Factor.html) takes the following options:

| [Extension](https://reference.wolfram.com/language/ref/Extension.html) | [None](https://reference.wolfram.com/language/ref/None.html) | coefficient field to be used |
| --- | --- | --- |
| [GaussianIntegers](https://reference.wolfram.com/language/ref/GaussianIntegers.html) | [False](https://reference.wolfram.com/language/ref/False.html) | whether to factor over Gaussian integers |
| [Modulus](https://reference.wolfram.com/language/ref/Modulus.html) | 0 | modulus to assume for integers |
| [Trig](https://reference.wolfram.com/language/ref/Trig.html) | [False](https://reference.wolfram.com/language/ref/False.html) | whether to do trigonometric as well as algebraic transformations  |

With the default setting [Extension](https://reference.wolfram.com/language/ref/Extension.html)->[None](https://reference.wolfram.com/language/ref/None.html), [Factor](https://reference.wolfram.com/language/ref/Factor.html)[*poly*] will treat algebraic number coefficients in `*poly*` like independent variables.

[Factor](https://reference.wolfram.com/language/ref/Factor.html)[*poly*,[Extension](https://reference.wolfram.com/language/ref/Extension.html)->[Automatic](https://reference.wolfram.com/language/ref/Automatic.html)] will extend the domain of coefficients to include any algebraic numbers that appear in `*poly*`.

[Factor](https://reference.wolfram.com/language/ref/Factor.html) automatically threads over lists, as well as equations, inequalities and logic functions.

## Examples

### Basic Examples

Factor univariate polynomials:

```wolfram
Factor[1+2x+x^2]
(* Output *)
(1+x)^2
```

```wolfram
Factor[x^10-1]
(* Output *)
(-1+x) (1+x) (1-x+x^2-x^3+x^4) (1+x+x^2+x^3+x^4)
```

Factor multivariate polynomials:

```wolfram
Factor[x^10-y^10]
(* Output *)
(x-y) (x+y) (x^4-x^3 y+x^2 y^2-x y^3+y^4) (x^4+x^3 y+x^2 y^2+x y^3+y^4)
```

Factor polynomials over the integers modulo 2:

```wolfram
Factor[x^10-1,Modulus->2]
(* Output *)
(1+x)^2 (1+x+x^2+x^3+x^4)^2
```

### Scope

#### Basic Uses

A univariate polynomial:

```wolfram
Factor[x^3-6x^2+11x-6]
(* Output *)
(-3+x) (-2+x) (-1+x)
```

A multivariate polynomial:

```wolfram
Factor[2x^3 y-2a^2x y-3a^2x^2+3a^4]
(* Output *)
(a-x) (a+x) (3 a^2-2 x y)
```

A rational function:

```wolfram
Factor[(x^3+2x^2)/(x^2-4y^2)-(x+2)/(x^2-4y^2)]
(* Output *)
((-1+x) (1+x) (2+x))/((x-2 y) (x+2 y))
```

A non-polynomial expression:

```wolfram
Factor[2^(3x)-1]
(* Output *)
(-1+2^x) (1+2^x+2^(2 x))
```

[Factor](https://reference.wolfram.com/language/ref/Factor.html) threads over lists:

```wolfram
Factor[{x^2-1,x^4-1,x^8-1}]
(* Output *)
{(-1+x) (1+x),(-1+x) (1+x) (1+x^2),(-1+x) (1+x) (1+x^2) (1+x^4)}
```

[Factor](https://reference.wolfram.com/language/ref/Factor.html) threads over equations and inequalities:

```wolfram
Factor[1<1+2 x+x^2+1/(1+x)<2]
(* Output *)
1<((2+x) (1+x+x^2))/(1+x)<2
```

#### Advanced Uses

Factor a polynomial over the Gaussian integers:

```wolfram
Factor[x^2+1,GaussianIntegers->True]
(* Output *)
(-ⅈ+x) (ⅈ+x)
```

Factor a polynomial over an algebraic extension:

```wolfram
Factor[x^4-2,Extension->Sqrt[2]]
(* Output *)
-(Sqrt[2]-x^2) (Sqrt[2]+x^2)
```

Factor a polynomial over the integers modulo 3:

```wolfram
Factor[x^3+1,Modulus->3]
(* Output *)
(1+x)^3
```

Factor polynomials over a finite field:

```wolfram
ℱ=FiniteField[17,3];
```

```wolfram
Factor[ℱ[3]x^2+ℱ[1771]]
(* Output *)
![image](img/image_001.png)
```

```wolfram
Factor[x^3+5x+19,Extension->ℱ]
(* Output *)
![image](img/image_003.png)
```

Factor a polynomial over an extension of a finite field:

```wolfram
ℱ=FiniteField[29,3];
𝒢=FiniteField[29,6];
ℰ=FiniteFieldEmbedding[ℱ,𝒢];
```

A polynomial irreducible over $\mathcal{F}$ factors after embedding $\mathcal{F}$ in a larger field $\mathcal{G}$:

```wolfram
Factor[ℱ[12]x^2+ℱ[34]x+ℱ[56]]
(* Output *)
![image](img/image_005.png)
```

```wolfram
Factor[ℱ[12]x^2+ℱ[34]x+ℱ[56],Extension->ℰ]
(* Output *)
![image](img/image_007.png)
```

Some non-polynomial expressions can be factored:

```wolfram
Factor[x^(2s)+2x^s+1]
(* Output *)
(1+x^s)^2
```

```wolfram
Factor[x^(2/3s)+2x^s+1]
(* Output *)
(1+x^(s/3)) (1-x^(s/3)+2 x^(2 s/3))
```

```wolfram
Factor[Exp[2s]+2Exp[s]+1]
(* Output *)
(1+ℯ^s)^2
```

Factor a polynomial of degree $2000$:

```wolfram
rpoly[n_]:=RandomInteger[{-2^10,2^10},{n+1}].x^Range[0,n]
SeedRandom[1234];
p=rpoly[1000]; q=rpoly[1000];
r=Expand[p q];
```

```wolfram
Factor[r]//Short//AbsoluteTiming
(* Output *)
{1.1077352,-((982+597 x+348 x^2-459 x^3+<<1505>>+161 x^996-155 x^997-695 x^998-86 x^999+460 x^1000) (<<1497>>+<<1>>))}
```

### Options

#### Extension

Factor over algebraic number fields:

```wolfram
Factor[1+x^4, Extension -> Sqrt[2]]
(* Output *)
-(-1+Sqrt[2] x-x^2) (1+Sqrt[2] x+x^2)
```

```wolfram
Factor[1+x^4,Extension->{Sqrt[2],I}]
(* Output *)
(1)/(4) (Sqrt[2]-(1+ⅈ) x) (Sqrt[2]-(1-ⅈ) x) (Sqrt[2]+(1-ⅈ) x) (Sqrt[2]+(1+ⅈ) x)
```

[Extension](https://reference.wolfram.com/language/ref/Extension.html)->[Automatic](https://reference.wolfram.com/language/ref/Automatic.html) automatically extends to a field that covers the coefficients:

```wolfram
Factor[2+2Sqrt[2]x+x^2]
(* Output *)
2+2 Sqrt[2] x+x^2
```

```wolfram
Factor[2+2Sqrt[2]x+x^2,Extension->Automatic]
(* Output *)
(Sqrt[2]+x)^2
```

Factor a polynomial with integer coefficients over a finite field:

```wolfram
ℱ=FiniteField[73,4];
```

```wolfram
Factor[x^4+21x+3,Extension->ℱ]
(* Output *)
![image](img/image_009.png)
```

Factor a polynomial with coefficients in a finite field:

```wolfram
ℱ=FiniteField[7,2];
```

```wolfram
Factor[ℱ[1]x^4+ℱ[234]x+ℱ[567]]
(* Output *)
![image](img/image_011.png)
```

Embedding $\mathcal{F}$ in a larger field $\mathcal{G}$ allows further factorization:

```wolfram
𝒢=FiniteField[7,4];
ℰ=FiniteFieldEmbedding[ℱ,𝒢];
```

```wolfram
Factor[ℱ[1]x^4+ℱ[234]x+ℱ[567],Extension->ℰ]
(* Output *)
![image](img/image_013.png)
```

#### GaussianIntegers

Factor over Gaussian integers:

```wolfram
Factor[1 + x^2, GaussianIntegers -> True]
(* Output *)
(-ⅈ+x) (ⅈ+x)
```

```wolfram
Factor[1 + x^2, Extension ->I]
(* Output *)
(-ⅈ+x) (ⅈ+x)
```

#### Modulus

Factor over finite fields:

```wolfram
Factor[5+x^12]
(* Output *)
5+x^12
```

```wolfram
Factor[5+x^12,Modulus->7]
(* Output *)
(2+x^3) (5+x^3) (4+x^6)
```

#### Trig

Factor a trigonometric expression:

```wolfram
Factor[Sin[x]+Sin[y],Trig->True]
(* Output *)
2 Cos[(x)/(2)-(y)/(2)] Sin[(x)/(2)+(y)/(2)]
```

### Applications

When modeling behavior with polynomials, it is important to determine when the polynomial evaluates to zero. For example, suppose the cost to produce a video game system is modeled by the following expression:

```wolfram
cost[x_]:=-.03x^3+4x^2-x+500
```

```wolfram
Plot[cost[x],{x,0,100}]
```

*([Graphics])*

Also suppose the revenue can be modeled by the equation:

```wolfram
revenue[x_]:=.03x^2-4x+200
```

```wolfram
Plot[revenue[x],{x,0,100}]
```

*([Graphics])*

If you wish to know the number of units you must sell before making a profit, calculate the difference:

```wolfram
profit[x_]:=revenue[x]-cost[x]
```

```wolfram
Plot[profit[x],{x,0,100}]
```

*([Graphics])*

Then solve to find where the profit function is zero using [Factor](https://reference.wolfram.com/language/ref/Factor.html):

```wolfram
Factor[profit[x]]
(* Output *)
0.030000000000000006 (-133.6415125480683+1. x) (74.82704894112294+1.3081792147349498 x+1. x^2)
```

This reveals to us there is a zero for profit at $x=133.642$:

```wolfram
Plot[profit[x],{x,130,140}]
```

*([Graphics])*

Find a number that is equal to its square:

```wolfram
x^2==x
(* Output *)
x^2==x
```

Subtract $x$ from both sides of the equation:

```wolfram
SubtractSides[%,x]
(* Output *)
-x+x^2==0
```

Use [Factor](https://reference.wolfram.com/language/ref/Factor.html) to find when a polynomial is zero:

```wolfram
Factor[%]
(* Output *)
(-1+x) x==0
```

The only numbers that are equal to their square are thus $x=0$ and $x=1$:

```wolfram
0^2
(* Output *)
0
```

```wolfram
1^2
(* Output *)
1
```

Compute the greatest common divisor of two polynomials:

```wolfram
Factor[-x^2+x]
(* Output *)
-(-1+x) x
```

```wolfram
Factor[x^4-1]
(* Output *)
(-1+x) (1+x) (1+x^2)
```

You can see they share a common factor of $x-1$. Confirm this result using [PolynomialGCD](https://reference.wolfram.com/language/ref/PolynomialGCD.html):

```wolfram
PolynomialGCD[-x^2+x,x^4-1]
(* Output *)
-1+x
```

### Properties & Relations

[Expand](https://reference.wolfram.com/language/ref/Expand.html) is effectively the inverse of [Factor](https://reference.wolfram.com/language/ref/Factor.html):

```wolfram
Factor[x^50-1]
(* Output *)
(-1+x) (1+x) (1-x+x^2-x^3+x^4) (1+x+x^2+x^3+x^4) (1-x^5+x^10-x^15+x^20) (1+x^5+x^10+x^15+x^20)
```

```wolfram
Expand[Factor[x^50-1]]
(* Output *)
-1+x^50
```

[FactorList](https://reference.wolfram.com/language/ref/FactorList.html) gives a list of factors:

```wolfram
FactorList[x^8+11x^7+43x^6+59x^5-35x^4-151x^3-63x^2+81x+54]
(* Output *)
{{1,1},{-1+x,2},{1+x,2},{2+x,1},{3+x,3}}
```

[FactorSquareFree](https://reference.wolfram.com/language/ref/FactorSquareFree.html) only pulls out multiple factors:

```wolfram
f=x^9+9x^8+21x^7-27x^6-153x^5-81x^4+239x^3+207x^2-108x-108;
```

```wolfram
Factor[f]
(* Output *)
(-2+x) (-1+x)^2 (1+x)^2 (2+x) (3+x)^3
```

```wolfram
FactorSquareFree[f]
(* Output *)
(3+x)^3 (-4+x^2) (-1+x^2)^2
```

### Neat Examples

The first factoring of $x^{n}-1$ where a `2 appears as a coefficient:

```wolfram
Factor[x^105-1]
(* Output *)
(-1+x) (1+x+x^2) (1+x+x^2+x^3+x^4) (1+x+x^2+x^3+x^4+x^5+x^6) (1-x+x^3-x^4+x^5-x^7+x^8) (1-x+x^3-x^4+x^6-x^8+x^9-x^11+x^12) (1-x+x^5-x^6+x^7-x^8+x^10-x^11+x^12-x^13+x^14-x^16+x^17-x^18+x^19-x^23+x^24) (1+x+x^2-x^5-x^6-2 x^7-x^8-x^9+x^12+x^13+x^14+x^15+x^16+x^17-x^20-x^22-x^24-x^26-x^28+x^31+x^32+x^33+x^34+x^35+x^36-x^39-x^40-2 x^41-x^42-x^43+x^46+x^47+x^48)
```

```wolfram
ListLinePlot[Table[Length[Factor[x^n-1]],{n,200}]]
```

*([Graphics])*

## Tech Notes ▪Transforming Algebraic Expressions ▪Putting Expressions into Different Forms ▪Structural Operations on Rational Expressions ▪Structural Operations on Polynomials ▪Polynomials Modulo Primes ▪Implementation notes: Algebra and Calculus

## Related Guides ▪Algebraic Transformations ▪Polynomial Factoring & Decomposition ▪Rational Functions ▪Precollege Education ▪Formula Manipulation ▪Finite Fields ▪Algebraic Number Theory ▪Polynomial Algebra

## Related Links [NKS|Online](http://www.wolframscience.com/nks/search/?q=Factor) ([A New Kind of Science](http://www.wolframscience.com/nks/))

## History Introduced in 1988 (1.0) | Updated in 1996 (3.0) ▪ 2007 (6.0) ▪ 2022 (13.2) ▪ 2023 (13.3)
