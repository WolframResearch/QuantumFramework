# PolynomialLCM | [SpanFromLeft]

> [PolynomialLCM](https://reference.wolfram.com/language/ref/PolynomialLCM.html)[*poly*_1,*poly*_2,…] — gives the least common multiple of the polynomials `*poly*_*i*`.
> [PolynomialLCM](https://reference.wolfram.com/language/ref/PolynomialLCM.html)[*poly*_1,*poly*_2,…,[Modulus](https://reference.wolfram.com/language/ref/Modulus.html)->*p*] — evaluates the LCM modulo the prime `*p*`.

## Details and Options

[PolynomialLCM](https://reference.wolfram.com/language/ref/PolynomialLCM.html)[*poly*_1,*poly*_2,…] will by default treat algebraic numbers that appear in the `*poly*_*i*` as independent variables.

[PolynomialLCM](https://reference.wolfram.com/language/ref/PolynomialLCM.html) takes the following options:

| [Extension](https://reference.wolfram.com/language/ref/Extension.html) | [None](https://reference.wolfram.com/language/ref/None.html) | generators for the algebraic number field to be used |
| --- | --- | --- |
| [Modulus](https://reference.wolfram.com/language/ref/Modulus.html) | 0 | modulus to assume for integers |
| [Trig](https://reference.wolfram.com/language/ref/Trig.html) | [False](https://reference.wolfram.com/language/ref/False.html) | whether to do trigonometric as well as algebraic transformations |

[PolynomialLCM](https://reference.wolfram.com/language/ref/PolynomialLCM.html)[*poly*_1,*poly*_2,…,[Extension](https://reference.wolfram.com/language/ref/Extension.html)->[Automatic](https://reference.wolfram.com/language/ref/Automatic.html)] extends the coefficient field to include algebraic numbers that appear in the `*poly*_*i*`.

## Examples

### Basic Examples

Compute the least common multiple (LCM) of polynomials:

```wolfram
PolynomialLCM[(1+x)^2(2+x)(4+x),(1+x)(2+x)(3+x)]
(* Output *)
(1+x)^2 (2+x) (3+x) (4+x)
```

Compute the least common multiple of several polynomials:

```wolfram
PolynomialLCM[(1+x)^2,(1+x)(2+x),(3+x)]
(* Output *)
(1+x)^2 (2+x) (3+x)
```

Compute the least common multiple of multivariate polynomials:

```wolfram
PolynomialLCM[(1+x)^2(2+y),(1+x)(2+y)(x+y)]
(* Output *)
(1+x)^2 (2+y) (x+y)
```

### Scope

#### Basic Uses

The LCM of univariate polynomials:

```wolfram
PolynomialLCM[x^4-4,x^4+4 x^2+4]
(* Output *)
(-2+x^2) (4+4 x^2+x^4)
```

The LCM of multivariate polynomials:

```wolfram
PolynomialLCM[x^2+2 x y+y^2,x^3+y^3]
(* Output *)
(x+y) (x^3+y^3)
```

The LCM of more than two polynomials:

```wolfram
PolynomialLCM[x^2-1,x^3-1,x^4-1,x^5-1,x^6-1,x^7-1]
(* Output *)
(-1+x) (1+x) (1+x^2) (1-x+x^2) (1+x+x^2) (1+x+x^2+x^3+x^4) (1+x+x^2+x^3+x^4+x^5+x^6)
```

The LCM of rational functions:

```wolfram
PolynomialLCM[(x-1)(x-2)/(x-4),(x-1)/((x-4)(x-6))]
(* Output *)
((-2+x) (-1+x))/(-4+x)
```

#### Advanced Uses

With [Extension](https://reference.wolfram.com/language/ref/Extension.html)->[Automatic](https://reference.wolfram.com/language/ref/Automatic.html), [PolynomialLCM](https://reference.wolfram.com/language/ref/PolynomialLCM.html) detects algebraically dependent coefficients:

```wolfram
PolynomialLCM[x^2-2,x-Sqrt[2],Extension->Automatic]
(* Output *)
-2+x^2
```

Compute the LCM over the integers modulo $3$:

```wolfram
PolynomialLCM[(x+2)^3,x^3+2x,Modulus->3]
(* Output *)
(1+x+x^2) (2 x+x^3)
```

Compute the LCM of polynomials over a finite field:

```wolfram
ℱ=FiniteField[17,3];
```

```wolfram
PolynomialLCM[ℱ[1]x^2+ℱ[246]x+ℱ[4436],ℱ[3]x^2+ℱ[1771]]
(* Output *)
![image](img/image_001.png)
```

With [Trig](https://reference.wolfram.com/language/ref/Trig.html)->[True](https://reference.wolfram.com/language/ref/True.html), [PolynomialLCM](https://reference.wolfram.com/language/ref/PolynomialLCM.html) recognizes identities between trigonometric functions:

```wolfram
PolynomialLCM[Sin[x/2]^2-1,Cos[x],Trig->True]
(* Output *)
Cos[x] (1+Cos[x])
```

The LCM of rational functions:

```wolfram
PolynomialLCM[(x-1)(x-2)/((x-3)(x-4)),(x-1)(x-5)/((x-3)(x-6))]
(* Output *)
((-5+x) (-2+x) (-1+x))/(-3+x)
```

### Options

#### Extension

By default, algebraic numbers are treated as independent variables:

```wolfram
PolynomialLCM[x^2-2,x-Sqrt[2]]
(* Output *)
(-Sqrt[2]+x) (-2+x^2)
```

With [Extension](https://reference.wolfram.com/language/ref/Extension.html)->[Automatic](https://reference.wolfram.com/language/ref/Automatic.html), [PolynomialLCM](https://reference.wolfram.com/language/ref/PolynomialLCM.html) detects algebraically dependent coefficients:

```wolfram
PolynomialLCM[x^2-2,x-Sqrt[2],Extension->Automatic]
(* Output *)
-2+x^2
```

#### Modulus

Compute the LCM over the integers modulo 2:

```wolfram
PolynomialLCM[(x+1)^3,x^3+x,Modulus->2]
(* Output *)
(1+x) (x+x^3)
```

#### Trig

By default, [PolynomialLCM](https://reference.wolfram.com/language/ref/PolynomialLCM.html) treats trigonometric functions as independent variables:

```wolfram
PolynomialLCM[Sin[2x],1-Cos[x]^2]
(* Output *)
(1-Cos[x]^2) Sin[2 x]
```

With [Trig](https://reference.wolfram.com/language/ref/Trig.html)->[True](https://reference.wolfram.com/language/ref/True.html), [PolynomialLCM](https://reference.wolfram.com/language/ref/PolynomialLCM.html) recognizes dependencies between trigonometric functions:

```wolfram
PolynomialLCM[Sin[2x],1-Cos[x]^2,Trig->True]
(* Output *)
-2 Cos[x] Sin[x]^2
```

### Applications

If $p$ divides $q$, then their least common multiple is equal to $q$:

```wolfram
p=x-1;
q=x^3-1;
```

```wolfram
PolynomialLCM[p,q]
(* Output *)
-1+x^3
```

If $p$ and $q$ are relatively prime, then their least common multiple is equal to $p q$:

```wolfram
p=x+1;
q=x^3-1;
```

```wolfram
PolynomialLCM[p,q]
(* Output *)
(1+x) (-1+x^3)
```

In general, the least common multiple of $p$ and $q$ is $p q$ divided by the greatest common divisor of $p$ and $q$:

```wolfram
p=x^6-1;
q=x^4-2x^2+1;
```

```wolfram
p q/PolynomialGCD[p,q]
(* Output *)
((1-2 x^2+x^4) (-1+x^6))/(-1+x^2)
```

```wolfram
PolynomialLCM[p,q]
(* Output *)
(1-2 x^2+x^4) (1+x^2+x^4)
```

Use [Together](https://reference.wolfram.com/language/ref/Together.html) to prove the equality:

```wolfram
Together[%%-%]
(* Output *)
0
```

Compute the LCM of the first five cyclotomic polynomials. Notice the coefficients are anti-palindromic:

```wolfram
PolynomialLCM@@Table[Cyclotomic[n,x],{n,0,5}]//Expand
(* Output *)
-1-2 x-3 x^2-3 x^3-2 x^4+2 x^6+3 x^7+3 x^8+2 x^9+x^10
```

This results from the fact that every cyclotomic polynomial is palindromic except the first:

```wolfram
Cyclotomic[105,x]
(* Output *)
1+x+x^2-x^5-x^6-2 x^7-x^8-x^9+x^12+x^13+x^14+x^15+x^16+x^17-x^20-x^22-x^24-x^26-x^28+x^31+x^32+x^33+x^34+x^35+x^36-x^39-x^40-2 x^41-x^42-x^43+x^46+x^47+x^48
```

The first cyclotomic polynomial is anti-palindromic:

```wolfram
Cyclotomic[1,x]
(* Output *)
-1+x
```

Thus when taking the product of palindromic polynomials with one anti-palindromic polynomial, we will always obtain an anti-palindromic polynomial:

```wolfram
PolynomialLCM@@Table[Cyclotomic[n,x],{n,0,15}]//Expand
(* Output *)
-1-x-3 x^2-5 x^3-9 x^4-14 x^5-23 x^6-33 x^7-49 x^8-68 x^9-94 x^10-124 x^11-164 x^12-207 x^13-261 x^14-319 x^15-386 x^16-455 x^17-532 x^18-606 x^19-684 x^20-755 x^21-823 x^22-878 x^23-925 x^24-952 x^25-965 x^26-955 x^27-926 x^28-872 x^29-799 x^30-702 x^31-588 x^32-456 x^33-312 x^34-158 x^35+158 x^37+312 x^38+456 x^39+588 x^40+702 x^41+799 x^42+872 x^43+926 x^44+955 x^45+965 x^46+952 x^47+925 x^48+878 x^49+823 x^50+755 x^51+684 x^52+606 x^53+532 x^54+455 x^55+386 x^56+319 x^57+261 x^58+207 x^59+164 x^60+124 x^61+94 x^62+68 x^63+49 x^64+33 x^65+23 x^66+14 x^67+9 x^68+5 x^69+3 x^70+x^71+x^72
```

### Properties & Relations

The LCM of polynomials is divisible by the polynomials; use [PolynomialMod](https://reference.wolfram.com/language/ref/PolynomialMod.html) to prove it:

```wolfram
f=x^7-2x^5-x^4+5x^3+4x^2-12x+5;
g=x^7-9x^5+x^4+17x^3-7x^2-6x+3;
m=PolynomialLCM[f,g]
(* Output *)
(5-2 x+x^4) (3-6 x-7 x^2+17 x^3+x^4-9 x^5+x^7)
```

```wolfram
{PolynomialMod[m,f],PolynomialMod[m,g]}
(* Output *)
{0,0}
```

[PolynomialGCD](https://reference.wolfram.com/language/ref/PolynomialGCD.html) finds the greatest common divisor of polynomials:

```wolfram
PolynomialGCD[f,g]
(* Output *)
1-2 x+x^3
```

```wolfram
m-(f g)/% //Together
(* Output *)
0
```

## Tech Notes ▪Algebraic Operations on Polynomials

## Related Guides ▪Polynomial Division ▪Finite Fields

## History Introduced in 1991 (2.0) | Updated in 1996 (3.0) ▪ 2022 (13.2) ▪ 2023 (13.3)
