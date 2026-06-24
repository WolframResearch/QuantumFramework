# PolynomialGCD | [SpanFromLeft]

> [PolynomialGCD](https://reference.wolfram.com/language/ref/PolynomialGCD.html)[*poly*_1,*poly*_2,…] — gives the greatest common divisor of the polynomials `*poly*_*i*`.
> [PolynomialGCD](https://reference.wolfram.com/language/ref/PolynomialGCD.html)[*poly*_1,*poly*_2,…,[Modulus](https://reference.wolfram.com/language/ref/Modulus.html)->*p*] — evaluates the GCD modulo the prime `*p*`.

## Details and Options

In [PolynomialGCD](https://reference.wolfram.com/language/ref/PolynomialGCD.html)[*poly*_1,*poly*_2,…], all symbolic parameters are treated as variables.

[PolynomialGCD](https://reference.wolfram.com/language/ref/PolynomialGCD.html)[*poly*_1,*poly*_2,…] will by default treat algebraic numbers that appear in the `*poly*_*i*` as independent variables.

[PolynomialGCD](https://reference.wolfram.com/language/ref/PolynomialGCD.html) takes the following options:

| [Extension](https://reference.wolfram.com/language/ref/Extension.html) | [None](https://reference.wolfram.com/language/ref/None.html) | generators for the algebraic number field to be used |
| --- | --- | --- |
| [Modulus](https://reference.wolfram.com/language/ref/Modulus.html) | 0 | modulus to assume for integers |
| [Trig](https://reference.wolfram.com/language/ref/Trig.html) | [False](https://reference.wolfram.com/language/ref/False.html) | whether to do trigonometric as well as algebraic transformations |

[PolynomialGCD](https://reference.wolfram.com/language/ref/PolynomialGCD.html)[*poly*_1,*poly*_2,…,[Extension](https://reference.wolfram.com/language/ref/Extension.html)->[Automatic](https://reference.wolfram.com/language/ref/Automatic.html)] extends the coefficient field to include algebraic numbers that appear in the `*poly*_*i*`.

## Examples

### Basic Examples

Compute the greatest common divisor of polynomials:

```wolfram
PolynomialGCD[(1+x)^2(2+x)(4+x),(1+x)(2+x)(3+x)]
(* Output *)
(1+x) (2+x)
```

Compute the GCD of multivariate polynomials:

```wolfram
PolynomialGCD[(x+y)^2(y+z),(x+y)(y+z)^3]
(* Output *)
(x+y) (y+z)
```

Show that polynomials are relatively prime:

```wolfram
PolynomialGCD[x^2+4x+4,x^2+2x+1]
(* Output *)
1
```

### Scope

#### Basic Uses

The GCD of univariate polynomials:

```wolfram
PolynomialGCD[x^4-4,x^4+4 x^2+4]
(* Output *)
2+x^2
```

The GCD of multivariate polynomials:

```wolfram
PolynomialGCD[x^2+2 x y+y^2,x^3+y^3]
(* Output *)
x+y
```

The GCD of more than two polynomials:

```wolfram
PolynomialGCD[x^2-1,x^3-1,x^4-1,x^5-1,x^6-1,x^7-1]
(* Output *)
-1+x
```

The GCD of polynomials with complex coefficients:

```wolfram
PolynomialGCD[x^2+1,x^3+3I x^2-3x-I]
(* Output *)
ⅈ+x
```

#### Advanced Uses

With [Extension](https://reference.wolfram.com/language/ref/Extension.html)->[Automatic](https://reference.wolfram.com/language/ref/Automatic.html), [PolynomialGCD](https://reference.wolfram.com/language/ref/PolynomialGCD.html) detects algebraically dependent coefficients:

```wolfram
PolynomialGCD[x^2-3,x-Sqrt[3],Extension->Automatic]
(* Output *)
Sqrt[3]-x
```

Compute the GCD of polynomials over the integers modulo $3$:

```wolfram
PolynomialGCD[(x^2-1)^2,x^3-1,Modulus->3]
(* Output *)
(2+x)^2
```

Compute the GCD of polynomials over a finite field:

```wolfram
ℱ=FiniteField[17,3];
```

```wolfram
PolynomialGCD[ℱ[1]x^2+ℱ[246]x+ℱ[4436],ℱ[3]x^2+ℱ[1771]]
(* Output *)
![image](img/image_001.png)
```

With [Trig](https://reference.wolfram.com/language/ref/Trig.html)->[True](https://reference.wolfram.com/language/ref/True.html), [PolynomialGCD](https://reference.wolfram.com/language/ref/PolynomialGCD.html) recognizes dependencies between trigonometric functions:

```wolfram
PolynomialGCD[Sin[2x],1-Cos[x]^2,Trig->True]
(* Output *)
Sin[x]
```

The GCD of rational functions:

```wolfram
PolynomialGCD[(x-1)(x-2)/((x-3)(x-4)),(x-1)(x-5)/((x-3)(x-6))]
(* Output *)
(-1+x)/((-6+x) (-4+x) (-3+x))
```

Compute the GCD of two polynomials of degree $20000$:

```wolfram
rpoly[n_]:=RandomInteger[{-2^10,2^10},{n+1}].x^Range[0,n]
SeedRandom[1234];
p=rpoly[10000]; q=rpoly[10000]; r=rpoly[10000];
pp=Expand[p r]; qq=Expand[q r];
```

```wolfram
(d=PolynomialGCD[pp,qq]);//AbsoluteTiming
(* Output *)
{0.185276,Null}
```

```wolfram
Exponent[d,x]
(* Output *)
10000
```

### Options

#### Extension

By default, algebraic numbers are treated as independent variables:

```wolfram
PolynomialGCD[x^2-2,x-Sqrt[2]]
(* Output *)
1
```

With [Extension](https://reference.wolfram.com/language/ref/Extension.html)->[Automatic](https://reference.wolfram.com/language/ref/Automatic.html), [PolynomialGCD](https://reference.wolfram.com/language/ref/PolynomialGCD.html) detects algebraically dependent coefficients:

```wolfram
PolynomialGCD[x^2-2,x-Sqrt[2],Extension->Automatic]
(* Output *)
Sqrt[2]-x
```

#### Modulus

Compute the GCD over the integers modulo 2:

```wolfram
PolynomialGCD[(x+1)^3,x^3+x,Modulus->2]
(* Output *)
(1+x)^2
```

#### Trig

By default, [PolynomialGCD](https://reference.wolfram.com/language/ref/PolynomialGCD.html) treats trigonometric functions as independent variables:

```wolfram
PolynomialGCD[Sin[2x],1-Cos[x]^2]
(* Output *)
1
```

With [Trig](https://reference.wolfram.com/language/ref/Trig.html)->[True](https://reference.wolfram.com/language/ref/True.html), [PolynomialGCD](https://reference.wolfram.com/language/ref/PolynomialGCD.html) recognizes dependencies between trigonometric functions:

```wolfram
PolynomialGCD[Sin[2x],1-Cos[x]^2,Trig->True]
(* Output *)
Sin[x]
```

### Applications

Find common roots of univariate polynomials:

```wolfram
f=x^7-2x^5-x^4+5x^3+4x^2-12x+5;
g=x^7-9x^5+x^4+17x^3-7x^2-6x+3;
d=PolynomialGCD[f,g]
(* Output *)
1-2 x+x^3
```

```wolfram
x/.Solve[d==0,x]
(* Output *)
{1,(1)/(2) (-1-Sqrt[5]),(1)/(2) (-1+Sqrt[5])}
```

```wolfram
Intersection[x/.Solve[f==0,x],x/.Solve[g==0,x]]
(* Output *)
{1,(1)/(2) (-1-Sqrt[5]),(1)/(2) (-1+Sqrt[5])}
```

Find multiple roots of univariate polynomials:

```wolfram
f=x^9-7x^8+19x^7-27x^6+35x^5-77x^4+145x^3-157x^2+88x-20;
d=PolynomialGCD[f,D[f,x]]
(* Output *)
-2+5 x-4 x^2+x^3
```

```wolfram
x/.Solve[d==0,x]
(* Output *)
{1,1,2}
```

```wolfram
N[x/.Solve[f==0,x]]
(* Output *)
{1.,1.,1.,2.,2.,-1.068757255260645-1.2688873677655133 ⅈ,-1.068757255260645+1.2688873677655133 ⅈ,1.068757255260645-0.8212240798160036 ⅈ,1.068757255260645+0.8212240798160036 ⅈ}
```

### Properties & Relations

The GCD of polynomials divides the polynomials; use [PolynomialMod](https://reference.wolfram.com/language/ref/PolynomialMod.html) to prove it:

```wolfram
f=x^7-2x^5-x^4+5x^3+4x^2-12x+5;
g=x^7-9x^5+x^4+17x^3-7x^2-6x+3;
d=PolynomialGCD[f,g]
(* Output *)
1-2 x+x^3
```

```wolfram
PolynomialMod[{f,g},d]
(* Output *)
{0,0}
```

[Cancel](https://reference.wolfram.com/language/ref/Cancel.html) divides the numerator and the denominator of a rational function by their GCD:

```wolfram
Cancel[f/g]
(* Output *)
(5-2 x+x^4)/(3-7 x^2+x^4)
```

[PolynomialLCM](https://reference.wolfram.com/language/ref/PolynomialLCM.html) finds the least common multiple of polynomials:

```wolfram
PolynomialLCM[f,g]
(* Output *)
(5-2 x+x^4) (3-6 x-7 x^2+17 x^3+x^4-9 x^5+x^7)
```

```wolfram
%-(f g)/d //Together
(* Output *)
0
```

[Resultant](https://reference.wolfram.com/language/ref/Resultant.html) of two polynomials is zero if and only if their GCD has a nonzero degree:

```wolfram
Resultant[x^2-4,x^2+4x+4,x]
(* Output *)
0
```

```wolfram
PolynomialGCD[x^2-4,x^2+4x+4]
(* Output *)
2+x
```

```wolfram
Resultant[3 x+9,6x^3-3x+12,x]
(* Output *)
-3807
```

```wolfram
PolynomialGCD[3 x+9,6x^3-3x+12]
(* Output *)
3
```

[Discriminant](https://reference.wolfram.com/language/ref/Discriminant.html) of a polynomial `*f*` is zero if and only if the degree of [GCD](https://reference.wolfram.com/language/ref/GCD.html)(*f*,*f*') is nonzero:

```wolfram
f=x^3-3x+2;
g=x^3-2x+1;
```

```wolfram
Discriminant[f,x]
(* Output *)
0
```

```wolfram
PolynomialGCD[f,D[f,x]]
(* Output *)
-1+x
```

```wolfram
Discriminant[g,x]
(* Output *)
5
```

```wolfram
PolynomialGCD[g,D[g,x],x]
(* Output *)
1
```

[Discriminant](https://reference.wolfram.com/language/ref/Discriminant.html) of a polynomial `*f*` is zero if and only if the polynomial has multiple roots:

```wolfram
x/.Solve[f==0,x]
(* Output *)
{-2,1,1}
```

```wolfram
x/.Solve[g==0,x]
(* Output *)
{1,(1)/(2) (-1-Sqrt[5]),(1)/(2) (-1+Sqrt[5])}
```

## Tech Notes ▪Algebraic Operations on Polynomials ▪Polynomials Modulo Primes

## Related Guides ▪Polynomial Division ▪Theorem Proving ▪Polynomial Algebra ▪Finite Fields

## History Introduced in 1991 (2.0) | Updated in 1996 (3.0) ▪ 2022 (13.2) ▪ 2023 (13.3)
