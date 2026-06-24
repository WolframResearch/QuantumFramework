# PolynomialMod | [SpanFromLeft]

> [PolynomialMod](https://reference.wolfram.com/language/ref/PolynomialMod.html)[*poly*,*m*] — gives the polynomial `*poly*` reduced modulo `*m*`.
> [PolynomialMod](https://reference.wolfram.com/language/ref/PolynomialMod.html)[*poly*,{*m*_1,*m*_2,…}] — reduces modulo all of the `*m*_*i*`.

## Details and Options

[PolynomialMod](https://reference.wolfram.com/language/ref/PolynomialMod.html)[*poly*,*m*] for integer `*m*` gives a polynomial in which all coefficients are reduced modulo `*m*`.

When `*m*` is a polynomial, [PolynomialMod](https://reference.wolfram.com/language/ref/PolynomialMod.html)[*poly*,*m*] reduces `*poly*` by subtracting polynomial multiples of `*m*`, to give a result with minimal degree and leading coefficient.

[PolynomialMod](https://reference.wolfram.com/language/ref/PolynomialMod.html) gives results according to a definite convention; other conventions could yield results differing by multiples of `*m*`.

Unlike [PolynomialRemainder](https://reference.wolfram.com/language/ref/PolynomialRemainder.html), [PolynomialMod](https://reference.wolfram.com/language/ref/PolynomialMod.html) never performs divisions in generating its results.

[PolynomialMod](https://reference.wolfram.com/language/ref/PolynomialMod.html) takes the following options:

| [CoefficientDomain](https://reference.wolfram.com/language/ref/CoefficientDomain.html) | Rationals | the type of objects assumed to be coefficients |
| --- | --- | --- |
| [Modulus](https://reference.wolfram.com/language/ref/Modulus.html) | 0 | modulus to assume for integers |
| [Trig](https://reference.wolfram.com/language/ref/Trig.html) | [False](https://reference.wolfram.com/language/ref/False.html) | whether to do trigonometric as well as algebraic transformations |

## Examples

### Basic Examples

Reduce a polynomial modulo 2:

```wolfram
PolynomialMod[3x^2+2x+1,2]
(* Output *)
1+x^2
```

Reduce a polynomial modulo another polynomial:

```wolfram
PolynomialMod[3x^2+2x+1,x^2+1]
(* Output *)
-2+2 x
```

### Scope

Reduce a polynomial modulo an integer:

```wolfram
PolynomialMod[3x^3+21x^2-7x+55,9]
(* Output *)
1+2 x+3 x^2+3 x^3
```

Reduce a multivariate polynomial modulo an integer:

```wolfram
PolynomialMod[35x^3+21x^2y^2-17x y^3+55z-123,19]
(* Output *)
10+16 x^3+2 x^2 y^2+2 x y^3+17 z
```

Reduce a polynomial modulo a polynomial:

```wolfram
PolynomialMod[3x^3+21x^2-7x+55,2x^2-7]
(* Output *)
(257)/(2)+(7 x)/(2)
```

The difference of the original polynomial and the result is divisible by the modulus:

```wolfram
Cancel[(3x^3+21x^2-7x+55-%)/(2x^2-7)]
(* Output *)
(3 (7+x))/(2)
```

Reduce a polynomial modulo a polynomial with complex coefficients:

```wolfram
PolynomialMod[3x^3+21x^2-7x+55,x^2+I]
(* Output *)
(55-21 ⅈ)-(7+3 ⅈ) x
```

Reduce a polynomial modulo a polynomial and an integer:

```wolfram
PolynomialMod[3x^3+21x^2-7x+55,{2x^2-7,9}]
(* Output *)
7+8 x
```

Reduce a polynomial modulo two polynomials and an integer:

```wolfram
PolynomialMod[3x^3+21x^2y^2-7x y^3+55,{2x^2-7,x y-3, 9}]
(* Output *)
1+7 x+x^3+4 y^2
```

### Options

#### CoefficientDomain

With the default `CoefficientDomain->[Rationals](https://reference.wolfram.com/language/ref/Rationals.html)`, integer coefficients can be inverted:

```wolfram
PolynomialMod[x^2,2x-1]
(* Output *)
(1)/(4)
```

With `CoefficientDomain->[Integers](https://reference.wolfram.com/language/ref/Integers.html)`, [PolynomialMod](https://reference.wolfram.com/language/ref/PolynomialMod.html) does not invert integer coefficients:

```wolfram
PolynomialMod[x^2,2x-1,CoefficientDomain->Integers]
(* Output *)
x^2
```

#### Modulus

Reduce a polynomial modulo a polynomial over the integers modulo 3:

```wolfram
PolynomialMod[x^2,2x-1,Modulus->3]
(* Output *)
1
```

### Applications

Reduce all coefficients of a polynomial modulo an integer:

```wolfram
PolynomialMod[1234(x+y)^7-(x-y)^5,7]
(* Output *)
6 x^5+2 x^7+5 x^4 y+4 x^3 y^2+3 x^2 y^3+2 x y^4+y^5+2 y^7
```

### Properties & Relations

For univariate rational polynomials, [PolynomialRemainder](https://reference.wolfram.com/language/ref/PolynomialRemainder.html) is the same as [PolynomialMod](https://reference.wolfram.com/language/ref/PolynomialMod.html):

```wolfram
PolynomialMod[2x^3+3x^2+4x+5,6x^2-7]
(* Output *)
(17)/(2)+(19 x)/(3)
```

```wolfram
PolynomialRemainder[2x^3+3x^2+4x+5,6x^2-7,x]
(* Output *)
(17)/(2)+(19 x)/(3)
```

[PolynomialRemainder](https://reference.wolfram.com/language/ref/PolynomialRemainder.html) considers all polynomials to be univariate in the specified variable:

```wolfram
PolynomialRemainder[x^2,x+a,x]
(* Output *)
a^2
```

For multivariate polynomials, [PolynomialMod](https://reference.wolfram.com/language/ref/PolynomialMod.html) picks its own variable order:

```wolfram
PolynomialMod[x^2,x+a]
(* Output *)
x^2
```

The main variable here is `a`:

```wolfram
PolynomialMod[a^2,x+a]
(* Output *)
x^2
```

[PolynomialRemainder](https://reference.wolfram.com/language/ref/PolynomialRemainder.html) considers parameters to be invertible:

```wolfram
PolynomialRemainder[x^2+a^2,a x-1,x]
(* Output *)
(1)/(a^2)+a^2
```

[PolynomialMod](https://reference.wolfram.com/language/ref/PolynomialMod.html) does not invert symbolic expressions:

```wolfram
PolynomialMod[x^2+a^2,a x-1]
(* Output *)
a^2+x^2
```

## Tech Notes ▪Algebraic Operations on Polynomials ▪Polynomials Modulo Primes

## Related Guides ▪Polynomial Algebra ▪Polynomial Division ▪Cryptographic Number Theory

## Related Links [NKS|Online](http://www.wolframscience.com/nks/search/?q=PolynomialMod) ([A New Kind of Science](http://www.wolframscience.com/nks/))

## History Introduced in 1991 (2.0) | Updated in 2003 (5.0)
