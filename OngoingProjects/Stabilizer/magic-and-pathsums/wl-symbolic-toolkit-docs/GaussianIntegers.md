# GaussianIntegers | [SpanFromLeft]

> [GaussianIntegers](https://reference.wolfram.com/language/ref/GaussianIntegers.html) — is an option for [FactorInteger](https://reference.wolfram.com/language/ref/FactorInteger.html), [PrimeQ](https://reference.wolfram.com/language/ref/PrimeQ.html), [Factor](https://reference.wolfram.com/language/ref/Factor.html), and related functions that specifies whether factorization should be done over Gaussian integers.

## Details

With [GaussianIntegers](https://reference.wolfram.com/language/ref/GaussianIntegers.html)->[False](https://reference.wolfram.com/language/ref/False.html), factorization is done over the ordinary ring of integers $\mathbb{Z}$.

With [GaussianIntegers](https://reference.wolfram.com/language/ref/GaussianIntegers.html)->[True](https://reference.wolfram.com/language/ref/True.html), factorization is done over the ring of integers with `*i*` adjoined $\mathbb{Z}[i]$.

The Gaussian primes used when [GaussianIntegers](https://reference.wolfram.com/language/ref/GaussianIntegers.html)->[True](https://reference.wolfram.com/language/ref/True.html) are chosen to have both real and imaginary parts positive.

The first entry in the list given by [FactorInteger](https://reference.wolfram.com/language/ref/FactorInteger.html) with [GaussianIntegers](https://reference.wolfram.com/language/ref/GaussianIntegers.html)->[True](https://reference.wolfram.com/language/ref/True.html) may be `-1 or `-[I](https://reference.wolfram.com/language/ref/I.html)`.

## Examples

### Basic Examples

Factor a polynomial over $\mathbb{Q}[i]$:

```wolfram
Factor[x^2+1,GaussianIntegers->True]
(* Output *)
(-ⅈ+x) (ⅈ+x)
```

Factor an integer over $\mathbb{Z}[i]$:

```wolfram
FactorInteger[2,GaussianIntegers->True]
(* Output *)
{{-ⅈ,1},{1+ⅈ,2}}
```

### Scope

By default polynomial factorization is performed over the rationals:

```wolfram
Factor[x^4-16]
(* Output *)
(-2+x) (2+x) (4+x^2)
```

This specifies that the factorization should be done over $\mathbb{Q}[i]$:

```wolfram
Factor[x^4-16,GaussianIntegers->True]
(* Output *)
(-2+x) (-2 ⅈ+x) (2 ⅈ+x) (2+x)
```

By default integer factorization is performed over the integers:

```wolfram
FactorInteger[12]
(* Output *)
{{2,2},{3,1}}
```

This specifies that the factorization should be done over the Gaussian integers:

```wolfram
FactorInteger[12,GaussianIntegers->True]
(* Output *)
{{-1,1},{1+ⅈ,4},{3,1}}
```

A number prime over the integers may not be prime over the Gaussian integers:

```wolfram
PrimeQ[17]
(* Output *)
True
```

```wolfram
PrimeQ[17,GaussianIntegers->True]
(* Output *)
False
```

### Properties & Relations

For [Factor](https://reference.wolfram.com/language/ref/Factor.html), [GaussianIntegers](https://reference.wolfram.com/language/ref/GaussianIntegers.html)->[True](https://reference.wolfram.com/language/ref/True.html) is equivalent to [Extension](https://reference.wolfram.com/language/ref/Extension.html)->[I](https://reference.wolfram.com/language/ref/I.html):

```wolfram
Factor[x^4+1,GaussianIntegers->True]
(* Output *)
(-ⅈ+x^2) (ⅈ+x^2)
```

```wolfram
Factor[x^4+1,Extension->I]
(* Output *)
(-ⅈ+x^2) (ⅈ+x^2)
```

## Tech Notes ▪Integer and Number Theoretic Functions ▪Polynomials over Algebraic Number Fields

## Related Guides ▪Complex Numbers ▪Polynomial Algebra ▪Polynomial Factoring & Decomposition ▪Number Theory ▪Algebraic Number Theory

## History Introduced in 1991 (2.0)
