# PrimitivePolynomialQ | [SpanFromLeft]

> [PrimitivePolynomialQ](https://reference.wolfram.com/language/ref/PrimitivePolynomialQ.html)[*poly*,*p*]  — tests whether `*poly*` is a primitive polynomial modulo a prime `*p*`.

## Details

The polynomial `*poly*` must be univariate.

## Examples

### Basic Examples

Test whether a polynomial is primitive modulo 13:

```wolfram
PrimitivePolynomialQ[x^4+3x^3+2x^2+x+7,13]
(* Output *)
True
```

This polynomial can be factored modulo 2, and therefore it is not primitive:

```wolfram
PrimitivePolynomialQ[x^2+1,2]
(* Output *)
False
```

### Scope

Test for primitivity of a univariate polynomial modulo a prime:

```wolfram
PrimitivePolynomialQ[x^2+5x-1,7]
(* Output *)
False
```

Polynomials can be given in non-expanded form:

```wolfram
PrimitivePolynomialQ[(x+1)(x+3)+2,7]
(* Output *)
True
```

Coefficients of the polynomial do not have to be integers:

```wolfram
PrimitivePolynomialQ[x^2/2+2x+5/2,7]
(* Output *)
True
```

### Properties & Relations

A polynomial must be irreducible in order to be primitive:

```wolfram
IrreduciblePolynomialQ[x^4+3x^3+2x^2+x+7,Modulus->13]
(* Output *)
True
```

Irreducibility is a necessary but not-sufficient condition for a polynomial to be primitive:

```wolfram
IrreduciblePolynomialQ[x^4+3x^3+2x^2+x+9,Modulus->13]
(* Output *)
True
```

```wolfram
PrimitivePolynomialQ[x^4+3x^3+2x^2+x+9,13]
(* Output *)
False
```

A trinomial whose order is a Mersenne prime exponent is primitive modulo 2 if and only if it is irreducible:

```wolfram
MersennePrimeExponent[10]
(* Output *)
89
```

```wolfram
poly=x^89+x^38+1;
```

```wolfram
PrimitivePolynomialQ[poly,2]
(* Output *)
True
```

```wolfram
IrreduciblePolynomialQ[poly,Modulus->2]
(* Output *)
True
```

Primitivity of a polynomial depends on the choice of prime:

```wolfram
poly=x^4+3x^3+2x^2+x+7;
```

```wolfram
PrimitivePolynomialQ[poly,17]
(* Output *)
True
```

```wolfram
PrimitivePolynomialQ[poly,23]
(* Output *)
False
```

## Related Guides ▪Polynomial Algebra ▪Number Theory

## History Introduced in 2017 (11.1)
