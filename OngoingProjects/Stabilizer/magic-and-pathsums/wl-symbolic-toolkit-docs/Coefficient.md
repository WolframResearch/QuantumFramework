# Coefficient | [SpanFromLeft]

> [Coefficient](https://reference.wolfram.com/language/ref/Coefficient.html)[*expr*,*form*] — gives the coefficient of `*form*` in the polynomial `*expr*`.
> [Coefficient](https://reference.wolfram.com/language/ref/Coefficient.html)[*expr*,*form*,*n*] — gives the coefficient of `*form*^*n*` in `*expr*`.

## Details and Options

[Coefficient](https://reference.wolfram.com/language/ref/Coefficient.html) picks only terms that contain the particular form specified. $x^{2}$ is not considered part of $x^{3}$.

`*form*` can be a product of powers.

[Coefficient](https://reference.wolfram.com/language/ref/Coefficient.html)[*expr*,*form*,0] picks out terms that are not proportional to `*form*`.

[Coefficient](https://reference.wolfram.com/language/ref/Coefficient.html) works whether or not `*expr*` is explicitly given in expanded form.

## Examples

### Basic Examples

Find coefficients of polynomials:

```wolfram
Coefficient[(x+1)^3,x, 2]
(* Output *)
3
```

```wolfram
Coefficient[(x+y)^4,x y^3]
(* Output *)
4
```

### Scope

Find a coefficient at `*x*`:

```wolfram
Coefficient[a x+b y+c,x]
(* Output *)
a
```

Find a coefficient at a power of `*x*`:

```wolfram
Coefficient[a x^3+b x^2+c x+d,x,2]
(* Output *)
b
```

Find the free term in a polynomial:

```wolfram
Coefficient[(x+2)^2+(x+3)^3,x,0]
(* Output *)
31
```

Find a coefficient at a multivariate monomial:

```wolfram
Coefficient[(x+y)(x+2y)(3x+4y+5),x y^2]
(* Output *)
18
```

### Options

#### Modulus

Find a coefficient over the integers modulo 2:

```wolfram
Coefficient[(x+1)^3,x, 2,Modulus->2]
(* Output *)
1
```

### Properties & Relations

[CoefficientList](https://reference.wolfram.com/language/ref/CoefficientList.html) gives a list of all polynomial coefficients:

```wolfram
f=(x+3)^5;
```

```wolfram
CoefficientList[f,x]
(* Output *)
{243,405,270,90,15,1}
```

The same list of coefficients obtained using [Coefficient](https://reference.wolfram.com/language/ref/Coefficient.html) and [Exponent](https://reference.wolfram.com/language/ref/Exponent.html):

```wolfram
Coefficient[f,x,#]&/@Range[0,Exponent[f,x]]
(* Output *)
{243,405,270,90,15,1}
```

For multivariate polynomials [CoefficientList](https://reference.wolfram.com/language/ref/CoefficientList.html) gives a tensor of the coefficients:

```wolfram
f=(3x+5y)^4;
```

```wolfram
cl=CoefficientList[f,{x, y}]
(* Output *)
{{0,0,0,0,625},{0,0,0,1500,0},{0,0,1350,0,0},{0,540,0,0,0},{81,0,0,0,0}}
```

[CoefficientArrays](https://reference.wolfram.com/language/ref/CoefficientArrays.html) gives the list of arrays of polynomial coefficients ordered by total degree:

```wolfram
ca=CoefficientArrays[f,{x,y}]
(* Output *)
{0,SparseArray[...],SparseArray[...],SparseArray[...],SparseArray[...]}
```

The coefficient of `*x* *y*^(3)`:

```wolfram
Coefficient[f,x y^3]
(* Output *)
1500
```

In `cl` the coefficient of `x^*a* y^*b*` is the element at position `*{**a*+1,*b*+1}`:

```wolfram
cl[[1+1,1+3]]
(* Output *)
1500
```

In `ca` the position of this coefficient is `*a*+*b*+1 followed by `*a*` `1s and `*b*` `2s (`1 and `2 indicate the first and second variables):

```wolfram
ca[[5,1,2,2,2]]
(* Output *)
1500
```

### Possible Issues

[Coefficient](https://reference.wolfram.com/language/ref/Coefficient.html) treats transcendental powers as being algebraically unrelated to algebraic powers:

```wolfram
Coefficient[x^s x,x^s]
(* Output *)
x
```

[Coefficient](https://reference.wolfram.com/language/ref/Coefficient.html) treats distinct transcendental powers as being algebraically unrelated to one another:

```wolfram
Coefficient[x^s x^t,x^s]
(* Output *)
x^t
```

## Tech Notes ▪Picking Out Pieces of Algebraic Expressions ▪Finding the Structure of a Polynomial

## Related Guides ▪Polynomial Algebra ▪Series Expansions ▪Formula Manipulation

## Related Links [MathWorld](https://mathworld.wolfram.com/CoefficientNotation.html)

## History Introduced in 1988 (1.0) | Updated in 1996 (3.0)
