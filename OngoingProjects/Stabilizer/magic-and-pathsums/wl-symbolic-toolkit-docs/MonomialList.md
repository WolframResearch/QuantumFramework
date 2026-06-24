# MonomialList | [SpanFromLeft]

> [MonomialList](https://reference.wolfram.com/language/ref/MonomialList.html)[*poly*] — gives the list of all monomials in the polynomial `*poly*`.
> [MonomialList](https://reference.wolfram.com/language/ref/MonomialList.html)[*poly*,{*x*_1,*x*_2,…}]  — gives the list of monomials with respect to the variables `*x*_*i*` in poly.
> [MonomialList](https://reference.wolfram.com/language/ref/MonomialList.html)[*poly*,{*x*_1,*x*_2,…},*order*] — puts the monomials in the specified order.

## Details and Options

[MonomialList](https://reference.wolfram.com/language/ref/MonomialList.html) works whether or not `*poly*` is explicitly given in expanded form.

[MonomialList](https://reference.wolfram.com/language/ref/MonomialList.html)[*poly*] is equivalent to [MonomialList](https://reference.wolfram.com/language/ref/MonomialList.html)[*poly*,[Variables](https://reference.wolfram.com/language/ref/Variables.html)[*poly*]].

Possible settings for `*order*` are `"Lexicographic"`, `"DegreeLexicographic"`, `"DegreeReverseLexicographic"`, `"NegativeLexicographic"`, `"NegativeDegreeLexicographic"`, `"NegativeDegreeReverseLexicographic"`, or an explicit weight matrix.

Monomials are sorted on the basis of their exponent vectors with respect to the variables `*x*_*i*`.

`"NegativeLexicographic"` corresponds to applying [Sort](https://reference.wolfram.com/language/ref/Sort.html) to the list of exponent vectors.

`"Lexicographic"` gives the reverse of `"NegativeLexicographic"`, and is the default for [MonomialList](https://reference.wolfram.com/language/ref/MonomialList.html).

`"DegreeLexicographic"` sorts first with respect to total degree, then by using the ordering defined by `"Lexicographic"`.

`"DegreeReverseLexicographic"` sorts first with respect to total degree, then in the negative lexicographic order by starting from the last variable.

`"NegativeDegreeLexicographic"` and `"NegativeDegreeReverseLexicographic"` sort from lower to higher total degree.

An explicit weight matrix `*w*` defines an ordering given by `"Lexicographic"` ordering of the `*w*.*v*_*i*`, where the `*v*_*i*` are the exponent vectors.

[MonomialList](https://reference.wolfram.com/language/ref/MonomialList.html)[*poly*,*vars*,[Modulus](https://reference.wolfram.com/language/ref/Modulus.html)->*m*] computes the coefficients modulo `*m*`.

[MonomialList](https://reference.wolfram.com/language/ref/MonomialList.html)[*poly*,[All](https://reference.wolfram.com/language/ref/All.html),*order*] is equivalent to [MonomialList](https://reference.wolfram.com/language/ref/MonomialList.html)[*poly*,[Variables](https://reference.wolfram.com/language/ref/Variables.html)[*poly*],*order*].

## Examples

### Basic Examples

Get the list of monomials:

```wolfram
MonomialList[(x+y)^3]
(* Output *)
{x^3,3 x^2 y,3 x y^2,y^3}
```

### Scope

Use `"DegreeLexicographic"` monomial ordering:

```wolfram
MonomialList[x^2 y^2+x^3,{x,y},"DegreeLexicographic"]
(* Output *)
{x^2 y^2,x^3}
```

Specify the same ordering using a weight matrix:

```wolfram
MonomialList[x^2 y^2+x^3,{x,y},{{1,1},{1,0}}]
(* Output *)
{x^2 y^2,x^3}
```

### Options

#### Modulus

Reduce the coefficients modulo 2:

```wolfram
MonomialList[(x+1)^5,x,Modulus->2]
(* Output *)
{x^5,x^4,x,1}
```

### Properties & Relations

[Plus](https://reference.wolfram.com/language/ref/Plus.html) or [Total](https://reference.wolfram.com/language/ref/Total.html) reconstructs the original polynomial:

```wolfram
MonomialList[a x^2+b x y+c y^2,{x,y}]
(* Output *)
{a x^2,b x y,c y^2}
```

```wolfram
Plus@@%
(* Output *)
a x^2+b x y+c y^2
```

[CoefficientRules](https://reference.wolfram.com/language/ref/CoefficientRules.html) gives a different representation:

```wolfram
CoefficientRules[a x^2+b x y+c y^2,{x,y}]
(* Output *)
{{2,0}->a,{1,1}->b,{0,2}->c}
```

Obtain `"NegativeDegreeReverseLexicographic"` from `"DegreeLexicographic"`:

```wolfram
poly=Sum[RandomInteger[{-10,10}]Times@@({x,y,z}^RandomInteger[10,3]),{10}]
(* Output *)
-x^3 y z-3 x^9 y^7 z+6 x^10 y^8 z-5 x y^10 z-x^2 y^5 z^4-9 x^10 y^10 z^4-2 x^8 z^6-9 x y^6 z^6-3 x^5 y^4 z^9+3 x^9 y^4 z^9
```

```wolfram
MonomialList[poly,{x,y,z},"NegativeDegreeReverseLexicographic"]===Reverse@MonomialList[poly,Reverse@{x,y,z},"DegreeLexicographic"]
(* Output *)
True
```

### Possible Issues

The list given by [Variables](https://reference.wolfram.com/language/ref/Variables.html)[*poly*] is not always sorted:

```wolfram
Variables[y+x z]
(* Output *)
{y,x,z}
```

```wolfram
MonomialList[y+x z]
(* Output *)
{y,x z}
```

```wolfram
MonomialList[y+x z,{x,y,z}]
(* Output *)
{x z,y}
```

## Tech Notes ▪Polynomial Orderings

## Related Guides ▪Polynomial Factoring & Decomposition ▪Polynomial Algebra

## History Introduced in 2008 (7.0)
