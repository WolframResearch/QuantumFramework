# CoefficientList | [SpanFromLeft]

> [CoefficientList](https://reference.wolfram.com/language/ref/CoefficientList.html)[*poly*,*var*] — gives a list of coefficients of powers of `*var*` in `*poly*`, starting with power `0
> [CoefficientList](https://reference.wolfram.com/language/ref/CoefficientList.html)[*poly*,{*var*_1,*var*_2,…}] — gives an array of coefficients of the `*var*_*i*`.
> [CoefficientList](https://reference.wolfram.com/language/ref/CoefficientList.html)[*poly*,{*var*_1,*var*_2,…},{*dim*_1,*dim*_2,…}] — gives an array of dimensions `{*dim*_1,*dim*_2,…}`, truncating or padding with zeros as needed.

## Details and Options

The dimensions of the array returned by [CoefficientList](https://reference.wolfram.com/language/ref/CoefficientList.html) are determined by the values of the [Exponent](https://reference.wolfram.com/language/ref/Exponent.html)[*poly*,*var*_*i*].

Terms that do not contain positive integer powers of a particular variable are included in the first element of the list for that variable.

[CoefficientList](https://reference.wolfram.com/language/ref/CoefficientList.html) always returns a full rectangular array. Combinations of powers that do not appear in `*poly*` give zeros in the array.

[CoefficientList](https://reference.wolfram.com/language/ref/CoefficientList.html)[0,*var*] gives `{}`.

[CoefficientList](https://reference.wolfram.com/language/ref/CoefficientList.html) works whether or not `*poly*` is explicitly given in expanded form.

## Examples

### Basic Examples

Find the coefficients in a polynomial:

```wolfram
CoefficientList[1+6x-x^4,x]
(* Output *)
{1,6,0,0,-1}
```

[CoefficientList](https://reference.wolfram.com/language/ref/CoefficientList.html) works even when the polynomial has not been expanded out:

```wolfram
CoefficientList[(1+x)^10 ,x]
(* Output *)
{1,10,45,120,210,252,210,120,45,10,1}
```

Matrix of coefficients for a quadratic function:

```wolfram
CoefficientList[1+a x^2+b x y + c y^2,{x,y}]
(* Output *)
{{1,0,c},{0,b,0},{a,0,0}}
```

### Scope

Univariate polynomial coefficient lists:

```wolfram
CoefficientList[(2x+3)^5,x]
(* Output *)
{243,810,1080,720,240,32}
```

```wolfram
CoefficientList[a x^4+b x^3+c x^2+d x+e,x]
(* Output *)
{e,d,c,b,a}
```

Multivariate polynomial coefficient lists:

```wolfram
CoefficientList[(3x+4 y+1)^3,{x,y}]
(* Output *)
{{1,12,48,64},{9,72,144,0},{27,108,0,0},{27,0,0,0}}
```

```wolfram
CoefficientList[(x+y+z+1)(2x+3y^2+4z^3+5),{x,y,z}]
(* Output *)
{{{5,5,0,4,4},{5,0,0,4,0},{3,3,0,0,0},{3,0,0,0,0}},{{7,2,0,4,0},{2,0,0,0,0},{3,0,0,0,0},{0,0,0,0,0}},{{2,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0}}}
```

### Options

#### Modulus

Coefficient list over the integers modulo 2:

```wolfram
CoefficientList[(x+1)^5,x,Modulus->2]
(* Output *)
{1,1,0,0,1,1}
```

### Properties & Relations

Use [Coefficient](https://reference.wolfram.com/language/ref/Coefficient.html) to get a coefficient at a specified power of the variable:

```wolfram
f=(2x+3)^3;
```

```wolfram
Coefficient[f,x,2]
(* Output *)
36
```

The list of coefficients can be obtained using [Coefficient](https://reference.wolfram.com/language/ref/Coefficient.html) and [Exponent](https://reference.wolfram.com/language/ref/Exponent.html):

```wolfram
Coefficient[f,x,#]&/@Range[0,Exponent[f,x]]
(* Output *)
{27,54,36,8}
```

```wolfram
CoefficientList[f,x]
(* Output *)
{27,54,36,8}
```

[FromDigits](https://reference.wolfram.com/language/ref/FromDigits.html) can reconstruct a univariate polynomial from the list of its coefficients:

```wolfram
CoefficientList[a +b x+c x^2,x]
(* Output *)
{a,b,c}
```

```wolfram
FromDigits[Reverse[%],x]
(* Output *)
a+b x+c x^2
```

Fold the operation for multivariate polynomials:

```wolfram
CoefficientList[(x+2y)^3,{x,y}]
(* Output *)
{{0,0,0,8},{0,0,12,0},{0,6,0,0},{1,0,0,0}}
```

```wolfram
Fold[FromDigits[Reverse[#1],#2]&,%,{x,y}]
(* Output *)
x^3+6 x^2 y+y^2 (12 x+8 y)
```

```wolfram
Expand[(x+2y)^3-%]
(* Output *)
0
```

Polynomial multiplication is convolution as performed by [ListConvolve](https://reference.wolfram.com/language/ref/ListConvolve.html):

```wolfram
CoefficientList[(a+b x+c x^2)(1+2x+3x^2+4x^3),x]
(* Output *)
{a,2 a+b,3 a+2 b+c,4 a+3 b+2 c,4 b+3 c,4 c}
```

```wolfram
ListConvolve[{a,b,c},{1,2,3,4},{1,-1},0]
(* Output *)
{a,2 a+b,3 a+2 b+c,4 a+3 b+2 c,4 b+3 c,4 c}
```

For multivariate polynomials, [CoefficientList](https://reference.wolfram.com/language/ref/CoefficientList.html) gives a tensor of the coefficients:

```wolfram
f=(3x+5y)^4;
```

```wolfram
cl=CoefficientList[f,{x, y}]
(* Output *)
{{0,0,0,0,625},{0,0,0,1500,0},{0,0,1350,0,0},{0,540,0,0,0},{81,0,0,0,0}}
```

[CoefficientArrays](https://reference.wolfram.com/language/ref/CoefficientArrays.html) gives the list of arrays of polynomial coefficients ordered by total degrees:

```wolfram
ca=CoefficientArrays[f,{x,y}]
(* Output *)
{0,SparseArray[...],SparseArray[...],SparseArray[...],SparseArray[...]}
```

The coefficient of $x y^{3}$:

```wolfram
Coefficient[f,x y^3]
(* Output *)
1500
```

In `cl`, the coefficient of `x^*a* y^*b*` is the element at position `*{**a*+1,*b*+1}`:

```wolfram
cl[[1+1,1+3]]
(* Output *)
1500
```

In `ca`, the position of this coefficient is `*a*+*b*+1 followed by `*a*` `1s and `*b*` `2s (`1 and `2 indicate the first and second variables):

```wolfram
ca[[5,1,2,2,2]]
(* Output *)
1500
```

## Tech Notes ▪Finding the Structure of a Polynomial

## Related Guides ▪Polynomial Systems ▪Series Expansions ▪Constructing Lists ▪Polynomial Algebra

## Related Links [NKS|Online](http://www.wolframscience.com/nks/search/?q=CoefficientList) ([A New Kind of Science](http://www.wolframscience.com/nks/))

## History Introduced in 1988 (1.0) | Updated in 1996 (3.0) ▪ 2003 (5.0) ▪ 2015 (10.3)
