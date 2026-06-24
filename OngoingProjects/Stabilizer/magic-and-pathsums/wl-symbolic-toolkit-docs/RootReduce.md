# RootReduce | [SpanFromLeft]

> [RootReduce](https://reference.wolfram.com/language/ref/RootReduce.html)[*expr*] — attempts to reduce `*expr*` to a single [Root](https://reference.wolfram.com/language/ref/Root.html) object.

## Details and Options

If `*expr*` consists only of integers and [Root](https://reference.wolfram.com/language/ref/Root.html) and [AlgebraicNumber](https://reference.wolfram.com/language/ref/AlgebraicNumber.html) objects combined using algebraic operations, then the result from [RootReduce](https://reference.wolfram.com/language/ref/RootReduce.html)[*expr*] will always be a single [Root](https://reference.wolfram.com/language/ref/Root.html) object.

Simple [Root](https://reference.wolfram.com/language/ref/Root.html) objects may in turn automatically evaluate to rational expressions or combinations of radicals.

[RootReduce](https://reference.wolfram.com/language/ref/RootReduce.html) automatically threads over lists, as well as equations, inequalities, and logic functions.

## Examples

### Basic Examples

Reduce to a single [Root](https://reference.wolfram.com/language/ref/Root.html) object:

```wolfram
RootReduce[Sqrt[2]+Sqrt[3]]
(* Output *)
Root
```

### Scope

Combinations of radical expressions:

```wolfram
RootReduce[2^(1/3)+Sqrt[5]]
(* Output *)
Root
```

Combinations of [Root](https://reference.wolfram.com/language/ref/Root.html) objects:

```wolfram
RootReduce[Root[#^5+11#+1&,1]Root[#^3+#+17&,1]]
(* Output *)
Root
```

Reduce any algebraic combination of radicals, [Root](https://reference.wolfram.com/language/ref/Root.html), and [AlgebraicNumber](https://reference.wolfram.com/language/ref/AlgebraicNumber.html) objects:

```wolfram
RootReduce[Sqrt[Root[#^3-7#-2&,1]+1]/AlgebraicNumber[3^(1/3),{1,2,3}]]
(* Output *)
Root
```

The result is always a [Root](https://reference.wolfram.com/language/ref/Root.html) object, a quadratic radical expression, or a rational number:

```wolfram
RootReduce[Sqrt[7]/Sqrt[5+2Sqrt[6]]]
(* Output *)
Root
```

```wolfram
RootReduce[(Sqrt[2]+Sqrt[3]+Sqrt[6]+3)/Sqrt[5+2Sqrt[6]]]
(* Output *)
1+Sqrt[3]
```

```wolfram
RootReduce[(Sqrt[18]+Sqrt[27])/Sqrt[5+2Sqrt[6]]]
(* Output *)
3
```

### Options

#### Method

By default, [RootReduce](https://reference.wolfram.com/language/ref/RootReduce.html) heuristically selects the method to use:

```wolfram
a=2 2^(1/3)+3 3^(1/7)+5 2^(1/3)3^(1/7);
RootReduce[a]//Timing
(* Output *)
{0.1248008,Root}
```

In this case conversion to [AlgebraicNumber](https://reference.wolfram.com/language/ref/AlgebraicNumber.html) objects in a common number field is used:

```wolfram
RootReduce[a,Method->"NumberField"]//Timing
(* Output *)
{0.1092007,Root}
```

The other available method recursively performs arithmetic operations:

```wolfram
RootReduce[a,Method->"Recursive"]//Timing
(* Output *)
{3.432022,Root}
```

Here the `"Recursive"` method is faster:

```wolfram
b=RootReduce[2^(1/3)+3^(1/3)+1];
c=b-2^(1/3)-3^(1/3);
RootReduce[c]//Timing
(* Output *)
{0.0700000000000181,1}
```

```wolfram
RootReduce[c,Method->"NumberField"]//Timing
(* Output *)
{0.7609999999999846,1}
```

```wolfram
RootReduce[c,Method->"Recursive"]//Timing
(* Output *)
{0.050000000000001044,1}
```

### Applications

The numeric test used by [Equal](https://reference.wolfram.com/language/ref/Equal.html) cannot prove the equality:

```wolfram
Sqrt[2]+Sqrt[3]+Sqrt[5]==Sqrt[10+2Sqrt[15]+4Sqrt[4+Sqrt[15]]]
(* Output *)
Sqrt[2]+Sqrt[3]+Sqrt[5]==Sqrt[10+2 Sqrt[15]+4 Sqrt[4+Sqrt[15]]]
```

[RootReduce](https://reference.wolfram.com/language/ref/RootReduce.html) proves that the two algebraic numbers are equal:

```wolfram
RootReduce[%]
(* Output *)
True
```

[FullSimplify](https://reference.wolfram.com/language/ref/FullSimplify.html) will use [RootReduce](https://reference.wolfram.com/language/ref/RootReduce.html):

```wolfram
FullSimplify[%%]
(* Output *)
True
```

### Properties & Relations

The results given by [RootReduce](https://reference.wolfram.com/language/ref/RootReduce.html) are canonical:

```wolfram
algs={Sqrt[2]+Sqrt[3]+Sqrt[5],Sqrt[10+2Sqrt[15]+4Sqrt[4+Sqrt[15]]],Sqrt[25+4Sqrt[15]+2Sqrt[46+8Sqrt[15]]]-Sqrt[5]}
(* Output *)
{Sqrt[2]+Sqrt[3]+Sqrt[5],Sqrt[10+2 Sqrt[15]+4 Sqrt[4+Sqrt[15]]],-Sqrt[5]+Sqrt[25+4 Sqrt[15]+2 Sqrt[46+8 Sqrt[15]]]}
```

```wolfram
RootReduce/@algs
(* Output *)
{Root,Root,Root}
```

In general the degree of the reduced polynomial will be the product of the degrees:

```wolfram
a=RootReduce[Root[#^3+#+11&,1]+Root[#^3+#+17&,1]]
(* Output *)
Root
```

```wolfram
MinimalPolynomial[a,x]
(* Output *)
21952-2352 x-2693 x^3+84 x^4+9 x^5+84 x^6+6 x^7+x^9
```

In exceptional cases the result can have a lower degree:

```wolfram
b=RootReduce[a-Root[#^3+#+11&,1]]
(* Output *)
Root
```

```wolfram
MinimalPolynomial[b,x]
(* Output *)
17+x+x^3
```

[Root](https://reference.wolfram.com/language/ref/Root.html) objects can be converted to [AlgebraicNumber](https://reference.wolfram.com/language/ref/AlgebraicNumber.html) objects:

```wolfram
ToNumberField[Root[#^5+#+11&,1]]
(* Output *)
AlgebraicNumber[Root,{0,1,0,0,0}]
```

[RootReduce](https://reference.wolfram.com/language/ref/RootReduce.html) converts from [AlgebraicNumber](https://reference.wolfram.com/language/ref/AlgebraicNumber.html) objects:

```wolfram
RootReduce[%]
(* Output *)
Root
```

## Tech Notes ▪Algebraic Numbers

## Related Guides ▪Algebraic Numbers ▪Algebraic Transformations ▪Algebraic Number Theory ▪Theorem Proving ▪Polynomial Algebra ▪Formula Manipulation ▪Number Recognition

## History Introduced in 1996 (3.0) | Updated in 2007 (6.0)
