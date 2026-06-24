# CoefficientRules | [SpanFromLeft]

> [CoefficientRules](https://reference.wolfram.com/language/ref/CoefficientRules.html)[*poly*,{*x*_1,*x*_2,…}]  — gives the list `{{*e*_11,*e*_12,…}->*c*_1,{*e*_21,…}->*c*_2,…}` of exponent vectors and coefficients for the monomials in `*poly*` with respect to the `*x*_*i*`.
> [CoefficientRules](https://reference.wolfram.com/language/ref/CoefficientRules.html)[*poly*,{*x*_1,*x*_2,…},*order*] — gives the result with the monomial ordering specified by `*order*`.

## Details and Options

[CoefficientRules](https://reference.wolfram.com/language/ref/CoefficientRules.html) works whether or not `*poly*` is explicitly given in expanded form.

[CoefficientRules](https://reference.wolfram.com/language/ref/CoefficientRules.html)[*poly*] is equivalent to [CoefficientRules](https://reference.wolfram.com/language/ref/CoefficientRules.html)[*poly*,[Variables](https://reference.wolfram.com/language/ref/Variables.html)[*poly*]].

Possible settings for `*order*` are the same as in [MonomialList](https://reference.wolfram.com/language/ref/MonomialList.html).

The default order is `"Lexicographic"`.

[CoefficientRules](https://reference.wolfram.com/language/ref/CoefficientRules.html)[*poly*,*vars*,[Modulus](https://reference.wolfram.com/language/ref/Modulus.html) ->*m*] computes the coefficients modulo `*m*`.

[CoefficientRules](https://reference.wolfram.com/language/ref/CoefficientRules.html)[*poly*,[All](https://reference.wolfram.com/language/ref/All.html),*order*] is the same as [CoefficientRules](https://reference.wolfram.com/language/ref/CoefficientRules.html)[*poly*,[Variables](https://reference.wolfram.com/language/ref/Variables.html)[*poly*],*order*].

## Examples

### Basic Examples

Get exponents and coefficients of monomials:

```wolfram
CoefficientRules[(x+y)^3]
(* Output *)
{{3,0}->1,{2,1}->3,{1,2}->3,{0,3}->1}
```

### Scope

Use `"DegreeReverseLexicographic"` monomial ordering:

```wolfram
CoefficientRules[a x y^2+b x^2 z,{x,y,z},"DegreeReverseLexicographic"]
(* Output *)
{{1,2,0}->a,{2,0,1}->b}
```

Specify the same ordering using weight matrix:

```wolfram
CoefficientRules[a x y^2+b x^2 z,{x,y,z},{{1,1,1},{0,0,-1},{0,-1,0}}]
(* Output *)
{{1,2,0}->a,{2,0,1}->b}
```

### Options

#### Modulus

Reduce the coefficients modulo 2:

```wolfram
CoefficientRules[(x+1)^5,x,Modulus->2]
(* Output *)
{{5}->1,{4}->1,{1}->1,{0}->1}
```

### Properties & Relations

[FromCoefficientRules](https://reference.wolfram.com/language/ref/FromCoefficientRules.html) reconstructs the original polynomial:

```wolfram
CoefficientRules[a x^2+b x y+c y^2,{x,y}]
(* Output *)
{{2,0}->a,{1,1}->b,{0,2}->c}
```

```wolfram
FromCoefficientRules[%,{x,y}]
(* Output *)
a x^2+b x y+c y^2
```

[MonomialList](https://reference.wolfram.com/language/ref/MonomialList.html) gives a different representation:

```wolfram
MonomialList[a x^2+b x y+c y^2,{x,y}]
(* Output *)
{a x^2,b x y,c y^2}
```

For two variables `"DegreeLexicographic"` and `"DegreeReverseLexicographic"` coincide:

```wolfram
poly=Sum[RandomInteger[{-10,10}]x^RandomInteger[10]y^RandomInteger[10],{20}]
(* Output *)
3+10 x+2 x^2-3 x^5 y+9 x^10 y-6 x^4 y^2-2 x y^4-10 y^5+10 x y^5-2 x^9 y^5+5 y^6-4 x^3 y^6-3 x^5 y^6+9 x y^7+3 x^5 y^7-x^6 y^7-x^10 y^8-9 x y^9+3 y^10
```

```wolfram
CoefficientRules[poly,{x,y},"DegreeLexicographic"]===CoefficientRules[poly,{x,y},"DegreeReverseLexicographic"]
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
CoefficientRules[y+x z]
(* Output *)
{{1,0,0}->1,{0,1,1}->1}
```

```wolfram
FromCoefficientRules[%,{y,x,z}]
(* Output *)
y+x z
```

### Neat Examples

Visualize monomial orderings in 2D:

```wolfram
MonomialOrderPlot[poly_,{x1_,x2_},ord_:"Lexicographic",o:OptionsPattern[]]:=
Graphics[{Blue,Arrow/@Partition[CoefficientRules[poly,{x1,x2},ord][[All,1]],2,1]},Sequence@@FilterRules[{o},Options[Graphics]]]
```

The standard built-in orderings:

```wolfram
orders={"Lexicographic","DegreeLexicographic","DegreeReverseLexicographic","NegativeLexicographic","NegativeDegreeLexicographic","NegativeDegreeReverseLexicographic"};
```

In 2D some orderings cannot be distinguished:

```wolfram
Table[MonomialOrderPlot[Sum[x^i y^j,{i,0,5},{j,0,5}],{x,y},o,PlotLabel->Style[o,Small]],{o,orders}]
(* Output *)
{[Graphics],[Graphics],[Graphics],[Graphics],[Graphics],[Graphics]}
```

Visualize monomial orderings in 3D:

```wolfram
MonomialOrderPlot3D[poly_,{x1_,x2_,x3_},ord_:"Lexicographic",o:OptionsPattern[]]:=
Graphics3D[{Arrowheads[.06],Arrow[Tube[#,.02]]&/@Partition[CoefficientRules[poly,{x1,x2,x3},ord][[All,1]],2,1]},Sequence@@FilterRules[{o},Options[Graphics3D]]]
```

```wolfram
orders={"Lexicographic","DegreeLexicographic","DegreeReverseLexicographic","NegativeLexicographic","NegativeDegreeLexicographic","NegativeDegreeReverseLexicographic"};
```

In 3D all orderings are distinct:

```wolfram
Table[MonomialOrderPlot3D[Sum[x^i y^j z^k,{i,0,2},{j,0,2},{k,0,2}],{x,y,z},o,PlotLabel->Style[o,Small]],{o,orders}]
(* Output *)
{[Graphics3D],[Graphics3D],[Graphics3D],[Graphics3D],[Graphics3D],[Graphics3D]}
```

## Tech Notes ▪Polynomial Orderings

## Related Guides ▪Polynomial Algebra

## History Introduced in 2008 (7.0)
