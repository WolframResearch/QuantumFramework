# Variables | [SpanFromLeft]

> [Variables](https://reference.wolfram.com/language/ref/Variables.html)[*poly*] — gives a list of all independent variables in a polynomial.

## Examples

### Basic Examples

Find a list of variables of a polynomial:

```wolfram
Variables[(x+y)^2+3z^2-y z+7]
(* Output *)
{x,y,z}
```

### Scope

A polynomial:

```wolfram
Variables[x^4+y^4+z^4-3x y z]
(* Output *)
{x,y,z}
```

A list of polynomials:

```wolfram
Variables[{x^2-a y,y^2-b z,z^2-c x}]
(* Output *)
{a,b,c,x,y,z}
```

A rational function:

```wolfram
Variables[(a-b)/(x+y)-2/z]
(* Output *)
{a,b,x,y,z}
```

Find variables in a radical expression:

```wolfram
Variables[Sqrt[x+y-z^2]+(-2 t)^(2/3)]
(* Output *)
{t,x,y,z}
```

### Options

#### Modulus

Find variables present after reducing coefficients modulo 2:

```wolfram
Variables[x+2y+3z,Modulus->2]
(* Output *)
{x,z}
```

### Properties & Relations

Use [CoefficientList](https://reference.wolfram.com/language/ref/CoefficientList.html) to find coefficients of polynomials:

```wolfram
f=x^2-2 x y+3 y^2+4 x+5y+6;
```

```wolfram
CoefficientList[f,Variables[f]]
(* Output *)
{{6,5,3},{4,-2,0},{1,0,0}}
```

Use [NonCommutativeVariables](https://reference.wolfram.com/language/ref/NonCommutativeVariables.html) to find variables in a noncommutative polynomial:

```wolfram
NonCommutativeVariables[(x+2y)**(3w**z+4z**z)]
(* Output *)
{w,x,y,z}
```

### Possible Issues

[Variables](https://reference.wolfram.com/language/ref/Variables.html) looks for variables only inside sums, products, and rational powers:

```wolfram
Variables[Sin[x]+Cos[x]]
(* Output *)
{Cos[x],Sin[x]}
```

```wolfram
Variables[E^x]
(* Output *)
{}
```

## Tech Notes ▪Finding the Structure of a Polynomial

## Related Guides ▪Polynomial Algebra

## History Introduced in 1988 (1.0)
