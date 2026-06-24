# Cancel | [SpanFromLeft]

> [Cancel](https://reference.wolfram.com/language/ref/Cancel.html)[*expr*] — cancels out common factors in the numerator and denominator of `*expr*`.

## Details and Options

[Cancel](https://reference.wolfram.com/language/ref/Cancel.html) cancels out the greatest common divisor of the numerator and denominator.

[Cancel](https://reference.wolfram.com/language/ref/Cancel.html) takes the following options:

| [Extension](https://reference.wolfram.com/language/ref/Extension.html) | [None](https://reference.wolfram.com/language/ref/None.html) | coefficient field to be used |
| --- | --- | --- |
| [Modulus](https://reference.wolfram.com/language/ref/Modulus.html) | 0 | modulus to assume for integers |
| [Trig](https://reference.wolfram.com/language/ref/Trig.html) | [False](https://reference.wolfram.com/language/ref/False.html) | whether to do trigonometric as well as algebraic transformations  |

[Cancel](https://reference.wolfram.com/language/ref/Cancel.html)[*expr*,[Extension](https://reference.wolfram.com/language/ref/Extension.html)->[Automatic](https://reference.wolfram.com/language/ref/Automatic.html)] allows operations to be performed on algebraic numbers in `*expr*`.

[Cancel](https://reference.wolfram.com/language/ref/Cancel.html)[*expr*,[Trig](https://reference.wolfram.com/language/ref/Trig.html)->[True](https://reference.wolfram.com/language/ref/True.html)] treats trigonometric functions as rational functions of exponentials and manipulates them accordingly.

[Cancel](https://reference.wolfram.com/language/ref/Cancel.html) automatically threads over lists, as well as equations, inequalities and logic functions.

## Examples

### Basic Examples

Cancel common factors:

```wolfram
(x^2-1)/(x-1)
(* Output *)
(-1+x^2)/(-1+x)
```

```wolfram
Cancel[%]
(* Output *)
1+x
```

### Scope

[Cancel](https://reference.wolfram.com/language/ref/Cancel.html) threads over sums:

```wolfram
Cancel[(x-y)/(x^2-y^2)+(x^3-27)/(x^2-9)+(x^3+1)/(x^2-x+1)]
(* Output *)
1+x+(9+3 x+x^2)/(3+x)+(1)/(x+y)
```

[Cancel](https://reference.wolfram.com/language/ref/Cancel.html) threads over Boolean combinations of equations and inequalities:

```wolfram
Cancel[(x-a)/(x^2-a^2)==0&&(x^2-2x+1)/(x-1)>=0]
(* Output *)
(1)/(a+x)==0&&-1+x>=0
```

Compute over a finite field:

```wolfram
ℱ=FiniteField[17,3];
```

```wolfram
Cancel[(ℱ[1]x^2+ℱ[246]x+ℱ[4436])/(ℱ[3]x^2+ℱ[1771])]
(* Output *)
![image](img/image_001.png)
```

### Options

#### Extension

By default, [Cancel](https://reference.wolfram.com/language/ref/Cancel.html) treats algebraic numbers as independent variables:

```wolfram
r=(x^2-2)/(Sqrt[2]+x);
```

```wolfram
Cancel[r]
(* Output *)
(-2+x^2)/(Sqrt[2]+x)
```

With [Extension](https://reference.wolfram.com/language/ref/Extension.html)->[Automatic](https://reference.wolfram.com/language/ref/Automatic.html), [Cancel](https://reference.wolfram.com/language/ref/Cancel.html) recognizes algebraically dependent coefficients:

```wolfram
Cancel[r,Extension->Automatic]
(* Output *)
-Sqrt[2]+x
```

Cancel common factors over a finite field:

```wolfram
ℱ=FiniteField[2,3];
```

```wolfram
Cancel[(x^4+1)/(x+1),Extension->ℱ]
(* Output *)
![image](img/image_003.png)
```

#### Modulus

Over the rational numbers the numerator and the denominator have no common factors:

```wolfram
r=(1+x^2)/(1+x);
```

```wolfram
Cancel[r]
(* Output *)
(1+x^2)/(1+x)
```

Over the integers modulo 2, the denominator divides the numerator:

```wolfram
Cancel[r,Modulus->2]
(* Output *)
1+x
```

#### Trig

By default, [Cancel](https://reference.wolfram.com/language/ref/Cancel.html) treats trigonometric functions as independent variables:

```wolfram
Cancel[Sin[2x]/Sin[x]]
(* Output *)
Csc[x] Sin[2 x]
```

With [Trig](https://reference.wolfram.com/language/ref/Trig.html)->[True](https://reference.wolfram.com/language/ref/True.html), [Cancel](https://reference.wolfram.com/language/ref/Cancel.html) recognizes dependencies between trigonometric functions:

```wolfram
Cancel[Sin[2x]/Sin[x],Trig->True]
(* Output *)
2 Cos[x]
```

### Properties & Relations

[Cancel](https://reference.wolfram.com/language/ref/Cancel.html) cancels common factors between numerators and denominators:

```wolfram
r=(x-1)/(x^2-1)+(x-2)/(x^2-4);
```

```wolfram
Cancel[r]
(* Output *)
(1)/(1+x)+(1)/(2+x)
```

[Together](https://reference.wolfram.com/language/ref/Together.html) combines terms over a common denominator and cancels common factors:

```wolfram
Together[r]
(* Output *)
(3+2 x)/((1+x) (2+x))
```

## Tech Notes ▪Putting Expressions into Different Forms ▪Structural Operations on Rational Expressions

## Related Guides ▪Algebraic Transformations ▪Rational Functions ▪Polynomial Division ▪Formula Manipulation ▪Finite Fields

## History Introduced in 1988 (1.0) | Updated in 1996 (3.0) ▪ 2007 (6.0) ▪ 2023 (13.3)
