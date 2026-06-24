# Collect | [SpanFromLeft]

> [Collect](https://reference.wolfram.com/language/ref/Collect.html)[*expr*,*x*] — collects together terms involving the same powers of objects matching `*x*`.
> [Collect](https://reference.wolfram.com/language/ref/Collect.html)[*expr*,{*x*_1,*x*_2,…}] — successively collects together terms that involve the same powers of objects matching `*x*_1, then `*x*_2, ….
> [Collect](https://reference.wolfram.com/language/ref/Collect.html)[*expr*,*var*,*h*] — applies `*h*` to the expression that forms the coefficient of each term obtained.

## Details and Options

[Collect](https://reference.wolfram.com/language/ref/Collect.html)[*expr*,*x*] effectively writes `*expr*` as a polynomial in `*x*` or a fractional power of `*x*`.

[Collect](https://reference.wolfram.com/language/ref/Collect.html)[*expr*,*x*,[Simplify](https://reference.wolfram.com/language/ref/Simplify.html)] can be used to simplify each coefficient separately.

[Collect](https://reference.wolfram.com/language/ref/Collect.html) automatically threads over lists in `*expr*`, as well as equations, inequalities and logic functions.

[Collect](https://reference.wolfram.com/language/ref/Collect.html) takes the following options:

| [Modulus](https://reference.wolfram.com/language/ref/Modulus.html) | 0 | modulus to assume for integers |
| --- | --- | --- |
| [Trig](https://reference.wolfram.com/language/ref/Trig.html) | [False](https://reference.wolfram.com/language/ref/False.html) | whether to do trigonometric as well as algebraic transformations |

## Examples

### Basic Examples

Combine like terms involving $x$:

```wolfram
Collect[b x^2+5x+7x^2+9a x+2,x]
(* Output *)
2+(5+9 a) x+(7+b) x^2
```

Collect terms involving $x$:

```wolfram
Collect[a x+b y+c x+d y,x]
(* Output *)
(a+c) x+b y+d y
```

Collect terms involving $y$:

```wolfram
Collect[a x+b y+c x + d y,y]
(* Output *)
a x+c x+(b+d) y
```

Collect each power of $x$:

```wolfram
Collect[(1+a+x)^3,x]
(* Output *)
1+3 a+3 a^2+a^3+(3+6 a+3 a^2) x+(3+3 a) x^2+x^3
```

Simplify each coefficient:

```wolfram
Collect[(1+a+x)^3,x,Simplify]
(* Output *)
(1+a)^3+3 (1+a)^2 x+3 (1+a) x^2+x^3
```

### Scope

#### Basic Uses

A polynomial:

```wolfram
Collect[a x^4+b x^4+2a^2x-3b x+x-7,x]
(* Output *)
-7+(1+2 a^2-3 b) x+(a+b) x^4
```

A Puiseux polynomial:

```wolfram
Collect[a Sqrt[x]+Sqrt[x]+x^(2/3)-c x+3x-2b x^(2/3)+5,x]
(* Output *)
5+(1+a) Sqrt[x]+(1-2 b) x^(2/3)+(3-c) x
```

Collect with respect to $x$ first, then collect with respect to $y$:

```wolfram
Collect[(x y+x z+y z+x+y)^3,{x,y}]
(* Output *)
y^3 (1+3 z+3 z^2+z^3)+x^3 (1+y^3+3 z+3 z^2+z^3+y^2 (3+3 z)+y (3+6 z+3 z^2))+x^2 (y^3 (3+3 z)+y^2 (6+12 z+6 z^2)+y (3+9 z+9 z^2+3 z^3))+x (y^3 (3+6 z+3 z^2)+y^2 (3+9 z+9 z^2+3 z^3))
```

#### Advanced Uses

Collect with respect to a pattern:

```wolfram
D[f[Sqrt[x^2+1]],{x,3}]
(* Output *)
(3 x^3 f^′[Sqrt[1+x^2]])/((1+x^2)^(5/2))-(3 x f^′[Sqrt[1+x^2]])/((1+x^2)^(3/2))-(3 x^3 f^′′[Sqrt[1+x^2]])/((1+x^2)^2)+(3 x f^′′[Sqrt[1+x^2]])/(1+x^2)+(x^3 f^(3)[Sqrt[1+x^2]])/((1+x^2)^(3/2))
```

Collect derivative terms:

```wolfram
Collect[%,Derivative[_][f][_],Together]
(* Output *)
-(3 x f^′[Sqrt[1+x^2]])/((1+x^2)^(5/2))+(3 x f^′′[Sqrt[1+x^2]])/((1+x^2)^2)+(x^3 f^(3)[Sqrt[1+x^2]])/((1+x^2)^(3/2))
```

Factor each coefficient after collection of terms:

```wolfram
Collect[(1+a+x)^4,x,Factor]
(* Output *)
(1+a)^4+4 (1+a)^3 x+6 (1+a)^2 x^2+4 (1+a) x^3+x^4
```

Collect terms over the integers modulo $3$:

```wolfram
Collect[(4+2a+x)^2,x,Modulus->3]
(* Output *)
1+a+a^2+(2+a) x+x^2
```

### Options

#### Modulus

Collect over the integers modulo 2:

```wolfram
Collect[(x+1)^2+(a x+b)^2,x,Modulus->2]
(* Output *)
1+b^2+(1+a^2) x^2
```

### Applications

When a polynomial has many variables, it can be put into many different forms. Here is a polynomial in three variables:

```wolfram
Expand[(y+2x+3x y+4x z+5y z)^3]
(* Output *)
8 x^3+12 x^2 y+36 x^3 y+6 x y^2+36 x^2 y^2+54 x^3 y^2+y^3+9 x y^3+27 x^2 y^3+27 x^3 y^3+48 x^3 z+108 x^2 y z+144 x^3 y z+72 x y^2 z+252 x^2 y^2 z+108 x^3 y^2 z+15 y^3 z+90 x y^3 z+135 x^2 y^3 z+96 x^3 z^2+288 x^2 y z^2+144 x^3 y z^2+270 x y^2 z^2+360 x^2 y^2 z^2+75 y^3 z^2+225 x y^3 z^2+64 x^3 z^3+240 x^2 y z^3+300 x y^2 z^3+125 y^3 z^3
```

[Collect](https://reference.wolfram.com/language/ref/Collect.html) reorganizes the polynomial so that $x$ is the dominant variable:

```wolfram
Collect[%,x]
(* Output *)
y^3+15 y^3 z+75 y^3 z^2+125 y^3 z^3+x^3 (8+36 y+54 y^2+27 y^3+48 z+144 y z+108 y^2 z+96 z^2+144 y z^2+64 z^3)+x^2 (12 y+36 y^2+27 y^3+108 y z+252 y^2 z+135 y^3 z+288 y z^2+360 y^2 z^2+240 y z^3)+x (6 y^2+9 y^3+72 y^2 z+90 y^3 z+270 y^2 z^2+225 y^3 z^2+300 y^2 z^3)
```

If $y$ is specified as the parameter, the terms are organized with $y$ as the leading variable:

```wolfram
Collect[%,y]
(* Output *)
8 x^3+48 x^3 z+96 x^3 z^2+64 x^3 z^3+y^3 (1+9 x+27 x^2+27 x^3+15 z+90 x z+135 x^2 z+75 z^2+225 x z^2+125 z^3)+y^2 (6 x+36 x^2+54 x^3+72 x z+252 x^2 z+108 x^3 z+270 x z^2+360 x^2 z^2+300 x z^3)+y (12 x^2+36 x^3+108 x^2 z+144 x^3 z+288 x^2 z^2+144 x^3 z^2+240 x^2 z^3)
```

When two variables are specified, the function will collect with respect to $x$ first, then collect with respect to $y$:

```wolfram
Collect[%,{x,y}]
(* Output *)
y^3 (1+15 z+75 z^2+125 z^3)+x^3 (8+27 y^3+48 z+96 z^2+64 z^3+y^2 (54+108 z)+y (36+144 z+144 z^2))+x^2 (y^3 (27+135 z)+y^2 (36+252 z+360 z^2)+y (12+108 z+288 z^2+240 z^3))+x (y^3 (9+90 z+225 z^2)+y^2 (6+72 z+270 z^2+300 z^3))
```

Or collect with respect to $y$ first, then collect with respect to $x$:

```wolfram
Collect[%,{y,x}]
(* Output *)
x^3 (8+48 z+96 z^2+64 z^3)+y^3 (1+27 x^3+15 z+75 z^2+125 z^3+x^2 (27+135 z)+x (9+90 z+225 z^2))+y (x^3 (36+144 z+144 z^2)+x^2 (12+108 z+288 z^2+240 z^3))+y^2 (x^3 (54+108 z)+x^2 (36+252 z+360 z^2)+x (6+72 z+270 z^2+300 z^3))
```

### Properties & Relations

[Expand](https://reference.wolfram.com/language/ref/Expand.html) is effectively the inverse of [Collect](https://reference.wolfram.com/language/ref/Collect.html):

```wolfram
f=Expand[(x+y+1)^5];
```

```wolfram
Collect[f,x]
(* Output *)
1+x^5+5 y+10 y^2+10 y^3+5 y^4+y^5+x^4 (5+5 y)+x^3 (10+20 y+10 y^2)+x^2 (10+30 y+30 y^2+10 y^3)+x (5+20 y+30 y^2+20 y^3+5 y^4)
```

```wolfram
Expand[%]===f
(* Output *)
True
```

The order of variables matters:

```wolfram
Collect[(1+x+y+z)^2,{x,y},Simplify]
(* Output *)
x^2+y^2+2 y (1+z)+(1+z)^2+x (2 y+2 (1+z))
```

```wolfram
Collect[(1+x+y+z)^2,{y,x},Simplify]
(* Output *)
x^2+y^2+2 x (1+z)+(1+z)^2+y (2 x+2 (1+z))
```

Use [NonCommutativeCollect](https://reference.wolfram.com/language/ref/NonCommutativeCollect.html) to collect terms in a noncommutative polynomial:

```wolfram
NonCommutativeCollect[x**x**y**x+2x**x**z**x+3x**y**x**x+4x**z**x**x,x]
(* Output *)
x**(3 y+4 z)**NonCommutativeMultiply+NonCommutativeMultiply**(y+2 z)**x
```

## Tech Notes ▪Putting Expressions into Different Forms ▪Structural Operations on Polynomials

## Related Guides ▪Algebraic Transformations ▪Formula Manipulation ▪Polynomial Algebra

## Related Links [NKS|Online](http://www.wolframscience.com/nks/search/?q=Collect) ([A New Kind of Science](http://www.wolframscience.com/nks/))

## History Introduced in 1988 (1.0) | Updated in 1996 (3.0) ▪ 2007 (6.0)
