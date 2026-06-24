# ToRadicals | [SpanFromLeft]

> [ToRadicals](https://reference.wolfram.com/language/ref/ToRadicals.html)[*expr*] — attempts to express all [Root](https://reference.wolfram.com/language/ref/Root.html) objects in `*expr*` in terms of radicals.

## Details and Options

[ToRadicals](https://reference.wolfram.com/language/ref/ToRadicals.html) can always give expressions in terms of radicals when the highest degree of the polynomial that appears in any [Root](https://reference.wolfram.com/language/ref/Root.html) object is four.

There are some cases in which expressions involving radicals can in principle be given, but [ToRadicals](https://reference.wolfram.com/language/ref/ToRadicals.html) cannot find them.

If [Root](https://reference.wolfram.com/language/ref/Root.html) objects in `*expr*` contain parameters, [ToRadicals](https://reference.wolfram.com/language/ref/ToRadicals.html)[*expr*] may yield a result that is not equal to `*expr*` for all values of the parameters.

[ToRadicals](https://reference.wolfram.com/language/ref/ToRadicals.html) automatically threads over lists, as well as equations, inequalities, and logic functions.

[ToRadicals](https://reference.wolfram.com/language/ref/ToRadicals.html) takes the following options:

| [Assumptions](https://reference.wolfram.com/language/ref/Assumptions.html) | [Automatic](https://reference.wolfram.com/language/ref/Automatic.html) | assumptions on parameters |
| --- | --- | --- |
| [Cubics](https://reference.wolfram.com/language/ref/Cubics.html) | [True](https://reference.wolfram.com/language/ref/True.html) | whether to use explicit radicals to solve all cubics |
| [Quartics](https://reference.wolfram.com/language/ref/Quartics.html) | [True](https://reference.wolfram.com/language/ref/True.html) | whether to use explicit radicals to solve all quartics |

## Examples

### Basic Examples

```wolfram
ToRadicals[Root[#^3+#+11&,1]+Root[#^5-2&,3]]
(* Output *)
(-1)^(4/5) 2^(1/5)-((2)/(3 (-99+Sqrt[9813])))^(1/3)+(((1)/(2) (-99+Sqrt[9813]))^(1/3))/(3^(2/3))
```

### Scope

All cubic [Root](https://reference.wolfram.com/language/ref/Root.html) objects can be converted into radicals:

```wolfram
ToRadicals[Root[#^3-5#^2-7#+9&,1]]
(* Output *)
(5)/(3)-(23^(2/3) (1+ⅈ Sqrt[3]))/(3 (7+3 ⅈ Sqrt[15])^(1/3))-(1)/(6) (1-ⅈ Sqrt[3]) (23 (7+3 ⅈ Sqrt[15]))^(1/3)
```

All quartic [Root](https://reference.wolfram.com/language/ref/Root.html) objects can be converted into radicals:

```wolfram
ToRadicals[Root[#^4+3#^3-5#^2-7#+9&,1]]
(* Output *)
-(3)/(4)-(1)/(4 Sqrt[(3)/(67+4 ((5555)/(2)-(3 Sqrt[82209])/(2))^(1/3)+2 2^(2/3) (5555+3 Sqrt[82209])^(1/3))])-(1)/(2) Sqrt((67)/(6)-(1)/(3) ((5555)/(2)-(3 Sqrt[82209])/(2))^(1/3)-(1)/(3) ((1)/(2) (5555+3 Sqrt[82209]))^(1/3)+(31)/(2) Sqrt[(3)/(67+4 ((5555)/(2)-(3 Sqrt[82209])/(2))^(1/3)+2 2^(2/3) (5555+3 Sqrt[82209])^(1/3))])
```

Some higher[Hyphen]degree [Root](https://reference.wolfram.com/language/ref/Root.html) objects can be represented in terms of radicals:

```wolfram
ToRadicals[Root[#^8-4#^7+14#^6-28#^5+82#^4-122#^3+221#^2-164#+399&,1]]
(* Output *)
(1)/(2) (1-Sqrt[-7+2 ⅈ Sqrt[2 (33+Sqrt[85])]])
```

```wolfram
Root[Cyclotomic[102,x],x,1]
(* Output *)
Root
```

```wolfram
ToRadicals[%]
(* Output *)
-(-1)^(2/51)
```

[ToRadicals](https://reference.wolfram.com/language/ref/ToRadicals.html) also works with [AlgebraicNumber](https://reference.wolfram.com/language/ref/AlgebraicNumber.html) objects:

```wolfram
ToRadicals[AlgebraicNumber[Root[#^3-11#+3&,1],{1,2,3}]]
(* Output *)
1+2 (-((1-ⅈ Sqrt[3]) ((1)/(2) (-27+ⅈ Sqrt[15243]))^(1/3))/(2 3^(2/3))-(11 (1+ⅈ Sqrt[3]))/(2^(2/3) (3 (-27+ⅈ Sqrt[15243]))^(1/3)))+3 (-((1-ⅈ Sqrt[3]) ((1)/(2) (-27+ⅈ Sqrt[15243]))^(1/3))/(2 3^(2/3))-(11 (1+ⅈ Sqrt[3]))/(2^(2/3) (3 (-27+ⅈ Sqrt[15243]))^(1/3)))^2
```

### Generalizations & Extensions

[ToRadicals](https://reference.wolfram.com/language/ref/ToRadicals.html) converts trigonometric functions of rational multiples of $\pi$:

```wolfram
ToRadicals[Sin[7Pi/16]]
(* Output *)
(1)/(2) Sqrt[2+Sqrt[2+Sqrt[2]]]
```

```wolfram
ToRadicals[Tan[Pi/19]]
(* Output *)
-(ⅈ (-1+(-1)^(2/19)))/(1+(-1)^(2/19))
```

### Options

#### Assumptions

The setting of [Assumptions](https://reference.wolfram.com/language/ref/Assumptions.html) affects conversion of parametric [Root](https://reference.wolfram.com/language/ref/Root.html) objects:

```wolfram
rt=Root[(#-a)(#^2-a)&,1]
(* Output *)
Root[a^2-a #1-a #1^2+#1^3&,1]
```

With the default [Automatic](https://reference.wolfram.com/language/ref/Automatic.html) setting, the result may not be equivalent to the [Root](https://reference.wolfram.com/language/ref/Root.html) object:

```wolfram
ToRadicals[rt]
(* Output *)
Sqrt[a]
```

When an assumption is given, the [Root](https://reference.wolfram.com/language/ref/Root.html) object is converted only if an equivalent result is found:

```wolfram
ToRadicals[rt,Assumptions->a>=0]
(* Output *)
-Sqrt[a]
```

```wolfram
ToRadicals[rt,Assumptions->a<0]
(* Output *)
a
```

```wolfram
ToRadicals[rt,Assumptions->True]
(* Output *)
Root[a^2-a #1-a #1^2+#1^3&,1]
```

With [Assumptions](https://reference.wolfram.com/language/ref/Assumptions.html)->[None](https://reference.wolfram.com/language/ref/None.html), parametric [Root](https://reference.wolfram.com/language/ref/Root.html) objects are not converted:

```wolfram
ToRadicals[Root[#^3-a&,1]+Root[#^3-2&,1],Assumptions->None]
(* Output *)
2^(1/3)+Root[-a+#1^3&,1]
```

#### Cubics

With [Cubics](https://reference.wolfram.com/language/ref/Cubics.html)->[False](https://reference.wolfram.com/language/ref/False.html) the general formulas for solving cubic equations are not used:

```wolfram
ToRadicals[Root[#^3-5#^2-7#+9&,1],Cubics->False]
(* Output *)
Root
```

Converting some cubic [Root](https://reference.wolfram.com/language/ref/Root.html) objects does not require the general formulas:

```wolfram
ToRadicals[Root[#^3-3&,1],Cubics->False]
(* Output *)
3^(1/3)
```

#### Quartics

With [Quartics](https://reference.wolfram.com/language/ref/Quartics.html)->[False](https://reference.wolfram.com/language/ref/False.html) the general formulas for solving quartic equations are not used:

```wolfram
ToRadicals[Root[#^4+3#^3-5#^2-7#+9&,1],Quartics->False]
(* Output *)
Root
```

Converting some quartic [Root](https://reference.wolfram.com/language/ref/Root.html) objects does not require the general formulas:

```wolfram
ToRadicals[Root[#^4-7#^2+3&,1],Quartics->False]
(* Output *)
-Sqrt[(7)/(2)+(Sqrt[37])/(2)]
```

### Properties & Relations

[RootReduce](https://reference.wolfram.com/language/ref/RootReduce.html)[ToRadicals[*r*]]==*r* for any algebraic number `*r*` given as a [Root](https://reference.wolfram.com/language/ref/Root.html) object:

```wolfram
r=Root[#^4+3#^3-5#^2-7#+9&,1];
```

```wolfram
RootReduce[ToRadicals[r]]==r
(* Output *)
True
```

By default [Reduce](https://reference.wolfram.com/language/ref/Reduce.html) will not produce radical solutions for general cubics:

```wolfram
Reduce[x^3+x+17==0,x]
(* Output *)
x==Root||x==Root||x==Root
```

Use [ToRadicals](https://reference.wolfram.com/language/ref/ToRadicals.html) to convert:

```wolfram
ToRadicals[%]
(* Output *)
x==-((2)/(3 (-153+Sqrt[23421])))^(1/3)+(((1)/(2) (-153+Sqrt[23421]))^(1/3))/(3^(2/3))||x==-((1+ⅈ Sqrt[3]) ((1)/(2) (-153+Sqrt[23421]))^(1/3))/(2 3^(2/3))+(1-ⅈ Sqrt[3])/(2^(2/3) (3 (-153+Sqrt[23421]))^(1/3))||x==-((1-ⅈ Sqrt[3]) ((1)/(2) (-153+Sqrt[23421]))^(1/3))/(2 3^(2/3))+(1+ⅈ Sqrt[3])/(2^(2/3) (3 (-153+Sqrt[23421]))^(1/3))
```

Alternatively set [Cubics](https://reference.wolfram.com/language/ref/Cubics.html)->[True](https://reference.wolfram.com/language/ref/True.html):

```wolfram
Reduce[x^3+x+17==0,x,Cubics->True]
(* Output *)
x==-((2)/(3 (-153+Sqrt[23421])))^(1/3)+(((1)/(2) (-153+Sqrt[23421]))^(1/3))/(3^(2/3))||x==-((1+ⅈ Sqrt[3]) ((1)/(2) (-153+Sqrt[23421]))^(1/3))/(2 3^(2/3))+(1-ⅈ Sqrt[3])/(2^(2/3) (3 (-153+Sqrt[23421]))^(1/3))||x==-((1-ⅈ Sqrt[3]) ((1)/(2) (-153+Sqrt[23421]))^(1/3))/(2 3^(2/3))+(1+ⅈ Sqrt[3])/(2^(2/3) (3 (-153+Sqrt[23421]))^(1/3))
```

### Possible Issues

In this case [ToRadicals](https://reference.wolfram.com/language/ref/ToRadicals.html) succeeds on the unreduced expression:

```wolfram
e=Sqrt[2]+Root[#^3+11#+3&,1];
```

```wolfram
ToRadicals[{e,RootReduce[e]}]
(* Output *)
{Sqrt[2]-11 ((2)/(3 (-27+Sqrt[16701])))^(1/3)+(((1)/(2) (-27+Sqrt[16701]))^(1/3))/(3^(2/3)),Root}
```

In this case [ToRadicals](https://reference.wolfram.com/language/ref/ToRadicals.html) succeeds on the reduced expression:

```wolfram
e=RootReduce[Sqrt[2]+Root[#^3+11#+3&,1]]-Sqrt[2];
```

```wolfram
ToRadicals[{e,RootReduce[e]}]
(* Output *)
{-Sqrt[2]+Root,-11 ((2)/(3 (-27+Sqrt[16701])))^(1/3)+(((1)/(2) (-27+Sqrt[16701]))^(1/3))/(3^(2/3))}
```

[ToRadicals](https://reference.wolfram.com/language/ref/ToRadicals.html) converts [Root](https://reference.wolfram.com/language/ref/Root.html) objects containing parameters:

```wolfram
ToRadicals[Root[#^3-a&,1]]
(* Output *)
(-1)^(2/3) a^(1/3)
```

The result may not be equal to the [Root](https://reference.wolfram.com/language/ref/Root.html) object for some values of the parameter:

```wolfram
{Root[#^3-a&,1],%}/.{{a->-1},{a->1}}
(* Output *)
{{-1,-1},{1,(-1)^(2/3)}}
```

## Tech Notes ▪Algebraic Numbers

## Related Guides ▪Algebraic Numbers ▪Algebraic Transformations ▪Algebraic Number Theory ▪Polynomial Algebra ▪Number Recognition ▪Polynomial Equations

## History Introduced in 1996 (3.0) | Updated in 2003 (5.0) ▪ 2007 (6.0) ▪ 2010 (8.0)
