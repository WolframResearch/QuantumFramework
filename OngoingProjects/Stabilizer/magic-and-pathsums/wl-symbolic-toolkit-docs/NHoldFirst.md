# NHoldFirst | [SpanFromLeft]

> [NHoldFirst](https://reference.wolfram.com/language/ref/NHoldFirst.html) — is an attribute which specifies that the first argument to a function should not be affected by [N](https://reference.wolfram.com/language/ref/N.html).

## Examples

### Basic Examples

Prevent [N](https://reference.wolfram.com/language/ref/N.html) from affecting the first argument of a function:

```wolfram
SetAttributes[f,NHoldFirst]
```

```wolfram
N[f[3,Pi]]
(* Output *)
f[3,3.141592653589793]
```

```wolfram
N[f[Pi,E],20]
(* Output *)
f[π,2.71828182845904523536028579257075851155]
```

### Scope

System symbols with the [NHoldFirst](https://reference.wolfram.com/language/ref/NHoldFirst.html) attribute:

```wolfram
ssymb=Cases[Map[ToExpression,Names["System`*"]],_Symbol];
Select[ssymb,MemberQ[Attributes[#],NHoldFirst]&]
(* Output *)
{AiryAiZero,AiryBiZero,BellB,EllipticTheta,EllipticThetaPrime,MathieuC,MathieuCharacteristicA,MathieuCharacteristicB,MathieuCPrime,MathieuS,MathieuSPrime,StieltjesGamma,ZetaZero}
```

The $k$$^{th}$ zero of the zeta function on the critical line with the imaginary part greater than $t+\pi$:

```wolfram
zz20=ZetaZero[20,t+Pi]
(* Output *)
ZetaZero[20,π+t]
```

[N](https://reference.wolfram.com/language/ref/N.html) does not affect the index $k$:

```wolfram
N[zz20]
(* Output *)
ZetaZero[20,3.141592653589793+t]
```

When the second argument is numeric, [N](https://reference.wolfram.com/language/ref/N.html) evaluates numerically:

```wolfram
N[zz20 /. t->1]
(* Output *)
0.5+77.1448400688748 ⅈ
```

### Applications

Define an inverse for $J_{k}(x)$:

```wolfram
e:binv[k_Integer,x_?InexactNumberQ]:= With[{r=Quiet[Check[z /. FindRoot[BesselJ[k,z]-x,{z,0},WorkingPrecision->Precision[x]],$Failed]]},r /;NumberQ[r]]
```

Set the [NHoldFirst](https://reference.wolfram.com/language/ref/NHoldFirst.html) attribute so that $k$ remains an integer:

```wolfram
SetAttributes[binv,NHoldFirst]
```

A symbolic representation of an inverse of $J_{1}(x)$:

```wolfram
binv[1,x]
(* Output *)
binv[1,x]
```

This remains unmodified with [N](https://reference.wolfram.com/language/ref/N.html):

```wolfram
N[%]
(* Output *)
binv[1,x]
```

With a numeric value of `x`, the function `binv` evaluates numerically:

```wolfram
N[binv[1,1/Pi]]
(* Output *)
0.674209497633699
```

### Properties & Relations

[HoldFirst](https://reference.wolfram.com/language/ref/HoldFirst.html) prevents evaluation while [NHoldFirst](https://reference.wolfram.com/language/ref/NHoldFirst.html) only prevents numerical evaluation:

```wolfram
SetAttributes[f1,HoldFirst]
```

```wolfram
f1[1+2,3+4]
(* Output *)
f1[1+2,7]
```

```wolfram
N[%]
(* Output *)
f1[1.+2.,7.]
```

```wolfram
SetAttributes[f2,NHoldFirst]
```

```wolfram
f2[1+2,3+4]
(* Output *)
f2[3,7]
```

```wolfram
N[%]
(* Output *)
f2[3,7.]
```

You can prevent both by setting both attributes:

```wolfram
SetAttributes[f,{HoldFirst,NHoldFirst}]
```

```wolfram
f[1+2,3+4]
(* Output *)
f[1+2,7]
```

```wolfram
N[%]
(* Output *)
f[1+2,7.]
```

## Tech Notes ▪Attributes ▪Controlling Numerical Evaluation

## Related Guides ▪Evaluation Control ▪Attributes

## History Introduced in 1996 (3.0)
