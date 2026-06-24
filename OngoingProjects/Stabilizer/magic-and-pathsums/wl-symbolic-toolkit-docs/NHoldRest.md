# NHoldRest | [SpanFromLeft]

> [NHoldRest](https://reference.wolfram.com/language/ref/NHoldRest.html) — is an attribute which specifies that all but the first argument to a function should not be affected by [N](https://reference.wolfram.com/language/ref/N.html).

## Examples

### Basic Examples

Prevent [N](https://reference.wolfram.com/language/ref/N.html) from affecting all but the first argument of a function:

```wolfram
SetAttributes[f,NHoldRest]
```

```wolfram
N[f[Pi,3]]
(* Output *)
f[3.141592653589793,3]
```

```wolfram
N[f[E,Pi],20]
(* Output *)
f[2.71828182845904523536028747135266249776,π]
```

### Scope

System symbols with the [NHoldRest](https://reference.wolfram.com/language/ref/NHoldRest.html) attribute:

```wolfram
ssymb=Cases[Map[ToExpression,Names["System`*"]],_Symbol];
Select[ssymb,MemberQ[Attributes[#],NHoldRest]&]
(* Output *)
{Drop,Extract,GraphicsComplex,HeldPart,Overscript,Part,Subscript,Subsuperscript,Superscript,Take,Underoverscript,Underscript}
```

Use [Part](https://reference.wolfram.com/language/ref/Part.html) symbolically:

```wolfram
e = Quiet[x[[1]]+x[[-1]] - 1]
(* Output *)
-1+x[[-1]]+x[[1]]
```

Using [N](https://reference.wolfram.com/language/ref/N.html) does not affect the specified parts:

```wolfram
ne = Quiet[N[e]]
(* Output *)
-1.+x[[-1]]+x[[1]]
```

The expression works when `x` is substituted with a list:

```wolfram
ne /. x->{2,4,8}
(* Output *)
9.
```

### Applications

Define a function that represents the real-valued root $x^{1/n}$ for positive $x$ and positive integer $n$:

```wolfram
nthroot[x_?InexactNumberQ /; Positive[x],n_Integer/;Positive[n]]:=FixedPoint[((n-1) #^n + x)/(n #^(n-1))&,N[1,Precision[x]]]
```

Prevent its second argument from being converted to real using [NHoldRest](https://reference.wolfram.com/language/ref/NHoldRest.html):

```wolfram
SetAttributes[nthroot, NHoldRest]
```

An exact representation of the cube root of 2:

```wolfram
crt=nthroot[2,3]
(* Output *)
nthroot[2,3]
```

Machine-number approximation:

```wolfram
N[crt]
(* Output *)
1.259921049894873
```

47-digit approximation:

```wolfram
N[crt,47]
(* Output *)
1.2599210498948731647672106072782283505702514647015079800819
```

### Properties & Relations

[Subscript](https://reference.wolfram.com/language/ref/Subscript.html) by default has the [NHoldRest](https://reference.wolfram.com/language/ref/NHoldRest.html) attribute:

```wolfram
Attributes[Subscript]
(* Output *)
{NHoldRest}
```

This means that subscripts generally do not change under [N](https://reference.wolfram.com/language/ref/N.html):

```wolfram
poly=Sum[(i+1)a_i x^i,{i,0,5}]
(* Output *)
a_0+2 x a_1+3 x^2 a_2+4 x^3 a_3+5 x^4 a_4+6 x^5 a_5
```

```wolfram
N[poly]
(* Output *)
a_0+2. x a_1+3. x^2 a_2+4. x^3 a_3+5. x^4 a_4+6. x^5 a_5
```

[HoldRest](https://reference.wolfram.com/language/ref/HoldRest.html) prevents evaluation while [NHoldRest](https://reference.wolfram.com/language/ref/NHoldRest.html) only prevents numerical evaluation:

```wolfram
SetAttributes[f1,HoldRest]
```

```wolfram
f1[1,2+3]
(* Output *)
f1[1,2+3]
```

```wolfram
N[%]
(* Output *)
f1[1.,2.+3.]
```

```wolfram
SetAttributes[f2,NHoldRest]
```

```wolfram
f2[1,2+3]
(* Output *)
f2[1,5]
```

```wolfram
N[%]
(* Output *)
f2[1.,5]
```

You can prevent both by setting both attributes:

```wolfram
SetAttributes[f,{HoldRest,NHoldRest}]
```

```wolfram
f[1,2+3]
(* Output *)
f[1,2+3]
```

```wolfram
N[%]
(* Output *)
f[1.,2+3]
```

## Tech Notes ▪Attributes ▪Controlling Numerical Evaluation

## Related Guides ▪Evaluation Control ▪Attributes

## History Introduced in 1996 (3.0)
