# Evaluate | [SpanFromLeft]

> [Evaluate](https://reference.wolfram.com/language/ref/Evaluate.html)[*expr*] — causes `*expr*` to be evaluated even if it appears as the argument of a function whose attributes specify that it should be held unevaluated.

## Details

You can use [Evaluate](https://reference.wolfram.com/language/ref/Evaluate.html) to override [HoldFirst](https://reference.wolfram.com/language/ref/HoldFirst.html) etc. attributes of built[Hyphen]in functions.

[Evaluate](https://reference.wolfram.com/language/ref/Evaluate.html) only overrides [HoldFirst](https://reference.wolfram.com/language/ref/HoldFirst.html) etc. attributes when it appears directly as the head of the function argument that would otherwise be held.

## Examples

### Basic Examples

Evaluate inside a [Hold](https://reference.wolfram.com/language/ref/Hold.html):

```wolfram
Hold[Evaluate[1+1],2+2]
(* Output *)
Hold[2,2+2]
```

### Scope

[Evaluate](https://reference.wolfram.com/language/ref/Evaluate.html) works for arguments of any symbol with attributes [HoldAll](https://reference.wolfram.com/language/ref/HoldAll.html), [HoldFirst](https://reference.wolfram.com/language/ref/HoldFirst.html), or [HoldRest](https://reference.wolfram.com/language/ref/HoldRest.html):

```wolfram
Attributes[Attributes]
(* Output *)
{HoldAll,Listable,Protected}
```

Since [Attributes](https://reference.wolfram.com/language/ref/Attributes.html) is [HoldAll](https://reference.wolfram.com/language/ref/HoldAll.html), use [Evaluate](https://reference.wolfram.com/language/ref/Evaluate.html) to find the properties of the value of `*x*`:

```wolfram
x=Plus;
```

```wolfram
{Attributes[x],Attributes[Evaluate[x]]}
(* Output *)
{{},{Flat,Listable,NumericFunction,OneIdentity,Orderless,Protected}}
```

### Applications

Unprotect a system symbol to make a definition for it:

```wolfram
protected=Unprotect[Sqrt]
(* Output *)
{"Sqrt"}
```

```wolfram
Sqrt[x_^2]:=x
```

Restore protection:

```wolfram
Protect[Evaluate[protected]]
(* Output *)
{"Sqrt"}
```

Force evaluation of the right-hand side of a delayed definition:

```wolfram
Expand[(1+x)^3]
(* Output *)
1+3 x+3 x^2+x^3
```

```wolfram
f[x_]:=Evaluate[%]
```

```wolfram
Definition[f]
(* Output *)
{{f[x_]:=1+3 x+3 x^2+x^3}}
```

Build a function from an expression:

```wolfram
ch=ChebyshevT[5,x]
(* Output *)
5 x-20 x^3+16 x^5
```

```wolfram
Function[x,Evaluate[ch]]
(* Output *)
Function[x,5 x-20 x^3+16 x^5]
```

```wolfram
%[10]
(* Output *)
1580050
```

### Properties & Relations

[Evaluate](https://reference.wolfram.com/language/ref/Evaluate.html) does not work inside functions with attribute [HoldAllComplete](https://reference.wolfram.com/language/ref/HoldAllComplete.html):

```wolfram
HoldComplete[Evaluate[1+2]]
(* Output *)
HoldComplete[Evaluate[1+2]]
```

Use [Unevaluated](https://reference.wolfram.com/language/ref/Unevaluated.html) to temporarily treat a function as if it were [HoldAll](https://reference.wolfram.com/language/ref/HoldAll.html):

```wolfram
Length[Unevaluated[1+2+3]]
(* Output *)
3
```

```wolfram
Length[1+2+3]
(* Output *)
0
```

[Evaluate](https://reference.wolfram.com/language/ref/Evaluate.html) does not work inside [Unevaluated](https://reference.wolfram.com/language/ref/Unevaluated.html):

```wolfram
Unevaluated[Evaluate[1+1]]
(* Output *)
Unevaluated[Evaluate[1+1]]
```

### Possible Issues

[Evaluate](https://reference.wolfram.com/language/ref/Evaluate.html) works only on the first level, directly inside a held function:

```wolfram
Hold[f[Evaluate[1+2]]]
(* Output *)
Hold[f[Evaluate[1+2]]]
```

## Tech Notes ▪Basic Plotting ▪Non[Hyphen]Standard Evaluation ▪Evaluation in Iteration Functions

## Related Guides ▪Evaluation Control ▪Expressions

## History Introduced in 1991 (2.0)
