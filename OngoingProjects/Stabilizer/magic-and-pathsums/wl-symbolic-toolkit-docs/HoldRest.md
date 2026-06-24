# HoldRest | [SpanFromLeft]

> [HoldRest](https://reference.wolfram.com/language/ref/HoldRest.html) — is an attribute which specifies that all but the first argument to a function are to be maintained in an unevaluated form.

## Examples

### Basic Examples

Give the function `f` the attribute [HoldAll](https://reference.wolfram.com/language/ref/HoldAll.html):

```wolfram
SetAttributes[f,HoldRest]
```

As a result, the second and subsequent arguments are not evaluated prior to entering the (in this case nonexistent) function body:

```wolfram
f[1+1,2+2,3+3]
(* Output *)
f[2,2+2,3+3]
```

Compare with a function `g` that has no attributes:

```wolfram
g[1+1,2+2,3+3]
(* Output *)
g[2,4,6]
```

### Scope

The first argument of [If](https://reference.wolfram.com/language/ref/If.html) is always evaluated, but the others are only evaluated if that corresponding branch is entered:

```wolfram
If[a||b||False,1+1,2+2]
(* Output *)
If[a||b,1+1,2+2]
```

This is because [If](https://reference.wolfram.com/language/ref/If.html) has the attribute [HoldRest](https://reference.wolfram.com/language/ref/HoldRest.html):

```wolfram
Attributes[If]
(* Output *)
{HoldRest,Protected}
```

Give the function `inout` the attribute [HoldRest](https://reference.wolfram.com/language/ref/HoldRest.html):

```wolfram
SetAttributes[inout, HoldRest]
```

Define a function of two arguments that returns them wrapped in [Hold](https://reference.wolfram.com/language/ref/Hold.html) as well as the second without wrapping:

```wolfram
inout[tag_,in_]:={Hold[tag],Hold[in], in}
```

The first argument is evaluated prior to entering the function body, while the second is evaluated inside the function body:

```wolfram
inout[1-1,1+1]
(* Output *)
{Hold[0],Hold[1+1],2}
```

### Applications

Implement your own conditional:

```wolfram
SetAttributes[if,HoldRest]
if[True,then_,else_]:=then;
if[False,then_,else_]:=else;
```

```wolfram
if[1<2,a,1/0]
(* Output *)
a
```

```wolfram
if[a==b,1+1,2+2]
(* Output *)
if[a==b,1+1,2+2]
```

### Properties & Relations

Use [Evaluate](https://reference.wolfram.com/language/ref/Evaluate.html) to force evaluation of an argument of a [HoldRest](https://reference.wolfram.com/language/ref/HoldRest.html) function:

```wolfram
SetAttributes[h,HoldRest]
```

```wolfram
h[1+1,Evaluate[2+2],3+3]
(* Output *)
h[2,4,3+3]
```

Suppress the evaluation of all but the first argument of a pure function:

```wolfram
Function[{e,f,g},Hold[e,f,g],{HoldRest}][1+2,2+3,3+4]
(* Output *)
Hold[3,2+3,3+4]
```

Sequence splitting still happens for [HoldRest](https://reference.wolfram.com/language/ref/HoldRest.html) functions:

```wolfram
SetAttributes[f,HoldRest]
```

```wolfram
f[1+1,Sequence[1+1,2+2]]
(* Output *)
f[2,1+1,2+2]
```

[NHoldRest](https://reference.wolfram.com/language/ref/NHoldRest.html) protects arguments from [N](https://reference.wolfram.com/language/ref/N.html), but evaluates them normally otherwise:

```wolfram
Attributes[Subscript]
(* Output *)
{NHoldRest}
```

```wolfram
N[a_1]
(* Output *)
a_1
```

```wolfram
a_1+2
(* Output *)
a_3
```

## Tech Notes ▪Attributes ▪Non[Hyphen]Standard Evaluation

## Related Guides ▪Evaluation Control ▪Attributes

[Graphics] | Related Workflows
▪Substitute Values of Variables in Functions That Hold Their Arguments

## History Introduced in 1988 (1.0)
