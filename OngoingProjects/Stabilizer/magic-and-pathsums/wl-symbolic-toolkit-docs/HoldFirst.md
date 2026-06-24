# HoldFirst | [SpanFromLeft]

> [HoldFirst](https://reference.wolfram.com/language/ref/HoldFirst.html) — is an attribute that specifies that the first argument to a function is to be maintained in an unevaluated form.

## Examples

### Basic Examples

Give the function `f` the attribute [HoldAll](https://reference.wolfram.com/language/ref/HoldAll.html):

```wolfram
SetAttributes[f,HoldFirst]
```

As a result, the first argument is not evaluated prior to entering the (in this case nonexistent) function body:

```wolfram
f[1+1,2+2,3+3]
(* Output *)
f[1+1,4,6]
```

Compare with a function `g` that has no attributes:

```wolfram
g[1+1,2+2,3+3]
(* Output *)
g[2,4,6]
```

### Scope

[SetAttributes](https://reference.wolfram.com/language/ref/SetAttributes.html) does not evaluate the symbols in its first argument:

```wolfram
f=c;
SetAttributes[f,Constant]
```

As a result, the preceding input set the attribute on `f`, not `c`:

```wolfram
Attributes[{f,c}]
(* Output *)
{{Constant},{}}
```

This is made possible by the fact that [SetAttributes](https://reference.wolfram.com/language/ref/SetAttributes.html) has the attribute [HoldFirst](https://reference.wolfram.com/language/ref/HoldFirst.html):

```wolfram
Attributes[SetAttributes]
(* Output *)
{HoldFirst,Protected}
```

Give the function `inout` the attribute [HoldFirst](https://reference.wolfram.com/language/ref/HoldFirst.html):

```wolfram
SetAttributes[inout, HoldFirst]
```

Define a function of two arguments that returns them wrapped in [Hold](https://reference.wolfram.com/language/ref/Hold.html) as well as the first without wrapping:

```wolfram
inout[in_,tag_]:={Hold[in], in,Hold[tag]}
```

The second argument is evaluated prior to entering the function body, while the first is evaluated inside the function body:

```wolfram
inout[1-1,1+1]
(* Output *)
{Hold[1-1],0,Hold[2]}
```

### Applications

Definitions for unevaluated expressions can implement call-by-name semantics:

```wolfram
SetAttributes[f,HoldFirst]
```

```wolfram
f[sym_,val_]:=(sym=val^2)
```

```wolfram
x=17;
```

```wolfram
f[x,5]
(* Output *)
25
```

The global variable has been modified:

```wolfram
x
(* Output *)
25
```

### Properties & Relations

Functions that operate on symbols often need the [HoldFirst](https://reference.wolfram.com/language/ref/HoldFirst.html) attribute:

```wolfram
Attributes[SetAttributes]
(* Output *)
{HoldFirst,Protected}
```

```wolfram
i=5;
SetAttributes[i,Constant]
```

Assignments do not evaluate their left-hand sides:

```wolfram
Attributes[Set]
(* Output *)
{HoldFirst,Protected,SequenceHold}
```

```wolfram
i=5;
```

```wolfram
i=6;
```

Use [Evaluate](https://reference.wolfram.com/language/ref/Evaluate.html) to force evaluation of an argument of a [HoldFirst](https://reference.wolfram.com/language/ref/HoldFirst.html) function:

```wolfram
SetAttributes[h,HoldFirst]
```

```wolfram
h[Evaluate[1+1],2+2,3+3]
(* Output *)
h[2,4,6]
```

Suppress the evaluation of the first argument of a pure function:

```wolfram
Function[{e,f},Hold[e,f],{HoldFirst}][1+2,2+3]
(* Output *)
Hold[1+2,5]
```

Sequence splitting still happens for [HoldFirst](https://reference.wolfram.com/language/ref/HoldFirst.html) functions:

```wolfram
SetAttributes[f,HoldFirst]
```

```wolfram
f[Sequence[1+1,2+2],3+3]
(* Output *)
f[1+1,4,6]
```

[NHoldFirst](https://reference.wolfram.com/language/ref/NHoldFirst.html) protects arguments from [N](https://reference.wolfram.com/language/ref/N.html) but evaluates them normally otherwise:

```wolfram
Attributes[EllipticTheta]
(* Output *)
{Listable,NHoldFirst,Protected}
```

```wolfram
N[EllipticTheta[1,a,b]]
(* Output *)
EllipticTheta[1,a,b]
```

```wolfram
EllipticTheta[1+1,a,b]
(* Output *)
EllipticTheta[2,a,b]
```

## Tech Notes ▪Attributes ▪Non[Hyphen]Standard Evaluation

## Related Guides ▪Evaluation Control ▪Attributes ▪Expressions ▪Defining Variables and Functions

[Graphics] | Related Workflows
▪Substitute Values of Variables in Functions That Hold Their Arguments

## History Introduced in 1988 (1.0)
