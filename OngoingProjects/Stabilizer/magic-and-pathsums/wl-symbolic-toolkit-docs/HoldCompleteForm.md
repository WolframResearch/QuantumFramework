# HoldCompleteForm

> [HoldCompleteForm](https://reference.wolfram.com/language/ref/HoldCompleteForm.html)[*expr*]  — prints as the expression `*expr*`, shielding `*expr*` completely from the standard Wolfram Language evaluation process.

## Details

[HoldCompleteForm](https://reference.wolfram.com/language/ref/HoldCompleteForm.html) allows you to see the output form of an expression without performing any evaluation of the expression.

[HoldCompleteForm](https://reference.wolfram.com/language/ref/HoldCompleteForm.html) has attribute [HoldAllComplete](https://reference.wolfram.com/language/ref/HoldAllComplete.html) and performs no operation on its arguments.

[HoldCompleteForm](https://reference.wolfram.com/language/ref/HoldCompleteForm.html) is removed by [ReleaseHold](https://reference.wolfram.com/language/ref/ReleaseHold.html).

[HoldCompleteForm](https://reference.wolfram.com/language/ref/HoldCompleteForm.html) can be inserted as a wrapper by such functions as [ToExpression](https://reference.wolfram.com/language/ref/ToExpression.html) and [Extract](https://reference.wolfram.com/language/ref/Extract.html).

Unlike [HoldForm](https://reference.wolfram.com/language/ref/HoldForm.html), [HoldCompleteForm](https://reference.wolfram.com/language/ref/HoldCompleteForm.html)[*expr*] remains unevaluated even if `*expr*` is of the form `*f*[*args*]` and upvalues for `*f*` have been defined.

## Examples

### Basic Examples

Addition in held form:

```wolfram
HoldCompleteForm[1+1]
(* Output *)
1+1
```

An unevaluated sequence:

```wolfram
HoldCompleteForm[Sequence[1+2,3+4]]
(* Output *)
Sequence[1+2,3+4]
```

Evaluate the expression by applying [ReleaseHold](https://reference.wolfram.com/language/ref/ReleaseHold.html):

```wolfram
ReleaseHold[%]
(* Output *)
Sequence[3,7]
```

### Scope

Display a sum of squares in unevaluated form:

```wolfram
∑_{i=1}^{10}(HoldCompleteForm[#1^2]&)[i]
(* Output *)
1^2+2^2+3^2+4^2+5^2+6^2+7^2+8^2+9^2+10^2
```

View the unevaluated form of an extracted part:

```wolfram
Extract[HoldComplete[{1,1+1,π}], {1,2},HoldCompleteForm]
(* Output *)
1+1
```

### Applications

Issue a message, inserting values with no evaluation:

```wolfram
Message[f::invl,HoldCompleteForm[Sequence[1,2,3]]]
(* Output *)
f
```

### Properties & Relations

[HoldCompleteForm](https://reference.wolfram.com/language/ref/HoldCompleteForm.html) displays without a wrapper:

```wolfram
HoldCompleteForm[1+1]
(* Output *)
1+1
```

See the [FullForm](https://reference.wolfram.com/language/ref/FullForm.html) expression:

```wolfram
FullForm[%]
(* Output *)
HoldCompleteForm[Plus[1,1]]
```

[HoldComplete](https://reference.wolfram.com/language/ref/HoldComplete.html) displays the held expression with a wrapper:

```wolfram
HoldComplete[1+1]
(* Output *)
HoldComplete[1+1]
```

[Evaluate](https://reference.wolfram.com/language/ref/Evaluate.html) cannot force evaluation of an argument of [HoldCompleteForm](https://reference.wolfram.com/language/ref/HoldCompleteForm.html):

```wolfram
HoldCompleteForm[Evaluate[1+1]]
(* Output *)
Evaluate[1+1]
```

Use [Evaluate](https://reference.wolfram.com/language/ref/Evaluate.html) to force evaluation of an argument of [HoldForm](https://reference.wolfram.com/language/ref/HoldForm.html):

```wolfram
HoldForm[Evaluate[1+1]]
(* Output *)
2
```

[Sequence](https://reference.wolfram.com/language/ref/Sequence.html) and [Splice](https://reference.wolfram.com/language/ref/Splice.html) splicing does not happen inside [HoldCompleteForm](https://reference.wolfram.com/language/ref/HoldCompleteForm.html):

```wolfram
HoldCompleteForm[Sequence[1+1]]
(* Output *)
Sequence[1+1]
```

```wolfram
HoldCompleteForm[Splice[{1+1},_]]
(* Output *)
Splice[{1+1},_]
```

Use [HoldForm](https://reference.wolfram.com/language/ref/HoldForm.html) to allow such transformations:

```wolfram
HoldForm[Sequence[1+1]]
(* Output *)
1+1
```

```wolfram
HoldForm[Splice[{1+1},_]]
(* Output *)
1+1
```

Upvalues do not work inside [HoldCompleteForm](https://reference.wolfram.com/language/ref/HoldCompleteForm.html):

```wolfram
HoldCompleteForm[h[x_]]^:=HoldCompleteForm[f[x]]
```

```wolfram
HoldCompleteForm[h[1+2]]
(* Output *)
h[1+2]
```

They do work inside [HoldForm](https://reference.wolfram.com/language/ref/HoldForm.html):

```wolfram
HoldForm[h[x_]]^:=HoldForm[f[x]]
```

```wolfram
HoldForm[h[1+2]]
(* Output *)
f[1+2]
```

[ReleaseHold](https://reference.wolfram.com/language/ref/ReleaseHold.html) removes one level of [HoldCompleteForm](https://reference.wolfram.com/language/ref/HoldCompleteForm.html):

```wolfram
ReleaseHold[HoldCompleteForm[1+2]]
(* Output *)
3
```

[HoldCompleteForm](https://reference.wolfram.com/language/ref/HoldCompleteForm.html) gives an object that is never evaluated:

```wolfram
HoldCompleteForm[1+2]
(* Output *)
1+2
```

Copy the output and paste it into an input cell. The `1+2 is still unevaluated:

```wolfram
1+2
(* Output *)
1+2
```

[Defer](https://reference.wolfram.com/language/ref/Defer.html) gives an object whose evaluation is merely deferred until it is explicitly given as Wolfram Language input:

```wolfram
Defer[1+2]
(* Output *)
1+2
```

Copy the output and paste it into an input cell. The `1+2 is evaluated:

```wolfram
1+2
(* Output *)
3
```

[HoldCompleteForm](https://reference.wolfram.com/language/ref/HoldCompleteForm.html) maintains expressions in unevaluated form, and all parts are inactive:

```wolfram
HoldCompleteForm[Sin[ArcTan[1]]]
(* Output *)
Sin[ArcTan[1]]
```

[Inactivate](https://reference.wolfram.com/language/ref/Inactivate.html) maintains symbols in inactive form and allows parts of expressions to be inactive:

```wolfram
Inactive[Sin][ArcTan[1]]
(* Output *)
Sin[(π)/(4)]
```

[Hold](https://reference.wolfram.com/language/ref/Hold.html) can be used to freeze the result of [ToExpression](https://reference.wolfram.com/language/ref/ToExpression.html) before it is evaluated:

```wolfram
ToExpression["{Sequence[a,b]}",InputForm,HoldCompleteForm]
(* Output *)
{Sequence[a,b]}
```

### Possible Issues

[HoldPattern](https://reference.wolfram.com/language/ref/HoldPattern.html)[*expr*] is equivalent to `*expr*` for pattern matching but maintains `*expr*` in an unevaluated form:

```wolfram
f[HoldPattern[x_+y_]]:=Hold[x,y]
```

```wolfram
f[Unevaluated[1+2]]
(* Output *)
Hold[1,2]
```

[HoldCompleteForm](https://reference.wolfram.com/language/ref/HoldCompleteForm.html)[*expr*] is not equivalent to `*expr*` for pattern matching:

```wolfram
g[HoldCompleteForm[x_+y_]]:=Hold[x,y]
```

```wolfram
g[Unevaluated[1+2]]
(* Output *)
g[Unevaluated[1+2]]
```

Only a literal [HoldCompleteForm](https://reference.wolfram.com/language/ref/HoldCompleteForm.html) expression is matched:

```wolfram
g[HoldCompleteForm[1+2]]
(* Output *)
Hold[1,2]
```

## Tech Notes ▪Evaluation ▪Non[Hyphen]Standard Evaluation ▪String[Hyphen]Oriented Output Formats

## Related Guides ▪Evaluation Control ▪Mathematical Typesetting

[Graphics] | Related Workflows
▪Handle Code Symbolically
▪Substitute Values of Variables in Functions That Hold Their Arguments

## Related Links [An Elementary Introduction to the Wolfram Language](https://www.wolfram.com/language/elementary-introduction/39-immediate-and-delayed-values.html)[: Immediate and Delayed Values](https://www.wolfram.com/language/elementary-introduction/39-immediate-and-delayed-values.html)

## History Introduced in 2025 (14.2)
