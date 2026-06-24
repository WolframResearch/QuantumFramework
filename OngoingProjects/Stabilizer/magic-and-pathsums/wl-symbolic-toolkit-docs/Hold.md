# Hold | [SpanFromLeft]

> [Hold](https://reference.wolfram.com/language/ref/Hold.html)[*expr*] — maintains `*expr*` in an unevaluated form.

## Details

[Hold](https://reference.wolfram.com/language/ref/Hold.html) allows you to use an expression that has not undergone normal evaluation.

[Hold](https://reference.wolfram.com/language/ref/Hold.html) has attribute [HoldAll](https://reference.wolfram.com/language/ref/HoldAll.html) and performs no operation on its arguments.

[Hold](https://reference.wolfram.com/language/ref/Hold.html) is removed by [ReleaseHold](https://reference.wolfram.com/language/ref/ReleaseHold.html).

[Hold](https://reference.wolfram.com/language/ref/Hold.html)[*e*_1,*e*_2,…] maintains a sequence of unevaluated expressions to which a function can be applied using [Apply](https://reference.wolfram.com/language/ref/Apply.html).

[Hold](https://reference.wolfram.com/language/ref/Hold.html) can be inserted as a wrapper by such functions as [ToExpression](https://reference.wolfram.com/language/ref/ToExpression.html) and [Extract](https://reference.wolfram.com/language/ref/Extract.html).

Even though `*expr*` itself is not evaluated, [Hold](https://reference.wolfram.com/language/ref/Hold.html)[*expr*] may still evaluate if `*expr*` is of the form `*f*[*args*]`, and upvalues for `*f*` have been defined.

## Examples

### Basic Examples

Hold an expression to prevent evaluation:

```wolfram
Hold[2+2]
(* Output *)
Hold[2+2]
```

Release the hold:

```wolfram
ReleaseHold[%]
(* Output *)
4
```

### Scope

Extract a part without allowing it to evaluate:

```wolfram
Extract[Hold[{1,1+1,π}], {1,2},Hold]
(* Output *)
Hold[1+1]
```

### Applications

Find the length of each expression in a held list without evaluation:

```wolfram
list=Hold[1+2,2 3 4 5,1/0,Quit[]];
```

```wolfram
Apply[List,Map[Hold,list]]
(* Output *)
{Hold[1+2],Hold[2 3 4 5],Hold[(1)/(0)],Hold[Quit[]]}
```

```wolfram
%/.Hold[e_]:>Length[Unevaluated[e]]
(* Output *)
{2,4,2,0}
```

Evaluate every sum (only) inside a held expression:

```wolfram
expr=Hold[{1+2,g[3+4,2 3],f[1+g[2+3]]}]
(* Output *)
Hold[{1+2,g[3+4,2 3],f[1+g[2+3]]}]
```

```wolfram
pos=Position[expr,_Plus]
(* Output *)
{{1,1},{1,2,1},{1,3,1,2,1},{1,3,1}}
```

```wolfram
val=Extract[expr,pos]
(* Output *)
{3,7,5,1+g[5]}
```

```wolfram
ReplacePart[expr,Thread[pos->val]]
(* Output *)
Hold[{3,g[7,2 3],f[1+g[5]]}]
```

### Properties & Relations

[Hold](https://reference.wolfram.com/language/ref/Hold.html) is a container with the attribute [HoldAll](https://reference.wolfram.com/language/ref/HoldAll.html):

```wolfram
Attributes[Hold]
(* Output *)
{HoldAll,Protected}
```

```wolfram
Hold[1+2]
(* Output *)
Hold[1+2]
```

Use [Evaluate](https://reference.wolfram.com/language/ref/Evaluate.html) to force evaluation of an argument of [Hold](https://reference.wolfram.com/language/ref/Hold.html):

```wolfram
Hold[Evaluate[1+1],2+2]
(* Output *)
Hold[2,2+2]
```

[Unevaluated](https://reference.wolfram.com/language/ref/Unevaluated.html) inside a held expression is not removed:

```wolfram
Hold[Unevaluated[1+1]]
(* Output *)
Hold[Unevaluated[1+1]]
```

Sequence splicing still happens inside [Hold](https://reference.wolfram.com/language/ref/Hold.html):

```wolfram
Hold[Sequence[1+1,2+2]]
(* Output *)
Hold[1+1,2+2]
```

Use the container [HoldComplete](https://reference.wolfram.com/language/ref/HoldComplete.html) to suppress even such transformations:

```wolfram
HoldComplete[Sequence[1+1,2+2]]
(* Output *)
HoldComplete[Sequence[1+1,2+2]]
```

Upvalues work inside [Hold](https://reference.wolfram.com/language/ref/Hold.html):

```wolfram
h/:Hold[h[x_]]:=f[x]
```

```wolfram
Hold[h[1+2]]
(* Output *)
f[3]
```

They do not work inside [HoldComplete](https://reference.wolfram.com/language/ref/HoldComplete.html):

```wolfram
h/:HoldComplete[h[x_]]:=f[x]
```

```wolfram
HoldComplete[h[1+2]]
(* Output *)
HoldComplete[h[1+2]]
```

Substitution works inside [Hold](https://reference.wolfram.com/language/ref/Hold.html):

```wolfram
Hold[f[1+2]]/.f[x_]:>g[x]
(* Output *)
Hold[g[1+2]]
```

Insert into a held expression:

```wolfram
Insert[Hold[x+x],y,{1,2}]
(* Output *)
Hold[x+y+x]
```

[ReleaseHold](https://reference.wolfram.com/language/ref/ReleaseHold.html) removes one level of [Hold](https://reference.wolfram.com/language/ref/Hold.html):

```wolfram
ReleaseHold[Hold[1+2]]
(* Output *)
3
```

[HoldForm](https://reference.wolfram.com/language/ref/HoldForm.html) is like [Hold](https://reference.wolfram.com/language/ref/Hold.html) but is normally not shown in the output:

```wolfram
HoldForm[1+2]
(* Output *)
1+2
```

```wolfram
FullForm[%]
(* Output *)
HoldForm[Plus[1,2]]
```

[Hold](https://reference.wolfram.com/language/ref/Hold.html) can be used to freeze the result of [ToExpression](https://reference.wolfram.com/language/ref/ToExpression.html) before it is evaluated:

```wolfram
ToExpression["1+1",InputForm,Hold]
(* Output *)
Hold[1+1]
```

## Tech Notes ▪Evaluation ▪Non[Hyphen]Standard Evaluation

## Related Guides ▪Evaluation Control ▪Expressions

[Graphics] | Related Workflows
▪Handle Code Symbolically
▪Substitute Values of Variables in Functions That Hold Their Arguments

## Related Links [An Elementary Introduction to the Wolfram Language](https://www.wolfram.com/language/elementary-introduction/39-immediate-and-delayed-values.html)[: Immediate and Delayed Values](https://www.wolfram.com/language/elementary-introduction/39-immediate-and-delayed-values.html)

## History Introduced in 1988 (1.0)
