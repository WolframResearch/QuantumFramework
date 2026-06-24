# HoldAllComplete | [SpanFromLeft]

> [HoldAllComplete](https://reference.wolfram.com/language/ref/HoldAllComplete.html) — is an attribute which specifies that all arguments to a function are not to be modified or looked at in any way in the process of evaluation.

## Details

By setting the attribute [HoldAllComplete](https://reference.wolfram.com/language/ref/HoldAllComplete.html), you can effectively shield the arguments of a function from all aspects of the standard Wolfram Language evaluation process.

[HoldAllComplete](https://reference.wolfram.com/language/ref/HoldAllComplete.html) not only prevents arguments from being evaluated, but also prevents [Sequence](https://reference.wolfram.com/language/ref/Sequence.html) objects from being flattened, [Unevaluated](https://reference.wolfram.com/language/ref/Unevaluated.html) wrappers from being stripped, and upvalues associated with arguments from being used.

[Evaluate](https://reference.wolfram.com/language/ref/Evaluate.html) cannot be used to override [HoldAllComplete](https://reference.wolfram.com/language/ref/HoldAllComplete.html).

## Examples

### Basic Examples

```wolfram
Attributes[HoldComplete]
(* Output *)
{HoldAllComplete,Protected}
```

```wolfram
HoldComplete[1+1,Evaluate[1+2],Sequence[3,4]]
(* Output *)
HoldComplete[1+1,Evaluate[1+2],Sequence[3,4]]
```

### Properties & Relations

[HoldComplete](https://reference.wolfram.com/language/ref/HoldComplete.html) is a standard container with attribute [HoldAllComplete](https://reference.wolfram.com/language/ref/HoldAllComplete.html):

```wolfram
Attributes[HoldComplete]
(* Output *)
{HoldAllComplete,Protected}
```

```wolfram
HoldComplete[Sequence[a,b],1+2]
(* Output *)
HoldComplete[Sequence[a,b],1+2]
```

No form of evaluation control affects an expression with attribute [HoldAllComplete](https://reference.wolfram.com/language/ref/HoldAllComplete.html):

```wolfram
HoldComplete[1+2]
(* Output *)
HoldComplete[1+2]
```

```wolfram
HoldComplete[Evaluate[1+2]]
(* Output *)
HoldComplete[Evaluate[1+2]]
```

```wolfram
HoldComplete[Sequence[a,b]]
(* Output *)
HoldComplete[Sequence[a,b]]
```

```wolfram
g/:HoldComplete[g[x_]]:=x
```

```wolfram
HoldComplete[g[1]]
(* Output *)
HoldComplete[g[1]]
```

Substitution still happens inside an expression with attribute [HoldAllComplete](https://reference.wolfram.com/language/ref/HoldAllComplete.html):

```wolfram
HoldComplete[f[1+2]]/.f[x_]:>g[x]
(* Output *)
HoldComplete[g[1+2]]
```

[Unevaluated](https://reference.wolfram.com/language/ref/Unevaluated.html) has the attribute [HoldAllComplete](https://reference.wolfram.com/language/ref/HoldAllComplete.html):

```wolfram
Length[Unevaluated[Sequence[a,b]]]
(* Output *)
2
```

### Possible Issues

[HoldAllComplete](https://reference.wolfram.com/language/ref/HoldAllComplete.html) affects only evaluation; input transformations are still applied:

```wolfram
FullForm[HoldComplete[a-b,a/b]]
(* Output *)
HoldComplete[Plus[a,Times[-1,b]],Times[a,Power[b,-1]]]
```

[HoldAllComplete](https://reference.wolfram.com/language/ref/HoldAllComplete.html) does not prevent formatting:

```wolfram
HoldComplete[Grid[{{1,2},{3,4}}]]
(* Output *)
HoldComplete[{{1, 2}, {3, 4}}]
```

Add [DisableFormatting](https://reference.wolfram.com/language/ref/DisableFormatting.html) to prevent formatting:

```wolfram
HoldComplete[DisableFormatting[Grid[{{1,2},{3,4}}]]]
(* Output *)
HoldComplete[Grid[{{1, 2}, {3, 4}}]]
```

### Neat Examples

A fast way to compute the Hofstadter-Conway sequence [[more info]](http://mathworld.wolfram.com/Hofstadter-Conway10000-DollarSequence.html):

```wolfram
SetAttributes[h,HoldAllComplete];
hc[m_Integer/;m>=1]:=
Module[{a=h@@Table[0,{m}],i},
a[[1]]=1;
Do[a[[i]]=a[[a[[i-1]]]]+a[[i-a[[i-1]]]],{i,2,m}];
a[[m]]
]
```

```wolfram
Timing[hc[5 10^5]]
(* Output *)
{5.047,500000}
```

## Tech Notes ▪Attributes ▪Non[Hyphen]Standard Evaluation

## Related Guides ▪Evaluation Control ▪Attributes

[Graphics] | Related Workflows
▪Substitute Values of Variables in Functions That Hold Their Arguments

## History Introduced in 1996 (3.0)
