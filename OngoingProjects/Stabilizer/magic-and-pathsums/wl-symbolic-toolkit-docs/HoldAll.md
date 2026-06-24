# HoldAll | [SpanFromLeft]

> [HoldAll](https://reference.wolfram.com/language/ref/HoldAll.html) — is an attribute that specifies that all arguments to a function are to be maintained in an unevaluated form.

## Details

You can use [Evaluate](https://reference.wolfram.com/language/ref/Evaluate.html) to evaluate the arguments of a [HoldAll](https://reference.wolfram.com/language/ref/HoldAll.html) function in a controlled way.

Even when a function has attribute [HoldAll](https://reference.wolfram.com/language/ref/HoldAll.html), [Sequence](https://reference.wolfram.com/language/ref/Sequence.html) objects that appear in its arguments are still by default flattened, [Unevaluated](https://reference.wolfram.com/language/ref/Unevaluated.html) wrappers are stripped, and upvalues associated with the arguments are used.

## Examples

### Basic Examples

Give the function `f` the attribute [HoldAll](https://reference.wolfram.com/language/ref/HoldAll.html):

```wolfram
SetAttributes[f,HoldAll]
```

As a result, arguments are not evaluated prior to entering the (in this case nonexistent) function body:

```wolfram
f[1+1,2+2,3+3]
(* Output *)
f[1+1,2+2,3+3]
```

Compare with a function `g` that has no attributes:

```wolfram
g[1+1,2+2,3+3]
(* Output *)
g[2,4,6]
```

### Scope

This is the full form of the evaluation result:

```wolfram
FullForm[1-1]
(* Output *)
0
```

This is the full form of the input, before evaluation:

```wolfram
HoldForm[FullForm[1-1]]
(* Output *)
Plus[1,-1]
```

This is made possible by the fact that [HoldForm](https://reference.wolfram.com/language/ref/HoldForm.html) has the [HoldAll](https://reference.wolfram.com/language/ref/HoldAll.html) attribute:

```wolfram
Attributes[HoldForm]
(* Output *)
{HoldAll,Protected}
```

Give the function `inout` the attribute [HoldAll](https://reference.wolfram.com/language/ref/HoldAll.html):

```wolfram
SetAttributes[inout, HoldAll]
```

Define a function of two arguments that returns them wrapped in [Hold](https://reference.wolfram.com/language/ref/Hold.html) as well as the first without wrapping:

```wolfram
inout[in_,tag_]:={Hold[in], in,Hold[tag]}
```

Neither argument is evaluated prior to entering the function body, but the first argument is evaluated inside the function:

```wolfram
inout[1-1,1+1]
(* Output *)
{Hold[1-1],0,Hold[1+1]}
```

### Applications

Many functions with scoping behavior have the [HoldAll](https://reference.wolfram.com/language/ref/HoldAll.html) attribute:

```wolfram
Attributes[Plot//Evaluate]
(* Output *)
{HoldAll,Protected,ReadProtected}
```

Plotting lists of functions will use separate styling for the different functions:

```wolfram
Plot[{Sin[x],Cos[x]},{x,0,2π}]
```

*([Graphics])*

If the list structure is not manifest, no separate styling is provided:

```wolfram
Plot[Through[{Sin,Cos}[x]],{x,0, 2π}]
```

*([Graphics])*

Use [Evaluate](https://reference.wolfram.com/language/ref/Evaluate.html) to make the list structure manifest:

```wolfram
Plot[Evaluate[Through[{Sin,Cos}[x]]],{x,0, 2π}]
```

*([Graphics])*

Different vector-valued functions in a list will still get separate styling:

```wolfram
Plot[{Through[{Sin,Cos}[x]],Through[{Tan,Cot}[x]]},{x,0,2π}]
```

*([Graphics])*

Use [HoldAll](https://reference.wolfram.com/language/ref/HoldAll.html) and [Unevaluated](https://reference.wolfram.com/language/ref/Unevaluated.html) to suppress evaluation of symbols wherever it would occur:

```wolfram
SetAttributes[symbolLength,HoldAll];
symbolLength[s_Symbol]:=StringLength[SymbolName[Unevaluated[s]]]
```

Find the length of a symbol's name even if it has a value:

```wolfram
xyzzy=1;
```

```wolfram
symbolLength[xyzzy]
(* Output *)
5
```

Implement your own control structure:

```wolfram
until::usage="until[cond,cmd] repeats evaluating cmd until cond is True.";
```

```wolfram
SetAttributes[until,HoldAll]
```

```wolfram
until[cond_,cmd_]:=While[True,cmd;If[cond,Return[]]]
```

```wolfram
i=1;
until[Prime[i]>10^6,i++];i
(* Output *)
78499
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

Functions that operate on symbols often need the [HoldAll](https://reference.wolfram.com/language/ref/HoldAll.html) attribute:

```wolfram
i=5;
Protect[i]
(* Output *)
{i}
```

Without the attribute, they would operate on the symbol's value:

```wolfram
Protect[Evaluate[i]]
(* Output *)
Protect
(* Output *)
Protect[5]
```

Control structures such as [Table](https://reference.wolfram.com/language/ref/Table.html) protect their arguments from evaluation:

```wolfram
i=5;
```

```wolfram
Table[i,{i,1,4}]
(* Output *)
{1,2,3,4}
```

Otherwise, global values might interfere with their operation:

```wolfram
Table[Evaluate[i],{i,1,4}]
(* Output *)
{5,5,5,5}
```

Use [Evaluate](https://reference.wolfram.com/language/ref/Evaluate.html) to force evaluation of an argument of a [HoldAll](https://reference.wolfram.com/language/ref/HoldAll.html) function:

```wolfram
Attributes[Unprotect]
(* Output *)
{HoldAll,Protected}
```

```wolfram
syms={Plus,Times};
```

```wolfram
Unprotect[Evaluate[syms]]
(* Output *)
{"Plus","Times"}
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

Use [Unevaluated](https://reference.wolfram.com/language/ref/Unevaluated.html) to temporarily treat a function as if it had the attribute [HoldAll](https://reference.wolfram.com/language/ref/HoldAll.html):

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

Suppress the evaluation of the arguments of a pure function:

```wolfram
Function[e,Hold[e],{HoldAll}][1+2]
(* Output *)
Hold[1+2]
```

```wolfram
Function[e,Hold[e]][1+2]
(* Output *)
Hold[3]
```

Sequence splicing still happens for functions with the attribute [HoldAll](https://reference.wolfram.com/language/ref/HoldAll.html):

```wolfram
SetAttributes[f,HoldAll]
f[Sequence[a,b,c]]
(* Output *)
f[a,b,c]
```

[Unevaluated](https://reference.wolfram.com/language/ref/Unevaluated.html) is stripped prior to entering the function body:

```wolfram
f[x_]:=Hold[x]
f[Unevaluated[1+1]]
(* Output *)
Hold[1+1]
```

```wolfram
_[___,g,___]^:=upvalue
f[Sequence[a,b,c],g]
(* Output *)
upvalue
```

The attribute [HoldAllComplete](https://reference.wolfram.com/language/ref/HoldAllComplete.html) suppresses all three behaviors:

```wolfram
Attributes[HoldComplete]
(* Output *)
{HoldAllComplete,Protected}
```

```wolfram
HoldComplete[Sequence[a,b,c],g]
(* Output *)
HoldComplete[Sequence[a,b,c],g]
```

Substitution works inside a held expression:

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

[NHoldAll](https://reference.wolfram.com/language/ref/NHoldAll.html) protects arguments from [N](https://reference.wolfram.com/language/ref/N.html) but evaluates them normally otherwise:

```wolfram
Attributes[Derivative]
(* Output *)
{NHoldAll,ReadProtected}
```

```wolfram
N[Derivative[1+1]]
(* Output *)
Derivative[2]
```

```wolfram
N[f[1+1]]
(* Output *)
f[2.]
```

[HoldPattern](https://reference.wolfram.com/language/ref/HoldPattern.html) protects patterns from evaluation but does not interfere with pattern matching:

```wolfram
{HoldPattern[_+_],_+_}
(* Output *)
{HoldPattern[_+_],2 _}
```

## Tech Notes ▪Attributes ▪Non[Hyphen]Standard Evaluation

## Related Guides ▪Evaluation Control ▪Expressions ▪Attributes

[Graphics] | Related Workflows
▪Substitute Values of Variables in Functions That Hold Their Arguments

## History Introduced in 1988 (1.0)
