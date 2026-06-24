# NHoldAll | [SpanFromLeft]

> [NHoldAll](https://reference.wolfram.com/language/ref/NHoldAll.html) — is an attribute which specifies that none of the arguments to a function should be affected by [N](https://reference.wolfram.com/language/ref/N.html).

## Details

[NHoldAll](https://reference.wolfram.com/language/ref/NHoldAll.html), [NHoldFirst](https://reference.wolfram.com/language/ref/NHoldFirst.html), and [NHoldRest](https://reference.wolfram.com/language/ref/NHoldRest.html) are useful in ensuring that arguments to functions are maintained as exact integers, rather than being converted by [N](https://reference.wolfram.com/language/ref/N.html) to approximate numbers.

## Examples

### Basic Examples

Prevent [N](https://reference.wolfram.com/language/ref/N.html) from affecting the arguments of a function:

```wolfram
SetAttributes[f,NHoldAll]
```

```wolfram
N[f[1+2,3+Pi]]
(* Output *)
f[3,3+π]
```

### Scope

System symbols with the [NHoldAll](https://reference.wolfram.com/language/ref/NHoldAll.html) attribute:

```wolfram
ssymb=Cases[Map[ToExpression,Names["System`*"]],_Symbol];
Select[ssymb,MemberQ[Attributes[#],NHoldAll]&]
(* Output *)
{AlgebraicNumber,C,CompiledFunction,Derivative,InverseFunction,RecurringDigitsForm,Root,Slot,SlotSequence}
```

The arguments of [Derivative](https://reference.wolfram.com/language/ref/Derivative.html) remain unchanged with [N](https://reference.wolfram.com/language/ref/N.html):

```wolfram
der=Derivative[1,0,1,0][f][1,2,3,x]
(* Output *)
f^(1,0,1,0)[1,2,3,x]
```

[N](https://reference.wolfram.com/language/ref/N.html) leaves the derivative order while changing the point of evaluation:

```wolfram
nder=N[der]
(* Output *)
f^(1,0,1,0)[1.,2.,3.,x]
```

```wolfram
nder /. f->Times
(* Output *)
2. x
```

Define a pure function:

```wolfram
f=#1^2+Pi^3(#1^2-#2)^2&;
```

The function with coefficients converted to numerical values:

```wolfram
nf=N[f]
(* Output *)
#1^2+31.006276680299816 (#1^2-1. #2)^2&
```

```wolfram
nf[1,2]
(* Output *)
32.00627668029982
```

The positional parameters remain unchanged with [N](https://reference.wolfram.com/language/ref/N.html) because [Slot](https://reference.wolfram.com/language/ref/Slot.html) has the [NHoldAll](https://reference.wolfram.com/language/ref/NHoldAll.html) attribute:

```wolfram
FullForm[f]
(* Output *)
Function[Plus[Power[Slot[1],2],Times[Power[Pi,3],Power[Plus[Power[Slot[1],2],Times[-1,Slot[2]]],2]]]]
```

### Applications

Use an indexed variable:

```wolfram
SetAttributes[x,NHoldAll]
```

With this attribute, the variables remain unchanged:

```wolfram
N[2x[1]+5x[2]^2]
(* Output *)
2. x[1]+5. x[2]^2
```

Define a data object that represents a polynomial $c_{1} x^{p_{1}}+...$ in a sparse form $\{\{\mathit{c}_{1},\mathit{p}_{1}\},\ldots \}$:

```wolfram
spoly[cp_][x_]:=Module[{c, n,p, y},
{p,c}=Transpose[SortBy[cp,Last]];
n=p[[-1]];
y=c[[-1]];
Do[y=c[[k]]+x^(n-p[[k]])*y;n=p[[k]],{k,Length[p]-1,1,-1}];
If[n> 0,y*=x^n];
y]
```

Make sure that [N](https://reference.wolfram.com/language/ref/N.html) only affects the coefficients, not the powers:

```wolfram
N[e:spoly[cp_],pa_]:= Module[{p,c,nc},
{c,p}=Transpose[cp];
nc=N[c,pa];
If[nc===c,e,spoly[Transpose[{nc,p}]]]]
```

Default [N](https://reference.wolfram.com/language/ref/N.html) evaluation of the argument needs to be prevented for the rule above to work:

```wolfram
SetAttributes[spoly,NHoldAll]
```

A representation of the polynomial $x+2 x^{2}+3 x^{4}+4 x^{8}$:

```wolfram
sp=spoly[{{1,1},{2,2},{3,4},{4,8}}]
(* Output *)
spoly[{{1,1},{2,2},{3,4},{4,8}}]
```

```wolfram
{sp[x],Expand[sp[x]]}
(* Output *)
{x (1+x (2+x (4+8 x))),x+2 x^2+4 x^3+8 x^4}
```

Get the representation with approximate real coefficients:

```wolfram
nsp= N[sp]
(* Output *)
spoly[{{1.,1},{2.,2},{3.,4},{4.,8}}]
```

Evaluate at $x=\pi$:

```wolfram
nsp[π]
(* Output *)
926.1786364489872
```

### Properties & Relations

[HoldAll](https://reference.wolfram.com/language/ref/HoldAll.html) prevents evaluation while [NHoldAll](https://reference.wolfram.com/language/ref/NHoldAll.html) only prevents numerical evaluation:

```wolfram
SetAttributes[f1,HoldAll]
```

```wolfram
f1[1+2,Pi+E]
(* Output *)
f1[1+2,π+ℯ]
```

```wolfram
N[%]
(* Output *)
f1[1.+2.,3.141592653589793+2.718281828459045]
```

```wolfram
SetAttributes[f2,NHoldAll]
```

```wolfram
f2[1+2,Pi+E]
(* Output *)
f2[3,ℯ+π]
```

```wolfram
N[%]
(* Output *)
f2[3,ℯ+π]
```

You can prevent both by setting both attributes:

```wolfram
SetAttributes[f,{HoldAll,NHoldAll}]
```

```wolfram
f[1+2,Pi+E]
(* Output *)
f[1+2,π+ℯ]
```

```wolfram
N[%]
(* Output *)
f[1+2,π+ℯ]
```

## Tech Notes ▪Attributes ▪Controlling Numerical Evaluation

## Related Guides ▪Evaluation Control ▪Attributes

## History Introduced in 1996 (3.0)
