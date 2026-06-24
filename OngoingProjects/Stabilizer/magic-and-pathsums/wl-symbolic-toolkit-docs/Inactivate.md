# Inactivate | [SpanFromLeft]

> [Inactivate](https://reference.wolfram.com/language/ref/Inactivate.html)[*expr*] — replaces all instances of `*f*` with [Inactive](https://reference.wolfram.com/language/ref/Inactive.html)[*f*] for symbols `*f*` used as heads in `*expr*`.
> [Inactivate](https://reference.wolfram.com/language/ref/Inactivate.html)[*expr*,*patt*] — inactivates all symbols in `*expr*` that match the pattern `*patt*`.

## Details and Options

[Inactivate](https://reference.wolfram.com/language/ref/Inactivate.html) has attribute [HoldFirst](https://reference.wolfram.com/language/ref/HoldFirst.html), and symbols in `*expr*` are inactivated before evaluation.

With the option setting [Heads](https://reference.wolfram.com/language/ref/Heads.html)->[False](https://reference.wolfram.com/language/ref/False.html), [Inactivate](https://reference.wolfram.com/language/ref/Inactivate.html) does not enter heads of expressions and inactivate their parts.

By default, certain semantically important heads are not inactivated. Common examples include [List](https://reference.wolfram.com/language/ref/List.html),  [Rule](https://reference.wolfram.com/language/ref/Rule.html) and [Blank](https://reference.wolfram.com/language/ref/Blank.html).

## Examples

### Basic Examples

Inactivate an expression:

```wolfram
Inactivate[Length[{a,b,c}]]
(* Output *)
Length[{a,b,c}]
```

Activate it:

```wolfram
Activate[%]
(* Output *)
3
```

Inactivate an expression with several terms:

```wolfram
expr=Inactivate[2+2+3^2]
(* Output *)
2+2+3^2
```

Activate different parts of the expression:

```wolfram
Activate[expr,Plus]
(* Output *)
4+3^2
```

```wolfram
Activate[expr,Power]
(* Output *)
2+2+9
```

```wolfram
Activate[expr]
(* Output *)
13
```

Inactivate a symbol in an expression:

```wolfram
Inactivate[Sum[k^2,{k,1,n}],Sum]
(* Output *)
k^2
```

Evaluate the expression:

```wolfram
Activate[%]
(* Output *)
(1)/(6) n (1+n) (1+2 n)
```

### Scope

Inactivate an expression:

```wolfram
expr=Inactivate[Cos[Pi Sin[Pi]]]
(* Output *)
Cos[π*Sin[π]]
```

Evaluate the expression:

```wolfram
Activate[expr]
(* Output *)
1
```

Inactivate the symbol `g` only:

```wolfram
Inactivate[f[g[h[x]],y,g[z]],g]
(* Output *)
f[g[h[x]],y,g[z]]
```

Activate `g`:

```wolfram
Activate[%]
(* Output *)
f[g[h[x]],y,g[z]]
```

Inactivate `g` and `h`:

```wolfram
Inactivate[f[g[h[x]],y,g[z]],g|h]
(* Output *)
f[g[h[x]],y,g[z]]
```

Activate `h` only:

```wolfram
Activate[%,h]
(* Output *)
f[g[h[x]],y,g[z]]
```

Inactivate all symbols except [Integrate](https://reference.wolfram.com/language/ref/Integrate.html) in an expression:

```wolfram
Inactivate[Integrate[Sin[x]+Cos[x],x],Except[Integrate]]
(* Output *)
∫(Sin[x]+Cos[x])ⅆx
```

Evaluate the expression:

```wolfram
Activate[%]
(* Output *)
-Cos[x]+Sin[x]
```

Prevent numeric functions from being inactivated:

```wolfram
Inactivate[2+Integrate[Sin[x],x],Except[_?(MemberQ[Attributes[#],NumericFunction]&)]]
(* Output *)
2+Sin[x]
```

Normally [Plus](https://reference.wolfram.com/language/ref/Plus.html) and [Sin](https://reference.wolfram.com/language/ref/Sin.html) would also have been inactivated:

```wolfram
Inactivate[2+Integrate[Sin[x],x]]
(* Output *)
2+Sin[x]
```

### Options

#### Heads

By default, [Inactivate](https://reference.wolfram.com/language/ref/Inactivate.html) will make heads inactive, even inside compound heads:

```wolfram
Inactivate[f[x][y]+g[z]]
(* Output *)
f[x][y]+g[z]
```

With the setting [Heads](https://reference.wolfram.com/language/ref/Heads.html)->[False](https://reference.wolfram.com/language/ref/False.html), inactivation does not proceed inside compound heads:

```wolfram
Inactivate[f[x][y]+g[z],Heads->False]
(* Output *)
f[x][y]+g[z]
```

### Applications

Define a trigonometric expression:

```wolfram
expr =Inactivate[Sin[π/3]+Cos[π/3]+Tan[π/3],Sin|Cos]
(* Output *)
Sqrt[3]+Cos[(π)/(3)]+Sin[(π)/(3)]
```

Activate different trigonometric functions:

```wolfram
Activate[expr,Sin]
(* Output *)
(3 Sqrt[3])/(2)+Cos[(π)/(3)]
```

```wolfram
Activate[expr,Cos]
(* Output *)
(1)/(2)+Sqrt[3]+Sin[(π)/(3)]
```

```wolfram
Activate[expr]
(* Output *)
(1)/(2)+(3 Sqrt[3])/(2)
```

Define $\partial_{x}\int_{0}^{x}(t+x)^{2}\mathrm{d}t$, leaving both the derivative and integral inactive:

```wolfram
inactive=Inactivate[D[Integrate[(t+x)^2,{t,0,x}], x],D|Integrate]
(* Output *)
Inactive
```

Differentiate the integral without evaluating the integral:

```wolfram
Activate[inactive,D]
(* Output *)
4 x^2+2 (t+x)
```

Activate the integral to compute the final result:

```wolfram
di = Activate[%]
(* Output *)
7 x^2
```

Integrate without performing the differentiation:

```wolfram
Activate[inactive,Integrate]
(* Output *)
Inactive
```

Activate the differentiation to compute the final result:

```wolfram
id = Activate[%]
(* Output *)
7 x^2
```

The results are mathematically the same:

```wolfram
Simplify[di-id]
(* Output *)
0
```

Show identities including Leibniz's rule for differentiating integrals:

```wolfram
Inactivate[∂_x∫_a[x]^b[x]f[x]ⅆx,D|Integrate]==∂_x∫_a[x]^b[x]f[x]ⅆx //
(* Output *)
(∂f(x))/(∂x)==f(b(x)) b^(′)(x)-f(a(x)) a^(′)(x)
```

Chain rule:

```wolfram
Inactivate[D[f[g[x]],x],D]==D[f[g[x]],x] //
(* Output *)
(∂f(g(x)))/(∂x)==g^(′)(x) f^(′)(g(x))
```

Indefinite integrals:

```wolfram
Inactivate[D[Integrate[f[x],x],x]==f[x],Integrate|D]//
(* Output *)
(∂f(x))/(∂x)==f(x)
```

```wolfram
Inactivate[Integrate[x^n,x],Integrate]==Integrate[x^n,x]//
(* Output *)
x^(n)==(x^(n+1))/(n+1)
```

```wolfram
Inactivate[Integrate[x^(-1),x],Integrate]==Integrate[x^(-1),x]//
(* Output *)
(1)/(x)==log(x)
```

```wolfram
Inactivate[Integrate[1/Sqrt[1-x^4],x],Integrate]==Integrate[1/Sqrt[1-x^4],x]//
(* Output *)
(1)/(Sqrt(1-x^(4)))==sin^(-1)(x)
```

Infinite sums and products:

```wolfram
Inactivate[Sum[1/n^4,{n,1,Infinity}],Sum]==Sum[1/n^4,{n,1,Infinity}]//
(* Output *)
(1)/(n^(4))==(π^(4))/(90)
```

```wolfram
Inactivate[Product[1-1/n^2,{n,2,Infinity}],Product]==Product[1-1/n^2,{n,2,Infinity}]//
(* Output *)
(1-(1)/(n^(2)))==(1)/(2)
```

### Properties & Relations

[Activate](https://reference.wolfram.com/language/ref/Activate.html) is the inverse of [Inactivate](https://reference.wolfram.com/language/ref/Inactivate.html):

```wolfram
Inactivate[f[x]]
(* Output *)
f[x]
```

```wolfram
Activate[%]
(* Output *)
f[x]
```

[Inactivate](https://reference.wolfram.com/language/ref/Inactivate.html) replaces specific symbols with their inactive forms:

```wolfram
Inactivate[f[x]+g[x]+h[x],f|g] // FullForm
(* Output *)
Plus[h[x],Inactive[f][x],Inactive[g][x]]
```

[Activate](https://reference.wolfram.com/language/ref/Activate.html) replaces all instances of inactive symbols with their active forms:

```wolfram
Activate[%]
(* Output *)
f[x]+g[x]+h[x]
```

[Inactivate](https://reference.wolfram.com/language/ref/Inactivate.html) maintains symbols in inactive form and allows parts of expressions to be inactive:

```wolfram
isin=Inactivate[Sin[ArcTan[1]],Sin]
(* Output *)
Sin[(π)/(4)]
```

```wolfram
Activate[isin]
(* Output *)
(1)/(Sqrt[2])
```

[Hold](https://reference.wolfram.com/language/ref/Hold.html) maintains expressions in unevaluated form, and all parts are inactive:

```wolfram
esin=Hold[Sin[ArcTan[1]]]
(* Output *)
Hold[Sin[ArcTan[1]]]
```

```wolfram
ReleaseHold[esin]
(* Output *)
(1)/(Sqrt[2])
```

Compare an inactive expression with the corresponding [FullForm](https://reference.wolfram.com/language/ref/FullForm.html):

```wolfram
Inactivate[Cos[x Sin[y]]]
(* Output *)
Cos[x*Sin[y]]
```

```wolfram
%//FullForm
(* Output *)
Inactive[Cos][Inactive[Times][x,Inactive[Sin][y]]]
```

```wolfram
FullForm[Cos[x Sin[y]]]
(* Output *)
Cos[Times[x,Sin[y]]]
```

[Inactivate](https://reference.wolfram.com/language/ref/Inactivate.html) is an idempotent operator:

```wolfram
Inactivate[Inactive[Sin][x],Sin]//InputForm
(* Output *)
Inactive[Sin][x]
```

Certain heads are not inactivated by default, including [List](https://reference.wolfram.com/language/ref/List.html), [Rule](https://reference.wolfram.com/language/ref/Rule.html) (`->`) and [Blank](https://reference.wolfram.com/language/ref/Blank.html) (`_`):

```wolfram
Inactivate[{f[x], a->b,x_}]
(* Output *)
{f[x],a->b,x_}
```

Using [Replace](https://reference.wolfram.com/language/ref/Replace.html) at all levels can inactivate all heads in an expression:

```wolfram
Replace[{f[x], a->b,x_},h_[args___]:>Inactive[h][args],{0,-1}]
(* Output *)
List[f[x],a->b,Pattern[x,Blank[]]]
```

### Possible Issues

As compound heads have no attributes, the use of [Inactivate](https://reference.wolfram.com/language/ref/Inactivate.html) can lead to evaluation leaks:

```wolfram
k=1;
Inactivate[Sum[f[k],{k,1,n}],Sum]
(* Output *)
f[1]
```

Normally, the [HoldAll](https://reference.wolfram.com/language/ref/HoldAll.html) attribute of [Sum](https://reference.wolfram.com/language/ref/Sum.html) would prevent the `k` from evaluating:

```wolfram
Sum[f[k],{k,1,n}]
(* Output *)
∑_{k=1}^{n}f[k]
```

### Neat Examples

Create a gallery of definite integrals:

```wolfram
defint=Inactivate[{∫_0^∞(1)/(x^4+x^2+2)ⅆx,∫_0^(π)/(2)Cos[Sin[x]^2]ⅆx,∫_0^∞AiryAi[x]^2ⅆx,∫_0^1Ceiling[x^2+Abs[3 x-1]]ⅆx,∫_0^2Floor[x^2]ⅆx,∫_0^1(Log[(1)/(2) (1+Sqrt[1+4 x])])/(x)ⅆx,∫_0^∞(∏_{k=0}^{5}Sin[(x)/(2 k+1)])/(x^6)ⅆx},Product|Integrate];
```

```wolfram
FormulaGallery[forms_List]:=Module[{vals=ParallelMap[Activate,forms]},Text[Grid[Table[{forms[[i]],"=",vals[[i]]},{i,Length[forms]}],Dividers->{{True,False,False,True},All},Alignment->{{Right,Center,Left},Baseline},Background->LightYellow,Spacings->{{6,{},6},1}]]];
```

```wolfram
FormulaGallery[defint]//
(* Output *)
| (1)/(x^(4)+x^(2)+2) | "=" | Sqrt((1)/(14 Sqrt(2))-(1)/(56)) π |
| --- | --- | --- |
| cos(sin^(2)(x)) | "=" | (1)/(2) π cos((1)/(2)) 0 |
| x^(2) | "=" | (1)/(3^(2/3) (1)/(3)^(2)) |
| 3 x-1+x^(2) | "=" | (1)/(2) (12-Sqrt(17)-Sqrt(21)) |
| x^(2) | "=" | 5-Sqrt(2)-Sqrt(3) |
| (log((1)/(2) (Sqrt(4 x+1)+1)))/(x) | "=" | (π^(2))/(15) |
| (sin((x)/(2 k+1)))/(x^(6)) | "=" | (π)/(20790) |
```

## Related Guides ▪Evaluation Control ▪Mathematical Typesetting ▪Tuning & Debugging

[Graphics] | Related Workflows
▪Handle Code Symbolically

## History Introduced in 2014 (10.0)
