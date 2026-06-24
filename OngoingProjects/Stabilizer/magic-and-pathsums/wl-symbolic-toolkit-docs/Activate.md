# Activate | [SpanFromLeft]

> [Activate](https://reference.wolfram.com/language/ref/Activate.html)[*expr*] — replaces all instances of [Inactive](https://reference.wolfram.com/language/ref/Inactive.html)[*f*] in `*expr*` with `*f*`.
> [Activate](https://reference.wolfram.com/language/ref/Activate.html)[*expr*,*patt*] — replaces only instances of [Inactive](https://reference.wolfram.com/language/ref/Inactive.html)[*f*] for which `*f*` matches the pattern `*patt*`.

## Details and Options

With the option setting [Heads](https://reference.wolfram.com/language/ref/Heads.html)->[False](https://reference.wolfram.com/language/ref/False.html), [Activate](https://reference.wolfram.com/language/ref/Activate.html) does not enter heads of expressions and activate their parts.

## Examples

### Basic Examples

Activate an [Inactive](https://reference.wolfram.com/language/ref/Inactive.html) expression:

```wolfram
Activate[Length[{a,b,c}]]
(* Output *)
3
```

Activate different parts of an inactive expression:

```wolfram
expr=Inactivate[2+2+3^2]
(* Output *)
2+2+3^2
```

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

### Scope

Define an inactive expression:

```wolfram
expr=Inactive[Sin][π/2]
(* Output *)
Sin[(π)/(2)]
```

Evaluate the expression using [Activate](https://reference.wolfram.com/language/ref/Activate.html):

```wolfram
Activate[expr]
(* Output *)
1
```

Create an inactive expression using [Inactivate](https://reference.wolfram.com/language/ref/Inactivate.html):

```wolfram
expr=Inactivate[Sin[Pi/2],Sin]
(* Output *)
Sin[(π)/(2)]
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

Activate `h`:

```wolfram
Activate[%,h]
(* Output *)
f[g[h[x]],y,g[z]]
```

Activate an inactive expression:

```wolfram
expr=Cos[π]+x;
```

```wolfram
Activate[expr]
(* Output *)
-1+(x^2)/(2)
```

Activate all symbols except [Integrate](https://reference.wolfram.com/language/ref/Integrate.html):

```wolfram
Activate[expr,Except[Integrate]]
(* Output *)
-1+x
```

Prevent numeric functions from being activated:

```wolfram
Activate[expr,Except[_?(MemberQ[Attributes[#],NumericFunction]&)]]
(* Output *)
Cos[π]+(x^2)/(2)
```

Formally differentiating a Laplace transform:

```wolfram
Inactive[LaplaceTransform][a t^2,t,s]
(* Output *)
LaplaceTransform[a t^2,t,s]
```

```wolfram
D[%,s]
(* Output *)
LaplaceTransform[-a t^3,t,s]
```

```wolfram
Activate[%]
(* Output *)
-(6 a)/(s^4)
```

```wolfram
D[LaplaceTransform[a t^2,t,s],s]
(* Output *)
-(6 a)/(s^4)
```

Similarly differentiating wrt `t` and `a`:

```wolfram
D[Inactive[LaplaceTransform][a t^2,t,s],t]
(* Output *)
0
```

```wolfram
D[Inactive[LaplaceTransform][a t^2,t,s],a]
(* Output *)
LaplaceTransform[t^2,t,s]
```

Inactive special function expression:

```wolfram
expr=Inactivate[Hypergeometric2F1[3,1,2,x],Hypergeometric2F1]
(* Output *)
Hypergeometric2F1[3,1,2,x]
```

Expression after automatic simplification:

```wolfram
Activate[expr]
(* Output *)
(2-x)/(2 (1-x)^2)
```

```wolfram
Hypergeometric2F1[3,1,2,x]
(* Output *)
(2-x)/(2 (1-x)^2)
```

### Options

#### Heads

An inactive [Derivative](https://reference.wolfram.com/language/ref/Derivative.html) expression:

```wolfram
inactive=Inactivate[Derivative[1][Cos][x]]
(* Output *)
Derivative[1][Cos][x]
```

Activate the expression:

```wolfram
Activate[inactive]
(* Output *)
-Sin[x]
```

Use the option setting [Heads](https://reference.wolfram.com/language/ref/Heads.html)->[False](https://reference.wolfram.com/language/ref/False.html) to avoid activating [Derivative](https://reference.wolfram.com/language/ref/Derivative.html):

```wolfram
Activate[inactive, Heads->False]
(* Output *)
Derivative[1][Cos][x]
```

### Applications

Define a trigonometric expression with two inactive terms:

```wolfram
expr =Inactivate[Sin[π/3]+Cos[π/3]+Tan[π/3],Sin|Cos]
(* Output *)
Sqrt[3]+Cos[(π)/(3)]+Sin[(π)/(3)]
```

Activate different parts of the expression:

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
inactive=Inactivate[D[Integrate[(t+x)^2,{t,0,x}],x],D|Integrate]
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
di=Activate[%]
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
id=Activate[%]
(* Output *)
7 x^2
```

The results are mathematically the same:

```wolfram
Simplify[di-id]
(* Output *)
0
```

Solution for the three-dimensional Laplace equation in inactive integral form:

```wolfram
V[x_, y_, z_]=Inactivate[Integrate[f[z+I x Cos[u]+I y Sin[u],u],{u,-Pi,Pi}],
   Integrate]
(* Output *)
f[z+ⅈ x Cos[u]+ⅈ y Sin[u],u]
```

Obtain a particular solution by specifying the function `*f*`:

```wolfram
f[a_,b_]:=3a^5+7b^4
```

```wolfram
V[x,y,z]
(* Output *)
(7 u^4+3 (z+ⅈ x Cos[u]+ⅈ y Sin[u])^5)
```

```wolfram
sol=Activate[%]//Simplify
(* Output *)
(14 π^5)/(5)+(3)/(4) π z (15 (x^2+y^2)^2-40 (x^2+y^2) z^2+8 z^4)
```

Visualize the solution:

```wolfram
Row[Table[Plot3D[sol/.{z-> j},{x,-3,3},{y,-3,3},Ticks->{Automatic,Automatic,None}],{j,-2,2}]]
(* Output *)
![image](img/image_001.png)
```

Verify the solution:

```wolfram
Laplacian[sol,{x,y,z}]//Simplify
(* Output *)
0
```

Formula for summation by parts:

```wolfram
sumparts = Inactivate[Sum[f[k]g[k],k]==g[k]Sum[f[k],k]-
     Sum[DifferenceDelta[g[k],k]DiscreteShift[Sum[f[k],k],k],k],
   Sum|DifferenceDelta|DiscreteShift]
(* Output *)
f[k] g[k]==g[k] f[k]-Inactive Inactive
```

Verify the formula in a special case:

```wolfram
f[k_]:=k
```

```wolfram
g[k_]:=HarmonicNumber[k]
```

```wolfram
Activate[sumparts]
(* Output *)
True
```

Evaluate the sum:

```wolfram
Inactive[Sum][k HarmonicNumber[k],k]==Sum[k HarmonicNumber[k],k]
(* Output *)
k HarmonicNumber[k]==-(1)/(4) (-1+k) k+(1)/(2) (-1+k) k HarmonicNumber[k]
```

Explore vector identities:

```wolfram
divcurl = Inactivate[Div[Curl[{f[x,y,z],g[x,y,z], h[x,y,z]},{x,y,z}],
    {x,y,z}],Div|Curl]
(* Output *)
({f[x,y,z],g[x,y,z],h[x,y,z]})
```

Activating [Curl](https://reference.wolfram.com/language/ref/Curl.html) is not very interesting:

```wolfram
Activate[divcurl,Curl]
(* Output *)
{-g^(0,0,1)[x,y,z]+h^(0,1,0)[x,y,z],f^(0,0,1)[x,y,z]-h^(1,0,0)[x,y,z],-f^(0,1,0)[x,y,z]+g^(1,0,0)[x,y,z]}
```

Activating [Div](https://reference.wolfram.com/language/ref/Div.html) demonstrates the relation $(v)=0$:

```wolfram
Activate[divcurl,Div]
(* Output *)
0
```

### Properties & Relations

[Inactive](https://reference.wolfram.com/language/ref/Inactive.html) expressions can be evaluated using [Activate](https://reference.wolfram.com/language/ref/Activate.html):

```wolfram
f[x_]:=x^2
```

```wolfram
expr=Inactive[f][3]
(* Output *)
f[3]
```

```wolfram
Activate[expr]
(* Output *)
9
```

[Activate](https://reference.wolfram.com/language/ref/Activate.html) is the inverse of [Inactivate](https://reference.wolfram.com/language/ref/Inactivate.html):

```wolfram
Inactivate[f[x],f]
(* Output *)
f[x]
```

```wolfram
Activate[%]
(* Output *)
f[x]
```

[Activate](https://reference.wolfram.com/language/ref/Activate.html) replaces all instances of inactive symbols in an expression:

```wolfram
Inactivate[f[x]+g[x]+h[x],f|g]
(* Output *)
h[x]+f[x]+g[x]
```

```wolfram
Activate[%]
(* Output *)
f[x]+g[x]+h[x]
```

[Activate](https://reference.wolfram.com/language/ref/Activate.html) evaluates inactive expressions and allows parts of expressions to be inactive:

```wolfram
isin=Inactivate[Sin[ArcTan[1]]]
(* Output *)
Sin[ArcTan[1]]
```

```wolfram
Activate[isin,Sin]
(* Output *)
Sin[ArcTan[1]]
```

```wolfram
Activate[%]
(* Output *)
(1)/(Sqrt[2])
```

[ReleaseHold](https://reference.wolfram.com/language/ref/ReleaseHold.html) evaluates expressions held in unevaluated form, and all parts are evaluated:

```wolfram
esin=Hold[Sin[ArcTan[1]]]
(* Output *)
Hold[Sin[ArcTan[1]]]
```

```wolfram
ReleaseHold[%]
(* Output *)
(1)/(Sqrt[2])
```

### Neat Examples

Create a gallery of infinite products:

```wolfram
infiniteproducts=Inactivate[{∏_{k=1}^{n}((k+1)^3 (k+5))/(k^2),∏_{k=0}^{n}k!,∏_{k=1}^{∞}(1+(1)/(k^2)),∏_{k}^{∞}(1-(4)/(3) Sin[(x)/(3^k)]^2),∏_{k}^{∞}(1+(1)/(Prime[k]^s)),∏_{k=1}^{n}((k+3)/(k+1))^k,∏_{k=1}^{n}(Sin[3 k+5])/(Cos[3 k+1])},Product];
```

```wolfram
FormulaGallery[forms_List]:=Module[{vals=ParallelMap[Activate,forms]},Text[Grid[Table[{forms[[i]],"=",vals[[i]]},{i,Length[forms]}],Dividers->{{True,False,False,True},All},Alignment->{{Right,Center,Left},Baseline},Background->LightYellow,Spacings->{{4,{},4},1}]]];
```

```wolfram
FormulaGallery[infiniteproducts]//
(* Output *)
| ((k+1)^(3) (k+5))/(k^(2)) | "=" | (1)/(120) (n+1)^(3) n+1 n+6 |
| --- | --- | --- |
| k! | "=" | n+2 |
| ((1)/(k^(2))+1) | "=" | (sinh(π))/(π) |
| (1-(4)/(3) sin^(2)(3^(-k) x)) | "=" | sinc(x) |
| ((k)^(-s)+1) | "=" | (s)/(2 s) |
| ((k+3)/(k+1))^(k) | "=" | (2 n+2 exp(ζ^((1,0))(-1,n+4)-ζ^((1,0))(-1,n+2)))/(n+4^(3)) |
| sin(3 k+5) sec(3 k+1) | "=" | ((1+ℯ^(2 ⅈ)) ⅈ^(-n) ℯ^(4 ⅈ (n+2)) ℯ^(-10 ⅈ))/((-1+ℯ^(10 ⅈ)) -ℯ^(-2 ⅈ)) |
```

## Related Guides ▪Evaluation Control ▪Tuning & Debugging

[Graphics] | Related Workflows
▪Handle Code Symbolically

## History Introduced in 2014 (10.0)
